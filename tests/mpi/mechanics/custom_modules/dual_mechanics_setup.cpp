#include "dual_mechanics_setup.h"
#include "heterogeneity.h"
#include "mechanics_comparison.h"
#include <iostream>
#include <vector>

/*================================================================================
 Global cell container and microenvironment pointers
 ================================================================================*/

Cell_Container* cell_container_omp = NULL;
Cell_Container* cell_container_mpi = NULL;
BioFVM::Microenvironment* microenv_omp = NULL;
BioFVM::Microenvironment* microenv_mpi = NULL;

extern BioFVM::Microenvironment microenvironment;  // Defined elsewhere

/*================================================================================
 Setup OMP-only microenvironment (serial, no MPI calls)
 ================================================================================*/

void setup_microenvironment_omp( void )
{
	// Create a new independent microenvironment for OMP testing
	microenv_omp = new BioFVM::Microenvironment;
	
	// Copy domain parameters from global microenvironment (which was set up by MPI version)
	// Domain bounds
	double x_min = microenvironment.mesh.bounding_box[0];
	double x_max = microenvironment.mesh.bounding_box[3];
	double y_min = microenvironment.mesh.bounding_box[1];
	double y_max = microenvironment.mesh.bounding_box[4];
	double z_min = microenvironment.mesh.bounding_box[2];
	double z_max = microenvironment.mesh.bounding_box[5];
	
	double dx = microenvironment.mesh.dx;
	double dy = microenvironment.mesh.dy;
	double dz = microenvironment.mesh.dz;
	
	// Initialize mesh
	microenv_omp->resize_uniform(x_min, x_max, y_min, y_max, z_min, z_max, dx, dy, dz);
	
	// Copy substrate definitions from MPI microenvironment
	microenv_omp->substrates = microenvironment.substrates;
	microenv_omp->substrate_diffusion_rates = microenvironment.substrate_diffusion_rates;
	microenv_omp->substrate_decay_rates = microenvironment.substrate_decay_rates;
	microenv_omp->bulk_supply_rate_function = microenvironment.bulk_supply_rate_function;
	
	std::cout << "[OMP] Initialized OMP-only microenvironment with " 
			  << microenv_omp->number_of_voxels() << " voxels\n";
}

/*================================================================================
 Setup OMP-only cell container
 ================================================================================*/

void setup_cell_container_omp( double mechanics_voxel_size )
{
	if( microenv_omp == NULL )
	{
		std::cerr << "ERROR: OMP microenvironment not initialized!\n";
		return;
	}
	
	// Create OMP cell container using the non-MPI signature (2 parameters only)
	cell_container_omp = create_cell_container_for_microenvironment( 
		*microenv_omp, 
		mechanics_voxel_size 
	);
	
	std::cout << "[OMP] Created OMP-only cell container\n";
}

/*================================================================================
 Distribute cells identically to both OMP and MPI containers
 
 Strategy:
 - Generate cell positions once (same as MPI version)
 - Add each cell to BOTH containers
 - This ensures identical initial state for comparison
 ================================================================================*/

void distribute_cells_to_both_containers( 
	BioFVM::Microenvironment* env_omp,
	BioFVM::Microenvironment* env_mpi,
	Cell_Container* container_omp,
	Cell_Container* container_mpi,
	mpi_Environment &world,
	mpi_Cartesian &cart_topo )
{
	std::string msg = "[DUAL] Distributing cells to both OMP and MPI containers...\n";
	if( world.rank == 0 ) std::cout << msg;
	
	// Generate identical cell positions
	double cell_radius = cell_defaults.phenotype.geometry.radius;
	double tumor_radius = parameters.doubles("tumor_radius");
	
	std::vector<std::vector<double>> positions;
	
	if( world.rank == 0 )
	{
		// Only rank 0 generates positions (same as MPI version)
		positions = create_cell_sphere_positions(cell_radius, tumor_radius);
		std::cout << "[DUAL] Generated " << positions.size() << " cell positions\n";
	}
	
	// Broadcast number of cells to all ranks
	int num_cells = positions.size();
	MPI_Bcast(&num_cells, 1, MPI_INT, 0, MPI_COMM_WORLD);
	
	// Broadcast all positions to all ranks
	if( num_cells > 0 )
	{
		std::vector<double> flat_positions;
		if( world.rank == 0 )
		{
			flat_positions.resize(num_cells * 3);
			for( int i = 0; i < num_cells; i++ )
			{
				flat_positions[3*i] = positions[i][0];
				flat_positions[3*i+1] = positions[i][1];
				flat_positions[3*i+2] = positions[i][2];
			}
		}
		else
		{
			flat_positions.resize(num_cells * 3);
		}
		
		MPI_Bcast(flat_positions.data(), num_cells * 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		
		if( world.rank != 0 )
		{
			positions.resize(num_cells);
			for( int i = 0; i < num_cells; i++ )
			{
				positions[i].resize(3);
				positions[i][0] = flat_positions[3*i];
				positions[i][1] = flat_positions[3*i+1];
				positions[i][2] = flat_positions[3*i+2];
			}
		}
	}
	
	// Get starting cell ID
	int start_cell_ID = (world.rank == 0) ? Basic_Agent::get_max_ID_in_parallel() : 0;
	MPI_Bcast(&start_cell_ID, 1, MPI_INT, 0, MPI_COMM_WORLD);
	
	// Add cells to BOTH containers (all ranks add all cells for OMP, MPI distributes via domain decomposition)
	int local_cell_count = 0;
	
	double p_mean = parameters.doubles("oncoprotein_mean");
	double p_sd = parameters.doubles("oncoprotein_sd");
	double p_min = parameters.doubles("oncoprotein_min");
	double p_max = parameters.doubles("oncoprotein_max");
	
	for( int i = 0; i < num_cells; i++ )
	{
		double x = positions[i][0];
		double y = positions[i][1];
		double z = positions[i][2];
		
		// Create cell
		Cell* pCell = create_cell(start_cell_ID + i);
		
		// Set position
		pCell->assign_position(x, y, z);
		
		// Set custom properties (heterogeneous)
		pCell->custom_data[0] = NormalRandom(p_mean, p_sd);
		if( pCell->custom_data[0] < p_min ) pCell->custom_data[0] = p_min;
		if( pCell->custom_data[0] > p_max ) pCell->custom_data[0] = p_max;
		
		// For OMP: add all cells to OMP container
		container_omp->add_cell_to_voxel( pCell, env_omp );
		local_cell_count++;
		
		// For MPI: check if cell is in this rank's domain
		// (MPI will handle distribution via domain decomposition)
		int voxel_index = env_mpi->get_closest_voxel_index(x, y, z);
		if( voxel_index >= 0 && voxel_index < (int)env_mpi->number_of_voxels() )
		{
			// Cell is in this rank's domain - it was already added by the MPI setup_tissue
			// Don't double-add; the MPI version already created this cell via create_cell()
		}
	}
	
	if( world.rank == 0 )
	{
		std::cout << "[DUAL] OMP: Added " << local_cell_count << " cells to OMP container\n";
	}
	
	// Update global cell ID counter
	Basic_Agent::set_max_ID_in_parallel(start_cell_ID + num_cells);
}

/*================================================================================
 Run one mechanics update iteration for both OMP and MPI
 ================================================================================*/

ComparisonResult run_dual_iteration( 
	double dt,
	mpi_Environment &world,
	mpi_Cartesian &cart_topo )
{
	// Capture state BEFORE update
	auto omp_snap_before = capture_all_cells(cell_container_omp);
	auto mpi_snap_before = capture_all_cells(cell_container_mpi);
	
	// Run OMP mechanics solver
	((Cell_Container *)microenv_omp->agent_container)->update_all_cells(dt);
	
	// Run MPI mechanics solver  
	((Cell_Container *)microenv_mpi->agent_container)->update_all_cells(dt, world, cart_topo);
	
	// Capture state AFTER update
	auto omp_snap_after = capture_all_cells(cell_container_omp);
	auto mpi_snap_after = capture_all_cells(cell_container_mpi);
	
	// Compare snapshots (use AFTER update snapshots)
	ComparisonResult result = compare_snapshots(omp_snap_after, mpi_snap_after, 1.0, 0.1);
	
	return result;
}

/*================================================================================
 Finalize comparison test and write reports
 ================================================================================*/

void finalize_comparison_test( 
	const std::vector<ComparisonResult>& all_results,
	int total_iterations,
	mpi_Environment &world )
{
	if( world.rank != 0 )
		return;  // Only rank 0 writes report
	
	// Create output directory if needed
	char output_path[1024];
	sprintf(output_path, "%s/comparison", PhysiCell_settings.folder.c_str());
	
	// Write comparison report
	char report_file[1024];
	sprintf(report_file, "%s/comparison_report.txt", output_path);
	write_comparison_report(report_file, all_results, total_iterations);
	
	std::cout << "\n[COMPARISON] Report written to: " << report_file << "\n";
	
	// Summary to console
	std::cout << "\n=================================================================\n";
	std::cout << "OMP vs MPI Mechanics Solver Comparison Summary\n";
	std::cout << "=================================================================\n";
	
	int pass_count = 0;
	for(const auto& result : all_results )
	{
		if( result.position_match && result.velocity_match )
			pass_count++;
	}
	
	std::cout << "Iterations passing tolerance (out of " << all_results.size() << "): " 
			  << pass_count << "\n";
	std::cout << "=================================================================\n\n";
}

/*================================================================================
 Helper: Check if iteration should be sampled (used in main loop)
 ================================================================================*/

bool should_sample_dual_iteration( int iteration )
{
	return should_sample_iteration(iteration);
}
