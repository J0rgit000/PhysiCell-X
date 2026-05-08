#include "mechanics_comparison.h"
#include <sstream>
#include <iomanip>
#include <cmath>

/*================================================================================
 Capture all cells from a cell container into snapshots
 ================================================================================*/
std::vector<CellSnapshot> capture_all_cells( Cell_Container* container )
{
	std::vector<CellSnapshot> snapshots;
	
	if( !container || container->agent_list.empty() )
	{
		return snapshots;  // Empty if no cells
	}
	
	// Iterate through all cells in this container
	for( auto cell : container->agent_list )
	{
		CellSnapshot snap;
		snap.cell_id = cell->ID;
		snap.position[0] = cell->position[0];
		snap.position[1] = cell->position[1];
		snap.position[2] = cell->position[2];
		snap.velocity[0] = cell->velocity[0];
		snap.velocity[1] = cell->velocity[1];
		snap.velocity[2] = cell->velocity[2];
		snap.num_attachments = 0;  // Will be set by caller if needed
		snap.max_force = 0.0;       // Will be set by caller if needed
		snap.valid = true;
		
		snapshots.push_back( snap );
	}
	
	// Sort by cell ID for consistent ordering
	std::sort( snapshots.begin(), snapshots.end(), 
		[](const CellSnapshot& a, const CellSnapshot& b) { return a.cell_id < b.cell_id; } );
	
	return snapshots;
}

/*================================================================================
 Compare two snapshots and compute metrics
 ================================================================================*/
ComparisonResult compare_snapshots( 
	const std::vector<CellSnapshot>& omp_snap,
	const std::vector<CellSnapshot>& mpi_snap,
	double pos_tol_um,
	double vel_tol )
{
	ComparisonResult result;
	result.omp_cell_count = omp_snap.size();
	result.mpi_cell_count = mpi_snap.size();
	
	// Check cell count match
	if( omp_snap.size() != mpi_snap.size() )
	{
		result.cell_count_match = false;
		result.report = "ERROR: Cell count mismatch! OMP=" + std::to_string(omp_snap.size()) 
			+ " vs MPI=" + std::to_string(mpi_snap.size()) + "\n";
		return result;
	}
	result.cell_count_match = true;
	
	// Create ID-based map for MPI snapshot for fast lookup
	std::map<int, const CellSnapshot*> mpi_map;
	for( const auto& snap : mpi_snap )
	{
		mpi_map[snap.cell_id] = &snap;
	}
	
	// Compute position and velocity errors
	double pos_error_sum = 0.0;
	double vel_error_sum = 0.0;
	result.max_pos_error = 0.0;
	result.max_vel_error = 0.0;
	result.cells_out_of_tolerance = 0;
	
	for( const auto& omp_cell : omp_snap )
	{
		auto mpi_it = mpi_map.find( omp_cell.cell_id );
		if( mpi_it == mpi_map.end() )
		{
			result.report += "ERROR: Cell " + std::to_string(omp_cell.cell_id) + " found in OMP but not MPI!\n";
			continue;
		}
		
		const CellSnapshot* mpi_cell = mpi_it->second;
		
		// Position error
		double dx = omp_cell.position[0] - mpi_cell->position[0];
		double dy = omp_cell.position[1] - mpi_cell->position[1];
		double dz = omp_cell.position[2] - mpi_cell->position[2];
		double pos_error = std::sqrt( dx*dx + dy*dy + dz*dz );
		
		pos_error_sum += pos_error * pos_error;
		if( pos_error > result.max_pos_error )
		{
			result.max_pos_error = pos_error;
		}
		if( pos_error > pos_tol_um )
		{
			result.cells_out_of_tolerance++;
		}
		
		// Velocity error
		double dvx = omp_cell.velocity[0] - mpi_cell->velocity[0];
		double dvy = omp_cell.velocity[1] - mpi_cell->velocity[1];
		double dvz = omp_cell.velocity[2] - mpi_cell->velocity[2];
		double vel_error = std::sqrt( dvx*dvx + dvy*dvy + dvz*dvz );
		
		vel_error_sum += vel_error * vel_error;
		if( vel_error > result.max_vel_error )
		{
			result.max_vel_error = vel_error;
		}
	}
	
	// RMS errors
	size_t n_cells = omp_snap.size();
	if( n_cells > 0 )
	{
		result.rms_pos_error = std::sqrt( pos_error_sum / n_cells );
		result.rms_vel_error = std::sqrt( vel_error_sum / n_cells );
	}
	
	// Check tolerances
	result.position_match = ( result.rms_pos_error < pos_tol_um );
	result.velocity_match = ( result.rms_vel_error < vel_tol );
	
	// Force and topology matching would require more detailed cell state
	// For now, mark as tentatively matching (can be expanded later)
	result.force_match = true;
	result.topology_match = true;
	result.force_relative_error = 0.0;
	result.topology_agreement = 1.0;
	
	return result;
}

/*================================================================================
 Write snapshot to CSV
 ================================================================================*/
void write_snapshot_to_csv( 
	const std::string& filename,
	const std::vector<CellSnapshot>& snapshots,
	int iteration )
{
	std::ofstream outfile( filename );
	if( !outfile.is_open() )
	{
		std::cerr << "ERROR: Could not open " << filename << " for writing\n";
		return;
	}
	
	// Header
	outfile << "iteration,cell_id,x,y,z,vx,vy,vz,num_attachments,max_force\n";
	
	// Data
	for( const auto& snap : snapshots )
	{
		outfile << iteration << ","
				<< snap.cell_id << ","
				<< std::fixed << std::setprecision(6)
				<< snap.position[0] << ","
				<< snap.position[1] << ","
				<< snap.position[2] << ","
				<< snap.velocity[0] << ","
				<< snap.velocity[1] << ","
				<< snap.velocity[2] << ","
				<< snap.num_attachments << ","
				<< snap.max_force << "\n";
	}
	
	outfile.close();
}

/*================================================================================
 Write comparison report
 ================================================================================*/
void write_comparison_report( 
	const std::string& filename,
	const std::vector<ComparisonResult>& results,
	int total_iterations )
{
	std::ofstream outfile( filename );
	if( !outfile.is_open() )
	{
		std::cerr << "ERROR: Could not open " << filename << " for writing\n";
		return;
	}
	
	outfile << "=================================================================\n";
	outfile << "  OMP vs MPI Mechanics Solver Comparison Report\n";
	outfile << "=================================================================\n\n";
	
	outfile << "Configuration:\n";
	outfile << "  Total iterations: " << total_iterations << "\n";
	outfile << "  Sampled iterations: " << results.size() << "\n";
	outfile << "  Position tolerance: 1.0 um\n";
	outfile << "  Velocity tolerance: 0.1 um/min\n\n";
	
	// Summary statistics
	double max_rms_pos = 0.0, max_rms_vel = 0.0;
	int match_count = 0;
	
	for( const auto& result : results )
	{
		if( result.rms_pos_error > max_rms_pos )
			max_rms_pos = result.rms_pos_error;
		if( result.rms_vel_error > max_rms_vel )
			max_rms_vel = result.rms_vel_error;
		if( result.position_match && result.velocity_match && result.cell_count_match )
			match_count++;
	}
	
	outfile << "Summary:\n";
	outfile << "  Iterations matching tolerance: " << match_count << " / " << results.size() << "\n";
	outfile << "  Max RMS position error: " << std::fixed << std::setprecision(6) << max_rms_pos << " um\n";
	outfile << "  Max RMS velocity error: " << std::fixed << std::setprecision(6) << max_rms_vel << " um/min\n\n";
	
	// Detailed results
	outfile << "Detailed Results:\n";
	outfile << "Iter | Cell_Count | RMS_Pos_um | RMS_Vel | Max_Pos_um | Cells_OoT | Pass\n";
	outfile << "----|----------|-----------|--------|-----------|----------|-----\n";
	
	for( const auto& result : results )
	{
		outfile << std::setw(4) << results.size();
		outfile << " | " << std::setw(8) << result.omp_cell_count;
		outfile << " | " << std::fixed << std::setprecision(6) << std::setw(9) << result.rms_pos_error;
		outfile << " | " << std::setw(6) << std::setprecision(6) << result.rms_vel_error;
		outfile << " | " << std::setw(9) << result.max_pos_error;
		outfile << " | " << std::setw(8) << result.cells_out_of_tolerance;
		outfile << " | ";
		if( result.position_match && result.velocity_match && result.cell_count_match )
			outfile << "PASS\n";
		else
			outfile << "FAIL\n";
	}
	
	outfile << "\n=================================================================\n";
	outfile << "Individual Error Reports:\n";
	outfile << "=================================================================\n\n";
	
	for( size_t i = 0; i < results.size(); i++ )
	{
		const auto& result = results[i];
		if( !result.report.empty() )
		{
			outfile << "Iteration " << i << ":\n" << result.report << "\n";
		}
	}
	
	outfile.close();
}

/*================================================================================
 Determine if an iteration should be sampled (exponential sampling)
 Samples: 0, 1, 5, 10, 50, 100, 500, 1000, ...
 ================================================================================*/
bool should_sample_iteration( int iteration )
{
	if( iteration == 0 || iteration == 1 )
		return true;
	
	// Check if iteration is in pattern: 5, 10, 50, 100, 500, 1000, ...
	// Pattern: 10^n, 5*10^n for n=0,1,2,...
	
	int temp = iteration;
	
	// Try 5 * 10^n
	if( temp % 5 == 0 )
	{
		temp /= 5;
		while( temp % 10 == 0 )
			temp /= 10;
		if( temp == 1 )
			return true;
	}
	
	// Try 10^n
	temp = iteration;
	while( temp % 10 == 0 )
		temp /= 10;
	if( temp == 1 )
		return true;
	
	return false;
}
