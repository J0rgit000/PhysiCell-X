#include <cmath>
#include <exception>
#include <iostream>
#include <string>

#include <omp.h>

#include "../../../modules/PhysiCell_settings.h"
#include "./custom_modules/custom.h"
#include "./custom_modules/snapshot_writer.h"

using namespace BioFVM;
using namespace PhysiCell;
using namespace DistPhy::mpi;

namespace
{

int mechanics_step_count( void )
{
	return static_cast<int>( std::llround( PhysiCell_settings.max_time / diffusion_dt ) );
}

void run_validation( const std::string& snapshot_path, mpi_Environment& world, mpi_Cartesian& cart_topo )
{
	PhysiCell_globals.current_time = 0.0;
	if( world.rank == 0 )
	{
		initialize_snapshot_csv( snapshot_path );
	}
	MPI_Barrier( cart_topo.mpi_cart_comm );
	append_mpi_snapshot( snapshot_path, 0, PhysiCell_globals.current_time, world, cart_topo );

	Cell_Container* container = static_cast<Cell_Container*>( microenvironment.agent_container );
	const int total_steps = mechanics_step_count();
	for( int iteration = 1; iteration <= total_steps; ++iteration )
	{
		container->update_all_cells( PhysiCell_globals.current_time, world, cart_topo );
		PhysiCell_globals.current_time += diffusion_dt;
		append_mpi_snapshot( snapshot_path, iteration, PhysiCell_globals.current_time, world, cart_topo );
	}
}

}

int main( int argc, char* argv[] )
{
	mpi_Environment world;
	world.Initialize();

	mpi_Cartesian cart_topo;
	cart_topo.Build_Cartesian_Topology( world );
	cart_topo.Find_Cartesian_Coordinates( world );
	cart_topo.Find_Left_Right_Neighbours( world );

	try
	{
		const std::string config_path = argc > 1 ? argv[1] : "./config/validation.xml";
		const std::string snapshot_path = argc > 2 ? argv[2] : "mpi_snapshots.csv";

		if( !load_PhysiCell_config_file( config_path, world ) )
		{
			if( world.rank == 0 )
			{
				std::cerr << "Failed to load config: " << config_path << std::endl;
			}
			world.Finalize();
			return 1;
		}

		omp_set_num_threads( PhysiCell_settings.omp_num_threads );

		setup_microenvironment( world, cart_topo );
		const double mechanics_voxel_size = 20.0;
		create_cell_container_for_microenvironment( microenvironment, mechanics_voxel_size, world, cart_topo );
		create_cell_types( world, cart_topo );
		setup_tissue( microenvironment, world, cart_topo );

		run_validation( snapshot_path, world, cart_topo );
	}
	catch( const std::exception& error )
	{
		if( world.rank == 0 )
		{
			std::cerr << "Mechanics MPI validation failed: " << error.what() << std::endl;
		}
		world.Finalize();
		return 1;
	}

	world.Finalize();
	return 0;
}
