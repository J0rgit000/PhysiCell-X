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

namespace
{

int mechanics_step_count( void )
{
	return static_cast<int>( std::llround( PhysiCell_settings.max_time / diffusion_dt ) );
}

void run_validation( const std::string& snapshot_path )
{
	PhysiCell_globals.current_time = 0.0;
	initialize_snapshot_csv( snapshot_path );
	append_serial_snapshot( snapshot_path, 0, PhysiCell_globals.current_time );

	Cell_Container* container = static_cast<Cell_Container*>( microenvironment.agent_container );
	const int total_steps = mechanics_step_count();
	for( int iteration = 1; iteration <= total_steps; ++iteration )
	{
		container->update_all_cells( PhysiCell_globals.current_time );
		PhysiCell_globals.current_time += diffusion_dt;
		append_serial_snapshot( snapshot_path, iteration, PhysiCell_globals.current_time );
	}
}

}

int main( int argc, char* argv[] )
{
	try
	{
		const std::string config_path = argc > 1 ? argv[1] : "./config/validation.xml";
		const std::string snapshot_path = argc > 2 ? argv[2] : "reference_snapshots.csv";

		if( !load_PhysiCell_config_file( config_path ) )
		{
			std::cerr << "Failed to load config: " << config_path << std::endl;
			return 1;
		}

		omp_set_num_threads( PhysiCell_settings.omp_num_threads );

		setup_microenvironment();
		const double mechanics_voxel_size = 20.0;
		create_cell_container_for_microenvironment( microenvironment, mechanics_voxel_size );
		create_cell_types();
		setup_tissue();

		run_validation( snapshot_path );
	}
	catch( const std::exception& error )
	{
		std::cerr << "Mechanics serial reference failed: " << error.what() << std::endl;
		return 1;
	}

	return 0;
}
