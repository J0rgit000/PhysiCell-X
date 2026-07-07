/*
###############################################################################
# If you use PhysiCell in your project, please cite PhysiCell and the version #
# number, such as below:                                                      #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1].    #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# See VERSION.txt or call get_PhysiCell_version() to get the current version  #
#     x.y.z. Call display_citations() to get detailed information on all cite-#
#     able software used in your PhysiCell application.                       #
#                                                                             #
# Because PhysiCell extensively uses BioFVM, we suggest you also cite BioFVM  #
#     as below:                                                               #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1],    #
# with BioFVM [2] to solve the transport equations.                           #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# [2] A Ghaffarizadeh, SH Friedman, and P Macklin, BioFVM: an efficient para- #
#     llelized diffusive transport solver for 3-D biological simulations,     #
#     Bioinformatics 32(8): 1256-8, 2016. DOI: 10.1093/bioinformatics/btv730  #
#                                                                             #
###############################################################################
#                                                                             #
# BSD 3-Clause License (see https://opensource.org/licenses/BSD-3-Clause)     #
#                                                                             #
# Copyright (c) 2015-2018, Paul Macklin and the PhysiCell Project             #
# All rights reserved.                                                        #
#                                                                             #
# Redistribution and use in source and binary forms, with or without          #
# modification, are permitted provided that the following conditions are met: #
#                                                                             #
# 1. Redistributions of source code must retain the above copyright notice,   #
# this list of conditions and the following disclaimer.                       #
#                                                                             #
# 2. Redistributions in binary form must reproduce the above copyright        #
# notice, this list of conditions and the following disclaimer in the         #
# documentation and/or other materials provided with the distribution.        #
#                                                                             #
# 3. Neither the name of the copyright holder nor the names of its            #
# contributors may be used to endorse or promote products derived from this   #
# software without specific prior written permission.                         #
#                                                                             #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" #
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE   #
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE  #
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE   #
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR         #
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF        #
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS    #
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN     #
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)     #
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE  #
# POSSIBILITY OF SUCH DAMAGE.                                                 #
#                                                                             #
###############################################################################
*/

/*================================================================================
+ If you use PhysiCell-X in your project, we would really appreciate if you can  +
+																				 +
+ [1] Cite the PhysiCell-X repository by giving its URL							 +
+																				 +
+ [2] Cite BioFVM-X: 															 +
+		Saxena, Gaurav, Miguel Ponce-de-Leon, Arnau Montagud, David Vicente Dorca,   +
+		and Alfonso Valencia. "BioFVM-X: An MPI+ OpenMP 3-D Simulator for Biological + 
+		Systems." In International Conference on Computational Methods in Systems    +
+		Biology, pp. 266-279. Springer, Cham, 2021. 							 +
=================================================================================*/

/*=======================================*/
/* Include mpi.h in the parallel version */
/*=======================================*/
/**
 *	\main main-prostate file
 *	\brief Main file of the prostate example

 *	\date 20/08/2025
 *	\author Arnau Montagud, Annika Meert, Gerard Pradas and Miguel Ponce de Leon, BSC-CNS
 */
#include <mpi.h>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <ctime>
#include <cmath>
#include <omp.h>
#include <fstream>
#ifdef EXTRAE_ENABLED
	#include "extrae.h"
#endif

#include "./core/PhysiCell.h"
#include "./modules/PhysiCell_standard_modules.h" 
#include "./addons/PhysiBoSS/src/maboss_intracellular.h"	

// put custom code modules here!

#include "./custom_modules/custom.h" 
#include "./custom_modules/drug_sensitivity.h"
#include "./custom_modules/boolean_model_interface.h"

using namespace BioFVM;
using namespace PhysiCell;

/*====================================================*/
/* Use DistPhy::mpi namespace in the parallel version */
/*====================================================*/
using namespace DistPhy::mpi;

int main( int argc, char* argv[] )
{
	/*=======================================================================================*/
	/* Create mpi_Environment object, initialize it, then create Cartesian Topology          */
	/*=======================================================================================*/
		
	mpi_Environment world;                         //Object contains size of communicator, rank of process
	world.Initialize();                            //Initialize using MPI_THREAD_FUNNELED, find comm. size and comm. rank
	mpi_Cartesian cart_topo;                       //Contains dims[3], ccoords[3] array and MPI_Comm mpi_cart_comm
	cart_topo.Build_Cartesian_Topology(world);     //Create 1-D X decomposition by setting dims[1]=size. 
	cart_topo.Find_Cartesian_Coordinates(world);   //Find Cartesian Topology coordinates of each process
	cart_topo.Find_Left_Right_Neighbours(world); 	 //Finds left/right immediate neighbour processes ranks and stores in X_LEFT/X_RIGHT

	// load and parse settings file(s)
	
	bool XML_status = false; 
	char copy_command [1024]; 
	char copy_command_2 [1024];

	#ifdef EXTRAE_ENABLED
    	Extrae_shutdown();
    #endif	

	if( argc > 1 )
	{ 
		XML_status = load_PhysiCell_config_file( argv[1], world ); 
	}
	else
	{ 
		XML_status = load_PhysiCell_config_file( "./config/PhysiCell_settings.xml", world ); 
	}
	
	if( !XML_status )
	{ 
		exit(-1); 
	}

	// copy config file to output directry 
	system( copy_command ); 
	
	if ( argc > 1 )
	{
		system( copy_command_2 );
	}
	
	// OpenMP setup
	//omp_set_num_threads(PhysiCell_settings.omp_num_threads);
	
	// PNRG setup 
	SeedRandom( parameters.ints("random_seed") ); // or specify a seed here
	std::string time_units = "min"; 

	/* Microenvironment setup */ 
	
	// store each column of drug sensitivity csv file in a custom vector variable
	// std::vector<std::pair<std::string, std::vector<double>>> csv_file = read_csv( parameters.strings("drug_sensitivity_file") );
	// for (std::pair<std::string, std::vector<double>> &element : csv_file ) {
	// 	if (element.first != "\"\""){
	// 		cell_defaults.custom_data.add_vector_variable( element.first, element.second );
	// 	}
	// }

	/* Microenvironment setup */ 
	setup_microenvironment(world, cart_topo);  

	
	/* PhysiCell setup */ 
 	
	// set mechanics voxel size, and match the data structure to BioFVM
	double mechanics_voxel_size = 20; 
	Cell_Container* cell_container = create_cell_container_for_microenvironment( microenvironment, mechanics_voxel_size, world, cart_topo );
	
	/* Users typically start modifying here. START USERMODS */ 
	
	create_cell_types( world, cart_topo );
	
	setup_tissue(microenvironment, world, cart_topo);

	/* Users typically stop modifying here. END USERMODS */ 
	
	// set MultiCellDS save options 

	set_save_biofvm_mesh_as_matlab( true ); 
	set_save_biofvm_data_as_matlab( true ); 
	set_save_biofvm_cell_data( true ); 
	set_save_biofvm_cell_data_as_custom_matlab( true );
	
	// save a simulation snapshot 
	
	char filename[1024];
	sprintf( filename , "%s/initial" , PhysiCell_settings.folder.c_str() ); 
	SVG_plot_mpi( filename , microenvironment, 0.0 , PhysiCell_globals.current_time, cell_coloring_function, substrate_coloring_function, world, cart_topo ); 
	// save_PhysiCell_to_MultiCellDS_xml_pugi( filename , microenvironment , PhysiCell_globals.current_time ); 


	sprintf( filename , "%s/states_initial.csv", PhysiCell_settings.folder.c_str());
	MaBoSSIntracellular::save( filename); // JM: No usa ni en alajuela ni en TNF y no sé si se debería paralelizar.

	
	// save a quick SVG cross section through z = 0, after setting its 
	// length bar to 200 microns 

	PhysiCell_SVG_options.length_bar = 200; 

	// for simplicity, set a pathology coloring function 
	
	std::vector<std::string> (*cell_coloring_function)(Cell*) = paint_by_number_cell_coloring; 
	// std::string (*ECM_coloring_function)(double, double, double) = my_coloring_function_for_stroma; 
	// std::vector<std::string> (*cell_coloring_function)(Cell*) = prolif_apoptosis_coloring;
	
	sprintf( filename , "%s/initial.svg" , PhysiCell_settings.folder.c_str() ); 
	SVG_plot_mpi( filename , microenvironment, 0.0 , PhysiCell_globals.current_time, cell_coloring_function, substrate_coloring_function, world, cart_topo );
	
	sprintf( filename , "%s/legend.svg" , PhysiCell_settings.folder.c_str() ); 
	create_plot_legend( filename , cell_coloring_function ); 

	if (IOProcessor(world)) {
		create_plot_legend( filename , cell_coloring_function ); 
		add_software_citation( "PhysiBoSS" , PhysiBoSS_Version , PhysiBoSS_DOI, PhysiBoSS_URL); 
		display_citations(); 
	}
	
	// set the performance timers 

	BioFVM::RUNTIME_TIC();
	BioFVM::TIC();
	
	std::ofstream report_file;
	
	// main loop 
	
	try 
	{		
		while( PhysiCell_globals.current_time < PhysiCell_settings.max_time + 0.1*diffusion_dt )
		{
			// save data if it's time. 
			// if( PhysiCell_globals.current_time > PhysiCell_globals.next_full_save_time - 0.5 * diffusion_dt )
			if( fabs( PhysiCell_globals.current_time - PhysiCell_globals.next_full_save_time ) < 0.01 * diffusion_dt )
			{
				display_simulation_status( std::cout, world, cart_topo );  
				if( PhysiCell_settings.enable_legacy_saves == true )
				{	
					log_output( PhysiCell_globals.current_time , PhysiCell_globals.full_output_index, microenvironment, report_file);
				}
				
				if( PhysiCell_settings.enable_full_saves == true )
				{	
					sprintf( filename , "%s/output%08u" , PhysiCell_settings.folder.c_str(),  PhysiCell_globals.full_output_index ); 
					
					save_PhysiCell_to_MultiCellDS_v2( filename , microenvironment , PhysiCell_globals.current_time , world, cart_topo); // Ultimos cambios echos aquí
					sprintf( filename , "%s/states_%08u.csv", PhysiCell_settings.folder.c_str(), PhysiCell_globals.full_output_index);
					MaBoSSIntracellular::save( filename); // JM: Aparición 2

				}
				
				PhysiCell_globals.full_output_index++; 
				PhysiCell_globals.next_full_save_time += PhysiCell_settings.full_save_interval;
			}
			
			// save SVG plot if it's time
			// if( PhysiCell_globals.current_time > PhysiCell_globals.next_full_save_time - 0.5 * diffusion_dt )
			if( fabs( PhysiCell_globals.current_time - PhysiCell_globals.next_SVG_save_time  ) < 0.01 * diffusion_dt )
			{
				if( PhysiCell_settings.enable_SVG_saves == true )
				{	
					sprintf( filename , "%s/snapshot%08u.svg" , PhysiCell_settings.folder.c_str() , PhysiCell_globals.SVG_output_index ); 
					SVG_plot_mpi( filename , microenvironment, 0.0 , PhysiCell_globals.current_time, cell_coloring_function,substrate_coloring_function, world, cart_topo);
					
					PhysiCell_globals.SVG_output_index++; 
					PhysiCell_globals.next_SVG_save_time  += PhysiCell_settings.SVG_save_interval;
				}
			}
			// Configure treatments
			// treatment_function();
	
			// update the microenvironment
			microenvironment.simulate_diffusion_decay( diffusion_dt, world, cart_topo );
			
			// run PhysiCell 
			((Cell_Container *)microenvironment.agent_container)->update_all_cells( PhysiCell_globals.current_time, world, cart_topo );

			/*
			  Custom add-ons could potentially go here. 
			*/

			PhysiCell_globals.current_time += diffusion_dt;
		}

		if( PhysiCell_settings.enable_legacy_saves == true )
		{			
			log_output(PhysiCell_globals.current_time, PhysiCell_globals.full_output_index, microenvironment, report_file);
			if (world.rank == 0)
		{
			report_file.close();
		}
		}
	}
	catch( const std::exception& e )
	{ // reference to the base of a polymorphic object
		std::cout << e.what(); // information from length_error printed
	}
	
	// save a final simulation snapshot 
	
	sprintf( filename , "%s/final" , PhysiCell_settings.folder.c_str() ); 
	save_PhysiCell_to_MultiCellDS_v2( filename , microenvironment , PhysiCell_globals.current_time , world, cart_topo );

	sprintf( filename , "%s/states_final.csv", PhysiCell_settings.folder.c_str());
	MaBoSSIntracellular::save( filename );
	
	sprintf( filename , "%s/final.svg" , PhysiCell_settings.folder.c_str() ); 
	SVG_plot_mpi( filename , microenvironment, 0.0 , PhysiCell_globals.current_time, cell_coloring_function, substrate_coloring_function, world, cart_topo );

	
	// timer 
	
	std::cout << std::endl << "Total simulation runtime: " << std::endl; 
	BioFVM::display_stopwatch_value( std::cout , BioFVM::runtime_stopwatch_value() ); 

	return 0; 
}
