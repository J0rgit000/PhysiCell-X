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
+																																							   +
+ [1] Cite the PhysiCell-X repository by giving its URL												   +
+																																							   +
+ [2] Cite BioFVM-X: 																													   +
+		Saxena, Gaurav, Miguel Ponce-de-Leon, Arnau Montagud, David Vicente Dorca,   +
+		and Alfonso Valencia. "BioFVM-X: An MPI+ OpenMP 3-D Simulator for Biological + 
+		Systems." In International Conference on Computational Methods in Systems    +
+		Biology, pp. 266-279. Springer, Cham, 2021. 																 +
=================================================================================*/

/*=======================================*/
/* Include mpi.h in the parallel version */
/*=======================================*/

#include <mpi.h>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <ctime>
#include <cmath>
#include <omp.h>
#include <fstream>
#include <string> 

/*=============================*/
/* Include ./DistPhy/DistPhy.h */
/*=============================*/

#include "./core/PhysiCell.h"
#include "./modules/PhysiCell_standard_modules.h" 

//Custom user modules 

#include "./custom_modules/heterogeneity.h" 

using namespace BioFVM;
using namespace PhysiCell;

/*============================*/
/* Use DistPhy::mpi namespace */
/*============================*/

using namespace DistPhy::mpi; 

int main( int argc, char* argv[] )
{
    
/*=======================================================================================*/
/* Create mpi_Environment object, initialize it, then create Cartesian Topology          */
/*=======================================================================================*/
	
	mpi_Environment world;                         //object contains size of communicator, rank of process
    world.Initialize();                            //Initialize using MPI_THREAD_FUNNELED, find comm. size and comm. rank
    mpi_Cartesian cart_topo;                       //Contains dims[3], coords[3] array and MPI_Comm mpi_cart_comm
    cart_topo.Build_Cartesian_Topology(world);     //Create 1-D X decomposition by setting dims[1]=size. 
    cart_topo.Find_Cartesian_Coordinates(world);   //Find Cartesian Topology coordinates of each process
    cart_topo.Find_Left_Right_Neighbours(world); 	 //Finds left/right immediate neighbour processes ranks and stores in X_LEFT/X_RIGHT
    
    bool XML_status = false;
    
/*====================================================================================*/    
/* Call parallel version of load_ function, all processes must load the complete file */
/*====================================================================================*/    

    
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

	//Setting OpenMP threads, this takes precedence over OMP_NUM_THREADS so be careful
	omp_set_num_threads(PhysiCell_settings.omp_num_threads);
	
	//PNRG setup 
	SeedRandom(); 
	
	//Time units setup 
	std::string time_units = "min"; 

/*================================================================================================*/
/* Objects of both mpi_Environment and mpi_Cartesian are needed in setup microenvironment(...) 		*/
/* These objects are passed to initialize_microenvironment(...) which calls resize_uniform(...)   */
/*================================================================================================*/
	
	setup_microenvironment(world, cart_topo); 

	//PhysiCell setup
 	
	//Set mechanics voxel size, must be >= Diffusion voxel size
	double mechanics_voxel_size = 20;
    
/*=========================================================*/
/* Calling the parallel version of Cell Container creation */
/*=========================================================*/
    
	Cell_Container* cell_container = create_cell_container_for_microenvironment( microenvironment, mechanics_voxel_size, world, cart_topo);
	
	create_cell_types(world, cart_topo);
    
	setup_tissue(microenvironment, world, cart_topo);      //Custom parallel function 
	
	/* Users typically start modifying here. START USERMODS */ 
	
	/*......................................................*/
	
	/* Users typically stop modifying here. END USERMODS 		*/ 
	
	//Set MultiCellDS save options, these MUST all be true 
	set_save_biofvm_mesh_as_matlab( true ); 
	set_save_biofvm_data_as_matlab( true ); 
	set_save_biofvm_cell_data( true ); 
	set_save_biofvm_cell_data_as_custom_matlab( true );
	
	//Save a simulation snapshot 
	char filename[1024];
	sprintf( filename , "%s/initial" , PhysiCell_settings.folder.c_str() ); 
	
	//Use the parallel version of the function for XML file writing
	save_PhysiCell_to_MultiCellDS_xml_pugi( filename , microenvironment , PhysiCell_globals.current_time, world, cart_topo ); 
	
	//Save a SVG cross section through z = 0, after setting its length bar to 200 microns 
	PhysiCell_SVG_options.length_bar = 200; 

	//For simplicity, set a pathology coloring function 
	std::vector<std::string> (*cell_coloring_function)(Cell*) = heterogeneity_coloring_function;
	std::string (*substrate_coloring_function)(double, double, double) = paint_by_density_percentage;
	
	//Use the parallel version of the function for SVG file plotting
	sprintf( filename , "%s/initial.svg" , PhysiCell_settings.folder.c_str() ); 
    SVG_plot_mpi( filename , microenvironment, 0.0 , PhysiCell_globals.current_time, cell_coloring_function, substrate_coloring_function, world, cart_topo );
	
	//Set the performance timers 
	BioFVM::RUNTIME_TIC();
	BioFVM::TIC();
	
	std::ofstream report_file;
	if( world.rank == 0 )
	{
        sprintf(filename , "%s/simulation_report.tsv" , PhysiCell_settings.folder.c_str() );
        report_file.open(filename);     // create the data log file 
        report_file << "timepoint";
        report_file << "\tbasic_agents\tcell_agents\talive\tdead\tonly_apoptotic\tnecrotic\tdead_but_ghost\tnecrotic_ghost"<<std::endl;;
		// report_file <<	"apoptotic_and_14\tapoptotic_not_14" <<std::endl;;
		std::cout   << "Apoptotic Code "<< PhysiCell_constants::apoptotic <<std::endl;;
		std::cout   << "necrotic_swelling Code "<< PhysiCell_constants::necrotic_swelling << std::endl;;
		std::cout   << "necrotic_lysed Code "<< PhysiCell_constants::necrotic_lysed << std::endl;;
		std::cout   << "necrotic Code "<< PhysiCell_constants::necrotic  << std::endl;;
	}
	//Main loop of the program
	try 
	{	
		while( PhysiCell_globals.current_time < PhysiCell_settings.max_time + 0.1*diffusion_dt )
		{
			//Save data if it is time. 
			if( fabs( PhysiCell_globals.current_time - PhysiCell_globals.next_full_save_time ) < 0.01 * diffusion_dt )
			{
				
				//Use the parallel version of the function to display total cells, time
				display_simulation_status( std::cout, world, cart_topo );
				
                double timepoint        = PhysiCell_globals.current_time;
				int basic_agents        = total_basic_agent_count(world, cart_topo);
				int cell_agents         = total_cell_agent_count(world, cart_topo);
                int alive               = total_live_cell_count(world, cart_topo);
                int dead                = total_dead_cell_count(world, cart_topo);
                int apoptotic           = total_apoptosis_cell_count(world, cart_topo);
                int necrotic            = total_necrosis_cell_count(world, cart_topo);
				int dead_but_ghost 		= total_dead_but_ghost(world, cart_topo);
				int necrotic_ghost 		= total_necrotic_ghost(world, cart_topo);
				// int apoptotic_dead_alive = total_apoptosis_and_codealive(world, cart_topo);
				// int apoptotic_dead_not_alive = total_apoptosis_not_codealive(world, cart_topo);
				if( world.rank == 0) 
				{
                    report_file << PhysiCell_globals.current_time;
                    report_file << "\t" << basic_agents <<"\t" << cell_agents << "\t" << alive;
					report_file << "\t" << dead << "\t" << apoptotic << "\t" << necrotic << "\t" << dead_but_ghost << "\t" << necrotic_ghost <<std::endl;;
					// report_file << "\t" << apoptotic_dead_alive << "\t" << apoptotic_dead_not_alive << std::endl;;

				}
				
				
				if( PhysiCell_settings.enable_full_saves == true )
				{	
					sprintf( filename , "%s/output%08u" , PhysiCell_settings.folder.c_str(),  PhysiCell_globals.full_output_index ); 
					
					//Use the parallel version of the function for XML file writing
					save_PhysiCell_to_MultiCellDS_xml_pugi( filename , microenvironment , PhysiCell_globals.current_time, world, cart_topo ); 
				}
				
				PhysiCell_globals.full_output_index++; 
				PhysiCell_globals.next_full_save_time += PhysiCell_settings.full_save_interval;
			}
			
			//Save SVG plot if it is time
			if( fabs( PhysiCell_globals.current_time - PhysiCell_globals.next_SVG_save_time  ) < 0.01 * diffusion_dt )
			{
				if( PhysiCell_settings.enable_SVG_saves == true )
				{	
										
					//Use the parallel version of the function for SVG file plotting
					sprintf( filename , "%s/snapshot%08u.svg" , PhysiCell_settings.folder.c_str() , PhysiCell_globals.SVG_output_index);
					SVG_plot_mpi( filename , microenvironment, 0.0 , PhysiCell_globals.current_time, cell_coloring_function,substrate_coloring_function, world, cart_topo);
					
					PhysiCell_globals.SVG_output_index++; 
					PhysiCell_globals.next_SVG_save_time  += PhysiCell_settings.SVG_save_interval;
				}
			}
			
			//Update the microenvironment i.e. solve diffusion equations
			microenvironment.simulate_diffusion_decay( diffusion_dt, world, cart_topo );
			
			//Run PhysiCell/PhysiCell-X to update all cells 
			((Cell_Container *)microenvironment.agent_container)->update_all_cells( PhysiCell_globals.current_time, world, cart_topo);
			
			PhysiCell_globals.current_time += diffusion_dt;
		}
			report_file.close();
		
	}
	catch( const std::exception& e )
	{ 
		//Reference to the base of a polymorphic object, information from length_error printed
		std::cout << e.what();  
	}
	
	//Save a final simulation snapshot 
	sprintf( filename , "%s/final" , PhysiCell_settings.folder.c_str() ); 
	save_PhysiCell_to_MultiCellDS_xml_pugi( filename , microenvironment , PhysiCell_globals.current_time, world, cart_topo ); 
	
	sprintf( filename , "%s/final.svg" , PhysiCell_settings.folder.c_str()); 
	SVG_plot_mpi( filename , microenvironment, 0.0 , PhysiCell_globals.current_time, cell_coloring_function, substrate_coloring_function, world, cart_topo );
	
	//Timer, printing is done ONLY through MPI process rank 0 i.e. the first MPI process 
	if(IOProcessor(world))
	{
		std::cout << std::endl << "Total simulation runtime: " << std::endl; 
		BioFVM::display_stopwatch_value( std::cout , BioFVM::runtime_stopwatch_value() ); 
	}
	
	//Gracefully shut-down MPI 
  	world.Finalize(); 

	return 0; 
}

