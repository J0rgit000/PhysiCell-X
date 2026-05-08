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

#include "./heterogeneity.h"
#include "../modules/PhysiCell_settings.h"
#include "../DistPhy/DistPhy_Utils.h"
#include "../DistPhy/DistPhy_Collective.h"

using namespace DistPhy::mpi;

void create_cell_types( mpi_Environment &world, mpi_Cartesian &cart_topo )
{
	// Use the same random seed so that future experiments have the 
	// same initial histogram of oncoprotein, even if threading means 
	// that future division and other events are still not identical 
	// for all runs 
	
	if (parameters.ints.find_index("random_seed") != -1)
	{
		SeedRandom(parameters.ints("random_seed"));
	}
	
	// Housekeeping 
	
	initialize_default_cell_definition();
	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment ); 
	
	cell_defaults.functions.volume_update_function = standard_volume_update_function;
	cell_defaults.functions.update_velocity = standard_update_cell_velocity;
	cell_defaults.functions.update_velocity_parallel = standard_update_cell_velocity;

	cell_defaults.functions.update_migration_bias = NULL; 
	cell_defaults.functions.update_phenotype = NULL; // update_cell_and_death_parameters_O2_based; 
	cell_defaults.functions.custom_cell_rule = NULL; 
	cell_defaults.functions.contact_function = NULL; 
	
	cell_defaults.functions.add_cell_basement_membrane_interactions = NULL; 
	cell_defaults.functions.calculate_distance_to_membrane = NULL; 

	initialize_cell_definitions_from_pugixml(world, cart_topo); 
	
	// Set default uptake and secretion 
	
	build_cell_definitions_maps(); 

	setup_signal_behavior_dictionaries(world, cart_topo); 


	/* 
	   Put any modifications to individual cell definitions here. 
	   
	   This is a good place to set custom functions. 
	*/ 
	
	static int oxygen_ID = microenvironment.find_density_index( "oxygen" );
	
	Cell_Definition* pCD = find_cell_definition( "cancer cell"); 
	pCD->functions.update_phenotype = tumor_cell_phenotype_with_oncoprotein; 
	pCD->functions.update_phenotype_parallel = tumor_cell_phenotype_with_oncoprotein;

	pCD->parameters.o2_proliferation_saturation = 38; 
	pCD->parameters.o2_reference = 38; 

	
	display_cell_definitions( std::cout, world, cart_topo );

	return; 
}

void phenotype_function( Cell* pCell, Phenotype& phenotype, double dt )
{ return; }

void custom_function( Cell* pCell, Phenotype& phenotype , double dt )
{ return; } 

void contact_function( Cell* pMe, Phenotype& phenoMe , Cell* pOther, Phenotype& phenoOther , double dt )
{ return; } 


void setup_microenvironment(mpi_Environment &world, mpi_Cartesian &cart_topo)
{
	
	//PhysiCell-X ONLY supports 3-D problems, hence MUST set 3-D option
	
	if( default_microenvironment_options.simulate_2D == true )
	{
		if(IOProcessor(world))
            std::cout << "Warning: overriding XML config option and setting to 3D!" << std::endl; 
		default_microenvironment_options.simulate_2D = false; 
	}
			
	initialize_microenvironment(world, cart_topo); 	

	return; 
}	

/*-------------------------------------------------------------------*/
/* Miguel Ponce-de-Leon's function for generating positions of cells */
/*-------------------------------------------------------------------*/

std::vector<std::vector<double>> create_cell_sphere_positions(double cell_radius, double sphere_radius)
{
	std::vector<std::vector<double>> cells;
	int xc=0,yc=0,zc=0;
	double x_spacing= cell_radius*sqrt(3);
  double y_spacing= cell_radius*2;
  double z_spacing= cell_radius*sqrt(3);
	
	std::vector<double> tempPoint(3,0.0);

	for(double z=-sphere_radius;z<sphere_radius;z+=z_spacing, zc++)
	{
		for(double x=-sphere_radius;x<sphere_radius;x+=x_spacing, xc++)
		{
			for(double y=-sphere_radius;y<sphere_radius;y+=y_spacing, yc++)
			{
				tempPoint[0]=x + (zc%2) * 0.5 * cell_radius;
				tempPoint[1]=y + (xc%2) * cell_radius;
				tempPoint[2]=z;

				if(sqrt(norm_squared(tempPoint))< sphere_radius)
				{ cells.push_back(tempPoint); }
			}

		}
	}
	return cells;

}


/*------------------------------------------------------------------------*/
/* Parallel version of setup_tissue(), replacing this function completely */
/* by Miguel's version of setup_tissue and then parallelizing             */
/*------------------------------------------------------------------------*/

void setup_tissue(Microenvironment &m, mpi_Environment &world, mpi_Cartesian &cart_topo)
{

	double Xmin = microenvironment.mesh.x_coordinates[0]; 
	double Ymin = microenvironment.mesh.y_coordinates[0]; 
	double Zmin = microenvironment.mesh.z_coordinates[0];

	int x_voxels = microenvironment.mesh.x_coordinates.size();
	int y_voxels = microenvironment.mesh.y_coordinates.size();
	int z_voxels = microenvironment.mesh.z_coordinates.size();

	double Xmax = microenvironment.mesh.x_coordinates[x_voxels - 1]; 
	double Ymax = microenvironment.mesh.y_coordinates[y_voxels - 1]; 
	double Zmax = microenvironment.mesh.z_coordinates[z_voxels - 1];

	std::vector<double> position(3,0.0); 
	
	double x_range = Xmax - Xmin; 
	double y_range = Ymax - Ymin;
	double z_range = Zmax - Zmin; 

	double relative_margin = 0.2;  
	double relative_outer_margin = 0.02; 
	
	bool use_csv_positions = false;
	pugi::xml_node node_ic = xml_find_node( physicell_config_root , "initial_conditions" ); 
	if( node_ic )
	{
		pugi::xml_node node_cp = xml_find_node( node_ic , "cell_positions" );
		if( node_cp && node_cp.attribute("enabled").as_bool() )
		{ use_csv_positions = true; }
	}

	for( int k=0; k < cell_definitions_by_index.size() ; k++ )
	{
		Cell* pC;
		Cell_Definition* pCD = cell_definitions_by_index[k]; 
		std::cout << "Placing cells of type " << pCD->name << " ... " << std::endl; 
		for( int n = 0 ; n < parameters.ints("number_of_cells") / world.size ; n++ )
		{
			std::vector<double> position = {0,0,0}; 
			position[0] = Xmin + UniformRandom()*x_range; 
			position[1] = Ymin + UniformRandom()*y_range; 
			position[2] = Zmin + UniformRandom()*z_range; 
			
			pC = create_cell( *pCD ); 
			pC->assign_position( position, world, cart_topo );
		}
	}
	std::cout << std::endl; 

    // Custom placement

	Cell_Definition* pCD = find_cell_definition( "cancer cell"); 
	double cell_radius = pCD->phenotype.geometry.radius; 
	double cell_spacing = 0.95 * 2.0 * cell_radius; 
	
	double tumor_radius = parameters.doubles( "tumor_radius" );
	
	//int i = parameters.doubles.find_index( "tumor_radius" ); 
	
	Cell* pCell = NULL; 
    
    std::vector<std::vector<double>> positions;		 
    std::vector<std::vector<double>> generated_positions_at_root;
    
    /*----------------------------------------------------------------------------------------------------*/
    /* Object of mpi_CellPositions must be declared for all processes because distribute_cell_positions() */
    /* function will pass 2 objects of the kind mpi_CellPositions and mpi_MyCells                         */
    /*----------------------------------------------------------------------------------------------------*/
    
    mpi_CellPositions cp;                //To store cell positions, cell IDs, no. of cell IDs at root only (for all processes)
    mpi_MyCells       mc;                //To store cell positions, cell IDs, no. of cells at each process.
	
    if(world.rank == 0) //Only the MPI Rank 0 process will generate positions
    {
        generated_positions_at_root = create_cell_sphere_positions(cell_radius,tumor_radius);   //Generate the cell positions
        
        int strt_cell_ID = Basic_Agent::get_max_ID_in_parallel();                               //IDs for new cells (positions) will start from the current highest ID
        
        
        cp.positions_to_rank_list(generated_positions_at_root, 
                                  m.mesh.bounding_box[0], m.mesh.bounding_box[3], m.mesh.bounding_box[1], m.mesh.bounding_box[4], m.mesh.bounding_box[2], m.mesh.bounding_box[5], 
                                  m.mesh.dx, m.mesh.dy, m.mesh.dz, 
                                  world, cart_topo, strt_cell_ID);
        
        Basic_Agent::set_max_ID_in_parallel(strt_cell_ID + generated_positions_at_root.size()); //Highest ID now is the starting ID + no. of generated coordinates ! 
    }
    
    distribute_cell_positions(cp, mc, world, cart_topo);                                        //Distribute cell positions to individual processes
	
    if(IOProcessor(world))
        std::cout << "creating " << generated_positions_at_root.size() << " closely-packed tumor cells ... " << std::endl;

	double x = 0.0; 
	double x_outer = tumor_radius; 
	double y = 0.0; 
	
		for( int i=0; i < mc.my_no_of_cell_IDs; i++ )
		{
		  
			pCell = create_cell(*pCD, mc.my_cell_IDs[i]); // tumor cell --> This has to be replaced by create_cell(mc.my_cell_IDs[i])
	  	  		
			pCell->assign_position(mc.my_cell_coords[3*i],mc.my_cell_coords[3*i+1],mc.my_cell_coords[3*i+2],world, cart_topo); //pCell->assign_position( positions[i] );
		} 
 
	if( use_csv_positions )
	{
		if( IOProcessor(world) )
		{
			std::cout << "cell_positions enabled in XML; loading positions from file." << std::endl;
		}

		double half_delta = m.mesh.dx/2.0;
		std::pair<double,double> x_range = {Xmin - half_delta, Xmax + half_delta };
		std::cout << "Rank: " << world.rank << " x_range: [" << x_range.first << " " << x_range.second << "]" << std::endl;
		load_cells_from_pugixml(world, cart_topo, x_range); 
	}

	Cell_Definition* pCD_onco = find_cell_definition( "cancer cell"); 
	if( pCD_onco != NULL )
	{
		double p_mean = parameters.doubles( "oncoprotein_mean" ); 
		double p_sd 	= parameters.doubles( "oncoprotein_sd" ); 
		double p_min 	= parameters.doubles( "oncoprotein_min" ); 
		double p_max 	= parameters.doubles( "oncoprotein_max" ); 

		for( int i=0; i < all_cells->size(); i++ )
		{
			Cell* pCell = (*all_cells)[i];
			if( pCell->type != pCD_onco->type )
			{ continue; }
			double p = NormalRandom( p_mean, p_sd );
			if( p < p_min ) { p = p_min; }
			if( p > p_max ) { p = p_max; }
			set_single_behavior( pCell, "custom:oncoprotein" , p );
		}
	}

	int local_cells = 0; 
	local_cells = (int) all_cells->size(); 
	
	double local_sum = 0.0, global_sum = 0.0; 
	double local_min = 9e9, global_min = 0.0; 
	double local_max = -9e9, global_max = 0.0; 
	int global_cells; 

/*------------------------------------------------------------------------------------------------------*/
/* The global_min/max are only for display purposes. We can just calculate local_min then do MPI_Reduce */
/* at the root process to display it. Since the mean is needed by all processes to calculate squared 	  */
/* sum of differences, I should do MPI_Allreduce(). Right now everything is globally distributed			  */
/* later we can decide to use MPI_Reduce selectively. 																									*/
/*------------------------------------------------------------------------------------------------------*/	
	
	for( int i=0; i < all_cells->size() ; i++ )
	{
		double r = (*all_cells)[i]->custom_data[0]; 
		local_sum += r;
		if( r < local_min )
		{ 
			local_min = r; 
		} 
		if( r > local_max )
		{ 
			local_max = r; 
		}
	}
	global_sum 		= distribute_global_sum(local_sum, cart_topo);
	global_cells 	= distribute_global_sum(local_cells, cart_topo); 
	global_max		= distribute_global_max(local_max, cart_topo); 
	global_min	  = distribute_global_min(local_min, cart_topo); 
	 
	double mean = global_sum / ( global_cells + 1e-15 ); 
	
	// Compute standard deviation 
	
	local_sum = 0.0;
	 
	for( int i=0; i < all_cells->size(); i++ )
	{
		local_sum +=  ( (*all_cells)[i]->custom_data[0] - mean )*( (*all_cells)[i]->custom_data[0] - mean ); 
	}
	global_sum = distribute_global_sum(local_sum, cart_topo);
		
	double standard_deviation = sqrt( global_sum / ( global_cells - 1.0 + 1e-15 ) ); 
	
	if(IOProcessor(world))
	{
		std::cout << std::endl << "Oncoprotein summary: " << std::endl<< "===================" << std::endl; 
		std::cout << "mean: " << mean << std::endl; 
		std::cout << "standard deviation: " << standard_deviation << std::endl; 
		std::cout << "[min max]: [" << global_min << " " << global_max << "]" << std::endl << std::endl;
	} 
	
	return; 
}

// Custom cell phenotype function to scale immunostimulatory factor with hypoxia 

void tumor_cell_phenotype_with_oncoprotein( Cell* pCell, Phenotype& phenotype, double dt )
{
	update_cell_and_death_parameters_O2_based(pCell,phenotype,dt);
	
	// If cell is dead, don't bother with future phenotype changes. 
	
	if( get_single_signal( pCell, "dead") > 0.5 )  // ← SIGNAL CHECK
    {
        pCell->functions.update_phenotype = NULL;
        return;
    }

	// Multiply proliferation rate by the oncoprotein 
	
	double cycle_rate = get_single_behavior( pCell, "cycle entry" );
	cycle_rate *= get_single_signal( pCell , "custom:oncoprotein" );
	set_single_behavior( pCell, "cycle entry" , cycle_rate );
	
	return; 
}
void tumor_cell_phenotype_with_oncoprotein( Cell* pCell, Phenotype& phenotype, double dt,mpi_Environment &world, mpi_Cartesian &cart_topo )
{
	update_cell_and_death_parameters_O2_based(pCell,phenotype,dt);
	
	// If cell is dead, don't bother with future phenotype changes. 
	
	if( phenotype.death.dead == true )
	{
		pCell->functions.update_phenotype = NULL;
		pCell->functions.update_phenotype_parallel = NULL;
		return; 
	}

	// Multiply proliferation rate by the oncoprotein 
	
	static int cycle_start_index = live.find_phase_index( PhysiCell_constants::live ); 
	static int cycle_end_index = live.find_phase_index( PhysiCell_constants::live ); 
	static int oncoprotein_i = pCell->custom_data.find_variable_index( "oncoprotein" ); 

	phenotype.cycle.data.transition_rate( cycle_start_index ,cycle_end_index ) *= pCell->custom_data[oncoprotein_i] ; 
	
	return; 
}


std::vector<std::string> heterogeneity_coloring_function( Cell* pCell )
{
	static int oncoprotein_i = pCell->custom_data.find_variable_index( "oncoprotein" ); 
	
	static double p_min = parameters.doubles( "oncoprotein_min" ); 
	static double p_max = parameters.doubles( "oncoprotein_max" ); 
	
	// Immune are black
	
	std::vector< std::string > output( 4, "black" ); 
	
	if( pCell->type == 1 )
	{ 
		return output; 
	} 
	
	// Live cells are green, but shaded by oncoprotein value
	 
	if( pCell->phenotype.death.dead == false )
	{
		int oncoprotein = (int) round( (1.0/(p_max-p_min)) * (pCell->custom_data[oncoprotein_i]-p_min) * 255.0 ); 
		char szTempString [128];
		sprintf( szTempString , "rgb(%u,%u,%u)", oncoprotein, oncoprotein, 255-oncoprotein );
		output[0].assign( szTempString );
		output[1].assign( szTempString );

		sprintf( szTempString , "rgb(%u,%u,%u)", (int)round(output[0][0]/p_max) , (int)round(output[0][1]/p_max) , (int)round(output[0][2]/p_max) );
		output[2].assign( szTempString );
		
		return output; 
	}

	// If not, dead colors 
	
	if (pCell->phenotype.cycle.current_phase().code == PhysiCell_constants::apoptotic )  // Apoptotic - Red
	{
		output[0] = "rgb(255,0,0)";
		output[2] = "rgb(125,0,0)";
	}
	
	// Necrotic - Brown
	
	if( pCell->phenotype.cycle.current_phase().code == PhysiCell_constants::necrotic_swelling || 
		pCell->phenotype.cycle.current_phase().code == PhysiCell_constants::necrotic_lysed || 
		pCell->phenotype.cycle.current_phase().code == PhysiCell_constants::necrotic )
	{
		output[0] = "rgb(250,138,38)";
		output[2] = "rgb(139,69,19)";
	}	
	
	return output; 
}
int total_basic_agent_count(mpi_Environment &world, mpi_Cartesian &cart_topo)
{
	int local_count = all_basic_agents.size();
	int global_count;
	MPI_Reduce(&local_count, &global_count, 1, MPI_INT, MPI_SUM, 0, cart_topo.mpi_cart_comm);
	return global_count;
}



int total_cell_agent_count(mpi_Environment &world, mpi_Cartesian &cart_topo)
{
	int local_count = (*all_cells).size();
	int global_count;
	MPI_Reduce(&local_count, &global_count, 1, MPI_INT, MPI_SUM, 0, cart_topo.mpi_cart_comm);
	return global_count;
}


// Cell count functions
int total_live_cell_count(mpi_Environment &world, mpi_Cartesian &cart_topo)
{
	int local_count = 0;
	#pragma omp parallel for reduction (+: local_count )
	for (int i = 0; i < (*all_cells).size(); i++)
	{
		Cell* pCell = (*all_cells)[i];
		if(pCell->phenotype.death.dead==false)
		{ local_count++;; }
	}
	int global_count;
	MPI_Reduce(&local_count, &global_count, 1, MPI_INT, MPI_SUM, 0, cart_topo.mpi_cart_comm);
	return global_count;
}

int total_dead_cell_count(mpi_Environment &world, mpi_Cartesian &cart_topo)
{
    int local_count = 0;
    #pragma omp parallel for reduction (+: local_count )
    for (int i = 0; i < (*all_cells).size(); i++)
    {
        Cell* pCell = (*all_cells)[i];
        if(pCell->phenotype.death.dead==true)
        { local_count++; }
    }
    int global_count;
    MPI_Reduce(&local_count, &global_count, 1, MPI_INT, MPI_SUM, 0, cart_topo.mpi_cart_comm);
    return global_count;
}

int total_necrosis_cell_count(mpi_Environment &world, mpi_Cartesian &cart_topo)
{
	int local_count = 0;
	#pragma omp parallel for reduction (+: local_count )
	for (int i = 0; i < (*all_cells).size(); i++)
	{
		Cell* pCell = (*all_cells)[i];
		if( pCell->phenotype.death.dead==false )
		{ continue; }
		if( pCell->phenotype.cycle.current_phase().code == PhysiCell_constants::necrotic_swelling || 
			pCell->phenotype.cycle.current_phase().code == PhysiCell_constants::necrotic_lysed || 
			pCell->phenotype.cycle.current_phase().code == PhysiCell_constants::necrotic )	
		{ local_count++; }
	}
	int global_count;
	MPI_Reduce(&local_count, &global_count, 1, MPI_INT, MPI_SUM, 0, cart_topo.mpi_cart_comm);
	return global_count;
}

int total_apoptosis_cell_count(mpi_Environment &world, mpi_Cartesian &cart_topo)
{	
	int local_count = 0;
	#pragma omp parallel for reduction (+: local_count )
	for (int i = 0; i < (*all_cells).size(); i++)
	{
		Cell* pCell = (*all_cells)[i];
		if( pCell->phenotype.death.dead==false )
		{ continue; }
		if( pCell->phenotype.cycle.current_phase().code == PhysiCell_constants::apoptotic)
		{ local_count++; }
	}
	int global_count;
	MPI_Reduce(&local_count, &global_count, 1, MPI_INT, MPI_SUM, 0, cart_topo.mpi_cart_comm);
	return global_count;
}


int total_dead_but_ghost(mpi_Environment &world, mpi_Cartesian &cart_topo)
{	
	int local_count = 0;
	#pragma omp parallel for reduction (+: local_count )
	for (int i = 0; i < (*all_cells).size(); i++)
	{
		Cell* pCell = (*all_cells)[i];
		if( pCell->phenotype.death.dead==false )
		{ continue; }
		if( pCell->phenotype.cycle.current_phase().code != PhysiCell_constants::apoptotic &&
			pCell->phenotype.cycle.current_phase().code != PhysiCell_constants::necrotic_swelling && 
			pCell->phenotype.cycle.current_phase().code != PhysiCell_constants::necrotic_lysed && 
			pCell->phenotype.cycle.current_phase().code != PhysiCell_constants::necrotic )
		{ local_count++; 
			// std::cout<<pCell->phenotype.cycle.current_phase().code<<" Line written by Rank:"<<world.rank<<std::endl;
		}
	}
	int global_count;
	MPI_Reduce(&local_count, &global_count, 1, MPI_INT, MPI_SUM, 0, cart_topo.mpi_cart_comm);
	return global_count;
}


int total_apoptosis_not_codealive(mpi_Environment &world, mpi_Cartesian &cart_topo)
{	
	int local_count = 0;
	#pragma omp parallel for reduction (+: local_count )
	for (int i = 0; i < (*all_cells).size(); i++)
	{
		Cell* pCell = (*all_cells)[i];
		if( pCell->phenotype.death.dead==false )
		{ continue; }
		if( pCell->phenotype.cycle.current_phase().code == PhysiCell_constants::apoptotic && pCell->phenotype.cycle.current_phase().code != 14 )
		{ local_count++; 
			// std::cout<<pCell->phenotype.cycle.current_phase().code<<" Line written by Rank:"<<world.rank<<std::endl;
		}
	}
	int global_count;
	MPI_Reduce(&local_count, &global_count, 1, MPI_INT, MPI_SUM, 0, cart_topo.mpi_cart_comm);
	return global_count;
}

int total_apoptosis_and_codealive(mpi_Environment &world, mpi_Cartesian &cart_topo)
{	
	int local_count = 0;
	#pragma omp parallel for reduction (+: local_count )
	for (int i = 0; i < (*all_cells).size(); i++)
	{
		Cell* pCell = (*all_cells)[i];
		if( pCell->phenotype.death.dead==false )
		{ continue; }
		if( pCell->phenotype.cycle.current_phase().code == PhysiCell_constants::apoptotic && pCell->phenotype.cycle.current_phase().code == 14 )
		{ local_count++; 
			// std::cout<<pCell->phenotype.cycle.current_phase().code<<" Line written by Rank:"<<world.rank<<std::endl;
		}
	}
	int global_count;
	MPI_Reduce(&local_count, &global_count, 1, MPI_INT, MPI_SUM, 0, cart_topo.mpi_cart_comm);
	return global_count;
}
int total_necrotic_ghost(mpi_Environment &world, mpi_Cartesian &cart_topo){
	int local_count = 0;
	#pragma omp parallel for reduction (+: local_count )
	for (int i = 0; i < (*all_cells).size(); i++)
	{
		Cell* pCell = (*all_cells)[i];
		if( pCell->phenotype.death.dead==false )
		{ continue; }
		if( pCell->phenotype.cycle.current_phase().code == PhysiCell_constants::necrotic_swelling && 
			pCell->phenotype.cycle.current_phase().code == PhysiCell_constants::necrotic_lysed && 
			pCell->phenotype.cycle.current_phase().code == PhysiCell_constants::necrotic)
			{
				if(pCell->phenotype.cycle.current_phase().code == 14)
					{local_count++; }
			}
		
	}
	int global_count;
	MPI_Reduce(&local_count, &global_count, 1, MPI_INT, MPI_SUM, 0, cart_topo.mpi_cart_comm);
	return global_count;


}
