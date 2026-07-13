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

#define _USE_MATH_DEFINES
#include <cmath>
#include <sstream>
#include "./custom.h"
#include "../BioFVM/BioFVM.h" 

/*namespace
{
bool is_necrotic_phase(const Cell* pCell)
{
	const int phase_code =
		pCell->phenotype.cycle.pCycle_Model->phases[pCell->phenotype.cycle.data.current_phase_index].code;
	return phase_code == PhysiCell_constants::necrotic_swelling ||
		   phase_code == PhysiCell_constants::necrotic_lysed ||
		   phase_code == PhysiCell_constants::necrotic;
}
} Comentado por no estar en la parte serial.*/

// declare cell definitions here 

// std::vector<bool> nodes;

bool is_necrotic_phase(const Cell* pCell)
{
	const int phase_code =
		pCell->phenotype.cycle.pCycle_Model->phases[pCell->phenotype.cycle.data.current_phase_index].code;
	return phase_code == PhysiCell_constants::necrotic_swelling ||
		   phase_code == PhysiCell_constants::necrotic_lysed ||
		   phase_code == PhysiCell_constants::necrotic;
}

void create_cell_types( mpi_Environment &world, mpi_Cartesian &cart_topo )
{
	// set the random seed 
	if (parameters.ints.find_index("random_seed") != -1)
	{
		SeedRandom(parameters.ints("random_seed"));
	}
	
	/* 
	   Put any modifications to default cell definition here if you 
	   want to have "inherited" by other cell types. 
	   
	   This is a good place to set default functions. 
	*/ 

	initialize_default_cell_definition(); 
	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment ); 

	cell_defaults.functions.volume_update_function = standard_volume_update_function;
	cell_defaults.functions.update_velocity = standard_update_cell_velocity;

	// cell_defaults.functions.pre_update_intracellular = pre_update_intracellular;  Estan comentandas en el serial.
	// cell_defaults.functions.post_update_intracellular = post_update_intracellular; 

	cell_defaults.functions.update_migration_bias = NULL;
	cell_defaults.functions.update_phenotype = NULL;  

	cell_defaults.functions.custom_cell_rule = NULL; 
	cell_defaults.functions.contact_function = NULL; 
	
	cell_defaults.functions.add_cell_basement_membrane_interactions = NULL; 
	cell_defaults.functions.calculate_distance_to_membrane = NULL; 
	
	// cell_defaults.custom_data.add_variable(parameters.strings("node_to_visualize"), "dimensionless", 0.0 ); //for paraview visualization

	/*
	   This parses the cell definitions in the XML config file. 
	*/
	
	initialize_cell_definitions_from_pugixml(world, cart_topo); 
	
	/* 
	   Put any modifications to individual cell definitions here. 
	   
	   This is a good place to set custom functions. 
	*/ 
	
	/*
	   This builds the map of cell definitions and summarizes the setup. 
	*/

	build_cell_definitions_maps(); 

	/*
	   This intializes cell signal and response dictionaries 
	*/

	setup_signal_behavior_dictionaries(world, cart_topo);
		/*
       Cell rule definitions 
	*/
	setup_cell_rules(world, cart_topo); 

	/* 
	   Put any modifications to individual cell definitions here. 
	   
	   This is a good place to set custom functions. 
	*/ 
	// > Parte serial, mirar si hay que paralelizar (hasta el <)
	cell_defaults.functions.pre_update_intracellular = pre_update_intracellular;
	cell_defaults.functions.post_update_intracellular = post_update_intracellular;
	
	cell_defaults.functions.update_phenotype = NULL; //phenotype_function; Funciones comentadas por no tener una definción clara
	cell_defaults.functions.custom_cell_rule = NULL; //custom_function; 
	cell_defaults.functions.contact_function = NULL; //contact_function; 

	Cell_Definition* pCD = find_cell_definition( "default");

	pCD->functions.pre_update_intracellular = pre_update_intracellular;
	pCD->functions.post_update_intracellular = post_update_intracellular;
	pCD->functions.custom_cell_rule = NULL; //custom_function; 
	pCD->functions.contact_function = NULL; //contact_function; 
	pCD->functions.update_velocity = standard_update_cell_velocity; 
	// <
	/*
	   This summarizes the setup. 
	*/
	
	display_cell_definitions( std::cout, world, cart_topo ); 


	return; 
}

double get_decay_rate(double half_life){
	// natural logarithm of 2 / half-life = k (reaction rate coefficient)
	double decay_rate = log(2)/(half_life);
	return decay_rate;
}

void setup_microenvironment( mpi_Environment &world, mpi_Cartesian &cart_topo )
{
	// set domain parameters 
	
	// put any custom code to set non-homogeneous initial conditions or 
	// extra Dirichlet nodes here. 
	
	// initialize BioFVM 

   /* if( default_microenvironment_options.simulate_2D == true )
	{
		if(IOProcessor(world))
            std::cout << "Warning: overriding XML config option and setting to 2D!" << std::endl; 
		default_microenvironment_options.simulate_2D = false; 
	}
		*/
	
	initialize_microenvironment(world, cart_topo); 	
	
	return; 
}

void setup_tissue(Microenvironment& microenvironment , mpi_Environment &world, mpi_Cartesian &cart_topo )
{
	// load cells from your CSV file
    std::pair<double, double> x_range = microenvironment.get_subdomain_x_limits();
    std::cout << "[Rank " << world.rank << "] Subdomain x-range: " << x_range.first << " to " << x_range.second << std::endl;
	load_cells_from_pugixml(world, cart_topo, x_range); 	
}

/*void pre_update_intracellular( Cell* pCell, Phenotype& phenotype, double dt )
{
	if (PhysiCell::PhysiCell_globals.current_time >= 100.0 
		&& pCell->phenotype.intracellular->get_parameter_value("$time_scale") == 0.0
	){
		pCell->phenotype.intracellular->set_parameter_value("$time_scale", 0.1);
	}

}

void post_update_intracellular( Cell* pCell, Phenotype& phenotype, double dt )
{
	color_node(pCell);
} Funciones comentadas por no estar en la parte serial.*/
/* old setup_tissue function, needs to be updated to new structure or removed
void setup_tissue_resistant( void )
{
	Cell* pC;

	std::vector<init_record> cells = read_init_file(parameters.strings("init_cells_filename"), ';', true);

	for (int i = 0; i < cells.size(); i++)
	{
		float x = cells[i].x;
		float y = cells[i].y;
		float z = cells[i].z;
		float radius = cells[i].radius;
		int phase = cells[i].phase;
		double elapsed_time = cells[i].elapsed_time;

		double random_num_1 = (double) rand()/RAND_MAX;
		double random_num_2 = (double) rand()/RAND_MAX;

		if (PhysiCell::parameters.ints("simulation_mode") == 0)
		{
			// single inhibition - just one drug is present 
			if (random_num_1 < PhysiCell::parameters.doubles("prop_drug_resistant_" + microenvironment.density_names[1]))
			{
				// cell is sensitive to the drug
				pC = create_cell(get_cell_definition(microenvironment.density_names[1] + "_resistant"));
			}
			else 
			{
				// cell is not sensitive to the drug
				pC = create_cell(get_cell_definition(microenvironment.density_names[1] + "_sensitive"));
			}
		}
		else if (PhysiCell::parameters.ints("simulation_mode") == 1)
		{
			// double inhibition - two drugs are present - we have 4 cell strains 
			if (random_num_1 < PhysiCell::parameters.doubles("prop_drug_resistant_" + microenvironment.density_names[1]))
			{
				if (random_num_2 < PhysiCell::parameters.doubles("prop_drug_resistant_" + microenvironment.density_names[2]))
				{
					// cell is resistant to both drugs
					pC = create_cell(get_cell_definition(microenvironment.density_names[1] + "_resistant"));
				}
				else 
				{
					// cell is only resistant to the first drug
					pC = create_cell(get_cell_definition(microenvironment.density_names[2] + "_sensitive"));
				}
			}
			else
			{
				if (random_num_2 < PhysiCell::parameters.doubles("prop_drug_resistant_" + microenvironment.density_names[2]))
				{
					// cell is only resistant to the second drug
					pC = create_cell(get_cell_definition(microenvironment.density_names[2] + "_resistant"));
				}
				else
				{
					// cell is resistant to no drug
					pC = create_cell(get_cell_definition(microenvironment.density_names[1] + "_sensitive"));
				}
				
			}
			
		}
		else
		{
			pC = create_cell(get_cell_definition("default"));
		}
 
		pC->assign_position( x, y, z );
		
		pC->phenotype.cycle.data.elapsed_time_in_phase = elapsed_time;	
		
		update_custom_variables(pC);
	}

	return; 
}
	*/
std::vector<std::string> my_coloring_function( Cell* pCell )
{
	std::vector< std::string > output( 4 , "rgb(0,0,0)" );
	
	if ( !pCell->phenotype.intracellular->get_boolean_variable_value( parameters.strings("node_to_visualize") ) )
	{
		output[0] = "rgb(255,0,0)";
		output[2] = "rgb(125,0,0)";
		
	}
	else{
		output[0] = "rgb(0, 255,0)";
		output[2] = "rgb(0, 125,0)";
	}
	
	return output;
}

void tumor_cell_phenotype_with_signaling( Cell* pCell, Phenotype& phenotype, double dt)
{
	update_cell_and_death_parameters_O2_based(pCell, phenotype, dt);

	// update motility state variable
	static int index_motility_state = pCell->custom_data.find_variable_index("motility_state");
	pCell->custom_data.variables.at(index_motility_state).value = int(pCell->phenotype.motility.is_motile);
	
	boolean_model_interface_main (pCell, phenotype, dt);
}

void color_node(Cell* pCell){
	std::string node_name = parameters.strings("node_to_visualize");
	pCell->custom_data[node_name] = pCell->phenotype.intracellular->get_boolean_variable_value(node_name);
}

int total_basic_agent_count(mpi_Environment &world, mpi_Cartesian &cart_topo)
{
	(void) world;
	const int local_count = static_cast<int>(all_basic_agents.size());
	return DistPhy::mpi::distribute_global_sum(local_count, cart_topo);
}

int total_cell_agent_count(mpi_Environment &world, mpi_Cartesian &cart_topo)
{
	(void) world;
	const int local_count = static_cast<int>((*all_cells).size());
	return DistPhy::mpi::distribute_global_sum(local_count, cart_topo);
}

int total_live_cell_count(mpi_Environment &world, mpi_Cartesian &cart_topo)
{
	(void) world;
	int local_count = 0;
	const int cell_count = static_cast<int>((*all_cells).size());

	#pragma omp parallel for reduction(+:local_count)
	for (int i = 0; i < cell_count; ++i)
	{
		Cell* pCell = (*all_cells)[i];
		if (!pCell->phenotype.death.dead)
		{
			++local_count;
		}
	}

	return DistPhy::mpi::distribute_global_sum(local_count, cart_topo);
}

int total_dead_cell_count(mpi_Environment &world, mpi_Cartesian &cart_topo)
{
	(void) world;
	int local_count = 0;
	const int cell_count = static_cast<int>((*all_cells).size());

	#pragma omp parallel for reduction(+:local_count)
	for (int i = 0; i < cell_count; ++i)
	{
		Cell* pCell = (*all_cells)[i];
		if (pCell->phenotype.death.dead)
		{
			++local_count;
		}
	}

	return DistPhy::mpi::distribute_global_sum(local_count, cart_topo);
}

int total_necrosis_cell_count(mpi_Environment &world, mpi_Cartesian &cart_topo)
{
	(void) world;
	int local_count = 0;
	const int cell_count = static_cast<int>((*all_cells).size());

	#pragma omp parallel for reduction(+:local_count)
	for (int i = 0; i < cell_count; ++i)
	{
		Cell* pCell = (*all_cells)[i];
		if (pCell->phenotype.death.dead && is_necrotic_phase(pCell))
		{
			++local_count;
		}
	}

	return DistPhy::mpi::distribute_global_sum(local_count, cart_topo);
}

int total_apoptosis_cell_count(mpi_Environment &world, mpi_Cartesian &cart_topo)
{
	(void) world;
	int local_count = 0;
	const int cell_count = static_cast<int>((*all_cells).size());

	#pragma omp parallel for reduction(+:local_count)
	for (int i = 0; i < cell_count; ++i)
	{
		Cell* pCell = (*all_cells)[i];
		if (pCell->phenotype.death.dead &&
			pCell->phenotype.cycle.current_phase().code == PhysiCell_constants::apoptotic)
		{
			++local_count;
		}
	}

	return DistPhy::mpi::distribute_global_sum(local_count, cart_topo);
}


std::vector<std::string> prolif_apoptosis_coloring( Cell* pCell )
{
	std::vector<std::string> output;
	if (pCell->phenotype.cycle.current_phase().code == PhysiCell_constants::apoptosis_death_model)
	{
		//apoptotic cells are colored red
		output = {"crimson", "black","darkred", "darkred"};
	}

	else if (pCell->phenotype.cycle.current_phase().code == PhysiCell_constants::necrosis_death_model)
	{
		//necrotic cells are colored brown
		output = {"peru", "black","saddlebrown", "saddlebrown"};
	}

	else if (PhysiCell::parameters.ints("simulation_mode") == 0) 
	{
		std::string drug_name = microenvironment.density_names[1];
		if (pCell->type_name == drug_name + "_sensitive")
		{
			//drug sensitive living cells are colored blue
			output = {"deepskyblue", "black", "darkblue", "darkblue"};
		} 
		else 
		{
			//drug resistant living cells are colored green
			output = {"limegreen", "black", "darkgreen", "darkgreen"};
		}
	}
	else if (PhysiCell::parameters.ints("simulation_mode") == 1) 
	{
		// // color living cells just in one color 
		// output = {"limegreen", "black", "darkgreen", "darkgreen"};

		// In case we want to color all 4 strains differently:
		// double inhibitions --> 4 cell strains 
		std::string drug1_name = microenvironment.density_names[1];
		std::string drug2_name = microenvironment.density_names[2]; 
		if (pCell->type_name == drug1_name + "_sensitive")
		{
			//cells that are sensitive to both drugs are colored blue
			output = {"deepskyblue", "black", "darkblue", "darkblue"};
		}
		else if (pCell->type_name == drug2_name + "_resistant")
		{
			// cells that are just sensitive to the first drug 
			output = {"limegreen", "black", "darkgreen", "darkgreen"};

		}
		else if (pCell->type_name == drug2_name + "_sensitive")
		{
			// cells are just sensitive to the second drug
			output = {"gold", "black", "orange", "orange"};
		}
		else 
		{
			// cells aren't sensitive to any drug
			output = {"mediumorchid", "black", "mediumpurple", "mediumpurple"};
		}
	}
	else 
	{
		// no drug simulation - living cells are colored green
		output = {"limegreen", "black", "darkgreen", "darkgreen"};
	}
	return output;

}