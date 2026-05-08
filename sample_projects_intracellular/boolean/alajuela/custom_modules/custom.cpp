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
# Copyright (c) 2015-2021, Paul Macklin and the PhysiCell Project             #
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

#include "./custom.h"

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

	cell_defaults.functions.update_migration_bias = NULL; 
	cell_defaults.functions.update_phenotype = NULL; // update_cell_and_death_parameters_O2_based; 
	cell_defaults.functions.custom_cell_rule = NULL; 
	cell_defaults.functions.contact_function = NULL; 
	
	cell_defaults.functions.add_cell_basement_membrane_interactions = NULL; 
	cell_defaults.functions.calculate_distance_to_membrane = NULL; 
	
	/*
	   This parses the cell definitions in the XML config file. 
	*/
	
	initialize_cell_definitions_from_pugixml(world, cart_topo); 

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

	setup_cell_rules( world, cart_topo); 

	/* 
	   Put any modifications to individual cell definitions here. 
	   
	   This is a good place to set custom functions. 
	*/ 




	Cell_Definition* pCD = find_cell_definition( "default");



	pCD->functions.update_velocity = standard_update_cell_velocity;
	pCD->functions.update_velocity_parallel = standard_update_cell_velocity; 
	pCD->functions.pre_update_intracellular =  pre_update_intracellular;
	pCD->functions.post_update_intracellular = post_update_intracellular;

	pCD = find_cell_definition( "blue");


	pCD->functions.update_velocity = standard_update_cell_velocity;
	pCD->functions.update_velocity_parallel = standard_update_cell_velocity;
	pCD->functions.pre_update_intracellular =  pre_update_intracellular;
	pCD->functions.post_update_intracellular = post_update_intracellular;


	pCD = find_cell_definition( "orange");


	pCD->functions.update_velocity = standard_update_cell_velocity;
	pCD->functions.update_velocity_parallel = standard_update_cell_velocity;
	pCD->functions.pre_update_intracellular =  pre_update_intracellular;
	pCD->functions.post_update_intracellular = post_update_intracellular;
	
	/*
	   This builds the map of cell definitions and summarizes the setup. 
	*/
		
	display_cell_definitions( std::cout, world, cart_topo ); 
	
	return; 
}

void setup_microenvironment( mpi_Environment &world, mpi_Cartesian &cart_topo  )
{
	initialize_microenvironment( world, cart_topo );
	
	double o2_conc;
	double drug_conc;
	
	std::vector<double> dirichlet_o2( 2 ); // number of sustrates, creates the vector of 3 positions

	o2_conc = parameters.doubles("o2_conc1");
	drug_conc = parameters.doubles("drug_conc");
		
	dirichlet_o2[0] = o2_conc;
	dirichlet_o2[1] = drug_conc;
	
	std::vector<std::vector<double>> positions;

	std::string csv_fname = parameters.strings("blood_source_file");
	positions = read_cells_positions(csv_fname, '\t', true);

	const double local_x_min = microenvironment.mesh.local_bounding_box[0];
	const double local_x_max = microenvironment.mesh.local_bounding_box[3];

	for (int i = 0; i < positions.size(); i++)
	{
		const double x = positions[i][0];
		const bool in_local_x_range =
			( x >= local_x_min ) &&
			( world.rank == world.size - 1 ? x <= local_x_max : x < local_x_max );

		if( in_local_x_range )
		{
			int local_voxel_index = microenvironment.mesh.nearest_lcl_voxel_index( positions[i] );
			microenvironment.add_dirichlet_node( local_voxel_index , dirichlet_o2 );
		}
		
		//microenvironment.set_substrate_dirichlet_activation( 1,  microenvironment.voxel_index(x,y,z),true );
		// microenvironment.add_dirichlet_node( microenvironment.voxel_index(data[0],data[1],data[2]) , dirichlet_drug );

		//microenvironment.set_substrate_dirichlet_activation( 2,  microenvironment.voxel_index(x,y,z),true );
		// microenvironment.add_dirichlet_node( microenvironment.voxel_index(positions[0],positions[1],positions[2]) , dirichlet_o2 );
	}
	
	return;
}

std::vector<std::vector<double>> create_cell_sphere_positions(double cell_radius, double sphere_radius)
{
	std::vector<std::vector<double>> cells;
	int xc=0,yc=0,zc=0;
	double x_spacing= cell_radius*sqrt(3);
	double y_spacing= cell_radius*2;
	double z_spacing= cell_radius*sqrt(3);
	
	std::vector<double> tempPoint(3,0.0);
	// std::vector<double> cylinder_center(3,0.0);
	
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


void setup_tissue(mpi_Environment &world, mpi_Cartesian &cart_topo)
{
	std::vector<std::vector<double>> positions;

	if( world.rank == 0 )
	{
		if ( parameters.bools("read_init") )
		{
			std::string csv_fname = parameters.strings("init_cells_filename");
			positions = read_cells_positions(csv_fname, ',', true);
		}
		else
		{
			double cell_radius = cell_defaults.phenotype.geometry.radius; 
			double tumor_radius =  parameters.doubles("tumor_radius");
			
			positions = create_cell_sphere_positions(cell_radius,tumor_radius);
		}
	}

	mpi_CellPositions cp;
	mpi_MyCells mc;
	int strt_cell_ID = Basic_Agent::get_max_ID_in_parallel();
	cp.positions_to_rank_list(
		positions,
		microenvironment.mesh.bounding_box[0], microenvironment.mesh.bounding_box[3],
		microenvironment.mesh.bounding_box[1], microenvironment.mesh.bounding_box[4],
		microenvironment.mesh.bounding_box[2], microenvironment.mesh.bounding_box[5],
		microenvironment.mesh.dx, microenvironment.mesh.dy, microenvironment.mesh.dz,
		world, cart_topo, strt_cell_ID
	);
	Basic_Agent::set_max_ID_in_parallel( strt_cell_ID + positions.size() );

	distribute_cell_positions( cp, mc, world, cart_topo );

	Cell* pCell = NULL; 
	for( int i = 0; i < mc.my_no_of_cell_IDs; i++ )
	{
		pCell = create_cell( get_cell_definition("default"), mc.my_cell_IDs[i] );
		pCell->assign_position(
			mc.my_cell_coords[3*i],
			mc.my_cell_coords[3*i+1],
			mc.my_cell_coords[3*i+2],
			world, cart_topo
		);
	}
	
	return;
}


void change_dirichlet_nodes ( void )
{

	double o2_conc;
	double drug_conc;
	
	std::vector<double> dirichlet_o2( 2 ); // number of sustrates, creates the vector of 3 positions

	o2_conc = parameters.doubles("o2_conc1");
	drug_conc = parameters.doubles("drug_conc");
		
	dirichlet_o2[0] = o2_conc;
	dirichlet_o2[1] = drug_conc;
	
	std::vector<std::vector<double>> positions;

	std::string csv_fname = parameters.strings("blood_source_file");
	positions = read_cells_positions(csv_fname, '\t', true);

	for (int i = 0; i < positions.size(); i++)
	{
		int x = (positions[i][0]);
		int y = ( positions[i][1]);
		int z = ( positions[i][2]);
		if (x >= microenvironment.mesh.local_bounding_box[0] && x <= microenvironment.mesh.local_bounding_box[3] )
		{
			microenvironment.set_substrate_dirichlet_activation( 1,  microenvironment.voxel_index(x,y,z),true );
		}
		
		// microenvironment.add_dirichlet_node( microenvironment.voxel_index(data[0],data[1],data[2]) , dirichlet_drug );

		//microenvironment.set_substrate_dirichlet_activation( 2,  microenvironment.voxel_index(x,y,z),true );
		// microenvironment.add_dirichlet_node( microenvironment.voxel_index(positions[0],positions[1],positions[2]) , dirichlet_o2 );
	}
	
	//std::cout << "voy" << std::endl;
	return;
}	

std::vector<std::vector<double>> read_cells_positions(std::string filename, char delimiter, bool header)
{
	// File pointer
	std::fstream fin;
	std::vector<std::vector<double>> positions;

	// Open an existing file
	fin.open(filename, std::ios::in);

	// Read the Data from the file
	// as String Vector
	std::vector<std::string> row;
	std::string line, word;

	if (header)
	{ getline(fin, line); }

	do
	{
		row.clear();

		// read an entire row and
		// store it in a string variable 'line'
		getline(fin, line);

		// used for breaking words
		std::stringstream s(line);

		while (getline(s, word, delimiter))
		{ 
			row.push_back(word); 
		}

		std::vector<double> tempPoint(3,0.0);
		tempPoint[0]= std::stof(row[0]);
		tempPoint[1]= std::stof(row[1]);
		tempPoint[2]= std::stof(row[2]);

		positions.push_back(tempPoint);
	} while (!fin.eof());

	return positions;
}


std::vector<std::string> my_coloring_function( Cell* pCell )
{ return paint_by_number_cell_coloring(pCell); }

std::vector<std::string> regular_colors( Cell* pCell )
{
	static int A_type = get_cell_definition( "green" ).type; 
	static int B_type = get_cell_definition( "blue" ).type; 
	static int C_type = get_cell_definition( "orange" ).type; 
	
	// start with flow cytometry coloring 
	
	std::vector<std::string> output = {"black" , "black" , "black" , "black"} ;

	// color live a
		
	if( pCell->type == A_type )
	{
		 output[0] = "rgb(34,139,34)";  
		 output[2] = "rgb(34,139,34)";  
	}
	
	// color live B

	if( pCell->type == B_type )
	{
		 output[0] = "rgb(0,0,255)";  
		 output[2] = "rgb(0,0,255)";  
	}
	
	// color live C

	if( pCell->type == C_type )
	{
		 output[0] = "rgb(255,140,0)";  
		 output[2] = "rgb(255,140,0)";  
	}
	
	if( pCell->phenotype.death.dead == true )
	{
		// Necrotic - Brown
		if( pCell->phenotype.cycle.current_phase().code == PhysiCell_constants::necrotic_swelling || 
			pCell->phenotype.cycle.current_phase().code == PhysiCell_constants::necrotic_lysed || 
			pCell->phenotype.cycle.current_phase().code == PhysiCell_constants::necrotic )
		{
			output[2] = "chocolate";
		}
		else
		{
			output[2] = "black"; 
		}
	}

	return output; 
}

void phenotype_function( Cell* pCell, Phenotype& phenotype, double dt )
{ return; }



void treatment_function (mpi_Environment &world, mpi_Cartesian &cart_topo) 
{
	if (PhysiCell::parameters.bools.find_index("treatment") != -1) 
	{
		int treatment_substrate_index = BioFVM::microenvironment.find_density_index(PhysiCell::parameters.strings("treatment_substrate"));

		if (PhysiCell::parameters.bools("treatment")){
		
			if (
				(((int)PhysiCell::PhysiCell_globals.current_time) % PhysiCell::parameters.ints("treatment_period")) == 0 
				&& !BioFVM::microenvironment.get_substrate_dirichlet_activation(treatment_substrate_index)
			)
			{
				if (IOProcessor(world))
					std::cout << PhysiCell::parameters.strings("treatment_substrate") << " activation at t=" << PhysiCell::PhysiCell_globals.current_time << std::endl;
				BioFVM::microenvironment.set_substrate_dirichlet_activation(treatment_substrate_index, true);	
			}

			if (
				(((int)PhysiCell::PhysiCell_globals.current_time) % PhysiCell::parameters.ints("treatment_period")) == PhysiCell::parameters.ints("treatment_duration") 
				&& BioFVM::microenvironment.get_substrate_dirichlet_activation(treatment_substrate_index)
			)
			{
				if (IOProcessor(world))
					std::cout << PhysiCell::parameters.strings("treatment_substrate") << " inactivation at t=" << PhysiCell::PhysiCell_globals.current_time << std::endl;
				BioFVM::microenvironment.set_substrate_dirichlet_activation(treatment_substrate_index, false);	
			}
			
		} else if ( BioFVM::microenvironment.get_substrate_dirichlet_activation(treatment_substrate_index) ){
			if (IOProcessor(world))
				std::cout << PhysiCell::parameters.strings("treatment_substrate") << " inactivation (NO TREATMENT) at t=" << PhysiCell::PhysiCell_globals.current_time << std::endl;
			BioFVM::microenvironment.set_substrate_dirichlet_activation(treatment_substrate_index, false);	
		}
	}
}
