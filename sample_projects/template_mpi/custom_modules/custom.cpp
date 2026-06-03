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
#include "../DistPhy/DistPhy_Utils.h"
#include "../DistPhy/DistPhy_Collective.h"

#include <vector>

using namespace DistPhy::mpi;

namespace
{
	void create_cell_types_base( mpi_Environment* world, mpi_Cartesian* cart_topo )
	{
		if( parameters.ints.find_index("random_seed") != -1 )
		{
			SeedRandom(parameters.ints("random_seed"));
		}

		initialize_default_cell_definition();
		cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment );

		cell_defaults.functions.volume_update_function = standard_volume_update_function;
		cell_defaults.functions.update_velocity = standard_update_cell_velocity;
		cell_defaults.functions.update_velocity_parallel = standard_update_cell_velocity;

		cell_defaults.functions.update_migration_bias = NULL;
		cell_defaults.functions.update_phenotype = NULL;
		cell_defaults.functions.custom_cell_rule = NULL;
		cell_defaults.functions.contact_function = NULL;

		cell_defaults.functions.add_cell_basement_membrane_interactions = NULL;
		cell_defaults.functions.calculate_distance_to_membrane = NULL;

		if( world && cart_topo )
		{
			initialize_cell_definitions_from_pugixml( *world, *cart_topo );
		}
		else
		{
			initialize_cell_definitions_from_pugixml();
		}

		build_cell_definitions_maps();

		if( world && cart_topo )
		{
			setup_signal_behavior_dictionaries( *world, *cart_topo );
		}
		else
		{
			setup_signal_behavior_dictionaries();
		}

		setup_cell_rules(world, cart_topo);

		cell_defaults.functions.update_phenotype = phenotype_function;
		cell_defaults.functions.custom_cell_rule = custom_function;
		cell_defaults.functions.contact_function = contact_function;

		if( world == nullptr || IOProcessor( *world ) )
		{
			display_cell_definitions( std::cout );
		}
	}

	void seed_positions_for_definition( Cell_Definition* pCD, Microenvironment &m, mpi_Environment &world, mpi_Cartesian &cart_topo )
	{
		double Xmin = m.mesh.bounding_box[0];
		double Ymin = m.mesh.bounding_box[1];
		double Zmin = m.mesh.bounding_box[2];

		double Xmax = m.mesh.bounding_box[3];
		double Ymax = m.mesh.bounding_box[4];
		double Zmax = m.mesh.bounding_box[5];

		if( default_microenvironment_options.simulate_2D == true )
		{
			Zmin = 0.0;
			Zmax = 0.0;
		}

		double Xrange = Xmax - Xmin;
		double Yrange = Ymax - Ymin;
		double Zrange = Zmax - Zmin;

		std::vector<std::vector<double>> generated_positions_at_root;

		mpi_CellPositions cp;
		mpi_MyCells mc;

		if( world.rank == 0 )
		{
			generated_positions_at_root.reserve( parameters.ints("number_of_cells") );

			for( int n = 0 ; n < parameters.ints("number_of_cells") ; n++ )
			{
				std::vector<double> position(3,0.0);
				position[0] = Xmin + UniformRandom()*Xrange;
				position[1] = Ymin + UniformRandom()*Yrange;
				position[2] = Zmin + UniformRandom()*Zrange;

				generated_positions_at_root.push_back( position );
			}

			int start_cell_ID = Basic_Agent::get_max_ID_in_parallel();
			cp.positions_to_rank_list( generated_positions_at_root,
				m.mesh.bounding_box[0], m.mesh.bounding_box[3],
				m.mesh.bounding_box[1], m.mesh.bounding_box[4],
				m.mesh.bounding_box[2], m.mesh.bounding_box[5],
				m.mesh.dx, m.mesh.dy, m.mesh.dz,
				world, cart_topo, start_cell_ID );

			Basic_Agent::set_max_ID_in_parallel( start_cell_ID + generated_positions_at_root.size() );
		}

		distribute_cell_positions( cp, mc, world, cart_topo );

		if( IOProcessor(world) )
		{
			std::cout << "Placing cells of type " << pCD->name << " ... " << std::endl;
		}

		for( int i=0; i < mc.my_no_of_cell_IDs; i++ )
		{
			Cell* pC = create_cell( *pCD, mc.my_cell_IDs[i] );
			pC->assign_position( mc.my_cell_coords[3*i], mc.my_cell_coords[3*i+1], mc.my_cell_coords[3*i+2], world, cart_topo );
		}
	}
}

void create_cell_types( void )
{
	create_cell_types_base( nullptr, nullptr );
	return;
}

void create_cell_types( mpi_Environment &world, mpi_Cartesian &cart_topo )
{
	create_cell_types_base( &world, &cart_topo );
	return;
}

void setup_microenvironment( void )
{
	initialize_microenvironment();
	return;
}

void setup_microenvironment( mpi_Environment &world, mpi_Cartesian &cart_topo )
{
	initialize_microenvironment( world, cart_topo );
	return;
}

void setup_tissue( void )
{
	double Xmin = microenvironment.mesh.bounding_box[0];
	double Ymin = microenvironment.mesh.bounding_box[1];
	double Zmin = microenvironment.mesh.bounding_box[2];

	double Xmax = microenvironment.mesh.bounding_box[3];
	double Ymax = microenvironment.mesh.bounding_box[4];
	double Zmax = microenvironment.mesh.bounding_box[5];

	if( default_microenvironment_options.simulate_2D == true )
	{
		Zmin = 0.0;
		Zmax = 0.0;
	}

	double Xrange = Xmax - Xmin;
	double Yrange = Ymax - Ymin;
	double Zrange = Zmax - Zmin;

	Cell* pC;

	for( int k=0; k < cell_definitions_by_index.size() ; k++ )
	{
		Cell_Definition* pCD = cell_definitions_by_index[k];
		std::cout << "Placing cells of type " << pCD->name << " ... " << std::endl;
		for( int n = 0 ; n < parameters.ints("number_of_cells") ; n++ )
		{
			std::vector<double> position = {0,0,0};
			position[0] = Xmin + UniformRandom()*Xrange;
			position[1] = Ymin + UniformRandom()*Yrange;
			position[2] = Zmin + UniformRandom()*Zrange;

			pC = create_cell( *pCD );
			pC->assign_position( position );
		}
	}
	std::cout << std::endl;

	load_cells_from_pugixml();
	set_parameters_from_distributions();

	return;
}

void setup_tissue( Microenvironment &m, mpi_Environment &world, mpi_Cartesian &cart_topo )
{
	for( int k=0; k < cell_definitions_by_index.size() ; k++ )
	{
		seed_positions_for_definition( cell_definitions_by_index[k], m, world, cart_topo );
	}

	if( IOProcessor(world) )
	{
		std::cout << std::endl;
	}

	set_parameters_from_distributions();

	return;
}

std::vector<std::string> my_coloring_function( Cell* pCell )
{ return paint_by_number_cell_coloring(pCell); }

void phenotype_function( Cell* pCell, Phenotype& phenotype, double dt )
{ return; }

void custom_function( Cell* pCell, Phenotype& phenotype , double dt )
{ return; } 

void contact_function( Cell* pMe, Phenotype& phenoMe , Cell* pOther, Phenotype& phenoOther , double dt )
{ return; } 
