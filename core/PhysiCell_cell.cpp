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

#include "PhysiCell_cell.h"
#include "PhysiCell_cell_container.h"
#include "PhysiCell_utilities.h"
#include "PhysiCell_constants.h"
#include "PhysiCell_signal_behavior.h"
#include "PhysiCell_rules.h"
#include "../BioFVM/BioFVM_vector.h"


#ifdef ADDON_PHYSIBOSS
#include "../addons/PhysiBoSS/src/maboss_intracellular.h"
#endif
#ifdef ADDON_ROADRUNNER
#include "../addons/libRoadrunner/src/librr_intracellular.h"
#endif
#ifdef ADDON_PHYSIDFBA
#include "../addons/dFBA/src/dfba_intracellular.h"
#endif

#include<limits.h>
#include <signal.h> //for segfault
#include <algorithm>
#include <iterator>

namespace PhysiCell
{

std::unordered_map<std::string,Cell_Definition*> cell_definitions_by_name;
std::unordered_map<int,Cell_Definition*> cell_definitions_by_type;
std::vector<Cell_Definition*> cell_definitions_by_index;


std::unordered_map<std::string,int> cell_definition_indices_by_name; 
std::unordered_map<int,int> cell_definition_indices_by_type;

// function pointer on how to choose cell orientation at division
// in case you want the legacy method 
std::vector<double> (*cell_division_orientation)(void) = UniformOnUnitSphere; // LegacyRandomOnUnitSphere;

Cell* standard_instantiate_cell()
{ return new Cell; }

Cell_Parameters::Cell_Parameters()
{
	o2_hypoxic_threshold = 15.0; // HIF-1alpha at half-max around 1.5-2%, and tumors often are below 2%
	o2_hypoxic_response = 8.0; // genomic / proteomic changes observed at 7-8 mmHg
	o2_hypoxic_saturation = 4.0; // maximum HIF-1alpha at 0.5% o2 (McKeown)

	o2_necrosis_threshold = 5.0;
	o2_necrosis_max = 2.5;

	o2_proliferation_threshold = 5.0; // assume no proliferation at same level as starting necrosis
	o2_proliferation_saturation = 160.0; // 5% = 38, 21% = 160 mmHg
	o2_reference = 160.0; // assume all was measured in normoxic 21% o2

	pReference_live_phenotype = NULL; // reference live (usually physioxic) phenotype

	// necrosis parameters

	max_necrosis_rate = 1.0 / (6.0 * 60.0); // assume cells survive 6 hours in very low oxygen
	necrosis_type = PhysiCell_constants::deterministic_necrosis;;

	return;
}

void Cell_Parameters::pack(std::vector<char>& snd_buffer, int& len_buffer, int &position){
	/* Packing members of class 'Cell_Parameters' except for pReference_live_phenotype pointer */
    /* There are 9 doubles and 1 int in this class except for the pointer */

    /*=====================================================================*/
    /* pReference_live_phenotype * IS NOT PACKED 													 */
    /*=====================================================================*/

   // Determine buffer size and pack using helper functions
	pack_buff(this->o2_hypoxic_threshold,      snd_buffer, len_buffer, position);
	pack_buff(this->o2_hypoxic_response,       snd_buffer, len_buffer, position);
	pack_buff(this->o2_hypoxic_saturation,     snd_buffer, len_buffer, position);
	pack_buff(this->o2_proliferation_saturation, snd_buffer, len_buffer, position);
	pack_buff(this->o2_proliferation_threshold,  snd_buffer, len_buffer, position);
	pack_buff(this->o2_reference,              snd_buffer, len_buffer, position);
	pack_buff(this->o2_necrosis_threshold,     snd_buffer, len_buffer, position);
	pack_buff(this->o2_necrosis_max,           snd_buffer, len_buffer, position);
	pack_buff(this->max_necrosis_rate,         snd_buffer, len_buffer, position);
	pack_buff(this->necrosis_type,             snd_buffer, len_buffer, position);
}

void Cell_Parameters::unpack(std::vector<char>& rcv_buffer, int& len_buffer, int &position){
	 // Parameters
    unpack_buff(this->o2_hypoxic_threshold, rcv_buffer, len_buffer, position);
    unpack_buff(this->o2_hypoxic_response, rcv_buffer, len_buffer, position);
    unpack_buff(this->o2_hypoxic_saturation, rcv_buffer, len_buffer, position);
    unpack_buff(this->o2_proliferation_saturation, rcv_buffer, len_buffer, position);
    unpack_buff(this->o2_proliferation_threshold, rcv_buffer, len_buffer, position);
    unpack_buff(this->o2_reference, rcv_buffer, len_buffer, position);
    unpack_buff(this->o2_necrosis_threshold, rcv_buffer, len_buffer, position);
    unpack_buff(this->o2_necrosis_max, rcv_buffer, len_buffer, position);
    unpack_buff(this->max_necrosis_rate, rcv_buffer, len_buffer, position);
    unpack_buff(this->necrosis_type, rcv_buffer, len_buffer, position);
}

Cell_Definition::Cell_Definition()
{
	// set the microenvironment pointer
	pMicroenvironment = NULL;
	if( BioFVM::get_default_microenvironment() != NULL )
	{ pMicroenvironment = BioFVM::get_default_microenvironment(); }

	// set up the default parameters
		// the default Cell_Parameters constructor should take care of this

	type = 0;
	name = "unnamed";

	is_movable = true;

	parameters.pReference_live_phenotype = &phenotype;

	// set up the default custom data
		// the default Custom_Cell_Data constructor should take care of this

	// set up the default functions
	functions.instantiate_cell = NULL;
	functions.volume_update_function = NULL; // standard_volume_update_function;
	functions.update_migration_bias = NULL;

	functions.update_phenotype = NULL;
	functions.custom_cell_rule = NULL;

	functions.update_velocity = NULL; 										// standard_update_cell_velocity;

	/*-------------------------------------------------------------------------------------*/
	/* For completion making update_velocity_parallel = NULL here as well 								 */
	/*-------------------------------------------------------------------------------------*/

	functions.update_velocity_parallel = NULL;
	functions.add_cell_basement_membrane_interactions = NULL; 
	functions.calculate_distance_to_membrane = NULL; 

	// bug fix July 31 2023
	functions.plot_agent_SVG = standard_agent_SVG;
	functions.plot_agent_legend = standard_agent_legend;
	// bug fix July 31 2023

	functions.set_orientation = NULL;

	// new March 2022 : make sure Cell_Interactions and Cell_Transformations 
	// 					are appropriately sized. Same on motiltiy. 
	phenotype.cell_interactions.sync_to_cell_definitions(); 
	phenotype.cell_transformations.sync_to_cell_definitions(); 
	phenotype.cycle.asymmetric_division.sync_to_cell_definitions();
	phenotype.motility.sync_to_current_microenvironment(); 
	phenotype.mechanics.sync_to_cell_definitions(); 

	cell_definitions_by_index.push_back( this ); 
	return;
}

Cell_Definition::Cell_Definition( Cell_Definition& cd )
{
	// set the microenvironment pointer
	pMicroenvironment = cd.pMicroenvironment;

	// set up the default parameters
		// the default Cell_Parameters constructor should take care of this

	type = cd.type;
	name = cd.name;

	parameters = cd.parameters;
	custom_data = cd.custom_data;
	functions = cd.functions;
	phenotype = cd.phenotype;

	// this is the whole reason we need ot make a copy constructor
	parameters.pReference_live_phenotype = &phenotype;

	cell_definitions_by_index.push_back( this );

	return;
}

Cell_Definition& Cell_Definition::operator=( const Cell_Definition& cd )
{
	// set the microenvironment pointer
	pMicroenvironment = cd.pMicroenvironment;

	// set up the default parameters
		// the default Cell_Parameters constructor should take care of this

	type = cd.type;
	name = cd.name;

	parameters = cd.parameters;
	custom_data = cd.custom_data;
	functions = cd.functions;
	phenotype = cd.phenotype;

	// this is the whole reason we need ot make a copy constructor
	parameters.pReference_live_phenotype = &phenotype;

	return *this;
}


Cell_Definition cell_defaults;

Cell_State::Cell_State()
{
	neighbors.resize(0);
	orientation.resize( 3 , 0.0 );

	simple_pressure = 0.0;
	mechanics_neighbor_candidates = 0;
	mechanics_neighbor_interactions = 0;
	
	attached_cells.clear();
	spring_attachments.clear(); 

	number_of_nuclei = 1; 

	//damage = 0.0; 
	total_attack_time = 0.0; 

	contact_with_basement_membrane = false; 

	return;
}

void Cell_State::pack(vector<char>& snd_buffer, int& len_buffer, int& position) {
	pack_buff( this->orientation, snd_buffer, len_buffer, position);
	pack_buff( this->simple_pressure, snd_buffer, len_buffer, position);
	pack_buff( this->number_of_nuclei, snd_buffer, len_buffer, position);
	pack_buff( this->total_attack_time, snd_buffer, len_buffer, position);
	pack_buff( this->contact_with_basement_membrane, snd_buffer, len_buffer, position);
}

void Cell_State::unpack(vector<char>& rcv_buffer, int& len_buffer, int& position) {
	unpack_buff(this->orientation, rcv_buffer, len_buffer, position);
    unpack_buff(this->simple_pressure, rcv_buffer, len_buffer, position);
    unpack_buff(this->number_of_nuclei, rcv_buffer, len_buffer, position);
    unpack_buff(this->total_attack_time, rcv_buffer, len_buffer, position);
    unpack_buff(this->contact_with_basement_membrane, rcv_buffer, len_buffer, position);
}
void Cell::update_motility_vector( double dt_ )
{
	if( phenotype.motility.is_motile == false )
	{
		phenotype.motility.motility_vector.assign( 3, 0.0 );
		return;
	}

	if( UniformRandom() < dt_ / phenotype.motility.persistence_time || phenotype.motility.persistence_time < dt_ )
	{
		std::vector<double> randvec(3,0.0);
		if( phenotype.motility.restrict_to_2D == true )
		{ randvec = UniformOnUnitCircle(); }
		else
		{ randvec = UniformOnUnitSphere(); }

		// if the update_bias_vector function is set, use it
		if( functions.update_migration_bias )
		{
			functions.update_migration_bias( this,phenotype,dt_ );
		}

		phenotype.motility.motility_vector = phenotype.motility.migration_bias_direction; // motiltiy = bias_vector
		phenotype.motility.motility_vector *= phenotype.motility.migration_bias; // motility = bias*bias_vector

		double one_minus_bias = 1.0 - phenotype.motility.migration_bias;

		axpy( &(phenotype.motility.motility_vector), one_minus_bias, randvec ); // motility = (1-bias)*randvec + bias*bias_vector

		normalize( &(phenotype.motility.motility_vector) );

		phenotype.motility.motility_vector *= phenotype.motility.migration_speed;
	}
	return;
}

void Cell::advance_bundled_phenotype_functions( double dt_ )
{
	standard_cell_transformations( this,this->phenotype,dt_ );

	if( PhysiCell_settings.rules_enabled )
	{ apply_ruleset( this ); }
	if( get_single_signal(this,"necrotic") > 0.5 )
	{
		double rupture = this->phenotype.volume.rupture_volume; 
		double volume = this->phenotype.volume.total; 
		if( volume > rupture )
		{
			std::cout << this->phenotype.volume.total << " vs " << this->phenotype.volume.rupture_volume << 
			" dead: " << get_single_signal( this, "dead") << 	std::endl; 
			std::cout << this->phenotype.cycle.current_phase_index() << " " 
			<< this->phenotype.cycle.pCycle_Model->name << std::endl; 
		}

	}

	// call the custom code to update the phenotype
	if( functions.update_phenotype )
	{	functions.update_phenotype( this , phenotype , dt_ ); }

	// update volume
	if( functions.volume_update_function )
	{
		functions.volume_update_function(this,phenotype,dt_);

		// The following line is needed in every volume
		// regulation method (it sets BioFVM total_volume)

		set_total_volume( phenotype.volume.total );
	}

	// update geometry
	phenotype.geometry.update( this, phenotype, dt_ );

	// update integrity 
	phenotype.cell_integrity.advance_damage( dt_ );

	// check for new death events
	if( phenotype.death.check_for_death( dt_ ) == true )
	{
		// if so, change the cycle model to the current death model
		phenotype.cycle.sync_to_cycle_model( phenotype.death.current_model() );

		// also, turn off motility.

		phenotype.motility.is_motile = false;
		phenotype.motility.motility_vector.assign( 3, 0.0 );
		functions.update_migration_bias = NULL;

		// turn off secretion, and reduce uptake by a factor of 10
		phenotype.secretion.set_all_secretion_to_zero();
		phenotype.secretion.scale_all_uptake_by_factor( 0.10 );

		// make sure to run the death entry function
		if( phenotype.cycle.current_phase().entry_function )
		{
			phenotype.cycle.current_phase().entry_function( this, phenotype, dt_ );
		}
	}

	// advance cycle model (for both cell cycle and death cycle models)
	phenotype.cycle.advance_cycle( this, phenotype, dt_ );
	if( phenotype.flagged_for_removal )
	{
		flag_for_removal();
		phenotype.flagged_for_removal = false;
	}
	if( phenotype.flagged_for_division )
	{
		flag_for_division();
		phenotype.flagged_for_division = false;
	}

	return;
}

void Cell::advance_bundled_phenotype_functions( double dt_ , mpi_Environment &world, mpi_Cartesian &cart_topo )
{
	// New March 2022
	// perform transformations 
	standard_cell_transformations( this,this->phenotype,dt_ ); //no mpi needed

	// New March 2023 in Version 1.12.0 
	// call the rules-based code to update the phenotype 
	if( PhysiCell_settings.rules_enabled )
	{ apply_ruleset( this ); } // no mpi needed
	if( get_single_signal(this,"necrotic") > 0.5 )
	{
		double rupture = this->phenotype.volume.rupture_volume; 
		double volume = this->phenotype.volume.total; 
		if( volume > rupture )
		{
			std::cout << this->phenotype.volume.total << " vs " << this->phenotype.volume.rupture_volume << 
			" dead: " << get_single_signal( this, "dead") << 	std::endl; 
			std::cout << this->phenotype.cycle.current_phase_index() << " " 
			<< this->phenotype.cycle.pCycle_Model->name << std::endl; 
		}

	}
	// call the custom code to update the phenotype
	if( functions.update_phenotype_parallel )
	{	functions.update_phenotype_parallel( this , phenotype , dt_ , world, cart_topo ); }
	else if( functions.update_phenotype )
	{	
		functions.update_phenotype( this , phenotype , dt_); 
		set_total_volume( phenotype.volume.total ); 
	}

	// update volume
	if( functions.volume_update_function )
	{
		functions.volume_update_function(this,phenotype,dt_);

		// The following line is needed in every volume
		// regulation method (it sets BioFVM total_volume)

		set_total_volume( phenotype.volume.total );
	}

	// update geometry
	phenotype.geometry.update( this, phenotype, dt_ );

	// update integrity 
	phenotype.cell_integrity.advance_damage( dt_ ); // no mpi needed

	// check for new death events
	if( phenotype.death.check_for_death( dt_ ) == true )
	{
		// if so, change the cycle model to the current death model
		phenotype.cycle.sync_to_cycle_model( phenotype.death.current_model() );

		// also, turn off motility.

		phenotype.motility.is_motile = false;
		phenotype.motility.motility_vector.assign( 3, 0.0 );
		functions.update_migration_bias = NULL;

		// turn off secretion, and reduce uptake by a factor of 10
		phenotype.secretion.set_all_secretion_to_zero();
		phenotype.secretion.scale_all_uptake_by_factor( 0.10 );

		// make sure to run the death entry function
		if( phenotype.cycle.current_phase().entry_function )
		{
			phenotype.cycle.current_phase().entry_function( this, phenotype, dt_ );
		}
	}

	// advance cycle model (for both cell cycle and death cycle models)
	phenotype.cycle.advance_cycle( this, phenotype, dt_ );
	if( phenotype.flagged_for_removal )
	{
		flag_for_removal();
		phenotype.flagged_for_removal = false;
	}
	if( phenotype.flagged_for_division )
	{
		flag_for_division();
		phenotype.flagged_for_division = false;
	}

	return;
}

Cell::Cell()
{
	// use the cell defaults;

	type = cell_defaults.type;
	type_name = cell_defaults.name;

	custom_data = cell_defaults.custom_data;
	parameters = cell_defaults.parameters;
	functions = cell_defaults.functions;

	phenotype = cell_defaults.phenotype;

	phenotype.molecular.sync_to_cell( this );

	// cell state should be fine by the default constructor

	current_mechanics_voxel_index=-1;

	updated_current_mechanics_voxel_index = 0;

	is_movable = true;
	is_out_of_domain = false;

	/*-----------------------------------*/
	/* Added by Gaurav Saxena to say that*/
	/* when cell created, both = 0			 */
	/*-----------------------------------*/

	crossed_to_left_subdomain = false;
	crossed_to_right_subdomain = false;


	displacement.resize(3,0.0); // state?

	assign_orientation();
	container = NULL;

	set_total_volume( phenotype.volume.total );

	return;
}

/*--------------------------------------------------------------------------------------------------*/
/* Parallel version of Cell class Constructor which takes as input the cell id in parallel settings */
/* It must call explicitly a new constructor of Basic_Agent which takes as input this cell id				*/
/*--------------------------------------------------------------------------------------------------*/

Cell::Cell(int p_ID):Basic_Agent(p_ID) //----> Correct syntax for calling parametrized base class constructor
{
	// use the cell defaults;

	//Basic_Agent(p_ID);---> this syntax is incorrect to call Base class parametrized constructor
	type = cell_defaults.type;
	type_name = cell_defaults.name;

	custom_data = cell_defaults.custom_data;
	parameters 	= cell_defaults.parameters;
	functions 	= cell_defaults.functions;

	phenotype 	= cell_defaults.phenotype;

	phenotype.molecular.sync_to_cell( this );

	// cell state should be fine by the default constructor

	current_mechanics_voxel_index =-1;
	updated_current_mechanics_voxel_index = 0;

	/*-----------------------------------------------------------------------------------------------*/
	/* Added by Gaurav Saxena to say that when cell created, both = 0. (I had forgotten to add this  */
	/* condition in this Cell::Cell(int p_ID):Basic_Agent(p_ID) constructor (See below) 						 */
	/*-----------------------------------------------------------------------------------------------*/

	crossed_to_left_subdomain = false;
	crossed_to_right_subdomain = false;

	is_movable 				= true;
	is_out_of_domain 	= false;
	displacement.resize(3,0.0); // state?

	assign_orientation();				//Just assigns a random unit vector to the cell.
	container = NULL;

	/* Gaurav Saxena added this line as this was in v1.7 */
	set_total_volume( phenotype.volume.total );

	return;
}

Cell::~Cell()
{
//	std::cout << std::endl << "=====-----------------------------=====" << std::endl; 
//	std::cout << "\tcell destructor " << this << " " << type_name << " at " << position << std::endl;
//	std::cout << "\t\tattached cells: " << this->state.attached_cells.size() << std::endl << std::endl; 
	if (all_cells != nullptr) {
	auto result = std::find( std::begin(*all_cells),std::end(*all_cells),this );
	if( result != std::end(*all_cells) )
	{
		std::cout << "Warning: Cell was never removed from data structure " << std::endl ; 
		std::cout << "I am of type " << this->type << " at " << this->position << std::endl; 

		int temp_index = -1; 
		bool found = false; 
		for( int n= 0 ; n < (*all_cells).size() ; n++ )
		{
			std::cout << this << " vs " << (*all_cells)[n] << std::endl; 
			if( (*all_cells)[n] == this )
			{ found = true; temp_index = n; } 
		}
		
		if( found )
		{
			// release any attached cells (as of 1.7.2 release)
			this->remove_all_attached_cells();
			// 1.11.0
			this->remove_all_spring_attachments();

			// new Dec 5, 2024: 
			this->remove_self_from_all_neighbors(); 
			// released internalized substrates (as of 1.5.x releases)
			this->release_internalized_substrates(); 

			// performance goal: don't delete in the middle -- very expensive reallocation
			// alternative: copy last element to index position, then shrink vector by 1 at the end O(constant)

			// move last item to index location  
			(*all_cells)[ (*all_cells).size()-1 ]->index=temp_index;
			(*all_cells)[temp_index] = (*all_cells)[ (*all_cells).size()-1 ];
			// shrink the vector
			(*all_cells).pop_back();	
			
			// deregister agent in from the agent container
			this->get_container()->remove_agent(this);
		}
	}
	}
	
	
	return; 
}


void Cell::flag_for_division( void )
{
	get_container()->flag_cell_for_division( this );
	return;
}

void Cell::flag_for_removal( void )
{
	get_container()->flag_cell_for_removal( this );
	return;
}

void Cell::start_death( int death_model_index )
{
	// set the death data struture to the indicated death model
	phenotype.death.trigger_death( death_model_index );
	// change the cycle model to the current death model
	phenotype.cycle.sync_to_cycle_model( phenotype.death.current_model() );

	// turn off secretion, and reduce uptake by a factor of 10
	phenotype.secretion.set_all_secretion_to_zero();
	phenotype.secretion.scale_all_uptake_by_factor( 0.10 );

	// turn off motility.
	phenotype.motility.is_motile = false;
	phenotype.motility.motility_vector.assign( 3, 0.0 );
	functions.update_migration_bias = NULL;

	// make sure to run the death entry function
	if( phenotype.cycle.current_phase().entry_function )
	{
		phenotype.cycle.current_phase().entry_function( this, phenotype, 0.0 );
	}

	return;
}

void Cell::assign_orientation()
{
	state.orientation.resize(3,0.0);
	if( functions.set_orientation != NULL )
	{
		functions.set_orientation(this, phenotype, 0.0 );
	}
	else
	{
		//assign a random unit vector
		double theta= UniformRandom()*6.28318530717959; //rand*2*pi
		double z= 2* UniformRandom()-1;
		double temp= sqrt(1-z*z);
		state.orientation[0]= temp * cos(theta);
		state.orientation[1]= temp * sin(theta);
		state.orientation[2]= z;
	}

	return;
}

Cell* Cell::divide( )
{
	// phenotype.flagged_for_division = false;
	// phenotype.flagged_for_removal = false;
	
	// make sure ot remove adhesions 
	remove_all_attached_cells();
	remove_all_spring_attachments(); 

	// version 1.10.3: 
	// conserved quantitites in custom data aer divided in half
	// so that each daughter cell gets half of the original ;
	for( int nn = 0 ; nn < custom_data.variables.size() ; nn++ )
	{
		if( custom_data.variables[nn].conserved_quantity == true )
		{ custom_data.variables[nn].value *= 0.5; }
	}
	for( int nn = 0 ; nn < custom_data.vector_variables.size() ; nn++ )
	{
		if( custom_data.vector_variables[nn].conserved_quantity == true )
		{ custom_data.vector_variables[nn].value *= 0.5; }
	}

	Cell* child = create_cell(functions.instantiate_cell);
	child->copy_data( this );
	child->copy_function_pointers(this);
	child->parameters = parameters;

	// evenly divide internalized substrates
	// if these are not actively tracked, they are zero anyway
	*internalized_substrates *= 0.5;
	*(child->internalized_substrates) = *internalized_substrates ;

	// The following is already performed by create_cell(). JULY 2017 ***
	// child->register_microenvironment( get_microenvironment() );

	// randomly place the new agent close to me, accounting for orientation and
	// polarity (if assigned)
	
	// May 30, 2020: 
	// Set cell_division_orientation = LegacyRandomOnUnitSphere to 
	// reproduce this code 
	/*
	double temp_angle = 6.28318530717959*UniformRandom();
	double temp_phi = 3.1415926535897932384626433832795*UniformRandom();

	double radius= phenotype.geometry.radius;
	std::vector<double> rand_vec (3, 0.0);

	rand_vec[0]= cos( temp_angle ) * sin( temp_phi );
	rand_vec[1]= sin( temp_angle ) * sin( temp_phi );
	rand_vec[2]= cos( temp_phi );
	rand_vec = rand_vec- phenotype.geometry.polarity*(rand_vec[0]*state.orientation[0]+
		rand_vec[1]*state.orientation[1]+rand_vec[2]*state.orientation[2])*state.orientation;

	if( norm(rand_vec) < 1e-16 )
	{
		std::cout<<"************ERROR********************"<<std::endl;
	}
	normalize( &rand_vec );
	rand_vec/= norm(rand_vec); */
	
	std::vector<double> rand_vec = cell_division_orientation(); 
	rand_vec = rand_vec- phenotype.geometry.polarity*(rand_vec[0]*state.orientation[0]+ 
		rand_vec[1]*state.orientation[1]+rand_vec[2]*state.orientation[2])*state.orientation;	
	rand_vec *= phenotype.geometry.radius;
	
	
	child->assign_position(position[0] + rand_vec[0],
						   position[1] + rand_vec[1],
						   position[2] + rand_vec[2]);
	//change my position to keep the center of mass intact 
	//and then see if I need to update my voxel index
	
	

	static double negative_one_half = -0.5;
	axpy( &position, negative_one_half , rand_vec );// position = position - 0.5*rand_vec;

	//If this cell has been moved outside of the boundaries, mark it as such.
	//(If the child cell is outside of the boundaries, that has been taken care of in the assign_position function.)
	if( !get_container()->underlying_mesh.is_position_valid(position[0], position[1], position[2])){
		is_out_of_domain = true;
		is_active = false;
		is_movable = false;
	}

	update_voxel_in_container();
	phenotype.volume.divide();
	child->phenotype.volume.divide();
	child->set_total_volume(child->phenotype.volume.total);
	set_total_volume(phenotype.volume.total);

	// child->set_phenotype( phenotype );
	child->phenotype = phenotype;
	
	if (child->phenotype.intracellular){
        child->phenotype.intracellular->start();
		child->phenotype.intracellular->inherit(this);
	}
	
// #ifdef ADDON_PHYSIDFBA
// 	child->fba_model = this->fba_model;
// #endif
	// changes for new phenotyp March 2022
	//state.damage = 0.0; 
	state.total_attack_time = 0; 
	//child->state.damage = 0.0; 
	child->state.total_attack_time = 0.0; 

	if( this->functions.cell_division_function )
	{ this->functions.cell_division_function( this, child); }

	return child;
}

/*-------------------------------------------------------------------------*/
/* Parallel version of the divide() function above												 */
/*-------------------------------------------------------------------------*/

Cell* Cell::divide(int p_ID, mpi_Environment &world, mpi_Cartesian &cart_topo)
{
	// phenotype.flagged_for_division = false;
	// phenotype.flagged_for_removal = false;
	
	// make sure ot remove adhesions 
	remove_all_attached_cells();
	remove_all_spring_attachments(); 

	for( int nn = 0 ; nn < custom_data.variables.size() ; nn++ )
	{
		if( custom_data.variables[nn].conserved_quantity == true )
		{ custom_data.variables[nn].value *= 0.5; }
	}
	for( int nn = 0 ; nn < custom_data.vector_variables.size() ; nn++ )
	{
		if( custom_data.vector_variables[nn].conserved_quantity == true )
		{ custom_data.vector_variables[nn].value *= 0.5; }
	}
	
	/* EXPERIMENTING HERE BY MAKING SYNTAX EXACTLY SAME AS IN SETUP_TISSUE() */
	/* AND IN UNPACKING FUNCTION 																						 */
	/* IF DOES NOT WORK THEN CHANGE BACK TO Cell* child = create_cell(p_ID); */	
	
	Cell* child = create_cell(functions.instantiate_cell,p_ID);					//Calls new version of create_cell function.
	//Cell* child = create_cell(p_ID);
	
	child->copy_data( this );
	child->copy_function_pointers(this);
	child->parameters = parameters;

	// evenly divide internalized substrates
	// if these are not actively tracked, they are zero anyway
	*internalized_substrates *= 0.5;
	*(child->internalized_substrates) = *internalized_substrates ;

	// The following is already performed by create_cell(). JULY 2017 ***
	// child->register_microenvironment( get_microenvironment() );

	// randomly place the new agent close to me, accounting for orientation and
	// polarity (if assigned)

	// double temp_angle = 6.28318530717959*UniformRandom();
// 	double temp_phi = 3.1415926535897932384626433832795*UniformRandom();
// 
// 	double radius= phenotype.geometry.radius;
// 	std::vector<double> rand_vec (3, 0.0);
// 
// 	rand_vec[0]= cos( temp_angle ) * sin( temp_phi );
// 	rand_vec[1]= sin( temp_angle ) * sin( temp_phi );
// 	rand_vec[2]= cos( temp_phi );
// 	rand_vec = rand_vec- phenotype.geometry.polarity*(rand_vec[0]*state.orientation[0]+
// 		rand_vec[1]*state.orientation[1]+rand_vec[2]*state.orientation[2])*state.orientation;
// 
// 	if( norm(rand_vec) < 1e-16 )
// 	{
// 		std::cout<<"************ERROR********************"<<std::endl;
// 	}
// 	normalize( &rand_vec );

	// rand_vec/= norm(rand_vec);<---- What is the use of it ? (Question: Gaurav Saxena)

	//std::cout<<" Daughter X Position="<<position[0]+ 0.5 * radius*rand_vec[0]<<std::endl;

	/*------------------------------------------------------------------------------------------*/
	/* The problem is that sometimes the new-born daughter cell CAN push out of the sub-domain  */
	/* boundary. Hence a rule needs to be imposed that the Daughter cell cannot cross the sub-	*/
	/* domain boundary.																																					*/
	/* (1) Get subdomain boundary limits. (2) See if x coordinate of child is outside subdomain */
	/* boundary limits. (3) If yes, pull it back, if no, do nothing 														*/
	/*------------------------------------------------------------------------------------------*/

	std::vector<double> child_position(3);
	std::vector<double> temp_subdomain_bndry = get_container()->underlying_mesh.get_subdomain_limits();
	double eps = 1e-5; // 1e-16 was the previous value, possibly causing errors like x - eps = x
	
	std::vector<double> rand_vec = cell_division_orientation(); 
	rand_vec = rand_vec- phenotype.geometry.polarity * 
						 (rand_vec[0]*state.orientation[0] + rand_vec[1]*state.orientation[1] + rand_vec[2]*state.orientation[2]) * 
						 state.orientation;	
	rand_vec *= phenotype.geometry.radius;
	

	child_position[0] = position[0] + rand_vec[0];
	child_position[1] = position[1] + rand_vec[1];
	child_position[2] = position[2] + rand_vec[2];

	/*---------------------------------------------------------------------------------*/
	/* The following 2 functions will adjust the position of the child x-coordinate IF */
	/* it crosses the sub-domains but NOT if it crosses the DOMAIN boundary 					 */
	/*---------------------------------------------------------------------------------*/

	if(get_container()->underlying_mesh.has_it_crossed_to_left_subdomain(child_position[0], world))
	{
		child_position[0] = temp_subdomain_bndry[0] + eps;
		//std::cout<<"Corrected Daughter X Position="<<child_position[0]<<std::endl;
	}
	if(get_container()->underlying_mesh.has_it_crossed_to_right_subdomain(child_position[0], world))
	{
		child_position[0] = temp_subdomain_bndry[3] - eps;
		//std::cout<<"Corrected Daughter X Position="<<child_position[0]<<std::endl;
	}
  
	/*-----------------------------------------------------------------------------------*/
  /* previously we were assigning position without checking if local subdomain crossed */
  /* Now we are checking in above. 																										 */
  /* BUT there assign_position(...) calls update_voxel_index(...) which calls 			   */
  /* correct_position_within_subdomain(...) <--- this function will NOT do anything to */
  /* the child position because we have already handled it above 											 */
  /* The parent is being handled explicitly below if it crosses to next sub-domain		 */
  /* */
	/*-----------------------------------------------------------------------------------*/


	child->assign_position(child_position[0], child_position[1], child_position[2], world, cart_topo);

	//Just trying if this changes anything, other wise DELETE next line
	//IT DOES REPAIR THE CYCLE_RATE PROBLEM BUT CELLS STILL BLOW UP
	//update_monitor_variables(child);
	//change my position to keep the center of mass intact and then see if I need to update my voxel index


	/*------------------------------------------------------------------------------------*/
	/* SIDENOTE: In naxpy() below we are passing &position. position is declared as 			*/
	/* std::vector<double> position. Now, "position" does not mean the address of the 		*/
	/* vector i.e. "position" is NOT a pointer. When we pass &position, it doesn't mean		*/
	/* that we are passing a pointer to a pointer. It JUST means we are passing the 			*/
	/* address of the vector - which is caught in std::vector<double> *y	i.e. a pointer	*/
	/* a type "vector". This is a big difference between arrays and vectors. 							*/
	/*------------------------------------------------------------------------------------*/

	/*------------------------------------------------------------------------------------*/
	/* TICKET: depending on the solution of the ticket, I may have to change the value of	*/
	/* negative_one_half = +0.5 and then call naxpy(). 																		*/
	/*------------------------------------------------------------------------------------*/

	static double negative_one_half = -0.5;
	axpy( &position, negative_one_half , rand_vec );// position = position - 0.5*rand_vec;

	// position[0] -= 0.5*radius*rand_vec[0];
	// position[1] -= 0.5*radius*rand_vec[1];
	// position[2] -= 0.5*radius*rand_vec[2];

	//If this cell has been moved outside of the boundaries, mark it as such.
	//(If the child cell is outside of the boundaries, that has been taken care of in the assign_position function.)

	if( !get_container()->underlying_mesh.is_position_valid(position[0], position[1], position[2]))
	{
		is_out_of_domain = true;
		is_active = false;
		is_movable = false;
		//std::cout<<"Cell ID"<<this->ID<<" is out of domain"<<std::endl; 
	}

	/*-----------------------------------------------------------------------------------*/
	/* Now the parent cell has also moved after division, so there is a possibility that */
	/* it may not be respecting sub-domain boundaries - therefore we need to call 			 */
	/* correct_position_within_subdomain() for the parent cell.													 */
	/*-----------------------------------------------------------------------------------*/

	get_container()->underlying_mesh.correct_position_within_subdomain(position, world, cart_topo);

	update_voxel_in_container(world, cart_topo);			//Though world, cart_topo are not needed in this, just to distinguish
	phenotype.volume.divide(); 												//Don't think anything to parallelize
	child->phenotype.volume.divide();
	child->set_total_volume(child->phenotype.volume.total);
	set_total_volume(phenotype.volume.total);

	// child->set_phenotype( phenotype );
	child->phenotype = phenotype;
		
	if (child->phenotype.intracellular) {
		child->phenotype.intracellular->start();
		child->phenotype.intracellular->inherit(this);
	}

  	
// #ifdef ADDON_PHYSIDFBA
// 	child->fba_model = this->fba_model;
// #endif

	// changes for new phenotyp March 2022
	//state.damage = 0.0; 
	state.total_attack_time = 0; 
	//child->state.damage = 0.0; 
	child->state.total_attack_time = 0.0;

	if( this->functions.cell_division_function )
	{ this->functions.cell_division_function( this, child); }
	
	return child;
}


bool Cell::assign_position(std::vector<double> new_position)
{
	return assign_position(new_position[0], new_position[1], new_position[2]);
}

void Cell::set_previous_velocity(double xV, double yV, double zV)
{
	previous_velocity[0] = xV;
	previous_velocity[1] = yV;
	previous_velocity[2] = zV;

	return;
}

bool Cell::assign_position(double x, double y, double z)
{
	
	position[0]=x;
	position[1]=y;
	position[2]=z;

	// update microenvironment current voxel index
	update_voxel_index();
	// update current_mechanics_voxel_index
	current_mechanics_voxel_index= get_container()->underlying_mesh.nearest_voxel_index( position );
	
	/*---------------------------------------------------------------------*/
	/* VALIDATION FIRST: Check if position is valid before registering.    */
	/* This prevents registering out-of-domain cells that would be clamped */
	/* by nearest_voxel_index(), leaving them in a valid voxel while also  */
	/* marked is_out_of_domain (inconsistent state).                       */
	/*---------------------------------------------------------------------*/
	
	if( !get_container()->underlying_mesh.is_position_valid(x,y,z) )
	{
		is_out_of_domain = true;
		is_active = false;
		is_movable = false;
		return false;
	}

	// Only register valid cells in the agent grid
	get_container()->register_agent(this);

	// Since it is most likely our first position, we update the max_cell_interactive_distance_in_voxel
	// which was not initialized at cell creation
	if( get_container()->max_cell_interactive_distance_in_voxel[get_current_mechanics_voxel_index()] < 
		phenotype.geometry.radius * phenotype.mechanics.relative_maximum_adhesion_distance )
	{
		// get_container()->max_cell_interactive_distance_in_voxel[get_current_mechanics_voxel_index()]= phenotype.geometry.radius*parameters.max_interaction_distance_factor;
		get_container()->max_cell_interactive_distance_in_voxel[get_current_mechanics_voxel_index()] = phenotype.geometry.radius
			* phenotype.mechanics.relative_maximum_adhesion_distance;
	}

	return true;
}

/*-------------------------------------------------------*/
/* Parallel version of assign_position() function				 */
/*-------------------------------------------------------*/

bool Cell::assign_position(double x, double y, double z, mpi_Environment &world, mpi_Cartesian &cart_topo)
{
	//std::cout << "Cell moves from " << "(" << position[0] << "," << position[1] << "," << position[2] << ") to (" << x << "," << y << "," << z << ")"  << std::endl;
	position[0]=x;
	position[1]=y;
	position[2]=z;

	// update microenvironment current voxel index - this becomes the local voxel index in parallel settings

	/*------------------------------------------------------------------------------------*/
	/* update_voxel_index(world, cart_topo) is in BioFVM_basic_agent.cpp									*/
	/* First it calls is_position_valid() for the "Diffusion" mesh. 											*/
	/* If position invalid then sets out_of_domain=true, current_voxel_index=-1 & returns	*/
	/* Otherwise calls nearest_voxel_local_index()<---new function (GS).									*/
	/*------------------------------------------------------------------------------------*/
	if( !get_container()->underlying_mesh.is_position_valid(x,y,z) )
	{
		is_out_of_domain = true;
		is_active = false;
		is_movable = false;
		return false;
	}

	update_voxel_index(world, cart_topo);

	//update current_mechanics_voxel_index
	//current_mechanics_voxel_index= get_container()->underlying_mesh.nearest_voxel_index( position);
	//Changed the call above to a parallel call.

	current_mechanics_voxel_index= get_container()->underlying_mesh.nearest_voxel_local_index( position, world, cart_topo);
	
	//std::cout<<"("<<position[0]<<","<<position[1]<<","<<position[2]<<")"<<" current_mechanics_voxel_index="<<current_mechanics_voxel_index<<std::endl;

	/*--------------------------------------------------------------------------------------------*/
	/* VALIDATION FIRST: Check if position is valid before registering in agent_grid.            */
	/* This prevents registering out-of-domain cells that would be clamped by nearest_voxel_*() */
	/* Note: update_voxel_index() above already called is_position_valid() internally, but we    */
	/* re-validate here explicitly to catch any position corrections and handle them correctly.  */
	/*--------------------------------------------------------------------------------------------*/

	
	// Only register valid cells in the agent grid
	get_container()->register_agent(this);

	// Since it is most likely our first position, we update the max_cell_interactive_distance_in_voxel
	// which was not initialized at cell creation
	if( get_container()->max_cell_interactive_distance_in_voxel[get_current_mechanics_voxel_index()] < 
		phenotype.geometry.radius * phenotype.mechanics.relative_maximum_adhesion_distance )
	{
		get_container()->max_cell_interactive_distance_in_voxel[get_current_mechanics_voxel_index()] = phenotype.geometry.radius
			* phenotype.mechanics.relative_maximum_adhesion_distance;
	}

	return true;
}

bool Cell::assign_position(std::vector<double> new_position, mpi_Environment &world, mpi_Cartesian &cart_topo) {
	return this->assign_position(new_position[0], new_position[1], new_position[2], world, cart_topo);
}

void Cell::set_total_volume(double volume)
{
	Basic_Agent::set_total_volume(volume);

	// If the new volume is significantly different than the
	// current total volume, adjust all the sub-volumes
	// proportionally.

	// if( fabs( phenotype.volume.total - volume ) < 1e-16 )
	if( fabs( phenotype.volume.total - volume ) > 1e-16 )
	{
		double ratio= volume/ phenotype.volume.total;
		phenotype.volume.multiply_by_ratio(ratio);
	}

	phenotype.geometry.update( this, phenotype, 0.0 );
	// phenotype.update_radius();
	//if( get_container()->max_cell_interactive_distance_in_voxel[get_current_mechanics_voxel_index()] <
	//	phenotype.geometry.radius * parameters.max_interaction_distance_factor )
	// Here the current mechanics voxel index may not be initialized, when position is still unknown. 
	if (get_current_mechanics_voxel_index() >= 0)
  {
		if( get_container()->max_cell_interactive_distance_in_voxel[get_current_mechanics_voxel_index()] <
				phenotype.geometry.radius * phenotype.mechanics.relative_maximum_adhesion_distance )
		{
			// get_container()->max_cell_interactive_distance_in_voxel[get_current_mechanics_voxel_index()]= phenotype.geometry.radius*parameters.max_interaction_distance_factor;
			get_container()->max_cell_interactive_distance_in_voxel[get_current_mechanics_voxel_index()] = phenotype.geometry.radius
			* phenotype.mechanics.relative_maximum_adhesion_distance;
		}
	}

	return;
}


void Cell::set_target_volume( double new_volume )
{

	// this function will keep the prior ratios (from targets)

	// first compute the actual raw totals on all these things
	double old_target_solid = phenotype.volume.target_solid_nuclear +
		phenotype.volume.target_solid_cytoplasmic;
	double old_target_total = old_target_solid / ( 1.0 - phenotype.volume.target_fluid_fraction );
	double old_target_fluid = phenotype.volume.target_fluid_fraction * old_target_total;

	// next whats the relative new size?
	double ratio = new_volume / (1e-16 + old_target_total );

	// scale the target solid cyto and target solid nuclear by this ratio
	phenotype.volume.target_solid_cytoplasmic *= ratio;
	phenotype.volume.target_solid_nuclear *= ratio;

	return;
}

void Cell::set_target_radius(double new_radius ) 
{
	static double four_thirds_pi =  4.188790204786391;

	// calculate the new target volume
	double new_volume = four_thirds_pi;
	new_volume *= new_radius;
	new_volume *= new_radius;
	new_volume *= new_radius;

	// now call the set_target_volume funciton
	this->set_target_volume( new_volume );
	return;
}

void Cell::set_radius(double new_radius )
{
	static double four_thirds_pi =  4.188790204786391;

	// calculate the new target volume
	double new_volume = four_thirds_pi;
	new_volume *= new_radius;
	new_volume *= new_radius;
	new_volume *= new_radius;

	this->set_total_volume( new_volume );
	return;
}



double& Cell::get_total_volume(void)
{
	static bool I_warned_you = false;
	if( I_warned_you == false )
	{
		std::cout << "Warning! Do not use " << __FUNCTION__ << "!" << std::endl
			<< "Use (some_cell).phenotype.volume.total instead!" << std::endl;
		I_warned_you = true;
	}
	return phenotype.volume.total;
}

void Cell::turn_off_reactions(double dt)
{
	is_active = false;

	for(int i=0;i< phenotype.secretion.uptake_rates.size();i++)
	{
		phenotype.secretion.uptake_rates[i] = 0.0;
		phenotype.secretion.secretion_rates[i] = 0.0;
		phenotype.secretion.net_export_rates[i] = 0.0;

	}
	set_internal_uptake_constants(dt);

	return;
}

Cell_Container * Cell::get_container()
{
	if(container == NULL)
	{
		container = (Cell_Container *)get_microenvironment()->agent_container;
	}

	return container;
}

void Cell::die() 
{
	delete_cell(this);
	return;
}

/*--------------------------------------------------------------------------------------*/
/* A function needed in the distributed parallel scenario which removes a cell that has */
/* crossed to the neighbouring sub-domain. 																							*/
/*--------------------------------------------------------------------------------------*/

void Cell::remove_crossed_cell()
{
	 vaporize_teleported_cell(this->index);
}

void Cell::update_position( double dt ) //done
{
	// BioFVM Basic_Agent::update_position(dt) returns without doing anything.
	// So we remove this to avoid any future surprises.
	//
	// Basic_Agent::update_position(dt);

	// use Adams-Bashforth
	static double d1;
	static double d2;
	static bool constants_defined = false;
	if( constants_defined == false )
	{
		d1 = dt;
		d1 *= 1.5;
		d2 = dt;
		d2 *= -0.5;
		constants_defined = true;
	}

	// new AUgust 2017
	if( default_microenvironment_options.simulate_2D == true )
	{ velocity[2] = 0.0; }

	std::vector<double> old_position(position);
	axpy( &position , d1 , velocity );
	axpy( &position , d2 , previous_velocity );
	// overwrite previous_velocity for future use
	// if(sqrt(dist(old_position, position))>3* phenotype.geometry.radius)
		// std::cout<<sqrt(dist(old_position, position))<<"old_position: "<<old_position<<", new position: "<< position<<", velocity: "<<velocity<<", previous_velocity: "<< previous_velocity<<std::endl;

	previous_velocity = velocity;

	velocity[0]=0; velocity[1]=0; velocity[2]=0;
	if(get_container()->underlying_mesh.is_position_valid(position[0],position[1],position[2]))
	{
		updated_current_mechanics_voxel_index=get_container()->underlying_mesh.nearest_voxel_index( position );
	}
	else
	{
		updated_current_mechanics_voxel_index=-1;

		is_out_of_domain = true;
		is_active = false;
		is_movable = false;
	}
	return;
}

/*----------------------------------------------*/
/* Parallel Implementation of update_position() */
/*----------------------------------------------*/

void Cell::update_position( double dt, mpi_Environment &world, mpi_Cartesian &cart_topo )
{
	/* use Adams-Bashforth */
	//if (std::isnan(position[0])) std::cout <<"Alert in update psoition!! " << ID << std::endl; 

	static double d1;
	static double d2;
	static bool constants_defined = false;

	if( constants_defined == false )
	{
		d1  = dt;
		d1 *= 1.5;
		d2  = dt;
		d2 *= -0.5;
		constants_defined = true;
	}


	if(default_microenvironment_options.simulate_2D == true )
		 velocity[2] = 0.0;

	
	std::vector<double> old_position(position);
	axpy( &position , d1 , velocity );
	axpy( &position , d2 , previous_velocity );

	//if (std::isnan(position[0])) std::cout <<"Alert in update psoition!! " << ID << std::endl; 
	/* overwrite previous_velocity for future use */

	previous_velocity = velocity;

	velocity[0]=0;
	velocity[1]=0;
	velocity[2]=0;

	crossed_to_left_subdomain  = get_container()->underlying_mesh.has_it_crossed_to_left_subdomain(position[0], world); 
	crossed_to_right_subdomain = get_container()->underlying_mesh.has_it_crossed_to_right_subdomain(position[0], world);

	/*------------------------------------------------------------------------------------------------*/
	/* According to Miguel/Arnau suggestion, have added a "subdomain_crossing_counter" which tells us */
	/* how many times a cell has crossed the sub-domain, this field exists for each cell 							*/
	/* The packing/unpacking of it is already taken care of in the custom_variables packing/unpacking */
	/* section of pack()/unpack() functions. It will automatically appear in XML										  */
	/*------------------------------------------------------------------------------------------------*/

	//int index = this->custom_data.find_variable_index( "subdomain_crossing_counter");


	bool crossed = crossed_to_left_subdomain || crossed_to_right_subdomain;

	if(crossed == true)
	{
		//this->custom_data[index]++ ; //Increase the value of subdomain_crossing_counter
		//std::cout<<std::endl;
		//std::cout<<"Rank="<<world.rank<<" Cell ID="<<ID<<" Old Position:("<<old_position[0]<<","<<old_position[1]<<","<<old_position[2]<<")"<<std::endl;
		//std::cout<<"Rank="<<world.rank<<" Cell ID="<<ID<<" New Position:("<<position[0]<<","<<position[1]<<","<<position[2]<<")"<<std::endl;
		//std::cout<<"Crossed to left ="<<crossed_to_left_subdomain<<" Crossed to right ="<<crossed_to_right_subdomain<<std::endl;
		
		/* There is a possibility that (1) Cell might cross to left/right sub-domain 		 */
		/* AND (2) because of the new Y/Z coordinate, it goes out of the domain			 		 */
		/* In this case, we need to set: crossed_to_left/right_subdomain as false 	 		 */
		/* We do NOT want such a cell to be packed/unpacked as it has gone out of domain */
		
		if(get_container()->underlying_mesh.is_position_valid(position[0],position[1],position[2]) == false)
		{
			updated_current_mechanics_voxel_index=-1;
			is_out_of_domain = true;
			is_active = false;
			is_movable = false;
			
			crossed_to_left_subdomain = false;
			crossed_to_right_subdomain = false; 
		}	
	}


	/*------------------------------------------------------------------------------*/
	/* If it has crossed to left OR right subdomain (from a non-boundary processes) */
	/* then we don't need to check if it is out-of-domain as it CANNOT be. 					*/
	/* If BOTH _left and _right above are false, then there is a possibility 				*/
	/* that process is boundary process and cell has gone out of physical boundary	*/
	/* In this case, the code below will execute 																		*/
	/* In case the cell is still within sub-domain, a new voxel index IS recomputed */
	/* IMPORTANT: DO NOT call the correct_position_within_subdomain() function from */
	/* here. THAT function is only meant to pull back the parent cell to the  			*/
	/* to the sub-domain IN CASE the parent cell crosses boundary due to center	of  */
	/* mass preservation.																														*/
	/*------------------------------------------------------------------------------*/

	if(crossed == false)
	{
		if(get_container()->underlying_mesh.is_position_valid(position[0],position[1],position[2]))
		{
			/* Need to call parallel version of nearest_voxel_index i.e. nearest_voxel_local_index */

			updated_current_mechanics_voxel_index=get_container()->underlying_mesh.nearest_voxel_local_index( position, world, cart_topo );
		}
		else
		{
			updated_current_mechanics_voxel_index=-1;
			is_out_of_domain = true;
			is_active = false;
			is_movable = false;
		}
	}
	return;
}

int Cell::get_current_mechanics_voxel_index()
{
	return current_mechanics_voxel_index;
}

void Cell::update_voxel_in_container()
{
	// call the method from BioFVM_basic_agent to update microenvironment's voxel index
	update_voxel_index();
	// int temp_current_voxel_index;
	// Check to see if we need to remove agents that are pushed out of boundary
	// if(!get_container()->underlying_mesh.is_position_valid(position[0],position[1],position[2]))

	if(updated_current_mechanics_voxel_index==-1)// updated_current_mechanics_voxel_index is updated in update_position
	{
		// check if this agent has a valid voxel index, if so, remove it from previous voxel
		if( get_current_mechanics_voxel_index() >= 0)
		{
			{get_container()->remove_agent_from_voxel(this, get_current_mechanics_voxel_index());}
		}
		{get_container()->add_agent_to_outer_voxel(this);}
		// std::cout<<"cell out of boundary..."<< __LINE__<<" "<<ID<<std::endl;
		current_mechanics_voxel_index=-1;
		is_out_of_domain=true;
		is_active=false;
		return;
	}

	// temp_current_voxel_index= get_current_mechanics_voxel_index();
	// updated_current_mechanics_voxel_index=get_container()->underlying_mesh.nearest_voxel_index( position );

	// update mesh indices (if needed)
	if(updated_current_mechanics_voxel_index!= get_current_mechanics_voxel_index())
	{
		{
			container->remove_agent_from_voxel(this, get_current_mechanics_voxel_index());
			container->add_agent_to_voxel(this, updated_current_mechanics_voxel_index);
		}
		current_mechanics_voxel_index=updated_current_mechanics_voxel_index;
	}

	return;
}

/*------------------------------------------------------------------*/
/* Parallel version of function: update_voxel_in_container() above  */
/*------------------------------------------------------------------*/

void Cell::update_voxel_in_container(mpi_Environment &world, mpi_Cartesian &cart_topo)
{

	update_voxel_index(world, cart_topo);							 //This is in BioFVM_basic_agent.cpp

	// Check to see if we need to remove agents that are pushed out of boundary

	if(updated_current_mechanics_voxel_index==-1)			//This is updated in update_position() and -1 indicates outside physical boundary
	{
		if( get_current_mechanics_voxel_index() >= 0)		//If agent had valid voxel index, need to remove agent from that voxel
			  get_container()->remove_agent_from_voxel(this, get_current_mechanics_voxel_index());


		get_container()->add_agent_to_outer_voxel(this); //std::cout<<"cell out of boundary..."<< __LINE__<<" "<<ID<<std::endl;
		current_mechanics_voxel_index=-1;
		is_out_of_domain=true;
		is_active=false;
		return;
	}

	// update mesh indices (if needed)
	if(updated_current_mechanics_voxel_index!= get_current_mechanics_voxel_index())
	{
			container->remove_agent_from_voxel(this, get_current_mechanics_voxel_index());
			container->add_agent_to_voxel(this, updated_current_mechanics_voxel_index);
			current_mechanics_voxel_index=updated_current_mechanics_voxel_index;
	}
	return;
}

void Cell::copy_data(Cell* copy_me)
{
	// phenotype=copyMe->phenotype; //it is taken care in set_phenotype
	type = copy_me->type;
	type_name = copy_me->type_name;

	custom_data = copy_me->custom_data;
	parameters = copy_me->parameters;

	velocity = copy_me->velocity;
	// expected_phenotype = copy_me-> expected_phenotype; //it is taken care in set_phenotype
	cell_source_sink_solver_temp1 = std::vector<double>(copy_me->cell_source_sink_solver_temp1);
	cell_source_sink_solver_temp2 = std::vector<double>(copy_me->cell_source_sink_solver_temp2);

	return;
}

void Cell::copy_function_pointers(Cell* copy_me)
{
	functions = copy_me->functions;
	return;
}

void Cell::add_potentials(Cell* other_agent)
{
	//if( this->ID == other_agent->ID )
	if( this == other_agent )
	{ return; }

	// 12 uniform neighbors at a close packing distance, after dividing out all constants
	static double simple_pressure_scale = 0.027288820670331; // 12 * (1 - sqrt(pi/(2*sqrt(3))))^2
	// 9.820170012151277; // 12 * ( 1 - sqrt(2*pi/sqrt(3)))^2

	double distance = 0;
	for( int i = 0 ; i < 3 ; i++ )
	{
		displacement[i] = position[i] - (*other_agent).position[i];
		distance += displacement[i] * displacement[i];
	}
	// Make sure that the distance is not zero

	distance = std::max(sqrt(distance), 0.00001);

	//Repulsive
	double R = phenotype.geometry.radius+ (*other_agent).phenotype.geometry.radius;

	//double RN = phenotype.geometry.nuclear_radius + (*other_agent).phenotype.geometry.nuclear_radius;
	double temp_r, c;
	if( distance > R )
	{
		temp_r=0;
	}
	// else if( distance < RN )
	// {
		// double M = 1.0;
		// c = 1.0 - RN/R;
		// c *= c;
		// c -= M;
		// temp_r = ( c*distance/RN  + M  );
	// }
	else
	{
		// temp_r = 1 - distance/R;
		temp_r = -distance; // -d
		temp_r /= R; // -d/R
		temp_r += 1.0; // 1-d/R
		temp_r *= temp_r; // (1-d/R)^2

		// add the relative pressure contribution
		state.simple_pressure += ( temp_r / simple_pressure_scale ); // New July 2017
	}

	// August 2017 - back to the original if both have same coefficient
	double effective_repulsion = sqrt( phenotype.mechanics.cell_cell_repulsion_strength * other_agent->phenotype.mechanics.cell_cell_repulsion_strength );
	temp_r *= effective_repulsion;

	// temp_r *= phenotype.mechanics.cell_cell_repulsion_strength; // original
	//////////////////////////////////////////////////////////////////

	// Adhesive
	//double max_interactive_distance = parameters.max_interaction_distance_factor * phenotype.geometry.radius +
	//	(*other_agent).parameters.max_interaction_distance_factor * (*other_agent).phenotype.geometry.radius;

	double max_interactive_distance = phenotype.mechanics.relative_maximum_adhesion_distance * phenotype.geometry.radius +
		(*other_agent).phenotype.mechanics.relative_maximum_adhesion_distance * (*other_agent).phenotype.geometry.radius;

	if(distance < max_interactive_distance )
	{
		// double temp_a = 1 - distance/max_interactive_distance;
		double temp_a = -distance; // -d
		temp_a /= max_interactive_distance; // -d/S
		temp_a += 1.0; // 1 - d/S
		temp_a *= temp_a; // (1-d/S)^2
		// temp_a *= phenotype.mechanics.cell_cell_adhesion_strength; // original

		// August 2017 - back to the original if both have same coefficient
		// May 2022 - back to oriinal if both affinities are 1
		int ii = find_cell_definition_index( this->type ); 
		int jj = find_cell_definition_index( other_agent->type ); 

		double adhesion_ii = phenotype.mechanics.cell_cell_adhesion_strength * phenotype.mechanics.cell_adhesion_affinities[jj]; 
		double adhesion_jj = other_agent->phenotype.mechanics.cell_cell_adhesion_strength * other_agent->phenotype.mechanics.cell_adhesion_affinities[ii]; 

		//double effective_adhesion = sqrt( phenotype.mechanics.cell_cell_adhesion_strength * other_agent->phenotype.mechanics.cell_cell_adhesion_strength );
		double effective_adhesion = sqrt( adhesion_ii*adhesion_jj );
		temp_a *= effective_adhesion;

		temp_r -= temp_a;
		state.mechanics_neighbor_interactions += 1;
		state.neighbors.push_back(other_agent);
	}
	/////////////////////////////////////////////////////////////////
	if( fabs(temp_r) < 1e-16 )
	{ return; }
	temp_r /= distance;
	// for( int i = 0 ; i < 3 ; i++ )
	// {
	//	velocity[i] += displacement[i] * temp_r;
	// }
	axpy( &velocity , temp_r , displacement );

	return;
}

/*----------------------------------------------------------------------------------*/
/* Parallel version of the function above: add_potentials 													*/
/* Instead of taking a Cell* argument, it takes an argument of type Moore_Cell_Info */
/* defined in PhysiCell_cell_container.h 																						*/
/*----------------------------------------------------------------------------------*/

void Cell::add_potentials(Moore_Cell_Info &other_agent, mpi_Environment &world, mpi_Cartesian &cart_topo)
{
	if( this->ID == other_agent.ID )
	//if( this == other_agent ) --> Commented out (Gaurav Saxena) because "other_agent" is NOT a pointer
	{
		return;
	}

	// 12 uniform neighbors at a close packing distance, after dividing out all constants
	static double simple_pressure_scale = 0.027288820670331; // 12 * (1 - sqrt(pi/(2*sqrt(3))))^2
																														// 9.820170012151277; // 12 * ( 1 - sqrt(2*pi/sqrt(3)))^2

	double distance = 0;
	for( int i = 0 ; i < 3 ; i++ )
	{
		displacement[i] = position[i] - other_agent.position[i];
		distance += displacement[i] * displacement[i];
	}
	// Make sure that the distance is not zero

	distance = std::max(sqrt(distance), 0.00001);

	//Repulsive
	double R = phenotype.geometry.radius+ other_agent.radius;

	double RN = phenotype.geometry.nuclear_radius + other_agent.nuclear_radius;
	double temp_r, c;
	if( distance > R )
	{
		temp_r=0;
	}
	// else if( distance < RN )
	// {
		// double M = 1.0;
		// c = 1.0 - RN/R;
		// c *= c;
		// c -= M;
		// temp_r = ( c*distance/RN  + M  );
	// }
	else
	{
		// temp_r = 1 - distance/R;
		temp_r = -distance; // -d
		temp_r /= R; // -d/R
		temp_r += 1.0; // 1-d/R
		temp_r *= temp_r; // (1-d/R)^2

		// add the relative pressure contribution
		state.simple_pressure += ( temp_r / simple_pressure_scale ); // New July 2017
	}

	// August 2017 - back to the original if both have same coefficient
	double effective_repulsion = sqrt( phenotype.mechanics.cell_cell_repulsion_strength * other_agent.cell_cell_repulsion_strength );
	temp_r *= effective_repulsion;

	// temp_r *= phenotype.mechanics.cell_cell_repulsion_strength; // original
	//////////////////////////////////////////////////////////////////

	// Adhesive
	//double max_interactive_distance = parameters.max_interaction_distance_factor * phenotype.geometry.radius +
	//	(*other_agent).parameters.max_interaction_distance_factor * (*other_agent).phenotype.geometry.radius;

	double max_interactive_distance = phenotype.mechanics.relative_maximum_adhesion_distance *
																		phenotype.geometry.radius +
		     														other_agent.relative_maximum_adhesion_distance *
		     														other_agent.radius;

	if(distance < max_interactive_distance )
	{
		// double temp_a = 1 - distance/max_interactive_distance;
		double temp_a = -distance; // -d
		temp_a /= max_interactive_distance; // -d/S
		temp_a += 1.0; // 1 - d/S
		temp_a *= temp_a; // (1-d/S)^2
		// temp_a *= phenotype.mechanics.cell_cell_adhesion_strength; // original

		// August 2017 - back to the original if both have same coefficient
		int ii = find_cell_definition_index( this->type ); 
		int jj = find_cell_definition_index( other_agent.type );

		double adhesion_ii = phenotype.mechanics.cell_cell_adhesion_strength * phenotype.mechanics.cell_adhesion_affinities[jj]; 
		double adhesion_jj = other_agent.cell_cell_adhesion_strength * phenotype.mechanics.cell_adhesion_affinities[ii]; 


		double effective_adhesion = sqrt( adhesion_ii*adhesion_jj ); 
		temp_a *= effective_adhesion;

		temp_r -= temp_a;
		state.mechanics_neighbor_interactions += 1;
	}
	/////////////////////////////////////////////////////////////////
	if( fabs(temp_r) < 1e-16 )
	{ return; }
	temp_r /= distance;
	// for( int i = 0 ; i < 3 ; i++ )
	// {
	//	velocity[i] += displacement[i] * temp_r;
	// }
	axpy( &velocity , temp_r , displacement );
	
	//state.neighbors.push_back(other_agent); // new 1.8.0 <-- THIS IS NOT POSSIBLE as "other_agent" is NOT a pointer
																						// more importantly "other_agent" contains information of cell that is in
																						//another process, so even if we get that cell pointer to this process
																						//it would be an invalid pointer in this process. We need to discuss this
																						//with Randy Heiland. 

	return;
}

Cell* create_cell( Cell* (*custom_instantiate)())
{
	Cell* pNew; 
	
	if (custom_instantiate) {
		pNew = custom_instantiate();
	} else {
		pNew = standard_instantiate_cell();
	}
	
	(*all_cells).push_back( pNew ); 
	pNew->index=(*all_cells).size()-1;
	
	// new usability enhancements in May 2017 
	
	if( BioFVM::get_default_microenvironment() )
	{
		pNew->register_microenvironment( BioFVM::get_default_microenvironment() );
	}

	// All the phenotype and other data structures are already set 
	// by virtue of the default Cell constructor. 
	
	pNew->set_total_volume( pNew->phenotype.volume.total ); 
	
	return pNew; 
}

/*---------------------------------------------------------------------*/
/* New function for parallel environment: Cell* create_cell(int p_ID)  */
/* Very similar body to the serial Cell* create_cell( void ) above 	 	 */
/* It will call a new Constructor pNew = new Cell(p_ID);							 */
/* This new Constructor will explicitly call a new Constructor of		 	 */
/* the Basic_Agent class : Basic_Agent(int p_ID)											 */
/*---------------------------------------------------------------------*/

Cell* create_cell(  int p_ID )
{
	Cell* pNew;

	pNew = new Cell(p_ID);
	(*all_cells).push_back( pNew );
	pNew->index=(*all_cells).size()-1;


	if( BioFVM::get_default_microenvironment() )
	{
		pNew->register_microenvironment( BioFVM::get_default_microenvironment() );
	}
	
	pNew->set_total_volume( pNew->phenotype.volume.total );

	return pNew;
}

Cell* create_cell( Cell* (*custom_instantiate)(),  int p_ID )
{
	Cell* pNew;

	if (custom_instantiate) {
		pNew = custom_instantiate();
	} else {
		pNew = standard_instantiate_cell();
	}

	// Set the ID on the instantiated cell (don't create a new one and discard the custom cell)
	pNew->ID = p_ID;
	
	(*all_cells).push_back( pNew );
	pNew->index=(*all_cells).size()-1;


	if( BioFVM::get_default_microenvironment() )
	{
		pNew->register_microenvironment( BioFVM::get_default_microenvironment() );
	}
	
	pNew->set_total_volume( pNew->phenotype.volume.total );

	return pNew;
}

// in the future, I might swap this with create_cell():
// In that "create_cell()" uses "create_cell( cell_defaults )"
Cell* create_cell( Cell_Definition& cd ) 
{
	Cell* pNew = create_cell(cd.functions.instantiate_cell);

	// use the cell defaults;
	pNew->type = cd.type;
	pNew->type_name = cd.name;

	pNew->custom_data = cd.custom_data;
	pNew->parameters = cd.parameters;
	pNew->functions = cd.functions;

	if( &cd != &cell_defaults )
	{ pNew->phenotype = cd.phenotype; }
	if (pNew->phenotype.intracellular)
		pNew->phenotype.intracellular->start();

	pNew->is_movable = cd.is_movable;
	pNew->is_out_of_domain = false;
	pNew->displacement.resize(3,0.0); // state?

	pNew->assign_orientation();

	pNew->set_total_volume( pNew->phenotype.volume.total );

	return pNew;
}

/*------------------------------------------------------------------------------------------*/
/* New function for parallel environment: Cell* create_cell(Cell_Definition &cd, int p_ID)  */
/* Very similar body to the serial Cell* create_cell( Cell_Definition &cd ) above 	 	 			*/
/* Within it, it will call the new function Cell *create_cell(int pID) which in turn 			  */
/* calls a new Constructor pNew = new Cell(p_ID);							 															*/
/* This new Constructor will explicitly call a new Constructor of		 	 											*/
/* the Basic_Agent class : Basic_Agent(int p_ID)											 											*/
/*------------------------------------------------------------------------------------------*/

Cell* create_cell( Cell_Definition& cd, int pID )
{
	Cell* pNew = create_cell(cd.functions.instantiate_cell);
	pNew->ID = pID;

	// use the cell defaults;
	pNew->type = cd.type;
	pNew->type_name = cd.name;

	pNew->custom_data = cd.custom_data;
	pNew->parameters = cd.parameters;
	pNew->functions = cd.functions;

	if( &cd != &cell_defaults )
	{ pNew->phenotype = cd.phenotype; }
	if (pNew->phenotype.intracellular)
		pNew->phenotype.intracellular->start();

	pNew->is_movable = cd.is_movable;
	pNew->is_out_of_domain = false;
	pNew->displacement.resize(3,0.0); // state?

	pNew->assign_orientation();

	pNew->set_total_volume( pNew->phenotype.volume.total );

	return pNew;
}

void Cell::convert_to_cell_definition( Cell_Definition& cd )
{
	Volume cell_volume = phenotype.volume;
	Geometry cell_geometry = phenotype.geometry;
	Molecular cell_molecular = phenotype.molecular;
	Custom_Cell_Data cell_custom_data = custom_data;
	// use the cell defaults;
	type = cd.type;
	type_name = cd.name;

	custom_data = cd.custom_data;
	parameters = cd.parameters;
	functions = cd.functions;

	phenotype = cd.phenotype;
	phenotype.volume = cell_volume; 
	phenotype.volume.target_solid_nuclear = cd.phenotype.volume.target_solid_nuclear;
	phenotype.volume.target_fluid_fraction = cd.phenotype.volume.target_fluid_fraction;
	phenotype.volume.target_cytoplasmic_to_nuclear_ratio = cd.phenotype.volume.target_cytoplasmic_to_nuclear_ratio;
	phenotype.volume.relative_rupture_volume = cd.phenotype.volume.relative_rupture_volume;

	phenotype.geometry = cell_geometry; // leave the geometry alone
	phenotype.molecular.internalized_total_substrates = cell_molecular.internalized_total_substrates;

	for( int nn = 0 ; nn < custom_data.variables.size() ; nn++ )
	{
		if( custom_data.variables[nn].conserved_quantity == true )
		{ custom_data.variables[nn].value = cell_custom_data.variables[nn].value; }
	}
	for( int nn = 0 ; nn < custom_data.vector_variables.size() ; nn++ )
	{
		if( custom_data.vector_variables[nn].conserved_quantity == true )
		{ custom_data.vector_variables[nn].value = cell_custom_data.vector_variables[nn].value; }
	}

	if (get_current_mechanics_voxel_index() >= 0)
    {
        if( get_container()->max_cell_interactive_distance_in_voxel[get_current_mechanics_voxel_index()] < 
            phenotype.geometry.radius * phenotype.mechanics.relative_maximum_adhesion_distance )
        {
            get_container()->max_cell_interactive_distance_in_voxel[get_current_mechanics_voxel_index()] = phenotype.geometry.radius
                * phenotype.mechanics.relative_maximum_adhesion_distance;
        }
	}
	//assign_orientation();

	//set_total_volume( phenotype.volume.total );

	return;
}

void delete_cell( int index )
{
//	std::cout << __FUNCTION__ << " " << (*all_cells)[index] 
//	<< " " << (*all_cells)[index]->type_name << std::endl; 
	
	Cell* pDeleteMe = (*all_cells)[index]; 
	
	// release any attached cells (as of 1.7.2 release)
	pDeleteMe->remove_all_attached_cells();
	// 1.11.0 
	pDeleteMe->remove_all_spring_attachments(); 

	pDeleteMe->remove_self_from_all_neighbors(); 

	// released internalized substrates (as of 1.5.x releases)
	pDeleteMe->release_internalized_substrates(); 
			

	// performance goal: don't delete in the middle -- very expensive reallocation
	// alternative: copy last element to index position, then shrink vector by 1 at the end O(constant)

	// move last item to index location  
	(*all_cells)[ (*all_cells).size()-1 ]->index=index;
	(*all_cells)[index] = (*all_cells)[ (*all_cells).size()-1 ];
	// shrink the vector
	(*all_cells).pop_back();	
	
	// deregister agent in from the agent container
	pDeleteMe->get_container()->remove_agent(pDeleteMe);
	// de-allocate (delete) the cell; 
	delete pDeleteMe; 


	return; 
}

void delete_cell_original( int index ) // before June 11, 2020
{
//	std::cout << __FUNCTION__ << " " << (*all_cells)[index] 
//	<< " " << (*all_cells)[index]->type_name << std::endl; 
	
	// release any attached cells (as of 1.7.2 release)
	(*all_cells)[index]->remove_all_attached_cells();
	// 1.11.0
	(*all_cells)[index]->remove_all_spring_attachments(); 
	
	// released internalized substrates (as of 1.5.x releases)
	(*all_cells)[index]->release_internalized_substrates(); 
	
	// deregister agent in from the agent container
	(*all_cells)[index]->get_container()->remove_agent((*all_cells)[index]);
	// de-allocate (delete) the cell; 
	delete (*all_cells)[index]; 

	// performance goal: don't delete in the middle -- very expensive reallocation
	// alternative: copy last element to index position, then shrink vector by 1 at the end O(constant)

	// move last item to index location  
	(*all_cells)[ (*all_cells).size()-1 ]->index=index;
	(*all_cells)[index] = (*all_cells)[ (*all_cells).size()-1 ];
	// shrink the vector
	(*all_cells).pop_back();	
	return; 
}

void delete_cell( Cell* pDelete )
{
	delete_cell(pDelete->index);
	return;
}

/*--------------------------------------------------------------------------------------------*/
/* The function vaporize_teleported_cell(int index) is only for distributed parallel scenario */
/* It is exactly like delete(int) function above except for it doesn't call 								  */
/* 'release_internalized_substrates()' as the cell does not die but only leaves the current		*/
/* sub-domain and enters a neighbouring sub-domain 																						*/
/*--------------------------------------------------------------------------------------------*/

void vaporize_teleported_cell(int index)
{
	
	/*------------------------------------------------------------------------------*/
	/* Changing the function as program is crashing in add_potentials(...) function */
	/* due to accessing illegal memory. 																						*/
	/*------------------------------------------------------------------------------*/

	Cell* pDeleteMe = (*all_cells)[index];
	
	pDeleteMe->remove_all_attached_cells(); 			//newly added statement 
	// 1.11.0
	pDeleteMe->remove_all_spring_attachments(); 
	// Ensure neighbors drop stale pointers before this cell is removed from all_cells.
	pDeleteMe->remove_self_from_all_neighbors();
	
	//NO NEED to release substrates as the cell is ONLY moving to another sub-domain and is not dying
	
	// performance goal: don't delete in the middle -- very expensive reallocation
	// alternative: copy last element to index position, then shrink vector by 1 at the end O(constant)
	// move last item to index location

	(*all_cells)[ (*all_cells).size()-1 ]->index=index;
	(*all_cells)[index] = (*all_cells)[ (*all_cells).size()-1 ];

	// shrink the vector
	(*all_cells).pop_back();
	
	// deregister agent in from the agent container
	/* get the container of the cell and remove the basic agent (basically cell) from it */
	pDeleteMe->get_container()->remove_agent(pDeleteMe); //moved this AFTER the above statements
	
	delete pDeleteMe; 	//This should be the last statement (see "void delete_cell( int index )" function above)

}

bool is_neighbor_voxel(Cell* pCell, std::vector<double> my_voxel_center, std::vector<double> other_voxel_center, int other_voxel_index)
{
	double max_interactive_distance = pCell->phenotype.mechanics.relative_maximum_adhesion_distance * pCell->phenotype.geometry.radius
		+ pCell->get_container()->max_cell_interactive_distance_in_voxel[other_voxel_index];

	int comparing_dimension = -1, comparing_dimension2 = -1;
	if(my_voxel_center[0] == other_voxel_center[0] && my_voxel_center[1] == other_voxel_center[1])
	{
		comparing_dimension = 2;
	}
	else if(my_voxel_center[0] == other_voxel_center[0] && my_voxel_center[2] == other_voxel_center[2])
	{
		comparing_dimension = 1;
	}
	else if(my_voxel_center[1] == other_voxel_center[1] && my_voxel_center[2] == other_voxel_center[2])
	{
		comparing_dimension = 0;
	}

	if(comparing_dimension != -1)
	{ //then it is an immediate neighbor (through side faces)
		double surface_coord= 0.5*(my_voxel_center[comparing_dimension] + other_voxel_center[comparing_dimension]);
		if(std::fabs(pCell->position[comparing_dimension] - surface_coord) > max_interactive_distance)
		{ return false; }
		return true;
	}
	comparing_dimension=-1;

	if(my_voxel_center[0] == other_voxel_center[0])
	{
		comparing_dimension = 1; comparing_dimension2 = 2;
	}
	else if(my_voxel_center[1] == other_voxel_center[1])
	{
		comparing_dimension=0; comparing_dimension2 = 2;
	}
	else if(my_voxel_center[2] == other_voxel_center[2])
	{
		comparing_dimension = 0; comparing_dimension2=1;
	}
	if(comparing_dimension != -1)
	{
		double line_coord1= 0.5*(my_voxel_center[comparing_dimension] + other_voxel_center[comparing_dimension]);
		double line_coord2= 0.5*(my_voxel_center[comparing_dimension2] + other_voxel_center[comparing_dimension2]);
		double distance_squared= std::pow( pCell->position[comparing_dimension] - line_coord1,2)+ std::pow( pCell->position[comparing_dimension2] - line_coord2,2);
		if(distance_squared > max_interactive_distance * max_interactive_distance)
		{ return false; }
		return true;
	}
	std::vector<double> corner_point= 0.5*(my_voxel_center+other_voxel_center);
	double distance_squared= (corner_point[0]-pCell->position[0])*(corner_point[0]-pCell->position[0])
		+(corner_point[1]-pCell->position[1])*(corner_point[1]-pCell->position[1])
		+(corner_point[2]-pCell->position[2]) * (corner_point[2]-pCell->position[2]);
	if(distance_squared > max_interactive_distance * max_interactive_distance)
	{ return false; }
	return true;
}

/*------------------------------------------------------------------------------------------*/
/* Parallel version of the is_neighbor_voxel() function above, the 4th argument here is the */
/* max_voxel_interactive_distance_of_other_voxel INSTEAD of the neighbour voxel index. 			*/
/* There is no need for mpi_Environment, mpi_Cartesian objects but maybe will be used later	*/
/*------------------------------------------------------------------------------------------*/

bool is_neighbor_voxel(Cell* pCell, std::vector<double> my_voxel_center, std::vector<double> other_voxel_center, double max_cell_interactive_distance_in_voxel, mpi_Environment &world, mpi_Cartesian &cart_topo)
{
	double max_interactive_distance = pCell->phenotype.mechanics.relative_maximum_adhesion_distance *
																		pCell->phenotype.geometry.radius +
																		max_cell_interactive_distance_in_voxel;

	int comparing_dimension = -1, comparing_dimension2 = -1;

	if(my_voxel_center[0] == other_voxel_center[0] && my_voxel_center[1] == other_voxel_center[1])
	{
		comparing_dimension = 2;
	}
	else if(my_voxel_center[0] == other_voxel_center[0] && my_voxel_center[2] == other_voxel_center[2])
	{
		comparing_dimension = 1;
	}
	else if(my_voxel_center[1] == other_voxel_center[1] && my_voxel_center[2] == other_voxel_center[2])
	{
		comparing_dimension = 0;
	}

	if(comparing_dimension != -1)
	{
		//then it is an immediate neighbor (through side faces)

		double surface_coord= 0.5*(my_voxel_center[comparing_dimension] + other_voxel_center[comparing_dimension]);
		if(std::fabs(pCell->position[comparing_dimension] - surface_coord) > max_interactive_distance)
		{
			return false;
		}
		return true;
	}

	comparing_dimension=-1;

	if(my_voxel_center[0] == other_voxel_center[0])
	{
		comparing_dimension = 1; comparing_dimension2 = 2;
	}
	else if(my_voxel_center[1] == other_voxel_center[1])
	{
		comparing_dimension=0; comparing_dimension2 = 2;
	}
	else if(my_voxel_center[2] == other_voxel_center[2])
	{
		comparing_dimension = 0; comparing_dimension2=1;
	}
	if(comparing_dimension != -1)
	{
		double line_coord1= 0.5*(my_voxel_center[comparing_dimension] + other_voxel_center[comparing_dimension]);
		double line_coord2= 0.5*(my_voxel_center[comparing_dimension2] + other_voxel_center[comparing_dimension2]);
		double distance_squared= std::pow( pCell->position[comparing_dimension] - line_coord1,2)+ std::pow( pCell->position[comparing_dimension2] - line_coord2,2);

		if(distance_squared > max_interactive_distance * max_interactive_distance)
		{
			return false;
		}
		return true;
	}

	std::vector<double> corner_point= 0.5*(my_voxel_center+other_voxel_center);
	double distance_squared = (corner_point[0]-pCell->position[0])* (corner_point[0]-pCell->position[0])+
														(corner_point[1]-pCell->position[1])* (corner_point[1]-pCell->position[1])+
														(corner_point[2]-pCell->position[2])* (corner_point[2]-pCell->position[2]);

	if(distance_squared > max_interactive_distance * max_interactive_distance)
	{
		return false;
	}

	return true;
}


std::vector<Cell*>& Cell::cells_in_my_container( void )
{
	return get_container()->agent_grid[get_current_mechanics_voxel_index()];
}

std::vector<Cell*> Cell::nearby_cells( void )
{ return find_nearby_cells( this ); }

std::vector<Cell*> Cell::nearby_interacting_cells( void )
{ return find_nearby_interacting_cells( this ); }

void Cell::ingest_cell( Cell* pCell_to_eat )
{
	// don't ingest self 
	if( pCell_to_eat == this )
	{ return; } 
	
	// don't ingest a cell that's already ingested 
	if( pCell_to_eat->phenotype.volume.total < 1e-15 )
	{ return; } 
		
	// make this thread safe 
	#pragma omp critical
	{
		bool volume_was_zero = false; 
		if( pCell_to_eat->phenotype.volume.total < 1e-15 )
		{
			volume_was_zero = true; 
			//std::cout << this << " " << this->type_name << " ingests " 
			//<< pCell_to_eat << " " << pCell_to_eat->type_name << std::endl; 
		}
		
		// mark it as dead 
		pCell_to_eat->phenotype.death.dead = true; 
		// set secretion and uptake to zero 
		pCell_to_eat->phenotype.secretion.set_all_secretion_to_zero( );  
		pCell_to_eat->phenotype.secretion.set_all_uptake_to_zero( );

		// deactivate all custom function 
		pCell_to_eat->functions.custom_cell_rule = NULL; 
		pCell_to_eat->functions.update_phenotype = NULL; 
		pCell_to_eat->functions.contact_function = NULL; 
		pCell_to_eat->functions.cell_division_function = NULL; 

		// should set volume fuction to NULL too! 
		pCell_to_eat->functions.volume_update_function = NULL; 

		// set cell as unmovable and non-secreting 
		pCell_to_eat->is_movable = false; 
		pCell_to_eat->is_active = false; 


		// absorb fluid volume (all into the cytoplasm) 
		phenotype.volume.cytoplasmic_fluid += pCell_to_eat->phenotype.volume.fluid; 
		pCell_to_eat->phenotype.volume.cytoplasmic_fluid = 0.0; 
		
		// absorb nuclear and cyto solid volume (into the cytoplasm) 
		phenotype.volume.cytoplasmic_solid += pCell_to_eat->phenotype.volume.cytoplasmic_solid; 
		pCell_to_eat->phenotype.volume.cytoplasmic_solid = 0.0; 
		
		phenotype.volume.cytoplasmic_solid += pCell_to_eat->phenotype.volume.nuclear_solid; 
		pCell_to_eat->phenotype.volume.nuclear_solid = 0.0; 
		
		// consistency calculations 
		
		phenotype.volume.fluid = phenotype.volume.nuclear_fluid + 
			phenotype.volume.cytoplasmic_fluid; 
		pCell_to_eat->phenotype.volume.fluid = 0.0; 
		
		phenotype.volume.solid = phenotype.volume.cytoplasmic_solid + 
			phenotype.volume.nuclear_solid; 
		pCell_to_eat->phenotype.volume.solid = 0.0; 
		
		// no change to nuclear volume (initially) 
		pCell_to_eat->phenotype.volume.nuclear = 0.0; 
		pCell_to_eat->phenotype.volume.nuclear_fluid = 0.0; 
		
		phenotype.volume.cytoplasmic = phenotype.volume.cytoplasmic_solid + 
			phenotype.volume.cytoplasmic_fluid; 
		pCell_to_eat->phenotype.volume.cytoplasmic = 0.0; 
		
		phenotype.volume.total = phenotype.volume.nuclear + 
			phenotype.volume.cytoplasmic; 
		pCell_to_eat->phenotype.volume.total = 0.0; 

		phenotype.volume.fluid_fraction = phenotype.volume.fluid / 
			(  phenotype.volume.total + 1e-16 ); 
		pCell_to_eat->phenotype.volume.fluid_fraction = 0.0; 

		phenotype.volume.cytoplasmic_to_nuclear_ratio = phenotype.volume.cytoplasmic_solid / 
			( phenotype.volume.nuclear_solid + 1e-16 );
			
		// update corresponding BioFVM parameters (self-consistency) 
		set_total_volume( phenotype.volume.total ); 
		pCell_to_eat->set_total_volume( 0.0 ); 
		
		// absorb the internalized substrates 
		
		// multiply by the fraction that is supposed to be ingested (for each substrate) 
		*(pCell_to_eat->internalized_substrates) *= 
			*(pCell_to_eat->fraction_transferred_when_ingested); // 
		
		*internalized_substrates += *(pCell_to_eat->internalized_substrates); 
		static int n_substrates = internalized_substrates->size(); 
		pCell_to_eat->internalized_substrates->assign( n_substrates , 0.0 );
		
		// conserved quantitites in custom data during phagocytosis
		// so that phagocyte cell absorbs the full amount from the engulfed cell;
		for( int nn = 0 ; nn < custom_data.variables.size() ; nn++ )
		{
			if( custom_data.variables[nn].conserved_quantity == true )
			{
				custom_data.variables[nn].value += 
				pCell_to_eat->custom_data.variables[nn].value; 			
			}
		}
		for( int nn = 0 ; nn < custom_data.vector_variables.size() ; nn++ )
		{
			if( custom_data.vector_variables[nn].conserved_quantity == true )
			{
				custom_data.vector_variables[nn].value += 
				pCell_to_eat->custom_data.vector_variables[nn].value; 
			}
		}
		
	}

	// things that have their own thread safety 
	pCell_to_eat->flag_for_removal();
	pCell_to_eat->remove_all_attached_cells();
	pCell_to_eat->remove_all_spring_attachments();
	
	return; 
}

void Cell::attack_cell( Interacting_Cell_Info* pCell_to_attack , double dt )
{
	if( pCell_to_attack->dead == true )
	{ return; } 

	// make this thread safe 
	#pragma omp critical
	{ 
		double new_damage = phenotype.cell_interactions.attack_damage_rate * dt;
		pCell_to_attack->attacked = true;
		pCell_to_attack->damage_suffered += new_damage; 
		pCell_to_attack->time_attacked += dt; 
		// Maintain coherency with local attack: update attacker's damage delivered
		phenotype.cell_interactions.total_damage_delivered += new_damage;
	}
	return; 
}

void Cell::was_attacked(double damage_suffered, double time_attacked) {
	#pragma omp critical
	{ 
		phenotype.cell_integrity.damage += damage_suffered;
		// Maintain coherency with local attack: update total attack time from ghost
		state.total_attack_time += time_attacked;
	}
}

void Cell::fuse_cell( Interacting_Cell_Info* pCell_to_fuse, mpi_Environment &world, mpi_Cartesian &cart_topo) 
{
	// don't ingest a cell that's already fused 
	if( pCell_to_fuse->phenotype_volume < 1e-15 )
	{ return; } 

	// make this thread safe 
	#pragma omp critical
	{
		bool volume_was_zero = false; 
		if( pCell_to_fuse->phenotype_volume < 1e-15 )
		{
			volume_was_zero = true; 
			std::cout << this << " " << this->type_name << " fuses " 
			<< pCell_to_fuse << " " << pCell_to_fuse->type << std::endl; 
		}

		// set new position at center of volume 
			// x_new = (vol_B * x_B + vol_S * x_S ) / (vol_B + vol_S )

		std::vector<double> new_position = position; // x_B
		new_position *= phenotype.volume.total; // vol_B * x_B 
		double total_volume = phenotype.volume.total; 
		total_volume += pCell_to_fuse->phenotype_volume ;  

		axpy( &new_position , pCell_to_fuse->phenotype_volume , pCell_to_fuse->position ); // vol_B*x_B + vol_S*x_S
		new_position /= total_volume; // (vol_B*x_B+vol_S*x_S)/(vol_B+vol_S);

		if (std::isnan(new_position[0])) std::cout << "Error in line " << __func__ << " " << __LINE__ << std::endl;

		static double xL = get_default_microenvironment()->mesh.bounding_box[0];		 
		static double xU = get_default_microenvironment()->mesh.bounding_box[3]; 

		static double yL = get_default_microenvironment()->mesh.bounding_box[1];		 
		static double yU = get_default_microenvironment()->mesh.bounding_box[4]; 

		static double zL = get_default_microenvironment()->mesh.bounding_box[2];		 
		static double zU = get_default_microenvironment()->mesh.bounding_box[5]; 

		if( new_position[0] < xL || new_position[0] > xU || 
		    new_position[1] < yL || new_position[1] > yU || 
			new_position[2] < zL || new_position[2] > zU )
		{
			std::cout << "cell fusion at " << new_position << " violates domain bounds" << std::endl; 
			std::cout << get_default_microenvironment()->mesh.bounding_box << std::endl << std::endl; 
		}
		position = new_position;
		//int aux = current_voxel_index; 
		//update_voxel_in_container(world, cart_topo);

		// set number of nuclei 

		state.number_of_nuclei += pCell_to_fuse->number_of_nuclei; 

		// absorb all the volume(s)

		// absorb fluid volume (all into the cytoplasm) 
		phenotype.volume.cytoplasmic_fluid += pCell_to_fuse->cytoplasmic_fluid; 

		phenotype.volume.nuclear_fluid += pCell_to_fuse->nuclear_fluid; 

		// absorb nuclear and cyto solid volume (into the cytoplasm) 
		phenotype.volume.cytoplasmic_solid += pCell_to_fuse->cytoplasmic_solid; 

		phenotype.volume.nuclear_solid += pCell_to_fuse->nuclear_solid; 

		// consistency calculations 

		phenotype.volume.fluid = phenotype.volume.nuclear_fluid + 
			phenotype.volume.cytoplasmic_fluid;  

		phenotype.volume.solid = phenotype.volume.cytoplasmic_solid + 
			phenotype.volume.nuclear_solid; 

		phenotype.volume.nuclear = phenotype.volume.nuclear_fluid + 
			phenotype.volume.nuclear_solid; 

		phenotype.volume.cytoplasmic = phenotype.volume.cytoplasmic_fluid + 
			phenotype.volume.cytoplasmic_solid;  

		phenotype.volume.total = phenotype.volume.nuclear + 
			phenotype.volume.cytoplasmic; 

		phenotype.volume.fluid_fraction = phenotype.volume.fluid / 
			(  phenotype.volume.total + 1e-16 ); 

		phenotype.volume.cytoplasmic_to_nuclear_ratio = phenotype.volume.cytoplasmic_solid / 
			( phenotype.volume.nuclear_solid + 1e-16 );

		// update corresponding BioFVM parameters (self-consistency) 
		set_total_volume( phenotype.volume.total ); 

		// absorb the internalized substrates 

		*internalized_substrates += pCell_to_fuse->internalized_substrates; 

		// set target volume(s)

		phenotype.volume.target_solid_cytoplasmic += pCell_to_fuse->target_solid_cytoplasmic;
		phenotype.volume.target_solid_nuclear += pCell_to_fuse->target_solid_nuclear;

		// trigger removal from the simulation 
		// pCell_to_eat->die(); // I don't think this is safe if it's in an OpenMP loop 

		// flag it for removal 
		// pCell_to_eat->flag_for_removal(); 
		// mark it as dead 
		pCell_to_fuse->dead = true;
		pCell_to_fuse->fused = true; 

	}

	return; 
}

void Cell::was_fused() {

	if( phenotype.death.dead == true )
	{ return; } 

	// absorb fluid volume (all into the cytoplasm) 
		phenotype.volume.cytoplasmic_fluid = 0.0; 

		phenotype.volume.nuclear_fluid = 0.0; 

		// absorb nuclear and cyto solid volume (into the cytoplasm) 
		phenotype.volume.cytoplasmic_solid = 0.0; 

		phenotype.volume.nuclear_solid = 0.0; 

		// consistency calculations 

		phenotype.volume.fluid = 0.0; 

		phenotype.volume.solid = 0.0; 

		phenotype.volume.nuclear = 0.0; 

		phenotype.volume.cytoplasmic = 0.0; 

		phenotype.volume.total = 0.0; 

		phenotype.volume.fluid_fraction = 0.0; 

		// update corresponding BioFVM parameters (self-consistency) 
		set_total_volume( 0.0 ); 

		//Delete internalized substrates
		static int n_substrates = internalized_substrates->size(); 
		internalized_substrates->assign( n_substrates , 0.0 ); 	

		phenotype.death.dead = true; 

		phenotype.secretion.set_all_secretion_to_zero( );  
		phenotype.secretion.set_all_uptake_to_zero( ); 

		// deactivate all custom function 
		functions.custom_cell_rule = NULL; 
		functions.update_phenotype = NULL; 
		functions.contact_function = NULL; 
		functions.volume_update_function = NULL; 

		// remove all adhesions 
		// pCell_to_eat->remove_all_attached_cells();

		// set cell as unmovable and non-secreting 
		is_movable = false; 
		is_active = false; 

		// things that have their own thread safety 
		flag_for_removal();
		remove_all_attached_cells();
}

void Cell::ingest_cell( Interacting_Cell_Info* pCell_to_eat )
{
	// don't ingest a cell that's already ingested 
	if( pCell_to_eat->phenotype_volume < 1e-15)
	{ return; } 
		
	// make this thread safe 
	#pragma omp critical
	{
		bool volume_was_zero = false; 
		if( pCell_to_eat->phenotype_volume < 1e-15 )
		{
			volume_was_zero = true; 
			std::cout << this << " " << this->type_name << " ingests " 
			<< pCell_to_eat << " " << pCell_to_eat->type << std::endl; 
		}
		// absorb all the volume(s)

		// absorb fluid volume (all into the cytoplasm) 
		phenotype.volume.cytoplasmic_fluid += pCell_to_eat->cytoplasmic_fluid + pCell_to_eat->nuclear_fluid; 
		
		// absorb nuclear and cyto solid volume (into the cytoplasm) 
		phenotype.volume.cytoplasmic_solid += pCell_to_eat->cytoplasmic_solid; 
		
		phenotype.volume.cytoplasmic_solid += pCell_to_eat->nuclear_solid; 
		
		// consistency calculations 
		
		phenotype.volume.fluid = phenotype.volume.nuclear_fluid + 
			phenotype.volume.cytoplasmic_fluid; 
		
		phenotype.volume.solid = phenotype.volume.cytoplasmic_solid + 
			phenotype.volume.nuclear_solid;  
		
		phenotype.volume.cytoplasmic = phenotype.volume.cytoplasmic_solid + 
			phenotype.volume.cytoplasmic_fluid; 
		
		phenotype.volume.total = phenotype.volume.nuclear + 
			phenotype.volume.cytoplasmic; 

		phenotype.volume.fluid_fraction = phenotype.volume.fluid / 
			(  phenotype.volume.total + 1e-16 ); 

		phenotype.volume.cytoplasmic_to_nuclear_ratio = phenotype.volume.cytoplasmic_solid / 
			( phenotype.volume.nuclear_solid + 1e-16 );
			
		// update corresponding BioFVM parameters (self-consistency) 
		set_total_volume( phenotype.volume.total ); 
		
		// absorb the internalized substrates 
		
		// multiply by the fraction that is supposed to be ingested (for each substrate) 
		pCell_to_eat->internalized_substrates *= 
			pCell_to_eat->fraction_transferred_when_ingested; // 
		
		*internalized_substrates += pCell_to_eat->internalized_substrates; 
		
		// trigger removal from the simulation 
		// pCell_to_eat->die(); // I don't think this is safe if it's in an OpenMP loop 
		
		// flag it for removal 
		// pCell_to_eat->flag_for_removal(); 
		// mark it as dead 
		pCell_to_eat->dead = true; 
		pCell_to_eat->ingested = true;
	
	}
	
	return; 
}

void Cell::was_ingested(){
	// absorb fluid volume (all into the cytoplasm) 
	#pragma omp critical
	{
		phenotype.volume.cytoplasmic_fluid = 0.0; 
		
		// absorb nuclear and cyto solid volume (into the cytoplasm) 
		phenotype.volume.cytoplasmic_solid = 0.0; 
		
		phenotype.volume.nuclear_solid = 0.0; 
		
		// consistency calculations 
		

		phenotype.volume.fluid = 0.0; 
		
		phenotype.volume.solid = 0.0; 
		
		// no change to nuclear volume (initially) 
		phenotype.volume.nuclear = 0.0; 
		phenotype.volume.nuclear_fluid = 0.0; 
		
		phenotype.volume.cytoplasmic = 0.0; 
		
		phenotype.volume.total = 0.0; 

		phenotype.volume.fluid_fraction = 0.0; 

			
		// update corresponding BioFVM parameters (self-consistency) 

		set_total_volume( 0.0 ); 
		
		// absorb the internalized substrates 
		static int n_substrates = internalized_substrates->size(); 
		internalized_substrates->assign( n_substrates , 0.0 ); 	
		
		// trigger removal from the simulation 
		// pCell_to_eat->die(); // I don't think this is safe if it's in an OpenMP loop 
		
		// flag it for removal 
		// pCell_to_eat->flag_for_removal(); 
		// mark it as dead 
		phenotype.death.dead = true; 
		// set secretion and uptake to zero 
		phenotype.secretion.set_all_secretion_to_zero( );  
		phenotype.secretion.set_all_uptake_to_zero( ); 
		
		// deactivate all custom function 
		functions.custom_cell_rule = NULL; 
		functions.update_phenotype = NULL; 
		functions.contact_function = NULL; 

		// remove all adhesions 
		// pCell_to_eat->remove_all_attached_cells();
		
		// set cell as unmovable and non-secreting 
		is_movable = false; 
		is_active = false; 
	}

	// things that have their own thread safety 
	flag_for_removal();
	remove_all_attached_cells();
}

void Cell::attack_cell( Cell* pCell_to_attack , double dt )
{
	// don't attack self 
	if( pCell_to_attack == this )
	{ return; } 
	
	// don't attack a dead or tiny cell 
	if( pCell_to_attack->phenotype.death.dead == true || pCell_to_attack->phenotype.volume.total < 1e-15 )
	{ return; }

	// make this thread safe 
	#pragma omp critical
	{ 
		// std::cout << this->type_name << " attacks " << pCell_to_attack->type_name << std::endl;
		// 
		double new_damage = phenotype.cell_interactions.attack_damage_rate * dt; 

		pCell_to_attack->phenotype.cell_integrity.damage += new_damage; 
		pCell_to_attack->state.total_attack_time += dt; 

		phenotype.cell_interactions.total_damage_delivered += new_damage;  
	}
	return; 
}

//Jose: check that the new position is inside the domain or traspassed to a corresponding sub-domain.
void Cell::fuse_cell( Cell* pCell_to_fuse )
{
	// don't ingest a cell that's already fused 
	if( pCell_to_fuse->phenotype.volume.total < 1e-15 || this == pCell_to_fuse )
	{ return; } 

	// make this thread safe 
	#pragma omp critical
	{

		// set new position at center of volume 
			// x_new = (vol_B * x_B + vol_S * x_S ) / (vol_B + vol_S )

		std::vector<double> new_position = position; // x_B
		new_position *= phenotype.volume.total; // vol_B * x_B 
		double total_volume = phenotype.volume.total; 
		total_volume += pCell_to_fuse->phenotype.volume.total ;  

		axpy( &new_position , pCell_to_fuse->phenotype.volume.total , pCell_to_fuse->position ); // vol_B*x_B + vol_S*x_S
		new_position /= total_volume; // (vol_B*x_B+vol_S*x_S)/(vol_B+vol_S);

		static double xL = get_default_microenvironment()->mesh.bounding_box[0];		 
		static double xU = get_default_microenvironment()->mesh.bounding_box[3]; 

		static double yL = get_default_microenvironment()->mesh.bounding_box[1];		 
		static double yU = get_default_microenvironment()->mesh.bounding_box[4]; 

		static double zL = get_default_microenvironment()->mesh.bounding_box[2];		 
		static double zU = get_default_microenvironment()->mesh.bounding_box[5]; 

		if( new_position[0] < xL || new_position[0] > xU || 
		    new_position[1] < yL || new_position[1] > yU || 
			new_position[2] < zL || new_position[2] > zU )
		{
			std::cout << "cell fusion at " << new_position << " violates domain bounds" << std::endl; 
			std::cout << get_default_microenvironment()->mesh.bounding_box << std::endl << std::endl; 
		}
		position = new_position; 
		update_voxel_in_container();

		// set number of nuclei 

		state.number_of_nuclei += pCell_to_fuse->state.number_of_nuclei; 

		// absorb all the volume(s)

		// absorb fluid volume (all into the cytoplasm) 
		phenotype.volume.cytoplasmic_fluid += pCell_to_fuse->phenotype.volume.cytoplasmic_fluid; 
		pCell_to_fuse->phenotype.volume.cytoplasmic_fluid = 0.0; 

		phenotype.volume.nuclear_fluid += pCell_to_fuse->phenotype.volume.nuclear_fluid; 
		pCell_to_fuse->phenotype.volume.nuclear_fluid = 0.0; 

		// absorb nuclear and cyto solid volume (into the cytoplasm) 
		phenotype.volume.cytoplasmic_solid += pCell_to_fuse->phenotype.volume.cytoplasmic_solid; 
		pCell_to_fuse->phenotype.volume.cytoplasmic_solid = 0.0; 

		phenotype.volume.nuclear_solid += pCell_to_fuse->phenotype.volume.nuclear_solid; 
		pCell_to_fuse->phenotype.volume.nuclear_solid = 0.0; 

		// consistency calculations 

		phenotype.volume.fluid = phenotype.volume.nuclear_fluid + 
			phenotype.volume.cytoplasmic_fluid; 
		pCell_to_fuse->phenotype.volume.fluid = 0.0; 

		phenotype.volume.solid = phenotype.volume.cytoplasmic_solid + 
			phenotype.volume.nuclear_solid; 
		pCell_to_fuse->phenotype.volume.solid = 0.0; 

		phenotype.volume.nuclear = phenotype.volume.nuclear_fluid + 
			phenotype.volume.nuclear_solid; 
		pCell_to_fuse->phenotype.volume.nuclear = 0.0; 

		phenotype.volume.cytoplasmic = phenotype.volume.cytoplasmic_fluid + 
			phenotype.volume.cytoplasmic_solid; 
		pCell_to_fuse->phenotype.volume.cytoplasmic = 0.0; 

		phenotype.volume.total = phenotype.volume.nuclear + 
			phenotype.volume.cytoplasmic; 
		pCell_to_fuse->phenotype.volume.total = 0.0; 

		phenotype.volume.fluid_fraction = phenotype.volume.fluid / 
			(  phenotype.volume.total + 1e-16 ); 
		pCell_to_fuse->phenotype.volume.fluid_fraction = 0.0; 

		phenotype.volume.cytoplasmic_to_nuclear_ratio = phenotype.volume.cytoplasmic_solid / 
			( phenotype.volume.nuclear_solid + 1e-16 );

		// update corresponding BioFVM parameters (self-consistency) 
		set_total_volume( phenotype.volume.total ); 
		pCell_to_fuse->set_total_volume( 0.0 ); 

		// absorb the internalized substrates 

		*internalized_substrates += *(pCell_to_fuse->internalized_substrates); 
		static int n_substrates = internalized_substrates->size(); 
		pCell_to_fuse->internalized_substrates->assign( n_substrates , 0.0 ); 	

		// set target volume(s)

		phenotype.volume.target_solid_cytoplasmic += pCell_to_fuse->phenotype.volume.target_solid_cytoplasmic;
		phenotype.volume.target_solid_nuclear += pCell_to_fuse->phenotype.volume.target_solid_nuclear;

		// trigger removal from the simulation 
		// pCell_to_eat->die(); // I don't think this is safe if it's in an OpenMP loop 

		// flag it for removal 
		// pCell_to_eat->flag_for_removal(); 
		// mark it as dead 
		pCell_to_fuse->phenotype.death.dead = true; 
		// set secretion and uptake to zero 
		pCell_to_fuse->phenotype.secretion.set_all_secretion_to_zero( );  
		pCell_to_fuse->phenotype.secretion.set_all_uptake_to_zero( ); 

		// deactivate all custom function 
		pCell_to_fuse->functions.custom_cell_rule = NULL; 
		pCell_to_fuse->functions.update_phenotype = NULL; 
		pCell_to_fuse->functions.contact_function = NULL; 
		pCell_to_fuse->functions.volume_update_function = NULL; 

		// remove all adhesions 
		// pCell_to_eat->remove_all_attached_cells();

		// set cell as unmovable and non-secreting 
		pCell_to_fuse->is_movable = false; 
		pCell_to_fuse->is_active = false; 

		

	}

	// things that have their own thread safety 
	pCell_to_fuse->flag_for_removal();
	pCell_to_fuse->remove_all_attached_cells();
	pCell_to_fuse->remove_all_spring_attachments();


	return; 
}

void Cell::lyse_cell( void )
{
	// don't lyse a cell that's already lysed 
	if( phenotype.volume.total < 1e-15 )
	{ return; } 	
	
	// flag for removal 
	flag_for_removal(); // should be safe now 
	
	// mark it as dead 
	phenotype.death.dead = true; 
	
	// set secretion and uptake to zero 
	phenotype.secretion.set_all_secretion_to_zero( );  
	phenotype.secretion.set_all_uptake_to_zero( ); 
	
	// deactivate all custom function 
	functions.custom_cell_rule = NULL; 
	functions.update_phenotype = NULL; 
	functions.contact_function = NULL; 
	
	// remove all adhesions 
	
	remove_all_attached_cells(); 
	
	// set volume to zero 
	set_total_volume( 0.0 ); 

	// set cell as unmovable and non-secreting 
	is_movable = false; 
	is_active = false; 	

	return; 
}


/*---------------------------------------------*/
/* Added by Gaurav Saxena to print 'most' of	 */
/* the fields of the Cell data structure,			 */
/* useful when comparing cell fields after 		 */
/* pack_and_send() function. 					*/
/*---------------------------------------------*/
void Cell::print_cell(mpi_Environment &world)
{
	/*-------------------------------------------------------------------------------------*/
	/* Trying to print all the data fields of the Cell class 															 */
	/* Since Cell is publicly derived from the Basic_Agent class, the protected and	public */
	/* member of Basic_Agent can be accessed in member functions of the Cell class but the */
	/* private members of the Basic_Agent class cannot be accessed by member functions of  */
	/* Derived class. For example: "volume" is not directly accessible but cell_source_sink*/
	/* solver_temp1[i] is directly accessible. 																						 */
	/*-------------------------------------------------------------------------------------*/

	/*----------------------------------------------------------------------------*/
	/* REWRITING THE PRINT_CELL() FUNCTION TO *NOT* PRINT ON THE COMMAND LINE   	*/
	/* Each process (1) will write its own file (2) File will contain cell data 	*/
	/* of cells BEFORE packing and AFTER unpacking. (3) Text will mention if cell */
	/* is moving to LEFT or RIGHT																									*/
	/*----------------------------------------------------------------------------*/

	std::string filename="CELLS_RANK_";
	filename = filename + std::to_string(world.rank);

	std::ofstream ofile;
	ofile.open(filename,std::ios::app);


	ofile<<"=================CELL DATA FOR CELL ID: "<<this->ID<<"======================"<<std::endl;

	if(crossed_to_left_subdomain == true)
		ofile<<"PACKED CELL MOVING TO LEFT"<<std::endl;

	if(crossed_to_right_subdomain == true)
		ofile<<"PACKED CELL MOVING TO RIGHT"<<std::endl;

	if(crossed_to_left_subdomain == false && crossed_to_right_subdomain == false)
		ofile<<"UNPACKED CELL"<<std::endl;

	ofile << "NEW POSITION ="<<"("<<position[0]<<","<<position[1]<<","<<position[2]<<")"<<std::endl;

	ofile<<"=> class Cell {"<<std::endl;
	ofile<<"TYPE_NAME:"<<type_name<<std::endl;


	ofile<<"=> class Cell { class Custom_Cell_Data {"<<std::endl;
	for(auto it = (custom_data.get_name_to_index_map()).cbegin(); it != (custom_data.get_name_to_index_map()).cend(); it++)
		ofile << "KEY="<<it->first<<" VALUE="<<it->second<<std::endl;


	ofile<<"=> class Cell { class Custom_Cell_Data { class Variable { "<<std::endl;
	for(int i=0; i<this->custom_data.variables.size(); i++)
		ofile<<"i:"<<i<<" name:"<<this->custom_data.variables[i].name<<" value:"<<this->custom_data.variables[i].value 
			<< " units:"<<this->custom_data.variables[i].units << " conserved quantity:"<<this->custom_data.variables[i].conserved_quantity <<std::endl;


	ofile<<"=> class Cell { class Custom_Cell_Data { class Vector_Variable { "<<std::endl;
	for(int i=0; i<this->custom_data.vector_variables.size();i++)
	{
		ofile<<"i:"<<i<<" name:"<<this->custom_data.vector_variables[i].name<<" units:"<<this->custom_data.vector_variables[i].units << " conserved quantity:"<<this->custom_data.vector_variables[i].conserved_quantity << std::endl;
		for(int j=0; j<this->custom_data.vector_variables[i].value.size();j++)
			ofile<<" value "<<j<<":"<<this->custom_data.vector_variables[i].value[j]<<std::endl;
	}


	ofile<<"=> class Cell { class Cell_Parameters { "<<std::endl;
	ofile<<"o2_hypoxic_threshold:"<<this->parameters.o2_hypoxic_threshold<<std::endl;
	ofile<<"o2_hypoxic_response:"<<this->parameters.o2_hypoxic_response<<std::endl;
	ofile<<"o2_hypoxic_saturation:"<<this->parameters.o2_hypoxic_saturation<<std::endl;
	ofile<<"o2_proliferation_saturation:"<<this->parameters.o2_proliferation_saturation<<std::endl;
	ofile<<"o2_proliferation_threshold:"<<this->parameters.o2_proliferation_threshold<<std::endl;
	ofile<<"o2_reference:"<<this->parameters.o2_reference<<std::endl;
	ofile<<"o2_necrosis_threshold:"<<this->parameters.o2_necrosis_threshold<<std::endl;
	ofile<<"o2_necrosis_max:"<<this->parameters.o2_necrosis_max<<std::endl;
	ofile<<"max_necrosis_rate:"<<this->parameters.max_necrosis_rate<<std::endl;
	ofile<<"necrosis_type:"<<this->parameters.necrosis_type<<std::endl;
	if(this->parameters.pReference_live_phenotype != NULL)
		ofile<<"pReference_live_phenotype pointer in Cell_Parameters is NOT NULL"<<std::endl;
	if(this->parameters.pReference_live_phenotype == &phenotype)
		ofile<<"pReference_live_phenotype pointer in Cell_Parameters is pointing to object phenotype of class Phenotype"<<std::endl;


	ofile<<"=> class Cell { class Functions { class Cycle_Model { "<<std::endl;
	std::vector<std::unordered_map<int,int>> &inverse_index_maps_cycle_model = this->functions.cycle_model.get_inverse_index_maps();
	ofile<<"Vector of unordered maps"<<std::endl;
	for(int i=0; i<inverse_index_maps_cycle_model.size(); i++)
	{
		for(auto it = inverse_index_maps_cycle_model[i].cbegin(); it != inverse_index_maps_cycle_model[i].cend(); it++)
			ofile<<"{"<<(*it).first<<":"<<(*it).second<<"}"<<std::endl;
	}
	ofile<<"name:"<<this->functions.cycle_model.name<<std::endl;
	ofile<<"code:"<<this->functions.cycle_model.code<<std::endl;

	ofile<<"=> class Cell { class Functions { class Cycle_Model { class Phase {"<<std::endl;
	for(int i=0; i<this->functions.cycle_model.phases.size();i++)
	{
		ofile<<"index:"	<<this->functions.cycle_model.phases[i].index<<std::endl;
		ofile<<"code:"	<<this->functions.cycle_model.phases[i].code<<std::endl;
		ofile<<"name:"	<<this->functions.cycle_model.phases[i].name<<std::endl;
		ofile<<"division_at_phase_exit:"	<<this->functions.cycle_model.phases[i].division_at_phase_exit<<std::endl;
		ofile<<"removal_at_phase_exit:"		<<this->functions.cycle_model.phases[i].removal_at_phase_exit<<std::endl;
	}

	ofile<<"=> class Cell { class Functions { class Cycle_Model { class Phase_Links {"<<std::endl;
	for(int i=0; i<this->functions.cycle_model.phase_links.size();i++)
		{
			ofile<<"Vector:"<<i<<std::endl;
			for(int j=0; j<this->functions.cycle_model.phase_links[i].size();j++)
			{
				ofile<<"Vector:"<<j<<std::endl;
				ofile<<"start_phase_index:"	<<this->functions.cycle_model.phase_links[i][j].start_phase_index	<<std::endl;
				ofile<<"end_phase_index:"		<<this->functions.cycle_model.phase_links[i][j].end_phase_index		<<std::endl;
				ofile<<"fixed_duration:"		<<this->functions.cycle_model.phase_links[i][j].fixed_duration		<<std::endl;
			}
		}

	ofile<<"=> class Cell { class Functions { class Cycle_Model {"<<std::endl;
	ofile<<"default phase index:"<<this->functions.cycle_model.default_phase_index<<std::endl;

	ofile<<"=> class Cell { class Functions { class Cycle_Model { class Cycle_Data { "<<std::endl;

	std::vector<std::unordered_map<int,int>> &inverse_index_maps_cycle_data_cmcd = this->functions.cycle_model.data.get_inverse_index_maps();
	ofile<<"Vector of unordered maps"<<std::endl;
	for(int i=0; i<inverse_index_maps_cycle_data_cmcd.size(); i++)
	{
		for(auto it = inverse_index_maps_cycle_data_cmcd[i].cbegin(); it != inverse_index_maps_cycle_data_cmcd[i].cend(); it++)
			ofile<<"{"<<(*it).first<<":"<<(*it).second<<"}"<<std::endl;
		ofile<<std::endl;
	}
	if(this->functions.cycle_model.data.pCycle_Model != NULL)
		ofile<<"pCycle_Model is NOT NULL"<<std::endl;
	ofile<<"time_units:"<<this->functions.cycle_model.data.time_units<<std::endl;
	for(int i=0; i<this->functions.cycle_model.data.transition_rates.size();i++)
		{
			ofile<<"Transition Rates Vector:"<<i<<std::endl;
			for(int j=0; j<this->functions.cycle_model.data.transition_rates[i].size();j++)
				ofile<<"rate "<<j<<":"<<this->functions.cycle_model.data.transition_rates[i][j]<<std::endl;
		}
	ofile<<"current_phase_index:"<<this->functions.cycle_model.data.current_phase_index<<std::endl;
	ofile<<"elapsed_time_in_phase:"<<this->functions.cycle_model.data.elapsed_time_in_phase<<std::endl;

	ofile<<"=> class Cell { class  State {"<<std::endl;
	ofile<<"Orientation[0]:"<< (this->state).orientation[0] << " Orientation[1]:" << (this->state).orientation[1] << " Orientation[2]:" << (this->state).orientation[2] << std::endl;
	ofile<<"Simple Pressure:"<<this->state.simple_pressure<<std::endl;


	ofile<<"=> class Cell { class Phenotype { "<<std::endl;
	ofile<<"flagged_for_division:"<<this->phenotype.flagged_for_division<<std::endl;
	ofile<<"flagged_for_removal:"<<this->phenotype.flagged_for_removal<<std::endl;

	ofile<<"=> class Cell { class Phenotype { class Cycle { class Cycle_Data {  "<<std::endl;
	std::vector<std::unordered_map<int,int>> &inverse_index_maps_cycle_data = this->phenotype.cycle.data.get_inverse_index_maps();
	ofile<<"Vector of unordered maps"<<std::endl;
	for(int i=0; i<inverse_index_maps_cycle_data.size(); i++)
	{
		for(auto it = inverse_index_maps_cycle_data[i].cbegin(); it != inverse_index_maps_cycle_data[i].cend(); it++)
			ofile<<"{"<<(*it).first<<":"<<(*it).second<<"}"<<std::endl;
		ofile<<std::endl;
	}
	if(this->phenotype.cycle.data.pCycle_Model != NULL)
		ofile<<"pCycle_Model is NOT NULL"<<std::endl;
	ofile<<"time_units:"<<this->phenotype.cycle.data.time_units<<std::endl;
 
	for(int i=0; i<this->phenotype.cycle.data.transition_rates.size();i++)
		{
			ofile<<"Transition Rates Vector:"<<i<<std::endl;
			for(int j=0; j<this->phenotype.cycle.data.transition_rates[i].size();j++)
				ofile<<"rate "<<j<<":"<<this->phenotype.cycle.data.transition_rates[i][j]<<std::endl;
		}
	
	ofile<<"current_phase_index:"<<this->phenotype.cycle.data.current_phase_index<<std::endl;
	ofile<<"elapsed_time_in_phase:"<<this->phenotype.cycle.data.elapsed_time_in_phase<<std::endl;

	ofile<<"=> class Cell { class Phenotype { class Death {"<<std::endl;
	for(int i=0; i < this->phenotype.death.rates.size(); i++)
		ofile<<"rates "<<i<<":"<<this->phenotype.death.rates[i]<<std::endl;
	ofile<<"Checking if std::vector<Cycle_Model*> models contains NULL entries"<<std::endl;
	for(int i=0; i< this->phenotype.death.models.size(); i++)
		if(this->phenotype.death.models[i] == NULL)
			ofile<<"models["<<i<<"]"<<" is NULL"<<std::endl;
		else
			ofile<<"models["<<i<<"]"<<" is NOT NULL"<<std::endl;

	ofile<<"=> class Cell { class Phenotype { class Death { class Death_Parameters { "<<std::endl;
	for(int i=0; i<this->phenotype.death.parameters.size(); i++)
		{
			ofile<<i<<":"<<std::endl;
			ofile<<"time_units:"<<this->phenotype.death.parameters[i].time_units<<std::endl;
			ofile<<"unlysed_fluid_change_rate:"<<this->phenotype.death.parameters[i].unlysed_fluid_change_rate<<std::endl;
			ofile<<"lysed_fluid_change_rate:"<<this->phenotype.death.parameters[i].lysed_fluid_change_rate<<std::endl;
			ofile<<"cytoplasmic_biomass_change_rate:"<<this->phenotype.death.parameters[i].cytoplasmic_biomass_change_rate<<std::endl;
			ofile<<"nuclear_biomass_change_rate:"<<this->phenotype.death.parameters[i].nuclear_biomass_change_rate<<std::endl;
			ofile<<"calcification rate:"<<this->phenotype.death.parameters[i].calcification_rate<<std::endl;
			ofile<<"relative_rupture_volume:"<<this->phenotype.death.parameters[i].relative_rupture_volume<<std::endl;
		}

	ofile<<"=> class Cell { class Phenotype { class Death {"<<std::endl;
	ofile<<"dead:"<<this->phenotype.death.dead<<std::endl;
	ofile<<"current_death_model_index:"<<this->phenotype.death.current_death_model_index<<std::endl;

	ofile<<"=> class Cell { class Phenotype { class Volume { "<<std::endl;
	ofile<<"Volume Data"<<std::endl;
	ofile<<"total:"<<this->phenotype.volume.total<<std::endl;
	ofile<<"solid:"<<this->phenotype.volume.solid<<std::endl;
	ofile<<"fluid:"<<this->phenotype.volume.fluid<<std::endl;
	ofile<<"fluid_fraction:"<<this->phenotype.volume.fluid_fraction<<std::endl;
	ofile<<"nuclear:"<<this->phenotype.volume.nuclear<<std::endl;
	ofile<<"nuclear_fluid:"<<this->phenotype.volume.nuclear_fluid<<std::endl;
	ofile<<"nuclear_solid:"<<this->phenotype.volume.nuclear_solid<<std::endl;
	ofile<<"cytoplasmic:"<<this->phenotype.volume.cytoplasmic<<std::endl;
	ofile<<"cytoplasmic_fluid:"<<this->phenotype.volume.cytoplasmic_fluid<<std::endl;
	ofile<<"cytoplasmic_solid:"<<this->phenotype.volume.cytoplasmic_solid<<std::endl;
	ofile<<"calcified_fraction:"<<this->phenotype.volume.calcified_fraction<<std::endl;
	ofile<<"cytoplasmic_to_nuclear_ratio:"<<this->phenotype.volume.cytoplasmic_to_nuclear_ratio<<std::endl;
	ofile<<"rupture_volume:"<<this->phenotype.volume.rupture_volume<<std::endl;
	ofile<<"cytoplasmic_biomass_change_rate:"<<this->phenotype.volume.cytoplasmic_biomass_change_rate<<std::endl;
	ofile<<"nuclear_biomass_change_rate:"<<this->phenotype.volume.nuclear_biomass_change_rate<<std::endl;
	ofile<<"fluid_change_rate:"<<this->phenotype.volume.fluid_change_rate<<std::endl;
	ofile<<"calcification_rate:"<<this->phenotype.volume.calcification_rate<<std::endl;
	ofile<<"target_solid_cytoplasmic:"<<this->phenotype.volume.target_solid_cytoplasmic<<std::endl;
	ofile<<"target_solid_nuclear:"<<this->phenotype.volume.target_solid_nuclear<<std::endl;
	ofile<<"target_fluid_fraction:"<<this->phenotype.volume.target_fluid_fraction<<std::endl;
	ofile<<"target_cytoplasmic_to_nuclear_ratio:"<<this->phenotype.volume.target_cytoplasmic_to_nuclear_ratio<<std::endl;
	ofile<<"relative_rupture_volume:"<<this->phenotype.volume.relative_rupture_volume<<std::endl;

	ofile<<"=> class Cell { class Phenotype { class  Geometry "<<std::endl;
	ofile<<"radius:"<<this->phenotype.geometry.radius<<std::endl;
	ofile<<"nuclear_radius:"<<this->phenotype.geometry.nuclear_radius<<std::endl;
	ofile<<"surface_area:"<<this->phenotype.geometry.surface_area<<std::endl;
	ofile<<"polarity:"<<this->phenotype.geometry.polarity<<std::endl;

	ofile<<"=> class Cell { class Phenotype { class  Mechanics "<<std::endl;
	ofile<<"cell_adhesion_strength:"<<this->phenotype.mechanics.cell_cell_adhesion_strength<<std::endl;
	ofile<<"cell_BM_adhesion_strength:"<<this->phenotype.mechanics.cell_BM_adhesion_strength<<std::endl;
	ofile<<"cell_cell_repulsion_strength:"<<this->phenotype.mechanics.cell_cell_repulsion_strength<<std::endl;
	ofile<<"cell_BM_repulsion_strength:"<<this->phenotype.mechanics.cell_BM_repulsion_strength<<std::endl;
	ofile<<"relative_maximum_adhesion_distance:"<<this->phenotype.mechanics.relative_maximum_adhesion_distance<<std::endl;

	ofile<<"=> class Cell { class  Phenotype { class Motility { "<<std::endl;
	ofile<<"is_motile:"<<this->phenotype.motility.is_motile<<std::endl;
	ofile<<"persistence_time:"<<this->phenotype.motility.persistence_time<<std::endl;
	ofile<<"migration_speed:"<<this->phenotype.motility.migration_speed<<std::endl;
	for(int i=0; i<this->phenotype.motility.migration_bias_direction.size(); i++)
		ofile<<"migration bias direction "<<i<<":"<<this->phenotype.motility.migration_bias_direction[i]<<std::endl;
	ofile<<"migration_bias:"<<this->phenotype.motility.migration_bias<<std::endl;
	ofile<<"restrict_to_2D:"<<this->phenotype.motility.restrict_to_2D<<std::endl;
	for(int i=0; i<this->phenotype.motility.motility_vector.size(); i++)
		ofile<<"motility vector "<<i<<":"<<this->phenotype.motility.motility_vector[i]<<std::endl;
	ofile<<"chemotaxis_index:"<<this->phenotype.motility.chemotaxis_index<<std::endl;
	ofile<<"chemotaxis_direction:"<<this->phenotype.motility.chemotaxis_direction<<std::endl;

	ofile<<"=> class Cell { class Phenotype { class Secretion { "<<std::endl;
	for(int i=0; i<this->phenotype.secretion.secretion_rates.size();i++)
			ofile<<"secretion_rates "<<i<<":"<<this->phenotype.secretion.secretion_rates[i]<<std::endl;
	for(int i=0; i<this->phenotype.secretion.uptake_rates.size();i++)
			ofile<<"uptake_rates "<<i<<":"<<this->phenotype.secretion.uptake_rates[i]<<std::endl;
	for(int i=0; i<this->phenotype.secretion.saturation_densities.size();i++)
			ofile<<"saturation_densities "<<i<<":"<<this->phenotype.secretion.saturation_densities[i]<<std::endl;
	for(int i=0; i<this->phenotype.secretion.net_export_rates.size();i++)
			ofile<<"net_export_rates "<<i<<":"<<this->phenotype.secretion.net_export_rates[i]<<std::endl;

	ofile<<"=> class Cell { class Phenotype { class Molecular { "<<std::endl;
	for(int i=0; i<this->phenotype.molecular.internalized_total_substrates.size();i++)
			ofile<<"internalized_total_substrates "<<i<<":"<<this->phenotype.molecular.internalized_total_substrates[i]<<std::endl;
	for(int i=0; i<this->phenotype.molecular.fraction_released_at_death.size();i++)
			ofile<<"fraction_released_at_death "<<i<<":"<<this->phenotype.molecular.fraction_released_at_death[i]<<std::endl;
	for(int i=0; i<this->phenotype.molecular.fraction_transferred_when_ingested.size();i++)
			ofile<<"fraction_transferred_when_ingested "<<i<<":"<<this->phenotype.molecular.fraction_transferred_when_ingested[i]<<std::endl;


if (this->phenotype.intracellular != NULL) {
		ofile<<"=> class Cell { class Phenotype { class Intracellular { "<<std::endl;
#ifdef ADDON_PHYSIBOSS
		if (this->phenotype.intracellular->intracellular_type.compare("maboss") == 0) {
			
			MaBoSSIntracellular* t_intracellular = static_cast<MaBoSSIntracellular*>(this->phenotype.intracellular); 

			ofile<<"bnd_file:"<<t_intracellular->get_bnd_filename()<<std::endl; 
			ofile<<"cfg_file:"<<t_intracellular->get_cfg_filename()<<std::endl; 
			ofile<<"time_step:"<<t_intracellular->time_step<<std::endl;
			ofile<<"discrete_time:"<<t_intracellular->discrete_time<<std::endl;
			ofile<<"time_tick:"<<t_intracellular->time_tick<<std::endl;
			ofile<<"scaling:"<<t_intracellular->scaling<<std::endl;
			
			for (auto t_mutation: t_intracellular->mutations)
				ofile<<"mutation:"<<t_mutation.first<<"="<<t_mutation.second<<std::endl;	
			
			for (auto t_initial_value: t_intracellular->initial_values)
				ofile<<"initial_value:"<<t_initial_value.first<<"="<<t_initial_value.second<<std::endl;
			
			for (auto t_parameter: t_intracellular->parameters)
				ofile<<"parameter:"<<t_parameter.first<<"="<<t_parameter.second<<std::endl;
			
			ofile<<"next_run:"<<t_intracellular->next_physiboss_run<<std::endl;
			ofile<<"time_to_update:"<<t_intracellular->maboss.get_time_to_update()<<std::endl;
			ofile<<"state:"<<t_intracellular->get_state()<<std::endl;
		}
#endif	
	}


	ofile<<"=> class Cell { "<<std::endl;
	ofile<<"is_out_of_domain:"<<is_out_of_domain<<std::endl;
	ofile<<"is_movable:"<<is_movable<<std::endl;
	
	
	for(int i=0; i<displacement.size(); i++)
		ofile<<"Displacement["<<i<<"]:"<<displacement[i]<<std::endl;

	ofile<<"=> class Cell inherits from Basic_Agent - its members now"<<std::endl;
	ofile<<"=> class Basic_Agent { "<<std::endl;
	ofile<<"Volume:"<<this->get_total_volume()<<std::endl;
	ofile<<"volume_is_changed:"<<this->get_is_volume_changed()<<std::endl;

	ofile<<"cell_source_sink_solver_temp1:"<<std::endl;
	for(int i=0; i<cell_source_sink_solver_temp1.size();i++)
	 	ofile<<i<<":"<<cell_source_sink_solver_temp1[i]<<std::endl;

	ofile<<"cell_source_sink_solver_temp2:"<<std::endl;
	for(int i=0; i<cell_source_sink_solver_temp2.size();i++)
	 	ofile<<i<<":"<<cell_source_sink_solver_temp2[i]<<std::endl;

	ofile<<"cell_source_sink_solver_temp_export1:"<<std::endl;
	for(int i=0; i<cell_source_sink_solver_temp_export1.size();i++)
	 	ofile<<i<<":"<<cell_source_sink_solver_temp_export1[i]<<std::endl;

	ofile<<"cell_source_sink_solver_temp_export2:"<<std::endl;
	for(int i=0; i<cell_source_sink_solver_temp_export2.size();i++)
	 	ofile<<i<<":"<<cell_source_sink_solver_temp_export2[i]<<std::endl;

	ofile<<"Prev_Vel[0]:"<<previous_velocity[0]<<" Prev_Vel[1]:"<<previous_velocity[1]<<" Prev_Vel[2]:"<<previous_velocity[2]<<std::endl;
	ofile<<"is_active:"<<is_active<<std::endl;

	ofile<<"total_extracellular_substrate_change"<<std::endl;
	for(int i=0; i<total_extracellular_substrate_change.size();i++)
	 	ofile<<i<<":"<<total_extracellular_substrate_change[i]<<std::endl;

	ofile<<"Secretion Rates"<<std::endl;
	for(int i=0; i<secretion_rates->size();i++)
	 	ofile<<i<<":"<<secretion_rates->data()[i]<<std::endl;

	ofile<<"Saturation Densities"<<std::endl;
	for(int i=0; i<saturation_densities->size();i++)
	 	ofile<<i<<":"<<saturation_densities->data()[i]<<std::endl;


	ofile<<"Uptake Rates"<<std::endl;
	for(int i=0; i<uptake_rates->size();i++)
	 	ofile<<i<<":"<<uptake_rates->data()[i]<<std::endl;

	ofile<<"Net Export Rates"<<std::endl;
	for(int i=0; i<net_export_rates->size();i++)
	 	ofile<<i<<":"<<net_export_rates->data()[i]<<std::endl;


	ofile<<"Internalized Substrates"<<std::endl;
	for(int i=0; i<internalized_substrates->size();i++)
	 	ofile<<i<<":"<<internalized_substrates->data()[i]<<std::endl;


	ofile<<"Fraction Released At Death"<<std::endl;
	for(int i=0; i<fraction_released_at_death->size();i++)
	 	ofile<<i<<":"<<fraction_released_at_death->data()[i]<<std::endl;


	ofile<<"Fraction Transferred When Ingested;"<<std::endl;
	for(int i=0; i<fraction_transferred_when_ingested->size();i++)
	 	ofile<<i<<":"<<fraction_transferred_when_ingested->data()[i]<<std::endl;

	ofile<<"Vel[0]:"<< velocity[0]<<" Vel[1]:"<<velocity[1]<<" Vel[2]:"<<velocity[2]<<std::endl;

	ofile.close();

}

void Cell::remove_self_from_all_neighbors( void )
{
	Cell* pCell = this; 
	// go through all neighbors (pN) of this (pC)

	for( int j = 0 ; j < pCell->state.neighbors.size(); j++ )
	{
	 	Cell* pN = pCell->state.neighbors[j]; 

		// for each pN, remove pC from list of neighbors 
			// find pC in neighbors 


			auto SearchResult = std::find( 
				pN->state.neighbors.begin(),pN->state.neighbors.end(),pCell );  		

			// if pC is indeed found, remove it  
			// erase pC from neighbors 
			if( SearchResult != pN->state.neighbors.end() )
			{
				// if the target is found, set the appropriate rate 
				pN->state.neighbors.erase( SearchResult ); 
			}
			else
			{ /* future error message */  }
	}

	return; 
}


bool cell_definitions_by_name_constructed = false;

void prebuild_cell_definition_index_maps( mpi_Environment &world, mpi_Cartesian &cart_topo) //1.14 revised
{
	// look in config file 
	extern pugi::xml_node physicell_config_root; 
	// find the start of cell definitions 
	pugi::xml_node node = physicell_config_root.child( "cell_definitions" ); 

	// find the first cell definition 
	node = node.child( "cell_definition" ); 

	// We won't declare and build the cell definitions just yet. 
	// All we want to do is know in advance the IDs and names 
	// of all cell definitions, so we can appropriately size 
	// Cell_Interactions and Cell_Transformations in the phenotype 
	// when we set up the cell definitions. 

	int n = 0; 
	while( node )
	{
		int ID = node.attribute( "ID" ).as_int();  
		std::string type_name = node.attribute( "name" ).value();   

		if (IOProcessor(world))
			std::cout << "Pre-processing type " << ID << " named " << type_name << std::endl; 

//		cell_definitions_by_name[ type_name ] = pCD; 
//		cell_definitions_by_type[ pCD->type ] = pCD; 

		cell_definition_indices_by_name[ type_name ] = n; 
		cell_definition_indices_by_type[ ID ] = n; 

		node = node.next_sibling( "cell_definition" ); 
		n++; 
	}	

	return; 
}

void build_cell_definitions_maps( void ) //1.14 revised
{
//	cell_definitions_by_name.
//	cell_definitions_by_index

	for( int n=0; n < cell_definitions_by_index.size() ; n++ )
	{
		Cell_Definition* pCD = cell_definitions_by_index[n];
		cell_definitions_by_name[ pCD->name ] = pCD;
		cell_definitions_by_type[ pCD->type ] = pCD;

		cell_definition_indices_by_name[ pCD->name ] = n; 
		cell_definition_indices_by_type[ pCD->type ] = n; 
	}

	cell_definitions_by_name_constructed = true;

	return;
}

void display_ptr_as_bool( void (*ptr)(Cell*,Phenotype&,double), std::ostream& os )
{
	if( ptr )
	{ os << "true"; return; }
	os << "false";
	return;
}

void display_ptr_as_bool( void (*ptr)(Cell*,Phenotype&,Cell*,Phenotype&,double), std::ostream& os )
{
	if( ptr )
	{ os << "true"; return; }
	os << "false"; 
	return;
}

void display_cell_definitions( std::ostream& os, mpi_Environment &world, mpi_Cartesian &cart_topo )
{
	if (IOProcessor(world)) {
		for( int n=0; n < cell_definitions_by_index.size() ; n++ )
		{
			Cell_Definition* pCD = cell_definitions_by_index[n]; 
			os << n << " :: type:" << pCD->type << " name: " << pCD->name << std::endl; 

			// summarize cycle model 
			if( pCD->phenotype.cycle.pCycle_Model != NULL )
			{
				os << "\t cycle model: " << pCD->phenotype.cycle.model().name  
					<< " (code=" << pCD->phenotype.cycle.model().code << ")" << std::endl; 
					
				// let's show the transition rates 
				Cycle_Model* pCM = &(pCD->phenotype.cycle.model() ); 
				Cycle_Data* pCMD = &(pCD->phenotype.cycle.data ); 
				for( int n=0 ; n < pCM->phases.size() ; n++ )
				{
					os << "\t\tPhase " << n << ": " << pCM->phases[n].name << std::endl; 
				}
				os << "\t\tCycle transitions: " << std::endl 
				<< "\t\t-----------------" << std::endl; 
				for( int n=0 ; n < pCM->phase_links.size() ; n++ )
				{
					for( int k=0; k < pCM->phase_links[n].size() ; k++ )
					{
						int start = pCM->phase_links[n][k].start_phase_index;
						int end = pCM->phase_links[n][k].end_phase_index; 
						os << "\t\t" << pCM->phases[start].name << " --> " 
							<< pCM->phases[end].name << " w mean duration " 
							<< 1.0 / pCMD->transition_rate( start,end) << " min" << std::endl; 
					}
				}			
				
			}
			else
			{ 	os << "\t cycle model not initialized" << std::endl; } 

			// summarize death models 
			os << "\t death models: " << std::endl; 
			for( int k=0 ; k < pCD->phenotype.death.models.size(); k++ )
			{
				os << "\t\t" << k << " : " << pCD->phenotype.death.models[k]->name 
				<< " (code=" << pCD->phenotype.death.models[k]->code << ")" 
				<< " with rate " << pCD->phenotype.death.rates[k] << " 1/min" << std::endl; 

				Cycle_Model* pCM = (pCD->phenotype.death.models[k] ); 
				Cycle_Data* pCMD = &(pCD->phenotype.death.models[k]->data ); 

				
				os << "\t\tdeath phase transitions: " << std::endl 
				<< "\t\t------------------------" << std::endl; 
				for( int n=0 ; n < pCM->phase_links.size() ; n++ )
				{
					for( int k=0; k < pCM->phase_links[n].size() ; k++ )
					{
						int start = pCM->phase_links[n][k].start_phase_index;
						int end = pCM->phase_links[n][k].end_phase_index; 
						os << "\t\t" << pCM->phases[start].name << " --> " 
							<< pCM->phases[end].name << " w mean duration " 
							<< 1.0 / pCMD->transition_rate( start,end) << " min" << std::endl; 
					}
				}			
				
				
				
			}
			
			// summarize functions 
			Cell_Functions* pCF = &(pCD->functions); 
			os << "\t key functions: " << std::endl; 
			os << "\t\t migration bias rule: "; display_ptr_as_bool( pCF->update_migration_bias , std::cout ); 
			os << std::endl; 
			os << "\t\t custom rule: "; display_ptr_as_bool( pCF->custom_cell_rule , std::cout ); 
			os << std::endl; 
			os << "\t\t phenotype rule: "; display_ptr_as_bool( pCF->update_phenotype , std::cout ); 
			os << std::endl; 
			os << "\t\t volume update function: "; display_ptr_as_bool( pCF->volume_update_function , std::cout ); 
			os << std::endl; 
			os << "\t\t mechanics function: "; display_ptr_as_bool( pCF->update_velocity , std::cout ); 
			os << std::endl;
			os << "\t\t contact function: "; display_ptr_as_bool( pCF->contact_function , std::cout ); 
			os << std::endl; 
			
			// summarize motility 
			
			Motility* pM = &(pCD->phenotype.motility); 
			std::string val = "true";
			if( pM->is_motile == false )
			{ val = "false"; } 
		
			std::string dimen = "3D"; 
			if( pM->restrict_to_2D == true )
			{ dimen = "2D"; } 

			os << "\tmotility (enabled: " << val << " in " << dimen << ")" << std::endl 
				<< "\t\tspeed: " << pM->migration_speed << " micron/min" << std::endl 
				<< "\t\tbias: " << pM->migration_bias << " " << std::endl 
				<< "\t\tpersistence time: " << pM->persistence_time << " min" << std::endl 
				<< "\t\tchemotaxis (enabled: ";
				
				val = "true" ;
				if( pCD->functions.update_migration_bias != chemotaxis_function )
				{ val = "false"; } 
			os << val << ")" << std::endl 
				<< "\t\t\talong " 
				<< pM->chemotaxis_direction << " * grad(" 
				<< microenvironment.density_names[ pM->chemotaxis_index ] << ") " << std::endl; 
				
			// secretion
			
			
			
			// mechanics
			
			// size 
		
			
			Custom_Cell_Data* pCCD = &(pCD->custom_data); 
			os << "\tcustom data: " << std::endl; 
			for( int k=0; k < pCCD->variables.size(); k++)
			{
				os << "\t\t" << pCCD->variables[k] << std::endl; 
			}
			os << "\tcustom vector data: " << std::endl; 
			for( int k=0; k < pCCD->vector_variables.size(); k++)
			{
				os << "\t\t" << pCCD->vector_variables[k] << std::endl; 
			}
			os << "\t\t\tNOTE: custom vector data will eventually be merged with custom data" << std::endl; 
				
		}
	}
	return; 
}

Cell_Definition* find_cell_definition( std::string search_string ) //1.14 ready
{
	// if the registry isn't built yet, then do it!
	if( cell_definitions_by_name_constructed == false )
	{
		build_cell_definitions_maps();
	}

	Cell_Definition* output = NULL;
	if( cell_definitions_by_name.count( search_string ) > 0 )
	{
		output = cell_definitions_by_name.find( search_string )->second;
	}

	if( output == NULL )
	{
		std::cout << "Warning! Cell_Definition for " << search_string << " not found!" << std::endl;
	}

	return output;
}

Cell_Definition* find_cell_definition( int search_type )
{
	// if the registry isn't built yet, then do it!
	if( cell_definitions_by_name_constructed == false )
	{
		build_cell_definitions_maps();
	}

	Cell_Definition* output = NULL;
	if( cell_definitions_by_type.count( search_type ) > 0 )
	{
		output = cell_definitions_by_type.find( search_type )->second;
	}

	if( output == NULL )
	{
		std::cout << "Warning! Cell_Definition for " << search_type << " not found!" << std::endl;
	}

	return output;
}

Cell_Definition& get_cell_definition( std::string search_string )
{
	// if the registry isn't built yet, then do it!
	if( cell_definitions_by_name_constructed == false )
	{
		build_cell_definitions_maps();
	}

	if( cell_definitions_by_name.count( search_string ) > 0 )
	{
		return *(cell_definitions_by_name.find( search_string )->second );
	}

	std::cout << "Warning! Cell_Definition for " << search_string << " not found!" << std::endl;
	std::cout << "Returning the default cell definition instead ... " << std::endl;

	return cell_defaults;
}

Cell_Definition& get_cell_definition( int search_type )
{
	// if the registry isn't built yet, then do it!
	if( cell_definitions_by_name_constructed == false )
	{
		build_cell_definitions_maps();
	}

	if( cell_definitions_by_type.count( search_type ) > 0 )
	{
		return *(cell_definitions_by_type.find( search_type )->second);
	}

	std::cout << "Warning! Cell_Definition for " << search_type << " not found!" << std::endl;
	std::cout << "Returning the default cell definition instead ... " << std::endl;

	return cell_defaults;
}


Cell_Definition* initialize_cell_definition_from_pugixml( pugi::xml_node cd_node, mpi_Environment &world, mpi_Cartesian &cart_topo) //1.14 revised
{
	Cell_Definition* pCD; 
	
	// if this is not "default" then create a new one 
	if( std::strcmp( cd_node.attribute( "name" ).value() , "default" ) != 0 
	  && std::strcmp( cd_node.attribute( "ID" ).value() , "0" ) != 0 )
	{ pCD = new Cell_Definition; }
	else
	{ pCD = &cell_defaults; }
	
	// set the name 
	pCD->name = cd_node.attribute("name").value();
	
	// set the ID 
	if( cd_node.attribute("ID" ) )
	{ pCD->type = cd_node.attribute("ID").as_int(); }
	else
	{ pCD->type = -1; } 

	// get the parent definition (if any) 
	Cell_Definition* pParent = NULL;
	if( cd_node.attribute( "parent_type" ) )
	{ pParent = find_cell_definition( cd_node.attribute( "parent_type" ).value() ); }
	// if it's not the default and no parent stated, inherit from default 
	bool use_default_as_parent_without_specifying = false; 
	if( pParent == NULL && pCD != &cell_defaults )
	{
		pParent = &cell_defaults; 
		use_default_as_parent_without_specifying = true; 
	}  

	// if we found something to inherit from, then do it! 
	if( pParent != NULL )
	{
		//std::cout << "\tCopying from type " << pParent->name << " ... " << std::endl; 
		*pCD = *pParent; 
		
		// but recover the name and ID (type)
		pCD->name = cd_node.attribute("name").value();
		pCD->type = cd_node.attribute("ID").as_int(); 
	} 
	
	/* bugfix on April 24, 2022 */ 
	// If we copied from cell_defaults and also wrote 
	// more properties to cell defaults, we have messed up 
	// some of the rates that are assumed to start at zero. 
	// So, let's overwrite with zeros. 
	if( use_default_as_parent_without_specifying )
	{

		pugi::xml_node node_options = xml_find_node( physicell_config_root , "options" ); 
		bool disable_bugfix = false; 
		if( node_options )
		{ disable_bugfix = xml_get_bool_value( node_options, "legacy_cell_defaults_copy" ); }

		if( disable_bugfix == false )
		{
			int number_of_substrates = microenvironment.density_names.size(); 
			int number_of_cell_defs = cell_definition_indices_by_name.size(); 

			// motility 
			pCD->phenotype.motility.is_motile = false; 
			pCD->phenotype.motility.chemotactic_sensitivities.assign(number_of_substrates,0.0); 
			pCD->functions.update_migration_bias = NULL; 

			// secretion  
			pCD->phenotype.secretion.secretion_rates.assign(number_of_substrates,0.0); 
			pCD->phenotype.secretion.uptake_rates.assign(number_of_substrates,0.0); 
			pCD->phenotype.secretion.net_export_rates.assign(number_of_substrates,0.0); 
			pCD->phenotype.secretion.saturation_densities.assign(number_of_substrates,0.0); 

			// interaction 
			pCD->phenotype.cell_interactions.apoptotic_phagocytosis_rate = 0.0; 
			pCD->phenotype.cell_interactions.necrotic_phagocytosis_rate = 0.0; 
			pCD->phenotype.cell_interactions.other_dead_phagocytosis_rate = 0.0; 
			pCD->phenotype.cell_interactions.live_phagocytosis_rates.assign(number_of_cell_defs,0.0); 
			pCD->phenotype.cell_interactions.attack_rates.assign(number_of_cell_defs,0.0); 
			pCD->phenotype.cell_interactions.attack_damage_rate = 1.0; 
			pCD->phenotype.cell_interactions.fusion_rates.assign(number_of_cell_defs,0.0); 

			// transformation 
			pCD->phenotype.cell_transformations.transformation_rates.assign(number_of_cell_defs,0.0); 

			// asymmetric division
			pCD->phenotype.cycle.asymmetric_division.asymmetric_division_probabilities.assign(number_of_cell_defs,0.0);
		}
		else 
		{
			std::cout << "[Rank " << world.rank << "] Warning! You have disabled a bugfix on cell definition inheritance" << std::endl 
			<< "\tBe VERY careful that you have manually specified every parameter value" << std::endl
			<< "\tfor every cell cell definition. Set legacy_cell_defaults_copy to false" << std::endl 
			<< "\tin the options section of your parameter file to re-enable the bug fix. " << std::endl << std::endl
			<< "\tSome good news: if you used the model builder, this has never affected your results." << std::endl; 
		}

	}

	// sync to microenvironment
	pCD->pMicroenvironment = NULL;
	if( BioFVM::get_default_microenvironment() != NULL )
	{ pCD->pMicroenvironment = BioFVM::get_default_microenvironment(); }

	// figure out if this ought to be 2D
	if( default_microenvironment_options.simulate_2D )
	{
		std::cout << "[Rank " << world.rank << "] Note: setting cell definition to 2D based on microenvironment domain settings ... "
		<< std::endl; 
		pCD->functions.set_orientation = up_orientation; 
		pCD->phenotype.geometry.polarity = 1.0; 
		pCD->phenotype.motility.restrict_to_2D = true; 
	}
	
	// make sure phenotype.secretions are correctly sized 
	
	// pCD->phenotype.secretion.sync_to_current_microenvironment();
	pCD->phenotype.secretion.sync_to_microenvironment( (pCD->pMicroenvironment) ); 
	pCD->phenotype.molecular.sync_to_microenvironment( (pCD->pMicroenvironment) );

	// make sure that interactions and transfomations are correctly sized 
	// this requires that prebuild_cell_definition_index_maps was already run 
	pCD->phenotype.cell_interactions.sync_to_cell_definitions(); 
	pCD->phenotype.cell_transformations.sync_to_cell_definitions();
	pCD->phenotype.cycle.asymmetric_division.sync_to_cell_definitions();
	pCD->phenotype.mechanics.sync_to_cell_definitions();
	
	
	// set the reference phenotype 
	pCD->parameters.pReference_live_phenotype = &(pCD->phenotype); 
	pugi::xml_node node = cd_node.child( "phenotype" ); 

	// set up the cell cycle 
	// make sure the standard cycle models are defined 
	create_standard_cycle_and_death_models();

	node = node.child( "cycle" ); 
	if( node )
	{
		int model; // = node.attribute("code").as_int() ; 
		if( strlen( node.attribute("code").as_string() ) > 0 )
		{ model = node.attribute("code").as_int(); }
		else
		{ model = find_cycle_model_code( node.attribute("name").as_string() ); } 
		if( model < 0 )
		{ 
			std::cout << "[Rank " << world.rank << "] Error. Unable to identify cycle model " 
				<< node.attribute("name").value() 
				<< " (" << node.attribute("code").value() << ")" << std::endl;
			exit(-1);
		}		
		
		// Set the model, but only if it was specified. 
		if( strlen( node.attribute("code").value() ) > 0 )
		{
			// set the model 
			if (model == PhysiCell_constants::advanced_Ki67_cycle_model)
            {
                pCD->functions.cycle_model = Ki67_advanced; 
            }
			else if (model == PhysiCell_constants::basic_Ki67_cycle_model)
            {
                pCD->functions.cycle_model = Ki67_basic; 
            }
            else if (model == PhysiCell_constants::flow_cytometry_cycle_model) 
            {
                pCD->functions.cycle_model = flow_cytometry_cycle_model;  
            }
            else if (model == PhysiCell_constants::live_apoptotic_cycle_model) // ?
            {
                pCD->functions.cycle_model = live;  // ?
                std::cout << "[Rank " << world.rank << "] Warning: live_apoptotic_cycle_model not directly supported." << std::endl		
                            << "         Substituting live cells model. Set death rates=0." << std::endl; 
            }
            else if (model == PhysiCell_constants::total_cells_cycle_model) 
            {
                pCD->functions.cycle_model = live; 
                std::cout << "[Rank " << world.rank << "] Warning: total_cells_cycle_model not directly supported." << std::endl		
                            << "         Substituting live cells model. Set death rates=0." << std::endl; 
            }
            else if (model == PhysiCell_constants::live_cells_cycle_model) 
            {
                pCD->functions.cycle_model = live; 
            }
            else if (model == PhysiCell_constants::flow_cytometry_separated_cycle_model) 
            {
                pCD->functions.cycle_model = flow_cytometry_separated_cycle_model; 
            }
            else if (model == PhysiCell_constants::cycling_quiescent_model) 
            {
                pCD->functions.cycle_model = cycling_quiescent; 
            }
            else
            {
                std::cout << "[Rank " << world.rank << "] Warning: Unknown cycle model " << std::endl;
                exit(-1); 
            }
			pCD->phenotype.cycle.sync_to_cycle_model( pCD->functions.cycle_model ); 
		}
		
		// now, if we inherited from another cell, AND 
		// if that parent type has the same cylce model, 
		// then overwrite with their transition rates 
		
		if( pParent != NULL )
		{
			if( pCD->phenotype.cycle.model().code == pParent->phenotype.cycle.model().code )
			{
				if (IOProcessor(world)) {
					std::cout << "copying data ... " << std::endl; 
					std::cout<< pParent->name << " to " << pCD->name << std::endl;
				} 
				pCD->phenotype.cycle.data = pParent->phenotype.cycle.data; 
			}
		}
		
		// set the transition rates 
		if( node.child( "phase_transition_rates" ) )
		{ node = node.child( "phase_transition_rates" ); }
		if( node.child( "transition_rates" ) )
		{
			node = node.child( "transition_rates" );
			if (IOProcessor(world)) 
				std::cout << "Warning: " << node.name() << " is deprecated. Use cycle.phase_transition_rates." 
				<< std::endl; 
		}
		if( node )
		{
			node = node.child( "rate");
			while( node )
			{
				// which rate 
				int start = node.attribute("start_index").as_int(); 
				int end = node.attribute("end_index").as_int(); 
				// fixed duration? 
				bool fixed = false; 
				if( node.attribute( "fixed_duration" ) )
				{ fixed = node.attribute("fixed_duration").as_bool(); }
				// actual value of transition rate 
				double value = xml_get_my_double_value( node ); 
				
				// set the transition rate 
				pCD->phenotype.cycle.data.transition_rate(start,end) = value; 
				// set it to fixed / non-fixed 
				pCD->phenotype.cycle.model().phase_link(start,end).fixed_duration = fixed; 
				
				node = node.next_sibling( "rate" ); 
			}
		}
		
		node = cd_node.child( "phenotype" );
		node = node.child( "cycle" ); 
		// Check for phase durations (as an alternative to transition rates)
		if( node.child( "phase_durations" ) )
		{ node = node.child( "phase_durations" ); }
		if( node )
		{
			node = node.child( "duration");
			while( node )
			{
				// which duration? 
				int start = node.attribute("index").as_int(); 
				// fixed duration? 
				bool fixed = false; 
				if( node.attribute( "fixed_duration" ) )
				{ fixed = node.attribute("fixed_duration").as_bool(); }
				// actual value of the duration 
				double value = xml_get_my_double_value( node ); 
				
				// set the transition rate 
				pCD->phenotype.cycle.data.exit_rate(start) = 1.0 / (value+1e-16); 
				// set it to fixed / non-fixed 
				pCD->phenotype.cycle.model().phase_links[start][0].fixed_duration = fixed; 
				
				node = node.next_sibling( "duration" ); 
			}
			
		}

		node = cd_node.child( "phenotype" );
		node = node.child( "cycle" );
		node = node.child( "standard_asymmetric_division" );
		if( node && node.attribute("enabled").as_bool() )
		{
			Asymmetric_Division *pAD = &(pCD->phenotype.cycle.asymmetric_division);

			// asymmetric division rates
			pugi::xml_node node_adp = node.child( "asymmetric_division_probability");
			while (node_adp)
			{
				// get the name of the target cell type
				std::string target_name = node_adp.attribute("name").value();
				// now find its index
				auto search = cell_definition_indices_by_name.find(target_name);
				// safety first!
				if( search != cell_definition_indices_by_name.end() )
				{
					// if the target is found, set the appropriate rate
					int target_index = search->second;

					double asymmetric_division_probability = xml_get_my_double_value(node_adp);
					pAD->asymmetric_division_probabilities[target_index] = asymmetric_division_probability;
				}
				else
				{
					
					std::cout << "[Rank " << world.rank << "] Error: When processing the " << pCD->name << " cell definition: " << std::endl
						<< "\tCould not find cell type " << target_name << " for asymmetric division." << std::endl
						<< "\tRemove this cell type from the asymmetric division probabilities!" << std::endl << std::endl;
					exit(-1);
				}
				node_adp = node_adp.next_sibling("asymmetric_division_probability");
			}
			if (IOProcessor(world))
				std::cout << "Asymmetric division probabilities for " << pCD->name << ": ";
			for (int i = 0; i < pAD->asymmetric_division_probabilities.size(); i++)
			{
				std::cout << pAD->asymmetric_division_probabilities[i] << " ";
			}
			std::cout << std::endl;
			pCD->functions.cell_division_function = standard_asymmetric_division_function;
			}
		}
		
	
	
	// here's what it ***should*** do: 
	// parse the model, get its code 
	// look for that model 
	// if the model is not yet there, then add it
	// otherwise, modify properties of that model 
	
	// set up the death models 
//	int death_model_index = 0; 
	node = cd_node.child( "phenotype" );
	node = node.child( "death" ); 
	if( node )
	{
		pugi::xml_node model_node = node.child( "model" );
		while( model_node )
		{
			node = model_node;

			int model; // = node.attribute("code").as_int() ; 
			if( strlen( node.attribute("code").as_string() ) > 0 )
			{ model = node.attribute("code").as_int(); }
			else
			{ model = find_cycle_model_code( node.attribute("name").as_string() ); } 
			if( model < 0 )
			{ 
				std::cout << "[Rank " << world.rank << "] Error. Unable to identify death model " 
					<< node.attribute("name").value() 
					<< " (" << node.attribute("code").value() << ")" << std::endl;
				exit(-1);
			}		

			
			// check: is that death model already there? 
			
			Death* pD = &( pCD->phenotype.death ); 
			int death_index = pD->find_death_model_index( model );
			bool death_model_already_exists = false; 
			if( pD->rates.size() > death_index )
			{
				if( pD->models[death_index]->code == model )
				{ death_model_already_exists = true; } 
			}
			
			// add the death model and its death rate 

			if( node.child( "death_rate" ) )
			{
				node = node.child( "death_rate" ); 
			}
			if( node.child( "rate" ) )
			{
				node = node.child( "rate" ); 
				std::cout << "Warning: " << node.name() << " is deprecated. Use death.model.death_rate." 
					<< std::endl; 
			}
			double rate = xml_get_my_double_value(node);
			node = node.parent();
			
			// get death model parameters 
			
			Death_Parameters death_params; 
			// if there is a parent and we already found this model, 
			// start with the inherited parameters 
			if( death_model_already_exists && pParent != NULL )
			{
				death_params = pParent->phenotype.death.parameters[death_index]; 
			}

			if( node.child("parameters" ) )
			{
				node = node.child( "parameters" );
			
				// only read these parameters if they are specified. 
				
				pugi::xml_node node_temp = node.child( "unlysed_fluid_change_rate" );
				if( node_temp )
				{ death_params.unlysed_fluid_change_rate = xml_get_my_double_value( node_temp ); }

				node_temp = node.child( "lysed_fluid_change_rate" );
				if( node_temp )
				{ death_params.lysed_fluid_change_rate = xml_get_my_double_value( node_temp ); }
			
				node_temp = node.child( "cytoplasmic_biomass_change_rate" );
				if( node_temp )
				{ death_params.cytoplasmic_biomass_change_rate = xml_get_my_double_value( node_temp ); }

				node_temp = node.child( "nuclear_biomass_change_rate" );
				if( node_temp )
				{ death_params.nuclear_biomass_change_rate = xml_get_my_double_value( node_temp ); }

				node_temp = node.child( "calcification_rate" );
				if( node_temp )
				{ death_params.calcification_rate = xml_get_my_double_value( node_temp ); }

				node_temp = node.child( "relative_rupture_volume" );
				if( node_temp )
				{ death_params.relative_rupture_volume = xml_get_my_double_value( node_temp ); }

				node_temp = node.child( "lysed_fluid_change_rate" );
				if( node_temp )
				{ death_params.lysed_fluid_change_rate = xml_get_my_double_value( node_temp ); }

	//			death_params.time_units = 
	//				get_string_attribute_value( node, "unlysed_fluid_change_rate", "units" ); 
				
				node = node.parent(); 
			}
					
			// set the model 
			// if the model already exists, just overwrite the parameters 
            if (model == PhysiCell_constants::apoptosis_death_model) 
            {
//					pCD->phenotype.death.add_death_model( rate , &apoptosis , apoptosis_parameters );
                if( death_model_already_exists == false )
                {
                    pCD->phenotype.death.add_death_model( rate , &apoptosis , death_params ); 
                    death_index = pD->find_death_model_index( model );
                }
                else
                {
                    pCD->phenotype.death.parameters[death_index] = death_params; 
                    pCD->phenotype.death.rates[death_index] = rate; 
                }
            }
            else if (model == PhysiCell_constants::necrosis_death_model) 
            {
                // set necrosis parameters 
//					pCD->phenotype.death.add_death_model( rate , &necrosis , necrosis_parameters );
                if( death_model_already_exists == false )
                {
                    pCD->phenotype.death.add_death_model( rate , &necrosis , death_params ); 
                    death_index = pD->find_death_model_index( model );
                }
                else
                {
                    pCD->phenotype.death.parameters[death_index] = death_params; 
                    pCD->phenotype.death.rates[death_index] = rate; 						
                }
            }
            else if (model == PhysiCell_constants::autophagy_death_model) 
            {
                std::cout << "[Rank " << world.rank << "] Warning: autophagy_death_model not yet supported." << std::endl		
                            << "         Skipping this model." << std::endl; 
            }
            else 
            {
                std::cout << "[Rank " << world.rank << "] Warning: Unknown death model " << std::endl;
                exit(-1); 
            }
			
			// now get transition rates within the death model 
			// set the rates 
			// node = node.child( "transition_rates" );
			
			
			pugi::xml_node node_death_transitions = node.child( "phase_transition_rates" ); 
			if( node.child( "transition_rates" ) )
			{
				node_death_transitions = node.child("transition_rates");
				if (IOProcessor(world))
				std::cout << "Warning: " << node_death_transitions.name() 
					<< " is deprecated. Use death.model.phase_transition_rates." 
					<< std::endl; 				
			}
			
			
			if( node_death_transitions )
			{
				pugi::xml_node node1 = node_death_transitions.child( "rate");
				while( node1 )
				{
					// which rate 
					int start = node1.attribute("start_index").as_int(); 
					int end = node1.attribute("end_index").as_int(); 
					// fixed duration? 
					bool fixed = false; 
					if( node1.attribute( "fixed_duration" ) )
					{ fixed = node1.attribute("fixed_duration").as_bool(); }
					// actual value of transition rate 
					double value = xml_get_my_double_value( node1 ); 
					
					// set the transition rate 
					pCD->phenotype.death.models[death_index]->transition_rate(start,end) = value; 
					// set it to fixed / non-fixed 
					pCD->phenotype.death.models[death_index]->phase_link(start,end).fixed_duration = fixed; 
					
					node1 = node1.next_sibling( "rate" ); 
				}
			}	

			if( node.child( "phase_durations" ) )
			{
				node = node.child("phase_durations"); // phase durations
				node = node.child( "duration" ); // duration
				while( node )
				{
					// which duration? 
					int start = node.attribute("index").as_int(); 
					// fixed duration? 
					bool fixed = false; 
					if( node.attribute( "fixed_duration" ) )
					{ fixed = node.attribute("fixed_duration").as_bool(); }
					// actual value of the duration 
					double value = xml_get_my_double_value( node ); 
					
					// set the transition rate 
					pCD->phenotype.death.models[death_index]->data.exit_rate(start) 
						= 1.0 / (value+1e-16); 
					// set it to fixed / non-fixed 
					pCD->phenotype.death.models[death_index]->phase_links[start][0].fixed_duration 
						= fixed; 
					
					node = node.next_sibling( "duration" ); 
				}
									
				node = node.parent(); // phase_durations 
				node = node.parent(); // model 
			}				
			
			model_node = model_node.next_sibling( "model" );
		}
		
	}
	
	// volume 
	node = cd_node.child( "phenotype" );
	node = node.child( "volume" );
	if( node )
	{
		Volume* pV = &(pCD->phenotype.volume);
		
		pugi::xml_node node_vol = node.child( "total" );
		if( node_vol )
		{ pV->total = xml_get_my_double_value( node_vol ); }

		node_vol = node.child( "fluid_fraction" );
		if( node_vol )
		{ pV->fluid_fraction = xml_get_my_double_value( node_vol ); }
		
		node_vol = node.child( "nuclear" );
		if( node_vol )
		{ pV->nuclear = xml_get_my_double_value( node_vol ); }

		node_vol = node.child( "fluid_change_rate" );
		if( node_vol )
		{ pV->fluid_change_rate = xml_get_my_double_value( node_vol ); }

		node_vol = node.child( "cytoplasmic_biomass_change_rate" );
		if( node_vol )
		{ pV->cytoplasmic_biomass_change_rate = xml_get_my_double_value( node_vol ); }

		node_vol = node.child( "nuclear_biomass_change_rate" );
		if( node_vol )
		{ pV->nuclear_biomass_change_rate = xml_get_my_double_value( node_vol ); }

		node_vol = node.child( "calcified_fraction" );
		if( node_vol )
		{ pV->calcified_fraction = xml_get_my_double_value( node_vol ); }

		node_vol = node.child( "calcification_rate" );
		if( node_vol )
		{ pV->calcification_rate = xml_get_my_double_value( node_vol ); }
	
		node_vol = node.child( "relative_rupture_volume" );
		if( node_vol )
		{ pV->relative_rupture_volume = xml_get_my_double_value( node_vol ); }

		// set all the parameters to be self-consistent 
		
		pV->fluid = pV->fluid_fraction * pV->total; 
		pV->solid = pV->total-pV->fluid; 

		pV->nuclear_fluid = pV->fluid_fraction * pV->nuclear; 
		pV->nuclear_solid = pV->nuclear - pV->nuclear_fluid;

		pV->cytoplasmic = pV->total - pV->nuclear;
		pV->cytoplasmic_fluid = pV->fluid_fraction*pV->cytoplasmic; 
		pV->cytoplasmic_solid = pV->cytoplasmic - pV->cytoplasmic_fluid; 
		

		pV->target_solid_cytoplasmic = pV->cytoplasmic_solid;
		pV->target_solid_nuclear = pV->nuclear_solid;
		pV->target_fluid_fraction = pV->fluid_fraction;
		
		pV->cytoplasmic_to_nuclear_ratio = pV->cytoplasmic / ( 1e-16 + pV->nuclear);
		pV->target_cytoplasmic_to_nuclear_ratio = pV->cytoplasmic_to_nuclear_ratio; 
		
		pV->rupture_volume = pV->relative_rupture_volume * pV->total; // in volume units 
		
		// update the geometry (radius, etc.) for consistency 

    	pCD->phenotype.geometry.update( NULL, pCD->phenotype, 0.0 ); 
	}
	
	// mechanics 
	node = cd_node.child( "phenotype" );
	node = node.child( "mechanics" ); 
	if( node )
	{
		Mechanics* pM = &(pCD->phenotype.mechanics);
		
		pugi::xml_node node_mech = node.child( "cell_cell_adhesion_strength" );
		if( node_mech )
		{ pM->cell_cell_adhesion_strength = xml_get_my_double_value( node_mech ); }	

		node_mech = node.child( "cell_BM_adhesion_strength" );
		if( node_mech )
		{ pM->cell_BM_adhesion_strength = xml_get_my_double_value( node_mech ); }	
		
		node_mech = node.child( "cell_cell_repulsion_strength" );
		if( node_mech )
		{ pM->cell_cell_repulsion_strength = xml_get_my_double_value( node_mech ); }
		
		node_mech = node.child( "cell_BM_repulsion_strength" );
		if( node_mech )
		{ pM->cell_BM_repulsion_strength = xml_get_my_double_value( node_mech ); }

		node_mech = node.child( "relative_maximum_adhesion_distance" );
		if( node_mech )
		{ pM->relative_maximum_adhesion_distance = xml_get_my_double_value( node_mech ); }	

		// cell adhesion affinities 
		node_mech = node.child( "cell_adhesion_affinities" );
		if( node_mech )
		{
			node_mech = node_mech.child("cell_adhesion_affinity"); 
			while( node_mech )
			{
				std::string target = node_mech.attribute("name").value(); 
				double value = xml_get_my_double_value( node_mech ) ; 

				// find the target
				// if found, assign taht affinity 
				int ind = find_cell_definition_index( target ); 
				if( ind > -1 )
				{ pM->cell_adhesion_affinities[ind] = value; }
				else
				{ std::cout << "what?!?" << std::endl; }

				node_mech = node_mech.next_sibling( "cell_adhesion_affinity"); 
			}
		}	

		node_mech = node.child( "options" );
		if( node_mech )
		{
			pugi::xml_node node_mech1 = node_mech.child( "set_relative_equilibrium_distance" ); 
			if( node_mech1 )
			{
				if( node_mech1.attribute("enabled").as_bool() )
				{
					double temp = xml_get_my_double_value( node_mech1 ); 
					pM->set_relative_equilibrium_distance( temp ); 
				}
			}

			node_mech1 = node_mech.child( "set_absolute_equilibrium_distance" ); 
			if( node_mech1 )
			{
				if( node_mech1.attribute("enabled").as_bool() )
				{
					double temp = xml_get_my_double_value( node_mech1 ); 
					pM->set_absolute_equilibrium_distance( pCD->phenotype , temp ); 
				}
			}
		}

        node_mech = node.child( "attachment_elastic_constant" );
		if( node_mech )
		{ pM->attachment_elastic_constant = xml_get_my_double_value( node_mech ); }

		if (IOProcessor(world))
			std::cout << "  --------- attachment_elastic_constant = " << pM->attachment_elastic_constant << std::endl;

        node_mech = node.child( "attachment_rate" );
		if( node_mech )
		{ pM->attachment_rate = xml_get_my_double_value( node_mech ); }	

        node_mech = node.child( "detachment_rate" );
		if( node_mech )
		{ pM->detachment_rate = xml_get_my_double_value( node_mech ); }	

		node_mech = node.child( "maximum_number_of_attachments" );
		if ( node_mech )
		{ pM->maximum_number_of_attachments = xml_get_my_int_value( node_mech ); }
	}
	
	// motility 
	node = cd_node.child( "phenotype" );
	node = node.child( "motility" ); 
	if( node )
	{
		Motility* pMot = &(pCD->phenotype.motility);
		
		pugi::xml_node node_mot = node.child( "speed" );
		if( node_mot )
		{ pMot->migration_speed = xml_get_my_double_value( node_mot ); }	

		node_mot = node.child( "migration_bias" );
		if( node_mot )
		{ pMot->migration_bias = xml_get_my_double_value( node_mot ); }	

		node_mot = node.child( "persistence_time" );
		if( node_mot )
		{ pMot->persistence_time = xml_get_my_double_value( node_mot ); }	

		node_mot = node.child( "options" );
		if( node_mot )
		{
			// enable motility? 
			pugi::xml_node node_mot1 = node_mot.child( "enabled" ); 
			if( node_mot1 )
            {
                pMot->is_motile = xml_get_my_bool_value( node_mot1 ); 
            } 
			
			// restrict to 2D? 
			node_mot1 = node_mot.child( "use_2D" ); 
			if( node_mot1 )
            {
                pMot->restrict_to_2D = xml_get_my_bool_value( node_mot1 ); 
		    }
			
			if( default_microenvironment_options.simulate_2D && pMot->restrict_to_2D == false )
			{
				if (IOProcessor(world))
				std::cout << "Note: Overriding to set cell motility for " << pCD->name << " to 2D based on " 
						  << "microenvironment domain settings ... " << std::endl; 				
				pMot->restrict_to_2D = true; 
			}
			
			// automated chemotaxis setup 
			node_mot1 = node_mot.child( "chemotaxis" ); 
			if( node_mot1 )
			{
				// enabled? if so, set the standard chemotaxis function
				if( xml_get_bool_value( node_mot1, "enabled" ) )
				{
					pCD->functions.update_migration_bias = chemotaxis_function;
				}	
				
				// search for the right chemo index 
				
				std::string substrate_name = xml_get_string_value( node_mot1 , "substrate" ); 
				pMot->chemotaxis_index = microenvironment.find_density_index( substrate_name ); 
				if( pMot->chemotaxis_index < 0)
				{
					std::cout << __FUNCTION__ << ": Error: parsing phenotype:motility:options:chemotaxis:  invalid substrate" << std::endl; 
					std::cout << substrate_name << " was not found in the microenvironment. Please check for typos!" << std::endl << std::endl; 
					exit(-1); 
				}
				
				std::string actual_name = microenvironment.density_names[ pMot->chemotaxis_index ]; 
				
				// error check 
				if( std::strcmp( substrate_name.c_str() , actual_name.c_str() ) != 0 )
				{
					std::cout << "Error: attempted to set chemotaxis to \"" 
						<< substrate_name << "\", which was not found in the microenvironment." << std::endl 
					<< "       Please double-check your substrate name in the config file." << std::endl << std::endl; 
					exit(-1); 
				}
				
				// set the direction 
				
				pMot->chemotaxis_direction = xml_get_int_value( node_mot1 , "direction" ); 
				
				//std::cout << pMot->chemotaxis_direction << " * grad( " << actual_name << " )" << std::endl; 

			}

			// automated advanced chemotaxis setup 
			node_mot1 = node_mot.child( "advanced_chemotaxis" ); 
			if( node_mot1 )
			{
				// enabled? if so, set the standard chemotaxis function
				if( xml_get_bool_value( node_mot1, "enabled" ) )
				{
					if( pCD->functions.update_migration_bias == chemotaxis_function )
					{
						if (IOProcessor(world))
						std::cout << "Warning: when processing motility for " << pCD->name << " cells: " << std::endl 
								 << "\tBoth chemotaxis and advanced_chemotaxis are enabled." << std::endl
						          << "\tThe settings for advanced_chemotaxis override those of chemotaxis." << std::endl; 
					}
					pCD->functions.update_migration_bias = advanced_chemotaxis_function;
					if( xml_get_bool_value( node_mot1, "normalize_each_gradient" ) )
					{ pCD->functions.update_migration_bias = advanced_chemotaxis_function_normalized; }
				}	

				// now process the chemotactic sensitivities 

				pugi::xml_node node_cs = node_mot1.child( "chemotactic_sensitivities"); 
				if( node_cs  )
				{
					node_cs = node_cs.child("chemotactic_sensitivity"); 

					while( node_cs )
					{
						std::string substrate_name = node_cs.attribute( "substrate").value(); 
						int index = microenvironment.find_density_index( substrate_name ); 
						std::string actual_name = ""; 
						if( index > -1 )
						{ actual_name = microenvironment.density_names[ index ]; }
			
						// error check 
						if( std::strcmp( substrate_name.c_str() , actual_name.c_str() ) != 0 )						
						{
							if (IOProcessor(world))
							std::cout << "Warning: when processing advanced chemotaxis for " << pCD->name << " cells: " << std::endl 
									  << "\tInvalid substrate " << substrate_name << " specified." << std::endl
						          	  << "\tIgnoring this invalid substrate in the chemotaxis function .. " << std::endl; 
						}
						else
						{ pCD->phenotype.motility.chemotactic_sensitivities[index] = xml_get_my_double_value(node_cs); }
						node_cs = node_cs.next_sibling( "chemotactic_sensitivity" ); 
					}

				}
				else
				{
					if (IOProcessor(world))
					std::cout << "Warning: when processing motility for " << pCD->name << " cells: " << std::endl 
								<< "\tAdvanced chemotaxis requries chemotactic_sensitivities." << std::endl
								<< "\tBut you have none. Your migration bias will be the zero vector." << std::endl; 
				}

			}



		}

		// display summary for diagnostic help 
		if( pCD->functions.update_migration_bias == chemotaxis_function && pMot->is_motile == true )
		{
			if (IOProcessor(world))
				std::cout << "Cells of type " << pCD->name << " use standard chemotaxis: " << std::endl 
				<< "\t d_bias (before normalization) = " << pMot->chemotaxis_direction << " * grad(" 
				<< microenvironment.density_names[pMot->chemotaxis_index] << ")" << std::endl; 
		}

		if( pCD->functions.update_migration_bias == advanced_chemotaxis_function && pMot->is_motile == true )
		{
			int number_of_substrates = microenvironment.density_names.size(); 
			if (IOProcessor(world))
				std::cout << "Cells of type " << pCD->name << " use advanced chemotaxis: " << std::endl 
				<< "\t d_bias (before normalization) = " 
				<< pMot->chemotactic_sensitivities[0] << " * grad(" << microenvironment.density_names[0] << ")"; 

			for( int n=1; n < number_of_substrates; n++ )
			{ std::cout << " + " << pMot->chemotactic_sensitivities[n] << " * grad(" << microenvironment.density_names[n] << ")"; }
			std::cout << std::endl; 
		}		

		if( pCD->functions.update_migration_bias == advanced_chemotaxis_function_normalized && pMot->is_motile == true )
		{
			int number_of_substrates = microenvironment.density_names.size(); 

			if (IOProcessor(world)) {
				std::cout << "Cells of type " << pCD->name << " use normalized advanced chemotaxis: " << std::endl 
				<< "\t d_bias (before normalization) = " 
				<< pMot->chemotactic_sensitivities[0] << " * grad(" << microenvironment.density_names[0] << ")" 
				<< " / ||grad(" << microenvironment.density_names[0] << ")||"; 

				for( int n=1; n < number_of_substrates; n++ )
				{
					std::cout << " + " << pMot->chemotactic_sensitivities[n] << " * grad(" << microenvironment.density_names[n] << ")"
					<< " / ||grad(" << microenvironment.density_names[n] << ")||"; 
				}
				std::cout << std::endl; 
			}
		}
	}	

	// secretion
	
	node = cd_node.child( "phenotype" );
	node = node.child( "secretion" ); 
	if( node )
	{
		Secretion* pS = &(pCD->phenotype.secretion);
		
		// find the first substrate 
		pugi::xml_node node_sec = node.child( "substrate" );
		while( node_sec )
		{
			// which substrate? 
			
			std::string substrate_name = node_sec.attribute( "name").value(); 
			int index = microenvironment.find_density_index( substrate_name ); 
			std::string actual_name = microenvironment.density_names[ index ]; 
			
			// error check 
			if( std::strcmp( substrate_name.c_str() , actual_name.c_str() ) != 0 )
			{
				std::cout << "[Rank " << world.rank << "] Error: attempted to set secretion/uptake/export for \"" 
					<< substrate_name << "\", which was not found in the microenvironment." << std::endl 
				<< "       Please double-check your substrate name in the config file." << std::endl << std::endl; 
				exit(-1); 
			}			
	
			// secretion rate
			pugi::xml_node node_sec1 = node_sec.child( "secretion_rate" ); 
			if( node_sec1 )
			{ pS->secretion_rates[index] = xml_get_my_double_value( node_sec1 ); }
			
			// secretion target 
			node_sec1 = node_sec.child( "secretion_target" ); 
			if( node_sec1 )
			{ pS->saturation_densities[index] = xml_get_my_double_value( node_sec1 ); }
	
			// uptake rate 
			node_sec1 = node_sec.child( "uptake_rate" ); 
			if( node_sec1 )
			{ pS->uptake_rates[index] = xml_get_my_double_value( node_sec1 ); }
			
			// net export rate 
			node_sec1 = node_sec.child( "net_export_rate" ); 
			if( node_sec1 )
			{ pS->net_export_rates[index] = xml_get_my_double_value( node_sec1 ); }
			
			node_sec = node_sec.next_sibling( "substrate" ); 
		}
	}	
	// cell interactions 

	node = cd_node.child( "phenotype" );
	node = node.child( "cell_interactions" ); 
	if( node )
	{
		Cell_Interactions* pCI = &(pCD->phenotype.cell_interactions);

		// dead_phagocytosis_rate
		pugi::xml_node node_dpr = node.child("dead_phagocytosis_rate");
		double dead_phagocytosis_rate = 0.0; 
		if( node_dpr )
		{
			// pCI->dead_phagocytosis_rate = xml_get_my_double_value(node_dpr); 
			dead_phagocytosis_rate = xml_get_my_double_value(node_dpr); 
		}
/*
		pCI->apoptotic_phagocytosis_rate = pCI->dead_phagocytosis_rate; 
		pCI->necrotic_phagocytosis_rate = pCI->dead_phagocytosis_rate; 
		pCI->other_dead_phagocytosis_rate = pCI->dead_phagocytosis_rate; 
*/
		pCI->apoptotic_phagocytosis_rate = dead_phagocytosis_rate; 
		pCI->necrotic_phagocytosis_rate = dead_phagocytosis_rate; 
		pCI->other_dead_phagocytosis_rate = dead_phagocytosis_rate; 

		// if specific apoptotic rate is specified, overwrite 
		pugi::xml_node node_apr = node.child("apoptotic_phagocytosis_rate"); 
		if( node_apr )
		{ pCI->apoptotic_phagocytosis_rate = xml_get_my_double_value(node_apr);  }

		// if specific necrotic rate is specified, overwrite 
		pugi::xml_node node_npr = node.child("necrotic_phagocytosis_rate"); 
		if( node_npr )
		{ pCI->necrotic_phagocytosis_rate = xml_get_my_double_value(node_npr);  }

		// if specific necrotic rate is specified, overwrite 
		pugi::xml_node node_odpr = node.child("other_dead_phagocytosis_rate"); 
		if( node_odpr )
		{ pCI->other_dead_phagocytosis_rate = xml_get_my_double_value(node_odpr);  }

		// live phagocytosis rates 
		pugi::xml_node node_lpcr = node.child( "live_phagocytosis_rates");
		if( node_lpcr )
		{ node_lpcr = node_lpcr.child( "phagocytosis_rate"); }
		while( node_lpcr )
		{
			// get the name of the target cell type
			std::string target_name = node_lpcr.attribute( "name").value(); 
			// now find its index 
			auto search = cell_definition_indices_by_name.find( target_name );
			// safety first! 
			if( search != cell_definition_indices_by_name.end() )
			{
				// if the target is found, set the appropriate rate 
				int target_index = search->second; 
				std::string target_name_check = search->first; 
				pCI->live_phagocytosis_rates[target_index] = xml_get_my_double_value(node_lpcr); 
			}
			else
			{
				if (IOProcessor(world))
					std::cout << "Warning: When processing the " << pCD->name << " cell definition: " << std::endl 
					<< "\tCould not find cell type " << target_name << " for phagocytosis." << std::endl
					<< "\tIgnoring this live phagocytosis rate!" << std::endl << std::endl; 
			}
			node_lpcr = node_lpcr.next_sibling( "phagocytosis_rate" ); 
		}

		// effector attack rates 
		pugi::xml_node node_ar = node.child( "attack_rates");
		if( node_ar )
		{ node_ar = node_ar.child( "attack_rate"); }
		while( node_ar )
		{
			// get the name of the target cell type
			std::string target_name = node_ar.attribute( "name").value(); 
			// now find its index 
			auto search = cell_definition_indices_by_name.find( target_name );
			// safety first! 
			if( search != cell_definition_indices_by_name.end() )
			{
				// if the target is found, set the appropriate rate 
				int target_index = search->second; 
				std::string target_name_check = search->first; 
				pCI->attack_rates[target_index] = xml_get_my_double_value(node_ar); 
			}
			else
			{
				if (IOProcessor(world))
					std::cout << "Warning: When processing the " << pCD->name << " cell definition: " << std::endl 
					<< "\tCould not find cell type " << target_name << " for cell attack." << std::endl
					<< "\tIgnoring this cell attack rate!" << std::endl << std::endl; 
			}
			node_ar = node_ar.next_sibling( "attack_rate" ); 
		}

		// damage_rate
		pugi::xml_node node_dr = node.child("attack_damage_rate");
		pCI->attack_damage_rate = xml_get_my_double_value(node_dr); 

		// attack_duration
		pugi::xml_node node_ad = node.child("attack_duration");
		pCI->attack_duration = xml_get_my_double_value(node_ad);

		// fusion_rates 
		pugi::xml_node node_fr = node.child( "fusion_rates");
		if( node_fr )
		{ node_fr = node_fr.child( "fusion_rate"); }
		while( node_fr )
		{
			// get the name of the target cell type
			std::string target_name = node_fr.attribute( "name").value(); 
			// now find its index 
			auto search = cell_definition_indices_by_name.find( target_name );
			// safety first! 
			if( search != cell_definition_indices_by_name.end() )
			{
				// if the target is found, set the appropriate rate 
				int target_index = search->second; 
				std::string target_name_check = search->first; 
				pCI->fusion_rates[target_index] = xml_get_my_double_value(node_fr); 
			}
			else
			{
				if (IOProcessor(world))
					std::cout << "Warning: When processing the " << pCD->name << " cell definition: " << std::endl 
					<< "\tCould not find cell type " << target_name << " for cell fusion." << std::endl
					<< "\tIgnoring this cell fusion rate!" << std::endl << std::endl; 
			}
			node_fr = node_fr.next_sibling( "fusion_rate" ); 
		}
	}	

	// cell_transformations>
    //            <transformation_rate

	node = cd_node.child( "phenotype" );
	node = node.child( "cell_transformations" ); 
	if( node )
	{
		Cell_Transformations * pCT = &(pCD->phenotype.cell_transformations);

		// transformation rates 
		pugi::xml_node node_tr = node.child( "transformation_rates");
		if( node_tr )
		{ node_tr = node_tr.child( "transformation_rate"); }
		while( node_tr )
		{
			// get the name of the target cell type
			std::string target_name = node_tr.attribute( "name").value(); 
			// now find its index 
			auto search = cell_definition_indices_by_name.find( target_name );
			// safety first! 
			if( search != cell_definition_indices_by_name.end() )
			{
				// if the target is found, set the appropriate rate 
				int target_index = search->second; 
				std::string target_name_check = search->first; 

				double transformation_rate = xml_get_my_double_value(node_tr);
				if( target_name == pCD->name && transformation_rate > 1e-16 )
				{
					if (IOProcessor(world))
						std::cout << "Warning: When processing the " << pCD->name << " cell definition: " << std::endl 
						<< "\tTransformation from " << pCD->name << " to " << target_name << " is not allowed." << std::endl
						<< "\tIgnoring this cell transformation rate!" << std::endl << std::endl; 
				}
				else
				{ pCT->transformation_rates[target_index] = transformation_rate; }
			}
			else
			{
				if (IOProcessor(world))
					std::cout << "Warning: When processing the " << pCD->name << " cell definition: " << std::endl 
					<< "\tCould not find cell type " << target_name << " for cell transformation." << std::endl
					<< "\tIgnoring this cell transformation rate!" << std::endl << std::endl; 
			}
			node_tr = node_tr.next_sibling( "transformation_rate" ); 
		}

	}

	// cell integrity 
	node = cd_node.child( "phenotype" );
	node = node.child( "cell_integrity" ); 
	if( node )
	{
		Cell_Integrity *pCI = &(pCD->phenotype.cell_integrity);

		pugi::xml_node node_dr = node.child("damage_rate");
		if( node_dr )
		{ pCI->damage_rate = xml_get_my_double_value( node_dr ); }

		pugi::xml_node node_drr = node.child("damage_repair_rate");
		if( node_drr )
		{ pCI->damage_repair_rate = xml_get_my_double_value( node_drr ); }
	}


	// intracellular
	
	node = cd_node.child( "phenotype" );
	node = node.child( "intracellular" ); 
	if( node )
	{
		std::string model_type = node.attribute( "type" ).value(); 
		
#ifdef ADDON_PHYSIBOSS
		if (model_type == "maboss") {
			if( pCD->phenotype.intracellular != NULL )
			{
				delete pCD->phenotype.intracellular;
				pCD->phenotype.intracellular = NULL;
			}

			MaBoSSIntracellular* pIntra = new MaBoSSIntracellular(node);
			pCD->phenotype.intracellular = pIntra->getIntracellularModel();
		}
#endif

#ifdef ADDON_ROADRUNNER
		if (model_type == "roadrunner") 
        {
			// If it has already be copied
			if (pParent != NULL && pParent->phenotype.intracellular != NULL) 
            {
                // std::cout << "------ " << __FUNCTION__ << ": copying another\n";
				pCD->phenotype.intracellular->initialize_intracellular_from_pugixml(node);
            }	
			// Otherwise we need to create a new one
			else 
            {
                std::cout << "\n------ " << __FUNCTION__ << ": creating new RoadRunnerIntracellular\n";
				RoadRunnerIntracellular* pIntra = new RoadRunnerIntracellular(node);
				pCD->phenotype.intracellular = pIntra->getIntracellularModel();
                pCD->phenotype.intracellular->validate_PhysiCell_tokens(pCD->phenotype);
                pCD->phenotype.intracellular->validate_SBML_species();
			}
		}
#endif

#ifdef ADDON_PHYSIDFBA
		if (model_type == "dfba") {
			// If it has already be copied
			if (pParent != NULL && pParent->phenotype.intracellular != NULL) {
				pCD->phenotype.intracellular->initialize_intracellular_from_pugixml(node);
			// Otherwise we need to create a new one
			} else {
				dFBAIntracellular* pIntra = new dFBAIntracellular(node);
				pCD->phenotype.intracellular = pIntra->getIntracellularModel();
			}
		}
#endif
	} else{

		pCD->phenotype.intracellular = NULL;

	}
	
	// set up custom data 
	node = cd_node.child( "custom_data" );
	pugi::xml_node node1 = node.first_child(); 
	while( node1 )
	{
		// name of teh custom data 
		std::string name = xml_get_my_name( node1 ); 
		
		// units 
		std::string units = node1.attribute( "units").value(); 
		
		std::vector<double> values; // ( length, 0.0 ); 

		bool conserved = node1.attribute( "conserved").as_bool(); 

		// get value(s)
		std::string str_values = xml_get_my_string_value( node1 ); 
		csv_to_vector( str_values.c_str() , values ); 
		
		// add variable if cell defaults  
		// if the custom data is not yet found, add it 
		// first, try scalar 
		if( values.size() == 1 )
		{
			// find the variable 
			int n = pCD->custom_data.find_variable_index( name ); 
			// if it exists, overwrite 
			if( n > -1 )
			{ pCD->custom_data.variables[n].value = values[0]; }
			// otherwise, add 
			else
			{ pCD->custom_data.add_variable( name, units, values[0] ); }

			n = pCD->custom_data.find_variable_index( name ); 
			pCD->custom_data.variables[n].conserved_quantity = conserved;
		}
		// or vector 
		else
		{ 
			// find the variable 
			int n = pCD->custom_data.find_vector_variable_index( name ); 
			// if it exists, overwrite 
			if( n > -1 )
			{ pCD->custom_data.vector_variables[n].value = values; }
			// otherwise, add 
			else
			{ pCD->custom_data.add_vector_variable( name, units, values ); }

			n = pCD->custom_data.find_vector_variable_index( name ); 
			pCD->custom_data.vector_variables[n].conserved_quantity = conserved;
		}
		
		node1 = node1.next_sibling(); 
	}
	
	return pCD;
}

void initialize_cell_definitions_from_pugixml( pugi::xml_node root , mpi_Environment world, mpi_Cartesian cart_topo) //1.14 revised
{
	pugi::xml_node node_options; 
	
	node_options = xml_find_node( root , "options" ); 
	if( node_options )
	{
		bool settings = 
			xml_get_bool_value( node_options, "virtual_wall_at_domain_edge" ); 
		if( settings )
		{
			std::cout << "virtual_wall_at_domain_edge: enabled" << std::endl; 
			cell_defaults.functions.add_cell_basement_membrane_interactions = standard_domain_edge_avoidance_interactions;
		}
	}
	
	prebuild_cell_definition_index_maps( world, cart_topo); 

	pugi::xml_node node = root.child( "cell_definitions" ); 
	
	node = node.child( "cell_definition" ); 
	
	while( node )
	{
		std::cout << "Processing " << node.attribute( "name" ).value() << " ... " << std::endl; 
		
		initialize_cell_definition_from_pugixml( node , world, cart_topo);	
		build_cell_definitions_maps(); 
		
		node = node.next_sibling( "cell_definition" ); 
	}
	
//	build_cell_definitions_maps(); 
//	display_cell_definitions( std::cout ); 
	
	return; 
}	

void initialize_cell_definitions_from_pugixml( mpi_Environment world, mpi_Cartesian cart_topo )
{
	initialize_cell_definitions_from_pugixml( physicell_config_root, world, cart_topo );
	return; 
}

int Cell_State::number_of_attached_cells( void )
{ return attached_cells.size(); } 

void Cell::attach_cell( Cell* pAddMe )
{
	#pragma omp critical
	{
		bool already_attached = false; 
		for( int i=0 ; i < state.attached_cells.size() ; i++ )
		{
			if( state.attached_cells[i] == pAddMe )
			{ already_attached = true; }
		}
		if( already_attached == false )
		{ state.attached_cells.push_back( pAddMe ); }
	}
	// pAddMe->attach_cell( this ); 
	return; 
}

void Cell::attach_cell_as_spring( Cell* pAddMe )
{
	#pragma omp critical
	{
		bool already_attached = false; 
		for( int i=0 ; i < state.spring_attachments.size() ; i++ )
		{
			if( state.spring_attachments[i] == pAddMe )
			{ already_attached = true; }
		}
		if( already_attached == false )
		{ state.spring_attachments.push_back( pAddMe ); }
	}
	// pAddMe->attach_cell( this ); 
	return; 
}

void Cell::detach_cell( Cell* pRemoveMe )
{
	#pragma omp critical
	{
		bool found = false; 
		int i = 0; 
		while( !found && i < state.attached_cells.size() )
		{
			// if pRemoveMe is in the cell's list, remove it
			if( state.attached_cells[i] == pRemoveMe )
			{
				int n = state.attached_cells.size(); 
				// copy last entry to current position 
				state.attached_cells[i] = state.attached_cells[n-1]; 
				// shrink by one 
				state.attached_cells.pop_back(); 
				found = true; 
			}
			i++; 
		}
	}
	return; 
}

void Cell::detach_cell_as_spring( Cell* pRemoveMe )
{
	#pragma omp critical
	{
		bool found = false; 
		int i = 0; 
		while( !found && i < state.spring_attachments.size() )
		{
			// if pRemoveMe is in the cell's list, remove it
			if( state.spring_attachments[i] == pRemoveMe )
			{
				int n = state.spring_attachments.size(); 
				// copy last entry to current position 
				state.spring_attachments[i] = state.spring_attachments[n-1]; 
				// shrink by one 
				state.spring_attachments.pop_back(); 
				found = true; 
			}
			i++; 
		}
	}
	return; 
}

void Cell::remove_all_attached_cells( void )
{
	{
		// remove self from any attached cell's list. 
		for( int i = 0; i < state.attached_cells.size() ; i++ )
		{
			state.attached_cells[i]->detach_cell( this ); 
		}
		// clear my list 
		state.attached_cells.clear(); 
	}
	return; 
}

void Cell::remove_all_spring_attachments( void )
{
	{
		// remove self from any attached cell's list. 
		for( int i = 0; i < state.spring_attachments.size() ; i++ )
		{
			state.spring_attachments[i]->detach_cell_as_spring( this ); 
		}
		// clear my list 
		state.spring_attachments.clear(); 
	}
	return; 
}

void attach_cells( Cell* pCell_1, Cell* pCell_2 )
{
	pCell_1->attach_cell( pCell_2 );
	pCell_2->attach_cell( pCell_1 );
	return; 
}

void attach_cells_as_spring( Cell* pCell_1, Cell* pCell_2 )
{
	pCell_1->attach_cell_as_spring( pCell_2 );
	pCell_2->attach_cell_as_spring( pCell_1 );
	return; 
}

void detach_cells( Cell* pCell_1 , Cell* pCell_2 )
{
	pCell_1->detach_cell( pCell_2 );
	pCell_2->detach_cell( pCell_1 );
	return; 
}

void detach_cells_as_spring( Cell* pCell_1 , Cell* pCell_2 )
{
	pCell_1->detach_cell_as_spring( pCell_2 );
	pCell_2->detach_cell_as_spring( pCell_1 );
	return; 
}

std::vector<Cell*> find_nearby_cells( Cell* pCell )
{
	std::vector<Cell*> neighbors = {}; 

	// First check the neighbors in my current voxel
	std::vector<Cell*>::iterator neighbor;
	std::vector<Cell*>::iterator end =
		pCell->get_container()->agent_grid[pCell->get_current_mechanics_voxel_index()].end();
	for( neighbor = pCell->get_container()->agent_grid[pCell->get_current_mechanics_voxel_index()].begin(); neighbor != end; ++neighbor)
	{ neighbors.push_back( *neighbor ); }

	std::vector<int>::iterator neighbor_voxel_index;
	std::vector<int>::iterator neighbor_voxel_index_end = 
		pCell->get_container()->underlying_mesh.moore_connected_voxel_indices[pCell->get_current_mechanics_voxel_index()].end();
	
	for( neighbor_voxel_index = 
		pCell->get_container()->underlying_mesh.moore_connected_voxel_indices[pCell->get_current_mechanics_voxel_index()].begin();
		neighbor_voxel_index != neighbor_voxel_index_end; 
		++neighbor_voxel_index )
	{
		if(!is_neighbor_voxel(pCell, pCell->get_container()->underlying_mesh.voxels[pCell->get_current_mechanics_voxel_index()].center, pCell->get_container()->underlying_mesh.voxels[*neighbor_voxel_index].center, *neighbor_voxel_index))
			continue;
		end = pCell->get_container()->agent_grid[*neighbor_voxel_index].end();
		for(neighbor = pCell->get_container()->agent_grid[*neighbor_voxel_index].begin();neighbor != end; ++neighbor)
		{ neighbors.push_back( *neighbor ); }
	}
	
	return neighbors; 
}

std::vector<Cell*> find_nearby_interacting_cells( Cell* pCell )
{
	std::vector<Cell*> neighbors = {}; 

	// First check the neighbors in my current voxel
	std::vector<Cell*>::iterator neighbor;
	std::vector<Cell*>::iterator end = pCell->get_container()->agent_grid[pCell->get_current_mechanics_voxel_index()].end();
	for( neighbor = pCell->get_container()->agent_grid[pCell->get_current_mechanics_voxel_index()].begin(); neighbor != end; ++neighbor)
	{
		std::vector<double> displacement = (*neighbor)->position - pCell->position; 
		double distance = norm( displacement ); 
		if( distance <= pCell->phenotype.mechanics.relative_maximum_adhesion_distance * pCell->phenotype.geometry.radius 
			+ (*neighbor)->phenotype.mechanics.relative_maximum_adhesion_distance * (*neighbor)->phenotype.geometry.radius 
			&& (*neighbor) != pCell )
		{ neighbors.push_back( *neighbor ); }
	}

	std::vector<int>::iterator neighbor_voxel_index;
	std::vector<int>::iterator neighbor_voxel_index_end = 
		pCell->get_container()->underlying_mesh.moore_connected_voxel_indices[pCell->get_current_mechanics_voxel_index()].end();
	
	for( neighbor_voxel_index = 
		pCell->get_container()->underlying_mesh.moore_connected_voxel_indices[pCell->get_current_mechanics_voxel_index()].begin();
		neighbor_voxel_index != neighbor_voxel_index_end; 
		++neighbor_voxel_index )
	{
		if(!is_neighbor_voxel(pCell, pCell->get_container()->underlying_mesh.voxels[pCell->get_current_mechanics_voxel_index()].center, pCell->get_container()->underlying_mesh.voxels[*neighbor_voxel_index].center, *neighbor_voxel_index))
			continue;
		end = pCell->get_container()->agent_grid[*neighbor_voxel_index].end();
		for(neighbor = pCell->get_container()->agent_grid[*neighbor_voxel_index].begin();neighbor != end; ++neighbor)
		{
			std::vector<double> displacement = (*neighbor)->position - pCell->position; 
			double distance = norm( displacement ); 
			if( distance <= pCell->phenotype.mechanics.relative_maximum_adhesion_distance * pCell->phenotype.geometry.radius 
				+ (*neighbor)->phenotype.mechanics.relative_maximum_adhesion_distance * (*neighbor)->phenotype.geometry.radius
				&& (*neighbor) != pCell	)
			{ neighbors.push_back( *neighbor ); }
		}
	}
	
	return neighbors; 
}


int find_cell_definition_index( std::string search_string )
{
	auto search = cell_definition_indices_by_name.find( search_string );
	// safety first! 
	if( search != cell_definition_indices_by_name.end() )
	{
		// if the target is found, set the appropriate rate 
		return search->second; 
	}
	return -1; 
}


int find_cell_definition_index( int search_type )
{	
	auto search = cell_definition_indices_by_type.find( search_type );
	// safety first! 
	if( search != cell_definition_indices_by_type.end() )
	{
		// if the target is found, set the appropriate rate 
		return search->second; 
	}
	return -1; 
}  

};
