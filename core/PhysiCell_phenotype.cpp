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

#include "./PhysiCell_phenotype.h"

#include "../BioFVM/BioFVM.h"
#include "./PhysiCell_constants.h"
#include "./PhysiCell_utilities.h"
#include "./MPI_helper.h"

#ifdef ADDON_PHYSIBOSS
#include "../addons/PhysiBoSS/src/maboss_intracellular.h"
#endif


using namespace BioFVM; 

namespace PhysiCell{
	
Phase::Phase() //UPDATED
{
	index = 0; 
	code = 0; 
	name = "unnamed"; 
	
	division_at_phase_exit = false; 
	removal_at_phase_exit = false; 
	
	entry_function = NULL; 
	return; 
}

void Phase::pack(std::vector<char>& snd_buffer, int& len_buffer, int& position){
	
	pack_buff( this->index, snd_buffer, len_buffer, position);
	pack_buff( this->code, snd_buffer, len_buffer, position);
	pack_buff( this->name, snd_buffer, len_buffer, position);
	pack_buff( this->division_at_phase_exit, snd_buffer, len_buffer, position);
	pack_buff( this->removal_at_phase_exit, snd_buffer, len_buffer, position);
}

void Phase::unpack(std::vector<char>& rcv_buffer, int& len_buffer, int& position){

	unpack_buff(this->index, rcv_buffer, len_buffer, position);
    unpack_buff(this->code, rcv_buffer, len_buffer, position);
    unpack_buff(this->name, rcv_buffer, len_buffer, position);
    unpack_buff(this->division_at_phase_exit, rcv_buffer,len_buffer,  position);
    unpack_buff(this->removal_at_phase_exit, rcv_buffer, len_buffer, position);
}
	
Phase_Link::Phase_Link() //UPDATED
{
	start_phase_index = 0; 
	end_phase_index = 0; 
	
	fixed_duration = false; 
	
	arrest_function = NULL; 
	exit_function = NULL; 
	return; 
}

void Phase_Link::pack(std::vector<char>& snd_buffer, int& len_buffer, int& position){
	pack_buff(this->start_phase_index, snd_buffer, len_buffer, position);
	pack_buff(this->end_phase_index, snd_buffer, len_buffer, position);
	pack_buff(this->fixed_duration, snd_buffer, len_buffer, position);
}

void Phase_Link::unpack(std::vector<char>& rcv_buffer, int& len_buffer, int& position){
	unpack_buff(this->start_phase_index, rcv_buffer, len_buffer, position);
	unpack_buff(this->end_phase_index, rcv_buffer, len_buffer, position);
	unpack_buff(this->fixed_duration, rcv_buffer, len_buffer, position);
}

Cycle_Data::Cycle_Data() //UPDATED
{
	inverse_index_maps.resize(0); 
	
	pCycle_Model = NULL; 

	time_units = "min"; 
	
	transition_rates.resize( 0 );

	current_phase_index = 0; 
	elapsed_time_in_phase = 0.0; 
	return; 
}

void Cycle_Data::sync_to_cycle_model( void ) //UPDATED
{
	// make sure the inverse map is the right size 
	int n = pCycle_Model->phases.size(); 
	inverse_index_maps.resize( n );
	
	// sync the inverse map to the cell cycle model by 
	// querying the phase_links 

	transition_rates.resize( n );
	
	// also make sure the transition_rates[] are the right size 
	
	for( int i=0 ; i < pCycle_Model->phase_links.size() ; i++ )
	{
		inverse_index_maps[i].clear(); 
		for( int j=0 ; j < pCycle_Model->phase_links[i].size() ; j++ )
		{
			inverse_index_maps[i][ pCycle_Model->phase_links[i][j].end_phase_index ] = j;
			transition_rates[i].resize( pCycle_Model->phase_links[i].size() ); 
		}
	}

	return; 
}

double& Cycle_Data::transition_rate( int start_phase_index , int end_phase_index ) //UPDATED
{
	return transition_rates[ start_phase_index ][ inverse_index_maps[start_phase_index][end_phase_index] ]; 
}

double& Cycle_Data::exit_rate(int phase_index ) //UPDATED
{
	return transition_rates[phase_index][0]; 
}

std::vector<std::unordered_map<int,int>> & Cycle_Data:: get_inverse_index_maps()
{
	return inverse_index_maps; 
}

void Cycle_Data::pack(std::vector<char>& snd_buffer, int& len_buffer, int& position) {

	pack_buff( static_cast<int>(this->inverse_index_maps.size()), snd_buffer, len_buffer, position);
	
    for(int i=0; i<(this->inverse_index_maps.size()); i++)
    {
		pack_buff( static_cast<int>(this->inverse_index_maps[i].size()), snd_buffer, len_buffer, position);

        for(auto it = this->inverse_index_maps[i].cbegin(); it != this->inverse_index_maps[i].cend(); it++)
        {
			pack_buff( it->first, snd_buffer, len_buffer, position);
			pack_buff( it->second, snd_buffer, len_buffer, position);
        }
    }

    /*=====================================================================*/
    /* Cycle_Model* pCycle_Model;  IS NOT PACKED 						   */
    /*=====================================================================*/

	pack_buff( this->time_units, snd_buffer, len_buffer, position);

	pack_buff( static_cast<int>(this->transition_rates.size()), snd_buffer, len_buffer, position);
    for(int i=0; i<this->transition_rates.size(); i++)
    {
		pack_buff(  this->transition_rates[i], snd_buffer, len_buffer, position);

    /* this->functions.cycle_data.data.transition_rates[i] - is INCORRECT, as MPI needs base address of a double array and NOT a vector address */
    //MPI_Pack(&(this->functions.cycle_model.data.transition_rates[i][0]), len_int, MPI_DOUBLE, &snd_buffer[0], len_buffer, &position, MPI_COMM_WORLD);
    }

    /* 1 int and 1 double are next to be packed in class Cycle_Data */
	pack_buff( this->current_phase_index, snd_buffer, len_buffer, position);
	pack_buff( this->elapsed_time_in_phase, snd_buffer, len_buffer, position);
}

void Cycle_Data::unpack(std::vector<char>& rcv_buffer, int& len_buffer, int& position) {

	// Cycle_Data
	auto& cycle_data_maps = this->get_inverse_index_maps();
	int data_map_size;
	unpack_buff(data_map_size, rcv_buffer, len_buffer, position);
	cycle_data_maps.resize(data_map_size);
	for (int i = 0; i < data_map_size; ++i) {
		int msize;
		cycle_data_maps[i].clear();
		unpack_buff(msize, rcv_buffer, len_buffer, position);
		for (int j = 0; j < msize; ++j) {
			int key, value;
			unpack_buff(key, rcv_buffer, len_buffer, position);
			unpack_buff(value, rcv_buffer, len_buffer, position);
			cycle_data_maps[i][key] = value;
		}
	}

	unpack_buff(this->time_units, rcv_buffer, len_buffer, position);

	int tr_size;
    unpack_buff(tr_size, rcv_buffer, len_buffer, position);
    this->transition_rates.resize(tr_size);
    for (int i = 0; i < tr_size; ++i) {
        unpack_buff(this->transition_rates[i], rcv_buffer, len_buffer, position);
    }

	unpack_buff(this->current_phase_index, rcv_buffer, len_buffer,position);
    unpack_buff(this->elapsed_time_in_phase, rcv_buffer,  len_buffer,position);
}
	
Cycle_Model::Cycle_Model() //UPDATED
{
	inverse_index_maps.resize( 0 );
	
	name = "unnamed";
	
	phases.resize(0); 
	phase_links.resize(0); 
	
	data.pCycle_Model = this; 
	
	code = PhysiCell_constants::custom_cycle_model; 
	default_phase_index = 0; 
	
	return; 
}	
	
int Cycle_Model::add_phase( int code, std::string name ) //UPDATED
{
	int n = phases.size();
	
	// resize the data structures 
	phases.resize( n+1 ); 
	phase_links.resize( n+1 );
	phase_links[n].resize(0);
	
	inverse_index_maps.resize( n+1 );
	inverse_index_maps[n].clear(); 
	
	// update phase n
	phases[n].code = code; 
	phases[n].index = n; 
	
	phases[n].name.assign( name ); 
	
	// make sure the cycle_data is also correctly sized 
	
	data.sync_to_cycle_model();
	
	return n; 
}
	
int Cycle_Model::add_phase_link( int start_index, int end_index , 
	bool (*arrest_function)(Cell* pCell, Phenotype& phenotype, double dt) ) //UPDATED
{
	// first, resize the phase links 
	int n = phase_links[start_index].size(); 
	phase_links[start_index].resize( n + 1 );
	
	// now, update the new phase links
	phase_links[start_index][n].start_phase_index  = start_index; 
	phase_links[start_index][n].end_phase_index = end_index; 
	phase_links[start_index][n].arrest_function = arrest_function; 
	
	// now, update the inverse index map 
	inverse_index_maps[start_index][end_index] = n; 
	
	// lastly, make sure the transition rates are the right size;
	
	data.sync_to_cycle_model(); 

	return n; 
}
	
int Cycle_Model::add_phase_link( int start_index, int end_index , double rate , 
	bool (*arrest_function)(Cell* pCell, Phenotype& phenotype, double dt) ) //UPDATED
{
	int n = add_phase_link( start_index , end_index , arrest_function );
	data.transition_rate( start_index , end_index ) = rate; 
	return n;
}

int Cycle_Model::find_phase_index( int code ) //UPDATED
{
	for( int i=0 ; i < phases.size() ; i++ )
	{
		if( phases[i].code == code )
		{ return i; }
	}
	return 0; 
}

int Cycle_Model::find_phase_index( std::string name ) //UPDATED
{
	for( int i=0 ; i < phases.size() ; i++ )
	{
		if( phases[i].name == name )
		{ return i; }
	}
	return 0; 
}
	
std::ostream& Cycle_Model::display( std::ostream& os ) //UPDATED
{
	os << "Cycle Model: " << name << " (PhysiCell code: " << code << ")" << std::endl; 
	os << "Phases and links: (* denotes phase with cell division)" << std::endl;
	for( int i=0; i < phases.size() ; i++ )
	{
		os << "Phase " << i << " (" << phases[i].name << ") "; 
		
		if( phases[i].division_at_phase_exit )
		{ os << "*"; }
		os << " links to: " << std::endl;
		for( int k=0 ; k < phase_links[i].size() ; k++ )
		{
			int j = phase_links[i][k].end_phase_index; 
			os << "\tPhase " << j << " (" << phases[j].name << ") with rate " << data.transition_rate(i,j) << " " << data.time_units << "^-1; " << std::endl; 
		}
		os << std::endl; 
	}
	
	return os; 
}
	
double& Cycle_Model::transition_rate( int start_index , int end_index ) //UPDATED
{
	return data.transition_rate( start_index , end_index ); 
}

Phase_Link& Cycle_Model::phase_link( int start_index, int end_index ) //UPDATED
{
	return phase_links[start_index][ inverse_index_maps[start_index][end_index] ]; 
}
	
void Cycle_Model::advance_model( Cell* pCell, Phenotype& phenotype, double dt ) //UPDATED
{
	int i = phenotype.cycle.data.current_phase_index; 
	
	phenotype.cycle.data.elapsed_time_in_phase += dt; 

	// Evaluate each linked phase: 
	// advance to that phase IF probabiltiy is in the range, 
	// and if the arrest function (if any) is false 
	
	int j; 
	for( int k=0 ; k < phase_links[i].size() ; k++ )
	{
		j = phase_links[i][k].end_phase_index; 
		
		// check for arrest. If arrested, skip to the next transition
		bool transition_arrested = false; 
		if( phase_links[i][k].arrest_function )
		{
			transition_arrested = phase_links[i][k].arrest_function( pCell,phenotype,dt ); 
		}
		if( !transition_arrested )
		{
			// check to see if we should transition 
			bool continue_transition = false; 
			if( phase_links[i][k].fixed_duration )
			{
				if( phenotype.cycle.data.elapsed_time_in_phase > ((1.0/phenotype.cycle.data.transition_rates[i][k]) - 0.5 * dt) )
				{
					continue_transition = true; 
				}
			}
			else
			{
				double prob = phenotype.cycle.data.transition_rates[i][k]*dt; 
				if( UniformRandom() < prob )
				{
					continue_transition = true; 
				}
			}
			
			// if we should transition, check if we're not supposed to divide or die 
			
			if( continue_transition )
			{
				// if the phase transition has an exit function, execute it
				if( phase_links[i][k].exit_function )
				{
					phase_links[i][k].exit_function( pCell,phenotype,dt ); 
				}
				
				// check if division or removal are required 
				if( phases[i].division_at_phase_exit )
				{
					// pCell->flag_for_division();
					phenotype.flagged_for_division = true; 
				}
				if( phases[i].removal_at_phase_exit )
				{
					// pCell->flag_for_removal(); 
					phenotype.flagged_for_removal = true; 
					return; 
				}
				// move to the next phase, and reset the elapsed time 
				phenotype.cycle.data.current_phase_index = j; 
				phenotype.cycle.data.elapsed_time_in_phase = 0.0; 
				
				// if the new phase has an entry function, execute it 
				if( phases[j].entry_function )
				{
					phases[j].entry_function( pCell,phenotype,dt );  
				}
				
				return; 
			}
			
		}
		
	}
	
	return; 
}

/*-----------------------------------------------*/
/* Adding a function get_inverse_index_maps() to */ 
/* return a reference to the private data 			 */
/* Added by Gaurav Saxena												 */
/*-----------------------------------------------*/

std::vector<std::unordered_map<int,int>> & Cycle_Model:: get_inverse_index_maps()
{
	return inverse_index_maps; 
}

void Cycle_Model::pack(std::vector<char>& snd_buffer, int& len_buffer, int &position){
	
	std::vector<std::unordered_map<int,int>> &inverse_index_maps_cycle_model = this->get_inverse_index_maps();
	
	pack_buff( static_cast<int>(inverse_index_maps_cycle_model.size()), snd_buffer, len_buffer, position);

    for(int i=0; i<inverse_index_maps_cycle_model.size(); i++)
    {
		pack_buff( static_cast<int>(inverse_index_maps_cycle_model[i].size()), snd_buffer, len_buffer, position);

        for(auto it = inverse_index_maps_cycle_model[i].cbegin(); it != inverse_index_maps_cycle_model[i].cend(); it++)
        {
            pack_buff( it->first, snd_buffer, len_buffer, position);
            pack_buff( it->second, snd_buffer, len_buffer, position);
        }
    }

	pack_buff( this->name, snd_buffer, len_buffer, position);

	pack_buff( this->code, snd_buffer, len_buffer, position);

	/* Packing of data members of Phases as its object is a member of Cycle_Model */
	/* std::vector<Phase> phases; is a member of it - its a vector */

	
	pack_buff( static_cast<int>(this->phases.size()), snd_buffer, len_buffer, position);

    for(int i=0; i<this->phases.size(); i++)
    {
		this->phases[i].pack(snd_buffer, len_buffer, position);
    }

	/* Packing members of Phase_Link - std::vector<std::vector<Phase_Link>> phase_links is a member of */
	/* class Cycle_Model which is a member of Functions - we're packing in DFS i.e. Inorder traversal	 */

	
	pack_buff( static_cast<int>(this->phase_links.size()), snd_buffer, len_buffer, position);

    for(int i=0; i<this->phase_links.size(); i++)
    {
		pack_buff(static_cast<int>(this->phase_links[i].size()), snd_buffer, len_buffer, position);
        for(int j=0; j<this->phase_links[i].size(); j++)
        {
			this->phase_links[i][j].pack(snd_buffer, len_buffer, position);
        }
    }

	pack_buff(this->default_phase_index, snd_buffer, len_buffer, position);

	//Pack Cycle_Data data
	this->data.pack(snd_buffer, len_buffer, position);
    
}

void Cycle_Model::unpack(std::vector<char>& rcv_buffer, int& len_buffer, int &position){
	
	 // Functions.cycle_model.inverse_index_maps
    int inv_map_size;
    unpack_buff(inv_map_size, rcv_buffer, len_buffer, position);
    auto& inverse_index_maps = this->get_inverse_index_maps();
    inverse_index_maps.resize(inv_map_size);
    for (int i = 0; i < inv_map_size; ++i) {
        int map_len;
        unpack_buff(map_len, rcv_buffer, len_buffer, position);
        for (int j = 0; j < map_len; ++j) {
            int key, value;
            unpack_buff(key, rcv_buffer, len_buffer, position);
            unpack_buff(value, rcv_buffer,len_buffer, position);
            inverse_index_maps[i][key] = value;
        }
    }

	unpack_buff(this->name, rcv_buffer, len_buffer, position);
    unpack_buff(this->code, rcv_buffer, len_buffer, position);

	//Phases
	int phase_count;
    unpack_buff(phase_count, rcv_buffer, len_buffer, position);
    this->phases.resize(phase_count);
    for (int i = 0; i < phase_count; ++i) {
		this->phases[i].unpack(rcv_buffer, len_buffer, position);
    }

	// Phase Links
	int link_size;
	unpack_buff(link_size, rcv_buffer, len_buffer, position);
	this->phase_links.resize(link_size);
	for (int i = 0; i < link_size; ++i) {
		int sub_size;
		unpack_buff(sub_size, rcv_buffer, len_buffer, position);
		this->phase_links[i].resize(sub_size);
		for (int j = 0; j < sub_size; ++j) {
			this->phase_links[i][j].unpack(rcv_buffer, len_buffer, position);
		}
	}

	unpack_buff(this->default_phase_index, rcv_buffer, len_buffer, position);

	this->data.unpack(rcv_buffer, len_buffer, position);
}
	
Phase& Cycle_Data::current_phase( void ) //UPDATED
{
	return pCycle_Model->phases[current_phase_index]; 
}

Death_Parameters::Death_Parameters() //UPDATED
{
	time_units = "min"; 
	
	// reference values: MCF-7 (1/min)
	unlysed_fluid_change_rate = 3.0/60.0; // apoptosis 
	lysed_fluid_change_rate = 0.05/60.0; // lysed necrotic cell 
	
	cytoplasmic_biomass_change_rate = 1.0/60.0; // apoptosis 
	nuclear_biomass_change_rate = 0.35/60.0; // apoptosis 
	
	calcification_rate = 0.0; // 0.0042 for necrotic cells 
	
	relative_rupture_volume = 2.0; 

	return; 
}

void Death_Parameters::pack(std::vector<char>& snd_buffer, int& len_buffer, int &position){
	pack_buff(this->time_units, snd_buffer, len_buffer, position);
    pack_buff(this->unlysed_fluid_change_rate,      snd_buffer, len_buffer, position);
	pack_buff(this->lysed_fluid_change_rate,        snd_buffer, len_buffer, position);
	pack_buff(this->cytoplasmic_biomass_change_rate,snd_buffer, len_buffer, position);
	pack_buff(this->nuclear_biomass_change_rate,    snd_buffer, len_buffer, position);
	pack_buff(this->calcification_rate,             snd_buffer, len_buffer, position);
	pack_buff(this->relative_rupture_volume,        snd_buffer, len_buffer, position);
}

void Death_Parameters::unpack(std::vector<char>& rcv_buffer, int& len_buffer, int &position){
	unpack_buff(this->time_units, rcv_buffer, len_buffer, position);
	unpack_buff(this->unlysed_fluid_change_rate,       rcv_buffer, len_buffer, position);
	unpack_buff(this->lysed_fluid_change_rate,         rcv_buffer, len_buffer, position);
	unpack_buff(this->cytoplasmic_biomass_change_rate, rcv_buffer, len_buffer, position);
	unpack_buff(this->nuclear_biomass_change_rate,     rcv_buffer, len_buffer, position);
	unpack_buff(this->calcification_rate,              rcv_buffer, len_buffer, position);
	unpack_buff(this->relative_rupture_volume,         rcv_buffer, len_buffer, position);
}
	
Death::Death() //UPDATED
{
	rates.resize( 0 ); 
	models.resize( 0 ); 
	parameters.resize( 0 ); 
	
	dead = false; 
	current_death_model_index = 0;
	
	return; 
}

int Death::add_death_model( double rate , Cycle_Model* pModel )
{
	rates.push_back( rate );
	models.push_back( pModel ); 
	
	parameters.resize( rates.size() ); 
	
	return rates.size() - 1; 
}

int Death::add_death_model( double rate, Cycle_Model* pModel, Death_Parameters& death_parameters)
{
	rates.push_back( rate );
	models.push_back( pModel ); 
	parameters.push_back( death_parameters ); 
	
	return rates.size() - 1; 
}

int Death::find_death_model_index( int code )
{
	for( int i=0 ; i < models.size() ; i++ )
	{
		if( models[i]->code == code )
		{ return i; }
	}
	return 0; 
}

int Death::find_death_model_index( std::string name )
{
	for( int i=0 ; i < models.size() ; i++ )
	{
		if( models[i]->name == name )
		{ return i; }
	}
	return 0; 
}
	
bool Death::check_for_death( double dt )
{
	// If the cell is already dead, exit. 
	if( dead == true )
	{
		return false;
	} 
	
	// If the cell is alive, evaluate all the 
	// death rates for each registered death type. 
	int i = 0; 
	while( !dead && i < rates.size() )
	{
		if( UniformRandom() < rates[i]*dt )
		{
			// update the Death data structure 
			dead = true; 
			current_death_model_index = i; 
			
			// and set the cycle model to this death model 
			
			return dead;
		}
		i++; 
	}
	
	return dead; 
}

void Death::trigger_death( int death_model_index )
{
	dead = true; 
	current_death_model_index = death_model_index; 
	
/*	
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
*/
		
	
	return; 
}

Cycle_Model& Death::current_model( void )
{
	return *models[current_death_model_index]; 
}

double& Death::apoptosis_rate(void)
{
	static int nApoptosis = find_death_model_index( PhysiCell_constants::apoptosis_death_model ); 
	return rates[nApoptosis];
}

double& Death::necrosis_rate(void)
{
	static int nNecrosis = find_death_model_index( PhysiCell_constants::necrosis_death_model ); 
	return rates[nNecrosis];
}

void Death::pack(std::vector<char>& snd_buffer, int& len_buffer, int &position){
	 /* Now packing data members of class Death as there is an object of class Death in Phenotype class */

	pack_buff( this->rates, snd_buffer, len_buffer, position);

    /* Next field is std::vector<Cycle_Model*> models ;  this is a vector of Cycle_Model * pointers */

    /*=====================================================================*/
    /* std::vector<Cycle_Model*> models ;  IS NOT PACKED  								 */
    /*=====================================================================*/

    /* Now packing members of class Death_Parameters as class Death contains */
    /* std::vector<Death_Parameters> parameters 														 */
	
	pack_buff( static_cast<int>(this->parameters.size()), snd_buffer, len_buffer, position);

    for(int i=0; i<this->parameters.size(); i++)
    {
		this->parameters[i].pack(snd_buffer, len_buffer, position);
	}

    /* Next 1 bool value (allocate sizeof(int)) and 1 int value are to be packed */

	pack_buff(this->dead ,snd_buffer, len_buffer, position);

	pack_buff(this->current_death_model_index ,snd_buffer, len_buffer, position);
}

void Death::unpack(std::vector<char>& rcv_buffer, int& len_buffer, int &position){
	// Unpack death.rates (std::vector<double>)
	unpack_buff(this->rates, rcv_buffer, len_buffer, position);

	// Unpack size of phenotype.death.parameters
	int death_parameters_size;
	unpack_buff(death_parameters_size, rcv_buffer, len_buffer, position);
	this->parameters.resize(death_parameters_size);

	// Unpack each Death_Parameters object
	for (size_t i = 0; i < death_parameters_size; i++)
	{
		this->parameters[i].unpack(rcv_buffer, len_buffer, position);
	}

	// Unpack bool phenotype.death.dead
	unpack_buff(this->dead, rcv_buffer, len_buffer, position);

	// Unpack int phenotype.death.current_death_model_index
	unpack_buff(this->current_death_model_index, rcv_buffer, len_buffer, position);

}


Cycle::Cycle()
{
	pCycle_Model = NULL; 
	return; 
}

void Cycle::advance_cycle( Cell* pCell, Phenotype& phenotype, double dt )
{
	pCycle_Model->advance_model( pCell, phenotype , dt ); 
	return; 
}

Cycle_Model& Cycle::model( void )
{
	return *pCycle_Model; 
}

Phase& Cycle::current_phase( void )
{
	return data.current_phase(); 
}

int& Cycle::current_phase_index( void )
{
	return data.current_phase_index; 
}

void Cycle::sync_to_cycle_model( Cycle_Model& cm )
{
	pCycle_Model = &cm; 
	data = cm.data; 
	return; 
}

void Cycle::pack(std::vector<char>& snd_buffer, int& len_buffer, int &position) {
	 /* Packing of members of class Cycle as the class Phenotype contains a data member of this class */
    /* Not packing this->phenotype.cycle.pCycle_Model pointer data member *

    /*=====================================================================*/
    /* Cycle_Model* pCycle_Model;  IS NOT PACKED 						   */
    /*=====================================================================*/

    /* Packing data fields of Cycle_Data as it is an object of Cycle class which is    */
    /* a data member of the Phenotype class which is a member of class Cell						*/
    /* First field is std::vector<std::unordered_map<int,int>> inverse_index_maps;			*/
	this->data.pack(snd_buffer, len_buffer, position);

	this->asymmetric_division.pack(snd_buffer, len_buffer, position);
	
}

void Cycle::unpack(std::vector<char>& rcv_buffer, int& len_buffer, int &position) {

	this->data.unpack(rcv_buffer, len_buffer, position);

	this->asymmetric_division.unpack(rcv_buffer, len_buffer, position);
	
}





Death_Parameters& Death::current_parameters( void )
{
	return parameters[ current_death_model_index ]; 
}
	
Volume::Volume()
{
	// reference parameter values for MCF-7, in cubic microns 
	fluid_fraction = 0.75;  

	total = 2494; 
	fluid = fluid_fraction * total; 
	solid = total-fluid; 
	
	nuclear = 540.0; 
	nuclear_fluid = fluid_fraction * nuclear; 
	nuclear_solid = nuclear - nuclear_fluid;

	cytoplasmic = total - nuclear;
	cytoplasmic_fluid = fluid_fraction*cytoplasmic; 
	cytoplasmic_solid = cytoplasmic - cytoplasmic_fluid; 
	
	// rates are in units of 1/min 
	cytoplasmic_biomass_change_rate = 0.27 / 60.0; 
	nuclear_biomass_change_rate = 0.33 / 60.0; 
	fluid_change_rate = 3.0 / 60.0;

	calcified_fraction = 0.0;
	
	calcification_rate = 0.0; 
	
	target_solid_cytoplasmic = cytoplasmic_solid;
	target_solid_nuclear = nuclear_solid;
	target_fluid_fraction = fluid_fraction;
	
	cytoplasmic_to_nuclear_ratio = cytoplasmic / ( 1e-16 + nuclear);
	target_cytoplasmic_to_nuclear_ratio = cytoplasmic_to_nuclear_ratio; 
	
	// the cell bursts at these volumes 
	relative_rupture_volume = 2.0; 
		// as fraction of volume at entry to the current phase
	rupture_volume = relative_rupture_volume * total; // in volume units 
	
	return; 
};


void Volume::multiply_by_ratio( double ratio )
{
	total *= ratio;
	solid *= ratio;
	fluid *= ratio;
	
	nuclear *= ratio;
	nuclear_fluid *= ratio;
	nuclear_solid *= ratio;
	
	cytoplasmic *= ratio;
	cytoplasmic_fluid *= ratio;
	cytoplasmic_solid *= ratio;
	
	rupture_volume *= ratio; 
	
	target_solid_nuclear *= ratio;
	target_solid_cytoplasmic *= ratio; 	
	
	return; 
}

void Volume::divide( void )
{
	multiply_by_ratio( 0.5 ); 
	return; 
}

void Volume::pack(std::vector<char>& snd_buffer, int& len_buffer, int &position){

	pack_buff(this->total ,snd_buffer, len_buffer, position);
	pack_buff(this->solid ,snd_buffer, len_buffer, position);
	pack_buff(this->fluid ,snd_buffer, len_buffer, position);

	pack_buff(this->fluid_fraction ,snd_buffer, len_buffer, position);
	pack_buff(this->nuclear,snd_buffer, len_buffer, position);
	pack_buff(this->nuclear_fluid ,snd_buffer, len_buffer, position);

    pack_buff(this->nuclear_solid,                            snd_buffer, len_buffer, position);
	pack_buff(this->cytoplasmic,                              snd_buffer, len_buffer, position);
	pack_buff(this->cytoplasmic_fluid,                        snd_buffer, len_buffer, position);

	pack_buff(this->cytoplasmic_solid,                        snd_buffer, len_buffer, position);
	pack_buff(this->calcified_fraction,                       snd_buffer, len_buffer, position);
	pack_buff(this->cytoplasmic_to_nuclear_ratio,            snd_buffer, len_buffer, position);

	pack_buff(this->rupture_volume,                           snd_buffer, len_buffer, position);
	pack_buff(this->cytoplasmic_biomass_change_rate,         snd_buffer, len_buffer, position);
	pack_buff(this->nuclear_biomass_change_rate,             snd_buffer, len_buffer, position);

	pack_buff(this->fluid_change_rate,                        snd_buffer, len_buffer, position);
	pack_buff(this->calcification_rate,                      snd_buffer, len_buffer, position);
	pack_buff(this->target_solid_cytoplasmic,                snd_buffer, len_buffer, position);

	pack_buff(this->target_solid_nuclear,                    snd_buffer, len_buffer, position);
	pack_buff(this->target_fluid_fraction,                   snd_buffer, len_buffer, position);
	pack_buff(this->target_cytoplasmic_to_nuclear_ratio,     snd_buffer, len_buffer, position);

	pack_buff(this->relative_rupture_volume,                 snd_buffer, len_buffer, position);
}

void Volume::unpack(std::vector<char>& rcv_buffer, int& len_buffer, int &position){
	// Unpack phenotype.volume fields (all doubles)
	unpack_buff(this->total,                              rcv_buffer, len_buffer, position);
	unpack_buff(this->solid,                              rcv_buffer, len_buffer, position);
	unpack_buff(this->fluid,                              rcv_buffer, len_buffer, position);

	unpack_buff(this->fluid_fraction,                     rcv_buffer, len_buffer, position);
	unpack_buff(this->nuclear,                            rcv_buffer, len_buffer, position);
	unpack_buff(this->nuclear_fluid,                      rcv_buffer, len_buffer, position);

	unpack_buff(this->nuclear_solid,                      rcv_buffer, len_buffer, position);
	unpack_buff(this->cytoplasmic,                        rcv_buffer, len_buffer, position);
	unpack_buff(this->cytoplasmic_fluid,                  rcv_buffer, len_buffer, position);

	unpack_buff(this->cytoplasmic_solid,                  rcv_buffer, len_buffer, position);
	unpack_buff(this->calcified_fraction,                 rcv_buffer, len_buffer, position);
	unpack_buff(this->cytoplasmic_to_nuclear_ratio,       rcv_buffer, len_buffer, position);

	unpack_buff(this->rupture_volume,                     rcv_buffer, len_buffer, position);
	unpack_buff(this->cytoplasmic_biomass_change_rate,    rcv_buffer, len_buffer, position);
	unpack_buff(this->nuclear_biomass_change_rate,        rcv_buffer, len_buffer, position);

	unpack_buff(this->fluid_change_rate,                  rcv_buffer, len_buffer, position);
	unpack_buff(this->calcification_rate,                 rcv_buffer, len_buffer, position);
	unpack_buff(this->target_solid_cytoplasmic,           rcv_buffer, len_buffer, position);

	unpack_buff(this->target_solid_nuclear,               rcv_buffer, len_buffer, position);
	unpack_buff(this->target_fluid_fraction,              rcv_buffer, len_buffer, position);
	unpack_buff(this->target_cytoplasmic_to_nuclear_ratio,rcv_buffer, len_buffer, position);

	unpack_buff(this->relative_rupture_volume,            rcv_buffer, len_buffer, position);
}

Geometry::Geometry()
{
	// reference values for MCF-7, based on 
	// volume = 2494 cubic microns
	// nuclear volume = 540 cubic microns 
	radius = 8.412710547954228; 
	nuclear_radius = 5.051670902881889; 
	surface_area = 889.3685284131693; 
	
	polarity = 0.0; 
	return; 
}

void Geometry::update_radius( Cell* pCell, Phenotype& phenotype, double dt )
{
	static double four_thirds_pi =  4.188790204786391;
	radius = phenotype.volume.total; 
	radius /= four_thirds_pi; 
	radius = pow( radius , 0.333333333333333333333333333333333333333 ); 
	return; 
}

void Geometry::update_nuclear_radius( Cell* pCell, Phenotype& phenotype, double dt )
{
	static double four_thirds_pi = 4.188790204786391;
	nuclear_radius = phenotype.volume.nuclear; 
	nuclear_radius /= four_thirds_pi; 
	nuclear_radius = pow( nuclear_radius , 0.333333333333333333333333333333333333333 ); 
	return; 
}

void Geometry::update_surface_area( Cell* pCell, Phenotype& phenotype, double dt )
{
	// 4pi / (4pi/3)^(2/3)
	static double the_constant = 4.835975862049409; 
	surface_area = pow( phenotype.volume.total , 0.666666666666667 );
	surface_area /= the_constant; 
	
	return; 
}

void Geometry::update( Cell* pCell, Phenotype& phenotype, double dt )
{
	update_radius(pCell,phenotype,dt); 
	update_nuclear_radius(pCell,phenotype,dt);
	
	// surface area = 4*pi*r^2 = (4/3)*pi*r^3 / (r/3)	
	surface_area = phenotype.volume.total; 
	surface_area /= radius; 
	surface_area *= 3.0; 
	return; 
}

void Geometry::pack(std::vector<char>& snd_buffer, int& len_buffer, int &position){
	 /* Next packing data members of class Geometry as this class has an object in class Phenotype */
    /* 4 doubles are to be packed */

    pack_buff(this->radius,         snd_buffer, len_buffer, position);
	pack_buff(this->nuclear_radius, snd_buffer, len_buffer, position);
	pack_buff(this->surface_area,   snd_buffer, len_buffer, position);
	pack_buff(this->polarity,       snd_buffer, len_buffer, position);
}

void Geometry::unpack(std::vector<char>& rcv_buffer, int& len_buffer, int &position){

	// Geometry: 4 doubles
	unpack_buff(this->radius,         rcv_buffer, len_buffer, position);
	unpack_buff(this->nuclear_radius, rcv_buffer, len_buffer, position);
	unpack_buff(this->surface_area,   rcv_buffer, len_buffer, position);
	unpack_buff(this->polarity,       rcv_buffer, len_buffer, position);
}
	
Mechanics::Mechanics()
{
	cell_cell_adhesion_strength = 0.4; 
	cell_BM_adhesion_strength = 4.0;
	
	cell_cell_repulsion_strength = 10.0; 
	cell_BM_repulsion_strength = 100.0; 

	cell_adhesion_affinities = {1}; 
	
	// this is a multiple of the cell (equivalent) radius
	relative_maximum_adhesion_distance = 1.25; 
	// maximum_adhesion_distance = 0.0; 

	/* for spring attachments */
	maximum_number_of_attachments = 12;
	attachment_elastic_constant = 0.01; 

	attachment_rate = 0; // 10.0 prior ot March 2023
	detachment_rate = 0; 

	/* to be deprecated */ 
	relative_maximum_attachment_distance = relative_maximum_adhesion_distance;
	relative_detachment_distance = relative_maximum_adhesion_distance;

	maximum_attachment_rate = 1.0; 
		
	return; 
}

void Mechanics::sync_to_cell_definitions()
{
	extern std::unordered_map<std::string,int> cell_definition_indices_by_name; 
	int number_of_cell_defs = cell_definition_indices_by_name.size(); 
	
	if( cell_adhesion_affinities.size() != number_of_cell_defs )
	{ cell_adhesion_affinities.resize( number_of_cell_defs, 1.0); }
	return; 
}

double& Mechanics::cell_adhesion_affinity( std::string type_name )
{
	extern std::unordered_map<std::string,int> cell_definition_indices_by_name; 
	int n = cell_definition_indices_by_name[type_name]; 
	return cell_adhesion_affinities[n]; 
}

void Mechanics::set_fully_heterotypic( void )
{
	extern std::unordered_map<std::string,int> cell_definition_indices_by_name; 
	int number_of_cell_defs = cell_definition_indices_by_name.size(); 	

	cell_adhesion_affinities.assign( number_of_cell_defs, 1.0);
	return; 
}

void Mechanics::set_fully_homotypic( Cell* pC )
{
	extern std::unordered_map<std::string,int> cell_definition_indices_by_name; 
	int number_of_cell_defs = cell_definition_indices_by_name.size(); 	

	cell_adhesion_affinities.assign( number_of_cell_defs, 0.0);

	// now find my type and set to 1 
//	cell_adhesion_affinity( pC->type_name ) = 1.0; 

	return; 
}


// new on July 29, 2018
// change the ratio without changing the repulsion strength or equilibrium spacing 
void Mechanics::set_relative_maximum_adhesion_distance( double new_value )
{
	// get old equilibrium spacing, based on equilibriation of pairwise adhesive/repulsive forces at that distance. 
	
		// relative equilibrium spacing (relative to mean cell radius)
	double s_relative = 2.0; 
	
	double temp1 = cell_cell_adhesion_strength; 
	temp1 /= cell_cell_repulsion_strength;
	temp1 = sqrt( temp1 ); 
	
	double temp2 = 1.0; 
	temp2 -= temp1; //  1 - sqrt( alpha_CCA / alpha_CCR );
	
	
	s_relative *= temp2; // 2*( 1 - sqrt( alpha_CCA / alpha_CCR ) ); 
	
	temp1 /= relative_maximum_adhesion_distance; // sqrt( alpha_CCA / alpha_CCR)/f;
	temp2 = 1.0; 
	temp2 -= temp1; // 1 - sqrt( alpha_CCA / alpha_CCR )/f;

	s_relative /= temp2; // 2*( 1 - sqrt( alpha_CCA / alpha_CCR ) ) / ( 1-1/f) ; 
	
	// now, adjust the relative max adhesion distance 
	
	relative_maximum_adhesion_distance = new_value; 
	
	// adjust the adhesive coefficient to preserve the old equilibrium distance

	temp1 = s_relative; 
	temp1 /= 2.0; 
	
	temp2 = 1.0;
	temp2 -= temp1; // 1 - s_relative/2.0 
	
	temp1 /= relative_maximum_adhesion_distance; // s_relative/(2*relative_maximum_adhesion_distance); 
	temp1 *= -1.0; // -s_relative/(2*relative_maximum_adhesion_distance); 
	temp1 += 1.0; // 1.0 -s_relative/(2*relative_maximum_adhesion_distance); 
	
	temp2 /= temp1; 
	temp2 *= temp2; 
	
	cell_cell_adhesion_strength = cell_cell_repulsion_strength;
	cell_cell_adhesion_strength *= temp2; 

	return; 
}		
		
// new on July 29, 2018
// set the cell-cell equilibrium spacing, accomplished by changing the 
// cell-cell adhesion strength, while leaving the cell-cell repulsion 
// strength and the maximum adhesion distance unchanged 
void Mechanics::set_relative_equilibrium_distance( double new_value )
{
	if( new_value > 2.0 )
	{
		std::cout << "**** Warning in function " << __FUNCTION__ << " in " << __FILE__ << " : " << std::endl 
			<< "\tAttempted to set equilibrium distance exceeding two cell radii." << std::endl
			<< "\tWe will cap the equilibrium distance at 2.0 cell radii." << std::endl 
			<< "****" << std::endl << std::endl; 
			
			new_value = 2.0; 
	}
 
	// adjust the adhesive coefficient to achieve the new (relative) equilibrium distance

	double temp1 = new_value; 
	temp1 /= 2.0; 
	
	double temp2 = 1.0;
	temp2 -= temp1; // 1 - s_relative/2.0 
	
	temp1 /= relative_maximum_adhesion_distance; // s_relative/(2*relative_maximum_adhesion_distance); 
	temp1 *= -1.0; // -s_relative/(2*relative_maximum_adhesion_distance); 
	temp1 += 1.0; // 1.0 -s_relative/(2*relative_maximum_adhesion_distance); 
	
	temp2 /= temp1; 
	temp2 *= temp2; 
	
	cell_cell_adhesion_strength = cell_cell_repulsion_strength;
	cell_cell_adhesion_strength *= temp2; 

	return; 
}

void Mechanics::set_absolute_equilibrium_distance( Phenotype& phenotype, double new_value )
{
	return set_relative_equilibrium_distance( new_value / phenotype.geometry.radius ); 
}

void Mechanics::pack(std::vector<char>& snd_buffer, int& len_buffer, int &position){
	/* Next packing object of class Mechanics as class Phenotype contains object of this class */
    /* 5 doubles are to be packed, now v1.8 has added 4 more doubles and 1 int */
	pack_buff(this->cell_cell_adhesion_strength,              snd_buffer, len_buffer, position);
	pack_buff(this->cell_BM_adhesion_strength,               snd_buffer, len_buffer, position);
	pack_buff(this->cell_cell_repulsion_strength,            snd_buffer, len_buffer, position);
	pack_buff(this->cell_BM_repulsion_strength,              snd_buffer, len_buffer, position);
	pack_buff(this->relative_maximum_adhesion_distance,      snd_buffer, len_buffer, position);
	
	pack_buff(this->attachment_rate,                         snd_buffer, len_buffer, position);
	pack_buff(this->detachment_rate,                        snd_buffer, len_buffer, position);
	
	pack_buff(this->relative_maximum_attachment_distance,    snd_buffer, len_buffer, position);
	pack_buff(this->relative_detachment_distance,            snd_buffer, len_buffer, position);
	pack_buff(this->attachment_elastic_constant,             snd_buffer, len_buffer, position);
	pack_buff(this->maximum_attachment_rate,                 snd_buffer, len_buffer, position);
	
	pack_buff(this->maximum_number_of_attachments,           snd_buffer, len_buffer, position);
   
    pack_buff(this->cell_adhesion_affinities, snd_buffer, len_buffer, position);
}

void Mechanics::unpack(std::vector<char>& rcv_buffer, int& len_buffer, int &position){
	// Mechanics: 5 original + 6 v1.8 new = 11 total (10 doubles + 1 int)
	unpack_buff(this->cell_cell_adhesion_strength,           rcv_buffer, len_buffer, position);
	unpack_buff(this->cell_BM_adhesion_strength,            rcv_buffer, len_buffer, position);
	unpack_buff(this->cell_cell_repulsion_strength,         rcv_buffer, len_buffer, position);
	unpack_buff(this->cell_BM_repulsion_strength,           rcv_buffer, len_buffer, position);
	unpack_buff(this->relative_maximum_adhesion_distance,   rcv_buffer, len_buffer, position);

	// New in v1.8
	unpack_buff(this->attachment_rate,                      rcv_buffer, len_buffer, position);
	unpack_buff(this->detachment_rate,                     rcv_buffer, len_buffer, position);
	unpack_buff(this->relative_maximum_attachment_distance, rcv_buffer, len_buffer, position);
	unpack_buff(this->relative_detachment_distance,         rcv_buffer, len_buffer, position);
	unpack_buff(this->attachment_elastic_constant,          rcv_buffer, len_buffer, position);
	unpack_buff(this->maximum_attachment_rate,              rcv_buffer, len_buffer, position);

	// 1 int
	unpack_buff(this->maximum_number_of_attachments,        rcv_buffer, len_buffer, position);

	// Unpack Phenotype::mechanics::cell_adhesion_affinities
	unpack_buff(this->cell_adhesion_affinities, rcv_buffer, len_buffer, position);
}
	
Motility::Motility()
{
	is_motile = false; 
	
	persistence_time = 1.0;
	migration_speed = 1.0;
	
	migration_bias_direction.resize( 3 , 0.0 ); 
	migration_bias = 0.0; 
		
	restrict_to_2D = false; 
	
	// update_migration_bias_direction = NULL; 
	
	motility_vector.resize( 3 , 0.0 ); 
	
	/*-----------------------------------------------------------*/
	/* Gaurav Saxena added the following 2 lines because of v1.7 */
	/*-----------------------------------------------------------*/
	
	chemotaxis_index = 0; 
	chemotaxis_direction = 1;
	
	sync_to_current_microenvironment(); 

	return; 
}

void Motility::sync_to_current_microenvironment( void )
{
	Microenvironment* pMicroenvironment = get_default_microenvironment(); 
	if( pMicroenvironment )
	{ sync_to_microenvironment( pMicroenvironment ); } 
	else
	{ chemotactic_sensitivities.resize( 1 , 0.0 ); }

	return; 
}

double& Motility::chemotactic_sensitivity( std::string name )
{
	int n = microenvironment.find_density_index(name); 
	return chemotactic_sensitivities[n]; 
}

void Motility::sync_to_microenvironment( Microenvironment* pNew_Microenvironment )
{
	chemotactic_sensitivities.resize( pNew_Microenvironment->number_of_densities() , 0.0 ); 
	return; 
}

void Motility::pack(std::vector<char>& snd_buffer, int& len_buffer, int &position){
	 /* Next Packing data members of class Motility */
    /* First bool value is to be packed - allocate sizeof(int) */

	pack_buff(this->is_motile, snd_buffer, len_buffer, position);

    /* Now pack 2 doubles */

    pack_buff(this->persistence_time, snd_buffer, len_buffer, position);
	pack_buff(this->migration_speed,  snd_buffer, len_buffer, position);

    pack_buff(this->migration_bias_direction, snd_buffer, len_buffer, position);

	pack_buff(this->restrict_to_2D, snd_buffer, len_buffer, position);

	pack_buff(this->motility_vector, snd_buffer, len_buffer, position);
    /*================================================================*/
    /* Version 1.7 introduces 2 new 'int' variables in class Motility */
    /* Now these will be packed 																			*/
    /*================================================================*/

    pack_buff(this->chemotaxis_index, snd_buffer, len_buffer, position);
	pack_buff(this->chemotaxis_direction, snd_buffer, len_buffer, position);
	
	//Pack vector<double> chemotactic_sensitivities
	pack_buff(this->chemotactic_sensitivities, snd_buffer, len_buffer, position);
}

void Motility::unpack(std::vector<char>& rcv_buffer, int& len_buffer, int &position){
	unpack_buff(this->is_motile, rcv_buffer, len_buffer, position);
	unpack_buff(this->persistence_time, rcv_buffer, len_buffer, position);
	unpack_buff(this->migration_speed,  rcv_buffer, len_buffer, position);
	unpack_buff(this->migration_bias_direction, rcv_buffer, len_buffer, position);
	unpack_buff(this->restrict_to_2D, rcv_buffer, len_buffer, position);
	unpack_buff(this->motility_vector, rcv_buffer, len_buffer, position);
	unpack_buff(this->chemotaxis_index, rcv_buffer, len_buffer, position);
	unpack_buff(this->chemotaxis_direction, rcv_buffer, len_buffer, position);
	unpack_buff(this->chemotactic_sensitivities, rcv_buffer, len_buffer, position);
}

Secretion::Secretion()
{
	pMicroenvironment = get_default_microenvironment(); 
	
	sync_to_current_microenvironment(); 
	return; 
}

void Secretion::sync_to_current_microenvironment( void )
{
	if( pMicroenvironment )
	{
		sync_to_microenvironment( pMicroenvironment ); 
	}
	else
	{
		secretion_rates.resize( 0 , 0.0 ); 
		uptake_rates.resize( 0 , 0.0 ); 
		saturation_densities.resize( 0 , 0.0 );
		net_export_rates.resize( 0, 0.0 );
	}
	return; 
}
	
void Secretion::sync_to_microenvironment( Microenvironment* pNew_Microenvironment )
{
	pMicroenvironment = pNew_Microenvironment;
	
	secretion_rates.resize( pMicroenvironment->number_of_densities() , 0.0 ); 
	uptake_rates.resize( pMicroenvironment->number_of_densities() , 0.0 ); 
	saturation_densities.resize( pMicroenvironment->number_of_densities() , 0.0 ); 
	net_export_rates.resize( pMicroenvironment->number_of_densities() , 0.0 );
	
	
	return; 
}

void Secretion::advance( Basic_Agent* pCell, Phenotype& phenotype , double dt )
{
	// if this phenotype is not associated with a cell, exit 
	if( pCell == NULL )
	{ return; }

	// if there is no microenvironment, attempt to sync. 
	if( pMicroenvironment == NULL )
	{
		// first, try the cell's microenvironment
		if( pCell->get_microenvironment() )
		{
			sync_to_microenvironment( pCell->get_microenvironment() ); 
		}
		// otherwise, try the default microenvironment
		else
		{
			sync_to_microenvironment( get_default_microenvironment() ); 
		}

		// if we've still failed, return. 
		if( pMicroenvironment == NULL ) 
		{
			return; 
		}
	}

	// make sure the associated cell has the correct rate vectors 
	if( pCell->secretion_rates != &secretion_rates )
	{
		delete pCell->secretion_rates; 
		delete pCell->uptake_rates; 
		delete pCell->saturation_densities;
		
		/*----------------------------------------------------------------------*/
		/* Gaurav Saxena added the following statement as it is present in v1.7 */
		/*----------------------------------------------------------------------*/

		delete pCell->net_export_rates; 
		
		pCell->secretion_rates = &secretion_rates; 
		pCell->uptake_rates = &uptake_rates; 
		pCell->saturation_densities = &saturation_densities;
		
		/*----------------------------------------------------------------------*/
		/* Gaurav Saxena added the following statement as it is present in v1.7 */
		/*----------------------------------------------------------------------*/		
		
		pCell->net_export_rates = &net_export_rates; 
		
		pCell->set_total_volume( phenotype.volume.total ); 
		pCell->set_internal_uptake_constants( dt );
	}

	// now, call the BioFVM secretion/uptake function 
	
	pCell->simulate_secretion_and_uptake( pMicroenvironment , dt ); 
	
	return; 
}

void Secretion::set_all_secretion_to_zero( void )
{
	for( int i=0; i < secretion_rates.size(); i++ )
	{ 
		secretion_rates[i] = 0.0; 
		net_export_rates[i] = 0.0;
	}
	return; 
}

void Secretion::set_all_uptake_to_zero( void )
{
	for( int i=0; i < uptake_rates.size(); i++ )
	{ uptake_rates[i] = 0.0; }
	return; 
}

void Secretion::scale_all_secretion_by_factor( double factor )
{
	for( int i=0; i < secretion_rates.size(); i++ )
	{ 
		secretion_rates[i]  *= factor; 
		net_export_rates[i] *= factor;
	}
	return; 
}

void Secretion::scale_all_uptake_by_factor( double factor )
{
	for( int i=0; i < uptake_rates.size(); i++ )
	{ uptake_rates[i] *= factor; }
	return; 
}

// ease of access
double& Secretion::secretion_rate( std::string name )
{
	int index = microenvironment.find_density_index(name); 
	return secretion_rates[index]; 
}

double& Secretion::uptake_rate( std::string name ) 
{
	int index = microenvironment.find_density_index(name); 
	return uptake_rates[index]; 
}

double& Secretion::saturation_density( std::string name ) 
{
	int index = microenvironment.find_density_index(name); 
	return saturation_densities[index]; 
}

double& Secretion::net_export_rate( std::string name )  
{
	int index = microenvironment.find_density_index(name); 
	return net_export_rates[index]; 
}

void Secretion::pack(std::vector<char>& snd_buffer, int& len_buffer, int &position){
	/* Now packing data members of class Secretion as its object is a data member of class Phenotype */
	pack_buff(this->secretion_rates, snd_buffer, len_buffer, position);
	pack_buff(this->uptake_rates, snd_buffer, len_buffer, position);
	pack_buff(this->saturation_densities, snd_buffer, len_buffer, position);
	pack_buff(this->net_export_rates, snd_buffer, len_buffer, position);
}

void Secretion::unpack(std::vector<char>& rcv_buffer, int& len_buffer, int &position){

	// Secretion
	unpack_buff(this->secretion_rates, rcv_buffer, len_buffer, position);
	unpack_buff(this->uptake_rates, rcv_buffer, len_buffer, position);
	unpack_buff(this->saturation_densities, rcv_buffer, len_buffer, position);
	unpack_buff(this->net_export_rates, rcv_buffer, len_buffer, position);
}


Molecular::Molecular()
{
	pMicroenvironment = get_default_microenvironment(); 
	sync_to_current_microenvironment(); 

	return; 
}

void Molecular::sync_to_current_microenvironment( void )
{
	if( pMicroenvironment )
	{
		sync_to_microenvironment( pMicroenvironment ); 
	}
	else
	{
		internalized_total_substrates.resize( 0 , 0.0 ); 
		fraction_released_at_death.resize( 0 , 0.0 ); 
		fraction_transferred_when_ingested.resize( 0, 1.0 ); 
	}
	return; 
}
	
void Molecular::sync_to_microenvironment( Microenvironment* pNew_Microenvironment )
{
	pMicroenvironment = pNew_Microenvironment;
	
	int number_of_densities = pMicroenvironment->number_of_densities() ; 

	internalized_total_substrates.resize( number_of_densities , 0.0 ); 
	fraction_released_at_death.resize( number_of_densities , 0.0 ); 
	fraction_transferred_when_ingested.resize( number_of_densities , 1.0 ); 
	
	return; 
}

void Molecular::sync_to_cell( Basic_Agent* pCell )
{
	if( pCell->internalized_substrates != &internalized_total_substrates )
	{
		delete pCell->internalized_substrates;
		pCell->internalized_substrates = &internalized_total_substrates;
	}
	
	if( pCell->fraction_released_at_death != &fraction_released_at_death )
	{
		delete pCell->fraction_released_at_death;
		pCell->fraction_released_at_death = &fraction_released_at_death; 
	}
	
	if( pCell->fraction_transferred_when_ingested != &fraction_transferred_when_ingested )
	{
		delete pCell->fraction_transferred_when_ingested; 
		pCell->fraction_transferred_when_ingested = &fraction_transferred_when_ingested; 
	}

	return; 
}

// ease of access 
double&  Molecular::internalized_total_substrate( std::string name )
{
	int index = microenvironment.find_density_index(name); 
	return internalized_total_substrates[index]; 
}

void Molecular::pack(std::vector<char>& snd_buffer, int& len_buffer, int &position){

	pack_buff(this->internalized_total_substrates, snd_buffer, len_buffer, position);

	pack_buff(this->fraction_released_at_death, snd_buffer, len_buffer, position);

	pack_buff(this->fraction_transferred_when_ingested, snd_buffer, len_buffer, position);
}

void Molecular::unpack(std::vector<char>& rcv_buffer, int& len_buffer, int &position){
	// Molecular
	unpack_buff(this->internalized_total_substrates, rcv_buffer, len_buffer, position);
	unpack_buff(this->fraction_released_at_death, rcv_buffer, len_buffer, position);
	unpack_buff(this->fraction_transferred_when_ingested, rcv_buffer, len_buffer, position);
}

/*
void Molecular::advance( Basic_Agent* pCell, Phenotype& phenotype , double dt )
{
	// if this phenotype is not associated with a cell, exit 
	if( pCell == NULL )
	{ return; }

	// if there is no microenvironment, attempt to sync. 
	if( pMicroenvironment == NULL )
	{
		// first, try the cell's microenvironment
		if( pCell->get_microenvironment() )
		{
			sync_to_microenvironment( pCell->get_microenvironment() ); 
		}
		// otherwise, try the default microenvironment
		else
		{
			sync_to_microenvironment( get_default_microenvironment() ); 
		}

		// if we've still failed, return. 
		if( pMicroenvironment == NULL ) 
		{
			return; 
		}
	}

	// make sure the associated cell has the correct rate vectors 
	if( pCell->internalized_substrates != &internalized_substrates )
	{
		// copy the data over 
		internalized_substrates = *(pCell->internalized_substrates);
		// remove the BioFVM copy 
		delete pCell->internalized_substrates; 
		// point BioFVM to this one  
		pCell->internalized_substrates = &internalized_substrates; 
	}

	// now, call the functions 
//	if( pCell->functions.internal_substrate_function )
//	{ pCell->functions.internal_substrate_function( pCell,phenotype,dt);  }
//	if( pCell->functions.molecular_model_function )
//	{ pCell->functions.molecular_model_function( pCell,phenotype,dt);  }


	return; 
}
*/


Cell_Functions::Cell_Functions()
{
	instantiate_cell = NULL;
	
	volume_update_function = NULL; 
	update_migration_bias = NULL; 
	
	update_phenotype = NULL; 
	custom_cell_rule = NULL; 
	
	update_velocity = NULL;
	
	/*------------------------------------------------------------------------------------*/
	/* For completion making update_velocity_parallel = NULL 															*/
	/*------------------------------------------------------------------------------------*/
	
	update_velocity_parallel = NULL; 
	 
	add_cell_basement_membrane_interactions = NULL; 
	calculate_distance_to_membrane = NULL; 
	
	set_orientation = NULL; 
	
	contact_function = NULL; 

/*	
	internal_substrate_function = NULL; 
	molecular_model_function = NULL; 
*/
	return; 
}

void Cell_Functions::pack(std::vector<char>& snd_buffer, int& len_buffer, int &position){

	//Only attribute is Cycle_Model cycle_model
	this->cycle_model.pack(snd_buffer, len_buffer, position);
}

void Cell_Functions::unpack(std::vector<char>& rcv_buffer, int& len_buffer, int &position){

	//Only attribute is Cycle_Model cycle_model
	this->cycle_model.unpack(rcv_buffer, len_buffer, position);
}

void Cell_Integrity::pack(std::vector<char>& snd_buffer, int& len_buffer, int& position){
	pack_buff(this->damage, snd_buffer, len_buffer, position);
	pack_buff(this->damage_rate, snd_buffer, len_buffer, position);
	pack_buff(this->damage_repair_rate, snd_buffer, len_buffer, position);
}

void Cell_Integrity::unpack(std::vector<char>& rcv_buffer, int& len_buffer, int& position){
	unpack_buff(this->damage, rcv_buffer, len_buffer, position);
	unpack_buff(this->damage_rate, rcv_buffer, len_buffer, position);
	unpack_buff(this->damage_repair_rate, rcv_buffer, len_buffer, position);
}

void Phenotype::sync_to_functions( Cell_Functions& functions )
{
	cycle.sync_to_cycle_model( functions.cycle_model );  
	
	return; 
}

Phenotype::Phenotype()
{
	flagged_for_division = false;
	flagged_for_removal = false; 
	
	intracellular = NULL;
	// sync the molecular stuff here automatically? 
	
	return; 
}


Phenotype::Phenotype(const Phenotype &p) 
{
	intracellular = NULL;
	*this = p;
}

Phenotype::~Phenotype() 
{
	if (intracellular != NULL)
		delete intracellular;
}

Phenotype& Phenotype::operator=(const Phenotype &p ) 
{ 
	if( this == &p )
	{ return *this; }

	Intracellular* new_intracellular = NULL;
	if( p.intracellular != NULL )
	{ new_intracellular = p.intracellular->clone(); }

		
	flagged_for_division = p.flagged_for_division;
	flagged_for_removal = p.flagged_for_removal;
	
	cycle = p.cycle;
	death = p.death;
	volume = p.volume;
	geometry = p.geometry;
	mechanics = p.mechanics;
	motility = p.motility;
	secretion = p.secretion;

	molecular = p.molecular;

	cell_integrity = p.cell_integrity; 
	
	delete intracellular;
	intracellular = new_intracellular;

	cell_interactions = p.cell_interactions; 
	cell_transformations = p.cell_transformations; 	
	
	return *this;
}

void Phenotype::pack(std::vector<char>& snd_buffer, int& len_buffer, int& position){

	pack_buff( this->flagged_for_division, snd_buffer, len_buffer, position);
    
	pack_buff( this->flagged_for_removal, snd_buffer, len_buffer, position);

	this->cycle.pack(snd_buffer, len_buffer, position);

	this->death.pack(snd_buffer, len_buffer, position);

	this->volume.pack(snd_buffer, len_buffer, position);

	this->geometry.pack(snd_buffer, len_buffer, position);

	this->mechanics.pack(snd_buffer, len_buffer, position);

	this->motility.pack(snd_buffer, len_buffer, position);

	this->secretion.pack(snd_buffer, len_buffer, position);

	this->molecular.pack(snd_buffer, len_buffer, position);

	this->cell_integrity.pack(snd_buffer, len_buffer, position);

	std::string temp_str;
	//Time to pack intracelluar
	 if (this->intracellular != NULL) 											//this is NULL and NOT nullptr
    {
        temp_str = this->intracellular->intracellular_type; 				//This was the original code
    }
    else
        temp_str = "not-maboss";
                
	pack_buff(temp_str, snd_buffer, len_buffer, position);
			
#ifdef ADDON_PHYSIBOSS

	if(this->intracellular != NULL)			//Added by Gaurav Saxena
		if (this->intracellular->intracellular_type.compare("maboss") == 0) 
		{
			MaBoSSIntracellular* t_intracellular = static_cast<MaBoSSIntracellular*>(this->intracellular); 
			t_intracellular->pack(snd_buffer, len_buffer, position);	
		}
				
#endif

	this->cell_interactions.pack(snd_buffer, len_buffer, position);

	this->cell_transformations.pack(snd_buffer, len_buffer, position);


}

void Phenotype::unpack(std::vector<char>& rcv_buffer, int& len_buffer, int& position){

	unpack_buff( this->flagged_for_division, rcv_buffer, len_buffer, position);
    
	unpack_buff( this->flagged_for_removal, rcv_buffer, len_buffer, position);

	this->cycle.unpack(rcv_buffer, len_buffer, position);

	this->death.unpack(rcv_buffer, len_buffer, position);

	this->volume.unpack(rcv_buffer, len_buffer, position);

	this->geometry.unpack(rcv_buffer, len_buffer, position);

	this->mechanics.unpack(rcv_buffer, len_buffer, position);

	this->motility.unpack(rcv_buffer, len_buffer, position);

	this->secretion.unpack(rcv_buffer, len_buffer, position);

	this->molecular.unpack(rcv_buffer, len_buffer, position);

	this->cell_integrity.unpack(rcv_buffer, len_buffer, position);

	std::string temp_str;
                
	unpack_buff(temp_str, rcv_buffer, len_buffer, position);

	if (temp_str == "not-maboss")
	{
		if (this->intracellular != NULL)
		{
			delete this->intracellular;
			this->intracellular = NULL;
		}
	}
	else
	{
#ifdef ADDON_PHYSIBOSS
		if (temp_str.compare("maboss") == 0)
		{
			if (this->intracellular == NULL || this->intracellular->intracellular_type.compare("maboss") != 0)
			{
				delete this->intracellular;
				this->intracellular = new MaBoSSIntracellular();
			}

			MaBoSSIntracellular* t_intracellular = static_cast<MaBoSSIntracellular*>(this->intracellular);
			t_intracellular->unpack(rcv_buffer, len_buffer, position);
		}
#endif
	}

	this->cell_interactions.unpack(rcv_buffer, len_buffer, position);

	this->cell_transformations.unpack(rcv_buffer, len_buffer, position);


}

Asymmetric_Division::Asymmetric_Division()
{
	asymmetric_division_probabilities = {0.0};
}

void Asymmetric_Division::sync_to_cell_definitions()
{
	extern std::unordered_map<std::string,int> cell_definition_indices_by_name; 
	int number_of_cell_defs = cell_definition_indices_by_name.size(); 
	
	if( asymmetric_division_probabilities.size() != number_of_cell_defs )
	{ asymmetric_division_probabilities.resize( number_of_cell_defs, 0.0); }
	
	return; 
}

double Asymmetric_Division::probabilities_total( void )
{
	double total = 0.0; 
	for( int i=0; i < asymmetric_division_probabilities.size(); i++ )
	{ total += asymmetric_division_probabilities[i]; }
	return total; 
}

void Asymmetric_Division::pack(std::vector<char>& snd_buffer, int& len_buffer, int &position) {

	pack_buff(this->asymmetric_division_probabilities, snd_buffer, len_buffer, position);
}

void Asymmetric_Division::unpack(std::vector<char>& rcv_buffer, int& len_buffer, int &position) {

	unpack_buff(this->asymmetric_division_probabilities, rcv_buffer, len_buffer, position);
}

// ease of access
double& Asymmetric_Division::asymmetric_division_probability( std::string type_name )
{
	extern std::unordered_map<std::string,int> cell_definition_indices_by_name; 
	int n = cell_definition_indices_by_name[type_name]; 
	return asymmetric_division_probabilities[n]; 
}

/*
class Bools
{
	public:
		std::vector<bool> values; 
		std::unordered_map<std::string,int> name_map; 
		std::string& name( int i ); 
		std::vector<std::string> units; 
		void resize( int n ); 
		int add( std::string name , std::string units , bool value ); 
		
		bool& operator[]( int i ); 
		bool& operator[]( std::string name ); 
		
		Bools(); 
}
*/

Bools::Bools()
{
	values.resize( 0 , true ); 
	name_map.clear(); 
	return; 
}

int Bools::size( void )
{ return values.size(); } 


void Phenotype::sync_to_microenvironment( Microenvironment* pMicroenvironment )
{
	secretion.sync_to_microenvironment( pMicroenvironment ); 
	molecular.sync_to_microenvironment( pMicroenvironment ); 

	return; 
}

Cell_Interactions::Cell_Interactions()
{
	// dead_phagocytosis_rate = 0.0; 

	apoptotic_phagocytosis_rate = 0.0; 
	necrotic_phagocytosis_rate = 0.0; 
	other_dead_phagocytosis_rate = 0.0; 

	live_phagocytosis_rates = {0.0}; 

	attack_damage_rate = 1.0; 
	attack_rates = {0.0}; 
	immunogenicities = {1}; 

	pAttackTarget = NULL; 
	total_damage_delivered = 0.0; 

	attack_duration = 30.0; // a typical attack duration for a T cell using perforin/granzyme is ~30 minutes

	fusion_rates = {0.0}; 
	
	return; 
}

void Cell_Interactions::sync_to_cell_definitions()
{
	extern std::unordered_map<std::string,int> cell_definition_indices_by_name; 
	int number_of_cell_defs = cell_definition_indices_by_name.size(); 

	if( live_phagocytosis_rates.size() != number_of_cell_defs )
	{
		live_phagocytosis_rates.resize( number_of_cell_defs, 0.0); 
		attack_rates.resize( number_of_cell_defs, 0.0); 
		fusion_rates.resize( number_of_cell_defs, 0.0); 
		immunogenicities.resize( number_of_cell_defs , 1.0 ); 
	}

	return; 
}

// ease of access 
double& Cell_Interactions::live_phagocytosis_rate( std::string type_name )
{
	extern std::unordered_map<std::string,int> cell_definition_indices_by_name; 
	int n = cell_definition_indices_by_name[type_name]; 
	return live_phagocytosis_rates[n]; 
}

double& Cell_Interactions::attack_rate( std::string type_name ) 
{
	extern std::unordered_map<std::string,int> cell_definition_indices_by_name; 
	int n = cell_definition_indices_by_name[type_name]; 
	return attack_rates[n]; 
}

double& Cell_Interactions::fusion_rate( std::string type_name )
{
	extern std::unordered_map<std::string,int> cell_definition_indices_by_name; 
	int n = cell_definition_indices_by_name[type_name]; 
	return fusion_rates[n]; 
}

double& Cell_Interactions::immunogenicity( std::string type_name )
{
	extern std::unordered_map<std::string,int> cell_definition_indices_by_name; 
	int n = cell_definition_indices_by_name[type_name]; 
	return immunogenicities[n]; 
}

void Cell_Interactions::pack(std::vector<char>& snd_buffer, int& len_buffer, int& position){
	pack_buff(this->apoptotic_phagocytosis_rate, snd_buffer, len_buffer, position);
	pack_buff(this->necrotic_phagocytosis_rate, snd_buffer, len_buffer, position);
	pack_buff(this->other_dead_phagocytosis_rate, snd_buffer, len_buffer, position);
	pack_buff(this->live_phagocytosis_rates, snd_buffer, len_buffer, position);
	pack_buff(this->attack_rates, snd_buffer, len_buffer, position);
	pack_buff(this->immunogenicities, snd_buffer, len_buffer, position);
	pack_buff(this->attack_damage_rate, snd_buffer, len_buffer, position);
	pack_buff(this->total_damage_delivered, snd_buffer, len_buffer, position);
	pack_buff(this->attack_duration, snd_buffer, len_buffer, position);
	pack_buff(this->fusion_rates, snd_buffer, len_buffer, position); 
}

void Cell_Interactions::unpack(std::vector<char>& rcv_buffer, int& len_buffer, int& position){
	unpack_buff(this->apoptotic_phagocytosis_rate, rcv_buffer, len_buffer, position);
	unpack_buff(this->necrotic_phagocytosis_rate, rcv_buffer, len_buffer, position);
	unpack_buff(this->other_dead_phagocytosis_rate, rcv_buffer, len_buffer, position);
	unpack_buff(this->live_phagocytosis_rates, rcv_buffer, len_buffer, position);
	unpack_buff(this->attack_rates, rcv_buffer, len_buffer, position);
	unpack_buff(this->immunogenicities, rcv_buffer, len_buffer, position);
	unpack_buff(this->attack_damage_rate, rcv_buffer, len_buffer, position);
	unpack_buff(this->total_damage_delivered, rcv_buffer, len_buffer, position);
	unpack_buff(this->attack_duration, rcv_buffer, len_buffer, position);
	unpack_buff(this->fusion_rates, rcv_buffer, len_buffer, position); 
}

Cell_Transformations::Cell_Transformations()
{
	transformation_rates = {0.0}; 

	return; 
}

void Cell_Transformations::sync_to_cell_definitions()
{
	extern std::unordered_map<std::string,int> cell_definition_indices_by_name; 	
	int number_of_cell_defs = cell_definition_indices_by_name.size();  

	if( transformation_rates.size() != number_of_cell_defs )
	{ transformation_rates.resize( number_of_cell_defs, 0.0); }

	return; 
}

// ease of access 
double& Cell_Transformations::transformation_rate( std::string type_name )
{
	extern std::unordered_map<std::string,int> cell_definition_indices_by_name; 
	int n = cell_definition_indices_by_name[type_name]; 
	return transformation_rates[n]; 
}

void Cell_Transformations::pack(std::vector<char>& snd_buffer, int& len_buffer, int& position){
	pack_buff(this->transformation_rates, snd_buffer, len_buffer, position);
}

void Cell_Transformations::unpack(std::vector<char>& snd_buffer, int& len_buffer, int& position){
	unpack_buff(this->transformation_rates, snd_buffer, len_buffer, position);
}

// beta functionality in 1.10.3 
Cell_Integrity::Cell_Integrity()
{
 	damage = 0;  
	damage_rate = 0.0; 
	damage_repair_rate = 0.0; 

	return; 
}

void Cell_Integrity::advance_damage( double dt )
{
	double temp1;
	double temp2; 
	static double tol = 1e-8; 

	// general damage 
	if( damage_rate > tol || damage_repair_rate > tol )
	{
		temp1 = dt; 
		temp2 = dt; 
		temp1 *= damage_rate;  
		temp2 *= damage_repair_rate; 
		temp2 += 1; 

		damage += temp1; 
		damage /= temp2; 
	}
	return; 
}

};


