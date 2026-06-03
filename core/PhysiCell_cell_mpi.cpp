#include "PhysiCell_cell.h"
#include "PhysiCell_cell_container.h"
#include "PhysiCell_utilities.h"
#include "PhysiCell_constants.h"
#include "MPI_helper.h"
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

#include <iostream>
#include <iomanip>  // for std::setprecision
#include <algorithm>
#include <cstdlib> // for rand()
#include <ctime>   // for seeding rand()


using namespace PhysiCell;

namespace
{

void sync_secretion_to_cell(Cell* pCell)
{
	if( pCell->secretion_rates != &(pCell->phenotype.secretion.secretion_rates) )
	{
		delete pCell->secretion_rates;
		delete pCell->uptake_rates;
		delete pCell->saturation_densities;
		delete pCell->net_export_rates;

		pCell->secretion_rates = &(pCell->phenotype.secretion.secretion_rates);
		pCell->uptake_rates = &(pCell->phenotype.secretion.uptake_rates);
		pCell->saturation_densities = &(pCell->phenotype.secretion.saturation_densities);
		pCell->net_export_rates = &(pCell->phenotype.secretion.net_export_rates);
	}
}

void restore_cell_from_definition(Cell* pCell, Cell_Definition* pCD)
{
	if( pCD == NULL )
	{
		return;
	}

	pCell->type = pCD->type;
	pCell->type_name = pCD->name;
	pCell->custom_data = pCD->custom_data;
	pCell->parameters = pCD->parameters;
	pCell->functions = pCD->functions;
	pCell->phenotype = pCD->phenotype;
	pCell->is_movable = pCD->is_movable;

	pCell->phenotype.cycle.sync_to_cycle_model( pCD->functions.cycle_model );

	BioFVM::Microenvironment* pMicroenvironment = pCell->get_microenvironment();
	if( pMicroenvironment != NULL )
	{
		pCell->phenotype.sync_to_microenvironment( pMicroenvironment );
	}

	pCell->phenotype.motility.sync_to_current_microenvironment();
	pCell->phenotype.mechanics.sync_to_cell_definitions();
	pCell->phenotype.cell_interactions.sync_to_cell_definitions();
	pCell->phenotype.cell_transformations.sync_to_cell_definitions();
	pCell->phenotype.cycle.asymmetric_division.sync_to_cell_definitions();
	pCell->phenotype.molecular.sync_to_cell( pCell );
	sync_secretion_to_cell( pCell );
}

struct Cycle_Model_Selector
{
	int kind;
	int death_model_index;
	int code;
	std::string name;
};

bool cycle_model_matches(Cycle_Model* pModel, const Cycle_Model_Selector& selector)
{
	if( pModel == NULL )
	{
		return false;
	}

	return pModel->code == selector.code && pModel->name == selector.name;
}

bool cycle_model_matches_loose(Cycle_Model* pModel, const Cycle_Model_Selector& selector)
{
	if( pModel == NULL )
	{
		return false;
	}

	return pModel->code == selector.code || pModel->name == selector.name;
}

Cycle_Model* definition_main_cycle_model(Cell_Definition* pCD)
{
	if( pCD == NULL )
	{
		return NULL;
	}

	return &(pCD->functions.cycle_model);
}

Cycle_Model* definition_death_cycle_model(Cell_Definition* pCD, int death_model_index)
{
	if( pCD == NULL )
	{
		return NULL;
	}

	if( death_model_index < 0 ||
		death_model_index >= static_cast<int>( pCD->phenotype.death.models.size() ) )
	{
		return NULL;
	}

	return pCD->phenotype.death.models[death_model_index];
}

Cycle_Model_Selector make_cycle_model_selector(Cell* pCell)
{
	Cycle_Model_Selector selector;
	selector.kind = 0;
	selector.death_model_index = -1;
	selector.code = -1;
	selector.name.clear();

	if( pCell == NULL )
	{
		return selector;
	}

	Cycle_Model* pModel = pCell->phenotype.cycle.pCycle_Model;
	if( pModel == NULL )
	{
		pModel = &(pCell->functions.cycle_model);
	}

	if( pModel == &(pCell->functions.cycle_model) )
	{
		selector.kind = 0;
	}
	else
	{
		selector.kind = 2;
		for( int i = 0 ; i < static_cast<int>( pCell->phenotype.death.models.size() ) ; i++ )
		{
			if( pCell->phenotype.death.models[i] == pModel )
			{
				selector.kind = 1;
				selector.death_model_index = i;
				break;
			}
		}
	}

	if( pModel != NULL )
	{
		selector.code = pModel->code;
		selector.name = pModel->name;
	}

	return selector;
}

void pack_cycle_model_selector_light(Cell* pCell, std::vector<char>& buffer, int& len_buffer, int& position)
{
	Cycle_Model_Selector selector = make_cycle_model_selector( pCell );
	pack_buff( selector.kind , buffer, len_buffer, position );
	pack_buff( selector.death_model_index , buffer, len_buffer, position );
	pack_buff( selector.code , buffer, len_buffer, position );
	pack_buff( selector.name , buffer, len_buffer, position );
}

Cycle_Model_Selector unpack_cycle_model_selector_light(std::vector<char>& buffer, int& len_buffer, int& position)
{
	Cycle_Model_Selector selector;
	selector.kind = 0;
	selector.death_model_index = -1;
	selector.code = -1;

	unpack_buff( selector.kind , buffer, len_buffer, position );
	unpack_buff( selector.death_model_index , buffer, len_buffer, position );
	unpack_buff( selector.code , buffer, len_buffer, position );
	unpack_buff( selector.name , buffer, len_buffer, position );

	return selector;
}

Cycle_Model* resolve_cycle_model_selector(Cell_Definition* pCD, const Cycle_Model_Selector& selector)
{
	if( pCD == NULL )
	{
		return NULL;
	}

	Cycle_Model* pMainModel = definition_main_cycle_model( pCD );

	if( selector.kind == 1 )
	{
		Cycle_Model* pDeathModel = definition_death_cycle_model( pCD , selector.death_model_index );
		if( pDeathModel != NULL )
		{
			return pDeathModel;
		}
	}

	if( selector.kind == 0 )
	{
		return pMainModel;
	}

	if( cycle_model_matches( pMainModel , selector ) )
	{
		return pMainModel;
	}

	for( int i = 0 ; i < static_cast<int>( pCD->phenotype.death.models.size() ) ; i++ )
	{
		Cycle_Model* pDeathModel = definition_death_cycle_model( pCD , i );
		if( cycle_model_matches( pDeathModel , selector ) )
		{
			return pDeathModel;
		}
	}

	if( cycle_model_matches_loose( pMainModel , selector ) )
	{
		return pMainModel;
	}

	for( int i = 0 ; i < static_cast<int>( pCD->phenotype.death.models.size() ) ; i++ )
	{
		Cycle_Model* pDeathModel = definition_death_cycle_model( pCD , i );
		if( cycle_model_matches_loose( pDeathModel , selector ) )
		{
			return pDeathModel;
		}
	}

	return pMainModel;
}

Cycle_Model* resolve_cycle_model_from_definition_state(Cell* pCell, Cell_Definition* pCD)
{
	if( pCD == NULL )
	{
		return NULL;
	}

	Cycle_Model* pSelectedModel = NULL;
	if( pCell != NULL )
	{
		pSelectedModel = pCell->phenotype.cycle.pCycle_Model;
	}

	Cycle_Model* pMainModel = definition_main_cycle_model( pCD );
	if( pSelectedModel == pMainModel )
	{
		return pMainModel;
	}

	for( int i = 0 ; i < static_cast<int>( pCD->phenotype.death.models.size() ) ; i++ )
	{
		Cycle_Model* pDeathModel = definition_death_cycle_model( pCD , i );
		if( pSelectedModel == pDeathModel )
		{
			return pDeathModel;
		}
	}

	if( pCell != NULL )
	{
		if( pCell->phenotype.death.dead )
		{
			Cycle_Model* pDeathModel =
				definition_death_cycle_model( pCD , pCell->phenotype.death.current_death_model_index );
			if( pDeathModel != NULL )
			{
				return pDeathModel;
			}
		}

		if( pSelectedModel != NULL )
		{
			Cycle_Model_Selector selector;
			selector.kind = 2;
			selector.death_model_index = pCell->phenotype.death.current_death_model_index;
			selector.code = pSelectedModel->code;
			selector.name = pSelectedModel->name;

			Cycle_Model* pResolvedModel = resolve_cycle_model_selector( pCD , selector );
			if( pResolvedModel != NULL )
			{
				return pResolvedModel;
			}
		}
	}

	return pMainModel;
}

void sync_cycle_to_model_preserving_state(Cell* pCell, Cycle_Model* pModel)
{
	if( pCell == NULL )
	{
		return;
	}

	if( pModel == NULL )
	{
		return;
	}

	std::vector<std::vector<double>> received_transition_rates = pCell->phenotype.cycle.data.transition_rates;
	int received_phase_index = pCell->phenotype.cycle.data.current_phase_index;
	double received_elapsed_time = pCell->phenotype.cycle.data.elapsed_time_in_phase;

	pCell->phenotype.cycle.sync_to_cycle_model( *pModel );

	int phase_count = static_cast<int>( pModel->phases.size() );
	int phase_link_count = static_cast<int>( pModel->phase_links.size() );
	int transition_row_count = static_cast<int>( pCell->phenotype.cycle.data.transition_rates.size() );
	int valid_phase_count = std::min( phase_count , std::min( phase_link_count , transition_row_count ) );
	if( received_phase_index >= 0 && received_phase_index < valid_phase_count )
	{
		pCell->phenotype.cycle.data.current_phase_index = received_phase_index;
		pCell->phenotype.cycle.data.elapsed_time_in_phase = std::max( 0.0 , received_elapsed_time );
	}
	else
	{
		if( valid_phase_count > 0 )
		{
			pCell->phenotype.cycle.data.current_phase_index =
				std::min( std::max( 0 , pModel->default_phase_index ) , valid_phase_count - 1 );
		}
		else
		{
			pCell->phenotype.cycle.data.current_phase_index = 0;
		}
		pCell->phenotype.cycle.data.elapsed_time_in_phase = 0.0;
	}

	int row_count = std::min(
		static_cast<int>( received_transition_rates.size() ),
		static_cast<int>( pCell->phenotype.cycle.data.transition_rates.size() ) );
	for( int i = 0 ; i < row_count ; i++ )
	{
		int col_count = std::min(
			static_cast<int>( received_transition_rates[i].size() ),
			static_cast<int>( pCell->phenotype.cycle.data.transition_rates[i].size() ) );
		std::copy(
			received_transition_rates[i].begin(),
			received_transition_rates[i].begin() + col_count,
			pCell->phenotype.cycle.data.transition_rates[i].begin() );
	}
}

void apply_cycle_model_safety_net(Cell* pCell, Cell_Definition* pCD)
{
	if( pCell == NULL )
	{
		return;
	}

	Cycle_Model* pModel = resolve_cycle_model_from_definition_state( pCell , pCD );
	sync_cycle_to_model_preserving_state( pCell , pModel );
}

void pack_custom_data_values_light(Custom_Cell_Data& custom_data, std::vector<char>& buffer, int& len_buffer, int& position)
{
	int scalar_count = static_cast<int>( custom_data.variables.size() );
	pack_buff( scalar_count , buffer, len_buffer, position );
	for( const auto& variable : custom_data.variables )
	{
		pack_buff( variable.value , buffer, len_buffer, position );
	}

	int vector_count = static_cast<int>( custom_data.vector_variables.size() );
	pack_buff( vector_count , buffer, len_buffer, position );
	for( const auto& variable : custom_data.vector_variables )
	{
		pack_buff( variable.value , buffer, len_buffer, position );
	}
}

void unpack_custom_data_values_light(Custom_Cell_Data& custom_data, std::vector<char>& buffer, int& len_buffer, int& position)
{
	int scalar_count = 0;
	unpack_buff( scalar_count , buffer, len_buffer, position );
	for( int i = 0 ; i < scalar_count ; i++ )
	{
		double value = 0.0;
		unpack_buff( value , buffer, len_buffer, position );
		if( i < static_cast<int>( custom_data.variables.size() ) )
		{
			custom_data.variables[i].value = value;
		}
	}

	int vector_count = 0;
	unpack_buff( vector_count , buffer, len_buffer, position );
	for( int i = 0 ; i < vector_count ; i++ )
	{
		std::vector<double> value;
		unpack_buff( value , buffer, len_buffer, position );
		if( i < static_cast<int>( custom_data.vector_variables.size() ) )
		{
			custom_data.vector_variables[i].value = value;
		}
	}
}

void pack_cycle_light(Cycle& cycle, std::vector<char>& buffer, int& len_buffer, int& position)
{
	int transition_count = static_cast<int>( cycle.data.transition_rates.size() );
	pack_buff( transition_count , buffer, len_buffer, position );
	for( const auto& rates : cycle.data.transition_rates )
	{
		pack_buff( rates , buffer, len_buffer, position );
	}

	pack_buff( cycle.data.current_phase_index , buffer, len_buffer, position );
	pack_buff( cycle.data.elapsed_time_in_phase , buffer, len_buffer, position );
	pack_buff( cycle.asymmetric_division.asymmetric_division_probabilities , buffer, len_buffer, position );
}

void unpack_cycle_light(Cycle& cycle, std::vector<char>& buffer, int& len_buffer, int& position)
{
	int transition_count = 0;
	unpack_buff( transition_count , buffer, len_buffer, position );
	for( int i = 0 ; i < transition_count ; i++ )
	{
		std::vector<double> rates;
		unpack_buff( rates , buffer, len_buffer, position );
		if( i < static_cast<int>( cycle.data.transition_rates.size() ) )
		{
			cycle.data.transition_rates[i] = rates;
		}
	}

	unpack_buff( cycle.data.current_phase_index , buffer, len_buffer, position );
	unpack_buff( cycle.data.elapsed_time_in_phase , buffer, len_buffer, position );
	unpack_buff( cycle.asymmetric_division.asymmetric_division_probabilities , buffer, len_buffer, position );
}

void pack_death_light(Death& death, std::vector<char>& buffer, int& len_buffer, int& position)
{
	pack_buff( death.rates , buffer, len_buffer, position );
	pack_buff( death.dead , buffer, len_buffer, position );
	pack_buff( death.current_death_model_index , buffer, len_buffer, position );
}

void unpack_death_light(Death& death, std::vector<char>& buffer, int& len_buffer, int& position)
{
	unpack_buff( death.rates , buffer, len_buffer, position );
	unpack_buff( death.dead , buffer, len_buffer, position );
	unpack_buff( death.current_death_model_index , buffer, len_buffer, position );
}

void pack_motility_light(Motility& motility, std::vector<char>& buffer, int& len_buffer, int& position)
{
	pack_buff( motility.is_motile , buffer, len_buffer, position );
	pack_buff( motility.persistence_time , buffer, len_buffer, position );
	pack_buff( motility.migration_speed , buffer, len_buffer, position );
	pack_buff( motility.migration_bias_direction , buffer, len_buffer, position );
	pack_buff( motility.migration_bias , buffer, len_buffer, position );
	pack_buff( motility.restrict_to_2D , buffer, len_buffer, position );
	pack_buff( motility.motility_vector , buffer, len_buffer, position );
	pack_buff( motility.chemotaxis_index , buffer, len_buffer, position );
	pack_buff( motility.chemotaxis_direction , buffer, len_buffer, position );
	pack_buff( motility.chemotactic_sensitivities , buffer, len_buffer, position );
}

void unpack_motility_light(Motility& motility, std::vector<char>& buffer, int& len_buffer, int& position)
{
	unpack_buff( motility.is_motile , buffer, len_buffer, position );
	unpack_buff( motility.persistence_time , buffer, len_buffer, position );
	unpack_buff( motility.migration_speed , buffer, len_buffer, position );
	unpack_buff( motility.migration_bias_direction , buffer, len_buffer, position );
	unpack_buff( motility.migration_bias , buffer, len_buffer, position );
	unpack_buff( motility.restrict_to_2D , buffer, len_buffer, position );
	unpack_buff( motility.motility_vector , buffer, len_buffer, position );
	unpack_buff( motility.chemotaxis_index , buffer, len_buffer, position );
	unpack_buff( motility.chemotaxis_direction , buffer, len_buffer, position );
	unpack_buff( motility.chemotactic_sensitivities , buffer, len_buffer, position );
}

void pack_intracellular_light(Phenotype& phenotype, std::vector<char>& buffer, int& len_buffer, int& position)
{
	std::string intracellular_type = "not-maboss";
	if( phenotype.intracellular != NULL )
	{
		intracellular_type = phenotype.intracellular->intracellular_type;
	}

	pack_buff( intracellular_type , buffer, len_buffer, position );

#ifdef ADDON_PHYSIBOSS
	if( phenotype.intracellular != NULL && intracellular_type == "maboss" )
	{
		MaBoSSIntracellular* intracellular = static_cast<MaBoSSIntracellular*>( phenotype.intracellular );

		int input_count = static_cast<int>( intracellular->listOfInputs.size() );
		pack_buff( input_count , buffer, len_buffer, position );
		for( const auto& input_pair : intracellular->listOfInputs )
		{
			pack_buff( input_pair.second.smoothed_value , buffer, len_buffer, position );
		}

		int output_count = static_cast<int>( intracellular->listOfOutputs.size() );
		pack_buff( output_count , buffer, len_buffer, position );
		for( const auto& output_pair : intracellular->listOfOutputs )
		{
			pack_buff( output_pair.second.probability , buffer, len_buffer, position );
			pack_buff( output_pair.second.initialized , buffer, len_buffer, position );
		}

		intracellular->maboss.pack( buffer, len_buffer, position );
		pack_buff( intracellular->next_physiboss_run , buffer, len_buffer, position );
	}
#endif
}

void unpack_intracellular_light(Phenotype& phenotype, std::vector<char>& buffer, int& len_buffer, int& position)
{
	std::string intracellular_type;
	unpack_buff( intracellular_type , buffer, len_buffer, position );

	if( intracellular_type == "not-maboss" )
	{
		if( phenotype.intracellular != NULL )
		{
			delete phenotype.intracellular;
			phenotype.intracellular = NULL;
		}
		return;
	}

#ifdef ADDON_PHYSIBOSS
	if( intracellular_type == "maboss" )
	{
		if( phenotype.intracellular == NULL || phenotype.intracellular->intracellular_type != "maboss" )
		{
			delete phenotype.intracellular;
			phenotype.intracellular = new MaBoSSIntracellular();
		}

		MaBoSSIntracellular* intracellular = static_cast<MaBoSSIntracellular*>( phenotype.intracellular );

		int input_count = 0;
		unpack_buff( input_count , buffer, len_buffer, position );
		auto input_it = intracellular->listOfInputs.begin();
		for( int i = 0 ; i < input_count ; i++ )
		{
			double smoothed_value = 0.0;
			unpack_buff( smoothed_value , buffer, len_buffer, position );
			if( input_it != intracellular->listOfInputs.end() )
			{
				input_it->second.smoothed_value = smoothed_value;
				++input_it;
			}
		}

		int output_count = 0;
		unpack_buff( output_count , buffer, len_buffer, position );
		auto output_it = intracellular->listOfOutputs.begin();
		for( int i = 0 ; i < output_count ; i++ )
		{
			double probability = 0.0;
			bool initialized = false;
			unpack_buff( probability , buffer, len_buffer, position );
			unpack_buff( initialized , buffer, len_buffer, position );
			if( output_it != intracellular->listOfOutputs.end() )
			{
				output_it->second.probability = probability;
				output_it->second.initialized = initialized;
				++output_it;
			}
		}

		intracellular->initialize_maboss();
		intracellular->maboss.unpack( buffer, len_buffer, position );
		unpack_buff( intracellular->next_physiboss_run , buffer, len_buffer, position );
		return;
	}
#endif
}

void pack_phenotype_light(Cell* pCell, std::vector<char>& buffer, int& len_buffer, int& position)
{
	Phenotype& phenotype = pCell->phenotype;

	pack_buff( phenotype.flagged_for_division , buffer, len_buffer, position );
	pack_buff( phenotype.flagged_for_removal , buffer, len_buffer, position );

	pack_death_light( phenotype.death , buffer, len_buffer, position );
	pack_cycle_model_selector_light( pCell , buffer, len_buffer, position );
	pack_cycle_light( phenotype.cycle , buffer, len_buffer, position );
	phenotype.volume.pack( buffer, len_buffer, position );
	phenotype.geometry.pack( buffer, len_buffer, position );
	pack_motility_light( phenotype.motility , buffer, len_buffer, position );
	phenotype.secretion.pack( buffer, len_buffer, position );
	phenotype.molecular.pack( buffer, len_buffer, position );
	phenotype.cell_integrity.pack( buffer, len_buffer, position );
	pack_intracellular_light( phenotype , buffer, len_buffer, position );
	phenotype.cell_interactions.pack( buffer, len_buffer, position );
	phenotype.cell_transformations.pack( buffer, len_buffer, position );
}

void unpack_phenotype_light(Cell* pCell, Cell_Definition* pCD, std::vector<char>& buffer, int& len_buffer, int& position)
{
	Phenotype& phenotype = pCell->phenotype;

	unpack_buff( phenotype.flagged_for_division , buffer, len_buffer, position );
	unpack_buff( phenotype.flagged_for_removal , buffer, len_buffer, position );

	unpack_death_light( phenotype.death , buffer, len_buffer, position );
	Cycle_Model_Selector selector = unpack_cycle_model_selector_light( buffer, len_buffer, position );
	Cycle_Model* pModel = resolve_cycle_model_selector( pCD , selector );
	if( pModel != NULL )
	{
		phenotype.cycle.sync_to_cycle_model( *pModel );
	}
	unpack_cycle_light( phenotype.cycle , buffer, len_buffer, position );
	phenotype.volume.unpack( buffer, len_buffer, position );
	phenotype.geometry.unpack( buffer, len_buffer, position );
	unpack_motility_light( phenotype.motility , buffer, len_buffer, position );
	phenotype.secretion.unpack( buffer, len_buffer, position );
	phenotype.molecular.unpack( buffer, len_buffer, position );
	phenotype.cell_integrity.unpack( buffer, len_buffer, position );
	unpack_intracellular_light( phenotype , buffer, len_buffer, position );
	phenotype.cell_interactions.unpack( buffer, len_buffer, position );
	phenotype.cell_transformations.unpack( buffer, len_buffer, position );
}

} // namespace

//To consider: USE MPI_Pack_size() to get the size of the packed buffer
void Cell::pack(std::vector<char>& snd_buffer, int& len_buffer, int &position){

    //Packing class cell attributes
	int len_int;
	std::string temp_str;
	int temp_int;		
	temp_int = this->ID;
	pack_buff(temp_int, snd_buffer, len_buffer, position);
	pack_buff(this->position, snd_buffer, len_buffer, position);
	pack_buff(this->type_name, snd_buffer, len_buffer, position);

    //Custom_Cell_Data custom_data


    this->custom_data.pack(snd_buffer, len_buffer, position);

    this->parameters.pack(snd_buffer, len_buffer, position);

	this->functions.pack(snd_buffer, len_buffer, position);

    this->state.pack(snd_buffer, len_buffer, position);

	pack_cycle_model_selector_light( this , snd_buffer, len_buffer, position );
    this->phenotype.pack(snd_buffer, len_buffer, position);

    pack_buff(this->is_out_of_domain, snd_buffer, len_buffer, position);
    
    pack_buff(this->is_movable, snd_buffer, len_buffer, position);

    pack_buff(this->displacement, snd_buffer, len_buffer, position);

    /* Now Packing data members of class Basic_Agent - the parent class of class Cell */
    /* First 2 private data members are to be packed: double and bool 								*/

    pack_buff(this->get_total_volume(), snd_buffer, len_buffer, position);

    pack_buff(this->get_is_volume_changed(), snd_buffer, len_buffer, position);

    /* Now 3 std::vector<double> are to be packed - but are protected, hence write helper functions */


    std::vector<double> & temp_cell_source_sink_solver_temp1 = this->get_cell_source_sink_solver_temp1();
    pack_buff(temp_cell_source_sink_solver_temp1, snd_buffer, len_buffer, position);

    std::vector<double> & temp_cell_source_sink_solver_temp2 = this->get_cell_source_sink_solver_temp2();
    pack_buff(temp_cell_source_sink_solver_temp2, snd_buffer, len_buffer, position);

    /*=============================================================================================*/
    /* v1.7 introduces 2 more new protected vector of doubles - helper functions have been written */
    /* Now they will be packed - similar to the way the other 3 protected vector doubles 					 */
    /*=============================================================================================*/

    std::vector<double> & temp_cell_source_sink_solver_temp_export1 = this->get_cell_source_sink_solver_temp_export1();
    pack_buff(temp_cell_source_sink_solver_temp_export1, snd_buffer, len_buffer, position);

    std::vector<double> & temp_cell_source_sink_solver_temp_export2 = this->get_cell_source_sink_solver_temp_export2();
    pack_buff(temp_cell_source_sink_solver_temp_export2, snd_buffer, len_buffer, position);

    std::vector<double> & temp_previous_velocity = this->get_previous_velocity();
    pack_buff(temp_previous_velocity, snd_buffer, len_buffer, position);

    std::vector<double> & temp_total_extracellular_substrate_change = this->get_total_extracellular_substrate_change();
    pack_buff(temp_total_extracellular_substrate_change, snd_buffer, len_buffer, position);

    /* Now need to pack multiple std::vector<double> *ptr ; type pointers to vectors */

    pack_buff(*(this->secretion_rates), snd_buffer, len_buffer, position);

    pack_buff(*(this->saturation_densities), snd_buffer, len_buffer, position);

    pack_buff(*(this->uptake_rates), snd_buffer, len_buffer, position);
    /*==============================================================================================*/
    /* A new public data member of type pointer to a double vector named 'net_export_rates' in v1.7 */
    /* Now this will be packed 																																			*/
    /*==============================================================================================*/

    pack_buff(*(this->net_export_rates), snd_buffer, len_buffer, position);

    pack_buff(*(this->internalized_substrates), snd_buffer, len_buffer, position);

    pack_buff(*(this->fraction_released_at_death), snd_buffer, len_buffer, position);

    pack_buff(*(this->fraction_transferred_when_ingested), snd_buffer, len_buffer, position);
    /* 1 single int type; to be packed next, remember ID and position are already packed */

    pack_buff(this->type, snd_buffer, len_buffer, position);

    pack_buff(this->velocity, snd_buffer, len_buffer, position);
	//pCell->print_cell(world);
}

void Cell::unpack(std::vector<char>& rcv_buffer, int& len_buffer, int& position) {
    int len_int;
    std::string temp_str;
    int temp_int;

    unpack_buff(this->ID, rcv_buffer, len_buffer, position);

    // vector<double>
    unpack_buff(this->position, rcv_buffer, len_buffer, position);

    // string
    unpack_buff(this->type_name, rcv_buffer, len_buffer, position);

   this->custom_data.unpack(rcv_buffer, len_buffer, position);
   
   this->parameters.unpack(rcv_buffer, len_buffer, position);

   this->functions.unpack(rcv_buffer, len_buffer, position);

   this->state.unpack(rcv_buffer, len_buffer, position);

   Cycle_Model_Selector cycle_selector = unpack_cycle_model_selector_light( rcv_buffer, len_buffer, position );

   this->phenotype.unpack(rcv_buffer, len_buffer, position);
   
	unpack_buff(this->is_out_of_domain, rcv_buffer, len_buffer, position);

	unpack_buff(this->is_movable, rcv_buffer, len_buffer, position);

	unpack_buff(this->displacement, rcv_buffer, len_buffer, position);

	double total_volume;
	unpack_buff(total_volume, rcv_buffer, len_buffer, position);
	this->set_total_volume(total_volume);

	bool is_volume_changed;
	unpack_buff(is_volume_changed, rcv_buffer, len_buffer, position);
	this->set_is_volume_changed(is_volume_changed);

	std::vector<double> vec_double;
	unpack_buff(vec_double, rcv_buffer, len_buffer, position);
	this->set_cell_source_sink_solver_temp1(vec_double);

	unpack_buff(vec_double, rcv_buffer, len_buffer, position);
	this->set_cell_source_sink_solver_temp2(vec_double);

	unpack_buff(vec_double, rcv_buffer, len_buffer, position);
	this->set_cell_source_sink_solver_temp_export1(vec_double);

	unpack_buff(vec_double, rcv_buffer, len_buffer, position);
	this->set_cell_source_sink_solver_temp_export2(vec_double);

	unpack_buff(vec_double, rcv_buffer, len_buffer, position);
	//this->set_previous_velocity(vec_double);
	previous_velocity = vec_double;

	unpack_buff(vec_double, rcv_buffer, len_buffer, position);
	this->set_total_extracellular_substrate_change(vec_double);

	unpack_buff(*(this->secretion_rates), rcv_buffer, len_buffer, position);
	unpack_buff(*(this->saturation_densities), rcv_buffer, len_buffer, position);
	unpack_buff(*(this->uptake_rates), rcv_buffer, len_buffer, position);

	unpack_buff(*(this->net_export_rates), rcv_buffer, len_buffer, position);
	unpack_buff(*(this->internalized_substrates), rcv_buffer, len_buffer, position);
	unpack_buff(*(this->fraction_released_at_death), rcv_buffer, len_buffer, position);
	unpack_buff(*(this->fraction_transferred_when_ingested), rcv_buffer, len_buffer, position);
	unpack_buff(this->type, rcv_buffer, len_buffer, position);
	unpack_buff(this->velocity, rcv_buffer, len_buffer, position);

    Cell_Definition* pCD = find_cell_definition(this->type);
	Cycle_Model* pModel = resolve_cycle_model_selector( pCD , cycle_selector );
	sync_cycle_to_model_preserving_state( this , pModel );

} 

void Cell::pack_light(std::vector<char>& snd_buffer, int& len_buffer, int& position)
{
	pack_buff( this->ID , snd_buffer, len_buffer, position );
	pack_buff( this->type , snd_buffer, len_buffer, position );
	pack_buff( this->position , snd_buffer, len_buffer, position );
    pack_buff( this->is_movable, snd_buffer, len_buffer, position );

	pack_custom_data_values_light( this->custom_data , snd_buffer, len_buffer, position );
	this->state.pack( snd_buffer, len_buffer, position );
	pack_phenotype_light( this , snd_buffer, len_buffer, position );

	pack_buff( this->is_out_of_domain , snd_buffer, len_buffer, position );
	pack_buff( this->displacement , snd_buffer, len_buffer, position );
	pack_buff( this->get_previous_velocity() , snd_buffer, len_buffer, position );
	pack_buff( this->velocity , snd_buffer, len_buffer, position );
}

void Cell::unpack_light(std::vector<char>& rcv_buffer, int& len_buffer, int& position)
{
	unpack_buff( this->ID , rcv_buffer, len_buffer, position );
	unpack_buff( this->type , rcv_buffer, len_buffer, position );

	Cell_Definition* pCD = find_cell_definition( this->type );
    if (pCD != NULL)
	    restore_cell_from_definition( this , pCD );

	unpack_buff( this->position , rcv_buffer, len_buffer, position );
    unpack_buff( this->is_movable, rcv_buffer, len_buffer, position );

	unpack_custom_data_values_light( this->custom_data , rcv_buffer, len_buffer, position );
	this->state.unpack( rcv_buffer, len_buffer, position );
	unpack_phenotype_light( this , pCD , rcv_buffer, len_buffer, position );

	unpack_buff( this->is_out_of_domain , rcv_buffer, len_buffer, position );
	unpack_buff( this->displacement , rcv_buffer, len_buffer, position );
	unpack_buff( this->get_previous_velocity() , rcv_buffer, len_buffer, position );
	unpack_buff( this->velocity , rcv_buffer, len_buffer, position );

	apply_cycle_model_safety_net( this , pCD );
	this->set_total_volume( this->phenotype.volume.total );
	this->set_is_volume_changed( true );
}


//note: it does not print all parametes, but enought to know if pack/unpack are coherent
void Cell::print_parameters(std::ofstream& outFile) {
    outFile << std::fixed << std::setprecision(2); // for nicer floating-point printing

    // Basic variables
    outFile << "ID: " << this->ID << "\n";

    outFile << "Position: ";
    for (auto v : this->position) outFile << v << " ";
    outFile << "\n";

    outFile << "Type name: " << this->type_name << "\n";

    // Custom Data - name_to_index map
    /*
    
    outFile << "Custom Data - name_to_index:\n";
    for (const auto& kv : this->custom_data.get_name_to_index_map()) {
        outFile << "  [" << kv.first << "] = " << kv.second << "\n";
    }*/

    // Custom Data Variables (scalar)
    outFile << "Custom Data Variables:\n";
    for (const auto& var : this->custom_data.variables) {
        outFile << "  Name: " << var.name 
                  << ", Value: " << var.value 
                  << ", Units: " << var.units 
                  << ", Conserved: " << (var.conserved_quantity ? "Yes" : "No") << "\n";
    }

    // Custom Data Variables (vector)
    outFile << "Custom Data Vector Variables:\n";
    for (const auto& var : this->custom_data.vector_variables) {
        outFile << "  Name: " << var.name 
                  << ", Units: " << var.units 
                  << ", Conserved: " << (var.conserved_quantity ? "Yes" : "No") 
                  << ", Values: ";
        for (auto v : var.value) outFile << v << " ";
        outFile << "\n";
    }

    // Parameters
    outFile << "Parameters:\n"
              << "  o2_hypoxic_threshold: " << this->parameters.o2_hypoxic_threshold << "\n"
              << "  o2_hypoxic_response: " << this->parameters.o2_hypoxic_response << "\n"
              << "  o2_hypoxic_saturation: " << this->parameters.o2_hypoxic_saturation << "\n"
              << "  o2_proliferation_saturation: " << this->parameters.o2_proliferation_saturation << "\n"
              << "  o2_proliferation_threshold: " << this->parameters.o2_proliferation_threshold << "\n"
              << "  o2_reference: " << this->parameters.o2_reference << "\n"
              << "  o2_necrosis_threshold: " << this->parameters.o2_necrosis_threshold << "\n"
              << "  o2_necrosis_max: " << this->parameters.o2_necrosis_max << "\n"
              << "  max_necrosis_rate: " << this->parameters.max_necrosis_rate << "\n"
              << "  necrosis_type: " << this->parameters.necrosis_type << "\n";

    // Functions - cycle_model
    outFile << "Cycle Model:\n"
              << "  Name: " << this->functions.cycle_model.name << "\n"
              << "  Code: " << this->functions.cycle_model.code << "\n"
              << "  Default Phase Index: " << this->functions.cycle_model.default_phase_index << "\n";

    // Phases
    outFile << "  Phases (" << this->functions.cycle_model.phases.size() << "):\n";
    for (const auto& phase : this->functions.cycle_model.phases) {
        outFile << "    Index: " << phase.index 
                  << ", Name: " << phase.name 
                  << ", Division at phase exit: " << phase.division_at_phase_exit 
                  << ", Removal at phase exit: " << phase.removal_at_phase_exit << "\n";
    }

    // Phase links
    outFile << "  Phase Links (" << this->functions.cycle_model.phase_links.size() << "):\n";
    for (const auto& link : this->functions.cycle_model.phase_links) {
		for (const auto& phase : link)
        outFile << "    " << phase.start_phase_index << " -> " << phase.end_phase_index << "\n";
    }

    // Cycle data
    outFile << "Cycle Data:\n"
              << "  Time units: " << this->phenotype.cycle.data.time_units << "\n"
              << "  Transition rates (" << this->phenotype.cycle.data.transition_rates.size() << "): ";
    for (auto r : this->phenotype.cycle.data.transition_rates) outFile << r << " ";
    outFile << "\n  Current phase index: " << this->phenotype.cycle.data.current_phase_index << "\n"
              << "  Elapsed time in phase: " << this->phenotype.cycle.data.elapsed_time_in_phase << "\n"
			  << "  Asymetric division: ";
	for (const auto& p : phenotype.cycle.asymmetric_division.asymmetric_division_probabilities) {
		outFile << p << " ";
	}
	outFile << "\n";

    /*
    // Inverse index maps of cycle data
    outFile << "  Inverse index maps (" << this->phenotype.cycle.data.get_inverse_index_maps().size() << "):\n";
    for (size_t i = 0; i < this->phenotype.cycle.data.get_inverse_index_maps().size(); i++) {
        const auto& map_i = this->phenotype.cycle.data.get_inverse_index_maps()[i];
        outFile << "    Map " << i << " (" << map_i.size() << " entries):\n";
        for (const auto& kv : map_i) {
            outFile << "      " << kv.first << " => " << kv.second << "\n";
        }
    }*/

    // Death data
    outFile << "Death Data:\n"
              << "  Rates (" << this->phenotype.death.rates.size() << "): ";
    for (auto r : this->phenotype.death.rates) outFile << r << " ";
    outFile << "\n  Parameters (" << this->phenotype.death.parameters.size() << "):\n";

    for (const auto& p : this->phenotype.death.parameters) {
        outFile << "    Time units: " << p.time_units << "\n"
                  << "    Unlysed fluid change rate: " << p.unlysed_fluid_change_rate << "\n"
                  << "    Lysed fluid change rate: " << p.lysed_fluid_change_rate << "\n"
                  << "    Cytoplasmic biomass change rate: " << p.cytoplasmic_biomass_change_rate << "\n"
                  << "    Nuclear biomass change rate: " << p.nuclear_biomass_change_rate << "\n"
                  << "    Calcification rate: " << p.calcification_rate << "\n"
                  << "    Relative rupture volume: " << p.relative_rupture_volume << "\n";
    }

    outFile << "  Dead: " << this->phenotype.death.dead << "\n"
              << "  Current death model index: " << this->phenotype.death.current_death_model_index << "\n";

    // Volume
    outFile << "Volume:\n"
              << "  Total: " << this->phenotype.volume.total << "\n"
              << "  Solid: " << this->phenotype.volume.solid << "\n"
              << "  Fluid: " << this->phenotype.volume.fluid << "\n";

    outFile << "Geometry:\n"
              << "  Radius: " << this->phenotype.geometry.radius << "\n";

    // Mechanics
    outFile << "Mechanics:\n"
              << "  Cell-cell adhesion strength: " << this->phenotype.mechanics.cell_cell_adhesion_strength << "\n"
              << "  Cell adhesion affinities: " << this->phenotype.mechanics.cell_adhesion_affinities << "\n";

    // Motility
    outFile << "Motility:\n"
              << "  Is motile: " << (this->phenotype.motility.is_motile ? "Yes" : "No") << "\n"
              << "  Persistence time: " << this->phenotype.motility.persistence_time << "\n"
              << "  Migration speed: " << this->phenotype.motility.migration_speed << "\n"
              << "  Migration bias direction: " << this->phenotype.motility.migration_bias_direction << "\n"
              << "  Restrict to 2D: " << (this->phenotype.motility.restrict_to_2D ? "Yes" : "No") << "\n"
              << "  Motility vector: " << this->phenotype.motility.motility_vector << "\n";

    // Intracellular
    outFile << "Intracellular type: ";
    if (this->phenotype.intracellular != nullptr) {
        outFile << this->phenotype.intracellular->intracellular_type << "\n";
#ifdef ADDON_PHYSIBOSS
        if (this->phenotype.intracellular->intracellular_type == "maboss") {
            outFile << "  MaBoSS Intracellular object data (custom):\n";
            //this->phenotype.intracellular->print();  // Assuming it has a print method
        }
#endif
    } else {
        outFile << "not-maboss\n";
    }

    // Other booleans and vectors
    outFile << "is_out_of_domain: " << (this->is_out_of_domain ? "true" : "false") << "\n";
    outFile << "is_movable: " << (this->is_movable ? "true" : "false") << "\n";

    outFile << "Displacement: ";
    for (auto v : this->displacement) outFile << v << " ";
    outFile << "\n";

    outFile << "Total volume: " << this->get_total_volume() << "\n";
    outFile << "Is volume changed: " << this->get_is_volume_changed() << "\n";

    // Vectors from getters
    outFile << "cell_source_sink_solver_temp_export1: ";
    for (auto v : this->get_cell_source_sink_solver_temp_export1()) outFile << v << " ";
    outFile << "\n";

    outFile << "cell_source_sink_solver_temp_export2: ";
    for (auto v : this->get_cell_source_sink_solver_temp_export2()) outFile << v << " ";
    outFile << "\n";

    outFile << "previous_velocity: ";
    for (auto v : this->get_previous_velocity()) outFile << v << " ";
    outFile << "\n";

    outFile << "total_extracellular_substrate_change: ";
    for (auto v : this->get_total_extracellular_substrate_change()) outFile << v << " ";
    outFile << "\n";

    // Secretion, saturation, uptake
    outFile << "Secretion rates: ";
    for (auto v : *(this->secretion_rates)) outFile << v << " ";
    outFile << "\n";

    outFile << "Saturation densities: ";
    for (auto v : *(this->saturation_densities)) outFile << v << " ";
    outFile << "\n";

    outFile << "Uptake rates: ";
    for (auto v : *(this->uptake_rates)) outFile << v << " ";
    outFile << "\n";

    // Type and velocity
    outFile << "Type: " << this->type << "\n";

    outFile << "Velocity: ";
    for (auto v : this->velocity) outFile << v << " ";
    outFile << "\n";

    // Custom data variables again, just to be exhaustive (optional)
    // ...

    outFile << "------ End of Cell parameters ------\n";
}

void Cell::initialize_random() {
    // Seed random number generator
    std::srand(static_cast<unsigned>(std::time(nullptr)));

    // Initialize ID
    this->ID = std::rand();

    // Initialize position (3D vector)
    this->position = {static_cast<double>(std::rand()) / RAND_MAX, 
                      static_cast<double>(std::rand()) / RAND_MAX, 
                      static_cast<double>(std::rand()) / RAND_MAX};

    // Initialize type_name
    this->type_name = "RandomType_" + std::to_string(std::rand() % 100);

    // Initialize custom_data name_to_index_map
    this->custom_data.get_name_to_index_map().clear();
    int map_size = std::rand() % 5 + 1; // Random size between 1 and 5
    for (int i = 0; i < map_size; ++i) {
        this->custom_data.get_name_to_index_map()["Key_" + std::to_string(i)] = std::rand();
    }

    // Initialize custom_data variables
    this->custom_data.variables.clear();
    int var_size = std::rand() % 5 + 1;
    for (int i = 0; i < var_size; ++i) {
        Variable var;
        var.name = "Var_" + std::to_string(i);
        var.value = static_cast<double>(std::rand()) / RAND_MAX;
        var.units = "unit_" + std::to_string(i);
        var.conserved_quantity = std::rand() % 2;
        this->custom_data.variables.push_back(var);
    }

    // Initialize custom_data vector_variables
    this->custom_data.vector_variables.clear();
    int vec_var_size = std::rand() % 5 + 1;
    for (int i = 0; i < vec_var_size; ++i) {
        Vector_Variable vec_var;
        vec_var.name = "VecVar_" + std::to_string(i);
        vec_var.units = "unit_" + std::to_string(i);
        vec_var.conserved_quantity = std::rand() % 2;
        int vec_size = std::rand() % 3 + 1;
        for (int j = 0; j < vec_size; ++j) {
            vec_var.value.push_back(static_cast<double>(std::rand()) / RAND_MAX);
        }
        this->custom_data.vector_variables.push_back(vec_var);
    }

    // Initialize parameters
    this->parameters.o2_hypoxic_threshold = static_cast<double>(std::rand()) / RAND_MAX;
    this->parameters.o2_hypoxic_response = static_cast<double>(std::rand()) / RAND_MAX;
    this->parameters.o2_hypoxic_saturation = static_cast<double>(std::rand()) / RAND_MAX;
    this->parameters.o2_proliferation_saturation = static_cast<double>(std::rand()) / RAND_MAX;
    this->parameters.o2_proliferation_threshold = static_cast<double>(std::rand()) / RAND_MAX;
    this->parameters.o2_reference = static_cast<double>(std::rand()) / RAND_MAX;
    this->parameters.o2_necrosis_threshold = static_cast<double>(std::rand()) / RAND_MAX;
    this->parameters.o2_necrosis_max = static_cast<double>(std::rand()) / RAND_MAX;
    this->parameters.max_necrosis_rate = static_cast<double>(std::rand()) / RAND_MAX;
    this->parameters.necrosis_type = std::rand();

    // Initialize functions.cycle_model
    this->functions.cycle_model.name = "CycleModel_" + std::to_string(std::rand() % 100);
    this->functions.cycle_model.code = std::rand();
    this->functions.cycle_model.default_phase_index = std::rand() % 10;
    this->functions.cycle_model.get_inverse_index_maps().clear();
    map_size = std::rand() % 5 + 1; // Random size between 1 and 5
    this->functions.cycle_model.get_inverse_index_maps().resize(map_size);
    for (int i = 0; i < map_size; ++i) {
        this->functions.cycle_model.get_inverse_index_maps()[i][i] = std::rand();
    }

    // Initialize phases
    this->functions.cycle_model.phases.clear();
    int phase_count = std::rand() % 5 + 1;
    for (int i = 0; i < phase_count; ++i) {
        Phase phase;
        phase.index = i;
        phase.code = std::rand();
        phase.name = "Phase_" + std::to_string(i);
        phase.division_at_phase_exit = std::rand() % 2;
        phase.removal_at_phase_exit = std::rand() % 2;
        this->functions.cycle_model.phases.push_back(phase);
    }

    // Initialize phase_links
    this->functions.cycle_model.phase_links.clear();
    int link_count = std::rand() % 3 + 1;
    for (int i = 0; i < link_count; ++i) {
        std::vector<Phase_Link> links;
        int sub_link_count = std::rand() % 3 + 1;
        for (int j = 0; j < sub_link_count; ++j) {
            Phase_Link link;
            link.start_phase_index = std::rand() % 10;
            link.end_phase_index = std::rand() % 10;
            link.fixed_duration = static_cast<double>(std::rand()) / RAND_MAX;
            links.push_back(link);
        }
        this->functions.cycle_model.phase_links.push_back(links);
    }

	std::vector<std::unordered_map<int,int>> &inverse_index_maps_cycle_model = this->functions.cycle_model.get_inverse_index_maps();
    inverse_index_maps_cycle_model.clear();
    int map_count = std::rand() % 5 + 1; // Random number of maps between 1 and 5
    for (int i = 0; i < map_count; ++i) {
        std::unordered_map<int, int> map;
        int entry_count = std::rand() % 5 + 1; // Random number of entries per map
        for (int j = 0; j < entry_count; ++j) {
            int key = std::rand() % 100;
            int value = std::rand() % 100;
            map[key] = value;
        }
        inverse_index_maps_cycle_model.push_back(map);
    }

    this->functions.cycle_model.data.time_units = "seconds";
    this->functions.cycle_model.data.transition_rates.clear();
    int transition_rate_count = std::rand() % 5 + 1; // Random number of transition rates
    for (int i = 0; i < transition_rate_count; ++i) {
        std::vector<double> transition_rate;
        for (int j = 0; j < 3; ++j) { // Assuming 3 transition rates
            transition_rate.push_back(static_cast<double>(std::rand()) / RAND_MAX);
        }
        this->functions.cycle_model.data.transition_rates.push_back(transition_rate);
    }



    this->functions.cycle_model.data.current_phase_index = std::rand();
    this->functions.cycle_model.data.elapsed_time_in_phase = static_cast<double>(std::rand()) / RAND_MAX;


    this->state.orientation = {static_cast<double>(std::rand()) / RAND_MAX, 
        static_cast<double>(std::rand()) / RAND_MAX, 
        static_cast<double>(std::rand()) / RAND_MAX};
    
    this->state.simple_pressure = static_cast<double>(std::rand()) / RAND_MAX;
    this->state.number_of_nuclei = std::rand() % 5; // Random number of nuclei between 0 and 4
    this->state.total_attack_time = static_cast<double>(std::rand()) / RAND_MAX;
    this->state.contact_with_basement_membrane = std::rand() % 2; // Random boolean

   // Initialize phenotype
    this->phenotype.flagged_for_division = std::rand() % 2;
    this->phenotype.flagged_for_removal = std::rand() % 2;

    this->phenotype.cycle.data.get_inverse_index_maps().clear(); // Ensure it's empty before filling
    map_size = std::rand() % 5 + 1; // Random size between 1 and 5
    for (int i = 0; i < map_size; ++i) {
        std::unordered_map<int, int> map;
        int entry_count = std::rand() % 5 + 1; // Random number of entries per map
        for (int j = 0; j < entry_count; ++j) {
            int key = std::rand() % 100;
            int value = std::rand() % 100;
            map[key] = value;
        }
        this->phenotype.cycle.data.get_inverse_index_maps().push_back(map);
    }

    this->phenotype.cycle.data.time_units = "seconds";
    this->phenotype.cycle.data.transition_rates.clear();
    transition_rate_count = std::rand() % 5 + 1; // Random number of transition rates
    for (int i = 0; i < transition_rate_count; ++i) {
        std::vector<double> transition_rate;
        for (int j = 0; j < 3; ++j) { // Assuming 10 transition rates
            transition_rate.push_back(static_cast<double>(std::rand()) / RAND_MAX);
        }
        this->phenotype.cycle.data.transition_rates.push_back(transition_rate);
    }
    this->phenotype.cycle.data.current_phase_index = std::rand() % 10; // Random phase index
    this->phenotype.cycle.data.elapsed_time_in_phase = static_cast<double>(std::rand()) / RAND_MAX;
    this->phenotype.cycle.asymmetric_division.asymmetric_division_probabilities.clear();
    int asym_div_prob_count = std::rand() % 5 + 1; // Random number of probabilities
    for (int i = 0; i < asym_div_prob_count; ++i) {
        this->phenotype.cycle.asymmetric_division.asymmetric_division_probabilities.push_back(static_cast<double>(std::rand()) / RAND_MAX);
    }
    
    this->phenotype.death.rates.clear();
    int death_rate_count = std::rand() % 5 + 1; // Random number of death rates
    for (int i = 0; i < death_rate_count; ++i) {
        this->phenotype.death.rates.push_back(static_cast<double>(std::rand()) / RAND_MAX);
    }
    this->phenotype.death.parameters.clear();
    int death_param_count = std::rand() % 5 + 1; // Random number of death parameters
    for (int i = 0; i < death_param_count; ++i) {
        Death_Parameters param;
        param.time_units = "seconds";
        param.unlysed_fluid_change_rate = static_cast<double>(std::rand()) / RAND_MAX;
        param.lysed_fluid_change_rate = static_cast<double>(std::rand()) / RAND_MAX;
        param.cytoplasmic_biomass_change_rate = static_cast<double>(std::rand()) / RAND_MAX;
        param.nuclear_biomass_change_rate = static_cast<double>(std::rand()) / RAND_MAX;
        param.calcification_rate = static_cast<double>(std::rand()) / RAND_MAX;
        param.relative_rupture_volume = static_cast<double>(std::rand()) / RAND_MAX;
        this->phenotype.death.parameters.push_back(param);
    }

    this->phenotype.death.dead = std::rand() % 2; // Random boolean
    this->phenotype.death.current_death_model_index = std::rand() % 10; // Random index
    this->phenotype.volume.total = static_cast<double>(std::rand()) / RAND_MAX;
    this->phenotype.volume.solid = static_cast<double>(std::rand()) / RAND_MAX;
    this->phenotype.volume.fluid = static_cast<double>(std::rand()) / RAND_MAX;
    this->phenotype.volume.fluid_fraction = static_cast<double>(std::rand()) / RAND_MAX;
    this->phenotype.volume.nuclear = static_cast<double>(std::rand()) / RAND_MAX;
    this->phenotype.volume.nuclear_fluid = static_cast<double>(std::rand()) / RAND_MAX;
    this->phenotype.volume.nuclear_solid = static_cast<double>(std::rand()) / RAND_MAX;
    this->phenotype.volume.cytoplasmic = static_cast<double>(std::rand()) / RAND_MAX;
    this->phenotype.volume.cytoplasmic_fluid = static_cast<double>(std::rand()) / RAND_MAX;
    this->phenotype.volume.cytoplasmic_solid = static_cast<double>(std::rand()) / RAND_MAX;
    this->phenotype.volume.calcified_fraction = static_cast<double>(std::rand()) / RAND_MAX;
	this->phenotype.volume.cytoplasmic_to_nuclear_ratio= static_cast<double>(std::rand()) / RAND_MAX;

	this->phenotype.volume.rupture_volume= static_cast<double>(std::rand()) / RAND_MAX;
	this->phenotype.volume.cytoplasmic_biomass_change_rate= static_cast<double>(std::rand()) / RAND_MAX;
	this->phenotype.volume.nuclear_biomass_change_rate= static_cast<double>(std::rand()) / RAND_MAX;

	this->phenotype.volume.fluid_change_rate= static_cast<double>(std::rand()) / RAND_MAX;
	this->phenotype.volume.calcification_rate= static_cast<double>(std::rand()) / RAND_MAX;
	this->phenotype.volume.target_solid_cytoplasmic= static_cast<double>(std::rand()) / RAND_MAX;

	this->phenotype.volume.target_solid_nuclear= static_cast<double>(std::rand()) / RAND_MAX;
	this->phenotype.volume.target_fluid_fraction= static_cast<double>(std::rand()) / RAND_MAX;
	this->phenotype.volume.target_cytoplasmic_to_nuclear_ratio= static_cast<double>(std::rand()) / RAND_MAX;

	this->phenotype.volume.relative_rupture_volume= static_cast<double>(std::rand()) / RAND_MAX;


    this->phenotype.geometry.radius= static_cast<double>(std::rand()) / RAND_MAX;
	this->phenotype.geometry.nuclear_radius= static_cast<double>(std::rand()) / RAND_MAX;
	this->phenotype.geometry.surface_area= static_cast<double>(std::rand()) / RAND_MAX;
	this->phenotype.geometry.polarity= static_cast<double>(std::rand()) / RAND_MAX;

    this->phenotype.mechanics.cell_cell_adhesion_strength= static_cast<double>(std::rand()) / RAND_MAX;
	this->phenotype.mechanics.cell_BM_adhesion_strength= static_cast<double>(std::rand()) / RAND_MAX;
	this->phenotype.mechanics.cell_cell_repulsion_strength= static_cast<double>(std::rand()) / RAND_MAX;
	this->phenotype.mechanics.cell_BM_repulsion_strength= static_cast<double>(std::rand()) / RAND_MAX;
	this->phenotype.mechanics.relative_maximum_adhesion_distance= static_cast<double>(std::rand()) / RAND_MAX;

    this->phenotype.mechanics.attachment_rate= static_cast<double>(std::rand()) / RAND_MAX;
	this->phenotype.mechanics.detachment_rate= static_cast<double>(std::rand()) / RAND_MAX;
	
	this->phenotype.mechanics.relative_maximum_attachment_distance= static_cast<double>(std::rand()) / RAND_MAX;
	this->phenotype.mechanics.relative_detachment_distance= static_cast<double>(std::rand()) / RAND_MAX;
	this->phenotype.mechanics.attachment_elastic_constant= static_cast<double>(std::rand()) / RAND_MAX;
	this->phenotype.mechanics.maximum_attachment_rate= static_cast<double>(std::rand()) / RAND_MAX;
	
	this->phenotype.mechanics.maximum_number_of_attachments= std::rand();

    this->phenotype.mechanics.cell_adhesion_affinities.clear();
    int adhesion_affinity_count = std::rand() % 5 + 1; //
    for (int i = 0; i < adhesion_affinity_count; ++i) {
        this->phenotype.mechanics.cell_adhesion_affinities.push_back(static_cast<double>(std::rand()) / RAND_MAX);
    }

	this->phenotype.motility.is_motile = std::rand() % 2; // Random boolean
    this->phenotype.motility.persistence_time = static_cast<double>(std::rand()) / RAND_MAX;
    this->phenotype.motility.migration_speed = static_cast<double>(std::rand()) / RAND_MAX;


    this->phenotype.motility.migration_bias_direction = {static_cast<double>(std::rand()) / RAND_MAX, 
        static_cast<double>(std::rand()) / RAND_MAX, 
        static_cast<double>(std::rand()) / RAND_MAX};
    this->phenotype.motility.restrict_to_2D = std::rand() % 2; // Random boolean
    this->phenotype.motility.motility_vector = {static_cast<double>(std::rand()) / RAND_MAX, 
        static_cast<double>(std::rand()) / RAND_MAX, 
        static_cast<double>(std::rand()) / RAND_MAX};
    this->phenotype.motility.chemotaxis_index = std::rand() % 10; // Random index
    this->phenotype.motility.chemotaxis_direction = std::rand() % 10;
    this->phenotype.motility.chemotactic_sensitivities.clear();
    int chemotactic_sensitivity_count = std::rand() % 5 + 1; // Random number of sensitivities
    for (int i = 0; i < chemotactic_sensitivity_count; ++i) {
        this->phenotype.motility.chemotactic_sensitivities.push_back(static_cast<double>(std::rand()) / RAND_MAX);
    }

    /* Now packing data members of class Secretion as its object is a data member of class Phenotype */


	// Assume some random size for the vectors, e.g., 5 to 10 elements
    int substrate_count = std::rand() % 6 + 5; // random number between 5 and 10

    // Secretion rates
    this->phenotype.secretion.secretion_rates.clear();
    this->phenotype.secretion.secretion_rates.resize(substrate_count);
    for (int i = 0; i < substrate_count; ++i)
        this->phenotype.secretion.secretion_rates[i] = static_cast<double>(std::rand()) / RAND_MAX;

    // Uptake rates
    this->phenotype.secretion.uptake_rates.clear();
    this->phenotype.secretion.uptake_rates.resize(substrate_count);
    for (int i = 0; i < substrate_count; ++i)
        this->phenotype.secretion.uptake_rates[i] = static_cast<double>(std::rand()) / RAND_MAX;

    // Saturation densities
    this->phenotype.secretion.saturation_densities.clear();
    this->phenotype.secretion.saturation_densities.resize(substrate_count);
    for (int i = 0; i < substrate_count; ++i)
        this->phenotype.secretion.saturation_densities[i] = static_cast<double>(std::rand()) / RAND_MAX;

    // Net export rates
    this->phenotype.secretion.net_export_rates.clear();
    this->phenotype.secretion.net_export_rates.resize(substrate_count);
    for (int i = 0; i < substrate_count; ++i)
        this->phenotype.secretion.net_export_rates[i] = static_cast<double>(std::rand()) / RAND_MAX;

    // Internalized substrates
    this->phenotype.molecular.internalized_total_substrates.clear();
    this->phenotype.molecular.internalized_total_substrates.resize(substrate_count);
    for (int i = 0; i < substrate_count; ++i)
        this->phenotype.molecular.internalized_total_substrates[i] = static_cast<double>(std::rand()) / RAND_MAX;

    // Fraction released at death
    this->phenotype.molecular.fraction_released_at_death.clear();
    this->phenotype.molecular.fraction_released_at_death.resize(substrate_count);
    for (int i = 0; i < substrate_count; ++i)
        this->phenotype.molecular.fraction_released_at_death[i] = static_cast<double>(std::rand()) / RAND_MAX;

    // Fraction transferred when ingested
    this->phenotype.molecular.fraction_transferred_when_ingested.clear();
    this->phenotype.molecular.fraction_transferred_when_ingested.resize(substrate_count);
    for (int i = 0; i < substrate_count; ++i)
        this->phenotype.molecular.fraction_transferred_when_ingested[i] = static_cast<double>(std::rand()) / RAND_MAX;

    //this->phenotype.intracellular = new Intracellular();
    //this->phenotype.intracellular->intracellular_type = "Random";

    // Random boolean for is_out_of_domain (true or false)
    this->is_out_of_domain = std::rand() % 2;

    // Random boolean for is_movable
    this->is_movable = std::rand() % 2;

    // Random 3D vector for displacement
    this->displacement = {
        static_cast<double>(std::rand()) / RAND_MAX,
        static_cast<double>(std::rand()) / RAND_MAX,
        static_cast<double>(std::rand()) / RAND_MAX
    };

    this->set_total_volume(static_cast<double>(std::rand()) / RAND_MAX);
    this->set_is_volume_changed(std::rand() % 2);

    int temp_vector_size = std::rand() % 6 + 5; // Random size between 5 and 10

    // Initialize temp1
    std::vector<double>& temp_cell_source_sink_solver_temp1 = this->get_cell_source_sink_solver_temp1();
    temp_cell_source_sink_solver_temp1.clear();
    temp_cell_source_sink_solver_temp1.resize(temp_vector_size);
    for (int i = 0; i < temp_vector_size; ++i) {
        temp_cell_source_sink_solver_temp1[i] = static_cast<double>(std::rand()) / RAND_MAX;
    }

    // Initialize temp2
    std::vector<double>& temp_cell_source_sink_solver_temp2 = this->get_cell_source_sink_solver_temp2();
    temp_cell_source_sink_solver_temp2.clear();
    temp_cell_source_sink_solver_temp2.resize(temp_vector_size);
    for (int i = 0; i < temp_vector_size; ++i) {
        temp_cell_source_sink_solver_temp2[i] = static_cast<double>(std::rand()) / RAND_MAX;
    }

    // Assume substrate-related vectors have N elements
    substrate_count = std::rand() % 6 + 5; // Or replace with known substrate count

    // temp_export1
    std::vector<double>& temp_cell_source_sink_solver_temp_export1 = this->get_cell_source_sink_solver_temp_export1();
    temp_cell_source_sink_solver_temp_export1.clear();
    temp_cell_source_sink_solver_temp_export1.resize(substrate_count);
    for (int i = 0; i < substrate_count; ++i) {
        temp_cell_source_sink_solver_temp_export1[i] = static_cast<double>(std::rand()) / RAND_MAX;
    }

    // temp_export2
    std::vector<double>& temp_cell_source_sink_solver_temp_export2 = this->get_cell_source_sink_solver_temp_export2();
    temp_cell_source_sink_solver_temp_export2.clear();
    temp_cell_source_sink_solver_temp_export2.resize(substrate_count);
    for (int i = 0; i < substrate_count; ++i) {
        temp_cell_source_sink_solver_temp_export2[i] = static_cast<double>(std::rand()) / RAND_MAX;
    }

    // previous_velocity (3D vector)
    std::vector<double>& temp_previous_velocity = this->get_previous_velocity();
    temp_previous_velocity.clear();
    temp_previous_velocity.resize(3);
    for (int i = 0; i < 3; ++i) {
        temp_previous_velocity[i] = static_cast<double>(std::rand()) / RAND_MAX;
    }

    // total_extracellular_substrate_change
    std::vector<double>& temp_total_extracellular_substrate_change = this->get_total_extracellular_substrate_change();
    temp_total_extracellular_substrate_change.clear();
    temp_total_extracellular_substrate_change.resize(substrate_count);
    for (int i = 0; i < substrate_count; ++i) {
        temp_total_extracellular_substrate_change[i] = static_cast<double>(std::rand()) / RAND_MAX;
    }

    substrate_count = std::rand() % 6 + 5; // Or use a known substrate count

    // Initialize secretion_rates
    if (!this->secretion_rates) {
        this->secretion_rates = new std::vector<double>();
    }
    this->secretion_rates->clear();
    this->secretion_rates->resize(substrate_count);
    for (int i = 0; i < substrate_count; ++i) {
        (*(this->secretion_rates))[i] = static_cast<double>(std::rand()) / RAND_MAX;
    }

    // Initialize saturation_densities
    if (!this->saturation_densities) {
        this->saturation_densities = new std::vector<double>();
    }
    this->saturation_densities->clear();
    this->saturation_densities->resize(substrate_count);
    for (int i = 0; i < substrate_count; ++i) {
        (*(this->saturation_densities))[i] = static_cast<double>(std::rand()) / RAND_MAX;
    }

    // Initialize uptake_rates
    if (!this->uptake_rates) {
        this->uptake_rates = new std::vector<double>();
    }
    this->uptake_rates->clear();
    this->uptake_rates->resize(substrate_count);
    for (int i = 0; i < substrate_count; ++i) {
        (*(this->uptake_rates))[i] = static_cast<double>(std::rand()) / RAND_MAX;
    }

    substrate_count = std::rand() % 6 + 5; // Or a known constant

    // net_export_rates
    if (!this->net_export_rates)
        this->net_export_rates = new std::vector<double>();
    this->net_export_rates->clear();
    this->net_export_rates->resize(substrate_count);
    for (int i = 0; i < substrate_count; ++i)
        (*this->net_export_rates)[i] = static_cast<double>(std::rand()) / RAND_MAX;

    // internalized_substrates
    if (!this->internalized_substrates)
        this->internalized_substrates = new std::vector<double>();
    this->internalized_substrates->clear();
    this->internalized_substrates->resize(substrate_count);
    for (int i = 0; i < substrate_count; ++i)
        (*this->internalized_substrates)[i] = static_cast<double>(std::rand()) / RAND_MAX;

    // fraction_released_at_death
    if (!this->fraction_released_at_death)
        this->fraction_released_at_death = new std::vector<double>();
    this->fraction_released_at_death->clear();
    this->fraction_released_at_death->resize(substrate_count);
    for (int i = 0; i < substrate_count; ++i)
        (*this->fraction_released_at_death)[i] = static_cast<double>(std::rand()) / RAND_MAX;

    // fraction_transferred_when_ingested
    if (!this->fraction_transferred_when_ingested)
        this->fraction_transferred_when_ingested = new std::vector<double>();
    this->fraction_transferred_when_ingested->clear();
    this->fraction_transferred_when_ingested->resize(substrate_count);
    for (int i = 0; i < substrate_count; ++i)
        (*this->fraction_transferred_when_ingested)[i] = static_cast<double>(std::rand()) / RAND_MAX;

    // type (random int type identifier, 0–9)
    this->type = std::rand() % 10;

    // velocity (3D vector)
    this->velocity = {
        static_cast<double>(std::rand()) / RAND_MAX,
        static_cast<double>(std::rand()) / RAND_MAX,
        static_cast<double>(std::rand()) / RAND_MAX
    };

}
