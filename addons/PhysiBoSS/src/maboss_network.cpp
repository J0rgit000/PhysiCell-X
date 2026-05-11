#include "maboss_network.h"
#include <fstream>
#include "../../../core/MPI_helper.h"



/* Default constructor */
void MaBoSSNetwork::init_maboss( std::string networkFile, std::string configFile)
{
	if (this->engine != NULL) {
		delete this->engine;
		this->engine = NULL;
	}

	if (this->network != NULL) {
		delete this->network;
		this->network = NULL;
	}
	
	if (this->config != NULL) {
		delete this->config;
		this->config = NULL;
	}

	this->nodesByName.clear();
	this->parametersByName.clear();
	this->state = NetworkState();
	this->output_mask = NetworkState();
	
	try{
		
		#pragma omp critical
		{
			std::ifstream f_bnd(networkFile.c_str());
			if (!f_bnd.good()) {
				std::cerr << "PhysiBoSS ERROR : Could not open the BND file " << networkFile.c_str() << std::endl;
				exit(1);
			}
			
			std::ifstream f_cfg(configFile.c_str());
			if (!f_cfg.good()) {
				std::cerr << "PhysiBoSS ERROR : Could not open the CFG file " << configFile.c_str() << std::endl;
				exit(1);
			}
			
			// Initialize MaBoSS Objects for a model
			this->network = new Network();
			this->network->parse(networkFile.c_str());

			this->config = new RunConfig();
			this->config->parse(this->network, configFile.c_str());
		}
		
		// Some models will have chosen to use the physical randon number generator 
		// This is a problem, as it will open /dev/urandom for each cell, and overload the number of file open
		// So for now we just don't use this, and choose by default mersen twister
		this->config->setParameter("use_physrandgen", false);
		this->config->setParameter("use_mtrandgen", true);
		
		IStateGroup::checkAndComplete(this->network);

		engine = new StochasticSimulationEngine(this->network, this->config, PhysiCell::UniformInt());
	
	} catch (BNException e) {
		std::cerr << "MaBoSS ERROR : " << e.getMessage() << std::endl;
		exit(1);
	}
	this->update_time_step = this->config->getMaxTime();
	
	// Building map of nodes for fast later access 
	for (auto node : this->network->getNodes()) {
		this->nodesByName[node->getLabel()] = node;
	}
	
	// Building map of parameters for fast later access
	for (auto parameter : this->network->getSymbolTable()->getSymbolsNames()) {
		if (parameter[0] == '$')
			this->parametersByName[parameter] = this->network->getSymbolTable()->getSymbol(parameter);
	}
	
	for (auto node : network->getNodes())
      if (!node->isInternal()) 
        output_mask.setNodeState(node, true);

}

void MaBoSSNetwork::mutate(std::map<std::string, double> mutations) 
{
	for (auto mutation : mutations) {
		if (nodesByName.find(mutation.first) != nodesByName.end())
			nodesByName[mutation.first]->mutate(mutation.second);
		else{
			std::cerr << "Mutation set for unknown node : can't find node " << mutation.first << std::endl;
			exit(1);
		}
	}
}

void MaBoSSNetwork::set_parameters(std::map<std::string, double> parameters) 
{	
	for (auto parameter: parameters) {
		set_parameter_value(parameter.first, parameter.second);
	}
}

double MaBoSSNetwork::get_parameter_value(std::string name) 
{
	return network->getSymbolTable()->getSymbolValue(parametersByName[name]);
}


void MaBoSSNetwork::set_parameter_value(std::string name, double value) 
{
	network->getSymbolTable()->setSymbolValue(parametersByName[name], value);
	network->getSymbolTable()->unsetSymbolExpressions();
}

/* Reset a vector of bools to the init state of the network */
void MaBoSSNetwork::restart_node_values()
{
	// NetworkState network_state;
	this->network->initStates(state, engine->random_generator);
	
	for (auto initial_value : initial_values) {
		if (nodesByName.find(initial_value.first) != nodesByName.end()) {
			state.setNodeState(nodesByName[initial_value.first], PhysiCell::UniformRandom() < initial_value.second);	
		} else {
			std::cerr << "Initial value set for unknown node : can't find node " << initial_value.first << std::endl;
			exit(1);
		}
	}
	
	this->set_time_to_update();
}

/* Run a MaBoSS simulation with the input values*/
void MaBoSSNetwork::run_simulation()
{	
	engine->setMaxTime(time_to_update/scaling);
	state = engine->run(state, NULL);
	this->set_time_to_update();

}

bool MaBoSSNetwork::has_node( std::string name ) {
	return nodesByName.find(name) != nodesByName.end();
}

void MaBoSSNetwork::set_node_value(std::string name, bool value) {
	if (has_node(name))
		state.setNodeState(nodesByName[name], value);
	else 
		std::cout << "Can't find node " << name  << "!!!!" << std::endl;
}

bool MaBoSSNetwork::get_node_value(std::string name) {
	if (has_node(name))
		return state.getNodeState(nodesByName[name]);
	else
		std::cout << "Can't find node " << name  << "!!!!" << std::endl;
		return true;
}

std::string MaBoSSNetwork::get_state() {
	return NetworkState(state.getState() & output_mask.getState()).getName(network);
}

/* Print current state of all the nodes of the network */
void MaBoSSNetwork::print_nodes()
{
	int i = 0;
	std::vector<Node*> nodes = this->network->getNodes();
	for ( auto node: nodes )
	{
		std::cout << node->getLabel() << "=" << state.getNodeState(node) << "; ";
		i++;
	}
	std::cout << std::endl;
}


void MaBoSSNetwork::pack(std::vector<char>& buffer, int& len_buffer, int& position){

	// MaBoSS's time to update
	pack_buff(time_to_update, buffer, len_buffer, position);
	// this->maboss.set_time_to_update(time_to_update);
	
	// This is the part which might only need to be sent : 
	// The present state of the model
	// This will be improved by sending vectors of unsigned long long
	std::vector<Node*> t_nodes = network->getNodes();
	int net_size = t_nodes.size();
	pack_buff(net_size,buffer, len_buffer, position);
	int temp_int;
	for (unsigned int i=0; i < t_nodes.size(); i++) {
		temp_int = state.getNodeState(t_nodes[i]) == true ? 1 : 0;
		pack_buff(temp_int, buffer, len_buffer, position);
	}
	
	/* Now this is an additional Data Structure which is being packed */
	/* corresponding unpacking code is also added 										*/
	
	// The present parameter table
	SymbolTable* symbol_table = getNetwork()->getSymbolTable();
	temp_int = symbol_table->getSymbolCount();
	pack_buff(temp_int, buffer, len_buffer, position);

	double temp_double;
	for (unsigned int i=0; i < symbol_table->getSymbolCount(); i++) {
		temp_double = symbol_table->getSymbolValue(symbol_table->getSymbol(symbol_table->getSymbolsNames()[i]));
		pack_buff(temp_double, buffer, len_buffer, position);
	}
}

void MaBoSSNetwork::unpack(std::vector<char>& buffer, int& len_buffer, int& position){

	int temp_int;
	double temp_double;
	
	//unpack time to update
	MPI_Unpack(&buffer[0], len_buffer, &position, &(temp_double), 1, MPI_DOUBLE, MPI_COMM_WORLD);
	time_to_update=temp_double;

	//Unpack nodes state
	MPI_Unpack(&buffer[0], len_buffer, &position, &temp_int, 1, MPI_INT, MPI_COMM_WORLD);
	std::vector<Node*> t_nodes = getNetwork()->getNodes();
	if (t_nodes.size() != temp_int) 
		cout <<"Warning at: " << __FILE__ << ":" << __LINE__ << " in function " << __func__ << " nodes received in not matching the exisiting nodes" <<endl;
	for (unsigned int i=0; i < temp_int; i++) {
		int t_node = 0;
		MPI_Unpack(&buffer[0], len_buffer, &position, &t_node, 1, MPI_INT, MPI_COMM_WORLD);
		
		state.setNodeState(t_nodes[i], t_node == 1?true:false);
	}
	
	//unpack symbol table
	MPI_Unpack(&buffer[0], len_buffer, &position, &temp_int, 1, MPI_INT, MPI_COMM_WORLD);
	SymbolTable* symbol_table = getNetwork()->getSymbolTable();
	for (unsigned int i=0; i < temp_int; i++) {
		double t_parameter = 0;
		MPI_Unpack(&buffer[0], len_buffer, &position, &t_parameter, 1, MPI_DOUBLE, MPI_COMM_WORLD);
		symbol_table->setSymbolValue(symbol_table->getSymbol(symbol_table->getSymbolsNames()[i]), t_parameter);
	}
}

void MaBoSSNetwork::print_state_after_unpack() const {
    std::cout << "=== MaBoSSNetwork===" << std::endl;

    // 1. Print time_to_update
    std::cout << "time_to_update: " << time_to_update << std::endl;

    // 2. Print node states
    std::vector<Node*> t_nodes = network->getNodes();
    std::cout << "Node states (" << t_nodes.size() << " nodes):" << std::endl;
    for (size_t i = 0; i < t_nodes.size(); ++i) {
        bool state_val = state.getNodeState(t_nodes[i]);
        std::cout << "  Node " << i << " (" << t_nodes[i]->getLabel() << "): "
                  << (state_val ? "ON" : "OFF") << std::endl;
    }

    // 3. Print symbol table values
    SymbolTable* symbol_table = network->getSymbolTable();
    const auto& symbol_names = symbol_table->getSymbolsNames();
    std::cout << "Symbol table (" << symbol_names.size() << " symbols):" << std::endl;
    for (const auto& name : symbol_names) {
        double value = symbol_table->getSymbolValue(symbol_table->getSymbol(name));
        std::cout << "  " << name << " = " << value << std::endl;
    }

    std::cout << "========================================" << std::endl;
}
