#include "./custom.h"

#include <cctype>
#include <fstream>
#include <sstream>
#include <stdexcept>

#include "../../../../modules/PhysiCell_pugixml.h"
#include "../../../../modules/PhysiCell_settings.h"

namespace
{

struct SeedCell
{
	int cell_id = -1;
	int type_id = 0;
	double x = 0.0;
	double y = 0.0;
	double z = 0.0;
};

std::string trim_copy( const std::string& input )
{
	size_t start = 0;
	while( start < input.size() && std::isspace( static_cast<unsigned char>( input[start] ) ) )
	{
		++start;
	}

	size_t end = input.size();
	while( end > start && std::isspace( static_cast<unsigned char>( input[end - 1] ) ) )
	{
		--end;
	}

	return input.substr( start, end - start );
}

std::vector<std::string> split_csv_line( const std::string& line )
{
	std::vector<std::string> tokens;
	std::stringstream stream( line );
	std::string token;

	while( std::getline( stream, token, ',' ) )
	{
		tokens.push_back( trim_copy( token ) );
	}

	return tokens;
}

std::string resolve_seed_cells_path( void )
{
	pugi::xml_node node = physicell_config_root.child( "initial_conditions" ).child( "cell_positions" );
	if( !node )
	{
		throw std::runtime_error( "Missing initial_conditions/cell_positions in mechanics validation config." );
	}

	if( node.attribute( "enabled" ).empty() || node.attribute( "enabled" ).as_bool() == false )
	{
		throw std::runtime_error( "Mechanics validation requires initial_conditions/cell_positions enabled=true." );
	}

	const std::string folder = xml_get_string_value( node, "folder" );
	const std::string filename = xml_get_string_value( node, "filename" );

	if( folder.empty() )
	{
		return filename;
	}

	return folder + "/" + filename;
}

std::vector<SeedCell> read_seed_cells_from_config( void )
{
	const std::string path = resolve_seed_cells_path();
	std::ifstream stream( path );
	if( !stream )
	{
		throw std::runtime_error( "Unable to open mechanics validation cell layout: " + path );
	}

	std::vector<SeedCell> cells;
	std::string line;
	while( std::getline( stream, line ) )
	{
		const std::string trimmed = trim_copy( line );
		if( trimmed.empty() || trimmed[0] == '#' )
		{
			continue;
		}

		const std::vector<std::string> fields = split_csv_line( trimmed );
		if( fields.size() != 5 )
		{
			throw std::runtime_error(
				"Mechanics validation CSV must have 5 columns: cell_id,x,y,z,type_id. Bad line: " + trimmed );
		}

		if( !fields[0].empty() && !std::isdigit( static_cast<unsigned char>( fields[0][0] ) )
			&& fields[0][0] != '-' )
		{
			continue;
		}

		SeedCell cell;
		cell.cell_id = std::stoi( fields[0] );
		cell.x = std::stod( fields[1] );
		cell.y = std::stod( fields[2] );
		cell.z = std::stod( fields[3] );
		cell.type_id = std::stoi( fields[4] );
		cells.push_back( cell );
	}

	if( cells.empty() )
	{
		throw std::runtime_error( "Mechanics validation CSV does not define any cells." );
	}

	return cells;
}

void make_definition_mechanics_only( Cell_Definition& definition, bool parallel )
{
	definition.functions.volume_update_function = standard_volume_update_function;
	definition.functions.update_velocity = standard_update_cell_velocity;
	if( parallel )
	{
		definition.functions.update_velocity_parallel = standard_update_cell_velocity;
	}
	definition.functions.update_migration_bias = NULL;
	definition.functions.update_phenotype = NULL;
	definition.functions.custom_cell_rule = NULL;
	definition.functions.contact_function = NULL;
	definition.functions.add_cell_basement_membrane_interactions = NULL;
	definition.functions.calculate_distance_to_membrane = NULL;

	for( size_t i = 0; i < definition.phenotype.cycle.data.transition_rates.size(); ++i )
	{
		for( size_t j = 0; j < definition.phenotype.cycle.data.transition_rates[i].size(); ++j )
		{
			definition.phenotype.cycle.data.transition_rates[i][j] = 0.0;
		}
	}

	for( size_t i = 0; i < definition.phenotype.death.rates.size(); ++i )
	{
		definition.phenotype.death.rates[i] = 0.0;
	}
}

void initialize_mechanics_cell_types( bool parallel, DistPhy::mpi::mpi_Environment* world = NULL,
	DistPhy::mpi::mpi_Cartesian* cart_topo = NULL )
{
	if( parameters.ints.find_index( "random_seed" ) != -1 )
	{
		SeedRandom( parameters.ints( "random_seed" ) );
	}

	initialize_default_cell_definition();
	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment );
	make_definition_mechanics_only( cell_defaults, parallel );

	if( parallel )
	{
		initialize_cell_definitions_from_pugixml( *world, *cart_topo );
		build_cell_definitions_maps();
		setup_signal_behavior_dictionaries( *world, *cart_topo );
	}
	else
	{
		DistPhy::mpi::mpi_Environment serial_world;
		DistPhy::mpi::mpi_Cartesian serial_cart_topo;
		serial_world.size = 1;
		serial_world.rank = 0;
		serial_world.init_comm = MPI_COMM_NULL;
		initialize_cell_definitions_from_pugixml( physicell_config_root, serial_world, serial_cart_topo );
		build_cell_definitions_maps();
		setup_signal_behavior_dictionaries( serial_world, serial_cart_topo );
	}

	build_cell_definitions_maps();

	Cell_Definition* definition = find_cell_definition( "cancer cell" );
	if( definition == NULL )
	{
		throw std::runtime_error( "Mechanics validation config must define a cell_definition named 'cancer cell'." );
	}

	make_definition_mechanics_only( *definition, parallel );
	build_cell_definitions_maps();
}

void create_seed_cell( const SeedCell& seed )
{
	Cell_Definition* definition = find_cell_definition( seed.type_id );
	if( definition == NULL )
	{
		throw std::runtime_error( "Unable to map mechanics validation cell type ID " + std::to_string( seed.type_id ) );
	}

	Cell* cell = create_cell( *definition, seed.cell_id );
	cell->assign_position( seed.x, seed.y, seed.z );
}

void create_seed_cell( const SeedCell& seed, DistPhy::mpi::mpi_Environment& world,
	DistPhy::mpi::mpi_Cartesian& cart_topo )
{
	Cell_Definition* definition = find_cell_definition( seed.type_id );
	if( definition == NULL )
	{
		throw std::runtime_error( "Unable to map mechanics validation cell type ID " + std::to_string( seed.type_id ) );
	}

	Cell* cell = create_cell( *definition, seed.cell_id );
	cell->assign_position( seed.x, seed.y, seed.z, world, cart_topo );
}

bool is_in_local_x_range( double x, const Microenvironment& microenvironment,
	DistPhy::mpi::mpi_Environment& world )
{
	const double local_x_min = microenvironment.mesh.local_bounding_box[0];
	const double local_x_max = microenvironment.mesh.local_bounding_box[3];
	return x >= local_x_min && ( world.rank == world.size - 1 ? x <= local_x_max : x < local_x_max );
}

}

void create_cell_types( void )
{
	initialize_mechanics_cell_types( false );
}

void create_cell_types( DistPhy::mpi::mpi_Environment& world, DistPhy::mpi::mpi_Cartesian& cart_topo )
{
	initialize_mechanics_cell_types( true, &world, &cart_topo );
}

void setup_microenvironment( void )
{
	if( default_microenvironment_options.simulate_2D == true )
	{
		std::cout << "Warning: overriding XML config option and setting to 3D!" << std::endl;
		default_microenvironment_options.simulate_2D = false;
	}

	initialize_microenvironment();
}

void setup_microenvironment( DistPhy::mpi::mpi_Environment& world, DistPhy::mpi::mpi_Cartesian& cart_topo )
{
	if( default_microenvironment_options.simulate_2D == true )
	{
		if( IOProcessor( world ) )
		{
			std::cout << "Warning: overriding XML config option and setting to 3D!" << std::endl;
		}
		default_microenvironment_options.simulate_2D = false;
	}

	initialize_microenvironment( world, cart_topo );
}

void setup_tissue( void )
{
	const std::vector<SeedCell> cells = read_seed_cells_from_config();
	for( const SeedCell& cell : cells )
	{
		create_seed_cell( cell );
	}
}

void setup_tissue( Microenvironment& microenvironment, DistPhy::mpi::mpi_Environment& world,
	DistPhy::mpi::mpi_Cartesian& cart_topo )
{
	const std::vector<SeedCell> cells = read_seed_cells_from_config();
	for( const SeedCell& cell : cells )
	{
		if( is_in_local_x_range( cell.x, microenvironment, world ) )
		{
			create_seed_cell( cell, world, cart_topo );
		}
	}
}
