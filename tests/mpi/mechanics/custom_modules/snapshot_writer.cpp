#include "./snapshot_writer.h"

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <vector>

#include <mpi.h>

using namespace PhysiCell;

namespace
{

struct SnapshotRecord
{
	int cell_id = -1;
	double x = 0.0;
	double y = 0.0;
	double z = 0.0;
	double vx = 0.0;
	double vy = 0.0;
	double vz = 0.0;
	double radius = 0.0;
	double relative_maximum_adhesion_distance = 0.0;
};

std::vector<SnapshotRecord> collect_local_snapshot_records( void )
{
	std::vector<SnapshotRecord> records;
	if( all_cells == NULL )
	{
		return records;
	}

	records.reserve( all_cells->size() );
	for( Cell* cell : *all_cells )
	{
		if( cell == NULL || cell->phenotype.death.dead )
		{
			continue;
		}

		SnapshotRecord record;
		record.cell_id = cell->ID;
		record.x = cell->position[0];
		record.y = cell->position[1];
		record.z = cell->position[2];
		record.vx = cell->velocity[0];
		record.vy = cell->velocity[1];
		record.vz = cell->velocity[2];
		record.radius = cell->phenotype.geometry.radius;
		record.relative_maximum_adhesion_distance = cell->phenotype.mechanics.relative_maximum_adhesion_distance;
		records.push_back( record );
	}

	return records;
}

std::vector<SnapshotRecord> collect_owned_mpi_snapshot_records( DistPhy::mpi::mpi_Environment& world )
{
	std::vector<SnapshotRecord> records;
	if( all_cells == NULL )
	{
		return records;
	}

	const double local_x_min = microenvironment.mesh.local_bounding_box[0];
	const double local_x_max = microenvironment.mesh.local_bounding_box[3];

	records.reserve( all_cells->size() );
	for( Cell* cell : *all_cells )
	{
		if( cell == NULL || cell->phenotype.death.dead )
		{
			continue;
		}

		const double x = cell->position[0];
		const bool in_local_x_range =
			x >= local_x_min &&
			( world.rank == world.size - 1 ? x <= local_x_max : x < local_x_max );
		if( !in_local_x_range )
		{
			continue;
		}

		SnapshotRecord record;
		record.cell_id = cell->ID;
		record.x = cell->position[0];
		record.y = cell->position[1];
		record.z = cell->position[2];
		record.vx = cell->velocity[0];
		record.vy = cell->velocity[1];
		record.vz = cell->velocity[2];
		record.radius = cell->phenotype.geometry.radius;
		record.relative_maximum_adhesion_distance = cell->phenotype.mechanics.relative_maximum_adhesion_distance;
		records.push_back( record );
	}

	return records;
}

std::string serialize_records( const std::vector<SnapshotRecord>& records )
{
	std::ostringstream stream;
	stream << std::setprecision( 17 );
	for( const SnapshotRecord& record : records )
	{
		stream << record.cell_id << ","
			<< record.x << ","
			<< record.y << ","
			<< record.z << ","
			<< record.vx << ","
			<< record.vy << ","
			<< record.vz << ","
			<< record.radius << ","
			<< record.relative_maximum_adhesion_distance << "\n";
	}
	return stream.str();
}

std::vector<SnapshotRecord> parse_records( const std::string& payload )
{
	std::vector<SnapshotRecord> records;
	std::stringstream lines( payload );
	std::string line;
	while( std::getline( lines, line ) )
	{
		if( line.empty() )
		{
			continue;
		}

		std::stringstream row( line );
		std::string token;
		std::vector<std::string> fields;
		while( std::getline( row, token, ',' ) )
		{
			fields.push_back( token );
		}

		if( fields.size() != 9 )
		{
			throw std::runtime_error( "Malformed snapshot payload row: " + line );
		}

		SnapshotRecord record;
		record.cell_id = std::stoi( fields[0] );
		record.x = std::stod( fields[1] );
		record.y = std::stod( fields[2] );
		record.z = std::stod( fields[3] );
		record.vx = std::stod( fields[4] );
		record.vy = std::stod( fields[5] );
		record.vz = std::stod( fields[6] );
		record.radius = std::stod( fields[7] );
		record.relative_maximum_adhesion_distance = std::stod( fields[8] );
		records.push_back( record );
	}

	return records;
}

std::map<int, SnapshotRecord> deduplicate_by_cell_id( const std::vector<SnapshotRecord>& records )
{
	std::map<int, SnapshotRecord> unique_records;
	for( const SnapshotRecord& record : records )
	{
		unique_records.emplace( record.cell_id, record );
	}
	return unique_records;
}

bool records_are_neighbors( const SnapshotRecord& left, const SnapshotRecord& right )
{
	const double dx = left.x - right.x;
	const double dy = left.y - right.y;
	const double dz = left.z - right.z;
	const double distance = std::sqrt( dx * dx + dy * dy + dz * dz );
	const double threshold =
		left.radius * left.relative_maximum_adhesion_distance
		+ right.radius * right.relative_maximum_adhesion_distance;
	return distance <= threshold;
}

std::map<int, std::set<int>> build_neighbor_sets( const std::map<int, SnapshotRecord>& records )
{
	std::map<int, std::set<int>> neighbors;
	for( const auto& entry : records )
	{
		neighbors[entry.first];
	}

	for( auto left = records.begin(); left != records.end(); ++left )
	{
		auto right = left;
		++right;
		for( ; right != records.end(); ++right )
		{
			if( records_are_neighbors( left->second, right->second ) )
			{
				neighbors[left->first].insert( right->first );
				neighbors[right->first].insert( left->first );
			}
		}
	}

	return neighbors;
}

std::string join_neighbor_ids( const std::set<int>& neighbors )
{
	std::ostringstream stream;
	bool first = true;
	for( int neighbor_id : neighbors )
	{
		if( !first )
		{
			stream << ";";
		}
		stream << neighbor_id;
		first = false;
	}
	return stream.str();
}

void append_snapshot_rows( const std::string& path, int iteration, double time,
	const std::map<int, SnapshotRecord>& records )
{
	const std::map<int, std::set<int>> neighbors = build_neighbor_sets( records );

	std::ofstream stream( path, std::ios::app );
	if( !stream )
	{
		throw std::runtime_error( "Unable to append mechanics snapshot file: " + path );
	}

	stream << std::setprecision( 17 );
	for( const auto& entry : records )
	{
		const SnapshotRecord& record = entry.second;
		stream << iteration << ","
			<< time << ","
			<< record.cell_id << ","
			<< record.x << ","
			<< record.y << ","
			<< record.z << ","
			<< record.vx << ","
			<< record.vy << ","
			<< record.vz << ","
			<< join_neighbor_ids( neighbors.at( record.cell_id ) ) << "\n";
	}
}

}

void initialize_snapshot_csv( const std::string& path )
{
	std::ofstream stream( path );
	if( !stream )
	{
		throw std::runtime_error( "Unable to create mechanics snapshot file: " + path );
	}

	stream << "iteration,time,cell_id,x,y,z,vx,vy,vz,neighbor_ids\n";
}

void append_serial_snapshot( const std::string& path, int iteration, double time )
{
	append_snapshot_rows( path, iteration, time, deduplicate_by_cell_id( collect_local_snapshot_records() ) );
}

void append_mpi_snapshot( const std::string& path, int iteration, double time,
	DistPhy::mpi::mpi_Environment& world, DistPhy::mpi::mpi_Cartesian& cart_topo )
{
	const std::string local_payload = serialize_records( collect_owned_mpi_snapshot_records( world ) );
	const int local_size = static_cast<int>( local_payload.size() );

	std::vector<int> recv_counts;
	if( world.rank == 0 )
	{
		recv_counts.resize( world.size, 0 );
	}

	MPI_Gather( &local_size, 1, MPI_INT, recv_counts.data(), 1, MPI_INT, 0, cart_topo.mpi_cart_comm );

	std::vector<int> displacements;
	std::vector<char> recv_buffer;
	if( world.rank == 0 )
	{
		displacements.resize( world.size, 0 );
		int total_size = 0;
		for( int i = 0; i < world.size; ++i )
		{
			displacements[i] = total_size;
			total_size += recv_counts[i];
		}
		recv_buffer.resize( total_size );
	}

	MPI_Gatherv(
		local_payload.empty() ? NULL : const_cast<char*>( local_payload.data() ),
		local_size,
		MPI_CHAR,
		recv_buffer.empty() ? NULL : recv_buffer.data(),
		recv_counts.empty() ? NULL : recv_counts.data(),
		displacements.empty() ? NULL : displacements.data(),
		MPI_CHAR,
		0,
		cart_topo.mpi_cart_comm );

	if( world.rank != 0 )
	{
		return;
	}

	std::vector<SnapshotRecord> gathered_records;
	for( int i = 0; i < world.size; ++i )
	{
		if( recv_counts[i] == 0 )
		{
			continue;
		}

		const std::string payload(
			recv_buffer.data() + displacements[i],
			static_cast<size_t>( recv_counts[i] ) );
		const std::vector<SnapshotRecord> rank_records = parse_records( payload );
		gathered_records.insert( gathered_records.end(), rank_records.begin(), rank_records.end() );
	}

	append_snapshot_rows( path, iteration, time, deduplicate_by_cell_id( gathered_records ) );
}
