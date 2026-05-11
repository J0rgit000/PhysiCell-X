#ifndef MPI_HELPER_H
#define MPI_HELPER_H

#include<mpi.h>
#include<vector>
#include<string>
#include<map>
using namespace std;

inline int mpi_packed_size(int count, MPI_Datatype datatype)
{
	int size = 0;
	MPI_Pack_size(count, datatype, MPI_COMM_WORLD, &size);
	return size;
}

inline void resize_for_pack(vector<char>& buffer, int& len_buffer, int position, int packed_size)
{
	len_buffer = position + packed_size;
	buffer.resize(len_buffer);
}

inline void finalize_pack_size(vector<char>& buffer, int& len_buffer, int position)
{
	len_buffer = position;
	buffer.resize(len_buffer);
}

//packing functions to avoid code replication and enhance clarity
//pack bool
inline void pack_buff(const bool& input, vector<char>& buffer, int& len_buffer, int& position){
	int temp_int = (input == true) ? 1 : 0;
	resize_for_pack(buffer, len_buffer, position, mpi_packed_size(1, MPI_INT));
	MPI_Pack(&(temp_int), 1, MPI_INT, &buffer[0], len_buffer, &position, MPI_COMM_WORLD); 
	finalize_pack_size(buffer, len_buffer, position);
}

//pack double
inline void pack_buff(const double& input, vector<char>& buffer, int& len_buffer, int& position){
	resize_for_pack(buffer, len_buffer, position, mpi_packed_size(1, MPI_DOUBLE));
	MPI_Pack(&(input), 1, MPI_DOUBLE, &buffer[0], len_buffer, &position, MPI_COMM_WORLD);
	finalize_pack_size(buffer, len_buffer, position);
}

//pack int
inline void pack_buff(const int& input, vector<char>& buffer, int& len_buffer, int& position){
	resize_for_pack(buffer, len_buffer, position, mpi_packed_size(1, MPI_INT));
	MPI_Pack(&(input), 1, MPI_INT, &buffer[0], len_buffer, &position, MPI_COMM_WORLD);
	finalize_pack_size(buffer, len_buffer, position);
}

//pack string
inline void pack_buff(const string& input, vector<char>& buffer, int& len_buffer, int& position){
	int len_str = static_cast<int>(input.length());
	int packed_size = mpi_packed_size(1, MPI_INT) + mpi_packed_size(len_str, MPI_CHAR);
	resize_for_pack(buffer, len_buffer, position, packed_size);
	MPI_Pack(&len_str, 1, MPI_INT, &buffer[0], len_buffer, &position, MPI_COMM_WORLD); 
	if (len_str > 0) {
		MPI_Pack(input.data(), len_str, MPI_CHAR, &buffer[0], len_buffer, &position, MPI_COMM_WORLD);
	}
	finalize_pack_size(buffer, len_buffer, position);
}

//pack vector<int>
inline void pack_buff(const vector<int>& input, vector<char>& buffer, int& len_buffer, int& position){
	int len_vector = static_cast<int>(input.size());
	int packed_size = mpi_packed_size(1, MPI_INT) + mpi_packed_size(len_vector, MPI_INT);
	resize_for_pack(buffer, len_buffer, position, packed_size);
	MPI_Pack(&(len_vector), 1, MPI_INT, &buffer[0], len_buffer, &position, MPI_COMM_WORLD);
	if (len_vector > 0) {
		MPI_Pack(input.data(), len_vector, MPI_INT, &buffer[0], len_buffer, &position, MPI_COMM_WORLD);
	}
	finalize_pack_size(buffer, len_buffer, position);
}

//pack vector<double>
inline void pack_buff(const vector<double>& input, vector<char>& buffer, int& len_buffer, int& position){
	int len_vector = static_cast<int>(input.size());
	int packed_size = mpi_packed_size(1, MPI_INT) + mpi_packed_size(len_vector, MPI_DOUBLE);
	resize_for_pack(buffer, len_buffer, position, packed_size);
	MPI_Pack(&(len_vector), 1, MPI_INT, &buffer[0], len_buffer, &position, MPI_COMM_WORLD);
	if (len_vector > 0) {
		MPI_Pack(input.data(), len_vector, MPI_DOUBLE, &buffer[0], len_buffer, &position, MPI_COMM_WORLD);
	}
	finalize_pack_size(buffer, len_buffer, position);
}

//Unpacking functions
//unpack bool
inline void unpack_buff(bool& output, const vector<char>& buffer, int len_buffer, int& position) {
    int temp_int;
    MPI_Unpack(buffer.data(), len_buffer, &position, &temp_int, 1, MPI_INT, MPI_COMM_WORLD);
    output = (temp_int != 0);
}

inline void unpack_buff(int& output, const vector<char>& buffer, int len_buffer, int& position) {
    MPI_Unpack(buffer.data(), len_buffer, &position, &output, 1, MPI_INT, MPI_COMM_WORLD);
}

//unpack double
inline void unpack_buff(double& output, const vector<char>& buffer, int len_buffer, int& position) {
    MPI_Unpack(buffer.data(), len_buffer, &position, &output, 1, MPI_DOUBLE, MPI_COMM_WORLD);
}

//unpack string
inline void unpack_buff(string& output, const vector<char>& buffer, int len_buffer, int& position) {
	int len_str = 0;
	MPI_Unpack(buffer.data(), len_buffer, &position, &len_str, 1, MPI_INT, MPI_COMM_WORLD);
	output.clear();
	if (len_str > 0) {
		output.resize(len_str);
		MPI_Unpack(buffer.data(), len_buffer, &position, &output[0], len_str, MPI_CHAR, MPI_COMM_WORLD);
	}
}

//unpack vector<int>
inline void unpack_buff(vector<int>& output, const vector<char>& buffer, int len_buffer, int& position) {
	int len_vector = 0;
	MPI_Unpack(buffer.data(), len_buffer, &position, &len_vector, 1, MPI_INT, MPI_COMM_WORLD);
	output.resize(len_vector);
	if (len_vector > 0) {
		MPI_Unpack(buffer.data(), len_buffer, &position, output.data(), len_vector, MPI_INT, MPI_COMM_WORLD);
	}
}

// unpack vector<double>
inline void unpack_buff(std::vector<double>& output, const std::vector<char>& buffer,int len_buffer, int& position) {
	int len_vector = 0;

	// Unpack the size of the vector
	MPI_Unpack(buffer.data(), len_buffer, &position, &len_vector, 1, MPI_INT, MPI_COMM_WORLD);

	output.resize(len_vector);
	if (len_vector > 0) {
		// Unpack the vector's data
		MPI_Unpack(buffer.data(), len_buffer, &position, output.data(), len_vector, MPI_DOUBLE, MPI_COMM_WORLD);
	}
}
#endif // MPI_HELPER_H
