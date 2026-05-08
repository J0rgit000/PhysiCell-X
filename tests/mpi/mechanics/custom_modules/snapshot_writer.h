#ifndef __MECHANICS_TEST_SNAPSHOT_WRITER_H__
#define __MECHANICS_TEST_SNAPSHOT_WRITER_H__

#include <string>

#include "../../../../core/PhysiCell.h"

void initialize_snapshot_csv( const std::string& path );
void append_serial_snapshot( const std::string& path, int iteration, double time );
void append_mpi_snapshot( const std::string& path, int iteration, double time,
	DistPhy::mpi::mpi_Environment& world, DistPhy::mpi::mpi_Cartesian& cart_topo );

#endif
