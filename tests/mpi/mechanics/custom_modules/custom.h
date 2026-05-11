#ifndef __MECHANICS_TEST_CUSTOM_H__
#define __MECHANICS_TEST_CUSTOM_H__

#include "../../../../core/PhysiCell.h"
#include "../../../../modules/PhysiCell_standard_modules.h"

using namespace BioFVM;
using namespace PhysiCell;

void create_cell_types( void );
void create_cell_types( DistPhy::mpi::mpi_Environment& world, DistPhy::mpi::mpi_Cartesian& cart_topo );

void setup_microenvironment( void );
void setup_microenvironment( DistPhy::mpi::mpi_Environment& world, DistPhy::mpi::mpi_Cartesian& cart_topo );

void setup_tissue( void );
void setup_tissue( Microenvironment& microenvironment, DistPhy::mpi::mpi_Environment& world,
	DistPhy::mpi::mpi_Cartesian& cart_topo );

#endif
