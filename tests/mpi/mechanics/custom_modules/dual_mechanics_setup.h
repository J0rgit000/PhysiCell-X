#ifndef __DUAL_MECHANICS_SETUP_H__
#define __DUAL_MECHANICS_SETUP_H__

#include <vector>
#include "../../../modules/PhysiCell_standard_modules.h"
#include "./mechanics_comparison.h"

using namespace BioFVM;
using namespace PhysiCell;
using namespace DistPhy::mpi;

/*================================================================================
 Global references to both OMP and MPI cell containers for dual execution
 ================================================================================*/

extern Cell_Container* cell_container_omp;      // OMP-only container
extern Cell_Container* cell_container_mpi;      // MPI-distributed container
extern BioFVM::Microenvironment* microenv_omp;  // OMP-only microenvironment
extern BioFVM::Microenvironment* microenv_mpi;  // MPI-distributed microenvironment

/*================================================================================
 Setup functions for dual OMP and MPI execution
 ================================================================================*/

/**
 * Initialize OMP-only microenvironment
 * (copies domain from config but creates non-distributed mesh)
 */
void setup_microenvironment_omp( void );

/**
 * Initialize OMP-only cell container
 */
void setup_cell_container_omp( double mechanics_voxel_size );

/**
 * Distribute cells to both OMP and MPI containers identically
 * Each cell is added to both containers at the same (x, y, z)
 */
void distribute_cells_to_both_containers( 
	BioFVM::Microenvironment* env_omp,
	BioFVM::Microenvironment* env_mpi,
	Cell_Container* container_omp,
	Cell_Container* container_mpi,
	mpi_Environment &world,
	mpi_Cartesian &cart_topo
);

/**
 * Run one mechanics update iteration for both OMP and MPI
 * Returns comparison result for this iteration
 */
ComparisonResult run_dual_iteration( 
	double dt,
	mpi_Environment &world,
	mpi_Cartesian &cart_topo
);

/**
 * Finalize and write comparison reports
 */
void finalize_comparison_test( 
	const std::vector<ComparisonResult>& all_results,
	int total_iterations,
	mpi_Environment &world
);

#endif
