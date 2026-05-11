#ifndef __MECHANICS_COMPARISON_H__
#define __MECHANICS_COMPARISON_H__

#include <vector>
#include <string>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <map>

using namespace BioFVM;
using namespace PhysiCell;

/*================================================================================
 Cell Snapshot Structure: captures state at a given iteration
 ================================================================================*/

struct CellSnapshot 
{
	int cell_id;
	double position[3];        // [x, y, z]
	double velocity[3];        // [vx, vy, vz]
	int num_attachments;       // adhesion count
	double max_force;          // max force magnitude this iteration
	bool valid;                // whether cell exists
};

/*================================================================================
 Comparison Result Structure: holds metrics comparing OMP vs MPI
 ================================================================================*/

struct ComparisonResult 
{
	// Flags
	bool position_match;           // RMS position error < threshold
	bool velocity_match;           // RMS velocity error < threshold
	bool force_match;              // Force statistics match
	bool topology_match;           // Cell neighbor graph matches
	bool cell_count_match;         // Same number of cells
	
	// Metrics
	double max_pos_error;          // Maximum position divergence across all cells
	double rms_pos_error;          // RMS position error
	double max_vel_error;          // Maximum velocity divergence
	double rms_vel_error;          // RMS velocity error
	double force_relative_error;   // Relative error in force magnitudes
	double topology_agreement;     // Fraction of matching neighbor pairs (0.0 to 1.0)
	int cells_out_of_tolerance;    // Count of cells exceeding position threshold
	int omp_cell_count;
	int mpi_cell_count;
	
	// Report
	std::string report;            // Human-readable summary
	
	// Constructor
	ComparisonResult() : 
		position_match(false), velocity_match(false), force_match(false), 
		topology_match(false), cell_count_match(false),
		max_pos_error(0.0), rms_pos_error(0.0), max_vel_error(0.0), rms_vel_error(0.0),
		force_relative_error(0.0), topology_agreement(0.0), 
		cells_out_of_tolerance(0), omp_cell_count(0), mpi_cell_count(0),
		report("") {}
};

/*================================================================================
 Function Declarations
 ================================================================================*/

/**
 * Capture all cells from a cell container
 * Returns vector of CellSnapshot sorted by cell_id
 */
std::vector<CellSnapshot> capture_all_cells( Cell_Container* container );

/**
 * Compare snapshots from OMP and MPI at a given iteration
 * Computes all comparison metrics
 */
ComparisonResult compare_snapshots( 
	const std::vector<CellSnapshot>& omp_snap,
	const std::vector<CellSnapshot>& mpi_snap,
	double pos_tol_um = 1.0,      // Position tolerance in micrometers
	double vel_tol = 0.1           // Velocity tolerance in micrometers/min
);

/**
 * Write snapshots to CSV file for inspection
 */
void write_snapshot_to_csv( 
	const std::string& filename,
	const std::vector<CellSnapshot>& snapshots,
	int iteration 
);

/**
 * Write comparison results to report file
 */
void write_comparison_report( 
	const std::string& filename,
	const std::vector<ComparisonResult>& results,
	int total_iterations 
);

/**
 * Determine if iteration should be sampled (exponential sampling)
 * Samples at: 0, 1, 5, 10, 50, 100, 500, 1000, ...
 */
bool should_sample_iteration( int iteration );

#endif
