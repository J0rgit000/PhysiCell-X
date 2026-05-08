#!/bin/bash
#SBATCH --nodes=2
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=1
#SBATCH --qos=gp_debug
#SBATCH -t 02:00:00
#SBATCH --account=cns119
#SBATCH -o output-%j
#SBATCH -e error-%j
#SBATCH --exclusive

set -euo pipefail

export OMP_DISPLAY_ENV=false
export OMP_PROC_BIND=spread
export OMP_PLACES=threads
export OMP_NUM_THREADS=1

module purge
module load gcc/13.2.0 openmpi/4.1.5-gcc ddt

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

CONFIG_PATH=${1:-"$SCRIPT_DIR/config/validation.xml"}
REFERENCE_CSV=${2:-"$SCRIPT_DIR/reference_snapshots.csv"}
MPI_CSV=${3:-"$SCRIPT_DIR/mpi_snapshots.csv"}
REPORT_PATH=${4:-"$SCRIPT_DIR/comparison_report.txt"}

srun --nodes=1 --ntasks=1 --cpus-per-task=1 "$SCRIPT_DIR/mechanics_ref" "$CONFIG_PATH" "$REFERENCE_CSV"
srun --nodes=2 --ntasks=2 --cpus-per-task=1 "$SCRIPT_DIR/mechanics_val" "$CONFIG_PATH" "$MPI_CSV"
python3 "$SCRIPT_DIR/compare_snapshots.py" "$REFERENCE_CSV" "$MPI_CSV" "$REPORT_PATH"
