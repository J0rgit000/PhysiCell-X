#!/usr/bin/env python3

import csv
import math
import sys
from pathlib import Path

POSITION_TOL = 1e-3
VELOCITY_TOL = 1e-5


def parse_neighbor_ids(raw_value: str):
    if not raw_value:
        return ()
    return tuple(sorted(int(value) for value in raw_value.split(";") if value))


def load_snapshots(path: Path):
    snapshots = {}
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            iteration = int(row["iteration"])
            cell_id = int(row["cell_id"])
            snapshots.setdefault(iteration, {})[cell_id] = {
                "time": float(row["time"]),
                "position": (
                    float(row["x"]),
                    float(row["y"]),
                    float(row["z"]),
                ),
                "velocity": (
                    float(row["vx"]),
                    float(row["vy"]),
                    float(row["vz"]),
                ),
                "neighbors": parse_neighbor_ids(row["neighbor_ids"]),
            }
    return snapshots


def vector_error(left, right):
    return math.sqrt(sum((a - b) ** 2 for a, b in zip(left, right)))


def compare(reference, mpi):
    report_lines = []
    ref_iterations = sorted(reference)
    mpi_iterations = sorted(mpi)

    if ref_iterations != mpi_iterations:
        report_lines.append("FAIL")
        report_lines.append("First failing iteration: n/a")
        report_lines.append("Failing cell IDs: n/a")
        report_lines.append(f"Iteration mismatch: reference={ref_iterations}, mpi={mpi_iterations}")
        return False, "\n".join(report_lines)

    worst_position_error = 0.0
    worst_velocity_error = 0.0
    first_failure_iteration = None
    failing_cells = []

    for iteration in ref_iterations:
        ref_cells = reference[iteration]
        mpi_cells = mpi[iteration]

        ref_ids = sorted(ref_cells)
        mpi_ids = sorted(mpi_cells)

        if len(ref_ids) != len(mpi_ids):
            report_lines.append("FAIL")
            report_lines.append(f"First failing iteration: {iteration}")
            report_lines.append("Failing cell IDs: n/a")
            report_lines.append(
                f"Iteration {iteration}: cell count mismatch reference={len(ref_ids)} mpi={len(mpi_ids)}"
            )
            return False, "\n".join(report_lines)

        if ref_ids != mpi_ids:
            report_lines.append("FAIL")
            report_lines.append(f"First failing iteration: {iteration}")
            report_lines.append("Failing cell IDs: n/a")
            report_lines.append(f"Iteration {iteration}: cell ID mismatch reference={ref_ids} mpi={mpi_ids}")
            return False, "\n".join(report_lines)

        for cell_id in ref_ids:
            ref_row = ref_cells[cell_id]
            mpi_row = mpi_cells[cell_id]

            if ref_row["neighbors"] != mpi_row["neighbors"]:
                if first_failure_iteration is None:
                    first_failure_iteration = iteration
                failing_cells.append(cell_id)
                report_lines.append("FAIL")
                report_lines.append(f"First failing iteration: {first_failure_iteration}")
                report_lines.append(f"Failing cell IDs: {sorted(set(failing_cells))}")
                report_lines.append(
                    f"Iteration {iteration}, cell {cell_id}: neighbor mismatch "
                    f"reference={ref_row['neighbors']} mpi={mpi_row['neighbors']}"
                )
                return False, "\n".join(report_lines)

            position_error = vector_error(ref_row["position"], mpi_row["position"])
            velocity_error = vector_error(ref_row["velocity"], mpi_row["velocity"])
            worst_position_error = max(worst_position_error, position_error)
            worst_velocity_error = max(worst_velocity_error, velocity_error)

            if position_error > POSITION_TOL or velocity_error > VELOCITY_TOL:
                if first_failure_iteration is None:
                    first_failure_iteration = iteration
                failing_cells.append(cell_id)
                report_lines.append("FAIL")
                report_lines.append(f"First failing iteration: {first_failure_iteration}")
                report_lines.append(f"Failing cell IDs: {sorted(set(failing_cells))}")
                report_lines.append(
                    f"Iteration {iteration}, cell {cell_id}: "
                    f"position_error={position_error:.6e}, velocity_error={velocity_error:.6e}"
                )
                report_lines.append(
                    f"Worst errors so far: position={worst_position_error:.6e}, velocity={worst_velocity_error:.6e}"
                )
                return False, "\n".join(report_lines)

    report_lines.append("PASS")
    report_lines.append(f"Iterations compared: {len(ref_iterations)}")
    report_lines.append(f"Cells per iteration: {len(reference[ref_iterations[0]]) if ref_iterations else 0}")
    report_lines.append(f"Worst position error: {worst_position_error:.6e}")
    report_lines.append(f"Worst velocity error: {worst_velocity_error:.6e}")
    return True, "\n".join(report_lines)


def main(argv):
    if len(argv) not in (3, 4):
        print("Usage: compare_snapshots.py reference.csv mpi.csv [report.txt]", file=sys.stderr)
        return 2

    reference_path = Path(argv[1])
    mpi_path = Path(argv[2])
    report_path = Path(argv[3]) if len(argv) == 4 else None

    reference = load_snapshots(reference_path)
    mpi = load_snapshots(mpi_path)
    success, report = compare(reference, mpi)

    if report_path is not None:
        report_path.write_text(report + "\n")

    print(report)
    return 0 if success else 1


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
