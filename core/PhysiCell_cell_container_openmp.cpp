#include "PhysiCell_cell_container.h"

#include <omp.h>
#include <vector>

namespace PhysiCell {


struct ThreadPackingBuffer {
    std::vector<char> buffer;
    int position = 0;
    int length = 0;
};


inline int calculate_voxel_pack_size(int no_of_cells) {
    // Header size: 2 ints + 3 doubles + 1 double = 2*4 + 4*8 = 40 bytes
    int header_size = 2 * sizeof(int) + 4 * sizeof(double);
    
    // Size per cell: 2 ints + 8 doubles = 2*4 + 8*8 = 72 bytes
    int cell_size = 2 * sizeof(int) + 8 * sizeof(double);
    
    return header_size + (no_of_cells * cell_size);
}


// Optimized function to pack voxel info using offset based packing
void pack_moore_voxel_to_buffer(
    const Cell_Container* container,
    unsigned int voxel_index,
    std::vector<char>& snd_buffer,
    int buffer_offset) {
    
    if (voxel_index < 0 || voxel_index >= container->underlying_mesh.voxels.size()) {
        return;
    }
    
    int position = buffer_offset;
    int global_mesh_index = container->underlying_mesh.voxels[voxel_index].global_mesh_index;
    int no_of_cells_in_vxl = container->agent_grid[voxel_index].size();
    const std::vector<double>& centers = container->underlying_mesh.voxels[voxel_index].center;
    double max_i_dist = container->max_cell_interactive_distance_in_voxel[voxel_index];

    // Pack voxel header
    MPI_Pack(&global_mesh_index, 1, MPI_INT, &snd_buffer[0], snd_buffer.size(), &position, MPI_COMM_WORLD);
    MPI_Pack(&no_of_cells_in_vxl, 1, MPI_INT, &snd_buffer[0], snd_buffer.size(), &position, MPI_COMM_WORLD);
    MPI_Pack(&centers[0], 3, MPI_DOUBLE, &snd_buffer[0], snd_buffer.size(), &position, MPI_COMM_WORLD);
    MPI_Pack(&max_i_dist, 1, MPI_DOUBLE, &snd_buffer[0], snd_buffer.size(), &position, MPI_COMM_WORLD);

    // Pack cells in voxel
    for (int vec_len = 0; vec_len < no_of_cells_in_vxl; vec_len++) {
        Cell* pCell = container->agent_grid[voxel_index][vec_len];

        MPI_Pack(&pCell->ID, 1, MPI_INT, &snd_buffer[0], snd_buffer.size(), &position, MPI_COMM_WORLD);
        MPI_Pack(&pCell->position[0], 3, MPI_DOUBLE, &snd_buffer[0], snd_buffer.size(), &position, MPI_COMM_WORLD);
        MPI_Pack(&pCell->phenotype.geometry.radius, 1, MPI_DOUBLE, &snd_buffer[0], snd_buffer.size(), &position, MPI_COMM_WORLD);
        MPI_Pack(&pCell->phenotype.geometry.nuclear_radius, 1, MPI_DOUBLE, &snd_buffer[0], snd_buffer.size(), &position, MPI_COMM_WORLD);
        MPI_Pack(&pCell->phenotype.mechanics.cell_cell_repulsion_strength, 1, MPI_DOUBLE, &snd_buffer[0], snd_buffer.size(), &position, MPI_COMM_WORLD);
        MPI_Pack(&pCell->phenotype.mechanics.relative_maximum_adhesion_distance, 1, MPI_DOUBLE, &snd_buffer[0], snd_buffer.size(), &position, MPI_COMM_WORLD);
        MPI_Pack(&pCell->phenotype.mechanics.cell_cell_adhesion_strength, 1, MPI_DOUBLE, &snd_buffer[0], snd_buffer.size(), &position, MPI_COMM_WORLD);
        MPI_Pack(&pCell->type, 1, MPI_INT, &snd_buffer[0], snd_buffer.size(), &position, MPI_COMM_WORLD);
    }
}

/*
 * pack_moore_info_parallel() - Do not use by now!
 * TODO: fix parallel bug
 * lines 163-173 (line 163): offset = position_left_thread is read outside 
 * the critical section. Multiple threads can read the same offset, then pack 
 * different voxels into overlapping buffer regions. Same bug on lines 190-200 (line 190) 
 * for the right buffer.
 * Suggested fix: precompute deterministic per-voxel prefix offsets before the parallel loop, 
 * like the serial implementation already does in PhysiCell_cell_container.cpp (line 837). 
 * Do not assign offsets through a shared counter inside the parallel loop.
 */
void Cell_Container::pack_moore_info_parallel(mpi_Environment &world, mpi_Cartesian &cart_topo) {

    position_left = 0;
    position_right = 0;
    snd_buf_left.resize(0);
    snd_buf_right.resize(0);

    int z_dim = underlying_mesh.z_coordinates.size();
    int y_dim = underlying_mesh.y_coordinates.size();
    int x_dim = underlying_mesh.x_coordinates.size();

    /*─────────────────────────────────────────────────────────────────────*/
    /*           FASE 1: PRE-CALCULAR TAMAÑOS DE BUFFERS                  */
    /*─────────────────────────────────────────────────────────────────────*/

    int total_size_left = 0;
    int total_size_right = 0;

    // Pre-calcular tamaño total para buffer LEFT (sin paralelización aquí,
    // los cálculos son rápidos)
    if (world.rank > 0) {
        for (int k = 0; k < z_dim; k++) {
            for (int j = 0; j < y_dim; j++) {
                unsigned int local_vxl_index = underlying_mesh.voxel_index(0, j, k);
                if (local_vxl_index >= 0 && local_vxl_index < agent_grid.size()) {
                    int no_of_cells = agent_grid[local_vxl_index].size();
                    total_size_left += calculate_voxel_pack_size(no_of_cells);
                }
            }
        }
    }

    // Pre-calcular tamaño total para buffer RIGHT
    if (world.rank < world.size - 1) {
        for (int k = 0; k < z_dim; k++) {
            for (int j = 0; j < y_dim; j++) {
                unsigned int local_vxl_index = underlying_mesh.voxel_index(x_dim - 1, j, k);
                if (local_vxl_index >= 0 && local_vxl_index < agent_grid.size()) {
                    int no_of_cells = agent_grid[local_vxl_index].size();
                    total_size_right += calculate_voxel_pack_size(no_of_cells);
                }
            }
        }
    }

    /*─────────────────────────────────────────────────────────────────────*/
    /*         FASE 2: RESIZE UNA SOLA VEZ (NO EN LOOP)                   */
    /*─────────────────────────────────────────────────────────────────────*/

    if (total_size_left > 0) {
        snd_buf_left.resize(total_size_left);
    }
    if (total_size_right > 0) {
        snd_buf_right.resize(total_size_right);
    }

    /*─────────────────────────────────────────────────────────────────────*/
    /*         FASE 3: EMPACAR DATOS CON PARALELIZACIÓN OPENMP            */
    /*─────────────────────────────────────────────────────────────────────*/

    int len_snd_buf_left = 0;
    int len_snd_buf_right = 0;
    int position_left_thread = 0;
    int position_right_thread = 0;

    // PACKING LEFT BOUNDARY (paralelizable)
    if (world.rank > 0) {
        #pragma omp parallel for collapse(2) reduction(+:len_snd_buf_left) 
        for (int k = 0; k < z_dim; k++) {
            for (int j = 0; j < y_dim; j++) {
                unsigned int local_vxl_index = underlying_mesh.voxel_index(0, j, k);
                if (local_vxl_index >= 0 && local_vxl_index < agent_grid.size()) {
                    int no_of_cells = agent_grid[local_vxl_index].size();
                    int voxel_size = calculate_voxel_pack_size(no_of_cells);
                    
                    // Calcular offset para este thread
                    int offset = position_left_thread;
                    
                    #pragma omp critical
                    {
                        position_left_thread += voxel_size;
                        len_snd_buf_left += voxel_size;
                    }
                    
                    // Empacar en offset específico (thread-safe)
                    pack_moore_voxel_to_buffer(this, local_vxl_index, snd_buf_left, offset);
                }
            }
        }
        position_left = len_snd_buf_left;
    }

    // PACKING RIGHT BOUNDARY (paralelizable)
    if (world.rank < world.size - 1) {
        #pragma omp parallel for collapse(2) reduction(+:len_snd_buf_right)
        for (int k = 0; k < z_dim; k++) {
            for (int j = 0; j < y_dim; j++) {
                unsigned int local_vxl_index = underlying_mesh.voxel_index(x_dim - 1, j, k);
                if (local_vxl_index >= 0 && local_vxl_index < agent_grid.size()) {
                    int no_of_cells = agent_grid[local_vxl_index].size();
                    int voxel_size = calculate_voxel_pack_size(no_of_cells);
                    
                    // Calcular offset para este thread
                    int offset = position_right_thread;
                    
                    #pragma omp critical
                    {
                        position_right_thread += voxel_size;
                        len_snd_buf_right += voxel_size;
                    }
                    
                    // Empacar en offset específico (thread-safe)
                    pack_moore_voxel_to_buffer(this, local_vxl_index, snd_buf_right, offset);
                }
            }
        }
        position_right = len_snd_buf_right;
    }

    /*─────────────────────────────────────────────────────────────────────*/
    /*         FASE 4: COMUNICACIÓN MPI (igual a original)                 */
    /*─────────────────────────────────────────────────────────────────────*/

    int size_of_data_recvd_from_right_process = 0;
    int size_of_data_recvd_from_left_process = 0;
    MPI_Request snd_req[2], rcv_req[2];

    // Send size info
    MPI_Irecv(&size_of_data_recvd_from_right_process, 1, MPI_INT, cart_topo.X_RIGHT, 1111, cart_topo.mpi_cart_comm, &rcv_req[0]);
    MPI_Isend(&position_left, 1, MPI_INT, cart_topo.X_LEFT, 1111, cart_topo.mpi_cart_comm, &snd_req[0]);

    MPI_Irecv(&size_of_data_recvd_from_left_process, 1, MPI_INT, cart_topo.X_LEFT, 2222, cart_topo.mpi_cart_comm, &rcv_req[1]);
    MPI_Isend(&position_right, 1, MPI_INT, cart_topo.X_RIGHT, 2222, cart_topo.mpi_cart_comm, &snd_req[1]);

    MPI_Waitall(2, snd_req, MPI_STATUSES_IGNORE);
    MPI_Waitall(2, rcv_req, MPI_STATUSES_IGNORE);

    // Resize receive buffers
    if (world.rank < world.size - 1)
        rcv_buf_right.resize(size_of_data_recvd_from_right_process);

    if (world.rank > 0)
        rcv_buf_left.resize(size_of_data_recvd_from_left_process);

    // Send actual data
    char* rcv_right_ptr = rcv_buf_right.empty() ? nullptr : rcv_buf_right.data();
    char* rcv_left_ptr = rcv_buf_left.empty() ? nullptr : rcv_buf_left.data();
    char* snd_left_ptr = snd_buf_left.empty() ? nullptr : snd_buf_left.data();
    char* snd_right_ptr = snd_buf_right.empty() ? nullptr : snd_buf_right.data();

    MPI_Irecv(rcv_right_ptr, size_of_data_recvd_from_right_process, MPI_PACKED, cart_topo.X_RIGHT, 3333, cart_topo.mpi_cart_comm, &rcv_req[0]);
    MPI_Isend(snd_left_ptr, len_snd_buf_left, MPI_PACKED, cart_topo.X_LEFT, 3333, cart_topo.mpi_cart_comm, &snd_req[0]);

    MPI_Irecv(rcv_left_ptr, size_of_data_recvd_from_left_process, MPI_PACKED, cart_topo.X_LEFT, 4444, cart_topo.mpi_cart_comm, &rcv_req[1]);
    MPI_Isend(snd_right_ptr, len_snd_buf_right, MPI_PACKED, cart_topo.X_RIGHT, 4444, cart_topo.mpi_cart_comm, &snd_req[1]);

    MPI_Waitall(2, snd_req, MPI_STATUSES_IGNORE);
    MPI_Waitall(2, rcv_req, MPI_STATUSES_IGNORE);

    /*─────────────────────────────────────────────────────────────────────*/
    /*         FASE 5: UNPACKING (igual a original)                       */
    /*─────────────────────────────────────────────────────────────────────*/

    position_right = 0;
    position_left = 0;

    int size_right = size_of_data_recvd_from_right_process;
    int size_left = size_of_data_recvd_from_left_process;

    um_mbfr.clear();
    um_mbfl.clear();

    // Unpack from right
    if (world.rank < world.size - 1) {
        for (int vxl_ctr = 0; vxl_ctr < y_dim * z_dim; vxl_ctr++) {
            Moore_Voxel_Info voxel(rcv_buf_right, size_right, position_right);
            um_mbfr[voxel.global_mesh_index] = std::move(voxel);
        }
    }

    // Unpack from left
    if (world.rank > 0) {
        for (int vxl_ctr = 0; vxl_ctr < y_dim * z_dim; vxl_ctr++) {
            Moore_Voxel_Info voxel(rcv_buf_left, size_left, position_left);
            um_mbfl[voxel.global_mesh_index] = std::move(voxel);
        }
    }
}

} // namespace PhysiCell
