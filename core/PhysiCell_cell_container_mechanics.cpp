#include "../BioFVM/BioFVM_agent_container.h"
#include "PhysiCell_constants.h"
#include "../BioFVM/BioFVM_vector.h"
#include "PhysiCell_cell.h"
#include "PhysiCell_standard_models.h"
#include "PhysiCell_cell_container.h"


#include <algorithm>
#include <iterator> 

using namespace BioFVM;


namespace PhysiCell{


// Constructor that unpacks a cell
Moore_Cell_Info::Moore_Cell_Info(std::vector<char>& buffer, int size, int &pos)
{
    MPI_Unpack(buffer.data(), size, &pos, &ID, 1, MPI_INT, MPI_COMM_WORLD);
    position.resize(3);
    MPI_Unpack(buffer.data(), size, &pos, position.data(), 3, MPI_DOUBLE, MPI_COMM_WORLD);
    MPI_Unpack(buffer.data(), size, &pos, &radius, 1, MPI_DOUBLE, MPI_COMM_WORLD);
    MPI_Unpack(buffer.data(), size, &pos, &nuclear_radius, 1, MPI_DOUBLE, MPI_COMM_WORLD);
    MPI_Unpack(buffer.data(), size, &pos, &cell_cell_repulsion_strength, 1, MPI_DOUBLE, MPI_COMM_WORLD);
    MPI_Unpack(buffer.data(), size, &pos, &relative_maximum_adhesion_distance, 1, MPI_DOUBLE, MPI_COMM_WORLD);
    MPI_Unpack(buffer.data(), size, &pos, &cell_cell_adhesion_strength, 1, MPI_DOUBLE, MPI_COMM_WORLD);
    MPI_Unpack(buffer.data(), size, &pos, &type, 1, MPI_INT, MPI_COMM_WORLD);
}

Moore_Voxel_Info::Moore_Voxel_Info(std::vector<char>& buffer, int size, int &pos)
{
    MPI_Unpack(buffer.data(), size, &pos, &global_mesh_index, 1, MPI_INT, MPI_COMM_WORLD);
    MPI_Unpack(buffer.data(), size, &pos, &no_of_cells_in_vxl, 1, MPI_INT, MPI_COMM_WORLD);
    center.resize(3);
    MPI_Unpack(buffer.data(), size, &pos, center.data(), 3, MPI_DOUBLE, MPI_COMM_WORLD);
    MPI_Unpack(buffer.data(), size, &pos, &max_cell_interactive_distance_in_voxel, 1, MPI_DOUBLE, MPI_COMM_WORLD);

    moore_cells.reserve(no_of_cells_in_vxl);
    for (int i = 0; i < no_of_cells_in_vxl; ++i)
    {
        moore_cells.emplace_back(buffer, size, pos);
    }
}

//TODO: 1: encapsulate the packing of each voxel's data into a separate function to avoid code repetition
std::vector<Moore_Cell_Info> Cell_Container::get_moore_neighbour_cells(Cell* pCell, mpi_Environment &world, mpi_Cartesian &cart_topo){
    std::vector<Moore_Cell_Info> neighbours;
    int x_dim = pCell->get_container()->underlying_mesh.x_coordinates.size(); 
    int y_dim = pCell->get_container()->underlying_mesh.y_coordinates.size();
    int z_dim = pCell->get_container()->underlying_mesh.z_coordinates.size();
    
    int local_vxl_index = pCell->get_current_mechanics_voxel_index();
    int position_voxel = (local_vxl_index / (z_dim * y_dim));

    if (position_voxel == 0  && world.rank > 0) { //left edge
        int yx_index = local_vxl_index %(y_dim*z_dim);
        std::vector<int> moore_list = pCell->get_container()->underlying_mesh.moore_connected_voxel_global_indices_left[yx_index];
        for(int i=0; i<moore_list.size(); i++)
        {
            Moore_Voxel_Info &mvi = um_mbfl.at(moore_list[i]);
            for(int cell_ctr=0; cell_ctr<mvi.moore_cells.size(); cell_ctr++) {
                neighbours.push_back(mvi.moore_cells[cell_ctr]);
            }
        } 
    } 
    if (position_voxel == (x_dim-1) && (world.rank < world.size-1)) { //right edge
        int yx_index = local_vxl_index %(y_dim*z_dim);
        std::vector<int> moore_list = pCell->get_container()->underlying_mesh.moore_connected_voxel_global_indices_right[yx_index];
        for(int i=0; i<moore_list.size(); i++)
        {
            Moore_Voxel_Info &mvi = um_mbfr.at(moore_list[i]);
                
            for(int cell_ctr=0; cell_ctr<mvi.moore_cells.size(); cell_ctr++) {
                neighbours.push_back(mvi.moore_cells[cell_ctr]);
            }
        } 
    }
    return neighbours;
}

//TODO: 1: encapsulate the packing of each voxel's data into a separate function to avoid code repetition
std::vector<Interacting_Cell_Info> Cell_Container::get_neighbour_interacting_cells(Cell* pCell, mpi_Environment &world, mpi_Cartesian &cart_topo){
    std::vector<Interacting_Cell_Info> neighbours;
    int x_dim = pCell->get_container()->underlying_mesh.x_coordinates.size(); 
    int y_dim = pCell->get_container()->underlying_mesh.y_coordinates.size();
    int z_dim = pCell->get_container()->underlying_mesh.z_coordinates.size();
    
    int local_vxl_index = pCell->get_current_mechanics_voxel_index();
    int position_voxel = (local_vxl_index / (z_dim * y_dim));

    if (position_voxel == 0  && world.rank > 0) { //left edge
        int yx_index = local_vxl_index %(y_dim*z_dim);
        std::vector<int> moore_list = pCell->get_container()->underlying_mesh.moore_connected_voxel_global_indices_left[yx_index];
        for(int i=0; i<moore_list.size(); i++)
        {
            Interacting_Voxel &ivi = pCell->get_container()->um_ivfl.at(moore_list[i]);
                
            for(int cell_ctr=0; cell_ctr<ivi.cells.size(); cell_ctr++) {
                neighbours.push_back(ivi.cells[cell_ctr]);
            }
        } 
    } 
    else if (position_voxel == (x_dim-1) && (world.rank < world.size-1)) { //right edge
        int yx_index = local_vxl_index %(y_dim*z_dim);
        std::vector<int> moore_list = pCell->get_container()->underlying_mesh.moore_connected_voxel_global_indices_right[yx_index];
        for(int i=0; i<moore_list.size(); i++)
        {
            Interacting_Voxel &ivi = pCell->get_container()->um_ivfr.at(moore_list[i]);
                
            for(int cell_ctr=0; cell_ctr<ivi.cells.size(); cell_ctr++) {
                neighbours.push_back(ivi.cells[cell_ctr]);
            }
        } 
    }
    return neighbours;
}

void Cell_Container::evaluate_cell_elastic_interactions( PhysiCell::Cell *pCell , double dt_mec,  mpi_Environment &world, mpi_Cartesian &cart_topo){

    dynamic_spring_attachments(pCell,pCell->phenotype,dt_mec); 
    for( int j=0; j < pCell->state.spring_attachments.size(); j++ )
    {
        Cell* pC1 = pCell->state.spring_attachments[j]; 
		standard_elastic_contact_function(pCell,pCell->phenotype,pC1,pC1->phenotype,dt_mec);  
    }

    //elastic interactions with other subdomains
    //check if voxel is left or right before!

    int max_elastic = pCell->phenotype.mechanics.maximum_number_of_attachments;
    int done_elastic = pCell->state.spring_attachments.size();
    vector<Moore_Cell_Info> neighbours = get_moore_neighbour_cells(pCell, world, cart_topo);
    
     double attachment_probability = pCell->phenotype.mechanics.attachment_rate * dt_mec; 
    int i = 0;
    while (done_elastic < max_elastic && i < neighbours.size()) {
        double probability_of_attachment = attachment_probability * pCell->phenotype.mechanics.cell_adhesion_affinities[neighbours[i].type];
        if (uniform_random() < probability_of_attachment) {
            standard_elastic_contact_function(pCell,pCell->phenotype, neighbours[i].position ,dt_mec, world, cart_topo);
			++done_elastic;
		}
        ++i;
        
    }
}

void Cell_Container::exchange_mechanics_halos(mpi_Environment &world, mpi_Cartesian &cart_topo) {
    // Exchange ghost cells for mechanical potentials and cell-cell interactions once per mechanics step.
    pack_moore_info(world, cart_topo);
    pack_cell_interact_info(world, cart_topo);
}

void Cell_Container::evaluate_cell_cell_interactions(double time_since_last_mechanics,  mpi_Environment &world, mpi_Cartesian &cart_topo) {
    // Boundary ghost data must be available prior to calling (exchange_mechanics_halos).

    #pragma omp parallel for 
    for( int i=0; i < (*all_cells).size(); i++ )
    {
        Cell* pC = (*all_cells)[i];
		if (!pC->is_out_of_domain)
        	standard_cell_cell_interactions( pC, pC->phenotype, time_since_last_mechanics , world, cart_topo );
    }

    unpack_cell_interact_info(world, cart_topo);

}

void Cell_Container::update_cell_potentials(double time_since_last_mechanics,  mpi_Environment &world, mpi_Cartesian &cart_topo){

    // Boundary data must be exchanged before calling this routine.
    #pragma omp parallel for schedule(guided)
    for( int i=0; i < (*all_cells).size(); i++ )
    {
        Cell* pC = (*all_cells)[i];
        
        if( pC->functions.contact_function && pC->is_out_of_domain == false )
			{ evaluate_interactions( pC,pC->phenotype,time_since_last_mechanics ); } //No parallelism necessary yet


        if( pC->functions.custom_cell_rule && pC->is_out_of_domain == false )
			{ pC->functions.custom_cell_rule( pC,pC->phenotype,time_since_last_mechanics ); } //No parallelism necessary yet

        if( pC->is_movable )
        {
            if(pC->functions.update_velocity_parallel && pC->is_out_of_domain == false)
            {
                pC->functions.update_velocity_parallel( pC, pC->phenotype, time_since_last_mechanics, world, cart_topo);
            }
            if( PhysiCell_settings.disable_automated_spring_adhesions == false ) 
            {
                evaluate_cell_elastic_interactions( pC , time_since_last_mechanics, world, cart_topo);
            }
        }
    }	
}

void Cell_Container::pack(std::vector<Cell*> *all_cells, mpi_Environment &world, mpi_Cartesian &cart_topo)
{

	int len_snd_buf_left  	 = 0;
	int len_snd_buf_right 	 = 0;

	int 	temp_int;
	double 	temp_double;
	std::string temp_str;
	std::unordered_map<std::string, int> :: iterator it;
	Cell *pCell;


	int len_int  = 0;
	int len_double  = 0;
	int len_str = 0;
	int len_vector = 0;

	position_left 			= 0;					//Must be initialized to 0
	position_right 			= 0;					//Must be initialized to 0
	no_cells_cross_left  	= 0;					//Must be initialized to 0
	no_cells_cross_right	= 0;					//Must be initialized to 0
	no_of_cells_from_right 	= 0;
	no_of_cells_from_left 	= 0;

	snd_buf_left.resize(0); 					//When we enter function again, this is reset
	snd_buf_right.resize(0);					//Reset the others also
	rcv_buf_left.resize(0);
	rcv_buf_right.resize(0);
    snd_pos_left.resize(0);
    snd_pos_right.resize(0);
    rcv_pos_left.resize(0);
    rcv_pos_right.resize(0);


	/*--------------------------------------------------------------------------------*/
	/* First count the number of cells crossing to left subdomain and right subdomain */
	/*--------------------------------------------------------------------------------*/

	for(int i=0; i<(*all_cells).size();i++)
	{
		pCell = (*all_cells)[i];
		if(pCell->crossed_to_left_subdomain == true)
			no_cells_cross_left++;
		if(pCell->crossed_to_right_subdomain == true)
			no_cells_cross_right++;
	}

    snd_pos_left.resize(no_cells_cross_left);
    snd_pos_right.resize(no_cells_cross_right);
	/*-----------------------------------------------------------------------------------*/
	/* There is no need to pack crossed_to_left_subdomain and crossed_to_right_subdomain */
	/* as we send these two fields as the First Communication between MPI processes. 	 */
	/* Also try to send the length of the buffer to be allocated when sending.			 */
	/*-----------------------------------------------------------------------------------*/

  //std::cout<<"Total cells crossing to left in Rank "<<world.rank<<":"<<no_cells_cross_left<<std::endl;
	//std::cout<<"Total cells crossing to right in Rank "<<world.rank<<":"<<no_cells_cross_right<<std::endl;

	/* IMPORTANT: CANNOT USE #pragma omp for HERE AS ALL THREADS WILL WRITE TO THE SAME SHARED BUFFER */
	int left = 0;
	int right = 0;
	int right_prev, left_prev;
	left_prev = 0;
	right_prev = 0;
	for(int i=0; i<(*all_cells).size();i++)
	{
		pCell = (*all_cells)[i];

		if(pCell->crossed_to_left_subdomain == true)
		{
            snd_pos_left[left] = position_left;
			//pCell->pack(snd_buf_left, len_snd_buf_left, position_left);
			pCell->pack_light(snd_buf_left, len_snd_buf_left, position_left);
			left++;
		}
		else if(pCell->crossed_to_right_subdomain == true)
		{
            snd_pos_right[right] = position_right;
			//pCell->pack(snd_buf_right, len_snd_buf_right, position_right);
			pCell->pack_light(snd_buf_right, len_snd_buf_right, position_right);
			right++;
		}
	}
 }

void Cell_Container::unpack(mpi_Environment &world, mpi_Cartesian &cart_topo)
{

		/*---------------------------------------------------------------------------*/
		/* position_left/right are members of class Cell_Container and contain the 	 */
		/* size of the buffers to be sent to the left and right neighbour, respect.	 */
		/* Since their purpose gets over after packing and sending/receiving, they	 */
		/* can be used for storing the position we are at in the recv buffers				 */
		/* (See beginning of Cell_Container::pack() function for 10 data members		 */
		/* relevant to packing, sending/receiving, unpacking data)  								 */
		/*---------------------------------------------------------------------------*/

		position_left  = 0;
		position_right = 0;

		/* Temporary variables to help in unpacking */

		
		int size_left, size_right;
		int cell_ID = 0;
		double cell_position[3];

		int 	 				temp_int;
		double 				temp_double;
		std::string 	temp_str;
		char 					*temp_char_array;
		double 				*temp_double_array;
		int 					temp_key_value[2];

		std::unordered_map<std::string, int> :: iterator it;
		std::vector<double> temp_double_vector;

		int len_int 		 		= 0;
		int len_double 	 		= 0;
		int len_str 		 		= 0;
		int len_vector	 		= 0;
		int len_vector_nest = 0;
		
		int right = 0;
		int left = 0;
		auto start_time = std::chrono::high_resolution_clock::now();
		if(no_of_cells_from_right > 0)
		{
			size_right = rcv_buf_right.size();

			/*-------------------------------------------------------------------------------------------------*/
			/* IMPORTANT: need to start a for(int i=0; i<no_of_cells_from_right; i++) here BUT cannot put this */
			/* loop BEFORE completing the unpacking, as the next cell would get filled with wrong values 			 */
			/* so complete unpacking of 1 cell first THEN put the loop 																				 */
			/*-------------------------------------------------------------------------------------------------*/
                //#pragma omp parallel for
				for(int loop_ctr = 0; loop_ctr < no_of_cells_from_right; loop_ctr++)
				{
					Cell *pCell = create_cell(cell_ID);
                    
                    //int  pos = rcv_pos_right[loop_ctr]; 

					//pCell->unpack(rcv_buf_right, size_right, position_right);
					pCell->unpack_light(rcv_buf_right, size_right, position_right);
					++right;

					pCell->assign_position(pCell->position[0], pCell->position[1], pCell->position[2], world, cart_topo);
					/*
					std::cout << "\t\t\t\t[Rank " << world.rank << " ] Unpacked cell ID: " << pCell->ID << " at position: (" 
					          << pCell->position[0] << "," << pCell->position[1] << "," << pCell->position[2] << ")  velocity: (" 
					          << pCell->velocity[0] << "," << pCell->velocity[1] << "," << pCell->velocity[2] << ") from RIGHT" 
					          << std::endl;*/
		 		}	 
		}

		if(no_of_cells_from_left > 0)
		{
			 size_left = rcv_buf_left.size();

			/*-------------------------------------------------------------------------------------------------*/
			/* IMPORTANT: need to start a for(int i=0; i<no_of_cells_from_left; i++) here BUT cannot put this */
			/* loop BEFORE completing the unpacking, as the next cell would get filled with wrong values 			 */
			/* so complete unpacking of 1 cell first THEN put the loop 																				 */
			/*-------------------------------------------------------------------------------------------------*/
                //#pragma omp parallel for
				for(int loop_ctr = 0; loop_ctr < no_of_cells_from_left; loop_ctr++)
				{
					Cell *pCell = create_cell(cell_ID);

                    int  pos = rcv_pos_left[loop_ctr];

					//pCell->unpack(rcv_buf_left, size_left, position_left);
					pCell->unpack_light(rcv_buf_left, size_left, position_left);
					
					pCell->assign_position(pCell->position[0], pCell->position[1], pCell->position[2], world, cart_topo);
					++left;
					/*
					std::cout << "\t\t\t\t[Rank " << world.rank << " ] Unpacked cell ID: " << pCell->ID << " at position: (" 
					          << pCell->position[0] << "," << pCell->position[1] << "," << pCell->position[2] << ")  velocity: (" 
					          << pCell->velocity[0] << "," << pCell->velocity[1] << "," << pCell->velocity[2] << ") from LEFT" 
					          << std::endl;*/
				}
		}
		auto end_time = std::chrono::high_resolution_clock::now();
		auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
		/*
		if (no_of_cells_from_left + no_of_cells_from_right > 0) {
			std::cout << "\t [Rank " << world.rank << " ] Unpacking completed in avg time per cell " << duration.count()/(no_of_cells_from_left+no_of_cells_from_right) << " microseconds. " << std::endl;
			std::cout << "\t [Rank " << world.rank << " ] Unpacking completed total time " << duration.count() << " microseconds for cells: " << no_of_cells_from_left + no_of_cells_from_right << std::endl;
		}*/
}

void Cell_Container::unpack_parallel(mpi_Environment &world, mpi_Cartesian &cart_topo)
{
	int len_rcv_buf_left = rcv_buf_left.size();
	int len_rcv_buf_right = rcv_buf_right.size();

	std::vector<Cell*> new_cells_from_left(no_of_cells_from_left, nullptr);
	std::vector<Cell*> new_cells_from_right(no_of_cells_from_right, nullptr);
	auto start_time = std::chrono::high_resolution_clock::now();
	#pragma omp parallel for schedule(static)
	for( int i = 0; i < no_of_cells_from_left; i++ )
	{
		// Use the received byte offsets to unpack independent cells concurrently.
		Cell* pCell = new Cell(0);
		int pos = rcv_pos_left[i];
		//pCell->unpack(rcv_buf_left, len_rcv_buf_left, pos);
		pCell->unpack_light(rcv_buf_left, len_rcv_buf_left, pos);
		new_cells_from_left[i] = pCell;
	}

	#pragma omp parallel for schedule(static)
	for( int i = 0; i < no_of_cells_from_right; i++ )
	{
		Cell* pCell = new Cell(0);
		int pos = rcv_pos_right[i];
		//pCell->unpack(rcv_buf_right, len_rcv_buf_right, pos);
		pCell->unpack_light(rcv_buf_right, len_rcv_buf_right, pos);
		new_cells_from_right[i] = pCell;
	}
	auto end_time = std::chrono::high_resolution_clock::now();
	auto duration_part = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
	/*
	if (no_of_cells_from_left + no_of_cells_from_right > 0) {
		std::cout << "\t [Rank " << world.rank << " ] Unpacking completed in avg time per cell " << duration_part.count()/(no_of_cells_from_left+no_of_cells_from_right) << " microseconds. " << std::endl;
		std::cout << "\t [Rank " << world.rank << " ] Unpacking completed total time " << duration_part.count() << " microseconds for cells: " << no_of_cells_from_left + no_of_cells_from_right << std::endl;
	}*/
	// Integration into global cell storage and mechanics voxels remains serialized
	// because all_cells and agent_grid are shared mutable containers.
	start_time = std::chrono::high_resolution_clock::now();
	for( int i = 0; i < no_of_cells_from_right; i++ )
	{
		Cell* pCell = new_cells_from_right[i];
		(*all_cells).push_back( pCell );
		pCell->index = (*all_cells).size() - 1;
		pCell->assign_position(pCell->position[0], pCell->position[1], pCell->position[2], world, cart_topo);
	}

	for( int i = 0; i < no_of_cells_from_left; i++ )
	{
		Cell* pCell = new_cells_from_left[i];
		(*all_cells).push_back( pCell );
		pCell->index = (*all_cells).size() - 1;
		pCell->assign_position(pCell->position[0], pCell->position[1], pCell->position[2], world, cart_topo);
	}
	end_time = std::chrono::high_resolution_clock::now();
	auto duration_integration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
	//if (no_of_cells_from_left + no_of_cells_from_right > 0) {
	//	std::cout << "\t [Rank " << world.rank << " ] Integration of unpacked cells completed in avg time per cell " << duration_integration.count()/(no_of_cells_from_left+no_of_cells_from_right) << " microseconds. " << std::endl;
	//	std::cout << "\t [Rank " << world.rank << " ] Integration of unpacked cells completed total time " << duration_integration.count() << " microseconds for cells: " << no_of_cells_from_left + no_of_cells_from_right << std::endl;
	//}
}

} // namespace PhysiCell
