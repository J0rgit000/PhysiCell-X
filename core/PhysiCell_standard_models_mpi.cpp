#include "PhysiCell_standard_models.h" 
#include "PhysiCell_cell.h" 

namespace PhysiCell{

void standard_elastic_contact_function( Cell* pC1, Phenotype& p1, std::vector<double> neighbour, double, mpi_Environment&, mpi_Cartesian& )
{
	if( pC1->position.size() != 3 || neighbour.size() != 3 )
	{ return; }

	std::vector<double> displacement = neighbour;
	displacement -= pC1->position;

	axpy( &(pC1->velocity) , p1.mechanics.attachment_elastic_constant , displacement );

	return;
}

//MPI parallel version
//TODO: Encapsulate left and right neighbor interactions
void standard_cell_cell_interactions( Cell* pCell, Phenotype& phenotype, double dt , mpi_Environment &world, mpi_Cartesian &cart_topo )
{
	if( phenotype.death.dead == true )
	{ return; }

	Cell* pC = NULL; 
	int type = -1; 
	double probability = 0.0; 

	// std::cout << "testing against " << pCell->state.neighbors.size() << " cells " << std::endl; 

	bool attacked = false; 
	for( int n=0; n < pCell->state.neighbors.size(); n++ )
	{
		pC = pCell->state.neighbors[n]; 
		type = pC->type; 
		if (std::isnan(pC->position[0])) std::cout << "Error in line " << __func__ << " " << __LINE__ << std::endl;
		if( pC->phenotype.death.dead == true )
		{
			// dead phagocytosis 
			probability = phenotype.cell_interactions.other_dead_phagocytosis_rate * dt; 
			if( UniformRandom() <= probability ) 
			{ pCell->ingest_cell(pC); 
			  //std::cout << "		Rank " << world.rank << " Cell with ID " << pC->ID << " is ingested " << std::endl;
			  if (std::isnan(pCell->position[0])) std::cout << "Error in line " << __func__ << " " << __LINE__ << std::endl;
			} 
		}
		else
		{
			// live phagocytosis
			probability = phenotype.cell_interactions.live_phagocytosis_rates[type] * dt;
			if (std::isnan(pC->position[0])) std::cout << "Error in line " << __func__ << " " << __LINE__ << std::endl;
			if (std::isnan(pCell->position[0])) std::cout << "Error in line " << __func__ << " " << __LINE__ << std::endl;
			if( UniformRandom() <= probability ) 
			{ pCell->ingest_cell(pC); 
			  //std::cout << "		Rank " << world.rank << " Cell with ID " << pC->ID << " is ingested " << std::endl;
			 if (std::isnan(pCell->position[0])) std::cout << "Error in line " << __func__ << " " << __LINE__ << std::endl;
			} 

			// attack 

			// assume you can only attack one cell at a time 
			probability = phenotype.cell_interactions.attack_rates[type] * dt;  
			if( UniformRandom() <= probability && attacked == false ) 
			{
				pCell->attack_cell(pC,dt);
				if (std::isnan(pCell->position[0])) std::cout << "Error in line " << __func__ << " " << __LINE__ << std::endl; 
				attacked = true;
				//std::cout << "		Rank " << world.rank << " Cell with ID " << pC->ID << " is attacked" << std::endl;
			} 

			// fusion 
			probability = phenotype.cell_interactions.fusion_rates[type] * dt;  
			if( UniformRandom() <= probability ) 
			{ pCell->fuse_cell(pC); 
				if (std::isnan(pCell->position[0])) std::cout << "Error in line " << __func__ << " " << __LINE__ << std::endl;
			  //std::cout << "		Rank " << world.rank << " Cell with ID " << pC->ID << " is fused " << std::endl;
			} 
		}

	

	}
//	std::cout << std::endl; 
	
	
	int x_dim = pCell->get_container()->underlying_mesh.x_coordinates.size(); 
	int y_dim = pCell->get_container()->underlying_mesh.y_coordinates.size();
	int z_dim = pCell->get_container()->underlying_mesh.z_coordinates.size();
	
	int local_vxl_inex = pCell->get_current_mechanics_voxel_index();
	//int position_voxel = (local_vxl_inex % (x_dim * y_dim)) % x_dim;  Layout warning!
	int position_voxel = (local_vxl_inex / (z_dim * y_dim));

	if (std::isnan(pCell->position[0])) std::cout << "Error in line " << __func__ << " " << __LINE__ << std::endl;

	if (position_voxel == 0) { //Left edge
		if (world.rank > 0) { //Not first rank 
			int yx_index = local_vxl_inex %(y_dim*z_dim);
			std::vector<int> moore_list = pCell->get_container()->underlying_mesh.moore_connected_voxel_global_indices_left[yx_index];
			for(int i=0; i<moore_list.size(); i++)
			{
				Interacting_Voxel &ivi = pCell->get_container()->um_ivfl.at(moore_list[i]);
					
				for(int cell_ctr=0; cell_ctr<ivi.cells.size(); cell_ctr++) {
					
					Interacting_Cell_Info * pC = &ivi.cells[cell_ctr];
					if (std::isnan(pC->position[0])) {
						std::cout << "Error in line " << __func__ << " " << __LINE__ << std::endl;
						continue;
					} 
					int type = pC->type; 
					if( pC->dead == true )
					{
						// dead phagocytosis 
						probability = 0; //phenotype.cell_interactions.dead_phagocytosis_rate * dt; 
						if( UniformRandom() <= probability ) 
						{ pCell->ingest_cell(pC); 
						  //std::cout << "		Rank " << world.rank << " remote Cell with ID " << pC->ID << " is ingested " << std::endl;
						} 
					}
					else
					{
						// live phagocytosis
						probability = phenotype.cell_interactions.live_phagocytosis_rates[type] * dt;  
						if( UniformRandom() <= probability ) 
						{ pCell->ingest_cell(pC); 
						  //std::cout << "		Rank " << world.rank << " remote Cell with ID " << pC->ID << " is ingested " << std::endl;
						} 

						// attack 

						// assume you can only attack one cell at a time 
						probability = phenotype.cell_interactions.attack_rates[type] * dt;  
						if( UniformRandom() <= probability && attacked == false ) 
						{
							pCell->attack_cell(pC,dt); 
							attacked = true;
							//std::cout << "		Rank " << world.rank << " remote Cell with ID " << pC->ID << " is attacked " << std::endl;
						} 

						// fusion 
						probability = phenotype.cell_interactions.fusion_rates[type] * dt;  
						if( UniformRandom() <= probability ) 
						{ pCell->fuse_cell(pC, world, cart_topo); 
						  //std::cout << "		Rank " << world.rank << " remote Cell with ID " << pC->ID << " is fused " << std::endl;
						} 
					}
				}				
			}
		}
	}
	if(position_voxel == (x_dim-1)) { //Right edge
		if(world.rank < world.size-1) { //

			int yx_index = local_vxl_inex %(y_dim*z_dim);
			std::vector<int> moore_list = pCell->get_container()->underlying_mesh.moore_connected_voxel_global_indices_right[yx_index];
			for(int i=0; i<moore_list.size(); i++)
			{
				Interacting_Voxel &ivi = pCell->get_container()->um_ivfr.at(moore_list[i]);
					
				for(int cell_ctr=0; cell_ctr<ivi.cells.size(); cell_ctr++) {
					
					Interacting_Cell_Info * pC = &ivi.cells[cell_ctr]; 
					int type = pC->type; 
					if( pC->dead == true )
					{
						// dead phagocytosis 
						probability = 0; // phenotype.cell_interactions.dead_phagocytosis_rate * dt; 
						if( UniformRandom() <= probability ) 
						{ 
						  pCell->ingest_cell(pC); 
						  //std::cout << "		Rank " << world.rank << " remote Cell with ID " << pC->ID << " is ingested " << std::endl;
						} 
					}
					else
					{
						// live phagocytosis
						probability = phenotype.cell_interactions.live_phagocytosis_rates[type] * dt;  
						if( UniformRandom() <= probability ) 
						{ 
						  pCell->ingest_cell(pC); 
						  //std::cout << "		Rank " << world.rank << " remote Cell with ID " << pC->ID << " is ingested " << std::endl;
						} 

						// attack 

						// assume you can only attack one cell at a time 
						probability = phenotype.cell_interactions.attack_rates[type] * dt;  
						if( UniformRandom() <= probability && attacked == false ) 
						{
							pCell->attack_cell(pC,dt); 
							attacked = true;
							//std::cout << "		Rank " << world.rank << " remote Cell with ID " << pC->ID << " is attacked" << std::endl;
						} 

						// fusion 
						probability = phenotype.cell_interactions.fusion_rates[type] * dt;  
						if( UniformRandom() <= probability ) 
						{ 
							pCell->fuse_cell(pC, world, cart_topo); 
						  	//std::cout << "		Rank " << world.rank << " remote Cell with ID " << pC->ID << " is fused " << std::endl;
						} 
					}
				}
			}
		}
	}
}

};
