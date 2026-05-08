/*
#############################################################################
# If you use BioFVM in your project, please cite BioFVM and the version     #
# number, such as below:                                                    #
#                                                                           #
# We solved the diffusion equations using BioFVM (Version 1.1.7) [1]        #
#                                                                           #
# [1] A. Ghaffarizadeh, S.H. Friedman, and P. Macklin, BioFVM: an efficient #
#    parallelized diffusive transport solver for 3-D biological simulations,#
#    Bioinformatics 32(8): 1256-8, 2016. DOI: 10.1093/bioinformatics/btv730 #
#                                                                           #
#############################################################################
#                                                                           #
# BSD 3-Clause License (see https://opensource.org/licenses/BSD-3-Clause)   #
#                                                                           #
# Copyright (c) 2015-2017, Paul Macklin and the BioFVM Project              #
# All rights reserved.                                                      #
#                                                                           #
# Redistribution and use in source and binary forms, with or without        #
# modification, are permitted provided that the following conditions are    #
# met:                                                                      #
#                                                                           #
# 1. Redistributions of source code must retain the above copyright notice, #
# this list of conditions and the following disclaimer.                     #
#                                                                           #
# 2. Redistributions in binary form must reproduce the above copyright      #
# notice, this list of conditions and the following disclaimer in the       #
# documentation and/or other materials provided with the distribution.      #
#                                                                           #
# 3. Neither the name of the copyright holder nor the names of its          #
# contributors may be used to endorse or promote products derived from this #
# software without specific prior written permission.                       #
#                                                                           #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS       #
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED #
# TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A           #
# PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER #
# OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,  #
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,       #
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR        #
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF    #
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING      #
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS        #
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.              #
#                                                                           #
#############################################################################
*/

#include "BioFVM_solvers.h" 
#include "BioFVM_vector.h" 

#include <iostream>
#include <omp.h>
#include <immintrin.h>
#include <cstdint>

namespace BioFVM{

// do I even need this? 
void diffusion_decay_solver__constant_coefficients_explicit( Microenvironment& M, double dt )
{
	static bool precomputations_and_constants_done = false; 
	if( !precomputations_and_constants_done )
	{
		std::cout	<< std::endl << "Using solver: " << __FUNCTION__ << std::endl 
					<< "     (constant diffusion coefficient with explicit stepping, implicit decay) ... " << std::endl << std::endl;  

		if( M.mesh.uniform_mesh == true )
		{
			std::cout << "Uniform mesh detected! Consider switching to a more efficient method, such as " << std::endl  
			<< "     diffusion_decay_solver__constant_coefficients_explicit_uniform_mesh" << std::endl  
			<< std::endl; 
		}

		precomputations_and_constants_done = true; 
	}

	return; 
}

void diffusion_decay_solver__constant_coefficients_explicit_uniform_mesh( Microenvironment& M, double dt )
{
	static bool precomputations_and_constants_done = false; 
	if( !precomputations_and_constants_done )
	{
		std::cout	<< std::endl << "Using solver: " << __FUNCTION__ << std::endl 
					<< "     (constant diffusion coefficient with explicit stepping, implicit decay, uniform mesh) ... " << std::endl << std::endl;  

		if( M.mesh.regular_mesh == false )
		{ std::cout << "Error. This code is only supported for regular meshes." << std::endl; }

		precomputations_and_constants_done = true; 
	}

	return; 
}

void diffusion_decay_solver__constant_coefficients_LOD_3D( Microenvironment& M, double dt )
{
	if( M.mesh.regular_mesh == false || M.mesh.Cartesian_mesh == false )
	{
		std::cout << "Error: This algorithm is written for regular Cartesian meshes. Try: other solvers!" << std::endl << std::endl; 
		return; 
	}

	if( !M.diffusion_solver_setup_done )
	{
		std::cout << std::endl << "Using method " << __FUNCTION__ << " (implicit 3-D LOD with Thomas Algorithm) ... "
		<< std::endl << std::endl;  

		M.n_subs = M.number_of_densities(); 

		M.thomas_denomx.resize( M.mesh.x_coordinates.size() , M.zero );
		M.thomas_cx.resize( M.mesh.x_coordinates.size() , M.zero );

		M.thomas_denomy.resize( M.mesh.y_coordinates.size() , M.zero );
		M.thomas_cy.resize( M.mesh.y_coordinates.size() , M.zero );

		M.thomas_denomz.resize( M.mesh.z_coordinates.size() , M.zero );
		M.thomas_cz.resize( M.mesh.z_coordinates.size() , M.zero );

		M.thomas_i_jump = M.mesh.y_coordinates.size() * M.mesh.z_coordinates.size() * M.n_subs; 
		M.thomas_j_jump = M.mesh.z_coordinates.size() * M.n_subs; 
		M.thomas_k_jump = M.n_subs; 

		M.thomas_constant1 = M.diffusion_coefficients; // dt*D/dx^2 
		M.thomas_constant1a = M.zero; // -dt*D/dx^2; 
		M.thomas_constant2 = M.decay_rates; // (1/3)*dt*lambda 
		M.thomas_constant3 = M.one; // 1 + 2*constant1 + constant2; 
		M.thomas_constant3a = M.one; // 1 + constant1 + constant2; 

		M.thomas_constant1 *= dt; 
		M.thomas_constant1 /= M.mesh.dx; 
		M.thomas_constant1 /= M.mesh.dx; 

		M.thomas_constant1a = M.thomas_constant1; 
		M.thomas_constant1a *= -1.0; 

		M.thomas_constant2 *= dt; 
		M.thomas_constant2 /= 3.0; // for the LOD splitting of the source 

		M.thomas_constant3 += M.thomas_constant1; 
		M.thomas_constant3 += M.thomas_constant1; 
		M.thomas_constant3 += M.thomas_constant2; 

		M.thomas_constant3a += M.thomas_constant1; 
		M.thomas_constant3a += M.thomas_constant2; 

		M.thomas_cx.assign( M.mesh.x_coordinates.size() , M.thomas_constant1a ); 
		M.thomas_denomx.assign( M.mesh.x_coordinates.size() , M.thomas_constant3 ); 
		M.thomas_denomx[0] = M.thomas_constant3a; 
		M.thomas_denomx[ M.mesh.x_coordinates.size()-1 ] = M.thomas_constant3a; 
		if( M.mesh.x_coordinates.size() == 1 )
		{ M.thomas_denomx[0] = M.one; M.thomas_denomx[0] += M.thomas_constant2; } 

		M.thomas_cx[0] /= M.thomas_denomx[0]; 
		for( unsigned int i=1 ; i < M.mesh.x_coordinates.size() ; i++ )
		{
			axpy( &M.thomas_denomx[i] , M.thomas_constant1 , M.thomas_cx[i-1] ); 
			M.thomas_cx[i] /= M.thomas_denomx[i]; 
		}

		M.thomas_cy.assign( M.mesh.y_coordinates.size() , M.thomas_constant1a ); 
		M.thomas_denomy.assign( M.mesh.y_coordinates.size() , M.thomas_constant3 ); 
		M.thomas_denomy[0] = M.thomas_constant3a; 
		M.thomas_denomy[ M.mesh.y_coordinates.size()-1 ] = M.thomas_constant3a; 
		if( M.mesh.y_coordinates.size() == 1 )
		{ M.thomas_denomy[0] = M.one; M.thomas_denomy[0] += M.thomas_constant2; } 

		M.thomas_cy[0] /= M.thomas_denomy[0]; 
		for( unsigned int j=1 ; j < M.mesh.y_coordinates.size() ; j++ )
		{
			axpy( &M.thomas_denomy[j] , M.thomas_constant1 , M.thomas_cy[j-1] ); 
			M.thomas_cy[j] /= M.thomas_denomy[j]; 
		}

		M.thomas_cz.assign( M.mesh.z_coordinates.size() , M.thomas_constant1a ); 
		M.thomas_denomz.assign( M.mesh.z_coordinates.size() , M.thomas_constant3 ); 
		M.thomas_denomz[0] = M.thomas_constant3a; 
		M.thomas_denomz[ M.mesh.z_coordinates.size()-1 ] = M.thomas_constant3a; 
		if( M.mesh.z_coordinates.size() == 1 )
		{ M.thomas_denomz[0] = M.one; M.thomas_denomz[0] += M.thomas_constant2; } 

		M.thomas_cz[0] /= M.thomas_denomz[0]; 
		for( unsigned int k=1 ; k < M.mesh.z_coordinates.size() ; k++ )
		{
			axpy( &M.thomas_denomz[k] , M.thomas_constant1 , M.thomas_cz[k-1] ); 
			M.thomas_cz[k] /= M.thomas_denomz[k]; 
		}

		M.diffusion_solver_setup_done = true; 
	}

	M.apply_dirichlet_conditions_fast();

	#pragma omp parallel for 
	for( unsigned int k=0; k < M.mesh.z_coordinates.size() ; k++ )
	{
		for( unsigned int j=0; j < M.mesh.y_coordinates.size() ; j++ )
		{
			unsigned int n = M.voxel_index(0,j,k) * M.n_subs; 
			for( int d=0; d < M.n_subs; d++ )
			{ (*M.p_density_vectors)[n+d] /= M.thomas_denomx[0][d]; }

			n += M.thomas_i_jump; 
			for( unsigned int i=1; i < M.mesh.x_coordinates.size() ; i++ )
			{
				for( int d=0; d < M.n_subs; d++ )
				{ (*M.p_density_vectors)[n+d] += M.thomas_constant1[d] * (*M.p_density_vectors)[n-M.thomas_i_jump+d]; }

				for( int d=0; d < M.n_subs; d++ )
				{ (*M.p_density_vectors)[n+d] /= M.thomas_denomx[i][d]; }

				n += M.thomas_i_jump; 
			}

			if( M.mesh.x_coordinates.size() > 1 )
			{
				n = M.voxel_index( static_cast<int>( M.mesh.x_coordinates.size() )-2 , j , k ) * M.n_subs; 
				for( int i = static_cast<int>( M.mesh.x_coordinates.size() )-2 ; i >= 0 ; i-- )
				{
					for( int d=0; d < M.n_subs; d++ )
					{ (*M.p_density_vectors)[n+d] -= M.thomas_cx[i][d] * (*M.p_density_vectors)[n+M.thomas_i_jump+d]; }

					n -= M.thomas_i_jump; 
				}
			}
		}
	}

	M.apply_dirichlet_conditions_fast();

	#pragma omp parallel for 
	for( unsigned int k=0; k < M.mesh.z_coordinates.size() ; k++ )
	{
		for( unsigned int i=0; i < M.mesh.x_coordinates.size() ; i++ )
		{
			unsigned int n = M.voxel_index(i,0,k) * M.n_subs; 
			for( int d=0; d < M.n_subs; d++ )
			{ (*M.p_density_vectors)[n+d] /= M.thomas_denomy[0][d]; }

			n += M.thomas_j_jump; 
			for( unsigned int j=1; j < M.mesh.y_coordinates.size() ; j++ )
			{
				for( int d=0; d < M.n_subs; d++ )
				{ (*M.p_density_vectors)[n+d] += M.thomas_constant1[d] * (*M.p_density_vectors)[n-M.thomas_j_jump+d]; }

				for( int d=0; d < M.n_subs; d++ )
				{ (*M.p_density_vectors)[n+d] /= M.thomas_denomy[j][d]; }

				n += M.thomas_j_jump; 
			}

			if( M.mesh.y_coordinates.size() > 1 )
			{
				n = M.voxel_index( i , static_cast<int>( M.mesh.y_coordinates.size() )-2 , k ) * M.n_subs; 
				for( int j = static_cast<int>( M.mesh.y_coordinates.size() )-2 ; j >= 0 ; j-- )
				{
					for( int d=0; d < M.n_subs; d++ )
					{ (*M.p_density_vectors)[n+d] -= M.thomas_cy[j][d] * (*M.p_density_vectors)[n+M.thomas_j_jump+d]; }

					n -= M.thomas_j_jump; 
				}
			}
		}
	}

	M.apply_dirichlet_conditions_fast();

	#pragma omp parallel for 
	for( unsigned int j=0; j < M.mesh.y_coordinates.size() ; j++ )
	{
		for( unsigned int i=0; i < M.mesh.x_coordinates.size() ; i++ )
		{
			unsigned int n = M.voxel_index(i,j,0) * M.n_subs; 
			for( int d=0; d < M.n_subs; d++ )
			{ (*M.p_density_vectors)[n+d] /= M.thomas_denomz[0][d]; }

			n += M.thomas_k_jump; 
			for( unsigned int k=1; k < M.mesh.z_coordinates.size() ; k++ )
			{
				for( int d=0; d < M.n_subs; d++ )
				{ (*M.p_density_vectors)[n+d] += M.thomas_constant1[d] * (*M.p_density_vectors)[n-M.thomas_k_jump+d]; }

				for( int d=0; d < M.n_subs; d++ )
				{ (*M.p_density_vectors)[n+d] /= M.thomas_denomz[k][d]; }

				n += M.thomas_k_jump; 
			}

			if( M.mesh.z_coordinates.size() > 1 )
			{
				n = M.voxel_index( i , j , static_cast<int>( M.mesh.z_coordinates.size() )-2 ) * M.n_subs; 
				for( int k = static_cast<int>( M.mesh.z_coordinates.size() )-2 ; k >= 0 ; k-- )
				{
					for( int d=0; d < M.n_subs; d++ )
					{ (*M.p_density_vectors)[n+d] -= M.thomas_cz[k][d] * (*M.p_density_vectors)[n+M.thomas_k_jump+d]; }

					n -= M.thomas_k_jump; 
				}
			}
		}
	}

	M.apply_dirichlet_conditions_fast();

	return; 
}

void diffusion_decay_solver__constant_coefficients_LOD_3D_BLOCKING(Microenvironment &M, double dt, mpi_Environment &world, mpi_Cartesian &cart_topo)
{
        
        //std::ofstream file(M.timing_csv, std::ios::app);
        uint32_t size = world.size;
        uint32_t rank = world.rank;
        uint granurality = world.size * 2;
        if (M.granurality >= 1)
            granurality = M.granurality;
        double t_strt_set, t_end_set;
        double t_strt_x, t_end_x;
        double t_strt_y, t_end_y;
        double t_strt_z, t_end_z;

        if (M.mesh.uniform_mesh == false || M.mesh.Cartesian_mesh == false)
        {
            std::cout << "Error: This algorithm is written for uniform Cartesian meshes. Try: other solvers!" << std::endl
                      << std::endl;
            return;
        }

        // define constants and pre-computed quantities

        if (!M.diffusion_solver_setup_done)
        {
            M.mesh.x_size = M.mesh.x_coordinates.size();
            M.mesh.y_size = M.mesh.y_coordinates.size();
            M.mesh.z_size = M.mesh.z_coordinates.size();
            M.mesh.n_substrates = M.number_of_densities();
            M.n_subs = M.number_of_densities(); 

            M.thomas_i_jump = M.mesh.y_size * M.mesh.z_size * M.mesh.n_substrates;
            M.thomas_j_jump = M.mesh.z_size * M.mesh.n_substrates;
            M.thomas_k_jump = M.mesh.n_substrates;

            uint32_t y_size = M.mesh.y_coordinates.size();
            uint32_t z_size = M.mesh.z_coordinates.size();
            uint32_t p_size = M.number_of_densities(); 

            uint32_t step_size = (z_size * y_size) / granurality;

            M.snd_data_size = step_size * p_size; // Number of data elements to be sent
            M.rcv_data_size = step_size * p_size; // All p_density_vectors elements have same size, use anyone

            M.snd_data_size_last = ((z_size * y_size) % granurality) * p_size; // Number of data elements to be sent
            M.rcv_data_size_last = ((z_size * y_size) % granurality) * p_size;
            M.last_iteration = ((z_size * y_size) % granurality) > 0;

            /*-------------------------------------------------------------*/
            /* x_coordinates are of size local_x_nodes                     */
            /* (see function resize() of class Cartesian Mesh in           */
            /* BioFVM_parallel.cpp.                                        */
            /* Each line of Voxels going from left to right forms          */
            /* a tridiagonal system of Equations                           */
            /* Now these lines are going to split in the X decomposition   */
            /*-------------------------------------------------------------*/

            M.thomas_denomx.resize(M.mesh.x_coordinates.size(), M.zero); // sizeof(x_coordinates) = local_x_nodes, denomx is the main diagonal elements
            M.thomas_cx.resize(M.mesh.x_coordinates.size(), M.zero);     // Both b and c of tridiagonal matrix are equal, hence just one array needed

            /*-------------------------------------------------------------*/
            /* y_coordinates are of size of local_y_nodes.                 */
            /* Each line of Voxels going                                   */
            /* from bottom to top forms a tridiagonal system of Equations  */
            /*-------------------------------------------------------------*/

            M.thomas_denomy.resize(M.mesh.y_coordinates.size(), M.zero);
            M.thomas_cy.resize(M.mesh.y_coordinates.size(), M.zero);

            /*-------------------------------------------------------------*/
            /* z_coordinates are of size of local_z_nodes.                 */
            /* Each line of Voxels going                                   */
            /* from front to back forms a tridiagonal system of Equations  */
            /*-------------------------------------------------------------*/

            M.thomas_denomz.resize(M.mesh.z_coordinates.size(), M.zero);
            M.thomas_cz.resize(M.mesh.z_coordinates.size(), M.zero);

            /*-------------------------------------------------------------*/
            /* For X-decomposition thomas_i_jump - 1 can be in the previous*/
            /* process and thomas_i_jump+1 can be in the next processs     */
            /* hence we can use thomas_j_jump and thomas_k_jump safely     */
            /* but we CANNOT use thomas_i_jump safely                      */
            /*-------------------------------------------------------------*/

            M.thomas_i_jump = M.number_of_densities() * M.mesh.z_coordinates.size() * M.mesh.y_coordinates.size();
            M.thomas_j_jump = M.number_of_densities() * M.mesh.z_coordinates.size();
            M.thomas_k_jump = M.number_of_densities(); // M.thomas_j_jump * M.mesh.y_coordinates.size();

            /*-------------------------------------------------------------*/
            /* This part below of defining constants SHOULD typically      */
            /* not change during parallelization.                          */
            /*-------------------------------------------------------------*/

            M.thomas_constant1 = M.diffusion_coefficients; // dt*D/dx^2
            M.thomas_constant1a = M.zero;                  // -dt*D/dx^2;
            M.thomas_constant2 = M.decay_rates;            // (1/3)* dt*lambda
            M.thomas_constant3 = M.one;                    // 1 + 2*constant1 + constant2;
            M.thomas_constant3a = M.one;                   // 1 + constant1 + constant2;

            M.thomas_constant1 *= dt;
            M.thomas_constant1 /= M.mesh.dx;
            M.thomas_constant1 /= M.mesh.dx;

            M.thomas_constant1a = M.thomas_constant1;
            M.thomas_constant1a *= -1.0;

            M.thomas_constant2 *= dt;
            M.thomas_constant2 /= 3.0; // for the LOD splitting of the source, division by 3 is for 3-D

            M.thomas_constant3 += M.thomas_constant1;
            M.thomas_constant3 += M.thomas_constant1;
            M.thomas_constant3 += M.thomas_constant2;

            M.thomas_constant3a += M.thomas_constant1;
            M.thomas_constant3a += M.thomas_constant2;

            // Thomas solver coefficients

            /*--------------------------------------------------------------------*/
            /* In 1-D X decomposition, y and z-lines are contiguous and typically */
            /* the assignments below for y,z should not be changed                */
            /*--------------------------------------------------------------------*/

            M.thomas_cx.assign(M.mesh.x_coordinates.size(), M.thomas_constant1a);    // Fill b and c elements with -D * dt/dx^2
            M.thomas_denomx.assign(M.mesh.x_coordinates.size(), M.thomas_constant3); // Fill diagonal elements with (1 + 1/3 * lambda * dt + 2*D*dt/dx^2)

            if (rank == 0)
                M.thomas_denomx[0] = M.thomas_constant3a; // First diagonal element is   (1 + 1/3 * lambda * dt + 1*D*dt/dx^2)

            if (rank == (size - 1))
                M.thomas_denomx[M.mesh.x_coordinates.size() - 1] = M.thomas_constant3a; // Last diagonal element  is   (1 + 1/3 * lambda * dt + 1*D*dt/dx^2)

            if (rank == 0)
                if (M.mesh.x_coordinates.size() == 1) // This is an extreme case, won't exist, still if it does
                {                                     // then this must be at rank 0
                    M.thomas_denomx[0] = M.one;
                    M.thomas_denomx[0] += M.thomas_constant2;
                }
            if (rank == 0)
                M.thomas_cx[0] /= M.thomas_denomx[0]; // The first c element of tridiagonal matrix is div by first diagonal el.

            // axpy(1st, 2nd, 3rd) => 1st = 1st + 2nd * 3rd
            // the value at  size-1 is not actually used
            // Since value of size-1 is not used, it means it is the value after the last Diagonal element
            // cout << "Rank " << rank << endl;
            for (uint32_t ser_ctr = 0; ser_ctr <= size - 1; ser_ctr++)
            {
                if (rank == ser_ctr)
                {
                    if (rank == 0 && rank <= size - 1) // If size=1, then this process does not send data
                    {

                        for (uint32_t i = 1; i <= M.mesh.x_coordinates.size() - 1; i++)
                        {
                            axpy(&M.thomas_denomx[i], M.thomas_constant1, M.thomas_cx[i - 1]);
                            M.thomas_cx[i] /= M.thomas_denomx[i]; // the value at  size-1 is not actually used
                        }
                    }
                    else
                    {
                        for (uint32_t i = 1; i <= M.mesh.x_coordinates.size() - 1; i++)
                        {
                            axpy(&M.thomas_denomx[i], M.thomas_constant1, M.thomas_cx[i - 1]);
                            M.thomas_cx[i] /= M.thomas_denomx[i]; // the value at  size-1 is not actually used
                        }
                    }

                    if (rank < (size - 1))
                    {
                        MPI_Request send_req;
                        MPI_Isend(&(M.thomas_cx[M.mesh.x_coordinates.size() - 1][0]), M.thomas_cx[M.mesh.x_coordinates.size() - 1].size(), MPI_DOUBLE, ser_ctr + 1, 1111, cart_topo.mpi_cart_comm, &send_req);
                    }
                }

                if (rank == (ser_ctr + 1) && (ser_ctr + 1) <= (size - 1))
                {

                    std::vector<double> temp_cx(M.thomas_cx[0].size());
                    MPI_Request recv_req;
                    MPI_Irecv(&temp_cx[0], temp_cx.size(), MPI_DOUBLE, ser_ctr, 1111, cart_topo.mpi_cart_comm, &recv_req);
                    MPI_Wait(&recv_req, MPI_STATUS_IGNORE);

                    axpy(&M.thomas_denomx[0], M.thomas_constant1, temp_cx); // CHECK IF &temp_cz[0] is OK, axpy() in BioFVM_vector.cpp
                    M.thomas_cx[0] /= M.thomas_denomx[0];                   // the value at  size-1 is not actually used
                }

                MPI_Barrier(cart_topo.mpi_cart_comm);
            }

            /*--------------------------------------------------------------------*/
            /* In 1-D X decomposition, z and y-lines are contiguous adn typically */
            /* the assignments below for z,y should not be changed                */
            /* Both the first voxel i.e. index 0 and last voxel i.e. index=       */
            /* y_coordinates.size()-1 are on the same process                     */
            /*--------------------------------------------------------------------*/

            M.thomas_cy.assign(M.mesh.y_coordinates.size(), M.thomas_constant1a);
            M.thomas_denomy.assign(M.mesh.y_coordinates.size(), M.thomas_constant3);
            M.thomas_denomy[0] = M.thomas_constant3a;
            M.thomas_denomy[M.mesh.y_coordinates.size() - 1] = M.thomas_constant3a;
            if (M.mesh.y_coordinates.size() == 1)
            {
                M.thomas_denomy[0] = M.one;
                M.thomas_denomy[0] += M.thomas_constant2;
            }
            M.thomas_cy[0] /= M.thomas_denomy[0];
            for (uint32_t i = 1; i <= M.mesh.y_coordinates.size() - 1; i++)
            {
                axpy(&M.thomas_denomy[i], M.thomas_constant1, M.thomas_cy[i - 1]);
                M.thomas_cy[i] /= M.thomas_denomy[i]; // the value at  size-1 is not actually used
            }

            M.thomas_cz.assign(M.mesh.z_coordinates.size(), M.thomas_constant1a);
            M.thomas_denomz.assign(M.mesh.z_coordinates.size(), M.thomas_constant3);
            M.thomas_denomz[0] = M.thomas_constant3a;
            M.thomas_denomz[M.mesh.z_coordinates.size() - 1] = M.thomas_constant3a;
            if (M.mesh.z_coordinates.size() == 1)
            {
                M.thomas_denomz[0] = M.one;
                M.thomas_denomz[0] += M.thomas_constant2;
            }
            M.thomas_cz[0] /= M.thomas_denomz[0];
            for (uint32_t i = 1; i <= M.mesh.z_coordinates.size() - 1; i++)
            {
                axpy(&M.thomas_denomz[i], M.thomas_constant1, M.thomas_cz[i - 1]);
                M.thomas_cz[i] /= M.thomas_denomz[i]; // the value at  size-1 is not actually used
            }

            if (world.rank == 0) { //Print message size and number of messages
                std::cout << "  Diffusion messages:" << granurality * (world.size-1) * 2;
                std::cout << " | Size in " << ((double)(M.thomas_i_jump * 8))/(1024.0*1024.0) <<" MB" << std::endl;   
            }
            M.diffusion_solver_setup_done = true;
            //if (rank == 0) file << "X-diffusion,Y-diffusion,Z-diffusion,Apply Dirichlet" << std::endl;
        }

        uint32_t n_req = granurality;
        if (M.last_iteration > 0) n_req+=1;
        MPI_Request send_req[n_req];
        MPI_Request recv_req[n_req];
        std::vector<double> block3d(M.thomas_i_jump); //The message to send is of the size Y_voxels * Z_voxels * Substrates

        M.apply_dirichlet_boundaries_conditions(rank, size);

        /*-----------------------------------------------------------------------------------*/
        /*                        FORWARD ELIMINATION - x DIRECTION/DECOMPOSITION            */
        /*-----------------------------------------------------------------------------------*/

        
        if (rank == 0)
        {
            for (uint32_t step = 0; step < granurality; ++step)
            {
                uint32_t initial_index = step * M.snd_data_size;
                #pragma omp parallel for
                for (uint32_t index = initial_index; index < initial_index + M.snd_data_size; index += M.mesh.n_substrates)
                {
                    uint32_t index_dec = index; 
                    for (uint32_t d = 0; d < M.thomas_k_jump; d++)
                    {
                        (*M.p_density_vectors)[index + d] /= M.thomas_denomx[0][d];
                    }

                    for (uint32_t i = 1; i < M.mesh.x_size; i++)
                    {
                        
                        uint32_t index_inc = index_dec + M.thomas_i_jump;
                        // axpy(&(*M.p_density_vectors)[n], M.thomas_constant1, (*M.p_density_vectors)[n - M.thomas_i_jump]);
                        for (uint32_t d = 0; d < M.thomas_k_jump; d++)
                        {
                            (*M.p_density_vectors)[index_inc + d] += M.thomas_constant1[d] * (*M.p_density_vectors)[index_dec + d];
                        }

                        //(*M.p_density_vectors)[n] /= M.thomas_denomx[i];
                        for (uint32_t d = 0; d < M.thomas_k_jump; d++)
                        {
                            (*M.p_density_vectors)[index_inc + d] /= M.thomas_denomx[i][d];
                        }
                        index_dec = index_inc;
                    }
                }

                if (size > 1) {
                    uint32_t x_end = M.mesh.x_size - 1;
                    uint32_t offset = step * M.snd_data_size;
                    MPI_Isend(&((*M.p_density_vectors)[x_end * M.thomas_i_jump + offset]), M.snd_data_size, MPI_DOUBLE, rank + 1, step, cart_topo.mpi_cart_comm, &send_req[step]);
                }
            }
            //Last iteration
            if (M.snd_data_size_last > 0) {
                uint32_t initial_index = granurality * M.snd_data_size;
                #pragma omp parallel for
                for (uint32_t index = initial_index; index < initial_index + M.snd_data_size_last; index += M.thomas_k_jump)
                {
                    uint32_t index_dec = index; 
                    for (uint32_t d = 0; d < M.thomas_k_jump; d++)
                    {
                        (*M.p_density_vectors)[index + d] /= M.thomas_denomx[0][d];
                    }

                    for (uint32_t i = 1; i < M.mesh.x_size; i++)
                    {
                        uint32_t index_inc = index_dec + M.thomas_i_jump;
                        // axpy(&(*(*M.p_density_vectors))[n], M.thomas_constant1, (*(*M.p_density_vectors))[n - M.thomas_i_jump]);
                        for (uint32_t d = 0; d < M.thomas_k_jump; d++)
                        {
                            (*M.p_density_vectors)[index_inc + d] += M.thomas_constant1[d] * (*M.p_density_vectors)[index_dec + d];
                        }
                        //(*(*M.p_density_vectors))[n] /= M.thomas_denomx[i];
                        for (uint32_t d = 0; d < M.thomas_k_jump; d++)
                        {
                            (*M.p_density_vectors)[index_inc + d] /= M.thomas_denomx[i][d];
                        }
                        index_dec = index_inc;
                    }
                }

                if (size > 1) {
                    uint32_t x_end = M.mesh.x_size - 1;
                    uint32_t offset = granurality * M.snd_data_size;
                    MPI_Isend(&((*M.p_density_vectors)[x_end * M.thomas_i_jump + offset]), M.snd_data_size_last, MPI_DOUBLE, rank + 1, granurality, cart_topo.mpi_cart_comm, &send_req[granurality]);
                    
                }
            }
        }
        else
        {
            if (rank >= 1 && rank <= (size - 1))
            {
                for (uint32_t step = 0; step < granurality; ++step)
                {
                    uint32_t initial_index = step * M.snd_data_size;
                    MPI_Irecv(&(block3d[initial_index]), M.rcv_data_size, MPI_DOUBLE, rank-1, step, cart_topo.mpi_cart_comm, &(recv_req[step]));
                }
                if (M.last_iteration)
                    MPI_Irecv(&(block3d[granurality*M.snd_data_size]), M.rcv_data_size_last, MPI_DOUBLE, rank-1, granurality, cart_topo.mpi_cart_comm, &(recv_req[granurality]));
                
                for (uint32_t step = 0; step < granurality; ++step)
                {
                    uint32_t initial_index = step * M.snd_data_size;
                    //MPI_Wait(&recv_req[step], MPI_STATUS_IGNORE);
                    MPI_Status status;
                    MPI_Wait(&recv_req[step], &status);

                    int received_count = 0;
                    MPI_Get_count(&status, MPI_DOUBLE, &received_count);

                    if (received_count != (int)M.rcv_data_size) {
                        std::fprintf(stderr, "RANK %d: TRUNCATION? step=%u expected=%zu received=%d source=%d tag=%d\n",
                                    rank, step, M.rcv_data_size, received_count, status.MPI_SOURCE, status.MPI_TAG);
                    }

                    #pragma omp parallel for
                    for (uint32_t index = initial_index; index < initial_index + M.snd_data_size; index += M.thomas_k_jump)
                    {
                        // axpy(&(*(*M.p_density_vectors))[n], M.thomas_constant1, block3d[k][j]);
                        uint32_t index_dec = index;
                        for (uint32_t d = 0; d < M.thomas_k_jump; d++)
                        {
                            (*M.p_density_vectors)[index + d] += M.thomas_constant1[d] * block3d[index + d];
                        }
                        //(*(*M.p_density_vectors))[n] /= M.thomas_denomx[0];
                        for (uint32_t d = 0; d < M.mesh.n_substrates; d++)
                        {
                            (*M.p_density_vectors)[index + d] /= M.thomas_denomx[0][d];
                        }
                        for (uint32_t i = 1; i < M.mesh.x_size; i++)
                        {
                            uint32_t index_inc = index_dec + M.thomas_i_jump;
                            // axpy(&(*(*M.p_density_vectors))[n], M.thomas_constant1, (*(*M.p_density_vectors))[n - M.thomas_i_jump]);
                            for (uint32_t d = 0; d < M.thomas_k_jump; d++)
                            {
                                (*M.p_density_vectors)[index_inc + d] += M.thomas_constant1[d] * (*M.p_density_vectors)[index_dec + d];
                            }
                            //(*(*M.p_density_vectors))[n] /= M.thomas_denomx[i];
                            for (uint32_t d = 0; d < M.thomas_k_jump; d++)
                            {
                                (*M.p_density_vectors)[index_inc + d] /= M.thomas_denomx[i][d];
                            }
                            index_dec = index_inc;
                        }
                    }
                    if (rank < (size - 1))
                    {
                        uint32_t x_end = M.mesh.x_size - 1;
                        MPI_Isend(&((*M.p_density_vectors)[(x_end * M.thomas_i_jump) + initial_index]), M.snd_data_size, MPI_DOUBLE, rank + 1, step, cart_topo.mpi_cart_comm, &send_req[step]);
                    }
                }
                if (M.snd_data_size_last > 0)
                {
                    uint32_t initial_index = granurality * M.snd_data_size;
                    MPI_Wait(&recv_req[granurality], MPI_STATUS_IGNORE); 
                    #pragma omp parallel for
                    for (uint32_t index = initial_index; index < initial_index + M.snd_data_size_last; index += M.thomas_k_jump)
                    {
                        // axpy(&(*(*M.p_density_vectors))[n], M.thomas_constant1, block3d[k][j]);
                        uint32_t index_dec = index;
                        for (uint32_t d = 0; d < M.thomas_k_jump; d++)
                        {
                            (*M.p_density_vectors)[index + d] += M.thomas_constant1[d] * block3d[index + d];
                        }
                        //(*(*M.p_density_vectors))[n] /= M.thomas_denomx[0];
                        for (uint32_t d = 0; d < M.thomas_k_jump; d++)
                        {
                            (*M.p_density_vectors)[index + d] /= M.thomas_denomx[0][d];
                        }

                        for (uint32_t i = 1; i < M.mesh.x_size; i++)
                        {
                            uint32_t index_inc = index_dec + M.thomas_i_jump;
                            // axpy(&(*(*M.p_density_vectors))[n], M.thomas_constant1, (*(*M.p_density_vectors))[n - M.thomas_i_jump]);
                            for (uint32_t d = 0; d < M.thomas_k_jump; d++)
                            {
                                (*M.p_density_vectors)[index_inc + d] += M.thomas_constant1[d] * (*M.p_density_vectors)[index_dec + d];
                            }
                            //(*(*M.p_density_vectors))[n] /= M.thomas_denomx[i];
                            for (uint32_t d = 0; d < M.thomas_k_jump; d++)
                            {
                                (*M.p_density_vectors)[index_inc + d] /= M.thomas_denomx[i][d];
                            }

                            index_dec = index_inc;
                        }
                        
                    }
                    // End of computation region
                    if (rank < (size - 1))
                    {
                        uint32_t x_end = M.mesh.x_size - 1;
                        MPI_Request aux;
                        MPI_Isend(&((*M.p_density_vectors)[x_end * M.thomas_i_jump + initial_index]), M.snd_data_size_last, MPI_DOUBLE, rank + 1, granurality, cart_topo.mpi_cart_comm, &send_req[granurality]);
                      
                    }
                }
            }
        }
        MPI_Barrier(cart_topo.mpi_cart_comm);
        /*-----------------------------------------------------------------------------------*/
        /*                         CODE FOR BACK SUBSITUTION                                 */
        /*-----------------------------------------------------------------------------------*/
        if (rank == (size - 1))
        {
            for (uint32_t step = 0; step < granurality; ++step)
            {
                uint32_t initial_index = ((M.mesh.x_size - 1)*M.thomas_i_jump) + (step * M.snd_data_size);
                #pragma omp parallel for
                for (uint32_t index = initial_index; index < initial_index + M.snd_data_size; index += M.mesh.n_substrates)
                {
                    uint32_t index_aux = index;
                    for (int32_t i = M.mesh.x_size - 2; i >= 0; i--)
                    {

                        uint32_t index_dec = index_aux - M.thomas_i_jump;
                        // naxpy(&(*(*M.p_density_vectors))[n], M.thomas_cx[i], (*(*M.p_density_vectors))[n + M.thomas_i_jump]);
                        for (uint32_t d = 0; d < M.thomas_k_jump; d++)
                        {
                            (*M.p_density_vectors)[index_dec + d] -= M.thomas_cx[i][d] * (*M.p_density_vectors)[index_aux + d];
                        }
                        index_aux = index_dec;
                    }
                }
                if (size > 1) {
                    MPI_Isend(&((*M.p_density_vectors)[step * M.snd_data_size]), M.snd_data_size, MPI_DOUBLE, rank - 1, step + granurality, cart_topo.mpi_cart_comm, &send_req[step]);
                }
            }

            //Last iteration
            if (M.snd_data_size_last > 0) {
                uint32_t initial_index = ((M.mesh.x_size - 1)*M.thomas_i_jump) + (granurality * M.snd_data_size);
                #pragma omp parallel for
                for (uint32_t index = initial_index; index < initial_index + M.snd_data_size_last; index += M.mesh.n_substrates)
                {
                    uint32_t index_aux = index;
                    for (int32_t i = M.mesh.x_coordinates.size() - 2; i >= 0; i--)
                    {

                        uint32_t index_dec = index_aux - M.thomas_i_jump;
                        // naxpy(&(*(*M.p_density_vectors))[n], M.thomas_cx[i], (*(*M.p_density_vectors))[n + M.thomas_i_jump]);
                        for (uint32_t d = 0; d < M.thomas_k_jump; d++)
                        {
                            (*M.p_density_vectors)[index_dec + d] -= M.thomas_cx[i][d] * (*M.p_density_vectors)[index_aux + d];
                        }
                        index_aux = index_dec;
                    }
                }
                if (size > 1) {
                    MPI_Request aux;
                    MPI_Isend(&((*M.p_density_vectors)[granurality * M.snd_data_size]), M.snd_data_size_last, MPI_DOUBLE, rank - 1, granurality + granurality, cart_topo.mpi_cart_comm, &send_req[granurality]);
                }
            
            }
        }
        else
        {
            for (uint32_t step = 0; step < granurality; ++step) {
                MPI_Irecv(&(block3d[step*M.snd_data_size]), M.rcv_data_size, MPI_DOUBLE, rank+1, step + granurality, cart_topo.mpi_cart_comm, &recv_req[step]);}
            if (M.last_iteration)
                MPI_Irecv(&(block3d[granurality*M.snd_data_size]), M.rcv_data_size_last, MPI_DOUBLE, rank+1, granurality + granurality, cart_topo.mpi_cart_comm, &recv_req[granurality]);
            
            for (uint32_t step = 0; step < granurality; ++step)
            {
                uint32_t initial_index = ((M.mesh.x_size - 1)*M.thomas_i_jump) + (step * M.snd_data_size);
                uint32_t index_3d_initial = (step * M.snd_data_size);
                //MPI_Wait(&recv_req[step], MPI_STATUS_IGNORE);
                MPI_Status status;
                MPI_Wait(&recv_req[step], &status);

                int received_count = 0;
                MPI_Get_count(&status, MPI_DOUBLE, &received_count);

                if (received_count != (int)M.rcv_data_size) {
                    std::fprintf(stderr, "RANK %d: TRUNCATION? step=%u expected=%zu received=%d source=%d tag=%d\n",
                                rank, step, M.rcv_data_size, received_count, status.MPI_SOURCE, status.MPI_TAG);
                }
                #pragma omp parallel for
                for (uint32_t offset = 0; offset < M.snd_data_size; offset += M.mesh.n_substrates)
                {
                    uint32_t index_aux = initial_index + offset;
                    uint32_t index_3d = index_3d_initial + offset;
                    for (uint32_t d = 0; d < M.thomas_k_jump; d++)
                    {
                        (*M.p_density_vectors)[index_aux + d] -= M.thomas_cx[M.mesh.x_size - 1][d] * block3d[index_3d + d];
                    }

                    for (int32_t i = M.mesh.x_size - 2; i >= 0; i--)
                    {

                        uint32_t index_dec = index_aux - M.thomas_i_jump;
                        // naxpy(&(*(*M.p_density_vectors))[n], M.thomas_cx[i], (*(*M.p_density_vectors))[n + M.thomas_i_jump]);
                        for (uint32_t d = 0; d < M.thomas_k_jump; d++)
                        {
                            (*M.p_density_vectors)[index_dec + d] -= M.thomas_cx[i][d] * (*M.p_density_vectors)[index_aux + d];
                        }
                        index_aux = index_dec;
                        
                    }
                }
                if (rank > 0)
                {
                    MPI_Request aux;
                    //MPI_Isend(&((*M.p_density_vectors)[step * M.snd_data_size]), M.snd_data_size, MPI_DOUBLE, rank - 1, step, cart_topo.mpi_cart_comm, &send_req[step]);
                    MPI_Send(&((*M.p_density_vectors)[step * M.snd_data_size]), M.snd_data_size, MPI_DOUBLE, rank - 1, step + granurality, cart_topo.mpi_cart_comm);
                }
            }
            if (M.snd_data_size_last > 0)
            {
                uint32_t initial_index = ((M.mesh.x_size - 1)*M.thomas_i_jump) + (granurality * M.snd_data_size);
                uint32_t index_3d_initial = (granurality * M.snd_data_size);
                MPI_Wait(&recv_req[granurality], MPI_STATUS_IGNORE);
                #pragma omp parallel for
                for (uint32_t offset = 0; offset < M.snd_data_size_last; offset += M.mesh.n_substrates)
                {
                    uint32_t index_aux = initial_index + offset;
                    uint32_t index_3d = index_3d_initial + offset;
                    //naxpy(&(*(*M.p_density_vectors))[n], M.thomas_cx[M.mesh.x_coordinates.size() - 1], block3d[k][j]);
                    for (uint32_t d = 0; d < M.thomas_k_jump; d++)
                    {
                        (*M.p_density_vectors)[index_aux + d] -= M.thomas_cx[M.mesh.x_size - 1][d] * block3d[index_3d + d];
                    }
                    for (int32_t i = M.mesh.x_size - 2; i >= 0; i--)
                    {
                        uint32_t index_dec = index_aux - M.thomas_i_jump;
                        // naxpy(&(*(*M.p_density_vectors))[n], M.thomas_cx[i], (*(*M.p_density_vectors))[n + M.thomas_i_jump]);
                        for (uint32_t d = 0; d < M.thomas_k_jump; d++)
                        {
                            (*M.p_density_vectors)[index_dec + d] -= M.thomas_cx[i][d] * (*M.p_density_vectors)[index_aux + d];
                        }
                        index_aux = index_dec;
                    }
                }
                if (rank > 0)
                {
                    MPI_Request aux;
                    MPI_Isend(&((*M.p_density_vectors)[granurality * M.snd_data_size]), M.snd_data_size_last, MPI_DOUBLE, rank - 1, granurality + granurality, cart_topo.mpi_cart_comm, &send_req[granurality]);
                }
            }
        }
        MPI_Barrier(cart_topo.mpi_cart_comm);

        M.apply_dirichlet_boundaries_conditions(rank, size);

        #pragma omp parallel for collapse(2)
        for (uint32_t i = 0; i < M.mesh.x_size; i++)
        {
            for (uint32_t k = 0; k < M.mesh.z_size; k++)
            {

                uint32_t index = i * M.thomas_i_jump + k * M.thomas_k_jump;
                //(*(*M.p_density_vectors))[n] /= M.thomas_denomy[0];
                for (uint32_t d = 0; d < M.thomas_k_jump; d++)
                {
                    (*M.p_density_vectors)[index + d] /= M.thomas_denomy[0][d];
                }

                for (uint32_t j = 1; j < M.mesh.y_size; j++)
                {

                    uint32_t index_inc = index + M.thomas_j_jump;
                    // axpy(&(*(*M.p_density_vectors))[n], M.thomas_constant1, (*(*M.p_density_vectors))[n - M.thomas_j_jump]);
                    for (uint32_t d = 0; d < M.thomas_k_jump; d++)
                    {
                        (*M.p_density_vectors)[index_inc + d] += M.thomas_constant1[d] * (*M.p_density_vectors)[index + d];
                    }
                    //(*(*M.p_density_vectors))[n] /= M.thomas_denomy[j];
                    for (uint32_t d = 0; d < M.thomas_k_jump; d++)
                    {
                        (*M.p_density_vectors)[index_inc + d] /= M.thomas_denomy[j][d];
                    }
                    index = index_inc;
                }

                // back substitution

                index = i * M.thomas_i_jump + k * M.thomas_k_jump + (M.thomas_j_jump * (M.mesh.y_size - 1));
                for (int32_t j = M.mesh.y_size - 2; j >= 0; j--)
                {

                    uint32_t index_dec = index - M.thomas_j_jump;
                    // naxpy(&(*(*M.p_density_vectors))[n], M.thomas_cy[j], (*(*M.p_density_vectors))[n + M.thomas_j_jump]);
                    for (uint32_t d = 0; d < M.thomas_k_jump; d++)
                    {
                        (*M.p_density_vectors)[index_dec + d] -= M.thomas_cy[j][d] * (*M.p_density_vectors)[index + d];
                    }
                    index = index_dec;
                }
            }
        }

        M.apply_dirichlet_boundaries_conditions(rank, size);

        #pragma omp parallel for collapse(2)
        for (uint32_t i = 0; i < M.mesh.x_size; i++)
        {
           for (uint32_t j = 0; j < M.mesh.y_size; j++)
            {

                uint32_t index = i * M.thomas_i_jump + j * M.thomas_j_jump;
                //(*(*M.p_density_vectors))[n] /= M.thomas_denomz[0];
                for (uint32_t d = 0; d < M.thomas_k_jump; d++)
                {
                    (*M.p_density_vectors)[index + d] /= M.thomas_denomz[0][d];
                }

                // should be an empty loop if mesh.z_coordinates.size() < 2
                for (uint32_t k = 1; k < M.mesh.z_size; k++)
                {
                    uint32_t index_inc = index + M.thomas_k_jump;
                    // axpy(&(*(*M.p_density_vectors))[n], M.thomas_constant1, (*(*M.p_density_vectors))[n - M.thomas_k_jump]);
                    for (uint32_t d = 0; d < M.thomas_k_jump; d++)
                    {
                        (*M.p_density_vectors)[index_inc + d] += M.thomas_constant1[d] * (*M.p_density_vectors)[index + d];
                    }
                    //(*(*M.p_density_vectors))[n] /= M.thomas_denomz[k];
                    for (uint32_t d = 0; d < M.thomas_k_jump; d++)
                    {
                        (*M.p_density_vectors)[index_inc + d] /= M.thomas_denomz[k][d];
                    }

                    index = index_inc;
                }
                // for parallelization need to break forward elimination and back substitution uint32_to
                // should be an empty loop if mesh.z_coordinates.size() < 2
                index = i * M.thomas_i_jump + j * M.thomas_j_jump + (M.thomas_k_jump * (M.mesh.z_size - 1));
                for (int32_t k = M.mesh.z_coordinates.size() - 2; k >= 0; k--)
                {
                    uint32_t index_dec = index - M.thomas_k_jump;
                    // naxpy(&(*(*M.p_density_vectors))[n], M.thomas_cz[k], (*(*M.p_density_vectors))[n + M.thomas_k_jump]);
                    for (uint32_t d = 0; d < M.thomas_k_jump; d++)
                    {
                        (*M.p_density_vectors)[index_dec + d] -= M.thomas_cz[k][d] * (*M.p_density_vectors)[index + d];
                    }
                    index = index_dec;
                }
            }
        }
        M.apply_dirichlet_boundaries_conditions(rank, size);
        return;
    }

    int gcd(int a, int b) {
        while (b != 0) {
            int temp = b;
            b = a % b;
            a = temp;
        }
        return a;
    }

    // Function to compute the least common multiple (LCM)
    int lcm(int a, int b) {
        return (a * b) / gcd(a, b);
    }

    void diffusion_decay_solver__constant_coefficients_LOD_3D_AVX256D(Microenvironment &M, double dt, mpi_Environment &world, mpi_Cartesian &cart_topo){
        int size = world.size;
        int rank = world.rank;
        uint granurality = 1;
        if (M.granurality >= 1)
            granurality = M.granurality;
        MPI_Request send_req[granurality], recv_req[granurality];
        //std::ofstream file(M.timing_csv, std::ios::app);
        int vl = 4;

        if (M.diffusion_solver_setup_done == false) {
            
            MPI_Request send_req[size];
            MPI_Request recv_req[size];
            
            M.mesh.x_size = M.mesh.x_coordinates.size();
            M.mesh.y_size = M.mesh.y_coordinates.size();
            M.mesh.z_size = M.mesh.z_coordinates.size();
            M.mesh.n_substrates = M.number_of_densities();
            M.n_subs = M.number_of_densities(); 

            M.thomas_i_jump = M.mesh.y_size * M.mesh.z_size * M.mesh.n_substrates;
            M.thomas_j_jump = M.mesh.z_size * M.mesh.n_substrates;
            M.thomas_k_jump = M.mesh.n_substrates;


            std::vector<double> zero(M.mesh.n_substrates, 0.0);
            std::vector<double> one(M.mesh.n_substrates, 1.0);

            int step_size = (M.mesh.z_size * M.mesh.y_size) / granurality;
            


            M.snd_data_size = step_size * M.mesh.n_substrates; // Number of data elements to be sent
            M.rcv_data_size = step_size * M.mesh.n_substrates; // All p_density_vectors elements have same size, use anyone

            M.snd_data_size_last = ((M.mesh.z_size * M.mesh.y_size) % granurality) * M.mesh.n_substrates; // Number of data elements to be sent
            M.rcv_data_size_last = ((M.mesh.z_size * M.mesh.y_size) % granurality) * M.mesh.n_substrates;

            //Thomas initialization
            M.thomas_denomx.resize(M.mesh.x_size, zero); // sizeof(x_coordinates) = local_x_nodes, denomx is the main diagonal elements
            M.thomas_cx.resize(M.mesh.x_size, zero);     // Both b and c of tridiagonal matrix are equal, hence just one array needed

            /*-------------------------------------------------------------*/
            /* y_coordinates are of size of local_y_nodes.                 */
            /* Each line of Voxels going                                   */
            /* from bottom to top forms a tridiagonal system of Equations  */
            /*-------------------------------------------------------------*/

            M.thomas_denomy.resize(M.mesh.y_size, zero);
            M.thomas_cy.resize(M.mesh.y_size, zero);

            /*-------------------------------------------------------------*/
            /* z_coordinates are of size of local_z_nodes.                 */
            /* Each line of Voxels going                                   */
            /* from front to back forms a tridiagonal system of Equations  */
            /*-------------------------------------------------------------*/

            M.thomas_denomz.resize(M.mesh.z_size, zero);
            M.thomas_cz.resize(M.mesh.z_size, zero);

            /*-------------------------------------------------------------*/
            /* For X-decomposition thomas_i_jump - 1 can be in the previous*/
            /* process and thomas_i_jump+1 can be in the next processs     */
            /* hence we can use thomas_j_jump and thomas_k_jump safely     */
            /* but we CANNOT use thomas_i_jump safely                      */
            /*-------------------------------------------------------------*/

            int i_jump = M.mesh.y_size*M.mesh.z_size*M.mesh.n_substrates;
            int j_jump = M.mesh.z_size*M.mesh.n_substrates;
            int k_jump = M.mesh.n_substrates; 

            /*-------------------------------------------------------------*/
            /* This part below of defining constants SHOULD typically      */
            /* not change during parallelization.                          */
            /*-------------------------------------------------------------*/

            M.thomas_constant1 = M.diffusion_coefficients; // dt*D/dx^2
            M.thomas_constant1a = M.zero;                  // -dt*D/dx^2;
            M.thomas_constant2 = M.decay_rates;            // (1/3)* dt*lambda
            M.thomas_constant3 = M.one;                    // 1 + 2*constant1 + constant2;
            M.thomas_constant3a = M.one;                   // 1 + constant1 + constant2;

            M.thomas_constant1 *= dt;
            M.thomas_constant1 /= M.mesh.dx; //dx
            M.thomas_constant1 /= M.mesh.dx; //dx

            M.thomas_constant1a = M.thomas_constant1;
            M.thomas_constant1a *= -1.0;

            M.thomas_constant2 *= dt;
            M.thomas_constant2 /= 3.0; // for the LOD splitting of the source, division by 3 is for 3-D

            M.thomas_constant3 += M.thomas_constant1;
            M.thomas_constant3 += M.thomas_constant1;
            M.thomas_constant3 += M.thomas_constant2;

            M.thomas_constant3a += M.thomas_constant1;
            M.thomas_constant3a += M.thomas_constant2;

            // Thomas solver coefficients

            /*--------------------------------------------------------------------*/
            /* In 1-D X decomposition, y and z-lines are contiguous and typically */
            /* the assignments below for y,z should not be changed                */
            /*--------------------------------------------------------------------*/

            M.thomas_cx.assign(M.mesh.x_size, M.thomas_constant1a);    // Fill b and c elements with -D * dt/dx^2
            M.thomas_denomx.assign(M.mesh.x_size, M.thomas_constant3); // Fill diagonal elements with (1 + 1/3 * lambda * dt + 2*D*dt/dx^2)

            if (rank == 0)
                M.thomas_denomx[0] = M.thomas_constant3a; // First diagonal element is   (1 + 1/3 * lambda * dt + 1*D*dt/dx^2)

            if (rank == (size - 1))
                M.thomas_denomx[M.mesh.x_size - 1] = M.thomas_constant3a; // Last diagonal element  is   (1 + 1/3 * lambda * dt + 1*D*dt/dx^2)

            if (rank == 0)
                if (M.mesh.x_size == 1) // This is an extreme case, won't exist, still if it does
                {                                     // then this must be at rank 0
                    M.thomas_denomx[0] = M.one;
                    M.thomas_denomx[0] += M.thomas_constant2;
                }
            if (rank == 0)
                M.thomas_cx[0] /= M.thomas_denomx[0]; // The first c element of tridiagonal matrix is div by first diagonal el.

            // axpy(1st, 2nd, 3rd) => 1st = 1st + 2nd * 3rd
            // the value at  size-1 is not actually used
            // Since value of size-1 is not used, it means it is the value after the last Diagonal element
            // cout << "Rank " << rank << endl;
            for (int ser_ctr = 0; ser_ctr <= size - 1; ser_ctr++)
            {
                if (rank == ser_ctr)
                {
                    if (rank == 0 && rank <= size - 1) // If size=1, then this process does not send data
                    {

                        for (int i = 1; i <= M.mesh.x_size - 1; i++)
                        {
                            axpy(&M.thomas_denomx[i], M.thomas_constant1, M.thomas_cx[i - 1]);
                            M.thomas_cx[i] /= M.thomas_denomx[i]; // the value at  size-1 is not actually used
                        }
                    }
                    else
                    {
                        for (int i = 1; i <= M.mesh.x_size - 1; i++)
                        {
                            axpy(&M.thomas_denomx[i], M.thomas_constant1, M.thomas_cx[i - 1]);
                            M.thomas_cx[i] /= M.thomas_denomx[i]; // the value at  size-1 is not actually used
                        }
                    }

                    if (rank < (size - 1))
                    {
                        MPI_Isend(&(M.thomas_cx[M.mesh.x_size - 1][0]), M.thomas_cx[M.mesh.x_size - 1].size(), MPI_DOUBLE, ser_ctr + 1, 1111, cart_topo.mpi_cart_comm, &send_req[0]);
                    }
                }

                if (rank == (ser_ctr + 1) && (ser_ctr + 1) <= (size - 1))
                {

                    std::vector<double> temp_cx(M.thomas_cx[0].size());

                    MPI_Irecv(&temp_cx[0], temp_cx.size(), MPI_DOUBLE, ser_ctr, 1111, cart_topo.mpi_cart_comm, &recv_req[0]);
                    MPI_Wait(&recv_req[0], MPI_STATUS_IGNORE);

                    axpy(&M.thomas_denomx[0], M.thomas_constant1, temp_cx); // CHECK IF &temp_cz[0] is OK, axpy() in BioFVM_vector.cpp
                    M.thomas_cx[0] /= M.thomas_denomx[0];                   // the value at  size-1 is not actually used
                }

                MPI_Barrier(cart_topo.mpi_cart_comm);
            }

            /*--------------------------------------------------------------------*/
            /* In 1-D X decomposition, z and y-lines are contiguous adn typically */
            /* the assignments below for z,y should not be changed                */
            /* Both the first voxel i.e. index 0 and last voxel i.e. index=       */
            /* y_coordinates.size()-1 are on the same process                     */
            /*--------------------------------------------------------------------*/

            M.thomas_cy.assign(M.mesh.y_size, M.thomas_constant1a);
            M.thomas_denomy.assign(M.mesh.y_size, M.thomas_constant3);
            M.thomas_denomy[0] = M.thomas_constant3a;
            M.thomas_denomy[M.mesh.y_size - 1] = M.thomas_constant3a;
            if (M.mesh.y_size == 1)
            {
                M.thomas_denomy[0] = M.one;
                M.thomas_denomy[0] += M.thomas_constant2;
            }
            M.thomas_cy[0] /= M.thomas_denomy[0];
            for (int i = 1; i <= M.mesh.y_size - 1; i++)
            {
                axpy(&M.thomas_denomy[i], M.thomas_constant1, M.thomas_cy[i - 1]);
                M.thomas_cy[i] /= M.thomas_denomy[i]; // the value at  size-1 is not actually used
            }

            M.thomas_cz.assign(M.mesh.z_size, M.thomas_constant1a);
            M.thomas_denomz.assign(M.mesh.z_size, M.thomas_constant3);
            M.thomas_denomz[0] = M.thomas_constant3a;
            M.thomas_denomz[M.mesh.z_size - 1] = M.thomas_constant3a;
            if (M.mesh.z_size == 1)
            {
                M.thomas_denomz[0] = M.one;
                M.thomas_denomz[0] += M.thomas_constant2;
            }
            M.thomas_cz[0] /= M.thomas_denomz[0];
            for (int i = 1; i <= M.mesh.z_size - 1; i++)
            {
                axpy(&M.thomas_denomz[i], M.thomas_constant1, M.thomas_cz[i - 1]);
                M.thomas_cz[i] /= M.thomas_denomz[i]; // the value at  size-1 is not actually used
            }
            M.diffusion_solver_setup_done = true;
            //if (rank == 0) file << "X-diffusion,Y-diffusion,Z-diffusion,Apply Dirichlet" << std::endl;
        }
        if (M.diffusion_solver_vectorized_setup_done == false) {
            //Vectorization initialization
            M.gvec_size = lcm(M.mesh.n_substrates, vl);

            //X-diffusion
            M.gthomas_constant1.resize(M.gvec_size, 0.0);
            auto dest_iter =  M.gthomas_constant1.begin();
            for (int j = 0; j < M.gvec_size; j+=M.mesh.n_substrates){
                copy(M.thomas_constant1.begin(), M.thomas_constant1.end(), dest_iter);
                dest_iter+=M.mesh.n_substrates;
            }

            M.gthomas_denomx.resize(M.mesh.x_size);
            M.gthomas_cx.resize(M.mesh.x_size);
            for (int i = 0; i < M.mesh.x_size; ++i){
                M.gthomas_denomx[i].resize(M.gvec_size, 0.0);
                M.gthomas_cx[i].resize(M.gvec_size, 0.0);
                auto dest_denomx = M.gthomas_denomx[i].begin();
                auto dest_cx = M.gthomas_cx[i].begin();
                for (int d = 0; d < M.gvec_size; d+=M.mesh.n_substrates){
                    copy(M.thomas_denomx[i].begin(), M.thomas_denomx[i].end(), dest_denomx);
                    copy(M.thomas_cx[i].begin(), M.thomas_cx[i].end(), dest_cx);
                    dest_denomx+=M.mesh.n_substrates;
                    dest_cx+=M.mesh.n_substrates;
                }
            }
            //Y-diffusion

            M.gthomas_denomy.resize(M.mesh.y_size);
            M.gthomas_cy.resize(M.mesh.y_size);
            for (int j = 0; j < M.mesh.y_size; ++j){
                M.gthomas_denomy[j].resize(M.gvec_size, 0.0);
                M.gthomas_cy[j].resize(M.gvec_size, 0.0);
                auto dest_denomy = M.gthomas_denomy[j].begin();
                auto dest_cy = M.gthomas_cy[j].begin();
                for (int d = 0; d < M.gvec_size; d+=M.mesh.n_substrates){
                    copy(M.thomas_denomy[j].begin(), M.thomas_denomy[j].end(), dest_denomy);
                    copy(M.thomas_cy[j].begin(), M.thomas_cy[j].end(), dest_cy);
                    dest_denomy+=M.mesh.n_substrates;
                    dest_cy+=M.mesh.n_substrates;
                }
            }
            //Z - diffusion

            M.gthomas_denomz.resize(M.mesh.z_size);
            M.gthomas_cz.resize(M.mesh.z_size);
            for (int k = 0; k < M.mesh.z_size; ++k){
                M.gthomas_denomz[k].resize(M.gvec_size, 0.0);
                M.gthomas_cz[k].resize(M.gvec_size, 0.0);
                auto dest_denomz = M.gthomas_denomz[k].begin();
                auto dest_cz = M.gthomas_cz[k].begin();
                for (int d = 0; d < M.gvec_size; d+=M.mesh.n_substrates){
                    copy(M.thomas_denomz[k].begin(), M.thomas_denomz[k].end(), dest_denomz);
                    copy(M.thomas_cz[k].begin(), M.thomas_cz[k].end(), dest_cz);
                    dest_denomz+=M.mesh.n_substrates;
                    dest_cz+=M.mesh.n_substrates;
                }
            }
            M.diffusion_solver_vectorized_setup_done = true;
            } 
    
        double block3d[M.thomas_i_jump]; //Aux structure of the size: Y*Z*Substrates

        M.apply_dirichlet_boundaries_conditions(rank, size);

        if (rank == 0)
        {
            for (int step = 0; step < granurality; ++step)
            {
                int initial_index = step * M.snd_data_size;
                int limit = (initial_index + M.snd_data_size);
                int limit_vec = limit -(M.snd_data_size%vl);
                #pragma omp parallel for
                for (int index = initial_index; index < limit_vec; index += vl)
                {
                    int index_dec = index;
                    int gd = index%M.gvec_size;

                    __m256d denomx1 = _mm256_loadu_pd(&M.gthomas_denomx[0][gd]);
                    __m256d density1 = _mm256_loadu_pd(&(*M.p_density_vectors)[index]);
                    __m256d aux1 = _mm256_div_pd(density1, denomx1);

                    _mm256_storeu_pd(&(*M.p_density_vectors)[index], aux1);

                    for (int i = 1; i < M.mesh.x_size; i++)
                    {
                        int index_inc = index_dec + M.thomas_i_jump;
                        __m256d constant1 = _mm256_loadu_pd(&M.gthomas_constant1[gd]);
                        __m256d density_curr1 = _mm256_loadu_pd(&(*M.p_density_vectors)[index_dec]);
                        __m256d density_inc1 = _mm256_loadu_pd(&(*M.p_density_vectors)[index_inc]);
                        __m256d denomy1 = _mm256_loadu_pd(&M.gthomas_denomx[i][gd]);
                    
                        density_curr1 = _mm256_fmadd_pd(constant1, density_curr1, density_inc1);
                
                        density_curr1 = _mm256_div_pd(density_curr1, denomy1);
                        _mm256_storeu_pd(&(*M.p_density_vectors)[index_inc], density_curr1);
                        
                        index_dec = index_inc;
                    }
                }

                //Epilogo vectorization
                for (int index = limit_vec; index < limit; ++index)
                {
                    int index_dec = index;
                    int d = index % M.mesh.n_substrates; 

                    (*M.p_density_vectors)[index] /= M.thomas_denomx[0][d];

                    for (int i = 1; i < M.mesh.x_size; i++)
                    {
                        int index_inc = index_dec + M.thomas_i_jump;
                        // axpy(&(*M.microenvironment)[n], M.thomas_constant1, (*M.microenvironment)[n - M.i_jump]);
                        (*M.p_density_vectors)[index_inc] += M.thomas_constant1[d] * (*M.p_density_vectors)[index_dec];
                        
                        //(*M.microenvironment)[n] /= M.thomas_denomx[i];
                        
                        (*M.p_density_vectors)[index_inc] /= M.thomas_denomx[i][d];
                        
                        index_dec = index_inc;
                    }
                }

                if (size > 1) {
                    int x_end = M.mesh.x_size - 1;
                    int offset = step * M.snd_data_size;
                    MPI_Status status;
                    MPI_Isend(&((*M.p_density_vectors)[x_end * M.thomas_i_jump + offset]), M.snd_data_size, MPI_DOUBLE, rank + 1, step, MPI_COMM_WORLD, &send_req[step]);
                }
            }
            //Last iteration
            if (M.snd_data_size_last != 0) {
                int initial_index = granurality * M.snd_data_size;
                #pragma omp parallel for
                for (int index = initial_index; index < initial_index + M.snd_data_size_last; index += M.mesh.n_substrates)
                {
                    int index_dec = index; 
                    for (int d = 0; d < M.mesh.n_substrates; d++)
                    {
                        (*M.p_density_vectors)[index + d] /= M.thomas_denomx[0][d];
                    }

                    for (int i = 1; i < M.mesh.x_size; i++)
                    {
                        int index_inc = index_dec + M.thomas_i_jump;
                        // axpy(&(*M.microenvironment)[n], M.thomas_constant1, (*M.microenvironment)[n - M.i_jump]);
                        for (int d = 0; d < M.mesh.n_substrates; d++)
                        {
                            (*M.p_density_vectors)[index_inc + d] += M.thomas_constant1[d] * (*M.p_density_vectors)[index_dec + d];
                        }

                        //(*M.microenvironment)[n] /= M.thomas_denomx[i];
                        for (int d = 0; d < M.mesh.n_substrates; d++)
                        {
                            (*M.p_density_vectors)[index_inc + d] /= M.thomas_denomx[i][d];
                        }
                        index_dec = index_inc;
                    }
                }

                if (size > 1) {
                    int x_end = M.mesh.x_size - 1;
                    int offset = granurality * M.snd_data_size;
                    MPI_Status status;
                    MPI_Isend(&((*M.p_density_vectors)[x_end * M.thomas_i_jump + offset]), M.snd_data_size_last, MPI_DOUBLE, rank + 1, granurality, MPI_COMM_WORLD, &send_req[granurality]);
                    
                }
            }
        }
        else
        {
            if (rank >= 1 && rank <= (size - 1))
            {
                for (int step = 0; step < granurality; ++step)
                {
                    int initial_index = step * M.snd_data_size;
                    MPI_Irecv(&(block3d[initial_index]), M.rcv_data_size, MPI_DOUBLE, rank-1, step, MPI_COMM_WORLD, &(recv_req[step]));
                }
                if (M.snd_data_size_last != 0)
                    MPI_Irecv(&(block3d[granurality*M.snd_data_size]), M.rcv_data_size_last, MPI_DOUBLE, rank-1, granurality, MPI_COMM_WORLD, &(recv_req[granurality]));
                for (int step = 0; step < granurality; ++step)
                {
                    int initial_index = step * M.snd_data_size;
                    int limit = (initial_index + M.snd_data_size);
                    int limit_vec = limit - (M.snd_data_size%vl);
                    MPI_Wait(&(recv_req[step]), MPI_STATUS_IGNORE);
                    #pragma omp parallel for
                    for (int index = initial_index; index < limit; index += vl)
                    {
                        // axpy(&(*M.microenvironment)[n], M.thomas_constant1, block3d[k][j]);
                        int index_dec = index;
                        int gd = index%M.gvec_size;
                        __m256d constant1 = _mm256_loadu_pd(&M.gthomas_constant1[gd]);
                        __m256d density_curr1 = _mm256_loadu_pd(&(*M.p_density_vectors)[index]);
                        __m256d density_inc1 = _mm256_loadu_pd(&block3d[index]);
                        __m256d denomy1 = _mm256_loadu_pd(&M.gthomas_denomx[0][gd]);
                

                        density_curr1 = _mm256_fmadd_pd(constant1, density_curr1, density_inc1);
            
                        //_mm256_storeu_pd(&microenvironment[index_inc + zd], density_curr);
            
                        //(*(*M.p_density_vectors))[n] /= M.thomas_denomy[j];
                        //Fer unrolling aqui
                
                        //__m256d density = _mm256_loadu_pd(&microenvironment[index_inc + zd]);
                        density_curr1 = _mm256_div_pd(density_curr1, denomy1);
                        _mm256_storeu_pd(&(*M.p_density_vectors)[index], density_curr1);

                        for (int i = 1; i < M.mesh.x_size; i++)
                        {

                            int index_inc = index_dec + M.thomas_i_jump;
                            __m256d constant1 = _mm256_loadu_pd(&M.gthomas_constant1[gd]);
                            __m256d density_curr1 = _mm256_loadu_pd(&(*M.p_density_vectors)[index_dec ]);
                            __m256d density_inc1 = _mm256_loadu_pd(&(*M.p_density_vectors)[index_inc]);
                            __m256d denomy1 = _mm256_loadu_pd(&M.gthomas_denomx[i][gd]);
                    

                            density_curr1 = _mm256_fmadd_pd(constant1, density_curr1, density_inc1);
                
                            //_mm256_storeu_pd(&microenvironment[index_inc + zd], density_curr);
                
                            //(*(*M.p_density_vectors))[n] /= M.thomas_denomy[j];

                            //__m256d density = _mm256_loadu_pd(&microenvironment[index_inc + zd]);
                            density_curr1 = _mm256_div_pd(density_curr1, denomy1);
                            _mm256_storeu_pd(&(*M.p_density_vectors)[index_inc], density_curr1);

                            index_dec = index_inc;
                        }
                        
                    }
                    //Epilogo vectorizacion
                    for (int index = limit_vec; index < limit; ++index)
                    {
                        int index_dec = index;
                        int d = index % M.mesh.n_substrates; 
                        
                        (*M.p_density_vectors)[index] += M.thomas_constant1[d] * block3d[index];
                        (*M.p_density_vectors)[index] /= M.thomas_denomx[0][d];
                        
                        for (int i = 1; i < M.mesh.x_size; i++)
                        {
                            int index_inc = index_dec + M.thomas_i_jump;
                            // axpy(&(*M.microenvironment)[n], M.thomas_constant1, (*M.microenvironment)[n - M.i_jump]);
                        
                            (*M.p_density_vectors)[index_inc] += M.thomas_constant1[d] * (*M.p_density_vectors)[index_dec];
                        
                            //(*M.microenvironment)[n] /= M.thomas_denomx[i];
                            (*M.p_density_vectors)[index_inc] /= M.thomas_denomx[i][d];
                            index_dec = index_inc;
                        }
                    }


                    if (rank < (size - 1))
                    {
                        int x_end = M.mesh.x_size - 1;
                        MPI_Isend(&((*M.p_density_vectors)[x_end * M.thomas_i_jump + initial_index]), M.snd_data_size, MPI_DOUBLE, rank + 1, step, MPI_COMM_WORLD, &send_req[step]);
                    }
                }
                if (M.snd_data_size_last != 0)
                {
                    int initial_index = granurality * M.snd_data_size;
                    MPI_Wait(&recv_req[granurality], MPI_STATUS_IGNORE); //Need to change
                    #pragma omp parallel for
                    for (int index = initial_index; index < initial_index + M.snd_data_size_last; index += M.mesh.n_substrates)
                    {
                        // axpy(&(*M.microenvironment)[n], M.thomas_constant1, block3d[k][j]);
                        int index_dec = index;
                        for (int d = 0; d < M.thomas_k_jump; d++)
                        {
                            (*M.p_density_vectors)[index + d] += M.thomas_constant1[d] * block3d[index + d];
                        }

                        //(*M.microenvironment)[n] /= M.thomas_denomx[0];
                        for (int d = 0; d < M.mesh.n_substrates; d++)
                        {
                            (*M.p_density_vectors)[index + d] /= M.thomas_denomx[0][d];
                        }

                        for (int i = 1; i < M.mesh.x_size; i++)
                        {

                            int index_inc = index_dec + M.thomas_i_jump;
                            // axpy(&(*M.microenvironment)[n], M.thomas_constant1, (*M.microenvironment)[n - M.i_jump]);
                            for (int d = 0; d < M.mesh.n_substrates; d++)
                            {
                                (*M.p_density_vectors)[index_inc + d] += M.thomas_constant1[d] * (*M.p_density_vectors)[index_dec + d];
                            }
                            //(*M.microenvironment)[n] /= M.thomas_denomx[i];
                            for (int d = 0; d < M.mesh.n_substrates; d++)
                            {
                                (*M.p_density_vectors)[index_inc + d] /= M.thomas_denomx[i][d];
                            }

                            index_dec = index_inc;
                        }
                        
                    }
                    // End of computation region
                    if (rank < (size - 1))
                    {
                        int x_end = M.mesh.x_size - 1;
                        MPI_Request aux;
                        MPI_Isend(&((*M.p_density_vectors)[x_end * M.thomas_i_jump + initial_index]), M.snd_data_size_last, MPI_DOUBLE, rank + 1, granurality, MPI_COMM_WORLD, &send_req[granurality]);
                      
                    }
                }
            }
        }
        
        /*-----------------------------------------------------------------------------------*/
        /*                         CODE FOR BACK SUBSITUTION                                 */
        /*-----------------------------------------------------------------------------------*/

        if (rank == (size - 1))
        {
            for (int step = 0; step < granurality; ++step)
            {
                int last_xplane = ((M.mesh.x_size - 1)*M.thomas_i_jump);
                int initial_index = last_xplane + (step * M.snd_data_size);
                int limit = initial_index + M.snd_data_size;
                int limit_vec = limit - (M.snd_data_size%vl);
                #pragma omp parallel for 
                for (int index = initial_index; index < limit_vec; index += vl)
                {
                    int index_aux = index;
                    int gd = (index - last_xplane)%M.gvec_size;
                    for (int i = M.mesh.x_size - 2; i >= 0; i--)
                    {
                        int index_dec = index_aux - M.thomas_i_jump;
                        __m256d cy1 = _mm256_loadu_pd(&M.gthomas_cx[i][gd]);
                        __m256d density_curr1 = _mm256_loadu_pd(&(*M.p_density_vectors)[index_aux]);
                        __m256d density_dec1 = _mm256_loadu_pd(&(*M.p_density_vectors)[index_dec]);

                        density_curr1 = _mm256_fnmadd_pd(cy1, density_curr1, density_dec1);

                        _mm256_storeu_pd(&(*M.p_density_vectors)[index_dec], density_curr1);
                        index_aux = index_dec;
                    }
                }

                //Epilogo Vectorizacion Back Last rank
                
                for (int index = limit_vec; index < limit; ++index){
                    int index_aux = index;
                    int d = (index - last_xplane) % M.mesh.n_substrates;
                    for (int i = M.mesh.x_size - 2; i >= 0; i--)
                    {
                        int index_dec = index_aux - M.thomas_i_jump;
                        // naxpy(&(*M.microenvironment)[n], M.thomas_cx[i], (*M.microenvironment)[n + M.i_jump]);
                        (*M.p_density_vectors)[index_dec] -= M.thomas_cx[i][d] * (*M.p_density_vectors)[index_aux];
                        
                        index_aux = index_dec;
                    }
                }

                if (size > 1) {
                    MPI_Request aux;
                    MPI_Isend(&((*M.p_density_vectors)[step * M.snd_data_size]), M.snd_data_size, MPI_DOUBLE, rank - 1, step, MPI_COMM_WORLD, &send_req[step]);
                }
            }

            //Last iteration
            if (M.snd_data_size_last != 0) {
                int initial_index = ((M.mesh.x_size - 1)*M.thomas_i_jump) + (granurality * M.snd_data_size);
                #pragma omp parallel for
                for (int index = initial_index; index < initial_index + M.snd_data_size_last; index += M.mesh.n_substrates)
                {
                    int index_aux = index;
                    for (int i = M.mesh.x_size - 2; i >= 0; i--)
                    {

                        int index_dec = index_aux - M.thomas_i_jump;
                        // naxpy(&(*M.microenvironment)[n], M.thomas_cx[i], (*M.microenvironment)[n + M.i_jump]);
                        for (int d = 0; d < M.thomas_k_jump; d++)
                        {
                            (*M.p_density_vectors)[index_dec + d] -= M.thomas_cx[i][d] * (*M.p_density_vectors)[index_aux + d];
                        }
                        index_aux = index_dec;
                    }
                }
                if (size > 1) {
                    MPI_Request aux;
                    MPI_Isend(&((*M.p_density_vectors)[granurality * M.snd_data_size]), M.snd_data_size_last, MPI_DOUBLE, rank - 1, granurality, MPI_COMM_WORLD, &send_req[granurality]);
                }
            
            }
        }
        else
        {
            for (int step = 0; step < granurality; ++step)
            {
                MPI_Irecv(&(block3d[step*M.snd_data_size]), M.rcv_data_size, MPI_DOUBLE, rank+1, step, MPI_COMM_WORLD, &recv_req[step]);
            }
            if (M.snd_data_size_last != 0)
                MPI_Irecv(&(block3d[granurality*M.snd_data_size]), M.rcv_data_size_last, MPI_DOUBLE, rank+1, granurality, MPI_COMM_WORLD, &recv_req[granurality]);

            
            for (int step = 0; step < granurality; ++step)
            {
                int last_xplane = ((M.mesh.x_size - 1)*M.thomas_i_jump);
                int initial_index = last_xplane + (step * M.snd_data_size);
                int limit = initial_index + M.snd_data_size;
                int limit_vec = limit - (M.snd_data_size%vl);
                MPI_Wait(&recv_req[step], MPI_STATUS_IGNORE);
                #pragma omp parallel for
                for (int index = initial_index; index < limit_vec; index += vl)
                {
                    int index_aux = index;
                    int index_3d = index - last_xplane;
                    int gd = index_3d%M.gvec_size;
                    __m256d cy1 = _mm256_loadu_pd(&M.gthomas_cx[M.mesh.x_size-1][gd]);
                    __m256d density_curr1 = _mm256_loadu_pd(&(*M.p_density_vectors)[index_aux]);
                    __m256d density_dec1 = _mm256_loadu_pd(&block3d[index_3d]);

                    density_curr1 = _mm256_fnmadd_pd(cy1, density_curr1, density_dec1);

                    _mm256_storeu_pd(&(*M.p_density_vectors)[index_aux], density_curr1);

                    for (int i = M.mesh.x_size - 2; i >= 0; i--)
                    {

                        int index_dec = index_aux - M.thomas_i_jump;
                        __m256d cy1 = _mm256_loadu_pd(&M.gthomas_cx[i][gd]);
                        __m256d density_curr1 = _mm256_loadu_pd(&(*M.p_density_vectors)[index_aux]);
                        __m256d density_dec1 = _mm256_loadu_pd(&(*M.p_density_vectors)[index_dec]);

                        density_curr1 = _mm256_fnmadd_pd(cy1, density_curr1, density_dec1);

                        _mm256_storeu_pd(&(*M.p_density_vectors)[index_dec], density_curr1);
                        index_aux = index_dec;
                    }
                }
                
                //Epilogo Vectorizacion
                for (int index = limit_vec; index < limit; ++index){
                    int index_aux = index;
                    int index_3d = index - last_xplane;
                    int d = (index - last_xplane) % M.mesh.n_substrates;

                    (*M.p_density_vectors)[index_aux] -= M.thomas_cx[M.mesh.x_size - 1][d] * block3d[index_3d];


                    for (int i = M.mesh.x_size - 2; i >= 0; i--)
                    {
                        int index_dec = index_aux - M.thomas_i_jump;
                        // naxpy(&(*M.microenvironment)[n], M.thomas_cx[i], (*M.microenvironment)[n + M.i_jump]);
                        (*M.p_density_vectors)[index_dec] -= M.thomas_cx[i][d] * (*M.p_density_vectors)[index_aux];
                        
                        index_aux = index_dec;
                    }
                }

                if (rank > 0)
                {
                    MPI_Request aux;
                    MPI_Isend(&((*M.p_density_vectors)[step * M.snd_data_size]), M.snd_data_size, MPI_DOUBLE, rank - 1, step, MPI_COMM_WORLD, &send_req[step]);
                    // cout << "Rank " << rank << " has send" << endl;
                }
            }
            if (M.snd_data_size_last != 0)
            {
                int initial_index = ((M.mesh.x_size - 1)*M.thomas_i_jump) + (granurality * M.snd_data_size);
                int index_3d_initial = (granurality * M.snd_data_size);
                MPI_Wait(&recv_req[granurality], MPI_STATUS_IGNORE);
                #pragma omp parallel for
                for (int offset = 0; offset < M.snd_data_size_last; offset += M.mesh.n_substrates)
                {
                    int index_aux = initial_index + offset;
                    //int index = j * M.j_jump + k * M.k_jump + (M.mesh.x_coordinates.size() - 1) * M.i_jump;
                    int index_3d = index_3d_initial + offset;
                    // naxpy(&(*M.microenvironment)[n], M.thomas_cx[M.mesh.x_coordinates.size() - 1], block3d[k][j]);
                    for (int d = 0; d < M.thomas_k_jump; d++)
                    {
                        (*M.p_density_vectors)[index_aux + d] -= M.thomas_cx[M.mesh.x_size - 1][d] * block3d[index_3d + d];
                    }

                    for (int i = M.mesh.x_size - 2; i >= 0; i--)
                    {

                        int index_dec = index_aux - M.thomas_i_jump;
                        // naxpy(&(*M.microenvironment)[n], M.thomas_cx[i], (*M.microenvironment)[n + M.i_jump]);
                        for (int d = 0; d < M.thomas_k_jump; d++)
                        {
                            (*M.p_density_vectors)[index_dec + d] -= M.thomas_cx[i][d] * (*M.p_density_vectors)[index_aux + d];
                        }
                        index_aux = index_dec;
                        
                    }
                }
                if (rank > 0)
                {
                    MPI_Request aux;
                    MPI_Isend(&((*M.p_density_vectors)[granurality * M.snd_data_size]), M.snd_data_size_last, MPI_DOUBLE, rank - 1, granurality, MPI_COMM_WORLD, &send_req[granurality]);
                }
            }
        }

    M.apply_dirichlet_boundaries_conditions(rank, size);

    #pragma omp parallel for
    for (int i = 0; i < M.mesh.x_size; i++)
    {
        //Forward Elimination
        //J = 0
        int gd = 0; //Density pointer
        int index = i * M.thomas_i_jump;
        int zd;
        for (zd = 0; zd + vl < M.thomas_j_jump; zd+=vl)
        {
            __m256d denomy1 = _mm256_loadu_pd(&M.gthomas_denomy[0][gd]);
            __m256d density1 = _mm256_loadu_pd(&(*M.p_density_vectors)[index + zd]);
            gd+=vl;
            if (gd == M.gvec_size) gd = 0;
            __m256d aux1 = _mm256_div_pd(density1, denomy1);
            _mm256_storeu_pd(&(*M.p_density_vectors)[index + zd], aux1);
        }
        //Epilogo
        int ep = M.thomas_j_jump - zd;
        int z_ini = zd / M.mesh.n_substrates;
        int d_ini = zd % M.mesh.n_substrates;
        index = index + z_ini * M.thomas_k_jump;
        for (int k = z_ini; k < M.mesh.z_size; ++k)
        {
            for (int d = d_ini; d < M.mesh.n_substrates; ++d){
                d_ini = 0;
                (*M.p_density_vectors)[index + d] /= M.thomas_denomy[0][d];
            }
            index+=M.thomas_k_jump;
        }
        //J = 1..(y_size-1)
        for (int j = 1; j < M.mesh.y_size; j++)
        {
            int index_base = i * M.thomas_i_jump +  (j-1)*M.thomas_j_jump;
            int index_inc =  index_base + M.thomas_j_jump;
            int zd;
            gd = 0;
            for (zd = 0; zd + vl < M.thomas_j_jump; zd+=vl)
            {
                __m256d constant1 = _mm256_loadu_pd(&M.gthomas_constant1[gd]);
                __m256d density_curr1 = _mm256_loadu_pd(&(*M.p_density_vectors)[index_base + zd]);
                __m256d density_inc1 = _mm256_loadu_pd(&(*M.p_density_vectors)[index_inc + zd]);
                __m256d denomy1 = _mm256_loadu_pd(&M.gthomas_denomy[j][gd]);
                gd+=vl;
                if (gd == M.gvec_size) gd = 0;
                density_curr1 = _mm256_fmadd_pd(constant1, density_curr1, density_inc1);
                density_curr1 = _mm256_div_pd(density_curr1, denomy1);
                _mm256_storeu_pd(&(*M.p_density_vectors)[index_inc + zd], density_curr1);
            }
            //Epilogo
            int ep = M.thomas_j_jump - zd;
            int z_ini = zd / M.mesh.n_substrates;
            int d_ini = zd % M.mesh.n_substrates;
            index_base = index_base + z_ini * M.thomas_k_jump;
            index_inc = index_inc + z_ini * M.thomas_k_jump;
            for (int k = z_ini; k < M.mesh.z_size; ++k)
            {
                for (int d = d_ini; d < M.mesh.n_substrates; ++d){
                    d_ini = 0;
                    (*M.p_density_vectors)[index_inc + d] += M.thomas_constant1[d] * (*M.p_density_vectors)[index_base + d];
                    (*M.p_density_vectors)[index_inc + d] /= M.thomas_denomy[j][d];
                }
                index_base+=M.thomas_k_jump;
                index_inc+=M.thomas_k_jump;
            }
        }
        // Back substitution
        for (int j = M.mesh.y_size - 2; j >= 0; j--)
        {
            int index_base = i * M.thomas_i_jump + (j+1) * M.thomas_j_jump;
            int index_dec = index_base - M.thomas_j_jump;
            int zd;
            gd = 0;
            for ( zd = 0; zd + vl < M.thomas_j_jump; zd+=vl)
            {
                __m256d cy1 = _mm256_loadu_pd(&M.gthomas_cy[j][gd]);
                __m256d density_curr1 = _mm256_loadu_pd(&(*M.p_density_vectors)[index_base + zd]);
                __m256d density_dec1 = _mm256_loadu_pd(&(*M.p_density_vectors)[index_dec+ zd]);
                gd+=vl;
                if (gd == M.gvec_size) gd = 0;

                density_curr1 = _mm256_fnmadd_pd(cy1, density_curr1, density_dec1);

                _mm256_storeu_pd(&(*M.p_density_vectors)[index_dec + zd], density_curr1);
                
            }

            //Epilogo
            int ep = M.thomas_j_jump - zd;
            int z_ini = zd / M.mesh.n_substrates;
            int d_ini = zd % M.mesh.n_substrates;
            index_base = index_base + z_ini * M.thomas_k_jump;
            index_dec = index_dec + z_ini * M.thomas_k_jump;
            for (int k = z_ini; k < M.mesh.z_size; ++k)
            {
                for (int d = d_ini; d < M.mesh.n_substrates; ++d){
                    d_ini = 0;
                    (*M.p_density_vectors)[index_dec + d] -= M.thomas_cy[j][d] * (*M.p_density_vectors)[index_base + d];
                }
                index_base+=M.thomas_k_jump;
                index_dec+=M.thomas_k_jump;
            }

        }
    }

    M.apply_dirichlet_boundaries_conditions(rank, size);

    #pragma omp parallel for collapse(2)
    for (int i = 0; i < M.mesh.x_size; i++)
        {
        for (int j = 0; j < M.mesh.y_size; j++)
            {
            int index = i * M.thomas_i_jump + j * M.thomas_j_jump;
            //(*(*M.p_density_vectors))[n] /= M.thomas_denomz[0];
            for (int d = 0; d < M.mesh.n_substrates; d++)
            {
                (*M.p_density_vectors)[index + d] /= M.thomas_denomz[0][d];
            }

            // should be an empty loop if mesh.z_coordinates.size() < 2
            for (int k = 1; k < M.mesh.z_size; k++)
            {

                int index_inc = index + M.thomas_k_jump;
                // axpy(&(*(*M.p_density_vectors))[n], M.thomas_constant1, (*(*M.p_density_vectors))[n - M.thomas_k_jump]);
                for (int d = 0; d < M.mesh.n_substrates; d++)
                {
                    (*M.p_density_vectors)[index_inc + d] += M.thomas_constant1[d] * (*M.p_density_vectors)[index + d]; //Jose:esta linea esta mal
                }
                //(*(*M.p_density_vectors))[n] /= M.thomas_denomz[k];
                for (int d = 0; d < M.mesh.n_substrates; d++)
                {
                    (*M.p_density_vectors)[index_inc + d] /= M.thomas_denomz[k][d];
                }

                index = index_inc;
            }

            index = i * M.thomas_i_jump + j * M.thomas_j_jump + (M.thomas_k_jump * (M.mesh.z_size - 1));
            for (int k = M.mesh.z_size - 2; k >= 0; k--)
            {

                int index_dec = index - M.thomas_k_jump;
                // naxpy(&(*(*M.p_density_vectors))[n], M.thomas_cz[k], (*(*M.p_density_vectors))[n + M.thomas_k_jump]);
                for (int d = 0; d < M.mesh.n_substrates; d++)
                {
                    (*M.p_density_vectors)[index_dec + d] -= M.thomas_cz[k][d] * (*M.p_density_vectors)[index + d];
                }
                index = index_dec;
            }
        }
    }
}   

//BioFVM-B: need to include communication 
void diffusion_decay_solver__constant_coefficients_LOD_2D( Microenvironment& M, double dt )
{
	if( M.mesh.regular_mesh == false )
	{
		std::cout << "Error: This algorithm is written for regular Cartesian meshes. Try: something else." << std::endl << std::endl; 
		return; 
	}
	
	// constants for the linear solver (Thomas algorithm) 
	
	if( !M.diffusion_solver_setup_done )
	{
		std::cout << std::endl << "Using method " << __FUNCTION__ << " (2D LOD with Thomas Algorithm) ... " << std::endl << std::endl;  
		
		M.thomas_denomx.resize( M.mesh.x_coordinates.size() , M.zero );
		M.thomas_cx.resize( M.mesh.x_coordinates.size() , M.zero );

		M.thomas_denomy.resize( M.mesh.y_coordinates.size() , M.zero );
		M.thomas_cy.resize( M.mesh.y_coordinates.size() , M.zero );
		
		// define constants and pre-computed quantities 

		M.thomas_i_jump = M.number_of_densities(); 
		M.thomas_j_jump = M.mesh.x_coordinates.size(); 

		M.thomas_constant1 =  M.diffusion_coefficients; //   dt*D/dx^2 
		M.thomas_constant1a = M.zero; // -dt*D/dx^2; 
		M.thomas_constant2 =  M.decay_rates; // (1/2)*dt*lambda 
		M.thomas_constant3 = M.one; // 1 + 2*constant1 + constant2; 
		M.thomas_constant3a = M.one; // 1 + constant1 + constant2; 
		
		M.thomas_constant1 *= dt; 
		M.thomas_constant1 /= M.mesh.dx; 
		M.thomas_constant1 /= M.mesh.dx; 

		M.thomas_constant1a = M.thomas_constant1; 
		M.thomas_constant1a *= -1.0; 

		M.thomas_constant2 *= dt; 
		M.thomas_constant2 *= 0.5; // for splitting via LOD

		M.thomas_constant3 += M.thomas_constant1; 
		M.thomas_constant3 += M.thomas_constant1; 
		M.thomas_constant3 += M.thomas_constant2; 

		M.thomas_constant3a += M.thomas_constant1; 
		M.thomas_constant3a += M.thomas_constant2; 
		
		// Thomas solver coefficients 

		M.thomas_cx.assign( M.mesh.x_coordinates.size() , M.thomas_constant1a ); 
		M.thomas_denomx.assign( M.mesh.x_coordinates.size()  , M.thomas_constant3 ); 
		M.thomas_denomx[0] = M.thomas_constant3a; 
		M.thomas_denomx[ M.mesh.x_coordinates.size()-1 ] = M.thomas_constant3a; 
		if( M.mesh.x_coordinates.size() == 1 )
		{ M.thomas_denomx[0] = M.one; M.thomas_denomx[0] += M.thomas_constant2; } 

		M.thomas_cx[0] /= M.thomas_denomx[0]; 
		for( unsigned int i=1 ; i <= M.mesh.x_coordinates.size()-1 ; i++ )
		{ 
			axpy( &M.thomas_denomx[i] , M.thomas_constant1 , M.thomas_cx[i-1] ); 
			M.thomas_cx[i] /= M.thomas_denomx[i]; // the value at  size-1 is not actually used  
		}

		M.thomas_cy.assign( M.mesh.y_coordinates.size() , M.thomas_constant1a ); 
		M.thomas_denomy.assign( M.mesh.y_coordinates.size()  , M.thomas_constant3 ); 
		M.thomas_denomy[0] = M.thomas_constant3a; 
		M.thomas_denomy[ M.mesh.y_coordinates.size()-1 ] = M.thomas_constant3a; 
		if( M.mesh.y_coordinates.size() == 1 )
		{ M.thomas_denomy[0] = M.one; M.thomas_denomy[0] += M.thomas_constant2; } 

		M.thomas_cy[0] /= M.thomas_denomy[0]; 
		for( unsigned int i=1 ; i <= M.mesh.y_coordinates.size()-1 ; i++ )
		{ 
			axpy( &M.thomas_denomy[i] , M.thomas_constant1 , M.thomas_cy[i-1] ); 
			M.thomas_cy[i] /= M.thomas_denomy[i]; // the value at  size-1 is not actually used  
		}

		M.diffusion_solver_setup_done = true; 
	}

	// set the pointer
	
	M.apply_dirichlet_conditions_fast();

	// x-diffusion 
	#pragma omp parallel for 
	for( unsigned int j=0; j < M.mesh.y_coordinates.size() ; j++ )
	{
		// Thomas solver, x-direction

		// remaining part of forward elimination, using pre-computed quantities 
		unsigned int n = M.voxel_index(0,j,0) * M.n_subs;
        for (int d = 0; d < M.n_subs; ++d)
		    (*M.p_density_vectors)[n+d] /= M.thomas_denomx[0][d]; 

		n += M.thomas_i_jump; 
		for( unsigned int i=1; i < M.mesh.x_coordinates.size() ; i++ )
		{
            //axpy( &(*M.p_density_vectors)[n] , M.thomas_constant1 , (*M.p_density_vectors)[n-M.thomas_i_jump] );
            for (int d = 0; d < M.n_subs; ++d)
                (*M.p_density_vectors)[n+d] += M.thomas_constant1[d] * (*M.p_density_vectors)[n-M.thomas_i_jump +d];
			
            for (int d = 0; d < M.n_subs; ++d)
                (*M.p_density_vectors)[n+d] /= M.thomas_denomx[i][d]; 
			
            n += M.thomas_i_jump; 
		}

		// back substitution 
		n = M.voxel_index( M.mesh.x_coordinates.size()-2 ,j,0) * M.n_subs; 

		for( int i = M.mesh.x_coordinates.size()-2 ; i >= 0 ; i-- )
		{
			//naxpy( &(*M.p_density_vectors)[n] , M.thomas_cx[i] , (*M.p_density_vectors)[n+M.thomas_i_jump] ); 
			for (int d = 0; d < M.n_subs; ++d)
                (*M.p_density_vectors)[n+d] -= M.thomas_cx[i][d] * (*M.p_density_vectors)[n+M.thomas_i_jump +d];
            n -= M.thomas_i_jump; 
		}
	}

	// y-diffusion 

	M.apply_dirichlet_conditions_fast();
	#pragma omp parallel for 
	for( unsigned int i=0; i < M.mesh.x_coordinates.size() ; i++ )
	{
		// Thomas solver, y-direction

		// remaining part of forward elimination, using pre-computed quantities 

		int n = M.voxel_index(i,0,0) * M.n_subs;
        for (int d = 0; d < M.n_subs; ++d)
		    (*M.p_density_vectors)[n + d] /= M.thomas_denomy[0][d]; 

		n += M.thomas_j_jump; 
		for( unsigned int j=1; j < M.mesh.y_coordinates.size() ; j++ )
		{
			//axpy( &(*M.p_density_vectors)[n] , M.thomas_constant1 , (*M.p_density_vectors)[n-M.thomas_j_jump] ); 
			for (int d = 0; d < M.n_subs; ++d)
                (*M.p_density_vectors)[n+d] += M.thomas_constant1[d] * (*M.p_density_vectors)[n-M.thomas_j_jump +d];
            for (int d = 0; d < M.n_subs; ++d)
                (*M.p_density_vectors)[n+d] /= M.thomas_denomy[j][d]; 
			n += M.thomas_j_jump; 
		}

		// back substitution 
		n = M.voxel_index( i,M.mesh.y_coordinates.size()-2, 0) * M.n_subs; 

		for( int j = M.mesh.y_coordinates.size()-2 ; j >= 0 ; j-- )
		{
			//naxpy( &(*M.p_density_vectors)[n] , M.thomas_cy[j] , (*M.p_density_vectors)[n+M.thomas_j_jump] ); 
			for (int d = 0; d < M.n_subs; ++d)
                (*M.p_density_vectors)[n+d] += M.thomas_cy[j][d] * (*M.p_density_vectors)[n+M.thomas_j_jump +d];
            n -= M.thomas_j_jump; 
		}
	}

	M.apply_dirichlet_conditions_fast();
	
	// reset gradient vectors 
//	M.reset_all_gradient_vectors(); 
	
	return; 
}

void diffusion_decay_explicit_uniform_rates( Microenvironment& M, double dt )
{
    std::cout << "Diffusion decay explicit uniform rates is desactivated" << std::endl;
    
	using std::vector; 
	using std::cout; 
	using std::endl; 

	// static int n_jump_i = 1; 
	// static int n_jump_j = M.mesh.x_coordinates.size(); 
	// static int n_jump_k = M.mesh.x_coordinates.size() * M.mesh.y_coordinates.size(); 

	if( !M.diffusion_solver_setup_done )
	{	
		M.thomas_i_jump = M.mesh.y_coordinates.size() * M.mesh.z_coordinates.size() * M.number_of_densities(); 
		M.thomas_j_jump = M.mesh.z_coordinates.size() * M.number_of_densities(); 
		M.thomas_k_jump = M.number_of_densities();
        M.n_subs = M.number_of_densities(); 
	
		M.diffusion_solver_setup_done = true; 
	}
	
	if( M.mesh.uniform_mesh == false )
	{
		cout << "Error: This algorithm is written for uniform Cartesian meshes. Try: something else" << endl << endl; 
		return; 
	}

	// double buffering to reduce memory copy / allocation overhead 

	static vector< double >* pNew = &(M.temporary_density_vectors1); //Alias for p_density_vectors
	static vector< double >* pOld = &(M.temporary_density_vectors2);

	// swap the buffers 

	vector< double >* pTemp = pNew; 
	pNew = pOld; 
	pOld = pTemp; 
	M.p_density_vectors = pNew; 

	// static bool reaction_diffusion_shortcuts_are_set = false; 

	static vector<double> constant1 = (1.0 / ( M.mesh.dx * M.mesh.dx )) * M.diffusion_coefficients; 
	static vector<double> constant2 = dt * constant1; 
	static vector<double> constant3 = M.one + dt * M.decay_rates;

	static vector<double> constant4 = M.one - dt * M.decay_rates;

	#pragma omp parallel for
	for( unsigned int i=0; i < M.number_of_voxels() ; i++ )
	{
        int index = i * M.n_subs;
		unsigned int number_of_neighbors = M.mesh.connected_voxel_indices[i].size(); 

		double d1 = -1.0 * number_of_neighbors; 

        for (int d = 0; d < M.n_subs; ++d) {
            (*pNew)[index +d] = (*pOld)[index +d] * constant4[d];  
        }
		 
		for( unsigned int j=0; j < number_of_neighbors ; j++ )
		{
			//axpy( &(*pNew)[i], constant2, (*pOld)[  M.mesh.connected_voxel_indices[i][j] ] ); 
            for (int d = 0; d < M.n_subs; ++d) {
                (*pNew)[index +d] += constant2[d] * (*pOld)[  M.mesh.connected_voxel_indices[i][j] * M.n_subs +d ];
            }
		}
		vector<double> temp = constant2; 
		temp *= d1; 
		//axpy( &(*pNew)[i] , temp , (*pOld)[i] ); 
        for (int d = 0; d < M.n_subs; ++d) {
                (*pNew)[index +d] += temp[d] * (*pOld)[index+d];
            }
	}
	
	// reset gradient vectors 
//	M.reset_all_gradient_vectors(); 

	return; 
}

void diffusion_decay_solver__constant_coefficients_LOD_1D( Microenvironment& M, double dt )
{
	if( M.mesh.regular_mesh == false )
	{
		std::cout << "Error: This algorithm is written for regular Cartesian meshes. Try: something else." << std::endl << std::endl; 
		return; 
	}

	// constants for the linear solver (Thomas algorithm) 
	
	if( !M.diffusion_solver_setup_done )
	{
		std::cout << std::endl << "Using method " << __FUNCTION__ << " (2D LOD with Thomas Algorithm) ... " << std::endl << std::endl;  
		
		M.thomas_denomx.resize( M.mesh.x_coordinates.size() , M.zero );
		M.thomas_cx.resize( M.mesh.x_coordinates.size() , M.zero );

		// define constants and pre-computed quantities 

		M.thomas_i_jump = M.mesh.y_coordinates.size() * M.n_subs; 
		M.thomas_j_jump = M.n_subs; 

		M.thomas_constant1 =  M.diffusion_coefficients; //   dt*D/dx^2 
		M.thomas_constant1a = M.zero; // -dt*D/dx^2; 
		M.thomas_constant2 =  M.decay_rates; // (1/2)*dt*lambda 
		M.thomas_constant3 = M.one; // 1 + 2*constant1 + constant2; 
		M.thomas_constant3a = M.one; // 1 + constant1 + constant2; 
		
		M.thomas_constant1 *= dt; 
		M.thomas_constant1 /= M.mesh.dx; 
		M.thomas_constant1 /= M.mesh.dx; 

		M.thomas_constant1a = M.thomas_constant1; 
		M.thomas_constant1a *= -1.0; 

		M.thomas_constant2 *= dt; 
		M.thomas_constant2 *= 1; // no splitting via LOD

		M.thomas_constant3 += M.thomas_constant1; 
		M.thomas_constant3 += M.thomas_constant1; 
		M.thomas_constant3 += M.thomas_constant2; 

		M.thomas_constant3a += M.thomas_constant1; 
		M.thomas_constant3a += M.thomas_constant2; 
		
		// Thomas solver coefficients 

		M.thomas_cx.assign( M.mesh.x_coordinates.size() , M.thomas_constant1a ); 
		M.thomas_denomx.assign( M.mesh.x_coordinates.size()  , M.thomas_constant3 ); 
		M.thomas_denomx[0] = M.thomas_constant3a; 
		M.thomas_denomx[ M.mesh.x_coordinates.size()-1 ] = M.thomas_constant3a; 
		if( M.mesh.x_coordinates.size() == 1 )
		{ M.thomas_denomx[0] = M.one; M.thomas_denomx[0] += M.thomas_constant2; } 

		M.thomas_cx[0] /= M.thomas_denomx[0]; 
		for( unsigned int i=1 ; i <= M.mesh.x_coordinates.size()-1 ; i++ )
		{ 
			axpy( &M.thomas_denomx[i] , M.thomas_constant1 , M.thomas_cx[i-1] ); 
			M.thomas_cx[i] /= M.thomas_denomx[i]; // the value at  size-1 is not actually used  
		}

		M.diffusion_solver_setup_done = true; 
	}

	// set the pointer
	
	M.apply_dirichlet_conditions_fast();

	// x-diffusion 
	#pragma omp parallel for 
	for( unsigned int j=0; j < M.mesh.y_coordinates.size() ; j++ )
	{
		// Thomas solver, x-direction

		// remaining part of forward elimination, using pre-computed quantities 
		unsigned int n = M.voxel_index(0,j,0) * M.n_subs;
		
        for (int d = 0; d < M.n_subs; ++d)
            (*M.p_density_vectors)[n +d] /= M.thomas_denomx[0][d]; 

		n += M.thomas_i_jump; 
		for( unsigned int i=1; i < M.mesh.x_coordinates.size() ; i++ )
		{
			//axpy( &(*M.p_density_vectors)[n] , M.thomas_constant1 , (*M.p_density_vectors)[n-M.thomas_i_jump] ); 
            for (int d = 0; d < M.n_subs; ++d)
                (*M.p_density_vectors)[n+d] += M.thomas_constant1[d] * (*M.p_density_vectors)[n-M.thomas_i_jump +d];
            for (int d = 0; d < M.n_subs; ++d)
			    (*M.p_density_vectors)[n+d] /= M.thomas_denomx[i][d]; 
			n += M.thomas_i_jump; 
		}

		// back substitution 
		n = M.voxel_index( M.mesh.x_coordinates.size()-2 ,j,0) * M.n_subs; 

		for( int i = M.mesh.x_coordinates.size()-2 ; i >= 0 ; i-- )
		{
			//naxpy( &(*M.p_density_vectors)[n] , M.thomas_cx[i] , (*M.p_density_vectors)[n+M.thomas_i_jump] ); 
            for (int d = 0; d < M.n_subs; ++d)
                (*M.p_density_vectors)[n+d] -= M.thomas_cx[i][d] * (*M.p_density_vectors)[n+M.thomas_i_jump+d];
			n -= M.thomas_i_jump; 
		}
	}

	M.apply_dirichlet_conditions_fast();
	
	// reset gradient vectors 
//	M.reset_all_gradient_vectors(); 
	
	return; 
}


};
