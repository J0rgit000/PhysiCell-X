/*
###############################################################################
# If you use PhysiCell in your project, please cite PhysiCell and the version #
# number, such as below:                                                      #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1].    #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# See VERSION.txt or call get_PhysiCell_version() to get the current version  #
#     x.y.z. Call display_citations() to get detailed information on all cite-#
#     able software used in your PhysiCell application.                       #
#                                                                             #
# Because PhysiCell extensively uses BioFVM, we suggest you also cite BioFVM  #
#     as below:                                                               #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1],    #
# with BioFVM [2] to solve the transport equations.                           #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# [2] A Ghaffarizadeh, SH Friedman, and P Macklin, BioFVM: an efficient para- #
#     llelized diffusive transport solver for 3-D biological simulations,     #
#     Bioinformatics 32(8): 1256-8, 2016. DOI: 10.1093/bioinformatics/btv730  #
#                                                                             #
###############################################################################
#                                                                             #
# BSD 3-Clause License (see https://opensource.org/licenses/BSD-3-Clause)     #
#                                                                             #
# Copyright (c) 2015-2018, Paul Macklin and the PhysiCell Project             #
# All rights reserved.                                                        #
#                                                                             #
# Redistribution and use in source and binary forms, with or without          #
# modification, are permitted provided that the following conditions are met: #
#                                                                             #
# 1. Redistributions of source code must retain the above copyright notice,   #
# this list of conditions and the following disclaimer.                       #
#                                                                             #
# 2. Redistributions in binary form must reproduce the above copyright        #
# notice, this list of conditions and the following disclaimer in the         #
# documentation and/or other materials provided with the distribution.        #
#                                                                             #
# 3. Neither the name of the copyright holder nor the names of its            #
# contributors may be used to endorse or promote products derived from this   #
# software without specific prior written permission.                         #
#                                                                             #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" #
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE   #
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE  #
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE   #
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR         #
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF        #
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS    #
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN     #
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)     #
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE  #
# POSSIBILITY OF SUCH DAMAGE.                                                 #
#                                                                             #
###############################################################################
*/

/*===================================================================================*
+ If you use PhysiCell-X in your project, we would really appreciate if you can 	 +
+																					 +
+ [1] Cite the PhysiCell-X repository by giving its URL								 +
+																					 +
+ [2] Cite BioFVM-X: 																 +
+		Saxena, Gaurav, Miguel Ponce-de-Leon, Arnau Montagud, David Vicente Dorca,   +
+		and Alfonso Valencia. "BioFVM-X: An MPI+ OpenMP 3-D Simulator for Biological + 
+		Systems." In International Conference on Computational Methods in Systems    +
+		Biology, pp. 266-279. Springer, Cham, 2021. 								 +
*====================================================================================*/

#include "../core/PhysiCell.h"
#include "../modules/PhysiCell_standard_modules.h" 
#include "./tnf_receptor_dynamics.h"
#include "./tnf_boolean_model_interface.h"

#include "../addons/PhysiBoSS/src/maboss_network.h"

using namespace BioFVM; 
using namespace PhysiCell;

struct init_record
{
	float x;
	float y;
	float z;
	float radius;
	int phase;
};

// setup functions to help us along 
void create_cell_types( mpi_Environment &world, mpi_Cartesian &cart_topo );
void setup_tissue( void );

/*======================================*/
/* Parallel prototype of setup_tissue() */
/*======================================*/
void setup_tissue(Microenvironment &m, mpi_Environment &world, mpi_Cartesian &cart_topo); 
// set up the intial cells 
void setup_tissue( void );

/*================================================*/
/* Parallel prototype of setup_microenvironment() */
/*================================================*/
void setup_microenvironment(mpi_Environment &world, mpi_Cartesian &cart_topo);
// set up the BioFVM microenvironment 
void setup_microenvironment( void ); 

// helper function to read init files
std::vector<init_record> read_init_file(std::string filename, char delimiter, bool header);

// helper function to calculate cell positions within a sphere
std::vector<std::vector<double>> create_cell_sphere_positions(double cell_radius, double sphere_radius);

// helper function that calculates phere volume
inline float sphere_volume_from_radius(float radius) {return 4/3 * PhysiCell_constants::pi * std::pow(radius, 3);}

/*===============================================*/
/* Parallel prototype of inject_density_sphere() */
/*===============================================*/
void inject_density_sphere( int density_index, double concentration, double membrane_lenght, 
							mpi_Environment &world, mpi_Cartesian &cart_topo );

// helper function to inject density surrounding a spheroid
void inject_density_sphere(int density_index, double concentration, double membrane_lenght);

/*===============================================*/
/* Parallel prototype of remove_density() */
/*===============================================*/
void remove_density( int density_index, mpi_Environment &world, mpi_Cartesian &cart_topo );

// helper function to remove a density
void remove_density( int density_index );

// custom pathology coloring function 
std::vector<std::string> my_coloring_function( Cell* );

/*===============================================*/
/* Parallel prototype of total_live_cell_count() */
/*===============================================*/
int total_basic_agent_count(mpi_Environment &world, mpi_Cartesian &cart_topo);
int total_cell_agent_count(mpi_Environment &world, mpi_Cartesian &cart_topo);
int total_live_cell_count(mpi_Environment &world, mpi_Cartesian &cart_topo);
int total_dead_cell_count(mpi_Environment &world, mpi_Cartesian &cart_topo);
int total_necrosis_cell_count(mpi_Environment &world, mpi_Cartesian &cart_topo);
int total_apoptosis_cell_count(mpi_Environment &world, mpi_Cartesian &cart_topo);

double get_total_tnf(mpi_Environment &world, mpi_Cartesian &cart_topo);

double total_custom_variable_live(std::string var_name, mpi_Environment &world, mpi_Cartesian &cart_topo);

double total_active_TNF_receptor(mpi_Environment &world, mpi_Cartesian &cart_topo);
double total_free_TNF_receptor(mpi_Environment &world, mpi_Cartesian &cart_topo);
double total_internalized_TNF_receptor(mpi_Environment &world, mpi_Cartesian &cart_topo);

double total_active_TNF_node(mpi_Environment &world, mpi_Cartesian &cart_topo);
double total_active_FADD_node(mpi_Environment &world, mpi_Cartesian &cart_topo);
double total_active_NFKb_node(mpi_Environment &world, mpi_Cartesian &cart_topo);
