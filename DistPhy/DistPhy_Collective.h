#ifndef __DistPhy_Collective_h__
#define __DistPhy_Collective_h__

#include "DistPhy_Cartesian.h"

namespace DistPhy
{
	namespace mpi
	{
				double distribute_global_sum(double x, mpi_Cartesian &cart_topo);
				int distribute_global_sum(int x, mpi_Cartesian &cart_topo);
				double distribute_global_max(double x, mpi_Cartesian &cart_topo); 
				double distribute_global_min(double x, mpi_Cartesian &cart_topo); 
	}
}

#endif