#include<iostream>
#include<cstdlib>                   //Contains header for exit() in C++

#include "DistPhy_Environment.h"    //See if we can possibly remove this


using namespace std;

namespace DistPhy
{
    namespace mpi
    {
        void mpi_Environment::Initialize()
        {
            //int req_thread_level = MPI_THREAD_FUNNELED; 
            int req_thread_level = MPI_THREAD_MULTIPLE;
            int provided_thread_level; 
            
            
            MPI_Init_thread(NULL, NULL, req_thread_level, &provided_thread_level);
            
            if(provided_thread_level != req_thread_level)
            {
                cout<<"The requested thread level cannot be provided, please provide a lower thread level";
                exit(1); 
            }
            
            init_comm = MPI_COMM_WORLD; 
            MPI_Comm_size(init_comm, &size);
            MPI_Comm_rank(init_comm, &rank);
            
            if(rank == 0)
            {
cout<<"  _____  _               _  _____     _ _     __   __        \n";
cout<<" |  __ \\| |             (_)/ ____|   | | |    \\ \\ / /       .-.--   \n";
cout<<" | |__) | |__  _   _ ___ _| |     ___| | |_____\\ V /     ,´,´.´   `.    \n";
cout<<" |  ___/| '_ \\| | | / __| | |    / _ \\ | |______> <      | | | BSC |    \n";
cout<<" | |    | | | | |_| \\__ \\ | |___|  __/ | |     / . \\     `.`.`. _ .´   \n";
cout<<" |_|    |_| |_|\\__, |___/_|\\_____\\___|_|_|    /_/ \\_\\       `·`··      \n";
cout<<"                __/ |                                        \n";
cout<<"               |___/                                         \n";  
}
                
        }
        
        int mpi_Environment::Size()
        {
            return size;
        }
        
        int mpi_Environment::Rank()
        {
            return rank; 
        }
        
        void mpi_Environment::Finalize()
        {
            MPI_Finalize(); 
        }
        
        //Non-class functions in namespace DistPhy::mpi
        int IOProcessor(mpi_Environment &world)
        {
            if(world.rank == 0)
                return 1;
            else
                return 0;
        }
    }
}
