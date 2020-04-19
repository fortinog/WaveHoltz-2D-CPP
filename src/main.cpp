#include <iostream>
#include <iomanip>
#include <cmath>
#include <chrono> 
#include <vector>
#include "Darray2.h"
#include "problem_setup.hpp"
#include "wave_solve.hpp"
#include <mpi.h>
#include "domain_decomposition.h"
using namespace std;
using namespace std::chrono;


int main(int argc, char * argv[])
{
    ProblemSetup setup;
    int num_proc;
    int rank;
    int dims[2] = {0,0};
    int rerank = 1;
    int N,M;
    int node_rank;
    int node_ID;
    int n_node;
    FILE *extFile;
    Subdomain loc_grid;
    MPI_Comm Cart_comm;
    MPI_Comm IO_Comm;
    MPI_Comm Inter_Node;
    MPI_Info info;

    // Initialize MPI program
    MPI_Init( &argc, &argv );
    MPI_Comm_size( MPI_COMM_WORLD, &num_proc );
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );

    // Dump out global x,y data
    if(rank == 0){
        extFile = fopen("x.txt", "w");
        for(int i=1;i<setup.N-1;i++){
            fprintf(extFile, "%18.10e\n", setup.x_L + i*setup.hx);
        }
        fclose(extFile);
        extFile = fopen("y.txt", "w");
        for(int i=1;i<setup.M-1;i++){
            fprintf(extFile, "%18.10e\n", setup.y_L + i*setup.hy);
        }
        fclose(extFile);

    }
    MPI_Barrier(MPI_COMM_WORLD);

    N = setup.N;
    M = setup.M;

    // Initialize new default communicator to account for architecture, as well as 
    // partition determine chunk of computational domain for each process
    Cart_comm  = mpi_decompose_domain(num_proc,N,M,rerank,&rank,dims,&loc_grid);

    // Group all processes in shared memory machine (node) for I/O
    MPI_Info_create(&info);
    MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, rank, info, &IO_Comm);

    // Grab rank on compute node, and create communicator only with rank 0 of each node
    MPI_Comm_rank(IO_Comm, &node_rank);
    MPI_Comm_split(MPI_COMM_WORLD, node_rank, rank, &Inter_Node);
    MPI_Comm_size(Inter_Node, &n_node);


    // Grab node name and then assign each node a unique number for this job
    if(node_rank == 0){
        char name[MPI_MAX_PROCESSOR_NAME];
        char names[n_node][MPI_MAX_PROCESSOR_NAME];
        int len;
        MPI_Get_processor_name(name,&len);

        // Create list of node names
        MPI_Allgather(name, MPI_MAX_PROCESSOR_NAME, MPI_CHAR, names, MPI_MAX_PROCESSOR_NAME, MPI_CHAR, Inter_Node);
        
        // Match our node name to the list of nodes to define a unique ID per node
        for(int i=0;i<n_node;i++){
            if(strcmp(names[i], name)==0){
                node_ID = i;
                break;
            }
        }
    }
    MPI_Barrier(Inter_Node);
    MPI_Comm_free(&Inter_Node);


    // Naming convention for files: w_NODE-ID_ITER-NUMBER.txt?

    // Do a forward wave solve
    Wave_Solve solver(loc_grid, node_ID, Cart_comm, IO_Comm);

    // Free communicators
    MPI_Comm_free(&Cart_comm);
    MPI_Comm_free(&IO_Comm);

    // Finalize program
    int finalize_retcode = MPI_Finalize();
    return finalize_retcode;
}