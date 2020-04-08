#include <iostream>
#include <iomanip>
#include <cmath>
#include <chrono> 
#include <vector>
#include "Darray.h"
#include "Darray2.h"
#include "problem_setup.hpp"
#include "wave_solve.hpp"
#include <mpi.h>
#include "domain_decomposition.h"
using namespace std;
using namespace std::chrono;

// void print_solution(int m_Ny, int m_Nx, Wave_Solve solver){

//     cout << "\n------------------------------\n";
//     // for (int j = 0; j<= m_Ny; j++){
//     //     for (int i = 0; i<= m_Nx; i++)
//     //         cout << solver.w(i,j) << " ";
//     //     cout << "\n";
//     // }
//     for (int i = 0; i<= m_Nx; i++){
//         for (int j = 0; j<= m_Ny; j++)
//             cout << solver.lap(i,j) << " ";
//         cout << "\n";
//     }
//     cout << "------------------------------\n";
            
// }

int main(int argc, char * argv[])
{
    // Wave_Solve solver;
    ProblemSetup setup;
    int num_proc;
    int rank;
    int dims[2] = {0,0};
    int rerank = 1;
    int N,M;
    Subdomain loc_grid;
    MPI_Comm Cart_comm;

    MPI_Init( &argc, &argv );
    MPI_Comm_size( MPI_COMM_WORLD, &num_proc );
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );

    N = setup.N;
    M = setup.M;

    // Initialize new communicator
    Cart_comm  = mpi_decompose_domain(num_proc,N,M,rerank,&rank,dims,&loc_grid);
    Wave_Solve solver(loc_grid, Cart_comm);
    // solver.Time_Step();
    // solver.Compute_Laplacian();

    // solver.Solve_PDE(Cart_comm);
    // print_solution(loc_grid.ny_loc+1, loc_grid.nx_loc+1, solver);

    // for(int i=0;i<num_proc;i++){
    //     if(rank == i){ 
    //         cout << "x_s : " << loc_grid.x_s << " x_e: " << loc_grid.x_e << "y_s : " << loc_grid.y_s << " y_e: " << loc_grid.y_e << "\n";
    //         cout << "nx_loc : " << loc_grid.nx_loc << " ny_loc: " <<  loc_grid.ny_loc << "\n";
    //         cout << "Left : " << loc_grid.left_neigh << " Right : " << loc_grid.right_neigh << " Down : " << loc_grid.down_neigh << " Up : " << loc_grid.up_neigh  << "\n";
    //         cout << "Dims : " << dims[0] << " Dims : " << dims[1] << "\n";
    //         print_solution(solver.size_ud-1, solver.size_lr-1, solver);
    //     }
    //     MPI_Barrier(Cart_comm);
    // }
    // cout << "The computed errors are: " << "\n";

    // // Compute error
    // if(rank == 0){
       // cout << scientific << right << setw(14)<< solver.Compute_Solution_Error(1, Cart_comm);
       //  cout << scientific << right << setw(14)<< solver.Compute_Solution_Error(2, Cart_comm);
       //  cout << scientific << right << setw(14)<< solver.Compute_Solution_Error(3, Cart_comm);
       // //  cout << "\n";     
    //    cout << scientific << right << setw(14)<< solver.Compute_Laplacian_Error(1, Cart_comm);
    //     cout << scientific << right << setw(14)<< solver.Compute_Laplacian_Error(2, Cart_comm);
    //     cout << scientific << right << setw(14)<< solver.Compute_Laplacian_Error(3, Cart_comm);
    //     cout << "\n";    
    // // }

    // solver.Finalize();
    MPI_Comm_free(&Cart_comm);
    int finalize_retcode = MPI_Finalize();
    if(0 == rank) fprintf(stderr, "Process, return_code\n");
    fprintf(stderr, "%i, %i\n", rank, finalize_retcode);
    // MPI_Finalize();

    return 0;
}

