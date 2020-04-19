#ifndef WAVE_SOLVE
#define WAVE_SOLVE

#include <string>   
#include <fstream>
#include <iostream>
#include <iomanip>
#include <math.h>
#include <memory>
#include "problem_setup.hpp"
#include "Darray1.h"
#include "Darray2.h"
#include "Iarray2.h"
#include "twilight_2d.h"
#include "domain_decomposition.h"
#include "H5Cpp.h"
using namespace std;


class Wave_Solve {
public:
	ProblemSetup setup;

	// MPI derived data types for communication
    MPI_Datatype sub_send_left;
    MPI_Datatype sub_recv_left;
    MPI_Datatype sub_send_right;
    MPI_Datatype sub_recv_right;
    MPI_Datatype sub_send_up;
    MPI_Datatype sub_recv_up;
    MPI_Datatype sub_send_down;
    MPI_Datatype sub_recv_down;

    // Requests for persistent communication
    MPI_Request send_req[4];
    MPI_Request recv_req[4];

	Iarray2 mask, ghost_list_left, ghost_list_right, ghost_list_up, ghost_list_down;
	Darray1 x, y, dist_list_left, dist_list_right, dist_list_up, dist_list_down;
	double t;
	double* IO_buf;

	int istart, iend, jstart, jend, N, M;
	int j_left_start, j_left_end, i_down_start, i_down_end;
	int j_right_start, j_right_end, i_up_start, i_up_end;

	int up_neigh, down_neigh, left_neigh, right_neigh;
	int size_lr, size_ud;

	Wave_Solve(Subdomain Local_Grid, int node_ID, MPI_Comm CART_COMM, MPI_Comm IO_Comm);
	virtual ~Wave_Solve() {if( IO_buf != 0 ) delete[] IO_buf;};
	void    Compute_Mask();
	void    Find_Ghost_Points();
	void    Fill_In_Ghost_Values();
	void    Enforce_BC(Darray2 &v);
	void    Compute_Laplacian(Darray2& w, Darray2& lap);
	void    Compute_Laplacian_NB(Darray2& w, Darray2& lap, MPI_Request* recv_req);
	void    Set_Initial_Data(Darray2& wm, Darray2& w, Darray2& lap);
	double  Compute_Laplacian_Error(const int flag, Darray2& lap, MPI_Comm CART_COMM);
	double  Compute_Solution_Error(const int flag, Darray2& w, MPI_Comm CART_COMM);
	void    Time_Step(Darray2& wm, Darray2& w, Darray2& wp, Darray2& lap);
	void    Taylor_Expand(Darray2& wm, Darray2& w, Darray2& lap);
	double  Compute_Energy(Darray2& wm, Darray2& w, Darray2& lap, MPI_Comm CART_COMM);
	double  forcing(const double x, const double y, const double t);
	void    Solve_PDE(Darray2& wm, Darray2& w, Darray2& wp, Darray2& lap, double* w_ptr, MPI_Comm CART_COMM);
	void    Communicate_Solution(MPI_Comm CART_COMM, double* w_ptr, MPI_Request* send_req, MPI_Request* recv_req);
	void    Setup_Subarrays(const int nolp);
	void    Print_Solution(char* FileName,double* IO_buf, Darray2& w, int* size_array, int proc_per_node, int node_rank, MPI_Comm IO_Comm);
	void    Finalize();
};
#endif