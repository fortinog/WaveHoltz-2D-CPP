/*
	This is the class that will (eventually) do a wave solve. For now
	we simply load the problem setup, define the grid, compute the 
	mask, and then compute the Laplacian
*/
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
	// Darray2 w, lap, wm, wp;
    MPI_Datatype sub_send_left;
    MPI_Datatype sub_recv_left;
    MPI_Datatype sub_send_right;
    MPI_Datatype sub_recv_right;
    MPI_Datatype sub_send_up;
    MPI_Datatype sub_recv_up;
    MPI_Datatype sub_send_down;
    MPI_Datatype sub_recv_down;


	Iarray2 mask;
	Darray1 x, y;
	double t;
	// double* w_ptr;

	int istart, iend, jstart, jend, N, M;
	int j_left_start, j_left_end, i_down_start, i_down_end;
	int j_right_start, j_right_end, i_up_start, i_up_end;

	int up_neigh, down_neigh, left_neigh, right_neigh;
	int size_lr, size_ud;

	Wave_Solve(Subdomain Local_Grid, MPI_Comm CART_COMM);
	virtual ~Wave_Solve() {};
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
	void    Finalize();
};
#endif