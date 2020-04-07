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
#include <math.h>
#include "problem_setup.hpp"
#include "Darray1.h"
#include "Darray2.h"
#include "Iarray2.h"
#include "twilight_2d.h"
#include <iomanip>
using namespace std;

class Wave_Solve {
public:
	ProblemSetup setup;
	Darray2 w, lap, wm, wp;
	Iarray2 mask;
	Darray1 x, y;
	double t;

	Wave_Solve();
	virtual ~Wave_Solve() {};
	void Compute_Mask();
	void Enforce_BC(Darray2 &v);
	void Compute_Laplacian();
	void Set_Initial_Data();
	double Compute_Laplacian_Error(const int flag);
	double Compute_Solution_Error(const int flag);
	void Refine_Grid();
	void Time_Step();
	void Taylor_Expand();
	void Clear_Data();
	double Compute_Energy();
	double forcing(const double x, const double y, const double t);
	void Solve_PDE();
};
#endif