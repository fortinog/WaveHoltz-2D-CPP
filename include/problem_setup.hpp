/*
	Here we supply some basic information about the problem, such as 
	number of gridpoints, number of levels, etc.
*/
#ifndef PROBLEMSETUP
#define PROBLEMSETUP

#include <string>   
#include <fstream>
#include <iostream>
#include <math.h>


struct ProblemSetup {
	int numberOfGridpoints;
	int N;
	int M;

	// Step size in either direction
	double hx;
	double hy;

	// Left/Right boundary points
	double x_L;
	double x_R;
	double y_L;
	double y_R;

	// // Value of solution at boundaries
	// double ua;
	// double ub;

	// Specify a type for MMS (if used)
	int twilightType;
	double kx, ky, kt;
	double x0, y0, t0;
	double cx[5],cy[5],ct[5];
	ProblemSetup();
	double Level_Set(const double x, const double y);
	double c2(const double x, const double y);
};
#endif