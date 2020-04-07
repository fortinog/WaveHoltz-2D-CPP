// Twilight class
#ifndef TWILIGHT
#define TWILIGHT

#define _USE_MATH_DEFINES
 
#include <cmath>
#include <iostream>
#include "problem_setup.hpp"

using namespace std;

class Twilight {
public:
	int type;
	double kx, x0, ky, y0, kt, t0;
	double cx[5], cy[5], ct[5];
	Twilight(){};
	Twilight(const ProblemSetup& setup);
	virtual ~Twilight() {};
	double trigTwilight(const int dx,const double x, const int dy,
						const double y, const int dt,const double t);
	double polyTwilight(const int dx,const double x, const int dy,
						const double y, const int dt,const double t);
};
#endif