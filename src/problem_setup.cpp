#include"problem_setup.hpp"

/************************************************************/
/*	    HERE IS WHERE WE EDIT THE PROBLEM PARAMETERS	    */
/************************************************************/

static const int	defaultN = 30;
static const int	defaultM = 30;

static const double	defaultXl					= 0.0;
static const double	defaultXr					= 1.0;
static const double	defaultYl					= 0.0;
static const double	defaultYr					= 1.0;

static const int	defaultTwilightType			= 1;
static const double	defaultKx					= M_PI;
static const double	defaultKy					= M_PI;
static const double	defaultKt					= M_PI;

static const double	defaultX0					= 0.0;
static const double	defaultY0					= 0.0;
static const double	defaultT0					= 0.0;

static const double	defaultCx[5]				= {1,0,0,0,0};
static const double	defaultCy[5]				= {0,0,0,0,1};
static const double	defaultCt[5]				= {1,0,0,0,0};

/************************************************************/
/*	         END USER SPECIFIED PROBLEM PARAMETERS	        */
/************************************************************/

// Constructor
ProblemSetup::ProblemSetup():
	N(defaultN),
	M(defaultM),
	x_L(defaultXl),
	x_R(defaultXr),
	y_L(defaultYl),
	y_R(defaultYr),
	twilightType(defaultTwilightType),
	kx(defaultKx),
	ky(defaultKy),
	kt(defaultKt),
	x0(defaultX0),
	y0(defaultY0),
	t0(defaultT0)
	{ 
		int i;
		for(i=0;i<5;i++){
			cx[i] = defaultCx[i];
			cy[i] = defaultCy[i];
			ct[i] = defaultCt[i];
		}

		// Include the physical boundary for now
		hx = (x_R-x_L)/double(N-1);
		hy = (y_R-y_L)/double(M-1);
	}

// The level set function defining the computational boundary of the domain
double ProblemSetup::Level_Set(const double x, const double y){
	// [-1,1]x[-1,1]
	double xval, yval, l_val;
	xval = std::min(x-(*this).x_L,(*this).x_R-x);
	yval = std::min(y-(*this).y_L,(*this).y_R-y);
	l_val = -std::min(xval,yval);
	return l_val;
	// The box [x_L,x_R]x[y_L,y_R]
	// return ((x_L<x && x<x_R) && (y_L<y && y<y_R))? 1:0;
}

// Speed of sound in the medium
double ProblemSetup::c2(const double x, const double y){
	// Normalized constant coefficient
	return 1.0;
}