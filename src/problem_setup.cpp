#include"problem_setup.hpp"

/************************************************************/
/*	    HERE IS WHERE WE EDIT THE PROBLEM PARAMETERS	    */
/************************************************************/

static const int	defaultN = 20;
static const int	defaultM = 20;

static const double	defaultXl					= -1.0;
static const double	defaultXr					= 1.0;
static const double	defaultYl					= -1.0;
static const double	defaultYr					= 1.0;

static const double	defaultCFL					= 0.5;
static const double	defaultFinalTime			= 1.0;

static const int	defaultTwilightType			= 1;
static const double	defaultKx					= 2*M_PI;
static const double	defaultKy					= 2*M_PI;
static const double	defaultKt					= 2*sqrt(2)*M_PI;

static const double	defaultX0					= 0.0;
static const double	defaultY0					= 0.0;
static const double	defaultT0					= 0.0*M_PI;

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
	CFL(defaultCFL),
	final_time(defaultFinalTime),
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
		// hx = (x_R-x_L)/double(N+1);
		// hy = (y_R-y_L)/double(M+1);
		hx = (x_R-x_L)/double(N-1);
		hy = (y_R-y_L)/double(M-1);

		// Calculate time step size
		dt     = CFL*std::min(hx,hy);
		nsteps = (int) ceil(final_time/dt);
		dt     = final_time/((double) nsteps);
		dt2    = pow(dt,2.0);
		idt2   = 1.0/dt2;
	}

// The level set function defining the computational boundary of the domain
double ProblemSetup::Level_Set(const double x, const double y){
	
    // [x_L,x_R]x[y_L,y_R]
	// double xval, yval, l_val;
	// xval = std::min(x-(*this).x_L,(*this).x_R-x);
	// yval = std::min(y-(*this).y_L,(*this).y_R-y);
	// l_val = -std::min(xval,yval);
	// return l_val;
    
	// double xval, yval, l_val;
	// xval = std::min(x+0.5,0.5-x);
	// yval = std::min(y+0.5,0.5-y);
	// l_val = -std::min(xval,yval);
	// return l_val;
    

    // Circle of radius 1
    double r = pow(x,2) + pow(y,2);
    double l_val = r-pow(0.8,2);
    return l_val;
}

// Speed of sound in the medium
double ProblemSetup::c2(const double x, const double y){
	// Normalized constant coefficient
	return 1.0;
}

// Routines for root finding
// AAL Start
// Use secant method to find the distance to the boundary
// xn+1 = xn - f(xn-1)(xn-1 - xn-2)/(f(xn-1) - f(xn-2))
double ProblemSetup::Dist_to_bdry(const double x, const double y,int dir){
    double d, tol,xold,xn,xnew,yold,yn,ynew;
    double num, denom, alpha;
    int k;
    tol = 10e-10; // Tolerance
    // Inital Guesses are midpoint of domain and Ghost pt
    if (dir == 1){
    	xold = 0.5;  // Midpoint of x interval
    	yold = y;  // y val
    }
    else if (dir == 2){
    	xold = x;    // x val
    	yold = 0.5;  // Midpoint of y interval
    }
    xn = x;
    yn = y;
    int Count = 0;
    while (abs(Level_Set(xn,yn)) > tol){
        num = Level_Set(xn,yn);
        denom = Level_Set(xn,yn) - Level_Set(xold,yold);
        alpha = num/denom;
        xnew = xn - alpha*(xn-xold);
        ynew = yn - alpha*(yn-yold);
        // Update everything
        xold = xn;
        yold = yn;
        xn = xnew;
        yn = ynew;
        Count++;
    }
    d = sqrt((xn-x)*(xn-x)+(yn-y)*(yn-y));
    return d;
}
// AAL Finish
