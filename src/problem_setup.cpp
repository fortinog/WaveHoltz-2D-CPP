#include"problem_setup.hpp"

/************************************************************/
/*	    HERE IS WHERE WE EDIT THE PROBLEM PARAMETERS	    */
/************************************************************/

static const int	defaultN = 60;
static const int	defaultM = 60;


static const int    defaultIter      = 100;
static const double	defaultTol	     = 1e-10;
static const int    defaultPrintFreq = 1;


static const double	defaultXl					= -1.0;
static const double	defaultXr					= 1.0;
static const double	defaultYl					= -1.0;
static const double	defaultYr					= 1.0;

static const double	defaultCFL					= 0.5;

static const int	defaultTwilightType			= 1;
static const double	defaultKx					= 0.5*M_PI;
static const double	defaultKy					= M_PI;
static const double	defaultKt					= 0.0;
static const double defaultomega                = 3.8;

static const double	defaultX0					= 0.0;
static const double	defaultY0					= 0.0;
static const double	defaultT0					= 0.5*M_PI;

static const double	defaultCx[5]				= {1,1,0,0,0};
static const double	defaultCy[5]				= {1,1,0,0,0};
static const double	defaultCt[5]				= {1,0,0,0,0};

/************************************************************/
/*	         END USER SPECIFIED PROBLEM PARAMETERS	        */
/************************************************************/

// Constructor
ProblemSetup::ProblemSetup():
	N(defaultN),
	M(defaultM),
	CFL(defaultCFL),
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
	t0(defaultT0),
    omega(defaultomega),
    cg_iter(defaultIter),
    cg_tol(defaultTol),
    print_freq(defaultPrintFreq)
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

		// Calculate time step size
        final_time = 2*M_PI/omega;
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
	// xval = std::min(x+0.0,1.0-x);
	// yval = std::min(y+0.0,1.0-y);
	// l_val = -std::min(xval,yval);
	// return l_val;
    

    // Circle of radius 1
   double r = pow(x,2) + pow(y,2);
   double l_val = r-pow(1.0,2);
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
double ProblemSetup::Dist_to_bdry(const double x_in, const double x_out, 
								  const double y_in, const double y_out,
								  int dir){
    double d, tol,xold,xn,xnew,yold,yn,ynew;
    double num, denom, alpha;
    int k;
	int maxIter = 30;

    tol = 10e-13; // Tolerance

    // Inital Guesses are midpoint of domain and Ghost pt
    if (dir == 1){
    	xold = x_out;  // exterior point
    	yold = y_in;  // y val
    }
    else if (dir == 2){
    	xold = x_in;    // x val
    	yold = y_out;  // exterior point
    }
    xn = x_in;
    yn = y_in;
    int Count = 0;
    while ((abs(Level_Set(xn,yn)) > tol) && (Count <= maxIter)){
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
    d = sqrt((xn-x_in)*(xn-x_in)+(yn-y_in)*(yn-y_in));

    return d;
}
// AAL Finish
