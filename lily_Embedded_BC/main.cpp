#include <iostream>
#include <iomanip>
#include <cmath>
#include <chrono>
#include <vector>
#include "stddef.h"
#include "silo.h"
#include "Darrays.h"
#include "parameters1D_setup.h"
#include "Wave_Solve1D.h"
#include "Test.h"
#include "parameters2D_setup.h"



using namespace std;
using namespace std::chrono;

int main()
{

    Darray1 x, y;
    //Darray1 w, wp, wm, lap;
    Darray2 w, wp, wm, lap;


    int Nx= 5;
    int Ny= 5;
    int n_steps = 5;

    double energy_old = 1e10;
    double energy;


    Parameters2D_setup setup(Nx, Ny, x, y);
    Wave_Solve1D solve(setup.dt, setup.hx, Nx);
    Test Test(Nx);

    //setup.coordinates();
    setup.grids(wm, w, wp, lap);
    setup.initial_value(wm,w,wp);
    setup.boundary_points(w);
    setup.ghost_points(w);




    for (int itr = 0 ; itr<= n_steps; itr++)
    {
        setup.ghost_points(w);

        solve.laplacian(w,lap);


        solve.advance(wm, w, wp, lap);


        setup.boundary_points(w);
        Test.print_matrix1D(-1,Nx+1,w);
    }

    return 0;
}

