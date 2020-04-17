#include <iostream>
#include <iomanip>
#include <cmath>
#include <chrono>
#include <vector>
#include "stddef.h"
#include "silo.h"
#include "Darrays.h"
#include "parameters_setup.h"
#include "Wave_Solve1D.h"
#include "Test.h"



using namespace std;
using namespace std::chrono;

int main()
{

    Darray1 x;
    Darray1 w, wp, wm, lap;
    double energy_old = 1e10;
    double energy;

    int Nx= 5;
    int n_steps = 5;

    Parameters_setup setup(Nx);
    Wave_Solve1D solve(setup.dt, setup.hx, Nx);
    Test test(Nx);
    setup.grids_1D(w, wp, wm, lap);
    setup.coordinates(x);


    setup.ghost_points1D(w, x);

    // Test initialization of ghost points
    //test.print_matrix1D(-1,Nx+1,x);


    //Test updating ghost points
    setup.Initial_value(w,x);

    for (int itr = 0 ; itr<= n_steps; itr++ )
    {
        setup.ghost_points1D(w,x);
        cout << "x_r: " << setup.x_r;
        cout << itr <<  " w(x): ";
        for(int i=-1;i<=Nx+1;i++) {
            cout << w(i) << " ";
        }




        cout << "\n";
    }

    return 0;
}

