#include <iostream>
#include <iomanip>
#include <cmath>
#include <chrono> 
#include <vector>
#include "Darray.h"
#include "Darray2.h"
#include "problem_setup.hpp"
#include "wave_solve.hpp"
using namespace std;
using namespace std::chrono;

int main()
{
    Darray2 u, v, res;
    double err;
    vector<int> size_vec(5);
    int start_size = 32;
    int N,M;
    int nsteps;
    double energy = 0.0;
    double energy_old = 0.0;

    Wave_Solve solver;
    nsteps = solver.setup.nsteps;
    solver.Solve_PDE();

    cout << "The computed errors are: " << "\n";
    // Compute error
    cout << scientific << right << setw(14)<< solver.Compute_Solution_Error(1);
    cout << scientific << right << setw(14)<< solver.Compute_Solution_Error(2);
    cout << scientific << right << setw(14)<< solver.Compute_Solution_Error(3);
    cout << "\n";
    return 0;
}

