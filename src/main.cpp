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
    Wave_Solve solver;
    solver.Solve_PDE();

    cout << "The computed errors are: " << "\n";
    
    // Compute error
    cout << scientific << right << setw(14)<< solver.Compute_Solution_Error(1);
    cout << scientific << right << setw(14)<< solver.Compute_Solution_Error(2);
    cout << scientific << right << setw(14)<< solver.Compute_Solution_Error(3);
    cout << "\n";
    return 0;
}

