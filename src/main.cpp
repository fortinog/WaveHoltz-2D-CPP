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

// Main program
int main()
{
    Darray2 u, v, res;
    double err;
    vector<int> size_vec(5);
    int start_size = 32;
    int N,M;

    // Set up sizes
    for(int i = 0;i<size_vec.size();i++) size_vec[i] = start_size*int(pow(2.0,i));

    cout << "Below are timings for addition. Increasing number of rows left to right, increasing number of columns top to bottom. \n\n";
    cout << right << setw(8) << " ";
    for(int i = 0;i<size_vec.size();i++) cout << right << setw(14) << size_vec[i];
    cout << "\n";
    
    for(int i = 0;i<size_vec.size();i++){
        N = size_vec[i];
        cout << right << setw(8) << N;
        for(int j = 0;j<size_vec.size();j++){
            M = size_vec[j];
            u.define(1,M,1,N);
            v.define(1,M,1,N);
            res.define(1,M,1,N);
            u.set_value(1.0);
            v.set_value(2.0);

            auto start = high_resolution_clock::now();
            for(int k=1;k<=M;k++) 
                for(int l=1;l<=N;l++) 
                    res(k,l) = u(k,l)+v(k,l);

            auto stop = high_resolution_clock::now(); 
            auto duration = duration_cast<microseconds>(stop - start); 

            cout << scientific << right << setw(14)
                 << duration.count();
        }
        cout << "\n";
    }

    cout << "\n\n";
    cout << "Below are timings for multiplication. Increasing number of rows left to right, increasing number of columns top to bottom. \n\n";
    cout << right << setw(8) << " ";
    for(int i = 0;i<size_vec.size();i++) cout << right << setw(14) << size_vec[i];
    cout << "\n";
    
    for(int i = 0;i<size_vec.size();i++){
        N = size_vec[i];
        cout << right << setw(8) << N;
        for(int j = 0;j<size_vec.size();j++){
            M = size_vec[j];
            u.define(1,M,1,N);
            v.define(1,M,1,N);
            res.define(1,M,1,N);
            u.set_value(1.0);
            v.set_value(2.0);

            auto start = high_resolution_clock::now();
            for(int k=1;k<=M;k++) 
                for(int l=1;l<=N;l++) 
                    res(k,l) = u(k,l)*v(k,l);

            auto stop = high_resolution_clock::now(); 
            auto duration = duration_cast<microseconds>(stop - start); 

            cout << scientific << right << setw(14)
                 << duration.count();
        }
        cout << "\n";
    }

    cout << "\n\n";
    cout << "Now we compute the Laplacian. \n";
    cout << right << setw(8) << "N" << right << setw(10)
            << "L_1"  << right << setw(12) << "L_2"
             << right << setw(16) << "L_inf";

    cout << "\n";
    Wave_Solve solver;
    for(int i=0;i<6;i++){
        cout << right << setw(8) << solver.setup.N;
        solver.Set_Initial_Data();
        solver.Compute_Laplacian();
        err = solver.Compute_Laplacian_Error(1);
        cout << scientific << right << setw(14)<< solver.Compute_Laplacian_Error(1);
        cout << scientific << right << setw(14)<< solver.Compute_Laplacian_Error(2);
        cout << scientific << right << setw(14)<< solver.Compute_Laplacian_Error(3);
        cout << "\n";
        solver.Refine_Grid();
    }

    solver.Compute_Mask();
    
    
    return 0;
}

