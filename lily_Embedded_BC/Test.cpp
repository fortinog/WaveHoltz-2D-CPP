//
// Created by Li-Yin Young on 4/17/20.
//

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

Test::Test(int nxs){
    double Nx = nxs;
}

void Test::print_matrix1D(int start, int end, Darray1& x){

    for(int i=start;i<=end;i++) {
        std::cout << i << " x: " << x(i) << " ";
    }
}

void Test::check_ghost_pts(){
    double ans_l1, ans_l2, ans_r1, ans_r2;

    ans_l1 = (setup.boundary_g(setup.x_l)-w(-1))/(setup.x_l-x(-1));
    ans_l2 = (w(-1)-w(0))/(x(-1)-x(0));
    ans_r1 = (setup.boundary_g(setup.x_r)-w(Nx))/(setup.x_r-x(Nx));
    ans_r2 = (w(Nx+1)-w(Nx))/(x(Nx+1)-x(Nx));

    if (ans_l1 != ans_l2)
    {cout << "\nWrong ghost points at left in " << itr << " step." << "\n";
        cout << "Left: " << ans_l1 << " " << ans_l2 << "\n";
        cout << w(-1) << " " << w(0);}

    if (ans_r1 != ans_r2)
    {cout << "\nWrong ghost points at right in " << itr << " step." << "\n";
        cout << "Right: " << ans_r1 << " " << ans_r2 << "\n";}
}