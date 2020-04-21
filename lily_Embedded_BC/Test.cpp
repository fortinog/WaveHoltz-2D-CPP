//
// Created by Li-Yin Young on 4/17/20.
//

#include <iostream>
#include <iomanip>
#include <cmath>
#include <chrono>
#include <vector>
#include "stddef.h"
#include "Test.h"



Test::Test(int nxs){
    Nx = nxs;

}

void Test::print_matrix1D(int start, int end, Darray1& x){

    std:: cout << "\n";
    for(int i=start;i<=end;i++) {
        std::cout << x(i) << " ";
    }

}

void Test::check_ghost_pts(Darray1& w, Darray1& x, int itr){
    double ans_l1, ans_l2, ans_r1, ans_r2;

    Parameters_setup setup(Nx,x);

    ans_l1 = (w(0)-w(-1))/(setup.x_l-x(-1));
    ans_l2 = (w(-1)-w(0))/(x(-1)-x(0));
    ans_r1 = (w(0)-w(Nx))/(setup.x_r-x(Nx));
    ans_r2 = (w(Nx+1)-w(Nx))/(x(Nx+1)-x(Nx));

    if (ans_l1 != ans_l2)
    {std::cout << "\nWrong ghost points at left in " << itr << " step." << "\n";
        std::cout << "Left: " << ans_l1 << " " << ans_l2 << "\n";
        std::cout << w(-1) << " " << w(0);}

    if (ans_r1 != ans_r2)
    {
        std::cout << "\nWrong ghost points at right in " << itr << " step." << "\n";
        std::cout << "Right: " << ans_r1 << " " << ans_r2 << "\n";
    }
}