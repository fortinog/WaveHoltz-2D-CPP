//
// Created by Li-Yin Young on 4/14/20.
//

#include "Wave_Solve1D.h"

Wave_Solve1D::Wave_Solve1D(double m_dt, double m_hx, int nxs) {

    dt = m_dt;
    hx = m_hx;
    Nx = nxs;

}



double Wave_Solve1D::boundary_g(double x){

    return 1;
    //return sin(x);
}

void Wave_Solve1D::laplacian(Darray1& w, Darray1& lap){

    double ihx2 = 1/hx/hx;;
    for (int i=0;i <= Nx;i++ ){
        lap(i) = 0.5 * ihx2 * (w(i+1) - 2*w(i) + w(i-1));
    }
}

void Wave_Solve1D::advance(Darray1& wm, Darray1& w, Darray1& wp, Darray1& lap){
    double dt2 = dt*dt;

    for(int i = 1 ; i <= Nx-1; i++)
        wp(i) = 2.0*w(i)-wm(i) + dt2*lap(i);

}

void Wave_Solve1D::swap_time(Darray1& wm, Darray1& w, Darray1& wp){


    wm.copy(w);

    for(int i = 1 ; i <= Nx-1; i++)
        w(i) = wp(i);

}

