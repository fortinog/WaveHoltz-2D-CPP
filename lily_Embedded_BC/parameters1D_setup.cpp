//
// Created by Li-Yin Young on 4/14/20.
//

#include "parameters1D_setup.h"

#include <cmath>
#include <iostream>



Parameters1D_setup::Parameters1D_setup(int nx, Darray1& domain){

    x = domain;

    Nx = nx;
    Ny = nx;

    x_beg = 0;
    x_end = 1.;
    y_beg = 0;
    y_end = 1.;

    alpha_1 = 0.2;
    alpha_2 = 0.2;

    hx = (x_end-x_beg)/double(Nx);
    hy = (y_end-y_beg)/double(Ny+1);


    final_time = 3.0;
    CFL = 0.3;

    coordinates(x);


}

void Parameters1D_setup::time_steps(){
    dt = CFL*std::min(hx,hy);
    nsteps = (int) ceil(final_time/dt);
    dt = final_time/((double) nsteps);
    dt2 = pow(dt,2.0);
    idt2 = 1.0/dt2;
}



void Parameters1D_setup::grids(Darray1& wm, Darray1& w, Darray1& wp, Darray1& lap){
    wm.define(1,-1,Nx+1);
    wm.set_value(0.0);

    w.define(1,-1,Nx+1);
    w.set_value(0.1);

    wp.define(1,-1,Nx+1);
    wp.set_value(0.0);

    lap.define(1,-1,Nx+1);
    lap.set_value(0.0);
}

void Parameters1D_setup::coordinates(Darray1& x){

    x.define(1,-1,Nx+1);
    for(int i=-1;i<=Nx+1;i++) x(i) = x_beg + i*hx;


}

void Parameters1D_setup::ghost_points(Darray1& w){

    w(-1) = w(0)/alpha_1 - w(1)*((1-alpha_1)/alpha_1);
    w(Nx+1) = w(Nx)/(alpha_2) - ((1-alpha_2)/alpha_2)*w(Nx-1);

}

void Parameters1D_setup::boundary_points(Darray1& w){

    x_l = x(1) - 0.2*hx;
    x_r = x(Nx-1) + 0.2*hx;

    w(0) = sin(x_l);
    w(Nx) = sin(x_r);
}

void Parameters1D_setup::initial_value(Darray1& w){

    w.set_value(0.1);

}

