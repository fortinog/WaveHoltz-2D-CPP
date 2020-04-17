//
// Created by Li-Yin Young on 4/14/20.
//

#include "parameters_setup.h"

#include <cmath>
#include <iostream>



Parameters_setup::Parameters_setup(int nx){
    Nx = 5;
    Ny = 5;

    xb = 0;
    xe = 1.;
    yb = 0;
    ye = 1.;

    hx = (xe-xb)/double(Nx);
    hy = (ye-yb)/double(Ny+1);


    final_time = 3.0;
    CFL = 0.3;


}

void Parameters_setup::time_steps(){
    dt = CFL*std::min(hx,hy);
    nsteps = (int) ceil(final_time/dt);
    dt = final_time/((double) nsteps);
    dt2 = pow(dt,2.0);
    idt2 = 1.0/dt2;
}

void Parameters_setup::grids_2D(Darray2& w, Darray2& wp, Darray2& wm, Darray2& lap){
    w.define(1,-1,Nx+1,0,Ny+1);
    w.set_value(0.0);
    wp.define(1,-1,Nx+1,0,Ny+1);
    wp.set_value(0.0);
    wm.define(1,-1,Nx+1,0,Ny+1);
    wm.set_value(0.0);
    lap.define(1,-1,Nx+1,0,Ny+1);
    lap.set_value(0.0);
}

void Parameters_setup::grids_1D(Darray1& w, Darray1& wp, Darray1& wm, Darray1& lap){
    w.define(1,-1,Nx+1);
    w.set_value(0.1);
    wp.define(1,0,Nx+1);
    wp.set_value(1.0);
    wm.define(1,0,Nx+1);

    wm.set_value(0.0);
    lap.define(1,-1,Nx+1);
    lap.set_value(0.0);
}

void Parameters_setup::coordinates(Darray1& x){

    x.define(1,-1,Nx+1);
    for(int i=-1;i<=Nx+1;i++) x(i) = xb + i*hx;


}

void Parameters_setup::ghost_points1D(Darray1& w, Darray1& x){

    alpha_1 = 0.2;
    alpha_2 = 0.2;

    x_l = x(0) - 0.2*hx;
    x_r = x(Nx) + 0.2*hx;

    w(-1) = boundary_g(x_l)/alpha_1 - w(0)*((1-alpha_1)/alpha_1);
    w(Nx+1) = boundary_g(x_r)/(alpha_2) - ((1-alpha_2)/alpha_2)*w(Nx);

}

double Parameters_setup::boundary_g(double x){

    return x;
    //return sin(x);
}

void Parameters_setup::Initial_value(Darray1& w, Darray1& x){

    for (int i=0; i<= Nx ; i++)
    {
        w(i) = boundary_g(x(i));
    }
}
