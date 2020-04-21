//
// Created by Li-Yin Young on 4/19/20.
//

#include "parameters2D_setup.h"
#include "parameters1D_setup.h"
#include <cmath>
Parameters2D_setup::Parameters2D_setup(int nx, int ny, Darray1& domain_x, Darray1& domain_y){

    x = domain_x;
    y = domain_y;


    Nx = nx;
    Ny = ny;

    x_beg = 0;
    x_end = 1.;
    y_beg = 0;
    y_end = 1.;

    alpha_1 = 0.2;
    alpha_2 = 0.2;

    CFL = 0.3;

    hx = (x_end-x_beg)/double(Nx);
    hy = (y_end-y_beg)/double(Ny+1);

    final_time = 3.0;
    CFL = 0.3;

    coordinates(x,y);
}

void Parameters2D_setup::grids(Darray2& w, Darray2& wp, Darray2& wm, Darray2& lap){
    w.define(1,-1,Nx+1,0,Ny+1);
    w.set_value(0.0);
    wp.define(1,-1,Nx+1,0,Ny+1);
    wp.set_value(0.0);
    wm.define(1,-1,Nx+1,0,Ny+1);
    wm.set_value(0.0);
    lap.define(1,-1,Nx+1,0,Ny+1);
    lap.set_value(0.0);
}

void Parameters2D_setup::initial_value(Darray2& wm, Darray2& w, Darray2& wp ){

    double y2;
    Darray2 ut(0,Nx,0,Ny);
    for(int j = 0; j<=Ny; j++){
        y2 = y(j)*y(j);
        for(int i = 0; i<= Nx; i++){
            wp(i,j) = w(i,j) = exp(-36.0*(x(i)*x(i)+y2));
            ut(i,j) = 0.0;
            wm(i,j) = w(i,j) - dt*ut(i,j);
        }
    }
}

void Parameters2D_setup::boundary_points(Darray1& w, double t){

    int i,j;

    x_left = x(1) - alpha_1*hx;
    x_right = x(Nx-1) + alpha_1*hx;

    y_bottom = y(1) - alpha_2*hy;
    y_top = y(Nx-1) + alpha_2*hy;

    // Left boundary
    for(j=1;j<=Ny;j++){
        w(i,j) = boundary_g(t,x_left,y(j));
    }

}

double Parameters2D_setup::boundary_g(double t, double xi, double yj){

    return sin(t*xi*yj);
}

void Parameters2D_setup::coordinates(Darray1& x, Darray1& y){

    x.define(1,-1,Nx+1);
    for(int i=-1;i<=Nx+1;i++) x(i) = x_beg + i*hx;

    x_left = x(1) - alpha_1*hx;
    x_right = x(Nx-1) + alpha_1*hx;

    x(0) = x_left;

    y_bottom = y(1) - alpha_2*hy;
    y_top = y(Nx-1) + alpha_2*hy;


}

