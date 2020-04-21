//
// Created by Li-Yin Young on 4/19/20.
//

#ifndef LILY_EMBEDDED_BC_PARAMETERS2D_SETUP_H
#define LILY_EMBEDDED_BC_PARAMETERS2D_SETUP_H

#include "Darrays.h"


class Parameters2D_setup {
public:
    Parameters2D_setup(int nx, int ny, Darray1& domain_x, Darray1& domain_y);

    //grids
    Darray1 x, y;
    int Nx, Ny;
    double x_beg, x_end, y_beg, y_end;
    double hx,hy;
    double alpha_1, x_left, x_right;
    double alpha_2, y_top, y_bottom;

    //time
    double final_time, dt, dt2, idt2;
    int nsteps;

    double CFL;


    void grids(Darray2& w, Darray2& wp, Darray2& wm, Darray2& lap);
    void initial_value(Darray2& wm, Darray2& w, Darray2& wp );

    void boundary_points(Darray1& w, double t);
    double boundary_g(double t, double x, double y);




private:
    void coordinates(Darray1& domain_x, Darray1& domain_y);


};


#endif //LILY_EMBEDDED_BC_PARAMETERS2D_SETUP_H
