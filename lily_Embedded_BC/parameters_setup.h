//
// Created by Li-Yin Young on 4/14/20.
//

#ifndef LILY_EMBEDDED_BC_PARAMETERS_SETUP_H
#define LILY_EMBEDDED_BC_PARAMETERS_SETUP_H
#include "Darrays.h"
#include "Wave_Solve1D.h"


class Parameters_setup {
public:
    Parameters_setup(int nx);


    //grids
    int Nx, Ny;
    double xb, xe, yb, ye;
    double hx,hy;
    double alpha_1, alpha_2, x_l, x_r;

    //time
    double final_time, dt, dt2, idt2;
    int nsteps;

    double CFL;


    void time_steps();
    void grids_2D(Darray2& w, Darray2& wp, Darray2& wm, Darray2& lap);
    void grids_1D(Darray1& w, Darray1& wp, Darray1& wm, Darray1& lap);
    void coordinates(Darray1& x);
    void ghost_points1D(Darray1& w, Darray1& x);
    double boundary_g(double x);
    void Initial_value(Darray1& w, Darray1& x);


};


#endif //LILY_EMBEDDED_BC_PARAMETERS_SETUP_H
