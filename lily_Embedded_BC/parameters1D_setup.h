//
// Created by Li-Yin Young on 4/14/20.
//

#ifndef LILY_EMBEDDED_BC_PARAMETERS1D_SETUP_H
#define LILY_EMBEDDED_BC_PARAMETERS1D_SETUP_H
#include "Darrays.h"
#include "Wave_Solve1D.h"


class Parameters1D_setup {
public:
    Parameters1D_setup(int nx, Darray1& domain);
    Darray1 x;


    //grids

    int Nx, Ny;
    double x_beg, x_end, y_beg, y_end;
    double hx,hy;
    double alpha_1, alpha_2, x_l, x_r;

    //time
    double final_time, dt, dt2, idt2;
    int nsteps;

    double CFL;


    void time_steps();

    void grids(Darray1& wm, Darray1& w, Darray1& wp, Darray1& lap);

    void ghost_points(Darray1& w);
    void boundary_points(Darray1& w);
    void initial_value(Darray1& w);


private:
    void coordinates(Darray1& x);


};


#endif //LILY_EMBEDDED_BC_PARAMETERS1D_SETUP_H
