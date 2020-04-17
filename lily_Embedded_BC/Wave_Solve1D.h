//
// Created by Li-Yin Young on 4/14/20.
//

#ifndef LILY_EMBEDDED_BC_WAVE_SOLVE1D_H
#define LILY_EMBEDDED_BC_WAVE_SOLVE1D_H


//#include "parameters_setup.h"
#include <cmath>
#include "Darrays.h"
#include "parameters_setup.h"

class Wave_Solve1D {
public:
    Wave_Solve1D(double m_dt, double m_hx, int nxs);

    int Nx;
    double dt, hx, ihx2;




    double boundary_g(double x);
    void ghost_points1D(Darray1& u, Darray1& x);
    void laplacian(Darray1& w, Darray1& lap);
    void advance(Darray1& up, Darray1& u, Darray1& um, Darray1& lap);


};


#endif //LILY_EMBEDDED_BC_WAVE_SOLVE1D_H
