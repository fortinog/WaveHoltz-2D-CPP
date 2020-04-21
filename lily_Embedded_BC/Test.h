//
// Created by Li-Yin Young on 4/17/20.
//

#ifndef LILY_EMBEDDED_BC_TEST_H
#define LILY_EMBEDDED_BC_TEST_H
#include "parameters1D_setup.h"


class Test{
public:

    double Nx;

    Test(int nxs);
    void print_matrix1D(int start, int end, Darray1& x);
    void check_ghost_pts(Darray1& w, Darray1& x, int itr);


};



#endif //LILY_EMBEDDED_BC_TEST_H
