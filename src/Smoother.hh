#ifndef SMOOTHER_H
#define SMOOTHER_H

#include <vector>

#include "Array.hh"
#include "Types.hh"

class Smoother
{
    public:
   

		void smooth_red_black_gauss_seidel_2d ( Array & u,    // modify this array
                                                Array & f,    // rhs
                                                int times,    // number of sweeps
                                                real h,       // spacing
												bool finest_grid
                                           );

		void update_neumann( Array & u, real h);
 

};

#endif

