
#include <cmath>
#include <iostream>

#include "MGSolver.hh"
#include "Array.hh"

#ifndef PI
#define PI (3.1415)
#endif

MGSolver::MGSolver ( int levels )
	: levels_ (levels)
	, v_grids_(levels, NULL)
	, r_grids_(levels, NULL)
	, h_intervals_(levels, 0)
{

	 //v_grids_(levels, 0);
	 //r_grids_(levels, 0);

	// allocate memory for v and r
	for (int i = 1; i <= levels; i++)
	{
		std::cout << i << std::endl;

		v_grids_[i-1] = new Array ( std::pow(2, i) - 1, std::pow(2, i) - 1);
		r_grids_[i-1] = new Array ( std::pow(2, i) - 1, std::pow(2, i) - 1);
		h_intervals_[i-1] = 1.0 / std::pow(2, i);

		std::cout << h_intervals_[i-1] << std::endl;
		v_grids_[i-1]->print();
		
	}	

}

void MGSolver::initialize_assignment_01 ()
{

	// initialize rhs r, by ADDING boundary values
	// this comes from:
	// - u'' = f  => - (u_i+1j - 2u_ij + u_i-1j)/h^2 = f => ADD boundaries

	// bc is sin(x*pi) sinh(y*pi) which is only != 0 at the top boundary (y = 1)

	Array * finest_grid = r_grids_.back();
	real h              = h_intervals_.back();

	for (int col = 0;
         col < finest_grid->getSize(DIM_1D); 
         col++)
	{
		finest_grid->operator()(col, finest_grid->getSize(DIM_2D)-1) = sin(PI * (real) (col+1) * h) * sinh(PI * h) / (h*h);
	}
	
	finest_grid->print();

}

void MGSolver::v_cycle     ( int pre_smooth,
                             int post_smooth,
                             int times        
                           )
{
	for (int i = 0; i < times; i++)
	{
		v_cycle_pvt (pre_smooth, post_smooth, levels_);
	}
}


void MGSolver::v_cycle_pvt ( int pre_smooth,
                             int post_smooth,
                             int level         // solve when level == 1
                           )
{

	(void) pre_smooth;
	(void) post_smooth;
	(void) level;
	// 1. do pre_smooth gauss seidel iterations on v

	// 2. 
}
