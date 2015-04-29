
#include <cmath>
#include <iostream>

#include "MGSolver.hh"
#include "Array.hh"

MGSolver::MGSolver ( int levels )
	: levels_ (levels)
	, v_grids_(levels, NULL)
	, r_grids_(levels, NULL)
{

	 //v_grids_(levels, 0);
	 //r_grids_(levels, 0);

	// allocate memory for v and r
	for (int i = 1; i <= levels; i++)
	{
		std::cout << i << std::endl;

		v_grids_[i-1] = new Array ( std::pow(2, i) - 1, std::pow(2, i) - 1);
		r_grids_[i-1] = new Array ( std::pow(2, i) - 1, std::pow(2, i) - 1);
		
		v_grids_[i-1]->print();

	}	


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
