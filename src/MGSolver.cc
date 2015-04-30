
#include <cmath>
#include <iostream>

#include "MGSolver.hh"
#include "Array.hh"
#include "Smoother.hh"

#ifndef PI
#define PI (3.1415)
#endif

MGSolver::MGSolver ( int levels, Smoother & smoother )
	: levels_ (levels)
	, v_grids_(levels, NULL)
	, r_grids_(levels, NULL)
	, h_intervals_(levels, 0)
	, smoother_(smoother)
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

		
#if 0
		std::cout << h_intervals_[i-1] << std::endl;
		v_grids_[i-1]->print();
#endif
		
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

	if (level == 1)
	{
		// SOLVE
		return;
	}
		
	// 1. do pre_smooth gauss seidel iterations on v
	smoother_.smooth_red_black_gauss_seidel_2d ( * v_grids_[level-1],
                                                 * r_grids_[level-1],
                                                 pre_smooth,
                                                 h_intervals_[level-1]);

	std::cout << std::endl;
	v_grids_[level-1]->print();	

	real residual = residual_2d ( * v_grids_[level-1],
                                  * r_grids_[level-1],
                                  h_intervals_[level-1]);

	std::cout << "Residual: " << residual << std::endl;

	// 2. 
}

real MGSolver::residual_2d ( Array & u,
                             Array & f,
							 real h
                           )
{

	real sum = 0.0;	
	real h_2_inv = 1.0 / (h*h);

	int width  = u.getSize(DIM_1D);
    int height = u.getSize(DIM_2D);

	// add up squares of the entries of the residual
	for (int j = 0; j < height; j++) {
	for (int i = 0; i < width; i++)
	{   
		// bottom left corner
		if (i == 0 && j == 0)
		{   
			sum += pow(f(i, j) - h_2_inv * ( 4 * u(i, j) - u(i, j+1)                         - u(i+1, j) ), 2.0);
		}   
		// top left corner
		else if (i == 0 && j == height-1)
		{   
			sum += pow(f(i, j) - h_2_inv * ( 4 * u(i, j)             - u(i, j-1)             - u(i+1, j) ), 2.0);
		}   
		// bottom right corner
		else if (i == width-1 && j == 0)
		{   
			sum += pow(f(i, j) - h_2_inv * ( 4 * u(i, j) - u(i, j+1)             - u(i-1, j)             ), 2.0);
		}   
		// top right corner
		else if (i == width-1 && j == height-1)
		{   
			sum += pow(f(i, j) - h_2_inv * ( 4 * u(i, j)             - u(i, j-1) - u(i-1, j)             ), 2.0);
		}   
		// left side
		else if (i == 0)
		{   
			sum += pow(f(i, j) - h_2_inv * ( 4 * u(i, j) - u(i, j+1) - u(i, j-1)             - u(i+1, j) ), 2.0);
		}   
		// right side
		else if (i == width-1)
		{   
			sum += pow(f(i, j) - h_2_inv * ( 4 * u(i, j) - u(i, j+1) - u(i, j-1) - u(i-1, j)             ), 2.0);
		}   
		// bottom
		else if (j == 0)
		{   
			sum += pow(f(i, j) - h_2_inv * ( 4 * u(i, j) - u(i, j+1)             - u(i-1, j) - u(i+1, j) ), 2.0);
		}   
		// top
		else if (j == height-1)
		{   
			sum += pow(f(i, j) - h_2_inv * ( 4 * u(i, j)             - u(i, j-1) - u(i-1, j) - u(i+1, j) ), 2.0);
		}   
		// inner domain
		else
		{   
			sum += pow(f(i, j) - h_2_inv * ( 4 * u(i, j) - u(i, j+1) - u(i, j-1) - u(i-1, j) - u(i+1, j) ), 2.0);
		}   
	}   
	}

	return sqrt(sum / (real) (width * height));
}
