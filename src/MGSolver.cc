
#include <cmath>
#include <iostream>
#include <cmath>
#include <math.h>

#include "MGSolver.hh"
#include "Array.hh"
#include "Smoother.hh"

#ifndef PI
#define PI (M_PI)
#endif

MGSolver::MGSolver ( int levels, Smoother & smoother )
	: levels_ (levels)
	, v_grids_(levels, NULL)
	, r_grids_(levels, NULL)
	, tmp_grids_(levels, NULL)
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
		tmp_grids_[i-1] = new Array ( std::pow(2, i) - 1, std::pow(2, i) - 1);
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
		std::cout << sin(PI * (real) (col+1) * h) * sinh(PI * h) << std::endl;
	}
	

}

void MGSolver::initialize_random ()
{
	v_grids_.back()->fillRandom();
}

void MGSolver::v_cycle     ( int pre_smooth,
                             int post_smooth,
                             int times        
                           )
{


	for (int i = 0; i < times; i++)
	{
		v_cycle_pvt (pre_smooth, post_smooth, levels_);

		std::cout << std::endl;

		std::cout << "Cycle no " << i + 1 << std::endl;
		
		v_grids_.back()->print();

		std::cout << std::endl;

		real residual = residual_2d ( * v_grids_.back(),
									  * r_grids_.back(),
									  h_intervals_.back());

		std::cout << "Residual (cylcle no " << i + 1 << "): " << residual << std::endl;
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
		//std::cout << v_grids_[level-1]->operator()(0, 0) << std::endl;
		v_grids_[level-1]->operator()(0, 0) = r_grids_[level-1]->operator()(0, 0) * h_intervals_[level-1] * h_intervals_[level-1] * 0.25;
	std::cout << "V GRID" << std::endl;
	v_grids_[level-1]->print();	
	std::cout << "R GRID" << std::endl;
	r_grids_[level-1]->print();	
	std::cout << "T GRID" << std::endl;
	tmp_grids_[level-1]->print();	
		
		return;
	}

	v_grids_[level-2]->fill(0.0);
	tmp_grids_[level-1]->fill(0.0);
		
	// 1. perform pre_smooth gauss seidel iterations on v
	smoother_.smooth_red_black_gauss_seidel_2d ( * v_grids_[level-1],
                                                 * r_grids_[level-1],
                                                 pre_smooth,
                                                 h_intervals_[level-1]);



	// 2. calculate coarser right hand side
	compose_right_hand_side ( * v_grids_[level-1],
                              * r_grids_[level-1],
                              * r_grids_[level-2],
                              level,
                              h_intervals_[level-1]);

	// tmp_grids_[level-2]->print();

	// 3.recursive call
	v_cycle_pvt ( pre_smooth, post_smooth, level-1 );


	// 4. correction
	error_correction ( * v_grids_[level-1],
                       * v_grids_[level-2],
                       level
                     );
	
	

	// 5. post smoothing
	smoother_.smooth_red_black_gauss_seidel_2d ( * v_grids_[level-1],
                                                 * r_grids_[level-1],
                                                 post_smooth,
                                                 h_intervals_[level-1]);
	std::cout << "V GRID" << std::endl;
	v_grids_[level-1]->print();	
	std::cout << "R GRID" << std::endl;
	r_grids_[level-1]->print();	
	std::cout << "T GRID" << std::endl;
	tmp_grids_[level-1]->print();	
	

}

void MGSolver::error_correction ( Array & u,
                                  Array & e_2h,
                                  int current_level
                                )
{
	Array & e_h = * tmp_grids_[current_level-1];
	e_h.fill(0.0);

	int width  = e_2h.getSize(DIM_1D);
    int height = e_2h.getSize(DIM_2D);

	// weighting for the restriction
	// stencil:
	//
	// w1 w2 w3
	// w4 w5 w6
	// w7 w8 w9

	real w1 = 0.25;
	real w2 = 0.5;
	real w3 = 0.25;
	real w4 = 0.5;
	real w5 = 1.0;
	real w6 = 0.5;
	real w7 = 0.25;
	real w8 = 0.5;
	real w9 = 0.25;
	

	// calculate I * e_2h (error to finer grid)
	for (int j = 0; j < height; j++) {
	for (int i = 0; i < width; i++)
	{   
		int mid_i = 2 * i + 1;	
		int mid_j = 2 * j + 1;	

		e_h(mid_i - 1, mid_j + 1) += w1 * e_2h(i, j);
		e_h(mid_i    , mid_j + 1) += w2 * e_2h(i, j);
		e_h(mid_i + 1, mid_j + 1) += w3 * e_2h(i, j);
		e_h(mid_i - 1, mid_j    ) += w4 * e_2h(i, j);
		e_h(mid_i    , mid_j    ) += w5 * e_2h(i, j);
		e_h(mid_i + 1, mid_j    ) += w6 * e_2h(i, j);
		e_h(mid_i - 1, mid_j - 1) += w7 * e_2h(i, j);
		e_h(mid_i    , mid_j - 1) += w8 * e_2h(i, j);
		e_h(mid_i + 1, mid_j - 1) += w9 * e_2h(i, j);

	}
	}

	// add to the solution

	width  = u.getSize(DIM_1D);
    height = u.getSize(DIM_2D);

	for (int j = 0; j < height; j++) {
	for (int i = 0; i < width; i++)
	{   
		u(i, j) += e_h(i, j);
	}
	}

}

void MGSolver::compose_right_hand_side ( Array & u, 
                                        Array & f,
                                        Array & r_2h,
                                        int current_level,
                                        real h
                                      )
{
	// store f - Au in temporary grid storage
	Array & res = * tmp_grids_[current_level-1];

	real h_2_inv = 1.0 / (h*h);

	int width  = u.getSize(DIM_1D);
    int height = u.getSize(DIM_2D);

	// calculate f - Au
	for (int j = 0; j < height; j++) {
	for (int i = 0; i < width; i++)
	{   
		// bottom left corner
		if (i == 0 && j == 0)
		{   
			res(i, j) = f(i, j) - h_2_inv * ( 4 * u(i, j) - u(i, j+1)                         - u(i+1, j) );
		}   
		// top left corner
		else if (i == 0 && j == height-1)
		{   
			res(i, j) = f(i, j) - h_2_inv * ( 4 * u(i, j)             - u(i, j-1)             - u(i+1, j) );
		}   
		// bottom right corner
		else if (i == width-1 && j == 0)
		{   
			res(i, j) = f(i, j) - h_2_inv * ( 4 * u(i, j) - u(i, j+1)             - u(i-1, j)             );
		}   
		// top right corner
		else if (i == width-1 && j == height-1)
		{   
			res(i, j) = f(i, j) - h_2_inv * ( 4 * u(i, j)             - u(i, j-1) - u(i-1, j)             );
		}   
		// left side
		else if (i == 0)
		{   
			res(i, j) = f(i, j) - h_2_inv * ( 4 * u(i, j) - u(i, j+1) - u(i, j-1)             - u(i+1, j) );
		}   
		// right side
		else if (i == width-1)
		{   
			res(i, j) = f(i, j) - h_2_inv * ( 4 * u(i, j) - u(i, j+1) - u(i, j-1) - u(i-1, j)             );
		}   
		// bottom
		else if (j == 0)
		{   
			res(i, j) = f(i, j) - h_2_inv * ( 4 * u(i, j) - u(i, j+1)             - u(i-1, j) - u(i+1, j) );
		}   
		// top
		else if (j == height-1)
		{   
			res(i, j) = f(i, j) - h_2_inv * ( 4 * u(i, j)             - u(i, j-1) - u(i-1, j) - u(i+1, j) );
		}   
		// inner domain
		else
		{   
			res(i, j) = f(i, j) - h_2_inv * ( 4 * u(i, j) - u(i, j+1) - u(i, j-1) - u(i-1, j) - u(i+1, j) );
		}   
	}   
	}


	// calculate new height
	width  = r_2h.getSize(DIM_1D);
    height = r_2h.getSize(DIM_2D);

	// weighting for the restriction
	// stencil:
	//
	// w1 w2 w3
	// w4 w5 w6
	// w7 w8 w9

	real w1 = 0.0625;
	real w2 = 0.125;
	real w3 = 0.0625;
	real w4 = 0.125;
	real w5 = 0.25;
	real w6 = 0.125;
	real w7 = 0.0625;
	real w8 = 0.125;
	real w9 = 0.0625;

	// restrict to coarser domain
	for (int j = 0; j < height; j++) {
	for (int i = 0; i < width; i++)
	{   
		int mid_i = 2*i + 1;
		int mid_j = 2*j + 1;
		
		r_2h(i, j) = w1 * u(mid_i - 1, mid_j + 1) +
		             w2 * u(mid_i    , mid_j + 1) +
		             w3 * u(mid_i + 1, mid_j + 1) +
		             w4 * u(mid_i - 1, mid_j    ) +
		             w5 * u(mid_i    , mid_j    ) +
		             w6 * u(mid_i + 1, mid_j    ) +
		             w7 * u(mid_i - 1, mid_j - 1) +
		             w8 * u(mid_i    , mid_j - 1) +
		             w9 * u(mid_i + 1, mid_j - 1);
	}   
	}
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
			sum += pow(f(i, j) - h_2_inv * ( 4.0 * u(i, j) - u(i, j+1)                         - u(i+1, j) ), 2.0);
		}   
		// top left corner
		else if (i == 0 && j == height-1)
		{   
			sum += pow(f(i, j) - h_2_inv * ( 4.0 * u(i, j)             - u(i, j-1)             - u(i+1, j) ), 2.0);
		}   
		// bottom right corner
		else if (i == width-1 && j == 0)
		{   
			sum += pow(f(i, j) - h_2_inv * ( 4.0 * u(i, j) - u(i, j+1)             - u(i-1, j)             ), 2.0);
		}   
		// top right corner
		else if (i == width-1 && j == height-1)
		{   
			sum += pow(f(i, j) - h_2_inv * ( 4.0 * u(i, j)             - u(i, j-1) - u(i-1, j)             ), 2.0);
		}   
		// left side
		else if (i == 0)
		{   
			sum += pow(f(i, j) - h_2_inv * ( 4.0 * u(i, j) - u(i, j+1) - u(i, j-1)             - u(i+1, j) ), 2.0);
		}   
		// right side
		else if (i == width-1)
		{   
			sum += pow(f(i, j) - h_2_inv * ( 4.0 * u(i, j) - u(i, j+1) - u(i, j-1) - u(i-1, j)             ), 2.0);
		}   
		// bottom
		else if (j == 0)
		{   
			sum += pow(f(i, j) - h_2_inv * ( 4.0 * u(i, j) - u(i, j+1)             - u(i-1, j) - u(i+1, j) ), 2.0);
		}   
		// top
		else if (j == height-1)
		{   
			sum += pow(f(i, j) - h_2_inv * ( 4.0 * u(i, j)             - u(i, j-1) - u(i-1, j) - u(i+1, j) ), 2.0);
		}   
		// inner domain
		else
		{   
			sum += pow(f(i, j) - h_2_inv * ( 4.0 * u(i, j) - u(i, j+1) - u(i, j-1) - u(i-1, j) - u(i+1, j) ), 2.0);
		}   
	}   
	}


	return sqrt(sum / (real) (width * height));
}
