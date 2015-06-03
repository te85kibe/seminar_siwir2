#include <cmath>
#include <iostream>
#include <cmath>
#include <math.h>
#include <fstream>
#include <omp.h>

#include "MGSolver.hh"
#include "Array.hh"
#include "Smoother.hh"

#ifndef PI
#define PI (M_PI)
#endif

#define PRINT_ERROR (1)
#define PRINT_RESIDUAL (1)

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

	// initialize solution
	solution_ = new Array(std::pow(2, levels) + 1, std::pow(2, levels) + 1);

	// allocate memory for v and r
	for (int i = 1; i <= levels; i++)
	{

		v_grids_[i-1] = new Array ( std::pow(2, i) + 1, std::pow(2, i) + 1);
		r_grids_[i-1] = new Array ( std::pow(2, i) + 1, std::pow(2, i) + 1);
		tmp_grids_[i-1] = new Array ( std::pow(2, i) + 1, std::pow(2, i) + 1);
		h_intervals_[i-1] = 2.0 / std::pow(2, i);

		
#if 0
		std::cout << h_intervals_[i-1] << std::endl;
		v_grids_[i-1]->print();
#endif
		
	}	

}

MGSolver::~MGSolver()
{
	delete solution_;
}

void MGSolver::initialize_assignment_01 ()
{

	// initialize bc on the unknown array boundary
	// bc is sin(x*pi) sinh(y*pi) which is only != 0 at the top boundary (y = 1)

	Array * finest_grid = v_grids_.back();
	real h              = h_intervals_.back();

	for (int col = 0;
         col < finest_grid->getSize(DIM_1D); 
         col++)
	{
		finest_grid->operator()(col, finest_grid->getSize(DIM_2D)-1) = sin(PI * (real) col * h) * sinh(PI);
	}



	// initialize solution
	for (int row = 0; row < solution_->getSize(DIM_2D); row++)
	{
		for (int col = 0; col < solution_->getSize(DIM_1D); col++)
		{
			solution_->operator()(col, row) = sin(PI * (real) col * h) * sinh(PI * (real) row * h);	
		}
	}
}

void MGSolver::initialize_seminar ()
{


	Array * finest_grid = v_grids_.back();
	real h              = h_intervals_.back();

	// initialize bc on top boundary

	for (int col = 0;
         col < finest_grid->getSize(DIM_1D); 
         col++)
	{
		real x = -1.0+h*col;
		real radius = sqrt(1.0+x*x);
		real phi = atan2(1.0,x);
		if(phi<0)
		{
			phi=phi+2*PI;
		}
		finest_grid->operator()(col, finest_grid->getSize(DIM_2D)-1) = std::pow(radius,0.5)*sin(0.5*phi);

	}
	// initialize bc on right boundary
	for (int row = 0;
         row < finest_grid->getSize(DIM_2D); 
         row++)
	{
		real y = -1.0+h*row;
		real radius = sqrt(1.0+y*y);
		real phi = atan2(y,1.0);
		if(phi<0)
		{
			phi=phi+2*PI;
		}

		finest_grid->operator()(finest_grid->getSize(DIM_1D)-1,row)  = std::pow(radius,0.5)*sin(0.5*phi);
	}
	
	// initialize bc on left boundary
	for (int row = 0;
         row < finest_grid->getSize(DIM_2D); 
         row++)
	{
		real y = -1.0+h*row;
		real radius = sqrt(1.0+y*y);
		real phi = atan2(y,-1.0);
		if(phi<0)
		{
			phi=phi+2*PI;
		}

		finest_grid->operator()(0,row)  = std::pow(radius,0.5)*sin(0.5*phi);
	}
	
	// initialize bc on bottom boundary
	for (int col = 0;
         col < finest_grid->getSize(DIM_1D); 
         col++)
	{
		real x = -1.0+h*col;
		real radius = sqrt(1.0+x*x);
		real phi = atan2(-1.0,x);
		if(phi<0)
		{
			phi=phi+2*PI;
		}
		finest_grid->operator()(col,0) = std::pow(radius,0.5)*sin(0.5*phi);

	}	

	// initialize solution
	for (int row = 0; row < solution_->getSize(DIM_2D); row++)
	{
		for (int col = 0; col < solution_->getSize(DIM_1D); col++)
		{
			real y = -1.0+(real)h*row;
			real x = -1.0+(real)h*col;
			real radius = sqrt(x*x+y*y);
			real phi = atan2(y,x);
			if(phi<0)
			{
				phi=phi+2*PI;
			}
			solution_->operator()(col, row) = sqrt(radius)*sin(0.5*phi);	
		}
	}

		
	// initialize rhs
	r_grids_.back()->fill(0.0);

}

void MGSolver::initialize_assignment_01_BONUS ()
{
	// init Dirichlet BC
	Array * finest_grid = v_grids_.back();
	real h				= h_intervals_.back();

	for (int col = 0;
         col < finest_grid->getSize(DIM_1D); 
         col++)
	{
		finest_grid->operator()(col, 0) 								= h*col * (1.0 - h*col);
		finest_grid->operator()(col, finest_grid->getSize(DIM_2D)-1) 	= h*col * (1.0 - h*col);               
	}

	// initialize solution
	// solution should be f(x, y) = - x^2 + x
	for (int row = 0; row < solution_->getSize(DIM_2D); row++)
	{
		for (int col = 0; col < solution_->getSize(DIM_1D); col++)
		{
			solution_->operator()(col, row) = - ((real) col * h) * ((real) col * h) + (real) col * h;
		}
	}
	


	// initialize rhs
	r_grids_.back()->fill(2.0);
	
}

void MGSolver::initialize_random ()
{
	Array * finest_grid = v_grids_.back();
	finest_grid->fillRandom();

	// TODO boundaries to zero for finest_grid

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

	// 


	if (level == 1)
	{
		// SOLVE
		//std::cout << v_grids_[level-1]->operator()(0, 0) << std::endl;
		v_grids_[level-1]->operator()(0, 0) = r_grids_[level-1]->operator()(0, 0) * h_intervals_[level-1] * h_intervals_[level-1] * 0.25;
	//std::cout << "V GRID" << std::endl;
	//v_grids_[level-1]->print();	
	//std::cout << "R GRID" << std::endl;
	//r_grids_[level-1]->print();	
	//std::cout << "T GRID" << std::endl;
	//tmp_grids_[level-1]->print();	
		
		return;
	}

	v_grids_[level-2]->fill(0.0);
	tmp_grids_[level-1]->fill(0.0);
		
	// 1. perform pre_smooth gauss seidel iterations on v
	//std::cout << "1. perform pre_smooth gauss seidel iterations on v" << std::endl;
	smoother_.smooth_red_black_gauss_seidel_2d ( * v_grids_[level-1],
                                                 * r_grids_[level-1],
                                                 pre_smooth,
                                                 h_intervals_[level-1],
												 level == (int)v_grids_.size());	// true if finest grid

#if 1

	// 2. calculate coarser right hand side
	//std::cout << "2. calculate coarser right hand side" << std::endl;
	compose_right_hand_side ( * v_grids_[level-1],
                              * r_grids_[level-1],
                              * r_grids_[level-2],
                              level,
                              h_intervals_[level-1]);

	// tmp_grids_[level-2]->print();

	// 3.recursive call
	v_cycle_pvt ( pre_smooth, post_smooth, level-1 );
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
                                                 h_intervals_[level-1],
												 level == (int)v_grids_.size());	// true if finest grid
#endif

#if 0
	std::cout << "V GRID" << std::endl;
	// v_grids_[level-1]->print();	
	std::cout << "R GRID" << std::endl;
	//r_grids_[level-1]->print();	
	std::cout << "T GRID" << std::endl;
	//tmp_grids_[level-1]->print();	
#endif
	

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


	for (int j = 1; j < height/2; j++) {
	for (int i = 1; i < width-1; i++)
	{   
		int mid_i = 2 * i;	
		int mid_j = 2 * j;	

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
	for (int i = 1; i < width/2; i++)
	{   
		int j=height/2;

		int mid_i = 2 * i;	
		int mid_j = 2 * j;	

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

	for (int j = height/2+1; j < height-1; j++) {
	for (int i = 1; i < width-1; i++)
	{   
		int mid_i = 2 * i;	
		int mid_j = 2 * j;	

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
//	}

	// add to the solution

	width  = u.getSize(DIM_1D);
   	height = u.getSize(DIM_2D);
	#pragma omp parallel for collapse(2) schedule(static)
	for (int j = 1; j < height/2; j++) {
	for (int i = 1; i < width-1; i++)
	{   
		u(i, j) += e_h(i, j);
	}
	}

	for (int i = 1; i < width/2; i++)
	{	
		int j=height/2;
		u(i, j) += e_h(i, j);
	}
	#pragma omp parallel for collapse(2) schedule(static)
	for (int j = height/2+1; j < height-1; j++) {
	for (int i = 1; i < width-1; i++)
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
	#pragma omp parallel for collapse(2) schedule(static)
	for (int j = 1; j < height/2; j++) {
	for (int i = 1; i < width-1; i++)
	{   
		res(i, j) = f(i, j) - h_2_inv * ( 4.0 * u(i, j) - u(i, j+1) - u(i, j-1) - u(i-1, j) - u(i+1, j) );
	}   
	}
	
	for (int i = 1; i < width/2; i++)
	{	
		int j=height/2;

		res(i, j) = f(i, j) - h_2_inv * ( 4.0 * u(i, j) - u(i, j+1) - u(i, j-1) - u(i-1, j) - u(i+1, j) );
	}   
	
	#pragma omp parallel for collapse(2) schedule(static)
	for (int j = height/2+1; j < height-1; j++) {
	for (int i = 1; i < width-1; i++)
	{   
		res(i, j) = f(i, j) - h_2_inv * ( 4.0 * u(i, j) - u(i, j+1) - u(i, j-1) - u(i-1, j) - u(i+1, j) );
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
	#pragma omp parallel for collapse(2) schedule(static)
	for (int j = 1; j < height/2; j++) {
	for (int i = 1; i < width-1; i++)
	{   
		int mid_i = 2*i;
		int mid_j = 2*j;
		
		r_2h(i, j) = w1 * res(mid_i - 1, mid_j + 1) +
		             w2 * res(mid_i    , mid_j + 1) +
		             w3 * res(mid_i + 1, mid_j + 1) +
		             w4 * res(mid_i - 1, mid_j    ) +
		             w5 * res(mid_i    , mid_j    ) +
		             w6 * res(mid_i + 1, mid_j    ) +
		             w7 * res(mid_i - 1, mid_j - 1) +
		             w8 * res(mid_i    , mid_j - 1) +
		             w9 * res(mid_i + 1, mid_j - 1);
	}   
	}
	for (int i = 1; i < width/2; i++)
	{	
		int j=height/2;
	  
		int mid_i = 2*i;
		int mid_j = 2*j;
		
		r_2h(i, j) = w1 * res(mid_i - 1, mid_j + 1) +
		             w2 * res(mid_i    , mid_j + 1) +
		             w3 * res(mid_i + 1, mid_j + 1) +
		             w4 * res(mid_i - 1, mid_j    ) +
		             w5 * res(mid_i    , mid_j    ) +
		             w6 * res(mid_i + 1, mid_j    ) +
		             w7 * res(mid_i - 1, mid_j - 1) +
		             w8 * res(mid_i    , mid_j - 1) +
		             w9 * res(mid_i + 1, mid_j - 1);
	  
	}
	#pragma omp parallel for collapse(2) schedule(static)
	for (int j = height/2+1; j < height-1; j++) {
		for (int i = 1; i < width-1; i++)
		{   
			int mid_i = 2*i;
			int mid_j = 2*j;

			r_2h(i, j) = w1 * res(mid_i - 1, mid_j + 1) +
				w2 * res(mid_i    , mid_j + 1) +
				w3 * res(mid_i + 1, mid_j + 1) +
				w4 * res(mid_i - 1, mid_j    ) +
				w5 * res(mid_i    , mid_j    ) +
				w6 * res(mid_i + 1, mid_j    ) +
				w7 * res(mid_i - 1, mid_j - 1) +
				w8 * res(mid_i    , mid_j - 1) +
				w9 * res(mid_i + 1, mid_j - 1);
		}   
	}


}

void MGSolver::restrict_2d ( Array & u,
                             Array & u_2h)
{

	int width  = u_2h.getSize(DIM_1D);
    int height = u_2h.getSize(DIM_2D);

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
	#pragma omp parallel for collapse(2) schedule(static)
	for (int j = 1; j < height/2; j++) {
	for (int i = 1; i < width-1; i++)
	{   
		int mid_i = 2*i;
		int mid_j = 2*j;
		
		u_2h(i, j) = w1 * u(mid_i - 1, mid_j + 1) +
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
	for (int i = 1; i < width/2; i++)
	{	
		int j=height/2;
		int mid_i = 2*i;
		int mid_j = 2*j;
		
		u_2h(i, j) = w1 * u(mid_i - 1, mid_j + 1) +
		             w2 * u(mid_i    , mid_j + 1) +
		             w3 * u(mid_i + 1, mid_j + 1) +
		             w4 * u(mid_i - 1, mid_j    ) +
		             w5 * u(mid_i    , mid_j    ) +
		             w6 * u(mid_i + 1, mid_j    ) +
		             w7 * u(mid_i - 1, mid_j - 1) +
		             w8 * u(mid_i    , mid_j - 1) +
		             w9 * u(mid_i + 1, mid_j - 1);
	}
	#pragma omp parallel for collapse(2) schedule(static)
	for (int j = height/2+1; j < height-1; j++) {
	for (int i = 1; i < width-1; i++)
	{   
		int mid_i = 2*i;
		int mid_j = 2*j;
		
		u_2h(i, j) = w1 * u(mid_i - 1, mid_j + 1) +
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
	for (int j = 1; j < height/2; j++) {
	for (int i = 1; i < width-1; i++)
	{   
		sum += pow(f(i, j) - h_2_inv * ( 4.0 * u(i, j) - u(i, j+1) - u(i, j-1) - u(i-1, j) - u(i+1, j) ), 2.0);
	}   
	}
	for (int i = 1; i < width/2; i++)
	{	
		int j=height/2;
		sum += pow(f(i, j) - h_2_inv * ( 4.0 * u(i, j) - u(i, j+1) - u(i, j-1) - u(i-1, j) - u(i+1, j) ), 2.0);
	}   
	
	for (int j = height/2+1; j < height-1; j++) {
	for (int i = 1; i < width-1; i++)
	{   
		sum += pow(f(i, j) - h_2_inv * ( 4.0 * u(i, j) - u(i, j+1) - u(i, j-1) - u(i-1, j) - u(i+1, j) ), 2.0);
	}   
	}




	return sqrt(sum / (real) ((width-2) * (height-2)));
}

real MGSolver::error_L2 ( Array & approximation,
							Array & solution,
							real h
						  )
{
	real sum = 0.0;	

	int width  = approximation.getSize(DIM_1D);
    int height = approximation.getSize(DIM_2D);

	// add up squares of the entries of the residual
	for (int j = 1; j < height/2; j++) {
	for (int i = 1; i < width-1; i++)
	{   
		sum += pow(approximation(i,j) - solution(i,j), 2.0);
	}   
	}
	for (int i = 1; i < width/2; i++)
	{	
		int j = height/2;
		sum += pow(approximation(i,j) - solution(i,j), 2.0);  
	}
	for (int j = height/2+1; j < height-1; j++) {
	for (int i = 1; i < width-1; i++)
	{   
		sum += pow(approximation(i,j) - solution(i,j), 2.0);
	}   
	}



	return sqrt(sum / (real) ((width-2) * (height-2)));

	
}

int MGSolver::saveToFile(std::string filename) const
{
	Array *u = v_grids_.back();
//	std::cout << "width: " << u->getSize(DIM_1D) << std::endl;

	std::ofstream gnuFile(filename);
	if (gnuFile.is_open())
	{
		for (int j = u->getSize(DIM_2D) - 1; j >= 0; j--)
		{
			for (int i = 0; i < u->getSize(DIM_1D); i++)
			{
				gnuFile << (double) -1.0+2.0*i/(u->getSize(DIM_1D)-1) << " " << (double) 1.0-2.0*j/(u->getSize(DIM_2D)-1) << " " << u->operator()(i,j) << "\n";
			}
		//	gnuFile << "\n";
		}
		gnuFile.close();
		return 0;
	}
	else
	{
		return 1;
	}
}
