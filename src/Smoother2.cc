
#if 1
#include "Smoother.hh"
#include "Types.hh"
#include "Array.hh"
#include "MGSolver.hh"
#include <omp.h>

void Smoother::smooth_red_black_gauss_seidel_2d ( Array & u,    // modify this array
                                                  Array & f,    // rhs
                                                  int times,     // number of sweeps
                                                  real h,       // spacing         
												  bool finest_grid
                                                )
{

#ifdef NEUMANN
	if (finest_grid)
	{
		update_neumann ( u, h );
	}
#endif

	(void) finest_grid;

	int width  = u.getSize(DIM_1D);
	int height = u.getSize(DIM_2D);

	real h_2     = h * h;
	real h_2_inv = 1.0 / h_2;
	real factor  = h_2 * 0.25;
	real omega = 1.4;

	#pragma omp parallel
	{ 
	for (int iter = 0; iter < times; iter++)
	{
		
		// red points
		#pragma omp for
		for (int j = 1; j < (height/2); j++)
		{
			for (int i = 1; i < width-1; i++)
			{
					// inner domain
					// i+j gerade
					if( ((i + j) % 2) == 0)
					//	u(i,j) = factor * (f(i,j) + h_2_inv * ( u(i-1, j) + u(i+1, j) + u(i, j+1) + u(i, j-1)));
						u(i,j) = omega * factor * (f(i,j) + h_2_inv * ( u(i-1, j) + u(i+1, j) + u(i, j+1) + u(i, j-1)))+(1-omega)*u(i,j);
					//	u(i+1,j+1) = omega * factor * (f(i+1,j+1) + h_2_inv * ( u(i, j+1) + u(i+2, j+1) + u(i+1, j+2) + u(i+1, j)))+(1-omega)*u(i+1,j+1);
			}
			//u(width-2,j) = omega * factor * (f(width-2,j) + h_2_inv * ( u(width-3, j) + u(width-1, j) + u(width-2, j+1) + u(width-2, j-1)))+(1-omega)*u(width-2,j);

		}
		/*for (int i =1; i < width-1;i+=2)
		{
			u(i,height/2-1) = omega * factor * (f(i,height/2-1) + h_2_inv * ( u(i-1, height/2-1) + u(i+1, height/2-1) + u(i, height/2) + u(i, height/2-2)))+(1-omega)*u(i,height/2-1);		
		}
		*/
		#pragma omp for
		for (int i = 1; i < (width/2); i++)
		{
			// inner domain
			// i+j gerade
			int j = height/2;
			if( ((i + j) % 2) == 0)
				//u(i,j) = factor * (f(i,j) + h_2_inv * ( u(i-1, j) + u(i+1, j) + u(i, j+1) + u(i, j-1)));
				u(i,j) = omega * factor * (f(i,j) + h_2_inv * ( u(i-1, j) + u(i+1, j) + u(i, j+1) + u(i, j-1)))+(1-omega)*u(i,j);
		}
		#pragma omp for
		for (int j = (height/2)+1; j < height-1; j++)
		{
			for (int i = 1; i < width-1; i++)
			{
					// inner domain
					// i+j gerade
					if( ((i + j) % 2) == 0)
						//u(i,j) = factor * (f(i,j) + h_2_inv * ( u(i-1, j) + u(i+1, j) + u(i, j+1) + u(i, j-1)));
						u(i,j) = omega * factor * (f(i,j) + h_2_inv * ( u(i-1, j) + u(i+1, j) + u(i, j+1) + u(i, j-1)))+(1-omega)*u(i,j);
			}
		}
		



		// black points
		#pragma omp for
		for (int j = 1; j < (height/2); j++)
		{
			for (int i = 1; i < width-1; i++)
			{
					// inner domain
					// i+j ungerade
					if( ((i + j) % 2) == 1)
						//u(i,j) = factor * (f(i,j) + h_2_inv * ( u(i-1, j) + u(i+1, j) + u(i, j+1) + u(i, j-1)));
						u(i,j) = omega * factor * (f(i,j) + h_2_inv * ( u(i-1, j) + u(i+1, j) + u(i, j+1) + u(i, j-1)))+(1-omega)*u(i,j);
			}
		}
		#pragma omp for
		for (int i = 1; i < (width/2); i++)
		{
			// inner domain
			// i+j ungerade
			int j = height/2;
			if( ((i + j) % 2) == 1)
				//u(i,j) = factor * (f(i,j) + h_2_inv * ( u(i-1, j) + u(i+1, j) + u(i, j+1) + u(i, j-1)));
				u(i,j) = omega * factor * (f(i,j) + h_2_inv * ( u(i-1, j) + u(i+1, j) + u(i, j+1) + u(i, j-1)))+(1-omega)*u(i,j);
		}
		#pragma omp for
		for (int j = (height/2)+1; j < height-1; j++)
		{
			for (int i = 1; i < width-1; i++)
			{
					// inner domain
					// i+j ungerade
					if( ((i + j) % 2) == 1)
						//u(i,j) = factor * (f(i,j) + h_2_inv * ( u(i-1, j) + u(i+1, j) + u(i, j+1) + u(i, j-1)));
						u(i,j) = omega * factor * (f(i,j) + h_2_inv * ( u(i-1, j) + u(i+1, j) + u(i, j+1) + u(i, j-1)))+(1-omega)*u(i,j);
			}
		}

	}
	}
}

void Smoother::update_neumann( Array & u,
                               real h ) 
{
    /*  
        To update neumann bc, set the boundary points to (u_inner - h)
        This comes from the first order difference: (u_boundary - u_inner) / h = -1
    */  

    for (int row = 0; row < u.getSize(DIM_2D); row++)
    {   
        u(0, row) = u(1, row) - h;
        u(u.getSize(DIM_1D)-1, row) = u(u.getSize(DIM_1D)-2, row) - h;
    }   
}


#endif
