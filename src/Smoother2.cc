
#include "Smoother.hh"
#include "Types.hh"
#include "Array.hh"

void Smoother::smooth_red_black_gauss_seidel_2d ( Array & u,    // modify this array
                                                  Array & f,    // rhs
                                                  int times,     // number of sweeps
                                                  real h        // spacing         
                                                )
{

	int width  = u.getSize(DIM_1D);
	int height = u.getSize(DIM_2D);

	real h_2     = h * h;
	real h_2_inv = 1.0 / h_2;
	real factor  = h_2 * 0.25;

	for (int iter = 0; iter < times; iter++)
	{
		// red points
		for (int j = 1; j < height-1; j++)
		{
			for (int i = 1; i < width-1; i++)
			{
					// inner domain
					// i+j gerade
					if( ((i + j) % 2) == 0)
						u(i,j) = factor * (f(i,j) + h_2_inv * ( u(i-1, j) + u(i+1, j) + u(i, j+1) + u(i, j-1)));
			}
		}

		// black points
		for (int j = 1; j < height-1; j++)
		{
			for (int i = 1; i < width-1; i++)
			{
					// inner domain
					// i+j ungerade
					if( ((i + j) % 2) == 1)
						u(i,j) = factor * (f(i,j) + h_2_inv * ( u(i-1, j) + u(i+1, j) + u(i, j+1) + u(i, j-1)));
			}
		}
	}
}
