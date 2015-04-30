
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
		for (int j = 0; j < height; j++)
		{
			int i;
			if (j % 2 == 0)
			{
				i = 0;
			}
			else
			{
				i = 1;
			}
				
			for (; i < width; i+=2)
			{
				// bottom left corner
				if (i == 0 && j == 0)
				{
					u(i,j) = factor * (f(i,j) + h_2_inv * ( u(i+1, j) + u(i, j+1)));
				}
				// top left corner
				else if (i == 0 && j == height-1)
				{
					u(i,j) = factor * (f(i,j) + h_2_inv * ( u(i+1, j) + u(i, j-1)));
				}
				// bottom right corner
				else if (i == width-1 && j == 0)
				{
					u(i,j) = factor * (f(i,j) + h_2_inv * ( u(i-1, j) + u(i, j+1)));
				}
				// top right corner
				else if (i == width-1 && j == height-1)
				{
					u(i,j) = factor * (f(i,j) + h_2_inv * ( u(i-1, j) + u(i, j-1)));
				}
				// left side
				else if (i == 0)
				{
					u(i,j) = factor * (f(i,j) + h_2_inv * ( u(i+1, j) + u(i, j+1) + u(i, j-1)));
				}
				// right side
				else if (i == width-1)
				{
					u(i,j) = factor * (f(i,j) + h_2_inv * ( u(i-1, j) + u(i, j+1) + u(i, j-1)));
				}
				// bottom
				else if (j == 0)
				{
					u(i,j) = factor * (f(i,j) + h_2_inv * ( u(i-1, j) + u(i+1, j) + u(i, j+1)));
				}
				// top
				else if (j == height-1)
				{
					u(i,j) = factor * (f(i,j) + h_2_inv * ( u(i-1, j) + u(i+1, j) + u(i, j-1)));
				}
				// inner domain
				else
				{
					u(i,j) = factor * (f(i,j) + h_2_inv * ( u(i-1, j) + u(i+1, j) + u(i, j+1) + u(i, j-1)));
				}
			}
		}
	
		// black points
		for (int j = 0; j < height; j++)
		{
			int i;
			if (j % 2 == 0)
			{
				i = 1;
			}
			else
			{
				i = 0;
			}
				
			for (; i < width; i+=2)
			{
				// bottom left corner
				if (i == 0 && j == 0)
				{
					u(i,j) = factor * (f(i,j) + h_2_inv * ( u(i+1, j) + u(i, j+1)));
				}
				// top left corner
				else if (i == 0 && j == height-1)
				{
					u(i,j) = factor * (f(i,j) + h_2_inv * ( u(i+1, j) + u(i, j-1)));
				}
				// bottom right corner
				else if (i == width-1 && j == 0)
				{
					u(i,j) = factor * (f(i,j) + h_2_inv * ( u(i-1, j) + u(i, j+1)));
				}
				// top right corner
				else if (i == width-1 && j == height-1)
				{
					u(i,j) = factor * (f(i,j) + h_2_inv * ( u(i-1, j) + u(i, j-1)));
				}
				// left side
				else if (i == 0)
				{
					u(i,j) = factor * (f(i,j) + h_2_inv * ( u(i+1, j) + u(i, j+1) + u(i, j-1)));
				}
				// right side
				else if (i == width-1)
				{
					u(i,j) = factor * (f(i,j) + h_2_inv * ( u(i-1, j) + u(i, j+1) + u(i, j-1)));
				}
				// bottom
				else if (j == 0)
				{
					u(i,j) = factor * (f(i,j) + h_2_inv * ( u(i-1, j) + u(i+1, j) + u(i, j+1)));
				}
				// top
				else if (j == height-1)
				{
					u(i,j) = factor * (f(i,j) + h_2_inv * ( u(i-1, j) + u(i+1, j) + u(i, j-1)));
				}
				// inner domain
				else
				{
					u(i,j) = factor * (f(i,j) + h_2_inv * ( u(i-1, j) + u(i+1, j) + u(i, j+1) + u(i, j-1)));
				}
			}
		}
	}
}

