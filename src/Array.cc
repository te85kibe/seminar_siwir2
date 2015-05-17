#include "Array.hh"
#include <cstdlib>
#include <iostream>
#include "Debug.hh"
#include <cmath>

//===================================================================================================================
//
//  Constructors
//
//===================================================================================================================


Array::Array( int xSize )
	: xSize_(xSize)
	, ySize_(1)
	, zSize_(1)
	, dimension_(DIM_1D)
{
   // TODO construct 1D array here
	CHECK_MSG(xSize > 0, "xSize must be positive");
	size_ = xSize;
	ar = new real[size_];
}

Array::Array( int xSize, int ySize )
	: xSize_(xSize)
	, ySize_(ySize)
	, zSize_(1)
	, dimension_(DIM_2D)
{
   // TODO construct 2D array here
	CHECK_MSG(xSize > 0 && ySize > 0, "xSize and ySize must be positive");
	size_ = xSize * ySize;
	ar = new real[size_];
}

Array::Array( int xSize, int ySize, int zSize )
	: xSize_(xSize)
	, ySize_(ySize)
	, zSize_(zSize)
	, dimension_(DIM_3D)
{
   // TODO construct 3D array here
	CHECK_MSG(xSize > 0 && ySize > 0 && zSize > 0, "xSize, ySize and zSize must be positive");
	size_ = xSize * ySize * zSize;
	ar = new real[size_];
}

Array::Array(const Array & s)
{
	dimension_ = s.getDimension();
	switch (dimension_)
	{
	case DIM_1D:
		xSize_ = s.getSize(DIM_1D);
		ySize_ = 1;
		zSize_ = 1;
		break;
	case DIM_2D:
		xSize_ = s.getSize(DIM_1D);
		ySize_ = s.getSize(DIM_2D);
		zSize_ = 1;
		break;
	case DIM_3D:
		xSize_ = s.getSize(DIM_1D);
		ySize_ = s.getSize(DIM_2D);
		zSize_ = s.getSize(DIM_3D);
		break;
	}
	size_ = s.getSize();
	ar = new real[size_];

	real *p_s = s.getArray();
	for (int i = 0; i < size_; i++)
	{
		ar[i] = p_s[i];
	}

}

Array& Array::operator=(const Array& s)
{
	real *local_ar = new real[s.getSize()];
	real *s_ar = s.getArray();

	for (int i = 0; i < s.getSize(); i++)
	{
		local_ar[i] = s_ar[i];
	}

	this->ar = local_ar;
	this->xSize_ = s.getSize(DIM_1D);
	this->ySize_ = s.getSize(DIM_2D);
	this->zSize_ = s.getSize(DIM_3D);
	this->dimension_ = s.getDimension();

	return *this;
}

//	Destructor
Array::~Array()
{
	if (NULL != ar)
	{
		delete[] ar;
	}
}




//===================================================================================================================
//
//  Convenience Functions
//
//===================================================================================================================


//initialize the whole array with a constant value
void Array::fill( real value )
{
   // TODO
   // you might want to use std::fill() here
	for (int i = 0; i < getSize(); i++)
	{
		ar[i] = value;
	}

}

void Array::fillRandom()
{
	for (int i = 0; i < getSize(); i++)
		{
			ar[i] = (real)rand() / RAND_MAX;;
		}
}

// returns dimension
int Array::getDimension() const
{
	return dimension_;
}

real * Array::getArray() const
{
	return ar;
}

// Print the whole array (for debugging purposes)
void Array::print()
{
   // TODO
   // For 2D Arrays the positive x-coordinate goes to the right
   //                   positive y-coordinate goes upwards
   //      -> the line with highest y-value should be printed first
	switch (dimension_)
	{
	case DIM_1D:
		for (int i = 0; i < xSize_; i++)
		{
			std::cout << Array::operator ()(i) << " ";
		}
		std::cout << std::endl;
		break;
	case DIM_2D:
		for (int i = ySize_ - 1; i >= 0; i--)
		{
			for (int j = 0; j < xSize_; j++)
			{
				std::cout << Array::operator ()(j,i) << " ";
			}
			std::cout << std::endl;
		}
		break;
	case DIM_3D:
		for (int k = 0; k < zSize_; k++)
		{
			std::cout << "z = " << k << std::endl;
			for (int i = ySize_ - 1; i >= 0; i--)
				{
					for (int j = 0; j < xSize_; j++)
					{
						std::cout << Array::operator ()(j,i,k) << " ";
					}
					std::cout << std::endl;
				}
			std::cout << std::endl;
		}
		break;

	}

}

int Array::getSize( int dimension ) const
{
   //TODO
	int ret = 0;
	switch (dimension)
	{
	case DIM_1D:
		ret = xSize_;
		break;
	case DIM_2D:
		ret = ySize_;
		break;
	case DIM_3D:
		ret = zSize_;
		break;
	}
	return ret;
}

//return total size of the array
int Array::getSize() const
{
   //TODO
   return size_;
}

real Array::getAbsMax()
{
	real max = fabs(ar[0]);

	for (int i = 1; i < size_; i++)
	{
		if (fabs(ar[i]) > max)
		{
			max = fabs(ar[i]);
		}
	}

	return max;
}

// Operator() 1D
real& Array::operator ()(int i)
{
   //TODO
   //static double dummy;
    CHECK_MSG(i >= 0 && i < size_, "Index i out of bounds");
    CHECK_MSG(dimension_ == DIM_1D, "Wrong dimension.");
   return ar[i];
}

// Operator() 2D
real& Array::operator ()(int i,int j)
{
   //TODO
   //static double dummy;
    CHECK_MSG(i >= 0 && i < xSize_, "Index i out of bounds");
    CHECK_MSG(j >= 0 && j < ySize_, "Index j out of bounds");
    CHECK_MSG(dimension_ == DIM_2D, "Wrong dimension.");
   return ar[xSize_*j + i]; 
}

// Operator() 3D
real& Array::operator ()(int i, int j, int k)
{
   //TODO
   //static double dummy;
    CHECK_MSG(i >= 0 && i < xSize_, "Index i out of bounds");
    CHECK_MSG(j >= 0 && j < ySize_, "Index j out of bounds");
    CHECK_MSG(k >= 0 && k < zSize_, "Index k out of bounds");
    CHECK_MSG(dimension_ == DIM_3D, "Wrong dimension.");
   return ar[(xSize_ * ySize_)*k + xSize_*j + i]; 
}

///////////
// CONST //
///////////

// Operator() 1D
const real& Array::operator ()(int i) const
{
   //TODO
   //static double dummy;
    CHECK_MSG(i >= 0 && i < size_, "Index i out of bounds");
    CHECK_MSG(dimension_ == DIM_1D, "Wrong dimension.");
   return ar[i];
}

// Operator() 2D
const real& Array::operator ()(int i,int j) const
{
   //TODO
   //static double dummy;
    CHECK_MSG(i >= 0 && i < xSize_, "Index i out of bounds");
    CHECK_MSG(j >= 0 && j < ySize_, "Index j out of bounds");
    CHECK_MSG(dimension_ == DIM_2D, "Wrong dimension.");
   return ar[xSize_*j + i];
}

// Operator() 3D
const real& Array::operator ()(int i, int j, int k) const
{
   //TODO
   //static double dummy;
    CHECK_MSG(i >= 0 && i < xSize_, "Index i out of bounds");
    CHECK_MSG(j >= 0 && j < ySize_, "Index j out of bounds");
    CHECK_MSG(k >= 0 && k < zSize_, "Index k out of bounds");
    CHECK_MSG(dimension_ == DIM_3D, "Wrong dimension.");
   return ar[(xSize_ * ySize_)*k + xSize_*j + i];
}


#if 0
void Array::normalize()
{
	// calc avg
	real avg;
	real sum = 0.0;

	for (int i = 0; i < size_; i++)
	{
		sum += ar[i];
	}

	avg = sum / size_;

	for (int i = 0; i < size_; i++)
	{
		ar[i] -= avg;
	}

}
#endif
