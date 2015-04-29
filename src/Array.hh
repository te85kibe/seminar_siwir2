#ifndef ARRAY_HH
#define ARRAY_HH

#include <iostream>
#include "Types.hh"
#include "Debug.hh"

#define DIM_1D (0x0)
#define DIM_2D (0x1)
#define DIM_3D (0x2)

//*******************************************************************************************************************
/*!  Array class for 1,2 and 3 dimensions
*
*    - all elements should be stored in a contiguous chunk of memory ( no vector<vector> ! )
*/
//*******************************************************************************************************************
class Array
{
public:
   // Constructors for 1D,2D and 3D
   Array( int xSize );
   Array( int xSize, int ySize );
   Array( int xSize, int ySize, int zSize );


   // Depending on your implementation you might need the following:
   ~Array();
   Array(const Array& s);
   Array& operator= (const Array& s);


   // Access Operators for 1D, 2D and 3D
   inline real & operator () ( int i );
   inline real & operator () ( int i ,int j );
   inline real & operator () ( int i, int j, int k );

   // for const Arrays the following access operators are required
   inline const real & operator () ( int i ) const;
   inline const real & operator () ( int i ,int j ) const;
   inline const real & operator () ( int i, int j, int k ) const;



   // initialize the whole array with a constant value
   void fill( real value );

   // initialize with random reals between 0.0 and 1.0
   void fillRandom();

   // return total size of the array
   int getSize() const;

   // return xSize for dimension==0, ySize for dimension==1 and zSize for dimension==2
   // other dimension values are not allowed
   int getSize(int dimension ) const;

   // added to skeleton
   // return dimension of the array
   int getDimension() const;

   // added to skeleton
   // return address of first element of the array
   real * getArray() const;

   // Print the whole array ( for debugging purposes )
   void print();

   // returns the maximum (absolute) value of the array
   real getAbsMax();

   // normalizes the array around zero (subtracting by avg)
   void normalize();

private:

   // Sizes
   int xSize_;
   int ySize_;
   int zSize_;
   int size_;

   // Dimension of the array
   int dimension_;

   // Array pointer holding the actual elements
   real *ar;


};


//===================================================================================================================
//
//  Inline Access Operators and Sizes
//
//===================================================================================================================



// Operator() 1D
inline real& Array::operator ()(int i)
{
   //TODO
   //static double dummy;
	CHECK_MSG(i >= 0 && i < size_, "Index i out of bounds");
	CHECK_MSG(dimension_ == DIM_1D, "Wrong dimension.");
   return ar[i];
}

// Operator() 2D
inline real& Array::operator ()(int i,int j)
{
   //TODO
   //static double dummy;
	CHECK_MSG(i >= 0 && i < xSize_, "Index i out of bounds");
	CHECK_MSG(j >= 0 && j < ySize_, "Index j out of bounds");
	CHECK_MSG(dimension_ == DIM_2D, "Wrong dimension.");
   return ar[xSize_*j + i];
}

// Operator() 3D
inline real& Array::operator ()(int i, int j, int k)
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
inline const real& Array::operator ()(int i) const
{
   //TODO
   //static double dummy;
	CHECK_MSG(i >= 0 && i < size_, "Index i out of bounds");
	CHECK_MSG(dimension_ == DIM_1D, "Wrong dimension.");
   return ar[i];
}

// Operator() 2D
inline const real& Array::operator ()(int i,int j) const
{
   //TODO
   //static double dummy;
	CHECK_MSG(i >= 0 && i < xSize_, "Index i out of bounds");
	CHECK_MSG(j >= 0 && j < ySize_, "Index j out of bounds");
	CHECK_MSG(dimension_ == DIM_2D, "Wrong dimension.");
   return ar[xSize_*j + i];
}

// Operator() 3D
inline const real& Array::operator ()(int i, int j, int k) const
{
   //TODO
   //static double dummy;
	CHECK_MSG(i >= 0 && i < xSize_, "Index i out of bounds");
	CHECK_MSG(j >= 0 && j < ySize_, "Index j out of bounds");
	CHECK_MSG(k >= 0 && k < zSize_, "Index k out of bounds");
	CHECK_MSG(dimension_ == DIM_3D, "Wrong dimension.");
   return ar[(xSize_ * ySize_)*k + xSize_*j + i];
}



#endif //ARRAY_HH

