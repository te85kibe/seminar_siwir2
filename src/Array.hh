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
   real & operator () ( int i );
   real & operator () ( int i ,int j );
   real & operator () ( int i, int j, int k );

   // for const Arrays the following access operators are required
   const real & operator () ( int i ) const;
   const real & operator () ( int i ,int j ) const;
   const real & operator () ( int i, int j, int k ) const;



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





#endif //ARRAY_HH

