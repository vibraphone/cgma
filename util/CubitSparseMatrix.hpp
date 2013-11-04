//- Class: CubitSparseMatrix
//-
//- Description: This file defines the CubitSparseMatrix class which is a
//- sparse NxM Matrix.
//-
//- Author: Matt Staten
//- Data: 4/15/2011
//- Checked by:

#ifndef CUBITSPARSEMATRIX_HPP
#define CUBITSPARSEMATRIX_HPP

#include <map>
#include <vector>
#include "CubitUtilConfigure.h"
#include <stddef.h>

class CUBIT_UTIL_EXPORT CubitSparseMatrix
{
public:
  
  CubitSparseMatrix();
  CubitSparseMatrix( int numRows,
                     int numCols,
                     std::vector<int> &is,
                     std::vector<int> &js,
                     std::vector<double> &es );

  ~CubitSparseMatrix();
    //- destructor
  
  void reset( int numRows,
              int numCols,
              std::vector<int> &is,
              std::vector<int> &js,
              std::vector<double> &es );
  void clear(); // Delete all data, free all memory and set sizes to 0
  void print( char *filename = 0 ) const;

  void add( int row, int col, double data );

  int num_non_zeros( void ) const;
  void num_non_zeros_per_row( double &ave, 
                              int &max, 
                              int &min ) const;
  int num_rows( void ) const { return numRows; };
  int num_cols( void ) const { return numCols; };

  // - Get an entry int the matrix, 0 <= idx < this->num_non_zeros()
  const std::map<int,double> * get_row( int row ) const;

  // Create a matrix containing the rows and cols of this that are true in
  // rows_to_include and cols_to_include.
  void sub_matrix( const std::vector<bool> &rows_to_include,
                   const std::vector<bool> &cols_to_include,
                   CubitSparseMatrix &submatrix );

  std::vector<double> operator* (const std::vector<double> &vec ) const;

     //- Add an identity matrix into this matrix.
  void plus_identity();

private:

  void delete_data( void );

  std::map< int, std::map< int, double> * > matrixData;

  int numRows;
  int numCols;
};

#endif



