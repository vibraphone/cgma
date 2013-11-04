//- Class: CubitSparseMatrix
//-
//- Description: This file defines the CubitSparseMatrix class which is a
//- sparse NxM Matrix.
//-
//- Author: Matt Staten
//- Data: 4/15/2011
//- Checked by:

#include "CubitSparseMatrix.hpp"
#include "CubitMessage.hpp"
#include "CubitDefines.h"
#include "CubitFileUtil.hpp"

#include <vector>
#include <set>
#include <stdexcept>
#include <stdio.h>

using std::vector;
using std::set;
using std::map;

//==============================================================================
// Description:  Constructor
// Notes:  
// Author: mlstate
// Date: 4/19/2011
//==============================================================================
CubitSparseMatrix::CubitSparseMatrix
(
  int num_rows,
  int num_cols,
  vector<int> &is,
  vector<int> &js,
  vector<double> &es
)
{
  numRows = numCols = 0;
  reset( num_rows, num_cols, is, js, es );
}

//==============================================================================
// Description: Constructor
// Notes:  
// Author: mlstate
// Date: 4/19/2011
//==============================================================================
CubitSparseMatrix::CubitSparseMatrix()
{
  numRows = numCols = 0;
}

//==============================================================================
// Description: 
// Notes:  
// Author: mlstate
// Date: 4/19/2011
//==============================================================================
void CubitSparseMatrix::reset
(
  int num_rows,
  int num_cols,
  vector<int> &is,
  vector<int> &js,
  vector<double> &es
)
{

  if( is.size() != js.size() ||
      is.size() != es.size() )
  {
    throw std::invalid_argument("The sizes of is, js, and es must be the same" );
  }

  delete_data();

  int num_data = is.size();

  numRows = num_rows;
  numCols = num_cols;
  for ( int i = 0; i < num_data; i++ )
  {
    add( is[i], js[i], es[i] );
  }
}

//==============================================================================
// Description: Delete all data, free all memory, and set size to 0.
// Notes:  
// Author: mlstate
// Date: 4/19/2011
//==============================================================================
void CubitSparseMatrix::clear()
{
  numRows = 0;
  numCols = 0;
  delete_data();
}

//==============================================================================
// Description: 
// Notes:  
// Author: mlstate
// Date: 4/19/2011
//==============================================================================
void CubitSparseMatrix::add( int row, int col, double data )
{
  if ( row+1 > numRows ||
       col+1 > numCols )
  {
    throw std::invalid_argument( "The dimensions of the 2 matrices must be the same." );
  }

  map<int,double> *submap = NULL;
  map<int,map<int,double>*>::iterator iter = matrixData.find( row );
  if ( iter == matrixData.end() )
  {
    submap = new map<int,double>;
    matrixData.insert( std::pair<int, map<int,double>*>(row, submap) );
  }
  else
  {
    submap = iter->second;
    map<int,double>::iterator iter2 = submap->find( col );
    if ( iter2 != submap->end() )
    {
      iter2->second += data;
      return;
    }
  }
  submap->insert( std::pair<int, double>( col, data ) );
}

//==============================================================================
// Description:  destructor
// Notes:  
// Author: mlstate
// Date: 4/19/2011
//==============================================================================
CubitSparseMatrix::~CubitSparseMatrix()
{
  delete_data();
}

//==============================================================================
// Description:  Delete all data in this matrix.
// Notes:  
// Author: mlstate
// Date: 4/19/2011
//==============================================================================
void CubitSparseMatrix::delete_data()
{
  while ( !matrixData.empty() )
  {
    map<int,double> *submap = matrixData.begin()->second;
    matrixData.begin()->second = NULL;
    matrixData.erase( matrixData.begin() );
    delete submap;
  }
}
  
  // Create a matrix containing the rows and cols of this that are true in
  // rows_to_include and cols_to_include.
//==============================================================================
// Description:  retrieve specified portions of this matrix.
// Notes:  
// Author: mlstate
// Date: 4/19/2011
//==============================================================================
void CubitSparseMatrix::sub_matrix
(
  const vector<bool> &rows_to_include,
  const vector<bool> &cols_to_include,
  CubitSparseMatrix &submatrix
)
{
  if ( (int) rows_to_include.size() != numRows )
  {
    throw std::invalid_argument( "The length of rows_to_include must the the same as the number of rows in the matrix." );
  }
  if ( (int) cols_to_include.size() != numCols )
  {
    throw std::invalid_argument( "The length of cols_to_include must the the same as the number of colums in the matrix." );
  }

  vector<int> is;
  vector<int> js;
  vector<double> es;
  int i;

  int *row_indices = new int[numRows];
  int *col_indices = new int[numCols];
  int num_rows = 0;
  int num_cols = 0;
  for ( i = 0; i < numRows; i++ )
  {
    if ( rows_to_include[i] )
    {
      row_indices[i] = num_rows;
      num_rows++;
    }
    else
    {
      row_indices[i] = 0;
    }
  }

  for ( i = 0; i < numCols; i++ )
  {
    if ( cols_to_include[i] )
    {
      col_indices[i] = num_cols;
      num_cols++;
    }
    else
    {
      col_indices[i] = 0;
    }
  }

  map<int, map<int,double>*>::iterator iter1 = matrixData.begin();
  while ( iter1 != matrixData.end() )
  {
    int row = iter1->first;
    map<int,double>*submap = iter1->second;
    iter1++;

    if ( !rows_to_include[row] )
      continue;

    map<int,double>::iterator iter2 = submap->begin();
    while ( iter2 != submap->end() )
    {
      int col = iter2->first;
      double data = iter2->second;
      iter2++;
      if ( !cols_to_include[col] )
        continue;

      is.push_back( row_indices[row] );
      js.push_back( col_indices[col] );
      es.push_back( data );
    }
  }

  delete [] col_indices;
  delete [] row_indices;

  submatrix.reset( num_rows, num_cols, is, js, es );
}

//==============================================================================
// Description:  Multiply this matrix by a vector.
// Notes:  
// Author: mlstate
// Date: 4/19/2011
//==============================================================================
vector<double> CubitSparseMatrix::operator* (const vector<double> &vec ) const
{
  if ( (int) vec.size() != numCols )
  {
    throw std::invalid_argument( "The length of the input vector must be the same as the number of colums in the matrix." );
  }

  vector<double> ans( numRows );

  int i;
  for ( i = 0; i < numRows; i++ )
  {
    ans[i] = 0.0;
  }

  map<int, map<int,double>*>::const_iterator iter1 = matrixData.begin();
  while ( iter1 != matrixData.end() )
  {
    int row = iter1->first;
    map<int,double>*submap = iter1->second;
    iter1++;

    map<int,double>::iterator iter2 = submap->begin();
    while ( iter2 != submap->end() )
    {
      int col = iter2->first;
      double data = iter2->second;
      iter2++;

      ans[row] += data * vec[col];
    }
  }
  return ans;
}

//==============================================================================
// Description:  retrieve a row of this matrix.
// Notes:  
// Author: mlstate
// Date: 4/19/2011
//==============================================================================
const map<int,double> *CubitSparseMatrix::get_row( int row ) const
{
  if ( row < 0 || row >= numRows )
  {
    throw std::invalid_argument( "The row index must be greater than or equal to 0, and less than the number of rows in the matrix." );
  }

  map<int, map<int,double>*>::const_iterator iter = matrixData.find( row );
  if ( iter == matrixData.end() )
    return NULL;
  return iter->second;
}

//==============================================================================
// Description: debug print of this matrix.
// Notes:  
// Author: mlstate
// Date: 4/19/2011
//==============================================================================
void CubitSparseMatrix::print( char *filename ) const
{
  printf( "CubitSparseMatrix::numRows: %d\n", numRows );
  printf( "CubitSparseMatrix::numCols: %d\n", numCols );
  int min, max;
  double ave;
  num_non_zeros_per_row( ave, max, min );

  printf( "CubitSparseMatrix::Average # Non Zeros per row: %lf\n", ave );
  printf( "CubitSparseMatrix::Minimum # Non Zeros per row: %d\n", min );
  printf( "CubitSparseMatrix::Maximum # Non Zeros per row: %d\n", max );

  CubitFile fp;
  if ( filename )
    fp.open( filename, "w" );

  map<int, map<int,double>*>::const_iterator iter1 = matrixData.begin();
  while ( iter1 != matrixData.end() )
  {
    int row = iter1->first;
    map<int,double>*submap = iter1->second;
    iter1++;

    map<int,double>::iterator iter2 = submap->begin();
    while ( iter2 != submap->end() )
    {
      int col = iter2->first;
      double data = iter2->second;
      iter2++;

      if ( fp.file() )
      {
        fprintf( fp.file(), "%d %d %lf\n",
                row, col, data );
      }
      else
      {
      printf( "CubitSparseMatrix::(%d,%d) : %lf\n", 
              row, col, data );
    }
  }
  }
}

//==============================================================================
// Description: debug print of this matrix.
// Notes:  
// Author: mlstate
// Date: 4/19/2011
//==============================================================================
void CubitSparseMatrix::num_non_zeros_per_row
(
  double &ave,
  int &max,
  int &min
) const
{
  max = 0;
  min = numCols;
  ave = 0.0;
  
  map<int, map<int,double>*>::const_iterator iter1 = matrixData.begin();
  while ( iter1 != matrixData.end() )
  {
    iter1->first;
    map<int,double>*submap = iter1->second;
    iter1++;

    int nz = submap->size();
    ave += nz;

    if ( nz > max ) max = nz;
    if ( nz < min ) min = nz;
  }
  if ( numCols > 0 )
    ave /= numCols;
}

//==============================================================================
// Description: Add an identity matrix to this matrix.
// Notes:  
// Author: mlstate
// Date: 4/19/2011
//==============================================================================
void CubitSparseMatrix::plus_identity()
{
  for ( int i = 0; i < numRows; i++ )
  {
    if ( i > numCols ) break;
    add( i, i, 1.0 );
  }
}

//==============================================================================
// Description: return the number of non-zeros in this matrix.
// Notes:  
// Author: mlstate
// Date: 4/19/2011
//==============================================================================
int CubitSparseMatrix::num_non_zeros( void ) const
{
  int num_nz = 0;
  map<int, map<int,double>*>::const_iterator iter1 = matrixData.begin();
  while ( iter1 != matrixData.end() )
  {
    iter1->first;
    map<int,double>*submap = iter1->second;
    if ( submap ) num_nz += submap->size();
    iter1++;
  }
  return num_nz;
}

