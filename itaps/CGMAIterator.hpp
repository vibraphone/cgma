/**
 * Copyright 2006 Sandia Corporation.  Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Coroporation, the U.S. Government
 * retains certain rights in this software.
 * 
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * 
 */
/* (2006) kraftche@cae.wisc.edu */
/* Renamed from TSTTGIterator to CGMAIterator - J.Kraftcheck - 2007/6/15 */

#ifndef CGMA_ITERATOR_HPP
#define CGMA_ITERATOR_HPP

class RefEntity;

#include <DLIList.hpp>
#include <vector>
#include <algorithm>

class CGMAIterator
{
public:

  /**\brief Constructor
   *
   *\param list  The list of entities to iterate over
   *\param size  Default query size - stored but not used by this class
   */
  CGMAIterator( DLIList<RefEntity*>& list, int size ) 
    : mArray(list.size()),
      mIter(mArray.begin()),
      mSize(size)
  {
    list.copy_to( &mArray[0] );
  }
  
  /**\brief Reset iterator */
  void reset() { mIter = mArray.begin(); }
  
  /**\brief Get and step 
   *
   *\param array  allocated array of pointers to RefEntitys.
   *              must be allocated by caller and contain at least
   *              array_size entries.
   *\param array_size  As input, the requested number of entries to
   *              retreive and step over.  As output, either the
   *              input value or a smaller value if there are insufficient
   *              entries remaining.
   */
  bool next( RefEntity** array, int& array_size )
  {
    if (array_size <= 0) {
      array_size = 0;
    }
    else {
      std::vector<RefEntity*>::const_iterator end = mArray.end();
      if ((end - mIter) > array_size)
        end = mIter + array_size;
      RefEntity** ptr = std::copy( mIter, end, array );
      array_size = ptr - array;
      mIter += array_size;
    }
    return mIter == mArray.end();
  }
 
  /**\brief Get the saved default request size 
   *        NOTE:  This is *not* the number of entities iterated over.
   *
   * Get saved default request size.  This value is saved in the
   * iterator class for convenience.  It has no affect on the
   * behavior of the class.  
   *
   * Note:  This is neither the number of entities iterated over 
   *        nor the remaining number of entities.
   */        
  int size() const { return mSize; }
  
  bool at_end() const
    { return mIter == mArray.end(); }

private:

  std::vector<RefEntity*> mArray;
  std::vector<RefEntity*>::const_iterator mIter;
  int mSize;
};

#endif

