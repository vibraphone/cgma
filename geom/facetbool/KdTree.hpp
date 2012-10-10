/*
 *
 *
 * Copyright (C) 2004 Sandia Corporation.  Under the terms of Contract DE-AC04-94AL85000
 * with Sandia Corporation, the U.S. Government retains certain rights in this software.
 *
 * This file is part of facetbool--contact via cubit@sandia.gov
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 *
 *
 */

#ifndef _KDTREE
#define _KDTREE
#include <vector>

#include "FBStructs.hpp"

class TreeStack {

public:

  TreeStack(int imin, int imax, int cut_dir, int iseq) {
   cuttingdir = cut_dir;
    min = imin;
    max = imax;
    sequence = iseq;
  }
  
  ~TreeStack() { }
  
  int min, max, cuttingdir, sequence;
};

class KDTree
{

public:

  KDTree();
  ~KDTree();
  void box_kdtree_intersect(FSBoundingBox& bbox, int *count, int *indexlist) const;
  int makeKDTree(int npoly, const FSBOXVECTOR& boxlist);

private:

  double epsilonkd;
  int numtris;
  FSBoundingBox* *treebox; 
  int *nextbranch;
  FSBOXVECTOR boxlistptr;
//  std::vector<FSBoundingBox* > boxlist;
  int *isequence;
  FSBoundingBox* getbox(int min, int max);
  void find_the_median(int k, int l, int r, double *array, int *ia);
  int getcuttingdirection(FSBoundingBox* box);
  inline void SWAP(int &x, int &y)
  {
  int temp;
    temp = x;
    x = y;
    y = temp;
  }   
  double rayxstart, rayystart, rayzstart, dx, dy, dz;
  bool rayintersectsbox(FSBoundingBox *box);
  inline double MAXX(double a, double b) {
    if ( a > b ) return a;
    else return b;
  }
  inline double MINN(double a, double b) {
    if ( a < b ) return a;
    else return b;
  }
};

#endif
