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

#ifndef _FACETEDBOOLEANTILER
#define _FACETEDBOOLEANTILER
#include "CubitDefines.h"

#include "FBDefines.hpp"
#include "FBStructs.hpp"

class FBTilerChainVert {

public:
  FBTilerChainVert(int v_0, int which_chain) {
    v0 = v_0;
    whichchain = which_chain;
  
  }
  ~FBTilerChainVert() { };
  int v0;
  int whichchain;
};

class FBTiler {

public:

  FBTiler(std::vector<FB_Coord *>& my_verts, int pd, int sd, int sequence,
          double a, double b, double c,
          std::vector<int> *tri_list);
  ~FBTiler();
  CubitStatus Tile_Poly(std::vector<int> *coordlist);
  
private:

  int p_dir;
  int s_dir;
  int parent;
  double xnorm, ynorm, znorm;
  std::vector<FB_Coord *> verts;
  std::vector<int> *my_tri_list;  
  int add_triangle(int v1, int v2, int v3);
//  bool reflex_angle(int v0, int v1, int v2, int v1chain);
  bool reflex_angle(int v0, int v1, int v2, std::vector<int> *coordlist);
  std::vector<FBTilerChainVert*> sortedchainlist;
  CubitStatus retriangulate(std::vector<int> *coordlist);
  void dud_from_coord_list(int val, std::vector<int> *ivec);
  inline void decrement_list_ptr(std::vector<int>::iterator& it,
         std::vector<int>::iterator itbegin,  
         std::vector<int>::iterator itend)
  {
    if ( it == itbegin ) it = itend-1;
    else it--;
  }

  inline void increment_list_ptr(std::vector<int>::iterator& it,
         std::vector<int>::iterator itbegin,  
         std::vector<int>::iterator itend)
  {
    it++;
    if ( it == itend ) it = itbegin;
  }

  inline int get_adjacency(int top, int thisone)
  {
//  This code works because LEFT is defined as 1, RIGHT as 2, and BOTH as 3.
//  Thus LEFT&LEFT = RIGHT&RIGHT = LEFT&BOTH = RIGHT&BOTH = 1;
//  LEFT&RIGHT = 0.

   if ( top & thisone ) return 2;
    return 1;
  } 
   
};

#endif
