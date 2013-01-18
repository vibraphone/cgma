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

#ifndef _FACETEDBOOLEANRETRIANGULATE
#define _FACETEDBOOLEANRETRIANGULATE
#include <math.h>
#include <vector>
#include "FBDefines.hpp"
#include "FBStructs.hpp"
#include "CubitDefines.h"

class FBRetriangulate {

public:
  FBRetriangulate(std::vector<FB_Coord *>& my_verts,
                  std::vector<FB_Triangle *>& my_tris,
                  std::vector<int>& my_newfacets,
                  std::vector<int>& my_newfacetsindex);

  FBRetriangulate(std::vector<FB_Coord *>& my_verts,
                  std::vector<FB_Triangle *>& my_tris,
                  std::vector<int>& my_newfacets);

  ~FBRetriangulate();
  CubitStatus retriangulate_this_tri(int sequence);
  
private:
  std::vector<FB_Coord *> verts;
  std::vector<FB_Triangle *> *tris;
  std::vector<int> *newfacets;
  std::vector<int> *newfacetsindex;  
  int p_dir, s_dir;
  enum WindingType { CCW = 0, CW = 1 };
  WindingType winding;
  FB_Triangle *my_tri;
  std::vector<VertexStuff* > vertstufflist;
  void make_vert_list();
  CubitStatus add_bdry_edges();  
  inline double get_dist(double x0, double y0, double z0,
                                 double x1, double y1, double z1) {
    return sqrt( (x1-x0)*(x1-x0) + (y1-y0)*(y1-y0) + (z1-z0)*(z1-z0) );
  };
  CubitStatus remove_min_max();
  void add_edge_up(int v0, int seq);
  void add_edge_down(int v0, int seq);
  bool test_for_crossing(int v0, int v1, int v2, int v3);
  void add_tri_edge(int v0, int v1, int v0_type, int v1_type);
  bool get_a_chain(std::vector<int> *chainlist);
  void sort_vertstufflist_edges();
  void get_next_vert_in_chain(int& vthis, int& vprev, int direction);
  void add_this_bdry_edge(int v0, int v1, int v0_type, int v1_type); 
  void classify_edges();
  void get_edge_indices(int v0, int v1, int v2, int parent, 
                        int &e0index, int &e1index, int &e2index);
                          
};

#endif
