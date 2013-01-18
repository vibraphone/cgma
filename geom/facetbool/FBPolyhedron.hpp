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

#ifndef _FBPOLYHEDRON
#define _FBPOLYHEDRON
#include <vector>
#include <map>
#include "FBStructs.hpp"
#include "IntegerHash.hpp"
#include "CubitDefines.h"
#include "KdTree.hpp"

class FBPolyhedron {

friend class FBIntersect;
friend class FBClassify;
friend class FBImprint;

public:
  FBPolyhedron();
  ~FBPolyhedron();
  CubitStatus makepoly(const std::vector<double>& coords,
                       const std::vector<int>& connections,
                       std::vector<int> *f_c_indices);
  bool boxintersection(FBPolyhedron *otherpoly);
  CubitStatus retriangulate(std::vector<int>& newfacets, 
                            std::vector<int>& newfacetsindex);
  CubitStatus retriangulate(std::vector<int>& newfacets);
bool edge_exists_in_tri(FB_Triangle& tri, int v0, int v1);

private:
  std::vector<FB_Coord *> verts;
  std::vector<FB_Triangle *> tris;
  std::vector<FB_Edge *> intersection_edges;
  std::multimap<unsigned int,unsigned int> edgmmap;
  std::vector<int> goodtris;
  double polyxmin, polyxmax, polyymin, polyymax, polyzmin, polyzmax;
  IntegerHash *hashobj;
  int makeahashvaluefrom_coord(double x, double y, double z);
  int addavertex(double x, double y, double z);
  unsigned int numpts, numtris;
  void putnewpoints(std::vector<double>& newpoints);
  void putedges(std::vector<int>& newedges); 
  bool edge_exists(int v0, int v1); 
  KDTree *kdtree;
  int original_numtris;
//  void putnewtriangles(std::vector<int>& newFacets);
  void add_new_triangle_data();
  void make_tri_plane_coeffs(FB_Triangle *tri);
  void make_tri_boundingbox(FB_Triangle *tri);
  void removeduddedtriangles();
  
    //debug functions
    //!Find the min and max angles in the given Triangle.
  bool min_max_angles_in_fb_triangle(FB_Triangle *tri, double& min_angle,
                                    double& max_angle );
    //!Find the min and max angles in the polyhedron (ie, min and max angle
    //! in any triangle).
  bool min_max_angles_in_polyhedron(double& min_angle, double& max_angle);
    //! Debug draw a triangle.
  void debug_draw_fb_triangle(FB_Triangle *tri);
    //! Debug draw the edges in the polyhedron that are marked boundary (ie,
    //! draw the edges that have a cubit_edge).
  void debug_draw_boundary_edges(int color);
    //! This a function that attempts to remove small angles from the
    //! polyhedron.
  bool remove_small_angles_in_triangle_range(int lower_index,
                                             int upper_index);
  
};

#endif
