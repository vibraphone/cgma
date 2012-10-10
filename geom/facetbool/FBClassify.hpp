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

#ifndef _FACETEDBOOLEANCLASSIFY
#define _FACETEDBOOLEANCLASSIFY
#include "FBDefines.hpp"
#include "FBStructs.hpp"
#include "FBPolyhedron.hpp"
#include "CubitDefines.h"
#include <vector>
#include <map>

const int NO_EDGE_NBR = -99999;

class FBClassify {

public:
  FBClassify();
  ~FBClassify();
  
  void SetPoly(FBPolyhedron *poly1, FBPolyhedron *poly2);
  CubitStatus Group(int which);
  CubitStatus CharacterizeGroups(int which, bool other_is_planar);
  void get_group(std::vector<int> **this_group,
                 std::vector<int> **this_group_characterization);

private:
  FBPolyhedron *polya, *polyb;
  std::vector<int> group, group_characterization;
  void fill_group(int itri, int ngroup);
  void perturb_the_ray(double &xbary, double &ybary, double &zbary);
  int pt_in_tri_2d(double xpt, double ypt,
                           double x0, double y0,
		           double x1, double y1,
		           double x2, double y2);
  int *e0, *e1, *e2;
  int number_of_groups;
  int classify(int itri, int which);
  int classify_against_plane(int itri, int which);
};
#endif
