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

#ifndef _FACETEDBOOLEANIMPRINT
#define _FACETEDBOOLEANIMPRINT
#include <math.h>
#include <vector>
#include "FBDefines.hpp"
#include "FBStructs.hpp"
#include "CubitDefines.h"

class FBPolyhedron;
class FBRetriangulate;

class FBImprint {

public:
  FBImprint();
  ~FBImprint();

  CubitStatus imprint_body_curve(const std::vector<double>& Bodycoords,
                       const std::vector<int>& Bodyconnections,
                       const std::vector<FB_Coord*>& FB_imprint_edge_coords,
                       const std::vector<FB_Edge*>& FB_imprint_edges,
                       const std::vector<FSBoundingBox*>& FB_imprint_edge_bboxes,
                       std::vector<int>* indices);

  CubitStatus update_surfs_and_curves(std::vector<double>& out_coords,
                                      std::vector<int>& out_connections,
                                      std::vector<int> *out_surf_index,
                                      std::vector<int> *out_curve_index
                                      );

private:
  std::vector<int> *f_c_indices;
  FBPolyhedron *poly;
  CubitStatus edges_tri_intersect(const std::vector<FB_Coord*>& FB_imprint_edge_coords,
                        const std::vector<FB_Edge*>& FB_imprint_edges,
                        const std::vector<FSBoundingBox*>& FB_imprint_edge_bboxes,
                        bool &new_edge_created);
                        
  CubitStatus single_edge_tri_intersect(double *edge_0,double *edge_1,
                                        bool &new_edge_created,
                                        FB_Triangle *tri,
                                        bool big_angle); 
  double imprint_res;
      
};

#endif
