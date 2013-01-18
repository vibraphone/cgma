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

#ifndef _FACETEDBOOLEANINTERSECT
#define _FACETEDBOOLEANINTERSECT
#include <math.h>
#include <vector>
#include "FBDefines.hpp"
#include "FBStructs.hpp"
#include "CubitDefines.h"
#include "FBClassify.hpp"

class FBPolyhedron;
class FBRetriangulate;

class FBIntersect {

public:
  FBIntersect();

  CubitStatus intersect(const std::vector<double>& Ticoords,
                        const std::vector<int>& Ticonnections,
                        const std::vector<double>& Tjcoords,
                        const std::vector<int>& Tjconnections,
                        std::vector<int>& duddedTiFacets, 
                        std::vector<int>& duddedTjFacets,
                        std::vector<int>& newTiFacets, 
                        std::vector<int>& newTjFacets,
                        std::vector<int>& newTiFacetsIndex,
                        std::vector<int>& newTjFacetsIndex,
                        std::vector<double>& newTiPoints, 
                        std::vector<double>& newTjPoints,
                        std::vector<int>& edgesTi, 
                        std::vector<int>& edgesTj
                        );


  CubitStatus intersect(const std::vector<double>& Ticoords,
                        const std::vector<int>& Ticonnections,
                        const std::vector<double>& Tjcoords,
                        const std::vector<int>& Tjconnections,
                        std::vector<int>& newTiFacets, 
                        std::vector<int>& newTjFacets,
                        std::vector<int> *indices1,
                        std::vector<int> *indices2                        
                        );
                     
  ~FBIntersect();
  void set_classify_flag(bool value);
  CubitStatus gather_by_boolean(std::vector<double>& out_coords,
                                std::vector<int>& out_connections,
                                std::vector<int> *out_surf_index,
                                std::vector<int> *out_curve_index,
                                std::vector<bool> *is_body_1,
                                const CubitFacetboolOp op
                                );

  CubitStatus update_surfs_and_curves(std::vector<double>& out_coords,
                                      std::vector<int>& out_connections,
                                      std::vector<int> *out_surf_index,
                                      std::vector<int> *out_curve_index,
                                      const int whichone
                                      );

  CubitStatus get_persistent_entity_info(bool *surfs_in,
                                         bool *curves_in,
                                         bool *surfs_out,
                                         bool *curves_out,
                                         const CubitFacetboolOp op,
                                         const int whichparent
                                         );
                                         
  
  void set_body1_planar();
  void set_body2_planar();
  void set_imprint();
   
private:
  double linecoeff[3];
  double linept[3];
  FBPolyhedron *poly1, *poly2; 
  bool do_edges_only;
  bool do_classify;
  bool do_imprint;
  bool body1_is_plane, body2_is_plane;
  bool nothing_intersected;
  std::vector<int> *f_c_indices1, *f_c_indices2;
  FBClassify *classify1, *classify2;
  CubitStatus pair_intersect();
  int get_vertex(FBPolyhedron *poly, int vtx, 
                 IntegerHash *hashobj,
                 std::vector<double>& out_coords,
                 int &num_sofar);
  int makeahashvaluefrom_coord(double x, double y, double z);
  CubitStatus tri_tri_intersect(FB_Triangle *tri1,
                                FB_Triangle *tri2);
  CubitStatus add_intersection_edges(FB_Triangle *tri1,
              FB_Triangle *tri2,
              double *tt,
              int *edge_vert_type);
  inline void get_point_from_parameter(double parameter,
              double *x, double *y, double *z);
  inline double get_distance_parameter(double *xc0,
              double *xc1,
              double d0, double d1);
  inline double get_distance_parameter_single(double *xc);
  int get_intersectionline_parameter_values(
                    double d0, 
                    double d1, 
                    double d2,
                    double *pt0,
                    double *pt1,
                    double *pt2,
                    double& t0,
                    double& t1,
                    int& vert_type_0,
                    int& vert_type_1);
  inline int determine_edge_vert_type(int vtype1, int vtype2)
  {
    if ( ( (vtype1 == VERTEX_0) && (vtype2 == VERTEX_1) ) ||
         ( (vtype1 == VERTEX_1) && (vtype2 == VERTEX_0) ) ) 
      return EDGE_0;
    else if ( ( (vtype1 == VERTEX_1) && (vtype2 == VERTEX_2) ) || 
              ( (vtype1 == VERTEX_2) && (vtype2 == VERTEX_1) ) ) 
      return EDGE_1;
    else if ( ( (vtype1 == VERTEX_2) && (vtype2 == VERTEX_0) ) ||
              ( (vtype1 == VERTEX_0) && (vtype2 == VERTEX_2) ) ) 
      return EDGE_2;
    else return INTERIOR_VERT;
  }  
  void  newplanecoefficients(FBPolyhedron *poly, FB_Triangle *tri);

  CubitStatus store_connectivity( std::vector<int>& out_connections,
                                  int vertnum1,
                                  int vertnum2,
                                  int vertnum3 );
};

#endif
