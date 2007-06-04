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

#include <vector>
#include "FBIntersect.hpp"
#include "FBPolyhedron.hpp"
#include "FBRetriangulate.hpp"
#include "CubitMessage.hpp"
#include "GfxDebug.hpp"
#include "GeometryDefines.h"

FBIntersect::FBIntersect()
{
  poly1 = poly2 = 0;
  do_classify = false;
  do_edges_only = false;
  do_imprint = false;
  classify1 = classify2 = 0;
  body1_is_plane = body2_is_plane = false;
  f_c_indices1 = f_c_indices2 = 0;
  nothing_intersected = false;
}

FBIntersect::~FBIntersect()
{
  if ( poly1 ) delete poly1;
  if ( poly2 ) delete poly2;
  if ( classify1 ) delete classify1;
  if ( classify2 ) delete classify2;
}

CubitStatus FBIntersect::intersect(const std::vector<double>& Ticoords,
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
                     )
{
CubitStatus status;
bool boxes_intersect;

  status = CUBIT_SUCCESS;   

  if ( Ticoords.size()%3 != 0 ) {
    PRINT_ERROR("Bad coordinates for first part fed to FBIntersect.\n");
    return CUBIT_FAILURE;
  }
  if ( Tjcoords.size()%3 != 0 ) {
    PRINT_ERROR("Bad coordinates for second part fed to FBIntersect.\n");
    return CUBIT_FAILURE;
  }
  if ( Ticonnections.size()%3 != 0 ) {
    PRINT_ERROR("Bad connection list for first part fed to FBIntersect.\n");
    return CUBIT_FAILURE;
  }
  if ( Tjconnections.size()%3 != 0 ) {
    PRINT_ERROR("Bad connection list for second part fed to FBIntersect.\n");
    return CUBIT_FAILURE;
  }    
     
  poly1 = new FBPolyhedron;
  status = poly1->makepoly(Ticoords,Ticonnections,f_c_indices1);
  if ( status != CUBIT_SUCCESS )
  {
      return status;
  }

  poly2 = new FBPolyhedron;
  status = poly2->makepoly(Tjcoords,Tjconnections,f_c_indices2);
  if ( status != CUBIT_SUCCESS )
  {
      return status;
  }
  
  boxes_intersect = poly1->boxintersection(poly2);
  if ( boxes_intersect == false ) {
    nothing_intersected = true;
    return status;
  }
  
  status = pair_intersect();  
  if ( status != CUBIT_SUCCESS )
  {
      return status;
  }

  poly1->putnewpoints(newTiPoints);
  poly2->putnewpoints(newTjPoints);

  poly1->putedges(edgesTi);
  poly2->putedges(edgesTj);
  
  status = poly1->retriangulate(newTiFacets, newTiFacetsIndex);
  if ( status != CUBIT_SUCCESS )
  {
      return status;
  }

  status = poly2->retriangulate(newTjFacets, newTjFacetsIndex);
  if ( status != CUBIT_SUCCESS )
  {
      return status;
  }

  unsigned int i;
  for ( i = 0; i < poly1->tris.size(); i++ ) {
    if ( poly1->tris[i]->dudded == true ) duddedTiFacets.push_back(i);
  }

  for ( i = 0; i < poly2->tris.size(); i++ ) {
    if ( poly2->tris[i]->dudded == true ) duddedTjFacets.push_back(i);
  }
  
//  newxxFacets is updated in FBRetriangulate::make_this_tri() after
//  a new triangle has been made.  Also newxxFacetsIndex is updated if
//  any new facets were made.  But this update points to the start of
//  newxxFacets.  Therefore we have to add on the last member to 
//  newxxFacetsIndex here. 
  newTiFacetsIndex.push_back(newTiFacets.size());
  newTjFacetsIndex.push_back(newTjFacets.size());    
     
//  If classification of the intersected surface parts is needed, we need to
//  do some things.
  bool other_is_planar;
  if ( do_classify == true ) {
    poly1->add_new_triangle_data();
    poly2->add_new_triangle_data();
    poly1->removeduddedtriangles();
    poly2->removeduddedtriangles();    
    classify1 = new FBClassify;
    classify1->SetPoly(poly1, poly2);    
    status = classify1->Group(1);
    if ( status != CUBIT_SUCCESS )
    {
      return status;
    }
    other_is_planar = body2_is_plane;
    status = classify1->CharacterizeGroups(1,other_is_planar);
    if ( status != CUBIT_SUCCESS )
    {
      return status;
    }
    classify2 = new FBClassify;
    classify2->SetPoly(poly1, poly2);
    status = classify2->Group(2);
    if ( status != CUBIT_SUCCESS )
    {
      return status;
    }
    other_is_planar = body1_is_plane;
    status = classify2->CharacterizeGroups(2,other_is_planar);
    if ( status != CUBIT_SUCCESS )
    {
      return status;
    }
  }  
  
  return status;

}

CubitStatus FBIntersect::intersect(const std::vector<double>& Ticoords,
                     const std::vector<int>& Ticonnections,
                     const std::vector<double>& Tjcoords,
                     const std::vector<int>& Tjconnections,
                     std::vector<int>& newTiFacets, 
                     std::vector<int>& newTjFacets, 
                     std::vector<int> *indices1,
                     std::vector<int> *indices2
                     )
{
CubitStatus status;
bool boxes_intersect;

  status = CUBIT_SUCCESS;   

  if ( Ticoords.size()%3 != 0 ) {
    PRINT_ERROR("Bad coordinates for first part fed to FBIntersect.\n");
    return CUBIT_FAILURE;
  }
  if ( Tjcoords.size()%3 != 0 ) {
    PRINT_ERROR("Bad coordinates for second part fed to FBIntersect.\n");
    return CUBIT_FAILURE;
  }
  if ( Ticonnections.size()%3 != 0 ) {
    PRINT_ERROR("Bad connection list for first part fed to FBIntersect.\n");
    return CUBIT_FAILURE;
  }
  if ( Tjconnections.size()%3 != 0 ) {
    PRINT_ERROR("Bad connection list for second part fed to FBIntersect.\n");
    return CUBIT_FAILURE;
  }    
  
  f_c_indices1 = indices1;
  f_c_indices2 = indices2;
   
  poly1 = new FBPolyhedron;
  status = poly1->makepoly(Ticoords,Ticonnections,f_c_indices1);
  if ( status != CUBIT_SUCCESS )
  {
      return status;
  }

  poly2 = new FBPolyhedron;
  status = poly2->makepoly(Tjcoords,Tjconnections,f_c_indices2);
  if ( status != CUBIT_SUCCESS )
  {
      return status;
  }

  double min_angle=0.0, max_angle=0.0;
  int mydebug = 0;
  if(mydebug){
    poly1->min_max_angles_in_polyhedron(min_angle, max_angle);
    
    PRINT_INFO("(0) Min angle in poly1 %f, max angle  %f\n", min_angle,
               max_angle);
    poly2->min_max_angles_in_polyhedron(min_angle, max_angle);
    PRINT_INFO("(0) Min angle in poly2 %f, max angle  %f\n", min_angle,
               max_angle);
  }
  boxes_intersect = poly1->boxintersection(poly2);
  if ( boxes_intersect == false ) {
    nothing_intersected = true;
      //don't return because in the case of union we need to unite the
      //bodies even if they do not intersect
      //return status;
  }
  else{
    status = pair_intersect();
    if ( status != CUBIT_SUCCESS )
    {
        return status;
    }
  }
  
  if(mydebug){
    poly1->min_max_angles_in_polyhedron(min_angle, max_angle);
    PRINT_INFO("(1) Min angle in poly1 %f, max angle  %f\n", min_angle,
               max_angle);
    poly2->min_max_angles_in_polyhedron(min_angle, max_angle);
    PRINT_INFO("(1) Min angle in poly2 %f, max angle  %f\n", min_angle,
               max_angle);
  }
  
  status = poly1->retriangulate(newTiFacets);
  if ( status != CUBIT_SUCCESS )
  {
      return status;
  }
  status = poly2->retriangulate(newTjFacets);
  if ( status != CUBIT_SUCCESS )
  {
      return status;
  }

  if(mydebug){
    GfxDebug::clear();
    poly1->debug_draw_boundary_edges(CUBIT_BLUE);
    GfxDebug::mouse_xforms();
    GfxDebug::clear();
    poly2->debug_draw_boundary_edges(CUBIT_RED);
    GfxDebug::mouse_xforms();
    poly1->min_max_angles_in_polyhedron(min_angle, max_angle);
    
    PRINT_INFO("(2) Min angle in poly1 %f, max angle  %f\n", min_angle,
               max_angle);
    poly2->min_max_angles_in_polyhedron(min_angle, max_angle);
    PRINT_INFO("(2) Min angle in poly2 %f, max angle  %f\n", min_angle,
               max_angle);
  }
//  If classification of the intersected surface parts is needed, we need to
//  do some things.
  bool other_is_planar;
  if ( do_classify == true ) {
    poly1->add_new_triangle_data();
    poly2->add_new_triangle_data();
    poly1->removeduddedtriangles();
    poly2->removeduddedtriangles();    
    classify1 = new FBClassify;
    classify1->SetPoly(poly1, poly2);    
    status = classify1->Group(1);
    if ( status != CUBIT_SUCCESS )
    {
      return status;
    }
    other_is_planar = body2_is_plane;
    status = classify1->CharacterizeGroups(1,other_is_planar);
    if ( status != CUBIT_SUCCESS )
    {
      return status;
    }
    classify2 = new FBClassify;
    classify2->SetPoly(poly1, poly2);
    status = classify2->Group(2);
    if ( status != CUBIT_SUCCESS )
    {
      return status;
    }
    other_is_planar = body1_is_plane;
    status = classify2->CharacterizeGroups(2,other_is_planar);
    if ( status != CUBIT_SUCCESS )
    {
      return status;
    }
  }  
  
  return status;

}

CubitStatus FBIntersect::pair_intersect()
{
CubitStatus status;
unsigned int i, j, k;
int numboxesfound, *boxlist;
double dtemp;

  boxlist = new int[poly2->tris.size()];
  
  status = CUBIT_SUCCESS;
  for ( i = 0; i < poly1->tris.size(); i++ ) {
     FB_Triangle *tri1 = poly1->tris[i];
     if ( (tri1->boundingbox.xmax < poly2->polyxmin) || 
          (tri1->boundingbox.xmin > poly2->polyxmax) ||
          (tri1->boundingbox.ymax < poly2->polyymin) || 
          (tri1->boundingbox.ymin > poly2->polyymax) ||
          (tri1->boundingbox.zmax < poly2->polyzmin) || 
          (tri1->boundingbox.zmin > poly2->polyzmax) ) continue;
     poly2->kdtree->box_kdtree_intersect(tri1->boundingbox,&numboxesfound,boxlist);

     for ( j = 0; j < (unsigned int)numboxesfound; j++ ) { 
        FB_Triangle *tri2 = poly2->tris[boxlist[j]];     
        if ( (tri1->boundingbox.xmax < tri2->boundingbox.xmin) || 
             (tri1->boundingbox.xmin > tri2->boundingbox.xmax) ||
             (tri1->boundingbox.ymax < tri2->boundingbox.ymin) || 
             (tri1->boundingbox.ymin > tri2->boundingbox.ymax) ||
             (tri1->boundingbox.zmax < tri2->boundingbox.zmin ) || 
             (tri1->boundingbox.zmin > tri2->boundingbox.zmax) ) continue;

     int mydebug = 0;
     if ( mydebug )
     {
         // Draw the 2 facets which are intersecting.
         GfxDebug::clear();
         poly1->debug_draw_boundary_edges(CUBIT_YELLOW);
         poly2->debug_draw_boundary_edges(CUBIT_PINK);
         poly1->debug_draw_fb_triangle(tri1);
         poly2->debug_draw_fb_triangle(tri2);
         GfxDebug::mouse_xforms();
     }

//  Find the line of intersection of the two triangle planes. 
        linecoeff[0] = tri1->c*tri2->b - tri1->b*tri2->c;
        linecoeff[1] = tri1->a*tri2->c - tri1->c*tri2->a;
        linecoeff[2] = tri1->b*tri2->a - tri1->a*tri2->b;

        if ( (fabs(linecoeff[0]) < EPSILON) && 
             (fabs(linecoeff[1]) < EPSILON) &&  
             (fabs(linecoeff[2]) < EPSILON) ) {
            if ( do_imprint == true ) continue;
            //  Don't intersect nearly-coplanar triangles.

            // calculate the distance between the triangles.  Just because they are coplanar
            // does not mean we have to intersect.  If the distance between them is larger than
            // GEOMETRY_RESABS, then just continue on.
            double dist = poly1->verts[tri1->v2]->coord[0]*tri2->a +
                          poly1->verts[tri1->v2]->coord[1]*tri2->b +
                          poly1->verts[tri1->v2]->coord[2]*tri2->c + tri2->d;
            if ( dist > GEOMETRY_RESABS )
            {
                continue;
            }

            //  coplanar triangles; tilt each vertex of the triangle up successively
            //  and then do the usual tri-tri intersection.
            double ta, tb, tc, td;
            double tx1[3],ty1[3],tz1[3];

            //  triangle 1
            ta = tri1->a; tb = tri1->b; tc = tri1->c; td = tri1->d;
            for ( k = 0; k < 3; k++ ) {
              tx1[k] = poly1->verts[tri1->v0]->coord[k];
              ty1[k] = poly1->verts[tri1->v1]->coord[k];
              tz1[k] = poly1->verts[tri1->v2]->coord[k];
            } 

            poly1->verts[tri1->v2]->coord[0] += ta;
            poly1->verts[tri1->v2]->coord[1] += tb;
            poly1->verts[tri1->v2]->coord[2] += tc;
            //  Compute new normal and linecoeff
            newplanecoefficients(poly1, tri1);
            linecoeff[0] = tri1->c*tri2->b - tri1->b*tri2->c;
            linecoeff[1] = tri1->a*tri2->c - tri1->c*tri2->a;
            linecoeff[2] = tri1->b*tri2->a - tri1->a*tri2->b;
            dtemp = sqrt(linecoeff[0]*linecoeff[0] +
                         linecoeff[1]*linecoeff[1] +
                         linecoeff[2]*linecoeff[2]);
            linecoeff[0] /= dtemp;
            linecoeff[1] /= dtemp;
            linecoeff[2] /= dtemp;

            status = tri_tri_intersect(tri1,tri2);
            if ( status != CUBIT_SUCCESS )
            {
                return status;
            }

            //  Restore values.
            poly1->verts[tri1->v2]->coord[0] = tz1[0];
            poly1->verts[tri1->v2]->coord[1] = tz1[1];
            poly1->verts[tri1->v2]->coord[2] = tz1[2];
            tri1->a = ta; tri1->b = tb; tri1->c = tc; tri1->d = td;

            poly1->verts[tri1->v1]->coord[0] += ta;
            poly1->verts[tri1->v1]->coord[1] += tb;
            poly1->verts[tri1->v1]->coord[2] += tc;
            //  Compute new normal and linecoeff
            newplanecoefficients(poly1, tri1);
            linecoeff[0] = tri1->c*tri2->b - tri1->b*tri2->c;
            linecoeff[1] = tri1->a*tri2->c - tri1->c*tri2->a;
            linecoeff[2] = tri1->b*tri2->a - tri1->a*tri2->b;
            dtemp = sqrt(linecoeff[0]*linecoeff[0] +
                         linecoeff[1]*linecoeff[1] +
                         linecoeff[2]*linecoeff[2]);
            linecoeff[0] /= dtemp;
            linecoeff[1] /= dtemp;
            linecoeff[2] /= dtemp;

            status = tri_tri_intersect(tri1,tri2);
            if ( status != CUBIT_SUCCESS )
            {
                return status;
            }

            //  Restore values.
            poly1->verts[tri1->v1]->coord[0] = ty1[0];
            poly1->verts[tri1->v1]->coord[1] = ty1[1];
            poly1->verts[tri1->v1]->coord[2] = ty1[2];
            tri1->a = ta; tri1->b = tb; tri1->c = tc; tri1->d = td;

            poly1->verts[tri1->v0]->coord[0] += ta;
            poly1->verts[tri1->v0]->coord[1] += tb;
            poly1->verts[tri1->v0]->coord[2] += tc;
            //  Compute new normal and linecoeff
            newplanecoefficients(poly1, tri1);
            linecoeff[0] = tri1->c*tri2->b - tri1->b*tri2->c;
            linecoeff[1] = tri1->a*tri2->c - tri1->c*tri2->a;
            linecoeff[2] = tri1->b*tri2->a - tri1->a*tri2->b;
            dtemp = sqrt(linecoeff[0]*linecoeff[0] +
                         linecoeff[1]*linecoeff[1] +
                         linecoeff[2]*linecoeff[2]);
            linecoeff[0] /= dtemp;
            linecoeff[1] /= dtemp;
            linecoeff[2] /= dtemp;

            status = tri_tri_intersect(tri1,tri2);
            if ( status != CUBIT_SUCCESS )
            {
                return status;
            }

            //  Restore values.
            poly1->verts[tri1->v0]->coord[0] = tx1[0];
            poly1->verts[tri1->v0]->coord[1] = tx1[1];
            poly1->verts[tri1->v0]->coord[2] = tx1[2];
            tri1->a = ta; tri1->b = tb; tri1->c = tc; tri1->d = td;

            //  triangle 2
            ta = tri2->a; tb = tri2->b; tc = tri2->c; td = tri2->d;
            for ( k = 0; k < 3; k++ ) {
              tx1[k] = poly2->verts[tri2->v0]->coord[k];
              ty1[k] = poly2->verts[tri2->v1]->coord[k];
              tz1[k] = poly2->verts[tri2->v2]->coord[k];
            } 

            poly2->verts[tri2->v2]->coord[0] += ta;
            poly2->verts[tri2->v2]->coord[1] += tb;
            poly2->verts[tri2->v2]->coord[2] += tc;
            //  Compute new normal and linecoeff
            newplanecoefficients(poly2, tri2);
            linecoeff[0] = tri1->c*tri2->b - tri1->b*tri2->c;
            linecoeff[1] = tri1->a*tri2->c - tri1->c*tri2->a;
            linecoeff[2] = tri1->b*tri2->a - tri1->a*tri2->b;
            dtemp = sqrt(linecoeff[0]*linecoeff[0] +
                         linecoeff[1]*linecoeff[1] +
                         linecoeff[2]*linecoeff[2]);
            linecoeff[0] /= dtemp;
            linecoeff[1] /= dtemp;
            linecoeff[2] /= dtemp;

            status = tri_tri_intersect(tri1,tri2);
            if ( status != CUBIT_SUCCESS )
            {
                return status;
            }

            //  Restore values.
            poly2->verts[tri2->v2]->coord[0] = tz1[0];
            poly2->verts[tri2->v2]->coord[1] = tz1[1];
            poly2->verts[tri2->v2]->coord[2] = tz1[2];
            tri2->a = ta; tri2->b = tb; tri2->c = tc; tri2->d = td;

            poly2->verts[tri2->v1]->coord[0] += ta;
            poly2->verts[tri2->v1]->coord[1] += tb;
            poly2->verts[tri2->v1]->coord[2] += tc;
            //  Compute new normal and linecoeff
            newplanecoefficients(poly2, tri2);
            linecoeff[0] = tri1->c*tri2->b - tri1->b*tri2->c;
            linecoeff[1] = tri1->a*tri2->c - tri1->c*tri2->a;
            linecoeff[2] = tri1->b*tri2->a - tri1->a*tri2->b;
            dtemp = sqrt(linecoeff[0]*linecoeff[0] +
                         linecoeff[1]*linecoeff[1] +
                         linecoeff[2]*linecoeff[2]);
            linecoeff[0] /= dtemp;
            linecoeff[1] /= dtemp;
            linecoeff[2] /= dtemp;

            status = tri_tri_intersect(tri1,tri2);
            if ( status != CUBIT_SUCCESS )
            {
                return status;
            }

            //  Restore values.
            poly2->verts[tri2->v1]->coord[0] = ty1[0];
            poly2->verts[tri2->v1]->coord[1] = ty1[1];
            poly2->verts[tri2->v1]->coord[2] = ty1[2];
            tri2->a = ta; tri2->b = tb; tri2->c = tc; tri2->d = td;

            poly2->verts[tri2->v0]->coord[0] += ta;
            poly2->verts[tri2->v0]->coord[1] += tb;
            poly2->verts[tri2->v0]->coord[2] += tc;
            //  Compute new normal and linecoeff
            newplanecoefficients(poly2, tri2);
            linecoeff[0] = tri1->c*tri2->b - tri1->b*tri2->c;
            linecoeff[1] = tri1->a*tri2->c - tri1->c*tri2->a;
            linecoeff[2] = tri1->b*tri2->a - tri1->a*tri2->b;
            dtemp = sqrt(linecoeff[0]*linecoeff[0] +
                         linecoeff[1]*linecoeff[1] +
                         linecoeff[2]*linecoeff[2]);
            linecoeff[0] /= dtemp;
            linecoeff[1] /= dtemp;
            linecoeff[2] /= dtemp;

            status = tri_tri_intersect(tri1,tri2);
            if ( status != CUBIT_SUCCESS )
            {
                return status;
            }

            //  Restore values.
            poly2->verts[tri2->v0]->coord[0] = tx1[0];
            poly2->verts[tri2->v0]->coord[1] = tx1[1];
            poly2->verts[tri2->v0]->coord[2] = tx1[2];
            tri2->a = ta; tri2->b = tb; tri2->c = tc; tri2->d = td;
            
            continue;
        }

        if ( do_imprint == true ) {
        //  Don't intersect nearly-coplanar triangles.
          if ( fabs(tri1->a*tri2->a + tri1->b*tri2->b +tri1->c*tri2->c) > 0.8 ) continue;
       }
        double dtemp;
        dtemp = sqrt(linecoeff[0]*linecoeff[0] +
                     linecoeff[1]*linecoeff[1] +
                     linecoeff[2]*linecoeff[2]);
        linecoeff[0] /= dtemp;
        linecoeff[1] /= dtemp;
        linecoeff[2] /= dtemp;

        status = tri_tri_intersect(tri1,tri2);
        if ( status != CUBIT_SUCCESS )
        {
            return status;
        }

     } // end of loop over poly2

  } // end of loop over poly1

    //delete entire array.
  delete [] boxlist;
  
  return status;
}

CubitStatus FBIntersect::tri_tri_intersect(FB_Triangle *tri1,
              FB_Triangle *tri2)
{
int ret1, ret2, i;
double xc10[3], xc11[3],xc12[3]; // coords for poly1 triangle
double xc20[3], xc21[3],xc22[3]; // coords for poly2 triangle
double d10, d11, d12, d20, d21, d22; // distance of vertex from plane of other triangle
double tt[4];
int edge_vert_type[4];
CubitStatus status;

  int mydebug = 0;
  status = CUBIT_SUCCESS;

  tt[0] = tt[1] = tt[2] = tt[3] = CUBIT_DBL_MAX;
//  Is tri1 entirely on one side of tri2?
  for ( i = 0; i < 3; i++ ) {
     xc10[i] = poly1->verts[tri1->v0]->coord[i];
     xc11[i] = poly1->verts[tri1->v1]->coord[i];
     xc12[i] = poly1->verts[tri1->v2]->coord[i];
   }

   //  distance of each tri1 vert to plane of tri2
   d10 = xc10[0]*tri2->a + xc10[1]*tri2->b + xc10[2]*tri2->c + tri2->d;
   d11 = xc11[0]*tri2->a + xc11[1]*tri2->b + xc11[2]*tri2->c + tri2->d;
   d12 = xc12[0]*tri2->a + xc12[1]*tri2->b + xc12[2]*tri2->c + tri2->d;
   if ( ( (d10 < -EPSILON2) && (d11 < -EPSILON2) && (d12 < -EPSILON2) ) ||
        ( (d10 > EPSILON2) && (d11 > EPSILON2) && (d12 > EPSILON2) ) ) 
      return CUBIT_SUCCESS;
//  Is tri2 entirely on one side of tri1?
  for ( i = 0; i < 3; i++ ) {
     xc20[i] = poly2->verts[tri2->v0]->coord[i];
     xc21[i] = poly2->verts[tri2->v1]->coord[i];
     xc22[i] = poly2->verts[tri2->v2]->coord[i];
   }

   //  distance of each tri2 vert to plane of tri1
   d20 = xc20[0]*tri1->a + xc20[1]*tri1->b + xc20[2]*tri1->c + tri1->d;
   d21 = xc21[0]*tri1->a + xc21[1]*tri1->b + xc21[2]*tri1->c + tri1->d;
   d22 = xc22[0]*tri1->a + xc22[1]*tri1->b + xc22[2]*tri1->c + tri1->d;
   if ( ( (d20 < -EPSILON2) && (d21 < -EPSILON2) && (d22 < -EPSILON2) ) ||
        ( (d20 > EPSILON2) && (d21 > EPSILON2) && (d22 > EPSILON2) ) )
     return CUBIT_SUCCESS;
//  Get a point on the line of intersection to serve as a reference point.
   double ta, tb, ts1, ts2, tdot11;
   ts1 = -tri1->d; ts2 = -tri2->d;

   tdot11 = tri1->a*tri2->a + tri1->b*tri2->b + tri1->c*tri2->c;
   ta = (ts2*tdot11 - ts1)/(tdot11*tdot11 - 1);
   tb = (ts1*tdot11 - ts2)/(tdot11*tdot11 - 1);
   linept[0] = ta*tri1->a + tb*tri2->a;
   linept[1] = ta*tri1->b + tb*tri2->b;
   linept[2] = ta*tri1->c + tb*tri2->c;
  
//  There are several cases for the distances.  
//  ret1 holds the number of intersections of the triangle with the
//  intersection line. 
//  Do tri1.
      ret1 = get_intersectionline_parameter_values(d10,d11,d12,
                                            xc10,xc11,xc12,
                   tt[0],tt[1],
                   edge_vert_type[0],edge_vert_type[1]);
//  Do tri2.
      ret2 = get_intersectionline_parameter_values(d20,d21,d22,
                                            xc20,xc21,xc22,
                   tt[2],tt[3],
                   edge_vert_type[2],edge_vert_type[3]);
//  If not two intersections for each triangle, no intersection edge exists.
    if ( (ret1 == 2) && (ret2 == 2) )
    {
        status = add_intersection_edges(tri1,tri2,tt,edge_vert_type);
        if ( status != CUBIT_SUCCESS )
        {
            return status;
        }
    }
    if(mydebug){
      double min_angle_1, max_angle_1;
      double min_angle_2, max_angle_2;
      poly1->min_max_angles_in_fb_triangle(tri1, min_angle_1, max_angle_1);
      poly2->min_max_angles_in_fb_triangle(tri2, min_angle_2, max_angle_2);
    
    
      if(min_angle_1 < 10 || min_angle_2 < 10){
        PRINT_INFO(" Tri 1 min angle = %f\n",min_angle_1);
        PRINT_INFO(" Tri 2 min angle = %f\n",min_angle_2);
      }
    }
    
    return status;
}

CubitStatus FBIntersect::add_intersection_edges(FB_Triangle *tri1,
                   FB_Triangle *tri2,
                   double *tt,
                   int *edge_vert_type
                   )
{
int v10, v11, v20, v21;
bool ifoundit, exists;
FB_Edge *edge;
double x1pt, y1pt, z1pt, x2pt, y2pt, z2pt;
int edge1_vert0_type, edge1_vert1_type, edge2_vert0_type, edge2_vert1_type;

  ifoundit = false;
  edge1_vert0_type = edge1_vert1_type = edge2_vert0_type = edge2_vert1_type = UNKNOWN;
 
//  Try to handle the epsilon cases by forcing tt[] values that are almost 
//  equal to be equal.

  if ( fabs(tt[0]-tt[2]) < EPSILON ) tt[2] = tt[0];
  if ( fabs(tt[0]-tt[3]) < EPSILON ) tt[3] = tt[0];
  if ( fabs(tt[1]-tt[2]) < EPSILON ) tt[2] = tt[1];
  if ( fabs(tt[1]-tt[3]) < EPSILON ) tt[3] = tt[1];
   
// cases 0 and 9, no overlap 
//  if ( (tt[1] < tt[2]) || (tt[3] < tt[0]) ) return CUBIT_SUCCESS;
  if ( tt[1] < tt[2] ) { return CUBIT_SUCCESS; }
  
  if ( tt[3] < tt[0] ) { return CUBIT_SUCCESS; }

//  These next four cases are by far the most common forms of overlap,
//  so check for them first.

//  The 1's (e.g.: 11111111) refer to the portion of the line of intersection
//  that is in tri1; the 2's for tri2.  The case numbers were arbirtarily assigned.

// case 4
//    11111111
//      2222        
  if ( (tt[2] > tt[0]) && (tt[3] < tt[1]) ) {
    ifoundit = true;    
    get_point_from_parameter(tt[2],&x1pt,&y1pt,&z1pt);
    get_point_from_parameter(tt[3],&x2pt,&y2pt,&z2pt);
    edge1_vert0_type = INTERIOR_VERT;
    edge1_vert1_type = INTERIOR_VERT;
    edge2_vert0_type = edge_vert_type[2];
    edge2_vert1_type = edge_vert_type[3];
    
// case 2
//    11111111
//        22222222     
  } else if ( (tt[2] < tt[1]) && (tt[2] > tt[0]) && (tt[3] > tt[1]) ) {
    ifoundit = true;
    get_point_from_parameter(tt[2],&x1pt,&y1pt,&z1pt);
    get_point_from_parameter(tt[1],&x2pt,&y2pt,&z2pt);
    edge1_vert0_type = INTERIOR_VERT;
    edge1_vert1_type = edge_vert_type[1];
    edge2_vert0_type = edge_vert_type[2];
    edge2_vert1_type = INTERIOR_VERT;
  
// case 7
//        11111111
//    22222222      
  } else if ( (tt[0] < tt[3]) && (tt[0] > tt[2]) && (tt[1] > tt[3]) ) {
    ifoundit = true;  
    get_point_from_parameter(tt[0],&x1pt,&y1pt,&z1pt);
    get_point_from_parameter(tt[3],&x2pt,&y2pt,&z2pt);
    edge1_vert0_type = edge_vert_type[0];
    edge1_vert1_type = INTERIOR_VERT;
    edge2_vert0_type = INTERIOR_VERT;
    edge2_vert1_type = edge_vert_type[3];
    
// case 12
//      1111
//    22222222      
  } else if ( (tt[0] > tt[2]) && (tt[1] < tt[3]) ) {
    ifoundit = true;  
    get_point_from_parameter(tt[0],&x1pt,&y1pt,&z1pt);
    get_point_from_parameter(tt[1],&x2pt,&y2pt,&z2pt);
    edge1_vert0_type = edge_vert_type[0];
    edge1_vert1_type = edge_vert_type[1];
    edge2_vert0_type = INTERIOR_VERT;
    edge2_vert1_type = INTERIOR_VERT;
    
// case 1 
//    11111111
//           22222222 
  } else if ( fabs(tt[1]-tt[2]) < EPSILON ) {
    ifoundit = true;
//    return CUBIT_SUCCESS;
    get_point_from_parameter(tt[1],&x1pt,&y1pt,&z1pt);
    get_point_from_parameter(tt[1],&x2pt,&y2pt,&z2pt);  
    edge1_vert0_type = edge_vert_type[1];
    edge1_vert1_type = edge_vert_type[1];
    edge2_vert0_type = edge_vert_type[2];
    edge2_vert1_type = edge_vert_type[2];    
   
// case 3
  } else if ( tt[0] < tt[2] ) {
// case 3
//    11111111
//       22222    
    if ( tt[1] == tt[3] ) {
      ifoundit = true;
      get_point_from_parameter(tt[2],&x1pt,&y1pt,&z1pt);
      get_point_from_parameter(tt[3],&x2pt,&y2pt,&z2pt);  
      edge1_vert0_type = INTERIOR_VERT;
      edge1_vert1_type = edge_vert_type[1];
      edge2_vert0_type = edge_vert_type[2];
      edge2_vert1_type = edge_vert_type[3];    
    }

// case 5, 6, or 11   
  } else if ( fabs(tt[0]-tt[2]) < EPSILON ) {

// case 5
//    11111111
//    22222222  
    if ( tt[1] == tt[3] ) { 
      ifoundit = true;
      get_point_from_parameter(tt[0],&x1pt,&y1pt,&z1pt);
      get_point_from_parameter(tt[1],&x2pt,&y2pt,&z2pt);
      edge1_vert0_type = edge_vert_type[0];
      edge1_vert1_type = edge_vert_type[1];
      edge2_vert0_type = edge_vert_type[2];
      edge2_vert1_type = edge_vert_type[3];    

// case 6
//    11111111
//    2222     
    } else if ( tt[1] > tt[3] ) { 
      ifoundit = true;
      get_point_from_parameter(tt[2],&x1pt,&y1pt,&z1pt);
      get_point_from_parameter(tt[3],&x2pt,&y2pt,&z2pt);
      edge1_vert0_type = edge_vert_type[0];
      edge1_vert1_type = INTERIOR_VERT;
      edge2_vert0_type = edge_vert_type[2];
      edge2_vert1_type = edge_vert_type[3];    

// case 11
//    1111
//    22222222      
    } else {
      ifoundit = true;
      get_point_from_parameter(tt[0],&x1pt,&y1pt,&z1pt);
      get_point_from_parameter(tt[1],&x2pt,&y2pt,&z2pt);
      edge1_vert0_type = edge_vert_type[0];
      edge1_vert1_type = edge_vert_type[1];
      edge2_vert0_type = edge_vert_type[2];
      edge2_vert1_type = INTERIOR_VERT;    
    }

// case 8 or 10      
  } else {

// case 8
//           11111111
//    22222222  
    if ( fabs(tt[0]-tt[3]) < EPSILON ) {
      ifoundit = true;
//      return CUBIT_SUCCESS;
      get_point_from_parameter(tt[0],&x1pt,&y1pt,&z1pt);
      get_point_from_parameter(tt[0],&x2pt,&y2pt,&z2pt);  
      edge1_vert0_type = edge_vert_type[0];
      edge1_vert1_type = edge_vert_type[0];
      edge2_vert0_type = edge_vert_type[3];
      edge2_vert1_type = edge_vert_type[3];    

// case 10
//        1111
//    22222222    
    } else if ( fabs(tt[1]-tt[3]) < EPSILON ) {
      ifoundit = true;
      get_point_from_parameter(tt[0],&x1pt,&y1pt,&z1pt);
      get_point_from_parameter(tt[1],&x2pt,&y2pt,&z2pt);       
      edge1_vert0_type = edge_vert_type[0];
      edge1_vert1_type = edge_vert_type[1];
      edge2_vert0_type = INTERIOR_VERT;
      edge2_vert1_type = edge_vert_type[3];    
    }
    
  }

  if ( ifoundit == false ) {
    PRINT_ERROR("unaccounted for case in add_intersection_edges: tt[] =  %le %le %le %le\n",
      tt[0],tt[1],tt[2],tt[3]);
    return CUBIT_FAILURE;
  } else {
    tri1->dudded = true;
    v10 = poly1->addavertex(x1pt,y1pt,z1pt);
    v11 = poly1->addavertex(x2pt,y2pt,z2pt);
    if ( v10 != v11 ) {
      exists = poly1->edge_exists_in_tri(*tri1,v10,v11);
      if ( exists == false ) {
        if ( edge1_vert0_type == INTERIOR_VERT )
          edge1_vert0_type = determine_edge_vert_type(edge_vert_type[0],
	                                              edge_vert_type[1]);
        if ( edge1_vert1_type == INTERIOR_VERT )
          edge1_vert1_type = determine_edge_vert_type(edge_vert_type[0],
	                                              edge_vert_type[1]);
        edge = new FB_Edge(v10,v11,edge1_vert0_type,edge1_vert1_type,true);
        tri1->edge_list.push_back(edge);
        if ( poly1->edge_exists(v10,v11) == false )
          poly1->intersection_edges.push_back(edge);
      }
    }  else {
      int edge_type = UNKNOWN;
      int vtype1 = UNKNOWN_VERT, vtype2 = UNKNOWN_VERT;
      int v_other1 = UNKNOWN_VERT, v_other2 = UNKNOWN_VERT;
//       edge_type = UNKNOWN;
      if ( edge1_vert0_type == EDGE_2 ) {
        v_other1 = tri1->v0;
        v_other2 = tri1->v2;
        edge_type = EDGE_2;
        vtype1 = VERTEX_0; 
        vtype2 = VERTEX_2;
      } else if ( edge1_vert0_type == EDGE_1 ) {
        v_other1 = tri1->v1;
        v_other2 = tri1->v2;
        edge_type = EDGE_1;
        vtype1 = VERTEX_1; 
        vtype2 = VERTEX_2;
      
      } else if ( edge1_vert0_type == EDGE_0 ) {
        v_other1 = tri1->v0;
        v_other2 = tri1->v1;
        edge_type = EDGE_0;
        vtype1 = VERTEX_0; 
        vtype2 = VERTEX_1;
      
      } 
      if ( edge_type != UNKNOWN ) {
        exists = poly1->edge_exists_in_tri(*tri1,v_other1,v10);
        if ( exists == false ) {    
          edge = new FB_Edge(v_other1,v10,vtype1,edge_type,false);
          tri1->edge_list.push_back(edge);
          if ( poly1->edge_exists(v_other1,v10) == false )
            poly1->intersection_edges.push_back(edge);
        }
        exists = poly1->edge_exists_in_tri(*tri1,v10,v_other2);
        if ( exists == false ) {    
          edge = new FB_Edge(v10,v_other2,edge_type,vtype2,false);
          tri1->edge_list.push_back(edge);
          if ( poly1->edge_exists(v10,v_other2) == false )
            poly1->intersection_edges.push_back(edge);
        } 
      }  
    }
    
    tri2->dudded = true;
    v20 = poly2->addavertex(x1pt,y1pt,z1pt);
    v21 = poly2->addavertex(x2pt,y2pt,z2pt);
    if ( v20 != v21) {
      exists = poly2->edge_exists_in_tri(*tri2,v20,v21);
      if ( exists == false ) {
        if ( edge2_vert0_type == INTERIOR_VERT )
          edge2_vert0_type = determine_edge_vert_type(edge_vert_type[2],
	                                              edge_vert_type[3]);
	if ( edge2_vert1_type == INTERIOR_VERT )
          edge2_vert1_type = determine_edge_vert_type(edge_vert_type[2],
	                                              edge_vert_type[3]);
        edge = new FB_Edge(v20,v21,edge2_vert0_type,edge2_vert1_type,true);
        tri2->edge_list.push_back(edge);
        if ( poly2->edge_exists(v20,v21) == false )
          poly2->intersection_edges.push_back(edge);
      }
    }  else {
      int edge_type = UNKNOWN;
      int vtype1 = UNKNOWN_VERT, vtype2 = UNKNOWN_VERT;
      int v_other1 = UNKNOWN_VERT, v_other2 = UNKNOWN_VERT;
//       edge_type = UNKNOWN;
      if ( edge2_vert0_type == EDGE_2 ) {
        v_other1 = tri2->v0;
        v_other2 = tri2->v2;
        edge_type = EDGE_2;
        vtype1 = VERTEX_0; 
        vtype2 = VERTEX_2;
      } else if ( edge2_vert0_type == EDGE_1 ) {
        v_other1 = tri2->v1;
        v_other2 = tri2->v2;
        edge_type = EDGE_1;
        vtype1 = VERTEX_1; 
        vtype2 = VERTEX_2;
      
      } else if ( edge2_vert0_type == EDGE_0 ) {
        v_other1 = tri2->v0;
        v_other2 = tri2->v1;
        edge_type = EDGE_0;
        vtype1 = VERTEX_0; 
        vtype2 = VERTEX_1;
      
      } 
      if ( edge_type != UNKNOWN ) {
        exists = poly2->edge_exists_in_tri(*tri2,v_other1,v20);
        if ( exists == false ) {    
          edge = new FB_Edge(v_other1,v20,vtype1,edge_type,false);
          tri2->edge_list.push_back(edge);
          if ( poly2->edge_exists(v_other1,v20) == false )
            poly2->intersection_edges.push_back(edge);
        }
        exists = poly2->edge_exists_in_tri(*tri2,v20,v_other2);
        if ( exists == false ) {    
          edge = new FB_Edge(v20,v_other2,edge_type,vtype2,false);
          tri2->edge_list.push_back(edge);
          if ( poly2->edge_exists(v20,v_other2) == false )
            poly2->intersection_edges.push_back(edge);
        }        
      }
    }
  }

  return CUBIT_SUCCESS;
  
}

void FBIntersect::get_point_from_parameter(double parameter,
                 double *x, double *y, double *z)
{
  *x = linept[0] + linecoeff[0]*parameter; 
  *y = linept[1] + linecoeff[1]*parameter; 
  *z = linept[2] + linecoeff[2]*parameter; 

}

double FBIntersect::get_distance_parameter(double *xc0,
                 double *xc1,
                 double d0, double d1)
{
double v1dot, v2dot;

//  dot the coordinates onto the line.
//  Assume that it has been checked already that d0 != d1.

  v1dot = linecoeff[0]*(xc0[0] - linept[0]) +
          linecoeff[1]*(xc0[1] - linept[1]) +
          linecoeff[2]*(xc0[2] - linept[2]);
  v2dot = linecoeff[0]*(xc1[0] - linept[0]) +
          linecoeff[1]*(xc1[1] - linept[1]) +
          linecoeff[2]*(xc1[2] - linept[2]);

  return v1dot + (v2dot - v1dot)*d0/(d0-d1);
}


double FBIntersect::get_distance_parameter_single(double *xc)
{
double coeff;
int coord;

  if ( (fabs(linecoeff[0]) >= fabs(linecoeff[1])) && 
       (fabs(linecoeff[0]) >= fabs(linecoeff[2])) ){
    coord = 0;
  } else if ( fabs(linecoeff[1]) >= fabs(linecoeff[2]) ) {
    coord = 1;
  }  else {
    coord = 2;
  }
  coeff = linecoeff[coord];
  return (xc[coord] - linept[coord])/coeff; 
}

int FBIntersect::get_intersectionline_parameter_values(
                    double d0, 
                    double d1, 
                    double d2,
                    double *pt0,
                    double *pt1,
                    double *pt2,
                    double& t0,
                    double& t1,
                    int& vert_type_0,
                    int& vert_type_1)
{
int ret = 0;
//  ret holds the nnumber of intersections of the line with the triangle.

  if ( fabs(d0) < EPSILON ) {
    if ( fabs(d1) < EPSILON ) {
      if ( fabs(d2) < EPSILON ) {
      // coplanar 
        ret = 0;
      } else {
      //  d0 = d1 = 0
        t0 = get_distance_parameter_single(pt0);
        t1 = get_distance_parameter_single(pt1);
   vert_type_0 = VERTEX_0;
   vert_type_1 = VERTEX_1;
        ret = 2;
      }
    } else if ( fabs(d2) < EPSILON ) {
    //  d0 = d2 = 0
      t0 = get_distance_parameter_single(pt0);
      t1 = get_distance_parameter_single(pt2);      
      vert_type_0 = VERTEX_0;
      vert_type_1 = VERTEX_2;     
      ret = 2;
    } else if ( d1*d2 < 0.0 ) {
    //  d0 = 0 and edge 12 crosses
      t0 = get_distance_parameter_single(pt0);
      t1 = get_distance_parameter(pt1,pt2,d1,d2);
      vert_type_0 = VERTEX_0;
      vert_type_1 = EDGE_1;     
      ret = 2;
    } else {
    //  d0 = 0 and no edges cross
      ret = 0;
    }
  } else if ( fabs(d1) < EPSILON ) {
    if ( fabs(d2) < EPSILON ) {
    //  d1 = d2 = 0
      t0 = get_distance_parameter_single(pt1);
      t1 = get_distance_parameter_single(pt2); 
      vert_type_0 = VERTEX_1;
      vert_type_1 = VERTEX_2;     
      ret = 2;   
    } else if ( d0*d2 < 0.0 ) {
    //  d1 = 0 and edge 20 crosses
      t0 = get_distance_parameter_single(pt1);    
      t1 = get_distance_parameter(pt0,pt2,d0,d2);
      vert_type_0 = VERTEX_1;
      vert_type_1 = EDGE_2;     
      ret = 2;
    } else {
    //  d1 = 0 and no edges cross
      ret = 0;
    }
  } else if ( fabs(d2) < EPSILON ) {
    if ( d0*d1 < 0.0 ) {
    //  d2 = 0 and edge 01 crosses
      t0 = get_distance_parameter_single(pt2);    
      t1 = get_distance_parameter(pt0,pt1,d0,d1);
      vert_type_0 = VERTEX_2;
      vert_type_1 = EDGE_0;     
      ret = 2;    
    } else {
    //  d2 = 0 and no edges cross
      ret = 0;
    }
  } else if ( d0*d1 < 0.0 ) {
    if ( d0*d2 < 0.0 ) {
    //  edges 01 and 02 cross
      t0 = get_distance_parameter(pt0,pt1,d0,d1);
      t1 = get_distance_parameter(pt0,pt2,d0,d2);
      vert_type_0 = EDGE_0;
      vert_type_1 = EDGE_2;     
      ret = 2;
    } else {
    //  edges 01 and 12 cross
      t0 = get_distance_parameter(pt0,pt1,d0,d1);
      t1 = get_distance_parameter(pt1,pt2,d1,d2);
      vert_type_0 = EDGE_0;
      vert_type_1 = EDGE_1;     
      ret = 2;
    }
  } else {
  //  edges 02 and 12 cross
      t0 = get_distance_parameter(pt0,pt2,d0,d2);
      t1 = get_distance_parameter(pt1,pt2,d1,d2);
      vert_type_0 = EDGE_2;
      vert_type_1 = EDGE_1;     
      ret = 2;
  }
  
  if ( ret == 2 ) {  //  Sort them if necessary.
    if ( t0 > t1 ) {
      int itemp;
      double dtemp;
      dtemp = t0;
      t0 = t1;
      t1 = dtemp;
      itemp = vert_type_0;
      vert_type_0 = vert_type_1;
      vert_type_1 = itemp;
    }
  
  } 
  
  return ret;
}

void FBIntersect::set_classify_flag(bool value)
{
  do_classify = value;
}

CubitStatus FBIntersect::get_persistent_entity_info(bool *surfs_in,
                                                    bool *curves_in,
                                                    bool *surfs_out,
                                                    bool *curves_out,
                                                    const CubitFacetboolOp op,
                                                    const int whichparent
                                       )
{
std::vector<int> *group, *groupcharacterization;
std::vector<int>::iterator it, ig;
int booltest1, booltest2, booltest1a, booltest2a, booltest_in, booltest_out;
FBPolyhedron *poly;

  if ( nothing_intersected == true ) return CUBIT_SUCCESS;
  if ( !poly1 || !poly2 ) {
    PRINT_ERROR("Error:  Objects for Booleans must first be created.\n");
    return CUBIT_FAILURE;
  }
  if ( !classify1 || !classify2 ) {
    PRINT_ERROR("Error:  Objects for Booleans must first be classified.\n");
    return CUBIT_FAILURE;
  } 
  if ( (whichparent < 1) || (whichparent > 2) ) {
    PRINT_ERROR("Error:  Requested nonexistent object.\n");
    return CUBIT_FAILURE;  
  }
  if ( op == CUBIT_FB_UNION ) {
    booltest1 = FB_ORIENTATION_INSIDE + FB_ORIENTATION_SAME;
    booltest2 = FB_ORIENTATION_INSIDE;
    booltest1a = FB_ORIENTATION_OUTSIDE + FB_ORIENTATION_SAME;
    booltest2a = FB_ORIENTATION_OUTSIDE;
  } else if ( op == CUBIT_FB_INTERSECTION ) {
    booltest1 = FB_ORIENTATION_INSIDE + FB_ORIENTATION_SAME;
    booltest2 = FB_ORIENTATION_INSIDE; 
    booltest1a = FB_ORIENTATION_OUTSIDE+ FB_ORIENTATION_OPPOSITE;
    booltest2a = FB_ORIENTATION_INSIDE;         
  } else if ( op == CUBIT_FB_SUBTRACTION ) {
    booltest1 = FB_ORIENTATION_OUTSIDE+ FB_ORIENTATION_OPPOSITE;
    booltest2 = FB_ORIENTATION_INSIDE;    
    booltest1a = FB_ORIENTATION_INSIDE + FB_ORIENTATION_SAME;
    booltest2a = FB_ORIENTATION_INSIDE; 
  } else {
    PRINT_ERROR("Error:  Unrecognized Boolean operation.\n");
    return CUBIT_FAILURE;
  } 
 
  if ( whichparent == 1 ) {
    booltest_in = booltest1;
    booltest_out = booltest1a;
    classify1->get_group(&group, &groupcharacterization);
    poly = poly1;
  } else {
    booltest_in = booltest2;
    booltest_out = booltest2a;
    classify2->get_group(&group, &groupcharacterization);
    poly = poly2;
  }
  unsigned int i;
  int isurf, icurve0, icurve1, icurve2;

  it = group->begin();
  ig = groupcharacterization->begin(); 

  i = 0;
  while ( it != group->end() ) {
    if ( ( *(ig + *it) & booltest_in ) != 0 ) {
//      vtx = poly1->tris[i]->v0;
      isurf = poly->tris[i]->cubitsurfaceindex;
      if ( isurf > 0 ) 
        surfs_in[isurf] = true;
      icurve0 = poly->tris[i]->cubitedge0index;
      if ( icurve0 > 0 ) 
        curves_in[icurve0] = true;
      icurve1 = poly->tris[i]->cubitedge1index;
      if ( icurve1 > 0 ) 
        curves_in[icurve1] = true;
       icurve2 = poly->tris[i]->cubitedge2index;
      if ( icurve2 > 0 ) 
        curves_in[icurve2] = true;
   }
    it++;
    i++;
  }  
  
  it = group->begin();
  ig = groupcharacterization->begin(); 

  i = 0;
  while ( it != group->end() ) {
    if ( ( *(ig + *it) & booltest_out ) != 0 ) {
      isurf = poly->tris[i]->cubitsurfaceindex;
      if ( isurf > 0 ) 
        surfs_out[isurf] = true;
      icurve0 = poly->tris[i]->cubitedge0index;
      if ( icurve0 > 0 ) 
        curves_out[icurve0] = true;
      icurve1 = poly->tris[i]->cubitedge1index;
      if ( icurve1 > 0 ) 
        curves_out[icurve1] = true;
       icurve2 = poly->tris[i]->cubitedge2index;
      if ( icurve2 > 0 ) 
        curves_out[icurve2] = true;
   }
    it++;
    i++;
  }  

  return CUBIT_SUCCESS;

}
                                         
CubitStatus FBIntersect::update_surfs_and_curves(std::vector<double>& out_coords,
                                    std::vector<int>& out_connections,
                                    std::vector<int> *out_surf_index,
                                    std::vector<int> *out_curve_index, 
                                    const int whichone)
{
unsigned int i;
FBPolyhedron *poly;
 
  if ( whichone == 1 ) poly = poly1;
  else poly = poly2;
  for ( i = 0; i < poly->verts.size(); i++ ) {
    out_coords.push_back(poly->verts[i]->coord[0]);
    out_coords.push_back(poly->verts[i]->coord[1]);
    out_coords.push_back(poly->verts[i]->coord[2]);
  } 
  for ( i = 0; i < poly->tris.size(); i++ ) {
    if ( poly->tris[i]->dudded == true ) continue;
    out_connections.push_back(poly->tris[i]->v0);
    out_connections.push_back(poly->tris[i]->v1);
    out_connections.push_back(poly->tris[i]->v2);
    out_surf_index->push_back(poly->tris[i]->cubitsurfaceindex);
    out_curve_index->push_back(poly->tris[i]->cubitedge0index);
    out_curve_index->push_back(poly->tris[i]->cubitedge1index);
    out_curve_index->push_back(poly->tris[i]->cubitedge2index); 
  }
  return CUBIT_SUCCESS;
}

CubitStatus FBIntersect::gather_by_boolean(std::vector<double>& out_coords,
                                    std::vector<int>& out_connections,
                                    std::vector<int> *out_surf_index,
                                    std::vector<int> *out_curve_index,
                                    std::vector<bool> *is_body_1,
                                    const CubitFacetboolOp op)
{
    CubitStatus status;
    std::vector<int> *group1,
                     *group2,
                     *groupcharacterization1, 
                     *groupcharacterization2;
    std::vector<int>::iterator it, ig;
    int booltest1, booltest2;

    if ( nothing_intersected == true && op !=CUBIT_FB_UNION)
    {
        return CUBIT_SUCCESS;
    }
    if ( !poly1 || !poly2 )
    {
        PRINT_ERROR("Error:  Objects for Booleans must first be created.\n");
        return CUBIT_FAILURE;
    }
    if ( !classify1 || !classify2 )
    {
        PRINT_ERROR("Error:  Objects for Booleans must first be classified.\n");
        return CUBIT_FAILURE;
    } 
  
    if ( op == CUBIT_FB_UNION )
    {
        booltest1 = FB_ORIENTATION_OUTSIDE + FB_ORIENTATION_SAME;
        booltest2 = FB_ORIENTATION_OUTSIDE;
    }
    else if ( op == CUBIT_FB_INTERSECTION )
    {
        booltest1 = FB_ORIENTATION_INSIDE + FB_ORIENTATION_SAME;
        booltest2 = FB_ORIENTATION_INSIDE;  
    }
    else if ( op == CUBIT_FB_SUBTRACTION )
    {
        booltest1 = FB_ORIENTATION_OUTSIDE+ FB_ORIENTATION_OPPOSITE;
        booltest2 = FB_ORIENTATION_INSIDE;    
    }
    else
    {
        PRINT_ERROR("Error:  Unrecognized Boolean operation.\n");
        return CUBIT_FAILURE;
    } 

    classify1->get_group(&group1, &groupcharacterization1);
    classify2->get_group(&group2, &groupcharacterization2);

    IntegerHash *hashobj = new IntegerHash(NUMHASHBINS,20);

    it = group1->begin();
    ig = groupcharacterization1->begin(); 
    unsigned int i;
    int vtx, vertnum1, vertnum2, vertnum3, verts_sofar;

    i = 0;
    verts_sofar = 0;
    while ( it != group1->end() )
    {
        if ( ( *(ig + *it) & booltest1 ) != 0 )
        {
            vtx = poly1->tris[i]->v0;
            vertnum1 = get_vertex(poly1, vtx, hashobj, out_coords, verts_sofar);

            vtx = poly1->tris[i]->v1;
            vertnum2 = get_vertex(poly1, vtx, hashobj, out_coords, verts_sofar);

            vtx = poly1->tris[i]->v2;
            vertnum3 = get_vertex(poly1, vtx, hashobj, out_coords, verts_sofar);

            if ( store_connectivity( out_connections,
                                     vertnum1,
                                     vertnum2,
                                     vertnum3 ) != CUBIT_SUCCESS )
            {
                return CUBIT_FAILURE;
            }
            if ( out_surf_index )
            {
                out_surf_index->push_back(poly1->tris[i]->cubitsurfaceindex);
            }
            if ( out_curve_index )
            {
                out_curve_index->push_back(poly1->tris[i]->cubitedge0index);
                out_curve_index->push_back(poly1->tris[i]->cubitedge1index);
                out_curve_index->push_back(poly1->tris[i]->cubitedge2index);
            }
            if ( is_body_1 ) is_body_1->push_back(true);     
        }
        it++;
        i++;
    }
    it = group2->begin();
    ig = groupcharacterization2->begin();
    i = 0;
    while ( it != group2->end() )
    {
        if ( ( *(ig + *it) & booltest2 ) != 0 )
        {
            vtx = poly2->tris[i]->v0;
            vertnum1 = get_vertex(poly2, vtx, hashobj, out_coords, verts_sofar);

            vtx = poly2->tris[i]->v1;
            vertnum2 = get_vertex(poly2, vtx, hashobj, out_coords, verts_sofar);

            vtx = poly2->tris[i]->v2;
            vertnum3 = get_vertex(poly2, vtx, hashobj, out_coords, verts_sofar);

            if ( op == CUBIT_FB_SUBTRACTION ) //  reverse the winding      
            {
                if ( store_connectivity( out_connections,
                                         vertnum1,
                                         vertnum3,
                                         vertnum2 ) != CUBIT_SUCCESS )
                {
                    return CUBIT_FAILURE;
                }
                if ( out_curve_index )
                {
                    out_curve_index->push_back(poly2->tris[i]->cubitedge2index);
                    out_curve_index->push_back(poly2->tris[i]->cubitedge1index);
                    out_curve_index->push_back(poly2->tris[i]->cubitedge0index);
                }
            }
            else
            {
                if ( store_connectivity( out_connections,
                                         vertnum1,
                                         vertnum2,
                                         vertnum3 ) != CUBIT_SUCCESS )
                {
                    return CUBIT_FAILURE;
                }
                if ( out_curve_index )
                {
                    out_curve_index->push_back(poly2->tris[i]->cubitedge0index);
                    out_curve_index->push_back(poly2->tris[i]->cubitedge1index);
                    out_curve_index->push_back(poly2->tris[i]->cubitedge2index); 
                }
            }
            if ( out_surf_index ) out_surf_index->push_back(poly2->tris[i]->cubitsurfaceindex);
            if ( is_body_1 ) is_body_1->push_back(false);
        }
        it++;
        i++;
    } 
//  GfxDebug::clear();
//  poly1->debug_draw_boundary_edges(CUBIT_MAGENTA);
//  GfxDebug::mouse_xforms();
//  GfxDebug::clear();
//  poly2->debug_draw_boundary_edges(CUBIT_GREEN);
    delete hashobj;   
 
    status = CUBIT_SUCCESS;  
    return status;
}

CubitStatus FBIntersect::store_connectivity
(
     std::vector<int>& out_connections,
     int vertnum1,
     int vertnum2,
     int vertnum3
)
{
    if ( vertnum1 == vertnum2 ||
         vertnum2 == vertnum3 ||
         vertnum1 == vertnum3 )
    {
        PRINT_ERROR( "Cannot continue without generating a degenerate facet.\n" );
        return CUBIT_FAILURE;
    }
    out_connections.push_back(vertnum1);   
    out_connections.push_back(vertnum2);   
    out_connections.push_back(vertnum3); 
    return CUBIT_SUCCESS;
}

int FBIntersect::get_vertex(FBPolyhedron *poly, int vtx, 
                            IntegerHash *hashobj,
                            std::vector<double>& out_coords,
                            int &num_sofar)
{
double xx, yy, zz, xval, yval, zval;
int i, hashvalue, *hasharrayptr, hasharraysize, hptr, ifoundit;

  xx = poly->verts[vtx]->coord[0];
  yy = poly->verts[vtx]->coord[1];
  zz = poly->verts[vtx]->coord[2];
  hashvalue = makeahashvaluefrom_coord(xx,yy,zz);
  hasharrayptr = hashobj->getHashBin(hashvalue,&hasharraysize);
  ifoundit = -1;
  for ( i = 0; i < hasharraysize; i++ ) {
    hptr = hasharrayptr[i];
    xval = out_coords[3*hptr];
    yval = out_coords[3*hptr+1];
    zval = out_coords[3*hptr+2];
    if ( ( fabs(xval-xx) < EPSILON ) && 
    ( fabs(yval-yy) < EPSILON ) &&
    ( fabs(zval-zz) < EPSILON ) ) {
      ifoundit = hasharrayptr[i];
      break;
    }
  }
  if ( ifoundit == -1 ) {
    ifoundit = num_sofar;
    hashobj->addtoHashList(hashvalue,num_sofar);
    out_coords.push_back(xx);
    out_coords.push_back(yy);
    out_coords.push_back(zz);    
    num_sofar++;
  }
    
  return ifoundit;
}
                                   
int FBIntersect::makeahashvaluefrom_coord(double x, double y, double z)
{
double mantissasum;

      if ( fabs(x) < 1.e-3 ) x = 0.0;
      if ( fabs(y) < 1.e-3 ) y = 0.0;
      if ( fabs(z) < 1.e-3 ) z = 0.0;
      mantissasum = (int)(10000.0*fabs(x) + 0.5) + 
                    (int)(10000.0*fabs(y) + 0.5) + 
          (int)(10000.0*fabs(z) + 0.5);
      
      return (int)(mantissasum) % NUMHASHBINS;
}
                                    
                                      
void FBIntersect::set_body1_planar()
{
  body1_is_plane = true;
}
                              
void FBIntersect::set_body2_planar()
{
  body2_is_plane = true;
}

void FBIntersect::set_imprint()
{
  do_imprint = true;
}

void FBIntersect::newplanecoefficients(FBPolyhedron *poly, FB_Triangle *tri)
{
FB_Coord *mycoord;
double x1, x2, x3, y1, y2, y3, z1, z2, z3, e1x, e1y, e1z, e2x, e2y, e2z;
double a, b, c, d, dtemp;

     mycoord = poly->verts[tri->v0];
     x1 = mycoord->coord[0];
     y1 = mycoord->coord[1];
     z1 = mycoord->coord[2];
     mycoord = poly->verts[tri->v1];
     x2 = mycoord->coord[0];
     y2 = mycoord->coord[1];
     z2 = mycoord->coord[2];
     mycoord = poly->verts[tri->v2];
     x3 = mycoord->coord[0];
     y3 = mycoord->coord[1];
     z3 = mycoord->coord[2];
     e1x = x1 - x2; e1y = y1 - y2; e1z = z1 - z2;
     e2x = x3 - x2; e2y = y3 - y2; e2z = z3 - z2;
     a = e1z*e2y - e2z*e1y;
     b = e1x*e2z - e2x*e1z;
     c = e1y*e2x - e2y*e1x;
     dtemp = sqrt(a*a + b*b + c*c);
     if ( dtemp > EPSILON2 ) {
       a /= dtemp;
       b /= dtemp;
       c /= dtemp;
     } else {
       PRINT_WARNING("small-area triangle\n");
     }
     d = -(a*x1 + b*y1 + c*z1);
     tri->a = a; tri->b = b; tri->c = c; tri->d = d;
          
}
