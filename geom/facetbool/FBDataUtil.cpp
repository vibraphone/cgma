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

#include "FBDataUtil.hpp"
#include "GeometryDefines.h"
#include "DLIList.hpp"

//===========================================================================
//Function Name: intersect_plane_with_boundingbox
//Description:   Find and return the number of intersections of a plane with
//               the edges of a bounding box.  There will be 6 or fewer.
//               These will be returned in the intersection_points array.
//Author: mlstate - rewritten from John Fowler's original version.
//Date: 3/2005
//===========================================================================
void FBDataUtil::intersect_plane_with_boundingbox(CubitBox &bbox, 
                                     const CubitVector &v0, 
                                     const CubitVector &v1,
                                     const CubitVector &v2,
                                     DLIList<CubitVector> &intersection_points )
{
CubitVector v3 = v0 - v1;
CubitVector v4 = v2 - v1;
CubitVector threepointnormal = v3*v4;
threepointnormal.normalize();
//CubitVector threepointcenter = (v0 + v1 + v2)/3.;
//  Need to get the size and location of the facetted cutting plane.  Do this by 
//  intersecting the 3-point plane with the bounding box edges.  Will get a max of
//  six intersections.  The size of the box for these intersections will determine
//  the cutting plane size; the centroid of the intersections will be the location
//  of the center of the cutting plane.

double a, b, c, d;  //  plane coefficients
  a = threepointnormal.x();
  b = threepointnormal.y();
  c = threepointnormal.z();
  d = -(a*v0.x() + b*v0.y() + c*v0.z());
double xbmin, ybmin, zbmin, xbmax, ybmax, zbmax;
double testpoint;

    CubitVector boxmin = bbox.minimum();
    xbmin = boxmin.x(); ybmin = boxmin.y(); zbmin = boxmin.z();
    CubitVector boxmax = bbox.maximum();
    xbmax = boxmax.x(); ybmax = boxmax.y(); zbmax = boxmax.z();

    if ( fabs(a) < GEOMETRY_RESABS &&
         fabs(b) < GEOMETRY_RESABS )
    {
        // Normal points in Z dir only.
        //  a = 0; b = 0
        if ( ((zbmin + d/c) < -GEOMETRY_RESABS) && ((zbmax + d/c) > GEOMETRY_RESABS) )
        {
            add_unique_point( intersection_points, CubitVector(xbmin, ybmin, -d/c) );
            add_unique_point( intersection_points, CubitVector(xbmax, ybmin, -d/c) );
            add_unique_point( intersection_points, CubitVector(xbmax, ybmax, -d/c) );
            add_unique_point( intersection_points, CubitVector(xbmin, ybmax, -d/c) );
        }
    }
    else if ( fabs(a) < GEOMETRY_RESABS &&
              fabs(c) < GEOMETRY_RESABS )
    {
        // Normal points in Y dir only.
        if ( ((ybmin + d/b) < -GEOMETRY_RESABS) && ((ybmax + d/b) > GEOMETRY_RESABS) )
        {
            add_unique_point( intersection_points, CubitVector(xbmin, -d/b, zbmin ) );
            add_unique_point( intersection_points, CubitVector(xbmax, -d/b, zbmin ) );
            add_unique_point( intersection_points, CubitVector(xbmax, -d/b, zbmax ) );
            add_unique_point( intersection_points, CubitVector(xbmin, -d/b, zbmax ) );
        }
    }
    else if ( fabs(b) < GEOMETRY_RESABS &&
              fabs(c) < GEOMETRY_RESABS )
    {
        // Normal points in X dir only.
        if ( ((xbmin + d/a) < -GEOMETRY_RESABS) && ((xbmax + d/a) > GEOMETRY_RESABS) )
        {
            add_unique_point( intersection_points, CubitVector(-d/a, ybmin, zbmin ) );
            add_unique_point( intersection_points, CubitVector(-d/a, ybmax, zbmin ) );
            add_unique_point( intersection_points, CubitVector(-d/a, ybmax, zbmax ) );
            add_unique_point( intersection_points, CubitVector(-d/a, ybmin, zbmax ) );
        }
    }
    else if ( fabs(a) < GEOMETRY_RESABS )
    {
        // Normal is in YZ plane.
        testpoint = -(c*zbmin + d)/b;
        if ( (ybmin < (testpoint + GEOMETRY_RESABS)) &&
             (ybmax > (testpoint - GEOMETRY_RESABS)) )
        {
            add_unique_point( intersection_points, CubitVector(xbmin, testpoint, zbmin) );
            add_unique_point( intersection_points, CubitVector(xbmax, testpoint, zbmin) );
        }
        testpoint = -(c*zbmax + d)/b;
        if ( (ybmin < (testpoint + GEOMETRY_RESABS)) &&
             (ybmax > (testpoint - GEOMETRY_RESABS)) )
        {
            add_unique_point( intersection_points, CubitVector(xbmin, testpoint, zbmax) );
            add_unique_point( intersection_points, CubitVector(xbmax, testpoint, zbmax) );
        }
        testpoint = -(b*ybmin + d)/c;
        if ( (zbmin < (testpoint + GEOMETRY_RESABS)) &&
             (zbmax > (testpoint - GEOMETRY_RESABS)) )
        {
            add_unique_point( intersection_points, CubitVector(xbmin, ybmin, testpoint) );
            add_unique_point( intersection_points, CubitVector(xbmax, ybmin, testpoint) );
        }
        testpoint = -(b*ybmax + d)/c;
        if ( (zbmin < (testpoint + GEOMETRY_RESABS)) &&
             (zbmax > (testpoint - GEOMETRY_RESABS)) )
        {
            add_unique_point( intersection_points, CubitVector(xbmin, ybmax, testpoint) );
            add_unique_point( intersection_points, CubitVector(xbmax, ybmax, testpoint) );
        }
    }
    else if ( fabs(b) < GEOMETRY_RESABS )
    {
        // Normal is in XZ plane
        testpoint = -(c*zbmin + d)/a;
        if ( (xbmin < (testpoint + GEOMETRY_RESABS)) &&
             (xbmax > (testpoint - GEOMETRY_RESABS)) )
        {
            add_unique_point( intersection_points, CubitVector(testpoint, ybmin, zbmin) );
            add_unique_point( intersection_points, CubitVector(testpoint, ybmax, zbmin) );
        }
        testpoint = -(c*zbmax + d)/a;
        if ( (xbmin < (testpoint + GEOMETRY_RESABS)) &&
             (xbmax > (testpoint - GEOMETRY_RESABS)) )
        {
            add_unique_point( intersection_points, CubitVector(testpoint, ybmin, zbmax) );
            add_unique_point( intersection_points, CubitVector(testpoint, ybmax, zbmax) );
        }
        testpoint = -(a*xbmin + d)/c;
        if ( (zbmin < (testpoint + GEOMETRY_RESABS)) &&
             (zbmax > (testpoint - GEOMETRY_RESABS)) )
        {
            add_unique_point( intersection_points, CubitVector(xbmin, ybmin, testpoint) );
            add_unique_point( intersection_points, CubitVector(xbmin, ybmax, testpoint) );
        }
        testpoint = -(a*xbmax + d)/c;
        if ( (zbmin < (testpoint + GEOMETRY_RESABS)) &&
             (zbmax > (testpoint - GEOMETRY_RESABS)) )
        {
            add_unique_point( intersection_points, CubitVector(xbmax, ybmin, testpoint) );
            add_unique_point( intersection_points, CubitVector(xbmax, ybmax, testpoint) );
        }
    }
    else if ( fabs(c) < GEOMETRY_RESABS )
    {
        // Normal is in XY plane
        testpoint = -(a*xbmin + d)/b;
        if ( (ybmin < (testpoint + GEOMETRY_RESABS)) &&
             (ybmax > (testpoint - GEOMETRY_RESABS)) )
        {
            add_unique_point( intersection_points, CubitVector(xbmin, testpoint, zbmin) );
            add_unique_point( intersection_points, CubitVector(xbmin, testpoint, zbmax) );
        }
        testpoint = -(a*xbmax + d)/b;
        if ( (ybmin < (testpoint + GEOMETRY_RESABS)) &&
             (ybmax > (testpoint - GEOMETRY_RESABS)) )
        {
            add_unique_point( intersection_points, CubitVector(xbmax, testpoint, zbmin) );
            add_unique_point( intersection_points, CubitVector(xbmax, testpoint, zbmax) );
        }
        testpoint = -(b*ybmin + d)/a;
        if ( (xbmin < (testpoint + GEOMETRY_RESABS)) &&
             (xbmax > (testpoint - GEOMETRY_RESABS)) )
        {
            add_unique_point( intersection_points, CubitVector(testpoint, ybmin, zbmin) );
            add_unique_point( intersection_points, CubitVector(testpoint, ybmin, zbmax) );
        }
        testpoint = -(b*ybmax + d)/a;
        if ( (xbmin < (testpoint + GEOMETRY_RESABS)) &&
             (xbmax > (testpoint - GEOMETRY_RESABS)) )
        {
            add_unique_point( intersection_points, CubitVector(testpoint, ybmax, zbmin) );
            add_unique_point( intersection_points, CubitVector(testpoint, ybmax, zbmax) );
        }
    }
    else
    {
        // The general case
        // a != 0; b != 0; c != 0
        testpoint = -(b*ybmin + c*zbmin + d)/a;
        if ( (testpoint > (xbmin-GEOMETRY_RESABS)) && (testpoint < (xbmax+GEOMETRY_RESABS)) )
        {
            add_unique_point( intersection_points, CubitVector(testpoint, ybmin, zbmin) );
        }
        if ( intersection_points.size() == 6 ) return;
        testpoint = -(b*ybmax + c*zbmin + d)/a;
        if ( (testpoint > (xbmin-GEOMETRY_RESABS)) && (testpoint < (xbmax+GEOMETRY_RESABS)) )
        {
            add_unique_point( intersection_points, CubitVector(testpoint, ybmax, zbmin) );
        }
        if ( intersection_points.size() == 6 ) return;
        testpoint = -(b*ybmax + c*zbmax + d)/a;
        if ( (testpoint > (xbmin-GEOMETRY_RESABS)) && (testpoint < (xbmax+GEOMETRY_RESABS)) )
        {
            add_unique_point( intersection_points, CubitVector(testpoint, ybmax, zbmax) );
        }
        if ( intersection_points.size() == 6 ) return;
        testpoint = -(b*ybmin + c*zbmax + d)/a;
        if ( (testpoint > (xbmin-GEOMETRY_RESABS)) && (testpoint < (xbmax+GEOMETRY_RESABS)) )
        {
            add_unique_point( intersection_points, CubitVector(testpoint, ybmin, zbmax) );
        }
        if ( intersection_points.size() == 6 ) return;
        testpoint = -(a*xbmin + c*zbmin + d)/b;
        if ( (testpoint > (ybmin-GEOMETRY_RESABS)) && (testpoint < (ybmax+GEOMETRY_RESABS)) )
        {
            add_unique_point( intersection_points, CubitVector(xbmin, testpoint, zbmin) );
        }
        if ( intersection_points.size() == 6 ) return;
        testpoint = -(a*xbmax + c*zbmin + d)/b;
        if ( (testpoint > (ybmin-GEOMETRY_RESABS)) && (testpoint < (ybmax+GEOMETRY_RESABS)) )
        {
            add_unique_point( intersection_points, CubitVector(xbmax, testpoint, zbmin) );
        }
        if ( intersection_points.size() == 6 ) return;
        testpoint = -(a*xbmax + c*zbmax + d)/b;
        if ( (testpoint > (ybmin-GEOMETRY_RESABS)) && (testpoint < (ybmax+GEOMETRY_RESABS)) )
        {
            add_unique_point( intersection_points, CubitVector(xbmax, testpoint, zbmax) );
        }
        if ( intersection_points.size() == 6 ) return;
        testpoint = -(a*xbmin + c*zbmax + d)/b;
        if ( (testpoint > (ybmin-GEOMETRY_RESABS)) && (testpoint < (ybmax+GEOMETRY_RESABS)) )
        {
            add_unique_point( intersection_points, CubitVector(xbmin, testpoint, zbmax) );
        }
        if ( intersection_points.size() == 6 ) return;
        testpoint = -(a*xbmin + b*ybmin + d)/c;
        if ( (testpoint > (zbmin-GEOMETRY_RESABS)) && (testpoint < (zbmax+GEOMETRY_RESABS)) )
        {
            add_unique_point( intersection_points, CubitVector(xbmin, ybmin, testpoint) );
        }
        if ( intersection_points.size() == 6 ) return;
        testpoint = -(a*xbmax + b*ybmin + d)/c;
        if ( (testpoint > (zbmin-GEOMETRY_RESABS)) && (testpoint < (zbmax+GEOMETRY_RESABS)) )
        {
            add_unique_point( intersection_points, CubitVector(xbmax, ybmin, testpoint) );
        }
        if ( intersection_points.size() == 6 ) return;
        testpoint = -(a*xbmax + b*ybmax + d)/c;
        if ( (testpoint > (zbmin-GEOMETRY_RESABS)) && (testpoint < (zbmax+GEOMETRY_RESABS)) )
        {
            add_unique_point( intersection_points, CubitVector(xbmax, ybmax, testpoint) );
        }
        if ( intersection_points.size() == 6 ) return;
        testpoint = -(a*xbmin + b*ybmax + d)/c;
        if ( (testpoint > (zbmin-GEOMETRY_RESABS)) && (testpoint < (zbmax+GEOMETRY_RESABS)) )
        {
            add_unique_point( intersection_points, CubitVector(xbmin, ybmax, testpoint) );
        }
    }
}

//===========================================================================
//Function Name: add_unique_point
//Description:   Given a point and a list of points, add the point to the list
//               unless it is already in the list within GEOMETRY_RESABS.
//Author:
//Date: 3/2005
//===========================================================================
void FBDataUtil::add_unique_point
(
    DLIList<CubitVector > &points,
    const CubitVector &pt
)
{
    int ipt;
    for ( ipt = 0; ipt < points.size(); ipt++ )
    {
        double dist = pt.distance_between( points[ipt] );
        if ( dist <= GEOMETRY_RESABS )
        {
            return;
        }
    }
    points.append( pt );
}

//===========================================================================
//Function Name: FBmake_xy_plane
//Description:   Makes a triangle plane centered on the origin in the x-y 
//               direction.  
//Author: jdfowle
//Date: 1/2004
//===========================================================================
CubitStatus FBDataUtil::FBmake_xy_plane(std::vector<double>& verts, 
                                        std::vector<int>& conns, 
                                        double xsize, 
                                        double ysize, 
                                        int numx, 
                                        int numy)
{
int i, j;
double xc, yc, zc, xinc, yinc;

  zc = 0.;
  xinc = xsize/(double)(numx-1);
  yinc = ysize/(double)(numy-1);
  xc = -0.5*xsize;
//  Make the coordinates
  for ( i = 0; i < numx; i++ ) {
    yc = -0.5*ysize;
    for ( j = 0; j < numy; j++ ) {    
       verts.push_back(xc);
       verts.push_back(yc);
       verts.push_back(zc);
       yc += yinc;
    }
    xc += xinc; 
  }
//  Make the connections
int i1, i2, i3, i4;

  for ( i = 0; i < numx-1; i++ ) {
  
    for ( j = 0; j < numy-1; j++ ) {
      i1 = numy*i + j;
      i2 = numy*(i+1) + j;
      i3 = numy*i + j + 1;
      i4 = numy*(i+1) + j + 1;
      conns.push_back(i1);
      conns.push_back(i2);
      conns.push_back(i4);
      conns.push_back(i1);
      conns.push_back(i4);
      conns.push_back(i3);          
    }
  }
  
  return CUBIT_SUCCESS;
}

//===========================================================================
//Function Name: rotate_FB_object
//Description:   Rotates a facetedboolean object so that the z-direction 
//               points in the direction of NormalDir
//Author: jdfowle
//Date: 1/2004
//===========================================================================
CubitStatus FBDataUtil::rotate_FB_object(std::vector<double>& verts,
                             CubitVector &NormalDir, CubitVector &CenterPt)
{
unsigned int i; 
double l1, l2, l3, m1, m2, m3, n1, n2, n3;
double tnx, tny, tnz, temp;
double tx, ty, tz;
double xpcen, ypcen, zpcen;

  tnx = NormalDir.x();
  tny = NormalDir.y();
  tnz = NormalDir.z();
  xpcen = CenterPt.x();
  ypcen = CenterPt.y();
  zpcen = CenterPt.z();
  
  n1 = tnx; n2 = tny; n3 = tnz;
  l1 = n3; l2 = 0.; l3 = -n1;
  m1 = -n1*n2; m2 = n1*n1 + n3*n3; m3 = -n2*n3;
  temp = sqrt(l1*l1 + l3*l3);
  if ( fabs(temp) > 1.e-6 ) {
    l1 /= temp; l3 /= temp;
    temp = sqrt(m1*m1 + m2*m2 + m3*m3);
    m1 /= temp; m2 /= temp; m3 /= temp;
  } else {
    l1 = 1.; l2 = l3 = 0.;
    m1 = m2 = 0.; m3 = 1.;
    n1 = n3 = 0.; n2 = -1.;
  }
//  Rotate the plane, whiuch is normal to the Z axis and through the origin,
//  so that it will be normal to threepointnormal.  Also translate it so
// that the center point, (0,0,0), moves to threepointcenter.
  for ( i = 0; i < verts.size(); i += 3 ) {
    tx = verts[i];
    ty = verts[i+1];
    tz = verts[i+2];
    verts[i] = tx*l1 + ty*m1 + tz*n1 + xpcen;
    verts[i+1] = tx*l2 + ty*m2 + tz*n2 + ypcen;
    verts[i+2] = tx*l3 + ty*m3 + tz*n3 + zpcen;
  } 
  
  return CUBIT_SUCCESS;
}

int FBDataUtil::makeahashvaluefrom_coord(double x, double y, double z, int numhashbins)
{
double sum;

      if ( fabs(x) < 1.e-3 ) x = 0.0;
      if ( fabs(y) < 1.e-3 ) y = 0.0;
      if ( fabs(z) < 1.e-3 ) z = 0.0;
      sum = (int)(10000.0*fabs(x) + 0.5) + 
                    (int)(10000.0*fabs(y) + 0.5) + 
                    (int)(10000.0*fabs(z) + 0.5);
      
      return (int)(sum) % numhashbins;
}
 
CubitStatus FBDataUtil::FBmake_cylinder(std::vector<double>& verts, 
                                        std::vector<int>& coords, 
                                        double radius, 
                                        double length, 
                                        int nr, 
                                        int nz)
{
  int i, j;
  CubitStatus status;   
  double cfac, rinc, linc;
  double x, y, z;
  int istart, iend, V3, pend;
  double zoffset, lpos, rpos, xrad, yrad;
  
  status = CUBIT_SUCCESS;
  
  rinc = 360.0/(double)nr;
  linc = length/(double)nz;
  cfac = CUBIT_PI/180.;
  
    istart = 0; iend = nz+1;
    V3 = (nz+1)*nr;
    pend = nz;
  
  
  //  Make the points.
  
  zoffset = 0.0;
  lpos = -0.5*length; 
  xrad = radius;
  yrad = radius;
  for ( i = istart; i < iend; i++ ) {
    rpos = 10.0;
    for ( j = 0; j < nr; j++ ) {
      x = xrad*cos(cfac*rpos);
      y = yrad*sin(cfac*rpos);
      z = lpos;
      verts.push_back(x);
      verts.push_back(y);
      verts.push_back(z);      
      rpos += rinc;   
    }
    lpos += linc;
    zoffset += linc;
  } 
  //  Add the two points on the axis at the ends.
  verts.push_back(0.);
  verts.push_back(0.);
  verts.push_back(-0.5*length);
  verts.push_back(0.);
  verts.push_back(0.);
  verts.push_back(0.5*length);
    
  //  Make the triangles.
  int vertnum;
  vertnum = 0;
  for ( i = 0; i < pend; i++ ) {
    for ( j = 0; j < nr-1; j++ ) {
//      facet_ptr = new CubitFacetData( points[vertnum+j],points[vertnum+j+1], points[vertnum+j+nr] );
    coords.push_back(vertnum+j);
    coords.push_back(vertnum+j+1);
    coords.push_back(vertnum+j+nr);
    
//      facet_ptr = new CubitFacetData( points[vertnum+j+1],points[vertnum+j+1+nr], points[vertnum+j+nr] );
    coords.push_back(vertnum+j+1);
    coords.push_back(vertnum+j+1+nr);
    coords.push_back(vertnum+j+nr);
 
    }
//    facet_ptr = new CubitFacetData( points[vertnum],points[vertnum+nr], points[vertnum+2*nr-1] );
    coords.push_back(vertnum);
    coords.push_back(vertnum+nr);
    coords.push_back(vertnum+2*nr-1);

//    facet_ptr = new CubitFacetData( points[vertnum+nr-1],points[vertnum], points[vertnum+2*nr-1] );
    coords.push_back(vertnum+nr-1);
    coords.push_back(vertnum);
    coords.push_back(vertnum+2*nr-1);
    
    vertnum += nr;
  }
  
  //  Endcap(s)
  for ( i = 0; i < nr-1; i++ ) { // top cap
//    facet_ptr = new CubitFacetData( points[vertnum+i],points[vertnum+i+1], points[V3+1] );
    coords.push_back(vertnum+i);
    coords.push_back(vertnum+i+1);
    coords.push_back(V3+1);

  }   
//  facet_ptr = new CubitFacetData( points[nr-1+vertnum],points[vertnum], points[V3+1] );
    coords.push_back(vertnum+nr-1);
    coords.push_back(vertnum);
    coords.push_back(V3+1);
  
  
  for ( i = 0; i < nr-1; i++ ) { // bottom cap
//    facet_ptr = new CubitFacetData( points[i+1],points[i], points[V3] );
    coords.push_back(i+1);
    coords.push_back(i);
    coords.push_back(V3);
 
  }   
//  facet_ptr = new CubitFacetData( points[0],points[nr-1], points[V3] );
    coords.push_back(0);
    coords.push_back(nr-1);
    coords.push_back(V3);
 
  return status;
  
}

double FBDataUtil::project_point_to_plane(double *point, double a, double b, 
                                          double c, double d, 
                                          double *projected_pt)
{
double signed_distance;

  signed_distance = point[0]*a + point[1]*b + point[2]*c + d;
  projected_pt[0] = point[0] - signed_distance*a;
  projected_pt[1] = point[1] - signed_distance*b;
  projected_pt[2] = point[2] - signed_distance*c;


  return signed_distance;

}


#define EPS_SEG_SEG 1.e-6
//  From Schneider and Eberly,"Geometric Tools for Computer Graphics",
//  Chap. 10.8, segment to segment
//  Note:  if the segments are parallel, *s and *t are undefined; *parallel
//  is true; *sunclipped  = 0.0, and *tunclipped is the corresponding value
//  for t.
double FBDataUtil::closest_seg_seg_dist(double *p0, double *d0, double *p1, 
                             double *d1, double *s, double *t, 
                             double *sunclipped, double *tunclipped, 
                             bool *parallel)
{
double ux, uy, uz;
double a, b, c, d, e, det;
double snum, sdenom, tnum, tdenom;

  *parallel = false;
  
  ux = p0[0] - p1[0];
  uy = p0[1] - p1[1];
  uz = p0[2] - p1[2];
  
  a = d0[0]*d0[0] + d0[1]*d0[1] + d0[2]*d0[2];
  b = d0[0]*d1[0] + d0[1]*d1[1] + d0[2]*d1[2];
  c = d1[0]*d1[0] + d1[1]*d1[1] + d1[2]*d1[2];
  d = d0[0]*ux + d0[1]*uy + d0[2]*uz;  
  e = d1[0]*ux + d1[1]*uy + d1[2]*uz;  
  det = a*c - b*b;

  if ( det < EPS_SEG_SEG ) {
    snum = 0.;
    tnum = e;
    tdenom = c;
    sdenom = 1.;
    *sunclipped = snum/sdenom;
    *tunclipped = tnum/tdenom;    
    double f = ux*ux + uy*uy + uz*uz;
    *parallel = true; 
    return sqrt(f - e*e/c); 
  } else {
    snum = b*e - c*d;
    tnum = a*e - b*d;
    sdenom = det;
    *sunclipped = snum/det;
    *tunclipped = tnum/det;    
  }


  if ( snum < 0.0 ) {
    snum = 0.0;
    tnum = e;
    tdenom = c;
  } else if ( snum > det ) {
    snum = det;
    tnum = e + b;
    tdenom = c;
  } else {
    tdenom = det;
  }
  
  if ( tnum < 0.0 ) {
    tnum = 0.0;
    if ( -d < 0.0 ) {
      snum = 0.0;
    } else if ( -d > a ) {
        snum = sdenom;
    } else {
        snum = -d;
        sdenom = a;
    }
  } else if ( tnum > tdenom ) {
      tnum = tdenom;
    if ( (-d + b) < 0.0 ) {
        snum = 0.0;
    } else if ( (-d + b) > a ) {
          snum = sdenom;
    } else {
          snum = -d + b;
          sdenom = a;
    }
  }

  *s = snum/sdenom;
  *t = tnum/tdenom;

  double vx, vy, vz;
  vx = p0[0] + (*s*d0[0]) - (p1[0] + (*t*d1[0]));
  vy = p0[1] + (*s*d0[1]) - (p1[1] + (*t*d1[1]));
  vz = p0[2] + (*s*d0[2]) - (p1[2] + (*t*d1[2]));

  return sqrt(vx*vx + vy*vy + vz*vz);
       
}
