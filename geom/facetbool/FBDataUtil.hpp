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

#ifndef _FBDATAUTIL
#define _FBDATAUTIL
#include <vector>
#include "CubitBox.hpp"
#include "CubitVector.hpp"
#include "CubitDefines.h"
#include "DLIList.hpp"

//===========================================================================
//  Arguments:
//  bbox (INPUT) the bounding box to intersect the plane with
//  v0 (INPUT) the coordinates of the first point in the plane
//  v1 (INPUT) the coordinates of the second point in the plane
//  v2 (INPUT) the coordinates of the third point in the plane
//  intersection_points (INPUT/OUTPUT) the array of points, up to six, that
//  are the intersections of the box edges with the plane
//  Returnn Value:  the number of intersection points
//===========================================================================
class FBDataUtil {

public:

static void intersect_plane_with_boundingbox(CubitBox &bbox, 
                                            const CubitVector &v0, 
                                            const CubitVector &v1,
                                            const CubitVector &v2,
                                            DLIList<CubitVector> &intersection_points); 

//===========================================================================
//  Arguments:
//  verts (INPUT/OUTPUT) the array of vertex coordinates (xyz,xyz,etc.)
//  coords (INPUT/OUTPUT) the array of connections
//  xsize (INPUT) the length in the x direction
//  ysize (INPUT) the length in the y direction
//  numx (INPUT) the number of points along the x direction, at least 2
//  numy (INPUT) the number of points along the y direction, at least 2
//  makes 2*(numx-1)*(numy-1) triangles in the x-y plane centered on the origin
//===========================================================================
static CubitStatus FBmake_xy_plane(std::vector<double>& verts, 
                                   std::vector<int>& coords, 
                                   double xsize, 
                                   double ysize, 
                                   int numx, 
                                   int numy);

//===========================================================================
//  Arguments:
//  verts (INPUT/OUTPUT) the array of vertex coordinates (xyz,xyz,etc.)
//  NormalDir (INPUT) the direction to rotate the z-axis into
//  CenterPt (INPUT) the location that (0,0,0) will be translated to
//===========================================================================
static CubitStatus rotate_FB_object(std::vector<double>& verts, 
                             CubitVector &NormalDir, 
                             CubitVector &CenterPt);

static int makeahashvaluefrom_coord(double x, double y, double z, 
                                                int numhashbins);

//===========================================================================
//  Arguments:
//  verts (INPUT/OUTPUT) the array of vertex coordinates (xyz,xyz,etc.)
//  coords (INPUT/OUTPUT) the array of connections
//  radius (INPUT) the radius
//  length (INPUT) the length in the z direction
//  nr (INPUT) the number of points along the radial direction
//  nz (INPUT) the number of points along the z direction
//  
//===========================================================================
static CubitStatus FBmake_cylinder(std::vector<double>& verts, 
                                   std::vector<int>& coords, 
                                   double radius, 
                                   double length, 
                                   int nr, 
                                   int nz);

//===========================================================================
//  Arguments:
//  point (INPUT) point coordinates
//  a, b, c, d (INPUT) plane coefficients
//  projected_pt (outPUT) projection of the point onto the plane
//  Returns signed distance of point to plane.
//===========================================================================
static double project_point_to_plane(double *point, double a, double b, 
                                     double c, double d, 
                                     double *projected_pt);

//===========================================================================
//  Finds closest distance between two 3-D line segments.  Also returns parametric
//  distance along each segment where the shortest line joining the two segments
//  intersects the segments.  And a boolean that says whether the segments are
//  parallel or not.  If parallel, *s and *t are undefined.
//  Arguments:
//  p0 (INPUT) beginning coordinates of first segment
//  d0 (INPUT) coordinates of end of first segment minus beginning
//  p1 (INPUT) beginning coordinates of second segment
//  d1 (INPUT) coordinates of end of second segment minus beginning
//  s (OUTPUT) Clipped parametric distance along first segment.  0. <= s <= 1.
//  t (OUTPUT) Clipped parametric distance along second segment.  0. <= s <= 1.
//  sunclipped (OUTPUT) parametric distance along first segment.  May be outside range.
//  tunclipped (OUTPUT) parametric distance along second segment.  May be outside range.
//  parallel (OUTPUT) true if segments are parallel; otherwise, false.
//  
//===========================================================================
static double closest_seg_seg_dist(double *p0, double *d0, 
                             double *p1, double *d1, double *s, double *t, 
                             double *sunclipped, double *tunclipped, 
                             bool *parallel);

//===========================================================================
//  Finds closest distance between two 3-D line segments.  Also returns parametric
//  distance along each segment where the shortest line joining the two segments
//  intersects the segments.  And a boolean that says whether the segments are
//  parallel or not.  If parallel, *s and *t are undefined.
//  Arguments:
//  p0 (INPUT) beginning coordinates of first segment
//  d0 (INPUT) coordinates of end of first segment minus beginning
//  p1 (INPUT) beginning coordinates of second segment
//  d1 (INPUT) coordinates of end of second segment minus beginning
//  s (OUTPUT) Clipped parametric distance along first segment.  0. <= s <= 1.
//  t (OUTPUT) Clipped parametric distance along second segment.  0. <= s <= 1.
//  sunclipped (OUTPUT) parametric distance along first segment.  May be outside range.
//  tunclipped (OUTPUT) parametric distance along second segment.  May be outside range.
//  parallel (OUTPUT) true if segments are parallel; otherwise, false.
//  
//===========================================================================
                             
                             
private:

//===========================================================================
//Function Name: add_unique_point
//Description:   Given a point and a list of points, add the point to the list
//               unless it is already in the list within GEOMETRY_RESABS.
//Author:
//Date: 3/2005
//===========================================================================
static void add_unique_point( DLIList<CubitVector > &points,
                              const CubitVector &pt );


};

#endif
