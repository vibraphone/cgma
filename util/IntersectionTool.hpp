// file: GeomIntersectionTool.hpp
// author: Michael Stephenson
//

#ifndef INTERSECTION_TOOL_HPP
#define INTERSECTION_TOOL_HPP

#include "CubitDefines.h"
#include "CubitUtilConfigure.h"

class CubitVector;
template <class X> class DLIList;

class CUBIT_UTIL_EXPORT IntersectionTool
{
public:
  IntersectionTool(const double &tol = 0.0001);
  virtual ~IntersectionTool() {}

  virtual CubitStatus closest_points_on_segments(CubitVector &p0,
                                                 CubitVector &p1,
                                                 CubitVector &p2,
                                                 CubitVector &p3,
                                                 CubitVector &point_1,
                                                 CubitVector &point_2,
                                                 double &sc, double &tc);
    //- Finds the two closest points on the line segments {p0,p1} and {p2,p3},
    //- and returns them in point_1 and point_2.  Will return failure
    //- if a divide by zero is attempted.

  double distance_point_line(const double point[3], const double start[3], 
                             const double end[3], double &t);
  int point_on_polyline(CubitVector& pt, DLIList<CubitVector*> &pt_list,
                         double *tol_in = NULL);

  double parametric_position(const double node[3],
                           const double pt1[3],
                           const double pt2[3]);
  
    ///
    /// Returns the parametric position of the node on the
    /// line segement between the points pt1 and pt2.  This
    /// value is between 0 and 1.  If the node is closest
    /// to pt1 it will return 0, if it is closest to pt2, it
    /// will return 1.  If the node is off the line, it
    /// returns the closest parametric value.
    /// If there is an error, like pt1 and pt2 are equal,
    /// then the function returns CUBIT_DBL_MAX.
    ///

  virtual CubitStatus initialize(){return CUBIT_SUCCESS;}

  // The following copyright applies to the following two functions...
  //
  // Copyright 2001, softSurfer (www.softsurfer.com)
  //
  // This code may be freely used and modified for any purpose
  // providing that this copyright notice is included with it.
  // SoftSurfer makes no warranty for this code, and cannot be held
  // liable for any real or imagined damage resulting from its use.
  // Users of this code must verify correctness for their application.

  static int intersect_triangle_with_ray( CubitVector &ray_origin, CubitVector &ray_direction,
	  const CubitVector *p0, const CubitVector *p1, const CubitVector *p2,
	  CubitVector* point, double &distance, int &edge_hit );
    //- Find intersection point of a ray and a triangle
    //    Return: -1 = triangle is degenerate (a segment or point)
    //             0 = disjoint (no intersect)
    //             1 = intersect at unique point
    //             2 = are in the same plane

  static int intersect_segment_with_ray( CubitVector &ray_origin, CubitVector &ray_direction,
	  const CubitVector *p0, const CubitVector *p1,
	  CubitVector* point, double &distance, int &point_hit, double tol=0.0 );
    //- Find intersection point of a ray and a facet edge
    //    Return: -1 = edge is degenerate (a point)
    //             0 = disjoint (no intersect)
    //             1 = intersect at unique point
    //             2 = are the same line (infinite points)

  static int intersect_point_with_ray( CubitVector &ray_origin, CubitVector &ray_direction, 
	  const CubitVector* point, double &distance, double tol=0.0);
    //- Find intersection of a ray and a point
    //	  Return: 0 = no intersection
    //			  1 = intersection


protected:
  double mTolerance;
  double mToleranceSquared;

  CubitBoolean ray_tri_test(const double start[3], const double dir[3],
                            const double tri1[3], const double tri2[3], 
                            const double tri3[3],
                            double &t, double &alpha, double &beta);

  CubitBoolean skew_line_test(const double start1[3], const double end1[3],
                              const double start2[3], const double end2[3],
                              double &t, double &u);

  void tolerance(const double &tol);
  double tolerance(void) const;
  
private:
};

inline IntersectionTool::IntersectionTool(const double &tol) 
    : mTolerance(tol)
{
  mToleranceSquared = tol * tol;
}

#endif // INTERSECTION_TOOL_HPP

