//-----------------------------------------------------------------------------------
// Class: PointLoopFacetor
// Description:
// Creates a set of CubitFacets which is a Delauney triagulation of the boundary
// points.  No interior points are added.  The points form a closed loop, and may
// consist of several loops.  The first loop must be the exterior loop, followed by
// any internal loops (or holes).
// Creator: David R. White
// Owner: David R. White
// Creation Date: 3/1/2003
//------------------------------------------------------------------------------------
#ifndef POINTLOOPFACETOR_HPP
#define POINTLOOPFACETOR_HPP

template <class X> class DLIList;
class CubitPoint;
class CubitFacet;

typedef DLIList<CubitPoint*> PointList;
typedef DLIList<PointList*> PointLoopList;
#include "CubitDefines.h"

class PointLoopFacetor
{
private:
  
public:
  PointLoopFacetor()
    {}
  ~PointLoopFacetor()
    {}
  static CubitStatus generate_facets( PointLoopList &boundary_loops,
                                      DLIList<CubitFacet*> &resulting_triangles);
    ///
    /// This function does the triangulation.
    ///

  static void write_xy(PointLoopList &boundary_loops);
  static void write_xyz(PointLoopList &boundary_loops);
  
};
#endif

