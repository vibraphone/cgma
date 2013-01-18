//- Class:       TDFacetBoundaryEdge
//- Description: Tool data for storing additional information at 
//-              the boundary of a facet set
//- Owner:       Steve Owen
//- Checked by:
//- Version:

#ifndef TD_FACET_BOUNDARY_EDGE_HPP


#define TD_FACET_BOUNDARY_EDGE_HPP

#include "CubitDefines.h"
#include "ToolData.hpp"
#include "DLIList.hpp"
#include "CubitVector.hpp"
#include "MemoryManager.hpp"
#include "CastTo.hpp"
class CubitFacetEdge;
class CubitFacet;
class CubitPoint;

typedef struct BoundaryEdgeData
{
  CubitVector bezierCtrlPts[3];
  CubitFacet *adjFacet;
  int surfID;
  int computed;
} BoundaryEdgeData;

class TDFacetBoundaryEdge : public ToolData
{
private:

  static MemoryManager memoryManager;
    //- memory management object

  CubitFacetEdge *edgePtr;
    //Because non-manifold geometry complicates things a little, if
    // we are using a represenation that is of a higher order than
    // linear.  We figure out which edges edgePtr should be C1 with
    // along a curve... The two edges (if exist) are store in nextEdgePtr
    // (if it is the edge across node(1)) or prevEdgePtr (if it is the
    // edge across node(0)).
  CubitFacetEdge *nextEdgePtr;
  CubitFacetEdge *prevEdgePtr;
  DLIList <BoundaryEdgeData *> edgeDataList;

public:

  TDFacetBoundaryEdge();
    //- constructor

  ~TDFacetBoundaryEdge();

  static int is_facet_boundary_edge(const ToolData* td)
     {return (CAST_TO(const_cast<ToolData*>(td), TDFacetBoundaryEdge) != NULL);}
  
  void add_surf(int new_id);

  CubitFacetEdge *get_edge()
    { return edgePtr; }
 
  void set_edge( CubitFacetEdge *edge_ptr )
    { edgePtr = edge_ptr; }

  CubitFacetEdge *get_next_edge()
    { return nextEdgePtr; }
 
  void set_next_edge( CubitFacetEdge *edge_ptr )
    { nextEdgePtr = edge_ptr; }

  CubitFacetEdge *get_prev_edge()
    { return prevEdgePtr; }
 
  void set_prev_edge( CubitFacetEdge *edge_ptr )
    { prevEdgePtr = edge_ptr; }
  
  
  SetDynamicMemoryAllocation(memoryManager)
    //- class specific new and delete operators
    
  static void set_memory_allocation_increment(int increment = 0)
    {memoryManager.set_memory_allocation_increment(increment);}
    //- set block memory size increment
  
  static void destroy_memory()
    {memoryManager.destroy_memory();}
    //- destroy all memory allocted to this object

  static CubitStatus add_facet_boundary_edge(CubitFacetEdge *edge_ptr);
  static TDFacetBoundaryEdge* get_facet_boundary_edge(CubitFacetEdge *edge_ptr);

  int control_points( CubitVector *ctrl_pts, int surf_id );
    //- get the control points (return the order)
  void control_points( CubitVector *ctrl_pts, int order, int surf_id );
    //- set the bezier control points on the edge.
    //- end points are assumed to be the first and last
    //- ctrl pt so they are not passed in.  Pass in order-1
    //- control points (max order = 4)
  CubitStatus control_points( CubitFacet *facet, CubitVector *ctrl_pts );
    // return the control pons on the edge of a facet based
    // on its edge use direction
  void get_control_points( CubitPoint *point_ptr, CubitVector ctrl_pts[3],
                           int surf_id );
  void set_control_points( CubitPoint *point_ptr, CubitVector ctrl_pts[3],
                           int surf_id );
    // get and set control points oriented with respect to the
    // point_ptr (ie.  point_ptr is the first control point on the edge)
    // Note: only gets and sets the middle three control points.  The
    // other 2 are the edge vertices.

  void add_surf_facet( DLIList<CubitFacet *> &facet_list );
    // add facet representing a new surface adjacent this edge

  void set_surf_id( CubitFacet *facet_ptr, int surf_id );
    // set id of surface defined by facet_ptr

  CubitBoolean is_internal_edge();
    // return whether this edge is internal to a surface

  CubitStatus init_control_points( double min_dot );
    // initialize control points on boundary edge

  CubitStatus merge_control_points();
    // find aerage of control points and set on edge

  CubitStatus compute_curve_tangent( int surf_id,
                                     CubitFacetEdge *edge,
                                     double min_dot,
                                     CubitVector &T0,
                                     CubitVector &T3 );
    // compute the tangents to the endpoints of a feature edge.

  CubitFacetEdge *next_feature_edge( int surf_id,
                                     CubitFacetEdge *this_edge,
                                     CubitPoint *p0 );
    // given a facet boundary edge and one of its nodes, find the
    // next edge on the same surface 

  CubitBoolean is_at_surf( int surf_id );
    // return whether the given edge has one of its adjacent 
    // faces on the given surface
};
    

#endif // TD_FACET_BOUNDARY_EDGE_HPP


