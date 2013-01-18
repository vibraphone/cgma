//-----------------------------------------------------------------------------
//
//  File: CubitFacetEdge.hpp
//
//  Purpose:  optionally used with the CubitFacet class if information
//            is required at the edges of facets     
//
//  Notes:    Note that this class does not contain any private data.
//            All data should be defined within the child classes inherited
//            from this class.  The current Cubit data class that inherits
//            from CubitFacetEdge is CubitFacetEdgeData.  This is done so that 
//            other applications using CubitFacetEdges can use their own
//            edge data, but take advantage of the CGM/Cubit functionality.
//            Please do not add private data to this class; instead add the
//            data to the children and access through virtual functions.
// 
//            Do not create a CubitFacetEdge directly.  For example, don't do: 
//              CubitFacetEdge *cfe = new CubitFacetEdge(...);
//            You should instead create the appropriate child class, and 
//            cast it to a CubitFacetEdge for use.  For example:
//              CubitFacetEdge *cfe = (CubitFacetEdge *) new CubitFacetEdgeData(...);
//                        
//-----------------------------------------------------------------------------

#ifndef CUBITFACETEDGE_HPP
#define CUBITFACETEDGE_HPP

// Include for CubitBoolean
#include "CubitDefines.h"
#include "CubitVector.hpp"
#include "DLIList.hpp"
#include "FacetEntity.hpp"

class CubitFacet;
class CubitPoint;
class CubitBox;

class CubitFacetEdge : public FacetEntity
{
private:

protected:

  CubitVector controlPoints[3];
    //- bezier points on the edge

  int bezierOrder;
    //- bezier order

  int markedFlag;
    //- marked flag

  IttyBit isFeature;
    //- set if this edge is a feature

  CubitBoolean isFlipped;
 
public:
  CubitFacetEdge();
    //- constructors
  virtual ~CubitFacetEdge();
    //- destructor
  virtual int id() = 0;
  virtual void set_id(int /*ent_id*/) {};

  int is_flipped()
    {
      if(isFlipped)
        return 1;
      return 0;
    }
  void is_flipped ( int flipped)
    {
      if(flipped)
        isFlipped = CUBIT_TRUE;
      else
        isFlipped = CUBIT_FALSE;
    }

  void toggle_is_flipped()
    {
      isFlipped = !isFlipped;
    }
  

  virtual CubitPoint *point( int index ) = 0;
    //- get one of its points
  virtual CubitFacet *adj_facet( int index) = 0;
    //- get one of its adjacent facets
  virtual int bezier_order() {return bezierOrder;};
  virtual int control_points( CubitVector *ctrl_pts );
    //- get the control points (return the order)
  virtual CubitVector *control_points() { return controlPoints; };
  virtual void control_points( CubitVector *ctrl_pts, int order );
    //- set the bezier control points on the edge.
    //- end points are assumed to be the first and last
    //- ctrl pt so they are not passed in.  Pass in order-1
    //- control points (max order = 4)
  virtual CubitStatus control_points( CubitFacet *facet, CubitVector *ctrl_pts );
    //- return the control pons on the edge of a facet based
    //- on its edge use direction
  virtual void get_control_points( CubitPoint *point_ptr, CubitVector ctrl_pts[3] );
  virtual void set_control_points( CubitPoint *point_ptr, CubitVector ctrl_pts[3] );
    //- get and set control points oriented with respect to the
    //- point_ptr (ie.  point_ptr is the first control point on the edge)
    //- Note: only gets and sets the middle three control points.  The
    //- other 2 are the edge vertices.
  void set_control_points( const double *ctrl_pt_array );
  virtual void set_flag( int my_flag ) {markedFlag = my_flag;};
  virtual int get_flag( ) {return markedFlag;};
  virtual void facets(DLIList<CubitFacet*> &facet_list ) = 0;
  virtual void edges(DLIList<CubitFacetEdge*> &edge_list ) = 0;
  virtual void points(DLIList<CubitPoint*> &point_list ) = 0;
  virtual int num_adj_facets() = 0;
  void tris(DLIList<CubitFacet*> &facet_list){ facets(facet_list); }
  void tris(int* /*tool_id*/,
	    DLIList<CubitFacet*> &facet_list){ facets(facet_list); }
  

  //- Implement in the child class if need be. 
  //- Assertion will occur if not implemented in the child class
  virtual void add_facet(CubitFacet * /*facet_ptr*/){ assert(0); }
  virtual CubitStatus remove_facet(CubitFacet * /*facet_ptr*/) = 0;
  virtual CubitPoint *start_node() = 0;
  virtual CubitPoint *end_node() = 0;
  virtual void flip() = 0;
 
  virtual int number_tris() { return num_adj_facets(); }
  virtual int number_faces() { return 0; }
 
  virtual void marked (int my_flag ) { set_flag(my_flag); }
  virtual int marked() { return get_flag(); }


  CubitStatus evaluate_position( const CubitVector &start_position,
                                 CubitVector *eval_point,
                                 CubitVector *eval_tangent);
    //- find closet point on non-linear edge
  CubitStatus evaluate( double &t, 
                        CubitVector *eval_point,
                        CubitVector *eval_tangent );
  CubitStatus evaluate_single(double &t,
                              CubitVector *outv);
  CubitStatus evaluate_single_tangent(double &t,
                                      CubitVector *outv);
  CubitStatus evaluate_2nd_derivative(double &t,
                                      CubitVector *outv);
    //- evaluate location -1 < t < 1
  CubitStatus closest_point(const CubitVector &point, 
                            CubitVector &closest_point );
    //- return closest point to linear segment

  CubitStatus intersect(CubitVector &aa, CubitVector &bb, // end point of segment
                        CubitVector &norm,  // normal of the common plane
                        CubitVector &qq,  // return the intersection point 
                        CubitBoolean &does_intersect );
    // intersect the edge with a segment.  Assumes segment and edge
    // are on the same plane (project to facet plane first)

  void boundary_edge_points( CubitPoint * &pt0, 
                             CubitPoint * &pt1,
                             int tool_id = 0);
    // return oriented points on a boundary facet edge

  double dist_to_edge( const CubitVector &this_point, 
                       CubitVector &close_point, 
                       CubitBoolean &outside_edge );
    // return distance from point to an edge

  CubitStatus proj_to_line( const CubitVector &this_point, 
                            CubitVector &proj_point );
    // project point to line defined by edge

  CubitStatus edge_tangent( const CubitVector &point_on_edge, 
                            CubitVector &tangent );
    // compute tangent vector of edge
  CubitStatus edge_curvature( const CubitVector &point_on_edge, 
                              CubitVector &curvature,
                              CubitFacetEdge *closest_edge );
    // compute curvature vector of edge
  double length();
  
  CubitVector position_from_fraction( double zero_to_one );
  CubitVector center() { return position_from_fraction(0.5); }
  
    // compute and return the edge length
  CubitPoint* other_point( CubitPoint *point_ptr );
    // return the other point on the facet edge

  void get_parents(DLIList<FacetEntity *> &facet_list);
    // return the adjacent facets to edge

  CubitFacet *other_facet( CubitFacet *facet_ptr );
    // return the other facet at the edge

  CubitFacet *other_facet_on_surf( CubitFacet *facet_ptr );
    // return the other facet at the edge that has the
    // same tool id

  int num_adj_facets_on_surf( int tool_id );
    // return number of adjacent facets with the indicated
    // tool id

  CubitFacet *adj_facet_on_surf( int tool_id );
    // return first facet adjacent the edge with the
    // indicated tool id

  CubitBoolean contains( CubitPoint *point_ptr );
    //- determines if point is contained in edge

  void set_as_feature() { isFeature = 1; }
  CubitBoolean is_feature( ){return (isFeature ? CUBIT_TRUE : CUBIT_FALSE); }
    // set and get the isFeature bit

  void debug_draw(int color = -1, int flush = 1, int draw_uv=0);
    // debug drawing

  CubitPoint *shared_point( CubitFacetEdge *edge_ptr );
    // get the common point
  void add_facets( );
    // add this edge to its adjacent facets

  CubitBox bounding_box();
    // return the bounding box of the edge

  static int intersect_2D_segments( double P0[2], double P1[2],
                             double P2[2], double P3[2],
                             double qq[4] );
  static int intersect_intervals( double u0, double u1,
                                  double v0, double v1,
                                  double w[2] );
  static CubitStatus order_edge_list(DLIList<CubitFacetEdge*> &edge_list,
                                     CubitPoint *start_point,
                                     CubitPoint *&end_point);
  static CubitPoint *find_start_point_for_edge_list(DLIList<CubitFacetEdge*> edge_list);

  double angle_between_facets();

  inline int less_than(CubitFacetEdge*& e1, CubitFacetEdge*& e2)
  {
    double len1 = e1->length();
    double len2 = e2->length();
  
    if (len1 == len2) return 0;
    return (len1 < len2) ? -1 : 1;
  }

};

inline void CubitFacetEdge::control_points(
  CubitVector *ctrl_pts, int order ) 
{
  assert(order > 0 && order <=4);
  bezierOrder = order;
  for(int i=0; i<order-1; i++){
    controlPoints[i] = ctrl_pts[i];
  }
}

//======================================================================
// Function: get_control_points (PUBLIC)
// Description: get control points oriented with respect to the
//              point_ptr (ie.  point_ptr is the first control point 
//              on the edge)
// Note: only gets and sets the middle three control points.  The
//       other 2 are the edge vertices.
// Author: sjowen
// Date: 05/01
//======================================================================
inline void CubitFacetEdge::get_control_points( CubitPoint *point_ptr, 
                                                CubitVector ctrl_pts[3] )
{
  DLIList<CubitPoint*> my_points;
  points(my_points);
  
  if (point_ptr == my_points.get())
  {
    ctrl_pts[0] = controlPoints[0];
    ctrl_pts[1] = controlPoints[1];
    ctrl_pts[2] = controlPoints[2];
  }
  else if(point_ptr == my_points.next())
  {
    ctrl_pts[0] = controlPoints[2];
    ctrl_pts[1] = controlPoints[1];
    ctrl_pts[2] = controlPoints[0];
  }
  else
  {
    assert(0);  // point_ptr does not match either point
  }
}

//======================================================================
// Function: set_control_points (PUBLIC)
// Description: set control points oriented with respect to the
//              point_ptr (ie.  point_ptr is the first control point 
//              on the edge)
// Note: only gets and sets the middle three control points.  The
//       other 2 are the edge vertices.
// Author: sjowen
// Date: 05/01
//======================================================================
inline void CubitFacetEdge::set_control_points( CubitPoint *point_ptr, 
                                                CubitVector ctrl_pts[3] )
{
  DLIList<CubitPoint*> my_points;
  points(my_points);

  if (point_ptr == my_points.get())
  {
    controlPoints[0] = ctrl_pts[0];
    controlPoints[1] = ctrl_pts[1];
    controlPoints[2] = ctrl_pts[2];
  }
  else if(point_ptr == my_points.next())
  {
    controlPoints[0] = ctrl_pts[2];
    controlPoints[1] = ctrl_pts[1];
    controlPoints[2] = ctrl_pts[0];
  }
  else
  {
    assert(0);  // point_ptr does not match either point
  }
}

#endif


