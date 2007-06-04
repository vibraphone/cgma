//- Class: FacetEvalTool
//- Description:  The FacetEvalTool is a general purpose tool that uses a set
//-               of facets to determine certain geometry computations.  This
//-               class can be used to define a surface or do fast bounding box
//-               calculations, or determine point inside or outside.
//- Assumptions:  It is assumed/required that the facets be continuous through
//-               the set.  In other words they must shared nodes and be conected
//-               through each other.  The algorithms assume, for efficient searching,
//-               that from one facet, you can get to any other facet in the set.
//- Owner: Steven J. Owen
//- Checked by: 

#ifndef SMOOTH_FACET_EVAL_TOOL_HPP
#define SMOOTH_FACET_EVAL_TOOL_HPP

#include "DLIList.hpp"
#include "AbstractTree.hpp"
#include "CubitFacet.hpp" //For inline function.
template <class Y> class KDDTree;

#define determ3(p1,q1,p2,q2,p3,q3) ((q3)*((p2)-(p1)) + (q2)*((p1)-(p3)) + (q1)*((p3)-(p2)))
#define sqr(a) ((a)*(a))
#define cube(a) (sqr(a) * (a))
#define quart(a) (sqr(a) * sqr(a))
#define blend(x) (-2.0*(x)*(x)*(x) + 3.0*(x)*(x))

class CubitBox;
class CubitFacet;
class CubitPoint;
class CubitVector;
class CubitTransformMatrix;

class FacetEvalTool
{
private:
  static double timeGridSearch;
  static double timeFacetProject;
  static int numEvals;

  int toolID;
  int interpOrder;
    //- interpolation order 0=linear, 1=gradient, 2=quadratic, 4=quartic Bezier
  DLIList<CubitFacet*> myFacetList;
    //- Facets that are used for evaluation.
  DLIList<CubitFacet*> markedFacetList;
    //- Facets that are used for evaluation.
  DLIList<CubitPoint*> myPointList;
    //- point list
  DLIList<CubitFacetEdge*> myEdgeList;
  DLIList<DLIList<CubitFacetEdge*>*> myLoopList;
    //- edges only used for bezier interp
  AbstractTree <CubitFacet*> *aTree;
    //- used for surfaces with many facets
  CubitBox *myBBox;
    //- Bounding box of facets
  double myArea;
    //- area of the facet surface
  double compareTol;
    //- for comparing bounding box
  CubitFacet *lastFacet;
    //- The last facet evaluated
  int isFlat;
    //- The surface represented by the facets is flat
  CubitBoolean isParameterized;
    //- Surface contains a parameterization (uv values initialized)
  double minDot;  
    //- computed from feature angle
  int output_id;

  void set_up_grid_search(double geom_factor);
  CubitStatus facets_from_search_grid( CubitVector &this_point,
                                       DLIList<CubitFacet *> &facet_list );
    //- set up a grid search data structure if we have lots of facets

  CubitStatus get_points_from_facets(DLIList<CubitFacet*> &facet_list,
                                     DLIList<CubitPoint*> &point_list );
    //- populates the point_list with points contained by the given facets.
  CubitStatus get_edges_from_facets(DLIList<CubitFacet*> &facet_list,
                                    DLIList<CubitFacetEdge*> &edge_list );
    //- populates the edge_list with edges contained by the given facets

  CubitStatus get_loops_from_facets(DLIList<CubitFacetEdge*> &all_edge_list,
                                    DLIList<DLIList<CubitFacetEdge*>*> &loop_list );
    //- populates the loop_list with edges contained by the given edge list
  void destroy_facets();
    //- Destroys the facets, and points.

  void debug_draw_edges(int color = -1 );
  void debug_draw_edge(CubitFacetEdge *edge, int color = -1 );
  void debug_draw_vec( CubitVector & vec, int color = -1 );
  void debug_draw_facet_normals( int color = -1 );
  void debug_draw_point_normals( int color = -1 );
  void debug_draw_bezier_edges( int color = -1 );
  void debug_draw_bezier_facets( int color = -1 );
  void debug_draw_bezier_facet( CubitFacet *facet, int color = -1 );
  void debug_draw_line( CubitVector &begin, CubitVector &end, int color = -1 );
  void debug_draw_eval_bezier_facet( CubitFacet *facet );
  void debug_draw_location( CubitVector &location, int color = -1 );
  void write_loops();
    //- drawing functions for debugging.

  CubitStatus project_to_facets(CubitVector &this_point,
                                CubitBoolean trim,
                                CubitBoolean *outside,
                                CubitVector *closest_point_ptr,
                                CubitVector *normal_ptr);

  CubitStatus init_gradient();
    //- initiaize the gradients (tangent planes at the points)

  CubitStatus init_quadrics();
    //- initialize quadric coefficients at the points

  CubitStatus eval_quadratic( CubitFacet *facet, 
                              int pt_idx, 
                              CubitVector &eval_pt,
                              CubitVector &qpoint,
                              CubitVector &qnorm );
    //- evaluate on facet based on quadratic approximation

  void on_circle( CubitVector &ptA,
                  CubitVector &normA,
                  CubitVector &ptB,
                  CubitVector &normB,
                  CubitVector &eval_pt,
                  CubitVector &pt_on_circle,
                  CubitVector &normal );
    //- compute a projection of a point onto a circle.  The circle
    //- is defined by two points and two normals.  Return the 
    //- projected point and its normal 

  void eval_quadric( CubitFacet *facet,
                     int pt_index,
                     CubitVector &eval_pt,
                     CubitVector &qpt );
    //- evaluate on facet based on quadratic approximation

  CubitStatus get_close_points( CubitPoint *point,
                                CubitPoint **close_point,
                                int &num_close,
                                int max_close,
                                int min_close );
    //- return an array of points that are close to the given point

  static CubitStatus eval_bezier_patch( CubitFacet *facet, CubitVector &areacoord,
                               CubitVector &pt );
  static CubitStatus eval_bezier_patch_normal( CubitFacet *facet, 
                                               CubitVector &areacoord,
                                               CubitVector &normal );
  static CubitStatus hodograph( CubitFacet *facet, 
                                CubitVector &areacoord,
                                CubitVector Nijk[10] );
  static CubitStatus project_to_patch( CubitFacet *facet,                                                
                                CubitVector &ac,
                                CubitVector &pt,
                                CubitVector &eval_pt,
                                CubitVector *eval_norm,
                                CubitBoolean &outside,
                                double compare_tol,
                                int edge_id /* = -1*/);
    //- Project a point to a bezier patch. Pass in the area
    //- of the point projected to the linear facet.  Function
    //- assumes that the point is contained within the patch -
    //- if not, it will project to one of its edges.

  static void ac_at_edge( CubitVector &fac,
                          CubitVector &eac,
                          int edge_id );
    //- determine the area coordinate of the facet at the edge

  static CubitBoolean is_at_vertex( CubitFacet *facet,
                                    CubitVector &pt,
                                    CubitVector &ac,
                                    double compare_tol,
                                    CubitVector &eval_pt,
                                    CubitVector *eval_norm_ptr );

  static CubitBoolean move_ac_inside( CubitVector &ac, double tol );
    //-  find the closest area coordinate to the boundary of the
    //-  patch if any of its components are < 0
    //-  Return if the ac was modified.

  void check_faceting();

  
    // For a given point, find which pairs of edges
    // should be C1 across that point.
  CubitStatus mark_edge_pairs(CubitPoint* point);

    //loop over the points and call the above function
  CubitStatus pair_edges();
  
public:
  
  FacetEvalTool( DLIList<CubitFacet*> &facet_list, DLIList<CubitPoint*> &point_list,
                       int interp_order/* = -1*/, double min_dot/* = 0.707106781185*/ );
  FacetEvalTool();
  //constructor
  ~FacetEvalTool();
    //destructor

  int tool_id()
    { return toolID; };

  CubitStatus reverse_facets();
  
  int get_output_id() { return output_id; }
  void set_output_id( int id ) { output_id = id; }

  CubitStatus save( FILE *fp );
    // save to a cubit file
  CubitStatus restore( FILE *fp, unsigned int endian,
                       int num_facets, int num_edges, 
                       int num_points, CubitFacet **facets, 
                       CubitFacetEdge **edges, CubitPoint **points );
  
    //- get and set the interpolation order
  CubitBox bounding_box();
  void set_bounding_box( CubitBox &box ) { *myBBox = box; }
  void reset_bounding_box();
    //- Returns the bounding box for the set of facets (based on the points
    //- used in the faceted set.

  double area();
    //- return the surface area of the facets

  void calculate_area();
    //- (re) calcualte area to make sure it is up to date.

  double get_min_dot() {return minDot;};
    //- return the mindot (feature angle)

  CubitStatus closest_point(CubitVector &this_point, 
                            CubitVector *closest_point_ptr,
                            CubitVector *normal_ptr = NULL);
    //- Finds the closest point from the vector (this_point) to the
    //- set of facets that lies on the set of facets.  If the point
    //- lies outside this set, the closest point will be on the plane
    //- of the closest facet.  The closest_point is set to be that point.
  
  CubitStatus closest_point_trimmed( CubitVector &this_point, 
                                     CubitVector *closest_point_ptr,
                                     CubitBoolean &lies_outside,
                                     CubitVector *normal_ptr = NULL);
    //- Finds the closest point from the vector (this_point) to the
    //- set of facets that lies on the set of facets.  If the point
    //- lies outside this set, the closest point will be on the edge
    //- of the closest facet.  
    //- This function also determines if the point is outside or
    //- inside the current set of facets.  You should call
    //- a bounding box test first before this...

  CubitFacet* closest_facet( const CubitVector& point );

  int is_flat();
    //- Determine if the set of facets are flat (all in the same plane)

  static CubitStatus eval_facet( CubitFacet *facet, 
                                 CubitVector &areacoord,
                                 CubitVector *eval_point,
                                 CubitVector *eval_normal );
    //- evaluate the facet at the areacoords

  CubitStatus eval_facet( CubitFacet *facet, 
                          CubitVector &pt,
                          CubitVector &areacoord, 
                          CubitVector &close_point,
                          CubitBoolean &outside_facet );
    //- evaluate a location and normal on (or near) a facet based on the
    //- the facet area coordinates

  static CubitStatus project_to_facet( CubitFacet *facet, 
                                       CubitVector &pt,
                                       CubitVector &areacoord,
                                       CubitVector &close_point,
                                       CubitBoolean &outside,
                                       double compare_tol);
    //- project a point to the given facet (assumes that it contains
    //- the point) areacoord is an initial guess - may be changed.

  static CubitStatus project_to_facets(DLIList <CubitFacet *> &facet_list,
                                 CubitFacet *&last_facet,
                                 int interp_order,
                                 double compare_tol,
                                 CubitVector &this_point,
                                 CubitBoolean trim,
                                 CubitBoolean *outside,
                                 CubitVector *closest_point_ptr,
                                 CubitVector *normal_ptr);
    //- project a point to the a list of facets using the interpolation
    //- order.  Option to trim to boundary

  static CubitStatus eval_edge( CubitFacet *facet, 
                         int vert0, int vert1,
                         CubitVector &pt_on_plane, 
                         CubitVector &close_point);
  static CubitStatus project_to_facetedge( CubitFacet *facet, 
                         int vert0, int vert1,
                         CubitVector &the_point,
                         CubitVector &pt_on_plane, 
                         CubitVector &close_point,
                         CubitBoolean &outside_facet, 
                         CubitBoolean must_be_on_edge = CUBIT_FALSE);
    //- project a point to the edge of a facet
    //- if must_be_on_edge is true, close_point will be on the facet edge

  static CubitStatus eval_point( CubitFacet *facet, 
                          int vertex_id, 
                          CubitVector &close_point );
    //- evaluate a location at a vertex of a facet

  static CubitStatus eval_facet_normal( CubitFacet *facet,
                                 CubitVector &areacoord,
                                 CubitVector &normal );
    //- evaluate the normal on a facet (use the interpOrder)

  static void project_to_facet_plane( CubitFacet *facet,
                               CubitVector &pt,
                               CubitVector &point_on_plane,
                               double &dist );
    //- project a point to the plane of a facet

  static void facet_area_coordinate( CubitFacet *facet,
                                     CubitVector &pt_on_plane,
                                     CubitVector &areacoord );
    //- define the area coordinates of a point on a plane of the facet

  static CubitBoolean is_outside( CubitFacet *facet, 
                           CubitVector &areacoord );
    //- determines if the point at the areacoord is outside 
    //- the range of facets

  DLIList<DLIList<CubitFacetEdge*>*> *loops() { return &myLoopList; };
    //- return the edge loop list

  int interp_order() { return interpOrder; };
    //- return the interpolation order of the facet eval tool

  void get_facets( DLIList<CubitFacet*> &facet_list ) { facet_list += myFacetList; };
  void get_edges( DLIList<CubitFacetEdge*> &edge_list) { edge_list += myEdgeList; };
  void get_points( DLIList<CubitPoint*> &point_list ) { point_list += myPointList; };
  void remove_facets( DLIList<CubitFacet*>& facet_list );
    //- return the facets and points defining the surface (append to list)

  void replace_facets( DLIList<CubitFacet *> &facet_list );
    // replace the facets on this FacetEvalTool

  int num_facets() const {return myFacetList.size();};
    //- get the number of facets
      
  int num_points() const {return myPointList.size();};
    //- get the number of points

  double compare_tol() { return compareTol; };
  void compare_tol(double tol) { compareTol = tol; };
    //- set/get the comparison tolerance.

  static double get_grid_search_time()
    { return timeGridSearch;}
  static double get_facet_project_time()
    { return timeFacetProject;}
    //- Returns static timing information for the two functions.

  void debug_draw_facets(int color = -1 );

  void transform_control_points( CubitTransformMatrix &rotmat );
    // transform all control points on the surface a specified 
    // transformation matrix

  CubitStatus init_bezier_surface( );
  static CubitStatus init_bezier_facet( CubitFacet *facet );
  static CubitStatus init_bezier_edge( CubitFacetEdge *edge, double min_dot );
  static CubitStatus init_boundary_bezier_edge( CubitFacetEdge *edge );
  static CubitStatus init_edge_control_points( CubitVector &P0,
                                               CubitVector &P3,
                                               CubitVector &N0,
                                               CubitVector &N3,
                                               CubitVector &T0,
                                               CubitVector &T3,
                                               CubitVector Pi[3] );
  static CubitStatus init_edge_control_points_single( 
                              CubitVector &P0,
                              CubitVector &P3,
                              CubitVector &T0,
                              CubitVector &T3,
                              CubitVector Pi[3] );                                      
  static CubitStatus init_facet_control_points( CubitVector N[6],
                                                CubitVector P[3][5],
                                                CubitVector G[6] );
    //- Initialize the bezier control points for 
    //- G1 bezier patch surface interpolation

  static CubitBoolean is_boundary_facet( CubitFacet *facet_ptr );
    //- decides if this is a boundary facet based on whether its
    //- edges or points have boundary tool datas attached.

  CubitStatus parameterize();
    //- define a parameterization

  static CubitStatus compute_curve_tangent( CubitFacetEdge *edge,
                                     double min_dot,
                                     CubitVector &T0,
                                     CubitVector &T3 );
    //- compute curve tangent vectors at an aedge on the boundary
  static CubitFacetEdge *next_boundary_edge( CubitFacetEdge *this_edge, CubitPoint *p0 );
    //- return the next edge on the boundary

};

#endif // SMOOTH_FACET_EVAL_TOOL_HPP



    

