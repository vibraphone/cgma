//- Class: CurveFacetEvalTool
//- Description:  The CurveFacetEvalTool is a general purpose tool that uses a set
//-               of facets to determine certain geometry computations.  This
//-               class can be used to define a curve or do fast bounding box
//-               calculations, or determine point inside or outside.
//- Assumptions:  It is assumed/required that the facets be continuous through
//-               the set.  In other words they must shared nodes and be conected
//-               through each other.  The algorithms assume, for efficient searching,
//-               that from one facet, you can get to any other facet in the set.
//- Owner: Steven J. Owen
//- Checked by: 

#ifndef CURVE_FACET_EVAL_TOOL_HPP
#define CURVE_FACET_EVAL_TOOL_HPP

#include "DLIList.hpp"
#include <map>

class CubitBox;
class CubitPoint;
class CubitVector;
class FacetEvalTool;
class CubitFacetEdge;
class FacetCurve;

class CurveFacetEvalTool
{
private:
  int toolID;
  int interpOrder;
    //- interpolation order 0=linear, 1=gradient, 2=quadratic, 
    //                      3=least_squares, 4=spline
  DLIList<CubitFacetEdge*> myEdgeList;
  DLIList<CubitPoint*> myPointList;
  CubitBox *myBBox;
  FacetEvalTool *surfFacetEvalTool;
  double facetLength;
  CubitSense curvSense;
  CubitBoolean goodCurveData;
  int output_id;

  void draw_edges(int color = -1);
  void draw_edge(CubitFacetEdge *edge, int color);
  void draw_line(CubitVector &begin, CubitVector &end, int color = -1);
  void draw_location(CubitVector &loc, int color = -1 );

  void destroy_facets();
    //- Destroys the facets, and points.

  CubitStatus get_segments_from_loops(
    DLIList<DLIList<CubitFacetEdge*>*> *facet_loop_list,
    CubitVector &start, CubitVector &end,
    DLIList<CubitFacetEdge*> &edge_list,
    DLIList<CubitPoint*> &point_list,
    CubitSense owrts );
  CubitStatus get_segments_from_loops(
    DLIList<DLIList<CubitFacetEdge*>*> *facet_loop_list,
    CubitPoint *start_pt, CubitPoint *end_pt,
    DLIList<CubitFacetEdge*> &edge_list,
    DLIList<CubitPoint*> &point_list,
    CubitSense owrts );
    // generate a list of edges and points from the loop list
    // for this curve

  CubitStatus get_segments_from_positions( 
    std::vector<CubitVector> &positions,
    CubitVector &start_pt,
    DLIList<CubitFacetEdge*> &edge_list, 
    DLIList<CubitPoint*> &point_list );
      // generate a list of edges and points from the loop list
    // for this curve
  
  /*
  CubitStatus project_to_facet_edge(CubitVector &this_point,
                                    CubitVector &closest_point,
                                    CubitVector *tangent_ptr,
                                    CubitVector *curvature_ptr,
                                    double *param,
                                    int *outside );

    // project to the facet edge
  CubitStatus project_to_linear_facet_edge(CubitVector &this_point,
                                           CubitVector &closest_point,
                                           CubitVector *tangent_ptr,
                                           CubitVector *curvature_ptr,
                                           double *param,
                                           int *outside );
    // same as above except assumes linear representation AND
    // surfFacetEvalTool does notexist
*/
  // evaluates a fractional length onto a Bezier edge 
  CubitStatus evaluate_bezier_edge(CubitFacetEdge *edge,
                                   int index0,
                                   int index1,
                                   double fraction,
                                   CubitVector &location_on_curve,
                                   double *tangent);

  // projects a point onto a Bezier edge 
  CubitStatus project_to_bezier_edge(CubitFacetEdge *edge,
                        CubitVector &point, CubitVector &projected_point,
                        double *tangent, double *tval);                                   
                                        
  double u_on_facet_edge( CubitFacetEdge *edge_at_pt, 
                          CubitVector pt );
    // return the u param on the facet curve given a point on the
    // curve and the facet edge it is on

//  CubitStatus project_to_edge_line( CubitFacetEdge *edge, 
//                                    CubitVector &this_point, 
//                                    CubitVector &pt_on_edge, 
//                                    double &dist_to_edge );
    // project to the line defined by the edge

  CubitStatus fix_point_edge_order();

  CubitSense find_curv_sense( CubitPoint *start_ptr );
  CubitSense find_curv_sense( CubitVector &start );

public:

  CurveFacetEvalTool();
  ~CurveFacetEvalTool();

  // Called immediately after the constructor to initialize the data.
  CubitStatus initialize( DLIList<CubitFacetEdge*> &edge_list,
                          DLIList<CubitPoint*> &point_list,
                          FacetEvalTool* surf_eval = NULL);
  CubitStatus initialize( FacetEvalTool *surf_eval_tool,
                          CubitPoint *start_point,
                          CubitPoint *end_point,
                          CubitSense orientation );
  CubitStatus initialize( FacetEvalTool *surf_eval_tool,
                          CubitVector &start,
                          CubitVector &end,                          
                          CubitSense orientation );
  CubitStatus initialize( FacetEvalTool *surf_eval_tool,
                          CubitVector &start,
                          std::vector<CubitVector> &positions );

  CubitStatus save(FILE *fp);
    // save to a cubit file
  CubitStatus restore(FILE *fp,
                      unsigned int endian,
                      int num_edges, 
                      int num_points,
                      CubitFacetEdge **edges, 
                      CubitPoint **points,
                      int num_fets,
                      FacetEvalTool **fet_list);

  int get_output_id() { return output_id; }
  void set_output_id( int id ) { output_id = id; }

  CubitBox bounding_box();
    //- Returns the bounding box for the set of facetedges (based on the points
    //- used in the faceted set.

  CubitStatus closest_point(CubitVector &this_point, 
                            CubitVector &closest_point_ptr,
                            CubitVector *normal_ptr = NULL,
                            CubitVector *curvature_ptr = NULL,
                            double *param = NULL);
    //- Finds the closest point from the vector (this_point) to the
    //- set of facets that lies on the set of facets.  If the point
    //- lies outside this set, the closest point will be on the plane
    //- of the closest facet.  The closest_point is set to be that point.

  // replace del_pnt with keep_pnt in point list
  CubitBoolean replace_point( CubitPoint *del_pnt, CubitPoint *keep_pnt );

  CubitBoolean replace_facets( DLIList< CubitFacetEdge *> &curv_edges );

  void remove_facets( DLIList<CubitFacetEdge*> &facet_edges);

  CubitStatus position_from_fraction( double fraction, // between 0 and 1
                                      CubitVector &location_on_curve );
    //- computes the location on the curve based on a fraction of the
    //- distance along the curve and the RefEdge sense

  double u_from_arc_length( double root_param, double arc_length );
  double length_from_u ( double root_param, double end_param );

  double length() { return facetLength; };
    //- return the length of the facet edges

  CubitSense sense() { return curvSense; };
    //- return the sense of the curve

  void get_facets(DLIList<CubitFacetEdge*>& facet_list)
    { facet_list +=  myEdgeList; }
    // append the edge facets for this curve onto facet_list
  void get_points(DLIList<CubitPoint*>& point_list)
    { point_list +=  myPointList; }
    // apend the points for this curve onto point_list

  FacetEvalTool *get_surf_eval_tool() { return surfFacetEvalTool; } 
  void set_facet_eval_tool(FacetEvalTool *surf_eval_tool_ptr) 
   { surfFacetEvalTool = surf_eval_tool_ptr;} 
    // get and set the associated surface facet eval tool 

  void set_length();
    // initialize the curve length
  
  CubitBoolean has_good_curve_data(){return goodCurveData;}
  
  void debug_draw_facet_edges(int color = -1 );

   
};

#endif // CURVE_FACET_EVAL_TOOL_HPP


