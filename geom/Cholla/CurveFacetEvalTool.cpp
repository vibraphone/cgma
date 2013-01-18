//- Class: CurveFacetEvalTool
//- Description:  The CurveFacetEvalTool is a tool to perform generic geometric 
//-               operations on a set of facets.
//- Assumptions:  All of the facets are connected topologically correct.
//-
//- Owner: Steve J. Owen
//- Checked by: 

#include "CubitVector.hpp"
#include "CubitPoint.hpp"
#include "CubitFacet.hpp"
#include "CubitFacetEdge.hpp"
#include "CubitFacetEdgeData.hpp"
#include "CubitBox.hpp"
#include "DLIList.hpp"
#include "CurveFacetEvalTool.hpp"
#include "FacetEvalTool.hpp"
#include "GeometryDefines.h"
#include "CastTo.hpp"
#include "GfxDebug.hpp"
#include "CubitFileIOWrapper.hpp"
#include "Cholla.h"
#include "FacetDataUtil.hpp"

//===========================================================================
//Function Name: initialize
//
//Member Type:  PUBLIC
//Description:  initialize a CurveFacetEvalTool - uses end locations of the curve
//              (does geometric comparisons to match facet points)
//===========================================================================
CubitStatus CurveFacetEvalTool::initialize(
  FacetEvalTool *surf_eval_tool,   // the adjacent surface facet eval tool
  CubitVector &start, CubitVector &end,  // begin and end vertex  
  CubitSense orientation_wrt_surface )  // direction wrt the surface  
{
//   static int counter = 0;
  output_id = -1;
  myBBox = NULL;

  surfFacetEvalTool = surf_eval_tool;  

  CubitStatus stat = get_segments_from_loops( surfFacetEvalTool->loops(), start, end, 
                           myEdgeList, myPointList, orientation_wrt_surface );
  goodCurveData = (stat == CUBIT_SUCCESS) ? CUBIT_TRUE : CUBIT_FALSE;
  if( !goodCurveData )
  {
    PRINT_ERROR( "Unable to initialize Curve Evaluation Tool for Faceted Geometry.\n" );
    return CUBIT_FAILURE;
  }

  interpOrder = surfFacetEvalTool->interp_order();
  curvSense = find_curv_sense( start );
  bounding_box();
  return CUBIT_SUCCESS;
}

CubitStatus CurveFacetEvalTool::initialize(
  FacetEvalTool *surf_eval_tool,   // the adjacent surface facet eval tool
  CubitVector &start, // begin vertex  
  std::vector<CubitVector> &positions )  
{
//   static int counter = 0;
  output_id = -1;
  myBBox = NULL;

  surfFacetEvalTool = surf_eval_tool;  
  
  CubitStatus stat = get_segments_from_positions( positions, start,
                           myEdgeList, myPointList ); 

  goodCurveData = (stat == CUBIT_SUCCESS) ? CUBIT_TRUE : CUBIT_FALSE;
  if( !goodCurveData )
  {
    PRINT_ERROR( "Unable to initialize Curve Evaluation Tool for Faceted Geometry.\n" );
    return CUBIT_FAILURE;
  }
  
  interpOrder = surfFacetEvalTool->interp_order();
  curvSense = find_curv_sense( start );
  bounding_box();
  return CUBIT_SUCCESS;
}

//===========================================================================
//Function Name: initialize
//
//Member Type:  PUBLIC
//Description:  initialize a CurveFacetEvalTool - uses end point pointers of the curve
//===========================================================================
CubitStatus CurveFacetEvalTool::initialize(
  FacetEvalTool *surf_eval_tool,             // the adjacent surface facet eval tool
  CubitPoint *start_pt, CubitPoint *end_pt,  // begin and end points  
  CubitSense orientation_wrt_surf )          // direction wrt the surface
{
//   static int counter = 0;
  output_id = -1;

  myBBox = NULL;
  surfFacetEvalTool = surf_eval_tool;
  CubitStatus stat = get_segments_from_loops( surfFacetEvalTool->loops(), 
                                              start_pt, end_pt, 
                                              myEdgeList, myPointList,
                                              orientation_wrt_surf );
  goodCurveData = (stat == CUBIT_SUCCESS) ? CUBIT_TRUE : CUBIT_FALSE;
  if( !goodCurveData )
  {
    PRINT_ERROR( "Unable to initialize Curve Evaluation Tool for Faceted Geometry.\n" );
    return CUBIT_FAILURE;
  }
  interpOrder = surfFacetEvalTool->interp_order();
  curvSense = find_curv_sense( start_pt );
  bounding_box();
  return CUBIT_SUCCESS;
}

//===========================================================================
//Function Name: CurveFacetEvalTool
//
//Member Type:  PUBLIC
//Description:  initialize a CurveFacetEvalTool - uses an existing set of
//              ordered facet edges to define curve
//===========================================================================
CubitStatus CurveFacetEvalTool::initialize(
  DLIList<CubitFacetEdge*> &edge_list, // the ordered facet edges on this curve
  DLIList<CubitPoint*> &point_list, FacetEvalTool* surf_eval)    // the ordered points on this curve
{
//   static int counter = 0;
  output_id = -1;
  myBBox = NULL;
  curvSense = CUBIT_FORWARD;
  surfFacetEvalTool = surf_eval;
  myEdgeList = edge_list;
  myPointList = point_list;
  set_length();
  
  CubitStatus status = fix_point_edge_order();

  if( CUBIT_SUCCESS != status )
    return status;

  interpOrder = 0;
  if(surf_eval)
    interpOrder = surf_eval->interp_order();
  bounding_box();
  goodCurveData = CUBIT_TRUE;
  return CUBIT_SUCCESS;
}

CurveFacetEvalTool::CurveFacetEvalTool()
{
  static int counter = 0;
  toolID = counter++;
  myBBox = NULL;
  goodCurveData = CUBIT_TRUE;
  
}
//===========================================================================
//Function Name: ~CurveFacetEvalTool
//
//Member Type:  PUBLIC
//Description:  destructor
//===========================================================================
CurveFacetEvalTool::~CurveFacetEvalTool()
{
  if ( myBBox ) delete myBBox;
  destroy_facets();
}

//===========================================================================
//Function Name: get_segments_from_loops
//
//Member Type:  PRIVATE
//Description:  extract the facet edges from the loops list that are 
//              used for the current curve.  Uses CubitVectors at curve
//              ends to determine curve -- has to do geometric comaprisons
//              so it may be less reliable than the overloaded function
// assumption: the sense (orientation) of the RefFace is the same as that of 
//             the facets and the facetedge loop list
//===========================================================================
CubitStatus CurveFacetEvalTool::get_segments_from_loops(
  DLIList<DLIList<CubitFacetEdge*>*> *facet_loop_list,  // lists of boundary edges
  CubitVector &start, CubitVector &end,        // begin and end vertex
  DLIList<CubitFacetEdge*> &edge_list,            // return the edges used for this curve
  DLIList<CubitPoint*> &point_list,              // return the points used for this curve
  CubitSense orientation_wrt_surf )
{
  CubitStatus stat = CUBIT_SUCCESS;
  int ii, jj;
  int mydebug = 0;
  
  CubitBoolean done = CUBIT_FALSE;
  DLIList<CubitFacetEdge*> *facet_loop;
  CubitFacetEdge *edge, *startedge = NULL;
  CubitPoint *pt0, *pt1;
  CubitVector loc_pt0, loc_pt1;
  double edge_length;

  facetLength = 0.0e0;

  int tool_id = surfFacetEvalTool->tool_id();
    //loop over the "loops"
  for (ii=0; ii<facet_loop_list->size() && !done && stat == CUBIT_SUCCESS; ii++) {
      //get the next facet loop
    facet_loop = facet_loop_list->get_and_step();   

      // now loop over the edges in that loop
    for (jj=0; jj<facet_loop->size() && !done; jj++) {
        //get the start edge... we will 
        //traverse the list forwards, if the surface is oriented forward
        //otherwise, traverse it backwards
      startedge =
        (orientation_wrt_surf == CUBIT_FORWARD) ? facet_loop->get_and_step() :
        facet_loop->get_and_back();
      
      if (startedge->get_flag() == 0) {
        if (orientation_wrt_surf == CUBIT_FORWARD) {
          startedge->boundary_edge_points( pt0, pt1, tool_id );
        }
        else {
          startedge->boundary_edge_points( pt1, pt0, tool_id );
        }
        loc_pt0 = pt0->coordinates();
        loc_pt1 = pt1->coordinates();
        if (loc_pt0.within_tolerance( start, GEOMETRY_RESABS )) {
          startedge->set_flag( 1 );
          edge = startedge;
          if (mydebug) draw_edge( edge, CUBIT_GREEN );          
          while (!done) {
            edge_list.append( edge );
            point_list.append( pt0 );
            edge_length = loc_pt0.distance_between( loc_pt1 );
            facetLength += edge_length;
              //if the other end of the edge is at the begginning of the
              //loop, then we're done.
            if (loc_pt1.within_tolerance( end, GEOMETRY_RESABS )) {
              done = CUBIT_TRUE;
              point_list.append( pt1 );
            }
            else {//otherwise...
              if (orientation_wrt_surf == CUBIT_FORWARD) {
                edge = facet_loop->get_and_step();
                edge->boundary_edge_points( pt0, pt1, tool_id );
              }
              else {
                edge = facet_loop->get_and_back();
                edge->boundary_edge_points( pt1, pt0, tool_id );
              }
              loc_pt0 = pt0->coordinates();
              loc_pt1 = pt1->coordinates();
              if (mydebug) draw_edge( edge, CUBIT_BLUE );              
                //if we got back to the startedge without the
                // other end of an edge getting to the start vertex,
                // we have a problem....
              if (edge == startedge)
              {
                stat = CUBIT_FAILURE;  // this shouldn't happen
                done = CUBIT_TRUE;
              }//end edge == startedge
            }//end else not within tolerance
          }//end while not done
        }//end if within tolerance
      }//end if not marked
    }//end loop over edges in facet loop
  }//end loop over facet loops
  if(done!=CUBIT_TRUE)
  {
    PRINT_ERROR("Can't define curve representation in mesh-based geometry\n");
    stat = CUBIT_FAILURE;
  }
  return stat;
  
}

//===========================================================================
//Function Name: get_segments_from_loops
//
//Member Type:  PRIVATE
//Description:  extract the facet edges from the loops list that are 
//              used for the current curve.  Uses the CubitPoints at
//              the ends of the curve to determine curve
// assumption: the sense (orientation) of the RefFace is the same as that of 
//             the facets and the facetedge loop list
//===========================================================================
CubitStatus CurveFacetEvalTool::get_segments_from_loops(
  DLIList<DLIList<CubitFacetEdge*>*> *facet_loop_list,  // lists of boundary edges
  CubitPoint *start_pt, CubitPoint *end_pt,   // begin and end point
  DLIList<CubitFacetEdge*> &edge_list,            // return the edges used for this curve
  DLIList<CubitPoint*> &point_list,               // return the points used for this curve
  CubitSense orientation_wrt_surf )
{
  CubitStatus stat = CUBIT_SUCCESS;
  int ii, jj;
  int mydebug = DEBUG_FLAG(181);  
  CubitBoolean done = CUBIT_FALSE;
  DLIList<CubitFacetEdge*> *facet_loop;
  CubitFacetEdge *edge, *startedge = NULL;
  CubitPoint *pt0, *pt1;
  CubitVector loc_pt0, loc_pt1;
  double edge_length;
  DLIList<CubitFacetEdge*> temp_edge_list;
  DLIList<CubitPoint*> temp_point_list;
  
  if (mydebug)
  {
    GfxDebug::draw_point( start_pt->x(),
                          start_pt->y(),
                          start_pt->z(),
                          CUBIT_RED );
    GfxDebug::draw_point( end_pt->x(),
                          end_pt->y(),
                          end_pt->z(),
                          CUBIT_GREEN );
    GfxDebug::flush();
  }

  int tool_id = surfFacetEvalTool->tool_id();
    //loop over the list of facet loops
  for (ii=0; ii<facet_loop_list->size() && !done && stat == CUBIT_SUCCESS; ii++) {
      //get the next facet loop
    facet_loop = facet_loop_list->get_and_step();
      //now loop over the edges in that facet loop
    for (jj=0; jj<facet_loop->size() && !done; jj++) {
        //get the start edge... we will 
        //traverse the list forwards, if the surface is oriented forward
        //otherwise, traverse it backwards
      startedge =
        (orientation_wrt_surf == CUBIT_FORWARD) ? facet_loop->get_and_step() :
        facet_loop->get_and_back();
        //if startedge isn't marked
      if (startedge->get_flag() == 0) {
          //get the points from the start edge in the appropriate order
        if (orientation_wrt_surf == CUBIT_FORWARD) {
          startedge->boundary_edge_points( pt0, pt1, tool_id );
        }
        else {
          startedge->boundary_edge_points( pt1, pt0, tool_id );
        }
          //convert to CubitVector
        loc_pt0 = pt0->coordinates();
        loc_pt1 = pt1->coordinates();
          //if the first point on the edge is the start_pt... otherwise step
        if (pt0 == start_pt) {
            //mark the startedge
          CubitFacetEdge* marked_edge = NULL;
          edge = startedge;
          if (mydebug) draw_edge( edge, CUBIT_GREEN );//draw the start edge
            //loop until we have reached an end or no more choices
          while (!done) {
              //if the first point is at the start, we are either
              // just starting, or we have looped back onto the
              // start point.  The latter _may_ be ok.  In case
              // the latter is the case, clean out the temporary
              // lists and (re)initialize the length.  Also, unmark
              // the edge that we marked the last time we started,
              // if we previously marked an edge.
            if (pt0 == start_pt){
              temp_point_list.clean_out();
              temp_edge_list.clean_out();
              facetLength = 0.0;
              if(marked_edge)
                marked_edge->set_flag(0);
              marked_edge=edge;
              edge->set_flag(1);
            }
              //add the current edge to the edge list
            temp_edge_list.append( edge );
              //add the current start point to the point list
            temp_point_list.append( pt0 );
              //measure the edge... essentially
            edge_length = loc_pt0.distance_between( loc_pt1 );
              //keep a running tally of the curve length
            facetLength += edge_length;
              
              //if the second point is at the end, we're done
            if (pt1 == end_pt ) {
              done = CUBIT_TRUE;
              temp_point_list.append( pt1 );
            }
              //otherwise, we need to step forward
            else {
              if (orientation_wrt_surf == CUBIT_FORWARD) {
                edge = facet_loop->get_and_step();
                edge->boundary_edge_points( pt0, pt1, tool_id );
              }
              else {
                edge = facet_loop->get_and_back();
                edge->boundary_edge_points( pt1, pt0, tool_id );
              }
                //convert point to CubitVector, as above
              loc_pt0 = pt0->coordinates();
              loc_pt1 = pt1->coordinates();
              if (mydebug) draw_edge( edge, CUBIT_BLUE );
              if (edge == startedge)
              {
                stat = CUBIT_FAILURE;  // this shouldn't happen
                done = CUBIT_TRUE;
              }
            }//end else (ie, if not pt1==end_pt)
          }//end while !done
        }//end if pt0 == start_pt
      }//end if start edge isn't marked
    }//end loop over edges in the loop
  }//end loop over facet loops

    //now put the temporary list entries onto the end of the list we return.
  point_list += temp_point_list;
  edge_list += temp_edge_list;
  
  if(done!=CUBIT_TRUE || stat == CUBIT_FAILURE)
  {
    stat = CUBIT_FAILURE;
    PRINT_WARNING("Can't define curve representation in mesh-based geometry\n");
    PRINT_INFO("  Hint: Try importing as a free mesh and examining nodes and attached\n");
    PRINT_INFO("        elements near the following locations:\n");
    PRINT_INFO("        Start curve = %f %f %f\n", start_pt->x(), start_pt->y(), start_pt->z());
    PRINT_INFO("        End curve = %f %f %f\n", end_pt->x(), end_pt->y(), end_pt->z());
    if (mydebug)
    {
      surfFacetEvalTool->debug_draw_facets( CUBIT_YELLOW );
      int i, j;
      for (i=0; i<facet_loop_list->size(); i++)
      {
        DLIList<CubitFacetEdge*> *my_facet_loop = facet_loop_list->get_and_step();
        for (j=0; j<my_facet_loop->size(); j++)
        {
          draw_edge( my_facet_loop->get_and_step(), CUBIT_RED );
        }
      }
      GfxDebug::mouse_xforms();
      mydebug = mydebug;
    }
  }
  return stat;
}

//===========================================================================
//Function Name: get_segments_from_positions
//
//Member Type:  PRIVATE
//Description:  extract the facet edges from the facets of the parent surface.
//              This function is used exclusively for hardlines.  
//===========================================================================
CubitStatus CurveFacetEvalTool::get_segments_from_positions(
  std::vector<CubitVector> &positions,   // postions of points on curve
  CubitVector &start,                    // begin point
  DLIList<CubitFacetEdge*> &edge_list,   // return the edges used for this curve
  DLIList<CubitPoint*> &point_list )     // return the points used for this curve  
{  
  DLIList<CubitPoint*> surface_pts;
  surfFacetEvalTool->get_points( surface_pts );  

  facetLength = 0;

  double smallest_dist = CUBIT_DBL_MAX;
  CubitPoint *start_pt = NULL;
  //find the closest point to 'start'
  for( int k=surface_pts.size(); k--; )
  {
    double tmp_dist = surface_pts.get()->coordinates().distance_between_squared( start );

    if( tmp_dist < smallest_dist )
    {
      smallest_dist = tmp_dist;
      start_pt = surface_pts.get();
    }
    surface_pts.get_and_step();
  }

  myPointList.append( start_pt );

  CubitPoint *current_pt = start_pt;
  CubitPoint *next_pt = NULL;
  for( int k=1; k<positions.size(); k++ )
  {
    CubitVector current_pos = positions[k];

    //find the CubitPoint attached to this point via an edge
    double smallest_dist = CUBIT_DBL_MAX;    
    CubitFacetEdge *next_edge;

    DLIList<CubitFacetEdge*> adj_edges;
    current_pt->edges( adj_edges );

    for( int m=adj_edges.size(); m--; )
    {
      CubitFacetEdge *tmp_edge = adj_edges.get_and_step();
      CubitPoint *other_pt = tmp_edge->other_point( current_pt );

      double tmp_dist = current_pos.distance_between_squared(  other_pt->coordinates() );
      if( tmp_dist < smallest_dist )
      {
        smallest_dist = tmp_dist;
        next_pt = other_pt;
        next_edge = tmp_edge;
      }
    }

    current_pt = next_pt;
    
    facetLength += next_edge->length();

    myEdgeList.append( next_edge );
    myPointList.append( current_pt );    
  }

  if( myPointList.size() != positions.size() ||
      myEdgeList.size() != positions.size()-1 )
    return CUBIT_FAILURE;

  return CUBIT_SUCCESS;
}

//===========================================================================
//Function Name: find_curv_sense
//
//Member Type:  PRIVATE
//Description:  Determines 'curvSense' depending on the CubitFacetEdges 
//              in 'myEdgeList'  Assumption is that all edges
//              are oriented the same in the curve
//===========================================================================
CubitSense CurveFacetEvalTool::find_curv_sense( CubitPoint *start_pt )
{
  myEdgeList.reset();
  CubitFacetEdge *edge_ptr = myEdgeList.get();
  CubitPoint *pt0 = edge_ptr->point(0);
  CubitPoint *pt1 = edge_ptr->point(1);
  CubitSense sense;

  if (start_pt == pt0)
  {
    sense = CUBIT_FORWARD;
  }
  else if (start_pt == pt1)
  {
    sense = CUBIT_REVERSED;
  }
  else
  {
    sense = CUBIT_UNKNOWN;
  }
  return sense;
}

//===========================================================================
//Function Name: find_curv_sense 
//
//Member Type:  PRIVATE
//Description:  Determines 'curvSense' depending on the CubitFacetEdges 
//              in 'myEdgeList' Assumption is that all edges
//              are oriented the same in the curve
//===========================================================================

CubitSense CurveFacetEvalTool::find_curv_sense( CubitVector &start )
{
  myEdgeList.reset();
  CubitFacetEdge* temp_edge = myEdgeList.get();
 
  CubitVector temp_vector = temp_edge->point(1)->coordinates() - 
                            temp_edge->point(0)->coordinates();
  
  double tol = temp_vector.length() / 100.0; 
  //if point(0) of first edge in list == start, it is CUBIT_FORWARD 
  if(start.within_tolerance(temp_edge->point(0)->coordinates(), tol ) )
    return CUBIT_FORWARD;
  //if point(1) of first edge in list == start, it is CUBIT_REVERSED 
  else if(start.within_tolerance(temp_edge->point(1)->coordinates(), tol ) )
    return CUBIT_REVERSED;
  else
    return CUBIT_UNKNOWN; 
}

//===========================================================================
//Function Name: bounding_box
//
//Member Type:  PUBLIC
//Description:  Calculates the bounding box of the curve.
//===========================================================================
CubitBox CurveFacetEvalTool::bounding_box()
{
  if ( !myBBox )
  {
    int ii;
    CubitPoint *point_ptr;
    double x, y, z;
    double x_min = CUBIT_DBL_MAX, x_max = -CUBIT_DBL_MAX;
    double y_min = CUBIT_DBL_MAX, y_max = -CUBIT_DBL_MAX;
    double z_min = CUBIT_DBL_MAX, z_max = -CUBIT_DBL_MAX;
    for ( ii = myPointList.size(); ii > 0; ii-- )
    {
      point_ptr = myPointList.get_and_step();
      x = point_ptr->x();
      y = point_ptr->y();
      z = point_ptr->z();
      if ( x < x_min )
        x_min = x;
      if ( y < y_min )
        y_min = y;
      if ( z < z_min )
        z_min = z;
      if ( x > x_max )
        x_max = x;
      if ( y > y_max )
        y_max = y;
      if ( z > z_max )
        z_max = z;
    }
    CubitVector min_v(x_min, y_min, z_min );
    CubitVector max_v(x_max, y_max, z_max );
    myBBox = new CubitBox( min_v, max_v );
  }
  return *(myBBox);
}

//===========================================================================
//Function Name: position_from_fraction
//
//Member Type:  PUBLIC
//Description:  evaluates the location a fraction of the way along the curve
//              based on the sum of the lengths of the facet edges on this curve
//===========================================================================
CubitStatus CurveFacetEvalTool::position_from_fraction( 
  double fraction, // between 0 and 1
  CubitVector &location_on_curve )
{
  CubitStatus status = CUBIT_SUCCESS;
  CubitBoolean done = CUBIT_FALSE;
  double length, curlength, targetlength, lastlength;
  CubitFacetEdge *edge = NULL;
  CubitPoint *pt0=NULL, *pt1=NULL;
  CubitVector loc_pt0, loc_pt1;
  CubitVector eval_location;
  myEdgeList.reset();

  myPointList.reset();
  pt0 = myPointList.get();
  if (fraction <= 0.0) {
    edge = myEdgeList.get_and_step();
    location_on_curve = pt0->coordinates();
    return status;
  }
  curlength = lastlength = 0.0e0;
  targetlength = fraction * facetLength;
  assert(myEdgeList.size() > 0);
  int ii = 0;
  for ( ; ii<myEdgeList.size() && !done; ii++) {
    edge = myEdgeList.get_and_step();
    loc_pt0 = edge->point( 0 )->coordinates();
    loc_pt1 = edge->point( 1 )->coordinates();
    length = loc_pt0.distance_between( loc_pt1 );
    curlength += length;
    if ( edge->point(0) == pt0 ){
      pt1 = edge->point(1);
    }
    else{
      pt1 = edge->point(0);
    }
    if (targetlength <= curlength) {
      done = CUBIT_TRUE;
      
      loc_pt0 = pt0->coordinates();
      loc_pt1 = pt1->coordinates();
      double local_fraction = (targetlength - lastlength) / length;
      eval_location = loc_pt0 + local_fraction * (loc_pt1 - loc_pt0); 
      if (interpOrder == 0 || !surfFacetEvalTool) {
        location_on_curve = eval_location;
      }
      else {
      //  evaluate Bezier edge
        int index0 = 0;        
        int index1 = 1;
        if ( pt1 == edge->point(0) ) {
            //switching the indices causes evaluate_bezier_edge to compute the 
            //wrong location, so I'm removing the switch (mbrewer)
          //index0 = 1;
          //index1 = 0;
          local_fraction = 1. - local_fraction;
        }
        status = evaluate_bezier_edge(edge,index0,index1,local_fraction,
                                        location_on_curve,0);
      }
    }
    else {
      lastlength = curlength;
    }
    if ( edge->point(0) == pt0 )
        pt0 = edge->point(1);
    else 
        pt0 = edge->point(0);
  }
  
  if (!done) {
    if(pt1 == NULL){
      PRINT_ERROR("Opposite point not found while evaluating a curve.\n");
      location_on_curve.set(0.0,0.0,0.0);
      return CUBIT_FAILURE;
    }
    location_on_curve = pt1->coordinates();
  }
  return status;
}

//===========================================================================
//Function Name: evaluate_bezier_edge
//
//Member Type:  PRIVATE
//Description:  evaluates a fractional length onto a Bezier edge 
//===========================================================================
CubitStatus CurveFacetEvalTool::evaluate_bezier_edge(CubitFacetEdge *edge,
                                        int index0,
                                        int index1,
                                        double fraction,
                                        CubitVector &location_on_curve,
                                        double* tangent)
{
CubitStatus status;
int numEdge, numVert, numLocs;
int edgeVert[2];
//int order;
double vert[6], edgeCtrlPts[9], location[3];
CubitVector *edge_ctrl_pts;
 double scaled_parameter;
 
  numEdge = 1;
  numVert = 2;
  numLocs = 1;
  edgeVert[0] = 0; edgeVert[1] = 1;
  vert[0] = edge->point(index0)->x();
  vert[1] = edge->point(index0)->y();
  vert[2] = edge->point(index0)->z();
  vert[3] = edge->point(index1)->x();
  vert[4] = edge->point(index1)->y();
  vert[5] = edge->point(index1)->z();
  
    edge_ctrl_pts = edge->control_points();
    edgeCtrlPts[0] = edge_ctrl_pts[0].x();  
    edgeCtrlPts[1] = edge_ctrl_pts[0].y();  
    edgeCtrlPts[2] = edge_ctrl_pts[0].z();  
    edgeCtrlPts[3] = edge_ctrl_pts[1].x();  
    edgeCtrlPts[4] = edge_ctrl_pts[1].y();  
    edgeCtrlPts[5] = edge_ctrl_pts[1].z();  
    edgeCtrlPts[6] = edge_ctrl_pts[2].x();  
    edgeCtrlPts[7] = edge_ctrl_pts[2].y();  
    edgeCtrlPts[8] = edge_ctrl_pts[2].z();  

      //evalBezierEdge takes a parameter in the range of -1 to 1.
      // evalBezierEdge calls functions which  convert back to
      // a range of 0 to 1,
      // but we must convert to -1 to 1, before sending it.
    scaled_parameter = (2.0 * fraction) - 1.0;
    
  //  static function in Cholla.cpp
    evalBezierEdge(numEdge,numVert,edgeVert,vert,edgeCtrlPts,numLocs,
                   &scaled_parameter, location,tangent);
  
  location_on_curve.x(location[0]); 
  location_on_curve.y(location[1]); 
  location_on_curve.z(location[2]); 
   
  status = CUBIT_SUCCESS;
  
  return status;
}
                                         
//===========================================================================
//Function Name: u_from_arc_length
//
//Member Type:  PUBLIC
//Description:  see notes in FacetCurve::u_from_arc_length
//===========================================================================
double CurveFacetEvalTool::u_from_arc_length ( double root_param,
                                               double arc_length )
{
  if( facetLength == 0 )
    return 0;

  double u = root_param + arc_length / facetLength;

  return u;

}

//===========================================================================
//Function Name: length_from_u
//
//Member Type:  PUBLIC
//Description:  return length along curve from u parameter
//===========================================================================
double CurveFacetEvalTool::length_from_u ( double root_param,
                                           double end_param )
{
  double length = (end_param - root_param) * facetLength;

  return length;
}

//===========================================================================
//Function Name: closest_point
//
//Member Type:  PUBLIC
//Description:  Finds the closest point from the vector (this_point) to the
//              set of facets that lies on the set of facets.  If the point
//              lies outside this set, the closest point will be on the plane
//              of the closest facet.  The closest_point is set to be that point.
//===========================================================================
CubitStatus CurveFacetEvalTool::closest_point(CubitVector &this_point,
                                              CubitVector &closest_point,
                                              CubitVector *tangent_ptr,
                                              CubitVector *curvature_ptr,
                                              double *param )
{
  CubitStatus stat = CUBIT_SUCCESS;
//  int outside;
  int mydebug = 0;
  if (mydebug)
  {
    draw_location(this_point, CUBIT_RED);
  }
/*
  if (surfFacetEvalTool == NULL)
  {
    stat = project_to_linear_facet_edge( this_point, closest_point,
                                         tangent_ptr, curvature_ptr,
                                         param, &outside );
  }
  else
  {
    stat = project_to_facet_edge( this_point, closest_point, 
                                  tangent_ptr, curvature_ptr, param, &outside );
  }
  */
  int i;
  CubitVector test_pt_v, closest_pt_v, pt0_v, pt1_v;
  double distance2, closest_distance2; //  distance squared
  CubitFacetEdge *edge, *best_edge = NULL;

  myEdgeList.reset();
  closest_distance2 = CUBIT_DBL_MAX;
  for ( i = myEdgeList.size(); i > 0; i-- ) {  
    edge = myEdgeList.get_and_step();
    pt0_v = edge->point(0)->coordinates();
    pt1_v = edge->point(1)->coordinates();
    test_pt_v = FacetDataUtil::squared_distance_to_segment(this_point,
                     pt0_v,pt1_v,distance2);
    if ( distance2 < closest_distance2 ) {
      closest_distance2 = distance2;
      closest_pt_v = test_pt_v;
      best_edge = edge;
      if ( distance2 < GEOMETRY_RESABS*GEOMETRY_RESABS ) break;
    }
  }


  if (mydebug)
  {
    myEdgeList.reset();
    for ( i = myEdgeList.size(); i > 0; i-- ) {  
      edge = myEdgeList.get_and_step();
      pt0_v = edge->point(0)->coordinates();
      pt1_v = edge->point(1)->coordinates();
      GfxDebug::draw_point(pt0_v, CUBIT_GREEN );
      GfxDebug::draw_point(pt1_v, CUBIT_RED );
      GfxDebug::flush();
      i = i;
    }
  }
  
  if ( interpOrder == 0 ) {
    closest_point = closest_pt_v;
    if (tangent_ptr) {
      CubitVector tangent;
      if (best_edge->edge_tangent( closest_point, tangent ) != CUBIT_SUCCESS) {
        return CUBIT_FAILURE;
      }
      if (curvSense == CUBIT_REVERSED)
        tangent = -tangent;
      *tangent_ptr = tangent;
    } 
    // evaluate the curvature if required
    if (curvature_ptr) 
    {
      if( myEdgeList.size() == 1 )
        (*curvature_ptr).set( 0, 0, 0); 
      else
      {
        CubitVector curvature;
        myEdgeList.move_to( best_edge );
        int index = myEdgeList.get_index();

        //"best_edge" could be last or first on curve 
        CubitFacetEdge* prev_edge = NULL; 
        myEdgeList.back();
        if( (index - 1) == myEdgeList.get_index() ) 
          prev_edge = myEdgeList.get();

        CubitFacetEdge* next_edge = NULL; 
        myEdgeList.step(2);
        if( (index + 1) == myEdgeList.get_index() ) 
          next_edge = myEdgeList.get();


        CubitFacetEdge *closest_edge;
        //determine which adjacent edge is closest to "best_point"
        if( prev_edge && next_edge )
        {
          CubitVector tmp_vec;
          double prev_dist, next_dist;
          tmp_vec = (prev_edge->point(0)->coordinates() + 
                     prev_edge->point(1)->coordinates() ) / 2;
          prev_dist = closest_point.distance_between( tmp_vec ); 

          tmp_vec = (next_edge->point(0)->coordinates() + 
                     next_edge->point(1)->coordinates() ) / 2;
          next_dist = closest_point.distance_between( tmp_vec ); 
        
          if( prev_dist < next_dist )
            closest_edge = prev_edge;
          else
            closest_edge = next_edge;
        }
        else if( prev_edge )
          closest_edge = prev_edge;
        else
          closest_edge = next_edge;

        if (best_edge->edge_curvature( closest_point, curvature, closest_edge ) != CUBIT_SUCCESS) {
          return CUBIT_FAILURE;
        }
        if (curvSense == CUBIT_REVERSED)
          curvature = -curvature;
        *curvature_ptr = curvature;
      }
    }
  } else if ( interpOrder == 4 ) { 
    double t_value;    
    double tangent[3];
    stat = project_to_bezier_edge(best_edge,this_point,
                                  closest_point,tangent,&t_value);
    if(curvSense == CUBIT_REVERSED){
      tangent[0]= -tangent[0];
      tangent[1]= -tangent[1];
      tangent[2]= -tangent[2];
    }
    
    if ( tangent_ptr ) {                                    
      tangent_ptr->x(tangent[0]);
      tangent_ptr->y(tangent[1]);
      tangent_ptr->z(tangent[2]);
    }

    if ( curvature_ptr ) {
    //  The idea here is to get the planes normal to two points on the curve
    //  on either side of the closest_point, and also to get the plane through 
    //  these three points.  The point of intersection of these three planes
    //  should be a good approximation of the center of curvature, from which
    //  the curvature vector is readily found.
    //  tan1 and tan2 are the normals to the two planes for the endpoints;
    //  norm3 is the normal for the plane through the three points. 
      CubitVector tan1, tan2, point_1,  point_2;
      CubitVector *tan1_ptr = &tan1;
      CubitVector *tan2_ptr = &tan2;
//       CubitVector *point_1_ptr = &point_1;
      CubitVector *point_2_ptr = &point_2;
      t_value -= 0.01;
      best_edge->evaluate_single_tangent(t_value,tan1_ptr);
      t_value += 0.02;
      best_edge->evaluate_single(t_value,tan2_ptr);
      best_edge->evaluate_single_tangent(t_value,point_2_ptr);
      tan1.normalize();
      tan2.normalize();
      //  Get the plane through the point and normal to the tangent for the two 
      //  points at the end.
      //  Get the plane for the three points.
      CubitVector arm1 = point_1 - closest_point;
      CubitVector arm2 = point_2 - closest_point;
      CubitVector norm3 = arm1*arm2;
      norm3.normalize();
      double denominator;
      
      denominator = tan1 % ( tan2 * norm3 );  

      if ( fabs(denominator) < GEOMETRY_RESABS ) { // planes are parallel; 
                                                    // curvature is infinite
        curvature_ptr->x(CUBIT_DBL_MAX);
        curvature_ptr->y(CUBIT_DBL_MAX);
        curvature_ptr->z(CUBIT_DBL_MAX);
      
      } else {
        CubitVector IntersectionPoint;
        IntersectionPoint = ( (point_1%tan1)*(norm3*tan2) + 
                              (closest_point%norm3)*(tan2*tan1) +
                              (point_2%tan2)*(tan1*norm3) )/denominator;
      
        curvature_ptr->x(closest_point.x() - IntersectionPoint.x());
        curvature_ptr->y(closest_point.y() - IntersectionPoint.y());
        curvature_ptr->z(closest_point.z() - IntersectionPoint.z());
        
      }
              
    }  
  } else {
    PRINT_ERROR("Error:  Only curves or order 0 or 4 are supported.\n");
    return CUBIT_FAILURE;    
  }
  

  
  if (param)
  {
    *param = u_on_facet_edge( best_edge, closest_point );
  }
   
  if (mydebug)
  {
    draw_location(closest_point, CUBIT_GREEN);
  }
  return stat;

}

//===========================================================================
//Function Name: project_to_bezier_edge
//
//Member Type:  PRIVATE
//Description:  projects a point onto a Bezier edge 
//===========================================================================
CubitStatus CurveFacetEvalTool::project_to_bezier_edge(CubitFacetEdge *edge,
                        CubitVector &point, CubitVector &projected_point,
                        double *tangent, double *tval)
{
CubitStatus stat;
int numEdge, numVert, edgeVert[2], numLocs;
double vert[6], edgeCtrlPts[9], xyz[3], xyzOnEdge[3];
CubitVector *edge_ctrl_pts;

  stat = CUBIT_SUCCESS;

  numEdge = 1;
  numVert = 2;
  numLocs = 1;
  edgeVert[0] = 0; edgeVert[1] = 1;
  vert[0] = edge->point(0)->x();
  vert[1] = edge->point(0)->y();
  vert[2] = edge->point(0)->z();
  vert[3] = edge->point(1)->x();
  vert[4] = edge->point(1)->y();
  vert[5] = edge->point(1)->z();
  
    edge_ctrl_pts = edge->control_points();
    edgeCtrlPts[0] = edge_ctrl_pts[0].x();  
    edgeCtrlPts[1] = edge_ctrl_pts[0].y();  
    edgeCtrlPts[2] = edge_ctrl_pts[0].z();  
    edgeCtrlPts[3] = edge_ctrl_pts[1].x();  
    edgeCtrlPts[4] = edge_ctrl_pts[1].y();  
    edgeCtrlPts[5] = edge_ctrl_pts[1].z();  
    edgeCtrlPts[6] = edge_ctrl_pts[2].x();  
    edgeCtrlPts[7] = edge_ctrl_pts[2].y();  
    edgeCtrlPts[8] = edge_ctrl_pts[2].z();  
  xyz[0] = point.x();
  xyz[1] = point.y(); 
  xyz[2] = point.z();
  projToBezierEdge( numEdge,numVert,edgeVert,vert,edgeCtrlPts,
                       numLocs,xyz,xyzOnEdge,tangent,tval );

  projected_point.x(xyzOnEdge[0]);
  projected_point.y(xyzOnEdge[1]);
  projected_point.z(xyzOnEdge[2]);
  
  return stat;
}
/*
//===========================================================================
//Function Name: project_to_facet_edge
//
//Member Type:  PRIVATE
//Description:  project the point to the facets defining this curve
//===========================================================================
CubitStatus CurveFacetEvalTool::project_to_facet_edge(CubitVector &this_point,
                                                      CubitVector &closest_point,
                                                      CubitVector *tangent_ptr,
                                                      CubitVector *curvature_ptr,
                                                      double *param,
                                                      int *outside )
{
  int trim = 0;
  int ncheck, ii, nincr=0;
  static int calls=0;
  static int nncheck=0;
  static int ntol=0;
  static int mydebug=0;
  CubitBoolean outside_facet, best_outside_facet;
  CubitVector boxmin, boxmax, p0, p1;
  CubitVector close_point, best_point, best_areacoord;
  CubitFacetEdge *best_edge, *edge;
  CubitFacet *best_facet, *facet;

  double tol = surfFacetEvalTool->compare_tol();
  double facet_tol = 1.0e-3 * surfFacetEvalTool->compare_tol();
  double mindist = 0.0e0;
  CubitVector pt_on_edge;
  CubitBoolean done = CUBIT_FALSE;
  while(!done) {

    // define a bounding box around the point

    CubitVector ptmin( this_point.x() - tol, 
                       this_point.y() - tol, 
                       this_point.z() - tol );
   
    CubitVector ptmax( this_point.x() + tol, 
                       this_point.y() + tol, 
                       this_point.z() + tol );

    mindist = CUBIT_DBL_MAX;
    ncheck = 0;
    best_outside_facet = CUBIT_TRUE;
    myEdgeList.reset();
    for ( ii = myEdgeList.size(); ii > 0 && !done; ii-- ) {
	    edge = myEdgeList.get_and_step();
      p0 = edge->point( 0 )->coordinates();
      p1 = edge->point( 1 )->coordinates();

      // Try to trivially reject this facet with a bounding box test

      boxmin.x( CUBIT_MIN( p0.x(), p1.x() ) );
      boxmax.x( CUBIT_MAX( p0.x(), p1.x() ) );
      if (ptmax.x() < boxmin.x() ||
                      ptmin.x() > boxmax.x()) {
        continue;
      }
      boxmin.y( CUBIT_MIN( p0.y(), p1.y() ) );
      boxmax.y( CUBIT_MAX( p0.y(), p1.y() ) );
      if (ptmax.y() < boxmin.y() ||
		      ptmin.y() > boxmax.y()) {
        continue;
      }
      boxmin.z( CUBIT_MIN( p0.z(), p1.z() ) );
      boxmax.z( CUBIT_MAX( p0.z(), p1.z() ) );
      if (ptmax.z() < boxmin.z() ||
		      ptmin.z() > boxmax.z()) {
        continue;
      }

      // Only facets that pass the bounding box test will get past here!

      // Project point to plane of the facet and determine its area coordinates

      ncheck++;
      CubitVector pt_on_plane;
      double dist_to_plane;
      facet = edge->adj_facet( 0 );
      if (!facet) {
        facet = edge->adj_facet( 1 );
      }
      project_to_edge_line( edge, this_point, pt_on_plane, dist_to_plane );
      CubitVector areacoord;
      pt_on_edge = pt_on_plane;
      surfFacetEvalTool->facet_area_coordinate( facet, pt_on_plane, areacoord );

      // If sign of areacoords are all positive then its inside the triangle
      // and we are done - go interpolate the point. (use an absolute
      // tolerance since the coordinates arenormalized)

      if (areacoord.x() > -GEOMETRY_RESABS && 
          areacoord.y() > -GEOMETRY_RESABS && 
          areacoord.z() > -GEOMETRY_RESABS) {
        if (interpOrder == 0 && dist_to_plane < facet_tol) {
          outside_facet = CUBIT_FALSE;
          close_point = this_point;
        }
        else {
          if (surfFacetEvalTool->eval_facet( facet, this_point, areacoord, 
                          close_point, outside_facet ) 
            != CUBIT_SUCCESS) {
            return CUBIT_FAILURE;
          }
        }
      }

      // otherwise find the closest vertex or edge to the projected point

      else if (areacoord.x() < GEOMETRY_RESABS) 
      {
        outside_facet = CUBIT_TRUE;
        if (areacoord.y() < GEOMETRY_RESABS) 
        {
          if (surfFacetEvalTool->eval_point( facet, 2, close_point ) 
            != CUBIT_SUCCESS) {
            return CUBIT_FAILURE;
          }
        }
        else if(areacoord.z() < GEOMETRY_RESABS) 
        {
          if (surfFacetEvalTool->eval_point( facet, 1, close_point ) 
            != CUBIT_SUCCESS) {
            return CUBIT_FAILURE;
          }
        }
        else 
        {
          if (surfFacetEvalTool->project_to_facetedge( facet, 1, 2, this_point, pt_on_plane, 
                         close_point, outside_facet ) !=CUBIT_SUCCESS) {
            return CUBIT_FAILURE;
          }
        }
      }
      else if (areacoord.y() < GEOMETRY_RESABS) 
      {
        outside_facet = CUBIT_TRUE;
        if (areacoord.z() < GEOMETRY_RESABS) 
        {
          if (surfFacetEvalTool->eval_point( facet, 0, close_point ) 
            != CUBIT_SUCCESS) {
            return CUBIT_FAILURE;
          }
        }
        else 
        {
          if (surfFacetEvalTool->project_to_facetedge( facet, 2, 0, this_point, pt_on_plane, 
                         close_point, outside_facet ) !=CUBIT_SUCCESS) {
            return CUBIT_FAILURE;
          }
        }
      }
      else 
      {
        outside_facet = CUBIT_TRUE;
        if (surfFacetEvalTool->project_to_facetedge( facet, 0, 1, this_point, pt_on_plane, 
                       close_point, outside_facet ) !=CUBIT_SUCCESS) {
          return CUBIT_FAILURE;
        }
      }
      
      // keep track of the minimum distance

      double dist = sqrt(sqr(close_point.x() - this_point.x()) +
                         sqr(close_point.y() - this_point.y()) +
                         sqr(close_point.z() - this_point.z()));
      if (dist < mindist) {
        if ((best_outside_facet == CUBIT_FALSE && outside_facet == CUBIT_TRUE)) {
          //int x=1;
        }
        else {
          mindist = dist;
          best_point = close_point;
          best_edge = edge;
          best_facet = facet;
          best_areacoord = areacoord;
          best_outside_facet = outside_facet;

          if (dist < facet_tol) {
            done = CUBIT_TRUE;
          }
        }
      }
    }

    // We are done if we found at least one triangle.  Otherwise
    // increase the tolerance and try again

    nincr++;
    if (ncheck > 0) {
      if (best_outside_facet && nincr < 10) {
        tol *= 2.0;
        ntol++;
      }
      else {
        done = CUBIT_TRUE;
      }
    }
    else {
      tol *= 2.0e0;
      ntol++;
    }
  }


  //If area coords are still negative...we are still not on a "best_facet" 
  //Find closest point to "best_edge" and use this point for area coordinates
  if (best_areacoord.x() <= -GEOMETRY_RESABS || 
      best_areacoord.y() <= -GEOMETRY_RESABS || 
      best_areacoord.z() <= -GEOMETRY_RESABS) 
  {
    surfFacetEvalTool->facet_area_coordinate( best_facet, best_point, best_areacoord );
  }


  // if the closest point is outside of a facet, then evaluate the point
  // on the facet using its area coordinates (otherwise it would be 
  // trimmed to an edge or point)
  if ( !trim && best_outside_facet && interpOrder != 4) {
    if (surfFacetEvalTool->eval_facet( best_facet, this_point, best_areacoord, 
      best_point, best_outside_facet ) 
      != CUBIT_SUCCESS) {
      return CUBIT_FAILURE;
    }
    
    // see if its really outside (it could just be on an edge where the
    // curvature is convex)

    best_outside_facet = surfFacetEvalTool->is_outside( best_facet, best_areacoord );
  }

  // evaluate the tangent if required

  CubitVector tangent_vec;
  if (tangent_ptr || curvature_ptr) {
    CubitVector facet_tangent = best_edge->point(1)->coordinates() - 
                                best_edge->point(0)->coordinates();
    if (curvSense == CUBIT_REVERSED)
      facet_tangent = -facet_tangent;
    facet_tangent.normalize();
    if (surfFacetEvalTool->interp_order() == 0)
    {
      tangent_vec = facet_tangent;
    }
    else
    {
      CubitVector normal;
      if (surfFacetEvalTool->eval_facet_normal( best_facet, best_areacoord, normal ) 
        != CUBIT_SUCCESS) {
        return CUBIT_FAILURE;
      }    
      tangent_vec = normal * (facet_tangent * normal);
    }
  }

  if( tangent_ptr )
    *tangent_ptr = tangent_vec;

  if (curvature_ptr) {
    //get adjacent facet edges to "best_edge"
    myEdgeList.move_to( best_edge );
    int index = myEdgeList.get_index();

    //"best_edge" could be last or first on curve 
    CubitFacetEdge* prev_edge = NULL; 
    myEdgeList.back();
    if( (index - 1) == myEdgeList.get_index() ) 
      prev_edge = myEdgeList.get();

    CubitFacetEdge* next_edge = NULL; 
    myEdgeList.step(2);
    if( (index + 1) == myEdgeList.get_index() ) 
      next_edge = myEdgeList.get();

    //now get 1 facet on each adj. edge
    CubitFacet *prev_facet = NULL;
    if( prev_edge )
    {
      prev_facet = prev_edge->adj_facet( 0 );
      if (!prev_facet)
        prev_facet = prev_edge->adj_facet( 1 );
    }

    CubitFacet *next_facet = NULL;
    if( next_edge )
    {
      next_facet = next_edge->adj_facet( 0 );
      if (!next_facet)
        next_facet = next_edge->adj_facet( 1 );
    }

    //get tangent at midpoint of prev_edge
    CubitVector prev_tangent_vec;
    if( prev_edge ) 
    {
      CubitVector facet_tangent = prev_edge->point(1)->coordinates() - 
                                  prev_edge->point(0)->coordinates();
      if (curvSense == CUBIT_REVERSED)
        facet_tangent = -facet_tangent;
      facet_tangent.normalize();

      if (surfFacetEvalTool->interp_order() == 0)
        prev_tangent_vec = facet_tangent;
      else
      {
        CubitVector prev_mid_point;
        prev_mid_point = (prev_edge->point(0)->coordinates() + 
                          prev_edge->point(1)->coordinates() ) / 2;

        CubitVector areacoord;
        surfFacetEvalTool->facet_area_coordinate( prev_facet, prev_mid_point, areacoord );
        CubitVector normal;
        if (surfFacetEvalTool->eval_facet_normal( prev_facet, areacoord, normal ) 
          != CUBIT_SUCCESS) 
        {
          return CUBIT_FAILURE;
        }    
        prev_tangent_vec = normal * (facet_tangent * normal);
      }
    }

    //get tangent at midpoint of next_edge
    CubitVector next_tangent_vec;
    if( next_edge ) 
    {
      CubitVector facet_tangent = next_edge->point(1)->coordinates() - 
                                  next_edge->point(0)->coordinates();
      if (curvSense == CUBIT_REVERSED)
        facet_tangent = -facet_tangent;
      facet_tangent.normalize();

      if (surfFacetEvalTool->interp_order() == 0)
        next_tangent_vec = facet_tangent;
      else
      {
        CubitVector next_mid_point; 
        next_mid_point = (next_edge->point(0)->coordinates() + 
                          next_edge->point(1)->coordinates() ) / 2;

        CubitVector areacoord;
        surfFacetEvalTool->facet_area_coordinate( next_facet, next_mid_point, areacoord );
        CubitVector normal;
        if (surfFacetEvalTool->eval_facet_normal( next_facet, areacoord, normal ) 
          != CUBIT_SUCCESS) 
        {
          return CUBIT_FAILURE;
        }    
        next_tangent_vec = normal * (facet_tangent * normal);
      }
    }

    //decide how to weight the tangents
    double prev_weight, next_weight; 
    if( !prev_edge ) 
    {
      prev_weight = 0;
      next_weight = 2;
    }
    else if( !next_edge )
    {
      prev_weight = 2;
      next_weight = 0;
    }
    else
    {
      double dist1 = best_edge->point(0)->coordinates().distance_between(
                     best_edge->point(1)->coordinates() );
      double dist2 = best_edge->point(0)->coordinates().distance_between(
                     best_point);
      next_weight = dist2/dist1;
      next_weight *= 2; 
      prev_weight = 2-next_weight; 
    }

    *curvature_ptr = (prev_weight*( tangent_vec - prev_tangent_vec ) + 
                     next_weight*( next_tangent_vec - tangent_vec ) );
  }

  closest_point = best_point;
 
  *outside = best_outside_facet;

  if (mydebug) {
    nncheck+= ncheck;
    calls++;
    if (calls%100==0){
      char message[100];
      sprintf(message,"calls = %d, ckecks = %d, ntol = %d\n",calls,nncheck,ntol);
      PRINT_INFO(message);
    }
  }

  // determine the parameter on the curve

  if (param)
  {
    *param = u_on_facet_edge( best_edge, closest_point );
  }
  
  return CUBIT_SUCCESS;
}

//===========================================================================
//Function Name: project_to_edge_line
//
//Member Type:  PRIVATE
//Description:  compute the area coordinate on the linear facet edge
//===========================================================================
CubitStatus CurveFacetEvalTool::project_to_edge_line( CubitFacetEdge *edge, 
                                                      CubitVector &this_point, 
                                                      CubitVector &pt_on_edge, 
                                                      double &dist_to_edge )
{
  CubitVector p0 = edge->point(0)->coordinates();
  CubitVector p1 = edge->point(1)->coordinates();
  CubitVector v0 = p1 - p0;
  CubitVector v1 = this_point - p0;
  double h = v1.normalize();
  v0.normalize();
  double costheta = v0 % v1;
  double l = h * costheta;
  dist_to_edge = sqrt( fabs((h*h) - (l*l)) );
  pt_on_edge = p0 + l*v0;
  
  return CUBIT_SUCCESS;
}

//===========================================================================
//Function Name: project_to_linear_facet_edge
//
//Member Type:  PRIVATE
//Description:  project the point to the facets defining this curve.
//              Note: this function assumes a piecewise linear representation
//              of the curve.  Does not rely on an adjacent surface
//              eval tool to project the point
//===========================================================================
CubitStatus CurveFacetEvalTool::project_to_linear_facet_edge(CubitVector &this_point,
                                                      CubitVector &closest_point,
                                                      CubitVector *tangent_ptr,
                                                      CubitVector *curvature_ptr,
                                                      double *param,
                                                      int *outside )
{
  int trim = 0;
  int ncheck, ii, nincr=0;
  static int ntol=0;
  CubitBoolean outside_facet, best_outside_facet;
  CubitVector boxmin, boxmax, p0, p1;
  CubitVector close_point, best_point;
  CubitFacetEdge *best_edge, *edge;

  // so we don't evaluate an edge more than once - mark the edges
  // as we evaluate them.  Put the evaluated edges on a used_edge_list
  // so we clear the marks off when we are done.  Note: this assumes
  // theat marks are initially cleared.
  
  DLIList<CubitFacetEdge *>used_edge_list;
  for(ii=0; ii<myEdgeList.size(); ii++)
    myEdgeList.get_and_step()->set_flag(0);
  CubitBoolean eval_all = CUBIT_FALSE;

  double tol = facetLength * 1.0e-3;
  double facet_tol = facetLength * GEOMETRY_RESABS;
  double mindist = CUBIT_DBL_MAX;
  double dist = 0.0e0;
  CubitBoolean done = CUBIT_FALSE;
  while(!done) {

    // define a bounding box around the point

    CubitVector ptmin( this_point.x() - tol, 
                       this_point.y() - tol, 
                       this_point.z() - tol );

    CubitVector ptmax( this_point.x() + tol, 
		       this_point.y() + tol, 
		       this_point.z() + tol );

    ncheck = 0;
    best_outside_facet = CUBIT_TRUE;
    myEdgeList.reset();
    for ( ii = myEdgeList.size(); ii > 0 && !done; ii-- ) 
    {
      edge = myEdgeList.get_and_step();
      
//      if( (edge->get_flag() == 1) ) 
//        continue;

      p0 = edge->point( 0 )->coordinates();
      p1 = edge->point( 1 )->coordinates();

      // Try to trivially reject this facet with a bounding box test

      if (!eval_all)
      {
        boxmin.x( CUBIT_MIN( p0.x(), p1.x() ) );
        boxmax.x( CUBIT_MAX( p0.x(), p1.x() ) );
	      if (ptmax.x() < boxmin.x() ||
		        ptmin.x() > boxmax.x()) {
          continue;
        }
        boxmin.y( CUBIT_MIN( p0.y(), p1.y() ) );
        boxmax.y( CUBIT_MAX( p0.y(), p1.y() ) );
        if (ptmax.y() < boxmin.y() ||
		        ptmin.y() > boxmax.y()) {
          continue;
        }
        boxmax.z( CUBIT_MAX( p0.z(), p1.z() ) );
        if (ptmax.z() < boxmin.z() ||
		        ptmin.z() > boxmax.z()) {
          continue;
        }
      }

      // Only edges that pass the bounding box test will get past here!

      // find distance to this edge
      dist = edge->dist_to_edge( this_point, close_point, outside_facet );

      if (dist <= mindist) {
        if ((best_outside_facet == CUBIT_FALSE && outside_facet == CUBIT_TRUE)) {
          //int x=1;
        }
        else {
          mindist = dist;
          best_point = close_point;
          best_edge = edge;
          best_outside_facet = outside_facet;

          if (dist < facet_tol) {
            done = CUBIT_TRUE;
          }
        }
      }
      ncheck++;
//      edge->set_flag(1);
//      used_edge_list.append(edge);
    }

    // We are done if we found at least one edge.  Otherwise
    // increase the tolerance and try again

    nincr++;
    if (ncheck > 0) {
      done = CUBIT_TRUE;
      if (best_outside_facet && nincr < 10) {
        done = CUBIT_FALSE;
      }
    }

    if (!done)
    {
      if (nincr < 10)
      {
        tol *= 2.0e0;
        ntol++;
      }
      else {
        eval_all = CUBIT_TRUE;
      }
    }
  }

  assert(best_edge != NULL);

  // if the closest point is outside of a facet, then evaluate the point
  // on the facet using its area coordinates (otherwise it would be 
  // trimmed to an edge or point)


  if ( !trim && best_outside_facet) {
//    if (best_edge->proj_to_line( this_point, best_point )!= CUBIT_SUCCESS) {
    if (best_edge->closest_point( this_point, best_point )!= CUBIT_SUCCESS) {
      return CUBIT_FAILURE;
    }
  }


  // evaluate the tangent if required
 
  if (tangent_ptr) {
    //get 2 adjacent edges
    CubitVector tangent;
    if (best_edge->edge_tangent( best_point, tangent ) 
      != CUBIT_SUCCESS) {
      return CUBIT_FAILURE;
    }
    if (curvSense == CUBIT_REVERSED)
      tangent = -tangent;
    *tangent_ptr = tangent;
  }

  // evaluate the curvature if required
  if (curvature_ptr) {
    CubitVector curvature;
    myEdgeList.move_to( best_edge );
    int index = myEdgeList.get_index();

    //"best_edge" could be last or first on curve 
    CubitFacetEdge* prev_edge = NULL; 
    myEdgeList.back();
    if( (index - 1) == myEdgeList.get_index() ) 
      prev_edge = myEdgeList.get();

    CubitFacetEdge* next_edge = NULL; 
    myEdgeList.step(2);
    if( (index + 1) == myEdgeList.get_index() ) 
      next_edge = myEdgeList.get();


    CubitFacetEdge *closest_edge;
    //determine which adjacent edge is closest to "best_point"
    if( prev_edge && next_edge )
    {
      CubitVector tmp_vec;
      double prev_dist, next_dist;
      tmp_vec = (prev_edge->point(0)->coordinates() + 
                 prev_edge->point(1)->coordinates() ) / 2;
      prev_dist = best_point.distance_between( tmp_vec ); 

      tmp_vec = (next_edge->point(0)->coordinates() + 
                 next_edge->point(1)->coordinates() ) / 2;
      next_dist = best_point.distance_between( tmp_vec ); 
      
      if( prev_dist < next_dist )
        closest_edge = prev_edge;
      else
        closest_edge = next_edge;
    }
    else if( prev_edge )
      closest_edge = prev_edge;
    else
      closest_edge = next_edge;

    if (best_edge->edge_curvature( best_point, curvature, closest_edge ) != CUBIT_SUCCESS) {
      return CUBIT_FAILURE;
    }
    if (curvSense == CUBIT_REVERSED)
      curvature = -curvature;
    *curvature_ptr = curvature;
  }

  closest_point = best_point;
  *outside = best_outside_facet;

  if (param)
  {
    *param = u_on_facet_edge( best_edge, closest_point );
  }

  // clear the marks from the used edges

  for (ii=0; ii<used_edge_list.size(); ii++)
  {
    edge = used_edge_list.get_and_step();
    edge->set_flag( 0 );
  }

  
  return CUBIT_SUCCESS;
}

*/
//===========================================================================
//Function Name: u_on_facet_edge
//
//Member Type:  PRIVATE
//Description:  return the u param on the facet curve given a point on the
//              curve and the facet edge it is on
//===========================================================================
double CurveFacetEvalTool::u_on_facet_edge( CubitFacetEdge *edge_at_pt, 
                                            CubitVector pt )
{
  int ii;
  double cum_len = 0.0;
  CubitBoolean done = CUBIT_FALSE;
  myEdgeList.reset();
  for (ii=0; ii<myEdgeList.size() && !done; ii++) {
    CubitFacetEdge *edge = myEdgeList.get_and_step();
    if (edge != edge_at_pt)
    {
      cum_len += edge->length();
    }
    else
    {
      CubitVector pt0;
      if (curvSense == CUBIT_REVERSED)
      {
        pt0 = edge->point( 1 )->coordinates();
      }
      else
      {
        pt0 = edge->point( 0 )->coordinates();
      }
      cum_len += pt.distance_between( pt0 );
      done  = CUBIT_TRUE;
    }
  }
  if (!done)
  {
    PRINT_DEBUG_122( "Error in CurveFacetEvalTool::u_on_facet_edge" );
  }
  double u = cum_len / facetLength;
  return u;
}

//===========================================================================
//Function Name: destroy_facets
//
//Member Type:  PRIVATE
//Description:  Deletes the points and facets.
//===========================================================================
void CurveFacetEvalTool::destroy_facets()
{

}
  

//===========================================================================
//Function Name: draw_edges
//
//Member Type:  PRIVATE
//Description:  draw the facet edges 
//===========================================================================
void CurveFacetEvalTool::draw_edges(int color)
{
  int ii;
  if ( color == -1 )
    color = CUBIT_BLUE;
  for ( ii = myEdgeList.size(); ii > 0; ii-- )
  {
    CubitFacetEdge *edge = myEdgeList.get_and_step();
    CubitPoint *begin_point = edge->point(0);
    CubitPoint *end_point = edge->point(1);
    GfxDebug::draw_line(begin_point->x(), 
                        begin_point->y(), 
                        begin_point->z(),
                        end_point->x(), 
                        end_point->y(), 
                        end_point->z(), 
                        color);
  }
  GfxDebug::flush();
}

//===========================================================================
//Function Name: draw_edge
//
//Member Type:  PRIVATE
//Description:  draw the facet edge 
//===========================================================================
void CurveFacetEvalTool::draw_edge(CubitFacetEdge *edge, int color)
{

  CubitPoint *begin_point = edge->point(0);
  CubitPoint *end_point = edge->point(1);
  GfxDebug::draw_line(begin_point->x(), 
                      begin_point->y(), 
                      begin_point->z(),
                      end_point->x(), 
                      end_point->y(), 
                      end_point->z(), 
                      color);

  GfxDebug::flush();
}

//===========================================================================
//Function Name: draw_line
//
//Member Type:  PRIVATE 
//===========================================================================
void CurveFacetEvalTool::draw_line(CubitVector &begin, CubitVector &end, int color)
{
  GfxDebug::draw_line(begin.x(), begin.y(), begin.z(),
                      end.x(), end.y(), end.z(), color);
  GfxDebug::flush();
}

//===========================================================================
//Function Name: draw_location
//
//Member Type:  PRIVATE 
//===========================================================================
void CurveFacetEvalTool::draw_location(CubitVector &loc, int color )
{
  if ( color == -1 )
    color = CUBIT_YELLOW;
  GfxDebug::draw_point(loc, color);
  GfxDebug::flush();
}

//===========================================================================
//Function Name: set_length
//
//Member Type:  PRIVATE
//Description:  compute the length of the facets in this curve
// assumption:  facet edges have been set up
//===========================================================================
void CurveFacetEvalTool::set_length()
{
  int ii;
  CubitFacetEdge *cf_edge;

  facetLength = 0.0e0;

  for (ii=0; ii<myEdgeList.size(); ii++)
  {
    cf_edge = myEdgeList.get_and_step();
    facetLength += cf_edge->length();
  }

}

//===========================================================================
//Function Name: save
//
//Member Type:  PUBLIC
//Description:  save the curve facet eval tool to a cubit file  
// Assumption:  contained edge facets have been previuosly saved.  This function
//              saves only the edge facet ids.
//===========================================================================
CubitStatus CurveFacetEvalTool::save( 
  FILE *fp )
{
  NCubitFile::CIOWrapper cio(fp);
  typedef NCubitFile::UnsignedInt32 int32;

    // write out "interpOrder" 
  cio.Write(reinterpret_cast<int32*>(&interpOrder), 1);

   // write the associated facet eval tool id.  If there is none then write -1 
  int surf_tool_id = -1;
  if (surfFacetEvalTool != NULL)
    surf_tool_id = surfFacetEvalTool->get_output_id();
  cio.Write(reinterpret_cast<int32*>(&surf_tool_id), 1);

    // convert "curvSense" in an int
  int sense;
  if( curvSense == CUBIT_UNKNOWN )
    sense = -1;
  else 
    sense = (curvSense == CUBIT_REVERSED) ? 1 : 0; 

    // write "curveSense" and "goodCurveData"
  cio.Write(reinterpret_cast<int32*>(&sense), 1);
  int32 is_good = goodCurveData ? 1 : 0;
  cio.Write(&is_good, 1);
 
    // write "facetLength"
  cio.Write( &facetLength, 1 );

    // write ids of facet edges in "myEdgeList"
  int ii;
  CubitFacetEdge *edge_ptr;
  int nedges = myEdgeList.size();
  int32* edge_id = new int32 [nedges];
  myEdgeList.reset();
  for (ii=0; ii<nedges; ii++)
  {
    edge_ptr = myEdgeList.get_and_step();
    edge_id[ii] = edge_ptr->id();
  }
  cio.Write(reinterpret_cast<int32*>(&nedges), 1);
  if (nedges > 0)
  {
    cio.Write(edge_id, nedges);
  }
  delete [] edge_id;

    // write ids of points in "myPointList" 
  CubitPoint *point_ptr;
  int npoints = myPointList.size();
  int32* point_id = new int32 [npoints];
  myPointList.reset();
  for (ii=0; ii<npoints; ii++)
  {
    point_ptr = myPointList.get_and_step();
    point_id[ii] = point_ptr->id();
  }
  cio.Write(reinterpret_cast<int32*>(&npoints), 1);
  if (npoints > 0)
  {
    cio.Write(point_id, npoints);
  }
  delete [] point_id;

  return CUBIT_SUCCESS;
}

//===========================================================================
//Function Name: restore
//Member Type:  PUBLIC
//Description:  restore curve eval tool from a CUB file  
//author: sjowen
//Date:1/28/2003
//===========================================================================
CubitStatus CurveFacetEvalTool::restore(FILE *fp,
                      unsigned int endian,
                      int num_edges, 
                      int num_points,
                      CubitFacetEdge **edges, 
                      CubitPoint **points,
                      int num_fets,
                      FacetEvalTool **fet_array)
{
  NCubitFile::CIOWrapper cio(endian, fp);
  typedef NCubitFile::UnsignedInt32 int32;

    // read stuff about this eval tool

  int ii;
  int int_data[4];
  cio.Read(reinterpret_cast<int32*>(int_data), 4);
  interpOrder = int_data[0];
  int surf_tool_id = int_data[1]; 

  if (int_data[2] == -1 )  
    curvSense = CUBIT_UNKNOWN;
  else
    curvSense= int_data[2] ? CUBIT_REVERSED : CUBIT_FORWARD;

  goodCurveData = (int_data[3]==0) ? CUBIT_FALSE : CUBIT_TRUE;

  if( surf_tool_id != -1 )
    surfFacetEvalTool = fet_array[ surf_tool_id ]; 
  else
    surfFacetEvalTool = NULL; 

  cio.Read(&facetLength, 1);

    // read the edges

  int nedges; 
  cio.Read(reinterpret_cast<int32*>(&nedges), 1);
  if (nedges > 0)
  {
    int32* edge_id = new int32 [nedges];
    cio.Read(edge_id, nedges);
    int id;
    for (ii=0; ii<nedges; ii++)
    {
      id = edge_id[ii];
      if (id <0 || id >= num_edges)
      {
        delete [] edge_id;
        return CUBIT_FAILURE;
      }
      myEdgeList.append(edges[id]);
    }
    delete [] edge_id;
  }
  
    // read the points

  int npoints; 
  cio.Read(reinterpret_cast<int32*>(&npoints), 1);
  int id;
  if (npoints > 0)
  {
    int32* point_id = new int32 [npoints];
    cio.Read(point_id, npoints);
    for (ii=0; ii<npoints; ii++)
    {
      id = point_id[ii];
      if (id <0 || id >= num_points)
      {
        delete [] point_id;
        return CUBIT_FAILURE;
      }
      myPointList.append(points[id]);
    }
    delete [] point_id;
  }

  bounding_box();

  return CUBIT_SUCCESS;
}

CubitStatus CurveFacetEvalTool::fix_point_edge_order()
{
  CubitFacetEdge* this_edge;
  CubitFacetEdge* next_edge;
  CubitPoint *this_pt1;
  CubitPoint *next_pt0, *next_pt1;
  
  if( myEdgeList.size() == 0 )
    return CUBIT_FAILURE;

  this_edge = myEdgeList.get_and_step();
  int ii;
  for( ii = myEdgeList.size() - 1; ii > 0; ii-- )
  {
    if( 0 != this_edge->num_adj_facets() )
       continue;
    
    next_edge = myEdgeList.get_and_step();

    this_pt1 = this_edge->point( 1 );
    next_pt0 = next_edge->point( 0 );
    next_pt1 = next_edge->point( 1 );

    if( this_pt1 != next_pt0 &&
        this_pt1 != next_pt1 )
    {
        //Swap direction of edge.  (mod. 3-7-06)  Now calling flip instead
        // of doing the flip manually.  The flip function also tracks
        // the orientation so that we can remember the original orientation
        // of the edge.
      this_edge->flip();
      
    }
    
    this_edge = next_edge;
  }
    //Handle last edge
  next_edge = myEdgeList.prev(2);
  if( 0 == this_edge->num_adj_facets() )
  {
    CubitPoint* this_pt0 = this_edge->point( 0 );
    next_pt0 = next_edge->point( 0 );
    next_pt1 = next_edge->point( 1 );

    if( this_pt0 != next_pt0 &&
        this_pt0 != next_pt1 )
    {
        //Swap direction of edge.(mod. 3-7-06)  Now calling flip instead
        // of doing the flip manually.  The flip function also tracks
        // the orientation so that we can remember the original orientation
        // of the edge.
      this_edge->flip();
      
    }
  }

  return CUBIT_SUCCESS;
}
//===========================================================================
//Function Name: debug_draw_facet_edges
//
//Member Type:  PUBLIC
//Descriptoin:  draw the facet edges for debug 
//===========================================================================
void CurveFacetEvalTool::debug_draw_facet_edges( int color )
{
  draw_edges(color);
}

CubitBoolean CurveFacetEvalTool::replace_point( CubitPoint *del_pnt, CubitPoint *keep_pnt )
{
  CubitPoint *point_ptr;
  int npoints = myPointList.size();
  myPointList.reset();
  CubitBoolean istat = CUBIT_FALSE;
  int i;
  for ( i = 0; i < npoints; i++ )
  {
    point_ptr = myPointList.get();
    
    if( point_ptr == del_pnt )
    {
      myPointList.remove();
      myPointList.insert( keep_pnt );
      istat = CUBIT_TRUE;
    }
    myPointList.step();
  }
  
  return istat;
}


CubitBoolean CurveFacetEvalTool::replace_facets( DLIList< CubitFacetEdge *> &curv_edges )
{

  // replace edges 
  this->myEdgeList = curv_edges;

  // replace points
  DLIList<CubitPoint *> point_list;
  int i;
  // insert start point of every facet_edge
  curv_edges.reset();
  for( i = 0; i < curv_edges.size(); i++ )
  {
    point_list.append( CAST_TO( curv_edges.get_and_step(), CubitFacetEdge )->point(0) );
  }
  // insert end point of last facet_edge
  curv_edges.step( curv_edges.size() - 1 );
  point_list.append( CAST_TO( curv_edges.get(), CubitFacetEdge )->point(1) );
  this->myPointList = point_list;

  return CUBIT_TRUE;
}


void CurveFacetEvalTool::remove_facets( DLIList<CubitFacetEdge*> &facet_edges)
{
  facet_edges = myEdgeList;
  myEdgeList.clean_out();
  myPointList.clean_out();
}
