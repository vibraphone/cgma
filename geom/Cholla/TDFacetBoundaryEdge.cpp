//- Class:       TDFacetBoundaryEdge
//- Description: Tool data for storing additional information at 
//-              the boundary of a facet set
//- Owner:       Steve Owen
//- Checked by:
//- Version:

#include "TDFacetBoundaryEdge.hpp"
#include "TDFacetBoundaryPoint.hpp"
#include "CubitFacetEdge.hpp"
#include "CastTo.hpp"
#include "CubitFacet.hpp"
#include "CubitPoint.hpp"
#include "FacetEvalTool.hpp"
#include "debug.hpp"
TDFacetBoundaryEdge::TDFacetBoundaryEdge()
{
  edgePtr = NULL;
  nextEdgePtr = NULL;
  prevEdgePtr = NULL;
}

TDFacetBoundaryEdge::~TDFacetBoundaryEdge()
{
  int ii;
  for (ii=0; ii<edgeDataList.size(); ii++)
  {
    BoundaryEdgeData *bed_ptr = edgeDataList.get_and_step();
    delete bed_ptr;
  }
}


//-------------------------------------------------------------------------
// Purpose       : create a new facet boundary edge
//
// Special Notes :
//
// Creator       : Steve Owen
//
// Creation Date : 05/01
//------------------------------------------------------------------------- 
CubitStatus TDFacetBoundaryEdge::add_facet_boundary_edge( 
  CubitFacetEdge *edge_ptr )
{
  ToolData *td;
  td = edge_ptr->get_TD(&TDFacetBoundaryEdge::is_facet_boundary_edge);
  if ( td == NULL )
  {
    TDFacetBoundaryEdge *td_gm = new TDFacetBoundaryEdge;
    edge_ptr->add_TD( td_gm);
    td_gm->set_edge( edge_ptr );
  }
  else
  {
    TDFacetBoundaryEdge *td_gm = CAST_TO(td, TDFacetBoundaryEdge);
    td_gm->set_edge( edge_ptr );
  }
  return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Purpose       : get the facet boundary edge from an edge
//
// Special Notes :
//
// Creator       : Steve Owen
//
// Creation Date : 05/01
//------------------------------------------------------------------------- 
TDFacetBoundaryEdge* TDFacetBoundaryEdge::get_facet_boundary_edge( 
  CubitFacetEdge *edge_ptr )
{
  ToolData *td;
  td = edge_ptr->get_TD(&TDFacetBoundaryEdge::is_facet_boundary_edge);
  if ( td != NULL )
  {
    TDFacetBoundaryEdge *td_gm = CAST_TO(td, TDFacetBoundaryEdge);
    return td_gm;
  }
  return (TDFacetBoundaryEdge*) NULL;
}

//-------------------------------------------------------------------------
// Purpose       : add facet representing a new surface adjacent this edge
//
// Special Notes :  pass in a list of facets and pick the ones that have
//                  this edge as one of its edges
//
// Creator       : Steve Owen
//
// Creation Date : 05/01
//------------------------------------------------------------------------- 
void TDFacetBoundaryEdge::add_surf_facet( DLIList<CubitFacet *> &facet_list)
{
  CubitFacet *facet_ptr;
  int found = 0;
  int ii,jj;
  for (ii=0; ii< facet_list.size(); ii++)
  {
    facet_ptr = facet_list.get_and_step();
    if (facet_ptr->edge( 0 ) == edgePtr ||
        facet_ptr->edge( 1 ) == edgePtr ||
        facet_ptr->edge( 2 ) == edgePtr)
    {
      // check if has already been added

      found = 0;
      BoundaryEdgeData *bed_ptr;
      for (jj=0; jj<edgeDataList.size() && !found; jj++)
      {
        bed_ptr = edgeDataList.get_and_step();
        if(bed_ptr->adjFacet == facet_ptr)
        {
          found = 1;
        }
      }
      if (!found)
      {
        // add a new bed and facet

        BoundaryEdgeData *new_bed_ptr = new BoundaryEdgeData;
        new_bed_ptr->adjFacet = facet_ptr;
        new_bed_ptr->surfID = -1;
        edgeDataList.append( new_bed_ptr );      
      }
    }
  }
}

//======================================================================
// Function: control_points (PUBLIC)
// Description: set the Bezier control points on the edge.  The first 
//              and last control points are assumed to be the end points 
//              of the edge, so they are not passed into this function.  
//              order-1 points should be passed in.  The order of the 
//              bezier is returned.
// Author: sjowen
// Date: 8/00
//======================================================================
void TDFacetBoundaryEdge::control_points(
  CubitVector *ctrl_pts, int order, int surf_id ) 
{
  assert(order<=4 && order >=0);
  int found = 0;
  BoundaryEdgeData *bed_ptr = NULL;
  for (int ii=0; ii<edgeDataList.size(); ii++)
  {
    bed_ptr = edgeDataList.get_and_step();
    if(bed_ptr->surfID == surf_id)
      found = 1;
  }
  assert( found );  // the surf_id does not match any in the list

  int jj;
  for(jj=0; jj<order-1; jj++){
    bed_ptr->bezierCtrlPts[jj] = ctrl_pts[jj];
  }
}

//======================================================================
// Function: control_points (PUBLIC)
// Description: return the Bezier control points (including the end points)
//        The order of the bezier is returned.
// Author: sjowen
// Date: 8/00
//======================================================================
int TDFacetBoundaryEdge::control_points( CubitVector *ctrl_pts, int surf_id ) 
{
  ctrl_pts[0] = edgePtr->point(0)->coordinates();

  int ii;
  int found = 0;
  BoundaryEdgeData *bed_ptr = NULL;
  for (ii=0; ii<edgeDataList.size(); ii++)
  {
    bed_ptr = edgeDataList.get_and_step();
    if(bed_ptr->surfID == surf_id)
      found = 1;
  }
  assert( found );  // the surf_id does not match any in the list

  for(int i=0; i<3; i++) {
    ctrl_pts[i+1] = bed_ptr->bezierCtrlPts[i];
  }
  ctrl_pts[4] = edgePtr->point(1)->coordinates();
  return 4;
}

//======================================================================
// Function: control_points (PUBLIC)
// Description: return the control points on the edge of a facet based
//        on its edge use direction
// Author: sjowen
// Date: 8/00
//======================================================================
CubitStatus TDFacetBoundaryEdge::control_points( 
  CubitFacet *facet, CubitVector *ctrl_pts )
{
  // find the edge on the facet

  int index = -1;
  CubitBoolean found = CUBIT_FALSE;
  for (int i=0; i<3 && !found; i++) {
    if (edgePtr == facet->edge(i)) {
      index = i;
      found = CUBIT_TRUE;
    }
  }
  if (!found) {
    return CUBIT_FAILURE;
  }

  // locate the facet in the boundary edge data

  found = CUBIT_FALSE;
  BoundaryEdgeData *bed_ptr = NULL;
  for (int ii=0; ii<edgeDataList.size(); ii++)
  {
    bed_ptr = edgeDataList.get_and_step();
    if(bed_ptr->adjFacet == facet)
      found = CUBIT_TRUE;
  }
  if (!found)  // the facet does not match any in the list
    return CUBIT_FAILURE;

  // retreive the control points

  switch (facet->edge_use(index)) {
  case 1:
    ctrl_pts[0] = edgePtr->point(0)->coordinates();
    ctrl_pts[1] = bed_ptr->bezierCtrlPts[0];
    ctrl_pts[2] = bed_ptr->bezierCtrlPts[1];
    ctrl_pts[3] = bed_ptr->bezierCtrlPts[2];
    ctrl_pts[4] = edgePtr->point(1)->coordinates();
    break;
  case -1:
    ctrl_pts[0] = edgePtr->point(1)->coordinates();
    ctrl_pts[1] = bed_ptr->bezierCtrlPts[2];
    ctrl_pts[2] = bed_ptr->bezierCtrlPts[1];
    ctrl_pts[3] = bed_ptr->bezierCtrlPts[0];
    ctrl_pts[4] = edgePtr->point(0)->coordinates();
    break;
  default:
    return CUBIT_FAILURE;
  }
  return CUBIT_SUCCESS;
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
void TDFacetBoundaryEdge::get_control_points( CubitPoint *point_ptr, 
                                             CubitVector ctrl_pts[3],
                                             int surf_id )
{
  CubitBoolean found = CUBIT_FALSE;
  BoundaryEdgeData *bed_ptr = NULL;
  for (int ii=0; ii<edgeDataList.size(); ii++)
  {
    bed_ptr = edgeDataList.get_and_step();
    if(bed_ptr->surfID == surf_id)
      found = CUBIT_TRUE;
  }
  assert( found );  // the surf_id does not match any in the list

  if (point_ptr == edgePtr->point(0))
  {
    ctrl_pts[0] = bed_ptr->bezierCtrlPts[0];
    ctrl_pts[1] = bed_ptr->bezierCtrlPts[1];
    ctrl_pts[2] = bed_ptr->bezierCtrlPts[2];
  }
  else if(point_ptr == edgePtr->point(1))
  {
    ctrl_pts[0] = bed_ptr->bezierCtrlPts[2];
    ctrl_pts[1] = bed_ptr->bezierCtrlPts[1];
    ctrl_pts[2] = bed_ptr->bezierCtrlPts[0];
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
void TDFacetBoundaryEdge::set_control_points( CubitPoint *point_ptr, 
                                              CubitVector ctrl_pts[3],
                                              int surf_id )
{
  CubitBoolean found = CUBIT_FALSE;
  BoundaryEdgeData *bed_ptr = NULL;
  for (int ii=0; ii<edgeDataList.size(); ii++)
  {
    bed_ptr = edgeDataList.get_and_step();
    if(bed_ptr->surfID == surf_id)
      found = CUBIT_TRUE;
  }
  assert( found );  // the surf_id does not match any in the list

  if (point_ptr == edgePtr->point(0))
  {
    bed_ptr->bezierCtrlPts[0] = ctrl_pts[0];
    bed_ptr->bezierCtrlPts[1] = ctrl_pts[1];
    bed_ptr->bezierCtrlPts[2] = ctrl_pts[2];
  }
  else if(point_ptr == edgePtr->point(1))
  {
    bed_ptr->bezierCtrlPts[0] = ctrl_pts[2];
    bed_ptr->bezierCtrlPts[1] = ctrl_pts[1];
    bed_ptr->bezierCtrlPts[2] = ctrl_pts[0];
  }
  else
  {
    assert(0);  // point_ptr does not match either point
  }
}

//-------------------------------------------------------------------------
// Purpose       : set the surface id of associated with one of the 
//                 facets adjacent to this edge
//
// Special Notes :
//
// Creator       : Steve Owen
//
// Creation Date : 05/01
//------------------------------------------------------------------------- 
void TDFacetBoundaryEdge::set_surf_id( CubitFacet *facet_ptr, int surf_id )
{
  int found = 0;
  int ii;
  BoundaryEdgeData *bed_ptr;
  for (ii=0; ii<edgeDataList.size() && !found; ii++)
  {
    bed_ptr = edgeDataList.get_and_step();
    if (bed_ptr->adjFacet == facet_ptr)
    {
      found = 1;
      bed_ptr->surfID = surf_id;
    }
  }
  assert(found);  // couldn't find the facet adjacent the surface
}

//-------------------------------------------------------------------------
// Purpose       : determine if this edge is internal to the surface
//                 (ie. surfID on both facets is the same)
//
// Special Notes :  Assumes surface ids have been set
//
// Creator       : Steve Owen
//
// Creation Date : 05/01
//------------------------------------------------------------------------- 
CubitBoolean TDFacetBoundaryEdge::is_internal_edge()
{
  if (edgeDataList.size() != 2)
    return CUBIT_FALSE;
  BoundaryEdgeData *bed0_ptr = edgeDataList.get_and_step();
  BoundaryEdgeData *bed1_ptr = edgeDataList.get_and_step();
  if (bed0_ptr->surfID == bed1_ptr->surfID)
  {
    return CUBIT_TRUE;
  }
  return CUBIT_FALSE;
}

//===========================================================================
//Function Name: init_control_points
//
//Member Type:  PRIVATE
//Descriptoin:  compute the control points for an edge
//===========================================================================
CubitStatus TDFacetBoundaryEdge::init_control_points( double min_dot )
{
  int ii;
  CubitStatus stat = CUBIT_SUCCESS;
  CubitVector N0, N3;
  for (ii=0; ii<edgeDataList.size(); ii++)
  {
    // get the normals specific to this facet.

    BoundaryEdgeData *bed_ptr = edgeDataList.get_and_step();
    TDFacetBoundaryPoint *td_fbp0 = 
      TDFacetBoundaryPoint::get_facet_boundary_point( edgePtr->point( 0 ) );
    TDFacetBoundaryPoint *td_fbp1 = 
      TDFacetBoundaryPoint::get_facet_boundary_point( edgePtr->point( 1 ) );
    if (!td_fbp0 || !td_fbp1)
      return CUBIT_FAILURE;
    stat = td_fbp0->get_normal( bed_ptr->adjFacet, N0 );
    if (stat != CUBIT_SUCCESS)
      return stat;
    td_fbp1->get_normal( bed_ptr->adjFacet, N3 );
    if (stat != CUBIT_SUCCESS)
      return stat;

    // determine the curve tangents

    CubitVector T0, T3;
    int tool_id = bed_ptr->adjFacet->tool_id();
    stat = compute_curve_tangent( tool_id, edgePtr, min_dot, T0, T3 );
    if (stat != CUBIT_SUCCESS)
      return stat;

    // init the control points

    CubitVector P0 = edgePtr->point(0)->coordinates();
    CubitVector P3 = edgePtr->point(1)->coordinates();
    stat = FacetEvalTool::init_edge_control_points( P0, P3, N0, N3, T0, T3,
                                                    bed_ptr->bezierCtrlPts );
    if (stat != CUBIT_SUCCESS)
      return stat;
  }
  return CUBIT_SUCCESS;
}

//===========================================================================
//Function Name: merge_control_points
//
//Member Type:  PRIVATE
//Descriptoin:  merge the control points - find average location for all
//===========================================================================
CubitStatus TDFacetBoundaryEdge::merge_control_points()
{
  int ii, jj;
  CubitVector ctrl_pts[3];
  for (ii=0; ii<edgeDataList.size(); ii++)
  {
    BoundaryEdgeData *bed_ptr = edgeDataList.get_and_step();
    for (jj=0; jj<3; jj++)
      ctrl_pts[jj] += bed_ptr->bezierCtrlPts[jj];
  }
  for (ii=0; ii<3; ii++)
    ctrl_pts[ii] = ctrl_pts[ii] / edgeDataList.size();
  edgePtr->control_points( ctrl_pts, 4 );
  return CUBIT_SUCCESS;
}

//===========================================================================
//Function Name: compute_curve_tangent
//
//Member Type:  PRIVATE
//Descriptoin:  compute the tangents to the endpoints of a feature edge.
//===========================================================================
CubitStatus TDFacetBoundaryEdge::compute_curve_tangent( 
  int tool_id,
  CubitFacetEdge *edge,
  double /*min_dot*/,
  CubitVector &T0,
  CubitVector &T3 )
{
  int mydebug = 0;
  
  CubitPoint *p0 = edge->point( 0 );
  CubitPoint *p1 = edge->point( 1 );
  CubitFacetEdge *prev_edge = NULL; //next_feature_edge( tool_id, edge, p0 );

    //If ther is a previous edge that this one should be C1 with, get
    // it, and make this one C1 with it.
  if(get_prev_edge()){
      prev_edge = get_prev_edge();
  }

  if (prev_edge == NULL)  // could be end of a hard line
  {
    T0 = p1->coordinates() - p0->coordinates();  
    T0.normalize();  
  }
  else
  {
      
    if(mydebug){
      dcolor(CUBIT_WHITE);
      dedraw(prev_edge);
      dcolor(CUBIT_MAGENTA);
      dedraw(edge);
      dview();
    }
    
    CubitPoint *p2 = prev_edge->other_point( p0 );
    T0 = (p0->coordinates() - p2->coordinates()) + 
      (p1->coordinates() - p0->coordinates());
    T0.normalize();
  }
  

  CubitFacetEdge *next_edge = NULL; //next_feature_edge( tool_id, edge, p1 );
    //If ther is a next edge that this one should be C1 with, get
    // it, and make this one C1 with it.
  if(get_next_edge()){
    next_edge = get_next_edge();
  }
  if (next_edge == NULL)  // could be end of a hard line
  {
    T3 = p1->coordinates() - p0->coordinates();
    T3.normalize();
  }
  else
  {
    if(mydebug){
      dcolor(CUBIT_YELLOW);
      dedraw(next_edge);
      dcolor(CUBIT_GREEN);
      dedraw(edge);
      dview();
    }
    CubitPoint *p2 = next_edge->other_point( p1 );
    T3 = (p2->coordinates() - p1->coordinates()) + 
      (p1->coordinates() - p0->coordinates());
    T3.normalize();
  }
  
  return CUBIT_SUCCESS;
}


//===========================================================================
//Function Name: next_feature_edge
//
//Member Type:  PRIVATE
//Descriptoin:  given a facet boundary edge and one of its nodes, find the
//              next edge on the same surface  
//===========================================================================
CubitFacetEdge *TDFacetBoundaryEdge::next_feature_edge( 
  int tool_id,
  CubitFacetEdge *this_edge,
  CubitPoint *p0 )
{
  CubitFacetEdge *next_edge = NULL;

  DLIList<CubitFacetEdge*> edge_list;
  p0->edges( edge_list );
  int ii;

  CubitFacetEdge *edge_ptr = NULL;
  TDFacetBoundaryEdge *td_fbe;
  for (ii=0; ii<edge_list.size() && next_edge == NULL; ii++)
  {
    edge_ptr = edge_list.get_and_step();
    if (edge_ptr != this_edge)
    {
      td_fbe = TDFacetBoundaryEdge::get_facet_boundary_edge( edge_ptr );
      if (td_fbe != NULL && td_fbe->is_at_surf( tool_id ))
      {
        next_edge = edge_ptr;
      }
    }
  }

  return next_edge;
}

//===========================================================================
//Function Name: is_at_surf
//
//Member Type:  PRIVATE
//Descriptoin:  return whether the given edge has one of its adjacent 
//              faces on the given surface
//===========================================================================
CubitBoolean TDFacetBoundaryEdge::is_at_surf( int surf_id )
{
  CubitBoolean found = CUBIT_FALSE;
  BoundaryEdgeData *bed_ptr;
  for (int ii=0; ii<edgeDataList.size() && !found; ii++)
  {
    bed_ptr = edgeDataList.get_and_step();
    if(bed_ptr->surfID == surf_id ||
       bed_ptr->adjFacet->tool_id() == surf_id)
      found = CUBIT_TRUE;
  }
  return found;
}
