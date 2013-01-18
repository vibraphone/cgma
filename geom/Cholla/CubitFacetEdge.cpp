//======================================================================= 
// 
//  File: CubitFacetEdge 
//  Description: used for edges of CubitFacets.  These are optional 
//               and used only if specific information must be stored at  
//               the edge common to two facets. 
//  Notes: See notes in CubitFacetEdge.hpp 
//  Owner: sjowen 
// 
//======================================================================= 
 
#include "CubitFacetEdge.hpp" 
#include "CubitFacet.hpp" 
#include "CubitPoint.hpp" 
#include "CubitVector.hpp" 
#include "GeometryDefines.h" 
#include "FacetEvalTool.hpp" 
#include "GfxDebug.hpp"
#include "IntersectionTool.hpp"
#include "CubitMessage.hpp"

 
//====================================================================== 
// Function: CubitFacetEdge (PUBLIC) 
// Description: constructor 
// Author: sjowen 
// Date: 8/00 
//====================================================================== 
CubitFacetEdge::CubitFacetEdge( )  
    : bezierOrder(0), markedFlag(0), isFeature(0), isFlipped(CUBIT_FALSE) 
{  
} 
 
 
//====================================================================== 
// Function: ~CubitFacetEdge (PUBLIC) 
// Description: destructor 
// Author: sjowen 
// Date: 8/00 
//====================================================================== 
CubitFacetEdge::~CubitFacetEdge() 
{ 
} 
 
//====================================================================== 
// Function: set_control_points (PUBLIC) 
// Description: set the control points from an array of doubles 
// Author: sjowen 
// Date: 8/00 
//====================================================================== 
void CubitFacetEdge::set_control_points( const double *ctrl_pt_array )  
{ 
  int ii;
  for (ii=0; ii<3; ii++)
  {
    controlPoints[ii].x( ctrl_pt_array[ii*3] );
    controlPoints[ii].y( ctrl_pt_array[ii*3+1] );
    controlPoints[ii].z( ctrl_pt_array[ii*3+2] );
  }
}

//====================================================================== 
// Function: control_points (PUBLIC) 
// Description: return the Bezier control points (including the end points) 
//        The order of the bezier is returned. 
// Author: sjowen 
// Date: 8/00 
//====================================================================== 
int CubitFacetEdge::control_points( CubitVector *ctrl_pts )  
{ 
  DLIList<CubitPoint*> my_points; 
  points(my_points); 
   
  ctrl_pts[0] = my_points.get()->coordinates(); 
  for(int i=0; i<bezierOrder-1; i++) { 
    ctrl_pts[i+1] = controlPoints[i]; 
  } 
  ctrl_pts[bezierOrder] = my_points.next()->coordinates(); 
  return bezierOrder; 
}
 
//====================================================================== 
// Function: control_points (PUBLIC) 
// Description: return the control points on the edge of a facet based 
//        on its edge use direction 
// Author: sjowen 
// Date: 8/00 
//====================================================================== 
CubitStatus CubitFacetEdge::control_points(  
  CubitFacet *facet, CubitVector *ctrl_pts ) 
{ 
  int index = -1; 
  CubitBoolean found = CUBIT_FALSE; 
  for (int i=0; i<3 && !found; i++) { 
    if (this == facet->edge(i)) { 
      index = i; 
      found = CUBIT_TRUE; 
    } 
  } 
  if (!found) { 
    return CUBIT_FAILURE; 
  } 
 
  DLIList<CubitPoint*> my_points; 
  points(my_points); 
   
  switch (facet->edge_use(index)) { 
  case 1: 
    ctrl_pts[0] = my_points.get()->coordinates(); 
    ctrl_pts[1] = controlPoints[0]; 
    ctrl_pts[2] = controlPoints[1]; 
    ctrl_pts[3] = controlPoints[2]; 
    ctrl_pts[4] = my_points.next()->coordinates(); 
    break; 
  case -1: 
    ctrl_pts[0] = my_points.next()->coordinates(); 
    ctrl_pts[1] = controlPoints[2]; 
    ctrl_pts[2] = controlPoints[1]; 
    ctrl_pts[3] = controlPoints[0]; 
    ctrl_pts[4] = my_points.get()->coordinates(); 
    break; 
  default: 
    return CUBIT_FAILURE; 
  } 
  return CUBIT_SUCCESS; 
} 
 
//=========================================================================== 
//Function Name: evaluate_position 
// 
//Member Type:  PUBLIC 
//Description:  evaluate the facet edge at a position 
//              eval_tangent is NULL if tangent not needed 
//=========================================================================== 
CubitStatus CubitFacetEdge::evaluate_position( const CubitVector &start_position, 
                                               CubitVector *eval_point, 
                                               CubitVector *eval_tangent)  
{ 
  CubitStatus stat = CUBIT_SUCCESS; 
 
  // find the adjacent facet 
 
  CubitFacet *facet_ptr = this->adj_facet( 0 ); 
 
  // If there is none or this is a linear representation -  
  // then project to the linear edge  
 
  if (!facet_ptr || facet_ptr->eval_order() == 0 || facet_ptr->is_flat()) 
  { 
    if (eval_point) 
    { 
      closest_point(start_position, *eval_point); 
    } 
    if (eval_tangent) 
    { 
      *eval_tangent = point(1)->coordinates() - point(0)->coordinates(); 
      (*eval_tangent).normalize(); 
    } 
  } 
  else 
  { 
    int vert0 = facet_ptr->point_index( point(0) ); 
    int vert1 = facet_ptr->point_index( point(1) ); 
    CubitVector pt_on_plane, close_point; 
    CubitVector start = start_position; 
    double dist_to_plane; 
    CubitBoolean outside_facet; 
    FacetEvalTool::project_to_facet_plane( facet_ptr, start,  
                                           pt_on_plane, dist_to_plane ); 
    stat = FacetEvalTool::project_to_facetedge( facet_ptr,  
                                                vert0, vert1, 
                                                start, 
                                                pt_on_plane,  
                                                close_point, 
                                                outside_facet ); 
    if (eval_point) 
    { 
      *eval_point = close_point; 
    } 
    if (eval_tangent) 
    { 
      CubitVector edvec = point(1)->coordinates() - point(0)->coordinates(); 
      edvec.normalize(); 
      CubitVector areacoord; 
      FacetEvalTool::facet_area_coordinate( facet_ptr, close_point, areacoord );  
      FacetEvalTool::eval_facet_normal(facet_ptr, areacoord, *eval_tangent); 
      CubitVector cross = edvec * *eval_tangent; 
      *eval_tangent = *eval_tangent * cross; 
      (*eval_tangent).normalize(); 
    } 
  } 
  return CUBIT_SUCCESS; 
} 
 
//=========================================================================== 
//Function Name: evaluate 
// 
//Member Type:  PUBLIC 
//Description:  evaluate the facet at area coordinates 
//              eval_normal is NULL if normal not needed 
//Note:         t is a value from -1 to 1.  t=0 is the midpoint 
//=========================================================================== 
CubitStatus CubitFacetEdge::evaluate( double &t,  
                                      CubitVector *eval_point, 
                                      CubitVector *eval_tangent )      
{ 
  CubitStatus stat = CUBIT_SUCCESS; 
 
  // project the position to the linear edge 
   
  double tt = (t + 1) * 0.5; 
  if (tt <= 0.0) tt = 0.0; 
  if (tt >= 1.0) tt = 1.0; 
  *eval_point = point(0)->coordinates() +  
    tt * (point(1)->coordinates() - point(0)->coordinates()); 
 
  // evaluate the point on the facet (if the order is higher than 0) 
 
  CubitFacet *facet_ptr = this->adj_facet( 0 ); 
  if (!facet_ptr || facet_ptr->is_flat()) 
  { 
    if (eval_tangent) 
    { 
      *eval_tangent = point(1)->coordinates() - point(0)->coordinates(); 
      (*eval_tangent).normalize(); 
    } 
  } 
  else 
  { 
    CubitVector areacoord; 
    FacetEvalTool::facet_area_coordinate( facet_ptr, *eval_point, areacoord ); 
    stat = facet_ptr->evaluate( areacoord, eval_point, eval_tangent ); 
    if (stat != CUBIT_SUCCESS) 
      return stat; 
    if (eval_tangent) 
    { 
      CubitVector edvec = point(1)->coordinates() - point(0)->coordinates(); 
      edvec.normalize(); 
      CubitVector cross = edvec * *eval_tangent; 
      *eval_tangent = *eval_tangent * cross; 
      (*eval_tangent).normalize(); 
    } 
  } 
  return stat; 
} 
 
//=========================================================================== 
//Function Name: evaluate_single 
// 
//Member Type:  PUBLIC 
//Description:  evaluate edge not associated with a facet. 
//Note:         t is a value from -1 to 1.
//=========================================================================== 
CubitStatus CubitFacetEdge::evaluate_single(double &t,
                            CubitVector *outv)
{
CubitVector P0, P1;
double t4, t3, t2, one_minus_t, one_minus_t2, one_minus_t3, one_minus_t4;
DLIList<CubitPoint*> my_points; 
  // project the position to the linear edge 
   
  double tt = (t + 1) * 0.5; 
  if (tt <= 0.0) tt = 0.0; 
  if (tt >= 1.0) tt = 1.0; 

  points(my_points);
 
  P0 = my_points.get()->coordinates();
  P1 = my_points.next()->coordinates();

  t2 = tt*tt;
  t3 = t2*tt;
  t4 = t3*tt;  
  one_minus_t = 1.-tt;
  one_minus_t2 = one_minus_t*one_minus_t;
  one_minus_t3 = one_minus_t2*one_minus_t;
  one_minus_t4 = one_minus_t3*one_minus_t;
  
  *outv = one_minus_t4*P0 + 
       4.*one_minus_t3*tt* controlPoints[0] + 
       6.*one_minus_t2*t2*controlPoints[1] + 
       4.*one_minus_t* t3*controlPoints[2] + 
                    t4*P1;

  return CUBIT_SUCCESS;
}
 
//=========================================================================== 
//Function Name: evaluate_single_tangent 
// 
//Member Type:  PUBLIC 
//Description:  evaluate tangent to edge not associated with a facet. 
//Note:         t is a value from -1 to 1.
//=========================================================================== 
CubitStatus CubitFacetEdge::evaluate_single_tangent(double &t,
                                                    CubitVector *outv)
{
CubitVector P0, P1;
double t4, t3, t2, one_minus_t, one_minus_t2, one_minus_t3, one_minus_t4;
DLIList<CubitPoint*> my_points; 
  // project the position to the linear edge 
   
  double tt = (t + 1) * 0.5; 
  if (tt <= 0.0) tt = 0.0; 
  if (tt >= 1.0) tt = 1.0; 

  points(my_points);
 
  P0 = my_points.get()->coordinates();
  P1 = my_points.next()->coordinates();

  t2 = tt*tt;
  t3 = t2*tt;
  t4 = t3*tt;  
  one_minus_t = 1.-tt;
  one_minus_t2 = one_minus_t*one_minus_t;
  one_minus_t3 = one_minus_t2*one_minus_t;
  one_minus_t4 = one_minus_t3*one_minus_t;
  
  *outv = -4.*one_minus_t3*P0 +
          4.*(one_minus_t3 -3.*tt*one_minus_t2)*controlPoints[0] +
          12.*(tt*one_minus_t2 - t2*one_minus_t)*controlPoints[1] +
          4.*(3.*t2*one_minus_t - t3)*controlPoints[2] +
          4.*t3*P1;

  return CUBIT_SUCCESS;
}
  
//=========================================================================== 
//Function Name: evaluate_2nd_derivative 
// 
//Member Type:  PUBLIC 
//Description:  evaluate tangent to edge not associated with a facet. 
//Note:         t is a value from -1 to 1.
//=========================================================================== 
CubitStatus CubitFacetEdge::evaluate_2nd_derivative(double &t,
                                                    CubitVector *outv)
{
CubitVector P0, P1, second_d;
DLIList<CubitPoint*> my_points; 
double val;
  // project the position to the linear edge 
   
  double tt = (t + 1) * 0.5; 
  if (tt <= 0.0) tt = 0.0; 
  if (tt >= 1.0) tt = 1.0; 

  points(my_points);
 
  P0 = my_points.get()->coordinates();
  P1 = my_points.next()->coordinates();

  val = 12.*(1.-tt)*(1.-tt)*P0.x() -
        24.*(2.*tt*tt - 3.*tt + 1.)*controlPoints[0].x() +
        12.*(6.*tt*tt - 6.*tt + 1.)*controlPoints[1].x() +
        24.*(tt - 2.*tt*tt)*controlPoints[2].x() +
        12.*tt*tt*P1.x();
  second_d.x(val);
  val = 12.*(1.-tt)*(1.-tt)*P0.y() -
        24.*(2.*tt*tt - 3.*tt + 1.)*controlPoints[0].y() +
        12.*(6.*tt*tt - 6.*tt + 1.)*controlPoints[1].y() +
        24.*(tt - 2.*tt*tt)*controlPoints[2].y() +
        12.*tt*tt*P1.y();
  second_d.y(val);
  val = 12.*(1.-tt)*(1.-tt)*P0.z() -
        24.*(2.*tt*tt - 3.*tt + 1.)*controlPoints[0].z() +
        12.*(6.*tt*tt - 6.*tt + 1.)*controlPoints[1].z() +
        24.*(tt - 2.*tt*tt)*controlPoints[2].z() +
        12.*tt*tt*P1.z();
  second_d.z(val);
  *outv = second_d;
  
  return CUBIT_SUCCESS;
}

//=========================================================================== 
//Function Name: closest_point 
// 
//Member Type:  PUBLIC 
//Description:  return the closest point on segment defined by the edge 
//=========================================================================== 
CubitStatus CubitFacetEdge::closest_point(const CubitVector &point,  
                                          CubitVector &closest_point ) 
{  
  //CubitStatus rv = CUBIT_SUCCESS; 
  CubitPoint *pt0 = this->point(0); 
  CubitPoint *pt1 = this->point(1); 
 
  // the edge vector 
 
  CubitVector e0 ( pt1->x() - pt0->x(), 
                   pt1->y() - pt0->y(), 
                   pt1->z() - pt0->z() ); 
  double elen = e0.normalize(); 
   
  // vector from vert0 to point 
 
  CubitVector v0 ( point.x() - pt0->x(), 
                   point.y() - pt0->y(), 
                   point.z() - pt0->z() ); 
   
  // project to edge 
 
  double len = v0%e0; 
  if (len <= 0.0) 
  { 
    closest_point = pt0->coordinates(); 
  } 
  else if( len >= elen ) 
  { 
    closest_point = pt1->coordinates(); 
  } 
  else 
  { 
    closest_point.x ( pt0->x() + e0.x() * len ); 
    closest_point.y ( pt0->y() + e0.y() * len ); 
    closest_point.z ( pt0->z() + e0.z() * len ); 
  } 
 
  return CUBIT_SUCCESS; 
}

 //=========================================================================== 
//Function Name: intersect 
// 
//Member Type:  PUBLIC 
//Description:  intersect the edge with a segment.  Assumes segment and edge
//              are on the same plane (project to facet plane first) 
//=========================================================================== 
CubitStatus CubitFacetEdge::intersect(
  CubitVector &aa, CubitVector &bb, // end point of segment
  CubitVector &norm,  // normal of the common plane
  CubitVector &qq,  // return the intersection point 
  CubitBoolean &does_intersect ) // return status of intersection
{  
 
  CubitPoint *pt0 = this->point(0); 
  CubitPoint *pt1 = this->point(1); 

  double P0[2], P1[2], AA[2], BB[2];
  CubitVector absnorm(fabs(norm.x()), fabs(norm.y()), fabs(norm.z()));
  if (absnorm.x() >= absnorm.y() && absnorm.x() >= absnorm.z())
  {
    P0[0] = pt0->coordinates().y();  P0[1] = pt0->coordinates().z();
    P1[0] = pt1->coordinates().y();  P1[1] = pt1->coordinates().z();
    AA[0] = aa.y();                  AA[1] = aa.z();
    BB[0] = bb.y();                  BB[1] = bb.z();
  }
  else if (absnorm.y() >= absnorm.x() && absnorm.y() >= absnorm.z())
  {
    P0[0] = pt0->coordinates().z();  P0[1] = pt0->coordinates().x();
    P1[0] = pt1->coordinates().z();  P1[1] = pt1->coordinates().x();
    AA[0] = aa.z();                  AA[1] = aa.x();
    BB[0] = bb.z();                  BB[1] = bb.x();
  }
  else
  {
    P0[0] = pt0->coordinates().x();  P0[1] = pt0->coordinates().y();
    P1[0] = pt1->coordinates().x();  P1[1] = pt1->coordinates().y();
    AA[0] = aa.x();                  AA[1] = aa.y();
    BB[0] = bb.x();                  BB[1] = bb.y();
  }

  double QQ[4], s;
  int ninter = intersect_2D_segments(P0, P1, AA, BB, QQ);

  if (ninter != 1)
  {
    does_intersect = CUBIT_FALSE;
    return CUBIT_SUCCESS;
  }
  does_intersect = CUBIT_TRUE;
 
  double dx = P1[0] - P0[0];
  double dy = P1[1] - P0[1];
  if (fabs(dx) > fabs(dy))
    s = (QQ[0] - P0[0]) / dx;
  else
    s = (QQ[1] - P0[1]) / dy;
  
  qq = pt0->coordinates() + s * (pt1->coordinates() - pt0->coordinates());
 
  return CUBIT_SUCCESS; 
}  
 
//====================================================================== 
// Function: boundary_edge_points (PUBLIC) 
// Description: return the oriented endpoints of a facet edge assuming 
//        the edge is at the frontier of the facets (only one adjacency) 
// Author: sjowen 
// Date: 8/00 
//====================================================================== 
void CubitFacetEdge::boundary_edge_points( CubitPoint * &pt0,  
                                           CubitPoint * &pt1, 
                                           int tool_id ) 
{ 
  if(num_adj_facets() == 0) 
  { 
    pt0 = point(0); 
    pt1 = point(1); 
    return; 
  } 
  CubitFacet *facet = NULL; 
  if (tool_id != 0) 
  { 
    int ii; 
    int found = 0; 
    DLIList <CubitFacet *> adj_facets; 
    facets(adj_facets); 
    for (ii=0; ii<adj_facets.size() && !found; ii++) 
    { 
      facet = adj_facets.get_and_step(); 
      if (facet->tool_id() == tool_id) 
        found = 1; 
    } 
    assert(found); 
  } 
  else 
  { 
    facet = adj_facet(0); 
    if (!facet) facet = adj_facet(1); 
    assert(facet != 0); 
  } 
   
  int index = facet->edge_index( this ); 
  int use = facet->edge_use( index ); 
  if (use == 1) { 
    pt0 = point(0); 
    pt1 = point(1); 
  } 
  else { 
    pt0 = point(1); 
    pt1 = point(0); 
  } 
} 
 
//====================================================================== 
// Function: dist_to_edge 
// Description: return the distance from the point to this edge 
// Author: sjowen 
// Date: 2/01 
// Corrected by JFowler 5/03
//====================================================================== 
double CubitFacetEdge::dist_to_edge( 
  const CubitVector &this_point,  
  CubitVector &close_point,  
  CubitBoolean &outside_edge ) 
{ 
  double dist = 0.0; 
  CubitVector p0 = point(0)->coordinates(); 
  CubitVector p1 = point(1)->coordinates(); 
  CubitVector edge_vec( p1, p0 ); 
  CubitVector point_vec( this_point, p0 ); 
  double edge_length;  
  edge_length = edge_vec.normalize(); 
  double dist_on_edge = edge_vec % point_vec; 
  if (dist_on_edge < 0.0e0) 
  { 
    close_point = p0; 
    outside_edge = CUBIT_TRUE; 
  } 
  else if (dist_on_edge > edge_length) 
  { 
    close_point = p1; 
    outside_edge = CUBIT_TRUE; 
  } 
  else 
  { 
    close_point = p0 - edge_vec * dist_on_edge; 
    outside_edge = CUBIT_FALSE; 
  } 
  dist = close_point.distance_between( this_point );  
  return dist; 
}
 
/*//====================================================================== 
// Function: dist_to_edge 
// Description: return the distance from the point to this edge 
// Author: sjowen 
// Date: 2/01 
//====================================================================== 
double CubitFacetEdge::dist_to_edge( 
  CubitVector &this_point,  
  CubitVector &close_point,  
  CubitBoolean &outside_edge ) 
{ 
  double dist = 0.0; 
  CubitVector p0 = point(0)->coordinates(); 
  CubitVector p1 = point(1)->coordinates(); 
  CubitVector edge_vec( p1, p0 ); 
  CubitVector point_vec( this_point, p0 ); 
  point_vec.normalize(); 
  double dist_on_edge = edge_vec % point_vec; 
  double edge_length; 
  edge_length = edge_vec.normalize(); 
  if (dist_on_edge < 0.0e0) 
  { 
    close_point = p0; 
    outside_edge = CUBIT_TRUE; 
  } 
  else if (dist_on_edge > edge_length) 
  { 
    close_point = p1; 
    outside_edge = CUBIT_TRUE; 
  } 
  else 
  { 
    close_point = p0 + edge_vec * dist_on_edge; 
    outside_edge = CUBIT_FALSE; 
  } 
  dist = close_point.distance_between( this_point );  
  return dist; 
} 
*/ 
//====================================================================== 
// Function: proj_to_line 
// Description: project the point to a line defined by the edge 
// Author: sjowen 
// Date: 2/01 
//====================================================================== 
CubitStatus  CubitFacetEdge::proj_to_line(  
  const CubitVector &this_point,  
  CubitVector &proj_point ) 
{ 
  CubitStatus stat = CUBIT_SUCCESS; 
  CubitVector p0 = point(0)->coordinates(); 
  CubitVector p1 = point(1)->coordinates(); 
  CubitVector edge_vec( p0,p1 ); 
  CubitVector point_vec( p0, this_point ); 
  edge_vec.normalize(); 
  double dist_on_edge = edge_vec % point_vec;   
  proj_point = p0 + (edge_vec * dist_on_edge); 
 
  return stat; 
} 
 
//====================================================================== 
// Function: edge_tangent 
// Description: return tangent at a point on the edge 
// Author: sjowen 
// Date: 2/01 
//====================================================================== 
CubitStatus  CubitFacetEdge::edge_tangent(  
  const CubitVector &/*point_on_edge*/,  
  CubitVector &tangent ) 
{ 
  CubitStatus stat = CUBIT_SUCCESS; 
  tangent = point(1)->coordinates() - 
            point(0)->coordinates(); 
  tangent.normalize(); 
  return stat; 
}

//====================================================================== 
// Function: edge_tangent 
// Description: return curvature at a point on the edge 
// Author: sjowen 
// Date: 2/01 
//====================================================================== 
CubitStatus  CubitFacetEdge::edge_curvature(  
  const CubitVector &/*point_on_edge*/,  
  CubitVector &curvature,
  CubitFacetEdge *closest_edge ) 
{ 
  CubitVector vec_ba, vec_ca, center_point;

  //if point(0) is middle point
  if( closest_edge->other_point( point(0) ) )  
  {
    center_point = point(0)->coordinates();
    vec_ba = closest_edge->point(0)->coordinates() - center_point; 
    vec_ca = point(1)->coordinates() - center_point; 
  }
  //if point(1) is middle point
  else if( closest_edge->other_point( point(1) ) )  
  {
    center_point = point(1)->coordinates();
    vec_ba = point(0)->coordinates() - center_point; 
    vec_ca = closest_edge->point(1)->coordinates() - center_point; 
  }
  else
    assert(0);
  
  // Squares of lengths of the edges incident to `a'.
  double ba_length = vec_ba.length_squared();
  double ca_length = vec_ca.length_squared();
  
  // Cross product of these edges.
  // (Take your chances with floating-point roundoff.)
  CubitVector cross_bc = vec_ba * vec_ca;
  
  // Calculate the denominator of the formulae.
  double temp_dbl = cross_bc % cross_bc;
  CubitVector circle_center(0.0,0.0,0.0);
  if(fabs(temp_dbl) > CUBIT_DBL_MIN){
    double denominator = 0.5 / (temp_dbl);
    assert(denominator != 0.0);
  
    // Calculate offset (from `a') of circumcenter.
    circle_center  = (ba_length * vec_ca - ca_length * vec_ba) * cross_bc;
    circle_center *= denominator;

    //store radius
    double radius = circle_center.length();
    circle_center.normalize();
    circle_center /= radius;
  } 
  curvature = circle_center; 

  return CUBIT_SUCCESS; 
} 
 
//====================================================================== 
// Function: length 
// Description: return length of an edge 
// Author: sjowen 
// Date: 2/01 
//====================================================================== 
double CubitFacetEdge::length() 
{ 
  CubitVector start = point(0)->coordinates(); 
  CubitVector end = point(1)->coordinates(); 
  return start.distance_between( end ); 
} 

//====================================================================== 
// Function: length 
// Description: return length of an edge 
// Author: jakraft 
// Date: 3/04
//====================================================================== 
CubitVector CubitFacetEdge::position_from_fraction( double f )
{
  return (1.0 - f) * point(0)->coordinates() +
                f  * point(1)->coordinates();
}

//====================================================================== 
// Function: other_point 
// Description: return the other point on the edge 
// Author: sjowen 
// Date: 4/01 
//====================================================================== 
CubitPoint* CubitFacetEdge::other_point( CubitPoint *point_ptr ) 
{ 
  if (point(0) == point_ptr) 
    return point(1); 
  if(point(1) == point_ptr) 
    return point(0); 
  return NULL; 
} 
 
//====================================================================== 
// Function: get_parents 
// Description: return entities attached to this edge 
// Author: sjowen 
// Date: 4/01 
//====================================================================== 
void CubitFacetEdge::get_parents(DLIList<FacetEntity *> &facet_list) 
{ 
  DLIList<CubitFacet *> cf_list; 
  facets( cf_list ); 
  for (int ii=0; ii<cf_list.size(); ii++) 
    facet_list.append(cf_list.get_and_step()); 
} 
 
//====================================================================== 
// Function: other_facet 
// Description: return the other facet on the edge (assumes 
//              max two facets per edge) 
// Author: sjowen 
// Date: 5/01 
//====================================================================== 
CubitFacet *CubitFacetEdge::other_facet( CubitFacet *facet_ptr ) 
{ 
  DLIList<CubitFacet *> cf_list; 
  facets( cf_list ); 
  assert(cf_list.size() < 3); 
  CubitFacet *adj_facet = NULL; 
  if (cf_list.size() > 0) 
  { 
    adj_facet = cf_list.get_and_step(); 
    if (adj_facet == facet_ptr) 
    { 
      if (cf_list.size() == 2) 
      { 
        adj_facet = cf_list.get(); 
      } 
    } 
  } 
  return adj_facet; 
} 
 
//====================================================================== 
// Function: other_facet_on_surf 
// Description: return the other facet on the edge that has the 
//              same tool id. (finds the first occurrence 
//              of the same tool id of the adjacent facets 
//              to this edge that is not facet_ptr) 
// Author: sjowen 
// Date: 5/01 
//====================================================================== 
CubitFacet *CubitFacetEdge::other_facet_on_surf( CubitFacet *facet_ptr ) 
{ 
  assert(facet_ptr != 0); 
  int tool_id = facet_ptr->tool_id(); 
  DLIList<CubitFacet *> cf_list; 
  facets( cf_list ); 
  CubitFacet *adj_facet = NULL; 
  int found = 0; 
  for (int ii=0; ii<cf_list.size() && !found; ii++) 
  { 
    adj_facet = cf_list.get_and_step(); 
    if (adj_facet != facet_ptr) 
    { 
      if (adj_facet->tool_id() == tool_id) 
        found = 1; 
    } 
  } 
  if (!found) 
    adj_facet = NULL; 
  return adj_facet; 
} 
 
//====================================================================== 
// Function: num_adj_facets_on_surf 
// Description: count the number of adjacent facets to this edge that 
//              have the specified tool id 
// Author: sjowen 
// Date: 5/01 
//====================================================================== 
int CubitFacetEdge::num_adj_facets_on_surf( int tool_id ) 
{ 
 
  DLIList<CubitFacet *> cf_list; 
  facets( cf_list ); 
  CubitFacet *adj_facet = NULL; 
  int nfacets = 0; 
  for (int ii=0; ii<cf_list.size(); ii++) 
  { 
    adj_facet = cf_list.get_and_step(); 
    if (adj_facet->tool_id() == tool_id) 
      nfacets++; 
  } 
  return nfacets; 
} 
 
//====================================================================== 
// Function: adj_facet_on_surf 
// Description: return the first facet on the adjacent facet list with  
//              the indicated tool id 
// Author: sjowen 
// Date: 5/01 
//====================================================================== 
CubitFacet *CubitFacetEdge::adj_facet_on_surf( int tool_id ) 
{ 
  DLIList<CubitFacet *> cf_list; 
  facets( cf_list ); 
  CubitFacet *adj_facet = NULL; 
  int found = 0; 
  for (int ii=0; ii<cf_list.size() && !found; ii++) 
  { 
    adj_facet = cf_list.get_and_step(); 
    if (adj_facet->tool_id() == tool_id) 
    { 
      found = 1; 
    } 
  } 
  if (!found) 
    adj_facet = NULL; 
  return adj_facet; 
} 
 
//====================================================================== 
// Function: contains 
// Description: determines if point is contained in edge 
// Author: sjowen 
// Date: 5/01 
//====================================================================== 
CubitBoolean CubitFacetEdge::contains( CubitPoint *point_ptr ) 
{ 
  if (point(0) == point_ptr || point(1) == point_ptr) 
    return CUBIT_TRUE; 
  return CUBIT_FALSE; 
} 
 
//====================================================================== 
// Function: draw 
// Description: draw the edge 
// Author: sjowen 
// Date: 5/01 
//====================================================================== 
void CubitFacetEdge::debug_draw(int color, int flush, int /*draw_uv*/) 
{ 
  if ( color == -1 ) 
    color = CUBIT_RED; 
  GfxDebug::draw_facet_edge(this, color);
  GfxDebug::draw_point( point(0)->coordinates(), color );
  GfxDebug::draw_point( point(1)->coordinates(), color );
  if ( flush ) 
    GfxDebug::flush(); 
} 

//======================================================================
// Function: shared_point (PUBLIC)
// Description: find the common point
// Author: sjowen
// Date: 10/02
//======================================================================
CubitPoint *CubitFacetEdge::shared_point( CubitFacetEdge *edge_ptr )
{
  CubitPoint *pA = this->point(0);
  CubitPoint *pB = this->point(1);
  CubitPoint *pC = edge_ptr->point(0);
  CubitPoint *pD = edge_ptr->point(1);
  CubitPoint *pShared = NULL;
  if (pA == pC || pA == pD)
  {
    pShared = pA;
  }
  else if (pB == pC || pB == pD)
  {
    pShared = pB;
  }
  return pShared;
}


//======================================================================
// Function: bounding_box (PUBLIC)
// Description: return the bounding box of the edge
// Author: sjowen
// Date: 10/02
//======================================================================
CubitBox CubitFacetEdge::bounding_box( )
{
  CubitPoint *p1 = point (0);
  CubitPoint *p2 = point (1);

  CubitVector bbox_min, bbox_max;
  bbox_min.x(CUBIT_MIN(p1->x(),p2->x()));
  bbox_min.y(CUBIT_MIN(p1->y(),p2->y()));
  bbox_min.z(CUBIT_MIN(p1->z(),p2->z()));
  bbox_max.x(CUBIT_MAX(p1->x(),p2->x()));
  bbox_max.y(CUBIT_MAX(p1->y(),p2->y()));
  bbox_max.z(CUBIT_MAX(p1->z(),p2->z()));
  CubitBox edge_box(bbox_min,bbox_max);
  return edge_box;
}

//======================================================================
//
// The following 2 function really belong in IntersectionTool in the util
// directory.  Because of the current distribution policy it was not easy
// to get it into there without breaking the build.  These should be
// removed at a later time and integrated with the IntersectionTool
//
//=======================================================================
// Intersection of two line segments in 2D. Intersect [P0,P1] with [P2,P3]  
// Return value is zero if there is no intersection, 1 if there is a unique 
// intersection, and 2 if the 2 segments overlap and the intersection set is 
// a segment itself.  The return value is the number of valid entries*2 in 
// the array qq. 
int CubitFacetEdge::intersect_2D_segments( double P0[2], double P1[2],
                                             double P2[2], double P3[2],
                                             double qq[4] )
{

  double D0[2], D1[2];
  D0[0] = P1[0] - P0[0];  D0[1] = P1[1] - P0[1];
  D1[0] = P3[0] - P2[0];  D1[1] = P3[1] - P2[1];

  // segments P0 + s * D0 for s in [0,1],
  //          P2 + t * D1 for t in [0,1]

  double sqr_epsilon = DBL_EPSILON * DBL_EPSILON;
  double E[2];
  E[0] = P2[0] - P0[0];
  E[1] = P2[1] - P0[1];

  double kross = D0[0] * D1[1] - D0[1] * D1[0];
  double sqr_kross = kross * kross;
  double sqr_len0 = D0[0] * D0[0] + D0[1] * D0[1];
  double sqr_len1 = D1[0] * D1[0] + D1[1] * D1[1];
  if (sqr_kross > sqr_epsilon  * sqr_len0 * sqr_len1)
  {
    // lines of the segment are not parallel

    double s = (E[0] * D1[1] - E[1] * D1[0]) / kross;
    if (s < 0.0 || s > 1.0)
    {
       // intersection of lines is not a point on segment P0 + s * D0
      return 0;
    }

    double t = (E[0] * D0[1] - E[1] * D0[0]) / kross;
    if (t < 0.0 || t > 1.0)
    {
      // intersection of lines is not a point on segment P1 + t * D1
      return 0;
    }

    // intersection of lines is a point on each segment

    qq[0] = P0[0] + s * D0[0];
    qq[1] = P0[1] + s * D0[1];
    return 1;
  }

  // lines of the segments are parallel

  double sqr_lenE = E[0] * E[0] + E[1] * E[1];
  kross = E[0] * D0[1] - E[1] * D0[1];
  sqr_kross = kross * kross;
  if (sqr_kross > sqr_epsilon  * sqr_len0 * sqr_lenE)
  {
    // lines of the segments are different
    return 0;
  }

  // lines of the segment are the same.  Need to test for overlap of segments

  double s0 = (D0[0] * E[0] + D0[1] * E[1]) / sqr_len0;
  double s1 = s0 + (D0[0] * D1[0] + D0[1] * D1[1]) / sqr_len0;
  double smin = CUBIT_MIN(s0, s1);
  double smax = CUBIT_MAX(s0, s1);

  double w[2];
  int imax = intersect_intervals(0.0, 1.0, smin, smax, w);
  for (int i=0; i<imax; i++)
  {
    qq[i*2]   = P0[0] + w[i] * D0[0];
    qq[i*2+1] = P0[1] + w[i] * D0[1];
  }
  return imax;
}

// The intersection of two intervals [u0,u1] and [v0,v1], where u0<u1 and v0<v1.
// Return value is zero if intervals do not intersect; 1 if they intersect at
// a single point, in which case w[0] contains that point; or 2 if they intersect 
// in an interval whose endpoints are stored in w[0] and w[1]
int CubitFacetEdge::intersect_intervals( double u0, double u1,
                                           double v0, double v1,
                                           double w[2] )
{
  if (u1 < v0 || u0 > v1)
    return 0;

  if (u1 > v0)
  {
    if (u0 < v1)
    {
      if (u0 < v0) 
        w[0] = v0;
      else 
        w[0] = u0;
      if (u1 > v1) 
        w[1] = v1;
      else 
        w[1] = u1;
      return 2;
    }
    else
    {
      w[0] = u0;
      return 1;
    }
  }
  else
  {
    w[0] = u1;
    return 1;
  }
  return 0;
}

//======================================================================
// Function: add_facets (PUBLIC)
// Description: add this edge to its adjacent facets
// Author: sjowen
// Date: 04/04
//======================================================================
void CubitFacetEdge::add_facets( )
{
  CubitPoint *p0 = point (0);
  CubitPoint *p1 = point (1);

  DLIList<CubitFacet *> adj_facets;
  p0->shared_facets(p1, adj_facets );

  CubitFacet *facet;
  int ii;
  for(ii=0; ii<adj_facets.size(); ii++)
  {
    facet = adj_facets.get_and_step();
    facet->add_edge( this );
  }
}

double CubitFacetEdge::angle_between_facets()
{
  CubitFacet *facet0 = adj_facet(0);
  CubitFacet *facet1 = adj_facet(1);

  // assumes facets are always oriented with outwards pointing normal

  CubitVector n0 = facet0->normal();
  CubitVector n1 = facet1->normal();

  // get orientation of edge with respect to facet0

  CubitPoint *p0 = point(0);
  CubitPoint *p1 = point(1);
  CubitPoint *pnext = facet0->next_node(p0);
  if (pnext != p1)
  {
    pnext = p0;
    p0 = p1;
    p1 = pnext;
  }
  CubitVector evec = p1->coordinates() - p0->coordinates();
  evec.normalize();
  CubitVector cross = n0 * n1;  cross.normalize();
  double angle;
  double edot = evec % cross;
  double ndot = n0 % n1;
  if (ndot >= 1.0)
    angle = 0.0;
  else if (ndot <= -1.0)
    angle = CUBIT_PI;
  else 
    angle = acos(ndot);
  if (edot <= 0.0)
  {
    angle = 2.0 * CUBIT_PI - angle;
  }
  return angle;
}

// order edges in list beginning at start_point
// report the endpoint
// return CUBIT_SUCCESS if all edges are connected and ordered successfully
// otherwise return CUBIT_FAILURE, in which case no changes are made
CubitStatus CubitFacetEdge::order_edge_list(DLIList<CubitFacetEdge*> &edge_list,
                                            CubitPoint *start_point,
                                            CubitPoint *&end_point)
{
  int i;
  assert(start_point);

  end_point = NULL;

  // invalid input
  if (0 == edge_list.size())
    return CUBIT_FAILURE;

  // simple case of a single edge - endpoitn
  if (1 == edge_list.size())
  {
    end_point = edge_list.get()->other_point(start_point);
    return end_point ? CUBIT_SUCCESS : CUBIT_FAILURE;
  }

  edge_list.reset();

  // note that a periodic/closed curve will fail
  // we could handle that case here if needed, but we may need more information
  // to know where to start and end the curve
  if (NULL == start_point)
    return CUBIT_FAILURE;

  // put edges in a set for faster searching
  std::set<CubitFacetEdge *> edge_set;
  for (i=0; i<edge_list.size(); i++)
    edge_set.insert(dynamic_cast<CubitFacetEdge*> (edge_list.step_and_get()));

  // a vector for the ordered list
  std::vector<CubitFacetEdge*> ordered_edges;

  // find connected edges from the start point
  CubitPoint *cur_pt = start_point;
  do
  {
    // get edges connected to the current point and find the next edge
    DLIList<CubitFacetEdge *> pt_edges;
    cur_pt->edges(pt_edges);

    std::set<CubitFacetEdge *>::iterator iter_found;
    CubitFacetEdge *cur_edge = NULL;
    for (i=0; i<pt_edges.size() && !cur_edge; i++)
    {
      CubitFacetEdge *tmp_edge = pt_edges.get_and_step();
      iter_found = edge_set.find(tmp_edge);
      if ( iter_found != edge_set.end() )
        cur_edge = tmp_edge;
    }

    // if we don't find a connection before we empty the set
    // then not all the edges are connected  -- return failure
    if (NULL == cur_edge)
      return CUBIT_FAILURE;

    // add the edge to the ordered list
    ordered_edges.push_back( cur_edge );
    edge_set.erase(iter_found);

    cur_pt = cur_edge->other_point(cur_pt);
  }
  while ( edge_set.size());

  if (ordered_edges.size() != edge_list.size())
    return CUBIT_FAILURE;

  // store the edges in the correct order
  edge_list.clean_out();

  std::vector<CubitFacetEdge*>::iterator iter;
  for (iter=ordered_edges.begin(); iter!=ordered_edges.end(); iter++)
    edge_list.append(*iter);

  // get the end point
  CubitFacetEdge *edge1 = edge_list[edge_list.size() - 1];
  CubitFacetEdge *edge2 = edge_list[edge_list.size() - 2];

  end_point = edge1->other_point( edge1->shared_point(edge2) );

  return CUBIT_SUCCESS;
}

CubitPoint *CubitFacetEdge::find_start_point_for_edge_list(DLIList<CubitFacetEdge*> edge_list)
{
  // look for an edge with a point only connected to the one edge in the set
  //
  // TODO - this algorithm could be made more efficient by only checking one endpoint
  //        per edge.  The current implementation should match the order points were
  //        checked in previous functionality.
  //        If speed becomes an issue it could be reimplemented
  //

  if (1 == edge_list.size())
  {
    return edge_list[0]->point(0);
  }

  int i;
  CubitPoint *start_point = NULL;
  for (i=0; i<edge_list.size() && (start_point == NULL); i++)
  {
    CubitFacetEdge *tmp_edge = edge_list.get_and_step();
    DLIList<CubitFacetEdge*> pt_edges;
    tmp_edge->point(0)->edges(pt_edges);

    pt_edges.intersect_unordered(edge_list);
    if (pt_edges.size() == 1)
    {
      start_point = tmp_edge->point(0);
    }
    else
    {
      pt_edges.clean_out();
      tmp_edge->point(1)->edges(pt_edges);
      pt_edges.intersect_unordered(edge_list);
      if (pt_edges.size() == 1)
      {
        start_point = tmp_edge->point(1);
      }
    }
  }
  return start_point;
}

//EOF
