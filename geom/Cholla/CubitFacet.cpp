#include "CubitFacet.hpp"
#include "CubitPoint.hpp"
#include "CubitVector.hpp"
#include "GeometryDefines.h"
#include "CubitPlane.hpp"
#include "CubitFacetEdge.hpp"
#include "CubitFacetData.hpp"
#include "FacetEvalTool.hpp"
#include "CastTo.hpp"
#include "FacetEntity.hpp"
#include "CubitPointData.hpp"
#include "CubitFacetEdgeData.hpp"
#include "TDFacetBoundaryPoint.hpp"
#include "GfxDebug.hpp"

#ifndef CUBIT_MAX 
#define CUBIT_MAX(a,b) (((a) > (b)) ? (a) : (b))
#endif
#ifndef CUBIT_MIN 
#define CUBIT_MIN(a,b) (((a) < (b)) ? (a) : (b))
#endif
#define min3(a,b,c) CUBIT_MIN((CUBIT_MIN((a),(b))),(c))
#define max3(a,b,c) CUBIT_MAX((CUBIT_MAX((a),(b))),(c))

//===========================================================================
//Function Name: CubitFacet
//
//Member Type:  PUBLIC
//Description:  constructor
//===========================================================================
CubitFacet::CubitFacet( )
  : cachedPlane(NULL), facetWeight(0.0), patchCtrlPts(NULL),
    markedFlag(0), isFlat(999), isBackwards(0), toolID(0)

{

}

//===========================================================================
//Function Name: ~CubitFacetData
//
//Member Type:  PUBLIC
//Description:  destructor
//===========================================================================
CubitFacet::~CubitFacet()
{
  if (cachedPlane )
    delete cachedPlane;
  cachedPlane = NULL;
  if (patchCtrlPts)
    delete [] patchCtrlPts;
}


//===========================================================================
//Function Name: normal
//
//Member Type:  PUBLIC
//Description:  return the facet normal
//===========================================================================
CubitVector CubitFacet::normal()
{
  CubitPlane fac_plane = plane();
  CubitVector the_normal = fac_plane.normal();
  if (isBackwards)
    the_normal = -the_normal;
  return the_normal;
}

//-------------------------------------------------------------------------
// Purpose       : set/get the bezier patch internal control points.
//
// Special Notes : Allocate if necessary
//
// Creator       : Steve Owen
//
// Creation Date : 06/28/00
//-------------------------------------------------------------------------
void CubitFacet::set_control_points( CubitVector points[6] )
{
  if (!patchCtrlPts) {
    patchCtrlPts = new CubitVector [6];
  }
  memcpy( patchCtrlPts, points, 6 *sizeof(CubitVector) );
}

void CubitFacet::set_control_points( const double *pt_array )
{
  if (!patchCtrlPts) {
    patchCtrlPts = new CubitVector [6];
  }
  int ii;
  for (ii=0; ii<6; ii++)
  {
    patchCtrlPts[ii].x(pt_array[3*ii]);
    patchCtrlPts[ii].y(pt_array[3*ii+1]);
    patchCtrlPts[ii].z(pt_array[3*ii+2]);
  }
}

void CubitFacet::get_control_points( CubitVector points[6] )
{
  assert(patchCtrlPts != 0);
  memcpy( points, patchCtrlPts, 6 *sizeof(CubitVector) );
}

const CubitPlane &CubitFacet::plane()
{
  if( ! cachedPlane )
  {
    CubitPoint *p0, *p1, *p2;
    points(p0, p1, p2);
    CubitVector v1 = p1->coordinates() - p0->coordinates();
    CubitVector v2 = p2->coordinates() - p0->coordinates();
    cachedPlane = new CubitPlane( v1 * v2, p0->coordinates() );
  }
  return *cachedPlane;
}

void CubitFacet::update_plane()
{
  if ( ! cachedPlane )
    return;
  
  CubitPoint *p0, *p1, *p2;
  points(p0, p1, p2);
  CubitVector v1 = p1->coordinates() - p0->coordinates();
  CubitVector v2 = p2->coordinates() - p0->coordinates();
  CubitVector normal = v1 * v2;
  if (is_backwards()) normal = -normal;
  cachedPlane->set(normal, p0->coordinates() );
}

//-------------------------------------------------------------------------
// Purpose       : determine if the facet is flat (all normal are the same)
//
// Special Notes :
//
// Creator       : Steve Owen
//
// Creation Date : 5/4/01
//-------------------------------------------------------------------------
int CubitFacet::is_flat()
{
  if (isFlat != 999)
  {
    return isFlat;
  }

  // check the edges first.  If on a boundary then assume not flat for now
  // This is to account for any in-plane curvature in the boundary edges

  int ii;
  CubitBoolean on_boundary = CUBIT_FALSE;
  for (ii=0; ii<3 && !on_boundary; ii++)
  {
    CubitPoint *point_ptr = point(ii);
    CubitFacetEdge *edge_ptr = edge(ii);
    TDFacetBoundaryPoint *td = TDFacetBoundaryPoint::get_facet_boundary_point(point_ptr);
    if (td != NULL || (edge_ptr && edge_ptr->num_adj_facets() <= 1))
    {
      on_boundary = CUBIT_TRUE;
    }
  }
  if (on_boundary)
  {
    isFlat = 0;
  }
  else
  {

    CubitPoint *p0, *p1, *p2;
    points(p0, p1, p2);

    CubitVector n0 = p0->normal( this );
    CubitVector n1 = p1->normal( this );
    CubitVector n2 = p2->normal( this );
    double dot0 = n0 % n1;
    double dot1 = n1 % n2;
    double dot2 = n2 % n0;

    double tol = 1.0 - GEOMETRY_RESABS;
    if (fabs(dot0) > tol && fabs(dot1) > tol && fabs(dot2) > tol)
      isFlat = 1;
    else
      isFlat = 0;
  }

  return isFlat;
}

//===========================================================================
//Function Name: evaluate_position
//
//Member Type:  PUBLIC
//Description:  evaluate the facet at a position (projects to facet)
//              eval_normal is NULL if normal not needed
//===========================================================================
CubitStatus CubitFacet::evaluate_position( const CubitVector &start_position,
                                           CubitVector *eval_point,
                                           CubitVector *eval_normal)
{
  CubitStatus stat = CUBIT_SUCCESS;

  if (is_flat())
  {
    if (eval_point != NULL)
      closest_point(start_position, *eval_point);
    if (eval_normal != NULL)
      *eval_normal = normal();
  }
  else  // project to the smooth facet
  {
    // project the position to the planar facet
    CubitVector close_point;
    stat = closest_point(start_position, close_point);

    // get the area coordinates of this point as a starting guess
    CubitVector area_coordinates;
    FacetEvalTool::facet_area_coordinate(this, close_point, area_coordinates);

    // now evaluate the smooth facet (this may alter the area coords)
    CubitVector proj_point;
    CubitBoolean outside;
    double tol = sqrt(area()) * 1.0e-3;
    stat = FacetEvalTool::project_to_facet( this, close_point, area_coordinates,
                                            proj_point, outside, tol );

    if (eval_point != NULL)
    {
      *eval_point = proj_point;
    }
    // compute the smooth normal if required
    if (eval_normal != NULL)
    {
      FacetEvalTool::eval_facet_normal(this, area_coordinates, *eval_normal);
    }
  }

  return stat;
}

//===========================================================================
//Function Name: evaluate
//
//Member Type:  PUBLIC
//Description:  evaluate the facet at area coordinates
//              eval_normal is NULL if normal not needed
//===========================================================================
CubitStatus CubitFacet::evaluate( CubitVector &areacoord,
                                  CubitVector *eval_point,
                                  CubitVector *eval_normal )
{
  if (isBackwards)
  {
    double temp = areacoord.y();
    areacoord.y( areacoord.z() );
    areacoord.z( temp );
  }
  return FacetEvalTool::eval_facet( this, areacoord, eval_point, eval_normal );
}

//===========================================================================
//Function Name: closest_point
//
//Member Type:  PUBLIC
//Description:  return the closest point on plane defined by the facet
//===========================================================================
CubitStatus CubitFacet::closest_point(const CubitVector &point,
                                          CubitVector &closest_point )
{
  CubitPlane fac_plane = plane();
  closest_point = fac_plane.project( point );
  return CUBIT_SUCCESS;
}

//===========================================================================
//Function Name: closest_point_trimmed
//
//Member Type:  PUBLIC
//Description:  return the closest point on the facet to the point (trimmed
//              to its boundary)
//===========================================================================
CubitStatus CubitFacet::closest_point_trimmed( const CubitVector &mypoint,
                                               CubitVector &closest_point,
                                               CubitPoint *&next_edge_p1,
                                               CubitPoint *&next_edge_p2)
{
  CubitVector p1 = point(0)->coordinates();
  CubitVector p2 = point(1)->coordinates();
  CubitVector p3 = point(2)->coordinates();
  //First get the edge vectors.
  p1.x();
  p2.x();

  CubitVector y1 = p2 - p1;
  CubitVector y2 = p3 - p2;
  CubitVector y3 = p1 - p3;
  //Now get the vectors from the point to the vertices of the facet.
  CubitVector w1 = mypoint - p1;
  CubitVector w2 = mypoint - p2;
  CubitVector w3 = mypoint - p3;
  //Now cross the edge vectors with the vectors to the point.  If the point
  //in questionis inside the facet then each of these vectors will be in in
  //the same direction as the normal of this facet.
  CubitVector x1 = y1*w1;
  CubitVector x2 = y2*w2;
  CubitVector x3 = y3*w3;

  //Now take the dot products to help us determine if it is in the triangle.
  CubitVector n = normal();
  double d1 = x1%n;
  double d2 = x2%n;
  double d3 = x3%n;

  //If this is true, then we just take the closest point to this facet.
  if ( d1 >= -GEOMETRY_RESABS && d2 >= -GEOMETRY_RESABS && d3 >= -GEOMETRY_RESABS )
  {
    CubitStatus rv = this->closest_point(mypoint, closest_point);
    if ( rv != CUBIT_SUCCESS )
    {
      PRINT_ERROR("Closest Point Trimmed Error.  Point in facet but can't"
                  " calc. point.\n");
      return CUBIT_FAILURE;
    }
    next_edge_p1 = NULL;
    next_edge_p2 = NULL;
    return CUBIT_SUCCESS;
  }
  CubitBoolean close_p1 = CUBIT_FALSE;
  CubitBoolean close_p2 = CUBIT_FALSE;
  CubitBoolean close_p3 = CUBIT_FALSE;
  double k1, k2, k3;
  //Now with the "d" values, determine which point or edge we
  //should project to.  We have to go through each of these and
  //determine which is the closest.
  if ( d1 < 0.0 )
  {
    double w1_dot_y1 = w1%y1;
    double y1_squared = y1.length_squared();
    if ( y1_squared <= CUBIT_DBL_MIN && y1_squared >= -CUBIT_DBL_MIN )
    {
      PRINT_ERROR("Length of facet edge too small.\n");
      return CUBIT_FAILURE;
    }
    k1 = w1_dot_y1/y1_squared;
    if ( k1 < 0.0 )
      close_p1 = CUBIT_TRUE;
    else if ( k1 > 1.0 )
      close_p2 = CUBIT_TRUE;
    else if ( k1 >= 0.0 && k1 <= 1.0 + GEOMETRY_RESABS )
    {
      //So we know that y1, is the closest edge.
      closest_point = p1 + k1*y1;
      next_edge_p1 = point(0);
      next_edge_p2 = point(1);
      return CUBIT_SUCCESS;
    }
  }
  if ( d2 < 0.0 )
  {
    double w2_dot_y2 = w2%y2;
    double y2_squared = y2.length_squared();
    if ( y2_squared <= CUBIT_DBL_MIN && y2_squared >= -CUBIT_DBL_MIN )
    {
      PRINT_ERROR("Length of facet edge too small.\n");
      return CUBIT_FAILURE;
    }
    k2 = w2_dot_y2/y2_squared;
    if ( k2 < 0.0 )
    {
      close_p2 = CUBIT_TRUE;
    }
    else if ( k2 > 1.0 )
    {
      close_p3 = CUBIT_TRUE;
    }
    else if ( k2 >= 0.0 && k2 <= 1.0 + GEOMETRY_RESABS )
    {
      //So we know that y2, is the closest edge.
      closest_point = p2 + k2*y2;
      next_edge_p1 = point(1);
      next_edge_p2 = point(2);
      return CUBIT_SUCCESS;
    }
  }
  if ( d3 < 0.0 )
  {
    double w3_dot_y3 = w3%y3;
    double y3_squared = y3.length_squared();
    if ( y3_squared <= CUBIT_DBL_MIN && y3_squared >= -CUBIT_DBL_MIN )
    {
      PRINT_ERROR("Length of facet edge too small.\n");
      return CUBIT_FAILURE;
    }
    k3 = w3_dot_y3/y3_squared;
    if ( k3 < 0.0 )
    {
      close_p3 = CUBIT_TRUE;
    }
    else if ( k3 > 1.0 )
    {
      close_p1 = CUBIT_TRUE;
    }
    else if ( k3 >= 0.0 && k3 <= 1.0 + GEOMETRY_RESABS )
    {
      //So we know that y3, is the closest edge.
      closest_point = p3 + k3*y3;
      next_edge_p1 = point(2);
      next_edge_p2 = point(0);
      return CUBIT_SUCCESS;
    }
  }
  //Now we have the distances, and which edges the point is closest to.
  if ( close_p1 && !close_p2 && !close_p3 )
  {
    closest_point = p1 ;
    next_edge_p1 = point(0);
    next_edge_p2 = NULL;
    return CUBIT_SUCCESS;
  }
  else if ( close_p2 && !close_p1 && !close_p3 )
  {
    closest_point = p2;
    next_edge_p1 = point(1);
    next_edge_p2 = NULL;
    return CUBIT_SUCCESS;
  }
  else if ( close_p3 && !close_p1 && !close_p2 )
  {
    closest_point = p3;
    next_edge_p1 = point(2);
    next_edge_p2 = NULL;
    return CUBIT_SUCCESS;
  }
  else if( close_p1 && close_p2 && !close_p3 )
  {
    if( w1.length_squared() < w2.length_squared() )
      closest_point = p1; 
    else
      closest_point = p2;
    return CUBIT_SUCCESS;
  }
  else if( close_p2 && close_p3 && !close_p1 )
  {
    if( w2.length_squared() < w3.length_squared() )
      closest_point = p2; 
    else
      closest_point = p3;
    return CUBIT_SUCCESS;
  }
  else if( close_p1 && close_p3 && !close_p2 )
  {
    if( w1.length_squared() < w3.length_squared() )
      closest_point = p1; 
    else
      closest_point = p3;
    return CUBIT_SUCCESS;
  }

  PRINT_ERROR("Problems finding the closest point to a facet.\n");
  return CUBIT_FAILURE;
}

//===========================================================================
//Function Name: contains
//
//Member Type:  PUBLIC
//Description:  determine if point is contained in facet
//===========================================================================
CubitBoolean CubitFacet::contains( CubitPoint *p1 )
{
  for ( int i = 0; i < 3; i++ )
    if ( point(i) == p1 )
      return CUBIT_TRUE;
  return CUBIT_FALSE;
}

//===========================================================================
//Function Name: shared_facet
//
//Member Type:  PUBLIC
//Description:  Find the other facet that shares these two points.
//              Note: assumes max two facets per edge
//===========================================================================
CubitFacet* CubitFacet::shared_facet( CubitPoint *p1, CubitPoint *p2 )
{
    //Find the other facet that shares these two points.
  int ii;
  DLIList<CubitFacet*> facet_list;
  p1->facets(facet_list);

  for ( ii = facet_list.size(); ii > 0; ii-- )
  {
    CubitFacet *t_facet = facet_list.get_and_step();
    if ( t_facet == this )
      continue;
    if ( t_facet->contains(p2) )
    {
      assert( t_facet->contains(p1) );
      return t_facet;
    }
  }
  return (CubitFacet*)NULL;
}

//===========================================================================
//Function Name: shared_facets
//
//Member Type:  PUBLIC
//Description:  Make a list of all facets adjacent this facet at the edge
//              defined by p1, p2
//===========================================================================
void CubitFacet::shared_facets( CubitPoint *p1, CubitPoint *p2,
                                DLIList <CubitFacet *> &adj_facet_list)
{
    //Find the other facets that share these two points.
  int ii;
  DLIList<CubitFacet*> facet_list;
  p1->facets(facet_list);

  for ( ii = facet_list.size(); ii > 0; ii-- )
  {
    CubitFacet *t_facet = facet_list.get_and_step();
    if ( t_facet == this )
      continue;
    if ( t_facet->contains(p2) )
    {
      assert( t_facet->contains(p1) );
      adj_facet_list.append( t_facet );
    }
  }
}


//===========================================================================
//Function Name: shared_facet_on_surf
//
//Member Type:  PUBLIC
//Description:  Find the other facet that shares these two points and that
//              is has the same tool id.
//              Note: assumes max two facets per edge
//===========================================================================
CubitFacet* CubitFacet::shared_facet_on_surf( CubitPoint *p1, CubitPoint *p2,
                                      int tool_id )
{
    //Find the other facet that shares these two points and that
    // is marked with flag.
  int ii;
  DLIList<CubitFacet*> facet_list;
  p1->facets(facet_list);

  for ( ii = facet_list.size(); ii > 0; ii-- )
  {
    CubitFacet *t_facet = facet_list.get_and_step();
    if ( t_facet == this )
      continue;
    if ( t_facet->contains(p2) )
    {
      if (tool_id == t_facet->tool_id())
      {
        assert( t_facet->contains(p1) );
        return t_facet;
      }
    }
  }
  return (CubitFacet*)NULL;
}

//===========================================================================
//Function Name: center
//
//Member Type:  PUBLIC
//Description:  return the centroid of this facet
//===========================================================================
CubitVector CubitFacet::center()
{
   CubitVector vec( point(0)->coordinates() );
   vec += point(1)->coordinates();
   vec += point(2)->coordinates();
   vec /= 3.0;
   return vec;
}

//===========================================================================
//Function Name: draw
//
//Member Type:  PUBLIC
//Description:  draw the facet
//===========================================================================
void CubitFacet::debug_draw(int color, int flush_it, int draw_uv )
{
  if ( color == -1 )
    color = CUBIT_RED;
  GfxDebug::draw_facet(this, color, draw_uv);
  if ( flush_it )
    GfxDebug::flush();
}

const int CubitFacet::point_edge_conn[30][2] = { {4, 8}, {8, 11}, {11, 13}, {13, 14},
                                            {3, 7}, {7, 10}, {10, 12},
                                            {2, 6}, {6, 9},
                                            {1, 5},
                                            {14, 12}, {12, 9}, {9, 5}, {5, 0},
                                            {13, 10}, {10, 6}, {6, 1},
                                            {11, 7}, {7, 2},
                                            {8, 3},
                                            {0, 1}, {1, 2}, {2, 3}, {3, 4},
                                            {5, 6}, {6, 7}, {7, 8},
                                            {9, 10}, {10, 11},
                                            {12, 13} };
const int CubitFacet::point_facet_conn[16][3] = {
  {0, 1, 5}, {1, 6, 5}, {1, 2, 6}, {2, 7, 6}, {2, 3, 7}, {3, 8, 7}, {3, 4, 8},
  {5, 6, 9}, {6, 10, 9}, {6, 7, 10}, {7, 11, 10}, {7, 8, 11},
  {9, 10, 12}, {10, 13, 12}, {10, 11, 13},
  {12, 13, 14} };
const double CubitFacet::my_points[15][3] = {
  {1.00, 0.00, 0.00},{0.75, 0.25, 0.00},{0.50, 0.50, 0.00},
  {0.25, 0.75, 0.00},{0.00, 1.00, 0.00},{0.75, 0.00, 0.25},
  {0.50, 0.25, 0.25},{0.25, 0.50, 0.25},{0.00, 0.75, 0.25},
  {0.50, 0.00, 0.50},{0.25, 0.25, 0.50},{0.00, 0.50, 0.50},
  {0.25, 0.00, 0.75},{0.00, 0.25, 0.75},{0.00, 0.00, 1.00} };



//===========================================================================
//Function Name: min_diagonal
//
//Member Type:  PUBLIC
//Description:  from the three points find the minimum diagonal.
//===========================================================================
double CubitFacet::min_diagonal()
{
    //from the three points find the minimum diagonal.
  CubitVector p1 = point(0)->coordinates();
  CubitVector p2 = point(1)->coordinates();
  CubitVector p3 = point(2)->coordinates();
  CubitVector temp1 = p2-p1;
  CubitVector mid_side_1 = p1 + temp1/2.0;
  CubitVector temp2 = p3-p2;
  CubitVector mid_side_2 = p2 + temp2/2.0;
  CubitVector temp3 = p1-p3;
  CubitVector mid_side_3 = p3 + temp3/2.0;
  CubitVector diagv_1 = p3 - mid_side_1;
  CubitVector diagv_2 = p2 - mid_side_3;
  CubitVector diagv_3 = p1 - mid_side_2;

  double diag_1 = diagv_1.length();
  double diag_2 = diagv_2.length();
  double diag_3 = diagv_3.length();
  if ( diag_1 >= diag_2 && diag_1 >= diag_3 )
    return diag_1;
  else if ( diag_2 > diag_1 && diag_2 > diag_3 )
    return diag_2;
  return diag_3;
}

//-------------------------------------------------------------------------
// Purpose       : Get the two points of the edge opposite the
//                 passed point.
//
// Special Notes :
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 03/07/00
//-------------------------------------------------------------------------
void CubitFacet::opposite_edge( CubitPoint* thepoint,
                                CubitPoint*& p1, CubitPoint *& p2 )
{
  p1 = p2 = 0;
  int i;

  for( i = 0; i < 3; i++ )
    if( point(i) == thepoint )
      break;

  if( i < 3 ) //point is on this triangle
  {
    p1 = point((i+1)%3);
    p2 = point((i+2)%3);
  }
}

//-------------------------------------------------------------------------
// Purpose       : Local modification functions.
//
// Special Notes :
//
// Creator       : 
//
// Creation Date : 03/25/00
//-------------------------------------------------------------------------
CubitPoint* CubitFacet::split_edge( CubitPoint* /* edge_pt1 */,
                                    CubitPoint* /* edge_pt2 */,
                                    const CubitVector& /* position */ )
{
  // this function should be implemented in child class since it creates
  // new facet and point entities
  assert(1);
  CubitPoint* new_point = NULL;
  return new_point;
}

//-------------------------------------------------------------------------
// Purpose       : insert a point into the facet
//
// Special Notes : create two new facets and return them
//
// Creator       :
//
// Creation Date :
//-------------------------------------------------------------------------
CubitPoint* CubitFacet::insert_point( const CubitVector& /* position */,
                                          CubitFacet*&   /* new_tri1 */,
                                          CubitFacet*&   /* new_tri2 */)
{
  // this function should be implemented in child class since it creates
  // new facet and point entities
  assert(1);
  CubitPoint* new_point = NULL;
  return new_point;
}


//-------------------------------------------------------------------------
// Purpose       : return the index of the other index not used by
//                 points p1 and p2
//
// Special Notes :
//
// Creator       :
//
// Creation Date :
//-------------------------------------------------------------------------
int CubitFacet::other_index( CubitPoint* pt1, CubitPoint* pt2 )
{
  int i;
  for( i = 0; (point(i) != pt1) && (i < 3); i++ );
  if( i == 3 ) return -1;

  int i1 = (i+1)%3;
  int i2 = (i+2)%3;

  if( point(i1) == pt2 ) return i2;
  else if( point(i2) == pt2 ) return i1;
  else return -1;
}

//-------------------------------------------------------------------------
// Purpose       : compute the angle at one of the points on the facet
//
// Special Notes : angle is in radians
//
// Creator       : Steve Owen
//
// Creation Date : 06/28/00
//-------------------------------------------------------------------------
double CubitFacet::angle( CubitPoint *pt )
{
  int ii, index = -1;
  CubitBoolean found=CUBIT_FALSE;
  for (ii=0; ii<3 && !found; ii++) {
    if (pt==this->point(ii)) {
      index = ii;
      found = CUBIT_TRUE;
    }
  }
  if (!found) {
    return 0.0e0;
  }
  CubitVector e0 = this->point((index+1)%3)->coordinates() - pt->coordinates();
  CubitVector e1 = this->point((index+2)%3)->coordinates() - pt->coordinates();
  e0.normalize();
  e1.normalize();
  double myangle;
  double cosangle = e0%e1;
  if (cosangle >= 1.0) {
    myangle = 0.0e0;
  }
  else if( cosangle <= -1.0 ) {
    myangle = CUBIT_PI;
  }
  else {
    myangle = acos(cosangle);
  }
  return myangle;
}

//-------------------------------------------------------------------------
// Purpose       : return the edge on a facet and its corresponding index.
//
// Special Notes : The index of the edge corresponds to the point index on
//                 the facet opposite the edge.
//
// Creator       : Steve Owen
//
// Creation Date : 06/28/00
//-------------------------------------------------------------------------
CubitFacetEdge *CubitFacet::edge_from_pts( CubitPoint *p1,
                                  CubitPoint *p2,
                                  int &edge_index )
{
  if ((point(0) == p1 && point(1) == p2) ||
    point(0) == p2 && point(1) == p1) {
    edge_index = 2;
    return edge(2);
  }
  if ((point(1) == p1 && point(2) == p2) ||
    point(1) == p2 && point(2) == p1) {
    edge_index = 0;
    return edge(0);
  }
  if ((point(2) == p1 && point(0) == p2) ||
    point(2) == p2 && point(0) == p1) {
    edge_index = 1;
    return edge(1);
  }
  edge_index = -1;
  return NULL;
}

//-------------------------------------------------------------------------
// Purpose       : return the index of an edge on a facet given two
//                 points on the facet
//
// Special Notes : The index of the edge corresponds to the point index on
//                 the facet opposite the edge.
//
// Creator       : Steve Owen
//
// Creation Date : 06/28/00
//-------------------------------------------------------------------------
int CubitFacet::edge_index( CubitPoint *p1, CubitPoint *p2, int &sense )
{
  if (point(0) == p1 && point(1) == p2)
  {
    sense = 1;
    return 2;
  }
  if (point(0) == p2 && point(1) == p1)
  {
    sense = -1;
    return 2;
  }
  if (point(1) == p1 && point(2) == p2)
  {
    sense = 1;
    return 0;
  }
  if (point(1) == p2 && point(2) == p1)
  {
    sense = -1;
    return 0;
  }
  if (point(2) == p1 && point(0) == p2)
  {
    sense = 1;
    return 1;
  }
  if (point(2) == p2 && point(0) == p1)
  {
    sense = -1;
    return 1;
  }
  sense = 0;
  return -1;
}

// overloaded function - based on edge pointer
int CubitFacet::edge_index( CubitFacetEdge *theedge )
{
  if (edge(0) == theedge) return 0;
  if (edge(1) == theedge) return 1;
  if (edge(2) == theedge) return 2;
  return -1;
}


//-------------------------------------------------------------------------
// Purpose       : return the index of a point on a facet
//
// Creator       : Steve Owen
//
// Creation Date : 06/28/00
//-------------------------------------------------------------------------
int CubitFacet::point_index( CubitPoint *pt )
{
  if (point(0) == pt) return 0;
  if (point(1) == pt) return 1;
  if (point(2) == pt) return 2;
  return -1;
}

//-------------------------------------------------------------------------
// Purpose       : get the points that define an edge of a facet.
//
// Special Notes : The index of the edge corresponds to the point index on
//                 the facet opposite the edge.  The points are returned
//                 based on the edgeUse order
//
// Creator       : Steve Owen
//
// Creation Date : 06/28/00
//-------------------------------------------------------------------------
void CubitFacet::get_edge_pts( int index, CubitPoint *&p1, CubitPoint *&p2 )
{
  assert(index >=0 && index <=2);

    p1 = point((index+1)%3);
    p2 = point((index+2)%3);
}

//-------------------------------------------------------------------------
// Purpose       : update the bounding box of the facet based on the control
//                 polygon of the facet.
//
// Creator       : Steve Owen
//
// Creation Date : 06/28/00
//-------------------------------------------------------------------------
void CubitFacet::update_bezier_bounding_box( )
{
  int i,j;
  CubitVector ctrl_pts[5], min, max;
  CubitFacetEdge *myedge;
  CubitBox bbox = bounding_box();
  min = bbox.minimum();
  max = bbox.maximum();
  for (i=0; i<3; i++) {
    myedge = edge(i);
    assert(myedge != 0);
    myedge->control_points( this, ctrl_pts );
    for (j=1; j<4; j++) {
      if (ctrl_pts[j].x() < min.x()) min.x( ctrl_pts[j].x() );
      if (ctrl_pts[j].y() < min.y()) min.y( ctrl_pts[j].y() );
      if (ctrl_pts[j].z() < min.z()) min.z( ctrl_pts[j].z() );
      if (ctrl_pts[j].x() > max.x()) max.x( ctrl_pts[j].x() );
      if (ctrl_pts[j].y() > max.y()) max.y( ctrl_pts[j].y() );
      if (ctrl_pts[j].z() > max.z()) max.z( ctrl_pts[j].z() );
    }
  }
  CubitVector patch_ctrl_pts[6];
  get_control_points( patch_ctrl_pts );
  assert(patch_ctrl_pts != 0);
  for (j=0; j<6; j++) {
    if (patch_ctrl_pts[j].x() < min.x()) min.x( patch_ctrl_pts[j].x() );
    if (patch_ctrl_pts[j].y() < min.y()) min.y( patch_ctrl_pts[j].y() );
    if (patch_ctrl_pts[j].z() < min.z()) min.z( patch_ctrl_pts[j].z() );
    if (patch_ctrl_pts[j].x() > max.x()) max.x( patch_ctrl_pts[j].x() );
    if (patch_ctrl_pts[j].y() > max.y()) max.y( patch_ctrl_pts[j].y() );
    if (patch_ctrl_pts[j].z() > max.z()) max.z( patch_ctrl_pts[j].z() );
  }
  bBox.reset(min,max);
}

//-------------------------------------------------------------------------
// Purpose       : reset the bounding box of the facet based on the control
//                 polygon of the facet.
//
// Creator       : Steve Owen
//
// Creation Date : 3/14/02
//-------------------------------------------------------------------------
void CubitFacet::reset_bounding_box( )
{

  CubitPoint *p1 = point (0);
  CubitPoint *p2 = point (1);
  CubitPoint *p3 = point (2);

  // define the bounding box

  if (!patchCtrlPts || is_flat())
  {
    CubitVector bbox_min, bbox_max;
    bbox_min.x(min3(p1->x(),p2->x(),p3->x()));
    bbox_min.y(min3(p1->y(),p2->y(),p3->y()));
    bbox_min.z(min3(p1->z(),p2->z(),p3->z()));
    bbox_max.x(max3(p1->x(),p2->x(),p3->x()));
    bbox_max.y(max3(p1->y(),p2->y(),p3->y()));
    bbox_max.z(max3(p1->z(),p2->z(),p3->z()));
    bBox.reset(bbox_min,bbox_max);
  }
  else
  {
    update_bezier_bounding_box( );
  }
}

//-------------------------------------------------------------------------
// Purpose       : get the next CCW (rhr) facetedge on the facet
//
// Special Notes :
//
// Creator       : Steve Owen
//
// Creation Date : 06/28/00
//-------------------------------------------------------------------------
CubitFacetEdge *CubitFacet::next_edge( CubitFacetEdge *theedge )
{
  int i;
  for (i=0; i<3; i++) {
    if (theedge == edge(i)) {
      return edge((i+1)%3);
    }
  }
  return NULL;
}

//-------------------------------------------------------------------------
// Purpose       : get the prev CCW (rhr) facetedge on the facet
//
// Special Notes :
//
// Creator       : Steve Owen
//
// Creation Date : 06/28/00
//-------------------------------------------------------------------------
CubitFacetEdge *CubitFacet::prev_edge( CubitFacetEdge *theedge )
{
  int i;
  for (i=0; i<3; i++) {
    if (theedge == edge(i)) {
      return edge((i+2)%3);
    }
  }
  return NULL;
}

//-------------------------------------------------------------------------
// Purpose       : reorient the facet
//
// Special Notes :
//
// Creator       : Steve Owen
//
// Creation Date : 06/28/00
//-------------------------------------------------------------------------
void CubitFacet::flip()
{
  // must be implemented in child class
  assert(0);

}

//-------------------------------------------------------------------------
// Purpose       : return the next edge on the triangle at the point
//
// Special Notes :
//
// Creator       : Steve Owen
//
// Creation Date : 4/31/01
//-------------------------------------------------------------------------
CubitFacetEdge *CubitFacet::next_edge_at_point( CubitFacetEdge *edge_ptr,
                                                CubitPoint *point_ptr )
{
  int eidx = edge_index( edge_ptr );
  int pidx = point_index( point_ptr );
  switch(eidx)
  {
    case 0:
      switch(pidx)
      {
        case 1: eidx = 2; break;
        case 2: eidx = 1; break;
        default: eidx = -1; break;
      }
      break;
    case 1:
      switch(pidx)
      {
        case 0: eidx = 2; break;
        case 2: eidx = 0; break;
        default: eidx = -1; break;
      }
      break;
    case 2:
      switch(pidx)
      {
        case 0: eidx = 1; break;
        case 1: eidx = 0; break;
        default: eidx = -1; break;
      }
      break;
    default:
      eidx = -1;
      break;
  }
  if (eidx == -1)
    return (CubitFacetEdge *)NULL;

  return edge( eidx );
}

//======================================================================
// Function: get_edge control_points (PUBLIC)
// Description: return the control points on the facet edges based
//        on its edge use directions
// Author: sjowen
// Date: 5/01
//======================================================================
CubitStatus CubitFacet::get_edge_control_points( CubitVector P[3][5] )
{
  CubitStatus stat;
  CubitFacetEdge *edge_ptr;
  CubitVector ctrl_pts[5];
  int ii, jj;
  for (ii=0; ii<3; ii++) {
    edge_ptr = edge(ii);
    stat = edge_ptr->control_points(this, ctrl_pts);
    if (stat!= CUBIT_SUCCESS)
      return stat;
    for (jj=0; jj<5; jj++)
    {
      P[ii][jj] = ctrl_pts[jj];
    }
  }
  return stat;
}

//======================================================================
// Function: area (PUBLIC)
// Description: compute the area of the facet
// Author: sjowen
// Date: 3/02
//======================================================================
double CubitFacet::area(  )
{
  CubitVector e0(point(0)->coordinates(), point(1)->coordinates());
  CubitVector e1(point(0)->coordinates(), point(2)->coordinates());
  double area = (e0 * e1).length() * 0.5;
  return area;
}

double CubitFacet::aspect_ratio()
{
  static const double normal_coeff = sqrt( 3. ) / 6.;

  // three vectors for each side 
  CubitVector a = point(1)->coordinates() - point(0)->coordinates();
  CubitVector b = point(2)->coordinates() - point(1)->coordinates();
  CubitVector c = point(0)->coordinates() - point(2)->coordinates();
    
  double a1 = a.length();
  double b1 = b.length();
  double c1 = c.length();
 
  double hm = a1 > b1 ? a1 : b1;
  hm = hm > c1 ? hm : c1;

  CubitVector ab = a * b;
  double denominator = ab.length();

  if( denominator < CUBIT_DBL_MIN ) 
    return (double)CUBIT_DBL_MAX;
  else
  {
    double aspect_ratio;
    aspect_ratio = normal_coeff * hm * (a1 + b1 + c1) / denominator;
    
    if( aspect_ratio > 0 )
      return (double) CUBIT_MIN( aspect_ratio, CUBIT_DBL_MAX );
    return (double) CUBIT_MAX( aspect_ratio, -CUBIT_DBL_MAX );
  }
}


//======================================================================
// Function: init_patch (PUBLIC)
// Description: computes the interior control points for the facet and
//              stores them with the facet.  Assumes edge control points
//              have already been computed
// Author: sjowen
// Date: 10/03
//======================================================================
CubitStatus CubitFacet::init_patch(  )
{
  return FacetEvalTool::init_bezier_facet( this );
}

//===========================================================================
//Function Name: add_edge
//
//Member Type:  PUBLIC
//Description:  add an existing edge to a facet.  The edge must contain
//              points that are already part of this facet
//===========================================================================
void CubitFacet::add_edge(CubitFacetEdge *edge)
{
  CubitPoint *p0 = edge->point(0);
  CubitPoint *p1 = edge->point(1);

  // find the points on this facet

  int p0_index = point_index( p0 );
  int p1_index = point_index( p1 );
  assert(p0_index >=0 && p1_index >= 0);
  assert(p0_index != p1_index);

  // add the edge based on the relative index of the points
  // set the edge use based on their order

  int edge_index;
  if ((p0_index+1)%3 == p1_index)
  {
    edge_index = (p1_index + 1)%3;
    assert(this->edge(edge_index) == NULL);
    this->edge( edge, edge_index );
    this->edge_use(1, edge_index);
  }
  else if((p0_index+2)%3 == p1_index)
  {
    edge_index = (p0_index + 1)%3;
    assert(this->edge(edge_index) == NULL);
    this->edge( edge, edge_index );
    this->edge_use(-1, edge_index);
  }
  else
  {
    assert(0);  // shouldn't happen
  }

  // add the facet to the edge

  edge->add_facet( this );
}


CubitPoint *CubitFacet::opposite_point( CubitFacetEdge *edge )
{
  int i;
  for( i = 0; i < 3; i++ )
    if( point(i) != edge->point(0) && point(i) != edge->point(1) )
      return point(i);

  return NULL;
}

CubitVector CubitFacet::update_normal( void )
{
   this->update_plane();
   return normal();
}

 void CubitFacet::unlink_from_children( void )
{
    
    this->point(0)->remove_facet(this);
    this->point(1)->remove_facet(this);
    this->point(2)->remove_facet(this);
    this->edge(0)->remove_facet(this);
    this->edge(1)->remove_facet(this);
    this->edge(2)->remove_facet(this);    
  
}

 CubitFacetEdge* CubitFacet::shared_edge( CubitFacet *cubit_facet )
 {
   for( int i=0; i<3; i++ )
     for( int j=0; j<3; j++ )
       if( edge(i) == cubit_facet->edge(j) )
         return edge(i);       

   return NULL;
 }
