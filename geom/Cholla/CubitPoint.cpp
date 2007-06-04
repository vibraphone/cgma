#include "CubitPoint.hpp"
#include "CubitFacet.hpp"
#include "CubitFacetEdge.hpp"
#include "TDFacetBoundaryPoint.hpp"
#include "CubitTransformMatrix.hpp"
#include "GfxDebug.hpp"
#include "CubitQuadFacet.hpp"
#include "GeometryDefines.h"

double CubitPoint::boxTol = GEOMETRY_RESABS;

//===========================================================================
// Function Name: CubitPoint
//
// Member Type:  PUBLIC
// Description:  constructor
// Author: 
// Date:
//===========================================================================
CubitPoint::CubitPoint() 
    : markedFlag(0), surfNormal(NULL), dCoef(0.0), 
      uVal(0.0), vVal(0.0), sizeVal(0.0),
      surfU(NULL), surfV(NULL),
      coefVector(NULL), isFeature(0)
{
}

//===========================================================================
// Function Name: ~CubitPoint
//
// Member Type:  PUBLIC
// Description:  destructor
// Author: 
// Date:
//===========================================================================
CubitPoint::~CubitPoint()
{
  if (surfNormal) {
    delete surfNormal;
  }
  if (surfU) {
    delete surfU;
  }
  if (surfV) {
    delete surfV;
  }
  if (coefVector) {
    delete coefVector;
  }
}

//===========================================================================
// Function Name: shared_facets
//
// Member Type:  PUBLIC
// Description:  return the facets sharing this point and another point
// Author: 
// Date:
//===========================================================================
void CubitPoint::shared_facets( 
  CubitPoint* other_pt, 
  CubitFacet *& f1, 
  CubitFacet*& f2 )
{
  f1 = f2 = 0;
  DLIList <CubitFacet *> attached_facets;
  facets( attached_facets );
  if( attached_facets.size() > 0 && other_pt )
  {
    for( int i = attached_facets.size(); i > 0; i-- )
    {
      CubitFacet* facet = attached_facets.get_and_step();
      if( facet->contains( other_pt ) )
      {
        //three facets??
        assert( !f2 );
        
        if( f1 ) f2 = facet;
        else f1 = facet;
      }
    }
  }
}
void CubitPoint::shared_facets( 
  CubitPoint* other_pt, DLIList<CubitFacet*>& result_set )
{
  DLIList <CubitFacet *> attached_facets;
  facets( attached_facets );
  if( attached_facets.size() > 0 && other_pt )
  {
    for( int i = attached_facets.size(); i > 0; i-- )
    {
      CubitFacet* facet = attached_facets.get_and_step();
      if( facet->contains( other_pt ) )
      {
        result_set.append( facet );
      }
    }
  }
}


//===========================================================================
// Function Name: shared_edge
//
// Member Type:  PUBLIC
// Description:  return the edge sharing this point and another point and 
//               assumes there is only one edge shared by two points
// Author: chynes
// Date: 6/6/2002
//===========================================================================
CubitFacetEdge * CubitPoint::shared_edge( CubitPoint* other_pt )
{
  
//   CubitFacetEdge *edge = NULL;
//   DLIList <CubitFacetEdge *> attached_edges;
//   edges( attached_edges );
//   if( attached_edges.size() > 0 && other_pt )
//   {
//     for( int i = attached_edges.size(); i > 0; i-- )
//     {
//       edge = attached_edges.get_and_step();
//       if( edge->contains( other_pt ) )
//       {
//         return edge;
//       }
//     }
//   }
//   return edge;
  return get_edge(other_pt);

}


//===========================================================================
// Function Name: normal
//
// Member Type:  PUBLIC
// Description:  return the facet normal with respect to the surface the
//               facet is on
// Author: 
// Date:
//===========================================================================
CubitVector CubitPoint::normal( CubitFacet *facet_ptr )
{
  TDFacetBoundaryPoint *td_bfp =
    TDFacetBoundaryPoint::get_facet_boundary_point( this );
  if (td_bfp == NULL)
  {
    return normal();
  }
  else
  {
    CubitVector norm;
    td_bfp->get_normal( facet_ptr, norm );
    return norm;
  }
}

//===========================================================================
// Function Name: normal
//
// Member Type:  PUBLIC
// Description:  return the facet normal with respect to the surface the
//               quad facet is on
// Author: 
// Date:
//===========================================================================
CubitVector CubitPoint::normal( CubitQuadFacet *qfacet_ptr )
{
  CubitFacet *facet_ptr = qfacet_ptr->get_tri_facet_at_point( this );
  return normal( facet_ptr );
}


//===========================================================================
// Function Name: normal
//
// Member Type:  PUBLIC
// Description:  return the facet normal with respect to the surface the
//               facet is on
// Author: 
// Date:
//===========================================================================
CubitVector CubitPoint::normal( CubitFacetEdge *edge_ptr )
{
  TDFacetBoundaryPoint *td_bfp =
    TDFacetBoundaryPoint::get_facet_boundary_point( this );
  if (td_bfp == NULL)
  {
    return normal();
  }
  else
  {
    CubitVector norm;
    td_bfp->get_normal( edge_ptr, norm );
    return norm;
  }
}

//===========================================================================
// Function Name: tangent
//
// Member Type:  PUBLIC
// Description:  return tangent with respect to the edge
// Notes      :  min_dot is the cosine of the feature angle.  Tangent vector
//               will use the feature angle to determine tangent
// Author:       sjowen
// Date:         10/28/2002
//===========================================================================
CubitVector CubitPoint::tangent( CubitFacetEdge *edge_ptr, 
                                 double min_dot )
{
  CubitPoint *p0 = edge_ptr->point( 0 );
  CubitPoint *p1 = edge_ptr->point( 1 );
  int ii;

  assert( p0 == this || p1 == this ); // the point isn't on the edge

    // if this isn't a feature edge, just return the tangent vector of 
    // the edge.  Otherwise compute the tangent based on neighboring 
    // feature edges

  CubitVector pt_tangent;
  if (!edge_ptr->is_feature())
  {
    CubitVector tmp_vec = coordinates();
    
    edge_ptr->edge_tangent( tmp_vec, pt_tangent );
  }
  else
  {
    // compute tangent for feature edge at previous
    if (p0 == this)
    {
      CubitFacetEdge *prev_edge;
      DLIList <CubitFacetEdge *>feature_edge_list;
      next_feature_edges( edge_ptr, feature_edge_list );

      // average the edges that meet the min_dot criteria
      CubitVector e1 = p1->coordinates() - p0->coordinates();         
      pt_tangent = e1;
      e1.normalize();
      for (ii=0; ii<feature_edge_list.size(); ii++)
      {
        prev_edge = feature_edge_list.get_and_step();
        CubitPoint *p2 = prev_edge->other_point( p0 );          
        CubitVector e0 = p0->coordinates() - p2->coordinates();     
        e0.normalize();
        if (e0 % e1 >= min_dot)
        {
          pt_tangent += (p0->coordinates() - p2->coordinates());
        }
      }
      if (feature_edge_list.size() == 0)
        pt_tangent = e1;
      else
        pt_tangent.normalize();
    }

    // compute tangent for feature edge at next
    else if (p1 == this)
    {
      CubitFacetEdge *next_edge;
      DLIList <CubitFacetEdge *>feature_edge_list;
      next_feature_edges( edge_ptr, feature_edge_list );

      // average the edges that meet the min_dot criteria
      CubitVector e1 = p1->coordinates() - p0->coordinates();         
      pt_tangent = e1;
      e1.normalize();
      for (ii=0; ii<feature_edge_list.size(); ii++)
      {
        next_edge = feature_edge_list.get_and_step();
        CubitPoint *p2 = next_edge->other_point( p1 );          
        CubitVector e0 = p2->coordinates() - p1->coordinates();     
        e0.normalize();
        if (e0 % e1 >= min_dot)
        {
          pt_tangent += (p2->coordinates() - p1->coordinates());
        }
      }
      if (feature_edge_list.size() == 0)
        pt_tangent = e1;
      else
        pt_tangent.normalize();
    }
  } 
  return pt_tangent;
}

//===========================================================================
//Function Name: next_feature_edges
//
//Member Type:  PRIVATE
//Descriptoin:  given a facet boundary edge and this point, get a list
//              of the next fetaure edges at this point
//===========================================================================
void CubitPoint::next_feature_edges( 
  CubitFacetEdge *this_edge_ptr,
  DLIList <CubitFacetEdge *> feature_edge_list )
{
  //CubitFacetEdge *next_edge_ptr = NULL;

  DLIList<CubitFacetEdge*> edge_list;
  edges( edge_list );
  int ii;

  CubitFacetEdge *edge_ptr = NULL;
  for (ii=0; ii<edge_list.size(); ii++)
  {
    edge_ptr = edge_list.get_and_step();
    if (edge_ptr != this_edge_ptr)
    {
      if (edge_ptr->is_feature())
      {
        feature_edge_list.append(edge_ptr);
      }
    }
  }
}

//===========================================================================
// Function Name: project_to_tangent_plane
//
// Member Type:  PUBLIC
// Descriptoin:  Project a point to the tangent plane defined at the CubitPoint
// Author: sjowen
// Date: 06/28/00
//===========================================================================
CubitVector CubitPoint::project_to_tangent_plane( CubitVector &pt )
{
  CubitVector surf_normal = normal();
  double dist = (surf_normal)%pt + d_coef();
  CubitVector point_on_plane( pt.x() - surf_normal.x() * dist,
                              pt.y() - surf_normal.y() * dist,
                              pt.z() - surf_normal.z() * dist );
  return point_on_plane;
}

//===========================================================================
// Function Name: adjacent_points
//
// Member Type:  PUBLIC
// Description:  return array of points sharing adjacent facets
// Author: sjowen
// Date: 06/28/00
//===========================================================================
void CubitPoint::adjacent_points( CubitPoint **adj_points,
                                  int &num_adj_points )
{
  int i, j, k, index = -1, nextindex = -1;
  CubitBoolean found;
  CubitFacet *facet;
  num_adj_points = 0;
  DLIList <CubitFacet *> attached_facets;
  facets( attached_facets );
  for(i=0; i<attached_facets.size(); i++) {
    facet = attached_facets.get_and_step();
    found = CUBIT_FALSE;
    for (j=0; j<3 && !found; j++) {
      if (facet->point(j) == this) {
        index = (j+1)%3;
        nextindex = (j+2)%3;
        found = CUBIT_TRUE;
      }
    }
    if (found) {
      found = CUBIT_FALSE;
      adj_points[num_adj_points++] = facet->point(index);
      for (k=0; k<num_adj_points-1 && !found; k++) {
        if(adj_points[k] == facet->point(nextindex)){
          found = CUBIT_TRUE;
        }
      }
      if (!found) {
        adj_points[num_adj_points++] = facet->point(nextindex);
      }
    }
  }
}
void CubitPoint::adjacent_points( DLIList<CubitPoint*>& result )
{
  DLIList <CubitFacet *> attached_facets;
  facets( attached_facets );
  for( int i = attached_facets.size(); i--; )
  {
    CubitPoint *pt1, *pt2;
    attached_facets.get_and_step()->opposite_edge( this, pt1, pt2 );
    result.append_unique(pt1);
    result.append_unique(pt2);  
  }
}

//===========================================================================
// Function Name: define_tangent_vectors
//
// Member Type:  PUBLIC
// Description:  define the surface tangent vectors at a point
// Author: sjowen
// Date: 06/28/00
//===========================================================================
void CubitPoint::define_tangent_vectors()
{
  // define orthogonal vectors to the normal that are tangent to the
  // the surface.  Note that du and dv are not defined in any global
  // parametric space - they are only defined locally.  Their directions
  // are defined arbitrarily in the tangent plane by taking the smallest
  // components of the normal vector and setting them to "1" for the
  // du vector and then solving for the other component so that the
  // dot product of du and the normal will be zero.  dv is just the cross
  // product of the normal and du

  CubitVector absnorm, duvec, dvvec;
  CubitVector surf_normal = normal();
  absnorm.x( fabs(surf_normal.x()) );
  absnorm.y( fabs(surf_normal.y()) );
  absnorm.z( fabs(surf_normal.z()) );
  if (absnorm.x() >= absnorm.y() && absnorm.x() >= absnorm.z()) {
    duvec.x( (-surf_normal.y() - surf_normal.z()) / surf_normal.x() );
    duvec.y( 1.0e0 );
    duvec.z( 1.0e0 );
  }
  else if (absnorm.y() >= absnorm.z() ) {
    duvec.x( 1.0e0 );
    duvec.y( (-surf_normal.x() - surf_normal.z()) / surf_normal.y() );
    duvec.z( 1.0e0 );
  }
  else {
    duvec.x( 1.0e0 );
    duvec.y( 1.0e0 );
    duvec.z( (-surf_normal.x() - surf_normal.y()) / surf_normal.z() );
  }
  duvec.normalize();
  dvvec = surf_normal * duvec;

//  //sjowen debug
//  CubitVector test = dvvec * surf_normal;
//  double dot = test % duvec;
//  if (dot < 0.999999999999) {
//    PRINT_ERROR("Error in define_tangent_vectors");
//  }

  du( duvec );
  dv( dvvec );
}

//===========================================================================
// Function Name: transform_to_local
//
// Member Type:  PUBLIC
// Description:  transform a vector from global to local system
// Author: sjowen
// Date: 06/28/00
//===========================================================================
void CubitPoint::transform_to_local( CubitVector &glob_vec,
                                     CubitVector &loc_vec )
{
  // Translate to local origin at point

  CubitVector vect = glob_vec - this->coordinates();

  // Multiply by transpose (inverse) of transformation vector */

  loc_vec.x( vect % du() );
  loc_vec.y( vect % dv() );
  loc_vec.z( vect % normal() );
}

//===========================================================================
// Function Name: transform_to_global
//
// Member Type:  PUBLIC
// Description:  transform a vector from local to global system
// Author: sjowen
// Date: 06/28/00
//===========================================================================
void CubitPoint::transform_to_global( CubitVector &loc_vec,
                                      CubitVector &glob_vec )
{
  // Multiply by transformation matrix

  CubitVector vect;
  CubitVector surf_u = du();
  CubitVector surf_v = dv();
  CubitVector surf_normal = normal();
  vect.x( loc_vec.x() * surf_u.x ()+
          loc_vec.y() * surf_v.x() +
          loc_vec.z() * surf_normal.x() );
  vect.y( loc_vec.x() * surf_u.y ()+
          loc_vec.y() * surf_v.y() +
          loc_vec.z() * surf_normal.y() );
  vect.z( loc_vec.x() * surf_u.z ()+
          loc_vec.y() * surf_v.z() +
          loc_vec.z() * surf_normal.z() );

  // Translate from origin

  glob_vec = vect + this->coordinates();
}

//===========================================================================
// Function Name: get_parents
//
// Member Type:  PUBLIC
// Description:  get the list of attached edges
// Author: sjowen
// Date: 06/28/00
//===========================================================================
void CubitPoint::get_parents( DLIList<FacetEntity*> &facet_list )
{
  DLIList<CubitFacetEdge*> edge_list;
  edges( edge_list );
  int ii;
  for (ii=0; ii<edge_list.size(); ii++)
    facet_list.append( edge_list.get_and_step() );
}

//===========================================================================
// Function Name: draw
//
// Member Type:  PUBLIC
// Description:  debug drawing
// Author: sjowen
// Date: 5/01
//===========================================================================
void CubitPoint::debug_draw( int color, int flush, int /*draw_uv*/ )
{
  if ( color == -1 )
    color = CUBIT_YELLOW;
  CubitVector vec = this->coordinates();
  GfxDebug::draw_point(vec, color);
  if (flush)
    GfxDebug::flush();
}

//===========================================================================
// Function Name: normal
//
// Member Type:  PUBLIC
// Description:  set the normal of the surface at the point. Allocate a 
//               new normal vector if necessary
// Author: 
// Date:
//===========================================================================
void CubitPoint::compute_avg_normal()
{
  int j;
  DLIList<CubitFacet*> adj_facet_list;
  facets(adj_facet_list);
  if (adj_facet_list.size() > 0) {
    CubitVector avg_normal(0.0e0, 0.0e0, 0.0e0);
    double totangle = 0.0e0;

      // weight the normal by the spanning angle at the point

    for (j = 0; j < adj_facet_list.size(); j++)
    {
      CubitFacet* facet = adj_facet_list.get_and_step();
      double angle = facet->angle( this );
      facet->weight( angle );
      totangle += angle;
    }
    for (j = 0; j < adj_facet_list.size(); j++)
    {
      CubitFacet* facet = adj_facet_list.get_and_step();
      CubitVector normal = facet->normal();
      normal.normalize();
      avg_normal += (facet->weight() / totangle) * normal;
    }
    avg_normal.normalize();
    if(!surfNormal) {
      surfNormal = new CubitVector ( avg_normal );
    }
    else
    {
      *surfNormal = avg_normal;
    }
    dCoef = -(this->coordinates()%avg_normal);
  }
}

//===========================================================================
// Function Name: facets_on_surf
//
// Member Type:  PUBLIC
// Description:  return facets adjacent this point that are on the given 
//               surface
//               surf_id is the FacetEvalTool ToolID
// Author:       sjowen
// Date:         6/26/01
//===========================================================================
CubitBox CubitPoint::bounding_box(  )
{
  CubitVector ptmin( coordinates().x() - boxTol,
                     coordinates().y() - boxTol,
                     coordinates().z() - boxTol );
  CubitVector ptmax( coordinates().x() + boxTol,
                     coordinates().y() + boxTol,
                     coordinates().z() + boxTol );
  CubitBox ptbox( ptmin, ptmax );
  return ptbox;
}

//===========================================================================
// Function Name: facets_on_surf
//
// Member Type:  PUBLIC
// Description:  return facets adjacent this point that are on the given 
//               surface
//               surf_id is the FacetEvalTool ToolID
// Author:       sjowen
// Date:         6/26/01
//===========================================================================
void CubitPoint::facets_on_surf( int surf_id, 
                                 DLIList<CubitFacet *> &facet_list,
                                 CubitBoolean &on_internal_boundary )
{
  DLIList<CubitFacet *>all_facets;
  facets( all_facets );
  int ii;
  CubitFacet *facet;
  for (ii=0; ii<all_facets.size(); ii++)
  {
    facet = all_facets.get_and_step();
    if (facet->tool_id() == surf_id)
    {
      facet_list.append( facet );
    }
  } 
  if (facet_list.size() == 0)
  {
    on_internal_boundary = CUBIT_FALSE;
  }
  else if (all_facets.size() == facet_list.size())
  {
    on_internal_boundary = CUBIT_FALSE;
  }
  else
  {
    on_internal_boundary = CUBIT_TRUE;
  }
}

//===========================================================================
// Function Name: get_uv
//
// Member Type:  PUBLIC
// Description:  return the u-v coordinates for the surface the facet is on
//               - assumes facet is adjacent to this point
// Author:       sjowen
// Date:         6/26/01
//===========================================================================
CubitStatus CubitPoint::get_uv( CubitFacet *facet, double &u, double &v )
{
  CubitStatus stat;
  TDFacetBoundaryPoint *td_bfp = 
    TDFacetBoundaryPoint::get_facet_boundary_point( this );
  if (!td_bfp)
  {
    u = uVal;
    v = vVal;
    stat = CUBIT_SUCCESS;
  }
  else
  {
    stat = td_bfp->get_uv( facet, u, v );
  }
  return stat;
}

//===========================================================================
// Function Name: get_uvs
//
// Member Type:  PUBLIC
// Description:  return the u-v coordinates and size for the surface the facet is on
//               - assumes facet is adjacent to this point
// Author:       chynes
// Date:         7/16/02
//===========================================================================
CubitStatus CubitPoint::get_uvs( CubitFacet *facet, double &u, double &v, double &s )
{
  CubitStatus stat;
  TDFacetBoundaryPoint *td_bfp = 
    TDFacetBoundaryPoint::get_facet_boundary_point( this );
  if (!td_bfp)
  {
    u = uVal;
    v = vVal;
	s = sizeVal;
    stat = CUBIT_SUCCESS;
  }
  else
  {
    stat = td_bfp->get_uvs( facet, u, v, s );
  }
  return stat;
}
//===========================================================================
// Function Name: merge_points
//
// Member Type:  PUBLIC
// Description:  merge two points
// Author: sjowen
// Date: 9/18/01
//===========================================================================
CubitStatus CubitPoint::merge_points( CubitPoint *  /*cp*/, CubitBoolean /* keep_point */ )
{
  // this virtual function must be defined in the inheriting class if you get here
  assert(0);
  return CUBIT_FAILURE;
}

//===========================================================================
// Function Name: get_edge
//
// Member Type:  PUBLIC
// Description:  return the CubitFacetEdge between the two points (if there is one)
// Author: sjowen
// Date: 9/25/01
//===========================================================================
CubitFacetEdge *CubitPoint::get_edge( CubitPoint *other_point )
{
  DLIList<CubitFacetEdge *>edge_list;
  edges(edge_list);
  int ii;
  CubitFacetEdge *edge_ptr;
  for (ii=0; ii<edge_list.size(); ii++)
  {
    edge_ptr = edge_list.get_and_step();
    if (edge_ptr->other_point( this ) == other_point)
      return edge_ptr;
  }
  return (CubitFacetEdge *)NULL;
}

//===========================================================================
// Function Name: transform
//
// Member Type:  PUBLIC
// Description:  transform the location of the point
// Author: sjowen
// Date: 3/16/02
//===========================================================================
void CubitPoint::transform(CubitTransformMatrix &tfmat)
{
  CubitVector loc;
  loc = tfmat * coordinates();
  set( loc );
}

//===========================================================================
// Function Name: rotate_normal
//
// Member Type:  PUBLIC
// Description:  transform the location of the point
// Author: sjowen
// Date: 3/16/02
//===========================================================================
void CubitPoint::rotate_normal(CubitTransformMatrix &rotmat)
{
  if (surfNormal)
  {
    *surfNormal = rotmat * (*surfNormal); 
  }
  TDFacetBoundaryPoint *td = TDFacetBoundaryPoint::get_facet_boundary_point( this );
  if (td)
  {
    td->rotate_normal(rotmat);
  }
}

//===========================================================================
// Function Name: check_inverted_facets
//
// Member Type:  PUBLIC
// Description:  check if moving a point will invert facets
// Author: jakraft
// Date: 05/09/03
//===========================================================================
CubitStatus CubitPoint::check_inverted_facets( const CubitVector& pos )
{
  DLIList<CubitFacet*> facets;
  this->facets(facets);
  while( facets.size() )
  {
    CubitFacet* facet = facets.pop();
    int index = facet->point_index(this);
    CubitVector corner = facet->point((index+1)%3)->coordinates();
    CubitVector opposite_edge = facet->point((index+2)%3)->coordinates();
    opposite_edge -= corner;
    CubitVector old_edge = corner - coordinates();
    CubitVector new_edge = corner - pos;
    old_edge *= opposite_edge;
    new_edge *= opposite_edge;
    if ( (old_edge % new_edge) <= 0.0 )
      return CUBIT_FAILURE;
  }
  return CUBIT_SUCCESS;
}

//EOF

