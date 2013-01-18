//
// File: CubitQuadFacet.cpp
//
// Owner: sjowen
//
#include "CubitDefines.h"
#include "CubitQuadFacet.hpp"
#include "CubitFacetEdge.hpp"
#include "CubitPoint.hpp"
#include "CubitFacet.hpp"
#include "CubitVector.hpp"
#include "FacetEvalTool.hpp"

#define DETERM(p1,q1,p2,q2,p3,q3) ((q3)*((p2)-(p1)) + (q2)*((p1)-(p3)) + (q1)*((p3)-(p2)))

//===========================================================================
//  Function: CubitQuadFacet
//  Purpose:  constructor
//  Date:     4/11/01
//  Author:   sjowen
//===========================================================================
CubitQuadFacet::CubitQuadFacet()
{
}

//===========================================================================
//  Function: CubitQuadFacet
//  Purpose:  destructor
//  Date:     4/11/01
//  Author:   sjowen
//===========================================================================
CubitQuadFacet::~CubitQuadFacet()
{
}

//===========================================================================
//  Function: evaluate
//  Purpose:  evaluate the quad at a specified u,v coordinate
//  Notes:    Uses the interpolation order previously set up for the
//            facet-based surface.  The u,v system is based on the point
//            order and is defined as follows:
//
//            [3] =  (-1,1)  *============*  [2] = (1,1)
//                           |            |  
//                           |            | 
//                           |            |
//                           |            |
//                           |            |
//            [0] = (-1,-1)  *============*  [1] = (1,-1)
//
//  Date:     4/11/01
//  Author:   sjowen
//===========================================================================
CubitStatus CubitQuadFacet::evaluate( double u, double v, 
                                      CubitVector *eval_point,
                                      CubitVector *eval_normal,
                                      CubitVector *eval_du,
                                      CubitVector *eval_dv)
{
  CubitVector normal;
  CubitStatus stat = CUBIT_SUCCESS;
  int tri_index;
  double A, B, C;
  map_to_tri_system( u, v, tri_index, A, B, C );
  CubitVector ac(A,B,C);
  CubitFacet *tri_facet = get_tri_facet( tri_index );
  CubitVector point;
  stat = tri_facet->evaluate( ac, &point, eval_normal );
  if (!stat)
    return stat;
  if (eval_point != NULL)
    *eval_point = point;
  if (eval_du != NULL || eval_dv != NULL)
  {
    CubitVector du, dv;
    stat = eval_derivatives( u, v, point, du, dv );
    normal = du * dv;
    normal.normalize();
    if (eval_du != NULL)
      *eval_du = du;
    if (eval_dv != NULL)
      *eval_dv = dv;
    if (eval_normal != NULL)
      *eval_normal = normal;
  }
  return stat;
}

//===========================================================================
//  Function: eval_derivatives
//  Purpose:  evaluate derivatives on quad with respect to the u,v system
//            on the quad
//  Date:     5/01
//  Author:   sjowen
//===========================================================================
CubitStatus CubitQuadFacet::eval_derivatives( double u, double v,
                                              CubitVector &point,
                                              CubitVector &du,
                                              CubitVector &dv )
{
  // do a forward, backward or central difference depending on where the
  // u, v is located on the quad

  CubitStatus stat = CUBIT_SUCCESS;
  double incr = 0.001;
  double min = -1.0;
  double max = 1.0;
  double uf, vf, ub, vb;
  CubitVector pf, pb;
  if (u < min + incr)
  {
    uf = u + incr;
    stat = evaluate(uf, v, &pf);
    du = ( pf - point ) / incr;
  }
  else if (u > max - incr)
  {
    ub = u - incr;
    stat = evaluate( ub, v, &pb);
    du = ( point - pb ) / incr;
  }
  else
  {
    uf = u + incr;
    ub = u - incr;
    stat = evaluate( uf, v, &pf );
    stat = evaluate( ub, v, &pb );
    du = ( pf - pb ) / (2.0 * incr);
  }
  if (v < min + incr)
  {
    vf = v + incr;
    stat = evaluate(u, vf, &pf);
    dv = ( pf - point ) / incr;
  }
  else if (v > max - incr)
  {
    vb = v - incr;
    stat = evaluate( u, vb, &pb);
    dv = ( point - pb ) / incr;
  }
  else
  {
    vf = v + incr;
    vb = v - incr;
    stat = evaluate( u, vf, &pf );
    stat = evaluate( u, vb, &pb );
    dv = ( pf - pb ) / (2.0 * incr);
  }
  return stat;
}

//===========================================================================
//  Function: debug_draw
//  Purpose:  draw the facets
//  Date:     4/11/01
//  Author:   sjowen
//===========================================================================
void CubitQuadFacet::debug_draw( int color, int flush_it, int draw_uv  )
{
  CubitFacet *tri_facet = get_tri_facet( 0 );
  tri_facet->debug_draw(color, flush_it, draw_uv);
  tri_facet = get_tri_facet( 1 );
  tri_facet->debug_draw(color, flush_it, draw_uv);
}

//===========================================================================
//  Function: map_to_tri_system
//  Purpose:  given a u,v coordinate, return the triangle and area coordinate
//  Date:     4/11/01
//  Author:   sjowen
//===========================================================================
void CubitQuadFacet::map_to_tri_system(double u, double v,
                                       int &tri_index,
                                       double &A, double &B, double &C)
{
  double trivert[3][2];
  double areacoords[2][3];
  int ii, jj, j0, j1;
  double small_ac[2];
  small_ac[0] = small_ac[1] = 1.0;

  // loop through both tris

  for (ii=0; ii<2; ii++)
  {
 
    // set up the (u,v) vertex coordinates for the tri

    for (jj=0; jj<3; jj++)
    {
      switch(tri_to_quad_index(ii,jj))
      {
      case 0: 
        trivert[jj][0] = -1.0;
        trivert[jj][1] = -1.0;
        break;
      case 1:
        trivert[jj][0] = 1.0;
        trivert[jj][1] = -1.0;
        break;
      case 2:
        trivert[jj][0] = 1.0;
        trivert[jj][1] = 1.0;
        break;
      case 3:
        trivert[jj][0] = -1.0;
        trivert[jj][1] = 1.0;
        break;
      }
    }

    // compute the area coordinates

    for (jj=0; jj<3; jj++)
    {
      j0 = (jj+1)%3;
      j1 = (jj+2)%3;
      areacoords[ii][jj] = DETERM( trivert[j0][0], trivert[j0][1],
                              trivert[j1][0], trivert[j1][1], u, v );
      if (small_ac[ii] > areacoords[ii][jj])
      {
        small_ac[ii] = areacoords[ii][jj];
      }
    }
  }
  tri_index = (small_ac[0] > small_ac[1]) ? 0 : 1;
  A = areacoords[tri_index][0];
  B = areacoords[tri_index][1];
  C = areacoords[tri_index][2];
  double sum = A + B + C;
  A /= sum;
  B /= sum;
  C /= sum;
}

//===========================================================================
//  Function: get_control_points
//  Purpose:  retreive the control points from the quad
//  Date:     4/11/01
//  Author:   sjowen
//===========================================================================
void CubitQuadFacet::get_control_points( double *ctrl_pts )
{
  // assumes 6 control points per triangle and 3 along the diagonal edge

  int ii, jj;
  int ipt = 0;
  int index = 0;
  int numtri = 2;
  int num_ctrl_pts_per_tri = 6;
  for (ii=0; ii<numtri; ii++)
  {
    CubitFacet *tri_facet = get_tri_facet( ii );
    CubitVector *tri_ctrl_pts = tri_facet->control_points();
    assert(tri_ctrl_pts != NULL);
    for (jj=0; jj<num_ctrl_pts_per_tri; jj++)
    {
      index = ipt * 3;
      ctrl_pts[index] = tri_ctrl_pts[jj].x();
      ctrl_pts[index+1] = tri_ctrl_pts[jj].y();
      ctrl_pts[index+2] = tri_ctrl_pts[jj].z();
      ipt++;
    }
  }

  int num_ctrl_pts_per_edge = 3;
  CubitVector *edge_ctrl_pts = NULL;
  CubitFacetEdge *edge_ptr = this->point(0)->shared_edge( this->point(2) );
  assert(edge_ptr != NULL);
  edge_ctrl_pts = edge_ptr->control_points();
  for (ii=0; ii<num_ctrl_pts_per_edge; ii++)
  {
    index = ipt * 3;
    ctrl_pts[index] = edge_ctrl_pts[ii].x();
    ctrl_pts[index+1] = edge_ctrl_pts[ii].y();
    ctrl_pts[index+2] = edge_ctrl_pts[ii].z();
    ipt++;
  }
}

//===========================================================================
//  Function: set_control_points
//  Purpose:  set the control points into the quad
//  Date:     4/11/01
//  Author:   sjowen
//===========================================================================
void CubitQuadFacet::set_control_points( double *ctrl_pts )
{
  // assumes 6 control points per triangle and 3 along the diagonal edge

  int ii;
  int index = 0;
  int numtri = 2;
  int num_ctrl_pts_per_tri = 6;
  for (ii=0; ii<numtri; ii++)
  {    
    CubitFacet *tri_facet = get_tri_facet( ii );
    index = 3 * ii * num_ctrl_pts_per_tri;
    tri_facet->set_control_points( &(ctrl_pts[index]) );
  }

  CubitFacetEdge *edge_ptr = this->point(0)->shared_edge( this->point(2) );
  assert(edge_ptr != NULL);
  index = 3 * numtri * num_ctrl_pts_per_tri;
  edge_ptr->set_control_points( &(ctrl_pts[index]) );
}

//===========================================================================
//  Function: points
//  Purpose:  get the points on the quad.  Virtual to FacetEntity
//  Date:     8/1/03
//  Author:   sjowen
//===========================================================================
void CubitQuadFacet::points(DLIList<CubitPoint*> &point_list ) 
{ 
  for ( int i = 0; i < 4; i++ ) 
    point_list.append(this->point(i)); 
}


//===========================================================================
//  Function: facets
//  Purpose:  append the two facets comprising this quad. Virtual to FacetEntity
//  Date:     8/1/03
//  Author:   sjowen
//===========================================================================
void CubitQuadFacet::facets(DLIList<CubitFacet*> &facet_list ) 
{ 
  facet_list.append( this->get_tri_facet(0) );
  facet_list.append( this->get_tri_facet(1) );
}


//===========================================================================
//  Function: edges
//  Purpose:  append the edges of this quad. Virtual to FacetEntity
//  Date:     8/1/03
//  Author:   sjowen
//===========================================================================
void CubitQuadFacet::edges(DLIList<CubitFacetEdge*> &edge_list ) 
{ 
  for ( int i = 0; i < 4; i++ ) 
    edge_list.append(this->edge(i)); 
} 

//===========================================================================
//  Function: init_patch
//  Purpose:   computes the interior control points for the quad facet and
//             stores them with the quad facet.  Assumes edge control points
//             on the four boundary edges have already been computed
//  Date:     10/03
//  Author:   sjowen
//===========================================================================
CubitStatus CubitQuadFacet::init_patch(  )
{
  // compute control point for the internal diagonal edge

  CubitFacetEdge *edge_ptr = this->point(0)->shared_edge( this->point(2) );
  assert(edge_ptr != NULL);
  CubitPoint *pt0 = edge_ptr->point(0);
  CubitPoint *pt1 = edge_ptr->point(1);
  CubitVector P0 = pt0->coordinates();
  CubitVector P1 = pt1->coordinates();
  CubitVector N0 = pt0->normal( edge_ptr );
  CubitVector N1 = pt1->normal( edge_ptr );
  CubitVector T0 = P1 - P0;
  T0.normalize();
  CubitVector T1 = T0;
  CubitVector Pi[3];
  CubitStatus stat = FacetEvalTool::init_edge_control_points( P0, P1, N0, N1, T0, T1, Pi);
  if (stat == CUBIT_FAILURE)
    return stat;
  edge_ptr->control_points( Pi, 4 );

  // compute control points on the triangles

  stat = this->get_tri_facet(0)->init_patch();
  if (stat != CUBIT_FAILURE)
    stat = this->get_tri_facet(1)->init_patch();
  return stat;
}



//EOF
