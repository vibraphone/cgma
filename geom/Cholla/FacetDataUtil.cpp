//- Class:       FacetDataUtil
//- Description: static library of general functions for querying and/or
//-              modifying facet entities
//- Owner:       Steve Owen
//- Checked by:
//- Version:

#include "FacetDataUtil.hpp"
#include "FacetEvalTool.hpp"
#include "CubitFacet.hpp"
#include "CubitFacetEdge.hpp"
#include "CubitPoint.hpp"
#include "CubitPointData.hpp"
#include "CubitFacetData.hpp"
#include "CubitFacetEdgeData.hpp"
#include "CubitQuadFacet.hpp"
#include "GeometryDefines.h"
#include "debug.hpp"
#include "GfxDebug.hpp"
#include "KDDTree.hpp"

//===========================================================================
//Function Name: edges_by_count
//Description:   find edges that are shared by "count" number of facets
//Author: bwhanks
//Date: 2003
//===========================================================================
void FacetDataUtil::edges_by_count(DLIList<CubitFacet*> &facets,
                                   unsigned int count,
                                   DLIList<CubitFacetEdge*> &edges)
{
  int i;
  int j;

  // get a list of all edges from the facets
  DLIList<CubitFacetEdge*> edge_list;
  facets.reset();
  for (i=facets.size(); i>0; i--)
  {
    CubitFacet* p_facet = facets.get_and_step();
    for (j=0; j<3; j++)
    {
      edge_list.append(p_facet->edge(j));
    }
  }

  // reset edge marks
  edge_list.reset();
  for (i=edge_list.size(); i>0; i--)
    edge_list.get_and_step()->marked(0);

  // mark with the hit count
  edge_list.reset();
  for (i=edge_list.size(); i>0; i--)
  {
    CubitFacetEdge* p_edge = edge_list.get_and_step();
    p_edge->marked(p_edge->marked() + 1);
  }

  // create the output list of edges hit the number of times passed in
  edge_list.reset();
  for (i=edge_list.size(); i>0; i--)
  {
    CubitFacetEdge* p_edge = edge_list.get_and_step();
    if (static_cast<unsigned>(p_edge->marked()) == count)
      edges.append(p_edge);
  }

  // reset edge marks
  edge_list.reset();
  for (i=edge_list.size(); i>0; i--)
    edge_list.get_and_step()->marked(0);
}

//===========================================================================
//Function Name: ordered_point_edge_bdry
//Description:   get an ordered list of points and edges around the boundary
//               of a set of facets. The output "chain" consists of point -
//               edge - point - edge - ... around the boundary of the input
//               facets.  Having the ordered points as well as the edges
//               helps provides sense information for the bounding edges.
//Author: bwhanks
//Date: 2003
//===========================================================================
void FacetDataUtil::ordered_point_edge_bdry(DLIList<CubitFacet*> &facets,
                                            DLIList<FacetEntity*> &point_edge_chain)
{
  assert(point_edge_chain.size() == 0);
  DLIList<CubitFacetEdge*> unordered_edges;

  int use_count = 1; // get boundary edges (used only once by facets in the facet list)
  FacetDataUtil::edges_by_count(facets, use_count, unordered_edges);

//  Adds the start_edge to the list of boundary_edges and then attempts to find the
//  next edge by getting all of the edges belonging to the second point on this
//  edge and searching for an unmarked (-1) boundary edge.  Does this until it can't
//  find another edge.


  // mark all edges connected to points of boundary edges with -2
  // NOTE: using negative marks, since as boundary edges are added to the ordered list they
  // get marked with the order in which they are found

  // TODO - this seems very inefficient
  int i;
  int j;
  int k;
  DLIList<CubitFacetEdge*> pt_edge_list;
  CubitFacetEdge* cur_edge;
  unordered_edges.reset();
  for (i=unordered_edges.size(); i>0; i--)
  {
    cur_edge = unordered_edges.get_and_step();
    for (j=0; j<2; j++)
    {
      pt_edge_list.clean_out();
      cur_edge->point(j)->edges(pt_edge_list);
      pt_edge_list.reset();
      for (k=pt_edge_list.size(); k>0; k--)
        pt_edge_list.get_and_step()->marked(-2);
    }
  }

  // mark all boundary edges with -1
  unordered_edges.reset();
  for (i=unordered_edges.size(); i>0; i--)
    unordered_edges.get_and_step()->marked(-1);




  CubitPoint *edge_pt1, *edge_pt2, *pt2, *originalpt1;
  CubitFacetEdge *start_edge;
  CubitFacetEdge *this_edge;
  CubitBoolean keepgoing;

  int i_found_some = 0;
  unordered_edges.reset();
  start_edge = unordered_edges.get();

  // find the orientation of the first edge with respect to one of the facets
  // then order the edges accordingly
  facets.reset();
  int e_index;
  int start_edge_sense = 0;
  for (i=facets.size(); i>0; i--)
  {
    CubitFacet* p_facet = facets.get_and_step();
    if ( (e_index = p_facet->edge_index(start_edge)) >= 0 )
    {
      start_edge_sense = p_facet->edge_use(e_index);
      break;
    }
  }


  start_edge->marked(i_found_some);
  i_found_some += 1;
  if (1 == start_edge_sense)
  {
    originalpt1 = start_edge->point(0);
    point_edge_chain.append(originalpt1);
    point_edge_chain.append(start_edge);
    pt2 = start_edge->point(1);
  }
  else
  {
    assert(-1 == start_edge_sense);
    originalpt1 = start_edge->point(1);
    point_edge_chain.append(originalpt1);
    point_edge_chain.append(start_edge);
    pt2 = start_edge->point(0);
  }


//  Look for an edge having pt2 as a point and also being on the boundary.
  pt_edge_list.clean_out();
  pt2->edges(pt_edge_list);

  keepgoing = CUBIT_TRUE;
  pt_edge_list.reset();
  while ( keepgoing == CUBIT_TRUE ) {
    keepgoing = CUBIT_FALSE;
    for ( i = pt_edge_list.size(); i > 0; i-- ) {
      this_edge = pt_edge_list.get_and_step();
      if ( this_edge->marked() == -1 ) {
        i_found_some += 1;
        this_edge->marked(i_found_some);
        point_edge_chain.append(pt2);
        point_edge_chain.append(this_edge);

        edge_pt1 = this_edge->point(0);
        edge_pt2 = this_edge->point(1);
        if (pt2 == edge_pt1)
        {
          pt2 = edge_pt2;
        }
        else
        {
          assert(pt2 == edge_pt2);
          pt2 = edge_pt1;
        }


	      keepgoing = CUBIT_TRUE;
	      pt_edge_list.clean_out();
        pt2->edges(pt_edge_list);
	      break;
      }
    }
  }
  
  // clear marks
  for (i=unordered_edges.size(); i>0; i--)
  {
    cur_edge = unordered_edges.get_and_step();
    for (j=0; j<2; j++)
    {
      pt_edge_list.clean_out();
      cur_edge->point(j)->edges(pt_edge_list);
      pt_edge_list.reset();
      for (k=pt_edge_list.size(); k>0; k--)
        pt_edge_list.get_and_step()->marked(0);
    }
  }

  assert(pt2 == originalpt1);
}

//===========================================================================
//Function Name: partial_chain
//Description:   given an ordered point - edge boundary and a start and end
//               point on the boundary, return the point - edge chain between
//               the two points, inclusive.
//Author: bwhanks
//Date: 2003
//===========================================================================
CubitStatus FacetDataUtil::partial_chain(DLIList<FacetEntity*> &point_edge_chain,
                                         FacetEntity* point1,
                                         FacetEntity* point2,
                                         DLIList<FacetEntity*> &chain_between)
{
  // make sure the edges_between list is empty
  assert(chain_between.size() == 0);
  assert(point1 != point2);

  // find the start point in the list
  int pt1_index = point_edge_chain.where_is_item(point1);

  if (-1 == pt1_index)
    return CUBIT_FAILURE;

  point_edge_chain.reset();
  point_edge_chain.step(pt1_index);

  CubitBoolean b_done = CUBIT_FALSE;
  while (!b_done)
  {
    if (point_edge_chain.get() == point2)
      b_done = CUBIT_TRUE;

    chain_between.append(point_edge_chain.get_and_step());
  }

  return CUBIT_SUCCESS;
}

//===========================================================================
//Function Name: get_facet_points
//Description:   return a unique list of the points from a given set of facets
//Author: bwhanks
//Date: 2003
//===========================================================================
void FacetDataUtil::get_facet_points(DLIList<CubitFacet*> &cubit_facets,
                                     DLIList<CubitPoint*> &facet_points)
{
  int i;
  DLIList<CubitPoint*> cubit_points(cubit_facets.size() * 3);
  for( i = cubit_facets.size(); i--; )
  {
    CubitFacet* facet = cubit_facets.step_and_get();
    for( int j = 0; j < 3; j++ )
    {
      CubitPoint* pt = dynamic_cast<CubitPoint*>(facet->point(j));
      assert(!!pt);
      pt->marked(0);
      cubit_points.append(pt);
    }
  }
  for( i = cubit_points.size(); i--; )
  {
    CubitPoint* pt = cubit_points.step_and_get();
    pt->marked( pt->marked() + 1);
  }
  for( i = cubit_points.size(); i--; )
  {
    CubitPoint* pt = cubit_points.step_and_get();
    pt->marked( pt->marked() - 1 );
    if( pt->marked() > 0 )
      cubit_points.change_to(0);
  }
  cubit_points.remove_all_with_value(0);

  facet_points = cubit_points;
}

//===========================================================================
//Function Name: get_boundary_points
//Description:   return the boundary points from a list of facets (unordered)
//               assumes edges exist on the facets
//Author: bwhanks
//Date: 2003
//===========================================================================
void FacetDataUtil::get_boundary_points(DLIList<CubitFacet*> &facet_list,
                                        DLIList<CubitPoint*> &point_list)
{
  int ii, jj;
  CubitPoint *pt;
  CubitFacetEdge *edge;
  CubitFacet *facet;
  DLIList<CubitFacetEdge *>edge_list;
  for(ii=0; ii<facet_list.size(); ii++)
  {
    facet = facet_list.get_and_step();
    for(jj=0; jj<3; jj++)
    {
      edge = facet->edge( jj );
      if (1 == edge->num_adj_facets())
      {
        edge_list.append(edge);
        edge->point( 0 )->marked( 0 );
        edge->point( 1 )->marked( 0 );
      }
    }
  }
  for(ii=0; ii<edge_list.size(); ii++)
  {
    edge = edge_list.get_and_step();
    pt = edge->point( 0 );
    if (pt->marked() == 0)
    {
      pt->marked( 1 );
      point_list.append( pt );
    }
    pt = edge->point( 1 );
    if (pt->marked() == 0)
    {
      pt->marked( 1 );
      point_list.append( pt );
    }
  }
  for(ii=0; ii<point_list.size(); ii++)
    point_list.get_and_step()->marked(0);

}

//===========================================================================
//Function Name: get_boundary_edges
//Description: return an unordered list of edges at the boundary of a
//             set of facets
//Author: sjowen
//Date: 1/19/2004
//===========================================================================
void FacetDataUtil::get_boundary_edges(DLIList<CubitFacet*> &facet_list,
                                       DLIList<CubitFacetEdge*> &edge_list)
{
  int ii, jj;
  CubitFacetEdge *edge;
  CubitFacet *facet;
  for(ii=0; ii<facet_list.size(); ii++)
  {
    facet = facet_list.get_and_step();
    for(jj=0; jj<3; jj++)
    {
      edge = facet->edge( jj );
      if (1 == edge->num_adj_facets())
      {
        edge_list.append(edge);
      }
    }
  }
}

//===========================================================================
//Function Name: get_points
//Description:   get a unique set of points from facets
//Author: sjowen
//Date: 9/11/03
//===========================================================================
void FacetDataUtil::get_points(DLIList<CubitFacet*> &facet_list,
                               DLIList<CubitPoint*> &point_list)
{
  CubitFacet *facet;
  CubitPoint *pt;
  int ii, jj;
  for(ii=0; ii<facet_list.size(); ii++)
  {
    facet = facet_list.get_and_step();
    facet->point( 0 )->marked( 0 );
    facet->point( 1 )->marked( 0 );
    facet->point( 2 )->marked( 0 );
  }

  for (ii=0; ii<facet_list.size(); ii++)
  {
    facet = facet_list.get_and_step();
    for(jj=0; jj<3; jj++)
    {
      pt = facet->point(jj);
      if (pt->marked() == 0)
      {
        pt->marked(1);
        point_list.append(pt);
      }
    }
  }

  for(ii=0; ii<facet_list.size(); ii++)
  {
    facet = facet_list.get_and_step();
    facet->point( 0 )->marked( 0 );
    facet->point( 1 )->marked( 0 );
    facet->point( 2 )->marked( 0 );
  }
}

//===========================================================================
//Function Name: copy_facets
//Description:  make a complete copy of facets, points and edges
//Note:  copies edges and points with their "feature" flag
//Author: sjowen
//Date: 9/11/03
//===========================================================================
void FacetDataUtil::copy_facets(DLIList<CubitFacet*> &old_facet_list,
                                DLIList<CubitFacet*> &new_facet_list,
                                DLIList<CubitPoint*> &new_point_list,
                                DLIList<CubitFacetEdge*> &new_edge_list)
{
  // get a unique set of points from the facets

  DLIList<CubitPoint *> old_point_list;
  get_points(old_facet_list, old_point_list);
  CubitPoint **point_array = new CubitPoint* [old_point_list.size()];

  //- copy the points

  int ii;
  old_point_list.reset();
  CubitPoint *new_point, *the_point;
  for(ii=0; ii<old_point_list.size(); ii++)
  {
    the_point = old_point_list.get_and_step();
    new_point = new CubitPointData( the_point->coordinates() );
    the_point->marked( ii );
    new_point_list.append( new_point );
    point_array[ii] = new_point;
    if (the_point->is_feature())
      new_point->set_as_feature();
  }

  //- copy the facets

  int jj, idx;
  CubitFacet *new_facet, *the_facet;
  CubitPoint *points[3];

  old_facet_list.reset();
  for (ii=0; ii<old_facet_list.size(); ii++)
  {
    the_facet = old_facet_list.get_and_step();
    for (jj=0; jj<3; jj++)
    {
      idx = the_facet->point(jj)->marked();
      points[jj] = point_array[idx];
    }
    new_facet = new CubitFacetData( points[0], points[1], points[2] );
    new_facet_list.append( new_facet );
  }

  //- copy the edges

  int idx0, idx1;
  CubitFacetEdge *new_edge;
  CubitFacetEdge *old_edge;
  DLIList<CubitFacetEdge *>old_edge_list;
  get_edges(old_facet_list, old_edge_list);
  for(ii=0; ii<old_edge_list.size(); ii++)
  {
    old_edge = old_edge_list.get_and_step();
    idx0 = old_edge->point(0)->marked();
    idx1 = old_edge->point(1)->marked();
    new_edge = new CubitFacetEdgeData( point_array[idx0], point_array[idx1] );
    if (old_edge->is_feature())
      new_edge->set_as_feature();
    new_edge_list.append( new_edge );
  }

  delete [] point_array;

}

//===========================================================================
//Function Name: get_edges
//Description:  Populates the edge list from the list of facets - creates the
//              edges if necessary
//Author: sjowen
//Date: 9/11/03
//===========================================================================
void FacetDataUtil::get_edges(
  DLIList<CubitFacet *> &facet_list,
  DLIList<CubitFacetEdge *> &edge_list )
{
  int i, j;
  CubitPoint *p0, *p1;
  CubitFacet *facet_ptr;
  CubitFacetEdge *edge_ptr;
  DLIList<CubitFacet *> adj_facet_list;

  // mark the edges and create any that are missing
  facet_list.reset();  //have to reset this list to that the CubitFacetEdgeData's
                         //get constructed correctly
  for ( i = 0; i < facet_list.size(); i++)
  {
    facet_ptr = facet_list.get_and_step();
    for (j=0; j<3; j++) {
      edge_ptr = facet_ptr->edge(j);
      if (!(edge_ptr))
      {
        facet_ptr->get_edge_pts(j, p0, p1);
        edge_ptr = (CubitFacetEdge *) new CubitFacetEdgeData( p0, p1 );
      }
      edge_ptr->set_flag( 0 );
    }
  }

  // create a unique list of edges

  for ( i = 0; i < facet_list.size(); i++)
  {
    facet_ptr = facet_list.get_and_step();
    for (j=0; j<3; j++)
    {
      edge_ptr = facet_ptr->edge(j);
      if (0 == edge_ptr->get_flag())
      {
        edge_ptr->set_flag( 1 );
        edge_list.append( edge_ptr );
      }
    }
  }

  // reset the flags on the edges

  for ( i = 0; i < facet_list.size(); i++)
  {
    facet_ptr = facet_list.get_and_step();
    for (j=0; j<3; j++)
    {
      edge_ptr = facet_ptr->edge(j);
      edge_ptr->set_flag( 0 );
    }
  }
}

//============================================================================
// Function: quality
// Author: sjowen
// Description: this is the S.H. Lo metric for triangles that also takes into
//               account a surface normal.
// Date: 2/2003
//============================================================================
double FacetDataUtil::quality(CubitVector &c1, CubitVector &c2, CubitVector &c3,
                              CubitVector &surf_normal)
{
#define TWO_ROOT_THREE  3.46410161514
   double area2, alpha;
   double length1, length2, length3;
   CubitVector edge1, edge2, edge3;

   // create edges from the vertices

   edge1 = c3 - c2;
   edge2 = c3 - c1;
   edge3 = c2 - c1;

   // compute twice the area

   CubitVector normal = edge3 * edge2;
   area2 = normal.length();

   length1 = edge1.length_squared();
   length2 = edge2.length_squared();
   length3 = edge3.length_squared();

   alpha = TWO_ROOT_THREE * area2 / (length1 + length2 + length3);

   // modify the alpha metric by the dot product of the triangle normal
   // with the surface normal.

   if (fabs(area2) < CUBIT_RESABS)
     alpha = 0.0;
   else
   {
     normal /= area2;
     double dot = normal % surf_normal;
     double penalty = pow(dot, 5);
     alpha *= penalty;
   }

   return alpha;
}

//================================================================================
// Description: compares length of two edges.
// Notes:       used for DLIList sort
// Author     : Steve Owen
// Date       : 9/27/2003
//================================================================================
int FacetDataUtil::edge_compare(CubitFacetEdge *&ea, CubitFacetEdge *&eb)
{
  double la = ea->length();
  double lb = eb->length();
  if (la < lb)
    return -1;
  else if (lb < la)
    return 1;
  return 0;
}

//================================================================================
// Description: collapse edges that are smaller than 10% the average length
// Notes:       non_manifold_only=true indicates that only edges that are attached
//              to more than two edges will be candidates for collapse.
// Author     : Steve Owen
// Date       : 9/27/2003
//================================================================================
CubitStatus FacetDataUtil::collapse_short_edges(
  DLIList<CubitFacet*> &facet_list,
  CubitBoolean non_manifold_only)
{
#define COLLAPSE_TOLERANCE 0.3

  CubitStatus rv = CUBIT_SUCCESS;
  int ncollapse = 0;

  // get the edges

  DLIList <CubitFacetEdge *>edge_list;
  DLIList <CubitPoint *> point_list;
  DLIList <CubitFacet *> one_facet;
  DLIList <CubitFacetEdge *> facet_edges;
  FacetDataUtil::get_edges( facet_list, edge_list );
  FacetDataUtil::get_points( facet_list, point_list );

  // determine average length and get the collapse threshold

  int ii, jj, kk;
  double len;
  double tot_len = 0.0;
  for(ii=0; ii<point_list.size(); ii++)
    point_list.get_and_step()->marked(0);
  CubitFacetEdge *edge;
  for(ii=0; ii<edge_list.size(); ii++)
  {
    edge = edge_list.get_and_step();
    len = edge->length();
    tot_len += len;
    edge->marked(0);
  }
  len = tot_len / edge_list.size();
  double collapse_tol = COLLAPSE_TOLERANCE * len;
  //check(edge_list, facet_list);

  // get a list of the edges that qualify for collapse

  DLIList<CubitFacetEdge *>short_edges;
  for(ii=0; ii<edge_list.size(); ii++)
  {
    edge = edge_list.get_and_step();
    if (non_manifold_only && edge->num_adj_facets() == 2)
      continue;
    len = edge->length();
    if (len < collapse_tol)
    {
      short_edges.append(edge);
    }
  }
  //sort them

  short_edges.sort(FacetDataUtil::edge_compare);

  // main loop

  int nedges = short_edges.size();
  for(ii=0; ii<nedges; ii++)
  {
    edge = short_edges.get_and_step();
    if (edge->marked() == 1)
      continue;
    len = edge->length();

    bool collapse = true;
    CubitPoint *pt[2];
    CubitPoint *pt1, *pt2;
    pt[0] = edge->point(0);
    pt[1] = edge->point(1);

    // compute a new candidate location for the merged points

    CubitVector new_location;
    CubitVector cpt[2];
    cpt[0] = pt[0]->coordinates();
    cpt[1] = pt[1]->coordinates();
    new_location = (cpt[0] + cpt[1]) * 0.5;

    for (jj=0; jj<2 && collapse; jj++)
    {
      DLIList<CubitFacet *> adjfacets;
      pt[jj]->facets( adjfacets );

      // check all facets adjacent this point that don't contain the edge
      // to make sure the resulting facet will be valid (we aren't inverting anything)

      for(kk=0; kk<adjfacets.size() && collapse; kk++)
      {
        CubitFacet *facet = adjfacets.get_and_step();
        pt1 = facet->next_node( pt[jj] );
        pt2 = facet->next_node( pt1 );
        int eidx = facet->edge_index( edge );
        if (eidx < 0)
          // don't check facets that have the current edge - they'll be deleted anyway
        {

          CubitVector cpt1 = pt1->coordinates();
          CubitVector cpt2 = pt2->coordinates();
          CubitVector norm = facet->normal();
          double q0 = FacetDataUtil::quality(cpt[jj], cpt1, cpt2, norm);
          double q1 = FacetDataUtil::quality(new_location, cpt1, cpt2, norm);
          if (!(q1 > 0.0 || q1 > q0))
          {
            collapse = false;
          }
        }
      }
    }

    if (collapse)
    {
      DLIList<CubitFacet *> new_facets;
      CubitPoint *collpt = pt[0];
      CubitPoint *delpt = pt[1];
      DLIList<CubitFacetEdge *> adjedges;
      CubitFacetEdge *adjedge;
      delpt->edges( adjedges );
      CubitFacet *facet;


      int mydebug = 0;
      if (mydebug)
      {
        dcolor(4);
        draw_edge( edge );
        dview();
      }

      collpt->set(new_location);

      // delete all facets adjacent to the delpt.

      DLIList<CubitFacet *> adjfacets;
      delpt->facets( adjfacets );

      for(jj=0; jj<adjfacets.size(); jj++)
      {
        facet = adjfacets.get_and_step();
        pt1 = facet->next_node(delpt);
        pt2 = facet->next_node(pt1);
        int eidx = facet->edge_index( edge );
        facet_list.move_to(facet);
        delete facet;
        facet = NULL;
        if (eidx >= 0)
        {
          //if this facet is adjecnt the edge, then just remove it from the list
          facet_list.extract();
        }
        else
        {

          // get or create edges as needed

          CubitFacetEdge *e[3];
          e[0] = collpt->get_edge( pt1 );
          e[1] = pt1->get_edge( pt2 );
          e[2] = pt2->get_edge( collpt );
          for (kk=0; kk<3; kk++)
          {
            if (!(e[kk]))
            {
              switch (kk)
              {
              case 0: e[kk] = (CubitFacetEdge *) new CubitFacetEdgeData( collpt, pt1 ); break;
              case 1: e[kk] = (CubitFacetEdge *) new CubitFacetEdgeData( pt1, pt2 ); break;
              case 2: e[kk] = (CubitFacetEdge *) new CubitFacetEdgeData( pt2, collpt ); break;
              }
              edge_list.append( e[kk] );
            }
          }

          // create a new facet with the points from the old facet and the collpt

          facet = new CubitFacetData( e[0], e[1], e[2] );
          new_facets.append(facet);

          // create edges on the facet

          facet_list.change_to( facet );
        }
      }

      for(jj=0; jj<adjedges.size(); jj++)
      {
        adjedge = adjedges.get_and_step();
        adjedge->marked(1);  // mark it for deletion later
        assert(adjedge->num_adj_facets() == 0);
      }
      assert(delpt->num_adj_facets() == 0);
      ncollapse++;

      //check(edge_list, facet_list);
    }

  }

  PRINT_INFO("Collapsed %d short edges in triangulation\n", ncollapse);

  // delete points and edges that aren't used
  if (ncollapse > 0)
  {
    for(ii=0; ii<edge_list.size(); ii++)
    {
      edge = edge_list.get_and_step();
      if (edge->marked())
      {
        assert(edge->num_adj_facets() == 0);
        delete edge;
      }
      else if(edge->num_adj_facets() < 2){
        PRINT_ERROR("Unexpected result while collapsing an edge.\n");
        return CUBIT_FAILURE;
          //assert(edge->num_adj_facets() >= 2);
      }
    }
    CubitPoint *point;
    for(ii=0; ii<point_list.size(); ii++)
    {
      point = point_list.get_and_step();
      if (point->num_adj_facets() == 0)
      {
        delete point;
      }
    }
  }

  return rv;
}

//===========================================================================
//  Function: check
//  Purpose:  debugging only
//  Date:     10/2003
//  Author:   sjowen
//===========================================================================
void FacetDataUtil::check(DLIList<CubitFacetEdge *> &edge_list,
                          DLIList<CubitFacet *> &facet_list)
{

  CubitBox box;
  CubitFacetEdge *edge;
  int ii;
  int nedges = 0;
  for(ii=0; ii<edge_list.size(); ii++)
  {
    edge = edge_list.get_and_step();
    if (edge->marked() == 1)
      continue;
    int nadj = edge->num_adj_facets();
    if (nadj <= 1)
    {
      CubitVector v0 = edge->point(0)->coordinates();
      CubitVector v1 = edge->point(1)->coordinates();
      CubitVector min(CUBIT_MIN(v0.x(), v1.x()), CUBIT_MIN(v0.y(), v1.y()), CUBIT_MIN(v0.z(), v1.z()));
      CubitVector max(CUBIT_MAX(v0.x(), v1.x()), CUBIT_MAX(v0.y(), v1.y()), CUBIT_MAX(v0.z(), v1.z()));
      CubitBox ebox(min, max);
      if (nedges == 0)
        box.reset(min, max);
      else
        box |= ebox;
//      dcolor(3);
//      dfldraw(facet_list);
      dcolor(4);
      dedraw(edge);
      dpoint(v0);
      dpoint(v1);
      nedges++;
    }
  }
  if (nedges > 0)
  {
    dzoom(box);
    dview();
  }
}

//===========================================================================
//  Function: draw_edge
//  Purpose:
//  Date:     10/2003
//  Author:   sjowen
//===========================================================================
void FacetDataUtil::draw_edge( CubitFacetEdge *edge )
{
  DLIList<CubitFacet *>adjfacets;
  edge->facets( adjfacets );
  int ii,jj;
  CubitFacet *facet;
  CubitBox box;
  for(jj=0; jj<2; jj++)
  {
    CubitPoint *p = edge->point(jj);
    adjfacets.clean_out();
    p->facets( adjfacets );
    for(ii=0; ii<adjfacets.size(); ii++)
    {
      facet = adjfacets.get_and_step();
      if (jj==0 && ii==0)
        box = facet->bounding_box();
      else
        box |= facet->bounding_box();
      dfdraw(facet);
    }
    dpoint(p->coordinates());
    dpoint(p->coordinates());
  }
  dzoom(box);
}

//=============================================================================
//Function:  is_point_in_polyhedron
//Description: point-in-polyhedron test; polyhedron must be water-tight,
//             manifold, triangle-tiled.  Casts a random ray in positive
//             x, y, and z.  Counts crossings.  Starts over with a recast
//             if a triangle is perpendicular to the ray or if the ray
//             hits a triangle edge or vertex.
//             Also tests for point on polyhedron.
//Author: John Fowler
//Date: 9/30/03
//=============================================================================
CubitStatus FacetDataUtil::is_point_in_polyhedron(
                                     DLIList<CubitFacet *> &tfacet_list,
                                     CubitVector &point_coords,
                                     CubitPointContainment &is_point_in)
{
unsigned int i;
CubitBox bbox;
CubitFacet *facet;
CubitVector bbox_min, bbox_max, ray;
CubitStatus status;
CubitPointContainment pt_status;
bool is_outside, recast, is_on_plane;
double xpt, ypt, zpt;
double rayx, rayy, rayz;
int number_of_recasts;

  point_coords.get_xyz(xpt, ypt, zpt);
  recast = true;
  number_of_recasts = 0;
  while ( (recast == true) && (number_of_recasts < 10) ) {
    recast = false;
    is_outside = true;
    random_positive_ray(ray); // a random positively-pointing ray
    ray.get_xyz(rayx,rayy,rayz);
    for ( i = tfacet_list.size(); i > 0; i-- ) {
      facet = tfacet_list.get_and_step();
      bbox = facet->bounding_box();
      bbox_min = bbox.minimum();
      bbox_max = bbox.maximum();
      //  Because the ray is positive in direction, discard bounding boxes
      //  that are entirely < the starting point.
      if ( (xpt-bbox_max.x()) > CUBIT_RESABS ) continue;
      if ( (ypt-bbox_max.y()) > CUBIT_RESABS ) continue;
      if ( (zpt-bbox_max.z()) > CUBIT_RESABS ) continue;
      if ( ray_intersects_boundingbox(point_coords,ray,bbox) == false )
        continue;
      CubitPlane plane = facet->plane();
      CubitVector normal = plane.normal();
      double xnorm, ynorm, znorm;
      xnorm = normal.x();
      ynorm = normal.y();
      znorm = normal.z();
      //  Intersect the ray with the facet plane.  If ray is perpendicular to
      //  facet plane and point is on it, recast the ray and try again.
      double denominator = rayx*xnorm + rayy*ynorm + rayz*znorm;
      double distanc = xnorm*xpt + ynorm*ypt + znorm*zpt + plane.coefficient();
      if ( fabs(denominator) < GEOMETRY_RESABS )
      {
        if ( fabs(distanc) < GEOMETRY_RESABS ) {
          recast = true;
	  break;
	} else continue;  // point is not on plane and ray is parallel to plane
      }
      double t, xintpt, yintpt, zintpt;
      t = -distanc;
      t /= denominator;
      if ( fabs(t) < GEOMETRY_RESABS ) is_on_plane = true;
      else if ( t < GEOMETRY_RESABS ) continue;
      else is_on_plane = false;
      xintpt = xpt + t*rayx;
      yintpt = ypt + t*rayy;
      zintpt = zpt + t*rayz;
      CubitPoint *pt;
      //  We need to project the triangle onto 2D.  Use the smaller components of
      //  the normal to detemine a good projection.
      double x0, y0, x1, y1, x2, y2;
      if ( (fabs(xnorm) >= fabs(ynorm)) && (fabs(xnorm) >= fabs(znorm)) ){  //  Use y,z
        pt = facet->point(0);
        x0 = pt->y(); y0 = pt->z();
        pt = facet->point(1);
        x1 = pt->y(); y1 = pt->z();
        pt = facet->point(2);
        x2 = pt->y(); y2 = pt->z();
        status = pt_in_tri_2d(yintpt,zintpt,x0,y0,x1,y1,x2,y2,pt_status);
      } else if (fabs(ynorm) >= fabs(znorm)) {  //  Use z,x
        pt = facet->point(0);
        x0 = pt->x(); y0 = pt->z();
        pt = facet->point(1);
        x1 = pt->x(); y1 = pt->z();
        pt = facet->point(2);
        x2 = pt->x(); y2 = pt->z();
        status = pt_in_tri_2d(xintpt,zintpt,x0,y0,x1,y1,x2,y2,pt_status);
      } else {  // Use x,y
        pt = facet->point(0);
        x0 = pt->x(); y0 = pt->y();
        pt = facet->point(1);
        x1 = pt->x(); y1 = pt->y();
        pt = facet->point(2);
        x2 = pt->x(); y2 = pt->y();
        status = pt_in_tri_2d(xintpt,yintpt,x0,y0,x1,y1,x2,y2,pt_status);
      }

      if ( status == CUBIT_FAILURE ) {
        recast = true;
	break;
      }
      if ( pt_status == CUBIT_PNT_OUTSIDE ) continue;
      //  Is the point on the triangle?
      if ( is_on_plane == true ) {
        is_point_in = CUBIT_PNT_BOUNDARY;
        return CUBIT_SUCCESS;
      }
      if ( pt_status == CUBIT_PNT_INSIDE ) is_outside = ! is_outside;
      else if ( pt_status == CUBIT_PNT_BOUNDARY ) {
        recast = true;
        break;
      }
    }
    if ( recast == true ) {
      tfacet_list.reset();
      number_of_recasts += 1;
    }
  }
  if ( recast == true ) {
    PRINT_ERROR("Number of recasts in point-in-polygon exceeded 10.\n");
    return CUBIT_FAILURE;
  }
  if ( is_outside == false ) is_point_in = CUBIT_PNT_INSIDE;
  else is_point_in = CUBIT_PNT_OUTSIDE;

  return CUBIT_SUCCESS;
}

//===========================================================================
//  Function: pt_in_tri_2d
//  Purpose:
//  Date:     10/2003
//  Author:   John Fowler
//===========================================================================
CubitStatus FacetDataUtil::pt_in_tri_2d(double xpt, double ypt,
                       double x0, double y0,
		       double x1, double y1,
		       double x2, double y2,
		       CubitPointContainment &is_point_in)
{
//  From Schneider & Eberly, "Geometric Tools for COmputer Graphics",
//  Chap. 13.3.1.  If triangle is needle-thin, CUBIT_FAILURE might be
//  returned, in wich case is_point_in is undefined.

double c0, c1, c2;
double e0x, e1x, e2x, e0y, e1y, e2y;
double n0x, n1x, n2x, n0y, n1y, n2y;
double denom0, denom1, denom2;

  e0x = x1 - x0; e0y = y1 - y0;
  e1x = x2 - x1; e1y = y2 - y1;
  e2x = x0 - x2; e2y = y0 - y2;
  n0x = e0y; n0y = -e0x;
  n1x = e1y; n1y = -e1x;
  n2x = e2y; n2y = -e2x;
  denom0 = n1x*e0x + n1y*e0y;
  if ( fabs(denom0) < CUBIT_RESABS ) {
    PRINT_ERROR("Failure in pt_in_tri_2d; needle-thin triangle encountered.\n");
    return CUBIT_FAILURE;
  }
  denom1 = n2x*e1x + n2y*e1y;
  if ( fabs(denom1) < CUBIT_RESABS ) {
    PRINT_ERROR("Failure in pt_in_tri_2d; needle-thin triangle encountered.\n");
    return CUBIT_FAILURE;
  }
  denom2 = n0x*e2x + n0y*e2y;
  if ( fabs(denom2) < CUBIT_RESABS ) {
    PRINT_ERROR("Failure in pt_in_tri_2d; needle-thin triangle encountered.\n");
    return CUBIT_FAILURE;
  }

  c0 = -( n1x*(xpt-x1) + n1y*(ypt-y1) )/denom0;
  c1 = -( n2x*(xpt-x2) + n2y*(ypt-y2) )/denom1;
  c2 = -( n0x*(xpt-x0) + n0y*(ypt-y0) )/denom2;

  if ( (c0 > 0.0) && (c1 > 0.0) && (c2 > 0.0) ) is_point_in = CUBIT_PNT_INSIDE;
  else if ( (c0 < 0.0) || (c1 < 0.0) || (c2 < 0.0) ) is_point_in = CUBIT_PNT_OUTSIDE;
  else is_point_in = CUBIT_PNT_BOUNDARY;

  return CUBIT_SUCCESS;

}

//===========================================================================
//  Function: random_positive_ray
//  Purpose:
//  Date:     10/2003
//  Author:   John Fowler
//===========================================================================
void FacetDataUtil::random_positive_ray(CubitVector& ray)
{
  double temp;
  double rayx = 0.0, rayy = 0.0, rayz = 0.0;

  temp = 0.0;
  while ( temp < 1.e-6 ) {
    rayx = (double(rand())/(RAND_MAX+1.0));
    rayy = (double(rand())/(RAND_MAX+1.0));
    rayz = (double(rand())/(RAND_MAX+1.0));
    temp = sqrt(rayx*rayx + rayy*rayy + rayz*rayz);
  }
  rayx /= temp;
  rayy /= temp;
  rayz /= temp;

  ray.set(rayx,rayy,rayz);

}

//===========================================================================
//  Function: ray_intersects_boundingbox
//  Purpose:
//  Date:     10/2003
//  Author:   John Fowler
//===========================================================================
bool FacetDataUtil::ray_intersects_boundingbox(CubitVector& point, CubitVector& ray, const CubitBox& bbox)
{
double t, xtest, ytest, ztest;
double xdir, ydir, zdir, xpt, ypt, zpt, xmin, ymin, zmin, xmax, ymax, zmax;
CubitVector bbox_min, bbox_max;

  point.get_xyz(xpt,ypt,zpt);
  ray.get_xyz(xdir,ydir,zdir);
  bbox_min = bbox.minimum();
  bbox_max = bbox.maximum();
  bbox_min.get_xyz(xmin,ymin,zmin);
  bbox_max.get_xyz(xmax,ymax,zmax);
  xmin -= GEOMETRY_RESABS;  //  So we don't miss any.
  ymin -= GEOMETRY_RESABS;
  zmin -= GEOMETRY_RESABS;
  xmax += GEOMETRY_RESABS;
  ymax += GEOMETRY_RESABS;
  zmax += GEOMETRY_RESABS;

 //  Notice that we are only interested in bboxes on the non-negative side of t.
  if ( fabs(xdir) > CUBIT_RESABS ) {
    // test xmin plane
    t = (xmin - xpt)/xdir;
    if ( t >= 0.0 ) {
      ytest = ypt + t*ydir;
      ztest = zpt + t*zdir;
      if ( (ytest >= ymin) && (ytest <= ymax) && (ztest >= zmin) && (ztest <= zmax) )
        return true;
    }
    // test xmax plane
    t = (xmax - xpt)/xdir;
    if ( t > 0.0 ) {
      ytest = ypt + t*ydir;
      ztest = zpt + t*zdir;
      if ( (ytest >= ymin) && (ytest <= ymax) && (ztest >= zmin) && (ztest <= zmax) )
        return true;
    }
  }
  if ( fabs(ydir) > CUBIT_RESABS ) {
    // test ymin plane
    t = (ymin - ypt)/ydir;
    if ( t >= 0.0 ) {
      xtest = xpt + t*xdir;
      ztest = zpt + t*zdir;
      if ( (xtest >= xmin) && (xtest <= xmax) && (ztest >= zmin) && (ztest <= zmax) )
        return true;
    }
    // test ymax plane

    t = (ymax - ypt)/ydir;
    if ( t > 0.0 ) {
      xtest = xpt + t*xdir;
      ztest = zpt + t*zdir;
      if ( (xtest >= xmin) && (xtest <= xmax) && (ztest >= zmin) && (ztest <= zmax) )
        return true;
    }
  }
  if ( fabs(zdir) > CUBIT_RESABS ) {
    // test zmin plane
    t = (zmin - zpt)/zdir;
    if ( t > 0.0 ) {
      xtest = xpt + t*xdir;
      ytest = ypt + t*ydir;
      if ( (xtest >= xmin) && (xtest <= xmax) && (ytest >= ymin) && (ytest <= ymax) )
        return true;
    }
    // test zmax plane
    t = (zmax - zpt)/zdir;
    if ( t > 0.0 ) {
      xtest = xpt + t*xdir;
      ytest = ypt + t*ydir;
      if ( (xtest >= xmin) && (xtest <= xmax) && (ytest >= ymin) && (ytest <= ymax) )
        return true;
    }
  }

  return false;
}

//===========================================================================
//  Function: write_facets
//  Purpose:  write a list of facets to a cubit facet file
//  Date:     11/28/2002
//  Author:   sjowen
//===========================================================================
CubitStatus FacetDataUtil::write_facets( const char *file_name, DLIList<CubitFacet *> &facet_list)
{
  FILE *fp = fopen(file_name, "w");
  if (!fp)
  {
    PRINT_ERROR("Couldn't open file %s for writing a facet file.\n", file_name);
    return CUBIT_FAILURE;
  }

  DLIList<CubitPoint*> point_list;
  get_points(facet_list, point_list);

  fprintf(fp, "%d %d\n", point_list.size(), facet_list.size());
  CubitPoint *pt;
  int ii;
  for (ii=1; ii<=point_list.size(); ii++)
  {
    pt = point_list.get_and_step();
    pt->set_id(ii);
    fprintf(fp, "%d %f %f %f\n", ii, pt->x(), pt->y(), pt->z() );
  }

  CubitFacet *facet;
  for (ii=1; ii<=facet_list.size(); ii++)
  {
    facet = facet_list.get_and_step();
    fprintf(fp, "%d %d %d %d\n", ii, facet->point(0)->id(),
             facet->point(1)->id(), facet->point(2)->id());
  }

  fclose(fp);
  return CUBIT_SUCCESS;

}

//=============================================================================
//Function:  split_into_shells (PUBLIC)
//Description: split this set of facets into multiple surfaces where there are
//             may be disjoint regions.
//Author: sjowen
//Date: 8/27/2003
//=============================================================================
CubitStatus FacetDataUtil::split_into_shells(
  DLIList<CubitFacet *> &tfacet_list,
  DLIList<CubitQuadFacet *> &qfacet_list,
  DLIList<DLIList<CubitFacet *> *> &shell_list,
  CubitBoolean &is_water_tight)
{
  CubitStatus stat = CUBIT_SUCCESS;
    //need to init this variable otherwise the caller must have.
  is_water_tight = CUBIT_TRUE;

  // combine the quads and tri lists

  DLIList <CubitFacet *> facet_list = tfacet_list;
  int ii;
  CubitQuadFacet *qfacet_ptr;
  for (ii=0; ii<qfacet_list.size(); ii++)
  {
    qfacet_ptr = qfacet_list.get_and_step();
    facet_list.append(qfacet_ptr->get_tri_facet( 0 ));
    facet_list.append(qfacet_ptr->get_tri_facet( 1 ));
  }

  // mark all the facets to begin with

  CubitFacet *facet_ptr;
  for (ii=0; ii<facet_list.size(); ii++)
  {
    facet_ptr = facet_list.get_and_step();
    facet_ptr->marked( 1 );
  }

  // populate the facet edge lists

  DLIList<CubitFacetEdge *> edge_list;
  FacetDataUtil::get_edges( facet_list, edge_list );

  // some debug stuff to draw the facets

  int mydebug=0;
  if (mydebug)
  {
    for (ii=0; ii<facet_list.size(); ii++)
    {
      facet_ptr = facet_list.get_and_step();
      GfxDebug::draw_facet(facet_ptr,CUBIT_YELLOW);
    }
    GfxDebug::flush();
    GfxDebug::mouse_xforms();
  }

  // Go through the surfaceElemList and pull facets off one by one as we
  // determine which surface it belongs to.  Continue until we have depleted
  // the list

  int jj;
  int num_surfs_created = 0;
  int num_facets_remaining = facet_list.size();
  while( num_facets_remaining > 0)
  {
    // create a new shell to hold the face info

    num_surfs_created++;
    DLIList<CubitFacet *> *shell_ptr = new DLIList<CubitFacet *>;

    // start with the first facet and create a list of all elements
    // attached to the facet

    CubitBoolean shell_is_water_tight = CUBIT_TRUE;
    CubitFacet *start_facet_ptr = facet_list.get_and_step();
    stat = get_adj_facets_on_shell( start_facet_ptr, shell_ptr,
                                    shell_is_water_tight, mydebug );
    if (stat != CUBIT_SUCCESS)
      return stat;
    if (!shell_is_water_tight)
      is_water_tight = CUBIT_FALSE;

    shell_list.append( shell_ptr );

    // remove the facets in this shell from the facet list

    for( jj = facet_list.size(); jj > 0; jj-- )
    {
      facet_ptr = facet_list.get();
      if (facet_ptr->marked() == 0)
      {
        facet_list.change_to( NULL );
      }
      facet_list.step();
    }
    facet_list.remove_all_with_value( NULL );
    num_facets_remaining = facet_list.size();
  }

  return CUBIT_SUCCESS;
}

//=============================================================================
//Function:  get_adj_facets_on_shell (PRIVATE)
//Description: non recursive function that creates a list of all facets connected
//             the passed in face that are on the same surface
//Author: sjowen
//Date: 12/22/00
//=============================================================================
CubitStatus FacetDataUtil::get_adj_facets_on_shell(
  CubitFacet *start_facet_ptr,
  DLIList<CubitFacet *> *shell_ptr,
  CubitBoolean &is_water_tight,
  int mydebug)
{
  int found = 0;
  int ii, jj;
  CubitStatus stat = CUBIT_SUCCESS;
  DLIList<CubitFacet*> temp_list;
  CubitFacet *facet_ptr = NULL;
  CubitFacet *adj_facet_ptr = NULL;
  DLIList<CubitFacetEdge *>edge_list;
  CubitFacetEdge *edge_ptr = NULL;
  DLIList<CubitFacet *>adj_facet_list;

  // add this facet to the list

  temp_list.append( start_facet_ptr );
  start_facet_ptr->marked( 0 );

  while (temp_list.size())
  {
    facet_ptr = temp_list.pop();
    if (facet_ptr->marked() == 0)
    {
      shell_ptr->append( facet_ptr );
      if (mydebug)
      {
        GfxDebug::draw_facet(facet_ptr, CUBIT_RED);
        GfxDebug::flush();
      }
      edge_list.clean_out();
      facet_ptr->edges( edge_list );
      for (ii=0; ii<edge_list.size(); ii++)
      {
        edge_ptr = edge_list.get_and_step();
        adj_facet_list.clean_out();
        edge_ptr->facets( adj_facet_list );
        found = 0;

        for (jj=0; jj<adj_facet_list.size() && !found; jj++)
        {
          adj_facet_ptr = adj_facet_list.get_and_step();
          if (adj_facet_ptr != facet_ptr)
          {

            // go to its neighbor if it is part of the surface

            if (adj_facet_ptr->marked() == 1)
            {
              temp_list.append( adj_facet_ptr );
              adj_facet_ptr->marked( 0 );
              found = 1;
            }
          }

        }
        if (is_water_tight && adj_facet_list.size() == 1)
        {
          is_water_tight = CUBIT_FALSE;
        }
      }
    }
  }

  return stat;
}

//=============================================================================
//Function:  stitch_facets (PRIVATE)
//Description: attempt to merge facets to form a watertight model
//Author: sjowen
//Date: 9/9/03
//=============================================================================
CubitStatus FacetDataUtil::stitch_facets(
  DLIList<DLIList<CubitFacet *> *> &shell_list,
  double tol,
  CubitBoolean &is_water_tight,
  CubitBoolean write_result)
{

  CubitStatus rv = CUBIT_SUCCESS;
  int npmerge = 0;
  int nemerge = 0;
  DLIList<CubitPoint*> unmerged_points;
  is_water_tight = CUBIT_FALSE;

  rv = merge_coincident_vertices( shell_list, tol, npmerge, nemerge, unmerged_points );
  if (rv != CUBIT_SUCCESS)
    return rv;

  int nnomerge = unmerged_points.size();
  if (nnomerge == 0)
  {
    is_water_tight = CUBIT_TRUE;
  }
  else
  {
    rv = merge_colinear_vertices( shell_list, tol, unmerged_points,
                                  npmerge, nemerge, nnomerge );
    if (rv != CUBIT_SUCCESS)
      return rv;
    if (nnomerge == 0)
    {
      is_water_tight = CUBIT_TRUE;
      int mydebug = 0;
      if (mydebug)  // make sure its really water-tight
      {
        DLIList<CubitFacet *> *shell_ptr;
        DLIList<CubitFacetEdge *>boundary_edges;
        for (int ii=0; ii<shell_list.size(); ii++)
        {
          shell_ptr = shell_list.get_and_step();
          FacetDataUtil::get_boundary_edges(*shell_ptr, boundary_edges);
        }
        if (boundary_edges.size() > 0)
        {
          PRINT_ERROR("Not Water-tight!\n");
        }
      }
    }
  }

  if (write_result)
  {
    if (npmerge > 0 || nemerge > 0)
    {
      PRINT_INFO("%d facet vertices and %d facet edges were successfully merged.\n",
                 npmerge, nemerge);
    }
    if (is_water_tight)
    {
      PRINT_INFO("Facets are water-tight.\n");
    }
    else
    {
      PRINT_INFO("%d facet vertices on model boundary detected that could not be merged.\n", nnomerge);
    }
  }

  return rv;
}

//=============================================================================
//Function:  merge_colinear_vertices (PRIVATE)
//Description: check if any vertices fall on a free facet edge.  If so - split
//             the  adjacent facet to merge
//Author: sjowen
//Date: 1/19/2004
//=============================================================================
CubitStatus FacetDataUtil::merge_colinear_vertices(
  DLIList<DLIList<CubitFacet *> *> &shell_list,
  double tol,
  DLIList<CubitPoint *> &merge_point_list,  // points to attempt to merge
  int &npmerge,    // return number of vertices merged
  int &nemerge,    // return number of edges merged
  int &nnomerge)   // return number of vertices in list NOT merged
{
  nnomerge = 0;
  int mydebug = 0;

  int ii, jj, kk, ll, mm, pt_shell_id, edge_shell_id;
  RTree <CubitFacetEdge*> r_tree(tol);
  shell_list.reset();
  DLIList<CubitFacet *> *shell_ptr = NULL;
  CubitPoint *pt;
  CubitFacetEdge *edge;
  DLIList<CubitFacetEdge *>boundary_edges;
  for (ii=0; ii<shell_list.size(); ii++)
  {
    shell_ptr = shell_list.get_and_step();

    // get the boundary edges from the shell - these are candidates for merging.
    // note that this won't merge internal facets

    DLIList<CubitFacetEdge *>shell_boundary_edges;
    FacetDataUtil::get_boundary_edges(*shell_ptr, shell_boundary_edges);

    // mark each of the edges with a shell id so we know when to merge shells
    // and add them to the kdtree for fast spatial searching
    edge_shell_id = ii+1;
    for(jj=0; jj<shell_boundary_edges.size(); jj++)
    {
      edge = shell_boundary_edges.get_and_step();
      edge->point(0)->marked(edge_shell_id);
      edge->point(1)->marked(edge_shell_id);
      edge->marked(edge_shell_id);
      r_tree.add(edge);
      boundary_edges.append(edge);
    }
    for(jj=0; jj<shell_ptr->size(); jj++)
      shell_ptr->get_and_step()->marked( edge_shell_id );
  }
  if (mydebug)
  {
    dview();
  }

  // find points in merge_point_list that are colinear with edges in boundary_edges

  CubitBox ptbox;
  CubitVector ptmin, ptmax, coord;
  DLIList<CubitFacetEdge *>close_edges;
  CubitFacetEdge *close_edge;
  CubitFacet *facet;
  CubitPoint *p0, *p1;
  DLIList<CubitPoint*>adj_pt_list;
  DLIList<CubitPoint*>del_points;
  for(ii=0; ii<merge_point_list.size(); ii++)
  {
    pt = merge_point_list.get_and_step();
    if (pt->marked() < 0)  // has already been merged
      continue;

    // find the closest edges in the rtree

    coord = pt->coordinates();
    ptmin.set( coord.x() - tol, coord.y() - tol, coord.z() - tol );
    ptmax.set( coord.x() + tol, coord.y() + tol, coord.z() + tol );
    ptbox.reset( ptmin, ptmax );
    close_edges.clean_out();
    r_tree.find(ptbox, close_edges);

    // We did find something - go try to merge

    CubitBoolean was_merged = CUBIT_FALSE;
    for (jj=0; jj<close_edges.size() && !was_merged; jj++)
    {
      // test the next edge on the list

      close_edge = close_edges.get_and_step();

      // make sure the edge does not contain the merge point

      if (close_edge->point(0) == pt || close_edge->point(1) == pt)
        continue;

      // check to see if the edge is within tolerance of the merge point

      CubitVector loc_on_edge;
      CubitBoolean is_outside_edge;
      double dist = close_edge->dist_to_edge(pt->coordinates(),
        loc_on_edge, is_outside_edge);

      // allow for some surface curvature. permit the point to be close to
      // the edge but not on.  May want to modify this factor if it not closing
      // the facets correctly.  If it's too big, it may get edges that aren't
      // really adjacent.
      double edge_tol = 0.1 * close_edge->length();
      if (is_outside_edge || dist > edge_tol)
        continue;

      // go merge the point with the edge
      was_merged = CUBIT_TRUE;

      if (mydebug)
      {
        DLIList<CubitFacet *> fl;
        pt->facets( fl );
        dcolor(CUBIT_YELLOW);
        dfldraw(fl);
        dcolor(CUBIT_BLUE);
        dpoint(pt->coordinates());
        CubitFacet *f = close_edge->adj_facet(0);
        dcolor(CUBIT_RED);
        dfdraw(f);
        dview();
        f = f;
      }

      // remove close_edge from the rtree
    
      r_tree.remove( close_edge );
      edge_shell_id = close_edge->marked();
      close_edge->marked(0);

      // split the edge

      assert(1 == close_edge->num_adj_facets());  // assumes we are splitting a boundary facet
      facet = close_edge->adj_facet(0);
      CubitFacetData *dfacet = dynamic_cast<CubitFacetData *>(facet);
      CubitPoint *new_pt = dfacet->split_edge(close_edge->point(0),
                                              close_edge->point(1),
                                              pt->coordinates());
      int facet_tool_id = facet->tool_id();
      new_pt->marked(edge_shell_id);

      // add any new facets to the shell and create missing edges

      for (ll=0; ll<edge_shell_id; ll++)
        shell_ptr = shell_list.get_and_step();
      DLIList<CubitFacet *>adj_facets;
      new_pt->facets( adj_facets );
      for (ll=0; ll<adj_facets.size(); ll++)
      {
        facet = adj_facets.get_and_step();
        if (!facet->marked())
        {
          facet->marked(edge_shell_id);
          shell_ptr->append(facet);
          for (mm=0; mm<3; mm++) {
            if (!(edge = facet->edge(mm)))
            {
              facet->get_edge_pts(mm, p0, p1);
              edge = (CubitFacetEdge *) new CubitFacetEdgeData( p0, p1 );
              edge->marked( 0 );
            }
          }
          facet->set_tool_id(facet_tool_id);
        }
      }

      // merge the points,
                                     
      merge_points( pt, new_pt, nemerge, &r_tree );
      npmerge++;

      // add any new edges to the rtree

      DLIList<CubitFacetEdge *> adj_edges;
      pt->edges( adj_edges );
      for (kk=0; kk<adj_edges.size(); kk++)
      {
        CubitFacetEdge *adj_edge = adj_edges.get_and_step();
        if (!adj_edge->marked() && adj_edge->num_adj_facets() == 1)
        {
          adj_edge->marked(pt->marked());
          r_tree.add(adj_edge);
        }
      }

      // see if shells need to merge and then merge

      pt_shell_id = pt->marked();
      if (pt_shell_id != edge_shell_id)
      {

        // get the shell containing the close point.  Nullify the
        // pointer in the shell list and move all of its facets
        // to the other shell.

        int delete_shell = edge_shell_id;
        DLIList<CubitFacet *> *delete_shell_ptr = NULL;
        shell_list.reset();
        for (ll=0; ll<delete_shell; ll++)
          delete_shell_ptr = shell_list.get_and_step();
        shell_list.back();
        shell_list.change_to(NULL);

        if(!delete_shell_ptr)
          return CUBIT_FAILURE;

        // mark all the points on the delete_shell as now being part of
        // the other shell.

        for(ll=0; ll<delete_shell_ptr->size(); ll++)
        {
          facet = delete_shell_ptr->get_and_step();
          facet->marked( pt_shell_id );
          facet->point( 0 )->marked( pt_shell_id );
          facet->point( 1 )->marked( pt_shell_id );
          facet->point( 2 )->marked( pt_shell_id );
          facet->edge( 0 )->marked( pt_shell_id );
          facet->edge( 1 )->marked( pt_shell_id );
          facet->edge( 2 )->marked( pt_shell_id );
        }

        // append the two lists together

        shell_list.reset();
        for (ll=0; ll<pt_shell_id; ll++)
          shell_ptr = shell_list.get_and_step();
        *shell_ptr += (*delete_shell_ptr);
        delete delete_shell_ptr;
      }

      // set the marked flag to negative to indicate that it has been
      // merged and it needs to be deleted.

      del_points.append(new_pt);
      new_pt->marked( -new_pt->marked() );
      pt->marked( -pt->marked() );
    } // close_edges

    if (!was_merged)
    {
      nnomerge++;
      if (mydebug)
      {
        dcolor(CUBIT_BLUE);
        dpdraw(pt);
        dcolor(CUBIT_RED);
        CubitVector mmin(-1e10, -1e10, -1e10);
        CubitVector mmax(1e10, 1e10, 1e10);
        ptbox.reset(mmin, mmax);
        r_tree.find(ptbox, close_edges);
        for(int zz=0; zz<close_edges.size(); zz++)
        {
          edge = close_edges.get_and_step();
          dedraw(edge);
          CubitVector loc_on_edge;
          CubitBoolean is_on_edge;
          double dist = edge->dist_to_edge(pt->coordinates(), loc_on_edge, is_on_edge);
          PRINT_INFO("edge %d dist = %f is_on_edge = %d\n", edge->id(), dist, is_on_edge);
        }
        dview();
        pt = pt;
      }
    }
  } // merge_point_list

  // compress the shell list

  shell_list.remove_all_with_value(NULL);

  // delete points that were merged

  for (ii=0; ii<del_points.size(); ii++)
  {
    pt = del_points.get_and_step();
    if (pt->marked() < 0)
    {
      assert( pt->num_adj_facets() == 0 );
      delete pt;
    }
  }

  return CUBIT_SUCCESS;
}

//=============================================================================
//Function:  merge_points (PRIVATE)
//Description: merge two cubit points where it has been determined that they
//             are coincident.  Also handle merging their adjacent edges where
//             appropriate
//Author: sjowen
//Date: 1/19/2004
//=============================================================================
CubitStatus FacetDataUtil::merge_points(
  CubitPoint *pt0, CubitPoint *pt1,
  int &nemerge, RTree <CubitFacetEdge*> *r_tree)  // number of edges we had to merge top do this
{
  int mydebug = 0;

  // merge the points
  pt0->merge_points( pt1, CUBIT_TRUE );
  //pt0->set_as_feature();

  if (mydebug)
  {
    DLIList<CubitFacet *>adj_facets;
    pt0->facets( adj_facets );
    dcolor(CUBIT_RED);
    dfldraw(adj_facets);
  }

  // merge edges

  int ll, mm;
  bool edge_merged;
  do {
    edge_merged = false;
    CubitFacetEdge *edge, *other_edge;

    DLIList<CubitFacetEdge *> adj_edges;
    pt0->edges( adj_edges );
    for(ll=0; ll<adj_edges.size() && !edge_merged; ll++)
    {
      edge = adj_edges.get_and_step();
      for(mm=0; mm<adj_edges.size() && !edge_merged; mm++)
      {
        other_edge = adj_edges.get_and_step();
        if (other_edge != edge &&
            edge->other_point(pt0) == other_edge->other_point(pt0))
        {
          CubitFacetEdgeData *dedge = dynamic_cast<CubitFacetEdgeData*>(edge);
          CubitFacetEdgeData *dother_edge = dynamic_cast<CubitFacetEdgeData*>(other_edge);
          // mbrewer (11/16/2005 for Bug 5049)
          if(r_tree){
//            r_tree->remove(dother_edge);
             CubitBoolean temp_bool = r_tree->remove(dother_edge);
             if(!temp_bool){
               PRINT_DEBUG_139("FacetDataUtil:  D_OTHER_EDGE did not exist in RTREE.\n");
             }
          }
          if (dedge->merge_edges( dother_edge ) == CUBIT_SUCCESS)
          {   
            nemerge++;
            edge_merged = true;
            dedge->set_as_feature();
          }
        }
      }
    }
  } while (edge_merged);
  return CUBIT_SUCCESS;
}

//=============================================================================
//Function:  merge_coincident_vertices (PRIVATE)
//Description: merge vertices (and connected facets and edges) if they are
//             within tolerance.
//Author: sjowen
//Date: 9/9/03
//=============================================================================
CubitStatus FacetDataUtil::merge_coincident_vertices(
  DLIList<DLIList<CubitFacet *> *> &shell_list,
  double tol,
  int &npmerge,    // return number of vertices merged
  int &nemerge,    // return number of edges merged
  DLIList<CubitPoint *> &unmerged_points)   // return the vertices on boundary NOT merged
{
  int mydebug = 0;
  npmerge = 0;
  nemerge = 0;
  int ii, jj, kk, ll, shell_id;
  KDDTree <CubitPoint*> kd_tree(tol, false);
  CubitPoint::set_box_tol( tol );
  shell_list.reset();
  DLIList<CubitFacet *> *shell_ptr = NULL;
  CubitPoint *pt;
  DLIList<CubitPoint *>boundary_points;
  DLIList<CubitPoint *> del_points;
  for (ii=0; ii<shell_list.size(); ii++)
  {
    shell_ptr = shell_list.get_and_step();

    // get the boundary points from the shell - these are candidates for merging.
    // note that this won't merge internal facets

    DLIList<CubitPoint *>shell_boundary_points;
    FacetDataUtil::get_boundary_points(*shell_ptr, shell_boundary_points);

    // mark each of the points with a shell id so we know when to merge shells
    // and add them to the kdtree for fast spatial searching

    shell_id = ii+1;
    for(jj=0; jj<shell_boundary_points.size(); jj++)
    {
      pt = shell_boundary_points.get_and_step();
      pt->marked(shell_id);
      kd_tree.add(pt);
      boundary_points.append(pt);
    }
  }
  kd_tree.balance();

  // find points that are coincident

  CubitBox ptbox;
  CubitVector ptmin, ptmax, coord;
  DLIList<CubitPoint *>close_points;
  CubitPoint *close_pt = NULL;
  CubitFacet *facet;
  DLIList<CubitPoint*>adj_pt_list;
  for(ii=0; ii<boundary_points.size(); ii++)
  {
    pt = boundary_points.get_and_step();
    if (pt->marked() < 0)  // has already been merged
      continue;

    // find the closest points in the kdtree

    coord = pt->coordinates();
    ptmin.set( coord.x() - tol, coord.y() - tol, coord.z() - tol );
    ptmax.set( coord.x() + tol, coord.y() + tol, coord.z() + tol );
    ptbox.reset( ptmin, ptmax );
    close_points.clean_out();
    kd_tree.find(ptbox, close_points);

    // if it didn't find anything to merge with, then we aren't water-tight

    CubitBoolean was_merged = CUBIT_FALSE;

    // We did find something - go try to merge

    for (jj=0; jj<close_points.size(); jj++)
    {
      close_pt = close_points.get_and_step();
      if (close_pt == pt)
        continue;
      if (close_pt->marked() < 0)  // has already been merged
        continue;
      assert(close_points.size() >= 1);

      // make sure this point is not already one of its neighbors
      // so we don't collapse a triangle

      CubitBoolean is_adjacent = CUBIT_FALSE;
      adj_pt_list.clean_out();
      close_pt->adjacent_points( adj_pt_list );
      for (kk=0; kk<adj_pt_list.size() && !is_adjacent; kk++)
      {
        if (adj_pt_list.get_and_step() == pt)
          is_adjacent = CUBIT_TRUE;
      }
      if (!is_adjacent)
      {
        // merge the points

        merge_points( pt, close_pt, nemerge );
        npmerge++;
        was_merged = CUBIT_TRUE;

        // see if shells need to merge and then merge

        shell_id = pt->marked();
        if (shell_id != close_pt->marked())
        {

          // get the shell containing the close point.  Nullify the
          // pointer in the shell list and move all of its facets
          // to the other shell.

          int delete_shell = close_pt->marked();
          DLIList<CubitFacet *> *delete_shell_ptr = NULL;
          shell_list.reset();
          for (ll=0; ll<delete_shell; ll++)
            delete_shell_ptr = shell_list.get_and_step();
          shell_list.back();
          shell_list.change_to(NULL);

          // mark all the points on the delete_shell as now being part of
          // the other shell.

          for(ll=0; ll<delete_shell_ptr->size(); ll++)
          {
            facet = delete_shell_ptr->get_and_step();
            facet->point( 0 )->marked( shell_id );
            facet->point( 1 )->marked( shell_id );
            facet->point( 2 )->marked( shell_id );
          }

          // append the two lists together

          shell_list.reset();
          for (ll=0; ll<shell_id; ll++)
            shell_ptr = shell_list.get_and_step();
          *shell_ptr += (*delete_shell_ptr);
        }

        // set the marked flag to negative to indicate that it has been
        // merged and it need to be deleted.

        if (close_pt->marked() > 0)
          close_pt->marked( -close_pt->marked() );
        del_points.append( close_pt );
      }
    }
    if (was_merged)
    {
      if (pt->marked() > 0)
        pt->marked( -pt->marked() );
    }
    else
    {
      // check to see if it was already merged
      if (close_points.size() == 1 && close_pt == pt)
      {
        CubitFacetEdge *edge_ptr;
        DLIList<CubitFacetEdge *>adj_edges;
        pt->edges(adj_edges);
        CubitBoolean on_boundary = CUBIT_FALSE;
        for(kk=0; kk<adj_edges.size() && !on_boundary; kk++)
        {
          edge_ptr = adj_edges.get_and_step();
          if (edge_ptr->num_adj_facets() == 1)
            on_boundary = CUBIT_TRUE;
        }
        if (on_boundary)
        {
            //PRINT_INFO("Merging 'boundary' points.\n");
            //if (pt->marked() > 0) pt->marked( -pt->marked() );
          unmerged_points.append(pt);
        }
        else
        {
          if (pt->marked() > 0) pt->marked( -pt->marked() );
        }
      }
      else
      {
        // otherwise save it as an unmerged point to be handled later
        unmerged_points.append(pt);
      }
    }
  }

  // compress the shell list

  shell_list.remove_all_with_value(NULL);

  // delete points that were merged

  for (ii=0; ii<del_points.size(); ii++)
  {
    pt = del_points.get_and_step();
    if (pt->marked() < 0)
    {
      assert( pt->num_adj_facets() == 0 );
      delete pt;
    }
  }

  if (mydebug)
  {
    CubitFacetEdge *edge;
    for (ii=0; ii<shell_list.size(); ii++)
    {
      shell_ptr = shell_list.get_and_step();
      dcolor(CUBIT_GREEN);
      dfldraw(*shell_ptr);
      dcolor(CUBIT_RED);
      for(jj=0; jj<shell_ptr->size(); jj++)
      {
        facet = shell_ptr->get_and_step();
        for(kk=0; kk<3; kk++)
        {
          edge = facet->edge(kk);
          DLIList<FacetEntity *> myfacet_list;
          edge->get_parents( myfacet_list );
          int np = myfacet_list.size();
          int na = edge->num_adj_facets();
          assert(np == na);
          if (myfacet_list.size() != 2)
          {
            dedraw( edge );
          }
        }
      }
    }
    dcolor(CUBIT_BLUE);
    dpldraw(unmerged_points);
    dview();
  }

  return CUBIT_SUCCESS;

}

//=============================================================================
//Function:  delete_facets (PUBLIC)
//Description: delete the facets and all associated edges and points from
//             a list of lists of facets (shell_list)
//Author: sjowen
//Date: 1/21/2004
//=============================================================================
void FacetDataUtil::delete_facets(DLIList<DLIList<CubitFacet*>*> &shell_list)
{
  int ii;
  for (ii=0; ii<shell_list.size(); ii++)
  {
    DLIList<CubitFacet*> *facet_list_ptr = shell_list.get_and_step();
    delete_facets( *facet_list_ptr );
  }
}

//=============================================================================
//Function:  delete_facets (PUBLIC)
//Description: delete the facets and all associated edges and points from
//             a list of facets
//Author: sjowen
//Date: 1/21/2004
//=============================================================================
void FacetDataUtil::delete_facets(DLIList<CubitFacet*> &facet_list)
{
  int ii;
  CubitFacet *facet_ptr;
  for (ii=0; ii<facet_list.size(); ii++)
  {
    facet_ptr = facet_list.get_and_step();
    delete_facet( facet_ptr );
  }
}

//=============================================================================
//Function:  delete_facet (PUBLIC)
//Description: delete a single facet and its underlying edges and points if
//             they are no longer attached to anything
//Author: sjowen
//Date: 1/21/2004
//=============================================================================
void FacetDataUtil::delete_facet(CubitFacet *facet_ptr)
{
  DLIList<CubitPoint *>point_list;
  DLIList<CubitFacetEdge *>edge_list;
  facet_ptr->points(point_list);
  facet_ptr->edges(edge_list);

  delete facet_ptr;

  CubitFacetEdge *edge_ptr;
  CubitPoint *point_ptr;
  int ii;

  for (ii=0; ii<edge_list.size(); ii++)
  {
    edge_ptr = edge_list.get_and_step();
    if (edge_ptr->num_adj_facets() == 0)
      delete edge_ptr;
  }

  for (ii=0; ii<3; ii++)
  {
    point_ptr = point_list.get_and_step();
    if (point_ptr->num_adj_facets() == 0)
      delete point_ptr;
  }

}


void FacetDataUtil::destruct_facet_no_delete(CubitFacet *facet_ptr)
{
  CubitFacetData* facet_d_ptr = dynamic_cast<CubitFacetData*>(facet_ptr);
  if(!facet_d_ptr){
    PRINT_ERROR("Can't work with Facet pointer that isn't a facet data object.\n");
    return;
  }
  
  DLIList<CubitPoint *>point_list;
  DLIList<CubitFacetEdge *>edge_list;
  facet_ptr->points(point_list);
  facet_ptr->edges(edge_list);
  
  facet_d_ptr->destruct_facet_internals();

  CubitFacetEdge *edge_ptr;
  CubitPoint *point_ptr;
  int ii;

  for (ii=0; ii<edge_list.size(); ii++)
  {
    edge_ptr = edge_list.get_and_step();
    if (edge_ptr->num_adj_facets() == 0)
      delete edge_ptr;
  }

  for (ii=0; ii<3; ii++)
  {
    point_ptr = point_list.get_and_step();
    if (point_ptr->num_adj_facets() == 0)
      delete point_ptr;
  }

}

//=============================================================================
//Function:  intersect_facet (PUBLIC)
//Description: determine intersection of a segment with a facet
//             returns CUBIT_PNT_UNKNOWN: if segment is in plane of facet
//                     CUBIT_PNT_OUTSIDE: if segment does not intersect
//                     CUBIT_PNT_INSIDE: if segment intersects inside facet
//                     CUBIT_PNT_BOUNDARY: if segment intersects a vertex or edge
//Author: sjowen
//Date: 1/30/2004
//=============================================================================
CubitPointContainment FacetDataUtil::intersect_facet(
  CubitVector &start, CubitVector &end, // start and end points of vector
  CubitFacet *facet_ptr,      // the facet to intersect
  CubitVector &qq,            // return the intersection point
  CubitVector &ac,    // area coordinates of qq if is in or on facet
  CubitBoolean bound) // if true, only check for intersections between the end points.
{

  CubitPlane fplane = facet_ptr->plane();
  double dstart = fplane.distance(start);
  double dend = fplane.distance(end);

  // points are both in the plane of the facet - can't handle this case

  if (fabs(dstart) < GEOMETRY_RESABS &&
      fabs(dend) < GEOMETRY_RESABS)
  {
    return CUBIT_PNT_UNKNOWN;
  }

  // one point is on the plane

  if (fabs(dstart) < GEOMETRY_RESABS)
  {
    qq = start;
  }
  else if (fabs(dend) < GEOMETRY_RESABS)
  {
    qq = end;
  }

  // points are both on the same side of the plane
  else if(dstart*dend > 0.0 &&
          (bound || fabs(dstart-dend) < CUBIT_RESABS) )
  {
    return CUBIT_PNT_OUTSIDE;
  }
  // points are on opposite sides of plane: if bound == false then compute intersection with plane

  else
  {
    CubitVector dir = end-start;
    dir.normalize();
    qq = fplane.intersect(start, dir);
  }

  FacetEvalTool::facet_area_coordinate( facet_ptr, qq, ac );
   
//mod mbrewer ... the original code would call a point
    // on the boundary if any area coordinate was near
    // zero, regardless of whether another area coordinate
    // was negative... making it outside.
//   if (fabs(ac.x()) < GEOMETRY_RESABS ||
//       fabs(ac.y()) < GEOMETRY_RESABS ||
//       fabs(ac.z()) < GEOMETRY_RESABS)
//   {
//     return CUBIT_PNT_BOUNDARY;
//   }
  if ( (fabs(ac.x()) < GEOMETRY_RESABS && (ac.y() > -GEOMETRY_RESABS &&
                                           ac.z() > -GEOMETRY_RESABS) )||
       (fabs(ac.y()) < GEOMETRY_RESABS && (ac.x() > -GEOMETRY_RESABS &&
                                           ac.z() > -GEOMETRY_RESABS) )||
       (fabs(ac.z()) < GEOMETRY_RESABS && (ac.x() > -GEOMETRY_RESABS &&
                                           ac.y() > -GEOMETRY_RESABS) ) ){
    return CUBIT_PNT_BOUNDARY;
  }   
  else if (ac.x() < 0.0 || ac.y() < 0.0 || ac.z() < 0.0)
  {
    return CUBIT_PNT_OUTSIDE;
  }

  return CUBIT_PNT_INSIDE;
}


//=============================================================================
//Function:  get_bbox_of_points (PUBLIC)
//Description: Find the bounding box of a list of CubitPoints
//Author: jdfowle
//Date: 12/15/03
//=============================================================================
CubitStatus FacetDataUtil::get_bbox_of_points(DLIList<CubitPoint*>& point_list, CubitBox& bbox)
{
double x, y, z, min[3], max[3];
int i;
CubitPoint *point;
  min[0] = min[1] = min[2] = CUBIT_DBL_MAX;
  max[0] = max[1] = max[2] = -CUBIT_DBL_MAX + 1.;

  for ( i = 0; i < point_list.size(); i++ ) {
    point = point_list.get_and_step();
    x = point->x();
    if ( min[0] > x ) min[0] = x;
    if ( max[0] < x ) max[0] = x;
    y = point->y();
    if ( min[1] > y ) min[1] = y;
    if ( max[1] < y ) max[1] = y;
    z = point->z();
    if ( min[2] > z ) min[2] = z;
    if ( max[2] < z ) max[2] = z;
  }
  bbox.reset(min,max);

  return CUBIT_SUCCESS;
}

//=============================================================================
//Function:  squared_distance_to_segment (PUBLIC)
//Description: get square of distance from point to closest point on a line segment;
//  returns closest point.  Taken from "Geometric Tools for Computer Graphics",
// by Schneider & Eberly, sec. 6.1.
//Author: jdfowle
//Date: 03/08/03
//=============================================================================
CubitVector FacetDataUtil::squared_distance_to_segment(CubitVector &p, CubitVector &p0, 
                                                      CubitVector &p1, double &distance2)
{
CubitVector YmPO, D;
double t;

  D = p1 - p0;
  YmPO = p - p0;
  t = D.x()*YmPO.x() + D.y()*YmPO.y() + D.z()*YmPO.z();
  
  if ( t < 0.0 ) {
    distance2 = YmPO.x()*YmPO.x() + YmPO.y()*YmPO.y() + YmPO.z()*YmPO.z();
    return p0;
  }
  
  double DdD;
  DdD = D.x()*D.x() + D.y()*D.y() + D.z()*D.z();
  if ( t >= DdD ) {
    CubitVector YmP1;
    YmP1 = p - p1;
    distance2 = YmP1.x()*YmP1.x() + YmP1.y()*YmP1.y() + YmP1.z()*YmP1.z();
    return p1;
  }
  
  //  closest point is interior to segment
  distance2 = YmPO.x()*YmPO.x() + YmPO.y()*YmPO.y() + YmPO.z()*YmPO.z() - t*t/DdD;
  return p0 + (t/DdD)*(p1 - p0);
}

//=============================================================================
//Function:  get_bbox_intersections (PUBLIC)
//Description: Get the intersection of the line defined by point1 and point2 with
//bbox.  Returns 0,1 or 2 for the number of intersections.  A line
//in one of the planes of the box will return 0.  Returns -1 for failure.
//Author mborden
//Date: 04/05/07
//=============================================================================
int FacetDataUtil::get_bbox_intersections(CubitVector& point1,
                                          CubitVector& point2,
                                          const CubitBox& bbox,
                                          CubitVector& intersection_1,
                                          CubitVector& intersection_2)
{
  int debug = 0;
  if( debug )
  {
    GfxDebug::draw_point( point1, CUBIT_RED );
    GfxDebug::draw_point( point2, CUBIT_BLUE );
    GfxDebug::flush();
  }
  
  double coords[6];
  coords[0] = bbox.min_x();
  coords[1] = bbox.max_x();
  coords[2] = bbox.min_y();
  coords[3] = bbox.max_y();
  coords[4] = bbox.min_z();
  coords[5] = bbox.max_z();

  DLIList<CubitVector*> intersections;
  
  int ii;
  for( ii = 0; ii < 3; ii++ )
  {
      //Define four points for each plane.
    double box_points[4][3];

      //ii = 0 -> x-planes
      //ii = 1 -> y_planes
      //ii = 2 -> z_planes

      //Only the coordinates for the plane we are in
      //change.  The other two are constant for the
      //jj loops below.
    box_points[0][(ii + 1) % 3] = coords[((ii*2)+2) % 6];
    box_points[1][(ii + 1) % 3] = coords[((ii*2)+3) % 6];
    box_points[2][(ii + 1) % 3] = coords[((ii*2)+3) % 6];
    box_points[3][(ii + 1) % 3] = coords[((ii*2)+2) % 6];

    box_points[0][(ii + 2) % 3] = coords[((ii*2)+4) % 6];
    box_points[1][(ii + 2) % 3] = coords[((ii*2)+4) % 6];
    box_points[2][(ii + 2) % 3] = coords[((ii*2)+5) % 6];
    box_points[3][(ii + 2) % 3] = coords[((ii*2)+5) % 6];
      
    int jj;
    for( jj = 0; jj < 2; jj++ )
    {
      CubitPoint* points[4];
      int kk;
      for( kk = 0; kk < 4; kk++ )
      {
        box_points[kk][ii] = coords[(ii*2)+jj];
        points[kk] = new CubitPointData( box_points[kk][0], box_points[kk][1], box_points[kk][2] );
      }
      
        //Create two facets for this plane to check for intersections.
      CubitFacet* facets[2];
      facets[0] = new CubitFacetData( points[0], points[1], points[3] );
      facets[1] = new CubitFacetData( points[1], points[2], points[3] );

      for( kk = 0; kk < 2; kk++ )
      {
        CubitVector intersection;
        CubitVector area_coord;

          //Make sure the points are not parrellel with the facet.
        CubitVector dir = point2 - point1;
        CubitVector normal = facets[kk]->normal();
        if( fabs(dir % normal) < CUBIT_RESABS )
            continue;
        
        CubitPointContainment contain = intersect_facet( point1, point2, facets[kk],
                                                         intersection, area_coord, CUBIT_FALSE );
        if( CUBIT_PNT_UNKNOWN == contain )
        {
            //The points are in a plane.  Return 0.
          delete facets[0];
          delete facets[1];
          int ll;
          for( ll = 0; ll < 4; ll++ )
              delete points[ll];
          for( ll = intersections.size(); ll > 0; ll-- )
              delete intersections.get_and_step();

          return 0;
        }
        if( CUBIT_PNT_BOUNDARY == contain ||
            CUBIT_PNT_INSIDE == contain )
        {
            //The point intersects the facet so it's inside the box's surface.
          CubitVector* new_intersection = new CubitVector;
          *new_intersection = intersection;
          intersections.append( new_intersection );

          if( debug )
          {
            GfxDebug::draw_point( *new_intersection, CUBIT_CYAN );
            GfxDebug::flush();
          }
          
          break;          
        }
      }
      
      delete facets[0];
      delete facets[1];
      for( kk = 0; kk < 4; kk++ )
          delete points[kk];      
    }
  }

    //Check for duplicate intersections.
  intersections.reset();
  for( ii = 0; ii < intersections.size(); ii++ )
  {
    CubitVector* base_vec = intersections.next(ii);
    if( NULL == base_vec )
        continue;
    
    int jj;
    for( jj = ii+1; jj < intersections.size(); jj++ )
    {
      CubitVector* compare_vec = intersections.next(jj);
      if( NULL != compare_vec )
      {
        if( base_vec->distance_between_squared( *compare_vec ) < GEOMETRY_RESABS * GEOMETRY_RESABS )
        {
          intersections.step(jj);
          delete intersections.get();
          intersections.change_to( NULL );
          intersections.reset();
        }
      }
    }
  }
  intersections.remove_all_with_value( NULL );

  
  if( intersections.size() > 2 )
  {
    assert( intersections.size() <= 2 );
    return -1;
  }
  else if( intersections.size() > 0 )
  {
    intersection_1 = *intersections.get();
    if( intersections.size() > 1 )
        intersection_2 = *intersections.next();
  }
  
    //Delete memory.
  for( ii = intersections.size(); ii > 0; ii-- )
      delete intersections.get_and_step();
  
  return intersections.size();
}


//=============================================================================
//Function:  mark_facets (PUBLIC)
//Description: mark facets and their children.  assumes facets have points and edges
//Author sjowen
//Date: 09/18/09
//=============================================================================
void FacetDataUtil::mark_facets( DLIList<FacetEntity *> &facet_list, int mark_value )
{
  int ifacet;
  FacetEntity *facet_ptr;
  CubitFacet *cfacet_ptr;
  for (ifacet = 0; ifacet<facet_list.size(); ifacet++)
  {
    facet_ptr = facet_list.get_and_step();
    cfacet_ptr = dynamic_cast<CubitFacet *> (facet_ptr);
    cfacet_ptr->marked(mark_value);
    for (int ii=0; ii<3; ii++)
    {
      cfacet_ptr->point(ii)->marked(mark_value);
      cfacet_ptr->edge(ii)->marked(mark_value);
    }
  }
  
}

//EOF

