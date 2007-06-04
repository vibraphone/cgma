//- Class:       ChollaCurve
//- Description: Temporary class for constructing the facet-based geometry
//-
//- Owner:       Steven J. Owen
//- Checked by:
//- Version:

#include "CubitVector.hpp"
#include "ChollaCurve.hpp"
#include "ChollaSurface.hpp"
#include "TDGeomFacet.hpp"
#include "CastTo.hpp"
#include "CubitFacet.hpp"
#include "CubitFacetData.hpp"
#include "CubitFacetEdge.hpp"
#include "CubitFacetEdgeData.hpp"
#include "debug.hpp"
#include "GfxDebug.hpp"

//===============================================================================
//Function:  ChollaCurve (PUBLIC) (constructor)
//===============================================================================
ChollaCurve::ChollaCurve( int block_id )
{
  static int count = 1;
  id = count++;
  myCurve = NULL;
  myEvalTool = NULL;
  startPoint = NULL;
  endPoint = NULL;
  blockID = block_id;
}
//===============================================================================
//Function:  ~ChollaCurve (PUBLIC) (destructor)
//===============================================================================
ChollaCurve::~ChollaCurve()
{
}

//===============================================================================
//Function:  remove_td_associativity (PUBLIC)
//===============================================================================
void ChollaCurve::remove_td_associativity( ChollaSurface *fsm_ptr )
{
  int i;
  TDGeomFacet *td;
  FacetEntity *edge_ptr;
  for(i=0; i<curveEdgeList.size(); i++)
  {
    edge_ptr = curveEdgeList.get_and_step();
    td = TDGeomFacet::get_geom_facet( edge_ptr );
    if (td)
    {
      td->remove_cholla_curve( this );
      td->remove_cholla_surf( fsm_ptr );
    }
  }
}


//=============================================================================
//Function:  get_ends (PUBLIC)
//Description: returns the end locations of the curve.  Determines their
//             location if not yet defined
//Author: sjowen
//Date: 12/4/00
//=============================================================================
CubitStatus ChollaCurve::get_ends( CubitVector &start, CubitVector &end )
{
  if (startPoint && endPoint)
  {
    start = startPoint->coordinates();
    end = endPoint->coordinates();
  }
  else
  {
    CubitStatus stat = determine_ends();
    if (stat == CUBIT_FAILURE)
      return stat;
    start = startPoint->coordinates();
    end = endPoint->coordinates();
  }
  return CUBIT_SUCCESS;
}

//=============================================================================
//Function:  get_ends (PUBLIC)
//Description: returns the end nodes of the curve.  Determines the
//             nodes if not yet defined
//Author: sjowen
//Date: 12/4/00
//=============================================================================
CubitStatus ChollaCurve::get_ends( CubitPoint *&start_ptr, CubitPoint *&end_ptr )
{
  if (startPoint && endPoint)
  {
    start_ptr = startPoint;
    end_ptr = endPoint;
  }
  else
  {
    CubitStatus stat = determine_ends();
    if (stat == CUBIT_FAILURE)
      return stat;
    start_ptr = startPoint;
    end_ptr = endPoint;
  }
  return CUBIT_SUCCESS;
}

//=============================================================================
//Function:  split_curve (PRIVATE)
//Description: split this curve into multiple ChollaCurve where there are
//             discontinuous strings of edges.  Define start and end nodes
//             for each curve while we are at it
//Author: sjowen
//Date: 12/4/00
//=============================================================================
CubitStatus ChollaCurve::split_curve(
  DLIList<ChollaCurve*> &facet_curve_list)
{
  DLIList<ChollaCurve*> new_curve_list;

  // Go through the curveEdgeList and pull edges off one by one as we
  // determine which curve it belongs to.  Continue until we have depleted
  // the list

  int periodic = 0;
  int start_size = curveEdgeList.size();
  int icount = 0;
  int mydebug = 0;

  curveEdgeList.reset();
  while( curveEdgeList.size() > 0)
  {

    // First, find an edge that has a start point on it

    CubitFacetEdge *start_edge_ptr = (CubitFacetEdge *)curveEdgeList.get_and_step();
    CubitPoint *point0_ptr = start_edge_ptr->point(0);
    CubitPoint *point1_ptr = start_edge_ptr->point(1);
    CubitPoint *start_point = NULL;
    if (periodic)
    {
      start_point = startPoint;
    }
    else
    {
      if (next_edge( point0_ptr, start_edge_ptr ) == NULL)
        start_point = point0_ptr;
      else if(next_edge( point1_ptr, start_edge_ptr ) == NULL)
        start_point = point1_ptr;
    }
    if (start_point != NULL || periodic)
    {

      // create a new curve to hold the edge info

      TDGeomFacet *td_gm_edge = TDGeomFacet::get_geom_facet(start_edge_ptr);
      int block_id = (td_gm_edge == NULL) ? -1 : td_gm_edge->get_block_id();
      ChollaCurve *fcm_ptr = new ChollaCurve( block_id );
      new_curve_list.append( fcm_ptr );
      fcm_ptr->set_start( start_point );
      start_point->set_as_feature();

      fcm_ptr->add_facet( start_edge_ptr );
      curveEdgeList.remove( start_edge_ptr );
      int iedgecount = 0;
      CubitFacetEdge *edge_ptr = start_edge_ptr;
      point0_ptr = start_point;
      CubitPoint *end_point = NULL;
      while(!end_point)
      {
        point1_ptr = edge_ptr->other_point( point0_ptr );
        if ((edge_ptr = next_edge( point1_ptr, edge_ptr )) == NULL)
        {
          end_point = point1_ptr;
        }
        else
        {
          iedgecount++;
          if (iedgecount > start_size)
          {
            PRINT_ERROR("ChollaCurve has start, but no end\n");
            return CUBIT_FAILURE;
          }

          fcm_ptr->add_facet( edge_ptr );
          curveEdgeList.remove( edge_ptr );
          if (periodic && point1_ptr == start_point)
            end_point = start_point;
          point0_ptr = point1_ptr;
        }
      }
      fcm_ptr->set_end( end_point );
      end_point->set_as_feature();
      start_size = curveEdgeList.size();
      icount = 0;
      periodic = 0;

      // make sure all the edges are oriented correctly

      int i;
      DLIList<FacetEntity *> flist = fcm_ptr->get_facet_list();
      flist.reset();
      DLIList<CubitFacetEdge *> elist;
      CAST_LIST( flist, elist, CubitFacetEdge );
      elist.reset();
      CubitPoint *cur_pt = start_point, *tmp_pt;
      for ( i = elist.size(); i > 0; i-- ) 
      {
        edge_ptr = elist.get_and_step();
        point0_ptr = edge_ptr->point(0);
        point1_ptr = edge_ptr->point(1);
        if (point0_ptr != cur_pt)
        {
          assert( cur_pt == point1_ptr );
          edge_ptr->flip();
          tmp_pt = point0_ptr;
          point0_ptr = point1_ptr;
          point1_ptr = tmp_pt;
          assert( point0_ptr == edge_ptr->point(0) &&
                  point1_ptr == edge_ptr->point(1) );
        }
        cur_pt = point1_ptr;
      }

      if (mydebug)
      {
        int i;
        DLIList<FacetEntity *> flist = fcm_ptr->get_facet_list();
        flist.reset();
        DLIList<CubitFacetEdge *> elist;
        CAST_LIST( flist, elist, CubitFacetEdge );
        elist.reset();
        for ( i = elist.size(); i > 0; i-- ) {  
          CubitFacetEdge *edge = elist.get_and_step();
          CubitVector pt0_v = edge->point(0)->coordinates();
          CubitVector pt1_v = edge->point(1)->coordinates();
          GfxDebug::draw_point(pt0_v, CUBIT_GREEN );
          GfxDebug::draw_point(pt1_v, CUBIT_RED );
          GfxDebug::draw_line( pt0_v, pt1_v, CUBIT_YELLOW );
          GfxDebug::flush();
          int view = 0;
          if (view)
            dview();
        }
      }
    }

    // if we have gone through all of the edges without finding an end,
    // then we have a periodic curve.  Choose an arbirary node to act as
    // the beginning and end

    if (curveEdgeList.size() > 0)
    {
      icount++;
      if (icount > start_size)
      {
        curveEdgeList.reset();
        CubitFacetEdge *edge = (CubitFacetEdge *)curveEdgeList.get();
        CubitPoint *point_ptr = edge->point(0);
        startPoint = point_ptr;
        endPoint = point_ptr;
        periodic = 1;
      }
    }
  }

  // add the new curves to the global curve list

  int ii, jj;
  for (ii=new_curve_list.size(); ii>0; ii--)
  {
    ChollaCurve *fcm_ptr = new_curve_list.get_and_step();

    facet_curve_list.append( fcm_ptr );

    // update the surface info

    for (jj=surfaceList.size(); jj>0; jj--)
    {
      ChollaSurface *fsm_ptr = surfaceList.get_and_step();
      fcm_ptr->add_surface( fsm_ptr );
      fsm_ptr->remove_curve( this );
      fsm_ptr->add_curve( fcm_ptr );
    }

    // update the geometric curve pointer

    fcm_ptr->assign_geometric_curve( NULL );

    // update the curve pointers in the edge tool data

    DLIList<FacetEntity*> facet_list = fcm_ptr->get_facet_list();
    for (jj=facet_list.size(); jj > 0; jj--)
    {
      FacetEntity *edge_ptr = facet_list.get_and_step();
      TDGeomFacet *td_gm_edge = TDGeomFacet::get_geom_facet(edge_ptr);
      td_gm_edge->remove_cholla_curve( this );
      td_gm_edge->add_cholla_curve( fcm_ptr );
    }
  }

  return CUBIT_SUCCESS;
}


//=============================================================================
//Function:  next_edge (PRIVATE)
//Description: determine the next edge from a given edge - return NULL if at
//             the end
//Author: sjowen
//Date: 12/4/00
//=============================================================================
CubitFacetEdge *ChollaCurve::next_edge( CubitPoint *node_ptr,
                                        CubitFacetEdge *edge_ptr )
{
  // check if this node has its hit flag set - we are at a feature break.

  TDGeomFacet *td_gm_node = TDGeomFacet::get_geom_facet(node_ptr);
  if (td_gm_node->get_hit_flag() == 1 || node_ptr->is_feature())
    return NULL;

  int jj, kk;

  // find the next edge

  CubitFacetEdge *next_edge_on_curve = NULL;
  DLIList<CubitFacetEdge*> edge_list;
  node_ptr->edges( edge_list );
  int num_adj_curves = 1;  // keep track of the number of curves at this node
  for (jj=0; jj<edge_list.size(); jj++)
  {
    CubitFacetEdge *node_edge_ptr = edge_list.get_and_step();
    if (node_edge_ptr != edge_ptr)
    {
      TDGeomFacet *td_gm_edge = TDGeomFacet::get_geom_facet(node_edge_ptr);
      if (td_gm_edge != NULL)
      {
        DLIList<ChollaCurve*> fcurve_list;
        td_gm_edge->get_cholla_curves( fcurve_list );
        if (fcurve_list.size() > 0)
        {

          // if 3 or more curves meet at this node, then force the curve to terminate here

          num_adj_curves++;
          if (num_adj_curves >= 3)
          {
            return NULL;
          }

          // otherwise try to match the curve to the edge to find the next edge

          for (kk=0; kk<fcurve_list.size(); kk++)
          {
            ChollaCurve *fcm_ptr = fcurve_list.get_and_step();
            if (fcm_ptr == this)
            {
              next_edge_on_curve = node_edge_ptr;
            }
          }
        }
      }
    }
  }
  return next_edge_on_curve;
}


//=============================================================================
//Function:  determine_ends (PRIVATE)
//Description: determine the end nodes of the curve
//             Assumes that there is one continuous string of edges
//             (may need to call split_curve first)
//Author: sjowen
//Date: 12/4/00
//=============================================================================
CubitStatus ChollaCurve::determine_ends( )
{
  int ii, jj, kk, inode;
  CubitFacetEdge *edge_ptr;
  CubitPoint *node0_ptr, *node1_ptr, *node_ptr;
  startPoint = endPoint = NULL;
  int done = 0;
  for(ii=0; ii<curveEdgeList.size() && !done; ii++)
  {
    edge_ptr = (CubitFacetEdge *)curveEdgeList.get_and_step();
    node0_ptr = edge_ptr->point(0);
    node1_ptr = edge_ptr->point(1);
    for (inode=0; inode<2 && !done; inode++)
    {
      node_ptr = (inode==0) ? node0_ptr : node1_ptr;
      DLIList<CubitFacetEdge*> edge_list;
      node_ptr->edges( edge_list );
      for (jj=0; jj<edge_list.size() && !done; jj++)
      {
        CubitFacetEdge *node_edge_ptr = edge_list.get_and_step();
        if (node_edge_ptr != edge_ptr)
        {
          TDGeomFacet *td_gm_edge = TDGeomFacet::get_geom_facet(node_edge_ptr);
          if (td_gm_edge != NULL)
          {
            int found = 0;
            DLIList<ChollaCurve*> fcurve_list;
            td_gm_edge->get_cholla_curves( fcurve_list );
            for (kk=0; kk<fcurve_list.size() && !found; kk++)
            {
              ChollaCurve *fcm_ptr = fcurve_list.get_and_step();
              if (fcm_ptr == this)
                found = 1;
            }
            if (!found)
            {
              if (startPoint == NULL)
              {
                startPoint = node_ptr;
              }
              else
              {
                endPoint = node_ptr;
                done = 1;
              }
            }
          }
        }
      }
    }

  }

  // check for periodic condition - just choose an arbitrary node to serve as
  // both start and end of the curve

  if (startPoint == NULL && endPoint == NULL)
  {
    curveEdgeList.reset();
    edge_ptr = (CubitFacetEdge *)curveEdgeList.get();
    node_ptr = edge_ptr->point(0);
    startPoint = node_ptr;
    endPoint = node_ptr;
  }
  else if (startPoint == NULL && endPoint != NULL ||
           startPoint != NULL && endPoint == NULL)
  {
    PRINT_ERROR("Could not determine start and end of curve in ChollaCurve\n");
    return CUBIT_FAILURE;
  }
  return CUBIT_SUCCESS;
}


//=============================================================================
//Function:  feature_angle (PRIVATE)
//Description: compute angles at nodes on the curve to see if we need to split
//             the curve.  Mark the node tooldata hitflag if the node will
//             break the curve (this is refernced in next_edge)
//Author: sjowen
//Date: 12/4/00
//=============================================================================
CubitStatus ChollaCurve::feature_angle(
  double min_dot )
{
  // first compute all of the edge vector and store with the edge tooldata

  int ii, jj;
  FacetEntity *facet_ptr;
  CubitFacetEdge *edge_ptr;
  CubitPoint *start_node;
  CubitPoint *end_node;
  CubitVector tangent;
  TDGeomFacet *td_gm;
  for (ii=0; ii<curveEdgeList.size(); ii++)
  {

    // compute the tangent vector of the edge and store it with its tooldata

    facet_ptr = curveEdgeList.get_and_step();
    edge_ptr = CAST_TO( facet_ptr, CubitFacetEdge );
    start_node = edge_ptr->point(0);
    end_node = edge_ptr->point(1);
    tangent = end_node->coordinates() - start_node->coordinates();
    tangent.normalize();
    td_gm = TDGeomFacet::get_geom_facet( edge_ptr );
    td_gm->set_normal( tangent );

    // initialize the nodes tooldata hit flags - set them all to -1

    td_gm = TDGeomFacet::get_geom_facet(start_node);
    td_gm->set_hit_flag(-1);
    td_gm = TDGeomFacet::get_geom_facet(end_node);
    td_gm->set_hit_flag(-1);
  }

  // now go through them again and compute the dot product between edges

  CubitVector tang0;
  CubitVector tang1;
  double dot;
  CubitPoint *node_ptr;
  CubitFacetEdge *next_edge_ptr;
  TDGeomFacet *td_gm_node;

  for (ii=0; ii<curveEdgeList.size(); ii++)
  {
    facet_ptr = curveEdgeList.get_and_step();
    edge_ptr = CAST_TO( facet_ptr, CubitFacetEdge );
    start_node = edge_ptr->point(0);
    end_node = edge_ptr->point(1);
    for (jj=0; jj<2; jj++)
    {
      node_ptr = (jj==0) ? start_node : end_node;
      td_gm_node = TDGeomFacet::get_geom_facet( node_ptr );
      if (td_gm_node->get_hit_flag() == -1)
      {
        next_edge_ptr = next_edge( node_ptr, edge_ptr );
        if (next_edge_ptr == NULL)
        {
          td_gm_node->set_hit_flag( 1 );
          node_ptr->set_as_feature();
        }
        else
        {
          td_gm = TDGeomFacet::get_geom_facet( edge_ptr );
          tang0 = td_gm->get_normal();
          td_gm = TDGeomFacet::get_geom_facet( next_edge_ptr );
          tang1 = td_gm->get_normal();

          // change the sign of the tangent vectors if the
          // sense of the edges are not the same

          if (node_ptr == start_node)
          {
            if (node_ptr != next_edge_ptr->point(1))
              tang0 = -tang0;
          }
          else
          {
            if (node_ptr != next_edge_ptr->point(0))
              tang0 = -tang0;
          }

          // compute the dot product between tangemt vectors

          dot = tang0 % tang1;

          // set the hit flag if there needs to be a feature break here

          if (dot <= min_dot)
          {
            td_gm_node->set_hit_flag( 1 );
            node_ptr->set_as_feature();
          }
          else
          {
            td_gm_node->set_hit_flag( 0 );
          }
        }
      }
    }
  }
  return CUBIT_SUCCESS;
}

static int icolor = 0;
void ChollaCurve::debug_draw()
{
  icolor++;
  icolor = icolor%15;
  dcolor(icolor);
  dldraw(curveEdgeList);
}

void ChollaCurve::print()
{
  FILE *fp = fopen("debug.curve", "a");
  fprintf(fp,"*** Curve %d ***\n", id);
  for (int ii=0; ii<curveEdgeList.size(); ii++)
  {
    FacetEntity *fe_ptr = curveEdgeList.get_and_step();
    CubitFacetEdge *cfe_ptr = CAST_TO(fe_ptr, CubitFacetEdge );
    CubitPoint *cp0_ptr = cfe_ptr->point(0);
    CubitPoint *cp1_ptr = cfe_ptr->point(1);
    fprintf(fp,"  Edge (%d)\n", cfe_ptr->id() );
    fprintf(fp,"     Point (%d)  %8.4f  %8.4f  %8.4f\n",
      cp0_ptr->id(), cp0_ptr->x(), cp0_ptr->y(), cp0_ptr->z());
    fprintf(fp,"     Point (%d)  %8.4f  %8.4f  %8.4f\n",
      cp1_ptr->id(), cp1_ptr->x(), cp1_ptr->y(), cp1_ptr->z());
  }
  fclose(fp);
}

//EOF

