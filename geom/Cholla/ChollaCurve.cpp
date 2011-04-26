//- Class:       ChollaCurve
//- Description: Temporary class for constructing the facet-based geometry
//-
//- Owner:       Steven J. Owen
//- Checked by:
//- Version:

#include "CubitVector.hpp"
#include "ChollaCurve.hpp"
#include "ChollaSurface.hpp"
#include "ChollaPoint.hpp"
#include "TDGeomFacet.hpp"
#include "CastTo.hpp"
#include "CubitFacet.hpp"
#include "CubitFacetData.hpp"
#include "CubitFacetEdge.hpp"
#include "CubitFacetEdgeData.hpp"
#include "debug.hpp"
#include "GfxDebug.hpp"
#include "CurveFacetEvalTool.hpp"
#include "ChollaEngine.hpp"

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
  myLength = MYLENGTH_UNINITIALIZED;
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
      
      // assign the edges to the new curve in the correct order and orientation
      
      CubitStatus rv = fcm_ptr->build_curve_from_edges( start_point,  periodic, start_size, start_edge_ptr, this );
      if (rv != CUBIT_SUCCESS)
        return rv;
      
      // remove the edges in the new curve from this curve
      
      int ii;
      DLIList<FacetEntity *> flist = fcm_ptr->get_facet_list();
      DLIList<CubitFacetEdge *> elist;
      CubitFacetEdge *edge_ptr;
      CAST_LIST( flist, elist, CubitFacetEdge );
      for ( ii = elist.size(); ii > 0; ii-- ) 
      {
        edge_ptr = elist.get_and_step();
        curveEdgeList.remove( edge_ptr );
      }
      start_size = curveEdgeList.size();
      icount = 0;
      periodic = 0;
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
//Function: build_curve_from_edges 
//Description: insert the ordered and oriented edges into this cholla curve
//Notes:  traverses starting at start_point and gathers facet edges until it 
//        runs into another curve.
//        start_point is an existing CubitPoint at either end of the curve
//        max_edges is the maximum number of edges on this curve.  should be 
//        known beforehand (used for error checking).
//
// ***this function used to be part of split_curve. ***
//Author: sjowen
//Return: 
//Date: 09/07/2009
//=============================================================================
CubitStatus ChollaCurve::build_curve_from_edges( CubitPoint *start_point,
                                                int periodic,
                                                int max_edges,
                                                CubitFacetEdge *start_edge_ptr,
                                                ChollaCurve *parent_curve )
{
  
  // find the first edge.  Match the chollacurve owner with this curve
  // do this only if the start_edge_ptr was not passed in
  
  DLIList<CubitFacetEdge *> point_edge_list;
  start_point->edges(point_edge_list);
  CubitFacetEdge *edge_ptr;
  if (start_edge_ptr == NULL)
  {
    for (int ii=0; ii<point_edge_list.size() && !start_edge_ptr; ii++)
    {
      edge_ptr = point_edge_list.get_and_step();
      TDGeomFacet *td_geom = TDGeomFacet::get_geom_facet( edge_ptr );
      
      // assumes that the TDGeomFacet info has already been set up for the edges
      assert(td_geom != NULL);
      
      DLIList<ChollaCurve *> cholla_curves;
      td_geom->get_cholla_curves(cholla_curves);
      
      // currently should be only one-to-one relationship
      // could also be edge on surface in which case no curves associated
      assert(cholla_curves.size() <= 1);
      if (cholla_curves.size())
      {
        if (cholla_curves.get() == this)
          start_edge_ptr = edge_ptr;
      }
    }
    assert(start_edge_ptr != NULL);  // didn't find an edge that marched this chollacurve
  }
  
  // create a new curve to hold the edge info
  
  this->set_start( start_point );
  start_point->set_as_feature();
  
  this->add_facet( start_edge_ptr );
  int iedgecount = 0;
  edge_ptr = start_edge_ptr;
  CubitPoint *point0_ptr = start_point, *point1_ptr;
  CubitPoint *end_point = NULL;
  while(!end_point)
  {
    point1_ptr = edge_ptr->other_point( point0_ptr );
    if ((edge_ptr = parent_curve->next_edge( point1_ptr, edge_ptr )) == NULL)
    {
      end_point = point1_ptr;
    }
    else
    {
      iedgecount++;
      if (iedgecount > max_edges)
      {
        PRINT_ERROR("ChollaCurve has start, but no end\n");
        return CUBIT_FAILURE;
      }
      
      this->add_facet( edge_ptr );
      if (periodic && point1_ptr == start_point)
        end_point = start_point;
      point0_ptr = point1_ptr;
    }
  }
  this->set_end( end_point );
  end_point->set_as_feature();
  
  // make sure all the edges are oriented correctly
  
  int i;
  DLIList<FacetEntity *> flist = this->get_facet_list();
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
  
  int mydebug = 0;
  if (mydebug)
  {
    int i;
    DLIList<FacetEntity *> flist = this->get_facet_list();
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
  return CUBIT_SUCCESS;
}

//=============================================================================
//Function: length
//Description: 
//Author: sjowen
//Date: 04/21/2009
//============================================================================
double ChollaCurve::length()
{
  if (myLength > MYLENGTH_UNINITIALIZED)
    return myLength;
  
  CubitFacetEdge *edge;
  FacetEntity *fent;
  myLength = 0.0;
  for (int iedge=0; iedge<curveEdgeList.size(); iedge++)
  {
    fent = curveEdgeList.get_and_step();
    edge = dynamic_cast<CubitFacetEdge *> (fent);
    myLength += edge->length();
  }
  
  return myLength;
}

//=============================================================================
//Function: find adjacent edges at a point that lie on the curve 
//Description: determine the next edge from a given edge - return NULL if at the end
//Author: william roshan quadros
//Return: returns false if no adj_edges can be found
//Date: 04/21/2009
//=============================================================================
bool ChollaCurve::adj_facet_edges( CubitPoint *node_ptr, CubitFacetEdge *&adj_edge1, CubitFacetEdge *&adj_edge2 )
{
  // initialize adj_edge1 and adj_edge2
  adj_edge1 = adj_edge2 = NULL;

  DLIList<CubitFacetEdge*> edge_list;
  node_ptr->edges( edge_list );
  int jj, kk;
  for (jj=0; jj<edge_list.size(); jj++)
  {
    CubitFacetEdge *node_edge_ptr = edge_list.get_and_step();
    TDGeomFacet *td_gm_edge = TDGeomFacet::get_geom_facet(node_edge_ptr);
    if (td_gm_edge != NULL)
    {
      DLIList<ChollaCurve*> fcurve_list;
      td_gm_edge->get_cholla_curves( fcurve_list );
      if (fcurve_list.size() > 0)
      { // match the curve to the edge to find the next edge
        for (kk=0; kk<fcurve_list.size(); kk++)
        {
          ChollaCurve *fcm_ptr = fcurve_list.get_and_step();
          if (fcm_ptr == this)
          {
            if( NULL == adj_edge1 )
              adj_edge1 = node_edge_ptr;
            else 
              if( NULL == adj_edge2 )
                adj_edge2 = node_edge_ptr;
              else
                assert( false ); // More than two adj_edges can't be incident on a curve
          }
        }
      }
    }    
  }
  if( NULL == adj_edge1 )
    return false;
  else
    return true;
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


  //  disassociate from cholla points
CubitStatus ChollaCurve::disassociate_from_points( void )
{  

  /*
  if( startPoint )
  {
    startPoint->remove_curve( this );
    startPoint = NULL;
  }
  if( endPoint )
  {
    endPoint->remove_curve( this );
    endPoint = NULL;
  }
  */

  int i; 
  for( i = 0; i < pointList.size(); i++ )
  {
    pointList.get_and_step()->remove_curve( this );
  }
  pointList.clean_out();

  return CUBIT_SUCCESS;
}
  
  // disassociate from cholla surface
CubitStatus ChollaCurve::disassociate_from_surfaces( void)
{
  int i;
  for( i = 0; i < surfaceList.size(); i++ )
  {
    surfaceList.get_and_step()->remove_curve( this );
  }
  surfaceList.clean_out();
  return CUBIT_SUCCESS;
}


CubitStatus ChollaCurve::replace_facet( FacetEntity *remove_edge, FacetEntity *replace_edge )
{
  curveEdgeList.move_to( remove_edge );
  curveEdgeList.insert( replace_edge );
  curveEdgeList.remove( remove_edge );
  myLength = MYLENGTH_UNINITIALIZED;
  return CUBIT_SUCCESS;
}

CubitStatus ChollaCurve::build_curve_facet_eval_tool( void )
{ 
  //debug this function as point list is not valid

  if( this->get_eval_tool() )
  {
    assert(false); //WARNING: Curve facet eval tool already exist!
  }

  CurveFacetEvalTool *curv_eval_tool_ptr = new CurveFacetEvalTool();

  // Step 1: Initialize facet_edge_list and point_list
  CubitStatus stat;
  DLIList<CubitPoint *> point_list;
  int i;
  // insert start point of every facet_edge
  curveEdgeList.reset();
  for( i = 0; i < curveEdgeList.size(); i++ )
  {
    point_list.append( CAST_TO( curveEdgeList.get_and_step(), CubitFacetEdge )->point(0) );
  }
  // insert end point of last facet_edge
  curveEdgeList.step( curveEdgeList.size() - 1 );
  point_list.append( CAST_TO( curveEdgeList.get(), CubitFacetEdge )->point(1) );

  DLIList<CubitFacetEdge *> edge_list;
  CAST_LIST( curveEdgeList, edge_list, CubitFacetEdge );
  stat = curv_eval_tool_ptr->initialize( edge_list, point_list );
  if( stat != CUBIT_SUCCESS )
  {
    return stat;
  }
  
 /*
 // Step 2: find sense of curve_facet_eval_tool  /// this is done internally in next Step in initialize()
  if( this->startPoint )
  {
    stat = curv_eval_tool_ptr->find_curv_sense( this->startPoint );
    if( stat != CUBIT_SUCCESS )
    {
      return stat;
    }  
  }
*/
  // Step 2: Initialize adj_surface_facet_eval_tool with orientation_wrt_surface
  if( surfaceList.size() )
  {
    CubitSense orientation_wrt_surface;
    if( CUBIT_SUCCESS == ChollaEngine::determine_curve_orientation( surfaceList.get(), this, orientation_wrt_surface ) )
    {
      if( this->startPoint && this->endPoint )
      {
        stat = curv_eval_tool_ptr->initialize( surfaceList.get()->get_eval_tool(),
                      this->startPoint,
                      this->endPoint,
                      orientation_wrt_surface);
      }
    }
    else
    {
      assert(false);
    }

    if( stat != CUBIT_SUCCESS )
    {
      assert( false );
      return stat;
    }      
  }
  else
  {
    assert(false); //WARNING: No adjacent cholla surface available
  }

  // Step 4: assign the new curv_eval_tool to cholla_curve
 assign_eval_tool( curv_eval_tool_ptr );

  return stat;
}

//=============================================================================
//Function:  is_in_volume (PUBLIC)
//Description:  return whether this curve is contained within the specified volume
//Author: sjowen
//Date: 9/11/2009
//=============================================================================
CubitBoolean ChollaCurve::is_in_volume( ChollaVolume *chvol_ptr )
{
  for (int ii=0; ii<surfaceList.size(); ii++)
  {
    ChollaSurface *chsurf_ptr = surfaceList.get_and_step();
    DLIList<ChollaVolume *> chvol_list;
    chsurf_ptr->get_volumes(chvol_list);
    for (int jj=0; jj<chvol_list.size(); jj++)
    {
      ChollaVolume *mychvol_ptr = chvol_list.get_and_step();
      if (mychvol_ptr == chvol_ptr)
        return CUBIT_TRUE;
    }
  }
  return CUBIT_FALSE;
}

//=============================================================================
//Function:  is_in_surface (PUBLIC)
//Description:  return whether this curve is contained within the specified surface
//Author: sjowen
//Date: 9/18/2009
//=============================================================================
CubitBoolean ChollaCurve::is_in_surface( ChollaSurface *chsurf_ptr )
{
  for (int ii=0; ii<surfaceList.size(); ii++)
  {
    ChollaSurface *mysurf_ptr = surfaceList.get_and_step();
    if (mysurf_ptr == chsurf_ptr)
    {
      return CUBIT_TRUE;
    }
  }
  return CUBIT_FALSE;
}


//=============================================================================
//Function:  get_facet_points (PUBLIC)
//Description:  return the list of facet points on this chollacurve
//Notes: inclusive = true will return end points as well, otherwise only
//       interior points will be returned
//Author: sjowen
//Date: 9/11/2009
//=============================================================================
void ChollaCurve::get_facet_points( DLIList<CubitPoint *> &point_list, CubitBoolean inclusive)
{
  FacetEntity *fe_ptr;
  CubitFacetEdge *edge_ptr;
  CubitPoint *pts[2];
  for (int ii=0; ii<curveEdgeList.size(); ii++)
  {
    fe_ptr = curveEdgeList.get_and_step();
    edge_ptr = dynamic_cast<CubitFacetEdge *> (fe_ptr);
    assert(edge_ptr != NULL);
    for (int jj=0; jj<2; jj++)
    {
      pts[jj] = edge_ptr->point(jj);
      if (inclusive)
      {
        point_list.append(pts[jj]);
      }
      else
      {
        if (pts[jj] != startPoint && pts[jj] != endPoint)
        {
          point_list.append(pts[jj]);
        }
      }
    }
  }
  point_list.uniquify_ordered();
}

//=============================================================================
//Function:  has_point (PUBLIC)
//Description:  return whether the curve contains the gicen point
//Notes: 
//=============================================================================
CubitBoolean ChollaCurve::has_point( ChollaPoint *chpt )
{
  for(int ii=0; ii<pointList.size(); ii++)
  {
    ChollaPoint *pt = pointList.get_and_step();
    if (pt == chpt)
      return CUBIT_TRUE;
  }
  return CUBIT_FALSE;
}

//=============================================================================
//Function:  verify_points (PUBLIC)
//Description:  verify that all points on this curve have this curve as an adjacency
//Notes: 
//=============================================================================
CubitStatus ChollaCurve::verify_points()
{
  for(int ii=0; ii<pointList.size(); ii++)
  {
    ChollaPoint *pt = pointList.get_and_step();
    if (!pt->is_in_curve(this))
      return CUBIT_FAILURE;
  }
  return CUBIT_SUCCESS;
}

//EOF

