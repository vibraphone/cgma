//- Class:       TDGeomFacet
//- Description: Tool data for building a geometry from a mesh.
//- Owner:       Steve Owen
//- Checked by:
//- Version:
#include "TDGeomFacet.hpp"
#include "CubitFacetEdge.hpp"
#include "CubitFacet.hpp"
#include "CastTo.hpp"
#include "ChollaVolume.hpp"
#include "ChollaSurface.hpp"
#include "ChollaCurve.hpp"
#include "ChollaPoint.hpp"
#include "CubitPoint.hpp"

TDGeomFacet::TDGeomFacet()
{
  blockId = -1;
  hitFlag = 0;
  partnerEdgeList = NULL;
  partnerPointList = NULL;
}

TDGeomFacet::~TDGeomFacet()
{
  if (partnerEdgeList != NULL)
    delete partnerEdgeList;
  if (partnerPointList != NULL)
    delete partnerPointList;
}

int TDGeomFacet::is_geom_facet(const ToolData* td)
{
  return (CAST_TO(const_cast<ToolData*>(td), TDGeomFacet) != NULL);
}

CubitStatus TDGeomFacet::add_geom_facet( FacetEntity *facet_ptr, int block_id )
{
  TDGeomFacet* td = (TDGeomFacet*) facet_ptr->get_TD( &TDGeomFacet::is_geom_facet );
  if ( td == NULL )
  {
    td = new TDGeomFacet;
    facet_ptr->add_TD( td );
  }
  td->set_block_id( block_id);
  return CUBIT_SUCCESS;
}

CubitStatus TDGeomFacet::add_geom_facet( CubitFacet *facet_ptr, int block_id)
{
  TDGeomFacet *td = (TDGeomFacet*) facet_ptr->get_TD(&TDGeomFacet::is_geom_facet);
  if ( td == NULL )
  {
    td = new TDGeomFacet;
    facet_ptr->add_TD( td );
  }
  td->set_block_id( block_id);
  return CUBIT_SUCCESS;
}
CubitStatus TDGeomFacet::add_geom_facet( CubitFacetEdge *edge_ptr, int block_id )
{
  TDGeomFacet *td = (TDGeomFacet*) edge_ptr->get_TD(&TDGeomFacet::is_geom_facet);
  if ( td == NULL )
  {
    td = new TDGeomFacet;
    edge_ptr->add_TD( td );
  }
  td->set_block_id( block_id);
  return CUBIT_SUCCESS;
}
CubitStatus TDGeomFacet::add_geom_facet( CubitPoint *point_ptr, int block_id )
{
  TDGeomFacet* td = (TDGeomFacet*) point_ptr->get_TD(&TDGeomFacet::is_geom_facet);
  if ( td == NULL )
  {
    td = new TDGeomFacet;
    point_ptr->add_TD( td );
  }
  td->set_block_id( block_id);

  return CUBIT_SUCCESS;
}

TDGeomFacet* TDGeomFacet::get_geom_facet( FacetEntity *facet_ptr )
{
  TDGeomFacet* td = (TDGeomFacet*) facet_ptr->get_TD(&TDGeomFacet::is_geom_facet);
  if ( td != NULL )
  {
    return td;
  }
  return (TDGeomFacet*) NULL;
}
TDGeomFacet* TDGeomFacet::get_geom_facet( CubitPoint *point_ptr )
{
  TDGeomFacet *td = (TDGeomFacet*) point_ptr->get_TD(&TDGeomFacet::is_geom_facet);
  if ( td != NULL )
  {
    return td;
  }
  return (TDGeomFacet*) NULL;
}
TDGeomFacet* TDGeomFacet::get_geom_facet( CubitFacetEdge *edge_ptr )
{
  TDGeomFacet *td = (TDGeomFacet*) edge_ptr->get_TD(&TDGeomFacet::is_geom_facet);
  if ( td != NULL )
  {
    return td;
  }
  return (TDGeomFacet*) NULL;
}
TDGeomFacet* TDGeomFacet::get_geom_facet( CubitFacet *facet_ptr )
{
  TDGeomFacet *td = (TDGeomFacet*) facet_ptr->get_TD(&TDGeomFacet::is_geom_facet);
  if ( td != NULL )
  {
    return td;
  }
  return (TDGeomFacet*) NULL;
}


int TDGeomFacet::get_block_id( FacetEntity *facet_ptr )
{
  TDGeomFacet *td = (TDGeomFacet*) facet_ptr->get_TD(&TDGeomFacet::is_geom_facet);
  if ( td != NULL )
  {
    return td->get_block_id();
  }
  return -1;
}

int TDGeomFacet::get_block_id( CubitFacet *facet_ptr )
{
  TDGeomFacet *td = (TDGeomFacet*) facet_ptr->get_TD(&TDGeomFacet::is_geom_facet);
  if ( td != NULL )
  {
    return td->get_block_id();
  }
  return -1;
}

int TDGeomFacet::get_block_id( CubitFacetEdge *edge_ptr )
{
  TDGeomFacet *td = (TDGeomFacet*) edge_ptr->get_TD(&TDGeomFacet::is_geom_facet);
  if ( td != NULL )
  {
    return td->get_block_id();
  }
  return -1;
}
void TDGeomFacet::add_cholla_owner( ChollaEntity *cholla_entity )
{
  ChollaSurface *cholla_surface = dynamic_cast<ChollaSurface *> (cholla_entity);
  if (cholla_surface != NULL)
    add_cholla_surf(cholla_surface);
  else
  {
    ChollaCurve *cholla_curve = dynamic_cast<ChollaCurve *> (cholla_entity);
    if (cholla_curve != NULL)
      add_cholla_curve( cholla_curve );
    else
    {
      ChollaPoint *cholla_point = dynamic_cast<ChollaPoint *> (cholla_entity);
      if (cholla_point != NULL)
        add_cholla_point( cholla_point );
      else
        assert(0); // not a recognized cholla entity
    }
  }  
  return;
}
void TDGeomFacet::add_cholla_surf( ChollaSurface *f_s_m )
{
  int ii;
  for ( ii = ChollaSurfaceList.size(); ii > 0; ii-- )
  {
    ChollaSurface *fsm_ptr = ChollaSurfaceList.get_and_step();
    if (f_s_m == fsm_ptr)
    {
      return;
    }
  }
  ChollaSurfaceList.append(f_s_m);
  return;
}
int TDGeomFacet::get_hit_flag( FacetEntity *facet_ptr )
{
  TDGeomFacet *td = (TDGeomFacet*) facet_ptr->get_TD(&TDGeomFacet::is_geom_facet);
  if ( td != NULL )
  {
    return td->get_hit_flag();
  }
  return -1;
}
void TDGeomFacet::set_hit_flag( FacetEntity *facet_ptr, int new_val )
{
  TDGeomFacet *td = (TDGeomFacet*) facet_ptr->get_TD(&TDGeomFacet::is_geom_facet);
  if ( td != NULL )
  {
    td->set_hit_flag(new_val);
  }
  return;
}

CubitBoolean TDGeomFacet::is_partner( CubitFacetEdge *edge_ptr )
{
  if (partnerEdgeList == NULL)
    return CUBIT_FALSE;
  for (int ii=0; ii<partnerEdgeList->size(); ii++)
  {
    if (partnerEdgeList->get_and_step() == edge_ptr)
      return CUBIT_TRUE;
  }
  return CUBIT_FALSE;
}
CubitBoolean TDGeomFacet::is_partner( CubitPoint *point_ptr )
{
  if (partnerPointList == NULL)
    return CUBIT_FALSE;
  for (int ii=0; ii<partnerPointList->size(); ii++)
  {
    if (partnerPointList->get_and_step() == point_ptr)
      return CUBIT_TRUE;
  }
  return CUBIT_FALSE;
}

void TDGeomFacet::reset_TD_as_new()
{
  blockId = -1;
  hitFlag = 0;
  partnerEdgeList = NULL;
  partnerPointList = NULL;
  ChollaSurfaceList.clean_out();
  ChollaCurveList.clean_out();
  myPoints.clean_out();
    //normal??????;
    //need to delete these?
  partnerEdgeList=NULL;
  partnerPointList=NULL;

}

int TDGeomFacet::geo_type()
{
  if (ChollaSurfaceList.size() > 0)
    return 2;
  else if (ChollaCurveList.size() > 0)
    return 1;
  else if (myPoints.size() > 0)
    return 0;
  return -1;
}

CubitBoolean TDGeomFacet::is_in_volume( ChollaVolume *chvol_ptr )
{
  for (int ii=0; ii<ChollaSurfaceList.size(); ii++)
  {
    ChollaSurface *chsurf_ptr = ChollaSurfaceList.get_and_step();
    if (chsurf_ptr->is_in_volume( chvol_ptr ))
      return CUBIT_TRUE;
  }
  for (int jj=0; jj<ChollaCurveList.size(); jj++)
  {
    ChollaCurve *chcurv_ptr = ChollaCurveList.get_and_step();
    if (chcurv_ptr->is_in_volume( chvol_ptr ))
      return CUBIT_TRUE;
  }
  return CUBIT_FALSE;
}



