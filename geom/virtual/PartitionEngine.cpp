//-------------------------------------------------------------------------
// Filename      : PartitionEngine.cpp
//
// Purpose       : 
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 04/23/02
//-------------------------------------------------------------------------

#include "PartitionEngine.hpp"
#include "CompositeEngine.hpp"

#include "SubCurve.hpp"
#include "SubSurface.hpp"
#include "SegmentedCurve.hpp"
#include "PartitionPoint.hpp"
#include "PartitionCoEdge.hpp"
#include "PartitionLoop.hpp"
#include "PartitionSurface.hpp"
#include "PartitionLump.hpp"
#include "PartitionBody.hpp"
#include "PartPTCurve.hpp"
#include "PartitionShell.hpp"
#include "BodySM.hpp"


#include "GMem.hpp"
#include "GeometryQueryEngine.hpp"
#include "GeometryModifyEngine.hpp"
#include "GeometryModifyTool.hpp"

#include "VGLoopTool.hpp"

#include "FacetProjectTool.hpp"
#include "CubitFacetData.hpp"
#include "CubitFacetEdgeData.hpp"
#include "CubitPointData.hpp"
#include "TDVGFacetOwner.hpp"
#include "TDVGFacetSplit.hpp"
#include "FacetEvalTool.hpp"
#include "FacetDataUtil.hpp"
#include "PartSurfFacetTool.hpp"

#include "CompositeSurface.hpp"
#include "CompositePoint.hpp"
#include "CompositeCurve.hpp"
#include "PartitionLumpImprint.hpp"
#include "GfxDebug.hpp"
#include "BridgeManager.hpp"

#include "CADefines.hpp"

#include <set>
#include "GeometryQueryTool.hpp"

#include "CubitTransformMatrix.hpp"
#include "CubitObservable.hpp"

const char* const PARTITION_GEOM_ATTRIB_NAME = "PARTITION_GEOM";

static CubitStatus get_edge_replacements(
                         std::vector<CubitFacetData*> &facet_list,
                         std::vector<CubitFacetEdgeData*> replacement_edges[3]);


PartitionEngine* PartitionEngine::instance_ = NULL; 


PartitionEngine::~PartitionEngine()
{
  GeometryQueryTool::instance()->unregister_intermediate_engine(this);
}

void PartitionEngine::delete_instance()
{
  if( NULL != instance_ )
  {
    delete instance_;
    instance_ = NULL;
  }
}    


//-------------------------------------------------------------------------
// Purpose       : Constructor
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 08/25/03
//-------------------------------------------------------------------------
PartitionEngine::PartitionEngine()
{
  CubitStatus r = GeometryQueryTool::instance()->register_intermediate_engine(this);
  assert(r == CUBIT_SUCCESS);
}

//-------------------------------------------------------------------------
// Purpose       : Get instance
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 02/23/03
//-------------------------------------------------------------------------
PartitionEngine& PartitionEngine::instance()
{
  if( instance_ == NULL )
  {
    instance_ = new PartitionEngine();
    assert( instance != NULL );
  }
  
  return *instance_;
}

//-------------------------------------------------------------------------
// Purpose       : Used by SubEntitySet to update ID map.
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 02/23/03
//-------------------------------------------------------------------------
CubitStatus PartitionEngine::add_to_id_map( SubEntitySet* set, int unique_id )
{
  if( uniqueIdMap.find(unique_id) != uniqueIdMap.end() )
    { assert(0); return CUBIT_FAILURE; }
  uniqueIdMap[unique_id] = set;
  return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Purpose       : Used by SubEntitySet to update ID map.
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 02/23/03
//-------------------------------------------------------------------------
CubitStatus PartitionEngine::remove_from_id_map( SubEntitySet* set, int id )
{
  std::map<int,SubEntitySet*>::iterator itor = uniqueIdMap.find(id);
  if( itor == uniqueIdMap.end() || itor->second != set )
    { assert(0); return CUBIT_FAILURE; }
  uniqueIdMap.erase(itor);
  return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Purpose       : Retreive an entity from the ID map
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 02/23/03
//-------------------------------------------------------------------------
SubEntitySet* PartitionEngine::get_from_id_map( int id )
{
  std::map<int,SubEntitySet*>::iterator itor = uniqueIdMap.find(id);
  return itor == uniqueIdMap.end() ? 0 : itor->second;
}

//-------------------------------------------------------------------------
// Purpose       : Replicate solid model topology in partition layer
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 02/23/03
//-------------------------------------------------------------------------
SubCurve* PartitionEngine::replace_curve( Curve* curve )
{
  SubCurve* result = dynamic_cast<SubCurve*>(curve);
  if( result )
    return result;
 
  if( dynamic_cast<SegmentedCurve*>(curve) )
    return 0;
  
  DLIList<TopologyBridge*> list;
  //curve->get_children_virt( list );
  //fix_up_query_results( list );
  curve->get_children(list, false, layer());
  PartitionPoint *start_pt = 0, *end_pt = 0;
  bool new_start = false, new_end = false;
  list.reset();
  TopologyBridge* bridge = list.get_and_step();
  start_pt = dynamic_cast<PartitionPoint*>(bridge);
  if( ! start_pt ) 
  {
    start_pt = replace_point( dynamic_cast<TBPoint*>(bridge) );
    if( ! start_pt )
      return 0;
    new_start = true;
  }
  if( list.size() == 1 )
  {
    end_pt = start_pt;
  }
  else if( list.size() == 2 )
  {
    bridge = list.get();
    end_pt = dynamic_cast<PartitionPoint*>(bridge);
    if( !end_pt )
    {
      end_pt = replace_point( dynamic_cast<TBPoint*>(bridge) );
      new_end = true;
      if( !end_pt )
      {
        if( new_start )
          restore_point( start_pt );
        return 0;
      }
    }
  }
  
  assert( start_pt && end_pt );
  
  result = new SubCurve( curve );
  if( ! result->start_point( start_pt ) ||
      ! result->end_point( end_pt ) )
  {
    result->start_point(0);
    result->end_point(0);
    delete result;
    if( new_start )
      restore_point( start_pt );
    if( new_end )
      restore_point( end_pt );
    return 0;
  }
  if( curve->owner() )
  {
    curve->owner()->swap_bridge(curve, result, false);
  }
  curve->owner(&(result->sub_entity_set()));
  
  list.clean_out();
  curve->get_parents_virt( list );
  //fix_up_query_results( list );
  
  for( int i = list.size(); i--; )
  {
    bridge = list.get_and_step();
    CoEdgeSM* coedge = dynamic_cast<CoEdgeSM*>(bridge);
    assert( coedge != NULL );
    PartitionCoEdge* pcoedge = new PartitionCoEdge( coedge );
    CubitStatus s = result->add( pcoedge );
    assert( s );
    if( coedge->owner() )
      coedge->owner()->swap_bridge(coedge, pcoedge, false);
    coedge->owner(&(pcoedge->sub_entity_set()));
  }
  
  
  PartitionBody* body = 0;
  if( result->start_point()->sub_entity_set().body() )
    body = result->start_point()->sub_entity_set().body();
  else if( result->end_point()->sub_entity_set().body() )
    body = result->end_point()->sub_entity_set().body();
  else
    body = make_body(result);
    
  if( body )
  {
    if( !result->start_point()->sub_entity_set().body() )
      body->add(result->start_point()->sub_entity_set());
    if( !result->end_point()->sub_entity_set().body() )
      body->add(result->end_point()->sub_entity_set());
    body->add(result->sub_entity_set());
    PartitionCoEdge* coedge = 0;
    while ( (coedge = result->next_coedge(coedge)) )
      if ( !coedge->sub_entity_set().body() )
        body->add(coedge->sub_entity_set());
  }

  // we must notify the graphics of the modify from "real" to virtual -- KGM
  // I realize that this is not where the other notifies are completed but there
  // is no knowledge of the change later on. 
  CubitObservable* observer = dynamic_cast<CubitObservable*>(result->topology_entity());
  if (observer)
  {
    observer->notify_all_observers(GEOMETRY_MODIFIED);
  }

  return result;
}

//-------------------------------------------------------------------------
// Purpose       : Replicate solid model topology in partition layer
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 02/23/03
//-------------------------------------------------------------------------
SubSurface* PartitionEngine::replace_surface( Surface* surface )
{
  int i, j;

  // check if the input surface is already a SubSurface
  SubSurface* result = dynamic_cast<SubSurface*>(surface);
  if( result )
    return result;
 
  // check if the input surface is already a PartitionSurface
  if( dynamic_cast<PartitionSurface*>(surface) )
    return 0;
  
  // create a new SubSurface, swap out the TBOwner if it
  // exists and reset the input surface TBOwner.
  DLIList<TopologyBridge*> loop_list, coe_list, curve_list;
  surface->get_children( loop_list, false, layer() );
  result = new SubSurface( surface );
  if( surface->owner() )
    surface->owner()->swap_bridge(surface, result, false);
  surface->owner(&(result->sub_entity_set()));
  
  // Add the loop onto the new surface
  loop_list.reset();
  for( i = loop_list.size(); i--; )
  {
    LoopSM* loop = dynamic_cast<LoopSM*>(loop_list.get_and_step());
    PartitionLoop* partloop = new PartitionLoop();
    result->add( partloop );
    if( loop->owner() )
      loop->owner()->swap_bridge( loop, partloop, false );
    loop->owner(0);
    
    // Take care of the curves in the loop
    coe_list.clean_out();
    loop->get_children( coe_list, false, layer() );
    
    coe_list.reset();
    for( j = coe_list.size(); j--; )
    {
      CoEdgeSM* coedge = dynamic_cast<CoEdgeSM*>(coe_list.get_and_step());
      if(!dynamic_cast<PartitionCoEdge*>(coedge) &&
         !dynamic_cast<SubEntitySet*>(coedge->owner()) )
      {
        curve_list.clean_out();
        coedge->get_children_virt( curve_list );
        assert( curve_list.size() == 1 );
      
        Curve* curve = dynamic_cast<Curve*>(curve_list.get());
        assert(! dynamic_cast<PartitionCurve*>(curve) );
        SubCurve* pcurve = replace_curve( curve );
        assert(pcurve != NULL);
      }
    }
    
    coe_list.clean_out();
    loop->get_children( coe_list, false, layer() );
    
    coe_list.reset();
    PartitionCoEdge *curr, *prev = 0;
    coe_list.reset();
    for( j = coe_list.size(); j--; )
    {
      curr = dynamic_cast<PartitionCoEdge*>(coe_list.get_and_step());
      assert( curr != NULL );
      partloop->insert_after( curr, prev );
      prev = curr;
    }
    
  }
  
  // If the facet data is incorrect put everything back
  if (!result->init_facet_data()) {
    restore_surface( result );
    return 0;
  }
  
  // finalize by making the body
  PartitionBody* body = make_body(result);
  if( body )
    body->add(result->sub_entity_set());
  
  return result;
}

//-------------------------------------------------------------------------
// Purpose       : Replicate solid model topology in partition layer
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 02/23/03
//-------------------------------------------------------------------------
PartitionLump* PartitionEngine::replace_lump( Lump* lump )
{
  int i, j;
  PartitionLump* result = dynamic_cast<PartitionLump*>(lump);
  if( result )
    return result;
 
  DLIList<TopologyBridge*> shell_list, surface_list;
  lump->get_children( shell_list, false, layer() );
  
  DLIList<PartitionShell*> new_shells;
  shell_list.reset();
  for( i = shell_list.size(); i--; )
  {
    ShellSM* shell = dynamic_cast<ShellSM*>(shell_list.get_and_step());
    PartitionShell* pshell = new PartitionShell();
    new_shells.append(pshell);
    if( shell->owner() )
      shell->owner()->swap_bridge( shell, pshell, false );
    shell->owner(0);
    
    surface_list.clean_out();
    shell->get_children( surface_list, false, layer() );
    
    surface_list.reset();
    for( j = surface_list.size(); j--; )
    {
      Surface* surface = dynamic_cast<Surface*>(surface_list.get_and_step());
      PartitionSurface* psurf = dynamic_cast<PartitionSurface*>(surface);
      if( !psurf )
      {
        psurf = replace_surface( surface );
        if (!psurf) 
        {
          while (new_shells.size())
            destroy_shell(new_shells.pop());
          return 0;
        }
      }
      
      Surface* real_surf = dynamic_cast<Surface*>(psurf->partitioned_entity());
      assert(!!real_surf);
      CubitSense sense = real_surf->get_shell_sense(shell);
      assert(sense == CUBIT_FORWARD || sense == CUBIT_REVERSED);
      pshell->add( psurf, sense );
    }
  }

  result = new PartitionLump( lump );
  if( lump->owner() )
    lump->owner()->swap_bridge(lump, result, false);
  lump->owner(&(result->sub_entity_set()));
  new_shells.reverse();
  while( new_shells.size() )
    result->add(new_shells.pop());
  
  PartitionBody* body = make_body(result);
  if( body )
    body->add( result->sub_entity_set() );
  
  return result;
}

//-------------------------------------------------------------------------
// Purpose       : Replicate solid model topology in partition layer
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 02/23/03
//-------------------------------------------------------------------------
PartitionPoint* PartitionEngine::replace_point( TBPoint* point )
{
  PartitionPoint* result = dynamic_cast<PartitionPoint*>(point);
  if( !result )
  {
    result = new PartitionPoint( point );
  }  
  
  if( result )
  {
    if( point->owner() )
      point->owner()->swap_bridge( point, result, false );
    point->owner( &(result->sub_entity_set()) );
  }
  
  return result;
}

//-------------------------------------------------------------------------
// Purpose       : Restore original solid model topology.  (Inverse of
//                 replace_point.)
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 02/23/03
//-------------------------------------------------------------------------
TBPoint* PartitionEngine::restore_point( PartitionPoint* point )
{
  if( point->num_curves() )
    return 0;
  
  TBPoint* result = point->real_point();
  if( result )
  {
    point->sub_entity_set().unwrap_attributes();
    point->sub_entity_set().remove_bridge(result);
    if( point->owner() )
      point->owner()->swap_bridge( point, result, false );
    delete point;
  }
  
  return result;
}

//-------------------------------------------------------------------------
// Purpose       : Restore original solid model topology.  (Inverse of
//                 replace_curve.)
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 02/23/03
//-------------------------------------------------------------------------
Curve* PartitionEngine::restore_curve( SubCurve* curve )
{
  if( curve->num_partitions() > 1 )
    return 0;
  
  PartitionCoEdge* coedge = 0;
  while( (coedge = curve->next_coedge( coedge )) )
    if( coedge->get_loop() )
      return 0;
  
  while( (coedge = curve->next_coedge(0)) )
  {
    curve->remove( coedge );
    CoEdgeSM* real_coedge = coedge->real_coedge();
    if( real_coedge )
    {
      coedge->sub_entity_set().remove_bridge( real_coedge );
      if( coedge->owner() )
        coedge->owner()->swap_bridge( coedge, real_coedge, false );
    }
    delete coedge;
  }
  
  PartitionPoint* start = curve->start_point();
  curve->start_point( 0 );
  if( start->num_curves() == 0 )
  {
    if( start->real_point() )
      restore_point( start );
    else
      delete start;
  }
  
  PartitionPoint* end = curve->end_point();
  curve->end_point(0);
  if( end->num_curves() == 0 )
  {
    if( end->real_point() )
      restore_point( end );
    else
      delete end;
  }
  
  Curve* result = curve->real_curve();
  curve->sub_entity_set().unwrap_attributes();
  curve->sub_entity_set().remove_bridge( result );
  if( curve->owner() )
    curve->owner()->swap_bridge( curve, result, false );
  delete curve;

  // we must notify the graphics of the modify from "real" to virtual -- KGM
  // I realize that this is not where the other notifies are completed but there
  // is no knowledge of the change later on. 
  CubitObservable* observer = dynamic_cast<CubitObservable*>(result->topology_entity());
  if (observer)
  {
    observer->notify_all_observers(GEOMETRY_MODIFIED);
  }

  return result;
}


//-------------------------------------------------------------------------
// Purpose       : Restore original solid model topology.  (Inverse of
//                 replace_surface.)
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 02/23/03
//-------------------------------------------------------------------------
Surface* PartitionEngine::restore_surface( SubSurface* surface )
{
  if( surface->sub_entity_set().has_lower_order() || surface->next_co_surface() )
    return 0;
  
  while( PartitionLoop* loop = surface->next_loop() )
    // while surface has loops...
  {
    surface->remove( loop );
    
    while( loop->first_coedge() )
    {
      PartitionCoEdge* coedge = loop->first_coedge();
      loop->remove( coedge );
      
      SubCurve* curve = dynamic_cast<SubCurve*>(coedge->get_curve());
      if( !curve || curve->sub_entity_set().has_lower_order() )
        continue;
      
      bool remove = true;
      coedge = 0;
      while( (coedge = curve->next_coedge( coedge )) )
      {
        if( coedge->get_loop() )
        {
          remove = false;
          break;
        }
      }
      
      if( remove )
      {
        restore_curve( curve );
      }
      
    }
    
    delete loop;
  }
  
  
  Surface* result = surface->partitioned_surface();
  surface->sub_entity_set().unwrap_attributes();
  surface->sub_entity_set().remove_bridge( result );
  if( surface->owner() )
    surface->owner()->swap_bridge( surface, result, false );
  delete surface;
  return result;
}

//-------------------------------------------------------------------------
// Purpose       : Destroy a single partition of a volume.
//
// Special Notes : Will restore original volume if this is the last/only
//                  partition.
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 02/11/03
//-------------------------------------------------------------------------
CubitStatus PartitionEngine::destroy_lump( PartitionLump* lump )
{
    /* only partition lump? */
  if ( ! lump->sub_entity_set().has_multiple_sub_entities() )
  {
    DLIList<PartitionEntity*> child_entities;
    lump->sub_entity_set().get_lower_order( child_entities );
    DLIList<PartitionSurface*> surf_list;
    CAST_LIST( child_entities, surf_list, PartitionSurface );
    
    while( surf_list.size() )
      if( ! remove_surface( surf_list.pop() ) )
        return CUBIT_FAILURE;
    
    return restore_lump(lump) ? CUBIT_SUCCESS : CUBIT_FAILURE;
  }
  
  PartitionShell* shell;
  while ( (shell = lump->next_shell(0)) )
  {
    lump->remove(shell);
    destroy_shell(shell);
  }
  delete lump;
  
  return CUBIT_SUCCESS;
}
      
    

//-------------------------------------------------------------------------
// Purpose       : Restore original solid model topology.  (Inverse of
//                 replace_lump.)
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 02/23/03
//-------------------------------------------------------------------------
Lump* PartitionEngine::restore_lump( PartitionLump* lump )
{
  if( lump->sub_entity_set().has_lower_order() )
    return 0;
  
  while( PartitionShell* shell = lump->next_shell() )
    // while surface has loops...
  {
    lump->remove( shell );
    destroy_shell( shell );
  }
  
  
  Lump* result = lump->real_lump();
  lump->sub_entity_set().unwrap_attributes();
  lump->sub_entity_set().remove_bridge( result );
  if( lump->owner() )
    lump->owner()->swap_bridge( lump, result, false );
  delete lump;
  return result;
}

bool PartitionEngine::is_partition(TBOwner *bridge_owner)
{
  bool ret = false;
  if(bridge_owner)
  {
    if(dynamic_cast<SubEntitySet*>(bridge_owner))
    {
      ret = true;
    }
  }
  return ret;
}

bool PartitionEngine::is_composite(TBOwner *bridge_owner)
{
  return false;
}

bool PartitionEngine::is_composite(TopologyBridge *bridge)
{
  return false;
}

//-------------------------------------------------------------------------
// Purpose       : Destroy a shell
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 03/11/04
//-------------------------------------------------------------------------
CubitStatus PartitionEngine::destroy_shell( PartitionShell* shell )
{
  if (shell->get_lump())
  {
    assert(0);
    return CUBIT_FAILURE;
  }
  
  while( PartitionCoSurf* cosurf = shell->next_co_surface() )
  {
    shell->remove( cosurf );

    SubSurface* psurf = dynamic_cast<SubSurface*>(cosurf->get_surface());
    if (psurf)
    {
      psurf->remove(cosurf);

        // If PartitionSurface is the same as the underlying surface
        // and is not part of any other PartitionLumps, remove it.
      if( !psurf->sub_entity_set().has_lower_order() &&
           psurf->next_co_surface() == 0 )
        restore_surface( psurf );
    }
    
    delete cosurf;
  }
    
  delete shell;
  return CUBIT_SUCCESS;
}



//-------------------------------------------------------------------------
// Purpose       : Partition a curve
//
// Special Notes : public interface to curve partitioning
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 02/23/03
//-------------------------------------------------------------------------
TBPoint* PartitionEngine::insert_point( Curve* curve, double u )
{
  PartitionCurve* pcurve = dynamic_cast<PartitionCurve*>(curve);
  bool replaced_curve = false;
  
  if( !pcurve )
  {
    if( CompositeCurve* comp = dynamic_cast<CompositeCurve*>(curve) )
    {
      return CompositeEngine::instance().insert_point( comp, u );
    }
  
    pcurve = replace_curve( curve );
    if( !pcurve )
      return 0;
    replaced_curve = true;
  }
  
  CubitVector coords;
  if( ! pcurve->position_from_u( u, coords ) )
  {
    if( replaced_curve )
      restore_curve( dynamic_cast<SubCurve*>(pcurve) );
    return 0;
  }
  
  PartitionPoint* npoint = new PartitionPoint( coords, pcurve );
  PartitionCurve* ncurve = insert_point(pcurve, npoint);
  
  if( !ncurve )
  {
    if( replaced_curve )
      restore_curve( dynamic_cast<SubCurve*>(pcurve) );
    delete npoint;
    return 0;
  }
  
  return npoint;
}


//-------------------------------------------------------------------------
// Purpose       : Partition a curve
//
// Special Notes : not public
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 02/23/03
//-------------------------------------------------------------------------
PartitionCurve* PartitionEngine::insert_point( PartitionCurve* pcurve, 
                                               PartitionPoint* npoint )
{
  const double TOLSQR = GEOMETRY_RESABS * GEOMETRY_RESABS;
  CubitVector coords = npoint->coordinates();
  double u = pcurve->u_from_position(coords);
  double u_min, u_max;
  pcurve->get_param_range( u_min, u_max );

  if(u_min > u_max)
  {
    double dtemp = u_min;
    u_min = u_max;
    u_max = dtemp;
  }

  if( u - u_min < CUBIT_RESABS || u_max - u < CUBIT_RESABS )
    return 0;
  
  if( (pcurve->start_point()->coordinates() - coords).length_squared() < TOLSQR 
   || (pcurve->end_point()->coordinates() - coords).length_squared() < TOLSQR )
    return 0;
  
  PartitionCurve* ncurve = pcurve->split( u );
  if( !ncurve ) 
    return 0;
  
  ncurve->end_point( pcurve->end_point() );
  pcurve->end_point( npoint );
  ncurve->start_point( npoint );
  pcurve->fix_facet_data( ncurve );
  
  PartitionCoEdge *pcoedge = 0, *ncoedge = 0;
  while( (pcoedge = pcurve->next_coedge( pcoedge )) )
  {
    ncoedge = new PartitionCoEdge( pcoedge );
    ncurve->add( ncoedge );
    if( pcoedge->get_loop() )
    {
      if( pcoedge->sense() == CUBIT_FORWARD )
        pcoedge->get_loop()->insert_after( ncoedge, pcoedge );
      else
        pcoedge->get_loop()->insert_before( ncoedge, pcoedge );
    }
  }
  
  if( pcurve->owner() )
    pcurve->owner()->notify_split( ncurve, pcurve );
  
  return ncurve;
}

//-------------------------------------------------------------------------
// Purpose       : un-partition a curve
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 02/23/03
//-------------------------------------------------------------------------
Curve* PartitionEngine::remove_point( PartitionPoint* point, 
                                      PartitionCurve* dead_curves[2] )
{
  if( dead_curves )
    dead_curves[0] = dead_curves[1] = 0;
  
  if( point->num_curves() > 2 )
    return 0;
      
  DLIList<TopologyBridge*> curve_bridges;
  point->get_parents_virt( curve_bridges );
  if( curve_bridges.size() > 2 )
    return 0;
  
  if( curve_bridges.size() == 1 )
  {
    PartPTCurve* point_curve = dynamic_cast<PartPTCurve*>(curve_bridges.get());
    if( ! point_curve )
      return 0;
    
    if( dead_curves )
      dead_curves[0] = point_curve;
    remove_point_curve( point );
    return 0;
  }

/*  
  TopologyBridge* real = point->partitioned_entity();
  if( dynamic_cast<TBPoint*>(real) )
  {
    Curve* curve1 = dynamic_cast<Curve*>(curve_bridges.get());
    Curve* curve2 = dynamic_cast<Curve*>(curve_bridges.next());
    PartitionCurve* pc1 = dynamic_cast<PartitionCurve*>(curve1);
    PartitionCurve* pc2 = dynamic_cast<PartitionCurve*>(curve2);
    assert( pc1 || pc2 );
    if( !pc1 )
      pc1 = replace_curve( curve1 );
    if( !pc2 )
      pc2 = replace_curve( curve2 );
    
    CompositeEngine::instance()->remove_point( dynamic_cast<TBPoint*>(real) );
    bool reverse = ( pc1->start_point() == pc2->start_point() ||
                     pc2->end_point() == pc2->end_point() );
    if( pc1->partitioned_entity() == 0 )
    {
      bool prepend = false;
      if( pc1->other_point( pc2->start_point() ) )
        prepend = true;
      pc2->sub_entity_set().merge( pc1->sub_entity_set(), reverse, prepend );
    }
    else if( pc2->partitioned_entity() == 0 )
    {
      bool prepend = false;
      if( pc2->other_point( pc1->start_point() ) )
        prepend = true;
      pc1->sub_entity_set().merge( pc2->sub_entity_set(), reverse, prepend );
    }
    else
      return 0;
  }
*/  
  if( point->num_curves() != 2 )
    return 0;
  
  PartitionCurve* survivor = point->next_curve(0);
  PartitionCurve* casualty = point->next_curve(survivor);
  if( survivor->start_point() == point &&
      casualty->end_point() == point )
  {
      // swap
    survivor = casualty;
    casualty = point->next_curve(0);
  }
  
  if( ! survivor->combine( casualty ) )
    return 0;
    
  PartitionPoint* other_pt = 
    casualty->start_point() == point ? 
    casualty->end_point() : casualty->start_point();
  casualty->start_point(0);
  casualty->end_point(0);
  
  if( survivor->start_point() == point )
    survivor->start_point( other_pt );
  else
    survivor->end_point( other_pt );
  
  delete point;
  
  while( PartitionCoEdge* coedge = casualty->next_coedge(0) )
  {
    if( coedge->get_loop() )
      coedge->get_loop()->remove( coedge );
    casualty->remove( coedge );
    delete coedge;
  }
  delete casualty;
  if( dead_curves )
    dead_curves[0] = casualty;
  
  Curve* result = survivor;
  SubCurve* subcurve = dynamic_cast<SubCurve*>(survivor);
  if( subcurve && !survivor->sub_entity_set().has_multiple_sub_entities() )
  {
    Curve* real_curve = restore_curve( subcurve );
    if (real_curve)
    { 
      result = real_curve;
      if (dead_curves)
        dead_curves[1] = survivor;
    } 
  }

  return result;
}

//-------------------------------------------------------------------------
// Purpose       : Remove a point-curve (point imprinted on a surface)
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 02/23/03
//-------------------------------------------------------------------------
CubitStatus PartitionEngine::remove_point_curve( PartitionPoint* point )
{
  if( point->num_curves() != 1 )
    return CUBIT_FAILURE;
  
  PartPTCurve* curve = dynamic_cast<PartPTCurve*>(point->next_curve(0));
  if( !curve )
    return CUBIT_FAILURE;
  
  PartitionSurface* surf = 0;  
  if( curve->num_coedges() )
  {
    assert( curve->num_coedges() == 1 );
    PartitionCoEdge* coedge = curve->next_coedge(0);
    PartitionLoop* loop = coedge->get_loop();
    if( loop )
    {
      surf = loop->get_surface();
      surf->remove( loop );
      loop->remove( coedge );
      assert( loop->num_coedges() == 0 );
      delete loop;
    }
    curve->remove( coedge );
    delete coedge;
  }
  
  curve->start_point(0);
  curve->end_point(0);
  delete curve;
  delete point;
  
  if( surf && surf->owner() )
    surf->owner()->notify_topology_modified( surf );
  
  return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Purpose       : Topological modifications for surface partitioning
//
// Special Notes : non-public, assumes surface facets have already been
//                 imprinted and the polyline of the imprint is already
//                 associated with the passed curve.
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 02/23/03
//-------------------------------------------------------------------------
PartitionSurface* PartitionEngine::insert_curve( PartitionSurface* surface, 
                                                 SegmentedCurve* curve )
{
  PartitionSurface* new_surface = 0;
  assert( &(curve->sub_entity_set()) == &(surface->sub_entity_set()) );
  
  PartitionPoint *start_point = curve->start_point();
  PartitionPoint *end_point   = curve->end_point();
  PartitionLoop  *start_loop = 0, *end_loop = 0;
  PartitionCoEdge *start_prev_coedge = 0, *start_next_coedge = 0;
  PartitionCoEdge *end_prev_coedge = 0, *end_next_coedge = 0;
  
    // find where to insert new curve in loop
  if ( !find_coedges( surface, curve, start_point, start_prev_coedge, start_next_coedge ) )
    return 0;
  if ( !find_coedges( surface, curve, end_point, end_prev_coedge, end_next_coedge )) 
    return 0;
    
  assert((!start_prev_coedge && !start_next_coedge) ||
         (start_prev_coedge && start_next_coedge && 
          start_prev_coedge->get_loop() == start_next_coedge->get_loop()));
  assert((!end_prev_coedge && !end_next_coedge) ||
         (end_prev_coedge && end_next_coedge && 
          end_prev_coedge->get_loop() == end_next_coedge->get_loop()));
  if ( start_prev_coedge )
    start_loop = start_prev_coedge->get_loop();
  if ( end_prev_coedge )
    end_loop = end_prev_coedge->get_loop();
  
  
    // remove any point-curves
  PartitionCurve* pt_curve;
  PartPTCurve* point_curve;
  for ( int i = 0; i < 2; i++ ) // do once for each end point
  {
    PartitionPoint* pt = i ? end_point : start_point;
    pt_curve = pt->next_curve(0);
    while ( pt_curve  )
    {
      point_curve = dynamic_cast<PartPTCurve*>(pt_curve);
      pt_curve = pt->next_curve(pt_curve);
      if (point_curve && point_curve->next_coedge(0)->get_loop()->get_surface() == surface)
      {
        point_curve->start_point(0);
        point_curve->end_point(0);

        PartitionCoEdge* coedge = point_curve->next_coedge(0);
        point_curve->remove(coedge);
        assert( !point_curve->next_coedge(0) );
        PartitionLoop* loop = coedge->get_loop();
        loop->remove(coedge);
        assert( !loop->first_coedge() );
        surface->remove(loop);
        delete loop;
        delete coedge;
        delete point_curve;
      }
    }
  }
  
  
  PartitionCoEdge* forward = new PartitionCoEdge( surface, CUBIT_FORWARD );
  PartitionCoEdge* reverse = new PartitionCoEdge( surface, CUBIT_REVERSED );
  curve->add(forward);
  curve->add(reverse);

    // new one-curve loop?
  if( !start_loop && !end_loop )
  {
    PartitionLoop* new_loop = new PartitionLoop();
    surface->add( new_loop );
    new_loop->insert_after( forward, 0 );
    
    if( curve->start_point() == curve->end_point() )
    {
      PartitionLoop* loop2 = new PartitionLoop();
      loop2->insert_after( reverse, 0 );
      
      if( VGLoopTool
          <PartitionSurface,PartitionLoop,PartitionCoEdge,PartitionCurve,PartitionPoint>
          ::loop_angle_metric( forward ) > 0 )
      {
        new_loop = loop2;
      }
      else
      {
        surface->remove(new_loop);
        surface->add( loop2 );
      }
      
      new_surface = split_surface( surface, new_loop->first_coedge() );
      new_surface->add( new_loop );
    }
    else
    {
      new_loop->insert_after( reverse, forward );
    }
  }
      
    // sipe
  else if( !start_loop || !end_loop )
  {
    PartitionCoEdge* prev;
    PartitionLoop* loop;
    if( start_loop )
    {
      loop = start_loop;
      prev = start_prev_coedge;
    }
    else
    {
      loop = end_loop;
      prev = end_prev_coedge;
    }
    
    if( forward->start_point() == prev->end_point() )
    {
      loop->insert_after( forward, prev );
      loop->insert_after( reverse, forward );
    }
    else
    {
      assert( reverse->start_point() == prev->end_point() );
      loop->insert_after( reverse, prev );
      loop->insert_after( forward, reverse );
    }
  }
  
    // join loops
  else if( start_loop != end_loop )
  {
    PartitionCoEdge* prev = start_prev_coedge;
    PartitionCoEdge* coedge = end_next_coedge;
    PartitionCoEdge* next = 0, *other_coedge;
    
    if( forward->start_point() == prev->end_point() )
    {
      start_loop->insert_after( forward, prev );
      prev = forward;
      other_coedge = reverse;
    }
    else 
    {
      assert( reverse->start_point() == prev->end_point() );
      start_loop->insert_after( reverse, prev );
      prev = reverse;
      other_coedge = forward;
    }
    
    while ( end_loop->first_coedge() )
    {
      next = end_loop->next_coedge( coedge );
      end_loop->remove( coedge );
      start_loop->insert_after( coedge, prev );
      prev = coedge;
      coedge = next;
    }
    
      // The other coedge for the curve we are restoring...
    start_loop->insert_after( other_coedge, prev );
    
    assert( end_loop->num_coedges() == 0 );
    surface->remove( end_loop );
    delete end_loop;
  }
  
    // split a loop (and create a new surface)
  else
  {
    assert( start_loop == end_loop );
    
      // If the curve we are adding results in a hole that
      // intersects the loop we are splitting at a single
      // vertex, we need to figure out which of the coedges
      // goes in the hole and which is added to the loop
      // we are splitting.
    if( forward->start_point() == forward->end_point() )
    {
      assert( start_next_coedge == end_next_coedge );
      assert( start_prev_coedge == end_prev_coedge );
      
      CubitVector prev_tan, forward_tan, reverse_tan, normal, junk;
      CubitVector point = forward->start_point()->coordinates();
      start_prev_coedge->get_curve()->closest_point( point, junk, &prev_tan );
      if( start_prev_coedge->sense() == CUBIT_FORWARD ) // yes, forward!!!
        prev_tan *= -1.0;
      forward->get_curve()->closest_point( point, junk, &forward_tan );
      reverse_tan = forward_tan * -1.0;
      surface->closest_point( point, 0, &normal );
      
      double angle1 = normal.vector_angle( prev_tan, forward_tan );
      double angle2 = normal.vector_angle( prev_tan, reverse_tan );
      
      if( angle2 < angle1 )
      {
        PartitionCoEdge* temp = reverse;
        reverse = forward;
        forward = temp;
      }
    }
    
    else
    {
      if( forward->end_point() == end_next_coedge->start_point() )
      {
        PartitionCoEdge* temp = reverse;
        reverse = forward;
        forward = temp;
      }
    }
    
    end_loop = new PartitionLoop();
    start_loop->insert_after( forward, end_prev_coedge );
    end_loop->insert_after( reverse, 0 );
    
    PartitionCoEdge* coedge = end_next_coedge;
    PartitionCoEdge* prev = reverse;
    while( coedge != start_next_coedge )
    {
      PartitionCoEdge* next = start_loop->next_coedge( coedge );
      start_loop->remove( coedge );
      end_loop->insert_after( coedge, prev );
      prev = coedge;
      coedge = next;
    }
    
    new_surface = split_surface( surface, reverse );
    new_surface->add( end_loop );
  }
  
  if( new_surface )
  {
    PartitionCoSurf* cos = 0;
    while( (cos = surface->next_co_surface( cos )) )
      cos->get_shell()->add( new_surface, cos->sense() );
    if( surface->owner() )
      surface->owner()->notify_split( new_surface, surface );
    return new_surface;
  }
  else
  {
    if( surface->owner() )
      surface->owner()->notify_topology_modified( surface );
    return surface;
  }
}


//-------------------------------------------------------------------------
// Purpose       : Find where to insert a curve in a loop
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 03/24/03
//-------------------------------------------------------------------------
CubitStatus PartitionEngine::find_coedges( PartitionSurface* surface,
                                           PartitionCurve* curve,
                                           PartitionPoint* point,
                                           PartitionCoEdge*& previous,
                                           PartitionCoEdge*& next )
{
  const char* const bad_loop_message = "Internal error: Invalid loop. (%s:%d)\n";
  
    // Find list of all coedges around passed point
    // and in passed surface.
  DLIList<PartitionCoEdge*> point_coedges;
  PartitionCoEdge* coedge;
  PartitionCurve* pt_curve = 0;
  while ( (pt_curve = point->next_curve(pt_curve)) )
  {
    coedge = 0;
    while ( (coedge = pt_curve->next_coedge(coedge)) )
    {
      if (coedge->get_loop() &&
          coedge->get_loop()->get_surface() == surface &&
          !dynamic_cast<PartPTCurve*>(coedge->get_curve()))
        point_coedges.append(coedge);
    }
  }
  
    // TBPoint is at end of a sipe/hardline
  if ( point_coedges.size() == 0 )
  {
    previous = next = 0;
    return CUBIT_SUCCESS;
  }
  
    // One coedge - closed curve
  if ( point_coedges.size() == 1 &&
       point_coedges.get()->start_point() ==
       point_coedges.get()->end_point() )
  {
    previous = next = point_coedges.get();
    return CUBIT_SUCCESS;
  }
  
    // Broken loop?
  if ( point_coedges.size() % 2 != 0 )
  {
    PRINT_ERROR(bad_loop_message,__FILE__,__LINE__);
    assert(0);
    return CUBIT_FAILURE;
  }
  
    // If only two coedges, then we're done
    // If closed curve, the answer may be the same coedge for both 
    // prev and next.  E.g. A torus cracked along the outside major
    // radius and a sipe that intersects that crack curve.  Need
    // to fall though to the more complex check below in that case.
  if ( point_coedges.size() == 2 &&
       point_coedges.get()->get_curve()->start_point() != 
         point_coedges.get()->get_curve()->end_point() &&
       point_coedges.next()->get_curve()->start_point() !=
         point_coedges.next()->get_curve()->end_point() )
  {
    previous = point_coedges.get();
    next = point_coedges.next();
    
    if ( previous->start_point() == point )
      std::swap(previous, next);
    
    return CUBIT_SUCCESS;
  }
  
    // Find previous/next coedges by using order of
    // facet edges about point in surface facetting.
  CubitPointData* facet_point = point->facet_point();
  DLIList<CubitFacetEdgeData*> curve_edges;
  curve->get_facet_data(curve_edges);
  
  CubitFacetEdgeData* edge;
  //bool edge_reversed;
  if ( point == curve->start_point() )
  {
    curve_edges.reset();
    edge = curve_edges.get();
  }
  else
  {
    assert(curve->end_point() == point);
    curve_edges.last();
    edge = curve_edges.get();
  }
    
  PartitionCurve *prev_curve, *next_curve;
  prev_curve = next_curve_around_point( surface, edge, facet_point, true );
  next_curve = next_curve_around_point( surface, edge, facet_point, false );
  
  if ( !prev_curve || !next_curve )
  {
    PRINT_ERROR(bad_loop_message,__FILE__,__LINE__);
    assert(0);
    return CUBIT_FAILURE;
  }
  
  coedge = 0;
  while ( (coedge = prev_curve->next_coedge(coedge)) )
  {
    if ( coedge->get_loop() &&  // parent not a partition surface?
         coedge->get_loop()->get_surface() == surface &&
         coedge->end_point() == point )
    {
      previous = coedge;
      break;
    }
  }
  coedge = 0;
  while ( (coedge = next_curve->next_coedge(coedge)) )
  {
    if ( coedge->get_loop() &&
         coedge->get_loop()->get_surface() == surface &&
         coedge->start_point() == point )
    {
      next = coedge;
      break;
    }
  }

  return CUBIT_SUCCESS;
}
    
PartitionCurve* PartitionEngine::next_curve_around_point(
                               PartitionSurface *const surface, 
                               CubitFacetEdgeData *const start_edge,
                               CubitPointData *const point, 
                               const bool backwards )
{
  //PartitionCoEdge* result = 0;
  CubitFacetEdge* edge = start_edge;
  
  while ( edge->num_adj_facets() == 2 )
  {
    CubitFacet* facet1 = edge->adj_facet(0);
    CubitFacet* facet2 = edge->adj_facet(1);
    assert (TDVGFacetOwner::get(facet1) == surface);
    assert (TDVGFacetOwner::get(facet2) == surface);
    
    bool edge_reversed = edge->point(1) == point;
    int facet_sense = edge_reversed == backwards ? -1 : 1;
    int facet1_index = facet1->edge_index(edge);
    int facet2_index = facet2->edge_index(edge);
    assert(facet1->edge_use(facet1_index) == -facet2->edge_use(facet2_index));
    
    if ( facet_sense == facet1->edge_use(facet1_index) )
      edge = facet1->next_edge_at_point( edge, point );
    else
      edge = facet2->next_edge_at_point( edge, point );
      
    assert(edge != start_edge);
    
    PartitionEntity* edge_owner = TDVGFacetOwner::get(edge);
    PartitionCurve* curve = dynamic_cast<PartitionCurve*>(edge_owner);
    if ( curve )
      return curve;
  }
  
  return 0;
}



//-------------------------------------------------------------------------
// Purpose       : Un-partition a surface
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 02/23/03
//-------------------------------------------------------------------------
Surface* PartitionEngine::remove_curve( PartitionCurve* passed_curve,
                                        PartitionSurface* dead_surfs[2] )
{
  if( dead_surfs )
    dead_surfs[0] = dead_surfs[1] = 0;
  
  PartPTCurve* point_curve = dynamic_cast<PartPTCurve*>(passed_curve);
  if( point_curve )
  {
    PartitionSurface* surf = 
      point_curve->next_coedge(0)->get_loop()->get_surface();
    if( !remove_point_curve( point_curve->start_point() ) )
      return 0;
    
    SubSurface* subsurf = dynamic_cast<SubSurface*>(surf);
    if( subsurf && !subsurf->next_co_surface() &&
        !subsurf->sub_entity_set().has_lower_order() )
    {
      Surface* result = restore_surface( subsurf );
      if( result && dead_surfs )
        dead_surfs[0] = subsurf;
      return result;
    }
    else
      return surf;
  }
  
  SegmentedCurve* curve = dynamic_cast<SegmentedCurve*>(passed_curve);

  if( !curve )
  {
    assert(0);
  }
/*
    assert( dynamic_cast<SubCurve*>(passed_curve) );
    
    PartitionCoEdge* coedge1 = passed_curve->next_coedge(0);
    PartitionCoEdge* coedge2 = passed_curve->next_coedge(coedge1);
    if( !coedge2 || passed_curve->next_coedge(coedge2) )
      return 0;
      
    DLIList<TopologyBridge*> surfaces;
    coedge1->find_parent_loop()->get_parents_virt( surfaces );
    fix_up_query_results( surfaces );
    assert( surfaces.size() == 1 );
    Surface* surface1 = dynamic_cast<Surface*>(surfaces.get());
    surfaces.clean_out();
    coedge2->find_parent_loop()->get_parents_virt( surfaces );
    fix_up_query_results( surfaces );
    assert( surfaces.size() == 1 );
    Surface* surface2 = dynamic_cast<Surface*>(surfaces.get());
    
    SubSurface* surf1 = dynamic_cast<SubSurface*>(surface1);
    SubSurface* surf2 = dynamic_cast<SubSurface*>(surface2);
    if( surface1 == surface2 )
    {
      if( !surf1 )
        surf1 = surf2 = replace_surface( surface1 );
    }
    else
    {
      if( !surf1 )
        surf1 = replace_surface( surface1 );
      if( !surf2 )
        surf2 = replace_surface( surface2 );
    }
    
      // create SegmentedCurves from SubCurves
    
    
    CompositeEngine::instance()->
      remove_curve( dynamic_cast<Curve*>(passed_curve->partitioned_entity()) );
  }
*/  
  DLIList<TopologyBridge*> tb_list;
  
    // get two coedges to remove
  curve->get_parents_virt( tb_list );
  if( tb_list.size() != 2 )
    return 0;
  PartitionCoEdge* coedge1 = dynamic_cast<PartitionCoEdge*>(tb_list.get_and_step());
  PartitionCoEdge* coedge2 = dynamic_cast<PartitionCoEdge*>(tb_list.get_and_step());
  assert( coedge1 && coedge2 );
  
    // get two loops and surfaces
  PartitionLoop* loop1 = coedge1->get_loop();
  PartitionLoop* loop2 = coedge2->get_loop();
  PartitionSurface* surf1 = loop1->get_surface();
  PartitionSurface* surf2 = loop2->get_surface();
  assert( surf1 && surf2 );
  
    // surfaces must be partitions of the same real surface
  if( surf1->partitioned_entity() != surf2->partitioned_entity() )
    return 0;
  
    // surfaces must have same parent shells
  PartitionCoSurf* cos = 0;
  while( (cos = surf1->next_co_surface( cos )) )
    if( ! surf2->find_first( cos->get_shell(), cos->sense() ) )
      return 0;
    
    // same loop : remove a sipe or split a loop
  if( loop1 == loop2 )
  {
    if( loop1->next_coedge( coedge1 ) != coedge2 &&
        loop1->prev_coedge( coedge1 ) != coedge2 )
    {
      loop2 = new PartitionLoop();
      PartitionCoEdge* coedge = loop1->next_coedge( coedge1 );
      PartitionCoEdge* prev = 0;
      while( coedge != coedge2 )
      {
        loop1->remove( coedge );
        loop2->insert_after( coedge, prev );
        prev = coedge;
        coedge = loop1->next_coedge( coedge1 );
      }
      
      surf1->add( loop2 );
    }
    
    loop1->remove( coedge1 );
    loop1->remove( coedge2 );
    if( loop1->first_coedge() == 0 )
    {
      surf1->remove( loop1 );
      delete loop1;
    }
  }
    // stitch loops
  else
  {
      // insert coedges
    while( loop2->next_coedge(coedge2) != coedge2 ) // all but the dead one
    {
      PartitionCoEdge* coedge = loop2->next_coedge( coedge2 );
      loop2->remove( coedge );
      loop1->insert_before( coedge, coedge1 );
    }
    loop1->remove( coedge1 );
    loop2->remove( coedge2 );
    assert( loop2->first_coedge() == 0 );
    surf2->remove( loop2 );
    if( loop1->num_coedges() == 0 )
    {
      surf1->remove( loop1 );
      delete loop1;
    }

    delete loop2;
  }
  
  curve->remove( coedge1 );
  curve->remove( coedge2 );
  delete coedge1;
  delete coedge2;
  PartitionPoint* start_pt = curve->start_point();
  PartitionPoint* end_pt = curve->end_point();
  curve->start_point(0);
  curve->end_point(0);
  delete curve;
  if ( !start_pt->next_curve() )
  {
    if ( start_pt->real_point() )
      restore_point(start_pt);
    else
      delete start_pt;
  }
  if ( end_pt != start_pt && !end_pt->next_curve() )
  {
    if ( end_pt->real_point() )
      restore_point(end_pt);
    else
      delete end_pt;
  }
  
    // need to combine surfaces?
  if( surf1 != surf2 )
  {
    surf1->combine( surf2 );
    while( PartitionLoop* loop = surf2->next_loop() )
    { 
      surf2->remove( loop );
      surf1->add( loop );
    }
  
    while( PartitionCoSurf* cos = surf2->next_co_surface() )
    {
      if( cos->get_shell() )
        cos->get_shell()->remove( cos );
      surf2->remove( cos );
      delete cos;
    }
    
    if( dead_surfs )
      dead_surfs[0] = surf2;
    delete surf2;
  }
  
    // can remove partition surface?
  SubSurface* subsurf = dynamic_cast<SubSurface*>(surf1);
  if( subsurf && !subsurf->next_co_surface() && 
      !subsurf->sub_entity_set().has_lower_order() )
  {
    if( dead_surfs )
    {
      if( dead_surfs[0] )
        dead_surfs[1] = surf1;
      else
        dead_surfs[0] = surf1;
    }
    return restore_surface( subsurf );
  }
  else
    return surf1;
}

//-------------------------------------------------------------------------
// Purpose       : Find the facet that:
//                   - contains the passed edge 
//                   - uses the passed edge in the specified orientation
//                   - is owned by the passed entity 
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 02/23/03
//-------------------------------------------------------------------------
static CubitFacetData* find_facet( CubitFacetEdge* edge, 
                                   bool sense,
                                   PartitionEntity* owner )
{
  DLIList<CubitFacet*> facets(2);
  edge->facets(facets);
  facets.reset();
  for ( int i = facets.size(); i--; ) {
    CubitFacet* facet = facets.get();
    int edge_index = facet->edge_index(edge);
    assert(edge_index >= 0);
    bool forward = edge->point(0) == facet->point((edge_index+1)%3);
    if ( forward != sense || TDVGFacetOwner::get(facet) != owner )
      facets.extract();
    else
      facets.step();
  }
  
  return facets.size() == 1 ? dynamic_cast<CubitFacetData*>(facets.get()) : 0;
}

//-------------------------------------------------------------------------
// Purpose       : Find a facet adjacent to the passed point and owned by
//                 the passed PartitionEntity
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 02/23/03
//-------------------------------------------------------------------------
static CubitFacetData* find_facet( CubitPoint* pt, PartitionEntity* owner )
{
  DLIList<CubitFacet*> facets;
  pt->facets(facets);
  facets.reset();
  for ( int i = facets.size(); i--; )
  {
    if ( TDVGFacetOwner::get(facets.get()) == owner )
      facets.step();
    else
      facets.extract();
  }
  
  return facets.size() == 1 ? dynamic_cast<CubitFacetData*>(facets.get()) : 0;
}

//-------------------------------------------------------------------------
// Purpose       : Find the facet adjacent to the passed facet on the
//                 specified edge, and owned by the passed PartitionEntity
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 02/23/03
//-------------------------------------------------------------------------
static CubitFacetData* other_facet( CubitFacetData* facet,
                                    int edge,
                                    PartitionEntity* owner )
{
  DLIList<CubitFacet*> facets;
  facet->shared_facets( facet->point((edge+1)%3), 
                        facet->point((edge+2)%3),
                        facets );
  for ( int i = facets.size(); i--; )
  {
    CubitFacet* facet = facets.get();
    if( TDVGFacetOwner::get(facet) != owner )
      facets.extract();
    else
      facets.step();
  }
  
  return facets.size() == 1 ? dynamic_cast<CubitFacetData*>(facets.get()) : 0;
}
/*
static void draw_facet_edges( CubitFacet* facet, int color )
{
  float x1, y1, z1, x2, y2, z2;
  for ( int i = 0; i < 3; i++ )
  {
    x1 = (float)facet->point(i)->coordinates().x();
    y1 = (float)facet->point(i)->coordinates().y();
    z1 = (float)facet->point(i)->coordinates().z();
    x2 = (float)facet->point((i+1)%3)->coordinates().x();
    y2 = (float)facet->point((i+1)%3)->coordinates().y();
    z2 = (float)facet->point((i+1)%3)->coordinates().z();
    GfxDebug::draw_line(x1,y1,z1,x2,y2,z2,color);
  }
  GfxDebug::flush();
}

static void draw_facet_point( CubitPoint* pt, int color )
{
  float x = (float)pt->coordinates().x();
  float y = (float)pt->coordinates().y();
  float z = (float)pt->coordinates().z();
  GfxDebug::draw_point( x, y, z, color );
  GfxDebug::flush();
}

static void draw_facet_points( CubitFacet* facet, int color )
{
  draw_facet_point(facet->point(0), color);
  draw_facet_point(facet->point(1), color);
  draw_facet_point(facet->point(2), color);
}
*/


//-------------------------------------------------------------------------
// Purpose       : Split surface using facet data, given a coedge along
//                 the split chain.
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 02/23/03
//-------------------------------------------------------------------------
PartitionSurface* PartitionEngine::split_surface( PartitionSurface* surface,
                                                  PartitionCoEdge* new_coedge )
{
  int i;

    // get surface facets and make sure marks are cleared
  DLIList<CubitFacetData*> faces;
  surface->get_facet_data(faces);
  for( i = faces.size(); i--; )
    faces.get_and_step()->marked(0);


    // get curve and curve's facet edges
  PartitionCurve* new_curve = dynamic_cast<PartitionCurve*>(new_coedge->get_curve());
  DLIList<CubitFacetEdgeData*> edges;
  new_curve->get_facet_data( edges );
  assert(edges.size());
  
    // get the first point in the curve's facet edges
  edges.reset();
  CubitFacetEdgeData* edge = edges.get();
  CubitPoint* pt = edge->point(0);
  if ( TDVGFacetOwner::get(pt) != new_curve->start_point() )
  {
    pt = edge->point(1);
    assert(TDVGFacetOwner::get(pt) == new_curve->start_point());
  }
  
    // get the facets on the surface on the inside size of
    // the coedge
  bool coedge_forward = new_coedge->sense() == CUBIT_FORWARD;
  DLIList<CubitFacetData*> front( faces.size() );
  for ( i = edges.size(); i--; )
  {
    edge = edges.get_and_step();
    bool edge_forward = edge->point(0) == pt;
    bool sense = edge_forward == coedge_forward;
    CubitFacetData* face = find_facet( edge, sense, surface );
    if( face->marked() == 0 )
    {
//draw_facet_edges(face,CUBIT_RED);
      face->marked(1);
      front.append( face );
    }
    pt = edge->other_point(pt);
  }    

    // mark all adjacent facets traversing only over
    // edges that are in the interior of the surface
  while( front.size() )
  {
    CubitFacetData* face = front.pop();
    for( i = 0; i < 3; i++ )
    {
      if( face->edge(i) && TDVGFacetOwner::get(face->edge(i)) )
        continue;
        
      CubitFacetData* other = other_facet( face, i, surface );
      if ( other->marked() == 0 )
      {
//draw_facet_edges(face,CUBIT_BLUE);
        other->marked(1);
        front.append(other);
      }
    }
  }
  
    // reuse front to hold all the facets we marked
  for( i = faces.size(); i--; )
  {
    CubitFacetData* face = faces.get_and_step();
    if( face->marked() )
      front.append( face );
  }
  
    // if we got back all the facets, then the surface
    // was not split!
  if( front.size() == faces.size() )
  {
    for( i = faces.size(); i--; )
      faces.get_and_step()->marked(0);
      
    return surface;
  }
    // split off a new PartitionSurface
  PartitionSurface* new_surf = surface->split( front );
  
    // get list of all loops in surface
  DLIList<PartitionLoop*> loops( surface->num_loops() );
  PartitionLoop* loop = 0;
  while( (loop = surface->next_loop( loop )) )
    loops.append( loop );
  
    // check which loops we need to move to the new surface
  while( loops.size() )
  {
    loop = loops.pop();
    PartitionCoEdge* coe = loop->first_coedge();

    if( !coe )
      continue;

    PartitionCurve* curve = coe->get_curve();
    
    CubitFacetData* face =0;
    if( PartPTCurve* ptc = dynamic_cast<PartPTCurve*>(curve) )
    {
      face = find_facet( ptc->start_point()->facet_point(), new_surf );
    }
    else
    {
      edges.clean_out();
      curve->get_facet_data( edges );
      if (edges.size() > 0) // normal facet
      {
        edges.reset();
        edge = edges.get();
        CubitPoint* start_pt = curve->start_point()->facet_point();
        bool forward = edge->point(0) == start_pt;
        assert( forward || edge->point(1) == start_pt );
        if ( coe->sense() == CUBIT_REVERSED )
          forward = !forward;
        face = find_facet( edge, forward, new_surf );
      }
      else // this is a hardpoint there is a curve with no facets (length)
           // and the start and end points are the same
      {
        TBPoint* hardpoint = curve->start_point()->real_point();
        if (hardpoint && curve->start_point() == curve->end_point() )
        {
          CubitVector hardpoint_coord = hardpoint->coordinates(); 
          CubitVector new_closest, old_closest;
          new_surf->closest_point_trimmed( hardpoint_coord, new_closest );
          surface->closest_point_trimmed( hardpoint_coord,  old_closest );
          // It appears that the new_surf is really the original surface
          if ( new_closest.length_squared() > old_closest.length_squared() )
            face = reinterpret_cast<CubitFacetData*>(1); // just set the point to non-zero (true)
        }
      }
    }
    
    if( face )
    {
      surface->remove( loop );
      new_surf->add( loop );
    }
  }
/*  
  PartitionCoSurf* cos = 0;
  while( (cos = surface->next_co_surface( cos )) )
    cos->get_shell()->add( new_surf, cos->sense() );
*/  
  for( i = faces.size(); i--; )
    faces.get_and_step()->marked(0);
  
  return new_surf;
}

    
//-------------------------------------------------------------------------
// Purpose       : Public interface to surface partitioning.
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 02/23/03
//-------------------------------------------------------------------------
Surface* PartitionEngine::insert_curve( Surface* surface,
                         DLIList<CubitVector*>& segment_points,
                         DLIList<Curve*>& new_curves,
                         const double *tolerance_length)
{
  DLIList<Surface*> input_surfs(1), output_surfs(2);
  input_surfs.append(surface);
  if ( ! insert_curve(input_surfs, segment_points, output_surfs, new_curves,
                tolerance_length) )
    return 0;

  return output_surfs.size() ? output_surfs.get() : surface;
} 

//-------------------------------------------------------------------------
// Purpose       : Public interface to surface partitioning
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 02/23/03
//-------------------------------------------------------------------------
CubitStatus PartitionEngine::insert_curve( 
                                    DLIList<Surface*>& input_surfaces,
                                    DLIList<CubitVector*>& segment_points,
                                    DLIList<Surface*>& new_surfaces,
                                    DLIList<Curve*>& new_curves,
                                    const double *tolerance_length,
                                    DLIList<Surface*>* surfs_to_reverse)
{
  int i;
  CubitStatus result = CUBIT_SUCCESS;
  
  for ( i = input_surfaces.size(); i--; )
    if ( dynamic_cast<CompositeSurface*>(input_surfaces.get_and_step()) )
      return CompositeEngine::instance().insert_curve(
        input_surfaces, segment_points, new_surfaces, new_curves );
  
  DLIList<SubSurface*> replaced_surfaces;
  DLIList<PartitionSurface*> surface_list, new_part_surfs;
  DLIList<PartitionCurve*> curve_list;
  
    // get partition surfaces
  for ( i = input_surfaces.size(); i--; )
  {
    int reverse_surf = 0;
    Surface* surf = input_surfaces.get_and_step();

    // We will need to reverse the surface and its facets if it
    // is in the passed-in list of surfaces to reverse.
    if(surfs_to_reverse && surfs_to_reverse->is_in_list(surf))
      reverse_surf = 1;
    PartitionSurface* partsurf = dynamic_cast<PartitionSurface*>(surf);
    if( partsurf )
    {
      surface_list.append(partsurf);
      if(reverse_surf)
        partsurf->reverse_sense();
    }
    else
    {
      SubSurface* subsurf = replace_surface(surf);
      if (!subsurf)
      {
        result = CUBIT_FAILURE;
        break;
      }
      
      replaced_surfaces.append( subsurf );
      surface_list.append(subsurf);
      if(reverse_surf)
        subsurf->reverse_sense();
    }
  }
  
  if (result)
  {
    result = insert_curve( surface_list, segment_points, new_part_surfs, curve_list,
      tolerance_length);
  }
  if (!result)
  {
    // we replaced the original surface we need to clean up also
    for ( i = replaced_surfaces.size(); i--;  )
    {
      SubSurface* surf = replaced_surfaces.step_and_get();
      restore_surface( surf ); // this also deletes the virtual surface
    }

    return result;
  }

  CAST_LIST_TO_PARENT(curve_list,new_curves);
  CAST_LIST_TO_PARENT(new_part_surfs,new_surfaces);
  
    // clean up any partition surfaces we created but didn't modify
  for ( i = replaced_surfaces.size(); i--;  )
  {
    SubSurface* surf = replaced_surfaces.step_and_get();
    if ( !surf->sub_entity_set().has_lower_order() &&
         !surf->next_co_surface(0) )
    {
      restore_surface( surf );
    }
    else
    {
      new_surfaces.append( surf );
    }
  }
  
  return result;
} 

void PartitionEngine::curves_from_surfaces( const DLIList<PartitionSurface*>& surfs,
                                            DLIList<PartitionCurve*>& curves )
{
  int i;
  
  for (i = curves.size(); i--; )
    curves.step_and_get()->mark = 1;
  
  for (i = surfs.size(); i--; )
  {
    PartitionSurface* surf = surfs.next(i);
    PartitionLoop* loop = 0;
    while ( (loop = surf->next_loop(loop)) )
    {
      PartitionCoEdge* coedge = loop->first_coedge();
      do {
        PartitionCurve* curve = coedge->get_curve();
        if (!curve->mark)
        {
          curve->mark = 1;
          curves.append(curve);
        }
      
        coedge = loop->next_coedge(coedge);
      } while (coedge != loop->first_coedge());
    }
  }
  
  for ( i = curves.size(); i--; )
    curves.step_and_get()->mark = 0;
  curves.reset();
}


CubitStatus PartitionEngine::insert_curve( DLIList<PartitionSurface*>& surface_list,
                                           DLIList<CubitVector*>& segment_points,
                                           DLIList<PartitionSurface*>& new_part_surfs,
                                           DLIList<PartitionCurve*>& curve_list,
                                           const double *tolerance_length)
{
  int i;
  const double TOL_SQR = GEOMETRY_RESABS*GEOMETRY_RESABS;
  
    // Get curves
  DLIList<PartitionCurve*> curves;
  curves_from_surfaces( surface_list, curves );
  
    // Check for intersections with boundary curves
  segment_points.reset();
  curves.reset();
  for ( i = segment_points.size(); i--;  )
  {
    CubitVector& pos = *segment_points.get_and_step();
    bool debug = false;
    if (debug)
    {
      GfxDebug::draw_point(pos, CUBIT_BLUE);
      GfxDebug::mouse_xforms();
    }
    for ( int j = curves.size(); j--; )
    {
      PartitionCurve* curve = curves.get_and_step();
      if ((curve->start_point()->coordinates() - pos).length_squared() < TOL_SQR)
      {
        pos = curve->start_point()->coordinates();
        break;
      }
      
      if ((curve->end_point()->coordinates() - pos).length_squared() < TOL_SQR)
      {
        pos = curve->end_point()->coordinates();
        break;
      }

      CubitVector closest;
      curve->closest_point_trimmed( pos, closest );
        
      if ( (pos - closest).length_squared() > TOL_SQR )
        continue;
      
      if ((curve->start_point()->coordinates() - closest).length_squared() < TOL_SQR)
      {
        pos = curve->start_point()->coordinates();
        break;
      }
      
      if ((curve->end_point()->coordinates() - closest).length_squared() < TOL_SQR)
      {
        pos = curve->end_point()->coordinates();
        break;
      }
/*      
      double u = curve->u_from_position(closest);
      TBPoint* point = insert_point( curve, u );
      if (!point)
      {
        PRINT_ERROR("Error splitting curve.  Aborting.\n");
        return CUBIT_FAILURE;
      }
      pos = point->coordinates();
*/
      pos = closest;
      break;
    }
  }
      
  
    // get surface facets
  DLIList<CubitFacetData*> surf_facets, facet_list;
  surface_list.reset();
  for ( i = surface_list.size(); i--; )
  {
    surf_facets.clean_out();
    surface_list.get_and_step()->get_facet_data( surf_facets );
    facet_list += surf_facets;
  }

    // Do facet projection
  DLIList<CubitFacetEdgeData*> polyline;
  DLIList<CubitPointData*> polyline_pts;
  project_to_surface( facet_list, segment_points, polyline, polyline_pts, tolerance_length);
  if (segment_points.size() == polyline_pts.size())
  { 
    // Move facet points to real geometry
    DLIList<CubitFacetEdge*> edge_list;
    DLIList<CubitFacet*> pt_facets;
    segment_points.reset();
    polyline_pts.reset();
    int count = 0;
    for ( i = segment_points.size(); i--; )
    {
      CubitPointData* facet_pt = polyline_pts.get_and_step();
      CubitVector position( *segment_points.get_and_step() );
      if ( !facet_pt )  
        continue;
    
      PartitionPoint* point = dynamic_cast<PartitionPoint*>(TDVGFacetOwner::get(facet_pt));
      if ( point ) 
      {
        count++;
        continue;
      }

      PartitionCurve* curve = 0;
      edge_list.clean_out();
      facet_pt->edges( edge_list );
      while ( edge_list.size() )
      {
        PartitionEntity* owner = TDVGFacetOwner::get(edge_list.pop());
        if ( owner )
        {
          assert(!curve || curve == owner);
          curve = dynamic_cast<PartitionCurve*>(owner);
        }
      }
      if ( curve )
      {
        Surface* surface = dynamic_cast<Surface*>(curve->partitioned_entity());
        if ( surface )
        {
          CubitVector surf_pos;
          surface->closest_point( position, &surf_pos );
          if (facet_pt->check_inverted_facets(surf_pos))
            facet_pt->set( surf_pos );
        }
        else
          curve->relax_to_geometry( facet_pt, &position );

        continue;
      }
    
      PartitionSurface* surf = 0;
      pt_facets.clean_out();
      facet_pt->facets( pt_facets );
      while ( pt_facets.size() )
      {
        PartitionEntity* owner = TDVGFacetOwner::get(pt_facets.pop());
        if ( owner )
        {
          assert(!surf || surf == owner);
          surf = dynamic_cast<PartitionSurface*>(owner);
        }
      }
      if ( surf && dynamic_cast<SubSurface*>(surf) )
      {
        surf->relax_to_geometry( facet_pt, &position );
      }
    }    
      
    if (count < polyline_pts.size())
    {
      for( i = 0; i < polyline_pts.size(); i++ )
      {
        if( polyline_pts[i] )
        {
          *(segment_points[i]) = polyline_pts[i]->coordinates();
        }
        else
        {
          return CUBIT_FAILURE;
        }
      }
    }
  }

  if ( !polyline.size() )
    return CUBIT_FAILURE;
    
  
    // construct the curves
  DLIList<CubitFacetEdgeData*> cubit_edges;
  DLIList<CubitFacetEdge*> pt_edges;
  polyline.reset();
  for ( i = polyline.size(); i > 0; )
  {
      // get list of FacetPoints for Curve
    cubit_edges.clean_out();
    CubitFacetEdgeData* edge = polyline.get_and_step();
    if (!edge) { i--; continue; } 
    
    assert(!TDVGFacetOwner::get(edge) && 
           edge->num_adj_facets() < 3 &&
           edge->num_adj_facets() > 0);
           
    PartitionEntity* owner = TDVGFacetOwner::get(edge->adj_facet(0));
    assert( edge->num_adj_facets() == 1 ||
            TDVGFacetOwner::get(edge->adj_facet(1)) == owner );
            
    cubit_edges.append( edge );
    
    while( polyline.get() )
    {
      CubitPoint* pt = edge->shared_point(polyline.get());
      edge = polyline.get();
      
      assert(!TDVGFacetOwner::get(edge) && 
             edge->num_adj_facets() < 3 &&
             edge->num_adj_facets() > 0 && pt);
           
      PartitionEntity* tmp = TDVGFacetOwner::get(edge->adj_facet(0));
      assert( edge->num_adj_facets() == 1 ||
            TDVGFacetOwner::get(edge->adj_facet(1)) == tmp );
            
      bool edge_owner = false;
      pt->edges(pt_edges);
      while (pt_edges.size())
        if (TDVGFacetOwner::get(pt_edges.pop()))
          edge_owner = true;
            
      if ( tmp != owner || TDVGFacetOwner::get(pt) || edge_owner )
        break;
        
      cubit_edges.append( polyline.get_and_step() );
    }
    if (NULL == polyline.get())
    {
      polyline.step();
      i -= 1;
    }
    i -= cubit_edges.size();

      // create curve and partition surface
    PartitionCurve* new_curve = insert_curve( cubit_edges, &segment_points );
    if ( new_curve )
      curve_list.append(new_curve);
  }
  
  
    // find new surfaces to return
  for ( i = curve_list.size(); i--; )
  {
    PartitionCurve* curve = curve_list.step_and_get();
    PartitionCoEdge* coedge = 0;
    while( (coedge = curve->next_coedge(coedge)) )
    {
      new_part_surfs.append( coedge->get_loop()->get_surface() );
    }
  }
  
  if ( curve_list.size() )
  {
      // do facet cleanup
    for ( i = curve_list.size(); i--; )
      curve_list.get_and_step()->do_facet_cleanup();
    //for ( i = new_part_surfs.size(); i--; )
    //  new_part_surfs.get_and_step()->do_facet_cleanup();
  }
    
  
  for ( i = new_part_surfs.size(); i--; )
    new_part_surfs.step_and_get()->mark = 1;
  for ( i = surface_list.size(); i--; )
    surface_list.step_and_get()->mark = 0;
  for ( i = new_part_surfs.size(); i--; )
    if ( new_part_surfs.step_and_get()->mark )
      new_part_surfs.get()->mark = 0;
    else 
      new_part_surfs.change_to(0);
  new_part_surfs.remove_all_with_value(0);

  return curve_list.size() ? CUBIT_SUCCESS : CUBIT_FAILURE;
} 

//-------------------------------------------------------------------------
// Purpose       : Given a polyline-imprint on the surface facetting,
//                 construct a new and insert it in the surface topology.
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 02/23/03
//-------------------------------------------------------------------------
SegmentedCurve* PartitionEngine::insert_curve( 
                         DLIList<CubitFacetEdgeData*>& segments,
                         DLIList<CubitVector*>* /*curve_geom*/ )
{
  int i;
  
    // find surface
  PartitionSurface* surface = 0;
  DLIList<CubitFacet*> facets;
  for ( i = segments.size(); i--; )
  {
    CubitFacetEdgeData* edge = segments.get_and_step();
    facets.clean_out();
    edge->facets(facets);
    assert(facets.size() == 2);
    for (int j = 0; j < facets.size(); j++ )
    {
      CubitFacet* facet = facets.get_and_step();
      PartitionEntity* fowner = TDVGFacetOwner::get(facet);
      PartitionSurface* surf = dynamic_cast<PartitionSurface*>(fowner);
      if ( !surf || (surface && surf != surface) )
      {
        PRINT_ERROR("Internal Error in PartitionEngine.  Corrupt facet data.\n");
        assert(0);
        return 0;
      }
      surface = surf;
    }
  }
  
    // if there is a real surface, relax facet points to surface
  Surface* real_surf = dynamic_cast<Surface*>(surface->partitioned_entity());
  
    // find start and end facet points
  CubitPointData* pts[2];
  CubitPoint* pt;
  
  if( segments.size() == 1 )
  {
    pts[0] = dynamic_cast<CubitPointData*>(segments.get()->point(0));
    pts[1] = dynamic_cast<CubitPointData*>(segments.get()->point(1));
  }
  else
  {
    segments.last();
    pt = segments.get()->shared_point( segments.prev() );
    pt = segments.get()->other_point( pt );
    pts[1] = dynamic_cast<CubitPointData*>(pt);
    segments.reset();
    pt = segments.get()->shared_point( segments.next() );
    pt = segments.get()->other_point( pt );
    pts[0] = dynamic_cast<CubitPointData*>(pt);
  }
  
    // create end points (partition curves if necessary)
  bool okay = true;
  DLIList<CubitFacetEdge*> edges;
  bool created_pts[2] = {false,false};
  PartitionPoint* end_pts[2] = {0,0};
  for ( i = 0; okay && i < 2; i++ )
  {
    end_pts[i] = dynamic_cast<PartitionPoint*>(TDVGFacetOwner::get(pts[i]));
    if( end_pts[i] ) continue;
    
    edges.clean_out();
    pts[i]->edges(edges);
    PartitionCurve* curve = 0;
    for ( int j = edges.size(); !curve && j--; )
      curve = dynamic_cast<PartitionCurve*>(TDVGFacetOwner::get(edges.get_and_step()));
    
    if( curve )
    {
      curve->relax_to_geometry( pts[i] );
      PartitionPoint* new_pt = new PartitionPoint( pts[i]->coordinates(), curve );
      if ( insert_point( curve, new_pt ) )
      {
        end_pts[i] = new_pt;
        created_pts[i] = true;
        assert( TDVGFacetOwner::get(pts[i]) == new_pt );
      }
      else
      {
        delete new_pt;
        okay = false;
      }
      continue;
    }
    
    facets.clean_out();
    pts[i]->facets(facets);
    for ( int k = facets.size(); k--; )
      if( TDVGFacetOwner::get(facets.get_and_step()) != surface )
      {
        PRINT_ERROR("Internal Error in PartitionEngine.  Corrupt facet data.\n");
        okay = false;
        break;
      }
      
    if( okay )
    {
      if( real_surf )
        surface->relax_to_geometry( pts[i] );
      
      end_pts[i] = new PartitionPoint( pts[i]->coordinates(), surface );
      end_pts[i]->facet_point( pts[i] );
    }
  }
  
    // create curve
  SegmentedCurve* curve = 0;
  if( okay )
  {
    DLIList<CubitVector*> vectors(segments.size()+1), *curve_geom_ptr;
//    if( curve_geom )
//    {
//      curve_geom_ptr = curve_geom;
//    }
//    else
    {
      segments.reset();
      pt = pts[0];
      vectors.append( new CubitVector( pt->coordinates() ) );
      for( i = segments.size(); i--; )
      {
        pt = segments.get_and_step()->other_point(pt);
        if ( real_surf && i > 0 )
        {
          surface->relax_to_geometry( pt );
        }
        assert(!!pt);
        vectors.append( new CubitVector( pt->coordinates() ) );
      }
      curve_geom_ptr = &vectors;
    }
    
    curve = new SegmentedCurve( surface, *curve_geom_ptr );
    
    while( vectors.size() ) delete vectors.pop();
    
    curve->start_point( end_pts[0] );
    curve->end_point( end_pts[1] );
    curve->set_facet_data(segments);
  
    if( insert_curve( surface, curve ) )
      return curve;
  }
  
  // If we got this far, then something went wrong.
  // Try to clean up.
    
  if( curve )
  {
    curve->start_point(0);
    curve->end_point(0);
    segments.clean_out();
    curve->set_facet_data(segments);
    delete curve;
  }
  
  for( i = 0; i < 2; i++ )
  {
    if( ! created_pts[i] ) continue;
    
    end_pts[i]->facet_point(0);
    if( dynamic_cast<Curve*>(end_pts[i]->partitioned_entity()) )
      remove_point( end_pts[i] );
    else
      delete end_pts[i];
  }
  
  return 0;
}

/*
//-------------------------------------------------------------------------
// Purpose       : Partition a curve using a facet point on the curve.
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 02/23/03
//-------------------------------------------------------------------------
PartitionPoint* PartitionEngine::insert_point( CubitPointData* pt )
{
  DLIList<CubitFacetEdge*> edges;
  pt->edges(edges);
  
  PartitionCurve* curve = 0;
  for ( int i = edges.size(); i--; )
  {
    CubitFacetEdge* edge = edges.step_and_get();
    PartitionCurve* c = dynamic_cast<PartitionCurve*>(TDVGFacetOwner::get(edge));
    if( c )
    {
      if( !curve )
        curve = c;
      else if( curve != c )
      {
        PRINT_ERROR("Internal Error in PartitionEngine.  Corrupt facet data.\n");
        return 0;
      }
    }
  }
      
  DLIList<CubitFacetEdgeData*> segments, other_segments;
  curve->get_facet_data(segments);
  segments.last();
  while( !segments.get()->contains( pt ) )
  {
    other_segments.append( segments.pop() );
    segments.last();
  }
  other_segments.reverse();
  
  PartitionPoint* new_point = new PartitionPoint( pt->coordinates(), curve );
  PartitionCurve* curve2 = insert_point( curve, new_point );
  if( !curve2 )
  {
    delete new_point;
    return 0;
  }
  
  curve->set_facet_data(segments);
  curve2->set_facet_data(other_segments);
  new_point->facet_point(pt);
  return new_point;
}
*/

//-------------------------------------------------------------------------
// Purpose       : Add point to surface facetting.
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 02/23/03
//-------------------------------------------------------------------------
CubitPointData* PartitionEngine::project_to_surface( 
                                 PartitionSurface* surface,
                                 const CubitVector& position )
{
  
  CubitVector pos;
  CubitFacet* facet = surface->closest_facet( position, pos );
  surface->closest_point( position, &pos );
  CubitFacetData* facet_d = dynamic_cast<CubitFacetData*>(facet);
  DLIList<CubitFacetEdge*> pt_edges;
  assert(!!facet_d);

  CubitVector p[3], v[3], d[3];
  int i;
  
    // Initialize data
  double dist_tol_sqr = CUBIT_DBL_MAX;
  for ( i = 0; i < 3; i++ ) {
    p[i] = facet_d->point(i)->coordinates();
    v[i] = facet_d->point((i+1)%3)->coordinates() - p[i];
    double dist_sqr = v[i].length_squared();
    if ( dist_sqr < dist_tol_sqr )
      dist_tol_sqr = dist_sqr;
    d[i] = pos - p[i];
  }
  dist_tol_sqr *= 0.2;
  if( dist_tol_sqr < GEOMETRY_RESABS*GEOMETRY_RESABS )
    dist_tol_sqr = GEOMETRY_RESABS*GEOMETRY_RESABS;
  
    // Check if we are within GEOMETRY_RESABS of a vertex of the
    // triangle or if one of the triangle vertices can be moved
    // to the input location.
  int closest_pt = -1;
  double closest_dist = dist_tol_sqr;
  for ( i = 0; i < 3; i++ ) {
    double dist_sqr = d[i].length_squared();
    bool has_owner = (0 != TDVGFacetOwner::get(facet->point(i)));
    if( !has_owner && dist_sqr < closest_dist ) {
      pt_edges.clean_out();
      facet->point(i)->edges( pt_edges);
      while( pt_edges.size() )
        if ( TDVGFacetOwner::get(pt_edges.pop()) )
          has_owner = true;
    }
    
    if ( (!has_owner && dist_sqr < closest_dist) ||
         dist_sqr < GEOMETRY_RESABS*GEOMETRY_RESABS ) {
      closest_pt = i;
      closest_dist = dist_sqr;
    }
  }
  
  if( closest_pt >= 0 )
  {
    facet->point(closest_pt)->set(pos);
    return dynamic_cast<CubitPointData*>(facet->point(closest_pt));
  }
  
    // Check if the input position is within GEOMETRY_RESABS of one
    // of the triangle edges or if the edge can be 'bent' such that
    // the split point is the input position.
  int closest_edge = -1;
  for ( i = 0; i < 3; i++ ) {
    double t = (v[i] % d[i]) / v[i].length_squared();
    double dist_sqr = (d[i] - t * v[i]).length_squared();
    CubitFacetEdge* edge = facet->point(i)->shared_edge(facet->point((i+1)%3));

    bool has_owner = edge && TDVGFacetOwner::get(edge);
    
    if( (!has_owner && dist_sqr < closest_dist) ||
         dist_sqr < GEOMETRY_RESABS*GEOMETRY_RESABS )
    {
      CubitFacet* other_facet = 0;
      for ( int j = 0; j < edge->num_adj_facets(); j++ )
      {
        if (edge->adj_facet(j) != facet && 
            TDVGFacetOwner::get(edge->adj_facet(j)) == surface)
        {
          other_facet = edge->adj_facet(j);
          break;
        }
      }
      
      if (other_facet)
      {
          // check if splitting the edge will invert some other facet
        int edge_index = other_facet->edge_index( edge );
        CubitVector other_pt = other_facet->point(edge_index)->coordinates();
        CubitVector prev_pt = other_facet->point((edge_index+1)%3)->coordinates();
        CubitVector next_pt = other_facet->point((edge_index+2)%3)->coordinates();
        prev_pt -= other_pt;
        next_pt -= other_pt;
        CubitVector norm = prev_pt * next_pt;
        CubitVector new_e = pos - other_pt;
        
        if ( norm.vector_angle_quick(prev_pt, new_e) > CUBIT_PI ||
             norm.vector_angle_quick(new_e, next_pt) > CUBIT_PI )
          continue;
        
          // don't create adjacent sliver facets
        double prev_t = prev_pt % new_e / prev_pt.length_squared();
        prev_pt *= prev_t;
        prev_pt += other_pt;
        if ( (prev_pt - pos).length_squared() < dist_sqr )
          continue;
        
        double next_t = next_pt % new_e / next_pt.length_squared();
        next_pt *= next_t;
        next_pt += other_pt;
        if ( (next_pt - pos).length_squared() < dist_sqr )
          continue;
      }
        
      
      closest_dist = dist_sqr;
      closest_edge = i;
    }
  }

  if( closest_edge >= 0 )
  {      
    CubitPoint* pt1 = facet->point(closest_edge);
    CubitPoint* pt2 = facet->point((closest_edge+1)%3);
    CubitFacetEdge* edge = pt1->shared_edge(pt2);
    CubitFacetEdge* new_edge;
    CubitFacet* new_facet;
    CubitPoint* new_point;
    PartSurfFacetTool::split_edge(edge, pos, 0, new_point, new_edge, new_facet );
    if (dynamic_cast<SubSurface*>(surface))
      surface->relax_to_geometry( new_point, &position );
    return dynamic_cast<CubitPointData*>(new_point);
  }
  
  CubitFacet *new_tri1 = 0, *new_tri2 = 0;
  CubitPoint* pt1 = facet->point(0);
  CubitPoint* pt2 = facet->point(1);
  CubitPoint* pt3 = facet->point(2);
  CubitPoint* new_pt = facet->insert_point( pos, new_tri1, new_tri2 );
  new CubitFacetEdgeData( new_pt, pt1 );
  new CubitFacetEdgeData( new_pt, pt2 );
  new CubitFacetEdgeData( new_pt, pt3 );
  surface->notify_split( facet, new_tri1 );
  surface->notify_split( facet, new_tri2 );

  if (dynamic_cast<SubSurface*>(surface))
    surface->relax_to_geometry( new_pt, &position );

  return dynamic_cast<CubitPointData*>(new_pt);
}

void PartitionEngine::delete_facet( CubitFacet* p_dead )
{
assert(!TDVGFacetOwner::get(p_dead));
assert(dynamic_cast<CubitFacetData*>(p_dead) != NULL);

  CubitFacetEdge* edges[3] = {p_dead->edge(0), p_dead->edge(1), p_dead->edge(2)};
  CubitPoint* pts[3] = {p_dead->point(0), p_dead->point(1), p_dead->point(2)};

  int k;
  for( k = 0; k < 3; k++ )
    if( !pts[k]->num_adj_facets() )
      delete pts[k];

  for( k = 0; k < 3; k++ )
    if( edges[k] && !edges[k]->num_adj_facets() )
      delete edges[k];
  
  delete p_dead;
}        

//-------------------------------------------------------------------------
// Purpose       : Create projection of a polyline in surface facetting.
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 02/23/03
//-------------------------------------------------------------------------
CubitStatus PartitionEngine::project_to_surface( 
                                 DLIList<CubitFacetData*>& cubit_facets,
                                 DLIList<CubitVector*>& polyline,
                                 DLIList<CubitFacetEdgeData*>& polyline_out,
                                 DLIList<CubitPointData*>& polyline_points,
                                 const double *tolerance_length)
{
  int i, j;
  unsigned u;
  PartitionSurface* p_adj_surf;
  PartitionCurve* p_curve = 0;
  PartitionEntity* p_owner;

  DLIList<PartitionSurface*> surfaces;
  for (i = cubit_facets.size(); i--; )
    surfaces.append_unique( dynamic_cast<PartitionSurface*>(TDVGFacetOwner::get(cubit_facets.step_and_get())));

    // get points

  DLIList<CubitPoint*> cubit_points(cubit_facets.size() * 3);
  DLIList<CubitFacet*> orig_facets;
  CAST_LIST_TO_PARENT(cubit_facets, orig_facets);
  FacetDataUtil::get_facet_points(orig_facets, cubit_points);

  // get bounding curves from the facets
  std::set<PartitionCurve*> boundary_curves;
  DLIList<CubitFacetEdge*> free_edges;
  FacetDataUtil::edges_by_count(orig_facets, 1, free_edges);
  free_edges.reset();
  for (i=free_edges.size(); i>0; i--)
  {
    p_owner = TDVGFacetOwner::get(free_edges.get_and_step());
    p_curve = CAST_TO(p_owner, PartitionCurve);
    assert(!!p_curve); // every bounding edge should be owned by a curve
    boundary_curves.insert(p_curve);
  }

  std::vector<double> facet_points(cubit_points.size() * 3);
  std::vector<int> simple_facets(cubit_facets.size() * 3);
  cubit_points.reset();
  for ( i = 0; i < cubit_points.size(); i++ )
  {
    CubitPoint* pt = cubit_points.get_and_step();
    pt->marked( i );
    
    int j = i * 3;
    facet_points[j  ] = pt->coordinates().x();
    facet_points[j+1] = pt->coordinates().y();
    facet_points[j+2] = pt->coordinates().z();
    if (DEBUG_FLAG(145))
    {
      GfxDebug::draw_point(pt->coordinates(), CUBIT_WHITE);
      GfxDebug::draw_label(i, (float)(pt->coordinates().x()),
                              (float)(pt->coordinates().y()),
                              (float)(pt->coordinates().z()),
                              CUBIT_WHITE);
    }
  }
  
  cubit_facets.reset();
  for ( i = 0; i < cubit_facets.size(); i++ )
  {
    CubitFacetData* facet = cubit_facets.get_and_step();
    int j = i * 3;
    simple_facets[j  ] = facet->point(0)->marked();
    simple_facets[j+1] = facet->point(1)->marked();
    simple_facets[j+2] = facet->point(2)->marked();
    if (DEBUG_FLAG(145))
    {
      GfxDebug::draw_facet( facet, CUBIT_CYAN );
      CubitVector c = facet->center();
      GfxDebug::draw_label( i, (float)(c.x()),
                               (float)(c.y()),
                               (float)(c.z()),
                               CUBIT_CYAN);
      PRINT_DEBUG_145("%2d : %2d %2d %2d\n", i,
                                             facet->point(0)->marked(),
                                             facet->point(1)->marked(),
                                             facet->point(2)->marked());
    }
  }
  
    // clear point marks
  for ( i = 0; i < cubit_points.size(); i++ )
    cubit_points.step_and_get()->marked(0);  
/*
    // reduce polyline resolution
  const double dist_tol = GEOMETRY_RESABS*GEOMETRY_RESABS;
  const double COS15 = 0.96592582628906831;
  const double angle_tol = COS15*COS15;
  DLIList<CubitVector*> polyline2;
  polyline.reset();
  CubitVector *prev = polyline.get_and_step();
  polyline2.append(prev);
  for( i = polyline.size() - 2; i--; )
  {
    CubitVector *pt = polyline.get_and_step();
    CubitVector *next = polyline.get();
    CubitVector v1(*pt - *prev);
    CubitVector v2(*next - *pt);
    double len1 = v1.length_squared();
    double len2 = v2.length_squared();
    double dot = v1 % v2;
    if( len1 > dist_tol && len2 > dist_tol &&
        ((dot*dot)/(len1*len2)) < angle_tol )
    {
      polyline2.append(pt);
      prev = pt;
    }
  }
  polyline2.append(polyline.get());
*/
    // Call tool to do projection

  std::vector<double> new_facet_points;
  std::vector<int> dead_simple_facets, new_simple_facets, 
                   facet_replacement, new_polyline, seg_pts;
  FacetProjectTool tool;
  if ( !tool.project( polyline, facet_points, simple_facets,
                      dead_simple_facets, new_simple_facets,
                      facet_replacement, new_facet_points, 
                      new_polyline, seg_pts, tolerance_length ) )
  {
    // output the map between input polyline positions and
    // facet points
    std::vector<int>::iterator iitor;
    for ( iitor = seg_pts.begin(); iitor != seg_pts.end(); ++iitor )
    {
      int pt_index = *iitor;
      if ( pt_index == -1 )
      {
        polyline_points.append(NULL);
      }
      else
      {
        CubitPointData* ptd = NULL;
        if(pt_index >= 0 && pt_index < cubit_points.size())
        {
          CubitPoint* pt = cubit_points.next( pt_index );
          ptd = dynamic_cast<CubitPointData*>(pt);
        }
        polyline_points.append(ptd);
      }
    }

    return CUBIT_FAILURE;
  }

    // Convert back to CubitFacet rep
   
    // make new points
  assert( new_facet_points.size() % 3 == 0 );
  std::vector<CubitPointData*> new_cubit_points(new_facet_points.size()/3);
  std::vector<CubitPointData*>::iterator cpitor = new_cubit_points.begin();
  std::vector<double>::iterator ditor = new_facet_points.begin();
  for ( ; cpitor != new_cubit_points.end(); ++cpitor )
  {
    double x = *ditor++;
    double y = *ditor++;
    double z = *ditor++;
    *cpitor = new CubitPointData(x,y,z);
    cubit_points.append(*cpitor);
  }

    // make new facets
  cubit_facets.reset();
  cubit_points.reset();
  int prev_index = 0;
  assert( new_simple_facets.size() % 3 == 0 );

  std::vector<CubitFacetData*> facet_list;
  std::vector<CubitFacetData*>::iterator fitor;
  std::vector<int>::iterator iitor = new_simple_facets.begin();

  int id=0; // dummy for CubitFacetData constructor

  for( u = 0; u < dead_simple_facets.size(); u++ )
  {
      // get number of replacement facets
    int index = facet_replacement[u];
    int count = (index - prev_index) / 3;
    assert(((index - prev_index) % 3) == 0);

      // make a new list to hold data
    facet_list.clear();
    facet_list.resize(count+1);
    fitor = facet_list.begin();

      // put dead facet at beginning of list
    CubitFacetData* p_dead = cubit_facets.next(dead_simple_facets[u]);
    *fitor++ = p_dead;

    // get the surface that owns the dead facet
    p_owner = TDVGFacetOwner::get(p_dead);
    PartitionSurface* p_dead_surf = CAST_TO(p_owner, PartitionSurface);
    assert(!!p_dead_surf);

      // construct and append replacement facets
    for ( j=0; j < count; j++ )
    {
      CubitPoint* pt1 = cubit_points.next( *iitor++ );
      CubitPoint* pt2 = cubit_points.next( *iitor++ );
      CubitPoint* pt3 = cubit_points.next( *iitor++ );
      *fitor++ = new CubitFacetData(pt1, pt2, pt3, &id);
    }

      // set up for next iteration
    assert(fitor == facet_list.end());
    prev_index = index;

      // do something with "facet_list"

    DLIList<CubitFacetData*> old_facet_list;

    // fill edge replacement lists
    std::vector<CubitFacetEdgeData*> edge_lists[3];
    get_edge_replacements(facet_list, edge_lists);

    // update edges/facets on curves and adjacent surfaces
    CubitFacetEdgeData* old_edge;
    for (j=0; j<3; j++)
    {
      DLIList<CubitFacetEdgeData*> new_edges;
      if (edge_lists[j].size()) // check for replacements
      {
        std::vector<CubitFacetEdgeData*>::iterator eitor = edge_lists[j].begin();
        old_edge = *eitor++;

        p_owner = TDVGFacetOwner::get(old_edge);
        p_curve = CAST_TO(p_owner, PartitionCurve);

        for (; eitor != edge_lists[j].end(); eitor++)
          new_edges.append(*eitor);

        assert(new_edges.size() > 1); // shouldn't ever get 1 to 1 replacement here

        // split facets on adjacent surfaces if the split crosses a curve
        // on the outside boundary of the facets
        if (boundary_curves.find(p_curve) != boundary_curves.end())
        {
          // get facets to be replaced on adjacent surf
          DLIList<CubitFacetData*> adj_facets;
          DLIList<CubitFacet*> tmp_list;
          old_edge->facets(tmp_list);
          CAST_LIST(tmp_list, adj_facets, CubitFacetData);

          if (adj_facets.size() > 1)
          {
            DLIList<CubitFacetData*> adj_replacements;
            int i_adj;
            for (i_adj=adj_facets.size(); i_adj>0; i_adj--)
            {
              CubitFacetData* p_adjacent = adj_facets.get_and_step();
              if (p_adjacent == p_dead)
                continue;

              // create replacement facets
              adj_replacements.clean_out();

              // get the point opposite the edge being replaced
              int e_index = p_adjacent->edge_index(old_edge);
              bool old_reversed = -1 == p_adjacent->edge_use(e_index);
              CubitPoint* prev_pt = NULL;
              //if old_reversed is false, the new_edge' reversed flag is true
              if (old_reversed)
                 prev_pt = old_edge->point(0);
              else
                 prev_pt = old_edge->point(1);
              CubitPoint* opposite_pt = p_adjacent->point(e_index);
              //CubitPoint* opposite_pt = p_adjacent->point(p_adjacent->edge_index(old_edge)); // opposite edge and point have same index
              new_edges.reset();
              int i_new;
              for (i_new=new_edges.size(); i_new>0; i_new--)
              {
                CubitFacetEdge* p_edge = new_edges.get_and_step();
                
                prev_pt = p_edge->other_point(prev_pt);
                assert(!!prev_pt);
                bool new_reversed = p_edge->point(0) == prev_pt;
                int reversed = (int)(new_reversed != old_reversed);
                
                CubitFacetData* p_facet = new CubitFacetData(p_edge->point(reversed),
                                                             p_edge->point(1-reversed),
                                                             opposite_pt, &id);
                adj_replacements.append(p_facet);
              }

              p_owner = TDVGFacetOwner::get(p_adjacent);
              p_adj_surf = CAST_TO(p_owner, PartitionSurface);
              assert(!!p_adj_surf);
              old_facet_list.clean_out();
              old_facet_list.append(p_adjacent);
              p_adj_surf->replace_facets(old_facet_list, adj_replacements);
              delete_facet( p_adjacent );
            }
          }
        }

        // replace the curve edges
        p_curve->replace_facet(old_edge, new_edges);
        
        // Don't delete edge yet.  It's still connected to
        // at least one facet.  Wait 'till the facet is deleted.
        //delete old_edge;
      }
    }

    // replace dead facet
    fitor = facet_list.begin();
    fitor++; // skip dead facet
    DLIList<CubitFacetData*> new_facets;
    for (; fitor != facet_list.end(); fitor++)
      new_facets.append(*fitor);
    old_facet_list.clean_out();
    old_facet_list.append(p_dead);
    p_dead_surf->replace_facets(old_facet_list, new_facets);
    delete_facet(p_dead);
  }

  assert(iitor == new_simple_facets.end());

  // output the edges that make up the projected curve
  cubit_points.reset();
  iitor = new_polyline.begin();
  while( *iitor )
  {
    int count = *iitor++;
    CubitPoint* prev_point = cubit_points.next( *iitor++ );
    for ( i = 1; i < count; i++ ) {
      CubitPoint* pt = cubit_points.next( *iitor++ );
      CubitFacetEdge* edge = prev_point->shared_edge(pt);
      CubitFacetEdgeData* edge_d = dynamic_cast<CubitFacetEdgeData*>(edge);
      assert(!!edge_d);
      if (TDVGFacetOwner::get(edge_d)) 
        polyline_out.append(0);
      else
      {
        if (edge_d->point(0) == pt)
          edge_d->flip();
        polyline_out.append(edge_d);
      }
      prev_point = pt;
    }
    polyline_out.append(0);
  }
  
  // output the map between input polyline positions and 
  // facet points
  for ( iitor = seg_pts.begin(); iitor != seg_pts.end(); ++iitor )
  {
    int pt_index = *iitor;
    if ( pt_index == -1 )
    {
      polyline_points.append(NULL);
    }
    else
    {
      assert(pt_index >= 0 && pt_index < cubit_points.size());
      CubitPoint* pt = cubit_points.next( pt_index );
      CubitPointData* ptd = dynamic_cast<CubitPointData*>(pt);
      assert(!!ptd);
      polyline_points.append(ptd);
    }
  }  
  

  return CUBIT_SUCCESS;
}
  
//-------------------------------------------------------------------------
// Purpose       : Create a point-curve
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 02/23/03
//-------------------------------------------------------------------------
TBPoint* PartitionEngine::insert_point_curve( Surface* surf,
                                            const CubitVector& position,
                                            Surface *&partitioned_surf )
{
  if ( CompositeSurface* csurf = dynamic_cast<CompositeSurface*>(surf) )
  { 
    TBPoint* result = 
      CompositeEngine::instance().insert_point_curve(csurf, position); 
    partitioned_surf = csurf;
    return result;
  }
  
  bool replaced_surface = false;
  PartitionSurface* psurf = dynamic_cast<PartitionSurface*>(surf);
  if ( ! psurf ) {
    psurf = replace_surface(surf);
    replaced_surface = true;
    if ( ! psurf ) return 0;
  }
  
  CubitPointData* facet_point = project_to_surface( psurf, position );
  if ( !facet_point ) {
    if( replaced_surface )
      restore_surface( dynamic_cast<SubSurface*>(psurf) );
    partitioned_surf = NULL;
    return 0;
  }
  
  PartitionPoint* result = insert_point_curve( psurf, facet_point );
  if ( !result ) {
    if( replaced_surface )
      restore_surface( dynamic_cast<SubSurface*>(psurf) );
    return 0;
  }
  
  PartitionCurve* curve = result->next_curve(0);
  if (replaced_surface && !dynamic_cast<PartPTCurve*>(curve))
  {
    TBPoint* real = result->real_point();
    partitioned_surf = restore_surface( dynamic_cast<SubSurface*>(psurf) );
    if (real)
    {
      result = dynamic_cast<PartitionPoint*>(real->owner());
      return result ? result : real;
    }
  }
  else 
    partitioned_surf = psurf; 

  return result;
}



//-------------------------------------------------------------------------
// Purpose       : Create point-curve
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 02/23/03
//-------------------------------------------------------------------------
PartitionPoint* PartitionEngine::insert_point_curve( PartitionSurface* surf,
                                                     CubitPointData* point )
{
    // If imprint at existing vertex, just return that vertex
  if ( TDVGFacetOwner::get(point) )
    return dynamic_cast<PartitionPoint*>(TDVGFacetOwner::get(point));
  
    // Check if imprint point is on a curve of the surface
  DLIList<CubitFacetEdge*> edges;
  point->edges(edges);
  PartitionEntity* owner = 0;
  int owner_count = 0;
  int owner_edges = 0;
  while (edges.size())
  {
    if (PartitionEntity* ent = TDVGFacetOwner::get(edges.pop()))
    {
      if (owner != ent)
        owner_count++;
      owner = ent;
      owner_edges++;
    }
  }
  
    // If imprint is on curve, split the curve
  if (owner)
  {
    PartitionCurve* curve = dynamic_cast<PartitionCurve*>(owner);

      // check valid facet configuration.  should be only one
      // owner, and owner should own exactly two of the adjacent
      // edges.  And owner should be a partition curve.
    if (!curve || owner_count != 1 || owner_edges != 2)
    {
      PRINT_ERROR("Internal error at %s:%d\n  Invalid state.\n", __FILE__, __LINE__ );
      assert(0);
      return NULL;
    }
    
    PartitionPoint* new_point = new PartitionPoint(point->coordinates(), curve);
    new_point->facet_point(point);
    curve = insert_point( curve, new_point );
    if (!curve)
    {
      delete new_point;
      return 0;
    }
    
    return new_point;
  }
      
    // If here, point is in surface interior.  Create point-curve.
  PartitionPoint* new_point = new PartitionPoint( point->coordinates(), surf );
  new_point->facet_point(point);
  PartPTCurve* new_curve = new PartPTCurve( surf );
  new_curve->start_point( new_point );
  new_curve->end_point( new_point );
  
  insert_point_curve( surf, new_curve );
  
  return new_point;
}

CubitStatus PartitionEngine::insert_point_curve( PartitionSurface* surf,
                                                 PartPTCurve* new_curve,
                                                 bool update_topology )
{
  PartitionCoEdge* new_coedge = new PartitionCoEdge( surf, CUBIT_FORWARD );
  new_curve->add( new_coedge );
  PartitionLoop* new_loop = new PartitionLoop();
  new_loop->insert_after( new_coedge, 0 );
  
  surf->add( new_loop );
  if ( surf->owner() && update_topology )
    surf->owner()->notify_topology_modified(surf);
  return CUBIT_SUCCESS;
}
  


//-------------------------------------------------------------------------
// Purpose       : Construct surface for volume partitioning
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 01/13/03
//-------------------------------------------------------------------------
Lump* PartitionEngine::insert_surface( Surface* surface, Lump* lump)
{
  GMem gmem;
  int i;
  CubitStatus status = surface->get_geometry_query_engine()->
    get_graphics( surface, &gmem );
  
  if( !status )
  {
    PRINT_ERROR("PartitionEngine::insert_surface: GeometryEngine::"
                "get_graphics failed for surface.\n");
    return 0;
  }

  int num_pts = gmem.pointListCount;
  
  CubitPointData** ptarray = new CubitPointData*[num_pts];
  GPoint *p_itor = gmem.point_list();
  GPoint *p_end = p_itor + num_pts;
  CubitPointData** array_itor = ptarray;
  for( ; p_itor < p_end; p_itor++, array_itor++ )
    *array_itor = new CubitPointData(CubitVector(p_itor->x,p_itor->y,p_itor->z));
  
  DLIList<CubitFacetData*> facets;
  int* f_itor = gmem.facet_list();
  int* f_end = f_itor + gmem.fListCount;
  i = 0;
  for( ; f_itor < f_end ; f_itor += (1 + *f_itor) )
  {
    if( *f_itor != 3 )
      PRINT_ERROR("Non-triangular facet in PartitionEngine::insert_surface\n");
    else
      facets.append( new CubitFacetData(ptarray[f_itor[1]],
                                        ptarray[f_itor[2]],
                                        ptarray[f_itor[3]], &i) );
  }
    
    // clean up any unused points
  for( i = 0; i < num_pts; i++ )
    if( ptarray[i]->num_adj_facets() == 0 )
      delete ptarray[i];
  delete [] ptarray;
  
  
    // create partition lump
  bool created_lump = false;
  PartitionLump* partlump = dynamic_cast<PartitionLump*>(lump);
  if( !partlump )
  {
    created_lump = true;
    partlump = replace_lump(lump);
    assert(partlump != NULL);
  }
  
    // Get vertex locations
  DLIList<CubitVector*> coords;
  DLIList<TopologyBridge*> loops, coedges, curves, points, all_points;
  surface->get_children_virt(loops);
  while( loops.size() )
  {
    loops.pop()->get_children( coedges, true, layer() - 1 );
    while ( coedges.size() )
    {
      coedges.pop()->get_children_virt(curves);
      assert(curves.size() == 1);
      TopologyBridge* curve = curves.pop();
      curve->get_children( points, true, layer() - 1 );
      assert ( points.size() < 3 );
      all_points.merge_unique( points );
      points.clean_out();
    }
  }
  while( all_points.size() )
    coords.append( new CubitVector( dynamic_cast<TBPoint*>(all_points.pop())->coordinates() ) );
  
  PartitionLumpImprint tool( partlump );
  DLIList<PartitionEntity*> new_entities;
  PartitionSurface* new_surf = tool.imprint( facets, coords, new_entities );
  if ( !new_surf )
  {
    if( created_lump ) restore_lump(partlump);
    return 0;
  }
  
  PartitionLump* result = insert_surface( partlump, new_surf );
  if( !result )
  {
    destroy_surface( new_surf );
    if( created_lump ) restore_lump( partlump );
    return 0;
  }

  return result;
}

Surface* PartitionEngine::insert_surface( DLIList<CubitFacet*>& facets, Lump* lump)
{
  
    // create partition lump and surface
  bool created_lump = false;
  PartitionLump* partlump = dynamic_cast<PartitionLump*>(lump);
  if( !partlump )
  {
    created_lump = true;
    partlump = replace_lump(lump);
    if (!partlump)
      return 0;
  }
  
  PartitionLumpImprint tool( partlump );
  DLIList<PartitionEntity*> new_list;
  PartitionSurface* new_surf = tool.imprint( facets, new_list );
  if ( !new_surf )
  {
    if( created_lump ) restore_lump(partlump);
    return 0;
  }
  
  PartitionLump* result = insert_surface( partlump, new_surf );
  if( !result )
  {
    destroy_surface( new_surf );
    if( created_lump ) restore_lump( partlump );
    return 0;
  }
  
  return new_surf;
}


  
//-------------------------------------------------------------------------
// Purpose       : Insert a surface into a volume
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 01/12/03
//-------------------------------------------------------------------------
PartitionLump* PartitionEngine::insert_surface( PartitionLump* lump, 
                                           PartitionSurface* surface )
{
  assert( &(lump->sub_entity_set()) == &(surface->sub_entity_set()) );
  
    // Find all affected shells.  If any curve is not part of an
    // existing shell in the volume, set all_curves_connected to
    // false.  If any curve does not intersect the volume, then
    // the volume does not need to be split into two (or the surface
    // forms a new void in the volume.)
  bool all_curves_connected = true;
  DLIList<PartitionShell*> modified_list;
  PartitionLoop* loop = 0;
  PartitionCoEdge *coedge, *curve_coedge;
  PartitionCurve* curve;
  PartitionShell* shell;
  PartitionCoSurf* cosurf;
  while( (loop = surface->next_loop(loop)) )
  {
    coedge = loop->first_coedge();
    do {
      curve = coedge->get_curve();
      curve_coedge = 0;
      shell = 0;
      while( (curve_coedge = curve->next_coedge(curve_coedge)) )
      {
        PartitionSurface* surf = curve_coedge->get_loop()->get_surface();
        cosurf = 0;
        while( (cosurf = surf->next_co_surface( cosurf )) )
        {
          if( cosurf->get_shell()->get_lump() == lump )
          {
            assert( !shell || cosurf->get_shell() == shell );
            shell = cosurf->get_shell();
          }
        }
      }
      
      if( shell )
        modified_list.append_unique(shell);
      else
        all_curves_connected = false;
      
      coedge = loop->next_coedge(coedge);
    } while( coedge != loop->first_coedge() );
  }

    // Construct CoSurfaces
  PartitionCoSurf* cosurf1 = new PartitionCoSurf(CUBIT_FORWARD);
  PartitionCoSurf* cosurf2 = new PartitionCoSurf(CUBIT_REVERSED);
  surface->add(cosurf1);
  surface->add(cosurf2);
  
    // If there were no affected shells, we need to make a new
    // shell for a possible void in the volume.
  PartitionShell* shell_to_split = 0;
  if( modified_list.size() == 0 )
  {
    shell_to_split = new PartitionShell();
    lump->add(shell_to_split);
    
      // If surface is closed, set all_curves_connected
      // to true to indicate that volume needs to be split.
      // (Creating a one-surface void in the volume.)
    loop = 0;
    all_curves_connected = true;
    while( (loop = surface->next_loop(loop)) )
    {
      coedge = loop->first_coedge();
      do
      {
        curve = coedge->get_curve();
        curve_coedge = 0;
        int coedge_count = 0;
        while( (curve_coedge = curve->next_coedge(curve_coedge)) )
        {
          if( curve_coedge->get_loop() == loop )
            coedge_count++;
          else 
            assert(curve_coedge->get_loop()->get_surface() != surface);
        }
    
        if( coedge_count == 1 )
          all_curves_connected = false;
        else
          assert(coedge_count == 2);
      
        coedge = loop->next_coedge(coedge);
      } while( coedge != loop->first_coedge() );
    }
  }
  
    // Otherwise combine all the affected shells
  else
  {
    modified_list.reverse();
    shell_to_split = modified_list.pop();
    while( modified_list.size() )
    {
      shell = modified_list.pop();
      while( PartitionCoSurf* cosurf = shell->next_co_surface() )
      {
        shell->remove(cosurf);
        shell_to_split->add(cosurf);
      }
      lump->remove(shell);
      delete shell;
    }
  }
  
  shell_to_split->add(cosurf1);
  shell_to_split->add(cosurf2);
  
    // If we don't need to split the volume, we are done
  if( !all_curves_connected )
    return lump;
    // We probably need to split the the lump...
  else    
    return split_lump( shell_to_split );
}

//-------------------------------------------------------------------------
// Purpose       : Test if a PartitionShell needs to be split, and
//                 if so partition the shell and owning lump
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 01/13/03
//-------------------------------------------------------------------------
PartitionLump* PartitionEngine::split_lump( PartitionShell* shell_to_split )
{
    // Split shell
  PartitionShell* new_shell = split_shell(shell_to_split);
  if( !new_shell )
  {
      // shell did not need to be split
    return shell_to_split->get_lump();
  }
   
    // Create the new lump and shell
  PartitionLump* new_lump = new PartitionLump( shell_to_split->get_lump() );
  new_lump->add( new_shell );
  
    // Check for any remaining shells that need to be moved
    // from the old lump to the new lump.
  PartitionShell* shell = 0;
  PartitionLump* lump = shell_to_split->get_lump();
  while( (shell = lump->next_shell(shell)) )
  {
    if( shell == shell_to_split )
      continue;
    
    bool some_contained = false;
    bool all_contained = true;
    bool all_tested = true;
    PartitionCoSurf* cosurf = 0;
    while( (cosurf = shell->next_co_surface(cosurf)) )
    {
      PartitionSurface* surf = cosurf->get_surface();
      CubitVector pt, tmp = surf->bounding_box().center();
      surf->closest_point_trimmed( tmp, pt );
      switch( new_shell->point_containment(pt) )
      {
        case CUBIT_PNT_INSIDE:
          some_contained = true;
          break;
        case CUBIT_PNT_OUTSIDE:
          all_contained = false;
         break;
        default:
          all_tested = false;
          break;
      }
    }
    
    
    if( !all_tested || (some_contained && !all_contained) )
    {
      PRINT_ERROR("Error splitting volume.  Could not determine\n"
                  "which volume shell belongs in.  Geometry may be\n"
                  "invalid.\n");
    }
    
    if(some_contained && all_contained)
    {
      lump->remove(shell);
      new_lump->add(shell);
    }
  }
  
  return new_lump;
}

/*
#include "GfxDebug.hpp"
static void draw_curve( PartitionCurve* curve, int color, bool flush = true )
{
  GMem graphics;
  int i;
  curve->get_geometry_query_engine()->get_graphics( curve, i, &graphics );
  GfxDebug::draw_polyline(graphics.point_list(),
                                          graphics.pointListCount,
                                          color);
  if(flush) GfxDebug::flush();
}

static void draw_surface( PartitionSurface* surface, int color, bool flush = true )
{
//  RefEntity* re = dynamic_cast<RefEntity*>(surface->owner());
//  if( re )
//    GfxDebug::draw_ref_entity( re, color, true, false );
//  else 
  for( PartitionLoop* loop = 0; loop = surface->next_loop(loop); )
  {
    PartitionCoEdge* coe = loop->first_coedge();
    do {
      draw_curve(coe->get_curve(), color, false);
      coe = loop->next_coedge(coe);
    } while(coe != loop->first_coedge());
  }
  if(flush) GfxDebug::flush();
}

static void draw_cosurface( PartitionCoSurf* cosurf, int color, bool flush = true )
{
  //draw_surface(cosurf->get_surface(), color, false);
  CubitVector p, n;
  n = cosurf->get_surface()->bounding_box().center();
  cosurf->get_surface()->closest_point_trimmed( n, p );
  cosurf->get_surface()->closest_point( p, 0, &n );
  double l = cosurf->get_surface()->bounding_box().diagonal().length() * 0.1;
  n.length(l);
  if(cosurf->sense() == CUBIT_REVERSED)
    n *= -1.0;
    GfxDebug::draw_vector(p, p+n, color);
  if(flush) GfxDebug::flush();
}
*/

//-------------------------------------------------------------------------
// Purpose       : Test if a PartitionShell needs to be split, and
//                 split it
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 01/13/03
//-------------------------------------------------------------------------
PartitionShell* PartitionEngine::split_shell( PartitionShell* shell_to_split )
{
    // Make sure all cosurface surface marks are cleared
  PartitionCoSurf* cosurf = 0;
  while( (cosurf = shell_to_split->next_co_surface( cosurf )) )
  {
    cosurf->mark = 0;
    cosurf->get_surface()->mark = 0;
  }
  
    // Identify non-manifold surfaces, marking them either
    // with a 2 if they can be used to split the volume
    // (if they are part of a connected patch for which the
    // bounadary of that patch intersects the volume boundary
    // at all curves) or a 3 if they are other non-manifold
    // surfaces.  This will get a bit tricky if there are
    // non-manifold surfaces hanging off of the patch of
    // split surfaces.
  
    // First for all non-manifold surfaces, if the surface has
    // a curve that is not shared with any other surface, mark
    // it with a 3.  Otherwise mark it with a 2.
  cosurf = 0;
  DLIList<PartitionSurface*> surf_stack;
  while ( (cosurf = shell_to_split->next_co_surface(cosurf)) )
  {
    PartitionSurface* surf = cosurf->get_surface();
      // If we haven't done this surface yet and it is non-manifold
    if ( !surf->mark && surf->find_next(cosurf) )
    {
      surf->mark = 2;
      bool no_free_curve = true;
      PartitionLoop* loop = 0;
      while( no_free_curve && (loop = surf->next_loop(loop)) ) 
      {
        PartitionCoEdge* coedge = loop->first_coedge();
        do 
        {
          PartitionCurve* curve = coedge->get_curve();
            // If the curve has more than one coedge, it 
            // is not a free curve (this also accounts for
            // the case where the curve is a non-manifold 
            // curve on the surface interioir -- e.g. a sipe)
          if ( !curve->next_coedge(coedge) && 
               curve->next_coedge(0) == coedge )
          {
            no_free_curve = false;
            break;
          }
          coedge = loop->next_coedge(coedge);
        } while( coedge != loop->first_coedge() );
      }
      
      if( !no_free_curve )
      {
        surf->mark = 3;
        surf_stack.append( surf );
      }
    }
  }
  
    // Now for each surface we marked with a three, traverse
    // and mark adjacent surfaces until we come to a curve
    // connected to more that two surfaces.
  while( surf_stack.size() ) 
  {
    PartitionSurface* surf = surf_stack.pop();
    PartitionLoop* loop = 0;
    while ( (loop = surf->next_loop(loop)) )
    {
      PartitionCoEdge* coedge = loop->first_coedge();
      do 
      {
        PartitionCurve* curve = coedge->get_curve();
        int split_count = 0;
        int boundary_count = 0;
        PartitionCoEdge* curve_coe = 0;
        while ( (curve_coe = curve->next_coedge(curve_coe) ) != NULL )
        {
          PartitionSurface* curve_surf = curve_coe->get_loop()->get_surface();
          switch ( curve_surf->mark ) 
          {
            case 0: boundary_count++; break;
            case 2: split_count++;    break;
          }
        }
        
        if ( split_count == 1 && !boundary_count )
        {
          curve_coe = 0;
          while ( (curve_coe = curve->next_coedge(curve_coe) ) != NULL )
          {
            PartitionSurface* curve_surf = curve_coe->get_loop()->get_surface();
            if ( curve_surf->mark == 2 )
            {
              curve_surf->mark = 3;
              surf_stack.append(curve_surf);
            }
          }
        }
        
        coedge = loop->next_coedge( coedge );
      } while( coedge != loop->first_coedge() );
    }
  }
  
    // Now build a new shell by traversing cofaces, marking
    // each with that will go in a new shell with a 1.
  
    // Start with any cosurf that does not have a free
    // non-manifold surface (marked with a 3).  We'll handle
    // free non-manifold surfaces later.
  DLIList<PartitionCoSurf*> cosurf_stack;
  cosurf = 0;
  while ( (cosurf = shell_to_split->next_co_surface(cosurf)) )
    if ( cosurf->get_surface()->mark != 3 )
      break;
  if ( cosurf )
  {
    cosurf->mark = 1;
    cosurf_stack.append( cosurf );
  }
  
    // Traverse over adjacent cosurfaces, marking them with a 1
  while (cosurf_stack.size())
  {
    cosurf = cosurf_stack.pop();
    PartitionSurface* surf = cosurf->get_surface();
    PartitionLoop* loop = 0;
    while ( (loop = surf->next_loop(loop)) )
    {
      PartitionCoEdge* coedge = loop->first_coedge();
      do
      {
        PartitionCurve* curve = coedge->get_curve();
        PartitionCoEdge* curve_coe = 0;
        PartitionCoSurf *boundary_cosurf = 0, *split_cosurf = 0;
        int split_cosurf_count = 0;
        while ( (curve_coe = curve->next_coedge(curve_coe)) )
        {
          if ( curve_coe == coedge )
            continue;
          
          bool same_coe_sense = curve_coe->sense() == coedge->sense();
          PartitionSurface* curve_surf = curve_coe->get_loop()->get_surface();
          PartitionCoSurf* curve_cosurf = 0;
          while ( (curve_cosurf = curve_surf->next_co_surface(curve_cosurf)) )
          {
            if ( curve_cosurf->get_shell() != shell_to_split )
              continue;
            
            bool same_cos_sense = curve_cosurf->sense() == cosurf->sense();
            if ( same_cos_sense == same_coe_sense )
              continue;
            
              // Always choose split surface first if we
              // found one
            if ( curve_cosurf->get_surface()->mark == 2 ) {
              split_cosurf_count++;
              split_cosurf = curve_cosurf;
            }
            
              // Skip other non-manifold surfaces.  We'll
              // handle those later.
            else if( curve_cosurf->get_surface()->mark != 3 )
              boundary_cosurf = curve_cosurf;
          }
        }
        
        PartitionCoSurf* next_cosurf = split_cosurf ? split_cosurf : boundary_cosurf;
        if ( !next_cosurf->mark && split_cosurf_count < 2 )
        {
          next_cosurf->mark = 1;
          cosurf_stack.append(next_cosurf);
        }
      
        coedge = loop->next_coedge(coedge);
      } while( coedge != loop->first_coedge() );
    } // end while (loop)
  } // end while (cosurf_stack.size())
    
  
    // build lists of cosurfaces, one for each shell and
    // one of other non-manifold surfaces
  DLIList<PartitionCoSurf*> marked_list, unmarked_list, other_list;
  while( (cosurf = shell_to_split->next_co_surface(0)) )
  {
    shell_to_split->remove(cosurf);
    if ( cosurf->get_surface()->mark == 3 )
      other_list.append( cosurf );
    else if( cosurf->mark )
      marked_list.append( cosurf );
    else
      unmarked_list.append( cosurf );
  }
  
    // If one of marked_list or unmarked_list is empty,
    // we can't split the shell yet.  Put cofaces back in
    // shell and exit.
  if ( !marked_list.size() || !unmarked_list.size() )
  {
    marked_list += unmarked_list;
    marked_list += other_list;
    marked_list.reverse();
    while ( marked_list.size() )
    {
      cosurf = marked_list.pop();
      cosurf->mark = 0;
      cosurf->get_surface()->mark = 0;
      shell_to_split->add( cosurf );
    }
    return 0;
  }
  
    // Put unmarked list back in old shell
  unmarked_list.reverse();
  while ( unmarked_list.size() )
  {
    cosurf = unmarked_list.pop();
    cosurf->get_surface()->mark = 0;
    shell_to_split->add(cosurf);
  }
  
    // Put marked list in new shell
  PartitionShell* new_shell = new PartitionShell;
  marked_list.reverse();
  while ( marked_list.size() )
  {
    cosurf = marked_list.pop();
    cosurf->mark = 0;
    cosurf->get_surface()->mark = 0;
    new_shell->add(cosurf);
  }
  
    // Now sort out other non-manifold surfaces
  
    // Clear marks and get list of surface from cosurfaces
  surf_stack.clean_out();
  while( other_list.size() )
  {
    cosurf = other_list.pop();
    PartitionSurface* surf = cosurf->get_surface();
    if ( surf->mark )
    {
      surf->mark = 0;
      surf_stack.append(surf);
    }
  }
  
  insert_nonmanifold_surfaces( surf_stack, shell_to_split, new_shell );
  return new_shell;
}
  
  
//-------------------------------------------------------------------------
// Purpose       : After a shell is split, determine which of the two
//                 resulting shells each non-manifold surface belongs in.
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 02/23/03
//-------------------------------------------------------------------------
void PartitionEngine::insert_nonmanifold_surfaces( 
                                   DLIList<PartitionSurface*>& surf_stack,
                                   PartitionShell* shell1,
                                   PartitionShell* shell2 )
{
  DLIList<PartitionSurface*> known_list(surf_stack.size()), 
                           unknown_list(surf_stack.size());
  
  PartitionCoSurf* cosurf;
  PartitionSurface* surf;
  PartitionLoop* loop;
  PartitionCoEdge *coedge, *curve_coe;
  PartitionCurve* curve;
  
    // Loop until we've placed all the surfaces in one
    // shell or the other.
  while ( surf_stack.size() )
  {
    bool did_some = false;
    
      // Put any surfaces for which we immediately
      // know the shell into the appropriate shell.
      // Put others in known_list or unknown_list
      // depending on if we can determine which shell
      // they go in using a geometric comparison.
    known_list.clean_out();
    unknown_list.clean_out();

      // Take all surfaces out of stack in this loop.
      // We might put some back in after the loop (thus
      // the outer loop.)
    while ( surf_stack.size() )
    {
      surf = surf_stack.pop();
      
      PartitionShell* known_shell = 0;
      bool found_shell = false;
      loop = 0;
      while ( (loop = surf->next_loop(loop)) )
      {
        coedge = loop->first_coedge();
        do
        {
          curve = coedge->get_curve();
          curve_coe = 0;
          while ( (curve_coe = curve->next_coedge(curve_coe)) )
          {
            PartitionSurface* surf = curve_coe->get_loop()->get_surface();
            cosurf = 0;
            while ( (cosurf = surf->next_co_surface( cosurf )) )
            {
              if( cosurf->get_shell() == shell1 ||
                  cosurf->get_shell() == shell2 )
              {
                found_shell = true;
                if ( known_shell && known_shell != cosurf->get_shell() )
                  known_shell = 0;
                else
                  known_shell = cosurf->get_shell();
              }
            }
          }
        
          coedge = loop->next_coedge(coedge);
        } while( coedge != loop->first_coedge() );
      } // end while(loop)
        
        // This surface does not intersect the shell at
        // any curve.  We can't do this one yet.
      if ( !found_shell )
      {
        unknown_list.append( surf );
        continue;
      }
      
        // This surface intersected both shells at some
        // curves, but did not have a curve that intersected
        // only one shell.  We can do this one geometricly
        // if we have to.
      if ( !known_shell )
      {
        known_list.append( surf );
        continue;
      }
      
        // If we got this far, then the surface had at least
        // one curve that intersected only one of the shells.
        // We know it goes in that shell.
      did_some = true;
      PartitionCoSurf* cosurf1 = surf->find_first((PartitionShell*)0);
      PartitionCoSurf* cosurf2 = surf->find_next(cosurf1);
      known_shell->add(cosurf1);
      known_shell->add(cosurf2);
    
    } // while(surf_stack)  -- the inside one
    
      // Unknown_list always goes back in surf_stack to
      // try again.
    surf_stack += unknown_list;
    
      // If we did some surfaces, then put the rest back
      // in surf_stack and try again.  If they intersect
      // one of the surfaces we did place in this iteration,
      // we can avoid needing to do geometric checks.
    if ( did_some )
    {
      surf_stack += known_list;
      continue;
    }
    
      // If known_list is empty, somethings wrong (we're
      // going to loop forever.)  Abort the loop and try
      // to recover as best we can.
    if( !known_list.size() )
      break;
    
      // choose a single surface in do a geometric comparison
      // for, and put the rest back in surf_stack
    surf = known_list.pop();
    surf_stack += known_list;
    
    bool in_shell = false;
    if ( ! inside_shell( shell2, surf, in_shell ) )
    {
        // if inside_shell failed, abort.
      surf_stack.append(surf);
      break;
    }
    
    PartitionShell* shell = in_shell ? shell2 : shell1;
    shell->add(surf->find_first((PartitionShell*)0));
    shell->add(surf->find_first((PartitionShell*)0));
  }
  
    // something went wrong
  if( surf_stack.size() )
  {
    PRINT_ERROR("Internal error splitting volume at %s:%d\n"
                "Topology may be invalid.  Please report this.\n",
                __FILE__, __LINE__);
    while( surf_stack.size() ) 
    {
      PartitionSurface* surf = surf_stack.pop();
      cosurf = 0;
      while( (cosurf = surf->next_co_surface(cosurf)) )
        if( !cosurf->get_shell() )
          shell1->add(cosurf);
    }
  }
}


CubitStatus PartitionEngine::inside_shell( PartitionShell* const shell,
                                           PartitionSurface* const surf,
                                           bool& result )
{
    // Find the curve and coedge at which the nonmanifold surface
    // intersects the shells
  PartitionLoop* loop = 0;
  PartitionCoEdge *nonman_coedge = 0;
  while ( !nonman_coedge && (loop = surf->next_loop(loop)) )
  {
    PartitionCoEdge* loop_coedge = loop->first_coedge();
    do 
    {
        // Check if this curve is the curve of intersection
        // Iterate through curve coedges.
      PartitionCurve* curve = loop_coedge->get_curve();
      PartitionCoEdge* curve_coedge = 0;
      while ( (curve_coedge = curve->next_coedge(curve_coedge)) )
      {
        
        PartitionSurface* coedge_surf = curve_coedge->get_loop()->get_surface();
        if( coedge_surf->find_first(shell) )
        {
          nonman_coedge = curve_coedge;
          break;
        }
      }
    
      loop_coedge = loop->next_coedge(loop_coedge);
    } while( !nonman_coedge && loop_coedge != loop->first_coedge() );
  }
  
  if ( !nonman_coedge ) // bad input!
    return CUBIT_FAILURE;

    // There must exist two surfaces in the shell that are manifold 
    // in the shell and that are adjacent to the curve
  PartitionCurve* common_curve = nonman_coedge->get_curve();
  PartitionCoSurf *cosurf1 = 0, *cosurf2 = 0, *cosurf;
  PartitionCoEdge *coedge1 = 0, *coedge2 = 0, *coedge = 0;
  while ( (coedge = common_curve->next_coedge(coedge)) )
  {
    PartitionSurface* surf = coedge->get_loop()->get_surface();
    if ( (cosurf = surf->find_first(shell)) && !surf->find_next(cosurf) )
    {
      if( coedge1 ) {
        coedge2 = coedge;
        cosurf2 = cosurf;
      } else {
        coedge1 = coedge;
        cosurf1 = cosurf;
      }
    }
  }
  if ( !coedge1 || !coedge2 )
    return CUBIT_FAILURE;
  
    // Evaluate normals at midpoint of curve
  CubitVector base, tangent, point;
  double u = (common_curve->start_param()+common_curve->end_param())/2.0;
  common_curve->position_from_u( u, base );
  common_curve->closest_point( base, point, &tangent );
  tangent.normalize();
  
  CubitVector normal1, normal2, normal;
  surf->closest_point( base, 0, &normal );
  cosurf1->get_surface()->closest_point( base, 0, &normal1 );
  cosurf2->get_surface()->closest_point( base, 0, &normal2 );
  
    // Try to handle tangencies
  bool fix1 = (normal1 * normal).length_squared() < CUBIT_RESABS;
  bool fix2 = (normal2 * normal).length_squared() < CUBIT_RESABS;
  if ( fix1 || fix2 )
  {
    CubitVector dir = tangent * normal;
    double len = dir.length();
    assert(len > GEOMETRY_RESABS);
    dir /= len;
    if ( nonman_coedge->sense() == CUBIT_FORWARD )
      dir = -dir;
      
    CubitVector diag = surf->bounding_box().diagonal();
    if ( fix1 )
    {
      CubitVector d = cosurf1->get_surface()->bounding_box().diagonal();
      if ( diag.x() < d.x() ) diag.x(d.x());
      if ( diag.y() < d.y() ) diag.y(d.y());
      if ( diag.z() < d.z() ) diag.z(d.z());
    }
    if( fix2 )
    {
      CubitVector d = cosurf2->get_surface()->bounding_box().diagonal();
      if ( diag.x() < d.x() ) diag.x(d.x());
      if ( diag.y() < d.y() ) diag.y(d.y());
      if ( diag.z() < d.z() ) diag.z(d.z());
    }
    
    double step = 1e-3 * fabs( diag % dir );
    if ( step < 2*GEOMETRY_RESABS )
      step = 2*GEOMETRY_RESABS;
    
    for ( int i = 0; i < 1000; i++ )
    {
      point = base + i * step * dir;
      surf->closest_point( base, 0, &normal );
      if( fix1 )
        cosurf1->get_surface()->closest_point( base, 0, &normal1 );
      if( fix2 )
        cosurf2->get_surface()->closest_point( base, 0, &normal2 );
      
      bool done1 = !fix1 || (normal1 * normal).length_squared() > CUBIT_RESABS;
      bool done2 = !fix2 || (normal2 * normal).length_squared() > CUBIT_RESABS;
      if ( done1 && done2 )
      {
        fix1 = fix2 = false;
        break;
      }
    }
  }
  
  if ( fix1 || fix2 )
  {
    PRINT_ERROR("Failed to adjust for surface tangencies.\n"
                "This is a BUG.  %s:%d\n", __FILE__, __LINE__ );
    return CUBIT_FAILURE;
  }
  
  if ( nonman_coedge->sense() == CUBIT_FORWARD )
    normal = -normal;
  if ( cosurf1->sense() == coedge1->sense() )
    normal1 = -normal1;
  if ( cosurf2->sense() == coedge2->sense() )
    normal2 = -normal2;
  
  result = tangent.vector_angle(normal1, normal ) <=
           tangent.vector_angle(normal1, normal2);
  return CUBIT_SUCCESS;
}

  




//-------------------------------------------------------------------------
// Purpose       : Remove a partition surface previously inserted in a
//                 lump via insert_surface(..).
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 01/13/03
//-------------------------------------------------------------------------
Lump* PartitionEngine::remove_surface( PartitionSurface* surface )
{
    // If surface belongs to non-partition volumes, we cannot proceed.
  DLIList<TopologyBridge*> shell_bridges;
  surface->get_parents_virt( shell_bridges );
  while( shell_bridges.size() )
    if( ! dynamic_cast<PartitionShell*>(shell_bridges.pop()) )
      return 0;
  
  PartitionCoSurf* cosurf1 = surface->next_co_surface(0);
  if( ! cosurf1 )
    return 0; // doesn't belong to any PartitionLumps;
    
  PartitionCoSurf* cosurf2 = surface->next_co_surface(cosurf1);
  if( !cosurf2 || surface->next_co_surface(cosurf2) )
    return 0; // must have exactly two parent cosurfaces.
    
  PartitionShell* shell1 = cosurf1->get_shell();
  PartitionShell* shell2 = cosurf2->get_shell();
  PartitionLump* lump1 = shell1->get_lump();
  PartitionLump* lump2 = shell2->get_lump();
  
  if( surface->partitioned_entity() != lump1->partitioned_entity() )
    return 0; // shouldn't have gotten this far, but just in case ...
  
    // remove and destroy surface
  shell1->remove( cosurf1 );
  shell2->remove( cosurf2 );
  surface->remove(cosurf1);
  surface->remove(cosurf2);
  delete cosurf1;
  delete cosurf2;
  destroy_surface(surface);
    
  
    // combine volumes?
  if( shell1 != shell2 )
  {
    assert(lump1 != lump2);
    while( PartitionCoSurf* cosurf = shell2->next_co_surface(0) )
    {
      shell2->remove(cosurf);
      shell1->add(cosurf);
    }
    lump2->remove(shell2);
    delete shell2;
    
    if( !shell1->next_co_surface(0) )
    {
      lump1->remove(shell1);
      delete shell1;
    }
    
    while( PartitionShell* shell = lump2->next_shell(0) )
    {
      lump2->remove(shell);
      lump1->add(shell);
    }
    delete lump2;
  }
    // one-surface shell
  else if( !shell1->next_co_surface(0) )
  {
    lump1->remove( shell1 );
    delete shell1;
  }
    // possibly split shells
  else 
  {
    PartitionShell* shell = split_shell( shell1 );
    if( shell )
      lump1->add(shell);
  }
    
  if( !lump1->sub_entity_set().has_lower_order() )
    return restore_lump( lump1 );
  else
    return lump1;
}


//-------------------------------------------------------------------------
// Purpose       : tear down and free a surface
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 01/14/03
//-------------------------------------------------------------------------
CubitStatus PartitionEngine::destroy_surface( PartitionSurface* surface )
{
  while( PartitionLoop* dead_loop = surface->next_loop(0) )
  {
    surface->remove(dead_loop);
    while( PartitionCoEdge* dead_coe = dead_loop->first_coedge() )
    {
      PartitionCurve* dead_curve = dead_coe->get_curve();
      dead_loop->remove(dead_coe);
      dead_curve->remove(dead_coe);
      delete dead_coe;
      
      if( dead_curve->partitioned_entity() == surface->partitioned_entity() )
      {
        if( dead_curve->next_coedge(0) == 0 )
        {
          PartitionPoint* start_pt = dead_curve->start_point();
          PartitionPoint* end_pt = dead_curve->end_point();
          delete dead_curve;
          if( start_pt->partitioned_entity() == surface->partitioned_entity() &&
              start_pt->next_curve() == 0 )
            delete start_pt;
          if( end_pt->partitioned_entity() == surface->partitioned_entity() &&
              end_pt->next_curve() == 0 )
            delete end_pt;
        }
      }
      else if( SubCurve* subcurve = dynamic_cast<SubCurve*>(dead_curve) )
      {
        restore_curve(subcurve);
      }
    }
    delete dead_loop;
  }
  delete surface;
  return CUBIT_SUCCESS;
}


//-------------------------------------------------------------------------
// Purpose       : Retreive entity from UID attrib data.
//
// Special Notes : Used when restoring saved geometry.
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 02/23/03
//-------------------------------------------------------------------------
PartitionEntity* PartitionEngine::entity_from_id( int set_id,
                                                  int entity_id,
                                                  SubEntitySet& default_set )
{
  if( set_id == 0 )
    return default_set.entity_from_id( entity_id );
  
  SubEntitySet* set = get_from_id_map(set_id);
  if( set )
    return set->entity_from_id( entity_id );
  
  return 0;
}
 
//-------------------------------------------------------------------------
// Purpose       : Restore curve partitions
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 02/23/03
//-------------------------------------------------------------------------
CubitStatus PartitionEngine::restore_from_attrib( Curve* curve )
{
  int i, id, dim, max_id = 0;
  DLIList<CubitSimpleAttrib*> attribs, point_attribs, curve_attribs;
  DLIList<CubitVector*> empty;
  DLIList<int> junk, point_conn;
  DLIList<PartitionPoint*> new_points;
  CubitStatus result = CUBIT_FAILURE;
  
  {
      // Get list of partition attributes.
    curve->get_simple_attribute(PARTITION_GEOM_ATTRIB_NAME,attribs);
    if( attribs.size() == 0 )
      return CUBIT_SUCCESS;

      // Group attributes into those defining partition points and
      // those defining the subcurves.  Also find max ID of all 
      // partition geometry on this curve
    attribs.reset();
    for( i = attribs.size(); i--; )
    {
      CubitSimpleAttrib* attrib = attribs.get_and_step();
      dim = SubEntitySet::get_geom_dimension(*attrib);
      if( dim == 0 )
        point_attribs.append(attrib);
      else if( dim == 1 )
        curve_attribs.append(attrib);
      else
        goto CLEANUP_CURVE_FROM_ATTRIB;

      id = SubEntitySet::get_geom_id( *attrib );
      if( id > max_id )
        max_id = id;
    }

      // Create first curve (duplicate of input curve)
    PartitionCurve* first_curve = replace_curve( curve );
    if( ! first_curve )
      goto CLEANUP_CURVE_FROM_ATTRIB;

      // Set up SubEntitySet such that new entities get IDs
      // larger than max_id as they are created.
    first_curve->sub_entity_set().renumerate( max_id + 1, false );

      // Create all the partition points
    point_attribs.reset();
    for( i = point_attribs.size(); i--; )
    {
      CubitSimpleAttrib* attrib = point_attribs.get_and_step();
      PartitionPoint* point = new PartitionPoint( *attrib, first_curve );
      if( !point )
        goto CLEANUP_CURVE_FROM_ATTRIB;

      curve->remove_simple_attribute_virt(attrib);
    }

      // Now partition the curve and assign the ID from each
      // curve attribute the the appropriate curve partition.
      //
      // We're going to do this in order along the curve, such
      // that the first subcurve will begin at the start_point() of
      // the original curve.
    PartitionPoint* prev_point = first_curve->start_point();
    PartitionCurve* next_curve = first_curve;
    while( curve_attribs.size() )
    {
        // Search for attribute corresponding to curve
        // beginning at prev_point
      PartitionPoint *start = 0, *end = 0;
      CubitSimpleAttrib* curve_attrib = 0;
      curve_attribs.last();
      for( i = curve_attribs.size(); i--; )
      {
          // Get data from attribute
        CubitSimpleAttrib* attrib = curve_attribs.step_and_get();
        point_conn.clean_out();
        if( !SubEntitySet::read_geometry(id,dim,empty,junk,point_conn,junk,*attrib) 
            || empty.size() || junk.size() || dim != 1 || point_conn.size() != 4)
          goto CLEANUP_CURVE_FROM_ATTRIB;

          // Get start point from IDs
        point_conn.reset();
        int sid = point_conn.get_and_step();
        int eid = point_conn.get_and_step();
        PartitionEntity* ent = entity_from_id(sid, eid, first_curve->sub_entity_set());
        start = dynamic_cast<PartitionPoint*>(ent);
        
        if( start == prev_point )
        {
            // Found the attrib we want. 
            // Get the end point and break the loop.
          sid = point_conn.get_and_step();
          eid = point_conn.get_and_step();
          ent = entity_from_id(sid, eid, first_curve->sub_entity_set());
          end = dynamic_cast<PartitionPoint*>(ent);

          curve_attribs.extract();
          curve_attrib = attrib;
          break;
        }
      }

        // Make sure we found an attrib to do
      if( !curve_attrib || !start || !end )
        goto CLEANUP_CURVE_FROM_ATTRIB;

        // Unless this is the last partition (no more splitting
        // required), split the curve at the end point.  Make
        // next_curve point to the (unpartitioned) remainder of
        // the curve.
      PartitionCurve* new_curve = next_curve;
      if( next_curve->end_point() != end )
      {
        next_curve = insert_point( new_curve, end );
        if( !next_curve )
          goto CLEANUP_CURVE_FROM_ATTRIB;
      }

        // Set the new curve's ID to whatever it is stored
        // as in the attribute.
      new_curve->sub_entity_set().set_id( new_curve, 
        SubEntitySet::get_geom_id( *curve_attrib ) );

        // Remove processed attributes.
      curve->remove_simple_attribute_virt( curve_attrib );
      
        // Next iteration
      prev_point = end;
    }

      // Compress SubEntitySet ID space
    next_curve->sub_entity_set().renumerate( max_id + 1, true );
    result = CUBIT_SUCCESS;
  }
  
CLEANUP_CURVE_FROM_ATTRIB:

    // free attribs
  while( attribs.size() )
    delete attribs.pop();
    
  for( i = new_points.size(); i--; )
  {
    PartitionPoint* pt = new_points.get_and_step();
    if( !pt->next_curve(0) )
      delete pt;
  }
  
  if (!result)
    PRINT_ERROR("Error restoring curve partitions from attributes.\n");
    
  return result;
}  


//-------------------------------------------------------------------------
// Purpose       : Restore surface partitions
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 02/23/03
//-------------------------------------------------------------------------
CubitStatus PartitionEngine::restore_from_attrib( Surface* surf )
{
  int i, j, k, id, dim, max_id = 0;
  DLIList<CubitSimpleAttrib*> attribs;
  DLIList<CubitFacetEdgeData*> polyline_edges;
  DLIList<CubitPointData*> polyline_pts;
  DLIList<CubitFacetData*> surf_facets;
  DLIList<CubitVector*> positions;
  DLIList<PartitionEntity*> new_geom;
  DLIList<PartitionLoop*> loops;
  DLIList<int> facets(0), connectivity, pt_owners;
  CubitStatus result = CUBIT_FAILURE;
  SubSurface* first_surf = 0;
  
  
  {
      // Get partition attributes from entity.  Return if there aren't any.
    surf->get_simple_attribute(PARTITION_GEOM_ATTRIB_NAME,attribs);
    if( attribs.size() == 0 )
      return CUBIT_SUCCESS;

    PRINT_DEBUG_86("Reading partition geometry on %s %p\n", 
      fix_type_name( typeid(*surf).name() ), surf );
      
      // Find max ID
    attribs.reset();
    for( i = attribs.size(); i--; )
    {
      id = SubEntitySet::get_geom_id( *attribs.get_and_step() );
      if( id > max_id )
        max_id = id;
    }

      // Create first partition surface (initially a duplicate of
      // the input surface)
    first_surf = replace_surface( surf );
    if( ! first_surf )
    {  
      PRINT_ERROR("Internal Error: PartitionEngine::replace_surface() failed.\n");
      goto CLEANUP_SURF_FROM_ATTRIB;
    }
    
      // Update SubEntitySet such that new entities get assigned
      // IDs greater than max_id when they are created.  This
      // leaves the range 0->max_id free so that we can reassign
      // the saved IDs after the entities are created.
    first_surf->sub_entity_set().renumerate( max_id + 1, false );

      // create all interior vertices
    attribs.reset();
    for( i = attribs.size(); i--; )
    {
        // Only want point attributes (dimension of zero)
      CubitSimpleAttrib* attrib = attribs.get();
      if( SubEntitySet::get_geom_dimension(*attrib) != 0 ) 
      {
        attribs.step();
        continue;
      }

        // Create the point at the location specified in the attrib.
      attribs.extract();
      surf->remove_simple_attribute_virt(attrib);
      PartitionPoint* new_pt = new PartitionPoint( *attrib, first_surf );
      CubitPointData* data = project_to_surface( first_surf, new_pt->coordinates() );
      new_pt->facet_point(data);
      new_geom.append(new_pt);
      delete attrib;
      PRINT_DEBUG_86("  Created Point %p (%d in subset) at (%f,%f,%f)\n", 
        new_pt, new_pt->sub_entity_set().get_id(new_pt), 
        new_pt->coordinates().x(), new_pt->coordinates().y(), 
        new_pt->coordinates().z() );
      if (DEBUG_FLAG(86))
      {
        GfxDebug::draw_point(new_pt->coordinates(), new_pt->sub_entity_set().get_id(new_pt));
        GfxDebug::flush();
      }
    }

      // create all interior curves
    attribs.reset();
    for( i = attribs.size(); i--; )
    {
        // Only want curve attributes (dimension of one)
      CubitSimpleAttrib* attrib = attribs.get();
      if( SubEntitySet::get_geom_dimension(*attrib) != 1 ) 
      {
        attribs.step();
        continue;
      }

        // Remove this attrib from the real geom, as we are
        // handling it
      attribs.extract();
      surf->remove_simple_attribute_virt(attrib);
      
      if ( SubEntitySet::get_segment_count(*attrib) == 0 )
      {
        PartPTCurve* pt_curve = PartPTCurve::construct(attrib, first_surf);
        if (!pt_curve)
        {
          PRINT_ERROR("point-curve construction failed -- corrupt data?\n");
          goto CLEANUP_SURF_FROM_ATTRIB;
        }
        new_geom.append( pt_curve );
        PRINT_DEBUG_86("  Created point-curve %p (%d in subset) with\n"
                       "    point %p (%d in subset) at (%f,%f,%f)\n", 
          pt_curve, 
          pt_curve->sub_entity_set().get_id(pt_curve), 
          pt_curve->start_point(),
          pt_curve->start_point()->sub_entity_set().get_id(pt_curve->start_point()),
          pt_curve->start_point()->coordinates().x(),
          pt_curve->start_point()->coordinates().y(), 
          pt_curve->start_point()->coordinates().z() );

        CubitPointData* pt = project_to_surface( first_surf, pt_curve->start_point()->coordinates() );
        if ( !pt )
        {
          PRINT_ERROR("projection of position onto surface facets failed.\n");
          goto CLEANUP_SURF_FROM_ATTRIB;
        }
          
        pt_curve->start_point()->facet_point( pt );
        continue;
      }
      
        // Create the curve
      SegmentedCurve* new_curve = SegmentedCurve::construct( attrib, first_surf );
      if( !new_curve ) 
      {
        PRINT_ERROR("segmented curve construction failed -- corrupt data?\n");
        goto CLEANUP_SURF_FROM_ATTRIB;
      }
      new_geom.append(new_curve);
      PRINT_DEBUG_86("  Created polyline curve %p (%d in subset) with %d segments.\n"
                     "    start point %p (%d in subset) at (%f,%f,%f)\n" 
                     "      end point %p (%d in subset) at (%f,%f,%f)\n", 
          new_curve, 
          new_curve->sub_entity_set().get_id(new_curve),
          new_curve->point_count() - 1, 
          new_curve->start_point(),
          new_curve->start_point()->sub_entity_set().get_id(new_curve->start_point()),
          new_curve->start_point()->coordinates().x(),
          new_curve->start_point()->coordinates().y(), 
          new_curve->start_point()->coordinates().z(),
          new_curve->end_point(),
          new_curve->end_point()->sub_entity_set().get_id(new_curve->end_point()),
          new_curve->end_point()->coordinates().x(),
          new_curve->end_point()->coordinates().y(), 
          new_curve->end_point()->coordinates().z() );
 
        // Get curve segments
      positions.clean_out();
      polyline_edges.clean_out();
      new_curve->get_segments( positions );
      
      if (DEBUG_FLAG(86))
      {
        GPoint* array = new GPoint[positions.size()];
        positions.reset();
        for (int dd = 0; dd < positions.size(); dd++ )
        {
          array[dd].x = (float)positions.next(dd)->x(); 
          array[dd].y = (float)positions.next(dd)->y();
          array[dd].z = (float)positions.next(dd)->z();
        }
        
        GfxDebug::draw_polyline(array, positions.size(),
                                new_curve->sub_entity_set().get_id(new_curve));
        GfxDebug::flush();
      }
      
        // Make sure segment end points are on curve end points
      positions.last();
      *positions.get() = new_curve->end_point()->coordinates();
      positions.reset();
      *positions.get() = new_curve->start_point()->coordinates();
      
        // Project curve into surface facets
      surf_facets.clean_out();
      first_surf->get_facet_data( surf_facets );
      CubitStatus s = project_to_surface( surf_facets, positions, 
                              polyline_edges, polyline_pts );
      while( positions.size() ) 
      {
        CubitVector* position = positions.pop();
        if(s)
        {
   //       assert(polyline_pts.size());
          CubitPointData* point = polyline_pts.pop();
          if (point && point->check_inverted_facets(*position))
            point->set(*position);
        }
        delete position;
      }
      delete attrib;
      if( !s || !polyline_edges.size() ) 
      {
        PRINT_ERROR("Imprint of polyline onto surface facets failed.\n");
        goto CLEANUP_SURF_FROM_ATTRIB;
      }
      
      
        // move facet points onto surface
      for ( j = 0; j < 2; j++ )
        for ( k = polyline_edges.size(); k--; )
          if ( CubitFacetEdgeData* edge = polyline_edges.get_and_step() )
            edge->point(j)->marked(1);
      
      for ( j = polyline_edges.size(); j--; )
      {
        CubitFacetEdgeData* edge = polyline_edges.get_and_step();
        if ( !edge ) continue;
        for ( k = 0; k < 2; k++ )
        {
          CubitPoint* pt = edge->point(k);
          if ( pt->marked() ) 
          {
            pt->marked(0);
            first_surf->relax_to_geometry(pt);
          }
        }
      }
        
      
        // Associate facet data with geometry
        
      CubitPoint *last_pt, *first_pt;
      polyline_edges.last();
      assert(!polyline_edges.get());
      polyline_edges.pop();
      
      if( polyline_edges.size() == 1 )
      {
        last_pt = polyline_edges.get()->point(1);
        first_pt = polyline_edges.get()->point(0);
      }
      else
      {
        CubitFacetEdgeData *edge, *neighbor;
        
        polyline_edges.last();
        edge = polyline_edges.get();
        neighbor = polyline_edges.prev();
        if( !edge || !neighbor )
        {
          PRINT_ERROR("Bad surface facets.\n");
          goto CLEANUP_SURF_FROM_ATTRIB;
        }
        last_pt = edge->shared_point(neighbor);
        last_pt = edge->other_point(last_pt);

        polyline_edges.reset();
        edge = polyline_edges.get();
        neighbor = polyline_edges.next();
        if( !edge || !neighbor )
        {
          PRINT_ERROR("Bad surface facets.\n");
          goto CLEANUP_SURF_FROM_ATTRIB;
        }
        first_pt = edge->shared_point(neighbor);
        first_pt = edge->other_point(first_pt);
      }

      PartitionEntity* owner1 = TDVGFacetOwner::get( first_pt );
      PartitionEntity* owner2 = TDVGFacetOwner::get( last_pt );

      if ( new_curve->start_point() == new_curve->end_point() )
      {
        if ( owner1 == owner2 && owner1 == new_curve->start_point() )
        {
          first_pt = polyline_edges.get()->other_point( first_pt );
          polyline_edges.last();
          last_pt = polyline_edges.get()->other_point( last_pt );
          double u1 = new_curve->u_from_position( first_pt->coordinates() );
          double u2 = new_curve->u_from_position( last_pt->coordinates() );
          if ( u2 < u1 )
            polyline_edges.reverse();
          polyline_edges.reset();
        }
        else
        {
          s = CUBIT_FAILURE;
        }
      }
      else if ( owner1 == new_curve->end_point() &&
                owner2 == new_curve->start_point() )
      {
        polyline_edges.reverse();
        polyline_edges.reset();
      }
      else if( owner1 != new_curve->start_point() ||
               owner2 != new_curve->end_point() )
      {
        s = CUBIT_FAILURE;
      }
     
      if (!s) 
      {
        PRINT_ERROR("Tolerance problems when restoring surface partitions.\n");
        goto CLEANUP_SURF_FROM_ATTRIB;
      }
      
      new_curve->set_facet_data( polyline_edges );
    }

      // while there are more surface partitions to create
    PartitionEntity* ent = 0;
    while( attribs.size() )
    {
        // read attribute
      CubitSimpleAttrib* surf_attr = attribs.pop();
      connectivity.clean_out();
      CubitStatus s = SubEntitySet::
        read_geometry( id, dim, positions, facets, connectivity, pt_owners, *surf_attr );
      delete surf_attr;
      
      if( !s || dim != 2 || positions.size() || facets.size() )
      {
        PRINT_ERROR("Corrupt/inconsistent subsurface data.\n");
        goto CLEANUP_SURF_FROM_ATTRIB;
      }
      
      PRINT_DEBUG_86("  Constructing SubSurface (%d in subset)\n", id);

        // construct surface loops
      int j = connectivity.size();
      PartitionCoEdge* curve_coedge = 0; // save one non-point-curve coedge
      while ( j > 0 )
      {
        int loop = connectivity.get_and_step();
        j--;

        PartitionLoop* new_loop = new PartitionLoop();
        loops.append( new_loop );
        PartitionCoEdge* prev_coedge = 0;
        PRINT_DEBUG_86("   Loop: ");

        while( loop-- )
        {
          int set_id = connectivity.get_and_step();
          int ent_id = connectivity.get_and_step();
          int sense = ent_id < 0 ? -1 : 1;
          ent_id *= sense;
          j -= 2;

          ent = entity_from_id( set_id, ent_id, first_surf->sub_entity_set() );
          PartitionCurve* curve = dynamic_cast<PartitionCurve*>(ent);
          PRINT_DEBUG_86(" %p,%d,%d%c", curve, set_id, ent_id*sense, loop?',':'\n');
          if(!curve)
          {
            PRINT_ERROR("Nonexistant curve specified in saved connectivity.\n");
            goto CLEANUP_SURF_FROM_ATTRIB;
          }

          CubitSense cubit_sense = sense == -1 ? CUBIT_REVERSED : CUBIT_FORWARD;
          PartitionCoEdge* coedge = 0;
          while ( (coedge = curve->next_coedge(coedge)) )
          {
            if ( coedge->sense() == cubit_sense &&
                 coedge->get_loop() &&
                 coedge->get_loop()->get_surface() == first_surf )
            {
              coedge->get_loop()->remove( coedge );
              break;
            }
          }
          if ( !coedge )
          {
            coedge = new PartitionCoEdge( first_surf, cubit_sense );
            curve->add(coedge);
          }
          
          if (!curve_coedge && curve->geometry_type() != POINT_CURVE_TYPE)
            curve_coedge = coedge;
            
          new_loop->insert_after( coedge, prev_coedge );
          prev_coedge = coedge;
        }
      }
      
      if (!curve_coedge)
      {
        PRINT_ERROR("Nonexistant co-edge specified in saved connectivity.\n");
        goto CLEANUP_SURF_FROM_ATTRIB;
      }
      
      PartitionSurface* new_surf = split_surface( first_surf, curve_coedge );
      if (new_surf == first_surf && attribs.size() > 0)
      {
        PRINT_ERROR("Surface splitting failed -- bad topology?.\n");
        goto CLEANUP_SURF_FROM_ATTRIB;
      }
            
      new_surf->sub_entity_set().set_id( new_surf, id );
      while ( loops.size() )
        new_surf->add( loops.pop() );
    }

    result = CUBIT_SUCCESS;
  }
  
CLEANUP_SURF_FROM_ATTRIB:

  while( attribs.size() )
    delete attribs.pop();
    
    // clean up any unused loops from failed surface creation
  while ( loops.size() )
  {
    PartitionLoop* loop = loops.pop();
    while ( loop->first_coedge() )
    {
      PartitionCoEdge* coedge = loop->first_coedge();
      loop->remove(coedge);
      if (coedge->get_curve())
        coedge->get_curve()->remove(coedge);
      delete coedge;
    }
    delete loop;
  }  
    
    // first_surf will contain dead loops - clean them up
  PartitionLoop* surf_loop = 0;
  while ( (surf_loop = first_surf->next_loop( surf_loop ) ) != NULL )
    if ( !surf_loop->first_coedge() )
      loops.append( surf_loop );
  while ( loops.size() )
  {
    surf_loop = loops.pop();
    first_surf->remove( surf_loop );
    delete surf_loop;
  }
    
  new_geom.reset();
  while( new_geom.size() )
  {
    PartitionEntity* ent = new_geom.pop();
    if( PartitionPoint* pt = dynamic_cast<PartitionPoint*>(ent) )
    {
      if( !pt->next_curve(0) )
        delete pt;
    }
    else if( PartitionCurve* curve = dynamic_cast<PartitionCurve*>(ent) )
    {
      if( !curve->next_coedge(0) )
        delete curve;
    }
  }
  
  first_surf->sub_entity_set().renumerate( max_id + 1, true );
  
  if (!result)
    PRINT_ERROR("Error restoring surface partitions from attributes.\n");
    
  return result;
}

//-------------------------------------------------------------------------
// Purpose       : Restore volume partitions from an attribute
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 01/22/03
//-------------------------------------------------------------------------
CubitStatus PartitionEngine::restore_from_attrib( Lump* lump )
{
  int i, id, dim, max_id = 0;
  DLIList<CubitSimpleAttrib*> attribs;
  DLIList<CubitVector*> positions;
  DLIList<int> facets, connectivity, pt_owners;
  DLIList<PartitionEntity*> new_geom;
  CubitStatus result = CUBIT_FAILURE;
  PartitionLump* first_lump = 0;
  PartitionLump* new_lump = 0;
  
  { // pseudo-'try'-block
  
      // get partition attributes from entity.  Return if none.
    lump->get_simple_attribute(PARTITION_GEOM_ATTRIB_NAME,attribs);
    if( attribs.size() == 0 )
      return CUBIT_SUCCESS;

      // find max ID
    attribs.reset();
    for( i = attribs.size(); i--; )
    {
      id = SubEntitySet::get_geom_id( *attribs.get_and_step() );
      if( id > max_id )
        max_id = id;
    }

      // Create first parition lump (initially a duplicate of
      // the input lump)
    first_lump = replace_lump( lump );
    if( ! first_lump )
      goto CLEANUP_LUMP_FROM_ATTRIB;
    
      // Update SubEntitySet such that new entities get assigned
      // IDs greater than max_id when they are created.  This
      // leaves the range 0->max_id free so that we can reassign
      // the saved IDs after the entities are created.
    first_lump->sub_entity_set().renumerate( max_id + 1, false );

      // create all interior vertices
    attribs.reset();
    for( i = attribs.size(); i--; )
    {
        // skip any non-point geometry
      CubitSimpleAttrib* attrib = attribs.get();
      if( SubEntitySet::get_geom_dimension(*attrib) != 0 ) 
      {
        attribs.step();
        continue;
      }

        // create the point as specified in the attrib
      attribs.extract();
      lump->remove_simple_attribute_virt(attrib);
      PartitionPoint* new_pt = new PartitionPoint( *attrib, first_lump );
      new_geom.append(new_pt);
      delete attrib;
    }

      // create all interior curves
    attribs.reset();
    for( i = attribs.size(); i--; )
    {
        // skip attribs for anything but curves
      CubitSimpleAttrib* attrib = attribs.get();
      if( SubEntitySet::get_geom_dimension(*attrib) != 1 ) 
      {
        attribs.step();
        continue;
      }

        // create curve from attrib
      attribs.extract();
      lump->remove_simple_attribute_virt(attrib);
      SegmentedCurve* new_curve = SegmentedCurve::construct( attrib, first_lump );
      if( !new_curve )
        goto CLEANUP_LUMP_FROM_ATTRIB;

      new_geom.append(new_curve);
      delete attrib;
    }

      // create all interior surfaces
    attribs.reset();
    for( i = attribs.size(); i--; )
    {
        // skip non-surface attributes
      CubitSimpleAttrib* attrib = attribs.get();
      if( SubEntitySet::get_geom_dimension(*attrib) != 2 ) 
      {
        attribs.step();
        continue;
      }

        // create surface from attribute
      attribs.extract();
      lump->remove_simple_attribute_virt(attrib);
      PartitionSurface* new_surf = PartitionSurface::construct( *attrib, first_lump );
      if( !new_surf )
        goto CLEANUP_LUMP_FROM_ATTRIB;

      new_geom.append(new_surf);
      delete attrib;
    }


      // while there are more volume partitions to create
    PartitionEntity* ent = 0;
    while( attribs.size() )
    {
        // Read data from attribute
      CubitSimpleAttrib* attr = attribs.pop();
      connectivity.clean_out();
      CubitStatus s = SubEntitySet::
        read_geometry( id, dim, positions, facets, connectivity, pt_owners, *attr );
      if( !s || dim != 3 || positions.size() || facets.size() )
      {
        new_lump = 0;
        goto CLEANUP_LUMP_FROM_ATTRIB;
      }
      
      new_lump = new PartitionLump( first_lump );

        // Get list of bounding co-surfaces from attrib data
      int j = connectivity.size();
      while( j > 0 )
      {
        int shell = connectivity.get_and_step();
        j--;
        PartitionShell* new_shell = new PartitionShell();
        new_lump->add(new_shell);
        while( shell-- )
        {
          int set_id = connectivity.get_and_step();
          int ent_id = connectivity.get_and_step();
          j -= 2;
          
          CubitSense sense = CUBIT_FORWARD;
          if ( ent_id < 0 )
          {
            ent_id = -ent_id;
            sense = CUBIT_REVERSED;
          }

          ent = entity_from_id( set_id, ent_id, first_lump->sub_entity_set() );
          PartitionSurface* surf = dynamic_cast<PartitionSurface*>(ent);
          if(!surf)
            goto CLEANUP_LUMP_FROM_ATTRIB;

          PartitionCoSurf* new_cosurf = new PartitionCoSurf(sense);
          surf->add(new_cosurf);
          new_shell->add(new_cosurf);
        }
      }

      new_lump->sub_entity_set().set_id( new_lump, id );
      new_lump = 0;  // we succeeded with this one, so don't try to 
                     // destroy it if we fail before creating the next.
    }
    
    result = CUBIT_SUCCESS;
  }
  
CLEANUP_LUMP_FROM_ATTRIB:

  if( new_lump ) // created a lump but didn't successfully complete it
  {
    PartitionShell* shell = 0;
    while ( shell )
    {
      new_lump->remove(shell);
      shell->remove_all_surfaces();
      PartitionShell* old_shell = shell;
      shell = new_lump->next_shell(shell);
      delete old_shell;
    }
  }
  
  DLIList<PartitionEntity*> lump_list;
  if( first_lump )
  {
      // Get list of other lumps created
    first_lump->sub_entity_set().get_sub_entities( lump_list );
    lump_list.move_to(first_lump);
    assert( lump_list.get() == first_lump );
    lump_list.extract();
    
      // Destroy first_lump (it is the union of all the partitions)
    DLIList<PartitionSurface*> surfaces;
    PartitionShell* shell = 0;
    while ( shell )
    {
      first_lump->remove(shell);
      shell->remove_all_surfaces( &surfaces );
      PartitionShell* old_shell = shell;
      shell = first_lump->next_shell(shell);
      delete old_shell;
      while( surfaces.size() )
      {
        SubSurface* surf = dynamic_cast<SubSurface*>(surfaces.pop());
        if ( surf && !surf->next_co_surface(0) && 
             !surf->sub_entity_set().has_lower_order() )
          restore_surface( surf );
      }
    }
    delete first_lump;
  }
  
    // deallocate any attribute objects still remaining
  while( attribs.size() )
    delete attribs.pop();
    
    // destroy any unused split geometry
  new_geom.reset();
  while( new_geom.size() )
  {
    PartitionEntity* ent = new_geom.pop();
    if( PartitionPoint* pt = dynamic_cast<PartitionPoint*>(ent) )
    {
      if( !pt->next_curve(0) )
        delete pt;
    }
    else if( PartitionCurve* curve = dynamic_cast<PartitionCurve*>(ent) )
    {
      if( !curve->next_coedge(0) )
        delete curve;
    }
    else if( PartitionSurface* surf = dynamic_cast<PartitionSurface*>(ent) )
    {
      if( !surf->next_co_surface(0) )
        destroy_surface( surf );
    }
  }
  
    // if there are any remaining partitions of the volume,
    // compress the ID space
  if ( lump_list.size() )
  {
    lump_list.get()->sub_entity_set().renumerate( max_id + 1, true );
  }
  
  if (!result)
    PRINT_ERROR("Error restoring volume partitions from attributes.\n");
    
  return result;
}  

      

//-------------------------------------------------------------------------
// Purpose       : Save partition geometry
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 02/23/03
//-------------------------------------------------------------------------
CubitStatus PartitionEngine::export_geometry( DLIList<TopologyBridge*>& list )
{
  int i;
  CubitStatus result = CUBIT_SUCCESS;

  if ( CGMApp::instance()->attrib_manager()->auto_update_flag(CA_PARTITION_VG) &&
       CGMApp::instance()->attrib_manager()->auto_write_flag(CA_PARTITION_VG) )
  {

      // Get all child topology to export

    DLIList<Curve*> curve_list;
    DLIList<Surface*> surface_list;
    DLIList<Lump*> lump_list;
    DLIList<BodySM*> body_list;
    DLIList<TopologyBridge*> bridge_list, temp_list, coe_curves(1), surf_curves;

    CAST_LIST( list, curve_list, Curve );
    CAST_LIST( list, surface_list, Surface );
    CAST_LIST( list, lump_list, Lump );
    CAST_LIST( list, body_list, BodySM );

    for( i = body_list.size(); i--; )
    {
      BodySM* body = body_list.get_and_step();
      body->get_children( bridge_list, true, layer() );
      while( bridge_list.size() )
        lump_list.append( dynamic_cast<Lump*>(bridge_list.pop()) );
    }

    for( i = lump_list.size(); i--; )
    {
      Lump* lump = lump_list.get_and_step();
      lump->get_children( bridge_list, true, layer() );
      while( bridge_list.size() )
      {
        bridge_list.pop()->get_children( temp_list, true, layer() );
        while( temp_list.size() )
          surface_list.append( dynamic_cast<Surface*>(temp_list.pop()) );
      }
    }

    for( i = surface_list.size(); i--; )
    {
      Surface* surf = surface_list.get_and_step();
      surf->get_children( bridge_list, true, layer() );
      while( bridge_list.size() )
      {
        bridge_list.pop()->get_children( temp_list, true, layer() );
        while( temp_list.size() )
        {
          temp_list.pop()->get_children_virt( coe_curves );
          assert(coe_curves.size() == 1);
          surf_curves.append(coe_curves.pop());
        }
      }
      surf_curves.uniquify_unordered();
      while( surf_curves.size() )
        curve_list.append( dynamic_cast<Curve*>(surf_curves.pop()) );
    }


      // save curve partitions

    DLIList<SubEntitySet*> set_list;
    for( i = curve_list.size(); i--; )
    {
      SubCurve* curve = dynamic_cast<SubCurve*>(curve_list.get_and_step());
      if( curve )
        set_list.append( &(curve->sub_entity_set()) );
    }
    set_list.uniquify_unordered();
    for( i = set_list.size(); i--; )
      if( ! set_list.get_and_step()->save_geometry() )
        result = CUBIT_FAILURE;


      // save surface partitions

    set_list.clean_out();
    for( i = surface_list.size(); i--; )
    {
      SubSurface* surf = dynamic_cast<SubSurface*>(surface_list.get_and_step());
      if( surf )
        set_list.append( &(surf->sub_entity_set()) );
    }
    set_list.uniquify_unordered();
    for( i = set_list.size(); i--; )
      if( ! set_list.get_and_step()->save_geometry() )
        result = CUBIT_FAILURE;


      // save volume partitions

    set_list.clean_out();
    for( i = lump_list.size(); i--; )
    {
      PartitionLump* lump = dynamic_cast<PartitionLump*>(lump_list.get_and_step());
      if( lump )
        set_list.append( &(lump->sub_entity_set()) );
    }
    set_list.uniquify_unordered();
    for( i = set_list.size(); i--; )
      if( ! set_list.get_and_step()->save_geometry() )
        result = CUBIT_FAILURE;
  }
    
    // replace partitions in passed list with the real
    // entities they partition
  
  PartitionEntity* p_ent = 0; 
  list.last() ;
  for ( i = list.size(); i--; )
  {
    if( (p_ent = dynamic_cast<PartitionEntity*>(list.step_and_get())) )
    {
      if( list.is_in_list(p_ent->partitioned_entity()) )
        list.change_to(0);
      else
        list.change_to(p_ent->partitioned_entity());
    }
  }
  list.remove_all_with_value(0);
  
  return result;
}


//-------------------------------------------------------------------------
// Purpose       : Restore form attributes
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 01/23/03
//-------------------------------------------------------------------------
CubitStatus PartitionEngine::import_geometry( DLIList<TopologyBridge*>& list )
{
  int i;
  CubitStatus result = CUBIT_SUCCESS;
  
  DLIList<Curve*> curve_list, temp_curves;
  DLIList<Surface*> surface_list, temp_surfaces;
  DLIList<Lump*> lump_list, temp_lumps;
  DLIList<BodySM*> body_list;
  
  CAST_LIST( list, curve_list, Curve );
  CAST_LIST( list, surface_list, Surface );
  CAST_LIST( list, lump_list, Lump );
  CAST_LIST( list, body_list, BodySM );
  
  for( i = body_list.size(); i--; )
  {
    temp_lumps.clean_out();
    body_list.get_and_step()->lumps( temp_lumps );
    lump_list += temp_lumps;
  }
  
  for( i = lump_list.size(); i--; )
  {
    temp_surfaces.clean_out();
    lump_list.get_and_step()->surfaces( temp_surfaces );
    surface_list += temp_surfaces;
  }
  
  for( i = surface_list.size(); i--; )
  {
    temp_curves.clean_out();
    surface_list.get_and_step()->curves( temp_curves );
    curve_list += temp_curves;
  }
  
  
  if ( CGMApp::instance()->attrib_manager()->auto_actuate_flag(CA_PARTITION_VG) &&
       CGMApp::instance()->attrib_manager()->auto_read_flag(CA_PARTITION_VG) )
  {
  
    curve_list.uniquify_unordered();
    for( i = curve_list.size(); i--; )
      if( ! restore_from_attrib( curve_list.get_and_step() ) )
        result = CUBIT_FAILURE;

    surface_list.uniquify_unordered();
    for( i = surface_list.size(); i--; )
      if( ! restore_from_attrib( surface_list.get_and_step() ) )
        result = CUBIT_FAILURE;

    lump_list.uniquify_unordered();
    for( i = lump_list.size(); i--; )
      if( ! restore_from_attrib( lump_list.get_and_step() ) )
        result = CUBIT_FAILURE;

      // update imported list
    list.last();
    DLIList<PartitionEntity*> entity_list;
    DLIList<TopologyBridge*> temp_list;
    for( i = list.size(); i--; )
    {
      TopologyBridge* bridge = list.step_and_get();
      SubEntitySet* set = dynamic_cast<SubEntitySet*>(bridge->owner());
      if( !set ) continue;

      entity_list.clean_out();
      temp_list.clean_out();
      set->get_sub_entities( entity_list );
      CAST_LIST( entity_list, temp_list, TopologyBridge );

      temp_list.reset();
      list.change_to( temp_list.get_and_step() );
      for( int j = temp_list.size(); j > 1; j-- )
        list.append( temp_list.get_and_step() );
    }
  }
  else
  {
    curve_list.uniquify_unordered();
    surface_list.uniquify_unordered();
    lump_list.uniquify_unordered();
    
    DLIList<TopologyBridge*> point_list(curve_list.size()*2), tmp_points(2);
    for ( i = curve_list.size(); i--; )
    {
      tmp_points.clean_out();
      curve_list.get()->get_children(tmp_points);
      point_list += tmp_points;
      SubEntitySet::strip_attributes( curve_list.get_and_step() );
    }
    point_list.uniquify_unordered();
    for ( i = point_list.size(); i--; )
      SubEntitySet::strip_attributes( point_list.get_and_step() );
    for ( i = surface_list.size(); i--; )
      SubEntitySet::strip_attributes( surface_list.get_and_step() );
    for ( i = lump_list.size(); i--; )
      SubEntitySet::strip_attributes( lump_list.get_and_step() );
  }
    
  return result;
}
/*
void PartitionEngine::destroy_facet( CubitFacetData* facet )
{
  int i;

  for( i = 0; i < 3; i++ )
  {
    CubitFacetEdgeData* edge = dynamic_cast<CubitFacetEdgeData*>(facet->edge(i));
    if( edge )
    {
      facet->edge(0, i);
      if( edge->number_tris() == 0 )
        delete edge;
    }
  }

  for( i = 0; i < 3; i++ )
  {
    CubitPointData* point = dynamic_cast<CubitPointData*>(facet->point(0));
    assert(point);
    point->remove_facet( facet );
    if( point->num_adj_facets() == 0 )
      delete point;
  }
  
  delete facet;
}
*/
//-------------------------------------------------------------------------
// Purpose       : Core functionality of make_body( X* ) methods.
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 02/23/03
//-------------------------------------------------------------------------
PartitionBody* PartitionEngine::make_body_internal( TopologyBridge* bridge )
{
    // query topology for BodySM
  DLIList<TopologyBridge*> list;
  BodySM* body;
  while( !(body = dynamic_cast<BodySM*>(bridge)) )
  {
    list.clean_out();
    bridge->get_parents_virt(list);
    if( !list.size() )
      return 0;
    bridge = list.get();
  }
  
  SubEntitySet* set = dynamic_cast<SubEntitySet*>(body->owner());
  if( set )
    return set->body();
  
  PartitionBody* result = new PartitionBody(body);
  if( body->owner() )
    body->owner()->swap_bridge( body, result, false );
  body->owner( &(result->sub_entity_set()) );
  return result;
}

//-------------------------------------------------------------------------
// Purpose       : Find or make an owning bodySM
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 02/23/03
//-------------------------------------------------------------------------
PartitionBody* PartitionEngine::make_body( PartitionLump* lump )
{
  if( lump->get_body() )
    return lump->get_body();
  
  return make_body_internal( lump );
}

//-------------------------------------------------------------------------
// Purpose       : Find or make an owning bodySM
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 02/23/03
//-------------------------------------------------------------------------
PartitionBody* PartitionEngine::make_body( PartitionSurface* surf )
{
  if( surf->sub_entity_set().body() )
    return surf->sub_entity_set().body();
    
  if( surf->next_co_surface(0) )
    return make_body(surf->next_co_surface(0)->get_shell()->get_lump());
  
  return make_body_internal(surf);
}

//-------------------------------------------------------------------------
// Purpose       : Find or make an owning bodySM
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 02/23/03
//-------------------------------------------------------------------------
PartitionBody* PartitionEngine::make_body( PartitionCurve* curve )
{
  if( curve->sub_entity_set().body() )
    return curve->sub_entity_set().body();
    
  PartitionCoEdge* coedge = 0;
  while( (coedge = curve->next_coedge(coedge)) )
    if( coedge->get_loop() )
      return make_body(coedge->get_loop()->get_surface());
  
  return make_body_internal(curve);
}
 
//-------------------------------------------------------------------------
// Purpose       : Find or make an owning bodySM
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 02/23/03
//-------------------------------------------------------------------------
PartitionBody* PartitionEngine::make_body( PartitionPoint* pt )
{
  if( pt->sub_entity_set().body() )
    return pt->sub_entity_set().body();
    
  if( pt->next_curve(0) )
    return make_body(pt->next_curve(0));
 
  return make_body_internal(pt);
}
static CubitStatus get_edge_replacements(
                         std::vector<CubitFacetData*> &facet_list,
                         std::vector<CubitFacetEdgeData*> replace_edge_lists[3])
{
  // first facet is the facet to replace
  // make sure more than one facet to replace the first
  assert(facet_list.size() > 2);

  // make sure the lists of replacement edges are empty
  assert(replace_edge_lists[0].size() == 0);
  assert(replace_edge_lists[1].size() == 0);
  assert(replace_edge_lists[2].size() == 0);


  CubitFacetData* p_dead = facet_list[0]; // facet to be replaced

  // look for edges on the facet that are owned by curves
  int i;
  //int replace_edge = 0;

  // get the edges
  CubitFacetEdgeData *edges[3];
  for (i=0; i<3; i++)
  {
    edges[i] = CAST_TO(p_dead->edge(i), CubitFacetEdgeData);
    assert(edges[i] != NULL);
  }

  PartitionEntity *p_owner;
  PartitionCurve *p_curve;

  // see if any are owned by a curve
  CubitBoolean b_boundary = CUBIT_FALSE;
  for (i=0; i<3 && !b_boundary; i++)
  {
    p_owner = TDVGFacetOwner::get(edges[i]);
    p_curve = CAST_TO(p_owner, PartitionCurve);
    if (p_curve)
      b_boundary = CUBIT_TRUE;
  }

  // if the facet has an edge on a curve, replace the edge with the edges of the
  // replacement facets and split the facet on the adjacent surface.
  if (b_boundary)
  {
    // make a list with just the new facets -- skip the first facet
    DLIList<CubitFacet*> new_facets(facet_list.size() - 1);
    for (unsigned u = 1; u < facet_list.size(); u++)
      new_facets.append(facet_list[u]);

    // get a point edge chain around the boundary of the replacement facets
    DLIList<FacetEntity*> point_edge_chain;
    FacetDataUtil::ordered_point_edge_bdry(new_facets, point_edge_chain);

    // for each edge on a curve get corresponding points and edges from the
    // replacement facets
    //std::vector<CubitFacetEdgeData*>::iterator eitor;
    for (i=0; i<3; i++)
    {
      p_owner = TDVGFacetOwner::get(edges[i]);
      p_curve = CAST_TO(p_owner, PartitionCurve);
      if (p_curve)
      {
        // fill the replacement list for this edge
        // add the edge to be replaced first
        //eitor = replace_edge_lists[i].begin();

        // get the points from the edge of the dead facet - in the correct order
        CubitPoint* pt1;
        CubitPoint* pt2;

        p_dead->get_edge_pts(i, pt1, pt2);

        // get the corresponding chain of points and edges from the replacement facets
        DLIList<FacetEntity*> replacement_chain;
        CubitStatus res = FacetDataUtil::partial_chain(point_edge_chain,
                                         pt1, pt2, replacement_chain);
        if (CUBIT_FAILURE == res)
          return CUBIT_FAILURE;

        // there should be an odd number of list entries, since it includes points
        // and edges, beginning and ending with a point
        assert( (replacement_chain.size() % 2) == 1 );

        // if only one replacement edge, it should be the same edge
        replacement_chain.reset();
        if (replacement_chain.size() == 3)
        {
          assert(replacement_chain.next(1) == edges[i]);
        }
        else
        {
          // add the replacement edges to the list - edge to replace first
          //*eitor++ = edges[i];
          replace_edge_lists[i].push_back(edges[i]);
          int j;
          CubitFacetEdgeData* p_edge;
          for (j=1; j<replacement_chain.size(); j+=2)
          {
            //replacement_chain.step();
            p_edge = CAST_TO(replacement_chain.next(j), CubitFacetEdgeData);
            assert(p_edge != NULL);
            //*eitor++ = p_edge;
            replace_edge_lists[i].push_back(p_edge);
          }
        }
      }
    }
  }

  return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Purpose       : clean up for deleted geometry
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 04/03/03
//-------------------------------------------------------------------------
CubitStatus PartitionEngine::delete_solid_model_entities( 
                                PartitionBody* body, BodySM*& real_body )
{
  real_body = body->real_body();
  body->destroy_all_children();
  return CUBIT_SUCCESS;
}

CubitStatus PartitionEngine::delete_solid_model_entities(
                                PartitionSurface* surf, 
                                Surface*& real_surf )
{
  DLIList<TopologyBridge*> parents(0);
  surf->get_parents_virt(parents);
  if ( parents.size() )
    return CUBIT_FAILURE;
  
  real_surf = 0;
  if ( !surf->sub_entity_set().has_multiple_sub_entities() )
    real_surf = dynamic_cast<Surface*>(surf->partitioned_entity());
  destroy_surface(surf);
  return CUBIT_SUCCESS;
}

CubitStatus PartitionEngine::delete_solid_model_entities(
                                PartitionCurve* curve,
                                Curve*& real_curve )
{
  if( curve->next_coedge(0) )
    return CUBIT_FAILURE;
    
  real_curve = 0;
  if ( !curve->sub_entity_set().has_multiple_sub_entities() )
    real_curve = dynamic_cast<Curve*>(curve->partitioned_entity());
  PartitionPoint* start = curve->start_point();
  PartitionPoint* end = curve->end_point();
  delete curve;
  if ( start->num_curves() == 0 )
    delete start;
  if ( end->num_curves() == 0 )
    delete end;
  return CUBIT_SUCCESS;
}

CubitStatus PartitionEngine::translate( PartitionEntity* ent, 
                                        const CubitVector& delta )
{
  GeometryEntity* geom = dynamic_cast<GeometryEntity*>(ent->partitioned_entity());
  if (!geom->get_geometry_query_engine()->translate( geom, delta ))
    return CUBIT_FAILURE;
  
  CubitTransformMatrix xform;
  xform.translate( delta );
  notify_transform_internal( geom, xform );
  return CUBIT_SUCCESS;
}
  
CubitStatus PartitionEngine::rotate( PartitionEntity* ent, 
                                     const CubitVector& axis, 
                                     double degrees )
{
  GeometryEntity* geom = dynamic_cast<GeometryEntity*>(ent->partitioned_entity());
  if (!geom->get_geometry_query_engine()->rotate( geom, axis, degrees ))
    return CUBIT_FAILURE;
  
  CubitTransformMatrix xform;
  xform.rotate( degrees, axis );
  notify_transform_internal( geom, xform );
  return CUBIT_SUCCESS;
}

CubitStatus PartitionEngine::scale( PartitionEntity* ent, 
                                    const CubitVector& factors )
{
  GeometryEntity* geom = dynamic_cast<GeometryEntity*>(ent->partitioned_entity());
  if (!geom->get_geometry_query_engine()->scale( geom, factors ))
    return CUBIT_FAILURE;
  
  CubitTransformMatrix xform;
  xform.scale_about_origin( factors );
  notify_transform_internal( geom, xform );
  return CUBIT_SUCCESS;
}

CubitStatus PartitionEngine::reflect( PartitionEntity* ent, 
                                      const CubitVector& axis )
{
  GeometryEntity* geom = dynamic_cast<GeometryEntity*>(ent->partitioned_entity());
  if (!geom->get_geometry_query_engine()->reflect( geom, axis ))
    return CUBIT_FAILURE;
  
  CubitTransformMatrix xform;
  xform.reflect( axis );
  notify_transform_internal( geom, xform );
  return CUBIT_SUCCESS;
}

CubitStatus PartitionEngine::translate( PartitionBody* ent, 
                                        const CubitVector& delta )
{
  BodySM* body = ent->real_body();
  if (!body->get_geometry_query_engine()->translate( body, delta ))
    return CUBIT_FAILURE;
  
  CubitTransformMatrix xform;
  xform.translate( delta );
  notify_transform_internal( ent, xform );
  return CUBIT_SUCCESS;
}

CubitStatus PartitionEngine::rotate( PartitionBody* ent, 
                                     const CubitVector& axis, 
                                     double degrees )
{
  BodySM* body = ent->real_body();
  if (!body->get_geometry_query_engine()->rotate( body, axis, degrees ))
    return CUBIT_FAILURE;
  
  CubitTransformMatrix xform;
  xform.rotate( degrees, axis );
  notify_transform_internal( ent, xform );
  return CUBIT_SUCCESS;
}

CubitStatus PartitionEngine::scale( PartitionBody* ent, 
                                    const CubitVector& factors )
{
  BodySM* body = ent->real_body();

  if( factors.x() != factors.y() ||
      factors.y() != factors.z() ||
      factors.x() != factors.z() )
  {
    GeometryModifyEngine *tmp_engine = GeometryModifyTool::instance()->get_engine( body );
    tmp_engine->scale( body, factors );  
  }
  else if(!body->get_geometry_query_engine()->scale( body, factors ))
    return CUBIT_FAILURE;
  
  CubitTransformMatrix xform;
  xform.scale_about_origin( factors );
  notify_transform_internal( ent, xform );
  return CUBIT_SUCCESS;
}

CubitStatus PartitionEngine::reflect( PartitionBody* ent, 
                                      const CubitVector& axis )
{
  BodySM* body = ent->real_body();
  if (!body->get_geometry_query_engine()->reflect( body, axis ))
    return CUBIT_FAILURE;
  
  CubitTransformMatrix xform;
  xform.reflect( axis );
  notify_transform_internal( ent, xform );
  return CUBIT_SUCCESS;
}

CubitStatus PartitionEngine::restore_transform( PartitionBody* ent )
{
  BodySM* body = ent->real_body();
  CubitTransformMatrix xform;
  body->get_transforms( xform );
  if (!body->get_geometry_query_engine()->restore_transform( body ))
    return CUBIT_FAILURE;
  
  notify_transform_internal( ent, xform );
  return CUBIT_SUCCESS;
}

CubitStatus PartitionEngine::notify_transform( TopologyBridge* ,
                                               const CubitTransformMatrix&  )
{
  return CUBIT_SUCCESS;
}

CubitStatus PartitionEngine::notify_transform_internal( TopologyBridge* ent,
                                               const CubitTransformMatrix& xform )
{
  int i;
  DLIList<Surface*> surfaces;
  DLIList<Curve*> curves;
  DLIList<TBPoint*> points;
  
  if (PartitionBody* body = dynamic_cast<PartitionBody*>(ent))
  {
    DLIList<PartitionEntity*> ent_list;
    body->get_all_children( ent_list );
    CAST_LIST( ent_list, surfaces, PartitionSurface );
    CAST_LIST( ent_list, curves, PartitionCurve );
    CAST_LIST( ent_list, points, PartitionPoint );
  }
  else
  {
    if(BodySM* bodysm = dynamic_cast<BodySM*>(ent))
    {
      DLIList<TopologyBridge*> lumps, shells, tmp_surfs;
      bodysm->get_children( lumps, true, layer() );
      while (lumps.size())
      {
        lumps.pop()->get_children( shells, true, layer() );
        while (shells.size())
        {
          tmp_surfs.clean_out();
          shells.pop()->get_children( tmp_surfs, true, layer() );
          while (tmp_surfs.size())
            surfaces.append( dynamic_cast<Surface*>(tmp_surfs.pop()) );
        }
      }
      surfaces.uniquify_unordered();
    }
    else if(Surface* surf = dynamic_cast<Surface*>(ent))
      surfaces.append( surf );
    else if(Curve* curv = dynamic_cast<Curve*>(ent))
      curves.append( curv );
    else if(TBPoint* point = dynamic_cast<TBPoint*>(ent))
      points.append( point );

    if (surfaces.size())
    {
      DLIList<TopologyBridge*> loops, coedges, tmp_curves;
      for (i = surfaces.size(); i--; )
      {
        surfaces.get_and_step()->get_children( loops, true, layer() );
        while( loops.size() )
        {
          loops.pop()->get_children( coedges, true, layer() );
          while (coedges.size())
          {
            coedges.pop()->get_children( tmp_curves, true, layer() );
            while (tmp_curves.size())
              curves.append( dynamic_cast<Curve*>(tmp_curves.pop()) );
          }
        }
      }
      curves.uniquify_unordered();
    }
    
    if (curves.size())
    {
      DLIList<TopologyBridge*> tmp_points;
      for (i = curves.size(); i--; )
      {
        curves.get_and_step()->get_children( tmp_points, true, layer() );
        while (tmp_points.size())
          points.append( dynamic_cast<TBPoint*>(tmp_points.pop()) );
      }
      points.uniquify_unordered();
    }
  }

  // see if the transformation is a reflection
  double det = xform.sub_matrix( 3, 3 ).determinant();
  bool reflection = det < 0.0;

  while (surfaces.size())
  {
    PartitionSurface* surf = dynamic_cast<PartitionSurface*>(surfaces.pop());
    if (surf)
    {
      surf->transform( xform );
    }
  }

  while (curves.size())
  {
    PartitionCurve* curv = dynamic_cast<PartitionCurve*>(curves.pop());
    if (curv)
    {
      curv->transform( xform );

      if (reflection)
      {
        // PartitionSurface transform above reverses loops and coedges on Partition
        // surfaces.  Need to reverse sense of PartitionCoEdges in loops on non
        // partition surfaces
        // reverse the coedges of this curve
        PartitionCoEdge* p_coedge = 0;
        while ((p_coedge = curv->next_coedge(p_coedge)) != NULL)
        {
          TopologyBridge* loop_bridge = p_coedge->find_parent_loop();
          if (loop_bridge)
          {
            if (0 == dynamic_cast<PartitionLoop*>(loop_bridge))
            {
              p_coedge->reverse_sense();
            }
          }
        }
      }
    }
  }
  
  while (points.size())
  {
    PartitionPoint* pnt = dynamic_cast<PartitionPoint*>(points.pop());
    if (pnt)
      pnt->transform( xform );
  }

  return CUBIT_SUCCESS;
}

void PartitionEngine::remove_attributes( DLIList<TopologyBridge*> &list )
{
  //for each SubEntitySet in map, change it's unique id to zero
  std::map<int,SubEntitySet*>::iterator itor = uniqueIdMap.begin();
  while( itor != uniqueIdMap.end() )
  {
    if(itor->second->get_entity())
    {
      SubEntitySet::strip_attributes( itor->second->get_entity() );
      itor->second->reset_unique_id();
    }
    itor++;
  }

  //clean out the map
  uniqueIdMap.clear();
}

void PartitionEngine::notify_deactivated (PartitionBody* body)
{
}
void PartitionEngine::notify_deactivated(PartitionLump* vol)
{
}
void PartitionEngine::add_to_deactivated_list(PartitionLump* vol)
{
}
void PartitionEngine::add_to_deactivated_list(PartitionBody* body)
{
}
void PartitionEngine::notify_deactivated (PartitionSurface* surface)
{
}
void PartitionEngine::add_to_deactivated_list (PartitionSurface* surface)
{
}
void PartitionEngine::notify_deactivated (PartitionCurve* curve)
{
}
void PartitionEngine::add_to_deactivated_list (PartitionCurve* curve)
{
}
void PartitionEngine::notify_deactivated (PartitionPoint* point)
{
}
void PartitionEngine::add_to_deactivated_list (PartitionPoint* point)
{
}
void PartitionEngine::clean_out_deactivated_geometry()
{
}
void PartitionEngine::remove_modified(DLIList<Surface*> &all_surfs,
    DLIList<Curve*> &all_curves, DLIList<TBPoint*> &all_pts)
{
}

void PartitionEngine::get_tbs_with_bridge_manager_as_owner( TopologyBridge *source_bridge, 
                                                            DLIList<TopologyBridge*> &tbs )
{
  SubEntitySet* set = dynamic_cast<SubEntitySet*>(source_bridge->owner());
  if( !set ) 
    return;

  DLIList<PartitionEntity*> entity_list;
  set->get_sub_entities( entity_list );
  DLIList<TopologyBridge*> tb_list;
  CAST_LIST( entity_list, tb_list, TopologyBridge );
  
  bool at_top = false;

  while( tb_list.size() )
  {
    TopologyBridge *tb = tb_list.pop();
    if( tb->bridge_manager() )    
      tbs.append( tb );          
    else
    {  
      TBOwner *owner = tb->owner();
      TopologyBridge *tmp_tb = CAST_TO( owner, TopologyBridge );
      if( tmp_tb )
        tb_list.append( tmp_tb );
      else
      {
        CompositePoint *comp_pt = CAST_TO( tb, CompositePoint );
        if( comp_pt )
          CompositeEngine::instance().get_tbs_with_bridge_manager_as_owner( tb, tb_list );
        else
        {
          CompositeCurve *comp_curve = CAST_TO( tb, CompositeCurve );
          if( comp_curve )
            CompositeEngine::instance().get_tbs_with_bridge_manager_as_owner( tb, tb_list );
        }
      }
    }    
  }

  tbs.uniquify_unordered();

  return;
}
  
