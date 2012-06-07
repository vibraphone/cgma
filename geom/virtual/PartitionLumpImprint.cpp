//-------------------------------------------------------------------------
// Filename      : PartitionLumpImprint.cpp
//
// Purpose       : Imprint a lump with a polyline-loop
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 02/05/03
//-------------------------------------------------------------------------

#include "PartitionLumpImprint.hpp"

#include "PartitionLump.hpp"
#include "PartitionSurface.hpp"
#include "PartitionLoop.hpp"
#include "PartitionCoEdge.hpp"
#include "PartitionCurve.hpp"
#include "PartitionPoint.hpp"

#include "CubitVector.hpp"
#include "PartitionEngine.hpp"
#include "SegmentedCurve.hpp"
#include "PartPTCurve.hpp"
#include "PartSurfFacetTool.hpp"

#include "CubitFacetData.hpp"
#include "CubitFacetEdgeData.hpp"
#include "CubitPointData.hpp"

const double DIST_TOL = 100.0 * GEOMETRY_RESABS;
const double DIST_TOL_SQR = DIST_TOL*DIST_TOL;

const int PART_LUMP_DEBUG = 88;
#define PART_LUMP_PRINT PRINT_DEBUG(PART_LUMP_DEBUG)
#include "GfxDebug.hpp"

//-------------------------------------------------------------------------
// Purpose       : Constructor - initialize some stuff
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 02/05/03
//-------------------------------------------------------------------------
PartitionLumpImprint::PartitionLumpImprint( PartitionLump* set_lump )
  : lump(set_lump)
{
  lump->get_all_children( entityList );
  entityList.reset();
  
  for( int i = 0; i < entityList.size(); i++ )
    rTree.add(entityList.get_and_step());
  entityList.append(lump);
  
  newSurface = 0;
}

//-------------------------------------------------------------------------
// Purpose       : Add an entity to internal structures
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 02/05/03
//-------------------------------------------------------------------------
bool PartitionLumpImprint::add( PartitionEntity* entity )
{
  if( !entityList.append_unique(entity) )
    return false;
    
  rTree.add(entity);
  return true;
}


//-------------------------------------------------------------------------
// Purpose       : Set up for imprint
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 02/05/03
//-------------------------------------------------------------------------
void PartitionLumpImprint::init( DLIList<CubitFacet*>* facets )
{
  int i, j, junk = 0;
  PART_LUMP_PRINT("Beginning imprint\n");
  
  boundaryEdges.clean_out();
  newEntities.clean_out();
  assert(newSurface == 0);

  if ( facets )
  {
      // clear all marks
    for ( i = facets->size(); i--; )
    {
      CubitFacet* facet = facets->get_and_step();
      facet->point(0)->marked(0);
      facet->point(1)->marked(0);
      facet->point(2)->marked(0);
    }

      // copy facets 
#ifdef BOYD17
    DLIList<CubitFacetData*> new_facets;
#endif
    DLIList<CubitPointData*> new_facet_pts;
    for ( i = facets->size(); i--; )
    {
      CubitFacet* facet = facets->get_and_step();
      CubitPointData *pts[3];
      for ( j = 0; j < 3; j++ )
      {
        CubitPoint* pt = facet->point(j);
        if ( !pt->marked() )
        {
          pts[j] = new CubitPointData(pt->coordinates());
          new_facet_pts.append(pts[j]);
          pt->marked(new_facet_pts.size());
        }
        else
        {
          new_facet_pts.reset();
          pts[j] = new_facet_pts.next(pt->marked()-1);
        }
      }
      facetList.append( new CubitFacetData(pts[0], pts[1], pts[2], &junk ) );
    }

      // clear marks on input facets
    for ( i = facets->size(); i--; )
    {
      CubitFacet* facet = facets->get_and_step();
      facet->point(0)->marked(0);
      facet->point(1)->marked(0);
      facet->point(2)->marked(0);
    }
  }

    // mark boundary edges with a 1, others with a 0
  for ( i = facetList.size(); i--; )
  {
    CubitFacetData* facet = facetList.get_and_step();
    for ( j = 0; j < 3; j++ )
      facet->edge(j)->marked( facet->edge(j)->marked() + 1 );
  }
  
  for ( i = facetList.size(); i--; )
  {
    CubitFacetData* facet = facetList.get_and_step();
    for ( j = 0; j < 3; j++ )
    {
      if ( facet->edge(j)->marked() == 1 )
        boundaryEdges.append( facet->edge(j) );
      facet->edge(j)->marked(0);
    }
  }
  
  for ( i = boundaryEdges.size(); i--; )
    boundaryEdges.get_and_step()->marked(1);
}

void PartitionLumpImprint::begin_loop( DLIList<CubitPoint*>& points )
{
  loopPoints = points;
  
    // associate each point with a geometric entity using rtree
  DLIList<PartitionEntity*> closest_list;
  CubitBox box;
  loopPoints.reset();
  for( int i = 0; i < loopPoints.size(); i++ )
  {
    CubitPoint* pt = loopPoints.next(i);
    closest_list.clean_out();
    CubitVector v = pt->coordinates();
    box.reset( CubitVector( v.x() - DIST_TOL,
                            v.y() - DIST_TOL,
                            v.z() - DIST_TOL ),
               CubitVector( v.x() + DIST_TOL,
                            v.y() + DIST_TOL,
                            v.z() + DIST_TOL ) );
                            
    rTree.find( box, closest_list );
    PartitionEntity* closest = find_closest( v, closest_list );
    if ( closest ) {
      set_point_owner(pt, closest);
    } else {
      set_point_owner(pt, lump);
    }
  }
}
    
//-------------------------------------------------------------------------
// Purpose       : Find closest entity in passed list to position
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 02/26/03
//-------------------------------------------------------------------------
PartitionEntity* PartitionLumpImprint::find_closest(
  const CubitVector& pt, 
  DLIList<PartitionEntity*>& closest_list,
  bool use_tolerance )
{ 
  PartitionPoint* ppoint;
  PartitionCurve* pcurve;
  PartitionSurface* psurf;
  closest_list.reset();
  double closest_dist_sqr = use_tolerance ? DIST_TOL_SQR : CUBIT_DBL_MAX;
  PartitionEntity* closest = 0;
  int closest_dim = 3;
  CubitVector c;

  for( int i = closest_list.size(); i--; )
  {
    PartitionEntity* entity = closest_list.get_and_step();
    int dim = 3;
    if( (ppoint = dynamic_cast<PartitionPoint*>(entity)) )
    {
      c = ppoint->coordinates();
      dim = 0;
    }
    else if( (pcurve = dynamic_cast<PartitionCurve*>(entity)) )
    {
      pcurve->closest_point_trimmed( pt, c );
      dim = 1;
    }
    else if( (psurf = dynamic_cast<PartitionSurface*>(entity)) )
    {
      psurf->closest_point_trimmed( pt, c );
      dim = 2;
    }
    else
      assert(0);


      // want the entity of smallest dimension within
      // tolerance of point.
    double dist_sqr = (pt - c).length_squared();
    if ( use_tolerance ) 
    {
      if( (dim == closest_dim && dist_sqr <  closest_dist_sqr) || 
          (dim <  closest_dim && dist_sqr <= DIST_TOL_SQR ) )
      {
        closest_dist_sqr = dist_sqr;
        closest = entity;
        closest_dim = dim;
      }
    } 
    else 
    {
      double diff = closest_dist_sqr - dist_sqr;
      if ( (dim  > closest_dim && diff >  DIST_TOL_SQR) ||
           (dim == closest_dim && diff >  0.0         ) ||
           (dim  < closest_dim && diff > -DIST_TOL_SQR) )
      {
        closest_dist_sqr = dist_sqr;
        closest = entity;
        closest_dim = dim;
      }
    }
  }

  return closest;
}
  
//-------------------------------------------------------------------------
// Purpose       : Clean up after imprint
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 02/05/03
//-------------------------------------------------------------------------
void PartitionLumpImprint::clean_up_loop()
{
  for( int i = entityList.size(); i--; )
    entityList.get_and_step()->mark = 0;
  entityList.reset();
  loopPoints.clean_out();
  pointAssoc.clear();
}

//-------------------------------------------------------------------------
// Purpose       : get the owner of a point in the loop of imprint points
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 02/05/03
//-------------------------------------------------------------------------
PartitionEntity* PartitionLumpImprint::point_owner( CubitPoint* pt ) 
{
  return pointAssoc[pt];
}

//-------------------------------------------------------------------------
// Purpose       : Set point owner
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 02/07/03
//-------------------------------------------------------------------------
void PartitionLumpImprint::set_point_owner( CubitPoint* pt, 
                                            PartitionEntity* owner )
{
  CubitVector closest;
  PartitionPoint* ppt = 0;
  PartitionCurve* pcv = 0;
  PartitionSurface* psf = 0;
  if( 0 != (ppt = dynamic_cast<PartitionPoint*>(owner)) ) {
    pt->set( ppt->coordinates() );
  } else if( 0 != (pcv = dynamic_cast<PartitionCurve*>(owner)) ) {
    pcv->closest_point( pt->coordinates(), closest );
    pt->set( closest );
  } else if( 0 != (psf = dynamic_cast<PartitionSurface*>(owner)) ) {
    psf->closest_point( pt->coordinates(), &closest );
    pt->set( closest );
  }

  if( DEBUG_FLAG(PART_LUMP_DEBUG) )
  {
    const char* type = 0;
    int color = 0;
    if( ppt ) {
      type = "Point";   color = CUBIT_BLUE;
    } else if( pcv ) {
      type = "Curve";   color = CUBIT_CYAN;
    } else if( psf ) {
      type = "Surface"; color = CUBIT_YELLOW;
    } else {
      type = "Lump";    color = CUBIT_RED;
    }
    
    int index = -1;
    if ( loopPoints.move_to(pt) )
      index = loopPoints.get_index();
    loopPoints.reset();

    TopologyBridge* tb = dynamic_cast<TopologyBridge*>(owner);
    RefEntity* re = dynamic_cast<RefEntity*>(tb->topology_entity());
    PART_LUMP_PRINT("%d. (%f,%f,%f) -> %s %p (%d)\n", index,
      pt->coordinates().x(), pt->coordinates().y(), pt->coordinates().z(), 
      type, static_cast<void*>(owner), re?re->id():0 );
    GfxDebug::draw_point( pt->coordinates(), color );
    GfxDebug::flush();
  }

  pointAssoc[pt] = owner;
}

//-------------------------------------------------------------------------
// Purpose       : Get points owned by passed entity
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 03/03/03
//-------------------------------------------------------------------------
void PartitionLumpImprint::get_owned_points( PartitionEntity* owner,
                                             DLIList<CubitPoint*>& results )
{
  std::map<CubitPoint*,PartitionEntity*>::iterator itor = pointAssoc.begin();
  for ( ; itor != pointAssoc.end(); ++itor )
    if( itor->second == owner )
      results.append(itor->first);
}


//-------------------------------------------------------------------------
// Purpose       : Entry point for imprinting
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 04/14/03
//-------------------------------------------------------------------------
PartitionSurface* PartitionLumpImprint::imprint( DLIList<CubitFacet*>& facets,
                                                 DLIList<PartitionEntity*>& new_list )
{
  init(&facets);
  DLIList<CubitVector*> empty_list;
  return imprint( empty_list, new_list );
}

PartitionSurface* PartitionLumpImprint::imprint( DLIList<CubitFacetData*>& facets,
                                                 DLIList<CubitVector*>& vertices,
                                                 DLIList<PartitionEntity*>& new_list )
{
  facetList = facets;
  init(0);
  return imprint( vertices, new_list );
}

PartitionSurface* PartitionLumpImprint::imprint( DLIList<CubitVector*>& vertices,
                                                 DLIList<PartitionEntity*>& new_list )
{
  int i, j;
  
  newSurface = new PartitionSurface( lump );
  
  
    // Make vertices where requested.
  DLIList<CubitPoint*> vtx_points, all_points;
  for ( j = 0; j < 3; j++ )
    for ( i = facetList.size(); i--; )
      facetList.get_and_step()->point(j)->marked(1);
  
  for ( j = 0; j < 3; j++ ) {
    for ( i = facetList.size(); i--; ) {
      if ( facetList.step_and_get()->point(j)->marked() ) {
        facetList.get()->point(j)->marked(0);
        all_points.append( facetList.get()->point(j) );
      }
    }
  }
  
  for ( i = vertices.size(); i--; )
  {
    double closest_sqr = CUBIT_DBL_MAX;
    CubitPoint* closest_pt = 0;
    CubitVector* vect = vertices.get_and_step();
    for ( j = all_points.size(); j--; )
    {
      CubitPoint* pt = all_points.step_and_get();
      double dist_sqr = (pt->coordinates() - *vect).length_squared();
      if ( dist_sqr < closest_sqr ) 
      {
        closest_sqr = dist_sqr;
        closest_pt = pt;
      }
    }
    vtx_points.append( closest_pt );
  }
  
  

    // For each closed chain of boundary edges
  DLIList<CubitFacetEdge*> pt_edges;
  DLIList<CubitPoint*> chain, chain_verts;
  int unused = boundaryEdges.size();
  while( unused )
  {
      // find first edge
    for ( i = boundaryEdges.size(); i--; )
      if ( boundaryEdges.step_and_get()->marked() )
        break;
    
    if ( !boundaryEdges.get()->marked() )
      break;
    
    chain.clean_out();
    chain_verts.clean_out();
    CubitFacetEdge* edge = boundaryEdges.get();
    edge->marked(0);
    unused--;
    CubitPoint* first_pt = edge->point(0);
    CubitPoint* point = edge->point(1);
    
    CubitFacet* facet = edge->adj_facet(0);
    int edge_index = facet->edge_index(edge);
    if ( facet->edge_use(edge_index) == -1 )
      std::swap(first_pt, point);
    
    while( point != first_pt )
    {
      chain.append( point );
      if ( vtx_points.is_in_list(point) )
        chain_verts.append(point);
        
      pt_edges.clean_out();
      point->edges(pt_edges);
      for ( j = pt_edges.size(); j--; )
        if ( pt_edges.step_and_get()->marked() )
          break;
      
      if ( !pt_edges.get()->marked() )
        break;
    
      edge = pt_edges.get();
      edge->marked(0);
      unused--;
      point = edge->other_point(point);
    }
    
    if ( point != first_pt )
      break;
    
    chain.append( first_pt );
    if ( vtx_points.is_in_list(first_pt) )
      chain_verts.append(first_pt);
      
      // remove edges from boundaryList that will be 
      // destroyed by merging with other facet edges.
    for ( i = boundaryEdges.size(); i--; )
      if( boundaryEdges.step_and_get()->marked() == 0 )
        boundaryEdges.change_to(0);
    boundaryEdges.remove_all_with_value(0);
    
    
    PartitionLoop* new_loop = imprint( chain, chain_verts );
    if (!new_loop)
    {
      unused = 1;
      break;
    }
    
    newSurface->add(new_loop);
  }
  
  if ( !unused )
  {
    new_list = newEntities;
    PartitionSurface* result = newSurface;
    result->set_facet_data( facetList );
    newSurface = 0;
    return result;
  }
  
  abort_imprint();
  return 0;
}
  

//-------------------------------------------------------------------------
// Purpose       : Imprint one loop onto the volume
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 02/05/03
//-------------------------------------------------------------------------
PartitionLoop* PartitionLumpImprint::imprint( DLIList<CubitPoint*>& loop ,
                                              DLIList<CubitPoint*>& vtx_points)
{
  begin_loop( loop );
  DLIList<PartitionCoEdge*> results;
  
  if( ! make_vertices( vtx_points ) )
  {
    PRINT_ERROR("PartitionLumpImprint::make_vertices() failed.\n");
    return 0;
  }
  
  if( ! do_imprint() )
  {
    PRINT_ERROR("PartitionLumpImprint::do_imprint() failed.\n");
    return 0;
  }
  
  if( ! make_volume_curves() )
  {
    PRINT_ERROR("PartitionLumpImprint::make_volume_curves() failed.\n");
    return 0;
   }
 
  if( ! get_curves(results) )
  {
    PRINT_ERROR("PartitionLumpImprint::get_curves() failed.\n");
    return 0;
  }

//for ( int i = 0; i < results.size(); i++ )
//  results.get_and_step()->draw_facets(CUBIT_RED);

  
  PartitionLoop* new_loop = new PartitionLoop;
  results.reverse();
  PartitionCoEdge* prev = 0;
  while ( results.size() )
  {
    PartitionCoEdge* coe = results.pop();
    new_loop->insert_after( coe, prev );
    prev = coe;
  }
  
  clean_up_loop();
  
  return new_loop;
}

//-------------------------------------------------------------------------
// Purpose       : Create vertices where asked to do so
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 02/05/03
//-------------------------------------------------------------------------
CubitStatus PartitionLumpImprint::make_vertices( DLIList<CubitPoint*>& vtx_points )
{
  for( int i = 0; i < vtx_points.size(); i++ )
  {
    CubitPoint* pt = vtx_points.get_and_step();
    PartitionEntity* entity = point_owner(pt);
    if( !entity )
      continue;
    
    
    if( dynamic_cast<Point*>(entity) ) {
     ; // already a point
    }
    else if( dynamic_cast<Curve*>(entity) ) {
      if( !partitionCurve( pt ) )
        return CUBIT_FAILURE;
    }
    else if( dynamic_cast<Surface*>(entity) ) {
      if( !makePointCurve( pt ) )
        return CUBIT_FAILURE;
    }
    else if( !makeFreePoint( pt ) ) {
      return CUBIT_FAILURE;
    }
  }
  return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Purpose       : Try to clean up failed imprint
//
// Special Notes : Always returns CUBIT_FAILURE
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 02/05/03
//-------------------------------------------------------------------------
CubitStatus PartitionLumpImprint::abort_imprint()
{
  PRINT_ERROR("Loop-Volume imprint failed in PartitionLumpImprint\n"
              "\tAttempting to restore state.\n");
  
  clean_up_loop();

  DLIList<PartitionCurve*> curves;
  DLIList<PartitionPoint*> points;
  CAST_LIST( newEntities, curves, PartitionCurve );
  CAST_LIST( newEntities, points, PartitionPoint );  
  
  while( curves.size() )
  {
    PartitionCurve* curve = curves.pop();

    // if num curves for the startpoint/endpoint equals one, this point will
    // get cleaned up when the curve is deleted or removed
    // so remove the point from the list of points to clean up
    if (1 == curve->start_point()->num_curves())
      points.remove(curve->start_point());
    if (1 == curve->end_point()->num_curves())
      points.remove(curve->end_point());

    if( &(curve->sub_entity_set()) == &(lump->sub_entity_set()) )
      delete curve;
    else
      PartitionEngine::instance().remove_curve(curve);
  }


  while( points.size() )
  {
    PartitionPoint* point = points.pop();
    if( point->next_curve() == 0 )
      delete point;
    else if( dynamic_cast<PartPTCurve*>(point->next_curve()) )
      PartitionEngine::instance().remove_point_curve(point);
    else if( dynamic_cast<Curve*>(point->partitioned_entity()) )
      PartitionEngine::instance().remove_point(point);
  }
  
  while( facetList.size() )
    PartitionEngine::delete_facet( facetList.pop() );
  
  return CUBIT_FAILURE;
}

//-------------------------------------------------------------------------
// Purpose       : Partition a curve and update local data
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 02/05/03
//-------------------------------------------------------------------------
CubitStatus PartitionLumpImprint::partitionCurve( CubitPoint* pt )
{
    // get entities
  PartitionEntity* ent = point_owner(pt);
  PartitionCurve* curve = dynamic_cast<PartitionCurve*>(ent);
  assert(!!curve);
  
  PART_LUMP_PRINT("Splitting curve %p (%d) at (%f,%f,%f)\n",
    static_cast<void*>(curve),
    dynamic_cast<RefEntity*>(curve->topology_entity()) ?
      dynamic_cast<RefEntity*>(curve->topology_entity())->id() : 0,
    pt->coordinates().x(), pt->coordinates().y(), pt->coordinates().z() );
  
    // get list of associated points that will need to be moved
    // to the new partition of the curve after the split.
  double u_end = curve->u_from_position( pt->coordinates() );
  DLIList<CubitPoint*> points_to_move;
  get_owned_points( curve, points_to_move );
  
  CubitPoint* curve_pt;
  for( int i = points_to_move.size(); i--; )
  {
    curve_pt = points_to_move.step_and_get();
    double u = curve->u_from_position( curve_pt->coordinates() );
    if( u < u_end || curve_pt == pt )
      points_to_move.change_to(0);
  }
  points_to_move.remove_all_with_value(0);
  
    // partition curve
  PartitionPoint* new_pt = new PartitionPoint( pt->coordinates(), curve );
  PartitionCurve* new_curve = 
    PartitionEngine::instance().insert_point( curve, new_pt );
  
    // partition failed
  if( !new_curve )
  {
    delete new_pt;
    return CUBIT_FAILURE;
  }
  
    // update internal lists
  add(new_pt);
  add(new_curve);
  newEntities.append(new_pt);
  pt->set( new_pt->coordinates() );
  set_point_owner( pt, new_pt );
  
    // update point associativity
  while( points_to_move.size() )
    set_point_owner( points_to_move.pop(), new_curve );
  
  return CUBIT_SUCCESS;
}


//-------------------------------------------------------------------------
// Purpose       : Partition a surface and update internal data structures
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 02/05/03
//-------------------------------------------------------------------------
CubitStatus PartitionLumpImprint::partitionSurface( int first_point_id,
                                                    int last_point_id,
                                                    PartitionSurface* surface )
{
  PartitionEntity* start_owner = point_owner(point(first_point_id));
  PartitionEntity* end_owner = point_owner(point(last_point_id));
  PartitionPoint* start_point = dynamic_cast<PartitionPoint*>(start_owner);
  PartitionPoint* end_point = dynamic_cast<PartitionPoint*>(end_owner);
  if ( (start_owner != lump && !start_point) ||
       (end_owner != lump && !end_point) ) {
    PRINT_ERROR("Facet boundary too coarse to resolve imprint.\n");
    return CUBIT_FAILURE;
  } 

  int i, next_point_id = (first_point_id + 1) % num_points();
  if( !surface )
  {
    assert(next_point_id != last_point_id);
    surface = dynamic_cast<PartitionSurface*>(point_owner(point(next_point_id)));
    assert(!!surface);
  }
  
  if ( PART_LUMP_DEBUG ) {
    RefEntity* re = dynamic_cast<RefEntity*>(surface->topology_entity());
    PART_LUMP_PRINT("Inserting curve into surface %p (%d)\n",
      static_cast<void*>(surface), re ? re->id() : 0 );
  }
  
  DLIList<CubitPoint*> segments;
  first_point_id %= num_points();
  int last_id = last_point_id % num_points();
  
  int start, stop;
  if ( start_point )
    start = first_point_id;
  else
    start = next_point_id;
  if ( end_point )
    stop = (last_id + 1) % num_points();
  else
    stop = last_id;
  
  if( start == last_id )
  {
    i = start;
    do {
      segments.append(point(i));
      i = (i+1) % num_points();
    } while( i != start );
    segments.append(point(start));
  }
  else
  {
    for ( i = start; i != stop; i = (i + 1) % num_points() )
      segments.append(point(i));
  }
  
  assert(surface&&start_owner&&end_owner && start_owner!=surface && end_owner!=surface );

  DLIList<CubitVector*> segment_points;
  segments.reset();
  for ( i = segments.size(); i--; )
    segment_points.append( new CubitVector( segments.get_and_step()->coordinates() ) );
  
  DLIList<PartitionSurface*> input_list(1), new_surfs;
  DLIList<PartitionCurve*> new_curves;
#ifdef BOYD17
  DLIList<PartitionPoint*> points;
#endif
  input_list.append(surface);
  CubitStatus result = PartitionEngine::instance().
    insert_curve( input_list, segment_points, new_surfs, new_curves );
  while( segment_points.size() )
    delete segment_points.pop();
    
  if ( !result )
  {
    PRINT_ERROR("Surface partitioning failed in PartitionLumpImprint\n");
    return CUBIT_FAILURE;
  }
  
    // Add all new geometry to 'entities'.  Append everything except
    // the surfaces to curves_and_points.
  DLIList<PartitionEntity*> entities;
  DLIList<PartitionPoint*> new_points;
    
    // Do Surfaces
  if( new_surfs.move_to( surface ) )
    new_surfs.extract();
  new_surfs.reset();
  for ( i = new_surfs.size(); i--; )
    entities.append(new_surfs.get_and_step());
    
    // Get new points
  for ( i = new_curves.size(); i--; )
  {
    new_curves.get()->start_point()->mark = 1;
    new_curves.get()->end_point()->mark = 1;
    new_curves.step();
  }
  if ( start_point )
    start_point->mark = 0;
  if ( end_point )
    end_point->mark = 0;
  for ( i = new_curves.size(); i--; )
  {
    if ( new_curves.get()->start_point()->mark )
    {
      PartitionPoint* pt = new_curves.get()->start_point();
      pt->mark = 0;
      new_points.append( pt );
      entities.append( pt );
    }
    if ( new_curves.get()->end_point()->mark )
    {
      PartitionPoint* pt = new_curves.get()->end_point();
      pt->mark = 0;
      new_points.append( pt );
      entities.append( pt );
    }
    entities.append( new_curves.get() );
    new_curves.step();
  }

    // If first and/or last segment points were already
    // associated with a vertex, don't search for a new
    // entity to associate them with.
  if( start_point ) {
    segments.reset();
    segments.remove();
  }
  if( end_point ) {
    segments.pop();
  }

  // Associate split segments with points
  for ( i = new_points.size(); i--; )
  {
    if (segments.size() == 0)
	break;

    PartitionPoint* vert = new_points.get_and_step();
    CubitPoint* pt = 0;
    segments.reset();
    int closest_index = -1;
    double closest_sqr = CUBIT_DBL_MAX;
    for ( int j = 0; j < segments.size(); j++ )
    {
      pt = segments.get_and_step();
      double dist_sqr = (pt->coordinates() - vert->coordinates()).length_squared();
      if ( dist_sqr < closest_sqr ) {
        closest_sqr = dist_sqr;
        closest_index = j;
      }
    }
    assert(closest_index >= 0);
    segments.reset();
    segments.step(closest_index);
    pt = segments.extract();
    assert(point_owner(pt) == surface);
    pt->set( vert->coordinates() );
    set_point_owner( pt, vert );
  }

  // Associate remaining segment points with curves
  PartitionEntity* entity;
  while( segments.size() )
  {
    CubitPoint* pt = segments.pop();
    PartitionCurve* closest = 0;
    double closest_sqr = CUBIT_DBL_MAX;
    CubitVector pos;
    for ( i = new_curves.size(); i--; )
    {
      PartitionCurve* curve = new_curves.get_and_step();
      curve->closest_point_trimmed( pt->coordinates(), pos );
      double dist_sqr = (pos - pt->coordinates()).length_squared();
      if ( dist_sqr < closest_sqr )
      {
        closest_sqr = dist_sqr;
        closest = curve;
      }
    }
    assert(closest != NULL);
    closest->closest_point( pt->coordinates(), pos );
    pt->set(pos);
    set_point_owner(pt, closest);
  }

  
    // append everything to newEntites and call add for each
  newEntities += entities;
  entities.reset();
  for ( i = entities.size(); i--; )
    add( entities.get_and_step() );
    
    // now add the orignal surface to entities  -- entities
    // becomes the list of entities that other polyline points
    // on the surface must lie on after the split
  entities.append(surface);
    
    // Associate any other points on the surface that 
    // haven't been used yet with the appropriate 
    // sub-region of the surface partitions.
  DLIList<CubitPoint*> surf_points;
  get_owned_points( surface, surf_points );
  for ( i = surf_points.size(); i--; )
  {
    CubitPoint* pt = surf_points.get_and_step();
    entity = find_closest( pt->coordinates(), entities, false );
    assert(!!entity);
    set_point_owner(pt, entity);
  }
  
  return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Purpose       : Create a point-curve on a surface
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 02/05/03
//-------------------------------------------------------------------------
CubitStatus PartitionLumpImprint::makePointCurve( CubitPoint* pt )
{
  PartitionEntity* ent = point_owner(pt);
  PartitionSurface* surf = dynamic_cast<PartitionSurface*>(ent);
  CubitPointData* cpd = 
    PartitionEngine::instance().project_to_surface(surf, pt->coordinates());
  
  PART_LUMP_PRINT("Creating vertex on surface %p (%d) at (%f,%f,%f)\n",
    static_cast<void*>(surf),
    dynamic_cast<RefEntity*>(surf->topology_entity()) ?
      dynamic_cast<RefEntity*>(surf->topology_entity())->id() : 0,
    pt->coordinates().x(), pt->coordinates().y(), pt->coordinates().z() );
    
  PartitionPoint* new_point = 0;
  if( cpd )
    new_point = PartitionEngine::instance().insert_point_curve( surf, cpd );
  
  if ( !new_point )
  {
   PRINT_ERROR("PartitionLumpImprint: point-curve creation failed.\n");
   return CUBIT_FAILURE;
  }
  
  add(new_point);
  newEntities.append(new_point);
  pt->set(new_point->coordinates());
  set_point_owner( pt, new_point );
  return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Purpose       : Create a free point in the volume
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 02/05/03
//-------------------------------------------------------------------------
CubitStatus PartitionLumpImprint::makeFreePoint( CubitPoint* pt )
{
  
  PART_LUMP_PRINT("Creating free vertex at (%f,%f,%f)\n",
    pt->coordinates().x(), pt->coordinates().y(), pt->coordinates().z() );

  assert( point_owner(pt) == lump );
  PartitionPoint* new_pt = new PartitionPoint( pt->coordinates(), lump );
  add(new_pt);
  newEntities.append(new_pt);
  set_point_owner( pt, new_pt );
  return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Purpose       : Create a free curve in the volume
//
// Special Notes : End points must already exist
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 02/05/03
//-------------------------------------------------------------------------
CubitStatus PartitionLumpImprint::makeFreeCurve( int start_id, int end_id )
{
  PartitionPoint* start_pt 
    = dynamic_cast<PartitionPoint*>(point_owner(point(start_id)));
  PartitionPoint* end_pt 
    = dynamic_cast<PartitionPoint*>(point_owner(point(end_id)));
  if( !start_pt || !end_pt )
    return CUBIT_FAILURE;
  
  PART_LUMP_PRINT("Creating free curve from point %p (%d) (%f,%f,%f)"
                  "  to point %p (%d) (%f,%f,%f)\n",
    static_cast<void*>(start_pt),
    dynamic_cast<RefEntity*>(start_pt->topology_entity()) ?
      dynamic_cast<RefEntity*>(start_pt->topology_entity())->id() : 0,
    start_pt->coordinates().x(), start_pt->coordinates().y(), start_pt->coordinates().z(),
    static_cast<void*>(end_pt),
    dynamic_cast<RefEntity*>(end_pt->topology_entity()) ?
      dynamic_cast<RefEntity*>(end_pt->topology_entity())->id() : 0,
    end_pt->coordinates().x(), end_pt->coordinates().y(), end_pt->coordinates().z());
    
  bool all_pts = false;
  if( start_id == end_id )
  {
    end_id = (end_id+1)%num_points();
    all_pts = true;
  }
  
  int i;
  DLIList<CubitPoint*> points;
  points.append( point(start_id) );
  end_id %= num_points();
  for( i = (start_id + 1) % num_points(); i != end_id; i = (i+1)%num_points() )
  {
    if( point_owner(point(i)) != lump )
      return CUBIT_FAILURE;
    points.append( point(i) );
  }
  points.append( point(all_pts ? start_id : end_id) );
  
  DLIList<CubitVector*> curve_segments;
  points.reset();
  for( i = points.size(); i--; )
    curve_segments.append( new CubitVector(points.get_and_step()->coordinates()));
  SegmentedCurve* curve = new SegmentedCurve( lump, curve_segments );
  while( curve_segments.size() )
    delete curve_segments.pop();
  curve->start_point(start_pt);
  curve->end_point(end_pt);
  add(curve);
  newEntities.append(curve);
  for( i = (start_id + 1) % num_points(); i != end_id; i = (i+1)%num_points() )
    set_point_owner( point(i), curve );

  return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Purpose       : Partition curves and surfaces on imprint loop
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 02/06/03
//-------------------------------------------------------------------------
CubitStatus PartitionLumpImprint::do_imprint() 
{
    // partition curves
  int i;
  for( i = 0; i <= num_points(); i++ )
  {
      // If the point owner isn't a curve, skip it.
    PartitionEntity* curr = point_owner(point(i));
    PartitionCurve* curve = dynamic_cast<PartitionCurve*>(curr);
    if( !curve )
      continue;
      
      // Split the curve at any point at which the polyline 
      // leaves the curve (except the end points, of course)
    PartitionEntity* prev = point_owner(point(i-1));
    PartitionEntity* next = point_owner(point(i+1));
    if( (prev != curr && 
         prev != curve->start_point() && 
         prev != curve->end_point()) ||
        (next != curr &&
         next != curve->start_point() &&
         next != curve->end_point()) )
    {
      if( !partitionCurve( point(i) ) )
        return CUBIT_FAILURE;
    }
  }            
      
    // partition surfaces
  
    // find position to start at
  for( i = 0; i < num_points(); i++ )
    if( point_owner(point(i)) != point_owner(point(i+1)) ||
        !dynamic_cast<PartitionSurface*>(point_owner(point(i))) )
      break;
  
    // loop through all points
  int index = i + 1;
  for( i = 0; i <= num_points(); )
  {
      // loop until we are at a point on a surface
    PartitionSurface* surf = 0;
    while( i <= num_points() &&
         !(surf = dynamic_cast<PartitionSurface*>(point_owner(point(index)))) )
    {
      i++;
      index++;
    }
    
    if( i > num_points() )
      break;
    
      // need the previous point too
    int prev_index = index - 1;
    
      // find remaining points on surface
    index++;
    i++;
    while( i <= num_points() && (point_owner(point(index)) == surf) )
    {
      index++;
      i++;
    }
    
      // loop stoped past last point on surface
      // store it and back up one.
    int next_index = index;
    index--;
    i--;
    
      // partition the surface
    if( ! partitionSurface( prev_index, next_index ) )
      return CUBIT_FAILURE;
  }
  
    // check for additional surface partitions where
    // the polyline has points on the boundary but
    // not the interior of the surface
  for( i = 0; i <= num_points(); i++ )
  {
    PartitionPoint* pt1 = dynamic_cast<PartitionPoint*>(point_owner(point(i)));
    PartitionPoint* pt2 = dynamic_cast<PartitionPoint*>(point_owner(point(i+1)));
    if( !pt1 || !pt2 || pt1->common_curve(pt2) )
      continue;

    CubitVector close, midpt = (pt1->coordinates() + pt2->coordinates()) / 2;
    PartitionSurface* surf = 0;
    PartitionCurve* curve = 0;
    double closest_dist = DIST_TOL_SQR;
    while( (curve = pt1->next_curve(curve)) )
    {
      PartitionCoEdge* coedge = 0;
      while( (coedge = curve->next_coedge(coedge)) )
      {
        PartitionSurface* surface = coedge->get_loop()->get_surface();
        if( surface == surf )
          continue;
          
        bool surf_contains_pt2 = false;
        PartitionCurve* curve2 = 0;
        while( (curve2 = pt2->next_curve(curve2)) && !surf_contains_pt2 )
        {
          PartitionCoEdge* coedge2 = 0;
          while( (coedge2 = curve2->next_coedge(coedge2)) )
          {
            if( coedge2->get_loop()->get_surface() == surface )
            {
              surf_contains_pt2 = true;
              break;
            }
          }
        }
        if( ! surf_contains_pt2 )
          continue;
          
        surface->closest_point( midpt, &close );
        double dist_sqr = (close - midpt).length_squared();
        if( dist_sqr <= DIST_TOL_SQR &&
            (!surf || dist_sqr < closest_dist) )
          surf = surface;
          closest_dist = dist_sqr;
      }
    }
      
    int end = i + 1;
    if( surf && !partitionSurface(i, end, surf) )
      return CUBIT_FAILURE;
    i = end - 1;
  }

  return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Purpose       : Construct free curves in lump interior to complete
//                 the imprint loop
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 02/06/03
//-------------------------------------------------------------------------
CubitStatus PartitionLumpImprint::make_volume_curves()
{
    // start at any point not owned by the volume
  int i;
  for( i = 0; i <= num_points() && point_owner(point(i)) == lump; i++ );
   
    // if everything is owned by the volume, make one vertex
  if( point_owner(point(i)) == lump && !makeFreePoint(point(i)) )
    return CUBIT_FAILURE;
  
    // create free curves from points in volume interior
  int index = i % num_points();
  for( i = 0; i <= num_points(); )
  {
      // advance until a point in the volume interior
    while( i <= num_points() && point_owner(point(index)) != lump )
    {
      i++;
      index = (index + 1) % num_points();
    }
    
    if( i > num_points() )
      break;
    
      // vertex to begin curve is at previous index
    int prev_index = (index + num_points() - 1) % num_points();
    
      // advance util we find a point not in the volume interior
      // (actually until a vertex possibly in the volume interior)
    while( i <= num_points() && point_owner(point(index)) == lump )
    {
      i++;
      index = (index + 1) % num_points();
    }
    
      // construct the curve
    if( !makeFreeCurve( prev_index, index ) )
      return CUBIT_FAILURE;
  }
  
    // create free curves where two points in are on vertices
    // but there does not already exist a curve connecting the
    // points
  for( i = 0; i <= num_points(); i++ )
  {
    PartitionPoint* pt1 = dynamic_cast<PartitionPoint*>(point_owner(point(i)));
    PartitionPoint* pt2 = dynamic_cast<PartitionPoint*>(point_owner(point(i+1)));
    if( pt1 && pt2 && !pt1->common_curve(pt2) && !makeFreeCurve(i,i+1) )
      return CUBIT_FAILURE;
  }
  
  return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Purpose       : Get the curves in the imprint loop
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 02/06/03
//-------------------------------------------------------------------------
CubitStatus PartitionLumpImprint::get_curves( DLIList<PartitionCoEdge*>& list )
{
  list.clean_out();
  DLIList<CubitFacetEdgeData*> facet_edges, curve_edges;
  CubitStatus result = CUBIT_SUCCESS;
  
    // search for a vertex to begin at
  int i;
  PartitionPoint* pt = 0;
  for( i = 0; i < num_points(); i++ )
    if( (pt = dynamic_cast<PartitionPoint*>(point_owner(point(i)))) )
      break;
  if( !pt )
    return CUBIT_FAILURE;
  CubitPoint* prev_facet_pt = point(i);
  
  int index = i + 1;
  PartitionCurve* curve = 0;
  for( i = 0; result && i < num_points(); i++, index++ )
  {
    CubitPoint* facet_pt = point(index);
    CubitFacetEdgeData* edge = 
      dynamic_cast<CubitFacetEdgeData*>(prev_facet_pt->shared_edge( facet_pt ));
    if ( !edge ) {
      result = CUBIT_FAILURE;
      break;
    }
    facet_edges.append(edge);
    prev_facet_pt = facet_pt;
    
    PartitionPoint* pt2;
    PartitionCurve* curve2;
    if( (pt2 = dynamic_cast<PartitionPoint*>(point_owner(point(index)))) )
    {
      if( curve )
      {
        if( curve->other_point(pt) != pt2 )
          result = CUBIT_FAILURE;
      }
      else
      {
          //find common curves
        curve = 0;
        PartitionCurve* tmp_curve = 0;
        while( (tmp_curve = pt->next_curve(tmp_curve) ) != NULL )
        {
          if ( tmp_curve->other_point(pt) != pt2 )
            continue;
          
          if( !curve )
          {
            curve = tmp_curve;
            continue;
          }
          
          CubitVector mid = 0.5 * (pt->coordinates() + pt2->coordinates());
          CubitVector cpt, tpt;
          curve->closest_point( mid, cpt );
          tmp_curve->closest_point( mid, tpt );
          cpt -= mid;
          tpt -= mid;
          if ( cpt.length_squared() > tpt.length_squared() )
            curve = tmp_curve;
        }

        if( !curve )
          result = CUBIT_FAILURE;
      }
      
      if ( !result )
        break;

      CubitSense sense;
      if ( curve->start_point() != curve->end_point() ) {
        sense = curve->start_point() == pt ? CUBIT_FORWARD : CUBIT_REVERSED;
      } else {
        CubitVector junk, tangent;
        facet_edges.last();
        CubitVector es = facet_edges.get()->point(0)->coordinates();
        CubitVector ee = facet_edges.get()->point(1)->coordinates();
        curve->closest_point( 0.5*(es+ee), junk, &tangent );
        sense = (ee-es)%tangent < 0.0 ? CUBIT_REVERSED : CUBIT_FORWARD;
      }

      PartitionCoEdge* new_coedge = new PartitionCoEdge( newSurface, sense );
      curve->add(new_coedge);
      list.append(new_coedge);
      curve = 0;
      pt = pt2;
      facet_edges.append(0);
    }
    else 
    {
      curve2 = dynamic_cast<PartitionCurve*>(point_owner(point(index)));
      
      if( !curve2 || (curve && curve != curve2) )
        result = CUBIT_FAILURE;
      else
        curve = curve2;
    }
  }
  
  if ( result )
  {
      // seam facet edges
    facet_edges.reset();
    list.reset();
    int used = facet_edges.size();
    for ( i = list.size(); i--; )
    {
      PartitionCurve* curve = list.get_and_step()->get_curve();
      curve_edges.clean_out();
      CubitFacetEdgeData* edge;
      while ( used-- && (edge = facet_edges.get_and_step()) )
        curve_edges.append(edge);
  
      assert(curve_edges.size());
      if ( !seam_curve( curve_edges, curve, facetList ) )
      {
        result = CUBIT_FAILURE;
        break;
      }
    }
  }
  
  if (result)
    return CUBIT_SUCCESS;
  
  while ( list.size() )
  {
    PartitionCoEdge* coedge = list.pop();
    coedge->get_curve()->remove(coedge);
    delete coedge;
  }
  return CUBIT_FAILURE;
}

CubitStatus PartitionLumpImprint::seam_curve ( DLIList<CubitFacetEdgeData*>& edges,
                                               PartitionCurve* curve,
                                               DLIList<CubitFacetData*>& facets )
{
  // Get start and end facet points for chain of edges
  CubitPoint *start_point_t, *end_point_t;
  if (edges.size() == 1)
  {
    start_point_t = edges.get()->point(0);
    end_point_t = edges.get()->point(1);
  }
  else
  {
    edges.last();
    end_point_t = edges.get()->shared_point( edges.prev() );
    end_point_t = edges.get()->other_point( end_point_t );
    edges.reset();
    start_point_t = edges.get()->shared_point( edges.next() );
    start_point_t = edges.get()->other_point( start_point_t );
  }
  
  CubitPointData* start_point = dynamic_cast<CubitPointData*>(start_point_t);
  CubitPointData*   end_point = dynamic_cast<CubitPointData*>(  end_point_t);
  assert(start_point && end_point);
  
  // If list if edges is backwards on curve, reverse it
  if (curve->start_point() == curve->end_point())
  {
    assert( start_point == end_point );
    edges.reset();
    CubitVector s = start_point->coordinates();
    CubitVector e = edges.get()->other_point(start_point)->coordinates();
    CubitVector j, t;
    curve->closest_point ( 0.5 * (s + e), j, &t );
    if ( (e - s) % t < 0.0 )
      edges.reverse();
  }
  else
  {
    if (curve->start_point()->facet_point() == end_point ||
        curve->end_point()->facet_point() == start_point)
    {
      std::swap(start_point, end_point);
      edges.reverse();
    }
    else if (curve->start_point()->facet_point() != start_point &&
             curve->end_point()->facet_point() != end_point)
    {
      double fwd = (curve->start_point()->coordinates() - start_point->coordinates()).length()
                 + (curve->end_point()->coordinates() - end_point->coordinates()).length();
      double rev = (curve->end_point()->coordinates() - start_point->coordinates()).length()
                 + (curve->start_point()->coordinates() - end_point->coordinates()).length();
      if (fwd > rev)
      {
        std::swap(start_point, end_point);
        edges.reverse();
      }
    }
  }
  
  // Merge facet points at start and end vertices
  if (curve->start_point()->facet_point() != start_point)
  {
    if (curve->start_point()->facet_point())
      curve->start_point()->facet_point()->merge_points(start_point);
    else
      curve->start_point()->facet_point(start_point);
  }
  if (curve->start_point() != curve->end_point() &&
      curve->end_point()->facet_point() != end_point)
  {
    if (curve->end_point()->facet_point())
      curve->end_point()->facet_point()->merge_points(end_point);
    else
      curve->end_point()->facet_point(end_point);
  }
  
  // Now seam interior points/edges
  return PartSurfFacetTool::seam_curve( edges, curve, facets );
}
