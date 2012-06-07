//-------------------------------------------------------------------------
// Filename      : PartitionCurve-new.cpp
//
// Purpose       : 
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 04/10/02
//-------------------------------------------------------------------------

#include "PartitionCurve.hpp"
#include "PartitionPoint.hpp"
#include "PartitionCoEdge.hpp"
#include "PartitionSurface.hpp"
#include "VirtualQueryEngine.hpp"
#include "PartSurfFacetTool.hpp"

#include "TDVGFacetOwner.hpp"
#include "CubitFacetEdgeData.hpp"
#include "CubitPointData.hpp"
#include "CubitFacetData.hpp"
#include "CubitTransformMatrix.hpp"
#include "Surface.hpp"

#include "GfxDebug.hpp"
#include "GMem.hpp"

/*
void print_point_list( PartitionPoint* point )
{
  PRINT_INFO("  Point %p: (%d curves)\n", static_cast<void*>(point),
    point ? point->num_curves() : 0 );
  if(!point) { PRINT_INFO("\n\n"); return; }
  
  PartitionCurve* curve = 0;
  while( curve = point->next_curve(curve) )
    PRINT_INFO("%10p",curve);
  PRINT_INFO("\n");
  curve = 0;

  while( curve = point->next_curve(curve) )
    PRINT_INFO("%10s",
      curve->start_point() == point && curve->end_point() == point ? "(both)" :
      curve->start_point() == point ? "(start)" :
      curve->end_point() == point ? "(end)" : "(none)" );
  PRINT_INFO("\n");
}
*/

PartitionCurve::PartitionCurve( )
  : firstCoEdge(0),
    startPoint(0), 
    endPoint(0), 
    startNext(0), 
    endNext(0)
{ 
}

PartitionCurve::~PartitionCurve()
{
  start_point(0);
  end_point(0);
  
  remove_all_coedges();
  remove_facet_data();
}

CubitStatus PartitionCurve::add( PartitionCoEdge* coedge )
{
  if( coedge->myCurve )
    return CUBIT_FAILURE;
  
  coedge->curveNext = firstCoEdge;
  firstCoEdge = coedge;
  coedge->myCurve = this;
  return CUBIT_SUCCESS;
}

CubitStatus PartitionCurve::remove( PartitionCoEdge* coedge )
{
  if( coedge->myCurve != this )
    return CUBIT_FAILURE;
  
  if( firstCoEdge == coedge )
  {
    firstCoEdge = coedge->curveNext;
  }
  else
  {
    PartitionCoEdge* prev = firstCoEdge;
    while( prev && prev->curveNext != coedge )
      prev = prev->curveNext;
    
    if( !prev )
    {
      assert(0);
      return CUBIT_FAILURE;
    }
    
    prev->curveNext = coedge->curveNext;
  }
  
  coedge->curveNext = 0;
  coedge->myCurve = 0;
  return CUBIT_SUCCESS;
}

int PartitionCurve::num_coedges() const
{
  int count = 0;
  for( PartitionCoEdge* coe = firstCoEdge; coe; coe = coe->curveNext )
    count++;
  return count;
}

void PartitionCurve::remove_all_coedges( )
{
  while( firstCoEdge && remove(firstCoEdge) ) {;}
  assert( !firstCoEdge );
}

CubitStatus PartitionCurve::start_point( PartitionPoint* point )
{
  if( point == startPoint )
    return CUBIT_SUCCESS;
  
  if( !remove_start_point() )
    { assert(0); return CUBIT_FAILURE; }

  if( !point )
    return CUBIT_SUCCESS;
  
  startPoint = point;
  if( point == endPoint )
  {
    startNext = endNext;
    endNext = 0;
  }
  else
  {
    startNext = point->firstCurve;
    point->firstCurve = this;
    point->curveCount++;
  }

  return CUBIT_SUCCESS;
}

CubitStatus PartitionCurve::end_point( PartitionPoint* point )
{
  if( point == endPoint )
    return CUBIT_SUCCESS;
  
  if( !remove_end_point() )
    { assert(0); return CUBIT_FAILURE; }
  
  if( !point )
    return CUBIT_SUCCESS;

  endPoint = point;
  if( point != startPoint )
  {
    endNext = point->firstCurve;
    point->firstCurve = this;
    point->curveCount++;
  }

  return CUBIT_SUCCESS;
}

CubitStatus PartitionCurve::remove_start_point()
{
  if( startPoint == endPoint )
    endNext = startNext;
  else if( startPoint && !remove_from_point( startPoint, startNext ) )
    return CUBIT_FAILURE;

  startPoint = 0;
  startNext = 0;

  return CUBIT_SUCCESS;
}

CubitStatus PartitionCurve::remove_end_point()
{
  assert( startPoint != endPoint || !endNext );
    // endNext should be null if points are the same

  if( endPoint && endPoint != startPoint &&
      !remove_from_point( endPoint, endNext ) )
    return CUBIT_FAILURE;
  
  endPoint = 0;
  endNext = 0;

  return CUBIT_SUCCESS;
}

CubitStatus PartitionCurve::remove_from_point( PartitionPoint* point, 
                                               PartitionCurve* next )
{
  point->curveCount--;
  if( point->firstCurve == this )
  {
    point->firstCurve = next;
  }
  else
  {
    PartitionCurve* prev = point->firstCurve; 
    while( prev )
    {
      PartitionCurve* temp = prev->next_curve(point);
      if( temp == this )
        break;
      prev = temp;
    }
    
    if( !prev )
    {
      assert(0);
      return CUBIT_FAILURE;
    }
    
    if( prev->startPoint == point )
      prev->startNext = next;
    else if( prev->endPoint == point )
      prev->endNext = next;
    else
      assert(0);
  }
  
  return CUBIT_SUCCESS;
}
    
    
    

bool PartitionCurve::is_nonmanifold( const PartitionSurface* surf ) const
{
  int count = 0;
  PartitionCoEdge* coedge = 0;
  while ( (coedge = next_coedge(coedge)) )
    if ( coedge->get_loop() && coedge->get_loop()->get_surface() == surf )
      count++;
  return count > 1;
}

bool PartitionCurve::is_in_surface( const PartitionSurface* surf,
                                    bool manifold_only ) const
{
  int count = 0;
  PartitionCoEdge* coedge = 0;
  while ( (coedge = next_coedge(coedge)) )
    if ( coedge->get_loop() && coedge->get_loop()->get_surface() == surf )
      count++;
  
  if (manifold_only)
    return count == 1;
  else 
    return count > 0;
}
  

CubitStatus PartitionCurve::move_to_geometry( CubitVector& position )
{
  const CubitVector copy(position);
  return closest_point( copy, position );
}

void PartitionCurve::get_parents_virt( DLIList<TopologyBridge*>& result )
{
  for( PartitionCoEdge* coe = firstCoEdge; coe; coe = coe->curveNext )
    result.append( coe );
}

void PartitionCurve::get_children_virt( DLIList<TopologyBridge*>& result )
{
  result.append( startPoint );
  if( startPoint != endPoint )
    result.append( endPoint );
}

GeometryQueryEngine* PartitionCurve::get_geometry_query_engine() const
{
  return VirtualQueryEngine::instance();
}

void PartitionCurve::print_debug_info( const char* prefix,
                                       bool entset ) const
{
  const bool print_children = true;
  
  if( prefix == 0 ) prefix = "";
  char* new_prefix = new char[strlen(prefix)+3];
  strcpy( new_prefix, prefix );
  strcat( new_prefix, "  ");
  PRINT_INFO("%sPartitionCurve %p\n", prefix, static_cast<void*>(this) );
  
  if(entset)
    sub_entity_set().print_debug_info( new_prefix );
  else
    print_partitioned_entity();
  
  PartitionCoEdge* coedge = 0;
  DLIList<Surface*> surface_list(1);
  while( (coedge = next_coedge(coedge)) )
  {
    surface_list.clean_out();
    coedge->surfaces(surface_list);
    Surface* surf_ptr = surface_list.size() == 1 ? surface_list.get() : 0;
    PRINT_INFO("%s  %s on Surface %p\n", prefix, 
      coedge->sense() == CUBIT_FORWARD ? "FORWARD" :
      coedge->sense() == CUBIT_REVERSED ? "REVERSE" : "UNKNOWN",
      static_cast<void*>(surf_ptr) );
  }
  
  PRINT_INFO("%s  start point: %p\n", prefix,
    static_cast<void*>(start_point()) );
  if( start_point() && print_children )
    start_point()->print_debug_info( new_prefix, false );
  
  if( end_point() && end_point() == start_point() )
    PRINT_INFO("%s  end point SAME as start point\n", prefix ); 
  else
    PRINT_INFO("%s  end point: %p\n", prefix, static_cast<void*>(end_point()) );
  if( end_point() && end_point() != start_point() && print_children )
    end_point()->print_debug_info( new_prefix, false );
  
  delete [] new_prefix;
}

void PartitionCurve::append_simple_attribute_virt(CubitSimpleAttrib* csa)
  { sub_entity_set().add_attribute( this, csa ); }
void PartitionCurve::remove_simple_attribute_virt(CubitSimpleAttrib* csa)
  { sub_entity_set().rem_attribute( this, csa ); }
void PartitionCurve::remove_all_simple_attribute_virt()
  { sub_entity_set().rem_all_attrib( this ); }
CubitStatus PartitionCurve::get_simple_attribute(DLIList<CubitSimpleAttrib*>& list)
{ 
  sub_entity_set().get_attributes( this, list ); 
  return CUBIT_SUCCESS;
}
CubitStatus PartitionCurve::get_simple_attribute(const CubitString& name,
                                       DLIList<CubitSimpleAttrib*>& list)
{ 
  sub_entity_set().get_attributes( this, name.c_str(), list ); 
  return CUBIT_SUCCESS;
}

void PartitionCurve::reverse_point_order()
{
  if( startPoint != endPoint )
  {
    PartitionPoint* tmp_pt = startPoint;
    startPoint = endPoint;
    endPoint = tmp_pt;
    
    PartitionCurve* tmp_cv = startNext;
    startNext = endNext;
    endNext = tmp_cv;
  }
}

void PartitionCurve::get_facet_data( DLIList<CubitFacetEdgeData*>& result_list ) const
  { result_list = facetEdges; }
  
void PartitionCurve::set_facet_data( const DLIList<CubitFacetEdgeData*>& new_list )
{
  remove_facet_data();
  facetEdges = new_list;
  for( int i = facetEdges.size(); i--; )
    TDVGFacetOwner::set( facetEdges.step_and_get(), this );
}

void PartitionCurve::remove_facet_data( )
{
  while( facetEdges.size() )
    TDVGFacetOwner::remove( facetEdges.pop() );
}  


void PartitionCurve::notify_split( FacetEntity* old_entity, FacetEntity* new_entity )
{
  CubitFacetEdgeData* old_ptr = dynamic_cast<CubitFacetEdgeData*>(old_entity);
  CubitFacetEdgeData* new_ptr = dynamic_cast<CubitFacetEdgeData*>(new_entity);
  assert(TDVGFacetOwner::get(old_entity) == this);
  assert(TDVGFacetOwner::get(new_entity) == 0);
  assert(old_ptr && new_ptr);
  
  facetEdges.reset();
  if (facetEdges.get() == old_ptr)
  {
    if (new_ptr->contains(start_point()->facet_point()))
      facetEdges.insert_first(new_ptr);
    else
      facetEdges.insert(new_ptr);
  }
  else if (facetEdges.move_to(old_ptr))
  {
    CubitFacetEdgeData* prev = facetEdges.prev();
    if (prev->shared_point(new_ptr))
      facetEdges.back();
    facetEdges.insert(new_ptr);
  }
  else
  {
    assert( 0 /*old_entity not in list*/);
  }

  facetEdges.reset();
  TDVGFacetOwner::set( new_ptr, this );
}


CubitStatus PartitionCurve::fix_facet_data( PartitionCurve* new_curve )
{
  if ( !facetEdges.size() )
    return CUBIT_SUCCESS;
    
    
  DLIList<CubitFacetEdgeData*> my_edges(facetEdges);
  assert( end_point() == new_curve->start_point() );
  CubitVector pos( end_point()->coordinates() );
  
  double min_edge_len_sqr = CUBIT_DBL_MAX;
  double min_edge_dst_sqr = CUBIT_DBL_MAX;
  CubitFacetEdgeData* closest_edge = 0;
  CubitVector point, vect, edge_pos;
  int i;
  
  for ( i = my_edges.size(); i--; ) {
    CubitFacetEdgeData* edge = my_edges.step_and_get();
    point = edge->point(0)->coordinates();
    vect = edge->point(1)->coordinates() - point;
    
    double len_sqr = vect.length_squared();
    if ( len_sqr < min_edge_len_sqr )
      min_edge_len_sqr = len_sqr;
    
    point -= pos;  
    double t = (vect % -point) / len_sqr;
    if ( t < 0.0 ) t = 0.0;
    else if( t > 1.0 ) t = 1.0;
    double dst_sqr = (point + t * vect).length_squared();

    if ( dst_sqr < min_edge_dst_sqr ) {
      min_edge_dst_sqr = dst_sqr;
      closest_edge = edge;
      edge_pos = edge->point(0)->coordinates() + t * vect;
    }
  }
  
  my_edges.reset();
  DLIList<CubitFacetEdgeData*> new_curve_edges;
  my_edges.last();
  while( my_edges.get() != closest_edge ) {
    new_curve_edges.append( my_edges.pop() );
    my_edges.last();
  }
  
  min_edge_len_sqr *= 0.2;
  double dst_tol_sqr = CUBIT_MAX( min_edge_len_sqr, GEOMETRY_RESABS*GEOMETRY_RESABS );
  
  bool edge_reversed = false;
  if ( my_edges.size() > 1 )
  {
    my_edges.last();
    CubitFacetEdgeData* other = my_edges.prev();
    edge_reversed = other->contains( closest_edge->point(1) );
  }
  else if( new_curve_edges.size() > 1 )
  {
    new_curve_edges.last();
    CubitFacetEdgeData* other = new_curve_edges.get();
    edge_reversed = other->contains( closest_edge->point(0) );
  }
  else
  {
    edge_reversed = (start_point()->facet_point() == closest_edge->point(1));
    assert(edge_reversed || start_point()->facet_point() == closest_edge->point(0));
  }
  
  
  CubitPoint* new_point = 0;
  double d1 = (pos - closest_edge->point(0)->coordinates()).length_squared();
  double d2 = (pos - closest_edge->point(1)->coordinates()).length_squared();
  if ( (!TDVGFacetOwner::get(closest_edge->point(0)) && d1 < dst_tol_sqr) ||
       (d1 < GEOMETRY_RESABS*GEOMETRY_RESABS) )
  {
    if ( !edge_reversed )
      new_curve_edges.append( my_edges.pop() );
    new_point = closest_edge->point(0);
  }
  else if( (!TDVGFacetOwner::get(closest_edge->point(1)) && d2 < dst_tol_sqr) ||
      (d2 < GEOMETRY_RESABS*GEOMETRY_RESABS) )
  {
    if( edge_reversed )
      new_curve_edges.append( my_edges.pop() );
    new_point = closest_edge->point(1);
  }
  else
  {
    CubitFacetEdge* new_edge;
    CubitFacet* new_facet;
    TDVGFacetOwner::remove(closest_edge);
    PartSurfFacetTool::split_edge(closest_edge, edge_pos, 0, new_point, new_edge, new_facet );
    TDVGFacetOwner::set(closest_edge, this);

    if ( edge_reversed ) {
      new_curve_edges.append( my_edges.pop() );
      my_edges.append( dynamic_cast<CubitFacetEdgeData*>(new_edge) );
    } else {
      new_curve_edges.append( dynamic_cast<CubitFacetEdgeData*>(new_edge) );
      //TDVGFacetOwner::remove(new_edge);
    }
  }

  CubitPointData* new_point_d = dynamic_cast<CubitPointData*>(new_point);
  new_point_d->set(end_point()->coordinates());
  end_point()->facet_point( new_point_d );
  set_facet_data(my_edges);
  new_curve_edges.reverse();
  new_curve->set_facet_data(new_curve_edges);
  
  PartitionCoEdge* coedge = 0;
  while ((coedge = next_coedge(coedge)) != NULL )
  {
    PartitionLoop* loop = coedge->get_loop();
    if (!loop) continue;
    
    PartitionSurface* surf = loop->get_surface();
    if (!surf->notify_moving_point( new_point, pos ))
      return CUBIT_FAILURE;
  }

  new_point->set(pos);
  return CUBIT_SUCCESS;
}

void PartitionCurve::remove_dead_facet( CubitFacetEdgeData* edge )
{
  facetEdges.remove(edge);
}
  
CubitStatus PartitionCurve::get_save_topology( DLIList<int>& end_points )
{
  for( int i = 0; i < 2; i++ )
  {
    PartitionPoint* point = i ? end_point() : start_point();
    int set_id, pt_id;
    if( &(point->sub_entity_set()) == &sub_entity_set() )
      set_id = 0;
    else 
      set_id = point->sub_entity_set().get_unique_id();
    pt_id = point->sub_entity_set().get_id(point);
    end_points.append(set_id);
    end_points.append(pt_id);
  }
  return CUBIT_SUCCESS;
}



void PartitionCurve::transform(const CubitTransformMatrix& xform)
{
  int i;
  
  DLIList<CubitPoint*> points(facetEdges.size());
  for ( i = facetEdges.size(); i--; )
  {
    CubitFacetEdgeData* edge = facetEdges.step_and_get();
    edge->point(0)->marked(1);
    edge->point(1)->marked(1);
  }
  for ( i = facetEdges.size(); i--; )
  {
    CubitFacetEdgeData* edge = facetEdges.step_and_get();
    for ( int j = 0; j < 2; j++ )
    {
      CubitPoint* pt = edge->point(j);
      if( pt->marked() )
      {
        pt->marked(0);
      
        if( TDVGFacetOwner::get(pt) == 0 )
          points.append(pt);
      }
    }
  }
   
  for( i = points.size(); i--; )
  {
    CubitPoint* pt = points.get_and_step();
    pt->set( xform * pt->coordinates() );
  }   
}

void PartitionCurve::replace_facet(CubitFacetEdgeData* dead_facet,
                                   DLIList<CubitFacetEdgeData*> &new_facets)
{
  assert(new_facets.size() > 1); // doesn't make sense to replace 1 to 1

  DLIList<CubitFacetEdgeData*> new_copy = new_facets;

  new_facets.reset();
  facetEdges.reset();
  int dead_index = facetEdges.where_is_item(dead_facet);
  assert(-1 != dead_index);

  // special case for a single facet on the edge
  // -  use curve endpoints to determine order
  if (facetEdges.size() == 1)
  {
    facetEdges.clean_out();

    new_copy.last();
    if (new_copy.get()->contains(start_point()->facet_point()) )
      new_copy.reverse();

    new_copy.reset();
    assert( new_copy.get()->contains(start_point()->facet_point()) );
    new_copy.last();
    assert( new_copy.get()->contains(end_point()->facet_point()) );

    facetEdges = new_copy;
  }
  else
  {
    DLIList<CubitFacetEdgeData*> ordered_edges;
    int i;

    // copy the edges before the dead edge into a new list
    for (i=0; i< dead_index; i++)
      ordered_edges.append(facetEdges.next(i));

    // if there are edges before the new edges, find out whether the first or last
    // new edge is adjacent to the one before the dead edge
    CubitFacetEdgeData *prev = NULL;
    CubitFacetEdgeData *next = NULL;
    CubitPoint *new_start_pt = NULL;
    CubitPoint *new_end_pt = NULL;
    if (dead_index > 0)
    {
      prev = facetEdges.next(dead_index-1);
      new_start_pt = prev->shared_point(dead_facet);
      assert(new_start_pt != NULL);
    }

    if (dead_index < (facetEdges.size() - 1) )
    {
      next = facetEdges.next(dead_index+1);
      new_end_pt = next->shared_point(dead_facet);
      assert(new_end_pt != NULL);
    }

    if (prev)
    {
      new_copy.last();
      if (new_copy.get()->contains(new_start_pt))
        new_copy.reverse();
    }
    else
    {
      // use the edge after the dead edge to find out whether to reverse the new
      // edges or not
      assert(next != NULL);
      new_copy.reset();
      if (new_copy.get()->contains(new_end_pt))
        new_copy.reverse();
    }

    if (prev)
    {
      new_copy.reset();
      assert( prev->shared_point(new_copy.get()) != NULL);
    }
    if (next)
    {
      new_copy.last();
      assert( next->shared_point(new_copy.get()) != NULL);
    }

    // copy the new edges into place in the ordered list
    new_copy.reset();
    for (i=new_copy.size(); i>0; i--)
      ordered_edges.append(new_copy.get_and_step());

    // add the rest of the edges to the list
    facetEdges.reset();
    for (i=dead_index+1; i<facetEdges.size(); i++)
      ordered_edges.append(facetEdges.next(i));


    facetEdges.clean_out();
    facetEdges = ordered_edges;
  }

  TDVGFacetOwner::remove(dead_facet);
  for( int j = 0; j < new_facets.size(); j++ )
    TDVGFacetOwner::set(new_facets.step_and_get(),this);

  // TODO - memory management - where should the facets get deleted
}

void PartitionCurve::draw_facets( int color)
{
  if( !facetEdges.size() )
  {
    PRINT_ERROR("No facet data.\n");
    return ;
  }
  
  int i;
  GPoint* pts = new GPoint[facetEdges.size() + 1];
  
  facetEdges.reset();
  CubitPoint* pt = 0;
  if( facetEdges.size() == 1 )
  {
    pt = facetEdges.get()->point(0);
  }
  else
  {
    pt = facetEdges.get()->shared_point(facetEdges.next());
    pt = facetEdges.get()->other_point(pt);
  }
  
  CubitVector v = pt->coordinates();
  pts[0].x = (float)v.x();
  pts[0].y = (float)v.y();
  pts[0].z = (float)v.z();

  for( i = 1; i <= facetEdges.size(); i++ )
  {
    pt = facetEdges.get_and_step()->other_point(pt);
    CubitVector v = pt->coordinates();
    pts[i].x = (float)v.x();
    pts[i].y = (float)v.y();
    pts[i].z = (float)v.z();
  }
  
  GfxDebug::draw_polyline( pts, facetEdges.size() + 1, color );
  for ( i = 0; i <= facetEdges.size(); i++ )
  {
    GfxDebug::draw_point( pts[i].x, pts[i].y, pts[i].z, color );
  }
  
  facetEdges.reset();
  for ( i = 0; i < facetEdges.size(); i++ )
  {
    CubitFacetEdge* edge = facetEdges.get_and_step();
    CubitVector mid = 0.5*(edge->point(0)->coordinates()+edge->point(1)->coordinates());
    GfxDebug::draw_label( i, (float)mid.x(), (float)mid.y(), (float)mid.z(), color );
  }
  GfxDebug::flush();
  
  delete [] pts;
}

void PartitionCurve::do_facet_cleanup()
{
  int i;
  const double tol = 10.0 * GEOMETRY_RESABS;
  const double tol_sqr = tol * tol;
  
  DLIList<CubitFacetEdgeData*> edge_list(facetEdges);
  CubitFacetEdgeData* prev_edge = 0;
  edge_list.reset();
  
  for ( i = edge_list.size(); i--; )
  {
    CubitFacetEdgeData* edge = edge_list.get_and_step();
    PartitionEntity* start_pt = TDVGFacetOwner::get(edge->point(0));
    PartitionEntity* end_pt = TDVGFacetOwner::get(edge->point(1));
    if ( start_pt && end_pt )
    {
      prev_edge = edge;
      continue;
    }
    
    CubitVector edge_vect(edge->point(1)->coordinates());
    edge_vect -= edge->point(0)->coordinates();
    if ( edge_vect.length_squared() > tol_sqr )
    {
      prev_edge = edge;
      continue;
    }
    
    CubitPoint *keep, *dead;
    if ( start_pt ) 
    {
      keep = edge->point(0);
      dead = edge->point(1);
    }
    else if ( end_pt )
    {
      keep = edge->point(1);
      dead = edge->point(0);
    }
    else
    {
      double prev = (prev_edge->point(0)->coordinates() -
        prev_edge->point(1)->coordinates()).length_squared();
      double next = (edge_list.get()->point(0)->coordinates() -
        edge_list.get()->point(1)->coordinates()).length_squared();
      if ( prev < next )
      {
        dead = prev_edge->shared_point(edge);
        keep = edge_list.get()->shared_point(edge);
      }
      else
      {
        dead = edge_list.get()->shared_point(edge);
        keep = prev_edge->shared_point(edge);
      }
      assert( dead && keep );
    }
    
    if ( PartSurfFacetTool::collapse_edge( keep, dead ) ) 
    {
      ;
    }
    else
    {
      prev_edge = edge;
    }
  }
}


CubitStatus PartitionCurve::move_to_geometry( CubitPoint* point )
{
  CubitVector new_pos;
  if (!closest_point( point->coordinates(), new_pos ))
    return CUBIT_FAILURE;
    
  point->set(new_pos);
  
  DLIList<PartitionSurface*> surfaces;
  PartitionCoEdge* coedge = 0;
  while ((coedge = next_coedge(coedge)) != NULL )
  {
    PartitionLoop* loop = coedge->get_loop();
    if (!loop)
      continue;
    
    PartitionSurface* surf = loop->get_surface();
    surfaces.append_unique(surf);
    
    coedge = coedge->next();
  }
  
  int i;
  surfaces.reset();
  for ( i = 0; i < surfaces.size(); i++ )
    if (!surfaces.get_and_step()->notify_moving_point( point, new_pos ))
      break;
      
  if (i < surfaces.size()) 
    return CUBIT_FAILURE;
  
  point->set(new_pos);
  return CUBIT_SUCCESS;
}
