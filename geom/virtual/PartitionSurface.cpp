#include "PartitionLump.hpp"
#include "PartitionSurface.hpp"
#include "PartitionLoop.hpp"
#include "PartitionShell.hpp"
#include "PartitionCurve.hpp"
#include "VirtualQueryEngine.hpp"
#include "PartitionEngine.hpp"
#include "PartSurfFacetTool.hpp"

#include "CubitFacetData.hpp"
#include "TDVGFacetOwner.hpp"
#include "GMem.hpp"
#include "CubitPoint.hpp"
#include "CubitPointData.hpp"
#include "CubitFacetEdgeData.hpp"
#include "CubitTransformMatrix.hpp"
#include "GfxDebug.hpp"

#include "RefVertex.hpp"

void PartitionSurface::draw_facets( int color ) const
{
  DLIList<CubitFacetData*> facets;
  get_facet_data(facets);
  int i;
  for ( i = facets.size(); i--; ) {
    CubitFacet* facet = facets.step_and_get();
    for ( int j = 0; j < 3; j++ )
      if ( facet->edge(j) )
        facet->edge(j)->marked(0);
  }
  for ( i = facets.size(); i--; ) {
    CubitFacet* facet = facets.step_and_get();
    for ( int j = 0; j < 3; j++ )
      if ( facet->edge(j) )
        facet->edge(j)->marked(facet->edge(j)->marked()+1);
  }
  
  
  for ( i = facets.size(); i--; ) {
    CubitFacetData* facet = facets.step_and_get();
    CubitVector pts[3] = {facet->point(0)->coordinates(),
                          facet->point(1)->coordinates(),
                          facet->point(2)->coordinates()};
    for ( int j = 0; j < 3; j++ )
    {
      CubitVector p1 = pts[j];
      CubitVector p2 = pts[(j+1)%3];
      int c = color == CUBIT_RED ? color - 1 : CUBIT_RED;
      CubitFacetEdge* edge = facet->edge((j+2)%3);
      if (!edge) continue;
      
      switch( edge->marked() ) {
          case 1: c = color + 1; break;
          case 2: c = color; break;
          default : c = color == CUBIT_RED ? color - 1 : CUBIT_RED;
       
        if (edge->marked()==1 && !TDVGFacetOwner::get(edge))
          PRINT_WARNING("Boundary edge in surface facetting not owned by a curve.\n");
      }
    
      GfxDebug::draw_line( (float)p1.x(), (float)p1.y(), (float)p1.z(), 
                      (float)p2.x(), (float)p2.y(), (float)p2.z(), c );
    }
    
    
    CubitVector center = ( pts[0] + pts[1] + pts[2] ) / 3.0;
    CubitVector normal = (pts[1] - pts[0]) * (pts[2] - pts[0]);
    double len = normal.length();
    if ( len > GEOMETRY_RESABS ) {
        normal /= sqrt(len);
//      if ( len > 1e-2 )
//        normal /= len;
//      else
//        normal.length(1.0);

        GfxDebug::draw_vector( center, center + normal, color );
    } else {
      GfxDebug::draw_point( center, color + 1 );
    }
  }
  GfxDebug::flush();

  for ( i = facets.size(); i--; ) {
    CubitFacet* facet = facets.step_and_get();
    for ( int j = 0; j < 3; j++ )
      if ( facet->edge(j) )
        facet->edge(j)->marked(0);
  }
}

void PartitionSurface::print_debug_info( const char* prefix ,
                                         bool print_sub_entity_set ) const
{
  if( !prefix ) prefix = "";
  PRINT_INFO("%sPartitionSurface %p\n", prefix, this );
  PartitionLoop* loop = 0;
  while( (loop = next_loop(loop)) )
  {
    PRINT_INFO("%s  Loop %p:\n", prefix, loop );
    PartitionCoEdge* coedge = loop->first_coedge();
    do
    {
      PRINT_INFO("%s    CoEdge %p %s -> Curve %p\n",
        prefix, coedge, coedge->sense() == CUBIT_FORWARD ? "FORWARD" :
                        coedge->sense() == CUBIT_REVERSED ? "REVERSE" :
                        "UNKNOWN", coedge->get_curve() );
      coedge = loop->next_coedge(coedge);
    } while( coedge != loop->first_coedge() );
  }

  char buffer[128];
  sprintf(buffer,"%s  ",prefix);
  if( print_sub_entity_set )
    sub_entity_set().print_debug_info(buffer);
} 

PartitionSurface::PartitionSurface()
  : firstLoop(0), firstCoSurf(0), geometry_sense(CUBIT_FORWARD)
{}

PartitionSurface::PartitionSurface( PartitionLump* lump )
  : firstLoop(0), firstCoSurf(0), geometry_sense(CUBIT_FORWARD)
{
  lump->sub_entity_set().add_lower_order( this );
}

PartitionSurface::~PartitionSurface()
{
  while( firstLoop )
    remove( firstLoop );
  
  while( firstCoSurf )
  {
    PartitionCoSurf* cos = firstCoSurf;
    if( cos->get_shell() )
      cos->get_shell()->remove( cos );
    remove( cos );
    delete cos;
    assert( firstCoSurf != cos );
  }
  
    // delete facets
  DLIList<CubitPoint*> facet_points(3);
  DLIList<CubitFacetEdge*> facet_edges(3);
  while( facetList.size() )
  {
    CubitFacet* facet = facetList.pop();
    facet_edges.clean_out();
    facet->edges(facet_edges);
    facet_points.clean_out();
    facet->points(facet_points);
    TDVGFacetOwner::remove(facet);
    delete facet;
    
    PartitionCurve* curve;
    while( facet_edges.size() )
    {
      CubitFacetEdge* edge = facet_edges.pop();
      if( edge && !edge->num_adj_facets() ) {
        if ( TDVGFacetOwner::get(edge) ) {
          curve = dynamic_cast<PartitionCurve*>(TDVGFacetOwner::get(edge));
          curve->remove_facet_data();
        }
        delete edge;
      }
    }
    
    PartitionPoint* ppoint;
    while( facet_points.size() )
    {
      CubitPoint* point = facet_points.pop();
      if( point && !point->num_adj_facets() ) {
        if ( TDVGFacetOwner::get(point) ) {
          ppoint = dynamic_cast<PartitionPoint*>(TDVGFacetOwner::get(point));
          if( ppoint )
            ppoint->facet_point(0);
        }
        delete point;
      }
    }
  }
}

CubitStatus PartitionSurface::add( PartitionLoop* loop )
{
  if( loop->mySurface )
    return CUBIT_FAILURE;
  
  loop->mySurface = this;
  loop->nextInSurface = firstLoop;
  firstLoop = loop;
  return CUBIT_SUCCESS;
}
  
  
CubitStatus PartitionSurface::remove( PartitionLoop* loop )
{
  if( loop->mySurface != this )
    return CUBIT_FAILURE;
  
  if( firstLoop == loop )
  {
    firstLoop = loop->nextInSurface;
  }
  else
  {
    PartitionLoop* prev = firstLoop;
    while( prev && prev->nextInSurface != loop )
      prev = prev->nextInSurface;
    
    if( !prev )
    {
      assert(0);
      return CUBIT_FAILURE;
    }
    
    prev->nextInSurface = loop->nextInSurface;
  }
  
  loop->mySurface = 0;
  loop->nextInSurface = 0;
  return CUBIT_SUCCESS;
}

int PartitionSurface::num_loops() const
{
  int count = 0;
  for( PartitionLoop* loop = firstLoop; loop; loop = loop->nextInSurface )
    count++;
  return count;
}

CubitStatus PartitionSurface::add( PartitionCoSurf* cosurf )
{
  if( cosurf->mySurface )
    return CUBIT_FAILURE;
  
  cosurf->mySurface = this;
  cosurf->surfaceNext = firstCoSurf;
  firstCoSurf = cosurf;
  return CUBIT_SUCCESS;
}

CubitStatus PartitionSurface::remove( PartitionCoSurf* cosurf )
{
  if( cosurf->mySurface != this )
    return CUBIT_FAILURE;
  
  if( cosurf == firstCoSurf )
  {
    firstCoSurf = cosurf->surfaceNext;
  }
  else
  {
    PartitionCoSurf* prev = firstCoSurf;
    while( prev && prev->surfaceNext != cosurf )
      prev = prev->surfaceNext;
    
    if( !prev ) { assert(prev != NULL); return CUBIT_FAILURE; }
    
    prev->surfaceNext = cosurf->surfaceNext;
  }
  
  cosurf->mySurface = 0;
  cosurf->surfaceNext = 0;
  return CUBIT_SUCCESS;
}

PartitionCoSurf* PartitionSurface::find_first( const PartitionShell* shell ) const
{
  PartitionCoSurf* cos = firstCoSurf;
  while( cos && cos->get_shell() != shell )
    cos = cos->surfaceNext;
  return cos;
}


PartitionCoSurf* PartitionSurface::find_first( const PartitionLump* lump ) const
{
  PartitionCoSurf* cos = firstCoSurf;
  while( cos && (!cos->get_shell() || cos->get_shell()->get_lump() != lump ) )
    cos = cos->surfaceNext;
  return cos;
}

PartitionCoSurf* PartitionSurface::find_next( const PartitionCoSurf* prev ) const
{
  if( prev->mySurface != this )
    return 0;
  PartitionCoSurf* cos = prev->surfaceNext;
  while( cos && cos->get_shell() != prev->get_shell() )
    cos = cos->surfaceNext;
  return cos;
}

PartitionCoSurf* PartitionSurface::find_first( const PartitionShell* shell,
                                               CubitSense sense ) const
{
  PartitionCoSurf* cos = firstCoSurf;
  while( cos && (cos->get_shell() != shell || cos->sense() != sense) )
    cos = cos->surfaceNext;
  return cos;
}

//-------------------------------------------------------------------------
// Purpose       : Get child point bridges
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 
//-------------------------------------------------------------------------
void PartitionSurface::get_points( DLIList<PartitionPoint*>& list ) const
{
  PartitionLoop* loop = 0;
  while ( (loop = next_loop(loop)) ) 
  {
    PartitionCoEdge* coedge = loop->first_coedge();
    do 
    {
      PartitionCurve* curve = coedge->get_curve();

      list.append( curve->start_point() );
      list.append( curve->end_point() );
      
      coedge = loop->next_coedge(coedge);
    } while( coedge != loop->first_coedge() );
  }
  
  list.reset();
  int i;
  for ( i = list.size(); i--; )
    list.get_and_step()->mark = 1;
  
  for ( i = list.size(); i--; )
  {
    list.back();
    if( list.get()->mark )
      list.get()->mark = 0;
    else  
      list.change_to(0);
  }
  list.remove_all_with_value(0);
}


    

GeometryQueryEngine* PartitionSurface::get_geometry_query_engine() const
{
  return VirtualQueryEngine::instance();
}

void PartitionSurface::get_children_virt( DLIList<TopologyBridge*>& list )
{
  for( PartitionLoop* loop = firstLoop; loop; loop = loop->nextInSurface )
    list.append( loop );
}

void PartitionSurface::get_parents_virt( DLIList<TopologyBridge*>& result_list )
{
  Surface* real_surf = dynamic_cast<Surface*>(partitioned_entity());
  if( ! real_surf )
  { 
    for( PartitionCoSurf* cos = firstCoSurf; cos; cos = cos->surfaceNext )
      if( cos->get_shell() )
        result_list.append_unique( cos->get_shell() );
  }
  else
  {
    int i;
    DLIList<TopologyBridge*> real_surf_shells, tmp_list;
    real_surf->get_parents_virt( real_surf_shells );
    
      // get real volumes from real shells
    DLIList<TopologyBridge*> real_surf_vols(real_surf_shells.size());
    real_surf_shells.reset();
    for( i = real_surf_shells.size(); i--; )
    {
      tmp_list.clean_out();
      real_surf_shells.get_and_step()->get_parents_virt( tmp_list );
      assert(tmp_list.size() == 1);
      real_surf_vols.append(tmp_list.get());
    }

      // replace real volumes with owning partition volumes (if any)
    DLIList<TopologyBridge*> vol_list(real_surf_vols.size());
    real_surf_vols.reset();
    for( i = real_surf_vols.size(); i--; )
    {
      TopologyBridge* vol_bridge = real_surf_vols.get_and_step();
      SubEntitySet* set = dynamic_cast<SubEntitySet*>(vol_bridge->owner());
      if( set )
      {
        tmp_list.clean_out();
        set->get_owners(tmp_list);
        vol_list += tmp_list;
      }
      else
        vol_list.append(vol_bridge);
    }
    
      // for each volume, get all child shells that are parents of this
    vol_list.reset();
    DLIList<TopologyBridge*> vol_shells;
    for( i = vol_list.size(); i--; )
    {
      vol_shells.clean_out();
      vol_list.get_and_step()->get_children( vol_shells, false, layer() );
      vol_shells.reset();
      for( int j = vol_shells.size(); j--; )
      {
        TopologyBridge* shell = vol_shells.get_and_step();
        tmp_list.clean_out();
        shell->get_children( tmp_list, false, layer() );
        if( tmp_list.is_in_list(this) )
          result_list.append(shell);
      }
    }
  }
}

CubitSense PartitionSurface::get_shell_sense( ShellSM* shell_ptr ) const
{
  if( PartitionShell* pshell = dynamic_cast<PartitionShell*>(shell_ptr) )
    return pshell->find_sense( this );
  
  Surface* real_surf = dynamic_cast<Surface*>(partitioned_entity());
  if( real_surf )
  { 
    DLIList<TopologyBridge*> list(2);
    real_surf->get_parents_virt( list );
    if( list.is_in_list( shell_ptr ) )
      return real_surf->get_shell_sense( shell_ptr );
  }
  
  return CUBIT_UNKNOWN;
}


void PartitionSurface::append_simple_attribute_virt(CubitSimpleAttrib* csa)
{ sub_entity_set().add_attribute( this, csa ); }
void PartitionSurface::remove_simple_attribute_virt(CubitSimpleAttrib* csa)
{ sub_entity_set().rem_attribute( this, csa ); }
void PartitionSurface::remove_all_simple_attribute_virt()
{ sub_entity_set().rem_all_attrib( this ); }
CubitStatus PartitionSurface::get_simple_attribute(DLIList<CubitSimpleAttrib*>& list)
{ 
  sub_entity_set().get_attributes( this, list ); 
  return CUBIT_SUCCESS;
}
CubitStatus PartitionSurface::get_simple_attribute(const CubitString& name,
                                       DLIList<CubitSimpleAttrib*>& list)
{ 
  sub_entity_set().get_attributes( this, name.c_str(), list ); 
  return CUBIT_SUCCESS;
}

void PartitionSurface::reverse_loops()
{
  PartitionLoop* loop = 0;
  while( (loop = next_loop(loop)) )
    loop->reverse();
}

CubitBox PartitionSurface::bounding_box() const
{
  int i, j;
  CubitFacet* facet;
  for ( i = 0; i< facetList.size(); i++ ) {
    facet = facetList.next(i);
    facet->point(0)->marked(1);
    facet->point(1)->marked(1);
    facet->point(2)->marked(1);
  }
  
  facet = facetList.get();
  CubitBox result(facet->point(0)->coordinates());
  facet->point(0)->marked(0);
  
  for ( i = 0; i< facetList.size(); i++ ) {
    facet = facetList.next(i);
    for ( j = 0; j < 3; j++ ) {
      if ( facet->point(j)->marked() ) {
        facet->point(j)->marked(0);
        result |= facet->point(j)->coordinates();
      }
    }
  }
  
  return result;
}

  
  

double PartitionSurface::measure()
{ 
  double result = 0.;
  for ( int i = facetList.size(); i--; )
    result += facetList.step_and_get()->area();
  return result;
}

GeometryType PartitionSurface::geometry_type()
{
  return FACET_SURFACE_TYPE;
}

void PartitionSurface::closest_point_trimmed( CubitVector f, CubitVector& r )
{ 
  closest_facet( f, r );
}

CubitStatus PartitionSurface::get_point_normal( CubitVector &, CubitVector & )
{
  return CUBIT_FAILURE;
}

CubitStatus PartitionSurface::move_to_geometry( CubitVector& pos )
{
  const CubitVector copy(pos);
  return closest_point( copy, &pos );
}

CubitStatus PartitionSurface::closest_point_uv_guess(  
    CubitVector const& location,
    double&, double&,
    CubitVector* closest_location,
    CubitVector* unit_normal )
{
  return closest_point(location, closest_location, unit_normal);
}


CubitStatus PartitionSurface::closest_point( const CubitVector& position,
                                             CubitVector* closest,
                                             CubitVector* normal,
                                             CubitVector* curvature1,
                                             CubitVector* curvature2 )
{
  if ( closest || normal )
  {
    CubitVector tmp_closest;
    CubitFacet* facet = closest_facet( position, tmp_closest );
    if ( closest )
      *closest = tmp_closest;
    if ( normal )
      *normal = facet->normal();
  }
  
  if( curvature1 )
    curvature1->set( 0.0, 0.0, 0.0 );
  if( curvature2 )
    curvature2->set( 0.0, 0.0, 0.0 );
  
  return CUBIT_SUCCESS;
}

CubitStatus PartitionSurface::principal_curvatures( 
                                       const CubitVector& ,
                                       double&, double& ,
                                       CubitVector* )
{
  return CUBIT_FAILURE;
}

CubitVector PartitionSurface::position_from_u_v( double , double  )
{
  return CubitVector(0.0,0.0,0.0);
}

CubitStatus PartitionSurface::u_v_from_position( const CubitVector& ,
                                                 double&, double &,
                                                 CubitVector* )
{
  return CUBIT_FAILURE;
}

CubitBoolean PartitionSurface::is_periodic()
{
  return CUBIT_FALSE;
}

CubitBoolean PartitionSurface::is_periodic_in_U( double& )
{
  return  CUBIT_FALSE;
}

CubitBoolean PartitionSurface::is_periodic_in_V( double& )
{
  return CUBIT_FALSE;
}

CubitBoolean PartitionSurface::is_singular_in_U( double )
{
  return CUBIT_FALSE;
}

CubitBoolean PartitionSurface::is_singular_in_V( double )
{
  return  CUBIT_FALSE;
}

CubitBoolean PartitionSurface::is_closed_in_U()
{
  return CUBIT_FALSE;
}

CubitBoolean PartitionSurface::is_closed_in_V()
{
  return  CUBIT_FALSE;
}

CubitStatus PartitionSurface::uv_derivitives( double , double , 
                                              CubitVector& ,
                                              CubitVector&  )
{
  return CUBIT_FAILURE;
}

CubitBoolean PartitionSurface::is_parametric()
{
  return CUBIT_FALSE;
}

CubitBoolean PartitionSurface::get_param_range_U( double& , double&  )
{
  return CUBIT_FALSE;
}

CubitBoolean PartitionSurface::get_param_range_V( double&, double& )
{
  return CUBIT_FALSE;
}

CubitBoolean PartitionSurface::is_position_on( CubitVector& position )
{
  const double tolsqr = GEOMETRY_RESABS * GEOMETRY_RESABS;
  CubitVector closest;
  closest_facet( position, closest );
  double distsqr = (position - closest).length_squared();
  return (CubitBoolean)(distsqr < tolsqr);
}

CubitPointContainment PartitionSurface::point_containment( const CubitVector& point )
{
  PartitionCurve* curve = 0;
  CubitPointContainment result = point_containment( point, curve );
  if ( curve ) {
    CubitVector curve_closest;
    curve->closest_point( point, curve_closest );
    if ( (point - curve_closest).length_squared() < GEOMETRY_RESABS*GEOMETRY_RESABS )
      result = CUBIT_PNT_BOUNDARY;
  }
  return result;
}

CubitPointContainment PartitionSurface::point_containment( const CubitVector& point,
                                                       PartitionCurve*& curve )
{
  CubitVector facet_pos;
  CubitFacet* facet = closest_facet(point, facet_pos);
  
  int i;
  double tol = 10*GEOMETRY_RESABS;
  int closest_edge = -1;
  CubitPoint* closest_point = 0;
  double closest_dist_sqr = CUBIT_DBL_MAX;
  CubitVector closest_pos;
  CubitFacetEdge* edge;
  
    // Find closest edge/point on facet
  for ( i = 0; i < 3; i++ ) 
  {
    edge = facet->edge(i);
    
    CubitVector start = edge->point(0)->coordinates();
    CubitVector end   = edge->point(1)->coordinates();
    CubitVector dir = end - start;
    double len_sqr = dir.length_squared();
    double t = (dir % (point - start)) / len_sqr;
    double len = sqrt(len_sqr);
    double dist_from_start = t * len;
    double dist_from_end   = (1.0 - t) * len;
    CubitPoint* facet_point = 0;
    CubitVector pos;
    if ( dist_from_start < tol ) 
    {
      pos = start;
      facet_point = edge->point(0);
    }
    else if( dist_from_end < tol ) 
    {
      pos = end;
      facet_point = edge->point(1);
    }
    else 
      pos = t * dir + start;
     
    double dist_sqr = (pos - point).length_squared();
    if ( dist_sqr < closest_dist_sqr ) 
    {
      closest_edge = i;
      closest_point = facet_point;
      closest_dist_sqr = dist_sqr;
      closest_pos = pos;
    }
  }
  
  edge = facet->edge(closest_edge);
  
    // If closest to an edge (not a point) ...
  if ( !closest_point ) 
  {
      // If the edge is not on the boundary of the surface facets...
    curve = dynamic_cast<PartitionCurve*>(TDVGFacetOwner::get(edge));
    if ( !curve || curve->is_nonmanifold(this) )
      return CUBIT_PNT_INSIDE;  // not a boundary edge
  
      // If within tolerance of boundary
    if ( (point - closest_pos).length_squared() < tol*tol )
      return CUBIT_PNT_BOUNDARY;
    
      // Check which side of the edge the point is on
    CubitVector dir = edge->point(1)->coordinates() - edge->point(0)->coordinates();
    if ( facet->edge_use(closest_edge) == -1 )
      dir = -dir;
    double dot = facet->normal() % ((closest_pos - point) * dir);
    if ( dot < 0.0 )
      return CUBIT_PNT_OUTSIDE;
    else
      return CUBIT_PNT_INSIDE;
  }
  
  
    // Get two boundary edges adjacent to point
  DLIList<CubitFacetEdge*> edge_list;
  closest_point->edges( edge_list );
  for ( i = edge_list.size(); i--; )
  {
    edge = edge_list.step_and_get();
    curve = dynamic_cast<PartitionCurve*>(TDVGFacetOwner::get(edge));
    if ( !curve || !curve->is_in_surface(this,true) )
      edge_list.change_to(0);
  }
  edge_list.remove_all_with_value(0);
  
  if (!edge_list.size())
    return CUBIT_PNT_INSIDE;

  if (edge_list.size() > 2) 
  {
    // This is a special case.  The point was closest to an
    // intersection of more than two curves on the surface boundary
    // (for example the point where a circular hole is tangent to
    //  the boundary of the surface.).  If the input position is closest
    //  to such a point it must be outside the surface (or on the boundary).
    // However, a little work is required to decide which curve to pass back.
    DLIList<PartitionCurve*> curve_list;
    while (edge_list.size())
    {
      PartitionEntity* owner = TDVGFacetOwner::get(edge_list.pop());
      curve_list.append_unique( dynamic_cast<PartitionCurve*>(owner) );
    }
    assert( !! TDVGFacetOwner::get(closest_point) );
    closest_dist_sqr = CUBIT_DBL_MAX;
    curve = 0;
    while (curve_list.size())
    {
      PartitionCurve* temp_curve = curve_list.pop();
      temp_curve->closest_point( point, closest_pos );
      if ( (point - closest_pos).length_squared() < closest_dist_sqr )
      {
        curve = temp_curve;
        closest_dist_sqr = (point - closest_pos).length_squared();
      }
    }
    
    if ( (point - closest_point->coordinates()).length_squared() < tol*tol )
      return CUBIT_PNT_BOUNDARY;
    else
      return CUBIT_PNT_OUTSIDE;
  }
    
     
  
  assert(edge_list.size() == 2);
  CubitFacetEdge *bdy_edges[2] = {edge_list.get(), edge_list.next()};
  
  
    // need to pass back a curve
  curve = dynamic_cast<PartitionCurve*>(TDVGFacetOwner::get(bdy_edges[0]));
  
    // If within tolerance of boundary
  if ( (point - closest_pos).length_squared() < tol*tol )
    return CUBIT_PNT_BOUNDARY;

    // Fill arrays for each of two boundary edges
  CubitFacet* bdy_facets[2] = {0,0};      // adjacent facet
  bool inside_facets[2] = {false,false};  // inside/outside of facet
  CubitPoint* start_pts[2], *end_pts[2];  // edge points
  for ( int j = 0; j < 2; j++ )
  {
    for ( i = 0; i < bdy_edges[j]->num_adj_facets(); i++ )
    {
      CubitFacet* facet = bdy_edges[j]->adj_facet(i);
      if ( TDVGFacetOwner::get(facet) == this ) 
      {
        bdy_facets[j] = facet;
        break;
      }
    }
    assert(!!bdy_facets[j]);
    
    int index = bdy_facets[j]->edge_index(bdy_edges[j]);
    if ( bdy_facets[j]->edge_use(index) < 0 ) {
      start_pts[j] = bdy_edges[j]->point(1);
      end_pts[j] = bdy_edges[j]->point(0);
    } else {
      start_pts[j] = bdy_edges[j]->point(0);
      end_pts[j] = bdy_edges[j]->point(1);
    }   
    CubitVector dir = end_pts[j]->coordinates() - start_pts[j]->coordinates(); 
    double dot = facet->normal() % ((closest_pos - point) * dir);
    inside_facets[j] = dot > 0.0;
  }
  
    // mean normal for facets
  CubitVector normal = bdy_facets[0]->normal();
  normal += bdy_facets[1]->normal();
  
    // cross product of edge vectors
    // if we had edges in wrong order, need to reverse
    // crossproduct
  CubitVector cross;
  CubitVector dir1 = end_pts[0]->coordinates() - start_pts[0]->coordinates();
  CubitVector dir2 = end_pts[1]->coordinates() - start_pts[1]->coordinates();
  if ( end_pts[0] == closest_point && start_pts[1] == closest_point )
  {
    cross = dir1 * dir2;
  }
  else
  {
    assert( start_pts[0] == closest_point && end_pts[1] == closest_point );
    cross = dir2 * dir1;
  }
  
  bool convex = (cross % normal) > 0.0;
  if ( convex )
    return inside_facets[0] && inside_facets[1] ? CUBIT_PNT_INSIDE : CUBIT_PNT_OUTSIDE;
  else
    return inside_facets[0] || inside_facets[1] ? CUBIT_PNT_INSIDE : CUBIT_PNT_OUTSIDE;
}




CubitPointContainment PartitionSurface::point_containment( double , double )
{
  return CUBIT_PNT_UNKNOWN;
}
/*
CubitPointContainment PartitionSurface::point_containment( 
            CubitVector& point, double u, double v )
{
  return point_containment( point );
}
*/
CubitSense PartitionSurface::get_geometry_sense()
{ 
  return geometry_sense;
}

void PartitionSurface::reverse_sense()
{
  reverse_loops();
  if( owner() )
    owner()->notify_reversed(this);

  if(geometry_sense == CUBIT_FORWARD)
    geometry_sense = CUBIT_REVERSED;
  else
    geometry_sense = CUBIT_FORWARD;

  int j;
  DLIList<CubitFacetData*> surf_facets;
  this->get_facet_data( surf_facets );
  for(j=surf_facets.size(); j--;)
    surf_facets.get_and_step()->flip();
}

PartitionSurface::PartitionSurface( PartitionSurface* split_from )
  : firstLoop(0), firstCoSurf(0), geometry_sense(CUBIT_FORWARD)
{
  split_from->sub_entity_set().add_lower_order( this );
}
  
PartitionSurface* PartitionSurface::copy()
{
  return new PartitionSurface(this);
}

PartitionSurface* PartitionSurface::split( DLIList<CubitFacetData*>& facets )
{
  int i;
  bool okay = true;
  
    // mark facets
  for( i = facets.size(); i--; )
    facets.get_and_step()->marked(0);
  for( i = facetList.size(); i--; )
    facetList.get_and_step()->marked(1);
    
    // make sure all facets in passed list are in my list,
    // and clear the marks
  for( i = facets.size(); i--; ) {
    CubitFacetData* facet = facets.get_and_step();
    if( facet->marked() )
      facet->marked(0);
    else
      okay = false;
  }
  
  if (!okay) {
    for( i = facetList.size(); i--; )
      facetList.get_and_step()->marked(0);
    assert(okay);
    return 0;
  }
  
  for( i = facetList.size(); i--; ) {
    if( facetList.step_and_get()->marked() )
      facetList.get()->marked(0);
    else
      facetList.change_to(0);
  }
  facetList.remove_all_with_value(0);
  
  PartitionSurface* result = copy();
  result->set_facet_data(facets);
  return result;
}

CubitStatus PartitionSurface::combine( PartitionSurface* other_surf )
{
  for ( int i = other_surf->facetList.size(); i--;  )
    TDVGFacetOwner::set(other_surf->facetList.step_and_get(),this);
  facetList += other_surf->facetList;
  other_surf->facetList.clean_out();
  return CUBIT_SUCCESS;
}

void PartitionSurface::get_facet_data( DLIList<CubitFacetData*>& result_list ) const
{
  result_list = facetList;
}

void PartitionSurface::set_facet_data( const DLIList<CubitFacetData*>&  new_list )
{
  int i;
  
  for ( i = facetList.size(); i--; )
    TDVGFacetOwner::remove(facetList.step_and_get());
  facetList = new_list;
  for ( i = facetList.size(); i--; )
    TDVGFacetOwner::set(facetList.step_and_get(),this);
}  
  
void PartitionSurface::notify_split( FacetEntity* old_tri, FacetEntity* new_tri )
{
  if ( dynamic_cast<CubitFacetEdge*>(old_tri) )
    return;
  
  CubitFacetData* old_ptr = dynamic_cast<CubitFacetData*>(old_tri);
  CubitFacetData* new_ptr = dynamic_cast<CubitFacetData*>(new_tri);
  assert(old_ptr && new_ptr && facetList.is_in_list(old_ptr));
  TDVGFacetOwner::set(new_ptr,this);
  facetList.append(new_ptr);
}

CubitStatus PartitionSurface::init_facet_data()
{
  Surface* real_surf = dynamic_cast<Surface*>(partitioned_entity());
  //assert(are_marks_cleared());

  CubitStatus res = CUBIT_SUCCESS;
  GMem gMem;
  
  // TODO - tolerance arguments are defaulted.  Do we need to specify?
  res = real_surf->get_geometry_query_engine()->
    get_graphics(real_surf, &gMem);
  if( !res )
    return CUBIT_FAILURE;

  // Work around ACIS bug w/ coincident facet points 
  // for blend surfaces.  Consolidate coincident points.
  int old_count = gMem.pointListCount;
  gMem.consolidate_points(GEOMETRY_RESABS);
  if( old_count != gMem.pointListCount ) {
    PRINT_WARNING("Possible invalid facetting for surface.  "
                  "Coincident points found.\n");
  }

  DLIList<CubitFacetData*> surf_facets;
  CubitPointData** point_array = new CubitPointData*[gMem.pointListCount];

  // create CubitFacetData from GMem
  int i;
   
  CubitVector ipos, opos;
  for ( i = 0; i < gMem.pointListCount; i++ ) {
    ipos.set( gMem.point_list()[i].x, gMem.point_list()[i].y, gMem.point_list()[i].z );
    if( real_surf->closest_point( ipos, &opos ) )
      ipos = opos;
    point_array[i] = new CubitPointData( ipos );
  }

  int junk = 0;
  bool fail_out = false;
  CubitPoint *p1, *p2, *p3;
  for ( i = 0; i < gMem.fListCount; i += 4 ) 
  {
    if(gMem.facet_list()[i] != 3)
      fail_out = true;
    else
    {
      p1 = point_array[ gMem.facet_list()[i+1] ];
      p2 = point_array[ gMem.facet_list()[i+2] ];
      p3 = point_array[ gMem.facet_list()[i+3] ];

      if(p1 == p2 || p1 == p3 || p2 == p3)
        fail_out = true;
    }
    if(fail_out == true)
    {
      PRINT_ERROR("Non-triangular facet encountered.  Aborting.\n");
      while (surf_facets.size())
        delete surf_facets.pop();
      for (i = 0; i < gMem.pointListCount; i++)
        delete point_array[i];
      delete [] point_array;
      return CUBIT_FAILURE;
    }
    else
    {
      CubitFacetData* facet = new CubitFacetData(p1, p2, p3, &junk );
      surf_facets.append(facet);
    }
  }
  
  delete [] point_array; 

  if(geometry_sense == CUBIT_REVERSED)
  {
    for(i=surf_facets.size(); i--;)
      surf_facets.get_and_step()->flip();
  }
  
/*  
    // Make sure facet orientation is consistant
  DLIList<CubitFacet*> facets(surf_facets.size());
  for ( i = surf_facets.size(); i--; )
  {
    CubitFacet* facet = surf_facets.step_and_get();
    facet->marked(1);
    facets.append(facet);
  }
  bool reversed_some = false;
  while ( facets.size() )
  {
    CubitFacet* facet = facets.pop();
    facet->marked(0);
    bool reversed[3] = {false,false,false};
    CubitFacet* neighbors[3] = {0,0,0};
    for ( i = 0; i < 3; i++ )
    {
      CubitFacetEdge* edge = facet->edge(i);
      CubitFacet* other = edge->other_facet(facet);
      if( other && other != facet )
      {
        neighbors[i] = other;
        int index = other->edge_index(edge);
        if ( facet->edge_use(i) == other->edge_use(index) )
          reversed[i] = true;
      }
    }
    
    if ( !(reversed[0] || reversed[1] || reversed[2]) )
      continue;

    facet->flip();
    reversed_some = true;
    reversed[0] = !reversed[0];
    reversed[1] = !reversed[1];
    reversed[2] = !reversed[2];
    
    for ( i = 0; i < 3; i++ )
    {
      if( neighbors[i] && reversed[i] && !neighbors[i]->marked() )
      {
        neighbors[i]->marked(1);
        facets.append( neighbors[i] );
      }
    }
  }
  
  if( reversed_some ) 
  {
    PRINT_WARNING("Possible invalid facetting for surface.  "
                  "Inconsistent facet orientation.\n");
  }
  
    // pick a facet and compare to surface normal to see
    // if we have them all backwards
  surf_facets.reset();
  CubitFacet* facet = surf_facets.get();  
  CubitVector pts[3] = { facet->point(0)->coordinates(),
                         facet->point(1)->coordinates(),
                         facet->point(2)->coordinates()};
  CubitVector facet_norm = (pts[1] - pts[0]) * (pts[2] - pts[0]);
  CubitVector facet_cent = (pts[0] + pts[1] + pts[2]) / 3.0;
  CubitVector surf_norm;
  real_surf->closest_point( facet_cent, 0, &surf_norm );
  if ( facet_norm % surf_norm < 0 )
  {
    PRINT_WARNING("Possible invalid facetting for surface.  "
                  "Backwards facets.\n");
    for ( i = surf_facets.size(); i--; )
      surf_facets.step_and_get()->flip();
  }    
*/
  
  PartSurfFacetTool tool(this);
  return tool.init_facet_data(surf_facets);
}


void PartitionSurface::replace_facets( DLIList<CubitFacetData*> &dead_facets,
                                       DLIList<CubitFacetData*> &new_facets)
{
  int i;
  facetList -= dead_facets;
  facetList += new_facets;
  
  for( i = dead_facets.size(); i--; )
    TDVGFacetOwner::remove( dead_facets.step_and_get() );

  for( i = new_facets.size(); i--; )
    TDVGFacetOwner::set( new_facets.step_and_get(), this );

  // TODO - memory management - where should the facets get deleted
}

//-------------------------------------------------------------------------
// Purpose       : Find closest facet and point on facet
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 03/28/03
//-------------------------------------------------------------------------
CubitFacet* PartitionSurface::closest_facet( const CubitVector& input_position,
                                             CubitVector& result_position)
{
  return PartSurfFacetTool::closest_facet( input_position, facetList, result_position );
}


//-------------------------------------------------------------------------
// Purpose       : Save geometry in an attribute
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 01/21/03
//-------------------------------------------------------------------------
CubitStatus PartitionSurface::save( CubitSimpleAttrib& attrib )
{
  assert(dynamic_cast<Lump*>(partitioned_entity()) != 0);
  
  int i;
  int id = sub_entity_set().get_id(this);
  if( id <= 0 ) return CUBIT_FAILURE;
  
    // get facets
  DLIList<CubitFacetData*> facets;
  get_facet_data( facets );
  
    // get list of points from facets
  DLIList<CubitPoint*> points(facets.size()*3), facet_points(3);
  facets.reset();
  for( i = facets.size(); i--; )
  {
    facet_points.clean_out();
    facets.get_and_step()->points(facet_points);
    points += facet_points;
  }
  
  for( i = points.size(); i--; )
    points.step_and_get()->marked(0);
  for( i = points.size(); i--; )
  {
    CubitPoint* pt = points.step_and_get();
    pt->marked(pt->marked()+1);
  }
  points.last();
  for( i = points.size(); i--; )
  {
    if( points.get()->marked() > 1 )
    {
      points.get()->marked( points.get()->marked() - 1 );
      points.change_to(0);
    }
    points.back();
  }
  points.remove_all_with_value(0);
  
    // construct position list
  DLIList<CubitVector*> pt_list(points.size());
  points.reset();
  for( i = 0; i < points.size(); i++ )
  {
    CubitPoint* pt = points.get_and_step();
    pt_list.append( new CubitVector(pt->coordinates()) );
    pt->marked(i);
  }
  
    // connect facet connectivity list
  DLIList<int> facetlist( facets.size() * 3 );
  facets.reset();
  for( i = facets.size(); i--; )
  {
    facet_points.clean_out();
    facets.get_and_step()->points(facet_points);
    facet_points.reset();
    for( int j = facet_points.size(); j--; )
      facetlist.append( facet_points.get_and_step()->marked() );
  }
  
  DLIList<int> facet_point_owners;
  DLIList<CubitFacetEdge*> pt_edges;
  points.reset();
  for ( i = points.size(); i--; )
  {
    CubitPoint* pt = points.get_and_step();
    PartitionEntity* owner = TDVGFacetOwner::get(pt);
    if (!owner)
    {
      pt->edges(pt_edges);
      while ( pt_edges.size() )
      {
        CubitFacetEdge* edge = pt_edges.pop();
        PartitionEntity* tmp = TDVGFacetOwner::get(edge);
        if ( tmp )
        {
          assert(!owner || owner == tmp);
          owner = tmp;
        }
      }
    }
    if ( !owner )
      owner = this;

    int set_id, ent_id;
    if( &(owner->sub_entity_set()) == &(sub_entity_set()) )
    {
      set_id = 0;
      ent_id = sub_entity_set().get_id( owner );
    }
    else
    {
      set_id = owner->sub_entity_set().get_unique_id();
      ent_id = owner->sub_entity_set().get_id(owner);
    }
    facet_point_owners.append( set_id );
    facet_point_owners.append( ent_id );
  }
  
    // clean up point marks
  for( i = points.size(); i--; )
    points.step_and_get()->marked(0);
  
  DLIList<int> topo_list;
  get_save_topology( topo_list );
  
  return sub_entity_set().save_geometry( id, 2, &pt_list, &facetlist, &topo_list, 
                                         &facet_point_owners, attrib );
}

//-------------------------------------------------------------------------
// Purpose       : Get topology connectivity list for use in saving
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 01/21/03
//-------------------------------------------------------------------------
CubitStatus PartitionSurface::get_save_topology( DLIList<int>& topo_list )
{
  PartitionLoop* loop = 0;
  while( (loop = next_loop(loop)) )
  {
    topo_list.append(loop->num_coedges());
    PartitionCoEdge* coedge = loop->first_coedge();
    do
    {
      int set_id, crv_id;
      PartitionCurve* curve = coedge->get_curve();
      if( &(curve->sub_entity_set()) == &(sub_entity_set()) )
      {
        set_id = 0;
        crv_id = sub_entity_set().get_id( curve );
      }
      else
      {
        set_id = curve->sub_entity_set().get_unique_id();
        crv_id = curve->sub_entity_set().get_id(curve);
      }

      if( coedge->sense() == CUBIT_REVERSED )
        crv_id *= -1;
      topo_list.append(set_id);
      topo_list.append(crv_id);
    } while( (coedge = loop->next_coedge(coedge)) != loop->first_coedge() );
  }

  return CUBIT_SUCCESS;
}
    

void PartitionSurface::interior_facet_points( DLIList<CubitPoint*>& list ) const
{
  
  DLIList<CubitPoint*> point_list(facetList.size() * 3);
  DLIList<CubitFacetEdge*> edge_list;
  int i;
  for( i = 0; i < facetList.size(); i++ )
  {
    CubitFacetData* facet = facetList.next(i);
    for( int j = 0; j < 3; j++ )
    {
      CubitPoint* pt = facet->point(j);
      if( TDVGFacetOwner::get(pt) )
        continue;  // not interior - on a vertex
      
      edge_list.clean_out();
      pt->edges( edge_list );
      bool skip = false;
      for( int k = edge_list.size(); k--; )
        if( TDVGFacetOwner::get(edge_list.step_and_get()) )
          skip = true; // not interior - on a curve
        
      if( !skip )
        point_list.append(pt);
    }
  }
  
    // uniquify list
  for( i = point_list.size(); i--; )
    point_list.step_and_get()->marked(1);
  for( i = point_list.size(); i--; )
    if( point_list.step_and_get()->marked() )
      point_list.get()->marked(0);
    else
      point_list.change_to(0);
  point_list.remove_all_with_value(0);
  list = point_list;
}  
  

void PartitionSurface::transform( const CubitTransformMatrix& xform )
{
  int i;
  
  DLIList<CubitPoint*> points;
  interior_facet_points(points);
  for(i = points.size(); i--; )
  {
    CubitPoint* pt = points.step_and_get();
    pt->set( xform * pt->coordinates() );
  }
  
  double det = xform.sub_matrix( 3, 3 ).determinant();
  if (det < 0.0)
  {
    for (i = 0; i < facetList.size(); i++)
      facetList.next(i)->flip();
  
    PartitionLoop* loop = 0;
    while ((loop = next_loop(loop)))
      loop->reverse();
  }
  
  for (i = 0; i < facetList.size(); i++ )
  {
    facetList.next(i)->update_plane();
  }
  
}
      
  

PartitionSurface* PartitionSurface::construct( CubitSimpleAttrib& attrib,
                                               PartitionLump* vol )
{
    // construct surface and read attrib data
  PartitionSurface* new_surf = new PartitionSurface;
  DLIList<CubitVector*> facet_points;
  DLIList<int> facets, curve_conn, point_owners;
  vol->sub_entity_set().
    add_lower_order( new_surf, attrib, 2, facet_points, facets, curve_conn, point_owners );
  
  if( facets.size() % 3 != 0 )
  {
    assert( !(facets.size() % 3 ) );
    delete new_surf;
    return 0;
  }
  
    // create loops and coedges
  DLIList<PartitionCurve*> curve_list;
  int i = curve_conn.size();
  curve_conn.reset();
  bool okay = true;
  while( okay && i > 0 )
  {
    int k = curve_conn.get_and_step();
    i--;
    
    if( k <= 0 || i < 2*k )
      { okay = false; break; }
    
    PartitionLoop* new_loop = new PartitionLoop();
    PartitionCoEdge* prev_coedge = 0;
    new_surf->add( new_loop );
    
    for( int j = 0; okay && j < k; j++ )
    {
      int set_id = curve_conn.get_and_step();
      int ent_id = curve_conn.get_and_step();
      i -= 2;
      
      CubitSense sense = CUBIT_FORWARD;
      if( ent_id < 0 )
      {
        sense = CUBIT_REVERSED;
        ent_id = -ent_id;
      }
      
      PartitionEntity* ent = PartitionEngine::instance()
        .entity_from_id( set_id, ent_id, vol->sub_entity_set() );
      PartitionCurve* curve = dynamic_cast<PartitionCurve*>(ent);
      
      if( !curve )
        { okay = false; break;  }
      
      curve_list.append( curve );
      PartitionCoEdge* coedge = new PartitionCoEdge( new_surf, sense );
      new_loop->insert_after( coedge, prev_coedge );
      prev_coedge = coedge;
      curve->add( coedge );
    }
  }
  
  if( !okay )
  {
    PartitionEngine::instance().destroy_surface( new_surf );
    return 0;
  }
  
    // construct facet points
  CubitPointData** points = new CubitPointData*[facet_points.size()];
  facet_points.reset();
  point_owners.reset();
  for( i = 0; i < facet_points.size(); i++ )
  {
    CubitVector* pt = facet_points.get_and_step();
    points[i] = new CubitPointData(*pt);
    delete pt;
    
    int set_id = point_owners.get_and_step();
    int ent_id = point_owners.get_and_step();
    PartitionEntity* owner = PartitionEngine::instance()
      .entity_from_id( set_id, ent_id, vol->sub_entity_set() );
    
    if ( !owner )
    {
      okay = false;
      break;
    }
    
    if ( owner != new_surf )
    {
      if ( !dynamic_cast<PartitionPoint*>(owner) &&
           !dynamic_cast<PartitionCurve*>(owner) )
      {
        okay = false;
        break;
      }
      TDVGFacetOwner::set(points[i], owner);
    }
  }
  
  if ( !okay )
  {
    for ( i = 0; i < facet_points.size(); i++ )
      delete points[i];
    delete [] points;
    PartitionEngine::instance().destroy_surface( new_surf );
    return 0;
  }
    
    // construct facets
  facets.reset();
  DLIList<CubitFacetData*> facet_list;
  int junk = 0;
  for( i = 0; i < facets.size(); i += 3 )
  {
    CubitFacetData* facet = new CubitFacetData(
                                    points[facets.next(0)],
                                    points[facets.next(1)],
                                    points[facets.next(2)], &junk );
    facet_list.append(facet);
    facets.step(3);
    
      // set edge owners to this surface as a kind of mark
      // (will use this a little later)
    for ( int j = 0; j < 3; j++ )
      TDVGFacetOwner::set( facet->edge(j), new_surf );
  } 
  
    // Seam facet points on vertices
  for ( i = 0; i < facet_points.size(); i++ )
  {
    PartitionEntity* ent = TDVGFacetOwner::get(points[i]);
    PartitionPoint* vertex = dynamic_cast<PartitionPoint*>(ent);
    if (!vertex) continue;
    
    if ( !vertex->facet_point() )
      vertex->facet_point(points[i]);
    else if ( vertex->facet_point() != points[i] )
    {
      vertex->facet_point()->merge_points(points[i]);
      points[i] = vertex->facet_point();
    }
  }
  
    // done with point list -- delete it
  delete [] points;

    // seam facet edges with curves
  DLIList<CubitFacetEdge*> pt_edges;
  DLIList<CubitFacetEdgeData*> curve_edges;
  while ( curve_list.size() )
  {
    curve_edges.clean_out();
    PartitionCurve* curve = curve_list.pop();
    CubitPoint* pt = curve->start_point()->facet_point();
    if (!pt) {okay = false; break;}

    do {
      CubitFacetEdge *edge = 0, *pt_edge = 0;
      pt_edges.clean_out();
      pt->edges(pt_edges);
      while ( pt_edges.size() )
      {
        CubitFacetEdge* tmp_edge = pt_edges.pop();
        if ( TDVGFacetOwner::get(tmp_edge) != new_surf )
          continue;
          
        PartitionEntity* owner = TDVGFacetOwner::get(tmp_edge->other_point(pt));
        if ( owner == curve )
        {
          edge = tmp_edge;
          TDVGFacetOwner::set(tmp_edge->other_point(pt), 0);
          break;
        }
        else if ( owner == curve->end_point() )
        {
          pt_edge = tmp_edge;
        }
      }
      
      if ( !edge ) edge = pt_edge;
        
      CubitFacetEdgeData* edge_d = dynamic_cast<CubitFacetEdgeData*>(edge);
      if (!edge_d) {okay = false; break;}
      
      curve_edges.append(edge_d);
      pt = edge_d->other_point(pt);
    } while (pt != curve->end_point()->facet_point());
    
    if( !okay || !PartSurfFacetTool::seam_curve( curve_edges, curve, facet_list ) ) 
    {
      okay = false;
      break;
    }
  }
    
    // clean up any owner pointers on edges that were being
    // used as a kind of mark
  for ( i = facet_list.size(); i--; )
  {
    CubitFacet* facet = facet_list.get_and_step();
    for ( int j = 0; j < 3; j++ )
      if ( TDVGFacetOwner::get(facet->edge(j)) == new_surf )
        TDVGFacetOwner::set(facet->edge(j), 0);
  }
    
    // set surface facets
  new_surf->set_facet_data( facet_list );
  
  if (!okay)
  {
    PartitionEngine::instance().destroy_surface( new_surf );
    return 0;
  }

  return new_surf;
}


CubitStatus PartitionSurface::notify_moving_point( CubitPoint* point, 
                                           const CubitVector& new_pos )
{
  DLIList<CubitFacetData*> old_facets, new_facets;
  if (!PartSurfFacetTool::fix_move_point(point, new_pos, facetList, old_facets, new_facets, this))
    return CUBIT_FAILURE;
  
  assert( !old_facets.size() == !new_facets.size() );
  if (old_facets.size())
    replace_facets( old_facets, new_facets );
  return CUBIT_SUCCESS;
}

void PartitionSurface::notify_destroyed( CubitFacetData* facet )
{
  facetList.move_to(facet);
  assert(facetList.size() && facetList.get() == facet);
  facetList.extract();
}

