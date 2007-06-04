#include "PartitionShell.hpp"
#include "PartitionLump.hpp"
#include "PartitionSurface.hpp"
#include "VirtualQueryEngine.hpp"
#include "CubitFacetData.hpp"

PartitionShell::PartitionShell( )
  : myLump(0), lumpNext(0), firstCoSurf(0)
  { }

PartitionShell::~PartitionShell()
{
  assert( !myLump );
  remove_all_surfaces();
}

CubitStatus PartitionShell::add( PartitionCoSurf* cosurf )
{
  if( cosurf->myShell || cosurf->shellNext )
  {
    assert(!cosurf->myShell && !cosurf->shellNext);
    return CUBIT_FAILURE;
  }
  
  cosurf->myShell = this;
  cosurf->shellNext = firstCoSurf;
  firstCoSurf = cosurf;
  return CUBIT_SUCCESS;
}

CubitStatus PartitionShell::remove( PartitionCoSurf* cosurf )
{
  if( cosurf->myShell != this )
  {
    assert( cosurf->myShell == this );
    return CUBIT_FAILURE;
  }
  
  if( cosurf == firstCoSurf )
  {
    firstCoSurf->shellNext = cosurf->shellNext;
    firstCoSurf = firstCoSurf->shellNext;
  }
  else
  {
    PartitionCoSurf* prev = firstCoSurf;
    while( prev && prev->shellNext != cosurf )
      prev = prev->shellNext;
    
    if( !prev ) { assert(0); return CUBIT_FAILURE; }
  
    prev->shellNext = cosurf->shellNext;
  }
  
  cosurf->myShell = 0;
  cosurf->shellNext = 0;
  return CUBIT_SUCCESS;
}

PartitionCoSurf* PartitionShell::add( PartitionSurface* surf, CubitSense sense )
{
  PartitionCoSurf* new_cos = new PartitionCoSurf( sense );
  if( !surf->add(new_cos) || !this->add(new_cos) )
  {
    delete new_cos;
    new_cos = 0;
  }
  return new_cos;
}

PartitionCoSurf* PartitionShell::find_first( const PartitionSurface* surf ) const
{
  PartitionCoSurf* cos = firstCoSurf;
  while( cos && cos->get_surface() != surf )
    cos = cos->shellNext;
  return cos;
}
PartitionCoSurf* PartitionShell::find_next( const PartitionCoSurf* prev ) const
{
  PartitionCoSurf* cos = prev->shellNext;
  while( cos && cos->get_surface() != prev->get_surface() )
    cos = cos->shellNext;
  return cos;
}
CubitSense PartitionShell::find_sense( const PartitionSurface* surf ) const
{
  PartitionCoSurf* cos = find_first( surf );
  if( ! cos )
    return CUBIT_UNKNOWN;
  
  CubitSense result = cos->sense();
  while( (cos = find_next( cos )) )
    if( cos->sense() != result )
      return CUBIT_UNKNOWN;
  
  return result;
}




void PartitionShell::remove_all_surfaces( DLIList<PartitionSurface*>* list )
{
  while( PartitionCoSurf* cosurf = firstCoSurf )
  {
    if( cosurf->get_surface() )
    {
      if( list ) 
        list->append( cosurf->get_surface() );
      cosurf->get_surface()->remove( cosurf );
    }
    remove( cosurf );
    delete cosurf;
  }
  
  if( list )
    list->uniquify_unordered();
}

void PartitionShell::get_parents_virt( DLIList<TopologyBridge*>& parents )
{
  if( myLump )
    parents.append( myLump );
}

void PartitionShell::get_children_virt( DLIList<TopologyBridge*>& children )
{
  PartitionCoSurf* cos = 0;
  while( (cos = next_co_surface( cos )) )
    if( cos->get_surface() )
      children.append_unique( cos->get_surface() );
}

int PartitionShell::layer() const { return get_lump()->layer(); }

GeometryQueryEngine* PartitionShell::get_geometry_query_engine() const
{
  return VirtualQueryEngine::instance();
}

void PartitionShell::append_simple_attribute_virt( CubitSimpleAttrib* )
  { }
void PartitionShell::remove_simple_attribute_virt( CubitSimpleAttrib* )
  { }
void PartitionShell::remove_all_simple_attribute_virt() 
  { }
CubitStatus PartitionShell::get_simple_attribute( DLIList<CubitSimpleAttrib*>& )
  { return CUBIT_FAILURE; }
CubitStatus PartitionShell::get_simple_attribute( const CubitString&,
                                                 DLIList<CubitSimpleAttrib*>& )
{ return CUBIT_FAILURE; }

void PartitionShell::print_debug_info( const char* /*prefix*/ ) const
{
}

//-------------------------------------------------------------------------
// Purpose       : Determine point containment.
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 02/23/03
//-------------------------------------------------------------------------
CubitPointContainment PartitionShell::point_containment( const CubitVector& pt )
{
    // Find closest cosurface to passed point
  PartitionCoSurf* cosurf = 0;
  PartitionCoSurf* closest = 0;
  double dist_to_closest = CUBIT_DBL_MAX;
  CubitVector pt_on_surface;
  while( (cosurf = next_co_surface(cosurf)) )
  {
    PartitionSurface* surf = cosurf->get_surface();
    
      // Ignore surfaces that occur more than once in the 
      // shell.
    bool skip = false;
    PartitionCoSurf* surf_cosurf = 0;
    while( (surf_cosurf = surf->next_co_surface(surf_cosurf)) )
    {
      if( surf_cosurf->get_shell() == this && surf_cosurf != cosurf )
      {
        skip = true;
        break;
      }
    }
    
    if( skip )
      continue;
    
    surf->closest_point_trimmed( pt, pt_on_surface );
    double dist_sqr = (pt - pt_on_surface).length_squared();
    if( dist_sqr < dist_to_closest )
    {
      closest = cosurf;
      dist_to_closest = dist_sqr;
    }
  }
  
  if( !closest )
    return CUBIT_PNT_UNKNOWN;  // already printed error messages above
  
  if( dist_to_closest < (GEOMETRY_RESABS*GEOMETRY_RESABS) ) // point is on shell
    return CUBIT_PNT_BOUNDARY;
  
  CubitVector normal;
  closest->get_surface()->closest_point_trimmed( pt, pt_on_surface );
  closest->get_surface()->closest_point( pt_on_surface, 0, &normal );
  if( closest->sense() == CUBIT_REVERSED )
    normal *= -1.0;
  
  pt_on_surface -= pt;
  if ( pt_on_surface.length_squared() < GEOMETRY_RESABS*GEOMETRY_RESABS )
    return CUBIT_PNT_BOUNDARY;
  
  return (normal % pt_on_surface) > 0.0 ? CUBIT_PNT_INSIDE : CUBIT_PNT_OUTSIDE;
}

bool PartitionShell::is_nonmanifold( PartitionSurface* surf ) const
{
  int count = 0;
  PartitionCoSurf* cosurf = 0;
  while ((cosurf = surf->next_co_surface( cosurf )))
    if (cosurf->get_shell() == this)
      count++;
  
  return count != 1;
}

CubitStatus PartitionShell::mass_properties( CubitVector& centroid, 
                                             double& volume )
{
  PartitionCoSurf* cosurf = 0;
  DLIList<CubitFacetData*> facets;
  CubitVector p1, p2, p3, normal;
  const CubitVector p0(0.0, 0.0, 0.0);
  centroid.set(0.0, 0.0, 0.0 );
  volume = 0.0;
  
  while ((cosurf = next_co_surface( cosurf )))
  {
    if (is_nonmanifold( cosurf->get_surface() ))
      continue;
    
    facets.clean_out();
    cosurf->get_surface()->get_facet_data( facets );
    
    for (int i = facets.size(); i--; )
    {
      CubitFacet* facet = facets.step_and_get();
      p1 = facet->point(0)->coordinates();
      p2 = facet->point(1)->coordinates();
      p3 = facet->point(2)->coordinates();
      normal = (p3 - p1) * (p2 - p1);
  
      double two_area = normal.length();
      if (two_area > CUBIT_RESABS )
      {
        if (cosurf->sense() == CUBIT_REVERSED)
          normal = -normal;
        
        normal /= two_area;
        
        double height = normal % (p0 - p1);
        double vol = two_area * height;
        
        volume += vol;
        centroid += vol * (p0 + p1 + p2 + p3);
      }
    }
  }
  
  if (volume > CUBIT_RESABS)
    centroid /= 4.0 * volume;
  volume /= 6.0;
  return CUBIT_SUCCESS;
}
  
