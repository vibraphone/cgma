//-------------------------------------------------------------------------
// Filename      : CompositeSurface.cc
//
// Purpose       : Implementation of the CompositeSurface class.
//
// Special Notes : 
// 
// Creator       : Jason Kraftcheck
//
// Creation Date : 03/10/98
//
// Owner         : Jason Kraftcheck
//-------------------------------------------------------------------------

#include <assert.h>

#include "VirtualQueryEngine.hpp"

#include "Surface.hpp"

#include "CompositeSurface.hpp"
#include "CompositeLoop.hpp"
#include "CompositeCurve.hpp"
#include "CompositePoint.hpp"
#include "CompositeShell.hpp"
#include "CompositeLump.hpp"

#include "CompositeEngine.hpp"
#include "CompSurfFacets.hpp"
#include "GfxDebug.hpp"

/*
#include "CpuTimer.hpp"
static double trimmed_time = 0.0;
static double contain_time = 0.0;
static double closest_time = 0.0;
static int total_calls = 0;
static int contain_trim_count = 0;
static int error_count = 0;
static CpuTimer timer;
*/

CompositeSurface::CompositeSurface( Surface* surface )
  : stitchPartner(0), firstCoSurf(0), firstLoop(0), hiddenSet(0), facetTool(0),
  HadBridgeRemoved(0)
{
  assert( surface != NULL );
  compGeom = new CompositeGeom(1);
  compGeom->append( surface, CUBIT_FORWARD );
  if( surface->owner() )
    surface->owner()->swap_bridge( surface, this, false );
  surface->owner(this);
}

CompositeSurface::CompositeSurface( CompositeGeom* geometry )
  : compGeom( geometry ),
    stitchPartner(0),
    firstCoSurf(0),
    firstLoop(0),
    hiddenSet(0),
    facetTool(0),
    HadBridgeRemoved(0)
{
  assert( geometry != NULL );
  for( int i = 0; i < compGeom->num_entities(); i++ )
  {
    GeometryEntity* entity = compGeom->entity(i);
    assert( !entity->owner() );
    entity->owner(this);
  }
}



//-------------------------------------------------------------------------
// Purpose       : Desctructor
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 
//-------------------------------------------------------------------------
CompositeSurface::~CompositeSurface()
{
  while( firstCoSurf )
  {
    CompositeCoSurf* cosurf = firstCoSurf;
    remove( cosurf );
    if( cosurf->get_shell() )
      cosurf->get_shell()->remove( cosurf );
    delete cosurf;
  }
  
  while( firstLoop )
    remove( firstLoop );
  
  for( int j = 0; j < num_surfs(); j++ )
    if( get_surface(j)->owner() == this )
      get_surface(j)->owner(0);
  
  if( stitchPartner )
  {
    stitchPartner->stitchPartner = 0;
    stitchPartner = 0;
  }
  
  delete hiddenSet;
  delete compGeom;
  delete facetTool;
  hiddenSet = (HiddenEntitySet*)0xbdbdbdbd;
  compGeom = (CompositeGeom*)0xbdbdbdbd;
}

//-------------------------------------------------------------------------
// Purpose       : Add a CompositeLoop to child list
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 03/06/02
//-------------------------------------------------------------------------
CubitStatus CompositeSurface::add( CompositeLoop* loop )
{
  if( loop->mySurface )
  {
    assert(0);
    return CUBIT_FAILURE;
  }
  
  loop->mySurface = this;
  loop->loopNext = firstLoop;
  firstLoop = loop;
  
  return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Purpose       : Remove a child loop
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 03/06/02
//-------------------------------------------------------------------------
CubitStatus CompositeSurface::remove( CompositeLoop* loop )
{
  if( loop->mySurface != this )
    return CUBIT_FAILURE;
 
  if( firstLoop == loop )
  {
    firstLoop = loop->loopNext;
  }
  else
  {
    CompositeLoop *prev = firstLoop,
                  *next = firstLoop->loopNext;
    
    while( next != loop )
    {
      assert( next != NULL );
      prev = next;
      next = next->loopNext;
    }
    
    prev->loopNext = next->loopNext;
  }
  
  loop->loopNext = 0;
  loop->mySurface = 0;
  return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Purpose       : Add a CoSurface to this Surface
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 08/07/02
//-------------------------------------------------------------------------
CubitStatus CompositeSurface::add( CompositeCoSurf* cosurf )
{
  if( cosurf->mySurface )
    return CUBIT_FAILURE;
  
  cosurf->mySurface = this;
  cosurf->surfaceNext = firstCoSurf;
  firstCoSurf = cosurf;
  return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Purpose       : Remove a CoSurface
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 08/07/02
//-------------------------------------------------------------------------
CubitStatus CompositeSurface::remove( CompositeCoSurf* cosurf )
{
  if( cosurf->mySurface != this )
    return CUBIT_FAILURE;
  
  if( cosurf == firstCoSurf )
    firstCoSurf = cosurf->surfaceNext;
  else
  {
    CompositeCoSurf* prev = firstCoSurf;
    while( prev && prev->surfaceNext != cosurf )
      prev = prev->surfaceNext;
    assert( prev != NULL );
    prev->surfaceNext = cosurf->surfaceNext;
  }
  
  cosurf->mySurface = 0;
  cosurf->surfaceNext = 0;
  return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Purpose       : Get a CoSurface attaching this surface to the passed 
//                 Shell or Lump
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 08/07/02
//-------------------------------------------------------------------------
CompositeCoSurf* CompositeSurface::find_first( CompositeShell* shell ) const
{
  CompositeCoSurf* cos = firstCoSurf;
  while( cos && cos->get_shell() != shell )
    cos = cos->next_in_surface();
  return cos;
}
CompositeCoSurf* CompositeSurface::find_first( CompositeLump* lump ) const
{
  CompositeCoSurf* cos = firstCoSurf;
  while( cos && (!cos->get_shell() || cos->get_shell()->get_lump() != lump ) )
    cos = cos->next_in_surface();
  return cos;
}
CompositeCoSurf* CompositeSurface::find_next( CompositeCoSurf* cosurf ) const
{
  CompositeCoSurf* cos = cosurf;
  while( cos && cos->get_shell() != cosurf->get_shell() )
    cos = cos->next_in_surface();
  return cos;
}


//-------------------------------------------------------------------------
// Purpose       : Split this CompositeSurface into two.
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 03/06/02
//-------------------------------------------------------------------------
CompositeSurface* CompositeSurface::split( VGArray<int>& indices_to_move )
{
  int i;
  
  for( i = 0; i < indices_to_move.size(); i++ )
    if( indices_to_move[i] < 0 || indices_to_move[i] >= num_surfs() )
      return 0;
  
  CompositeGeom* new_geom = compGeom->split( indices_to_move );
  if( !new_geom )
    return 0;
  
  for( i = 0; i < new_geom->num_entities(); i++ )
    new_geom->entity(i)->owner( 0 );
    
  delete facetTool;
  facetTool = 0;
  
  return new CompositeSurface( new_geom );
}

  
//-------------------------------------------------------------------------
// Purpose       : Combine composite surfaces
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 03/06/02
//-------------------------------------------------------------------------
CubitStatus CompositeSurface::combine( CompositeSurface* dead_surf )
{
  int old_size = compGeom->num_entities();

  // Merge the "surfaces_to_ignore" list.
  surfacesToIgnore.merge_unique(dead_surf->surfacesToIgnore);

  compGeom->merge( *(dead_surf->compGeom) );
  if( dead_surf->hiddenSet != 0 )
    hidden_entities().merge( dead_surf->hiddenSet );
  for( int i = old_size; i < compGeom->num_entities(); i++ )
  {
    TopologyBridge* bridge = compGeom->entity(i);
    assert( bridge->owner() == dead_surf );
    bridge->owner( this );
  }
  
  delete facetTool;
  facetTool = 0;
  
  return CUBIT_SUCCESS;
}

void CompositeSurface::get_ignored_surfs(DLIList<Surface*> &surfs)
{
  int i;

  if(surfacesToIgnore.size() > 0)
  {
    for(i=0; i<num_surfs(); i++)
    {
      Surface *srf = get_surface(i);
      if(surfacesToIgnore.is_in_list(srf))
        surfs.append(srf);
    }
  }
}

//-------------------------------------------------------------------------
// Purpose       : Return the bounding box
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 05/28/99
//-------------------------------------------------------------------------
CubitBox CompositeSurface::bounding_box() const
{
  return compGeom->bounding_box();
}

//-------------------------------------------------------------------------
// Purpose       : TB Queries
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 03/06/02
//-------------------------------------------------------------------------
void CompositeSurface::get_parents_virt( DLIList<TopologyBridge*>& list )
{
  DLIList<TopologyBridge*> parents, parents2;
  for( int i = 0; i < num_surfs(); i++ )
  {
    parents.clean_out();
    get_surface(i)->get_parents( parents );
    parents.reset();
    for ( int j = parents.size(); j--; )
    {
      TopologyBridge* shell = parents.get_and_step();
      shell->get_parents( parents2 );
      assert (parents2.size() == 1);
      if (0 == dynamic_cast<CompositeLump*>(parents2.pop()->owner()))
        list.append_unique( shell );
    }
  }
  
  CompositeCoSurf* cosurf = 0;
  while ((cosurf = next_co_surface( cosurf )))
    list.append_unique( cosurf->get_shell() );
}
void CompositeSurface::get_children_virt( DLIList<TopologyBridge*>& list )
{
  for( CompositeLoop* loop = firstLoop; loop; loop = loop->loopNext )
    list.append( loop );
}


//-------------------------------------------------------------------------
// Purpose       : Attach a CubitSimpleAttribute
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 03/10/98
//-------------------------------------------------------------------------
void CompositeSurface::append_simple_attribute_virt(
						    CubitSimpleAttrib* simple_attrib_ptr )
{
  compGeom->add_attribute( simple_attrib_ptr );
}

//-------------------------------------------------------------------------
// Purpose       : Remove an attached CubitSimpleAttrib
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 03/10/98
//-------------------------------------------------------------------------
void CompositeSurface::remove_simple_attribute_virt(
						    CubitSimpleAttrib* simple_attrib_ptr )
{
  compGeom->rem_attribute( simple_attrib_ptr );
}


//-------------------------------------------------------------------------
// Purpose       : Remove an all attached CubitSimpleAttrib
//
// Special Notes : 
//
// Creator       : Greg Nielson
//
// Creation Date : 07/10/98
//-------------------------------------------------------------------------
void CompositeSurface::remove_all_simple_attribute_virt()
{
  compGeom->rem_all_attributes( );
}


//-------------------------------------------------------------------------
// Purpose       : Return the attached CubitSimpleAttribs.
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 03/10/98
//-------------------------------------------------------------------------
CubitStatus CompositeSurface::get_simple_attribute(
						   DLIList<CubitSimpleAttrib*>& attrib_list )
{
  compGeom->get_attributes( attrib_list );
  return CUBIT_SUCCESS;
}
CubitStatus CompositeSurface::get_simple_attribute(
					const CubitString& name, DLIList<CubitSimpleAttrib*>& attrib_list )
{
  compGeom->get_attributes( name.c_str(), attrib_list );
  return CUBIT_SUCCESS;
}


//-------------------------------------------------------------------------
// Purpose       : Methods from TBOwner
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 03/07/02
//-------------------------------------------------------------------------
CubitStatus CompositeSurface::remove_bridge( TopologyBridge* bridge )
{
  int i = compGeom->index_of(bridge);
  if( i < 0 )
    return CUBIT_FAILURE;
  
  delete facetTool;
  facetTool = 0;
  
  assert( bridge->owner() == this );
  bridge->owner(0);
  if (!compGeom->remove(i, true))
    return CUBIT_FAILURE;
  
  if (compGeom->num_entities() == 0)
    CompositeEngine::instance().notify_deactivated(this);
  HadBridgeRemoved = 1;
  return CUBIT_SUCCESS;
}

Surface* CompositeSurface::remove_surface( int index )
{
  Surface* result = get_surface(index);
  if ( !result || !compGeom->remove(index,false) ) 
    return 0;

  result->owner(0);
  return result;
}

  

CubitStatus CompositeSurface::swap_bridge( TopologyBridge* o,
                                           TopologyBridge* n, 
                                           bool reversed )
{
  if( n->owner() )
    return CUBIT_FAILURE;

  int i = compGeom->index_of(o);
  GeometryEntity* ge = dynamic_cast<GeometryEntity*>(n);
  if( i >= 0 && ge != 0 )
  {
    o->owner(0);
    n->owner(this);
    if ( !compGeom->swap( i, ge ) )
      return CUBIT_FAILURE;
    
    if (reversed)
      compGeom->reverse_sense(i);
    return CUBIT_SUCCESS;
  }
  else
    return CUBIT_FAILURE;
}
CubitBoolean CompositeSurface::contains_bridge( TopologyBridge* bridge ) const
{
  return (CubitBoolean)(compGeom->index_of(bridge) >= 0);
}


    


//-------------------------------------------------------------------------
// Purpose       : Return a pointer to VirtualQueryEngine
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 03/10/98
//-------------------------------------------------------------------------
GeometryQueryEngine* CompositeSurface::get_geometry_query_engine() const
{
  return VirtualQueryEngine::instance();
}


//-------------------------------------------------------------------------
// Purpose       : Return the surface area.
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 03/10/98
//-------------------------------------------------------------------------
double CompositeSurface::measure()
{
  return compGeom->measure();
}


//-------------------------------------------------------------------------
// Purpose       : Find the closest point to the passed point.
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 03/10/98
//-------------------------------------------------------------------------
void CompositeSurface::closest_point_trimmed( CubitVector from_point,
                                              CubitVector& point_on_surf )
{
  int index = closest_underlying_surface( from_point );
  get_surface(index)->closest_point_trimmed( from_point, point_on_surf );
}



CubitStatus CompositeSurface::get_point_normal( CubitVector& origin,
                                                CubitVector& normal )
{
  int count = num_surfs();
  CubitVector* vect_list = new CubitVector[count];
  double RESABS_SQUARED = CUBIT_RESABS * CUBIT_RESABS;
  normal.set(0.0,0.0,0.0);
	
  if (count == 1)
  {
    return get_surface(0)->get_point_normal(origin,normal);
  }
  
  for( int i = 0; i < count; i++ )
  {
    Surface* surf = get_surface(i);
    if( surf->get_point_normal(origin,vect_list[i]) == CUBIT_FAILURE )
    {
      delete [] vect_list;
      return CUBIT_FAILURE;
    }

    if( compGeom->sense(i) == CUBIT_REVERSED )
      vect_list[i] *= -1.0;
    normal += vect_list[i];
  }
  //If we reach this point, then all of the underlying surfaces are planar.
  //Next check if they are coplanar.
  if( normal.length_squared() < RESABS_SQUARED )
  {
    delete [] vect_list;
    return CUBIT_FAILURE;
  }
  normal.normalize();
  for( int j = 0; j < count; j++ )
  {
    if( fabs( 1.0 - (normal % ~vect_list[j]) ) > CUBIT_RESABS )
    {
      delete [] vect_list;
      return CUBIT_FAILURE;
    }
  }
	
  delete [] vect_list;
  CubitVector zero( 0.0, 0.0, 0.0 );
  closest_point( zero, &origin );
  return CUBIT_SUCCESS;
}	

// The CompositeSurface class has a variable 
// "surfacesToIgnore" which specifies which surfaces
// within the composite should be ignored during evaluation.
// The "facetTool" variable has a corresponding list of flags
// and this function syncs the flags in the "facetTool" variable
// with those that are in the "surfacesToIgnore" variable.
void CompositeSurface::update_facets_to_ignore()
{
   if(facetTool)
   {
      int i;
      DLIList<int> surfaces_to_ignore;
      int num_surfs_in_composite = num_surfs();
      for (i=0; i<num_surfs_in_composite; i++)
      {
         Surface *cur_surf = get_surface(i);
         surfacesToIgnore.reset();
         for(int j=surfacesToIgnore.size(); j--;)
         {
            if(cur_surf == surfacesToIgnore.get_and_step())
            {
               surfaces_to_ignore.append( i );
               j=0;
            }
         }
      }
      
      //do it all at once
      facetTool->set_ignore_flag( surfaces_to_ignore, 1 );
   }
}


// This function tells the composite to ignore one of its
// surfaces during evaluation.
void CompositeSurface::ignore_surface(int surface_id)
{
   update_facet_tool();
   if(facetTool)
   {
      int i;
      int num_surfs_in_composite = num_surfs();
      for (i=0; i<num_surfs_in_composite; i++)
      {
         Surface *cur_surf = get_surface(i);
         if(cur_surf->get_saved_id() == surface_id)
         {
            surfacesToIgnore.append_unique(cur_surf);
            i = num_surfs_in_composite;
            update_facets_to_ignore();
         }
      }
   }
}

void CompositeSurface::ignore_surface(Surface *surf)
{
   update_facet_tool();
   if(facetTool)
   {
      int i;
      int num_surfs_in_composite = num_surfs();
      for (i=0; i<num_surfs_in_composite; i++)
      {
         Surface *cur_surf = get_surface(i);
         if(cur_surf == surf)
         {
            surfacesToIgnore.append_unique(cur_surf);
            i = num_surfs_in_composite;
            update_facets_to_ignore();
         }
      }
   }
}

// This function tells the composite to unset the ignore
// flag for one of its surfaces.
void CompositeSurface::unignore_surface(int surface_id)
{
   update_facet_tool();
   if(facetTool)
   {
      int i;
      int num_surfs_in_composite = num_surfs();
      for (i=0; i<num_surfs_in_composite; i++)
      {
         Surface *cur_surf = get_surface(i);
         if(cur_surf->get_saved_id() == surface_id)
         {
            surfacesToIgnore.remove(cur_surf);
            update_facets_to_ignore();
            i = num_surfs_in_composite;
         }
      }
   }
}

void CompositeSurface::update_facet_tool()
{
  if( ! facetTool )
  {
    std::vector<Surface*> surf_vect(num_surfs());
    for ( int i = 0; i < num_surfs(); i++ )
      surf_vect[i] =  get_surface(i);
    facetTool = new CompSurfFacets();
    if ( ! facetTool->setup( surf_vect ) )
    {
      delete facetTool;
      facetTool = 0;
    }
    else
    {
      // Make sure to update the facetTool to reflect
      // any surfaces we think we need to ignore.
      update_facets_to_ignore();
    }
  }
}

CubitStatus CompositeSurface::closest_point_uv_guess(  
    CubitVector const& location,
    double &u, double &v,
    CubitVector* closest_location,
    CubitVector* unit_normal )
{
  if ( num_surfs() == 1)
    return get_surface(0)->
      closest_point_uv_guess(location, u, v, closest_location, unit_normal);
  else
    return closest_point(location, closest_location, unit_normal);
}

CubitStatus CompositeSurface::evaluate( double u, double v,
                                        CubitVector *position,                                   
                                        CubitVector *normal,
                                        CubitVector *curvature1,
                                        CubitVector *curvature2 )
{
  if( position || normal || (curvature1 && curvature2) )
  {
    if ( num_surfs() == 1)
      return get_surface(0)->evaluate(u, v, position, normal, curvature1, curvature2 ); 
    else
      return CUBIT_FAILURE; 
  }
  else
    return CUBIT_FAILURE;
}



/*
//-------------------------------------------------------------------------
// Purpose       : Find the closest point to the passed point, on the surface.
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 03/10/98
//-------------------------------------------------------------------------
CubitStatus CompositeSurface::closest_point(
  CubitVector const& location,
  CubitVector* closest_location,
  CubitVector* unit_normal,
  CubitVector* curvature1,
  CubitVector* curvature2 )
{
  if ( num_surfs() == 1 )
    return get_surface(0)->closest_point( location, 
      closest_location, unit_normal, curvature1, curvature2 );
  
  update_facet_tool();
  if ( facetTool )
  {
    CubitVector facet_closest, surf_closest;
    int index = facetTool->closest_index( location, &facet_closest );
    CubitStatus result = get_surface(index)->closest_point( location,
      &surf_closest, unit_normal, curvature1, curvature2 );
    
    if (!result)
      return result;
      
    CubitVector facet_delta = facet_closest - location;
    CubitVector  surf_delta =  surf_closest - location;
    double facet_len = facet_delta.length();
    double  surf_len =  surf_delta.length();
    
    if ( facet_len > CUBIT_RESABS && surf_len > CUBIT_RESABS )
    {
      facet_delta /= facet_len;
      surf_delta /= surf_len;
      const double cos_angle = facet_delta % surf_delta;
      if (0.985 < cos_angle)  // angle greater than about 10 degrees
      {
        result = get_surface(index)->closest_point( facet_closest,
                     &surf_closest, unit_normal, curvature1, curvature2 );
      
        if (!result)
          return result;
      }
    }
    
    if (closest_location)
      *closest_location = surf_closest;
    
    if (unit_normal && sense(index) == CUBIT_REVERSED)
      *unit_normal = -*unit_normal;
      
    return result;
  }
  
  
  double shortest_dist_sqr, current_dist_sqr;
  CubitVector closest_point, current_point;
  int closest_surf, current_surf;
	
  //initialize CompositeEntity data structures
  closest_surf = compGeom->closest_box( location );
	
  closest_trimmed( closest_surf, location, closest_point );
  shortest_dist_sqr = (location - closest_point).length_squared();
	
  while( (current_surf = compGeom->next_box_within_dist( shortest_dist_sqr ) ) >= 0 )
  {
    closest_trimmed( current_surf, location, current_point );
    current_dist_sqr = (location - current_point).length_squared();
	
    if( current_dist_sqr < shortest_dist_sqr )
    {
      closest_surf = current_surf;
      closest_point = current_point;
      shortest_dist_sqr = current_dist_sqr;
    }
  }

  if( closest_location ) *closest_location = closest_point;
  if( unit_normal || curvature1 || curvature2 )
  {
    get_surface( closest_surf )->closest_point( closest_point, NULL,
                              unit_normal, curvature1, curvature2 );
    if( unit_normal && compGeom->sense(closest_surf) == CUBIT_REVERSED )
      *unit_normal *= -1;
  }

  return CUBIT_SUCCESS;
}
*/
//-------------------------------------------------------------------------
// Purpose       : Find the closest point to the passed point, on the surface.
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 03/10/98
//-------------------------------------------------------------------------
CubitStatus CompositeSurface::closest_point( CubitVector const& location,
                                             CubitVector* closest_location,
                                             CubitVector* unit_normal,
                                             CubitVector* curvature1,
                                             CubitVector* curvature2 )
{
  if ( num_surfs() == 1 )
    return get_surface(0)->closest_point( location, closest_location,
                                          unit_normal,
                                          curvature1, curvature2 );
  
  update_facet_tool();
  
  if ( facetTool )
  {
    CubitStatus result;
    CubitVector facet_closest;

      // look for multiple surfaces if normal is requested
    if (unit_normal)
    {
      DLIList<int> index_list;
      int num_found = facetTool->closest_index( location, index_list,
                                                &facet_closest );

      CubitVector normal(0.0, 0.0, 0.0);
      int i;
      for (i = 0; i < num_found; i++) 
      {
        int index = index_list[i];
        if(index > -1)
        {
          Surface* surf = get_surface(index);
          
          result = surf->closest_point( facet_closest, closest_location,
                                        &normal, curvature1, curvature2 );
          
          if (get_sense(index) == CUBIT_REVERSED)
            *unit_normal += (-normal);
          else
            *unit_normal += normal;
        }
      }
      unit_normal->normalize();
    }
    else
    {  
      int index = facetTool->closest_index( location, &facet_closest );
      if(index > -1)
      {
        Surface* surf = get_surface(index);
      
        result = surf->closest_point( facet_closest, closest_location,
                                      unit_normal, curvature1, curvature2 );
    
  //       if (unit_normal && get_sense(index) == CUBIT_REVERSED)
  //         *unit_normal = -*unit_normal;
      }
      else
        result = CUBIT_FAILURE;
    }
    
    return result;

      // this code is never accessed
    Surface* surf = NULL;
    int index = -1;
    double u, v;
    result = surf->u_v_from_position( facet_closest, u, v );
    if (!result) return CUBIT_FAILURE;
    
    CubitVector surf_closest;
    result = surf->closest_point_uv_guess( location, u, v, &surf_closest, 
                                           unit_normal );
    if (!result) return CUBIT_FAILURE;
    
    if (unit_normal && get_sense(index) == CUBIT_REVERSED)
      *unit_normal = -*unit_normal;
    
    if (curvature1 || curvature2)
    {
      result = surf->closest_point( surf_closest, 0, 0, 
                                    curvature1, curvature2 );
      if (!result) return CUBIT_FAILURE;
    }
    
    if(closest_location)
      *closest_location = surf_closest;
  
    return CUBIT_SUCCESS;
  }
  
  double shortest_dist_sqr, current_dist_sqr;
  CubitVector closest_point, current_point;
  int closest_surf, current_surf;
	
  //initialize CompositeEntity data structures
  closest_surf = compGeom->closest_box( location );
	
  closest_trimmed( closest_surf, location, closest_point );
  shortest_dist_sqr = (location - closest_point).length_squared();
	
  while( (current_surf = compGeom->next_box_within_dist( shortest_dist_sqr ) ) >= 0 )
  {
    closest_trimmed( current_surf, location, current_point );
    current_dist_sqr = (location - current_point).length_squared();
	
    if( current_dist_sqr < shortest_dist_sqr )
    {
      closest_surf = current_surf;
      closest_point = current_point;
      shortest_dist_sqr = current_dist_sqr;
    }
  }

  if( closest_location ) *closest_location = closest_point;
  if( unit_normal || curvature1 || curvature2 )
  {
    get_surface( closest_surf )->closest_point( closest_point, NULL,
                              unit_normal, curvature1, curvature2 );
    if( unit_normal && compGeom->sense(closest_surf) == CUBIT_REVERSED )
      *unit_normal *= -1;
  }

  return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Purpose       : Find closest underlying surface
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 02/18/03
//-------------------------------------------------------------------------
int CompositeSurface::closest_underlying_surface( const CubitVector& pos )
{
  if( num_surfs() == 1 )
    return 0;
  
  update_facet_tool();
  if( facetTool )
    return facetTool->closest_index( pos );
 
  double shortest_dist_sqr, current_dist_sqr;
  CubitVector closest_point;
  int closest_surf, current_surf;
	
  //initialize CompositeEntity data structures
  closest_surf = compGeom->closest_box( pos );
	
  closest_trimmed( closest_surf, pos, closest_point );
  shortest_dist_sqr = (pos - closest_point).length_squared();
	
  while( (current_surf = compGeom->next_box_within_dist( shortest_dist_sqr ) ) >= 0 )
  {
    closest_trimmed( current_surf, pos, closest_point );
    current_dist_sqr = (pos - closest_point).length_squared();
	
    if( current_dist_sqr < shortest_dist_sqr )
    {
      closest_surf = current_surf;
      shortest_dist_sqr = current_dist_sqr;
    }
  }

  return closest_surf;
}

//-------------------------------------------------------------------------
// Purpose       : Evaluate closest_point_trimmed un an underliying surface
//
// Special Notes : honors use_gme_cpt, and updates stats
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 11/10/99
//-------------------------------------------------------------------------
void CompositeSurface::print_cpt_stats()
{
/*
  PRINT_INFO("Total Calls    %10d\n",          total_calls);
  PRINT_INFO("ClosestTrimmed %10.0f %10.5f\n", trimmed_time, trimmed_time / total_calls );
  PRINT_INFO("Closest        %10.0f %10.5f\n", closest_time, closest_time / total_calls );
  PRINT_INFO("Containment    %10.0f %10.5f\n", contain_time, contain_time / total_calls );
  double average_outside = (double)contain_trim_count / total_calls;
  PRINT_INFO("Outside Count  %10d   %10.5f\n", contain_trim_count, average_outside );
  double containment = contain_time + 
                       average_outside * trimmed_time + 
                       (1.0-average_outside) * closest_time;
  PRINT_INFO("Contain Est.   %10.0f %10.5f\n", containment, containment / total_calls );
  PRINT_INFO("Error Count    %10d   %10.5f\n", error_count, (double)error_count / total_calls );
*/
}
void CompositeSurface::reset_cpt_stats()
{
/*
  trimmed_time = contain_time = closest_time = 0.0;
  total_calls = contain_trim_count = error_count = 0;
*/
}
CubitStatus CompositeSurface::closest_trimmed( int index,
	                           const CubitVector& position, 
                             CubitVector& result )
{
  get_surface(index)->closest_point_trimmed( position, result );
  return CUBIT_SUCCESS;
/*
  total_calls++;
  CubitVector close, copy(position);
  Surface* surf = get_surface(index);
  
  timer.cpu_secs();
  surf->closest_point_trimmed( position, result );
  trimmed_time += timer.cpu_secs();
  
  surf->closest_point( position, &close );
  closest_time += timer.cpu_secs();
  
  CubitPointContainment contain = surf->point_containment( copy );
  contain_time += timer.cpu_secs();
  
  if ( contain == CUBIT_PNT_OUTSIDE )
    contain_trim_count++;
  else if( (result - close).length_squared() > GEOMETRY_RESABS*GEOMETRY_RESABS)
    error_count++;
  
  return CUBIT_SUCCESS;
*/
//  if( use_gme_cpt )
//	{
//		get_surface( index )->closest_point_trimmed( position, result );
//		return CUBIT_SUCCESS;
//  }
/*  
  DLIList<TopologyBridge*> bridge_list;
	Surface* surf_ptr = get_surface(index);
  CubitVector surf_pt, normal, curve_pt;
  if( !surf_ptr->closest_point( position, &surf_pt, &normal ) )
    return CUBIT_FAILURE;

  CoEdgeSM *closest_coedge, *other_coedge = 0;
  cptInfo.setup(surf_ptr);
  cptInfo.closest_coedge( position, closest_coedge, other_coedge, curve_pt );
  if ( !closest_coedge )
    return CUBIT_FAILURE;

  CubitVector coe_normal, cross, tangent1, tangent2, junk;
  bool inside;


  if ( !other_coedge )
  {
    bridge_list.clean_out();
    closest_coedge->get_children_virt( bridge_list );
    Curve* curve_ptr = dynamic_cast<Curve*>(bridge_list.get());
    assert( !!curve_ptr );
    double u = curve_ptr->u_from_position( curve_pt );

    if( !curve_ptr->G1_discontinuous( u, &tangent1, &tangent2 ) )
    {
      curve_ptr->closest_point( curve_pt, junk, &tangent1 );
      bool inside = is_inside( tangent1, curve_pt, surf_pt, normal );
      result = inside ? surf_pt : curve_pt;
      return CUBIT_SUCCESS;
    }

    if( closest_coedge->sense() == CUBIT_REVERSED )
    {
      tangent1 = -tangent1;
      tangent2 = -tangent2;
    }
  }
  else
  {
    bridge_list.clean_out();
    closest_coedge->get_children_virt( bridge_list );
    Curve* curve1 = dynamic_cast<Curve*>(bridge_list.get());
    bridge_list.clean_out();
    other_coedge->get_children_virt( bridge_list );
    Curve* curve2 = dynamic_cast<Curve*>(bridge_list.get());
    assert(curve1 && curve2);

    curve1->closest_point( curve_pt, junk, &tangent1 );
    curve2->closest_point( curve_pt, junk, &tangent2 );
    if( closest_coedge->sense() == CUBIT_REVERSED ) tangent1 = -tangent1;
    if( other_coedge->sense() == CUBIT_REVERSED ) tangent2 = -tangent2;
  }
    
  surf_ptr->closest_point( curve_pt, 0, &coe_normal );
  cross = tangent1 * tangent2;
  bool inside1 = is_inside( tangent1, curve_pt, surf_pt, normal );
  bool inside2 = is_inside( tangent2, curve_pt, surf_pt, normal );
  
  if ( (cross % coe_normal) > 0.0 )
    inside = inside1 && inside2;
  else
    inside = inside1 || inside2;

  result = inside ? surf_pt : curve_pt;
  return CUBIT_SUCCESS;
*/
}
/*
bool CompositeSurface::is_inside( const CubitVector& tangent,
                                  const CubitVector& curve_pt,
                                  const CubitVector& surf_pt,
                                  const CubitVector& normal )
{
  CubitVector cross = tangent * ( surf_pt - curve_pt );
  return cross % normal >= 0.0;
}
*/


//-------------------------------------------------------------------------
// Purpose       : Get the magnitudes of the principal curvatures.
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 03/10/98
//-------------------------------------------------------------------------
CubitStatus CompositeSurface::principal_curvatures(
  CubitVector const& location,
  double& curvature_1,
  double& curvature_2,
  CubitVector* closest_location )
{
  if (num_surfs() == 1)
    return get_surface(0)->
      principal_curvatures(location, curvature_1, curvature_2, closest_location);
  
  CubitVector curvature1, curvature2;
  CubitStatus s = closest_point( location, closest_location, NULL, 
                                 &curvature1, &curvature2 );
  if( s == CUBIT_FAILURE ) return CUBIT_FAILURE;
	
  curvature_1 = curvature1.length();
  curvature_2 = curvature2.length();
  return CUBIT_SUCCESS;
}


//-------------------------------------------------------------------------
// Purpose       : Evaluate the parameter values to get a position on the
//				   surface.
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 03/10/98
//-------------------------------------------------------------------------
CubitVector CompositeSurface::position_from_u_v( double u, double v )
{
  if (num_surfs() == 1)
    return get_surface(0)->position_from_u_v(u,v);
    
  PRINT_ERROR("CompositeSurface::position_from_u_v for non-paramtric surface.\n");
  CubitVector nulvect( 0., 0., 0.);
  return nulvect;
}


//-------------------------------------------------------------------------
// Purpose       : Determine values of u and v which evaluate to the passed
//                 position on the surface.
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 03/10/98
//-------------------------------------------------------------------------
CubitStatus CompositeSurface::u_v_from_position( 
  CubitVector const& pos ,
  double& u,
  double& v,
  CubitVector* closest )
{
  if (num_surfs() == 1)
    return get_surface(0)->u_v_from_position(pos, u, v, closest);

  PRINT_ERROR("CompositeSurface::u_v_from_position for non-paramtric surface.\n");
	u = v = 0.0;
	return CUBIT_FAILURE;
}


//-------------------------------------------------------------------------
// Purpose       : Is the prameterization of the surface periodic in nature?
//
// Special Notes : Always false, because the surface is not parametric.
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 03/10/98
//-------------------------------------------------------------------------
CubitBoolean CompositeSurface::is_periodic()
{
  if (num_surfs() == 1)
    return get_surface(0)->is_periodic();
    
  return CUBIT_FALSE;
}

//-------------------------------------------------------------------------
// Purpose       : Is the prameterization of the surface periodic in nature?
//
// Special Notes : Always false, because the surface is not parametric.
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 03/10/98
//-------------------------------------------------------------------------
CubitBoolean CompositeSurface::is_periodic_in_U( double& period )
{
  if (num_surfs() == 1)
    return get_surface(0)->is_periodic_in_U(period);
    
	period = 0.0;
	return CUBIT_FALSE;
}
CubitBoolean CompositeSurface::is_periodic_in_V( double& period )
{
  if (num_surfs() == 1)
    return get_surface(0)->is_periodic_in_V(period);
    
	period = 0.0;
	return CUBIT_FALSE;
}

//-------------------------------------------------------------------------
// Purpose       : Is the parameterization singular at the passed parameter
//                 value.
//
// Special Notes : Never: surface is not parametric.
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 03/10/98
//-------------------------------------------------------------------------
CubitBoolean CompositeSurface::is_singular_in_U( double param )
{
  if (num_surfs() == 1)
    return get_surface(0)->is_singular_in_U(param);
    
	return CUBIT_FALSE;
}
CubitBoolean CompositeSurface::is_singular_in_V( double param )
{
   if (num_surfs() == 1)
    return get_surface(0)->is_singular_in_V(param);
    
 return CUBIT_FALSE;
}


//-------------------------------------------------------------------------
// Purpose       : Is the surface closed along a param.
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 03/10/98
//-------------------------------------------------------------------------
CubitBoolean CompositeSurface::is_closed_in_U()
{
  if (num_surfs() == 1)
    return get_surface(0)->is_closed_in_U();
    
  return CUBIT_FALSE;
}
CubitBoolean CompositeSurface::is_closed_in_V()
{
  if (num_surfs() == 1)
    return get_surface(0)->is_closed_in_V();
    
  return CUBIT_FALSE;
}


//-------------------------------------------------------------------------
// Purpose       : Find u and v derivitives at the specified location.
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 03/10/98
//-------------------------------------------------------------------------
CubitStatus CompositeSurface::uv_derivitives( double u,
                                              double v,
                                              CubitVector& du,
                                              CubitVector& dv )
{
  if (num_surfs() == 1)
    return get_surface(0)->uv_derivitives(u, v, du, dv);
  
  PRINT_ERROR("CompositeSurface::uv_derivitives for non-paramtric surface.\n");
	return CUBIT_FAILURE;
}


//-------------------------------------------------------------------------
// Purpose       : Is the surface parametric.
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 03/10/98
//-------------------------------------------------------------------------
CubitBoolean CompositeSurface::is_parametric()
{
  if (num_surfs() == 1)
    return get_surface(0)->is_parametric();

  return CUBIT_FALSE;
}


//-------------------------------------------------------------------------
// Purpose       : Get the range of a parameter.
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 03/10/98
//-------------------------------------------------------------------------
CubitBoolean CompositeSurface::get_param_range_U( double& lower,
                                                  double& upper )
{
  if (num_surfs() == 1)
    return get_surface(0)->get_param_range_U(lower, upper);
  
  lower = upper = 0;
	return CUBIT_FALSE;
}
CubitBoolean CompositeSurface::get_param_range_V( double& lower,
                                                  double& upper )
{
  if (num_surfs() == 1)
    return get_surface(0)->get_param_range_V(lower, upper);
  
  lower = upper = 0;
  return CUBIT_FALSE;
}


//-------------------------------------------------------------------------
// Purpose       : Check if the passed position is on the surface.
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 03/10/98
//-------------------------------------------------------------------------
CubitBoolean CompositeSurface::is_position_on(CubitVector& position)
{
  for( int i = 0; i < compGeom->num_entities(); i++ )
  {
    Surface* surf = get_surface(i);
    if( surf->is_position_on( position ) ) return CUBIT_TRUE;
  }
	
  return CUBIT_FALSE;
}

//-------------------------------------------------------------------------
// Purpose       : Check if the passed position is outside, inside or
//                 on the boundary of the surface.
//
// Special Notes : 
//
// Creator       : Steve Storm
//
// Creation Date : 12/01/00
//-------------------------------------------------------------------------
CubitPointContainment CompositeSurface::point_containment( const CubitVector &point )
{
  bool boundary = false;
  for( int i = 0; i < num_surfs(); i++ )
  {
    CubitPointContainment cpc = get_surface(i)->point_containment( point );
    switch( cpc )
    {
      case CUBIT_PNT_OUTSIDE:
      //case CUBIT_PNT_OFF: 
        break;
      case CUBIT_PNT_INSIDE:
      //case CUBIT_PNT_ON:
        return cpc;
      case CUBIT_PNT_BOUNDARY:
        boundary = true;
        break;
      case CUBIT_PNT_UNKNOWN:
      default:
        return CUBIT_PNT_UNKNOWN;
    }
  }
  
  if( boundary )
    return CUBIT_PNT_BOUNDARY;
  else
    return CUBIT_PNT_OUTSIDE;
}
/*
CubitPointContainment CompositeSurface::point_containment( CubitVector &point, 
                                                           double ,
                                                           double  )
{
   return point_containment( point );
}
*/
CubitPointContainment CompositeSurface::point_containment( double u, 
                                                           double v )
{
   if (num_surfs() == 1)
    return get_surface(0)->point_containment(u,v);
   // Set this up when uv parameters are defined for composite surfaces.
   return CUBIT_PNT_UNKNOWN;
}

//-------------------------------------------------------------------------
// Purpose       : Relative sense of Surface wrt geometry underneath.
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 03/10/98
//-------------------------------------------------------------------------
CubitSense CompositeSurface::get_geometry_sense()
{
  return CUBIT_FORWARD;
}
void CompositeSurface::reverse_sense()
{
  compGeom->reverse();
}		


void CompositeSurface::get_curves( DLIList<CompositeCurve*>& result )
{
  for( CompositeLoop* loop = firstLoop; loop; loop = loop->loopNext )
  {
    CompositeCoEdge* coedge = loop->first_coedge();
    while( coedge )
    {
      result.append_unique( coedge->get_curve() );
      coedge = loop->next_coedge( coedge );
    }
  }
}


bool CompositeSurface::has_hidden_entities() const
{
  return hiddenSet && !hiddenSet->is_empty();
}

GeometryType CompositeSurface::geometry_type() 
  { return UNDEFINED_SURFACE_TYPE; }

void CompositeSurface::get_hidden_curves( DLIList<Curve*>& curves )
{
  if( hiddenSet )
    hiddenSet->hidden_curves( curves );
}

void CompositeSurface::print_debug_info( const char* line_prefix,
                                         bool brief ) 
{
  if( line_prefix == 0 ) line_prefix = "";
  CompositeLoop* loop = 0;
  
  if( brief )
  {
    int count = 0;
    while ( (loop = next_loop(loop) ) != NULL )
      count++;
#ifdef TOPOLOGY_BRIDGE_IDS
    PRINT_INFO("%sCompositeSurface %d : %d loops ", line_prefix, get_id(), count );
    if ( num_surfs() == 1 )
      PRINT_INFO("%s %d\n", fix_type_name(typeid(*get_surface(0)).name()), get_surface(0)->get_id());
    else
      PRINT_INFO("%d surfaces.\n", num_surfs());

#else
    PRINT_INFO("%sCompositeSurface %p : %d loops ", line_prefix, this, count );
    if ( num_surfs() == 1 )
      PRINT_INFO("%s %d\n", fix_type_name(typeid(*get_surface(0)).name()), get_surface(0)->get_saved_id());
   //   PRINT_INFO("%s %p\n", fix_type_name(typeid(*get_surface(0)).name()), get_surface(0));
    else
      PRINT_INFO("%d surfaces.\n", num_surfs());
#endif
    return;
  }
  
  char* new_prefix = new char[strlen(line_prefix)+3];
  strcpy( new_prefix, line_prefix );
  strcat( new_prefix, "  " );
#ifdef TOPOLOGY_BRIDGE_IDS
  PRINT_INFO("%sCompositeSurface %d\n", line_prefix, get_id() );
#else
  PRINT_INFO("%sCompositeSurface %d\n", line_prefix, this->get_saved_id() );
 // PRINT_INFO("%sCompositeSurface %p\n", line_prefix, this );
#endif
  compGeom->print_debug_info( new_prefix );

  // Print out info about any surfaces we are ingoring
  // during evaluation.
  if(surfacesToIgnore.size() > 0)
  {
      PRINT_INFO("%sSurfaces which are ignored:\n", new_prefix);
      for(int k=surfacesToIgnore.size(); k--;)
      {
         PRINT_INFO("%sSurface: %d\n", 
                    new_prefix, surfacesToIgnore.get_and_step()->get_saved_id());
      }
  }
  if( hiddenSet ) hiddenSet->print_debug_info( new_prefix );
  else PRINT_INFO("%s  No Hidden Entities.\n", line_prefix );
  while( (loop = next_loop(loop) ) != NULL )
    loop->print_debug_info( new_prefix );
  delete [] new_prefix;
  
  update_facet_tool();
  if ( facetTool )
  {
    facetTool->debug_draw_facets();
//    bool* reversed = new bool[num_surfs()];
//    for (int i = 0; i < num_surfs(); i++ )
//      reversed[i] = get_sense(i) == CUBIT_REVERSED;
//    facetTool->consolidate_points(reversed, GEOMETRY_RESABS);
//    delete [] reversed;
  }
}

//-------------------------------------------------------------------------
// Purpose       : Get sense in Shell
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 08/27/02
//-------------------------------------------------------------------------
CubitSense CompositeSurface::get_shell_sense( ShellSM* shellsm_ptr ) const
{
  CompositeShell* shell = dynamic_cast<CompositeShell*>(shellsm_ptr);
  if( shell )
    return shell->find_sense(this);
  
  DLIList<TopologyBridge*> parents;
  for( int i = 0; i < num_surfs(); i++ )
  {
    parents.clean_out();
    get_surface(i)->get_parents( parents );
    if( parents.is_in_list( shellsm_ptr ) )
    {
      CubitSense result = get_surface(i)->get_shell_sense(shellsm_ptr);
      if( get_sense(i) == CUBIT_REVERSED )
      {
        if( result == CUBIT_FORWARD )
          result = CUBIT_REVERSED;
        else if( result == CUBIT_REVERSED )
          result = CUBIT_FORWARD;
      }
      return result;
    }
  }
  return CUBIT_UNKNOWN;
}

void CompositeSurface::notify_reversed( TopologyBridge* bridge )
{
  int index = compGeom->index_of(bridge);
  if( index >= 0 )
    compGeom->reverse_sense(index);
}

CubitStatus CompositeSurface::stitch( CompositeSurface* partner )
{
  if( this == partner || this->stitchPartner || partner->stitchPartner )
  {
    assert(0);
    return CUBIT_FAILURE;
  }
  
  this->stitchPartner = partner;
  partner->stitchPartner = this;
  return CUBIT_SUCCESS;
}

CompositeSurface* CompositeSurface::unstitch()
{
  CompositeSurface* result = this->stitchPartner;
  if( result )
    this->stitchPartner = result->stitchPartner = 0;
  return result;
}

//-------------------------------------------------------------------------
// Purpose       : Update for change in underlying topology
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 12/19/02
//-------------------------------------------------------------------------
void CompositeSurface::notify_topology_modified( TopologyBridge* bridge )
{
  DLIList<CompositeCurve*> new_curves;
  Surface* surf = dynamic_cast<Surface*>(bridge);
  assert( surf && index_of(surf) >= 0 );
  update_modified();
  update_modified( surf, new_curves );
  
//  for( int i = new_curves.size(); i--; )
//    if( !CompositeEngine::instance().restore_curve(new_curves.get_and_step()) )
//      assert(0);
}

void CompositeSurface::update_modified( Surface* surf,
                                        DLIList<CompositeCurve*>& new_curves )
{
  int i;
  //int i = index_of(surf);
  //assert(i >= 0);
  //CubitSense rel_sense = get_sense(i);
  
    // find any new coedges in the surface
  DLIList<TopologyBridge*> bridge_list;
  DLIList<LoopSM*> loops;
  DLIList<CoEdgeSM*> coedges;
  
  surf->get_children_virt( bridge_list );
  CAST_LIST( bridge_list, loops, LoopSM );
  assert( bridge_list.size() == loops.size() );
  bridge_list.clean_out();
  
  for( i = loops.size(); i--; )
  {
    loops.get_and_step()->get_children_virt( bridge_list );
    while( bridge_list.size() )
    {
      CoEdgeSM* coedge = dynamic_cast<CoEdgeSM*>(bridge_list.pop());
      assert(0 != coedge);
      coedges.append(coedge);
    }
  }
  
  for( i = coedges.size(); i--; )
  {
    CoEdgeSM* coedge = coedges.get_and_step();
    bridge_list.clean_out();
    coedge->get_children_virt( bridge_list );
    assert( bridge_list.size() == 1 );
    Curve* curve = dynamic_cast<Curve*>(bridge_list.get());
    assert( 0 != curve );
    
    CompositeCoEdge* ccoedge = dynamic_cast<CompositeCoEdge*>(coedge->owner());
    if (ccoedge)
    {
        // If replace-curve was already done for the curve 
        // when processing the other surface, the composite
        // coedge will have been created already.  Add it to
        // the hidden set.
      if( !ccoedge->owner() && ccoedge->get_curve()->owner() == &hidden_entities())
        hidden_entities().hide( ccoedge );
      
        // If the coedge is a composite, the curve must be one
        // already as well.  Done with this coedge.
      continue;
    }
      
      // Replace curve with composite, and hide composite curve 
      // and any new child points.
    CompositeCurve* ccurve = CompositeEngine::instance().replace_curve(curve);
    assert(0 != ccurve);
    new_curves.append(ccurve);
    hidden_entities().hide( ccurve );
    CompositePoint* start = ccurve->start_point();
    if( ! start->owner() )
      hidden_entities().hide(start);
    CompositePoint* end = ccurve->start_point();
    if( ! end->owner() )
      hidden_entities().hide(end);
    
      // CompositeCoEdge was created by replace_curve(..)
      // Add it to the HiddenEntitySet
    ccoedge = dynamic_cast<CompositeCoEdge*>(coedge->owner());
    assert(ccoedge && !ccoedge->owner());
    hidden_entities().hide( ccoedge );
  }
}


  
  
   
//-------------------------------------------------------------------------
// Purpose       : Update for split in underlying surface
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 12/19/02
//-------------------------------------------------------------------------
void CompositeSurface::notify_split( TopologyBridge* new_bridge,
                                     TopologyBridge* old_bridge )
{
  assert(!new_bridge->owner());
  
  Surface* old_surf = dynamic_cast<Surface*>(old_bridge);
  Surface* new_surf = dynamic_cast<Surface*>(new_bridge);
  assert( old_surf && new_surf );
  int old_surf_index = index_of(old_surf);
  assert(old_surf_index >= 0);
  
  compGeom->append(new_surf, compGeom->sense(old_surf_index));
  new_surf->owner(this);
  
  DLIList<CompositeCurve*> new_curves;
  update_modified();
  update_modified( old_surf, new_curves );
  update_modified( new_surf, new_curves );
  
//  for( int i = new_curves.size(); i--; )
//    if( !CompositeEngine::instance().restore_curve(new_curves.get_and_step()) )
//      assert(0);
}

bool CompositeSurface::is_dead_coedge( CompositeCoEdge* coedge )
{
  if (coedge->num_coedges() > 0)
    return false;
  
  CompositeCurve* curve = coedge->get_curve();
  if (!curve)
    return true;
  
  if (curve->num_curves() == 0) // point-curve
  {
    CompositePoint* comp = curve->start_point();
    assert(comp == curve->end_point());
    return !(comp->get_point());
  }
  
  return false;
}

void CompositeSurface::update_modified( )
{
    // search for dead CoEdge-Curve pairs
  DLIList<CoEdgeSM*> coedge_list;
  if ( hiddenSet )
    hiddenSet->hidden_coedges( coedge_list );
  
  while (coedge_list.size())
  {
    CompositeCoEdge* coedge = dynamic_cast<CompositeCoEdge*>(coedge_list.pop());
    if ( is_dead_coedge(coedge) )
      remove_dead_coedge(coedge);
  }
  
  CompositeLoop* loop = next_loop();
  while (loop)
  {
    CompositeLoop* next = next_loop(loop);
    CompositeCoEdge* coedge = loop->first_coedge();
    while (coedge && is_dead_coedge(coedge))
    {
      remove_dead_coedge(coedge);
      coedge = loop->first_coedge();
    }
    
    if (coedge)
    {
      coedge = coedge->next();
      while ( coedge != loop->first_coedge() )
      {
        CompositeCoEdge* next_coe = coedge->next();
        if (is_dead_coedge(coedge))
          remove_dead_coedge(coedge);
        coedge = next_coe;
      }
    }
    
    if (loop->first_coedge() == NULL)
    {
      remove(loop);
      delete loop;
    }
    
    loop = next;
  }
  
}

void CompositeSurface::remove_dead_coedge( CompositeCoEdge* coedge )
{
  assert(is_dead_coedge(coedge));

  CompositeCurve* curve = coedge->get_curve();
  assert(curve->num_curves() == 0);
  curve->remove(coedge);
  delete coedge;
  if (curve->next_coedge(NULL))
    return;

  CompositePoint* start = curve->start_point();
  CompositePoint* end = curve->end_point();
  curve->start_point(0);
  curve->end_point(0);
  delete curve;

  if ( start->next_curve(NULL) == NULL )
  {
    if ( start->get_point() )
      CompositeEngine::instance().restore_point_in_curve(start);
    else 
      delete start;
  }

  if ( end != start && end->next_curve(NULL) == NULL )
  {
    if ( end->get_point() )
      CompositeEngine::instance().restore_point_in_curve(end);
    else
      delete end;
  }
}

//-------------------------------------------------------------------------
// Purpose       : Get graphics from cached facet data
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 08/18/03
//-------------------------------------------------------------------------
CubitStatus CompositeSurface::get_graphics( GMem& gmem )
{
  if (!facetTool)
    update_facet_tool();
  
  if (!facetTool)
    return CUBIT_FAILURE;
  
  facetTool->graphics( GEOMETRY_RESABS, gmem);
  return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Purpose       : Update for transform
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 10/17/03
//-------------------------------------------------------------------------
void CompositeSurface::notify_transformed()
{
  if (facetTool)
  {
    delete facetTool;
    facetTool = 0;
  }
}



CubitStatus CompositeSurface::get_projected_distance_on_surface( CubitVector *pos1,
                                                                 CubitVector *pos2, 
                                                                 double &distance )
{
  if ( num_surfs() == 1)
      return get_surface(0)->get_projected_distance_on_surface( pos1, pos2, distance );
    else
      return CUBIT_FAILURE; 
  return CUBIT_FAILURE;
}
CubitStatus CompositeSurface::get_sphere_params
(
  CubitVector &center,
  double &radius
) const
{
  PRINT_ERROR("Currently, Cubit is unable to determine sphere parameters for CompositeSurfaces.\n");
  return CUBIT_FAILURE;
}

CubitStatus CompositeSurface::get_cone_params
(
   CubitVector &center,
   CubitVector &normal,
   CubitVector &major_axis,
   double &radius_ratio,
   double &sine_angle,
   double &cos_angle
) const
{
  PRINT_ERROR("Currently, Cubit is unable to determine cone parameters for CompositeSurfaces.\n");
  return CUBIT_FAILURE;
}

CubitStatus CompositeSurface::get_torus_params
(
  CubitVector &center,
  CubitVector &normal,
  double &major_radius,
  double &minor_radius
) const
{
  PRINT_ERROR("Currently, Cubit is unable to determine torus parameters for CompositeSurface.\n");
  return CUBIT_FAILURE;
}

CubitStatus CompositeSurface::get_nurb_params
(
  bool &rational,
  int &degree_u,
  int &degree_v,
  int &num_cntrl_pts_u,
  int &num_cntrl_pts_v,
  DLIList<CubitVector> &cntrl_pts,
  DLIList<double> &cntrl_pt_weights,
  DLIList<double> &u_knots,
  DLIList<double> &v_knots
) const
{
  PRINT_ERROR("Currently, Cubit is unable to determine nurbs parameters for CompositeSurface.\n");
  return CUBIT_FAILURE;
}
