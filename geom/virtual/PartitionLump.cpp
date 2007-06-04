#include "PartitionLump.hpp"
#include "PartitionShell.hpp"
#include "PartitionSurface.hpp"
#include "VirtualQueryEngine.hpp"
#include "PartitionCurve.hpp"
#include "PartitionBody.hpp"
#include "CubitFacetData.hpp"

PartitionLump::PartitionLump( Lump* real_lump )
  : listHead(0)
{
  assert( dynamic_cast<SubEntitySet*>(real_lump->owner()) == 0 );
  new SubEntitySet( real_lump, this );
}

PartitionLump::PartitionLump( PartitionLump* split_from )
  : listHead(0)
{
  split_from->sub_entity_set().add_partition( this );
}

PartitionLump::~PartitionLump()
{
  remove_all_shells();
}

Lump* PartitionLump::real_lump() const
{
  return dynamic_cast<Lump*>(partitioned_entity());
}


CubitStatus PartitionLump::add( PartitionShell* shell )
{
  assert( shell->myLump == 0 );
  shell->lumpNext = listHead;
  shell->myLump = this;
  listHead = shell;
  return CUBIT_SUCCESS;
}

CubitStatus PartitionLump::remove( PartitionShell* shell )
{
  if( shell->myLump != this )
    return CUBIT_FAILURE;
  
  if( listHead == shell )
  {
    listHead = shell->lumpNext;
  }
  else
  {
    PartitionShell* sh = listHead;
    while( sh && sh->lumpNext != shell )
      sh = sh->lumpNext;
    assert( sh && sh->lumpNext == shell );
    sh->lumpNext = shell->lumpNext;
  }
  
  shell->myLump = 0;
  shell->lumpNext = 0;
  return CUBIT_SUCCESS;
}

void PartitionLump::remove_all_shells()
{
  while( listHead )
  {
    CubitStatus s = remove( listHead );
    assert( s );
  }
}

void PartitionLump::get_parents_virt( DLIList<TopologyBridge*>& list )
{
  if( get_body() )
    list.append( get_body() );
  else
  {
    real_lump()->get_parents_virt( list );
//    PartitionEngine::fix_up_query_results( list );
  }
}

void PartitionLump::get_children_virt( DLIList<TopologyBridge*>& list )
{
  PartitionShell* shell = first_shell();
  while( shell )
  {
    list.append( shell );
    shell = next_shell( shell );
  }
}

GeometryQueryEngine* PartitionLump::get_geometry_query_engine() const
{
  return VirtualQueryEngine::instance();
}



void PartitionLump::append_simple_attribute_virt(CubitSimpleAttrib* csa)
{ sub_entity_set().add_attribute( this, csa ); }
void PartitionLump::remove_simple_attribute_virt(CubitSimpleAttrib* csa)
{ sub_entity_set().rem_attribute( this, csa ); }
void PartitionLump::remove_all_simple_attribute_virt()
{ sub_entity_set().rem_all_attrib( this ); }
CubitStatus PartitionLump::get_simple_attribute(DLIList<CubitSimpleAttrib*>& list)
{ 
  sub_entity_set().get_attributes( this, list ); 
  return CUBIT_SUCCESS;
}
CubitStatus PartitionLump::get_simple_attribute(const CubitString& name,
                                       DLIList<CubitSimpleAttrib*>& list)
{ 
  sub_entity_set().get_attributes( this, name.c_str(), list ); 
  return CUBIT_SUCCESS;
}

void PartitionLump::notify_split( FacetEntity* , FacetEntity* )
  { assert(0); }

void PartitionLump::reverse_sense()
{
  PartitionShell* shell = 0;
  while( (shell = next_shell(shell)) )
  {
    PartitionCoSurf* cosurf = 0;
    while( (cosurf = shell->next_co_surface( cosurf )) )
    {
      bool reversed = (cosurf->sense() == CUBIT_REVERSED);
      cosurf->sense( reversed ? CUBIT_FORWARD : CUBIT_REVERSED );
    }
  }
}

void PartitionLump::transform(const CubitTransformMatrix&) {;}    

CubitBox PartitionLump::bounding_box() const
{
  return real_lump() ? real_lump()->bounding_box() : CubitBox();
}

double PartitionLump::measure() 
{
    // Calculate volume of this lump from surface facets.
    // Volume is calculated as the sum of the signed
    // volumes of the tetrahedrons formed by this point
    // with triagle.
  CubitVector p0 = bounding_box().center();
  
  double result = 0.0;
  PartitionShell* shell = 0;
  PartitionCoSurf* cosurf = 0;
  DLIList<CubitFacetData*> facets;
  CubitVector p1, normal;
  
    // calculate area around this point
  while ( (shell = next_shell(shell)) )
  {
      // mark surfaces with a count of how many times they 
      // occur in the shell
    cosurf = 0;
    while ( (cosurf = shell->next_co_surface(cosurf)) )
      cosurf->get_surface()->mark = 0;
    cosurf = 0;
    while ( (cosurf = shell->next_co_surface(cosurf)) )
      cosurf->get_surface()->mark++;
    
      // calculate partial area for each surface
    cosurf = 0;
    while ( (cosurf = shell->next_co_surface(cosurf)) )
    {
      PartitionSurface* surf = cosurf->get_surface();
      
        // skip non-manifold surfaces
      if( surf->mark > 1 ) continue;
      
      facets.clean_out();
      surf->get_facet_data( facets );
      for  ( int i = facets.size(); i--; )
      {
          // Tetrahedron volume is Ah/3
          // Calculate 2Ah and add to result.
          // Divide result by 6 when all done.
        
        CubitFacet* facet = facets.step_and_get();
        p1 = facet->point(0)->coordinates();
        normal = (facet->point(2)->coordinates() - p1) 
               * (facet->point(1)->coordinates() - p1);
        
          // Triangle area is 1/2 length of edge product
        double two_area = normal.length();
        
        if ( two_area > CUBIT_RESABS )
        {
            // Calculating signed area - need to 
            // reverse facet normal for reversed
            // cosurfaces.
          if ( cosurf->sense() == CUBIT_REVERSED )
            normal = -normal;
          
            // Make normal a unit vector.  Already
            // calculated length, so reuse that 
            // value.
          normal /= two_area;
          
            // This is where the signed part of the
            // calculation comes in.  If normal is 
            // in opposite direction as the vector
            // (p0-p1), then the height is negative.
          double height = normal % (p0 - p1);
          
            // Add the signed height value times 
            // twice the triangle area to result.
          result += two_area * height;
        }
      }
    }
    
      // clear marks
    cosurf = 0;
    while ( (cosurf = shell->next_co_surface(cosurf)) )
      cosurf->get_surface()->mark = 0;
  }
  
    // Calculated 2Ah for each tetrahedron, so divide
    // by 6 to get Ah/3.
  return result / 6.0;
}

PartitionBody* PartitionLump::get_body() const
  { return sub_entity_set().body(); }

TopologyBridge* PartitionLump::find_parent_body() const
{
  if( get_body() ) 
    return get_body();
  
  Lump* lump = real_lump();
  if(!lump) return 0;
  
  DLIList<TopologyBridge*> list;
  lump->get_parents_virt(list);
  return list.size() ? list.get() : 0;
}

CubitStatus PartitionLump::save( CubitSimpleAttrib& attrib )
{
  DLIList<int> surf_list;
  
  PartitionShell* shell = 0;
  while( (shell = next_shell(shell)) )
  {
    PartitionCoSurf* cosurf = 0;
    int cosurf_count = 0;
    while( (cosurf = shell->next_co_surface(cosurf)) )
      cosurf_count++;
    
    surf_list.append(cosurf_count);
    cosurf = 0;
    while( (cosurf = shell->next_co_surface(cosurf)) )
    {
      PartitionSurface* surf = cosurf->get_surface();
      int set_id = 0;
      int surf_id = surf->sub_entity_set().get_id(surf);
      if( &(surf->sub_entity_set()) != &sub_entity_set() )
        set_id = surf->sub_entity_set().get_unique_id();
    
      if ( cosurf->sense() == CUBIT_REVERSED )
        surf_id = -surf_id;
        
      surf_list.append( set_id );
      surf_list.append( surf_id );
    }
  }
  
  int id = sub_entity_set().get_id(this);
  return sub_entity_set().save_geometry( id, 3, 0, 0, &surf_list, 0, attrib );
}

void PartitionLump::get_all_children( DLIList<PartitionEntity*>& list )
{
  PartitionShell* shell = 0;
  while( (shell = next_shell(shell)) )
  {
    PartitionCoSurf* cosurf = 0;
    while( (cosurf = shell->next_co_surface(cosurf)) )
    {
      PartitionSurface* surface = cosurf->get_surface();
      list.append( surface );
      PartitionLoop* loop = 0;
      while( (loop = surface->next_loop(loop)) )
      {
        PartitionCoEdge* coedge = loop->first_coedge();
        do {
          PartitionCurve* curve = coedge->get_curve();
          list.append( curve );
          list.append( curve->start_point() );
          list.append( curve->end_point() );
          coedge = loop->next_coedge( coedge );
        } while( coedge != loop->first_coedge() );
      }
    }
  }
  
  int i;
  for( i = list.size(); i--; )
    list.step_and_get()->mark = 0;
  for( i = list.size(); i--; )
    list.step_and_get()->mark++;
  for( i = list.size(); i--; )
  {
    list.step_and_get()->mark--;
    if( list.get()->mark != 0 )
      list.change_to(0);
  }
  list.remove_all_with_value(0);
}

//-------------------------------------------------------------------------
// Purpose       : Mass properties
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 05/10/04
//-------------------------------------------------------------------------
CubitStatus PartitionLump::mass_properties( CubitVector& centroid, double& volume )
{
  PartitionShell* shell = 0;
  CubitVector s_cent;
  double s_vol;
  volume = 0.0;
  centroid.set( 0.0, 0.0, 0.0 );
  
  while ((shell = next_shell( shell )))
  {
    if (CUBIT_SUCCESS != shell->mass_properties( s_cent, s_vol ))
      return CUBIT_FAILURE;
    centroid += s_vol * s_cent;
    volume += s_vol;
  }
  
  if (volume > CUBIT_RESABS)
    centroid /= volume;
  else
    centroid.set( 0.0, 0.0, 0.0 );
  return CUBIT_SUCCESS;
}

CubitPointContainment PartitionLump::point_containment( const CubitVector& pos )
{
  PartitionCoSurf* closest_surf = NULL;
  double closest_dist = CUBIT_DBL_MAX;
  CubitVector closest, normal;
  
  PartitionShell* shell = 0;
  while ((shell = next_shell( shell )))
  {
    PartitionCoSurf* cosurf = 0;
    while ((cosurf = shell->next_co_surface( cosurf )))
    {
      PartitionSurface* surf = cosurf->get_surface();
      if (!shell->is_nonmanifold( surf ))
      {
        surf->closest_point( pos, &closest );
        double dist = (pos - closest).length_squared();
        if (dist < closest_dist)
        {
          closest_dist = dist;
          closest_surf = cosurf;
        }
      }
    }
  }
  
  if (!closest_surf)
    return CUBIT_PNT_UNKNOWN;
  
  closest_surf->get_surface()->closest_point( pos, &closest, &normal );
  if ((closest - pos).length_squared() < GEOMETRY_RESABS*GEOMETRY_RESABS)
    return CUBIT_PNT_BOUNDARY;
  
  if (closest_surf->sense() == CUBIT_REVERSED)
    normal = -normal;
  
  double dot = normal % (closest - pos);
  return dot < CUBIT_RESABS ? CUBIT_PNT_OUTSIDE :
         dot > CUBIT_RESABS ? CUBIT_PNT_INSIDE  :
                              CUBIT_PNT_UNKNOWN;
}

  
