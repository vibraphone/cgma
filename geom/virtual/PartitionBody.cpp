//-------------------------------------------------------------------------
// Filename      : PartitionBody.cpp
//
// Purpose       : 
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 02/14/03
//-------------------------------------------------------------------------

#include "PartitionBody.hpp"
#include "CubitTransformMatrix.hpp"
#include "CubitFacetData.hpp"
#include "CubitPointData.hpp"
#include "PartitionPoint.hpp"
#include "PartitionSurface.hpp"
#include "SegmentedCurve.hpp"
#include "PartitionCoEdge.hpp"
#include "PartitionLump.hpp"
#include "VirtualQueryEngine.hpp"

//-------------------------------------------------------------------------
// Purpose       : Constructor
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 02/14/03
//-------------------------------------------------------------------------
PartitionBody::PartitionBody( BodySM* body )
  : childList(0)
{ 
  SubEntitySet* set = new SubEntitySet( body,this ); 
  set->bodyNext = 0;
  set->bodyPtr = this;
}


//-------------------------------------------------------------------------
// Purpose       : Desturctor
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 02/14/03
//-------------------------------------------------------------------------
PartitionBody::~PartitionBody()
{ 
  assert(! sub_entity_set().bodyNext);
  sub_entity_set().bodyPtr = 0;
  
  while( SubEntitySet* dead = childList )
  {
    childList = dead->bodyNext;
    dead->bodyNext = 0;
    dead->bodyPtr = 0;
  }
}

//-------------------------------------------------------------------------
// Purpose       : Get transform matrix
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 02/14/03
//-------------------------------------------------------------------------
CubitStatus PartitionBody::get_transforms( CubitTransformMatrix& xforms )
{
  return real_body()->get_transforms( xforms );
}

//-------------------------------------------------------------------------
// Purpose       : Attributes
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 02/14/03
//-------------------------------------------------------------------------
void PartitionBody::append_simple_attribute_virt(CubitSimpleAttrib* csa)
  { real_body()->append_simple_attribute_virt(csa); }
void PartitionBody::remove_simple_attribute_virt(CubitSimpleAttrib* csa)
  { real_body()->remove_simple_attribute_virt(csa); }
void PartitionBody::remove_all_simple_attribute_virt()
  { real_body()->remove_all_simple_attribute_virt(); }
CubitStatus PartitionBody::get_simple_attribute(DLIList<CubitSimpleAttrib*>& list)
  { return real_body()->get_simple_attribute(list); }
CubitStatus PartitionBody::get_simple_attribute(const CubitString& name,
                                       DLIList<CubitSimpleAttrib*>& list)
  { return real_body()->get_simple_attribute(name,list); }

//-------------------------------------------------------------------------
// Purpose       : Which layer
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 02/14/03
//-------------------------------------------------------------------------
int PartitionBody::layer() const
  { return sub_entity_set().get_owner_layer(); }

//-------------------------------------------------------------------------
// Purpose       : get children
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 02/14/03
//-------------------------------------------------------------------------
void PartitionBody::get_parents_virt( DLIList<TopologyBridge*>& )
  { }
void PartitionBody::get_children_virt( DLIList<TopologyBridge*>& list )
  { real_body()->get_children_virt(list); }

//-------------------------------------------------------------------------
// Purpose       : get GQE
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 02/14/03
//-------------------------------------------------------------------------
GeometryQueryEngine* PartitionBody::get_geometry_query_engine() const
  { return VirtualQueryEngine::instance(); }
  
//-------------------------------------------------------------------------
// Purpose       : misc. junk from PartitionEntity
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 02/14/03
//-------------------------------------------------------------------------
void PartitionBody::reverse_sense() { assert(0); }
void PartitionBody::notify_split(FacetEntity*, FacetEntity* ) { assert(0); }
CubitStatus PartitionBody::save(CubitSimpleAttrib&) { return CUBIT_FAILURE; }


//-------------------------------------------------------------------------
// Purpose       : get bounding box
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 02/14/03
//-------------------------------------------------------------------------
CubitBox PartitionBody::bounding_box() const
{
  DLIList<TopologyBridge*> lumps;
  real_body()->get_children_virt(lumps);
  
  CubitBox result = dynamic_cast<Lump*>(lumps.step_and_get())->bounding_box();
  for( int i = 1; i < lumps.size(); i++ )
    result |= dynamic_cast<Lump*>(lumps.step_and_get())->bounding_box();
  
  return result;
}


//-------------------------------------------------------------------------
// Purpose       : Get BodySM
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 02/14/03
//-------------------------------------------------------------------------
BodySM* PartitionBody::real_body() const
  { return dynamic_cast<BodySM*>(partitioned_entity()); }


//-------------------------------------------------------------------------
// Purpose       : Add an entity to the list we need to transform
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 02/14/03
//-------------------------------------------------------------------------
void PartitionBody::add( SubEntitySet& set )
{
  assert(!set.bodyPtr);
  set.bodyPtr = this;
  set.bodyNext = childList;
  childList = &set;
}

//-------------------------------------------------------------------------
// Purpose       : Remove an entity from the list of entities to transform
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 02/14/03
//-------------------------------------------------------------------------
void PartitionBody::remove( SubEntitySet& set )
{
  assert(set.bodyPtr == this);
  if( childList == &set )
  {
    childList = set.bodyNext;
  }
  else
  {
    SubEntitySet* ptr = childList;
    while( ptr->bodyNext != &set )
    {
      ptr = ptr->bodyNext;
      assert(!!ptr);
    }
    ptr->bodyNext = set.bodyNext;
  }
  set.bodyNext = 0;
  set.bodyPtr = 0;
}

  
//-------------------------------------------------------------------------
// Purpose       : Remove all child partition geometry (body deleted)
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 04/04/03
//-------------------------------------------------------------------------
void PartitionBody::destroy_all_children()
{
  DLIList<PartitionLump*> lumps;
  DLIList<PartitionSurface*> surfs;
  DLIList<PartitionCoEdge*> coedges;
  DLIList<PartitionCurve*> curves;
  DLIList<PartitionPoint*> points;
  
  DLIList<PartitionEntity*> sub_geom, split_geom;
  for ( SubEntitySet* ptr = childList; ptr; ptr = ptr->next_in_body() )
  {
    ptr->get_sub_entities( sub_geom );
    ptr->get_lower_order( split_geom );
    sub_geom += split_geom;
    split_geom.clean_out();
    while( sub_geom.size() )
    { 
      PartitionEntity* ent = sub_geom.pop();
      if ( PartitionPoint* point = dynamic_cast<PartitionPoint*>(ent) )
        points.append(point);
      else if( PartitionCoEdge* coedge = dynamic_cast<PartitionCoEdge*>(ent) )
        coedges.append(coedge);
      else if( PartitionCurve* curve = dynamic_cast<PartitionCurve*>(ent) )
        curves.append(curve);
      else if( PartitionSurface* surf = dynamic_cast<PartitionSurface*>(ent) )
        surfs.append(surf);
      else if ( PartitionLump* lump = dynamic_cast<PartitionLump*>(ent) )
        lumps.append(lump);
      else
        assert(0);
    }
  }
  
  while( lumps.size() )
  {
    PartitionLump* lump = lumps.pop();
    while( PartitionShell* shell = lump->first_shell() )
    {
      lump->remove(shell);
      shell->remove_all_surfaces();
      delete shell;
    }
    delete lump;
  }

  while( surfs.size() )
  {
    PartitionSurface* surf = surfs.pop();
    while( PartitionLoop* loop = surf->next_loop(NULL) )
    {
      surf->remove(loop);
      loop->remove_all_coedges();
      delete loop;
    }
    delete surf;
  }

  while( coedges.size() )
  {
    PartitionCoEdge* coedge = coedges.pop();
    if (coedge->get_curve() != NULL)
       coedge->get_curve()->remove( coedge );
    delete coedge;
  }

  while( curves.size() )
  {
    PartitionCurve* curve = curves.pop();
    assert( !curve->next_coedge(NULL) );
    delete curve;
  }

  while( points.size() )
  {
    PartitionPoint* point = points.pop();
    assert( !point->num_curves() );
    delete point;
  }

  // All child SubEntitySets should have been deleted when they
  // became empty, and should be removed from the list in this
  // body as they are destroyed.  However, when the last SubEntitySet
  // in this body is destroyed, this body will be destroyed as well.
  // So don't do this check!
  //assert(!childList);
}

//-------------------------------------------------------------------------
// Purpose       : 
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 05/10/04
//-------------------------------------------------------------------------
CubitStatus PartitionBody::mass_properties( CubitVector& result, double& volume )
{
  DLIList<Lump*> lump_list;
  lumps( lump_list );
  
  DLIList<PartitionLump*> part_list;
  CAST_LIST( lump_list, part_list, PartitionLump );
  if (part_list.size() < lump_list.size())
    return real_body()->mass_properties( result, volume );
  
  CubitVector centroid(0.0, 0.0, 0.0), tmp_centroid;
  volume = 0.0;
  double tmp_volume;
  for (int i = part_list.size(); i--; )
  {
    if (CUBIT_FAILURE == 
        part_list.get_and_step()->mass_properties( tmp_centroid, tmp_volume ))
      return CUBIT_FAILURE;
    
    centroid += tmp_volume * tmp_centroid;
    volume += tmp_volume;
  }
  
  if (volume > CUBIT_RESABS)
  {
    result = centroid / volume;
  }
  else
  {
    result.set( 0.0, 0.0, 0.0 );
    volume = 0.0;
  }
  return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Purpose       : 
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 05/10/04
//-------------------------------------------------------------------------
CubitPointContainment PartitionBody::point_containment( const CubitVector& pos )
{
  DLIList<Lump*> lump_list;
  lumps( lump_list );
  
  DLIList<PartitionLump*> part_list;
  CAST_LIST( lump_list, part_list, PartitionLump );
  
  if (part_list.size() < lump_list.size())
  {
    return real_body()->point_containment( pos );
  }
  
  bool inside = false, on = false;
  part_list.reset();
  for (int i = part_list.size(); i--; )
  {
    switch( part_list.get_and_step()->point_containment( pos ) )
    {
      case CUBIT_PNT_INSIDE:
        inside = true;
        break;
      case CUBIT_PNT_BOUNDARY:
        on = true;
        break;
      case CUBIT_PNT_OUTSIDE:
        break;
      default:
        return CUBIT_PNT_UNKNOWN;
    }
  }
  
  if (inside)
    return CUBIT_PNT_INSIDE;
  else if(on)
    return CUBIT_PNT_BOUNDARY;
  else
    return CUBIT_PNT_OUTSIDE; 
}

//-------------------------------------------------------------------------
// Purpose       : 
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 05/26/04
//-------------------------------------------------------------------------
void PartitionBody::get_all_children( DLIList<PartitionEntity*>& list )
{
  DLIList<PartitionEntity*> tmp;
  for (SubEntitySet* ptr = childList; ptr; ptr = ptr->bodyNext )
  {
    tmp.clean_out();
    ptr->get_sub_entities( tmp );
    list += tmp;
    
    tmp.clean_out();
    ptr->get_lower_order( tmp );
    list += tmp;
  }
}
