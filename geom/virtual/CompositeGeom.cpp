//-------------------------------------------------------------------------
// Filename      : CompositeGeom.cpp
//
// Purpose       : Object used by CompositeSurface and CompositeCurve
//                 to manage underlying entities.
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 12/19/01
//-------------------------------------------------------------------------

#include "VGDefines.h"
#include "GeometryEntity.hpp"
#include "CompositeGeom.hpp"
#include "CompositeAttrib.hpp"
#include "DLIList.hpp"

  // used to print debug info:
#include "Lump.hpp"
#include "Surface.hpp"
#include "Curve.hpp"
#include "Point.hpp"

const char* const COMPOSITE_DATA_ATTRIB_NAME = "COMPOSITE_ATTRIB";
  

//-------------------------------------------------------------------------
// Purpose       : Constructor
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 12/19/01
//-------------------------------------------------------------------------
CompositeGeom::CompositeGeom( int size )
  : entityList(size), 
    currentIndex(0),
    firstIndex(-1),
    needToUpdate( true ),
    listHead(0)
{
    // initialization above both allocated space for size
    // entities, and set the current count in the list to
    // size.  We want the initial memory, but need to set
    // the count back to zero.
  entityList.size(0);
}

//-------------------------------------------------------------------------
// Purpose       : Destructor
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 12/19/01
//-------------------------------------------------------------------------
CompositeGeom::~CompositeGeom()
{
  while(listHead)
  {
    CompositeAttrib* dead = listHead;
    listHead = listHead->next;
    delete dead;
  }
}

//-------------------------------------------------------------------------
// Purpose       : Find the index of an entity
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 12/19/01
//-------------------------------------------------------------------------
int CompositeGeom::index_of( TopologyBridge* ptr ) const
{
  int i;
  for( i = entityList.size() - 1; i >= 0; i-- )
    if( entityList[i].entity == ptr )
      break;

  return i;
}

//-------------------------------------------------------------------------
// Purpose       : Insert an entry
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 12/19/01
//-------------------------------------------------------------------------
CubitStatus CompositeGeom::insert( int index, GeometryEntity* geom_ptr, 
                                   CubitSense sense )
{
  if( index < 0 )
  {
    assert( index >= 0 );
    index = 0;
  }
  else if( index > entityList.size() )
  {
    assert( index <= entityList.size() );
    index = entityList.size();
  }
  
  CompositeEntry ent;
  ent.entity  = geom_ptr;
  ent.sense   = sense;
  ent.dist_sqr = ent.measure = 0.;

  //force 0th surface to be one that has the composite attrib on it.
  DLIList<CubitSimpleAttrib*> list;
  geom_ptr->get_simple_attribute(COMPOSITE_DATA_ATTRIB_NAME,list);
  if( list.size() )
    index = 0;

  entityList.insert( ent, index );

  update_cached_data();
   
  return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Purpose       : Remove an entry
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 12/19/01
//-------------------------------------------------------------------------
CubitStatus CompositeGeom::remove( int index, bool dead )
{
  if( (index < 0) || (index >= entityList.size()) )
  {
    assert( index >= 0 && index < entityList.size() );
    return CUBIT_FAILURE;
  }
  
  if (!dead)
    clean_up_attribs(entityList[index].entity);
  
  entityList.remove( index );
  update_cached_data();
  
  return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Purpose       : Split this composite into two, where all entries
//                 before the passed index remain in this entity and all
//                 entries after the index become part of a new composite.
//
// Special Notes : The entry at the specified index remains in this.
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 12/19/01
//-------------------------------------------------------------------------
CompositeGeom* CompositeGeom::split( int index )
{
  int i;
  if( (index < 0) || (index >= entityList.size()) )
  {
    assert( index >= 0 && index < entityList.size() );
    return 0;
  }
  
    // find EntityName attribute
  CompositeAttrib* name_attrib = listHead;
  while (name_attrib && name_attrib->name() != "ENTITY_NAME")
    name_attrib = name_attrib->next;
      
  int first = index + 1;
  CompositeGeom* new_geom = new CompositeGeom( entityList.size() - first );
  
  for( i = first; i < entityList.size(); i++ )
    new_geom->append( entityList[i].entity, entityList[i].sense );
  
  entityList.size( first );
  update_cached_data();
  
    // copy entityname attrib to new entity
  if ( name_attrib )
  {
    assert(!new_geom->listHead);
    new_geom->listHead = new CompositeAttrib(*name_attrib);
  }
 
  return new_geom;
}
CompositeGeom* CompositeGeom::split( VGArray<int>& index_array )
{
  int i, j;
  
  if( index_array.size() == 0 )
    return 0;
  
  for( i = 0; i < index_array.size(); i++ )
    if( index_array[i] < 0 || index_array[i] >= entityList.size() )
    {
      assert(0);
      return 0;
    }
  
  
    // find EntityName attribute
  CompositeAttrib* name_attrib = listHead;
  while (name_attrib && name_attrib->name() != "ENTITY_NAME")
    name_attrib = name_attrib->next;
  
  CompositeGeom* new_geom = new CompositeGeom( index_array.size() );
  for( i = 0; i < index_array.size(); i++ )
  {
    int index = index_array[i];
    assert( entityList[index].entity != NULL );
    new_geom->append( entityList[index].entity, entityList[index].sense );
    entityList[index].entity = 0;
  }
    
  for( i = 0; i < entityList.size() && entityList[i].entity; i++ );
  for( j = i + 1; j < entityList.size(); j++ )
  {
    if( entityList[j].entity )
    {
      entityList[i].entity = entityList[j].entity;
      entityList[i].sense  = entityList[j].sense;
      entityList[j].entity = 0;
      i++;
    }
  }
  entityList.size( i );
  
    // copy entityname attrib to new entity
  if ( name_attrib )
  {
    assert(!new_geom->listHead);
    new_geom->listHead = new CompositeAttrib(*name_attrib);
  }
  
  update_cached_data();
  return new_geom;
}

//-------------------------------------------------------------------------
// Purpose       : change geometry entity
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 01/11/02
//-------------------------------------------------------------------------
CubitStatus CompositeGeom::swap( int index, GeometryEntity* new_geom )
{
  if( (index < 0) || (index >= entityList.size()) )
  {
    assert( index >= 0 && index < entityList.size() );
    return CUBIT_FAILURE;
  }
  
  entityList[index].entity = new_geom;
  update_cached_data();
  return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Purpose       : Return union of bounding boxes of all entities.
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 12/19/01
//-------------------------------------------------------------------------
CubitBox CompositeGeom::bounding_box()
{
  if( entityList.size() == 0 )
  {
    return CubitBox();
  }
  
  if( needToUpdate )
    update_data();
    
  CubitBox box( entityList[0].bbox );
  for( int i = 1; i < entityList.size(); i++ )
    box |= entityList[i].bbox;
  return box;
}

//-------------------------------------------------------------------------
// Purpose       : Return sum of measure() for all entities.
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 12/19/01
//-------------------------------------------------------------------------
double CompositeGeom::measure() 
{
  if( entityList.size() == 0 )
    return 0.0;
  
  if( needToUpdate )
    update_data();

  return entityList[entityList.size()-1].measure;
}

//-------------------------------------------------------------------------
// Purpose       : Update cached data for each entity
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 12/19/01
//-------------------------------------------------------------------------
void CompositeGeom::update_data()
{
  needToUpdate = false;
  double sum = 0.0;
  for( int i = 0; i < entityList.size(); i++ )
  {
    sum += entityList[i].entity->measure();
    entityList[i].measure = sum;
    entityList[i].bbox    = entityList[i].entity->bounding_box();
  }
}

//-------------------------------------------------------------------------
// Purpose       : Setup search based on distance between passed position
//                 and the bounding box of each entity, and return the
//                 index of the entity with the bounding box closest to
//                 the passed position.
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 12/19/01
//-------------------------------------------------------------------------
int CompositeGeom::closest_box( const CubitVector& position )
{
  if( entityList.size() <= 0 ) return -1;
  
  if( needToUpdate ) update_data();
  
  int min_index = 0;
  double min_dist = entityList[0].dist_sqr 
    = entityList[0].bbox.distance_squared( position );
    
  for( int i = 1; i < entityList.size(); i++  )
  {
    CompositeEntry& ent = entityList[i];
    double dist_sqr = ent.dist_sqr = ent.bbox.distance_squared( position );
    if( dist_sqr < min_dist )
      min_index = i;
  }

  firstIndex = -1;
  currentIndex = min_index;
  return currentIndex;
}

//-------------------------------------------------------------------------
// Purpose       : Interate through entities with bounding boxes
//                 within the specified distance (squared) from the
//                 position passed to closest_box().
//
// Special Notes : -1 is returned ONCE to indicate the end of the list.
//                 Behavior is undefined if closest_box() has never been
//                 called.
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 12/19/01
//-------------------------------------------------------------------------
int CompositeGeom::next_box_within_dist( double dist_squared )
{
  if( entityList.size() > 0 ) 
  {
    if( firstIndex < 0 )
    {
      firstIndex = currentIndex;
    }
    else if( firstIndex == currentIndex )
    {
      return -1;
    }
    
    while( (currentIndex = (currentIndex + 1) % entityList.size()) 
                != firstIndex )
    {
      if( entityList[currentIndex].dist_sqr < dist_squared )
      {
        return currentIndex;
      }
    }
  }
    
  return -1;
}

//-------------------------------------------------------------------------
// Purpose       : Reverse sense of composite
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 
//-------------------------------------------------------------------------
CubitStatus CompositeGeom::reverse()
{
  reverse_order();
  reverse_rel_senses();
  return CUBIT_SUCCESS;
}
CubitStatus CompositeGeom::reverse_order()
{
  int i, half = entityList.size() / 2;
  for( i = 0; i < half; i++ )
  {
    int j = entityList.size() - (i+1);
    GeometryEntity* temp_entity = entityList[i].entity;
    entityList[i].entity = entityList[j].entity;
    entityList[j].entity = temp_entity;
    CubitSense temp_sense = entityList[i].sense;
    entityList[i].sense = entityList[j].sense;
    entityList[j].sense = temp_sense;
  }
  update_cached_data();
  return CUBIT_SUCCESS;
}
CubitStatus CompositeGeom::reverse_rel_senses()
{
  for( int i = 0; i < entityList.size(); i++ )
  {
    entityList[i].sense = 
      entityList[i].sense == CUBIT_FORWARD  ? CUBIT_REVERSED :
      entityList[i].sense == CUBIT_REVERSED ? CUBIT_FORWARD  :
      CUBIT_UNKNOWN;
  }
  
  return CUBIT_SUCCESS;
}


//-------------------------------------------------------------------------
// Purpose       : Combine 
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 03/04/02
//-------------------------------------------------------------------------
CubitStatus CompositeGeom::merge( CompositeGeom& dead, bool prepend )
{
  int i;

  if (entityList.size() == 1)
  {
    DLIList<CubitSimpleAttrib*> list;
    entityList[0].entity->get_simple_attribute(list);
    list.reset();
    for (i = list.size(); i--; )
    {
      CubitSimpleAttrib* csa = list.get_and_step();
      if (csa->character_type() != COMPOSITE_DATA_ATTRIB_NAME)
        listHead = new CompositeAttrib(csa, listHead);
      delete csa;
    }
  }
  
    // find EntityName attribute
  CompositeAttrib* this_name = listHead;
  while (this_name && this_name->name() != "ENTITY_NAME")
    this_name = this_name->next;
  
    // merge entity name attributes
  CompositeAttrib* dead_name = dead.listHead;
  if (dead_name)
  {
    if (dead_name->name() == "ENTITY_NAME")
    {
      dead_name = dead.listHead;
      dead.listHead = dead_name->next;
      dead_name->next = 0;
    }
    else 
    {
      while(dead_name->next && dead_name->next->name() != "ENTITY_NAME")
        dead_name = dead_name->next;
      if(dead_name->next)
      {
        CompositeAttrib* prev = dead_name;
        prev->next = dead_name = dead_name->next;
        dead_name->next = 0;
      }
    }
  }
  

  int insert ;
  if ( prepend )
  {
    insert = 0;
    entityList.size_end( entityList.size() + dead.entityList.size() );
  }
  else
  {
    insert = entityList.size();
    entityList.size( entityList.size() + dead.entityList.size() );
  }
  
  for( i = 0; i < dead.entityList.size(); i++ )
  {
    entityList[insert].entity = dead.entityList[i].entity;
    entityList[insert].sense  = dead.entityList[i].sense;
    insert++;
  }

  dead.entityList.size(0);
  update_cached_data();

  return CUBIT_SUCCESS;
}


//-------------------------------------------------------------------------
// Purpose       : Store an attribute
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 06/18/02
//-------------------------------------------------------------------------
void CompositeGeom::add_attribute( CubitSimpleAttrib* csa )
{
  if (entityList.size() == 1)
    entityList[0].entity->append_simple_attribute_virt(csa);
  else
    listHead = new CompositeAttrib( csa, listHead );
}

//-------------------------------------------------------------------------
// Purpose       : Remove an attribute
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 06/18/02
//-------------------------------------------------------------------------
void CompositeGeom::rem_attribute( CubitSimpleAttrib* csa ) 
{
  while (listHead && listHead->equals(csa))
  {
    CompositeAttrib* dead = listHead;
    listHead = dead->next;
    delete dead;
  }
  
  if (listHead)
  {
    CompositeAttrib* attrib = listHead;
    while (attrib->next)
    {
      if(attrib->next->equals(csa))
      {
        CompositeAttrib* dead = attrib->next;
        attrib->next = dead->next;
        delete dead;
      }
      else
      {
        attrib = attrib->next;
      }
    }
  }
  
  if (entityList.size() == 1)
    entityList[0].entity->remove_simple_attribute_virt(csa);
}

//-------------------------------------------------------------------------
// Purpose       : Remove all attributes
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 06/18/02
//-------------------------------------------------------------------------
void CompositeGeom::rem_all_attributes() 
{
  while(listHead)
  {
    CompositeAttrib* dead = listHead;
    listHead = listHead->next;
    delete dead;
  }
  
  if (entityList.size() == 1)
    entityList[0].entity->remove_all_simple_attribute_virt();
}


//-------------------------------------------------------------------------
// Purpose       : Get all attributes
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 06/18/02
//-------------------------------------------------------------------------
void CompositeGeom::get_attributes( DLIList<CubitSimpleAttrib*>& list ) 
{
    // special case: single-entity 'composite'
  if (entityList.size() == 1)
  {
    TopologyBridge* entity = entityList[0].entity;
    entity->get_simple_attribute(list);
    
      // handle 8.1 attribs on single-entity 'composites'
    for (int i = list.size(); i--; )
    {
      CubitSimpleAttrib* attrib = list.step_and_get();
      attrib->int_data_list()->reset();
      if (attrib->character_type() == COMPOSITE_DATA_ATTRIB_NAME &&
          *attrib->int_data_list()->get() == 1)
      {
        entity->remove_simple_attribute_virt(attrib);
        attrib->string_data_list()->reset();
        attrib->int_data_list()->reset();
        delete attrib->string_data_list()->remove();
        delete attrib->int_data_list()->remove();
        entity->append_simple_attribute_virt(attrib);
      }
    }
  }
      
    
  for (CompositeAttrib* ptr = listHead; ptr; ptr = ptr->next)
    list.append(ptr->csa());
}

//-------------------------------------------------------------------------
// Purpose       : Get named attributes
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 03/03/03
//-------------------------------------------------------------------------
void CompositeGeom::get_attributes( const char* name,
                                    DLIList<CubitSimpleAttrib*>& list ) 
{
  if (entityList.size() == 1)
  {
      // handle 8.1 attribs on single-entity 'composites'
    list.clean_out();
    entityList[0].entity->get_simple_attribute(COMPOSITE_DATA_ATTRIB_NAME,list);
    while (list.size())
    {
      CubitSimpleAttrib* attrib = list.pop();
      attrib->int_data_list()->reset();
      if (*attrib->int_data_list()->get() == 1)
      {
        entityList[0].entity->remove_simple_attribute_virt(attrib);
        attrib->string_data_list()->reset();
        attrib->int_data_list()->reset();
        delete attrib->string_data_list()->remove();
        delete attrib->int_data_list()->remove();
        entityList[0].entity->append_simple_attribute_virt(attrib);
      }
      delete attrib;
    }
    
    entityList[0].entity->get_simple_attribute(name, list);
  }
    
  for (CompositeAttrib* ptr = listHead; ptr; ptr = ptr->next)
    if (ptr->name() == name)
      list.append(ptr->csa());
}

//-------------------------------------------------------------------------
// Purpose       : Write info about object for debugging
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 
//-------------------------------------------------------------------------
void CompositeGeom::print_debug_info( const char* line_prefix )
{
  if( needToUpdate )
    update_data();
  if( line_prefix == 0 )
    line_prefix = "";
  
  PRINT_INFO("%sCompositeGeom @ %p : \n", line_prefix, this );
  for( int i = 0; i < entityList.size(); i++ )
  {
    GeometryEntity* ptr = entityList[i].entity;
#ifdef TOPOLOGY_BRIDGE_IDS
    PRINT_INFO("%s  %15s %d %7s\n", line_prefix, 
      ptr ? fix_type_name(typeid(*ptr).name()) : "GeometryEntity",
      ptr ? ptr->get_id() : 0, 
      entityList[i].sense == CUBIT_FORWARD ? "Forward" :
      entityList[i].sense == CUBIT_REVERSED ? "Reverse" :
      "Unknown");
#else    
    /*
    PRINT_INFO("%s  %15s %p %7s\n", line_prefix, 
      ptr ? fix_type_name(typeid(*ptr).name()) : "GeometryEntity",
      ptr, 
      entityList[i].sense == CUBIT_FORWARD ? "Forward" :
      entityList[i].sense == CUBIT_REVERSED ? "Reverse" :
      "Unknown");
      */
    PRINT_INFO("%s  %15s %d %7s\n", line_prefix, 
      ptr ? fix_type_name(typeid(*ptr).name()) : "GeometryEntity",
      ptr ? ptr->get_saved_id() : 0, 
      entityList[i].sense == CUBIT_FORWARD ? "Forward" :
      entityList[i].sense == CUBIT_REVERSED ? "Reverse" :
      "Unknown");
            
#endif
  }
}

//-------------------------------------------------------------------------
// Purpose       : Reverse geometric sense of composite
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 
//-------------------------------------------------------------------------
void CompositeGeom::reverse_sense( int index )
{
  assert(index < entityList.size() );
  CubitSense old = entityList[index].sense;
  entityList[index].sense = (old == CUBIT_REVERSED) ? CUBIT_FORWARD : CUBIT_REVERSED;
}

//-------------------------------------------------------------------------
// Purpose       : Read and remove attributes from underlying entities
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 06/30/03
//-------------------------------------------------------------------------
void CompositeGeom::read_attributes( GeometryEntity* geom_ptr )
{
  DLIList<CubitSimpleAttrib*> list;
  int i;

    // remove any attributes from previous read
  rem_all_attributes();

  if (geom_ptr)
  {
      // Special case for point-curves (no real curves to write
      // attirbutes to.)  Write to passed entity instead.
    assert(entityList.size() == 0);
    geom_ptr->get_simple_attribute(COMPOSITE_DATA_ATTRIB_NAME,list);
    
    list.reset();
    for (i = list.size(); i--; )
    {
      CubitSimpleAttrib* attrib = list.get_and_step();
      attrib->int_data_list()->reset();
      assert(attrib->int_data_list()->size());
      if (*attrib->int_data_list()->get() == entityList.size())
      {
        geom_ptr->remove_simple_attribute_virt(attrib);
        delete attrib->int_data_list()->remove();
        attrib->string_data_list()->reset();
        delete attrib->string_data_list()->remove();
        listHead = new CompositeAttrib(attrib,listHead);
        delete attrib;
      }
    }
    
    return;
  }

  int index_of_entity_with_attribs = -1;

  for (i = 0; i < entityList.size(); i++)
  {
    list.clean_out();
    entityList[i].entity->get_simple_attribute(COMPOSITE_DATA_ATTRIB_NAME,list);

    if( list.size() )
      index_of_entity_with_attribs = i;
  
    list.reset();
    for (int j = list.size(); j--; )
    {
      CubitSimpleAttrib* attrib = list.get_and_step();
      attrib->int_data_list()->reset();
      assert(attrib->int_data_list()->size());
      if (*attrib->int_data_list()->get() == entityList.size())
      {
        // Take the attributes off of the current entity and put them on the first entity
        // in this list.  I believe this is ok to do because the attributes should apply to
        // the whole composite surface and not just the underlying entity they are on 
        // (the one exception to this might be UNIQUE_ID but I haven't seen any problems
        // with this yet).  The reason for doing this is that there is some code (I believe
        // in uncomposite() that assumes any attributes will be on the first entity
        // in the list.  Previous code actually moved the entity to the beginning
        // of the list but this reordering of the list does not fly with composite
        // curves because there is code depending on the curves in the list being
        // ordered so that they connect end to end in a contiguous manner (the 
        // faceting code, for one, relies on this).  BWC 1/7/07.
        entityList[i].entity->remove_simple_attribute_virt(attrib);
        entityList[0].entity->append_simple_attribute_virt(attrib);

        delete attrib->int_data_list()->remove();
        attrib->string_data_list()->reset();
        delete attrib->string_data_list()->remove();
        listHead = new CompositeAttrib(attrib,listHead);
        delete attrib;
      }
    }
  }
}

//-------------------------------------------------------------------------
// Purpose       : Save attributes on first underlying entity
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 06/30/03
//-------------------------------------------------------------------------
void CompositeGeom::write_attributes( GeometryEntity* geom_ptr )
{
  DLIList<CubitSimpleAttrib*> list;

  if (geom_ptr)
  {
      // Special case for point-curves (no real curves to write
      // attirbutes to.)  Write to passed entity instead.
    assert(entityList.size() == 0);
 
      // clean up any attributes from the previous write
    geom_ptr->get_simple_attribute(COMPOSITE_DATA_ATTRIB_NAME,list);
    while (list.size())
    {
      CubitSimpleAttrib* csa = list.pop();
      geom_ptr->remove_simple_attribute_virt(csa);
      delete csa;
    }
  }
  else
  {
    geom_ptr = entityList[0].entity;
 
      // clean up any attributes from the previous write
    for (int i = 0; i < entityList.size(); i++)
    {
      entityList[i].entity->get_simple_attribute(COMPOSITE_DATA_ATTRIB_NAME,list);
      while (list.size())
      {
        CubitSimpleAttrib* csa = list.pop();
        entityList[i].entity->remove_simple_attribute_virt(csa);
        delete csa;
      } 
    }
  }
  
  
  CubitString name = COMPOSITE_DATA_ATTRIB_NAME;
  int count = entityList.size();
  CubitSimpleAttrib attrib;
  
  for (CompositeAttrib* ptr = listHead; ptr; ptr = ptr->next)
  {
    attrib.string_data_list()->append(&name);
    attrib.int_data_list()->append(&count);
    
    ptr->append_to_csa(&attrib);
    geom_ptr->append_simple_attribute_virt(&attrib);
    
    attrib.string_data_list()->clean_out();
    attrib.int_data_list()->clean_out();
    attrib.double_data_list()->clean_out();
  }
}

//-------------------------------------------------------------------------
// Purpose       : Clean out composite attributes.
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 07/01/03
//-------------------------------------------------------------------------
void CompositeGeom::clean_up_attribs( GeometryEntity* ent )
{
  DLIList<CubitSimpleAttrib*> list;
  ent->get_simple_attribute(COMPOSITE_DATA_ATTRIB_NAME,list);
  while (list.size())
  {
    CubitSimpleAttrib* csa = list.pop();
    ent->remove_simple_attribute_virt(csa);
    delete csa;
  }
}
