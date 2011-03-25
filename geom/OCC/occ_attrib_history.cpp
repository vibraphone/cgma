
#include "occ_attrib_history.hpp"
#include "CubitMessage.hpp"
#include "OCCHistory.hpp"
#include "CubitSimpleAttrib.hpp"
#include "OCCAttribSet.hpp"

// ********** BEGIN ACIS INCLUDES             **********
#include "TopoDS_Shape.hxx"	
#include "BRepBuilderAPI_ModifyShape.hxx"
#include "BRepAlgoAPI_BooleanOperation.hxx"
#include "LocOpe_SplitShape.hxx"
#include "TopAbs_ShapeEnum.hxx"
// ********** END ACIS INCLUDES               **********

// ********** BEGIN MACRO DEFINITIONS         **********
#define THIS() OCC_ATTRIB_HISTORY
#define THIS_LIB NONE
#define PARENT() NONE
#define PARENT_LIB NONE

#define OCC_ATTRIB_HISTORY_NAME "history"

ATTRIB_DEF("history_attribute")

int OCC_ATTRIB_HISTORY::num_attribs = 0;
bool ATTRIB_HISTORY::addOrRemoveFromList = true;
std::set<OCC_ATTRIB_HISTORY*> OCC_ATTRIB_HISTORY::allHistoryAttribs; 

//When acis backs up this attribute, it gets to this macro right after exiting 
//the constructor.  Here we remove the backed-up history attribute from our list 
//since Acis takes care of destroying all backed-up attributes.
ATTCOPY_DEF("history_attribute")
  OCC_ATTRIB_HISTORY *history_attrib = const_cast<OCC_ATTRIB_HISTORY*>(rollback);
  std::set<OCC_ATTRIB_HISTORY*>::iterator iter;
  iter = allHistoryAttribs.find( history_attrib ); 
  if( iter != allHistoryAttribs.end() )
    allHistoryAttribs.erase( iter );

LOSE_DEF

DTOR_DEF
  num_attribs--;
  std::set<OCC_ATTRIB_HISTORY*>::iterator iter;
  if( addOrRemoveFromList )
  {
    iter = allHistoryAttribs.find( this );
    if( iter != allHistoryAttribs.end() )
      allHistoryAttribs.erase( iter );
  }
  if( num_attribs == 0 )
    PRINT_INFO("num_attribs = %d\n", num_attribs ); 

DEBUG_DEF


SAVE_DEF
// Don't save

RESTORE_DEF
// Don't restore

COPY_DEF
  //This attribute is being copied from TopoDS_Shape "from".  At this point
  //it would be nice to be able to map the original TopoDS_Shape 
  //(from->owner())
  //to the copy (this->owner()) but you can't since all the copying isn't 
  //complete.  So right now we just get the 'from' attribute.  Later in
  //FIX_POINTER_DEF after everything is hooked up correctly, we can get
  //the this->owner().
  fromHistoryAttrib = (OCC_ATTRIB_HISTORY*)from;
  //get the index of the copy in the array in fix_pointer_def 

SCAN_DEF

FIX_POINTER_DEF
  //change from index to actual value
  //ENTITY *thisENTITY = ( ENTITY* ) read_array( array, index );
  
  //we get in here in the last step of the copy...at this point 
  //everything should be hooked up nicely so that we can determine
  //the entity from which this attribute's owner was copied
  //AcisHistory::AcisEvent event;
  //event.eventType = AcisHistory::COPY;
  //append the original entity to the first list
  //int unique_id = fromHistoryAttrib->get_unique_id();
  //event.entities.push_back( unique_id ); 

  this->add_acis_history( fromHistoryAttrib->occHistory );
  this->trackingIds = fromHistoryAttrib->trackingIds;

  //Append the copied ENTITY to the 'other' list
  //event.other_entities.push_back( this->get_unique_id() ); 

  //add the event
  //fromHistoryAttrib->acisHistory->add_event( event ); 

TERMINATE_DEF
// Don't do anything special
// ********** END MACRO DEFINITIONS           **********

//bool ATTRIB_CUBIT_OWNER::copyingAttribs = true;
//bool ATTRIB_CUBIT_OWNER::splitCopy = false;



OCC_ATTRIB_HISTORY::OCC_ATTRIB_HISTORY( TopoDS_Shape *entity, 
                                        OCCHistory *occ_history )
{
  occHistory = NULL;
  num_attribs++;
  if( occ_history )
    add_occ_history( occ_history );

  if( addOrRemoveFromList )
    allHistoryAttribs.insert( this ); 

  DLIList<CubitString*> strings;
  DLIList<double> doubles;
  DLIList<int> ints;
  CubitSimpleAttrib* attb = new CubitSimpleAttrib(&strings, &doubles, &ints);
  OCCAttribSet::append_attribute(attb, *entity);
}

void ATTRIB_HISTORY::add_occ_history( OCCHistory *occ_history )
{
  occHistory = occ_history; 
}

std::set<int> ATTRIB_HISTORY::get_tracking_ids()
{
  return trackingIds;
}

void ATTRIB_HISTORY::split_owner( ENTITY *entity)
{
  /*
  PRINT_INFO("I'm in split_owner!!!!\n");
  if( is_BODY( this->entity() )) 
    PRINT_INFO("  It's a body\n");
  if( is_LUMP( this->entity() )) 
    PRINT_INFO("  It's a volume\n");
  if( is_FACE( this->entity() )) 
    PRINT_INFO("  It's a surface\n");
  if( is_EDGE( this->entity() )) 
    PRINT_INFO("  It's an edge\n");
  if( is_VERTEX( this->entity() )) 
    PRINT_INFO("  It's a vertex\n"); */

  //add an ATTRIB_HISTORY to the new entity
  ATTRIB_HISTORY *other_history_attrib = get_history_attrib( entity, true, acisHistory );
  other_history_attrib->trackingIds = this->trackingIds;

  if( trackingIds.size() > 1 )
    return;
  //Create an AcisEvent of type CUT 
  AcisHistory::AcisEvent event;
  event.eventType = CGMHistory::CUT;

  event.entities.push_back( *(trackingIds.begin()) ); 

  assert( acisHistory );
  
  //add the event to the history object
  acisHistory->add_event( event ); 

  //Create an AcisEvent of type Subdivision
  //AcisHistory::AcisEvent event;
  //event.eventType = AcisHistory::SUBDIVISION;
  //event.entities.push_back( this->get_unique_id() );

  //event.other_entities.push_back( other_history_attrib->get_unique_id() ); 
  //event.other_entities.push_back( this->get_unique_id() );

  //assert( acisHistory );
  
  //add the event to the history object
  //acisHistory->add_event( event ); 
}
  
void ATTRIB_HISTORY::merge_owner( ENTITY *entity, logical delete_this)
{
  if( delete_this )
    return;
/*
  PRINT_INFO("I'm in merge owner!!!!\n");
  if( is_BODY( this->entity() )) 
    PRINT_INFO("  It's a body\n");
  if( is_LUMP( this->entity() )) 
    PRINT_INFO("  It's a volume\n");
  if( is_FACE( this->entity() )) 
    PRINT_INFO("  It's a surface\n");
  if( is_EDGE( this->entity() )) 
    PRINT_INFO("  It's an edge\n");
  if( is_VERTEX( this->entity() )) 
    PRINT_INFO("  It's a vertex\n");
 */

  //if the entity we're merging this's->owner() with doesn't have
  //a history attrib on it, we don't care about this
  ATTRIB_HISTORY *other_history_attrib = get_history_attrib( entity, false ); 
  if( !other_history_attrib )
    return;

  //Create an AcisEvent of type Subdivision
  //AcisHistory::AcisEvent event;
  //event.eventType = AcisHistory::ABSORPTION;

  //'entity' got absorbed into 'this->entity()'
  //event.entities.push_back( other_history_attrib->get_unique_id() ); 
  //event.entities.push_back( this->get_unique_id() );
  //event.other_entities.push_back( this->get_unique_id() );

  std::set<int> other_tracking_ids = other_history_attrib->trackingIds;
  std::set<int>::iterator iter;
  for(iter = other_tracking_ids.begin(); iter != other_tracking_ids.end(); iter++ )
    trackingIds.insert( *iter ); 

  assert( acisHistory );
  
  //add the event to the history object
//  acisHistory->add_event( event ); 
}
  
void ATTRIB_HISTORY::trans_owner( SPAtransf const& )
{
}
void ATTRIB_HISTORY::to_tolerant_owner( ENTITY *tol_ent )
{
}

void ATTRIB_HISTORY::copy_owner( ENTITY *copy_ent )
{

  ATTRIB_HISTORY *other_history_attrib = get_history_attrib( copy_ent, true, acisHistory );
  other_history_attrib->trackingIds = this->trackingIds;
}

void ATTRIB_HISTORY::replace_owner( ENTITY *other_entity, logical replace_owner )
{
  if( replace_owner )
    return;
/*
  PRINT_INFO("I'm in replace owner!!!!\n");
  if( is_BODY( this->entity() )) 
    PRINT_INFO("  It's a body\n");
  if( is_LUMP( this->entity() )) 
    PRINT_INFO("  It's a volume\n");
  if( is_FACE( this->entity() )) 
    PRINT_INFO("  It's a surface\n");
  if( is_EDGE( this->entity() )) 
    PRINT_INFO("  It's an edge\n");
  if( is_VERTEX( this->entity() )) 
    PRINT_INFO("  It's a vertex\n");
*/
  //if this has alrady been merged w/ something else, it will have 
  //2 tracking ids in it's list.  I don't believe there's a need to 
  //track it then if it's already been merged
  
  if( trackingIds.size() > 1 )
    return;
  //Create an AcisEvent of type geometry_changed 
  AcisHistory::AcisEvent event;
  event.eventType = CGMHistory::GEOMETRY_CHANGED;

  event.entities.push_back( *(trackingIds.begin()) ); 

  assert( acisHistory );
  
  //add the event to the history object
  acisHistory->add_event( event ); 
}

void ATTRIB_HISTORY::lop_change_owner()
{
  /*
  PRINT_INFO("I'm in lop_change_owner owner!!!!\n");
  if( is_BODY( this->entity() )) 
    PRINT_INFO("  It's a body\n");
  if( is_LUMP( this->entity() )) 
    PRINT_INFO("  It's a volume\n");
  if( is_FACE( this->entity() )) 
    PRINT_INFO("  It's a surface\n");
  if( is_EDGE( this->entity() )) 
    PRINT_INFO("  It's an edge\n");
  if( is_VERTEX( this->entity() )) 
    PRINT_INFO("  It's a vertex\n");
*/
  //if this has alrady been merged w/ something else, it will have 
  //2 tracking ids in it's list.  I don't believe there's a need to 
  //track it then if it's already been merged
  
  if( trackingIds.size() > 1 )
    return;

  //Create an AcisEvent of type Subdivision
  AcisHistory::AcisEvent event;
  event.eventType = CGMHistory::GEOMETRY_CHANGED;

  event.entities.push_back( *(trackingIds.begin()) ); 

  assert( acisHistory );
  
  //add the event to the history object
  acisHistory->add_event( event ); 
}

void ATTRIB_HISTORY::replace_owner_geometry( ENTITY *new_geom )
{
  /*
  PRINT_INFO("I'm in replace_owner_geometry !!!!\n");
  if( is_BODY( this->entity() )) 
    PRINT_INFO("  It's a body\n");
  if( is_LUMP( this->entity() )) 
    PRINT_INFO("  It's a volume\n");
  if( is_FACE( this->entity() )) 
    PRINT_INFO("  It's a surface\n");
  if( is_EDGE( this->entity() )) 
    PRINT_INFO("  It's an edge\n");
  if( is_VERTEX( this->entity() )) 
    PRINT_INFO("  It's a vertex\n");
*/
  //if this has alrady been merged w/ something else, it will have 
  //2 tracking ids in it's list.  I don't believe there's a need to 
  //track it then if it's already been merged
  
  if( trackingIds.size() > 1 )
    return;
  //Create an AcisEvent of type Subdivision
  AcisHistory::AcisEvent event;
  event.eventType = CGMHistory::GEOMETRY_CHANGED;

  event.entities.push_back( *(trackingIds.begin()) ); 

  assert( acisHistory );
  
  //add the event to the history object
  acisHistory->add_event( event ); 
}

void ATTRIB_HISTORY::reverse_owner()
{
  PRINT_INFO("I'm in reverse_owner !!!!\n");
}

void ATTRIB_HISTORY::warp_owner( law *warp )
{
  PRINT_INFO("I'm in warp_owner !!!!\n");
}

ATTRIB_HISTORY* ATTRIB_HISTORY::get_history_attrib( ENTITY *acis_entity,
                                                    bool create_if_necessary,
                                                    AcisHistory *acis_history )
{

  if (acis_entity == NULL)
  {
    PRINT_ERROR("Trying to set a history attribute on a "
                " NULL ACIS ENTITY.\n");
    return NULL;
  }

    // Search for the attribute
  ATTRIB_HISTORY *history_attribute =
    (ATTRIB_HISTORY*)find_attrib(acis_entity,
                                 ATTRIB_SNL_TYPE, 
                                 ATTRIB_HISTORY_TYPE);
  
  //if we're not creating anything, return whatever we found
  if( false == create_if_necessary )
    return history_attribute;

    // If found, reset it's value
  if (history_attribute != NULL)
  {

  }
    // Otherwise, create a new attribute for acis_entity
  else
  {
    API_BEGIN;
    history_attribute = new ATTRIB_HISTORY( acis_entity, acis_history );
    API_END;
  }

  return history_attribute;
}


void ATTRIB_HISTORY::add_tracking_id( int id ) 
{
  trackingIds.insert( id );  
}

void ATTRIB_HISTORY::remove_all_attribs()
{
  //prevents 
  //addOrRemoveFromList = false;

  std::set<ATTRIB_HISTORY*>::iterator iter = allHistoryAttribs.begin();
  for( ; iter != allHistoryAttribs.end(); iter++ )
  {
    ATTRIB_HISTORY *history_attrib = *iter;
    history_attrib->unhook();
    history_attrib->lose();
  }

  allHistoryAttribs.clear();
  //addOrRemoveFromList = true;
}
