
#include "OCCHistory.hpp"
#include "OCCQueryEngine.hpp"
#include "CubitMessage.hpp"
#include "RefEntity.hpp"
#include "occ_attrib_history.hpp"
#include "CGMHistory.hpp"
#include "GeometryQueryTool.hpp"

/*Here's a brief description of how the history is created in the 
 * OCC port:
 * 1.  Create an OCCHistory object.
 *
 * 2. TopologyBridges are given to the port for modification.  From 
 * the TopologyBridges, get all associated parent and child TopologyBridges.  
 * From each TopologyBridge, get a) the corresponding RefEntity
 * and b) the corresponding TopoDS_Shape(s).  
 *
 * 3.  On every OCC TopoDS_Shape, attach an attrib_history attribute. 
 * The attrib_history basically holds a unique 'tracking' integer.
 *
 * 4.  Populate the refentity_tracking_id_map in the OCCHistory object, with
 * pairs, where each pair consists of a RefEntity pointer and it's associated
 * tracking id, in #3.
 *
 * 5.  OCC performs the operation(s), making callbacks into the 
 * occ_attrib_history's overloaded functions.  It does the following:
 * -For split_owner callbacks, where A is split into A and B, 
 *  put a history attrib on B and add A's tracking integer to B.
 * -For merge_owner callbacks, where A and B are merged into A, 
 *  add B's tracking id to A.  Now A will have 2 tracking ids 
 *  associated with it.
 *  -In the copying callback macros, if B has been copied from A, give B
 *  the tracking id that was in A.
 *  -For the other callbacks, replace_owner, lop_change_owner, and 
 *  replace_owner_geometry, add an OCCEvent to the OCCHistory object's
 *  event list.
 *
 * 6.  Once the operation has ended and TopologyBridges have been created
 * from the new TopoDS_Shape's (populate_topology_bridges has been called), for 
 * each TopologyBridge, create a pair consisting of
 * a set of tracking ids and a TopologyBridge pointer. 
 *
 * 7.  Find pairs in #6 that have any tracking id in common in their sets. 
 * If you find multiple pairs who's sets in size, you have a subdivision-absorption
 * event (n->m).  If you only find a single pair that has a set with 
 * multiple tracking ids in it, you have an absorption event (n->1). If you find
 * multiple pairs that all have a single tracking id in their set, you have
 * a subdivision event(1->n).  And if you only find a single pair with a single
 * tracking id in it's set, you have a geometry modified event.
 *
 * 8. Everytime you gather one of these pair group, in #7, you create a PortEvent 
 * which ties the original RefEntity pointer(s) to the new TopologyBridge(s).
 */

OCCHistory::OCCHistory()
{
  trackingId = -1;
  highLevelCopying = false;
}

OCCHistory::~OCCHistory()
{
 //OCC_ATTRIB_HISTORY::remove_all_attribs(); 
}

void OCCHistory::add_event( const OCCEvent &event )
{
  eventList.push_back( event );
}

void OCCHistory::add_refentity_tracking_id_pair( RefEntity *ref_ent, int tracking_id ) 
{
  refentity_tracking_id_map.push_back( 
    std::pair<RefEntity*, int>( ref_ent, tracking_id));
}


void OCCHistory::print()
{
  /*
  std::vector<CGMEvent>::iterator iter = cgmEventList.begin();
  for(; iter != cgmEventList.end(); iter++ )
  {
    CGMEvent cgm_event = *iter;
    PRINT_INFO("Event type = %s RefEntities ", event_names[ cgm_event.eventType ]);

    int i;
    for( i=0; i<cgm_event.RefEnts.size(); i++ )
    {
      RefEntity *ref_ent = cgm_event.RefEnts[i];
      PRINT_INFO(" %s %d ", ref_ent->class_name(), ref_ent->id() );
    }
    PRINT_INFO(" ---> ");
    for( i=0; i<cgm_event.TopologyBridges.size(); i++ )
    {
      TopologyBridge *tb = cgm_event.TopologyBridges[i];
      PRINT_INFO(" %p", tb );
    }
    PRINT_INFO("\n"); 
  } */
}

int OCCHistory::generate_next_tracking_id()
{
  trackingId++;
  return trackingId;
}

void OCCHistory::add_to_tracking_ids_tb_map( std::set<int> &tracking_ids, TopologyBridge *tb ) 
{
  tracking_ids_to_tb_map.push_back( std::pair< std::set<int>, TopologyBridge*>( tracking_ids, tb ) );
}

void OCCHistory::create_cgm_history_objects()
{
  //look for subdivision-absorption events.
  //Look in tracking_ids_to_tb_map for a pair whose vector 
  //has multiple integers in it.  If you find another pair with the 
  //same integers in it, you've found a subdivision-absorption event.

  std::vector< std::pair< std::set<int>, TopologyBridge*> > unrelated_pairs;
  std::list< std::pair< std::set<int>, TopologyBridge*> >::iterator iter;
 
  if( highLevelCopying )
  {
    while( tracking_ids_to_tb_map.size() )
    {
      iter = tracking_ids_to_tb_map.begin();

      std::pair< std::set<int>, TopologyBridge*> pair1 = *iter; 

      std::vector< std::pair< std::set<int>, TopologyBridge*> > found_pairs; 
      found_pairs.push_back( pair1 );
      
      create_port_event( CGMHistory::COPY, found_pairs );
      tracking_ids_to_tb_map.erase( iter++ );
    }
  }
  else
  {

    while( tracking_ids_to_tb_map.size() )
    {
      iter = tracking_ids_to_tb_map.begin();

      std::pair< std::set<int>, TopologyBridge*> pair1 = *iter; 
      std::set<int> set1 = pair1.first;

      //find all pairs that have a vector that contains any integer in set1,
      //and at the same time recursively do the same for the found pairs 
      //get back a 'found_pairs' list
      
      std::vector< std::pair< std::set<int>, TopologyBridge*> > found_pairs; 
      found_pairs.push_back( pair1 );
      tracking_ids_to_tb_map.erase( iter++ );
      collect_relevant_pairs( set1, found_pairs ); 

      //if the found_pairs just contains the original pair
      //it's a simple absorption
      std::vector< std::pair< std::set<int>, TopologyBridge*> >::iterator vec_iter;
      std::pair< std::set<int>, TopologyBridge*> tmp_pair; 
      if( found_pairs.size() == 1 )
      {
        //if the single set contains multiple integers, 
        //it's a simple absorption

        vec_iter = found_pairs.begin();
        tmp_pair = *vec_iter;
        std::set<int> tmp_set = tmp_pair.first;
        if( tmp_set.size() > 1 )
          create_port_event( CGMHistory::ABSORPTION, found_pairs ); 
        else
          unrelated_pairs.push_back( pair1 );
      }
      else
      {
        //if the found_pairs all have exactly one integer in their sets, 
        //it's a simple subdivision
        vec_iter = found_pairs.begin();
        tmp_pair = *vec_iter;
        std::set<int> tmp_set = tmp_pair.first;
        int size_of_set = tmp_set.size(); 
        vec_iter++; 
        bool varying_sizes = false;
        for(; vec_iter != found_pairs.end(); vec_iter++ )
        {
          tmp_pair = *vec_iter;
          tmp_set = tmp_pair.first;
          if( size_of_set != tmp_set.size() )
          {
            varying_sizes = true;
            break;
          }
        }

        if( varying_sizes )
          create_port_event( CGMHistory::SUBDIVISION_ABSORPTION, found_pairs ); 
        else
          create_port_event( CGMHistory::SUBDIVISION, found_pairs ); 
      }
    }
  }

  std::vector< std::pair< std::set<int>, TopologyBridge*> >::iterator vec_iter1;
  vec_iter1 = unrelated_pairs.begin();
  for(; vec_iter1 != unrelated_pairs.end(); vec_iter1++ ) 
  {
    std::set<int> set1 = (*vec_iter1).first;
    if( set1.size() == 0 )
      continue;
    int tracking_id = *set1.begin();

    std::vector<OCCEvent>::iterator vec_iter2 = eventList.begin();
    for(; vec_iter2 != eventList.end(); vec_iter2++ )
    {
      OCCEvent event = (*vec_iter2);
      if( event.eventType == CGMHistory::GEOMETRY_CHANGED ||
          event.eventType == CGMHistory::CUT ) 
      {
        int other_tracking_id = event.entities[0];

        if( tracking_id == other_tracking_id )
        {
          TopologyBridge *tb = (*vec_iter1).second;
          std::pair<RefEntity*, int> tmp_pair = refentity_tracking_id_map[tracking_id];
          RefEntity* ref_ent = tmp_pair.first; 
          CGMHistory::PortEvent port_event( event.eventType );
          port_event.RefEnts.push_back( ref_ent );
          port_event.TopologyBridges.push_back( tb );
          GeometryQueryTool::instance()->history().add_port_event( port_event ); 
          break;
        }
      }
    }
  }
}


void OCCHistory::collect_relevant_pairs( std::set<int> &seed_set,
           std::vector< std::pair< std::set<int>, TopologyBridge*> > &found_pairs ) 
{
   
  std::list< std::pair< std::set<int>, TopologyBridge*> >::iterator iter;
  iter = tracking_ids_to_tb_map.begin();

  for(; iter != tracking_ids_to_tb_map.end(); iter++ ) 
  {
    std::pair< std::set<int>, TopologyBridge*> my_pair1 = *iter; 
    std::set<int> set1 = my_pair1.first;

    if( set1.empty() )
      continue;
    
    std::set<int>::iterator iter1, iter2;
    //does set1 contain any integers that are in the seed set?
    //if so, remove the pair from the list and add it to 'found_pairs'
    iter1 = set1.begin();
    bool pair_found = false;
    for(; iter1 != set1.end(); iter1++ )
    {
      int int1 = *iter1;
      
      for( iter2 = seed_set.begin(); iter2 != seed_set.end(); iter2++ )
      {
        int int2 = *iter2; 

        if( int2 == int1 )
        {
          //clear out the set
          std::set<int> &tmp_set = (*iter).first; 
          tmp_set.clear();

          found_pairs.push_back( my_pair1 );
          //tracking_ids_to_tb_map.erase( iter++ );
          collect_relevant_pairs( set1, found_pairs );
          pair_found = true;
          break;
        }
      }
      if( pair_found )
        break;
    }
//    if( pair_found == false )
//      iter++;
  }
}


void OCCHistory::create_port_event( CGMHistory::EventType event_type, 
           std::vector< std::pair< std::set<int>, TopologyBridge*> > &found_pairs ) 
{
  CGMHistory::PortEvent port_event( event_type );
  
  //add all TBs to the event's TB list
  std::vector< std::pair< std::set<int>, TopologyBridge*> >::iterator vec_iter;
  vec_iter = found_pairs.begin();
  std::pair< std::set<int>, TopologyBridge*> tmp_pair; 
  DLIList<RefEntity*> ref_ent_list;
  DLIList<TopologyBridge*> tb_list;
  for(; vec_iter != found_pairs.end(); vec_iter++ )
  {
    tmp_pair = *vec_iter;
    TopologyBridge *tb = tmp_pair.second; 
    tb_list.append( tb );

    std::set<int> tmp_set = tmp_pair.first;
    std::set<int>::iterator iter=tmp_set.begin();
    for(; iter != tmp_set.end(); iter++ )
    {
      int index = *iter;
      std::pair<RefEntity*, int> pair2 = refentity_tracking_id_map[index];
      RefEntity* ref_ent = pair2.first; 
      ref_ent_list.append( ref_ent );
    }
  }
    
  //uniquify the lists...then append them to the CGMEvent's vector
  ref_ent_list.uniquify_unordered(); 
  tb_list.uniquify_unordered(); 
  
  int i;
  for(i=ref_ent_list.size(); i--; )
    port_event.RefEnts.push_back( ref_ent_list.get_and_step() );
  for(i=tb_list.size(); i--; )
    port_event.TopologyBridges.push_back( tb_list.get_and_step() );

  GeometryQueryTool::instance()->history().add_port_event( port_event ); 
}


void OCCHistory::add_TopoDS_Shape( TopoDS_Shape *entity )
{
  all_ents.push_back( entity );
}

void OCCHistory::remove_attributes()
{
/*
  int i;
  for( i=0; i<all_ents.size(); i++ )
  {
    OCC_ATTRIB_HISTORY *att = OCC_ATTRIB_HISTORY::get_history_attrib( all_ents[i], false ); 
    if( att )
      att->lose();
  }
*/
}

/*
void OCCHistory::collect_events_and_find_resultant_uniqueids( int source_id, 
                                     DLIList<OCCEvent*> &events, 
                                     DLIList<int> &result_unique_ids )
{
  //get an OCCEvent w/ 'entity' in  'other_entities' list
  std::vector<OCCEvent>::iterator iter;
  bool found_entity = false;
  for(iter=eventList.begin(); iter != eventList.end(); iter++ )
  {
    OCCEvent *tmp_event = &(*iter);

    //ignore events that are already in the list
    if( events.is_in_list( tmp_event ) )
      continue;
 
    //search 'entities' list 
    if( tmp_event->eventType == SUBDIVISION || 
        tmp_event->eventType == ABSORPTION || 
        tmp_event->eventType == COPY ) 
    {
      int k;
      for(k=0; k<tmp_event->entities.size(); k++ )
      {
        if( source_id == tmp_event->entities[k] )
        {
          found_entity = true;
          //we found the ENTITY...now get ents in 'entities' list
          //and look for them 
          events.append( tmp_event );
          int i;
          for( i=0; i<tmp_event->other_entities.size(); i++ )
            collect_events_and_find_resultant_uniqueids( tmp_event->other_entities[i], 
                events, result_unique_ids );
        }
      }
    }
    else if( tmp_event->eventType == GEOMETRY_CHANGED ) 
    {
      if( source_id == tmp_event->entities[0] )
      {
        events.append( tmp_event );
        collect_events_and_find_resultant_uniqueids( tmp_event->entities[0], 
            events, result_unique_ids);
      }
    }
  }

  //so you can't go any further, traversing the events...see if this one is in
  //the results list...if so, append it
  if( false == found_entity ) 
  {
    if( resultUniqueIds.find( source_id ) != resultUniqueIds.end() ) 
      result_unique_ids.append( source_id );
  }
} */


/*
void OCCHistory::collect_events_and_find_original_uniqueids( int source_id, 
                                     DLIList<OCCEvent*> &events, 
                                     DLIList<int> &original_unique_ids)
{
 
  //get an OCCEvent w/ 'entity' in  'other_entities' list
  std::vector<OCCEvent>::iterator iter;
  bool found_entity = false;
  for(iter=eventList.begin(); iter != eventList.end(); iter++ )
  {
    OCCEvent *tmp_event = &(*iter);

    //ignore events that are already in the list
    if( events.is_in_list( tmp_event ) )
      continue;
 
    //search 'other_entities' list 
    if( tmp_event->eventType == SUBDIVISION || 
        tmp_event->eventType == ABSORPTION || 
        tmp_event->eventType == COPY ) 
    {
      int k;
      for(k=0; k<tmp_event->other_entities.size(); k++ )
      {
        if( source_id == tmp_event->other_entities[k] )
        {
          found_entity = true;
          //we found the ENTITY...now get ents in 'entities' list
          //and look for them 
          events.append( tmp_event );
          int i;
          for( i=0; i<tmp_event->entities.size(); i++ )
            collect_events_and_find_original_uniqueids( tmp_event->entities[i], events, original_unique_ids );
        }
      }
    }
    else if( tmp_event->eventType == GEOMETRY_CHANGED ) 
    {
      if( source_id == tmp_event->entities[0] )
      {
        events.append( tmp_event );
        collect_events_and_find_original_uniqueids( tmp_event->entities[0], events, original_unique_ids );
      }
    }
  }

  if( false == found_entity )
    original_unique_ids.append( source_id );
} */


/*
void OCCHistory::map_back_to_ref_entity( TopologyBridge *topology_bridge, 
                                          RefEntity *ref_entity, 
                                          OCCHistory::EventType event_type )
{

  ENTITY *entity = OCCQueryEngine::instance()->get_ENTITY_of_entity( topology_bridge ); 

  ATTRIB_HISTORY *tmp_attrib = ATTRIB_HISTORY::get_history_attrib( entity );
  if( !tmp_attrib )
    return;

  int source_id = tmp_attrib->get_unique_id();
  
  DLIList<OCCEvent*> events;
  DLIList<int> original_unique_ids;

  //starting from the result ENTITY, work your waying up the events, 
  // to the original ENTITY
  collect_events_and_find_original_uniqueids( source_id, events, original_unique_ids );
  
  //PRINT_INFO("source entity = %p  ", entity );

  //PRINT_INFO("target_ents = ");
  int i;
  //for( i=original_ents.size(); i--; )
  //  PRINT_INFO("%p ", original_ents.get_and_step() );
  //PRINT_INFO("\n");


  EventType my_event_type = NO_EVENT;
  DLIList<RefEntity*> original_ref_ents; 

  //from the supposed original ENTITY, 
  //find the corresponding RefEntity in refentity_ENTITY_multimap 
  for( i=original_unique_ids.size(); i--; )
  {
    int original_id = original_unique_ids.get_and_step();

    std::multimap<RefEntity*, int>::iterator iter = refentity_uniqueid_multimap.begin();
    for( ; iter != refentity_uniqueid_multimap.end(); iter++ )
    {
      int tmp_id = (*iter).second;
      if( original_id == tmp_id ) 
      {
        //determine what type of event you have
        int j;
        for( j=events.size(); j--; )
        {
          OCCEvent *event = events.get_and_step();

          //ignore copy events for now
          if( event->eventType != COPY && 
              event->eventType != NO_EVENT )
          {

            //split events take priority over 
            if( event->eventType > my_event_type  )
              my_event_type = event->eventType;
          }
        }
        RefEntity *original_ref_ent = (*iter).first;
        original_ref_ents.append( original_ref_ent );
      }
    }
  }

  if( my_event_type != NO_EVENT )
  {
    PRINT_INFO("Event type = %s RefEntity = %s %d Result ENTITY = %p\n", 
      event_names[my_event_type], original_ref_ents.get()->class_name(),  
      original_ref_ents.get()->id(), entity ); 

    CGMEvent tmp_event;
    tmp_event.eventType = my_event_type;
    for( i=original_ref_ents.size(); i--; )
      tmp_event.RefEnts.push_back( original_ref_ents.get_and_step() );

    tmp_event.TopologyBridges.push_back( topology_bridge ); 
    
    cgmEventList.push_back( tmp_event );
  }

}*/

/*Here's a brief description of how the history is created in the 
 * OCC port:
 * 1.  Create an OCCHistory object.
 *
 * 2.  Create a RefEnity to ENTITY multimap in the OCCHistory object:
 * TopologyBridges are given to the port for modification.  From 
 * the TopologyBridges, get all associated parent and child TopologyBridges.  
 * From each TopologyBridge, get a) the corresponding RefEntity and b) the
 * corresponding TopoDS_Shape(s).  Put them in 'refentity_ENTITY_multimap', 
 * where the RefEntity is the key and the TopoDS_Shape(s) is the data.
 *
 * 3.  Attach an occ_attrib_history attribute to every OCC TopoDS_Shape.
 *
 * 4.  OCC preforms the operation(s), making callbacks into the 
 * occ_attrib_history's overloaded functions.  Within these functions, OCCEvents
 * are created and added to the OCCHistory object's event list.
 *
 * 5.  Once the operation has ended and TopologyBridges have been created
 * from the new TopoDS_Shape, get all TopoDS_Shape off the resulting bodies and
 * add them to the OCCHistory's 'resultTopoDS_Shape's' list.  
 *
 * 6.  Remove (ignore) any SUBDIVISION events that don't result in multiple 
 * TopoDS_Shape's in the result body.
 *
 * 7.  For each TopologyBridge, trace your way back through the OCCEvents until 
 * you can't go any further.  With this 'leaf' entity, try to find a 
 * corresponding RefEntity in the 'refentity_ENTITY_multimap'.  
 *
 * 8.  If you found a corresponding RefEntity, create a CGMEvent.
 *
 * 9.  Consolidate CGMEvents of type SUBIDIVISION and ABSORPTION.
 *
 */



