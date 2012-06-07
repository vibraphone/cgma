
#include "CGMHistory.hpp"
#include "TopologyBridge.hpp"
#include "CubitMessage.hpp"
#include "RefEntity.hpp"

const char* event_names[] = 
{
  "TOP_LEVEL_ENTITY_CREATED",
  "TOP_LEVEL_ENTITY_DELETED",
  "ENTITY_CREATED",
  "ENTITY_DELETED",
  "TOPOLOGY_CHANGED",
  "GEOMETRY_CHANGED",
  "GEOMETRY_TRANSFORMED",
  "ENTITY_ID_CHANGED",
  "ENTITY_NAME_CHANGED",
  "ENTITY_COLOR_CHANGED",
  "SUBDIVISION",
  "ABSORPTION",
  "SUBDIVISION_ABSORPTION",
  "MERGE",
  "COPY",
  "CUT"
};

CGMHistory::CGMHistory()
{
  mTracking = false;
  //mTracking = true;
}

CGMHistory::~CGMHistory()
{
}

CGMHistory::Event::Event(CGMHistory::EventType type, RefEntity* refentity)
{
  eventType = type;
  entities.push_back(refentity);
}

CGMHistory::Event::Event(CGMHistory::EventType type, const std::vector<RefEntity*>& refentity_list)
{
  eventType = type;
  entities = refentity_list;
}

CGMHistory::Event::~Event()
{
}
      
CGMHistory::EventType CGMHistory::Event::get_event_type() const
{
  return eventType;
}
      
const std::vector<RefEntity*>& CGMHistory::Event::get_entities() const
{
  return entities;
}

void CGMHistory::start_tracking()
{
  mTracking = true;
}

void CGMHistory::end_tracking()
{
  mTracking = false;
}

bool CGMHistory::is_tracking() const
{
  return mTracking;
}

int CGMHistory::get_number_of_events() const
{
  return eventList.size();
}
    
const CGMHistory::Event* CGMHistory::get_event( int index ) const
{
  if(index < eventList.size())
    return &eventList[index];
  return NULL;
}

void CGMHistory::add_event( const CGMHistory::Event &event )
{
  if(mTracking)
    eventList.push_back(event);
}

void CGMHistory::add_port_event( const CGMHistory::PortEvent &event )
{
  if(mTracking)
    portEventList.push_back(event);
}

void CGMHistory::compress()
{
  // TODO
}

void CGMHistory::clear()
{
  eventList.clear();
  portEventList.clear();
}

void CGMHistory::print_port_events() 
{
  
  std::vector<PortEvent>::iterator iter = portEventList.begin();
  for(; iter != portEventList.end(); iter++ )
  {
    PortEvent port_event = *iter;
    PRINT_INFO("Event type = %s RefEntities ", event_names[ port_event.eventType ]);

    int i;
    for( i=0; i<port_event.RefEnts.size(); i++ )
    {
      RefEntity *ref_ent = port_event.RefEnts[i];
      PRINT_INFO(" %s %d ", ref_ent->class_name(), ref_ent->id() );
    }
    PRINT_INFO(" ---> ");
    for( i=0; i<port_event.TopologyBridges.size(); i++ )
    {
      TopologyBridge *tb = port_event.TopologyBridges[i];
      PRINT_INFO(" %p", static_cast<void*>(tb) );
    }
    PRINT_INFO("\n"); 
  }
}

CGMHistory::PortEvent::PortEvent(CGMHistory::EventType type )
{
  eventType = type;
//  if( eventType != GEOMETRY_TRANSFORMED )
//    TopologyBridges = new std::vector<TopologyBridge*>;
//  else
//    TopologyBridges = NULL;
}

CGMHistory::PortEvent::~PortEvent()
{
//  if( TopologyBridges )
//    delete TopologyBridges;
}

CGMHistory::EventType CGMHistory::PortEvent::get_event_type() const
{
  return eventType;
}

