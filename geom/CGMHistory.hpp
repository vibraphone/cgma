

#ifndef CGMHistory_hpp
#define CGMHistory_hpp

#include <vector>
#include "CubitTransformMatrix.hpp"
#include "CubitGeomConfigure.h"
class RefEntity;
class TopologyBridge;

class CUBIT_GEOM_EXPORT CGMHistory
{
  public:
    CGMHistory();
    ~CGMHistory();
            
    enum EventType
    {
      // Top level creation means an entity with no parent is created
      // this could be bodies, surfaces, curves or vertices
      TOP_LEVEL_ENTITY_CREATED,
      
      // Top level deleted means an entity no longer has top level status
      // it may have been deleted or it may have been consumed into a higher
      // dimensional entity
      TOP_LEVEL_ENTITY_DELETED,
      
      // An entity has been created
      ENTITY_CREATED,
      
      // An entity has been deleted
      ENTITY_DELETED,
      
      // An entity's topology changes, meaning it has gained new child entities
      // or lost old child entities.  This can result from solid modeling
      // operations or from operations within cgm such as merging.
      TOPOLOGY_CHANGED, 

      // The geometry of an entity is changed such that it might have been
      // stretched, warped or other such changes where the topology may be the
      // same.  If child entities changed, they'll have their own event.
      GEOMETRY_CHANGED,

      // A transformation happend to this entity, all child entities are
      // transformed as well.  In this case, a GEOMETRY_CHANGED does not occur.
      GEOMETRY_TRANSFORMED,

      // The id of an entity changes
      ENTITY_ID_CHANGED,
      
      // The name of an entity changes
      ENTITY_NAME_CHANGED,
      
      // The name of an entity changes
      ENTITY_COLOR_CHANGED,

      // An entity is subdivided into multiple entities
      // this is a supplemental event such that other events completely
      // describe the topology changes, but this event specifies where the
      // subdivision occurred
      SUBDIVISION,

      // Multiple entities are absorbed into one entity
      // this is a supplemental event such that other events completely
      // describe the topology changes, but this event specifies where the
      // absorption occurred
      ABSORPTION,

      // Multiple entities modified by a combination of subdivision and
      // absorption, and intermediate entities don't exist by the time
      // this history object is given to someone
      // this is a supplemental event such that other events completely
      // describe the topology changes, but this event specifies where the
      // subdivision/absorption occurred
      SUBDIVISION_ABSORPTION,

      // An entity is merged into another
      // this is a supplemental event such that other events completely
      // describe the topology changes, but this event specifies where the
      // merge occurred
      MERGE,

      // An entity is copied from another
      // you may also get this when an unmerge happens in cgm
      // this is a supplemental event such that other events completely
      // describe the topology changes, but this event specifies where the
      // copy occurred
      COPY,
      
      // Cousin to SUBDIVISION.  Think of it as a 1-to-n subdivision where
      // only 1 of the n actually survives, so it's really and 1-to-1
      // modification.
      CUT
    };

    //Event is a record of a change that occurred during an operation
    class CUBIT_GEOM_EXPORT Event
    {
    public:
      Event(EventType type, RefEntity* refentity);
      Event(EventType type, const std::vector<RefEntity*>& refentity_list);
      ~Event();

      EventType get_event_type() const;
      const std::vector<RefEntity*>& get_entities() const;

    private:

      // the type of event this is
      EventType eventType;
      // the entities this event applies to
      std::vector<RefEntity*> entities;

      // extra data associated with event
      union 
      {
        std::vector<RefEntity*> *other_entities;
        CubitTransformMatrix *matrix;
        // TODO add data types for other events
      };
    };

    
    class CUBIT_GEOM_EXPORT PortEvent
    {
      public:
        PortEvent(EventType type, std::vector<RefEntity*> &source_entities,
                              std::vector<TopologyBridge*> &result_entities );
        PortEvent( EventType type );
        ~PortEvent();

        EventType get_event_type() const;
        const std::vector<RefEntity*>& get_entities() const;

        EventType eventType;
        // the entities this event applies to
        std::vector<RefEntity*> RefEnts; 
        std::vector<TopologyBridge*> TopologyBridges;

        // extra data associated with event
        //union 
       // {
        //  std::vector<TopologyBridge*> *TopologyBridges;
        //  CubitTransformMatrix *matrix;
          // TODO add data types for other events
       // };
    };



    // get the number of events in this history
    int get_number_of_events() const;
    // get an event by index
    const Event* get_event( int index ) const;

    void print_port_events(); 

    // add an event to this history
    void add_event( const Event &new_event );
    // compress the events in this history
    // for example, if an event for volume 1 created and an event for volume 1
    // deleted exists, both are removed

    void add_port_event( const PortEvent &event );

    void compress();
    
    // cleans out the history
    void clear();

    // start tracking events
    void start_tracking();
    // stop tracking events
    void end_tracking();
    // ask if tracking
    bool is_tracking() const;

  private:
    // the list of events
    std::vector<Event> eventList;
    std::vector<PortEvent> portEventList;
    bool mTracking;
};

#endif

