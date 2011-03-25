

#ifndef OCCHistory_hpp
#define OCCHistory_hpp

#include <map>
#include <vector>
#include <set>
#include <list>
#include "DLIList.hpp"
#include "TopologyBridge.hpp"
#include "CGMHistory.hpp"

class RefEntity;
class TopoDS_Shape;

class OCCHistory
{
  public:
    OCCHistory();
    ~OCCHistory();

    //Event is a record of a change that occurred during an operation
    struct OCCEvent
    {
      // the type of event this is
      CGMHistory::EventType eventType;

      // the entities this event applies to
      std::vector<int> entities;

      // extra data associated with event
      std::vector<int> other_entities;
    };

    struct CGMEvent
    {
      // the type of event this is
      CGMHistory::EventType eventType;

      // the entities this event applies to
      std::vector<RefEntity*> RefEnts;

      // extra data associated with event
      std::vector<TopologyBridge*> TopologyBridges;
    };

    // add an event to this history
    void add_event( const OCCEvent &event );

    // add a RefEntity, TopoDS_Shape pair
    void add_refentity_tracking_id_pair( RefEntity *ref_ent, int tracking_id ); 
    void add_to_tracking_ids_tb_map( std::set<int> &tracking_ids, TopologyBridge *tb ); 

    void collect_relevant_pairs( std::set<int> &seed_set,
           std::vector< std::pair< std::set<int>, TopologyBridge*> > &found_pairs );

    int generate_next_tracking_id();

    void print();

    void create_cgm_history_objects();

    void create_port_event( CGMHistory::EventType event_type, 
           std::vector< std::pair< std::set<int>, TopologyBridge*> > &found_pairs );

    void add_TopoDS_Shape( TopoDS_Shape *ent );
    void remove_attributes();

    bool highLevelCopying; 

  private:
    // the list of events
    std::vector<OCCEvent> eventList;
    std::vector< std::pair<RefEntity*, int> > refentity_tracking_id_map;
    std::list< std::pair< std::set<int>, TopologyBridge* > > tracking_ids_to_tb_map;
    std::vector<CGMEvent> cgmEventList;
    int trackingId;
    std::vector<TopoDS_Shape*> all_ents;
};

#endif

