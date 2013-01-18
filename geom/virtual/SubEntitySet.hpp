#ifndef SUB_ENTITY_SET_HPP
#define SUB_ENTITY_SET_HPP

#include "TBOwnerSet.hpp"

class PartitionEntity;
class TopologyBridge;
class CubitSimpleAttrib;
class CubitVector;
class PartitionBody;

class SubEntitySet : public TBOwnerSet
{

  public:
  
    SubEntitySet( TopologyBridge* real_entity,
                  PartitionEntity* first_partition );
    virtual ~SubEntitySet();
    
    TopologyBridge* get_entity() const
      { return myEntity; }
    
    void add_partition( PartitionEntity* partition,
                        PartitionEntity* insert_after = 0 );
                      
    void add_lower_order( PartitionEntity* partition );
    void add_lower_order( PartitionEntity* entity, 
                          const CubitSimpleAttrib& attrib,
                          int dimension,
                          DLIList<CubitVector*>& points_from_attrib,
                          DLIList<int>& facets_from_attrib,
                          DLIList<int>& children_from_attrib,
                          DLIList<int>& facet_point_owners_from_attrib );
    
    void remove( PartitionEntity* partition );
    
      
    bool has_lower_order() const
      { return lowerOrderHead != 0; }
    
    bool has_multiple_sub_entities() const;
    
    PartitionBody* body() const
      { return bodyPtr; }
    
    void get_owners( DLIList<TopologyBridge*>& owner_list ) const;
    
    CubitStatus remove_bridge( TopologyBridge* bridge );
    CubitStatus swap_bridge( TopologyBridge* remove, 
                             TopologyBridge* add,
                             bool reversed );
    void notify_reversed( TopologyBridge* bridge );
    CubitStatus bridge_destroyed( TopologyBridge* bridge );
  
    void print_debug_info( const char* prefix = 0 ) const;
  
    void add_attribute( PartitionEntity* entity, const CubitSimpleAttrib& csa );
    void rem_attribute( PartitionEntity* entity, const CubitSimpleAttrib& csa );
    void get_attributes( PartitionEntity* entity, 
                         DLIList<CubitSimpleAttrib>& list );
    void get_attributes( PartitionEntity* entity, const char* name,
                         DLIList<CubitSimpleAttrib>& list );
    void rem_all_attrib( PartitionEntity* entity );
    void unwrap_attributes();
  
    void get_sub_entities( DLIList<PartitionEntity*>& result_set ) const;
    void get_lower_order( DLIList<PartitionEntity*>& result_set ) const;

    int get_id( PartitionEntity* entity ) const;
    PartitionEntity* entity_from_id( int id ) const;
    
      // The following two methods are used by PartitionEngine
      // when restoring geometry.  Calling them at the wrong
      // time will result in the loss of attribute data.
    void set_id( PartitionEntity* entity, int id );
    void renumerate( int lowest_value, bool only_higher_ids  );
      // renumerate sets new SubEntitySet ids on owned entities.
      // The lowest id of any updated entity will be the passed value.
      // If only_higher_ids is true, only entities with IDs equal to
      // or higher than the passed value are updated.
      
    int get_owner_layer() const { return layerNumber; }
    
    
    CubitStatus save_geometry();
    
    int get_unique_id();
    void reset_unique_id() { uniqueId = 0; }
    
    CubitStatus save_geometry( int id, int dimension,
                               DLIList<CubitVector*>* point_list,
                               DLIList<int>* point_connectivity,
                               DLIList<int>* topo_connectivity,
                               DLIList<int>* point_owners,
                               CubitSimpleAttrib& attrib );

    static CubitStatus read_geometry( int& id, int& dimension,
                               DLIList<CubitVector*>& point_list,
                               DLIList<int>& point_connectivity,
                               DLIList<int>& topo_connectivity,
                               DLIList<int>& point_owners,
                               const CubitSimpleAttrib& attrib );
    
    static int get_geom_dimension( const CubitSimpleAttrib& attrib );
      // return the dimension of the geometry stored in the attrib
    
    static int get_geom_id( const CubitSimpleAttrib& attrib );
    
    static int get_segment_count( const CubitSimpleAttrib& attrib );
    
    static void remove_non_geom_attribs( DLIList<CubitSimpleAttrib>& list );
      // given a list of all the geometry attributes on an entity,
      // remove and destroy any that are not partition geometry.
    
    static void strip_attributes( TopologyBridge* bridge );
      // remove any attributes related to partition geometry
      // from the passed bridge.

    void strip_attributes(); 
      // remove all attributes on paritioned entities,
      // sub entities and lower-order entities
      
    inline SubEntitySet* next_in_body() const
      { return bodyNext; }
    
  private:
    
    friend class PartitionBody;
    SubEntitySet* bodyNext;
    PartitionBody* bodyPtr;
      
       // don't allow assignment (make assignment private)
    SubEntitySet& operator=( const SubEntitySet& ) { return *this;}
    bool operator==( const SubEntitySet& c ) const { return &c == this; }
  
    CubitStatus wrap_attribute( CubitSimpleAttrib& csa, int id ) const;
    int unwrap_attribute( CubitSimpleAttrib& csa ) const;
    bool is_attribute( const CubitSimpleAttrib& csa, int id = 0 ) const;
      
    TopologyBridge* myEntity;
    
    PartitionEntity* subEntityHead;
    PartitionEntity* lowerOrderHead;
    
    int lastId;
    
    int layerNumber;
    
    int uniqueId;
};

#endif
