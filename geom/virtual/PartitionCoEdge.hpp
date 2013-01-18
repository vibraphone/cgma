//-------------------------------------------------------------------------
// Filename      : PartitionCoEdge.hpp
//
// Purpose       : 
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 04/10/02
//-------------------------------------------------------------------------

#ifndef PARTITION_CO_EDGE_HPP
#define PARTITION_CO_EDGE_HPP

#include "CoEdgeSM.hpp"
#include "PartitionEntity.hpp"

class PartitionLoop;
class PartitionCurve;
class PartitionSurface;
class PartitionPoint;

class PartitionCoEdge : public CoEdgeSM, public PartitionEntity
{

  friend class PartitionLoop;
  friend class PartitionCurve;

  public:
  
    PartitionCoEdge(PartitionSurface* surf, CubitSense sense);
    PartitionCoEdge(PartitionCoEdge* split_from);
    PartitionCoEdge(CoEdgeSM* coedge);
  
    ~PartitionCoEdge();
    
    PartitionCurve* get_curve() const
      { return myCurve; }
    PartitionLoop* get_loop() const
      { return myLoop; }
      
    PartitionCoEdge* next() const
      { return loopNext; }

    TopologyBridge* find_parent_loop() const;
    CoEdgeSM* real_coedge() const;
    
    CubitSense sense() 
      { return mySense; }
    
    PartitionPoint* start_point() const;
    PartitionPoint*   end_point() const;
      
    void reverse_sense();
    
    virtual CubitBox bounding_box() const;
    virtual CubitStatus save( CubitSimpleAttrib& ) 
      { assert(0); return CUBIT_FAILURE; }
    virtual void transform( const CubitTransformMatrix& ) {;}
    
    void append_simple_attribute_virt( const CubitSimpleAttrib& );
    void remove_simple_attribute_virt( const CubitSimpleAttrib& );
    void remove_all_simple_attribute_virt();
    CubitStatus get_simple_attribute( DLIList<CubitSimpleAttrib>& );
    CubitStatus get_simple_attribute( const CubitString& name,
                                      DLIList<CubitSimpleAttrib>& );
    
    
    void get_parents_virt( DLIList<TopologyBridge*>& );
    void get_children_virt( DLIList<TopologyBridge*>& );
    int layer() const { return sub_entity_set().get_owner_layer(); }
    GeometryQueryEngine* get_geometry_query_engine() const;
    
    void notify_split( FacetEntity*, FacetEntity* );
    
  private:
  
    PartitionCoEdge( const PartitionCoEdge& );
  
    CubitSense mySense;
    
    PartitionLoop* myLoop;
    PartitionCoEdge* loopPrev;
    PartitionCoEdge* loopNext;
    
    PartitionCurve* myCurve;
    PartitionCoEdge* curveNext;
    
};


#endif
