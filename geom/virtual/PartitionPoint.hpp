#ifndef PARTITION_POINT_HPP
#define PARTITION_POINT_HPP

#include "Point.hpp"
#include "PartitionEntity.hpp"

class CubitPointData;
class PartitionCurve;

class PartitionPoint : public TBPoint, public PartitionEntity
{

friend class PartitionCurve;

public:

  PartitionPoint( const CubitVector& position, PartitionEntity* owner );
  PartitionPoint( const CubitSimpleAttrib& attrib, PartitionEntity* owner );
  PartitionPoint( TBPoint* real_point );
  
  ~PartitionPoint();
  
  virtual CubitVector coordinates() const;
  virtual CubitBox bounding_box() const;
  virtual CubitStatus move_to_geometry( CubitVector& );
 
  TBPoint* real_point() const;
  
  CubitStatus move( CubitVector& delta );
  
  void append_simple_attribute_virt( const CubitSimpleAttrib& );
  void remove_simple_attribute_virt( const CubitSimpleAttrib& );
  void remove_all_simple_attribute_virt();
  CubitStatus get_simple_attribute( DLIList<CubitSimpleAttrib>& );
  CubitStatus get_simple_attribute( const CubitString& name,
                                    DLIList<CubitSimpleAttrib>& );
  
  void get_parents_virt( DLIList<TopologyBridge*>& parents );
  void get_children_virt( DLIList<TopologyBridge*>& children );
  int layer() const { return sub_entity_set().get_owner_layer(); }
  GeometryQueryEngine* get_geometry_query_engine() const;
  
  int num_curves() const;
  PartitionCurve* next_curve( PartitionCurve* prev = 0 ) const;
  
  PartitionCurve* common_curve( PartitionPoint* other ) const;
  
  virtual void print_debug_info( const char* prefix = 0,
                                 bool print_subentity_set = true ) const;
  
  virtual void reverse_sense();
  
  virtual void notify_split( FacetEntity*, FacetEntity* );
  CubitPointData* facet_point() const { return facetPoint; }
  void facet_point( CubitPointData* set );
  
  virtual CubitStatus save( CubitSimpleAttrib& attrib );
  

  virtual void transform( const CubitTransformMatrix& );
private:

  PartitionCurve* firstCurve;
  int curveCount;

  CubitVector myPosition;
  CubitPointData* facetPoint;
};

inline int PartitionPoint::num_curves() const
  { return curveCount; }

#endif
