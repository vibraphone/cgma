//-------------------------------------------------------------------------
// Filename      : PartitionLump.hpp
//
// Purpose       : 
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 08/15/02
//-------------------------------------------------------------------------
#ifndef PARTITION_LUMP_HPP
#define PARTITION_LUMP_HPP

#include "DLIList.hpp"
#include "Lump.hpp"
#include "PartitionEntity.hpp"
#include "PartitionShell.hpp"

class PartitionBody;
class PartitionSurface;
class PartitionCurve;
class PartitionPoint;

class PartitionLump : public Lump, public PartitionEntity
{

public:
  friend class PartitionBody;
  
  PartitionLump( Lump* real_lump );
  PartitionLump( PartitionLump* split_from );
  virtual ~PartitionLump();
  
  Lump* real_lump() const;
  
  PartitionShell* first_shell() const;
  PartitionShell* next_shell( PartitionShell* after_this = 0 ) const;
  
  CubitStatus add( PartitionShell* shell );
  CubitStatus remove( PartitionShell* shell );
  void remove_all_shells();
  
  PartitionBody* get_body() const;
  TopologyBridge* find_parent_body() const;
  
  virtual CubitStatus save( CubitSimpleAttrib& );
  
  void get_parents_virt( DLIList<TopologyBridge*>& );
  void get_children_virt( DLIList<TopologyBridge*>& );
  int layer() const { return sub_entity_set().get_owner_layer(); }
  GeometryQueryEngine* get_geometry_query_engine() const;
  
  void append_simple_attribute_virt( CubitSimpleAttrib* );
  void remove_simple_attribute_virt( CubitSimpleAttrib* );
  void remove_all_simple_attribute_virt();
  CubitStatus get_simple_attribute( DLIList<CubitSimpleAttrib*>& );
  CubitStatus get_simple_attribute( const CubitString& name,
                                    DLIList<CubitSimpleAttrib*>& );
  
  CubitBox bounding_box() const;
  double measure();
  
  void reverse_sense();
  void transform(const CubitTransformMatrix&);
  
  //void print_debug_info( const char* prefix = 0 ) const;
  
  void notify_split( FacetEntity*, FacetEntity* );
  
  void get_all_children( DLIList<PartitionEntity*>& result );
  
  CubitStatus mass_properties( CubitVector& volume_centrioid,
                               double& volume );
  
  CubitPointContainment point_containment( const CubitVector& pos );
 
private:

  PartitionShell* listHead;
};

inline PartitionShell* PartitionLump::first_shell() const
  { return listHead; }

inline PartitionShell* PartitionLump::next_shell( PartitionShell* prev ) const
  { return !prev ? listHead : (prev->myLump == this) ? prev->lumpNext : 0; }

#endif



   
