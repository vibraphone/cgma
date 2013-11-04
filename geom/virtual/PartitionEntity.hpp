//-------------------------------------------------------------------------
// Filename      : PartitionEntity.hpp
//
// Purpose       : Base class for partition entities
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 04/21/02
//-------------------------------------------------------------------------

#ifndef PARTITION_ENTITY_HPP
#define PARTITION_ENTITY_HPP

#include "SubEntitySet.hpp"
#include "CubitBox.hpp"
class FacetEntity;
class CubitTransformMatrix;
class CubitPoint;

class PartitionEntity
{

  public:
  
    virtual ~PartitionEntity();
    
    TopologyBridge* partitioned_entity() const;
    
    SubEntitySet& sub_entity_set() const;
    
    virtual CubitStatus move_to_geometry( CubitVector& position );
    
    CubitStatus relax_to_geometry( CubitPoint* facet_point,
                                   const CubitVector* input_position = 0 );
    
    int mark;
    
    virtual void print_debug_info( const char* prefix = 0,
                                   bool print_subentity_set = true ) const;
    void print_partitioned_entity( const char* prefix = 0 ) const;
    
    virtual void transform(const CubitTransformMatrix& xform) = 0;
    
    virtual void reverse_sense() = 0;
    
    virtual void notify_split( FacetEntity* old_entity, 
                               FacetEntity* new_entity ) = 0;
    
    virtual CubitStatus save( CubitSimpleAttrib& attrib ) = 0;
    
    virtual CubitBox bounding_box() const = 0;  // for abstracttree
    
  protected:
  
    PartitionEntity( );
    
  private:
      
       // don't allow assignment (make assignment private)
    PartitionEntity( const PartitionEntity& );
    void operator=( const PartitionEntity& );
    
    friend class SubEntitySet;
    
    SubEntitySet* entitySet;
    PartitionEntity* entitySetNext;
      // These variables are set/cleared by SubEntitySet.
      // They are never modified directly by this class.
      
    int entitySetId;
      // An ID in the SubEntitySet used for save and restore.
};

inline TopologyBridge* PartitionEntity::partitioned_entity() const
{
  return entitySet->get_entity();
}

inline SubEntitySet& PartitionEntity::sub_entity_set() const
{
  return *entitySet;
}

#endif
