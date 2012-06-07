//-------------------------------------------------------------------------
// Filename      : PartitionEntity.cpp
//
// Purpose       : 
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 04/21/02
//-------------------------------------------------------------------------

#include "VGDefines.h"
#include "PartitionEntity.hpp"
#include "CubitMessage.hpp"
#include "CubitPoint.hpp"

PartitionEntity::~PartitionEntity()
{
  entitySet->remove(this);
}

PartitionEntity::PartitionEntity( )
  : mark(0), entitySet(0), entitySetNext(0), entitySetId(0)
{
}

void PartitionEntity::print_debug_info( const char* prefix, 
                                        bool print_subent_set ) const
{
  if( !prefix ) prefix = "";
  PRINT_INFO("%sPartitionEntity %p\n", prefix, static_cast<void*>(this) );
  if( print_subent_set )
    sub_entity_set().print_debug_info(prefix);
  else
    print_partitioned_entity( prefix );
}

void PartitionEntity::print_partitioned_entity( const char* prefix ) const
{
  if(!prefix) prefix = "";
  PRINT_INFO("%sPartitioned Entity: %s %p\n",
    prefix, 
    partitioned_entity() ? 
    fix_type_name(typeid(*partitioned_entity()).name()) : "TopologyBridge",
    static_cast<void*>(partitioned_entity()) );
}

CubitStatus PartitionEntity::move_to_geometry( CubitVector& )
  { return CUBIT_SUCCESS; }

CubitStatus PartitionEntity::relax_to_geometry( CubitPoint* facet_point,
                                                const CubitVector* input_pos )
{
  if (input_pos)
  {
    CubitVector closest(*input_pos);
    if (move_to_geometry(closest) && 
        facet_point->check_inverted_facets(closest))
    {
      facet_point->set(closest);
      return CUBIT_SUCCESS;
    }
  }
  
  CubitVector facet_pos(facet_point->coordinates());
  if (move_to_geometry(facet_pos) &&
      facet_point->check_inverted_facets(facet_pos))
  {
    facet_point->set(facet_pos);
    return CUBIT_SUCCESS;
  }
  
  return CUBIT_FAILURE;
}
