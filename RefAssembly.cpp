#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "RefAssembly.hpp"
#include "RefEntityName.hpp"
#include "DLIList.hpp"
#include "CollectionEntity.hpp"
#include "CastTo.hpp"
#include "RefEntityFactory.hpp"
#include "AssemblySM.hpp"

// The next two includes are used in the match_intervals function
// If that were split out into a MeshingController class, we could
// remove these includes.
#include "RefFace.hpp"

RefAssembly::RefAssembly(CollectionEntity* collection_entity_ptr)
#ifdef PROE
 : RefGroup( 1 )
#endif  
{
  recursionMark = 0;
  entityId = RefEntityFactory::instance()->next_ref_assembly_id();

    // assign default name - dont assign default name get it from proe
    //assign_entity_name();
  if( collection_entity_ptr != NULL )
  {
    CubitString assembly_name;
    
      // Get the associated name, if it exists
    assembly_name = collection_entity_ptr->name() + CubitString(entityId);
    RefEntityName::instance()->add_refentity_name(this, assembly_name);
  }
}

RefAssembly::~RefAssembly ()
{
     // remove this entity from the global list
   RefEntityFactory::instance()->remove(this);
   remove_all_ref_entities();
}


CubitStatus RefAssembly::remove()
{
 
   delete collectionEntityPtr;
   
   return CUBIT_SUCCESS;
}

GeometryQueryEngine* RefAssembly::get_geometry_query_engine() const
{
   return collectionEntityPtr->get_geometry_query_engine();   
}

AssemblySM* RefAssembly::get_assembly_sm_ptr() const
{
  return dynamic_cast<AssemblySM*>(collectionEntityPtr);
}
  
