#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "RefPart.hpp"
#include "RefEntityName.hpp"
#include "DLIList.hpp"
#include "CollectionEntity.hpp"
#include "CastTo.hpp"
#include "RefEntityFactory.hpp"
#include "PartSM.hpp"

RefPart::RefPart(CollectionEntity* collection_entity_ptr)
#ifdef PROE
 : RefGroup( 2 )
#endif
{
  recursionMark = 0;
  entityId = RefEntityFactory::instance()->next_ref_part_id();

  // assign default name
  //assign_entity_name(); // doesnt do anything now.
  
  // get the name of the part from the ce and add the
  // part id to it.
  
  CubitString sm_name;
  CubitString part_name;
  if( collection_entity_ptr != NULL ){
         
     // Get the associated name, if it exists
     sm_name = collection_entity_ptr->name();    
     //     collection_entity_ptr->remove_name_attribute( sm_name );
     
     part_name = sm_name + CubitString( entityId );
          
     RefEntityName::instance()->add_refentity_name(this, part_name);
  }
}

RefPart::~RefPart ()
{
     // remove this entity from the global list
   RefEntityFactory::instance()->remove(this);
   remove_all_ref_entities();
}


CubitStatus RefPart::remove()
{   
   // delete the PartXXX object.
   delete collectionEntityPtr;
   
   return CUBIT_SUCCESS;
}

GeometryQueryEngine* RefPart::get_geometry_query_engine() const
{
   return collectionEntityPtr->get_geometry_query_engine();   
}

PartSM* RefPart::get_part_sm_ptr() const
{
  return dynamic_cast<PartSM*>(collectionEntityPtr);
}
