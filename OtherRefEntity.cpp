//-------------------------------------------------------------------------
// Filename      : OtherRefEntity.cc
//
// Purpose       : Implementation of the OtherRefEntity class.
//                 
//		   
//
// Special Notes :
//
// Creator       : Eric Nielsen
//
// Creation Date : 03/13/98
//
// Owner         : Eric Nielsen
//-------------------------------------------------------------------------

#include "OtherRefEntity.hpp"
#include "OtherEntity.hpp"

OtherRefEntity::OtherRefEntity()
{
   // Default constructor
}

OtherRefEntity::~OtherRefEntity()
{
   // destructor
}

//-------------------------------------------------------------------------
// Purpose       : Sets the associated OtherEntity
//
// Special Notes :
//
// Creator       : Eric Nielsen
//
// Creation Date : 03/13/98
//-------------------------------------------------------------------------
CubitStatus OtherRefEntity::set_other_entity_ptr(
                           OtherEntity* other_entity_ptr ) 
{
   // Make sure that a CUBIT OWNER attribute is attached to the solid
   // model entity underlying the OE.  This attribute contains a 
   // (back)pointer to this VGI GroupingEntity.

   CubitStatus status = other_entity_ptr->set_owner_attribute( this );

   if (status == CUBIT_SUCCESS)
   {
      otherEntityPtr = other_entity_ptr ;
   }
   
   return status ;
}
