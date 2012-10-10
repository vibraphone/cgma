//-------------------------------------------------------------------------
// Filename      : OtherEntity.cc
//
// Purpose       : 
//
// Special Notes :
//
// Creator       : Eric Nielsen
//
// Creation Date : 03/12/98
//
// Owner         : Eric Nielsen
//-------------------------------------------------------------------------

// ********** BEGIN STANDARD INCLUDES      **********
// ********** END STANDARD INCLUDES        **********

// ********** BEGIN MOTIF INCLUDES         **********
// ********** END MOTIF INCLUDES           **********

// ********** BEGIN OPEN INVENTOR INCLUDES **********
// ********** END OPEN INVENTOR INCLUDES   **********

// ********** BEGIN ACIS INCLUDES          **********
// ********** END ACIS INCLUDES            **********

// ********** BEGIN CUBIT INCLUDES         **********

#include "CubitMessage.hpp"
#include "OtherEntity.hpp"

#include "GeometryQueryEngine.hpp"

// ********** END CUBIT INCLUDES           **********

// ********** BEGIN STATIC DECLARATIONS    **********
// ********** END STATIC DECLARATIONS      **********

// ********** BEGIN PUBLIC FUNCTIONS       **********

//-------------------------------------------------------------------------
// Purpose       : Default constructor. Does not do anything right now.
//
// Special Notes :
//
// Creator       : Eric Nielsen
//
// Creation Date : 03/12/98
//-------------------------------------------------------------------------

OtherEntity::OtherEntity()
{
}

//-------------------------------------------------------------------------
// Purpose       : Destructor. Does not do anything right now.
//
// Special Notes :
//
// Creator       : Eric Nielsen
//
// Creation Date : 03/12/98
//-------------------------------------------------------------------------

OtherEntity::~OtherEntity()
{
}


int OtherEntity::validate(const CubitString&)
{
  // No way to validate a generic OtherEntity.
  // Override in subclasses.
  return 0;
}

// ********** END PUBLIC FUNCTIONS         **********

// ********** BEGIN PROTECTED FUNCTIONS    **********

//-------------------------------------------------------------------------
// Purpose       : This function merges the owner attributes of two 
//                 OtherEntity. At the end of this function, the owner 
//                 of the deadOE will be the same as owner of this O. 
//                 The owner attribute of deadOE is changed to the owner 
//                 attribute of this OE.
//
// Special Notes :
//
// Creator       : Eric Nielsen
//
// Creation Date : 03/16/98
//-------------------------------------------------------------------------

CubitStatus OtherEntity::merge_owner_attribute(OtherEntity* deadOE)
{
   // set the deadOE's owner the OtherRefEntity which
   // is the owner of this OE>

   OtherRefEntity* ownerORE = this->get_owner_attribute();
   CubitStatus result = deadOE->set_owner_attribute(ownerORE) ;

   if( result == CUBIT_FAILURE )
   {
      PRINT_ERROR("ERROR: In OtherEntity::merge_owner_attribute()\n") ;
      PRINT_ERROR("       Cannot set the owner attribute of OtherEntity.\n");
      PRINT_ERROR("  THIS IS A BUG - PLEASE REPORT IT!\n");
      assert(0);
      return CUBIT_FAILURE;
   }

   return result ;
}
   
// ********** END PROTECTED FUNCTIONS      **********

// ********** BEGIN PRIVATE FUNCTIONS      **********
// ********** END PRIVATE FUNCTIONS        **********

// ********** BEGIN HELPER CLASSES         **********
// ********** END HELPER CLASSES           **********

// ********** BEGIN EXTERN FUNCTIONS       **********
// ********** END EXTERN FUNCTIONS         **********

// ********** BEGIN STATIC FUNCTIONS       **********
// ********** END STATIC FUNCTIONS         **********

