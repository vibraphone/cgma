//-------------------------------------------------------------------------
// Copyright Notice
//
// Copyright (c) 1996 
// by Malcolm J. Panthaki, DBA, and the University of New Mexico.
//-------------------------------------------------------------------------
//
//-------------------------------------------------------------------------
// Filename      : OtherSolidModelEntity.cc
//
// Purpose       : 
//
// Special Notes :
//
// Creator       : Xuechen Liu
//
// Creation Date : 08/06/96
//
// Owner         : Raikanta Sahu
//-------------------------------------------------------------------------

// ********** BEGIN STANDARD INCLUDES      **********
// ********** END STANDARD INCLUDES        **********

// ********** BEGIN CUBIT INCLUDES         **********
#include "OtherSolidModelEntity.hpp"
#include "CubitMessage.hpp"
// ********** END CUBIT INCLUDES           **********

// ********** BEGIN STATIC DECLARATIONS    **********
// ********** END STATIC DECLARATIONS      **********

// ********** BEGIN PUBLIC FUNCTIONS       **********

//-------------------------------------------------------------------------
// Purpose       : The default constructor.
//
// Special Notes :
//
// Creator       : Xuechen Liu
//
// Creation Date : 08/06/96
//-------------------------------------------------------------------------

OtherSolidModelEntity::OtherSolidModelEntity()
{
}

//-------------------------------------------------------------------------
// Purpose       : The destructor.
//
// Special Notes :
//
// Creator       : Raikanta Sahu
//
// Creation Date : 09/06/96
//-------------------------------------------------------------------------

OtherSolidModelEntity::~OtherSolidModelEntity()
{
}

OtherSolidModelEntity* OtherSolidModelEntity::copy()
{
  PRINT_ERROR("OtherSolidModelEntity::copy() called on an entity"
              " that doesn't implement it.\n");
  assert(CUBIT_TRUE);
  return (OtherSolidModelEntity *)NULL;
}

// ********** END PUBLIC FUNCTIONS         **********

// ********** BEGIN PROTECTED FUNCTIONS    **********
// ********** END PROTECTED FUNCTIONS      **********

// ********** BEGIN PRIVATE FUNCTIONS      **********
// ********** END PRIVATE FUNCTIONS        **********

// ********** BEGIN HELPER CLASSES         **********
// ********** END HELPER CLASSES           **********

// ********** BEGIN EXTERN FUNCTIONS       **********
// ********** END EXTERN FUNCTIONS         **********

// ********** BEGIN STATIC FUNCTIONS       **********
// ********** END STATIC FUNCTIONS         **********

