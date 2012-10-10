//-------------------------------------------------------------------------
// Copyright Notice
//
// Copyright (c) 1996 
// by Malcolm J. Panthaki, DBA, and the University of New Mexico.
//-------------------------------------------------------------------------
//
//-------------------------------------------------------------------------
// Filename      : AssemblySM.cc
//
// Purpose       : To implement the AssemblySM solid model class.
//
// Special Notes :
//
// Creator       : Madhan Narayanan
//
// Creation Date : 07/18/97
//
// Owner         : Madhan Narayanan
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

#include "AssemblySM.hpp"
#include "RefAssembly.hpp"
// ********** END CUBIT INCLUDES           **********

// ********** BEGIN STATIC DECLARATIONS    **********
// ********** END STATIC DECLARATIONS      **********

// ********** BEGIN PUBLIC FUNCTIONS       **********

//-------------------------------------------------------------------------
// Purpose       : The default constructor. Does not do anything right now.
//
// Special Notes :
//
// Creator       : Stephen J. Verzi
//
// Creation Date : 02/26/97
//-------------------------------------------------------------------------
AssemblySM::AssemblySM()
{
}

//-------------------------------------------------------------------------
// Purpose       : The destructor. Does not do anything right now.
//
// Special Notes :
//
// Creator       : Stephen J. Verzi
//
// Creation Date : 02/26/97
//-------------------------------------------------------------------------
AssemblySM::~AssemblySM()
{
}

//-------------------------------------------------------------------------
// Purpose       : Get type of TopologyEntity this GeometryEntity
//                 should be attached to.
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 06/14/01
//-------------------------------------------------------------------------
const type_info& AssemblySM::topology_entity_type_info() const
{ return typeid(RefAssembly); }

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

