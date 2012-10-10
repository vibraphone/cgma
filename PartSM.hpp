//-------------------------------------------------------------------------
// Copyright Notice
//
// Copyright (c) 1996 
// by Malcolm J. Panthaki, DBA, and the University of New Mexico.
//-------------------------------------------------------------------------
//
//-------------------------------------------------------------------------
// Filename      : PartSM.hpp
//
// Purpose       : To declare the PartSM solid model class.
//
// Special Notes :
//
// Creator       : Madhan Narayanan
//
// Creation Date : 07/18/97
//
// Owner         : Madhan Narayanan
//-------------------------------------------------------------------------

#ifndef PARTSM_HPP
#define PARTSM_HPP

// ********** BEGIN STANDARD INCLUDES      **********
// ********** END STANDARD INCLUDES        **********

// ********** BEGIN MOTIF INCLUDES         **********
// ********** END MOTIF INCLUDES           **********

// ********** BEGIN OPEN INVENTOR INCLUDES **********
// ********** END OPEN INVENTOR INCLUDES   **********

// ********** BEGIN ACIS INCLUDES          **********
// ********** END ACIS INCLUDES            **********

// ********** BEGIN CUBIT INCLUDES         **********

#include "CubitDefines.h"

#include "CollectionEntity.hpp"

// ********** END CUBIT INCLUDES           **********

// ********** BEGIN FORWARD DECLARATIONS   **********
// ********** END FORWARD DECLARATIONS     **********

class PartSM : public CollectionEntity
{
   public :

      PartSM() ;
      //- The default constructor

      virtual ~PartSM() ;
      //- The destructor
  
  virtual const type_info& topology_entity_type_info() const;
   
   protected: 

   private:
} ;


// ********** BEGIN INLINE FUNCTIONS       **********
// ********** END INLINE FUNCTIONS         **********

// ********** BEGIN FRIEND FUNCTIONS       **********
// ********** END FRIEND FUNCTIONS         **********

// ********** BEGIN EXTERN FUNCTIONS       **********
// ********** END EXTERN FUNCTIONS         **********

#endif

