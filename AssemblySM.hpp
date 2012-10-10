//-------------------------------------------------------------------------
// Filename      : AssemblySM.hpp
//
// Purpose       : To declare the AssemblySM solid model class.
//
// Special Notes :
//
// Creator       : Madhan Narayanan
//
// Creation Date : 07/18/97
//
// Owner         : Madhan Narayanan
//-------------------------------------------------------------------------

#ifndef ASSEMBLYSM_HPP
#define ASSEMBLYSM_HPP

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

class AssemblySM : public CollectionEntity
{
   public :

      AssemblySM() ;
      //- The default constructor

      virtual ~AssemblySM() ;
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

