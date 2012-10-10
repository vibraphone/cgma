//-------------------------------------------------------------------------
// Filename      : Csys.hpp
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

#ifndef CSYS_HPP
#define CSYS_HPP

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
#include "OtherEntity.hpp"
#include "CubitVector.hpp"


// ********** END CUBIT INCLUDES           **********

// ********** BEGIN FORWARD DECLARATIONS   **********
// ********** END FORWARD DECLARATIONS     **********

class Csys : public OtherEntity
{
   public :

      Csys() ;
      //- The default constructor

      virtual ~Csys() ;
      //- The destructor
      
      virtual CubitVector origin() const = 0;  // pure virtual
      //R- origin of the coordinate system
      
      
      

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

