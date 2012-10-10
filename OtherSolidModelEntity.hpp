//-------------------------------------------------------------------------
// Filename      : OtherSolidModelEntity.hpp
//
// Purpose       : This class is the base class for all the other solid 
//                 model entities in the Geometry subsystem. 
//
// Special Notes :
//
// Creator       : Xuechen Liu
//
// Creation Date : 08/06/96
//
// Owner         : Raikanta Sahu
//-------------------------------------------------------------------------

#ifndef OTHER_SOLID_MODEL_ENTITY_HPP
#define OTHER_SOLID_MODEL_ENTITY_HPP

// ********** BEGIN STANDARD INCLUDES      **********
#include <assert.h>
// ********** END STANDARD INCLUDES        **********

// ********** BEGIN CUBIT INCLUDES         **********
#include "TopologyBridge.hpp"
#include "CubitDefines.h"
#include "DLIList.hpp"
// ********** END CUBIT INCLUDES           **********

// ********** BEGIN FORWARD DECLARATIONS   **********
class TopologyEntity;
// ********** END FORWARD DECLARATIONS     **********

class OtherSolidModelEntity : public TopologyBridge
{
public:
  
  OtherSolidModelEntity(); 
    // Default constructor.
  
  virtual ~OtherSolidModelEntity(); 
    // Destructor.
  
  virtual OtherSolidModelEntity* copy();
    //R OtherSolidModelEntity*
    //R- Pointer to the copied entity - created on the heap.
    // Makes a copy of this OSME entity *and* the solid model entities
    // it contains and returns a pointer to the new OSME object.
  
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
