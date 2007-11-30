//-------------------------------------------------------------------------
// Filename      : OCCCoEdge.cpp
//
// Purpose       : 
//
// Special Notes :
//
// Creator       : Steven J. Owen
//
// Creation Date : 07/18/00
//
// Owner         : Steven J. Owen
//-------------------------------------------------------------------------

// ********** BEGIN STANDARD INCLUDES      **********
// ********** END STANDARD INCLUDES        **********
#include "config.h"
// ********** BEGIN CUBIT INCLUDES         **********
#include "CubitDefines.h"
#include "CastTo.hpp"
#include "OCCCoEdge.hpp"
#include "OCCLoop.hpp"
#include "OCCQueryEngine.hpp"
#include "CubitUtil.hpp"

#include "OCCBody.hpp"
#include "OCCLump.hpp"
#include "OCCShell.hpp"
#include "OCCSurface.hpp"
#include "OCCCurve.hpp"
// ********** END CUBIT INCLUDES           **********

// ********** BEGIN FORWARD DECLARATIONS   **********
// ********** END FORWARD DECLARATIONS     **********

// ********** BEGIN STATIC DECLARATIONS    **********
// ********** END STATIC DECLARATIONS      **********

// ********** BEGIN PUBLIC FUNCTIONS       **********


OCCCoEdge::OCCCoEdge( Curve *curv_ptr, 
	  	      LoopSM *loop_ptr, CubitSense sense )
	 	    : myCurve(curv_ptr), 
 	  	      myLoop(loop_ptr),edgeSense(sense)
{
}

//-------------------------------------------------------------------------
// Purpose       : The destructor
//
// Special Notes :
//
// Creator       : Steve Owen
//
// Creation Date : 07/18/00
//-------------------------------------------------------------------------
OCCCoEdge::~OCCCoEdge()
{
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

