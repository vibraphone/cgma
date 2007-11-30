//-------------------------------------------------------------------------
// Filename      : OCCCoEdge.hpp
//
// Purpose       : 
//
// Special Notes :
//
// Creator       : Steven J. Owen
//
// Creation Date : 07/23/00
//
// Owner         : Steven J. Owen
//-------------------------------------------------------------------------

#ifndef OCCCOEDGE_HPP
#define OCCCOEDGE_HPP

// ********** BEGIN STANDARD INCLUDES      **********
// ********** END STANDARD INCLUDES        **********

// ********** BEGIN CUBIT INCLUDES         **********
#include "CubitDefines.h"
#include "CubitEntity.hpp"
#include "CoEdgeSM.hpp"
#include <TopoDS_Edge.hxx>
// ********** END CUBIT INCLUDES           **********

// ********** BEGIN FORWARD DECLARATIONS   **********
class TopologyEntity;
class LoopSM;

class OCCBody;
class OCCLump;
class OCCShell;
class OCCSurface;
class OCCLoop;
class OCCCurve;
class OCCPoint;

// ********** END FORWARD DECLARATIONS     **********

class OCCCoEdge 
{
public:
  
  OCCCoEdge(Curve *curv_ptr, LoopSM *loop_ptr, CubitSense sense);
    //- A constructor
  
  virtual ~OCCCoEdge() ;
    //- The destructor
    
  inline CubitSense sense(){return edgeSense;}
    //- returns the sense of the underlying coedge wrt the underlying edge

  inline Curve *curve()
    {return myCurve;}
    //- get the curve associated with this coedge
    
  inline LoopSM* loop() const { return myLoop; }

protected: 
  
private:
  LoopSM *myLoop;
  Curve *myCurve;
  CubitSense edgeSense;
};


// ********** BEGIN INLINE FUNCTIONS       **********
// ********** END INLINE FUNCTIONS         **********

// ********** BEGIN FRIEND FUNCTIONS       **********
// ********** END FRIEND FUNCTIONS         **********

// ********** BEGIN EXTERN FUNCTIONS       **********
// ********** END EXTERN FUNCTIONS         **********

#endif

