//-------------------------------------------------------------------------
// Filename      : PartitionCoSurf.hpp
//
// Purpose       : 
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 08/09/02
//-------------------------------------------------------------------------

#ifndef PARTITION_CO_SURF
#define PARTITION_CO_SURF

#include "CubitDefines.h"

class PartitionSurface;
class PartitionShell;

class PartitionCoSurf
{
  
  friend class PartitionSurface;
  friend class PartitionShell;
  
  public:
  
  PartitionCoSurf( CubitSense sense = CUBIT_FORWARD );
  ~PartitionCoSurf();
  
  PartitionSurface* get_surface() const;
  PartitionShell*   get_shell()   const;
  
  CubitSense sense() const;
  void sense( CubitSense set );
  
  PartitionCoSurf* next_in_surface() const;
  PartitionCoSurf* next_in_shell() const;
  
  short mark;
  
  private:
  
  CubitSense mySense;
  
  PartitionSurface* mySurface;
  PartitionCoSurf* surfaceNext;
  
  PartitionShell* myShell;
  PartitionCoSurf* shellNext;
};

inline PartitionCoSurf::PartitionCoSurf( CubitSense sense )
  : mark(0),
    mySense( sense ),
    mySurface(0),
    surfaceNext(0),
    myShell(0),
    shellNext(0)
  {}

inline PartitionSurface* PartitionCoSurf::get_surface() const
  { return mySurface; }

inline PartitionShell* PartitionCoSurf::get_shell() const
  { return myShell; }

inline CubitSense PartitionCoSurf::sense() const
  { return mySense; }

inline void PartitionCoSurf::sense( CubitSense set )
 { mySense = set; }

#endif
