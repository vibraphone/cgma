//-------------------------------------------------------------------------
// Filename      : CompositeCoSurf.hpp
//
// Purpose       : CoSurface for composite topology
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 08/07/02
//-------------------------------------------------------------------------

#ifndef COMPOSITE_CO_SURFACE
#define COMPOSITE_CO_SURFACE

#include "CubitDefines.h"

class CompositeSurface;
class CompositeShell;

class CompositeCoSurf
{

  friend class CompositeSurface;
  friend class CompositeShell;

  public:

    CompositeCoSurf( CubitSense sense = CUBIT_FORWARD );
    ~CompositeCoSurf();

    CompositeSurface* get_surface() const;
    CompositeShell* get_shell() const;
    
    CubitSense sense() const;
    void sense( CubitSense set );
    
    CompositeCoSurf* next_in_surface() const;
    CompositeCoSurf* next_in_shell() const;
  
    void print_debug_info( const char* line_prefix = 0, bool brief = false );
    
  private:
  
    CubitSense mySense;
    
    CompositeSurface* mySurface;
    CompositeCoSurf* surfaceNext;
    
    CompositeShell* myShell;
    CompositeCoSurf* shellNext;
};

inline CompositeCoSurf::CompositeCoSurf( CubitSense sense ) 
  : mySense(sense),
    mySurface(0),
    surfaceNext(0),
    myShell(0),
    shellNext(0)
  {}

inline CompositeSurface* CompositeCoSurf::get_surface() const
  { return mySurface; }

inline CompositeShell* CompositeCoSurf::get_shell() const
  { return myShell; }

inline CompositeCoSurf* CompositeCoSurf::next_in_surface() const
  { return surfaceNext; }

inline CompositeCoSurf* CompositeCoSurf::next_in_shell() const
  { return shellNext; }

inline CubitSense CompositeCoSurf::sense() const
  { return mySense; }

inline void CompositeCoSurf::sense( CubitSense set ) 
  { mySense = set; }



#endif


    
