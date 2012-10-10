#include "CompositeCoSurf.hpp"
#include "CompositeShell.hpp"
#include "CompositeSurface.hpp"

CompositeCoSurf::~CompositeCoSurf()
{
  if( myShell )
  {
    myShell->remove( this );
    myShell = 0;
  }
  if( mySurface )
  {
    mySurface->remove( this );
    mySurface = 0;
  }
  
  mySense = CUBIT_UNKNOWN;
  assert( !shellNext && !surfaceNext );
}

//-------------------------------------------------------------------------
// Purpose       : Print debug output
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 09/16/04
//-------------------------------------------------------------------------
void CompositeCoSurf::print_debug_info( const char* prefix, bool brief )
{
  if( prefix == 0 ) prefix = "";
  
  const char* sense = mySense == CUBIT_FORWARD ? "Forward" :
                      mySense == CUBIT_REVERSED ? "Reverse" : "UNKNOWN";
                      
  PRINT_INFO("%sCompCoSurf %p %s ", prefix, this, sense );
   
  if( !mySurface )
    PRINT_INFO("NULL SURFACE\n");
  else if( brief )
    PRINT_INFO("surface %p\n", mySurface );
  else
    { PRINT_INFO("\n  ");  mySurface->print_debug_info(prefix, true); }
}
