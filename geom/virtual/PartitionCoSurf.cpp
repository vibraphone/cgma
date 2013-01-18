#include "PartitionCoSurf.hpp"
#include "PartitionShell.hpp"
#include "PartitionSurface.hpp"

PartitionCoSurf::~PartitionCoSurf()
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

