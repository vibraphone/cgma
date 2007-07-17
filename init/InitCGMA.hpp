#ifndef INIT_CGMA_HPP
#define INIT_CGMA_HPP

#include "CubitDefines.h"

class InitCGMA
{
public:

  static CubitStatus initialize_cgma();
  
  static CubitStatus initialize_engine( const char* name );
};

#endif
