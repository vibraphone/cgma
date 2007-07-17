#include "InitCGMA.hpp"
#include "CGMApp.hpp"
#include "VirtualQueryEngine.hpp"
#include "FacetQueryEngine.hpp"
#include "FacetModifyEngine.hpp"

#include <ctype.h>

#ifdef CUBIT_CGM
#include <new>

class AcisModifyEngine
{
public:
  AcisModifyEngine();
  static AcisModifyEngine *instance_;
};

class AcisQueryEngine
{
public:
  AcisQueryEngine();
  static AcisQueryEngine *instance_;
};

class dummym
{
  char dummyvar[4096];
};

#elif defined(HAVE_ACIS) 
#  include "AcisQueryEngine.hpp"
#  include "AcisModifyEngine.hpp"
#endif


static bool streq_nocase( const char* s, const char* t )
{
  for (; *s; ++s, ++t) 
    if (tolower(*s) != tolower(*t))
      return false;
  return !*t;
}


CubitStatus InitCGMA::initialize_cgma()
{
  CGMApp::instance()->startup( 0, NULL );
  return CUBIT_SUCCESS;
}


CubitStatus InitCGMA::initialize_engine( const char* name )
{
  if (streq_nocase(name,"ACIS")) {
#ifdef CUBIT_CGM
    if (!AcisQueryEngine::instance_)
      AcisQueryEngine::instance_ = new (reinterpret_cast<AcisQueryEngine*>(new dummym)) AcisQueryEngine;
    if (!AcisModifyEngine::instance_)
      AcisModifyEngine::instance_ = new (reinterpret_cast<AcisModifyEngine*>(new dummym)) AcisModifyEngine;
#elif defined (HAVE_ACIS)
    AcisQueryEngine::instance();
    AcisModifyEngine::instance();
#else
    return CUBIT_FAILURE;
#endif
  }
  
  else if (streq_nocase(name,"VIRTUAL") || streq_nocase(name,"vg")) {
    VirtualQueryEngine::instance();
  }
  
  else if (streq_nocase(name,"facet")) {
    FacetQueryEngine::instance();
    FacetModifyEngine::instance();
  }
  
  else {
    return CUBIT_FAILURE;
  }
  
  return CUBIT_SUCCESS;
}

