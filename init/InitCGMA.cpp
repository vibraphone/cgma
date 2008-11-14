#include "InitCGMA.hpp"
#include "CGMApp.hpp"
#include "VirtualQueryEngine.hpp"
#include "FacetQueryEngine.hpp"
#include "FacetModifyEngine.hpp"
#include "GeometryQueryTool.hpp"
#include "GeometryModifyTool.hpp"

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

#if defined(HAVE_OCC) 
#  include "OCCQueryEngine.hpp"
#  include "OCCModifyEngine.hpp"
#endif


static bool streq_nocase( const char* s, const char* t )
{
  for (; *s; ++s, ++t) 
    if (tolower(*s) != tolower(*t))
      return false;
  return !*t;
}


CubitStatus InitCGMA::initialize_cgma( const char* default_engine_name )
{
  CGMApp::instance()->startup( 0, NULL );
  GeometryModifyEngine* default_engine = 0;

#if defined(CUBIT_CGM) || defined(HAVE_ACIS)
  GeometryModifyEngine* acis_engine_ptr;
  #ifdef CUBIT_CGM
    if (!AcisQueryEngine::instance_)
      AcisQueryEngine::instance_ = new (reinterpret_cast<AcisQueryEngine*>(new dummym)) AcisQueryEngine;
    if (!AcisModifyEngine::instance_)
      AcisModifyEngine::instance_ = new (reinterpret_cast<AcisModifyEngine*>(new dummym)) AcisModifyEngine;
    acis_engine_ptr = reinterpret_cast<GeometryModifyEngine*>
                         (AcisModifyEngine::instance_);
  #else
    AcisQueryEngine::instance();
    AcisModifyEngine::instance();
    acis_engine_ptr = AcisModifyEngine::instance();
  #endif
  if (default_engine_name && streq_nocase("ACIS",default_engine_name))
    default_engine = acis_engine_ptr;
#endif

#ifdef HAVE_OCC  
  OCCQueryEngine::instance();
  OCCModifyEngine::instance();
  if (default_engine_name && streq_nocase("OCC",default_engine_name))
    default_engine = OCCModifyEngine::instance();
#endif  

  FacetQueryEngine::instance();
  FacetModifyEngine::instance();
  VirtualQueryEngine::instance();

  if (default_engine_name && streq_nocase("FACET",default_engine_name))
    default_engine = FacetModifyEngine::instance();

  if(default_engine_name) {
    if (!default_engine) {
      PRINT_ERROR("Invalid or unsupported engine: '%s'\n", default_engine_name);
      return CUBIT_FAILURE;
    }
    
    CubitStatus rval;
    rval = GeometryModifyTool::instance()->set_default_gme(default_engine);
    if (CUBIT_SUCCESS != rval)
      return rval;
    rval = GeometryQueryTool::instance()->set_default_gqe(default_engine->get_gqe());
    if (CUBIT_SUCCESS != rval)
      return rval;
  }
  
  return CUBIT_SUCCESS;
}
