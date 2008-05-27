/**
 * \file makept.cpp
 *
 * \brief makept, another simple C++ driver for CGM
 *
 * This program acts as a simple driver for CGM.  It reads in a geometry,
 * and performs varies checks for bodies, surfaces, curves and vertices.
 */
#include "GeometryModifyTool.hpp"
#include "GeometryModifyEngine.hpp"
#include "GeometryQueryTool.hpp"
#include "ModelQueryEngine.hpp"
#include "Body.hpp"
#include "RefFace.hpp"
#include "CubitObserver.hpp"

#include <typeinfo>
#include <assert.h>

#ifdef HAVE_ACIS
#include "AcisModifyEngine.hpp"
#include "AcisQueryEngine.hpp"
#endif

#ifdef HAVE_OCC
#include "OCCModifyEngine.hpp"
#include "OCCQueryEngine.hpp"
#endif

int test_sheet_query( GeometryModifyEngine* engine );

// main program - initialize, then send to proper function
int main (int argc, char **argv)
{
  CubitObserver::init_static_observers();
  CGMApp::instance()->startup( argc, argv );

  GeometryQueryEngine* const gqe_list[] = {
#ifdef HAVE_ACIS
    AcisQueryEngine::instance(),
#endif
#ifdef HAVE_OCC
    OCCQueryEngine::instance(),
#endif
    NULL
  };


  GeometryModifyEngine* const gme_list[] = {
#ifdef HAVE_ACIS
    AcisModifyEngine::instance(),
#endif
#ifdef HAVE_OCC
    OCCModifyEngine::instance(),
#endif
    NULL
  };

  int exit_val = 0;
  for (int i = 0; gme_list[i]; ++i)
    exit_val += test_sheet_query( gme_list[i] );

  return exit_val;
}


int test_sheet_query( GeometryModifyEngine* engine )
{
  CubitStatus rval = GeometryModifyTool::instance()->set_default_gme( engine );
  assert(rval);
  
  Body *sphere = 0, *sheet = 0;
  
  sphere = GeometryModifyTool::instance()->sphere( 1.0 );
  if (!sphere) {
    fprintf(stderr,"Sphere creation failed for engine (%s)\n", 
      typeid(*engine).name());
    return 1;
  }

  sheet = GeometryModifyTool::instance()->planar_sheet( 
    CubitVector(0,0,0), CubitVector(1,0,0), 
    CubitVector(1,1,0), CubitVector(0,1,0) );
  if (!sheet) {
    fprintf(stderr,"Planar sheet creation failed for engine (%s)\n", 
      typeid(*engine).name());
    return 1;
  }
  
  if (sphere->is_sheet_body()) {
    fprintf(stderr,"Solid sphere reported as sheet body (%s).\n", 
      typeid(*engine).name());
    return 1;
  }
  
  if (!sheet->is_sheet_body()) {
    fprintf(stderr,"Planar sheet reported as non-sheet body (%s).\n", 
      typeid(*engine).name());
    return 1;
  }
  
  printf("Success (%s)\n",typeid(*engine).name());
  return 0;
}

    
