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
#include "InitCGMA.hpp"
#include "DLIList.hpp"

#include <typeinfo>

int test_sheet_query( GeometryModifyEngine* engine );

// main program - initialize, then send to proper function
int main (int argc, char **argv)
{
  CubitStatus result = InitCGMA::initialize_cgma();
  if (CUBIT_SUCCESS != result) return 1;
  
  DLIList<GeometryModifyEngine*> gme_list;
  GeometryModifyTool::instance()->get_gme_list(gme_list);

  int exit_val = 0;
  for (int i = 0; i < gme_list.size(); ++i)
    exit_val += test_sheet_query( gme_list[i] );

  return exit_val;
}


int test_sheet_query( GeometryModifyEngine* engine )
{
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

    
