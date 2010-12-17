#include "MergeTool.hpp"
#include "InitCGMA.hpp"
#include "GeometryModifyTool.hpp"
#include "GeometryQueryTool.hpp"
#include "Body.hpp"
#include "CGMApp.hpp"
#include "CubitAttribManager.hpp"

#include <stdio.h>
#include <stdlib.h>

#ifdef TEST_ACIS
#  define ENGINE "ACIS"
#  define FORMAT "ACIS_SAT"
#elif defined (TEST_OCC)
#  define ENGINE "OCC"
#  define FORMAT "OCC"
#elif defined (TEST_FACET)
#  define ENGINE "FACET"
#  define FORMAT "FACET"
#else
#  error "Which engine to test?"
#endif
#define FILE_NAME "merge_test." FORMAT

#define ASSERT(A) if (!(A)) failed(#A,__FILE__,__LINE__)
void failed( const char* A, const char* FILE, int LINE )
{
  printf( "Condition failed at %s:%d : %s\n", FILE, LINE, A );
  abort();
}

RefFace* shared_face( Body* bod1, Body* bod2 );

int main()
{
  CubitStatus s = InitCGMA::initialize_cgma( ENGINE );
  ASSERT(s);
  CGMApp::instance()->attrib_manager()->auto_flag( CUBIT_TRUE );

  Body *brick1, *brick2;
  brick1 = GeometryModifyTool::instance()->brick( 1, 1, 1 );
  brick2 = GeometryModifyTool::instance()->brick( 1, 1, 1 );
  ASSERT(brick1 && brick2);
  
  s = GeometryQueryTool::instance()->translate( brick2, CubitVector(1,0,0) );
  ASSERT(s);
  
  DLIList<Body*> merge_list;
  merge_list.append( brick1 );
  merge_list.append( brick2 );
  s = MergeTool::instance()->merge_bodies( merge_list );
  ASSERT(s);
  
  ASSERT( shared_face( brick1, brick2 ) );
  
  int junk;
  DLIList<RefEntity*> export_list;
  export_list.append( brick1 );
  export_list.append( brick2 );
  s = GeometryQueryTool::instance()->export_solid_model( export_list,
                                                         FILE_NAME,
                                                         FORMAT,
                                                         junk,
                                                         CubitString(__FILE__) );
  ASSERT(s);
  
  
  DLIList<RefEntity*> import_list;
  s = GeometryQueryTool::instance()->import_solid_model( FILE_NAME,
                                                         FORMAT,
                                                         NULL,
                                                         CUBIT_FALSE,
                                                         CUBIT_TRUE,
                                                         CUBIT_TRUE,
                                                         CUBIT_TRUE,
                                                         CUBIT_TRUE,
                                                         CUBIT_TRUE,
                                                         &import_list );
  ASSERT(s);
  remove( FILE_NAME );
  
  DLIList<Body*> import_bodies;
  CAST_LIST( import_list, import_bodies, Body );
  ASSERT( 2 == import_bodies.size() );
  
  Body* new_bod1 = import_bodies.get_and_step();
  Body* new_bod2 = import_bodies.get_and_step();
  ASSERT( brick1 != new_bod1 && brick1 != new_bod2 );
  ASSERT( brick2 != new_bod1 && brick2 != new_bod2 );
  ASSERT( shared_face( new_bod1, new_bod2 ) );
  
  return 0;
}

RefFace* shared_face( Body* b1, Body* b2 )
{
  DLIList<RefFace*> faces1, faces2;
  b1->ref_faces( faces1 );
  b2->ref_faces( faces2 );
  faces1.intersect( faces2 );
  if (faces1.size() == 1)
    return faces1.get();
  else
    return 0;
}

