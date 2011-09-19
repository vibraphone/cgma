#include "InitCGMA.hpp"
#include "GeometryModifyTool.hpp"
#include "GeometryQueryTool.hpp"
#include "Body.hpp"
#include "CGMApp.hpp"
#include "CubitAttribManager.hpp"
#include "CAEntityId.hpp"

#include <iostream>
#include <string>
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

#define ASSERT(A) if (!(A)) failed(#A,__FILE__,__LINE__)
void failed( const char* A, const char* FILE, int LINE )
{
  printf( "Condition failed at %s:%d : %s\n", FILE, LINE, A );
  abort();
}

int main()
{
  CubitStatus s = InitCGMA::initialize_cgma( ENGINE );
  ASSERT(s);
  
  // actuate CA_BODIES attribute and turn on auto flag
  CGMApp::instance()->attrib_manager()->register_attrib_type(CA_ENTITY_ID, "id", "ENTITY_ID",
							     CAEntityId_creator, CUBIT_TRUE,
							     CUBIT_TRUE, CUBIT_TRUE, CUBIT_TRUE,
							     CUBIT_TRUE, CUBIT_FALSE);
  CGMApp::instance()->attrib_manager()->auto_flag(CUBIT_TRUE);

  // make 2 bricks
  int n = 2;
  Body** bricks = new Body*[n];
  DLIList<RefEntity*> export_list;
  for (int i = 0; i < n; i++) {
    bricks[i] = GeometryModifyTool::instance()->brick( 1, 1, 1 );
    ASSERT(bricks[i]);
    s = GeometryQueryTool::instance()->translate( bricks[i], CubitVector(i,0,0) );
    ASSERT(s);
    export_list.append( bricks[i] );
  }

  // export as buffer
  char* p_buffer;
  int n_bsize = 0;
  s = GeometryQueryTool::instance()->export_solid_model(export_list, p_buffer,
                                                        n_bsize, false);
  void *temp_buffer = malloc(n_bsize);
  p_buffer = (char *) temp_buffer;
  s = GeometryQueryTool::instance()->export_solid_model(export_list, p_buffer,
                                                        n_bsize, true);
  ASSERT(s);

  // delete geometry
  GeometryQueryTool::instance()->delete_geometry();

  //check that the two single volume bodys' attributes are exported as SINGLELUMP%

  std::string search ("SINGLELUMP%") ;
  std::string buffer (p_buffer, n_bsize);
  size_t found = 0;
  found = buffer.find(search);
  if(found==std::string::npos)
    assert(0);

  // import it again
  DLIList<RefEntity*> import_list;
  s = GeometryQueryTool::instance()->import_solid_model(NULL, p_buffer, n_bsize);
  ASSERT(s);
  delete p_buffer;

  // check imported entity has tool data actuated by attributes
  DLIList<RefEntity*> body_entity_list;
  s = GeometryQueryTool::instance()->ref_entity_list("body", body_entity_list, CUBIT_FALSE);
  int n_body = body_entity_list.size();
  if (n_body != 2) {
    std::cerr << "Error: # of imported bodies should be 2." << std::endl;
    return 1;
  }

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

