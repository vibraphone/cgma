#include "InitCGMA.hpp"
#include "GeometryModifyTool.hpp"
#include "GeometryQueryTool.hpp"
#include "Body.hpp"
#include "CGMApp.hpp"
#include "CubitAttribManager.hpp"
#include "CADefines.hpp"
#include "TDParallel.hpp"
#include "CABodies.hpp"

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
  CGMApp::instance()->attrib_manager()->register_attrib_type(CA_BODIES, "bodies", "BODIES",
							     CABodies_creator, CUBIT_TRUE,
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

    // add tool data to bodies
    RefEntity* entity = dynamic_cast<RefEntity*> (bricks[i]);
    DLIList<int> shared_procs;
    shared_procs.append(i);
    TDParallel *td_par = (TDParallel *) entity->get_TD(&TDParallel::is_parallel);
    if (td_par == NULL) td_par = new TDParallel(entity, NULL, &shared_procs);

    // add tool data to volumes too
    DLIList<RefVolume*> volumes;
    (dynamic_cast<TopologyEntity*> (entity))->ref_volumes(volumes);
    int n_vol = volumes.size();
    volumes.reset();
    for (int j = 0; j < n_vol; j++) {
      RefEntity *vol = volumes.get_and_step();
      td_par = (TDParallel *) vol->get_TD(&TDParallel::is_parallel);
      if (td_par == NULL) td_par = new TDParallel(vol, NULL, &shared_procs);
    }
  }

  // export as file
  int junk;
  s = GeometryQueryTool::instance()->export_solid_model( export_list, "bricks2.occ",
                                                         FORMAT, junk, CubitString(__FILE__) );
  ASSERT(s);

  // delete geometry
  GeometryQueryTool::instance()->delete_geometry();

  // import it again
  DLIList<RefEntity*> import_list;
  s = GeometryQueryTool::instance()->import_solid_model( "bricks2.occ",
                                                         FORMAT, NULL, CUBIT_FALSE,
                                                         CUBIT_TRUE, CUBIT_TRUE,
                                                         CUBIT_TRUE, CUBIT_TRUE,
                                                         CUBIT_TRUE, &import_list
                                                         );
  ASSERT(s);

  // check imported entity has tool data actuated by attributes
  DLIList<RefEntity*> body_entity_list;
  s = GeometryQueryTool::instance()->ref_entity_list("body", body_entity_list, CUBIT_FALSE);
  int n_body = body_entity_list.size();
  body_entity_list.reset();

  for (int i = 0; i < n_body; i++) {
    RefEntity* entity = body_entity_list.get_and_step();
    TDParallel *td_par = (TDParallel *) entity->get_TD(&TDParallel::is_parallel);

    // check bodies
    if (td_par == NULL) {
      std::cout << "Error: Body doesn't have tool data." << std::endl;

      // check child volumes
      DLIList<RefEntity*> volumes;
      entity->get_child_ref_entities(volumes);

      // check if the first Volume has tool data
      volumes.reset();
      RefEntity *vol = volumes.get();
      td_par = (TDParallel *) vol->get_TD(&TDParallel::is_parallel);

      if (td_par == NULL) {
	std::cout << "Error: Volume doesn't have tool data." << std::endl;
	return 1;
      }
    }
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

