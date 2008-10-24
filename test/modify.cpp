/**
 * \file makept.cpp
 *
 * \brief makept, another simple C++ driver for CGM
 *
 * This program acts as a simple driver for CGM.  It reads in a geometry,
 * and performs varies checks for bodies, surfaces, curves and vertices.
 */
#include "config.h"
#include "CpuTimer.hpp"
#include "GeometryModifyTool.hpp"
#include "GeometryQueryTool.hpp"
#include "OCCQueryEngine.hpp"
#include "CubitUtil.hpp"
#include "CubitMessage.hpp"
#include "CubitDefines.h"
#include "RefEntity.hpp"
#include "Body.hpp"
#include "RefVolume.hpp"
#include "RefFace.hpp"
#include "RefEdge.hpp"
#include "RefVertex.hpp"
#include "CubitObserver.hpp"
#include "CastTo.hpp"
#include "OCCModifyEngine.hpp"
#include "AppUtil.hpp"
#include "RefEntityFactory.hpp"
#include "RefEdge.hpp"
#include "BodySM.hpp"
#include "Lump.hpp"
#include "OCCLump.hpp"
#include "OCCBody.hpp"
#include "OCCSurface.hpp"
#include "OCCCurve.hpp"
#include "OCCShell.hpp"
#include "TopoDS_Shape.hxx"
#include "RefEntityName.hpp"
#include "RefEntityFactory.hpp"

#ifndef SRCDIR
# define SRCDIR .
#endif

#define STRINGIFY_(X) #X
#define STRINGIFY(X) STRINGIFY_(X)
#define SRCPATH STRINGIFY(SRCDIR) "/"

// forward declare some functions used and defined later
CubitStatus read_geometry(int, char **, bool local = false);
CubitStatus make_Point();
// macro for printing a separator line
#define PRINT_SEPARATOR   PRINT_INFO("=======================================\n");


// main program - initialize, then send to proper function
int main (int argc, char **argv)
{

  CubitObserver::init_static_observers();
    // Initialize the GeometryTool
  
  CGMApp::instance()->startup( argc, argv );
  OCCQueryEngine::instance();
  OCCModifyEngine::instance();

    // If there aren't any file arguments, print usage and exit
  //if (argc == 1) {
  //  PRINT_INFO("Usage: mergechk <geom_file> [<geom_file> ...]\n");
  //  exit(0);
  //}
  
  CubitStatus status = CUBIT_SUCCESS;



  //Do make point.
  status = make_Point();
  if (status == CUBIT_FAILURE) 
     PRINT_INFO("Operation Failed");

  int ret_val = ( CubitMessage::instance()->error_count() );
  if ( ret_val > 0 )
  {
    PRINT_ERROR("Errors found during Mergechk session.\n");
  }
  return ret_val;
  
}

/// attribs module: list, modify attributes in a give model or models
/// 
/// Arguments: file name(s) of geometry files in which to look
///
CubitStatus read_geometry(int num_files, char **argv, bool local) 
{
  CubitStatus status = CUBIT_SUCCESS;
  GeometryQueryTool *gti = GeometryQueryTool::instance();
  assert(gti);
  int i;
  
  PRINT_SEPARATOR;

  for (i = 0; i < num_files; i++) {
    std::string filename( local ? "./" : SRCPATH );
    filename += argv[i];
    status = gti->import_solid_model(filename.c_str(), "OCC");
    if (status != CUBIT_SUCCESS) {
      PRINT_ERROR("Problems reading geometry file %s.\n", filename.c_str());
      abort();
    }
  }
  PRINT_SEPARATOR;

  return CUBIT_SUCCESS;
}

CubitStatus make_Point()
{
  RefEntityFactory* ref = RefEntityFactory::instance();
  GeometryQueryTool *gti = GeometryQueryTool::instance();
  GeometryModifyTool *gmti = GeometryModifyTool::instance();

  OCCQueryEngine::instance();
  OCCModifyEngine* ome = OCCModifyEngine::instance();

  Body* body = gmti->brick(10, 10, 10);
  CubitVector v(15,0,0);
  BodySM* bodysm = body->get_body_sm_ptr();
  DLIList<OCCSurface*> occ_surfaces;
  CAST_TO(bodysm, OCCBody)->get_all_surfaces(occ_surfaces);  
  DLIList<RefFace*> ref_faces;
  DLIList<RefEdge*> ref_edges;
  body->ref_faces(ref_faces);

  DLIList<RefFace*> faces_to_stitch;
  RefFace* test_face = NULL;
  for(int i = 0 ; i < ref_faces.size(); i++)
  {
    RefFace* refface = gmti->make_RefFace(ref_faces.get_and_step());
    if(refface)
	faces_to_stitch.append(refface);
  }

  DLIList<Surface*> surface_list;
  test_face = gmti->make_RefFace(ref_faces.get());
  surface_list.append(test_face->get_surface_ptr());
  CubitVector v_test(0,5,5);
  CubitVector normal;
  normal = test_face->normal_at(v_test); //(-1,0,0)
  ome->flip_normals(surface_list);
  normal = test_face->normal_at(v_test); //(1,0, 0)
  DLIList<Body*> bodies;
  test_face->bodies(bodies);
  gti->delete_Body(bodies);
  surface_list.clean_out();
  bodies.clean_out(); 
  gti->delete_Body(body);

  DLIList<BodySM*> bodysm_list;
  DLIList<RefFace*> face_list;
  DLIList<RefVertex*> vertices;
  for(int i = 0; i < faces_to_stitch.size(); i++)
  {
    //move each refface by (15,0,0)
    RefFace* refface = faces_to_stitch.get_and_step();
    refface->ref_vertices(vertices);
    Body* body = refface->ref_volume()->get_body_ptr();
    bodysm = body->get_body_sm_ptr();
    Surface* surface = refface->get_surface_ptr();
    surface_list.append(surface);
    bodysm_list.append(bodysm);
    gti->translate(body, v);
    body->ref_faces(face_list);
    occ_surfaces.clean_out();
    CAST_TO(bodysm, OCCBody)->get_all_surfaces(occ_surfaces);
  }

  //create solid from surfaces 
  GeometryModifyEngine *gme = gmti->get_engine(occ_surfaces.get());
  BodySM* stitched_body = NULL;
  DLIList<Body*> new_bodies;
  gmti->create_solid_bodies_from_surfs(face_list, new_bodies);
  //Created volume(s): 9 
  //Destroyed volume(s): 2 to 7

  CubitStatus rsl = CUBIT_SUCCESS;
  DLIList<RefEntity*> ref_entity_list;
  int num_ents_exported=0;
  const CubitString cubit_version="10.2";
  const char * filename = "stitch.occ";
  const char * filetype = "OCC";

  rsl = gti->export_solid_model(ref_entity_list, filename, filetype,
                                 num_ents_exported, cubit_version);

  DLIList<RefEntity*>  free_entities;
  gti->bodies(bodies);
  gti->get_free_ref_entities(free_entities);
  //delete all entities
  gti->delete_Body(bodies);

  for (int j = free_entities.size(); j--;)
    {
      gti->delete_RefEntity( free_entities.get_and_step());
    }

  // Read in the geometry from files specified on the command line
  char *argv = "stitch.occ";
  CubitStatus status = read_geometry(1, &argv, true);
  if (status == CUBIT_FAILURE) exit(1);
  //Read in 1 volume.

  bodies.clean_out();
  free_entities.clean_out();
  gti->bodies(bodies); 
  gti->get_free_ref_entities(free_entities);

  //delete all entities
  gti->delete_Body(bodies);

  for (int j = free_entities.size(); j--;)
    {
      gti->delete_RefEntity( free_entities.get_and_step());
    }

  //test for cylinder making
  Body* from_body = gmti->cylinder(10, 4, 3, 2);
  Body* from_body2 = gmti->cylinder(8, 4, 2, 0);
  Body* tool_body = gmti->cylinder(10, 1, 1, 1);
  double d;
  d = from_body->measure(); //d = 219.91
  d = from_body2->measure();//d = 67.02
  d = tool_body->measure(); //d = 31.41

  // test for sphere making
  Body* test_body = gmti->sphere(5);
  d = test_body->measure(); //d = 523.598 

  //test for prism making
  test_body = gmti->prism(10, 4, 4,2);
  d =  test_body->measure(); //d = 320

  //test for pyramid making
  test_body = gmti->pyramid(10, 4, 5, 2, 3);
  d =  test_body->measure(); //d = 320

  //test for torus making
  test_body =  gmti->torus(10,5);
  d =  test_body->measure(); //d = 4934.8

  //test for planar_sheet making
  CubitVector p1(0, 0, 0);
  CubitVector p2(0, 0, 1);
  CubitVector p3(0, 1, 1);
  CubitVector p4(0, 1, 0);
  test_body =  gmti->planar_sheet(p1, p2, p3, p4);
  d =  test_body->measure(); //d = 1

  bodies.clean_out();
  gti->bodies(bodies);
  //delete all entities
  gti->delete_Body(bodies); 

  //test for subtract
  from_body = gmti->brick(10, 10, 10);
  int width = 10; //we can also test for width < 10
  from_body2 = gmti->brick(4, width, 4);
  tool_body = gmti->brick(10, 10, 10);  
  CubitVector v_move(1,0,-1);
  CubitVector v_movei(0,0,1);
  gti->translate(from_body2,v_move);
  gti->translate(tool_body, v_movei);
  DLIList<Body*> from_bodies;
  from_bodies.append(from_body);
  new_bodies.clean_out();
  rsl = gmti->subtract(from_body2, from_bodies, new_bodies, 
                       CUBIT_TRUE, CUBIT_FALSE);
  //Updated volume(s): 19
  //Destroyed volume(s): 20
  //new bodies has one body, new body has 10 ref-faces, 5 of them are remaining
  //with old id, 5 of them are new faces.

  from_bodies=new_bodies;
  new_bodies.clean_out();
  rsl = gmti->subtract(tool_body,from_bodies, new_bodies,
                       CUBIT_TRUE, CUBIT_FALSE);
  //Created volume(s): 22, 23
  //Destroyed volume(s): 19, 21
  d = new_bodies.step_and_get()->measure();
  v = new_bodies.get()->center_point();
  int n = new_bodies.get()->num_ref_faces();
  // n = 6
  //new bodies has 2 bodies, one has a volume = 10 and the other has a 
  //volume = 50; each of them has 6 ref_faces, of which 3 are new and 3 are
  //remaining (unchanged or modified).

  bodies.clean_out();
  gti->bodies(bodies);
  //delete all entities
  gti->delete_Body(bodies); 
  
  free_entities.clean_out();
  gti->get_free_ref_entities(free_entities);
  assert(free_entities.size() == 0);
  //there shouldn't be any free_entites.

  from_body = gmti->brick(10, 10, 10);
  tool_body = gmti->brick(11, 1, 1);
  CubitVector v_move4(0,1,-1);
  gti->translate(from_body,v_move4);
  Body* cp_from_body = gmti->copy_body(from_body);
  Body* cp_from_body2 = gmti->copy_body(cp_from_body);

  //test edge imprint on body
  ref_edges.clean_out();
  tool_body->ref_edges(ref_edges);
  from_bodies.clean_out();
  from_bodies.append(cp_from_body);
  new_bodies.clean_out();
  CubitStatus stat = gmti->imprint(from_bodies, ref_edges, new_bodies, CUBIT_FALSE, CUBIT_TRUE );
  //Updated volume(s): 26

  //test edge imprint on surface
  CubitVector vv(5,1,4.0);
  face_list.clean_out();
  cp_from_body2->ref_faces(face_list);
  int size = face_list.size();
  DLIList<RefFace*> unimprint_faces;
  for(int i = 0; i < size; i++)
  {
    CubitVector v = face_list.get()->center_point();
    if(!v.about_equal(vv))
      unimprint_faces.append(face_list.remove());
    else
      face_list.step();
  }
  assert(face_list.size() == 1);
  ref_edges.step_and_get()->ref_vertices(vertices);
  stat = gmti->imprint(unimprint_faces, ref_edges, new_bodies, CUBIT_FALSE);
  //Updated volume(s): 27

  new_bodies.clean_out();
  face_list.clean_out();
  cp_from_body2->ref_faces(face_list);
  for(int i = 0; i < size; i++)
  {
    CubitVector v = face_list.get()->center_point();
    if(!v.about_equal(vv))
      face_list.remove();
    else
      face_list.step();
  }

  //imprint a point on an edge, split it
  CubitVector pp1(10,1,8);
  CubitVector pp2(10,5,9);
  CubitVector pp3(10,1,6);
  DLIList<CubitVector*> vectors;
  vectors.append(&pp1);
  vectors.append(&pp2);
  vectors.append(&pp3);
  from_bodies.clean_out();
  from_bodies.append(cp_from_body2);
  new_bodies.clean_out(); 
  stat = gmti->imprint(from_bodies, vectors, new_bodies, CUBIT_FALSE);
  //Updated volume(s): 27

  n = new_bodies.get()->num_ref_edges();//n = 17

  new_bodies.clean_out();
  stat = gmti->imprint(face_list, ref_edges, new_bodies, CUBIT_FALSE);
  //Updated volume(s): 27

  n = new_bodies.get()->num_ref_edges();//n = 21

  //test for multi-cut imprint for subtract.
  from_bodies.clean_out();
  from_bodies.append(from_body);
  new_bodies.clean_out();
  rsl = gmti->subtract(tool_body, from_bodies, new_bodies,
                       CUBIT_TRUE, CUBIT_FALSE);
  //Updated volume(s): 24
  //Destroyed volume(s): 25

  n = new_bodies.get()->num_ref_faces();
  //n = 8
  n = new_bodies.get()->num_ref_edges();
  //n = 18

  bodies.clean_out();
  gti->bodies(bodies);
  //delete all entities
  gti->delete_Body(bodies);

  free_entities.clean_out();
  gti->get_free_ref_entities(free_entities);
  assert(free_entities.size() == 0);
  //there shouldn't be any free_entites.

  //test for shell body subtract.
  tool_body = gmti->brick(1, 1, 1);
  //just need two surfaces.
  DLIList<RefFace*> reffaces;
  tool_body->ref_faces(reffaces);
  CubitVector v1(1, 0.5, 0.5);
  CubitVector v2(0.5, 0.5, 1);
  DLIList<Surface*> surfaces;
  for (int i = 0; i < tool_body->num_ref_faces(); i++)
  {
    RefFace* face = reffaces.get_and_step();
    CubitVector v = face->center_point();
    if (v.about_equal(v1) || v.about_equal(v2))
      surfaces.append(face->get_surface_ptr());
  } 
  assert(surfaces.size() == 2);
  DLIList<BodySM*> body_list;
  for (int i = 0; i < surfaces.size(); i++)
  {
    Surface* surface = surfaces.get_and_step();
    surface = ome->make_Surface(surface);
    body_list.append(CAST_TO(surface,OCCSurface)->my_body());
  } 
  bodysm = NULL;

  //test stitch surfaces operation
  ome->stitch_surfs(body_list, bodysm);

  bodies.clean_out();
  bodies.append(tool_body);
  gti->delete_Body(bodies);

  //test flip_normal for a shell body.
  from_body2 = gti->make_Body(bodysm); 
  test_body = gmti->copy_body(from_body2);
  ref_faces.clean_out();
  test_body->ref_faces(ref_faces);
  normal = ref_faces.get()->normal_at(v1); //(1,0,0)
  normal = ref_faces.step_and_get()->normal_at(v2); //(0,0,1)
  surfaces.clean_out();
  surfaces.append(ref_faces.step_and_get()->get_surface_ptr());
  ome->flip_normals(surfaces);
  normal = ref_faces.get()->normal_at(v1); //(-1,0,0)
  normal = ref_faces.step_and_get()->normal_at(v2); //(0,0,1)
  gti->delete_Body(test_body);
  surfaces.clean_out();
   

  tool_body  = gmti->brick(4, 4, 4);
  CubitVector v_move3(0,1,0);
  gti->translate(tool_body,v_move3);
  Body* copy_tool_body = gmti->copy_body(tool_body);
  Body* copy_tool_body2 = gmti->copy_body(tool_body);

  from_bodies.clean_out();
  new_bodies.clean_out();
  from_bodies.append(copy_tool_body);
  from_bodies.append(from_body2);
  stat = gmti->imprint(from_bodies, new_bodies, CUBIT_FALSE); 
  //Updated volume(s): 32 one body gets cut once with 7 faces, the other is not

  //test body imprinted by curves.
  ref_edges.clean_out();
  from_body2->ref_edges(ref_edges);
  from_bodies.clean_out();
  new_bodies.clean_out();
  from_bodies.append(copy_tool_body2);
  stat = gmti->imprint(from_bodies, ref_edges, new_bodies, CUBIT_FALSE, CUBIT_TRUE );
  //Updated volume(s): 33, no cuts performed

  //test body cutting a shell, one surface got cut as the result. 
  CubitVector v_move6(1,-1,0);
  gti->translate(tool_body,v_move6);
  from_bodies.clean_out();
  from_bodies.append(from_body2);
  new_bodies.clean_out();
  rsl = gmti->subtract(tool_body, from_bodies, new_bodies,
                       CUBIT_TRUE, CUBIT_TRUE);
  //Created volume(s): 34

  d = new_bodies.step_and_get()->measure();
  v = new_bodies.get()->center_point();
  n = new_bodies.get()->num_ref_faces();
  // n = 1
 
  from_bodies.clean_out();
  from_bodies.append(tool_body);
  
  //test a shell cutting a body, failed operation with an Error message.
  rsl = gmti->subtract(from_body2, from_bodies, new_bodies,
                       CUBIT_TRUE, CUBIT_TRUE);
  //WARNING: Surfaces or Shells can't be used to cut a body.
  //ERROR: Subtract FAILED

  //test solid solid imprint
  tool_body  = gmti->brick(4, 4, 4);
  CubitVector v_move5(0,0.5,0);
  gti->translate(tool_body,v_move5);
  from_body = gmti->brick(1,1,1);
  from_bodies.clean_out();
  new_bodies.clean_out();
  from_bodies.append(from_body);
  from_bodies.append(tool_body);
  stat =  gmti->imprint(from_bodies, new_bodies, CUBIT_FALSE);
  //Updated volume(s): 35, 36
  //one body gets 4 cuts, the other 3.
  //one body gets 10 faces, the other 9.
  
  //test imprint projected edges
  bodies.clean_out();
  gti->bodies(bodies);
  //delete all entities
  gti->delete_Body(bodies);

  free_entities.clean_out();
  gti->get_free_ref_entities(free_entities);  
  assert(free_entities.size() == 0);
 
  OCCQueryEngine* oqe = OCCQueryEngine::instance();
  DLIList <OCCBody* > *occ_bodies = oqe->BodyList;

  from_body  = gmti->brick(4, 4, 4);
  CubitVector v_move9(0,3,0);
  gti->translate(from_body,v_move9); 
  ref_faces.clean_out();
  from_body->ref_faces(ref_faces);
  tool_body = gmti->brick(11, 1, 1);  
  ref_edges.clean_out();
  tool_body->ref_edges(ref_edges);
  new_bodies.clean_out();
  gmti->imprint_projected_edges(ref_faces,ref_edges, new_bodies, CUBIT_TRUE,
       CUBIT_FALSE); 
  //Created volume(s): 39
  if(new_bodies.size())
    n = new_bodies.get()->num_ref_faces();//n = 8, new_bodies.size() == 1

  //delete all bodies.
  bodies.clean_out();
  gti->bodies(bodies);
  //delete all entities
  gti->delete_Body(bodies);

  free_entities.clean_out();
  gti->get_free_ref_entities(free_entities);
  assert(free_entities.size() == 0);
  //test body-body intersect.
  //1. from body is the commom body, no update
  tool_body  = gmti->brick(4, 4, 4);
  from_body = gmti->brick(1,1,1);
  from_bodies.clean_out();
  new_bodies.clean_out();
  from_bodies.append(from_body);
  stat =  gmti->intersect(tool_body, from_bodies, new_bodies, CUBIT_FALSE);
  //Updated volume(s): 41
  //Destroyed volume(s): 40
  d = new_bodies.get()->measure(); //d = 1
  
  //2. common body is part of from body, update the correponding face
  tool_body  = gmti->brick(4, 4, 4);
  CubitVector v_move7(0,0.5,0);
  gti->translate(tool_body,v_move7);
  from_bodies.clean_out();
  from_bodies.append(from_body);
  new_bodies.clean_out();
  stat =  gmti->intersect(tool_body, from_bodies, new_bodies, CUBIT_FALSE);
  d = new_bodies.get()->measure(); //d = 0.5
  //Updated volume(s): 41
  //Destroyed volume(s): 42

  //3. there's no common body, from body is deleted or kept depending on
  //keep-old flag. 
  tool_body  = gmti->brick(4, 4, 4);
  gti->translate(tool_body,v_move7);
  gti->translate(tool_body,v_move7);
  from_bodies.clean_out();
  new_bodies.clean_out();
  from_bodies.append(from_body);
  //nothing changed 
  stat =  gmti->intersect(tool_body, from_bodies, new_bodies, CUBIT_TRUE);
  //"The 1 body did not have common part with the tool_body."

  //from_body get deleted
  stat =  gmti->intersect(tool_body, from_bodies, new_bodies, CUBIT_FALSE);
  //"The 1 body did not have common part with the tool_body."
  //Destroyed volume(s): 41, 43

  bodies.clean_out();
  gti->bodies(bodies);
  //delete all entities
  gti->delete_Body(bodies);

  free_entities.clean_out();
  gti->get_free_ref_entities(free_entities);
  assert(free_entities.size() == 0);

  //test chop operation
  tool_body  = gmti->brick(4, 4, 4);
  from_body = gmti->brick(1,1,1);
  gti->translate(tool_body,v_move7);
  from_bodies.clean_out();
  from_bodies.append(from_body);
  from_bodies.append(tool_body);
  new_bodies.clean_out();

  //new_bodies = intersect bodies;bodies = outside bodies; body is dummy 
  stat =  gmti->chop(from_bodies, new_bodies, bodies, body, CUBIT_FALSE);
  d = new_bodies.get()->measure(); //d = 0.5
  //Created volume(s): 46
  //Destroyed volume(s): 44
  //Updated volume(s): 45 

  bodies.clean_out();
  gti->bodies(bodies);
  //delete all entities
  gti->delete_Body(bodies);

  free_entities.clean_out();
  gti->get_free_ref_entities(free_entities);
  assert(free_entities.size() == 0);
  //test chop 2
  from_body = gmti->brick(4, 4, 4);
  tool_body = gmti->brick(1,1,1);
  gti->translate(tool_body,v_move7);
  from_bodies.clean_out();
  from_bodies.append(from_body);
  from_bodies.append(tool_body);
  new_bodies.clean_out();
  bodies.clean_out();
  stat =  gmti->chop(from_bodies, new_bodies, bodies, body, CUBIT_FALSE);
  //Created volume(s): 49
  //Destroyed volume(s): 48
  //Updated volume(s): 47
  d = bodies.get()->measure();//d = 63
  
  bodies.clean_out();
  gti->bodies(bodies);
  //delete all entities
  gti->delete_Body(bodies);

  free_entities.clean_out();
  gti->get_free_ref_entities(free_entities);
  assert(free_entities.size() == 0);

  //test unite 1
  tool_body = gmti->brick(1,1,1);
  gti->translate(tool_body,v_move7);
  gti->translate(tool_body,v_move7);
  from_body = gmti->brick(1,1,1);
  from_bodies.append(from_body);
  from_bodies.append(tool_body);
  new_bodies.clean_out();
  stat = gmti->unite(from_bodies, new_bodies, CUBIT_FALSE);
  //Updated volume(s): 51
  //Destroyed volume(s): 50
  d = new_bodies.get()->measure(); //d = 2
  n = new_bodies.get()->num_ref_faces(); //d = 10  
 
  bodies.clean_out();
  gti->bodies(bodies); //bodies.size() = 1
  free_entities.clean_out();
  gti->get_free_ref_entities(free_entities); //free_entities.size() = 0
  assert(free_entities.size() == 0);
  //delete all entities
  gti->delete_Body(bodies);
  free_entities.clean_out();
  gti->get_free_ref_entities(free_entities); //free_entities.size() = 0
  assert(free_entities.size() == 0);
  //test unite 2
  tool_body = gmti->brick(4, 4,4);
  gti->translate(tool_body,v_move7);
  from_body = gmti->brick(1,1,1);
  from_bodies.clean_out();
  from_bodies.append(from_body);
  from_bodies.append(tool_body);
  new_bodies.clean_out();
  stat = gmti->unite(from_bodies, new_bodies, CUBIT_FALSE);
  //Updated volume(s): 53
  //Destroyed volume(s): 52
  d = new_bodies.get()->measure(); //d = 64.5
  n = new_bodies.get()->num_ref_faces(); 
  //d = 13. because fusion keeps the imprintings.

  bodies.clean_out();
  gti->bodies(bodies); //bodies.size() = 1
  free_entities.clean_out();
  gti->get_free_ref_entities(free_entities); //free_entities.size() = 0
  assert(free_entities.size() == 0);
  //delete all entities
  gti->delete_Body(bodies);
  free_entities.clean_out();
  gti->get_free_ref_entities(free_entities); //free_entities.size() = 0
  assert(free_entities.size() == 0);

  //test making thick body.
  from_body = gmti->cylinder(10, 4, 4, 4); 
  tool_body = gmti->cylinder(5,1,1,1);
  //  v = from_body->center_point();
  //v = tool_body->center_point();
  CubitVector v_move8(0,0,10);
  gti->translate(tool_body,v_move8);  
  from_bodies.clean_out();
  from_bodies.append(from_body);
  from_bodies.append(tool_body);
  new_bodies.clean_out();
  stat = gmti->unite(from_bodies, new_bodies, CUBIT_FALSE);
  //Updated volume(s): 54
  //Destroyed volume(s): 55
  d = new_bodies.get()->measure(); //d = 518.3627
  n = new_bodies.get()->num_ref_faces(); //n = 5
  //find the top most surface as the opening of the thick body.
  ref_faces.clean_out();
  new_bodies.get()->ref_faces(ref_faces);
  CubitVector center(0,0,15);
  for(int i = 0; i < n; i++)
  {
    if(ref_faces.step_and_get()->is_planar() && 
       ref_faces.get()->center_point() == center )
      break;
  }
  RefFace* sweep_face = gmti->make_RefFace(ref_faces.get());
  DLIList<RefFace*> faces_to_remove;
  faces_to_remove.append(ref_faces.get());
  from_bodies = new_bodies;
  new_bodies.clean_out();
  stat = gmti->hollow(from_bodies, faces_to_remove, new_bodies, -.2);
  //Updated volume(s): 54
  n = new_bodies.get()->num_ref_faces(); //n = 10
  d = new_bodies.get()->measure(); //d = 72.3618

  DLIList<RefEntity*> refentities;
  refentities.append(sweep_face);
  RefFace* draft_face = gmti->make_RefFace(sweep_face);
  RefFace* perp_face = gmti->make_RefFace(sweep_face);
  RefFace* rotate_face = gmti->make_RefFace(sweep_face);
  RefFace* rotate_face2 = gmti->make_RefFace(sweep_face);
  gmti->sweep_translational(refentities, v_move8, 0, 1, CUBIT_FALSE, CUBIT_FALSE);  
  body = CAST_TO(refentities.get(), Body);
  d = body->measure();
  //d = 31.4159
  refentities.clean_out();
  refentities.append(draft_face);
  gmti->sweep_translational(refentities, v_move8, 0.087, 1, CUBIT_FALSE, CUBIT_FALSE); 
  body = CAST_TO(refentities.get(), Body);
  d = body->measure();
  //d = 66.3676  theoretical calculation is 66.7833, error 0.62%
  DLIList<RefEdge*> edges;
  body->ref_edges(edges);
  refentities.clean_out();
  refentities.append(edges.get());
  gmti->sweep_translational(refentities, v_move8, 0.087, 1, CUBIT_FALSE, CUBIT_FALSE);
  body = CAST_TO(refentities.get(), Body);
  d = body->measure();
  //d = area = 90.1292 theoretica calculation is 90.5754, error 0.49%
  refentities.clean_out();
  refentities.append(perp_face);
  gmti->sweep_perpendicular(refentities, 10, 0.087, 1, CUBIT_FALSE, CUBIT_FALSE);
  body = CAST_TO(refentities.get(), Body);
  d = body->measure();
  //d = 66.3676  theoretical calculation is 66.7833, error 0.62%
  
  //Testing for sweep_rotational function
  //1. Negative testing: surface rotates along an intersection axis, fails.
  refentities.clean_out();
  refentities.append(rotate_face); 
  DLIList<RefEdge*> rotate_edges;
  rotate_face->ref_edges(rotate_edges);
  CubitVector sweep_axis(1,0,0);   
  CubitVector point(0, 2, 15);
  double angle = 1.57;
  gmti->sweep_rotational(refentities, center, sweep_axis, angle, 0, 0, 1,
                         CUBIT_FALSE, CUBIT_TRUE, CUBIT_FALSE);

  //2. surface rotates along a non-intersection axis
  refentities.append(rotate_face);
  gmti->sweep_rotational(refentities, point, sweep_axis, angle, 0, 0, 1,
                         CUBIT_FALSE, CUBIT_TRUE, CUBIT_FALSE);
  body = CAST_TO(refentities.get(), Body);
  d = body->measure();
  //d = 9.8646

  //3. closed curve rotates along an intersecting axis, fails
  refentities.clean_out();
  refentities.append(rotate_edges.get());
  gmti->sweep_rotational(refentities, center, sweep_axis, angle, 0, 0, 1,
                         CUBIT_FALSE, CUBIT_TRUE, CUBIT_FALSE); 

  //4. closed curve rotates along a non-intersecting axis with make_solid
  //   option set to be true.
  refentities.append(rotate_edges.get());
  gmti->sweep_rotational(refentities, point, sweep_axis, angle, 0, 0, 1,
                         CUBIT_FALSE, CUBIT_TRUE, CUBIT_FALSE);
  body = CAST_TO(refentities.get(), Body);
  d = body->measure();
  //d = 9.8646

  //5. closed curve rotates along a non-intersecting axis with make_solid
  //   option set to be false.
  refentities.clean_out();
  refentities.append(rotate_edges.get());
  gmti->sweep_rotational(refentities, point, sweep_axis, -angle, 0, 0, 1,
                         CUBIT_FALSE, CUBIT_FALSE, CUBIT_FALSE);
  body = CAST_TO(refentities.get(), Body);
  d = body->measure();
  //d = 19.7292

  //6. open curve rotates along an intersecting axis with make_solid
  //   option set to be true. Failed
  CubitVector pt1(0,1,15);
  CubitVector pt2(0,-1,15);
  CubitVector pt3(-1,0,15);
  RefVertex* vt1 = gmti->make_RefVertex(pt1);
  RefVertex* vt2 = gmti->make_RefVertex(pt2);
  RefVertex* vt3 = gmti->make_RefVertex(pt3);
  RefEdge* edge1 = gmti->make_RefEdge(vt1, vt2, rotate_face, vt3);
  CubitVector apoint(0,-0.5,15);
  refentities.clean_out();
  refentities.append(edge1);
  gmti->sweep_rotational(refentities, apoint, sweep_axis, angle, 0, 0, 1,
                         CUBIT_FALSE, CUBIT_TRUE, CUBIT_FALSE);

  //7. open curve rotates along an intersecting axis at end points with 
  //   make_solid option set to be true.
  refentities.append(edge1);
  CubitVector rotate_axis(0, 1, 0);
  //  gmti->sweep_rotational(refentities, center, rotate_axis, angle, 0, 0, 1,
                         CUBIT_FALSE, CUBIT_TRUE, CUBIT_FALSE);
    //body = CAST_TO(refentities.get(), Body);
    //d = body->measure();
  //d = 19.7292

  return CUBIT_SUCCESS;
}
