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

// forward declare some functions used and defined later
CubitStatus read_geometry(int, char **);
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
CubitStatus read_geometry(int num_files, char **argv) 
{
  CubitStatus status = CUBIT_SUCCESS;
  GeometryQueryTool *gti = GeometryQueryTool::instance();
  assert(gti);
  int i;
  
  PRINT_SEPARATOR;

  for (i = 0; i < num_files; i++) {
    status = gti->import_solid_model(argv[i], "OCC");
    if (status != CUBIT_SUCCESS) {
      PRINT_ERROR("Problems reading geometry file %s.\n", argv[i]);
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
  BodySM* bodysm = body->get_body_sm_ptr();
  DLIList<OCCSurface*> occ_surfaces;
  CAST_TO(bodysm, OCCBody)->get_all_surfaces(occ_surfaces);  
  DLIList<RefFace*> ref_faces;
  DLIList<RefEdge*> ref_edges;
  body->ref_faces(ref_faces);

  DLIList<RefFace*> faces_to_stitch;
  for(int i = 0 ; i < ref_faces.size(); i++)
  {
    RefFace* refface = gmti->make_RefFace(ref_faces.get_and_step());
    if(refface)
	faces_to_stitch.append(refface);
  }

  gti->delete_Body(body);

  DLIList<BodySM*> bodysm_list;
  DLIList<RefFace*> face_list;
  DLIList<Surface*> surface_list;
  DLIList<RefVertex*> vertices;
  CubitVector v(15,0,0);
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
  //ome->stitch_surfs(bodysm_list, bodysm);
  //Lump* lump = ome->make_Lump(surface_list);
  //bodysm = CAST_TO(lump, OCCLump)->get_body();
  //gti->make_Body(bodysm);

  CubitStatus rsl = CUBIT_SUCCESS;
  DLIList<RefEntity*> ref_entity_list;
  int num_ents_exported=0;
  const CubitString cubit_version="10.2";
  const char * filename = "stitch.occ";
  const char * filetype = "OCC";

  rsl = gti->export_solid_model(ref_entity_list, filename, filetype,
                                 num_ents_exported, cubit_version);

  DLIList<Body*> bodies;
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
  char *argv = "./stitch.occ";
  CubitStatus status = read_geometry(1, &argv);
  if (status == CUBIT_FAILURE) exit(1);

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
  //new bodies has one body, new body has 10 ref-faces, 5 of them are remaining
  //with old id, 5 of them are new faces.
  from_bodies=new_bodies;
  new_bodies.clean_out();
  rsl = gmti->subtract(tool_body,from_bodies, new_bodies,
                       CUBIT_TRUE, CUBIT_FALSE);
  d = new_bodies.step_and_get()->measure();
  v = new_bodies.get()->center_point();
  int n = new_bodies.get()->num_ref_faces();
  // n = 6
  //new bodies has 2 bodies, one has a volume = 10 and the other has a 
  //volume = 50; each of them has 6 ref_faces, of which 3 are new and 3 are
  //remaining (unchanged or modified).

/*
  from_bodies.clean_out();
  from_bodies += new_bodies;
  new_bodies.clean_out();
  CubitVector v_move2(0, -2, -2);
  Body* body_new = from_bodies.step_and_get();
  d = body_new->measure();
  v = body_new->center_point();
  gti->translate(body_new,v_move2);
  rsl = gmti->subtract(tool_body, from_bodies, new_bodies,
                       CUBIT_TRUE, CUBIT_FALSE);
  d = new_bodies.step_and_get()->measure();
  n = new_bodies.get()->num_ref_faces();
  // n = 8
  n = new_bodies.get()->num_ref_edges();
  // n = 22

  bodies.clean_out();
  gti->bodies(bodies);
  //delete all entities
  gti->delete_Body(bodies); 
  
  //test for multi-cut imprint for subtract.
  from_body = gmti->brick(10, 10, 10);
  tool_body = gmti->brick(11, 1, 1);
  CubitVector v_move4(0,1,-1);
  gti->translate(from_body,v_move4);
  from_bodies.clean_out();
  from_bodies.append(from_body);
  new_bodies.clean_out();
  rsl = gmti->subtract(tool_body, from_bodies, new_bodies,
                       CUBIT_TRUE, CUBIT_TRUE); 
  n = new_bodies.get()->num_ref_faces();
  //n = 8
  n = new_bodies.get()->num_ref_edges();
  //n = 18
*/
  bodies.clean_out();
  gti->bodies(bodies);
  //delete all entities
  gti->delete_Body(bodies);

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

  from_body2 = gti->make_Body(bodysm); 

  tool_body  = gmti->brick(4, 4, 4);
  BodySM* copy_bodysm = ome->copy_body(tool_body->get_body_sm_ptr());
  CubitVector v_move3(0,1,0);
  gti->translate(tool_body,v_move3);
/*  from_bodies.clean_out();
  from_bodies.append(from_body2);
  new_bodies.clean_out();

  //test face body imprint
  TopoDS_Shape* tool_shape = CAST_TO(copy_bodysm,OCCBody)->get_TopoDS_Shape();  
  TopoDS_Shape* from_shape = CAST_TO(body_list[0],OCCBody)->my_sheet_surface()->get_TopoDS_Face(); 
  ome->imprint_toposhapes(tool_shape, from_shape);

  //test shell body imprint.
  from_shape = CAST_TO(bodysm,OCCBody)->shell()->get_TopoDS_Shell();
  ome->imprint_toposhapes(tool_shape, from_shape);

  //test body cutting a shell, one surface got cut as the result. 
  CubitVector v_move6(1,-1,0);
  gti->translate(tool_body,v_move6);
  rsl = gmti->subtract(tool_body, from_bodies, new_bodies,
                       CUBIT_TRUE, CUBIT_TRUE);
  d = new_bodies.step_and_get()->measure();
  v = new_bodies.get()->center_point();
  n = new_bodies.get()->num_ref_faces();
  // n = 1
 
  from_bodies.clean_out();
  from_bodies.append(tool_body);
  
  //test a shell cutting a body, failed operation with a warning message.
  rsl = gmti->subtract(from_body2, from_bodies, new_bodies,
                       CUBIT_TRUE, CUBIT_TRUE);

  //test solid solid imprint
  tool_body  = gmti->brick(4, 4, 4);
  CubitVector v_move5(0,0.5,0);
  gti->translate(tool_body,v_move5);
  from_body = gmti->brick(1,1,1);
  from_shape = CAST_TO(from_body->get_body_sm_ptr(), OCCBody)->get_TopoDS_Shape();
  tool_shape = CAST_TO(tool_body->get_body_sm_ptr(),OCCBody)->get_TopoDS_Shape();
  ome->imprint_toposhapes(tool_shape, from_shape);
  ome->imprint_toposhapes(from_shape, tool_shape);
  return stat;
*/
}
