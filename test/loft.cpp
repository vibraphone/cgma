/**
 * \file makept.cpp
 *
 * \brief makept, another simple C++ driver for CGM
 *
 * This program acts as a simple driver for CGM.  It reads in a geometry,
 * and performs varies checks for bodies, surfaces, curves and vertices.
 */

#undef NDEBUG
#include <cassert>

#include "GeometryModifyEngine.hpp"
#include "GeometryModifyTool.hpp"
#include "GeometryQueryTool.hpp"
#include "CubitMessage.hpp"
#include "Body.hpp"
#include "RefVolume.hpp"
#include "RefFace.hpp"
#include "RefEdge.hpp"
#include "RefVertex.hpp"
#include "InitCGMA.hpp"
#include "OCCBody.hpp"
#include "OCCSurface.hpp"
#include "OCCCurve.hpp"
#include "OCCDrawTool.hpp"
#include "CubitCompat.hpp"

#ifndef SRCDIR
# define SRCDIR .
#endif

#define STRINGIFY_(X) #X
#define STRINGIFY(X) STRINGIFY_(X)
#define SRCPATH STRINGIFY(SRCDIR) "/"

// forward declare some functions used and defined later
CubitStatus read_geometry(int, const char **, bool local = false);
CubitStatus make_Point();
// macro for printing a separator line
#define PRINT_SEPARATOR   PRINT_INFO("=======================================\n");


// main program - initialize, then send to proper function
int main (int argc, char **argv)
{
  CubitStatus status = InitCGMA::initialize_cgma("OCC");
  if (CUBIT_SUCCESS != status) return 1;

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
CubitStatus read_geometry(int num_files, const char **argv, bool local) 
{
  CubitStatus status = CUBIT_SUCCESS;
  GeometryQueryTool *gti = GeometryQueryTool::instance();
  assert(gti);
  int i;
  
  PRINT_SEPARATOR;

  for (i = 0; i < num_files; i++) {
    std::string filename( local ? "./" : SRCPATH );
    filename += argv[i];
    status = CubitCompat_import_solid_model(filename.c_str(), "OCC");
    if (status != CUBIT_SUCCESS) {
      PRINT_ERROR("Problems reading geometry file %s.\n", filename.c_str());
    }
  }
  PRINT_SEPARATOR;

  return CUBIT_SUCCESS;
}

CubitStatus make_Point()
{
  GeometryQueryTool *gti = GeometryQueryTool::instance();
  GeometryModifyTool *gmti = GeometryModifyTool::instance();

  //test for loft two surfaces:
  Body* body = gmti->brick(10, 10, 10);
  CubitVector vector1(23,13,11);
  gti->translate(body, vector1);
  Body* body2 = gmti->brick(10, 10, 10);

  //use the surface1 and surface8 to make a loft
  DLIList<RefFace*> face_list;
  body->ref_faces(face_list);
  DLIList<RefFace*> loft_surfaces;
  face_list.last();
  loft_surfaces.append(face_list.get());
  face_list.clean_out();
  body2->ref_faces(face_list);
  face_list.last();
  face_list.back();
  loft_surfaces.append(face_list.get());

  //Do loft
  Body* new_body = NULL;
  gmti->loft_surfaces_to_body(loft_surfaces.step_and_get(), 0.0, 
                           loft_surfaces.step_and_get(), 0.0, new_body, CUBIT_FALSE,
                           CUBIT_FALSE, CUBIT_FALSE, CUBIT_FALSE, CUBIT_FALSE); 
  double volume = new_body->measure();
  assert(fabs(volume - 2100) < 0.000001);

  DLIList<RefEntity*> ref_entity_list;
  const CubitString cubit_version="12.2";
  int num_ents_exported=0;
  const char *filename = "loft.brep";
  const char *filetype = "OCC";
  CubitStatus rsl = CubitCompat_export_solid_model(ref_entity_list, filename, filetype,
                                 num_ents_exported, cubit_version);

  DLIList<Body*> bodies;
  gti->bodies(bodies);

  //delete all entities
  gti->delete_Body(bodies);
  
  CubitVector p1(0, 0 ,0);
  CubitVector p2(1, 0 ,0);
  CubitVector p3(0, 1 ,0);
  CubitVector p4(1, 1 ,0);

  gmti->make_RefVertex(p1,5);
  gmti->make_RefVertex(p2,5);
  gmti->make_RefVertex(p3,5);
  gmti->make_RefVertex(p4,5);
  DLIList<RefEntity*>  free_entities;
  gti->get_free_ref_entities(free_entities);

  RefVertex* vertex1 = CAST_TO(free_entities.step_and_get(),RefVertex);
  RefVertex* vertex2 = CAST_TO(free_entities.step_and_get(),RefVertex);
  RefVertex* vertex3 = CAST_TO(free_entities.step_and_get(),RefVertex);
  RefVertex* vertex4 = CAST_TO(free_entities.step_and_get(),RefVertex);
  RefEdge* new_edge1 = gmti->make_RefEdge(STRAIGHT_CURVE_TYPE, vertex1, vertex2);
  RefEdge* new_edge2 = gmti->make_RefEdge(STRAIGHT_CURVE_TYPE, vertex3, vertex4);
  DLIList<RefEdge*> edges;
  edges.append(new_edge1);
  edges.append(new_edge2);
  DLIList<RefEdge*> guides;
  gmti->create_skin_surface(edges, new_body, guides); 
  CubitVector center = new_body->center_point();
  CubitVector comp (0.5, 0.5, 0);
  assert(center.distance_between(comp) < 0.00001);
  CubitBoolean is_sheet = new_body->is_sheet_body();
  assert(is_sheet == CUBIT_TRUE);  
  return CUBIT_SUCCESS;
}

