/**
 * \file makept.cpp
 *
 * \brief makept, another simple C++ driver for CGM
 *
 * This program acts as a simple driver for CGM.  It reads in a geometry,
 * performs imprints between all the bodies, merges them, and writes information
 * on the results.  It also performs pairwise intersections between the
 * bodies to check for overlaps.  Results are written to stardard output.
 *
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
#include "OCCBody.hpp"
#include "OCCSurface.hpp"

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
  GeometryQueryTool *gti = GeometryQueryTool::instance();
  GeometryModifyTool *gmti = GeometryModifyTool::instance();

  OCCQueryEngine::instance();
  OCCModifyEngine::instance();

  // Read in the geometry from files specified on the command line
  char *argv = "./66_shaver3.brep";
  CubitStatus status = read_geometry(1, &argv);
  if (status == CUBIT_FAILURE) exit(1);

  argv = "./62_shaver1.brep";
  status = read_geometry(1, &argv);
  if (status == CUBIT_FAILURE) exit(1);

  argv = "./72_shaver6.brep";
  status = read_geometry(1, &argv);
  if (status == CUBIT_FAILURE) exit(1);
  
  CubitVector vector1(10,10,10);
  CubitVector vector2(15,15,15);
  DLIList<RefEntity*> free_entities;
  gti->get_free_ref_entities(free_entities);
 
  for(int i = 1; i <= free_entities.size(); i++)
  {
     RefEntity * entity = free_entities.get_and_step();
     gti->translate((BasicTopologyEntity*)entity, i*vector1); 
  }

  // Read in the geometry from files specified on the command line
  gmti->make_RefVertex(vector1,5);
  gmti->make_RefVertex(vector2,5);

  CubitStatus rsl = CUBIT_SUCCESS;
  DLIList<RefEntity*> ref_entity_list;
  int num_ents_exported=0;
  const CubitString cubit_version="10.2";
  const char * filename = "point.occ";
  const char * filetype = "OCC";
  /*
  rsl = gti->export_solid_model(ref_entity_list, filename, filetype, 
                                 num_ents_exported, cubit_version);

  */
  //delete all entities
  DLIList<Body*> bodies;
  gti->bodies(bodies);
  gti->get_free_ref_entities(free_entities);
 
  RefVertex* vertex1 = CAST_TO(free_entities.get_and_step(),RefVertex);
  RefVertex* vertex2 = CAST_TO(free_entities.get(),RefVertex);
  CubitBoolean is_equal = gti->
		about_spatially_equal(vertex1,vertex2);
 
  CubitVector vi, vii;
  double d;
  gti->entity_entity_distance(vertex1,vertex2,vi, vii,d);
 
  d = bodies.get()->measure(); 
  vi = bodies.get()->center_point();
  
  CubitBox box = bodies.get()->bounding_box();

  gti->entity_entity_distance(gti->get_first_ref_volume(), vertex2,vi, vii,d);

  BodySM* body = bodies.get()->get_body_sm_ptr();
  CubitVector axis(10,0,0);
  gti->reflect(bodies, axis);
  vi = bodies.get()->center_point();
  gti->scale(bodies.get(),2);
  vi = bodies.get()->center_point();
  gti->translate(bodies.get(),axis);
  vi = bodies.get()->center_point();
  gti->rotate(bodies.get(), axis, 30);
  vi = bodies.get()->center_point();
  OCCBody* occ_body = CAST_TO(body, OCCBody);
  occ_body->mass_properties(vi, d);
  vi = occ_body->get_bounding_box().center(); 

  //check for surfaces
  DLIList<OCCSurface*> surfaces;

  CAST_TO(body, OCCBody)->get_all_surfaces(surfaces);
  OCCSurface* surface = surfaces.get();
  GeometryType type = surface->geometry_type();
  box = surface->bounding_box();
  vi = surface->center_point();
  CubitVector normal;
  surface->get_point_normal(vi, normal);
  CubitVector* closest_location = new CubitVector;
  CubitVector* unit_normal_ptr = new CubitVector;
  CubitVector* curvature1_ptr = new CubitVector;
  CubitVector* curvature2_ptr = new CubitVector;
  surface->closest_point(vi, closest_location, unit_normal_ptr, curvature1_ptr,
                         curvature2_ptr);

  double u = 0.5;      
  double v = 0.5;
  vi = surface->position_from_u_v(u,v);
  CubitBoolean periodic = surface->is_periodic();
  double p = 0; //period
  periodic = surface->is_periodic_in_U(p); 
  periodic = surface->is_periodic_in_V(p);
  CubitBoolean sigular = surface->is_singular_in_U(1);
  double lower, upper;
  surface->get_param_range_U(lower, upper);
  sigular = surface->is_singular_in_U(upper); 
  sigular = surface->is_singular_in_U(lower);
   
  surface->get_param_range_V(lower, upper);
  sigular = surface->is_singular_in_V(upper);
  sigular = surface->is_singular_in_V(lower);

  CubitBoolean closed = surface->is_closed_in_U();
  closed = surface->is_closed_in_V();

  CubitPointContainment pc = surface->point_containment(1,1);
  gti->delete_Body(bodies);

  for (int j = free_entities.size(); j--;)
    {
      gti->delete_RefEntity( free_entities.get_and_step());
    }

  return rsl;
}
