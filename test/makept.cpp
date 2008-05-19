/**
 * \file makept.cpp
 *
 * \brief makept, another simple C++ driver for CGM
 *
 * This program acts as a simple driver for CGM.  It reads in a geometry,
 * and performs varies checks for bodies, surfaces, curves and vertices.
 */
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
#include "OCCCurve.hpp"

#ifndef SRCDIR
# define SRCDIR .
#endif

#define STRINGIFY_(X) #X
#define STRINGIFY(X) STRINGIFY_(X)
#define SRCPATH STRINGIFY(SRCDIR) "/"

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
    std::string filename( SRCPATH );
    filename += argv[i];
    status = gti->import_solid_model(filename.c_str(), "OCC");
    if (status != CUBIT_SUCCESS) {
      PRINT_ERROR("Problems reading geometry file %s.\n", filename);
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
  
  // test create a Compound body.
  DLIList<Body*> test_bodies;
  gti->bodies(test_bodies);
 
  DLIList<RefVolume*> ref_volume_list;
  for(int i = 0; i < test_bodies.size(); i++)
    test_bodies.get_and_step()->ref_volumes(ref_volume_list);
  
  Body* CompBody = gmti->make_Body(ref_volume_list);

  BodySM* CompBodySM = CompBody->get_body_sm_ptr();

  OCCBody *occ_CompBody = CAST_TO(CompBodySM, OCCBody);

  test_bodies.clean_out();
  gti->bodies(test_bodies);

  CubitVector vi, vii;
  vi = CompBody->center_point(); 

  CubitVector axis(10,0,0);

  gti->translate(CompBody,axis);
  vi = CompBody->center_point();
  // After parellel move, center point moved by x (10)

  CubitVector vector1(10,10,10);
  CubitVector vector2(10,-10,10);
  DLIList<RefEntity*> free_entities;

  // Make two vertices.
  gmti->make_RefVertex(vector1,5);
  gmti->make_RefVertex(vector2,5);

  gti->get_free_ref_entities(free_entities);

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
 
  //check for vertex
  DLIList<Body*> bodies;
  gti->bodies(bodies);
  free_entities.clean_out();// get_free_ref_entities directly append
  //without checking for duplicates, so clean_out first. 
  gti->get_free_ref_entities(free_entities);
 
  RefVertex* vertex1 = CAST_TO(free_entities.get_and_step(),RefVertex);
  RefVertex* vertex2 = CAST_TO(free_entities.get(),RefVertex);
  CubitBoolean is_equal = gti->
		about_spatially_equal(vertex1,vertex2);
  //vertex1,vertex2 are not spatially equal. 
  
  double d;
  gti->entity_entity_distance(vertex1,vertex2,vi, vii,d);
  // distance (d) between vertex1,vertex2. vi (20, 20, 20) is vertex1 
  //translated by (10,10,10) and vii(30, 10, 30) is vertex2 translated
  //by 2*(10,10,10).
 
  //check for body
  d = bodies.get()->measure(); 
  // first body's volume.

  vi = bodies.get()->center_point();
  //first body's bounding box's center point.
  
  CubitBox box = bodies.get()->bounding_box();
  //first body's bounding box.

  gti->entity_entity_distance(gti->get_first_ref_volume(), vertex2,vi, vii,d);
  //first body and vertex2 's minimum distance(d) and locations for the minimum.

  BodySM* body = bodies.get()->get_body_sm_ptr();
  OCCBody* occ_body = CAST_TO(body, OCCBody);

  gti->reflect(bodies, axis);
  vi = bodies.get()->center_point();
  // After reflection, only x value should change.

  gti->scale(bodies.get(),2);
  vi = bodies.get()->center_point();
  // After scale, center point moved by 2 times 

  gti->translate(bodies.get(),axis);
  vi = bodies.get()->center_point();
  // After parellel move, center point moved by x (10)

  gti->rotate(bodies.get(), axis, 3.14/6);
  vi = bodies.get()->center_point();
  // After rotation, center point changed in y and z value.

  occ_body->mass_properties(vi, d);
  //true center and volume, volume should be 8 times of the original one.

  vi = occ_body->get_bounding_box().center(); 
  // bounding box center.

  //check for surface
  DLIList<OCCSurface*> surfaces;
  CAST_TO(body, OCCBody)->get_all_surfaces(surfaces);
  OCCSurface* surface = surfaces.step_and_get();
  GeometryType type = surface->geometry_type();
  // CONE_SURFACE_TYPE

  box = surface->bounding_box();
  // bounding box

  DLIList<RefFace*> ref_faces;
  gti->ref_faces(ref_faces);
  //RefFace* ref_face = ref_faces.step_and_get();
  RefFace* ref_face = ref_faces.get();

  //make a new refface out of existing refface.
  CubitBoolean extended_from = CUBIT_TRUE;
  RefFace* new_face = gmti->make_RefFace(ref_face, extended_from);

  rsl = gti->export_solid_model(ref_entity_list, filename, filetype,
                                 num_ents_exported, cubit_version);

  DLIList<DLIList<RefEdge*>*> ref_edge_loops;
  new_face->ref_edge_loops(ref_edge_loops);

  DLIList<RefEdge*>* ref_edge_list;
  ref_edge_list = ref_edge_loops.get();

  RefVertex* start = NULL;
  RefVertex* end = NULL;
  for (int i = 0; i < ref_edge_list->size(); i++)
  {
    RefEdge * edge = ref_edge_list->get_and_step();
    double d = edge->measure();
    start = edge->start_vertex();
    end = edge->end_vertex();
  }

  RefFace* new_face2 = gmti->make_RefFace(PLANE_SURFACE_TYPE, 
                         *ref_edge_list, ref_face, CUBIT_TRUE);

  bodies.clean_out();
  gti->bodies(bodies);
  //translate the new face by (40,40,40)
  for(int i = 1; i <= bodies.size(); i++)
  {
     bodies.step();
     if( i != 4)
        continue;
     Body * entity = bodies.get();
     gti->translate(entity, i*vector1);
  }

  RefVolume* volume = NULL;
  if ( new_face->get_surface_ptr()->is_closed_in_U())
    volume = new_face->ref_volume(); 
  else
    vi = new_face->center_point();
  //center point should moved by (40,40,40) compared with the original one below

  vi = ref_face->center_point();
  // center point

  CubitVector normal;
  normal = ref_face->normal_at(vi);
  // surface normal at center point.

  CubitVector closest_location ;
  CubitVector unit_normal_ptr ;
  CubitVector curvature1_ptr ;
  CubitVector curvature2_ptr ; 

  ref_face->find_closest_point_trimmed(vi, closest_location);  
  // Found surface projection location for vi.

  double curvature1, curvature2;
  ref_face->get_principal_curvatures(closest_location, curvature1, curvature2);
  // get principal curvatures at center point.

  double area = ref_face->area();
  // area of the face

  double u = 2.5;      
  double v = -20000;
  vi = ref_face->position_from_u_v(u,v);
  // get location on surface for it's u,v

  ref_face->u_v_from_position(vi, u, v);
  // get (u,v) from this vi.

  CubitBoolean periodic = ref_face->is_periodic();
  // found if surface is periodic.

  double p = 0; //period
  periodic = ref_face->is_periodic_in_U(p); 
  // found if surface is periodic in U and its period.

  periodic = ref_face->is_periodic_in_V(p);
  // found if surface is periodic in V and its period.

  //All OCC entities are parametric.

  double lower, upper;
  ref_face->get_param_range_U(lower, upper);
  // get surface U direction boundaries. here it's (2.48, 3.91)

  ref_face->get_param_range_V(lower, upper);
  // get surface V direction boundaries. here it's (-20010,-19998)

  CubitBoolean closed = surface->is_closed_in_U();
  // check if surface is closed in U direction.

  closed = surface->is_closed_in_V();
  // check if surface is closed in V direction.

  CubitPointContainment pc = ref_face->point_containment(1,-20000);
  // this (u,v) location should be outside of the surface.

  CubitPointContainment pc2 = ref_face->point_containment(3,-20000);
  // this (u,v) location should be inside of the surface.

  ref_edge_loops.clean_out();
  int num_loops = ref_face->ref_edge_loops(ref_edge_loops);
  DLIList<RefEdge*> *ref_edges1;
  ref_edges1 = ref_edge_loops.get();
  RefEdge* edge1 = ref_edges1->get();
  RefEdge* edge2 = ref_edges1->step_and_get();
  double angle = edge1->angle_between(edge2, ref_face);

  //test for curve
  CubitVector c_point, tangent, center;

  //make all kinds of curves.
  CubitVector center_pnt(0,0,0);
  RefEdge* new_edge_1 = gmti->make_RefEdge(SPLINE_CURVE_TYPE, vertex1,
					vertex2, &center_pnt);
  //Gives invalid curve type error.

  //straight line
  RefEdge* new_edge_2 = gmti->make_RefEdge(STRAIGHT_CURVE_TYPE, vertex1,
                                        vertex2, &center_pnt);
  d = new_edge_2->measure();
  new_edge_2->closest_point_trimmed(vi, c_point);

  //arc curve
  RefEdge* new_edge_3 = gmti->make_RefEdge(ARC_CURVE_TYPE, vertex1,
                                        vertex2, &center_pnt);
  d = new_edge_3->measure();
  new_edge_3->closest_point_trimmed(vi, c_point);

  //ellipse curve
  RefEdge* new_edge_4 = gmti->make_RefEdge(ELLIPSE_CURVE_TYPE, vertex1,
                                        vertex2, &center_pnt);
  d = new_edge_4->measure();
  new_edge_4->closest_point_trimmed(vi, c_point);

  RefEdge* new_edge_5 = gmti->make_RefEdge(ELLIPSE_CURVE_TYPE, vertex1,
                                        vertex2, &center_pnt, CUBIT_REVERSED);
  d = new_edge_5->measure();
  new_edge_5->closest_point_trimmed(vi, c_point);

  //PARABOLA_CURVE_TYPE
  RefEdge* new_edge_6 = gmti->make_RefEdge(PARABOLA_CURVE_TYPE, vertex1,
                                        vertex2, &center_pnt);
  d = new_edge_6->measure();
  new_edge_6->closest_point_trimmed(vi, c_point);

  //HYPERBOLA_CURVE_TYPE
  RefEdge* new_edge_7 = gmti->make_RefEdge(HYPERBOLA_CURVE_TYPE, vertex1,
                                        vertex2, &center_pnt);
  d = new_edge_7->measure();
  new_edge_7->closest_point_trimmed(vi, c_point);

  //delete all free vertices and edges
  for (int j = free_entities.size(); j--;)
  {
     gti->delete_RefEntity( free_entities.get_and_step());
  }

  DLIList<RefEdge *> ref_edges;
  gti->ref_edges(ref_edges);

  //make a new refedge out of existing refedge.
  RefEdge* ref_edge = ref_edges.step_and_get();

  RefEdge* new_edge = gmti->make_RefEdge(ref_edge);

  free_entities.clean_out();
  gti->get_free_ref_entities(free_entities);

  //translate the new curve by (10,10,10)
  RefEntity * entity = free_entities.get();
  gti->translate((BasicTopologyEntity*)entity, vector1);

  box = new_edge->get_curve_ptr()->bounding_box();
  // test free curve translation

  box = ref_edge->get_curve_ptr()->bounding_box();

  //general query
  DLIList<OCCCurve*> curves;
  CAST_TO(body, OCCBody)->get_all_curves(curves);

  OCCCurve *curve = curves.step_and_get();

  type = curve->geometry_type();
  // ARC_CURVE_TYPE

  box = curve-> bounding_box();
  // bounding box

  d = ref_edge->measure();
  // length of the curve.  

  ref_edge->get_param_range(lower,upper);
  // paremeter range.

  d = ref_edge->length_from_u(lower,upper);
  // double check whole curve length.

  d = ref_edge->length_from_u(lower/2+upper/2, upper);
  // half curve length.

  periodic = ref_edge->is_periodic(p);
  // if curve is periodic and its period (p). here it's not.

  ref_edge->position_from_u(lower/2+upper/2, vi);
  // middle point.

  u = ref_edge->u_from_position(vi); 
  // middle point's u value.

  double radius;
  ref_edge->closest_point(vi, c_point, &tangent, & curvature1_ptr);
  // Closed point on middle point.

  ref_edge->tangent(vi, tangent);
  // tangent at middle point.

  ref_edge->get_point_direction(c_point, tangent);
  // double check tangent

  ref_edge->get_center_radius(center, radius);
  // center and radius for arc curves

  pc = ref_edge->point_containment(c_point);
  // middle point should be on the curve.

  //delete all entities
  gti->delete_Body(bodies);

  for (int j = free_entities.size(); j--;)
    {
      gti->delete_RefEntity( free_entities.get_and_step());
    }

  return rsl;
}
