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

  //test for creating prisms
  Body* prism1 = gmti->prism(10, 7, 5, 2); 
  Body* prism2 = gmti->prism(10, 7, 5, 5);
  Body* prism3 = gmti->prism(10, 8, 5, 2);
  Body* prism4 = gmti->prism(10, 8, 5, 5);
  CubitStatus rsl = CUBIT_SUCCESS;
  DLIList<RefEntity*> ref_entity_list;
  int num_ents_exported=0;
  const CubitString cubit_version="10.2";
  const char * filename = "prism.occ";
  const char * filetype = "OCC";

  rsl = CubitCompat_export_solid_model(ref_entity_list, filename, filetype,
                                 num_ents_exported, cubit_version);

  CubitBox box1 = prism1->bounding_box();  
  CubitVector min(-5.0, -1.949, -5.0);
  CubitVector max(4.504, 1.949, 5.0);
  CubitBox test_box(min, max);
  assert(box1 >= test_box);
  test_box *= 1.02;
  assert( box1 <= test_box );
  
  box1 = prism2->bounding_box();
  min.set(-5.0, -4.874, -5.0);
  max.set(4.504, 4.874, 5.0);
  test_box.reset(min, max);
  assert(box1 >= test_box);
  test_box *= 1.02;
  assert( box1 <= test_box );

  box1 = prism3->bounding_box();
  min.set(-4.619, -1.847, -5.0);
  max.set(4.619, 1.847, 5.0);
  test_box.reset(min, max);
  assert(box1 >= test_box);
  test_box *= 1.02;
  assert( box1 <= test_box );

  box1 = prism4->bounding_box();
  min.set(-4.619, -4.619, -5.0);
  max.set(4.619, 4.619, 5.0);
  test_box.reset(min, max);
  assert(box1 >= test_box);
  test_box *= 1.02;
  assert( box1 <= test_box );

  DLIList<Body*> bodies;
  gti->bodies(bodies);
  gti->delete_Body(bodies);

  //test for creating pyramids
  Body* pyramid1 = gmti->pyramid(10, 7, 5, 2, 2);
  Body* pyramid2 = gmti->pyramid(10, 7, 5, 5, 2);
  Body* pyramid3 = gmti->pyramid(10, 8, 5, 2, 2);
  Body* pyramid4 = gmti->pyramid(10, 8, 5, 5, 2);
  Body* pyramid5 = gmti->pyramid(10, 4, 5, 2, 0);
  Body* pyramid6 = gmti->pyramid(10, 4, 5, 5, 0);
  ref_entity_list.clean_out();
  num_ents_exported=0;
  filename = "pyramid.occ";

  rsl = CubitCompat_export_solid_model(ref_entity_list, filename, filetype,
                                 num_ents_exported, cubit_version);

  box1 = pyramid1->bounding_box();
  min.set(-5.0, -1.949, -5.0);
  max.set(4.504, 1.949, 5.0);
  test_box.reset(min, max);
  assert(box1 >= test_box); 
  test_box *= 1.02;
  assert( box1 <= test_box );

  box1 = pyramid2->bounding_box();
  min.set(-5.0, -4.874, -5.0);
  max.set(4.504, 4.874, 5.0);
  test_box.reset(min, max);
  assert(box1 >= test_box);
  test_box *= 1.02;
  assert( box1 <= test_box );

  box1 = pyramid3->bounding_box();
  min.set(-4.619, -1.847, -5.0);
  max.set(4.619, 1.847, 5.0);
  test_box.reset(min, max);
  assert(box1 >= test_box);
  test_box *= 1.02;
  assert( box1 <= test_box );

  box1 = pyramid4->bounding_box();
  min.set(-4.619, -4.619, -5.0);
  max.set(4.619, 4.619, 5.0);
  test_box.reset(min, max);
  assert(box1 >= test_box);
  test_box *= 1.02;
  assert( box1 <= test_box ); 

  box1 = pyramid5->bounding_box();
  min.set(-3.535, -1.414, -5.0);
  max.set(3.535, 1.414, 5.0);
  test_box.reset(min, max);
  assert(box1 >= test_box);
  test_box *= 1.02;
  assert( box1 <= test_box );

  box1 = pyramid6->bounding_box(); 
  min.set(-3.535, -3.535, -5.0);
  max.set(3.535, 3.535, 5.0);
  test_box.reset(min, max);
  assert(box1 >= test_box);
  test_box *= 1.02;
  assert( box1 <= test_box );

  bodies.clean_out();
  gti->bodies(bodies);
  gti->delete_Body(bodies);

  //Create sphere
  RefEntity* sphereEnt= GeometryModifyTool::instance()->sphere(1.5);
  sphereEnt->entity_name("sphere");

  box1 = sphereEnt->bounding_box();
  min.set(-1.5,-1.5,-1.5);
  max.set(1.5,1.5,1.5);
  test_box.reset(min, max);
  assert(box1 >= test_box);
  test_box *= 1.02;
  assert( box1 <= test_box );

  TopoDS_Compound* objOCC;
  Body* tmpBd = GeometryQueryTool::instance()->get_first_body();
  DLIList<RefVertex*> vertices;
  tmpBd->ref_vertices(vertices);
  DLIList<RefEdge*> ref_edges;
  tmpBd->ref_edges(ref_edges);
  DLIList<RefFace*> sphere_faces;
  tmpBd->ref_faces(sphere_faces);
  for(int i = 0; i < ref_edges.size(); i++) 
  {
    vertices.clean_out();
    ref_edges.get()->ref_vertices(vertices);
    ref_edges.get_and_step()->measure();
  }
  BodySM* tmpBdSM = tmpBd->get_body_sm_ptr();
  objOCC = ( (OCCBody*) tmpBdSM )->get_TopoDS_Shape(); //Opencascade Object
  //OCCDrawTool::instance()->draw_TopoDS_Shape(objOCC, 200);

  bodies.clean_out();
  gti->bodies(bodies);
  DLIList<RefEntity*>  free_entities;
  gti->get_free_ref_entities(free_entities);
  gti->delete_Body(bodies);
     
  for (int j = free_entities.size(); j--;)
  {
     gti->delete_RefEntity( free_entities.get_and_step());
  }
  // Read in the geometry from files specified on the command line
  const char *argv = "66_shaver3.brep";
  CubitStatus status = read_geometry(1, &argv);
  if (status == CUBIT_FAILURE) exit(1);

  const char *argv2 = "62_shaver1.brep";
  status = read_geometry(1, &argv2);
  if (status == CUBIT_FAILURE) exit(1);

  const char *argv3 = "72_shaver6.brep";
  status = read_geometry(1, &argv3);
  if (status == CUBIT_FAILURE) exit(1);
   
  DLIList<Body*> test_bodies;
  gti->bodies(test_bodies);
  CubitVector vi, vii;
  vi = test_bodies.get()->center_point(); 

  CubitVector axis(10,0,0);

  gti->translate(test_bodies.get(),axis);
  vii = test_bodies.get()->center_point();
  assert(vii - vi == axis);
  // After parellel move, center point moved by x (10)

  CubitVector vector1(10,10,10);
  CubitVector vector2(-10,-10,10);
  CubitVector vector3(10, -10, 10);
  free_entities.clean_out();

  // Make two vertices.
  gmti->make_RefVertex(vector1,5);
  gmti->make_RefVertex(vector2,5);
  gmti->make_RefVertex(vector3,5);
  gti->get_free_ref_entities(free_entities);

  ref_entity_list.clean_out();
  num_ents_exported=0;
  filename = "point.occ";
  filetype = "OCC";
  
  rsl = CubitCompat_export_solid_model(ref_entity_list, filename, filetype, 
                                num_ents_exported, cubit_version);
 
  //check for vertex
  bodies.clean_out();
  gti->bodies(bodies);
  free_entities.clean_out();// get_free_ref_entities directly append
  //without checking for duplicates, so clean_out first. 
  gti->get_free_ref_entities(free_entities);
 
  RefVertex* vertex1 = CAST_TO(free_entities.step_and_get(),RefVertex);
  RefVertex* vertex2 = CAST_TO(free_entities.step_and_get(),RefVertex);
  RefVertex* vertex3 = CAST_TO(free_entities.step_and_get(),RefVertex); 
  CubitBoolean is_equal = gti->
		about_spatially_equal(vertex1,vertex2);
  assert(is_equal == CUBIT_FAILURE);
 
  double d;
  gti->entity_entity_distance(vertex1,vertex2,vi, vii,d);
  assert(d > 28.284 && d < 28.2843);
  // distance (d) between vertex1,vertex2.  
 
  //check for body
  d = bodies.get()->measure(); 
  assert( d > 2237.75 && d < 2237.752);
  // first body's volume.

  vi = bodies.get()->center_point();
  //first body's bounding box's center point.
  min.set(6.67, 0, 14.58);
  assert(vi.distance_between(min) < 0.01 );
 
  CubitBox box = bodies.get()->bounding_box();
  min.set(-11.5894,-11.58944, 7.79849);
  max.set(24.9303,11.58944,21.3615);
  test_box.reset(min,max);
  assert(box >= test_box);
  test_box *= 1.01;
  assert(box <= test_box);
  //first body's bounding box.

  gti->entity_entity_distance(gti->get_first_ref_volume(), vertex2,vi, vii,d);
  //first body and vertex2 's minimum distance(d) and locations for the minimum.
  assert(d <3.64214  && d > 3.64213);

  // test create a Compound body.
  test_bodies.clean_out();
  gti->bodies(test_bodies);

  bodies.clean_out();
  gmti->unite(test_bodies, bodies );
  assert(bodies.size() == 1);

  Body* CompBody = bodies.get();
  BodySM* body = CompBody->get_body_sm_ptr();
  OCCBody* occ_body = CAST_TO(body, OCCBody);
  occ_body->mass_properties(vi, d);

  bodies.last();
  vi = bodies.get()->center_point();
  gti->reflect(bodies, axis);
  vii = bodies.pop()->center_point();
  assert(vi.x() == -vii.x() && vi.y() == vii.y() && vi.z() == vii.z());
  // After reflection, only x value should change.

  vi = CompBody->center_point();
  gti->scale(CompBody,2);
  vii = CompBody->center_point();
  assert(vii == 2*vi);
  // After scale, center point moved by 2 times 

  vi = CompBody->center_point();
  gti->translate(CompBody,axis);
  vii =CompBody->center_point();
  assert(axis.about_equal(vii - vi) );
  // After parellel move, center point moved by x (10)

  vi = CompBody->center_point();
  gti->rotate(CompBody, axis, 3.14/6);
  vii = CompBody->center_point();
  assert(vi.x() == vii.x() && vi.y() != vii.y() && vi.z() != vii.z());
  // After rotation, center point changed in y and z value.

  double d_;
  occ_body->mass_properties(vii, d_);
  assert(fabs(d_ - 8*d) < 0.0001);
  //true center and volume, volume should be 8 times of the original one.

  vi = occ_body->get_bounding_box().center(); 
  // bounding box center.
  assert(vi != vii);

  //check for surface
  DLIList<OCCSurface*> surfaces;
  occ_body->get_all_surfaces(surfaces);
  OCCSurface* surface = surfaces.step_and_get();
  GeometryType type = surface->geometry_type();
  // CONE_SURFACE_TYPE

  box = surface->bounding_box();
  // bounding box

  DLIList<RefFace*> ref_faces;
  gti->ref_faces(ref_faces);
  //  RefFace* ref_face = ref_faces.step_and_get();
  RefFace* ref_face = ref_faces.get();

  //make a new refface out of existing refface.
  CubitBoolean extended_from = CUBIT_TRUE;
  RefFace* new_face = gmti->make_RefFace(ref_face, extended_from);

  ref_entity_list.clean_out();
  ref_entity_list.append(new_face);
  rsl = CubitCompat_export_solid_model(ref_entity_list, filename, filetype,
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
    d = edge->measure();
    start = edge->start_vertex();
    end = edge->end_vertex();
  }

  RefFace* new_face2;
  new_face2 = gmti->make_RefFace(PLANE_SURFACE_TYPE, 
                         *ref_edge_list, CUBIT_TRUE, new_face, CUBIT_TRUE);

  bodies.clean_out();
  gti->bodies(bodies);
  //translate the new face by (20,20,20)
  for(int i = 1; i <= bodies.size(); i++)
  {
     bodies.step();
     if( i != 2)
        continue;
     Body * entity = bodies.get();
     gti->translate(entity, i*vector1);
  }

  RefVolume* volume = NULL;
  if ( new_face->get_surface_ptr()->is_closed_in_U())
    volume = new_face->ref_volume(); 
  else
    vi = new_face->center_point();
  //center point should moved by (20,20,20) compared with the original one below

  vii = ref_face->center_point();
  assert(fabs(vi.distance_between(vii)-vector1.length()*2) < 0.0001);

  CubitVector normal;
  normal = ref_face->normal_at(vii);
  vi.set(0, -0.0091, 0.999958);
  assert(fabs(vi.distance_between( normal )) < 0.0001);
  // surface normal at center point.

  CubitVector closest_location ;
  CubitVector unit_normal_ptr ;
  CubitVector curvature1_ptr ;
  CubitVector curvature2_ptr ; 

  ref_face->find_closest_point_trimmed(vii, closest_location);  
  vi.set(-37.69135, -3.23057, -91.71684);
  assert(fabs(vi.distance_between(closest_location)) < 0.0001); 
  // Found surface projection location for vi.

  double curvature1, curvature2;
  ref_face->get_principal_curvatures(closest_location, curvature1, curvature2);
  assert(curvature1 == 0  && curvature2 == 0);
  // get principal curvatures at center point.

  double area = ref_face->area();
  assert(fabs(area-3) < 0.0001);
  // area of the face

  double u = -2.5;      
  double v = 24.3;
  vi = ref_face->position_from_u_v(u,v);
  vii.set(-37.5, -3.2805, -91.7173);
  assert(fabs(vi.distance_between(vii)) < 0.0001);
  // get location on surface for it's u,v

  ref_face->u_v_from_position(vi, u, v);
  assert(fabs(u + 2.5) < 0.0001 && fabs(v -24.3) < 0.00001);
  // get (u,v) from this vi.

  CubitBoolean periodic = ref_face->is_periodic();
  assert(periodic == CUBIT_FALSE);
  // found if surface is periodic.

  double p = 0; //period
  periodic = ref_face->is_periodic_in_U(p); 
  assert(periodic == CUBIT_FALSE  );
  // found if surface is periodic in U and its period.

  periodic = ref_face->is_periodic_in_V(p);
  assert(periodic == CUBIT_FALSE );

  //All OCC entities are parametric.

  double lower, upper;
  ref_face->get_param_range_U(lower, upper);
  assert(fabs(lower- -5.30865) < 0.0001 && fabs(upper - 0.69135) < 0.0001);

  ref_face->get_param_range_V(lower, upper);
  assert(fabs(lower - 24) < 0.01 && fabs(upper - 24.5) < 0.01 );
  // get surface V direction boundaries. here it's (24,24.5)

  CubitBoolean closed = ref_face->is_closed_in_U();
  assert(closed == CUBIT_FALSE);

  closed = ref_face->is_closed_in_V();
  assert(closed == CUBIT_FALSE);

  CubitPointContainment pc = ref_face->point_containment(7,25);
  assert(pc == CUBIT_PNT_OUTSIDE);

  CubitPointContainment pc2 = ref_face->point_containment(-3,24.2);
  assert(pc2 == CUBIT_PNT_INSIDE);

  ref_edge_loops.clean_out();
  int num_loops;
  num_loops = ref_face->ref_edge_loops(ref_edge_loops);
  DLIList<RefEdge*> *ref_edges1;
  ref_edges1 = ref_edge_loops.get();
  RefEdge* edge1 = ref_edges1->pop();
  RefEdge* edge2 = ref_edges1->pop();
  double angle = edge2->angle_between(edge1, ref_face);
  assert(fabs(angle - 1.5707) < 0.001);

  //test for curve
  CubitVector c_point, tangent, center;

  //make all kinds of curves.
  CubitVector center_pnt(0,0,10);
  DLIList<CubitVector*> list;
  CubitVector center_pnt1(5,8,10);
  CubitVector center_pnt2(1,2,10);
  CubitVector center_pnt3(-2,-3.5,10);
  list.append(&center_pnt1);
  list.append(&center_pnt2);
  list.append(&center_pnt);
  list.append(&center_pnt3);
  RefEdge* new_edge_1 = gmti->make_RefEdge(SPLINE_CURVE_TYPE, vertex1,
                                          vertex2, list);
  d = new_edge_1->measure();
  assert(fabs(d - 28.5)<0.01);

  //Spline with points and tangents
  CubitVector* first_v = new CubitVector(vertex1->coordinates());
  CubitVector* last_v = new CubitVector(vertex2->coordinates());
  list.insert_first(first_v);
  list.append(last_v);

  DLIList<CubitVector*> tangents;
  for (int i = 0; i< 6; i++)
    tangents.append(NULL);

  GeometryModifyEngine* gme = gmti->get_engine(new_edge_1);
  Curve* new_curve = gme->make_Curve(list, tangents);
  RefEdge * new_edge_11 = gti->make_free_RefEdge(new_curve);
  d = new_edge_11->measure();
  assert(fabs(d - 29.5517) <0.001);

  CubitVector v11(-1,-1,0);
  CubitVector v12(3,1,0);
  tangents.clean_out();
  tangents.append(&v11);
  for (int i = 0; i< 4; i++)
    tangents.append(NULL);
  tangents.append(&v12);

  new_curve = gme->make_Curve(list, tangents);
  RefEdge * new_edge_12 = gti->make_free_RefEdge(new_curve);
  d = new_edge_12->measure(); 
  assert(fabs(d - 34.33967) < 0.0001);

  num_ents_exported=0;
  filename = "BsplineCurve.occ";
  ref_entity_list.clean_out();
  ref_entity_list.append(new_edge_11);
  ref_entity_list.append(new_edge_12);
  rsl = CubitCompat_export_solid_model(ref_entity_list, filename, filetype,
                                 num_ents_exported, cubit_version);

  //straight line
  RefEdge* new_edge_2 = gmti->make_RefEdge(STRAIGHT_CURVE_TYPE, vertex1,
                                        vertex2, &center_pnt);
  d = new_edge_2->measure();
  assert(fabs(d - 28.284) <0.001);

  vi.set(-9.374537, 7.2082, 51.8845);
  new_edge_2->closest_point_trimmed(vi, c_point);
  vii.set(-1.083148, -1.083148, 10);
  assert(vii.distance_between( c_point ) < 0.0001);

  //arc curve
  RefEdge* new_edge_3 = gmti->make_RefEdge(ARC_CURVE_TYPE, vertex1,
                                        vertex3, &center_pnt);
  d = new_edge_3->measure();
  assert(fabs(d - 31.4159) < 0.0001);
  new_edge_3->closest_point_trimmed(vi, c_point);
  vii.set(0.62764, 3.486959, 10);
  assert(vii.distance_between( c_point ) < 0.0001);

  //ellipse curve
  RefEdge* new_edge_4 = gmti->make_RefEdge(ELLIPSE_CURVE_TYPE, vertex1,
                                        vertex3, &center_pnt);
  d = new_edge_4->measure();
  assert(fabs(d - 22.21) < 0.01); 
  new_edge_4->closest_point_trimmed(vi, c_point);
  vii.set(10, 10, 10);
  assert(vii.distance_between( c_point ) < 0.0001);

  RefEdge* new_edge_5 = gmti->make_RefEdge(ELLIPSE_CURVE_TYPE, vertex1,
                                        vertex3, &center_pnt, CUBIT_REVERSED);
  d = new_edge_5->measure();
  assert(fabs(d-66.643) < 0.001);
  new_edge_5->closest_point_trimmed(vi, c_point);
  vii.set(-11.2111, 8.62037, 10);
  assert(vii.distance_between( c_point ) < 0.0001);

  //PARABOLA_CURVE_TYPE
  RefEdge* new_edge_6 = gmti->make_RefEdge(PARABOLA_CURVE_TYPE, vertex1,
                                        vertex3, &center_pnt);
  d = new_edge_6->measure();
  assert(fabs(d-29.56546) < 0.0001);
  new_edge_6->closest_point_trimmed(vi, c_point);
  vii.set(0.58077, 2.41, 10);
  assert(vii.distance_between( c_point ) < 0.0001);

  //HYPERBOLA_CURVE_TYPE
  RefEdge* new_edge_7 = gmti->make_RefEdge(HYPERBOLA_CURVE_TYPE, vertex1,
                                        vertex3, &center_pnt);
  d = new_edge_7->measure();
  assert(fabs(d-21.6815) < 0.0001);

  new_edge_7->closest_point_trimmed(vi, c_point);
  vii.set( 6.5835, 2.8855, 10);
  assert(vii.distance_between( c_point ) < 0.0001);

  //delete all free vertices and edges
  for (int j = free_entities.size(); j--;)
  {
     gti->delete_RefEntity( free_entities.get_and_step());
  }

  ref_edges.clean_out();
  gti->ref_edges(ref_edges);

  //make a new refedge out of existing refedge.
  RefEdge* ref_edge = ref_edges.step_and_get();

  RefEdge* new_edge;
  new_edge = gmti->make_RefEdge(ref_edge);

  free_entities.clean_out();
  gti->get_free_ref_entities(free_entities);

  //translate the new curve by (10,10,10)
  RefEntity * entity = free_entities.get();
  box = entity->bounding_box();
  gti->translate((BasicTopologyEntity*)entity, vector1);
  test_box = entity->bounding_box();
  assert(fabs(test_box.minimum().distance_between( box.minimum()) -vector1.length())<0.001 &&
         fabs(test_box.maximum().distance_between( box.maximum()) -vector1.length())< 0.001);
  //general query
  DLIList<OCCCurve*> curves;
  CAST_TO(body, OCCBody)->get_all_curves(curves);

  OCCCurve *curve = curves.step_and_get();

  type = curve->geometry_type();
  assert(type == ARC_CURVE_TYPE);
  // ARC_CURVE_TYPE

  box = curve-> bounding_box();
  min.set(26.75282, 16.7499, 20.3036);
  max.set(33.24717, 23.24399, 20.36297);
  assert(min.distance_between(box.minimum()) < 0.001 &&
         max.distance_between(box.maximum()) < 0.001);
  // bounding box

  d = ref_edge->measure();
  assert(fabs(d - 21.6815) < 0.001);

  ref_edge->get_param_range(lower,upper);
  assert(lower + 1.06127 <= 0.0001 && fabs(upper-1.06127) < 0.001);
  // paremeter range.

  d = ref_edge->length_from_u(lower,upper);
  assert(fabs(d - 21.6815) < 0.001);

  d = ref_edge->length_from_u(lower/2+upper/2, upper);
  assert(fabs(d-10.8407) < 0.0001);
  // half curve length.

  periodic = ref_edge->is_periodic(p);
  assert(periodic == CUBIT_FALSE);
  // if curve is periodic and its period (p). here it's not.

  ref_edge->position_from_u(lower/2+upper/2, vi);
  // middle point.

  u = ref_edge->u_from_position(vi); 
  assert(fabs(u ) < 0.001);
  // middle point's u value.

  double radius;
  ref_edge->closest_point(vi, c_point, &tangent, & curvature1_ptr);
  vii.set(0, -1, 0);
  assert(tangent.distance_between(vii) < 0.001);
  // Closed point on middle point.

  ref_edge->tangent(vi, tangent);
  assert(tangent.distance_between(vii) < 0.001);
  // tangent at middle point.

  ref_edge->get_point_direction(c_point, tangent);
  assert(tangent.distance_between(vii) < 0.001);
  // double check tangent

  ref_edge->get_center_radius(center, radius);
  // center and radius for arc curves

  pc = ref_edge->point_containment(c_point);
  assert(pc == CUBIT_PNT_INSIDE);
  // middle point should be on the curve.

  //delete all entities
  bodies.clean_out();
  gti->bodies(bodies);
  gti->delete_Body(bodies);

  free_entities.clean_out();
  gti->get_free_ref_entities(free_entities);
  for (int j = free_entities.size(); j--;)
    {
      gti->delete_RefEntity( free_entities.get_and_step());
    }

  return rsl;
}
