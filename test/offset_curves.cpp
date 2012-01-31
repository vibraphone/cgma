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
  if ( ret_val > 5 )
  {
    PRINT_ERROR("Errors found during Mergechk session.\n");
  }
  return ret_val-5;
  
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

  CubitVector vector1(10,10,10);
  CubitVector vector2(-10,-10,10);
  CubitVector vector3(10, -10, 10);

  gmti->make_RefVertex(vector1,5);
  gmti->make_RefVertex(vector2,5);
  gmti->make_RefVertex(vector3,5);
  DLIList<RefEntity*>  free_entities;
  gti->get_free_ref_entities(free_entities);

  DLIList<RefEntity*> ref_entity_list;
  const CubitString cubit_version="10.2";
  int num_ents_exported=0;
  const char *filename = "curves.brep";
  const char *filetype = "OCC";
 
  RefVertex* vertex1 = CAST_TO(free_entities.step_and_get(),RefVertex);
  RefVertex* vertex2 = CAST_TO(free_entities.step_and_get(),RefVertex);
  RefVertex* vertex3 = CAST_TO(free_entities.step_and_get(),RefVertex); 
 
  //make all kinds of curves.
  CubitVector center_pnt(0,0,10);
  DLIList<CubitVector*> list;
  CubitVector center_pnt1(5,8,10);
  CubitVector center_pnt2(1,2,10);
  CubitVector center_pnt3(-2,-3.5,10);
  list.append(&vector1);
  list.append(&center_pnt1);
  list.append(&center_pnt2);
  list.append(&center_pnt);
  list.append(&center_pnt3);
  list.append(&vector2);
  RefEdge* new_edge_1 = gmti->make_RefEdge(SPLINE_CURVE_TYPE, vertex1,
                                          vertex2, list);
  double d = new_edge_1->measure();
  assert(fabs(d - 29.55)<0.01);

  //Added test for offset curve.
  DLIList<RefEdge*> offset_edges;
  offset_edges.append(new_edge_1);
  CubitVector v12(3,1,0);
  CubitStatus rsl = gmti->offset_curves(offset_edges, 1,  v12);
  assert(rsl == CUBIT_FAILURE);
 
  //straight line
  RefEdge* new_edge_2 = gmti->make_RefEdge(STRAIGHT_CURVE_TYPE, vertex1,
                                        vertex2, &center_pnt);
  d = new_edge_2->measure();
  assert(fabs(d - 28.284) <0.001);

  //Added test for offset curve.
  CubitVector v_null(0,0,0);
  offset_edges.clean_out();
  offset_edges.append(new_edge_2);
  CubitVector offset_v(1,-1,0);
  rsl = gmti->offset_curves(offset_edges, 1,  offset_v);
  DLIList<RefEdge*> ref_edges;
  gti->ref_edges(ref_edges);
  ref_edges.last();
  RefEdge* offset_edge = ref_edges.get();
  d = offset_edge->measure(); 
  assert(fabs(d - 28.284) <0.001);

  RefVertex* start_v = offset_edge->start_vertex();
  CubitVector loc = start_v->coordinates();
  CubitVector closest_p(10.70710, 9.292893, 10);
  assert(fabs(loc.distance_between(closest_p)) < 0.0001);

  rsl = gmti->offset_curves(offset_edges, 1,  offset_v);
  ref_edges.clean_out();
  gti->ref_edges(ref_edges);
  ref_edges.last();
  offset_edge = ref_edges.get();
  d = offset_edge->measure();
  assert(fabs(d - 28.284) <0.001);

  start_v = offset_edge->start_vertex();
  loc = start_v->coordinates();
  assert(fabs(loc.distance_between(closest_p)) < 0.0001);

  rsl = gmti->offset_curves(offset_edges, 1, v_null);
  assert(rsl == CUBIT_FAILURE);

  //arc curve
  RefEdge* new_edge_3 = gmti->make_RefEdge(ARC_CURVE_TYPE, vertex1,
                                        vertex3, &center_pnt);
  d = new_edge_3->measure();
  assert(fabs(d - 31.4159) < 0.0001);

  CubitVector p_center(20, 0, 10);
  RefEdge* new_edge_31 = gmti->make_RefEdge(ARC_CURVE_TYPE, vertex3,
                                        vertex1, &p_center);
  d = new_edge_31->measure();
  assert(fabs(d - 31.4159) < 0.0001);

  offset_edges.clean_out();
  offset_edges.append(new_edge_31);
  rsl = gmti->offset_curves(offset_edges, -1,  v12);
  ref_edges.clean_out();
  gti->ref_edges(ref_edges); 
  ref_edges.last();
  offset_edge = ref_edges.get();
  d = offset_edge->measure();
  assert(fabs(d-28.2743) < 0.0001);
  CubitVector curve_center;
  curve_center = offset_edge->center_point();
  CubitVector center_p(19, 0, 10);
  assert(fabs(curve_center.distance_between(center_p)) < 0.0001);

  offset_edges.clean_out();
  offset_edges.append(new_edge_3);
  rsl = gmti->offset_curves(offset_edges, 1,  v12);
  ref_edges.clean_out();
  gti->ref_edges(ref_edges);
  ref_edges.last();
  offset_edge = ref_edges.get();
  d = offset_edge->measure();
  assert(fabs(d-34.5575) < 0.0001);
  curve_center = offset_edge->center_point();
  CubitVector center_p2(-1, 0, 10);
  assert(fabs(curve_center.distance_between(center_p2)) < 0.0001);

  offset_edges.clean_out();
  offset_edges.append(new_edge_31);
  offset_edges.append(new_edge_3);
  rsl = gmti->offset_curves(offset_edges, 1,  v12);
  ref_edges.clean_out();
  gti->ref_edges(ref_edges);
  ref_edges.last();
  offset_edge = ref_edges.get();
  d = offset_edge->measure();
  assert(fabs(d-34.5575) < 0.0001);

  ref_edges.back();
  offset_edge = ref_edges.get();
  d = offset_edge->measure();
  assert(fabs(d-34.5575) < 0.0001);
  
  //ellipse curve
  RefEdge* new_edge_4 = gmti->make_RefEdge(ELLIPSE_CURVE_TYPE, vertex1,
                                        vertex3, &center_pnt);
  d = new_edge_4->measure();
  assert(fabs(d - 22.21) < 0.01); 

  offset_edges.clean_out();
  offset_edges.append(new_edge_4);
  rsl = gmti->offset_curves(offset_edges, 1.0, v_null);
  assert(rsl == CUBIT_FAILURE);

  //PARABOLA_CURVE_TYPE
  RefEdge* new_edge_6 = gmti->make_RefEdge(PARABOLA_CURVE_TYPE, vertex1,
                                        vertex3, &center_pnt);
  d = new_edge_6->measure();
  assert(fabs(d-29.56546) < 0.0001);

  offset_edges.clean_out();
  offset_edges.append(new_edge_6);
  rsl = gmti->offset_curves(offset_edges, 2.0, v_null );
  assert(rsl == CUBIT_FAILURE);

  //HYPERBOLA_CURVE_TYPE
  RefEdge* new_edge_7 = gmti->make_RefEdge(HYPERBOLA_CURVE_TYPE, vertex1,
                                        vertex3, &center_pnt);
  d = new_edge_7->measure();
  assert(fabs(d-21.6815) < 0.0001);

  offset_edges.clean_out();
  offset_edges.append(new_edge_7);
  rsl = gmti->offset_curves(offset_edges, -2.0, v_null);
  assert(rsl == CUBIT_FAILURE);

  rsl = CubitCompat_export_solid_model(ref_entity_list, filename, filetype,
                                 num_ents_exported, cubit_version);

  free_entities.clean_out();
  gti->get_free_ref_entities(free_entities);

  //delete all free vertices and edges
  for (int j = free_entities.size(); j--;)
  {
     gti->delete_RefEntity( free_entities.get_and_step());
  }

  return rsl;
}
