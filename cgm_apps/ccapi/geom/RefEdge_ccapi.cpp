#include "RefEdge_ccapi.h"

#include "CubitDefines.h"
#include "RefVolume.hpp"
#include "RefFace.hpp"
#include "RefEdge.hpp"
#include "RefVertex.hpp"

#include "DLCoEdgeList.hpp"
#include "DLCubitVectorList.hpp"

#include "CubitVectorStruct.h"

#include "copy_defines.h"

void *RefEdge_get_address(void *this_ref_edge,
                          enum EntityType inputEntityType) 
{
  RefEdge *temp_ref_edge = (RefEdge *) this_ref_edge;
  return temp_ref_edge->get_address(inputEntityType);
}
      
enum CubitStatus RefEdge_get_point_direction(void *this_ref_edge,
                                               /* CubitVector & */ struct CubitVectorStruct *origin, 
                                               /* CubitVector & */ struct CubitVectorStruct *direction ) 
{
  RefEdge *temp_ref_edge = (RefEdge *) this_ref_edge;
  CubitVector temp_origin(*origin);
  CubitVector temp_direction(*direction);
  
  return temp_ref_edge->get_point_direction(temp_origin, temp_direction);
}

enum CubitStatus RefEdge_get_center_radius(void *this_ref_edge,
                                             /* CubitVector & */ struct CubitVectorStruct *center,
                                             /* double& */ double *radius ) 
{
  RefEdge *temp_ref_edge = (RefEdge *) this_ref_edge;
  CubitVector temp_center(*center);

  return temp_ref_edge->get_center_radius(temp_center, *radius);
}
  
enum EntityType RefEdge_entity_type(void *) 
{ return RefEdge_TYPE;}

enum EntityType RefEdge_get_child_ref_entity_type(void *) 
{ return RefVertex_TYPE; }
      
enum EntityType RefEdge_get_parent_ref_entity_type(void *) 
{ return RefFace_TYPE; }
      
enum EntityType RefEdge_get_topology_bridge_type(void *) 
{ return Curve_TYPE ;}
  
enum EntityType RefEdge_get_grouping_entity_type(void *) 
{ return Chain_TYPE; }

enum EntityType RefEdge_get_sense_entity_type(void *) 
{ return CoEdge_TYPE; }

  /* RefVertex * */ void *RefEdge_startRefVertex(void *this_ref_edge) 
{
  RefEdge *temp_ref_edge = (RefEdge *) this_ref_edge;
  return temp_ref_edge->startRefVertex();
}

  /* RefVertex * */ void *RefEdge_endRefVertex(void *this_ref_edge) 
{
  RefEdge *temp_ref_edge = (RefEdge *) this_ref_edge;
  return temp_ref_edge->endRefVertex();
}
    
void RefEdge_switch_vertices(void *this_ref_edge) 
{
  RefEdge *temp_ref_edge = (RefEdge *) this_ref_edge;
  temp_ref_edge->switch_vertices();
}

  /* Chain * */ void *RefEdge_get_first_chain_ptr(void *this_ref_edge) 
{
  RefEdge *temp_ref_edge = (RefEdge *) this_ref_edge;
  return temp_ref_edge->get_first_chain_ptr();
}
  
enum CubitStatus RefEdge_get_co_edges(void *this_ref_edge,
                                        /* DLCoEdgeList & */ void ***co_edges_found_list, int *co_edges_found_list_size,
                                        /* RefFace * */ void *input_ref_face_ptr) 
{
  RefEdge *temp_ref_edge = (RefEdge *) this_ref_edge;

  RefFace *temp_input_ref_face_ptr = (RefFace *) input_ref_face_ptr;
  DLCoEdgeList temp_co_edges_found_list;
  
  CubitStatus result = temp_ref_edge->get_co_edges(temp_co_edges_found_list, temp_input_ref_face_ptr);

  COPY_LIST_TO_ARRAY(temp_co_edges_found_list, *co_edges_found_list, *co_edges_found_list_size);

  return result;
}

double RefEdge_angle_between(void *this_ref_edge,
                               /* RefEdge * */ void *other_edge_ptr, /* RefFace * */ void *face_ptr ) 
{
  RefEdge *temp_ref_edge = (RefEdge *) this_ref_edge;

  RefEdge *temp_other_edge_ptr = (RefEdge *) other_edge_ptr;
  RefFace *temp_face_ptr = (RefFace *) face_ptr;
  
  return temp_ref_edge->angle_between(temp_other_edge_ptr, temp_face_ptr);
}
  
enum CubitSense RefEdge_sense(void *this_ref_edge,
                                /* RefFace * */ void *face ) 
{
  RefEdge *temp_ref_edge = (RefEdge *) this_ref_edge;
  RefFace *temp_face = (RefFace *) face;

  return temp_ref_edge->sense(temp_face);
}
    
int RefEdge_dimension(void *this_ref_edge) 
{
  RefEdge *temp_ref_edge = (RefEdge *) this_ref_edge;
  return temp_ref_edge->dimension();
}

  /* RefVertex * */ void *RefEdge_common_ref_vertex_1(void *this_ref_edge,
                                                        /* RefEdge * */ void *otherRefEdgePtr ) 
{
  RefEdge *temp_ref_edge = (RefEdge *) this_ref_edge;
  RefEdge *temp_otherRefEdgePtr = (RefEdge *) otherRefEdgePtr;
  
  return temp_ref_edge->common_ref_vertex(temp_otherRefEdgePtr);
}
  /* RefVertex * */ void *RefEdge_common_ref_vertex_2(void *this_ref_edge,
                                                        /* RefEdge * */ void *next_ref_edge,
                                                        /* RefFace * */ void *ref_face_ptr ) 
{
  RefEdge *temp_ref_edge = (RefEdge *) this_ref_edge;
  RefEdge *temp_next_ref_edge = (RefEdge *) next_ref_edge;
  RefFace *temp_ref_face_ptr = (RefFace *) ref_face_ptr;

  return temp_ref_edge->common_ref_vertex(temp_next_ref_edge, temp_ref_face_ptr);
}

enum CubitBoolean RefEdge_common_vertices(void *this_ref_edge,
                                            /* RefEdge * */ void *otherRefEdgePtr, 
                                            /* DLRefVertexList & */ void ***common_verts, int *common_verts_size) 
{
  RefEdge *temp_ref_edge = (RefEdge *) this_ref_edge;
  RefEdge *temp_otherRefEdgePtr = (RefEdge *) otherRefEdgePtr;

  DLRefVertexList temp_common_verts;
  
  enum CubitBoolean result = temp_ref_edge->common_vertices(temp_otherRefEdgePtr, temp_common_verts);

  COPY_LIST_TO_ARRAY(temp_common_verts, *common_verts, *common_verts_size);

  return result;
}

  /* RefEdge * */ void *RefEdge_get_other_curve(void *this_ref_edge,
                                                  /* RefVertex * */ void *common_vertex,
                                                  /* RefFace * */ void *ref_face_ptr) 
{
  RefEdge *temp_ref_edge = (RefEdge *) this_ref_edge;

  RefVertex *temp_common_vertex = (RefVertex *) common_vertex;
  
  RefFace *temp_ref_face_ptr = (RefFace *) ref_face_ptr;
  
  return temp_ref_edge->get_other_curve(temp_common_vertex, temp_ref_face_ptr);
}
  
enum CubitStatus RefEdge_get_two_co_edges(void *this_ref_edge,
                                            /* RefEdge * */ void *next_ref_edge,
                                            /* RefFace * */ void *ref_face_ptr,
                                            /* CoEdge *& */ void **co_edge_this,
                                            /* CoEdge *& */ void **co_edge_next ) 
{
  RefEdge *temp_ref_edge = (RefEdge *) this_ref_edge;
  RefEdge *temp_next_ref_edge = (RefEdge *) next_ref_edge;
  
  RefFace *temp_ref_face_ptr = (RefFace *) ref_face_ptr;
  
  CoEdge *temp_co_edge_this;
  CoEdge *temp_co_edge_next;
  
  enum CubitStatus result = temp_ref_edge->get_two_co_edges(temp_next_ref_edge, temp_ref_face_ptr,
                                                            temp_co_edge_this, temp_co_edge_next);

  *co_edge_this = temp_co_edge_this;
  *co_edge_next = temp_co_edge_next;

  return result;
}

  /* RefFace * */ void *RefEdge_other_face(void *this_ref_edge,
                                             /* RefFace * */ void *not_face,
                                             /* RefVolume * */ void *ref_volume) 
{
  RefEdge *temp_ref_edge = (RefEdge *) this_ref_edge;

  RefFace *temp_not_face = (RefFace *) not_face;
  RefVolume *temp_ref_volume = (RefVolume *) ref_volume;
  
  return temp_ref_edge->other_face(temp_not_face, temp_ref_volume);
}

  /* RefVertex * */ void *RefEdge_other_vertex(void *this_ref_edge,
                                                 /* RefVertex * */ void *refVertexPtr) 
{
  RefEdge *temp_ref_edge = (RefEdge *) this_ref_edge;
  RefVertex *temp_refVertexPtr = (RefVertex *) refVertexPtr;
  
  return temp_ref_edge->other_vertex(temp_refVertexPtr);
}

  /* RefVertex * */ void *RefEdge_closest_vertex(void *this_ref_edge,
                                                 const /* CubitVector & */ struct CubitVectorStruct *point3) 
{
  RefEdge *temp_ref_edge = (RefEdge *) this_ref_edge;
  return temp_ref_edge->closest_vertex(*point3);
}
  
  /* Curve * */ void *RefEdge_get_curve_ptr(void *this_ref_edge)  
{
  RefEdge *temp_ref_edge = (RefEdge *) this_ref_edge;
  return temp_ref_edge->get_curve_ptr();
}
      
struct CubitVectorStruct RefEdge_start_coordinates(void *this_ref_edge) 
{
  RefEdge *temp_ref_edge = (RefEdge *) this_ref_edge;
  return temp_ref_edge->start_coordinates();
}
struct CubitVectorStruct RefEdge_end_coordinates(void *this_ref_edge) 
{
  RefEdge *temp_ref_edge = (RefEdge *) this_ref_edge;
  return temp_ref_edge->end_coordinates();
}
      
void RefEdge_move_to_curve (void *this_ref_edge,
                              /* CubitVector & */ struct CubitVectorStruct *vector ) 
{
  RefEdge *temp_ref_edge = (RefEdge *) this_ref_edge;
  CubitVector temp_vector(*vector);
  
  temp_ref_edge->move_to_curve (temp_vector);

  *vector = temp_vector;
}

enum CubitStatus RefEdge_get_interior_extrema(void *this_ref_edge,
                                                /* DLCubitVectorList& */ struct CubitVectorStruct **interior_points,
                                              int *interior_points_size,
                                                /* CubitSense& */ enum CubitSense *return_sense) 
{
  RefEdge *temp_ref_edge = (RefEdge *) this_ref_edge;
  DLCubitVectorList temp_interior_points;

  
  enum CubitStatus result = temp_ref_edge->get_interior_extrema(temp_interior_points, *return_sense);
  COPY_LIST_TO_STRUCTARRAY(temp_interior_points, *interior_points, *interior_points_size,
                           CubitVectorStruct);

  return result;
}

enum CubitStatus RefEdge_closest_point_trimmed(void *this_ref_edge,
                                                 /* CubitVector const& */ struct CubitVectorStruct *location, 
                                                 /* CubitVector & */ struct CubitVectorStruct *closest_location) 
{
  RefEdge *temp_ref_edge = (RefEdge *) this_ref_edge;
  CubitVector temp_closest_location;
  
  enum CubitStatus result = temp_ref_edge->closest_point_trimmed(*location, temp_closest_location);

  *closest_location = temp_closest_location;

  return result;
}
  
enum CubitStatus RefEdge_closest_point(void *this_ref_edge,
                                         /* CubitVector const& */ struct CubitVectorStruct *location, 
                                         /* CubitVector & */ struct CubitVectorStruct *closest_location,
                                         /* CubitVector * */ struct CubitVectorStruct *tangent_ptr,
                                         /* CubitVector * */ struct CubitVectorStruct *curvature_ptr) 
{
  RefEdge *temp_ref_edge = (RefEdge *) this_ref_edge;

  CubitVector temp_closest_location;
  CubitVector *temp_tangent_ptr = (tangent_ptr ? new CubitVector(*tangent_ptr) :
                                   NULL);
  CubitVector *temp_curvature_ptr = (curvature_ptr ? new CubitVector(*curvature_ptr) :
                                   NULL);
  
  
  enum CubitStatus result = temp_ref_edge->closest_point(*location, temp_closest_location,
                                                         temp_tangent_ptr, temp_curvature_ptr);

  *closest_location = temp_closest_location;
  
  if (temp_tangent_ptr) {
    *tangent_ptr = *temp_tangent_ptr;
    delete temp_tangent_ptr;
  }
  if (temp_curvature_ptr) {
    *curvature_ptr = *temp_curvature_ptr;
    delete temp_curvature_ptr;
  }

  return result;
}
      
void RefEdge_tangent_1(void *this_ref_edge,
                       const /* CubitVector & */ struct CubitVectorStruct *point,
                         /* CubitVector & */ struct CubitVectorStruct *tangent_vec ) 
{
  RefEdge *temp_ref_edge = (RefEdge *) this_ref_edge;

  CubitVector temp_tangent_vec, temp_point(*point);
  
  temp_ref_edge->tangent(temp_point, temp_tangent_vec);

  *tangent_vec = temp_tangent_vec;
}

enum CubitStatus RefEdge_tangent_2(void *this_ref_edge,
                                   const /* CubitVector & */ struct CubitVectorStruct *point,
                                     /* CubitVector & */ struct CubitVectorStruct *tangent_vec,
                                     /* RefFace * */ void *ref_face_ptr) 
{
  RefEdge *temp_ref_edge = (RefEdge *) this_ref_edge;
  RefFace *temp_ref_face_ptr = (RefFace *) ref_face_ptr;
  
  CubitVector temp_tangent_vec;

  enum CubitStatus result = temp_ref_edge->tangent(*point, temp_tangent_vec, temp_ref_face_ptr);

  *tangent_vec = temp_tangent_vec;
  
  return result;
  
}
  
enum CubitStatus RefEdge_tangent_3(void *this_ref_edge,
                                   const /* CubitVector & */ struct CubitVectorStruct *point,
                                     /* CubitVector & */ struct CubitVectorStruct *tangent_vec,
                                     /* RefEdge * */ void *next_ref_edge,
                                     /* RefFace * */ void *ref_face_ptr ) 
{
  RefEdge *temp_ref_edge = (RefEdge *) this_ref_edge;
  RefEdge *temp_next_ref_edge = (RefEdge *) next_ref_edge;
  RefFace *temp_ref_face_ptr = (RefFace *) ref_face_ptr;

  CubitVector temp_tangent_vec;
  
  enum CubitStatus result = temp_ref_edge->tangent(*point, temp_tangent_vec,
                                                   temp_next_ref_edge, temp_ref_face_ptr);
  *tangent_vec = temp_tangent_vec;

  return result;
  
}
  
double RefEdge_measure(void *this_ref_edge) 
{
  RefEdge *temp_ref_edge = (RefEdge *) this_ref_edge;
  return temp_ref_edge->measure();
}

/* CubitString */ const char *measure_label(void *this_ref_edge) 
{
  RefEdge *temp_ref_edge = (RefEdge *) this_ref_edge;
  return temp_ref_edge->measure_label().c_str();
}

double RefEdge_get_arc_length_1(void *this_ref_edge) 
{
  RefEdge *temp_ref_edge = (RefEdge *) this_ref_edge;
  return temp_ref_edge->get_arc_length();
}

double RefEdge_get_arc_length_2(void *this_ref_edge,
                                const /* CubitVector & */ struct CubitVectorStruct *point1,
                                const /* CubitVector & */ struct CubitVectorStruct *point2 ) 
{
  RefEdge *temp_ref_edge = (RefEdge *) this_ref_edge;
  return temp_ref_edge->get_arc_length(*point1, *point2);
}

double RefEdge_get_arc_length_3(void *this_ref_edge,
                                const /* CubitVector & */ struct CubitVectorStruct *point1, const int whichEnd ) 
{
  RefEdge *temp_ref_edge = (RefEdge *) this_ref_edge;
  return temp_ref_edge->get_arc_length(*point1, whichEnd);
}

double RefEdge_get_chord_length(void *this_ref_edge) 
{
  RefEdge *temp_ref_edge = (RefEdge *) this_ref_edge;
  return temp_ref_edge->get_chord_length();
}
  
/* CubitVector */ struct CubitVectorStruct RefEdge_center_point(void *this_ref_edge) 
{
  RefEdge *temp_ref_edge = (RefEdge *) this_ref_edge;
  return temp_ref_edge->center_point();
}

enum CubitStatus RefEdge_mid_point_1(void *this_ref_edge,
                                     const /* CubitVector & */ struct CubitVectorStruct *point1,
                                     const /* CubitVector & */ struct CubitVectorStruct *point2,
                                       /* CubitVector & */ struct CubitVectorStruct *midPoint ) 
{
  RefEdge *temp_ref_edge = (RefEdge *) this_ref_edge;
  CubitVector temp_midPoint;
  
  enum CubitStatus result = temp_ref_edge->mid_point(*point1, *point2, temp_midPoint);

  *midPoint = temp_midPoint;

  return result;
}

enum CubitStatus RefEdge_mid_point_2(void *this_ref_edge,
                                       /* CubitVector & */ struct CubitVectorStruct *mid_point) 
{
  RefEdge *temp_ref_edge = (RefEdge *) this_ref_edge;
  CubitVector temp_mid_point;

  enum CubitStatus result = temp_ref_edge->mid_point(temp_mid_point);

  *mid_point = temp_mid_point;

  return result;
}

enum CubitStatus RefEdge_position_from_fraction(void *this_ref_edge,
                                                const double fraction_along_curve,
                                                  /* CubitVector & */ struct CubitVectorStruct *output_position ) 
{
  RefEdge *temp_ref_edge = (RefEdge *) this_ref_edge;
  CubitVector temp_output_position;
  
  enum CubitStatus result = temp_ref_edge->position_from_fraction(fraction_along_curve, temp_output_position);

  *output_position = temp_output_position;
  
  return result;
}

double RefEdge_start_param(void *this_ref_edge) 
{
  RefEdge *temp_ref_edge = (RefEdge *) this_ref_edge;
  return temp_ref_edge->start_param();
}
  
double RefEdge_end_param(void *this_ref_edge) 
{
  RefEdge *temp_ref_edge = (RefEdge *) this_ref_edge;
  return temp_ref_edge->end_param();
}
  
enum CubitBoolean RefEdge_get_param_range(void *this_ref_edge,
                                            /* double& */ double *start_param,
                                            /* double& */ double *end_param ) 
{
  RefEdge *temp_ref_edge = (RefEdge *) this_ref_edge;
  return temp_ref_edge->get_param_range(*start_param, *end_param);
}

double RefEdge_u_from_position (void *this_ref_edge,
                                const /* CubitVector & */ struct CubitVectorStruct *input_position) 
{
  RefEdge *temp_ref_edge = (RefEdge *) this_ref_edge;
  return temp_ref_edge->u_from_position (*input_position);
}

enum CubitStatus RefEdge_position_from_u (void *this_ref_edge,
                                          double u_value,
                                            /* CubitVector & */ struct CubitVectorStruct *output_position) 
{
  RefEdge *temp_ref_edge = (RefEdge *) this_ref_edge;
  CubitVector temp_output_position;
  
  enum CubitStatus result = temp_ref_edge->position_from_u (u_value, temp_output_position);

  *output_position = temp_output_position;

  return result;
}

double RefEdge_u_from_arc_length (void *this_ref_edge,
                                  double root_param,
                                  double arc_length ) 
{
  RefEdge *temp_ref_edge = (RefEdge *) this_ref_edge;
  return temp_ref_edge->u_from_arc_length (root_param, arc_length);
}
      
enum CubitStatus RefEdge_point_from_arc_length (void *this_ref_edge,
                                                const /* CubitVector & */ struct CubitVectorStruct *root_point, 
                                                double const arc_length,
                                                  /* CubitVector & */ struct CubitVectorStruct *new_point ) 
{
  RefEdge *temp_ref_edge = (RefEdge *) this_ref_edge;
  CubitVector temp_new_point;
  
  enum CubitStatus result = temp_ref_edge->point_from_arc_length (*root_point, arc_length, temp_new_point);

  *new_point = temp_new_point;
  
  return result;
}
  
double RefEdge_length_from_u(void *this_ref_edge,
                             double parameter1,
                             double parameter2 ) 
{
  RefEdge *temp_ref_edge = (RefEdge *) this_ref_edge;
  return temp_ref_edge->length_from_u(parameter1, parameter2);
}

enum CubitBoolean RefEdge_is_periodic_1(void *this_ref_edge) 
{
  RefEdge *temp_ref_edge = (RefEdge *) this_ref_edge;
  return temp_ref_edge->is_periodic();
}
  
enum CubitBoolean RefEdge_is_periodic_2(void *this_ref_edge,
                                          /* double& */ double *period) 
{
  RefEdge *temp_ref_edge = (RefEdge *) this_ref_edge;
  return temp_ref_edge->is_periodic(*period);
}

int RefEdge_get_mark(void *this_ref_edge) 
{
  RefEdge *temp_ref_edge = (RefEdge *) this_ref_edge;
  return temp_ref_edge->get_mark();
}
void RefEdge_set_mark(void *this_ref_edge,
                      int set_value ) 
{
  RefEdge *temp_ref_edge = (RefEdge *) this_ref_edge;
  temp_ref_edge->set_mark(set_value);
}

enum CubitStatus RefEdge_relative_sense(void *this_ref_edge,
                                          /* RefEdge * */ void *ref_edge_ptr_2,
                                        double tolerance_factor,
                                        enum CubitSense *sense,
                                          /* CubitBoolean & */ enum CubitBoolean *spatially_equal,
                                          /* CubitBoolean & */ enum CubitBoolean *tangent_warning,
                                        enum CubitBoolean force_merge) 
{
  RefEdge *temp_ref_edge = (RefEdge *) this_ref_edge;
  RefEdge *temp_ref_edge_ptr_2 = (RefEdge *) ref_edge_ptr_2;
  
  return temp_ref_edge->relative_sense(temp_ref_edge_ptr_2, tolerance_factor, sense,
                                       *spatially_equal, *tangent_warning, force_merge);
}

enum CubitBoolean RefEdge_about_spatially_equal(void *this_ref_edge,
                                                  /* RefEdge * */ void *ref_edge_ptr_2,
                                                double tolerance_factor,
                                                enum CubitSense *sensePtr,
                                                enum CubitBoolean notify_refEntity,
                                                enum CubitBoolean force_merge) 
{
  RefEdge *temp_ref_edge = (RefEdge *) this_ref_edge;
  RefEdge *temp_ref_edge_ptr_2 = (RefEdge *) ref_edge_ptr_2;
  
  return temp_ref_edge->about_spatially_equal(temp_ref_edge_ptr_2, tolerance_factor,
                                              sensePtr, notify_refEntity, force_merge);
}


int RefEdge_validate(void *this_ref_edge) 
{
  RefEdge *temp_ref_edge = (RefEdge *) this_ref_edge;
  return temp_ref_edge->validate();
}

enum CubitBoolean RefEdge_is_tolerant(void *this_ref_edge) 
{
  RefEdge *temp_ref_edge = (RefEdge *) this_ref_edge;
  return temp_ref_edge->is_tolerant();
}

  /* RefVertex * */ void *RefEdge_get_startRefVertex(void *this_ref_edge) 
{
  RefEdge *temp_ref_edge = (RefEdge *) this_ref_edge;
  return temp_ref_edge->get_startRefVertex();
}
  /* RefVertex * */ void *RefEdge_get_endRefVertex(void *this_ref_edge) 
{
  RefEdge *temp_ref_edge = (RefEdge *) this_ref_edge;
  return temp_ref_edge->get_endRefVertex();
}
