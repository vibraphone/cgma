#include "RefVertex_ccapi.h"

#include "Point.hpp"
#include "RefEdge.hpp"
#include "RefFace.hpp"
#include "RefVolume.hpp"
#include "RefVertex.hpp"

#include "copy_defines.h"

  /* Point * */ void *RefVertex_get_point_ptr(void *this_ref_vertex)
{
  RefVertex *temp_vertex = (RefVertex *) this_ref_vertex;
  return temp_vertex->get_point_ptr();
}

enum EntityType RefVertex_entity_type(void *)  
{ return RefVertex_TYPE; }

void *RefVertex_get_address(void *this_ref_vertex, enum EntityType inputEntityType)
{
  RefVertex *temp_vertex = (RefVertex *) this_ref_vertex;
  return temp_vertex->get_address(inputEntityType);
}

enum EntityType RefVertex_get_child_ref_entity_type(void *) 
{ return InvalidEntity_TYPE; }

enum EntityType RefVertex_get_topology_bridge_type(void *)
{ return Point_TYPE;}

struct CubitVectorStruct RefVertex_coordinates (void *this_ref_vertex)
{
  RefVertex *temp_vertex = (RefVertex *) this_ref_vertex;
  return temp_vertex->coordinates();
}

struct CubitVectorStruct RefVertex_center_point(void *this_ref_vertex)
{
  RefVertex *temp_vertex = (RefVertex *) this_ref_vertex;
  return temp_vertex->center_point();
}

int RefVertex_dimension(void *)
{ return 0; }

void RefVertex_ref_edges_of_refvolume (void *this_ref_vertex, /* DLRefEdgeList& */ void ***edge_list, int *edge_list_size, 
                                         /* RefVolume * */ void *ref_volume_ptr )
{
  RefVertex *temp_vertex = (RefVertex *) this_ref_vertex;

  DLRefEdgeList temp_edge_list;
  RefVolume *temp_ref_volume_ptr = (RefVolume *) ref_volume_ptr;
  
  temp_vertex->ref_edges_of_refvolume (temp_edge_list, temp_ref_volume_ptr);

  COPY_LIST_TO_ARRAY(temp_edge_list, *edge_list, *edge_list_size);
}

  /* RefEdge * */ void *RefVertex_common_ref_edge(void *this_ref_vertex, /* RefVertex * */ void *other_vertex, 
                                                    /* RefFace * */ void *owning_face)
{
  RefVertex *temp_vertex = (RefVertex *) this_ref_vertex;
  RefVertex *temp_other_vertex = (RefVertex *) other_vertex;
  RefFace *temp_owning_face = (RefFace *) owning_face;
  
  return temp_vertex->common_ref_edge(temp_other_vertex, temp_owning_face);
}

enum CubitBoolean RefVertex_about_spatially_equal(void *this_ref_vertex, /* RefVertex * */ void *ref_vertex_ptr_2,
                                                  double tolerance_factor,
                                                  enum CubitBoolean notify_ref_entity,
                                                  enum CubitBoolean force_merge)
{
  RefVertex *temp_vertex = (RefVertex *) this_ref_vertex;
  RefVertex *temp_ref_vertex_ptr_2 = (RefVertex *) ref_vertex_ptr_2;

  return temp_vertex->about_spatially_equal(temp_ref_vertex_ptr_2, tolerance_factor, notify_ref_entity, force_merge);
}

enum CubitBoolean RefVertex_move(void *this_ref_vertex,  struct CubitVectorStruct delta )
{
  RefVertex *temp_vertex = (RefVertex *) this_ref_vertex;
  CubitVector temp_delta(delta);
  
  return temp_vertex->move(temp_delta);
}

