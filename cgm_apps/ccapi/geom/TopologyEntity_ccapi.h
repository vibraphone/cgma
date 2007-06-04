#ifndef TOPOLOGY_ENTITY_CCAPI
#define TOPOLOGY_ENTITY_CCAPI

#include "CubitDefines.h"
#include "EntityType.h"

#ifdef __cplusplus
extern "C" {
#endif

void *TopologyEntity_get_address(enum EntityType inputEntityType);

void TopologyEntity_update(enum EventType /*eventType*/);

/* GeometricModelingEngine* */ void * TopologyEntity_get_geometric_modeling_engine();

enum CubitBoolean TopologyEntity_is_directly_related( /* TopologyEntity* */ void * entity );

    /* static */ enum CubitStatus TopologyEntity_get_related_entities_1(/* DLModelEntityList & */ void ***in_list, int *in_list_size,
                                        enum EntityType target_type,
                                                         /* DLModelEntityList & */ void ***out_list, int *out_list_size,
                                        enum CubitBoolean unique);

    /* static */ enum CubitStatus TopologyEntity_get_related_entities_2(/* DLBodyList & */ void ***in_list, int *in_list_size,
                                      enum EntityType target_type,
                                                         /* DLRefEntityList & */ void ***out_list, int *out_list_size,
                                      enum CubitBoolean unique);

  enum CubitStatus TopologyEntity_bodies(/* DLBodyList& */ void ***body_list, int *body_list_size);
  enum CubitStatus TopologyEntity_shells(/* DLShellList& */ void ***shell_list, int *shell_list_size);
  enum CubitStatus TopologyEntity_loops( /* DLLoopList& */ void ***loop_list, int *loop_list_size);
  enum CubitStatus TopologyEntity_chains(/* DLChainList& */ void ***chain_list, int *chain_list_size);
  enum CubitStatus TopologyEntity_ref_volumes(/* DLRefVolumeList& */ void ***ref_volume_list, int *ref_volume_list_size);
  enum CubitStatus TopologyEntity_ref_faces(/* DLRefFaceList& */ void ***ref_face_list, int *ref_face_list_size);
  enum CubitStatus TopologyEntity_ref_edges(/* DLRefEdgeList& */ void ***ref_edge_list, int *ref_edge_list_size);
  enum CubitStatus TopologyEntity_ref_vertices(/* DLRefVertexList& */ void ***ref_vertex_list, int *ref_vertex_list_size);
  enum CubitStatus TopologyEntity_co_volumes(/* DLCoVolumeList& */ void ***co_volume_list, int *co_volume_list_size);
  enum CubitStatus TopologyEntity_co_faces(/* DLCoFaceList& */ void ***co_face_list, int *co_face_list_size);
  enum CubitStatus TopologyEntity_co_edges(/* DLCoEdgeList& */ void ***co_edge_list, int *co_edge_list_size);
  enum CubitStatus TopologyEntity_co_vertices(/* DLCoVertexList& */ void ***co_vertex_list, int *co_vertex_list_size);

/* RefVertex **/ void *TopologyEntity_ref_vertex();
/* RefEdge **/ void *TopologyEntity_ref_edge();
/* RefFace **/ void *TopologyEntity_ref_face();
/* RefVolume **/ void *TopologyEntity_ref_volume();
/* Body **/ void *TopologyEntity_body();

/* BridgeManager* */ void *TopologyEntity_bridge_manager();

  enum EntityType TopologyEntity_get_topology_bridge_type();
  


#ifdef __cplusplus
}
#endif

#endif
