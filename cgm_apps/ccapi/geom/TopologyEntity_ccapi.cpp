#include "TopologyEntity_ccapi.h"
#include "CubitDefines.h"
#include "copy_defines.h"
#include "EntityType.h"

#include "TopologyEntity.hpp"
#include "TopologyBridge.hpp"

#include "Body.hpp"
#include "DLModelEntityList.hpp"
#include "DLRefEntityList.hpp"
#include "DLBodyList.hpp"
#include "DLRefVolumeList.hpp"
#include "DLRefFaceList.hpp"
#include "DLRefEdgeList.hpp"
#include "DLRefVertexList.hpp"
#include "DLShellList.hpp"
#include "DLLoopList.hpp"
#include "DLChainList.hpp"
#include "DLCoVolumeList.hpp"
#include "DLCoFaceList.hpp"
#include "DLCoEdgeList.hpp"
#include "DLCoVertexList.hpp"

void *TopologyEntity_get_address(void *this_topology_entity, enum EntityType inputEntityType)
{

  return ((TopologyEntity *)this_topology_entity)->get_address(inputEntityType);
}

void TopologyEntity_update(void *this_topology_entity, enum EventType eventType)
{
  ((TopologyEntity *)this_topology_entity)->update(eventType);
}

/* GeometricModelingEngine* */ void * TopologyEntity_get_geometric_modeling_engine(void *this_topology_entity)
{
  return ((TopologyEntity *)this_topology_entity)->get_geometric_modeling_engine();
}

enum CubitBoolean TopologyEntity_is_directly_related(void *this_topology_entity,  /* TopologyEntity* */ void * entity )
{
  TopologyEntity *temp_entity = (TopologyEntity *) entity;
  
  return ((TopologyEntity *)this_topology_entity)->is_directly_related(temp_entity);
}

  /* static */ enum CubitStatus TopologyEntity_get_related_entities_1(/* DLModelEntityList & */ void ***in_list, int *in_list_size,
                                                                      enum EntityType target_type,
                                                                        /* DLModelEntityList & */ void ***out_list, int *out_list_size,
                                                                      enum CubitBoolean unique)
{
  DLModelEntityList temp_in_list;
  DLModelEntityList temp_out_list;

  COPY_ARRAY_TO_LIST(*in_list, *in_list_size, temp_in_list);
  
  enum CubitStatus result = TopologyEntity::get_related_entities(temp_in_list, target_type, temp_out_list, unique);

  COPY_LIST_TO_ARRAY(temp_out_list, *out_list, *out_list_size);
  
  return result;
}

    /* static */ enum CubitStatus TopologyEntity_get_related_entities_2(/* DLBodyList & */ void ***in_list, int *in_list_size,
                                      enum EntityType target_type,
                                                         /* DLRefEntityList & */ void ***out_list, int *out_list_size,
                                      enum CubitBoolean unique)
{
  DLBodyList temp_in_list;
  DLRefEntityList temp_out_list;

  COPY_ARRAY_TO_LIST(*in_list, *in_list_size, temp_in_list);
  
  enum CubitStatus result = TopologyEntity::get_related_entities(temp_in_list, target_type, temp_out_list, unique);

  COPY_LIST_TO_ARRAY(temp_out_list, *out_list, *out_list_size);
  
  return result;
}

  enum CubitStatus TopologyEntity_bodies(void *this_topology_entity, /* DLBodyList& */ void ***body_list, int *body_list_size)
{
  DLBodyList temp_body_list;
  
  enum CubitStatus result = ((TopologyEntity *)this_topology_entity)->bodies(temp_body_list);

    COPY_LIST_TO_ARRAY(temp_body_list, *body_list, *body_list_size);
    
    return result;
}
  enum CubitStatus TopologyEntity_shells(void *this_topology_entity, /* DLShellList& */ void ***shell_list, int *shell_list_size)
{
  DLShellList temp_shell_list;
  
  enum CubitStatus result = ((TopologyEntity *)this_topology_entity)->shells(temp_shell_list);

    COPY_LIST_TO_ARRAY(temp_shell_list, *shell_list, *shell_list_size);
    
    return result;
}
  enum CubitStatus TopologyEntity_loops(void *this_topology_entity,  /* DLLoopList& */ void ***loop_list, int *loop_list_size)
{
  DLLoopList temp_loop_list;
  
  enum CubitStatus result = ((TopologyEntity *)this_topology_entity)->loops(temp_loop_list);

    COPY_LIST_TO_ARRAY(temp_loop_list, *loop_list, *loop_list_size);
    
    return result;
}
  enum CubitStatus TopologyEntity_chains(void *this_topology_entity, /* DLChainList& */ void ***chain_list, int *chain_list_size)
{
  DLChainList temp_chain_list;
  
  enum CubitStatus result = ((TopologyEntity *)this_topology_entity)->chains(temp_chain_list);

    COPY_LIST_TO_ARRAY(temp_chain_list, *chain_list, *chain_list_size);
    
    return result;
}
  enum CubitStatus TopologyEntity_ref_volumes(void *this_topology_entity, /* DLRefVolumeList& */ void ***ref_volume_list, int *ref_volume_list_size)
{
  DLRefVolumeList temp_ref_volume_list;
  
  enum CubitStatus result = ((TopologyEntity *)this_topology_entity)->ref_volumes(temp_ref_volume_list);

    COPY_LIST_TO_ARRAY(temp_ref_volume_list, *ref_volume_list, *ref_volume_list_size);
    
    return result;
}
  enum CubitStatus TopologyEntity_ref_faces(void *this_topology_entity, /* DLRefFaceList& */ void ***ref_face_list, int *ref_face_list_size)
{
  DLRefFaceList temp_ref_face_list;
  
  enum CubitStatus result = ((TopologyEntity *)this_topology_entity)->ref_faces(temp_ref_face_list);

    COPY_LIST_TO_ARRAY(temp_ref_face_list, *ref_face_list, *ref_face_list_size);
    
    return result;
}
  enum CubitStatus TopologyEntity_ref_edges(void *this_topology_entity, /* DLRefEdgeList& */ void ***ref_edge_list, int *ref_edge_list_size)
{
  DLRefEdgeList temp_ref_edge_list;
  
  enum CubitStatus result = ((TopologyEntity *)this_topology_entity)->ref_edges(temp_ref_edge_list);

    COPY_LIST_TO_ARRAY(temp_ref_edge_list, *ref_edge_list, *ref_edge_list_size);
    
    return result;
}
  enum CubitStatus TopologyEntity_ref_vertices(void *this_topology_entity, /* DLRefVertexList& */ void ***ref_vertex_list, int *ref_vertex_list_size)
{
  DLRefVertexList temp_ref_vertex_list;
  
  enum CubitStatus result = ((TopologyEntity *)this_topology_entity)->ref_vertices(temp_ref_vertex_list);

    COPY_LIST_TO_ARRAY(temp_ref_vertex_list, *ref_vertex_list, *ref_vertex_list_size);
    
    return result;
}
  enum CubitStatus TopologyEntity_co_volumes(void *this_topology_entity, /* DLCoVolumeList& */ void ***co_volume_list, int *co_volume_list_size)
{
  DLCoVolumeList temp_co_volume_list;
  
  enum CubitStatus result = ((TopologyEntity *)this_topology_entity)->co_volumes(temp_co_volume_list);

    COPY_LIST_TO_ARRAY(temp_co_volume_list, *co_volume_list, *co_volume_list_size);
    
    return result;
}
  enum CubitStatus TopologyEntity_co_faces(void *this_topology_entity, /* DLCoFaceList& */ void ***co_face_list, int *co_face_list_size)
{
  DLCoFaceList temp_co_face_list;
  
  enum CubitStatus result = ((TopologyEntity *)this_topology_entity)->co_faces(temp_co_face_list);

    COPY_LIST_TO_ARRAY(temp_co_face_list, *co_face_list, *co_face_list_size);
    
    return result;
}
  enum CubitStatus TopologyEntity_co_edges(void *this_topology_entity, /* DLCoEdgeList& */ void ***co_edge_list, int *co_edge_list_size)
{
  DLCoEdgeList temp_co_edge_list;
  
  enum CubitStatus result = ((TopologyEntity *)this_topology_entity)->co_edges(temp_co_edge_list);

    COPY_LIST_TO_ARRAY(temp_co_edge_list, *co_edge_list, *co_edge_list_size);
    
    return result;
}
  enum CubitStatus TopologyEntity_co_vertices(void *this_topology_entity, /* DLCoVertexList& */ void ***co_vertex_list, int *co_vertex_list_size)
{
  DLCoVertexList temp_co_vertex_list;
  
  enum CubitStatus result = ((TopologyEntity *)this_topology_entity)->co_vertices(temp_co_vertex_list);

    COPY_LIST_TO_ARRAY(temp_co_vertex_list, *co_vertex_list, *co_vertex_list_size);
    
    return result;
}

/* RefVertex **/ void *TopologyEntity_ref_vertex(void *this_topology_entity)
{
  return ((TopologyEntity *)this_topology_entity)->ref_vertex();
  
}
/* RefEdge **/ void *TopologyEntity_ref_edge(void *this_topology_entity)
{
  return ((TopologyEntity *)this_topology_entity)->ref_edge();
  
}
/* RefFace **/ void *TopologyEntity_ref_face(void *this_topology_entity)
{
  return ((TopologyEntity *)this_topology_entity)->ref_face();
  
}
/* RefVolume **/ void *TopologyEntity_ref_volume(void *this_topology_entity)
{
  return ((TopologyEntity *)this_topology_entity)->ref_volume();
  
}
/* Body **/ void *TopologyEntity_body(void *this_topology_entity)
{
  return ((TopologyEntity *)this_topology_entity)->body();
  
}

/* BridgeManager* */ void *TopologyEntity_bridge_manager(void *this_topology_entity)
{
  return ((TopologyEntity *)this_topology_entity)->bridge_manager();
  
}

  enum EntityType TopologyEntity_get_topology_bridge_type(void *this_topology_entity)
{
  return ((TopologyEntity *)this_topology_entity)->get_topology_bridge_type();
  
}
  
