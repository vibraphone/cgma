#include "MergeTool_ccapi.h"

#include "MergeTool.hpp"

#include "Body.hpp"
#include "RefVolume.hpp"
#include "RefFace.hpp"
#include "RefEdge.hpp"
#include "RefVertex.hpp"
#include "RefEntity.hpp"
#include "DLBodyList.hpp"
#include "DLRefVolumeList.hpp"
#include "DLRefFaceList.hpp"
#include "DLRefEdgeList.hpp"
#include "DLRefVertexList.hpp"

#include "copy_defines.h"

#define MTI MergeTool::instance()

static /* MergeTool* */ void *MergeTool_instance()
{
  return MTI->instance();
}

  /*  Returns a static pointer to unique instance of this class. */
   
enum CubitBoolean MergeTool_contains_merged_entities_1( /* DLBodyList & */ void ***bodies, int *bodies_size )
{

  DLBodyList temp_bodies  ;
  COPY_ARRAY_TO_LIST(*bodies, *bodies_size, temp_bodies);
  
  return MTI->contains_merged_entities(temp_bodies);
}

enum CubitBoolean MergeTool_contains_merged_entities_2( /* DLRefEntityList & */ void ***ref_entities, int *ref_entities_size)
{
  DLRefEntityList temp_ref_entities;
  COPY_ARRAY_TO_LIST(*ref_entities, *ref_entities_size, temp_ref_entities);
  
  return MTI->contains_merged_entities(temp_ref_entities);
}

enum CubitBoolean MergeTool_entity_merged( /* TopologyEntity * */ void *entity )
{
  TopologyEntity *temp_entity = (TopologyEntity *) entity;
  return MTI->entity_merged(temp_entity);
}

enum CubitStatus MergeTool_merge_all_bodies()
{
  return MTI->merge_all_bodies();
}

enum CubitStatus MergeTool_merge_bodies( /* DLBodyList& */ void ***refbody_list, int *refbody_list_size )
{
  DLBodyList temp_refbody_list;
  COPY_ARRAY_TO_LIST(*refbody_list, *refbody_list_size, temp_refbody_list);
  
  return MTI->merge_bodies(temp_refbody_list);
}

enum CubitStatus MergeTool_merge_volumes( /* DLRefVolumeList& */ void ***vol_list, int *vol_list_size,
                                          enum CubitBoolean print_info)
{
  DLRefVolumeList temp_vol_list;
  COPY_ARRAY_TO_LIST(*vol_list, *vol_list_size, temp_vol_list);
  
  return MTI->merge_volumes(temp_vol_list, print_info);
}

enum CubitStatus MergeTool_merge_all_reffaces()
{
  return MTI->merge_all_reffaces();
}

enum CubitStatus MergeTool_merge_reffaces( /* DLRefFaceList& */ void ***refface_list, int *refface_list_size,
                                           enum CubitBoolean force_merge,
                                           enum CubitBoolean print_info)
{
  DLRefFaceList temp_refface_list;
  COPY_ARRAY_TO_LIST(*refface_list, *refface_list_size, temp_refface_list);
  
  return MTI->merge_reffaces(temp_refface_list, force_merge, print_info);
}


enum CubitStatus MergeTool_merge_all_refedges()
{
  return MTI->merge_all_refedges();
}

enum CubitStatus MergeTool_merge_refedges( /* DLRefEdgeList& */ void ***refedge_list, int *refedge_list_size,
                                           enum CubitBoolean force_merge,
                                           enum CubitBoolean should_clean_out,
                                           enum CubitBoolean print_info)
{
  DLRefEdgeList temp_refedge_list;
  COPY_ARRAY_TO_LIST(*refedge_list, *refedge_list_size, temp_refedge_list);
  
  return MTI->merge_refedges(temp_refedge_list, force_merge, should_clean_out, print_info);
}


    
enum CubitStatus MergeTool_merge_all_refvertices()
{
  return MTI->merge_all_refvertices();
}

enum CubitStatus MergeTool_merge_refvertices( /* DLRefVertexList& */ void ***refvertex_list, int *refvertex_list_size,
                                              enum CubitBoolean force_merge,
                                              enum CubitBoolean print_info)
{
  DLRefVertexList temp_refvertex_list;
  COPY_ARRAY_TO_LIST(*refvertex_list, *refvertex_list_size, temp_refvertex_list);
  
  return MTI->merge_refvertices(temp_refvertex_list, force_merge, print_info);
}

enum CubitStatus MergeTool_merge_entities( /* DLRefEntityList& */ void ***entity_list, int *entity_list_size,
                                           enum CubitBoolean force_merge,
                                           enum CubitBoolean should_clean_out,
                                           enum CubitBoolean print_info)
{
  DLRefEntityList temp_entity_list;
  COPY_ARRAY_TO_LIST(*entity_list, *entity_list_size, temp_entity_list);
  
  return MTI->merge_entities(temp_entity_list, force_merge, should_clean_out, print_info);
}


    
enum CubitStatus MergeTool_unmerge( /* DLRefEntityList& */ void ***)
{
  DLRefEntityList temp_list;
  
  return MTI->unmerge(temp_list);
}

static enum CubitBoolean MergeTool_merge_has_occured()
{
  return MTI->merge_has_occured();
}

static void MergeTool_set_merge_occurance( enum CubitBoolean t_or_f )
{
  MTI->set_merge_occurance(t_or_f);
}

static void MergeTool_group_results( enum CubitBoolean t_or_f )
{
  MTI->group_results(t_or_f);
}

void MergeTool_compare_notify(/* RefEntity * */ void *entity, enum EventType event)
{
  RefEntity *temp_entity = (RefEntity *) entity;
  MTI->compare_notify(temp_entity, event);
}

