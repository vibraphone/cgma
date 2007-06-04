#ifndef MERGETOOL_CCAPI_H
#define MERGETOOL_CCAPI_H

#include "CubitDefines.h"

#ifdef __cplusplus
extern "C" {
#endif

  static /* MergeTool* */ void *MergeTool_instance();
  /*  Returns a static pointer to unique instance of this class. */
   
  enum CubitBoolean MergeTool_contains_merged_entities_1( /* DLBodyList & */ void ***bodies, int *bodies_size );
  /* - Tests the entities in the list of bodies to see if they  */
  /* - have been merged with other entities and result in */
  /* - of multiple volumes.  This returns CUBIT_TRUE at the first */
  /* - lower entity that shares two volumes. */
   
  enum CubitBoolean MergeTool_contains_merged_entities_2( /* DLRefEntityList & */ void ***ref_entities, int *ref_entities_size);
  /* - Tests the entities in the list to see if they or their descendents */
  /* - have been merged with other entities */
   
  enum CubitBoolean MergeTool_entity_merged( /* TopologyEntity * */ void *entity );
  /* - Takes a RefEntity, makes sure it is not a volume, and */
  /* - checks to see if it has more than one volume. */
   
   enum CubitStatus MergeTool_merge_all_bodies();
  /* - Compare all RefFaces, RefEdges, and RefVertices in the */
  /* - model and merge the matches */
   
  enum CubitStatus MergeTool_merge_bodies( /* DLBodyList& */ void *** refbody_list, int *refbody_list_size );
  /* - Compare all RefFaces, RefEdges, and RefVertices in  */
  /* - refbody_list and merge the matches */
   
  enum CubitStatus MergeTool_merge_volumes( /* DLRefVolumeList& */ void *** vol_list, int *vol_list_size,
                              enum CubitBoolean print_info);
  /* - Compare all RefFaces, RefEdges, and RefVertices in  */
  /* - vol_list and merge the matches */
   
   enum CubitStatus MergeTool_merge_all_reffaces    ();
  /* - Compare all RefFaces in the model and merge the matches */
   
  enum CubitStatus MergeTool_merge_reffaces( /* DLRefFaceList& */ void *** refface_list, int *refface_list_size,
                               enum CubitBoolean force_merge,
                               enum CubitBoolean print_info);
  /* - Compare all input RefFaces and merge the matches */
   
   enum CubitStatus MergeTool_merge_all_refedges();
  /* - Compare all RefEdges in the model and merge the matches */
   
  enum CubitStatus MergeTool_merge_refedges( /* DLRefEdgeList& */ void *** refedge_list, int *refedge_list_size,
                               enum CubitBoolean force_merge,
                               enum CubitBoolean should_clean_out,
                               enum CubitBoolean print_info);
  /* - Compare all input RefEdges and merge the matches */
  /* - BE CAREFUL with the should_clean_out flag.  If you set */
  /* - it to false, then YOU (the caller) are responsible for */
  /* - cleaning out the deactivated geometry. */
   
   enum CubitStatus MergeTool_merge_all_refvertices();
  /* - Compare all RefVertices in the model and merge the matches */
   
  enum CubitStatus MergeTool_merge_refvertices( /* DLRefVertexList& */ void *** refvertex_list, int *refvertex_list_size,
                                  enum CubitBoolean force_merge,
                                  enum CubitBoolean print_info);
  /* - Compare all input RefVertices and merge the matches */
   
  enum CubitStatus MergeTool_merge_entities( /* DLRefEntityList& */ void *** entity_list, int *entity_list_size,
                               enum CubitBoolean force_merge,
                               enum CubitBoolean should_clean_out,
                               enum CubitBoolean print_info);
  /* - merge the entities in list; asserts if they're not all the same type */
   
  enum CubitStatus MergeTool_unmerge( /* DLRefEntityList& */ void *** );
  /* - Takes a list of entities.  NOTE this list should have been */
  /* - screened so that only surfaces, curves, and vertices are */
  /* - in it.  If not, it will assert. */
   
/*    enum CubitStatus has_merged_parents(); */
  /* R enum CubitStatus */
  /* R- CUBIT_SUCCESS/FAILURE */
  /* - This function determines if this entity has merged parents */
   
   static enum CubitBoolean MergeTool_merge_has_occured();
   static void MergeTool_set_merge_occurance( enum CubitBoolean t_or_f );
  /* - Sets/Gets the mergeCalled flag.  This signifies */
  /* - that currently a merge operation has been called for  */
  /* - the model. */

   static void MergeTool_group_results( enum CubitBoolean t_or_f );
  /* - tells us to group the results or not. */
  
  void MergeTool_compare_notify(/* RefEntity * */ void *entity, enum EventType event);
  /* - notifies MergeTool about comparisons found and put on ref entities */
  

#ifdef __cplusplus
}
#endif

#endif
