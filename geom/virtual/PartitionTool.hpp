//-------------------------------------------------------------------------
// Filename      : PartitionTool.hpp
//
// Purpose       : Functions for splitting geometry
//
// Special Notes :
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 12/13/98
//-------------------------------------------------------------------------

#ifndef PARTITION_TOOL_HPP
#define PARTITION_TOOL_HPP

class SenseEntity;
class PartitionEntity;
class PartitionCurve;
class PartitionSurface;
class Loop;
class Shell;

class RefEdge;
class RefVertex;
class Curve;
template <class X> class DLIList;
class PartitionCurve;

class CubitFacet;

#include "CubitDefines.h"
#include "RefEdge.hpp"

class PartitionTool
{
public:

  static PartitionTool* instance();
/*
  virtual CubitStatus partition( RefEdge* edge_ptr,
                                 RefVertex* split_point,
                                 RefEdge*& first_new_edge,
                                 RefEdge*& second_new_edge,
                                 CubitBoolean ignore_mesh = CUBIT_FALSE );
*/
  virtual RefVertex* partition( RefEdge* edge_ptr,
                                 const CubitVector& split_point,
                                 RefEdge*& first_new_edge,
                                 RefEdge*& second_new_edge,
                                 CubitBoolean ignore_mesh = CUBIT_FALSE );
  virtual RefVertex* partition( RefEdge* edge_ptr,
                                 const CubitVector& split_point,
                                 DLIList<RefEdge*>& result_edges,
                                 CubitBoolean ignore_mesh = CUBIT_FALSE );
  //R CubitStatus
  //R- CUBIT_SUCCESS/CUBIT_FAILURE
  //I edge_ptr
  //I- A pointer to the RefEdge to Partition.
  //I split_point
  //I- A pointer to the Vertex at the split location.
  //O first_new_edge / second_new_edge
  //O- Pointers to the new RefEdges created.
  //- Partition a RefEdge.

  virtual RefFace* insert_edge( RefFace* face_ptr,
                                DLIList<CubitVector*>& segments,
                                CubitBoolean       is_meshed,
                                DLIList<RefEdge*>& created_curves,
                                int level_of_recursion = 0,
                                const double *tolerance_length = NULL);

  virtual CubitStatus insert_edge( DLIList<RefFace*>& input_faces,
                                   DLIList<CubitVector*>& segments,
                                   DLIList<RefFace*>& result_faces,
                                   DLIList<RefEdge*>& created_edges,
                                   CubitBoolean do_split_curves = CUBIT_FALSE);


  virtual CubitStatus partition( RefFace* face_ptr,
                                 RefEdge* edge_ptr,
                                 RefFace*& new_face1,
                                 RefFace*& new_face2,
                                 CubitBoolean );

  virtual CubitStatus partition( RefFace* face_ptr,
                                 DLIList<RefEdge*>& split_edges,
                                 RefFace*& first_new_face,
                                 RefFace*& second_new_face,
                                 CubitBoolean ignore_mesh = CUBIT_FALSE );
  //R CubitStatus
  //R- CUBIT_SUCCESS/CUBIT_FAILURE
  //I face_ptr
  //I- A pointer to the RefFace to split.
  //I split_edges
  //I- The list of edges defining split location.
  //O first_new_face / second_new_face
  //O- The new RefFaces created.
  //- Partition a RefFace.

  virtual CubitStatus partition( RefVolume* vol_ptr,
                                 DLIList<RefFace*>& split_faces,
                                 RefVolume*& first_new_volume,
                                 RefVolume*& second_new_volume,
                                 CubitBoolean ignore_mesh = CUBIT_FALSE );
  virtual CubitStatus partition( RefVolume* vol_ptr,
                                 DLIList<CubitFacet*>& split_faces,
                                 RefVolume*& first_new_volume,
                                 RefVolume*& second_new_volume,
                                 DLIList<RefFace*>& new_surfaces,
                                 CubitBoolean ignore_mesh = CUBIT_FALSE );
  //R CubitStatus
  //I vol_ptr
  //I- A pointer to the RefVolume to split.
  //I split_faces
  //I- The list of faces defining the split location
  //O first_new_volume/second_new_volume
  //O- The new RefVolumes created
  //- Partition a RefVolume.

  CubitStatus destroy_volume_partition( RefVolume* partition );

  //CubitStatus partition( RefEdge* edge_ptr,
  //                       DLIList<RefVertex*>& split_points,
  //                       DLIList<RefEdge*>& new_edges,
  //                       CubitBoolean ignore_mesh = CUBIT_FALSE );
  CubitStatus partition( RefEdge* edge_ptr,
                         DLIList<CubitVector*>& split_points,
                         DLIList<RefEdge*>& new_edges,
                         CubitBoolean ignore_mesh = CUBIT_FALSE );
  CubitStatus partition( RefEdge* edge_ptr,
                         DLIList<CubitVector*>& split_points,
                         DLIList<RefVertex*>& new_vertices,
                         DLIList<RefEdge*>& new_edges,
                         CubitBoolean ignore_mesh = CUBIT_FALSE );
  //R CubitStatus
  //R- CUBIT_SUCCESS/CUBIT_FAILURE
  //I edge_ptr
  //I- A pointer to the RefEdge to Partition.
  //I split_points
  //I- A list of Vertices at which to split the curve.
  //O new_edges
  //O- The new RefEdges created.
  //- Partition a RefEdge.

  CubitStatus partition_face_by_curves( RefFace* face_ptr,
                         const DLIList<Curve*>& split_curves,
                         DLIList<RefFace*>& result_faces,
                         CubitBoolean do_split_curves = CUBIT_FALSE,
                         DLIList<RefEdge*>* new_edges = NULL,
                         CubitBoolean ignore_mesh = CUBIT_FALSE );
  //- Similar to the next function but takes curves as input

  CubitStatus partition( RefFace* face_ptr,
                         const DLIList<RefEdge*>& split_edges,
                         DLIList<RefFace*>& result_faces,
                         CubitBoolean do_split_curves = CUBIT_FALSE,
                         DLIList<RefEdge*>* new_edges = NULL,
                         CubitBoolean ignore_mesh = CUBIT_FALSE );
  //R CubitStatus
  //R- CUBIT_SUCCESS/CUBIT_FAILURE
  //I face_ptr
  //I- A pointer to the RefFace to partition
  //I split_edges
  //I- A list of RefEdges with which to partition the RefFace
  //O result_faces
  //O- The RefFaces resulting from the splitting operation
  //I split_curves
  //I- If true, spatially compare all RefVertices in the split_edge
  //I- list with edges in that list and on the RefFace, and partition
  //I- RefEdges as necessary.
  //O new_edges
  //O- If split_curves in CUBIT_TRUE, and this is non-null, a list
  //O- of new RefEdges created by partitioning will be passed back.
  //I ignore_mesh
  //I- If true, operator will not fail if the passed RefFace is meshed,
  //I- or if split_curves in true and any of the RefEdges are meshed.
  //I- NOTE: If this is passed as CUBIT_TRUE, it is the responsibility
  //I-       of the caller to remove the mesh from hidden entities.
  //- Do multiple partitioning of a RefFace.

  RefEdge* unpartition( RefEdge* first_partition_edge,
                        RefEdge* second_partition_edge );
  //R RefEdge*
  //R- The RefEdge restored to the model.
  //I first_partition_edge / second_partition_edge
  //I- Pointers to RefEdges with PartitionCurve geometry.
  //O free_vertex
  //O- A pointer to the freed vertex.
  //I keep_dead
  //I- If true, the passed-in RefEdges are not destroyed, but the
  //I- attached PartitionCurve is.
  //- Undo a partitioning operation.

  RefFace* unpartition( RefFace* first_partition_face,
                        RefFace* second_partition_face );
  //R RefFace*
  //R- The RefFace restored to the model.
  //I first_partition_face / second_partition_face
  //I- Pointers to RefFaces with PartitionSurface geometry.
  //O free_edges
  //O- A list of edges freed by this operation.
  //I keep_dead
  //I- If true, the passed-in RefFaces are not destroyed, but the
  //I- attached PartitionSurface is.
  //- Undo a partitioning operation.

   RefVolume* unpartition( RefVolume* first_partition_volume,
                           RefVolume* second_partition_volume );
  //R RefVolume*
  //R- The RefVolume restored to the model.
  //I first_partition_volume / second_partition_volume
  //I- Pointers to RefVolumes with PartitionLump geometry.
  //O free_faces
  //O- A list of faces freed by this operation.
  //- Undo a partitioning operation.

  RefVolume* unpartition( DLIList<RefVolume*>& partition_volumes );


  virtual CubitStatus unpartitionAll( DLIList<RefEdge*>& partition_edges,
                              DLIList<RefEdge*>& restored_edges );
  //R CubitStatus
  //R- CUBIT_SUCCESS/CUBIT_FAILURE
  //I partition_edges
  //I- A list of RefEdges with PartitionCurve geometry.
  //O restored_edges
  //O- A list of RefEdges restored to the model.
  //- Remove a bunch of partitions.
  //-
  //- Note: Any vertices freed by this operation will be deleted.

  virtual CubitStatus unpartitionAll( DLIList<RefFace*>& partition_faces,
                              DLIList<RefFace*>& restored_faces );
  //R CubitStatus
  //R- CUBIT_SUCCESS/CUBIT_FAILURE
  //I partition_faces
  //I- A list of RefFaces with PartitionSurfaces geometry.
  //O restored_faces
  //O- A list of RefFaces restored to the model.
  //- Remove a bunch of partitions.
  //-
  //- Note: any RefEdges freed by this operation and thier child
  //-       RefVertices will be deleted.


  virtual CubitStatus unpartitionAll( DLIList<RefVolume*>& partition_vols,
                              DLIList<RefVolume*>& restored_vols );

  RefVertex* make_point_curve( RefFace* ref_face,
                               const CubitVector& position );

  CubitStatus make_point_curves( RefFace* ref_face,
                                 DLIList<CubitVector> &position,
                                 DLIList<RefVertex*> &new_vertices);

  virtual ~PartitionTool();

  virtual CubitStatus can_remove( RefVertex* vertex );
  virtual CubitStatus can_remove( RefEdge* edge );
  virtual CubitStatus can_remove( RefFace* face );

protected:

  static PartitionTool* instance_;

  PartitionTool();

  DLIList<Surface*>* group_merged_surfaces( DLIList<RefFace*>& face_list,
                                            int& num_result_sets );


private:

  void notify_partition( DLIList<RefEntity*> &partitioning_entities,
                         BasicTopologyEntity *first_partitioned_entity,
                         BasicTopologyEntity *second_partitioned_entity,
                         BasicTopologyEntity *old_entity );
    // notify observers that the old_entity was partitioned



  CubitStatus commonGroupingEntities( BasicTopologyEntity* first_bte_ptr,
                                      BasicTopologyEntity* secnd_bte_ptr,
                                      GroupingEntity*& first_gpe_ptr,
                                      GroupingEntity*& second_gpe_ptr ) const;
    //R CubitStatus
    //R- CUBIT_SUCCESS/CUBIT_FAILURE
    //I first_bte_ptr/second_bte_ptr
    //I- Pointers to two BTEs of the same type.
    //O first_gpe_ptr/second_gpe_ptr
    //O- Child grouping entities of first_bte_ptr and second_bte_ptr,
    //O- respectively.
    //- This function searches for a grouping entity from each
    //- of the BTEs for which those GroupingEntites have atleast one
    //- common child.
    //- If there is only one pair of grouping entities which meet this
    //- criteria, they are passed back, and CUBIT_SUCCESS is returned.
    //- If there is more than one pair, an arbitrary pair is passed
    //- back, and CUBIT_FAILURE is returned.
    //- If no common GpEs are found, both pointers are passed back as
    //- NULL, and CUBIT_FAILURE is returned.

};

//-------------------------------------------------------------------------
// Purpose       : Singleton class
//
// Special Notes :
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 12/13/98
//-------------------------------------------------------------------------
inline PartitionTool* PartitionTool::instance()
{
  if( instance_ == NULL ) instance_ = new PartitionTool();
  return instance_;
}

#endif

