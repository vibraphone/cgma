#ifndef REF_VERTEX_CCAPI_H
#define REF_VERTEX_CCAPI_H

#include "CubitDefines.h"
#include "EntityType.h"
#include "CubitVectorStruct.h"

#ifdef __cplusplus
extern "C" {
#endif

    /* Point * */ void *RefVertex_get_point_ptr(void *this_ref_vertex);
    /* R Point* */
    /* R- A pointer to the Point to which the current  */
    /* R- vertex points.  */
    /* - This function returns a pointer to the Point */
    /* - which the current vertex points to. */

  enum EntityType RefVertex_entity_type(void *this_ref_vertex);
    /* - return the type of this RefEntity/TopologyEntity */

  void *RefVertex_get_address(void *this_ref_vertex, enum EntityType inputEntityType);
    /* R void* */
    /* R- Returned void pointer */
    /* I inputenum EntityType */
    /* I- The input type to get the address of. */
    /* - This function returns a void pointer that points to the */
    /* - "appropriate" portion of this object.  The appropriate */
    /* - portion is determined by the input enum EntityType variable. */
    /* - Returns NULL if the input type and the type of "this" */
    /* - are not related by inheritance. */
    /* - Note: The RTTI capabilities encoded in these functions */
    /* -       are designed to work with any form of multiple  */
    /* -       inheritance, as well.  Multiple inheritance is what */
    /* -       necessitates having this function defined in every  */
    /* -       class in the hierarchy. */
    /* - Note: This function can also be used to merely check if */
    /* -       an object of one type is related to another type */
    /* -       through inheritance.  If a non-NULL pointer is */
    /* -       returned, then this is true. */

  enum EntityType RefVertex_get_child_ref_entity_type(void *this_ref_vertex);
    /* R enum EntityType */
    /* R- A type value. */
    /* - This function returns the type of the child RefEntity of   */
    /* - RefVertex -- and there isn't one! */
      
  enum EntityType RefVertex_get_parent_ref_entity_type(void *this_ref_vertex);
    /* R enum EntityType */
    /* R- A type value. */
    /* - This function returns the type of the parent RefEntity of   */
    /* - RefVertex, which is the type value of the RefEdge class. */
      
  enum EntityType RefVertex_get_grouping_entity_type(void *this_ref_vertex);
    /* R enum EntityType */
    /* R- The type of the corresponding GroupingEntity type. */
    /* - This function returns the type of the corresponding */
    /* - GroupingEntity . */
  
  enum EntityType RefVertex_get_sense_entity_type(void *this_ref_vertex);
    /* R enum EntityType */
    /* R- The type of SenseEntity associated with this object */
    /* - This function returns the type of SenseEntity associated with  */
    /* - this BasicTopologyEntity.  */

  enum EntityType RefVertex_get_topology_bridge_type(void *this_ref_vertex);
    /* R enum EntityType */
    /* R- The type of GeometryEntity associated with this object */
    /* - This function returns the type of GeometryEntity associated with  */
    /* - this BasicTopologyEntity.  */

  struct CubitVectorStruct RefVertex_coordinates (void *this_ref_vertex);
    /* R struct CubitVectorStruct */
    /* R- Contains the coordinate values {x y z} of this RefVertex */
    /* - Returns the spatial coordinates of this RefVertex. */

  struct CubitVectorStruct RefVertex_center_point(void *this_ref_vertex);
    /* - Return the "center" of this RefVertex (its coordinates) */

  int RefVertex_dimension(void *this_ref_vertex);
    /* - returns the geometric dimension of RefVertex'es.  */

    /* int number_mesh_elements(int dimension = -1, */
    /*                         ElementType type = INVALID_ELEMENT_TYPE); */
    /* - Returns number of mesh entities owned by this geometry entity */
    /* - of the specified dimension. If dimension == -1, the highest  */
    /* - dimension elements are returned. */

  void RefVertex_ref_edges_of_refvolume (void *this_ref_vertex, /* DLRefEdgeList& */ void ***edge_list, int *edge_list_size, 
                                           /* RefVolume * */ void *ref_volume_ptr );
    /* - Returns the subset of the list of ref_edges attached to this  */
    /* - RefVertex owned by the input RefVolume. */
    /* - */

  /* RefEdge * */ void *RefVertex_common_ref_edge(void *this_ref_vertex, /* RefVertex * */ void *other_vertex, 
                                       /* RefFace * */ void *owning_face);
    /* - returns an edge sharing the other vertex and owned by owning_face  */
    /* - (if non-NULL) */

  enum CubitBoolean RefVertex_about_spatially_equal(void *this_ref_vertex, /* RefVertex * */ void *ref_vertex_ptr_2,
                                                    double tolerance_factor,
                                                    enum CubitBoolean notify_ref_entity,
                                                    enum CubitBoolean force_merge);
    /* R-CUBIT_TRUE/CUBIT_FALSE */
    /* I RefVertex*, double, enum CubitBoolean */
    /* I- Second RefVertex to compare, Tolerance factor to for GEOMETRY_RESABS, */
    /* I- and flag for notifying compared RefEntities. */
    /* O enum CubitBoolean */
    /* O- if the two vertices are spatially equal within the GEOMETRY_RESABS* */
    /* -the tolerance_factor, then CUBIT_TRUE will be returned.  Otherwise */
    /* -CUBIT_FALSE is returned. */

  enum CubitBoolean RefVertex_move(void *this_ref_vertex,  struct CubitVectorStruct delta );
    /* - moves the vertex by the delta vector.  This is pretty */
    /* - scary but may save a model if you move it very slightly. */

 
#ifdef __cplusplus
}
#endif

#endif
