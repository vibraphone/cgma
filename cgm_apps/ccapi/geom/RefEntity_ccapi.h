#ifndef REFENTITY_CCAPI_H
#define REFENTITY_CCAPI_H

#include "CubitDefines.h"
#include "EntityType.h"
#include "CubitVectorStruct.h"

#ifdef __cplusplus
extern "C" {
#endif

    /*  Name-related functions */
  /* RefEntity* */ void *RefEntity_get_by_name(void *this_ref_entity,
                                                        /* const CubitString& */ const char *name);
  
    /* CubitString */ const char *RefEntity_entity_name_1(void *this_ref_entity);
  
  enum CubitStatus RefEntity_entity_name_2(void *this_ref_entity, /* CubitString */ const char *name);
  
  void RefEntity_entity_names(void *this_ref_entity, /* DLCubitStringList & */ void ***names, int *names_size);
  
  int RefEntity_num_names(void *this_ref_entity);

  enum CubitStatus RefEntity_generate_default_name (void *this_ref_entity,  /* CubitString & */ char **name );

  enum CubitStatus RefEntity_assign_default_name(void *this_ref_entity,  enum CubitBoolean user_setting);

  enum CubitStatus RefEntity_remove_entity_name(void *this_ref_entity, /* CubitString const & */ const char *name);

  enum CubitStatus RefEntity_remove_entity_names(void *this_ref_entity);

  void RefEntity_merge_entity_names(void *this_ref_entity, /* RefEntity* */ void *dead_entity);

  void RefEntity_switch_entity_names(void *this_ref_entity, /* RefEntity * */ void *other_entity);

    /* - Functions related to naming of RefEntity objects */
  
  void RefEntity_marked_1(void *this_ref_entity, int value);

  int  RefEntity_marked_2(void *this_ref_entity);

    /* - Generic flag for temporary use by algorithms. */
    /* - Value is volatile and may change unexpectedly. */

  void RefEntity_is_mergeable_1(void *this_ref_entity, int val);

  int  RefEntity_is_mergeable_2(void *this_ref_entity);

    /* - Set/Get the mergeable status for this entity. */

  int RefEntity_children_mergeable(void *this_ref_entity);

    /* - Get whether all child entities are mergeable. */
  
  int RefEntity_dimension(void *this_ref_entity); 
    /*  Returns the geometric dimension of the entity.  */
    /*  vertex == 0, edge == 1, etc. */
  
  enum EntityType RefEntity_get_child_ref_entity_type(void *this_ref_entity);
    /* R EntityType */
    /* R- A type value. */
    /* - This function returns the type of the child RefEntity for  */
    /* - "this" type of RefEntity.  If there is no child RefEntity type, */
    /* - for example for RefVertex, then InvalidEntity_TYPE is returned. */
  
  enum EntityType RefEntity_get_parent_ref_entity_type(void *this_ref_entity);
    /* R EntityType */
    /* R- A type value. */
    /* - This function returns the type of the parent RefEntity for  */
    /* - "this" type of RefEntity.  If there is no parent RefEntity type, */
    /* - for example for RefVolume, then InvalidEntity_TYPE is returned. */
  
  void RefEntity_get_child_ref_entities(void *this_ref_entity,
                                          /* DLRefEntityList & */ void ***entity_list, int *entity_list_size); 
    /* - Appends all immediate (child) RefEntities owned by this RefEntity to  */
    /* - entity_list. (The query goes down just one dimension.) */
  
  void RefEntity_get_all_child_ref_entities(void *this_ref_entity,
                                              /* DLRefEntityList & */ void ***entity_list, int *entity_list_size);

    /* - Appends all child RefEntities owned by this entity to entity_list. */
    /* - (The query recurses all the way down to RefEntities of dimension 0). */
  
  void RefEntity_gather_bdy_entities(void *this_ref_entity,
                                              /* DLRefEntityList & */ void ***entity_list, int *entity_list_size, 
                                              /* DLRefEntityList & */ void ***bdy_list, int *bdy_list_size);

    /*  Gather the boundary entities of the entity_list into the bdy_list.  */
    /*  Entities appear once in bdy_list, and will not appear in bdy_list */
    /*  if they are already in the entity_list. */
    /*  Uses listMark. */

  void RefEntity_get_parent_ref_entities(void *this_ref_entity,
                                           /* DLRefEntityList & */ void ***entity_list,
                                         int *entity_list_size);

    /* - Appends all RefEntities that own this (parent RefEntities) to  */
    /* - entity_list. */
    /* - (The query goes up just one dimension. For example, if this is a */
    /* - vertex, the resulting list contains only RefEdges). */
  
  void RefEntity_change_to_parent_ref_entities(void *this_ref_entity,
                                                        /* DLRefEntityList & */ void ***ancestors,
                                               int *ancestors_size);

    /* - Modify the input list to contain the list of RefEntities that are */
    /* - the parents of each of the RefEntities in the original list. */
  
    /* RefEntity * */ void *RefEntity_join_1(void *this_ref_entity,  /* RefEntity* */ void *ref_entity_2,
                                           /* DLRefEntityList & */ void ***join_set, int *join_set_size);

/* RefEntity * */ void *RefEntity_join_2(void *this_ref_entity,
                                                /* DLRefEntityList & */ void ***ref_entities, int *ref_entities_size, 
                                                /* DLRefEntityList & */ void ***join_set, int *join_set_size);

    /*  RefEntity* join( RefEntity* ref_entity_2 ); */
    /* - Computes the geometric "join" of elements (elements on the list or  */
    /* - this and ref_entity_2).   */
    /* - Definition "Join" = The lowest dimensional entitity that */
    /* - is a higher dimensional ancestor of them all.  */
    /* - Note join could be one of the entities itself, NULL,  */
    /* - or multiple elements. */
    /* - E.g. The join of a vertex and a curve containing the vertex is  */
    /* - the curve. The join of two curves of a split cylinder is  */
    /* - both containing surfaces. */
    /* - The join of two entities in separate, unmerged volumes is null. */
    /* - Returns the first element of the join_set, or NULL if set is empty. */
  
    /* RefEntity * */ void *RefEntity_meet_1(void *this_ref_entity,  /* RefEntity* */ void *ref_entity_2,
                                           /* DLRefEntityList & */ void ***join_set, int *join_set_size);

/* RefEntity * */ void *RefEntity_meet_2(void *this_ref_entity,
                                                /* DLRefEntityList & */ void ***ref_entities, int *ref_entities_size, 
                                                /* DLRefEntityList & */ void ***join_set, int *join_set_size);

    /* - like join, except returns the lower order entities common to the input  */
    /* - entities */

  int RefEntity_valence(void *this_ref_entity, /* RefEntity * */ void *parent);

    /* - the valence of this entity with respect to parent (absolute */
    /* - valence if parent is null) */

  enum CubitBoolean RefEntity_is_child(void *this_ref_entity, /* RefEntity * */ void *entity);

  enum CubitBoolean RefEntity_is_parent(void *this_ref_entity, /* RefEntity * */ void *entity);

    /* - Return TRUE if this is entity, or a direct child (parent) of entity. */
  
/*    void common_draw_label(int color, const int label_style); */
/*      //- Common code used for drawing labels. Actually, the entire */
/*      //- draw_label function could be implemented here if the */
/*      //- is_*_labeling_on() functions were or better implemented. */

  struct CubitVectorStruct RefEntity_center_point(void *this_ref_entity);
    /* - Return the approximate (spatial) center of this RefEntity */
  
  double RefEntity_measure(void *this_ref_entity);

    /* - A generic geometric extent function. */
    /* - Returns volume for RefVolumes, area for RefFaces, length for RefEdge, */
    /* - and 1.0 for RefVertices */
    /* - A RefGroup calculates the maximum dimension of its contained */
    /* - entities and returns the sum of the measures() of all entities */
    /* - of that dimension. */
    /* - Default return value is 0.0 for all other entities. */
  
  /* CubitString */ const char *RefEntity_measure_label(void *this_ref_entity);

    /* - Returns the type of measure: (volume, area, length, or N/A) */
  
  int RefEntity_validate(void *this_ref_entity);

    /* - Perform checks to see if entity valid. */
  
enum EntityType RefEntity_entity_type(void *this_ref_entity);
    /* - Return the type of the entity */
  
  void* RefEntity_get_address(void *this_ref_entity, enum EntityType inputEntityType);

    /* R void* */
    /* R- Returned void pointer */
    /* I inputEntityType */
    /* I- The input type to get the address of. */
    /* - This function returns a void pointer that points to the */
    /* - "appropriate" portion of this object.  The appropriate */
    /* - portion is determined by the input EntityType variable. */
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
  
  void RefEntity_notify(void *this_ref_entity, /* RefEntity * */ void *partner,
                        enum EventType event);

    /* R void */
    /* I partner */
    /* I- The merge partner for this object */
    /* I event */
    /* I- The type of event */
    /* - This function takes actions depending on the type of event it */
    /* - is notified of. */
    /* -   COMPARISON_FOUND: */
    /* -     Make temporary TDCompare objects and attach to "this" and */
    /* -     the "partner" object. */
  
  void RefEntity_add_compare_data(void *this_ref_entity, /* RefEntity * */ void *partner);

    /* R void */
    /* I partner */
    /* I- The compare partner for this object */
    /* - This function makes the connection between the two RefEntities, */
    /* - this and partner. At the end of this function the two entities */
    /* - would know who they compare with. */
  
  void RefEntity_remove_compare_data(void *this_ref_entity);

    /* R void */
    /* - This function clears the compare related temporary data. */

/* ========  Change Code by RY of Cat,  5/7/99 11:08:12 AM  ======== */
  void RefEntity_get_related_entity_list(void *this_ref_entity, enum EntityType related_entity_type,
                                           /* DLRefEntityList & */ void ***entity_list, int *entity_list_size);

    /* - to parse group in <ref_entity> commands */
/* ========  Change End by RY of Cat,  5/7/99 11:08:12 AM  ======== */

  void RefEntity_set_id(void *this_ref_entity, int i);

    /* - set the id of this entity to i */

  void RefEntity_draw(void *this_ref_entity, int color);

#ifdef __cplusplus
}
#endif

#endif

