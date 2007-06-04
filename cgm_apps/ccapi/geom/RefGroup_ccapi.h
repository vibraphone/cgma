#ifndef REFGROUP_CCAPI_H
#define REFGROUP_CCAPI_H

#include "CubitDefines.h"
#include "CubitVectorStruct.h"
#include "CubitBoxStruct.h"
#include "EntityType.h"

#ifdef __cplusplus
extern "C" {
#endif

  int RefGroup_maximum_dimension(void *this_group);
    //- This routine returns the maximum dimension of its owned subentities.
    //- The 'only' kludge is that if one of its subentities is a RefGroup, it
    //- must call 'maximum_dimension' on that entity instead of dimension().

  enum EntityType RefGroup_entity_type(void *this_group);
    //- return the type for this geometryEntity

  enum CubitStatus RefGroup_add_ref_entity_1(void *this_group, 
                                               /* RefEntity * */ void *ref_entity);

  enum CubitStatus RefGroup_add_ref_entity_2(void *this_group, 
                                               /* DLRefEntityList& */ void ***entity_list,
                                             int *entity_list_size);

  enum CubitStatus RefGroup_remove_ref_entity(void *this_group, 
                                                /* RefEntity * */ void *ref_entity,
                                              const enum CubitBoolean from_observable);
  
    //- add or remove one or more ref entities to/from this group

  void RefGroup_get_sub_entities(void *this_group,
                                   /* DLRefEntityList & */ void ***entity_list,
                                 int *entity_list_size);
    //- appends all ref entities owned by this entity and recurses
    //- down to dimension 0

  void RefGroup_is_mergeable_1(void *this_group, int val);

  int RefGroup_is_mergeable_2(void *this_group);
  
  int RefGroup_remove_all_ref_entities(void *this_group);
    //- remove all geometry entities in the group.  
    //- Returns number of items removed.

  void RefGroup_get_parent_ref_entities(void *this_group,
                                          /* DLRefEntityList& */ void ***,
                                        int *);
    //- appends all ref entities that own this to entity_list.
    //- Goes up just one dimension.

  void RefGroup_get_child_ref_entities(void *this_group,
                                         /* DLRefEntityList& */ void ***entity_list,
                                       int *entity_list_size);
    //- appends all immediate ref entities owned by this entity on entity_list
    //- Goes down just one dimension.
  
  void RefGroup_get_child_entities(void *this_group,
                                     /* DLCubitEntityList& */ void ***cub_entity_list,
                                   int *cub_entity_list_size);
    //- appends all immediate ref entities owned by this entity on entity_list
    //- Goes down just one dimension.

  void RefGroup_expand_group(void *this_group,
                               /* DLRefEntityList& */ void ***entity_list,
                             int *entity_list_size); 
    //- appends all the ref entities owned by this group.  It will go down
    //- until there are no ref-groups in the entity_list.

  struct CubitBoxStruct RefGroup_bounding_box(void *this_group);

  struct CubitVectorStruct RefGroup_center_point(void *this_group);
  
  void* RefGroup_get_address(void *this_group,
                             enum EntityType inputEntityType);
    //R void*
    //R- Returned void pointer
    //I inputEntityType
    //I- The input type to get the address of.
    //- This function returns a void pointer that points to the
    //- "appropriate" portion of this object.  The appropriate
    //- portion is determined by the input EntityType variable.
    //- Returns NULL if the input type and the type of "this"
    //- are not related by inheritance.
    //- Note: The RTTI capabilities encoded in these functions
    //-       are designed to work with any form of multiple 
    //-       inheritance, as well.  Multiple inheritance is what
    //-       necessitates having this function defined in every 
    //-       class in the hierarchy.
    //- Note: This function can also be used to merely check if
    //-       an object of one type is related to another type
    //-       through inheritance.  If a non-NULL pointer is
    //-       returned, then this is true.  

  int RefGroup_subtract(void *this_group,
                          /* RefGroup * */ void *group_to_subtract,
                          /* RefGroup * */ void *target_group);
    //- subtract group_to_subtract from this group
  
  int RefGroup_intersect(void *this_group,
                           /* RefGroup * */ void *other_group,
                           /* RefGroup * */ void *target_group);
    //- intersect other_group with this one
  
  int RefGroup_unite(void *this_group,
                       /* RefGroup * */ void *other_group,
                     /* RefGroup * */ void *target_group);
    //- unite other_group with this one

  int RefGroup_validate(void *this_group);
    //- Return number of problems detected, 0 if none.
  
  void RefGroup_draw (void *this_group, int color);
    //- draw the group's contained entities

  static enum CubitStatus RefGroup_delete_group (void *this_group,
                                                 /* RefGroup * */ void *group_ptr,
                                                 enum CubitBoolean propagate);
    //- deletes a specified group from the global group list
    //- if the boolean "propagate" is true, the groups sub_groups are 
    //- deleted also; as well as the sub_group's groups, etc..

  static void RefGroup_delete_all_groups (void *this_group);
    //- deletes all the groups from the model (except for 'picked_group' and 'drawn_group')

  static void RefGroup_get_contained_groups (void *this_group,
                                             /* RefGroup * */ void *group_ptr,
                                               /* DLRefGroupList & */ void ***contained_groups,
                                             int *contained_groups_size);
    //- gets the groups owned by group_ptr, as well as any other groups
    //- owned by these groups, etc..  Current group (group_ptr) is also 
    //- added to the list.

  static void RefGroup_get_groups_within_1(void *this_group,
                                             /* CubitEntity* */ void *cubit_entity_ptr, 
                                             /* DLRefGroupList & */ void ***groups_within,
                                           int *groups_within_size,
                                           const enum CubitBoolean recursive);
  static void RefGroup_get_groups_within_2(void *this_group,
                                             /* RefEntity* */ void *ref_entity_ptr, 
                                             /* DLRefGroupList & */ void ***groups_within,
                                           int *groups_within_size,
                                           const enum CubitBoolean recursive);
  static void RefGroup_get_groups_within_3(void *this_group,
                                             /* RefGroup* */ void *ref_group_ptr, 
                                             /* DLRefGroupList & */ void ***groups_within,
                                           int *groups_within_size,
                                           const enum CubitBoolean recursive);
    //- Finds those groups which contain the input entity, at any level
    //-  (i.e., if group 2 contains group 3 which contains the entity, then
    //-  group 2 and group 3 will be returned in the output list).

  enum CubitStatus RefGroup_notify_observer(void *this_group,
                                            /* CubitObservable * */ void *observable,
                                            const /* CubitEvent & */ void *observer_event,
                                            enum CubitBoolean from_observable);
    //- handle notify observer function

#ifdef __cplusplus
}
#endif

#endif
