#include "EntityType.h"
#include "CubitDefines.h"

#include "RefEntity_ccapi.h"
#include "RefEntity.hpp"
#include "DLRefEntityList.hpp"

#include "copy_defines.h"

  /* RefEntity* */ void *RefEntity_get_by_name(void *,
                                                        /* const CubitString& */ const char *name) 
{
  return RefEntity::get_by_name(CubitString(name));
}

    /* CubitString */ const char *RefEntity_entity_name_1(void *this_ref_entity) 
{
  RefEntity *temp_ref_entity = (RefEntity *) this_ref_entity;
  return temp_ref_entity->entity_name().c_str();
}

  enum CubitStatus RefEntity_entity_name_2(void *this_ref_entity, /* CubitString */ const char *name) 
{
  RefEntity *temp_ref_entity = (RefEntity *) this_ref_entity;
  return temp_ref_entity->entity_name(name);
}

  void RefEntity_entity_names(void *this_ref_entity, /* DLCubitStringList & */ void ***names, int *names_size) 
{
  RefEntity *temp_ref_entity = (RefEntity *) this_ref_entity;
  DLCubitStringList temp_names;
  
  temp_ref_entity->entity_names(temp_names);

  COPY_LIST_TO_ARRAY(temp_names, *names, *names_size);
  
}

  int RefEntity_num_names(void *this_ref_entity) 
{
  RefEntity *temp_ref_entity = (RefEntity *) this_ref_entity;
  return temp_ref_entity->num_names();
}

  enum CubitStatus RefEntity_generate_default_name (void *this_ref_entity,  /* CubitString & */ char **name ) 
{
  RefEntity *temp_ref_entity = (RefEntity *) this_ref_entity;
  CubitString temp_name;
  
  CubitStatus result = temp_ref_entity->generate_default_name(temp_name);

  *name = new char[temp_name.length()];
  strcpy(*name, temp_name.c_str());
  
  return result;
}

  enum CubitStatus RefEntity_assign_default_name(void *this_ref_entity,  enum CubitBoolean user_setting) 
{
  RefEntity *temp_ref_entity = (RefEntity *) this_ref_entity;
  return temp_ref_entity->assign_default_name(user_setting);
}

  enum CubitStatus RefEntity_remove_entity_name(void *this_ref_entity, /* CubitString const & */ const char *name) 
{
  RefEntity *temp_ref_entity = (RefEntity *) this_ref_entity;
  CubitString temp_name(name);
  
  return temp_ref_entity->remove_entity_name(temp_name);
}

  enum CubitStatus RefEntity_remove_entity_names(void *this_ref_entity) 
{
  RefEntity *temp_ref_entity = (RefEntity *) this_ref_entity;
  return temp_ref_entity->remove_entity_names();
}

  void RefEntity_merge_entity_names(void *this_ref_entity, /* RefEntity * */ void *dead_entity) 
{
  RefEntity *temp_ref_entity = (RefEntity *) this_ref_entity;
  RefEntity *temp_dead_entity = (RefEntity *) dead_entity;
  
  temp_ref_entity->merge_entity_names(temp_dead_entity);
}

  void RefEntity_switch_entity_names(void *this_ref_entity, /* RefEntity * */ void *other_entity) 
{
  RefEntity *temp_ref_entity = (RefEntity *) this_ref_entity;
  RefEntity *temp_other_entity = (RefEntity *) other_entity;

  temp_ref_entity->switch_entity_names(temp_other_entity);
}
  
  void RefEntity_marked_1(void *this_ref_entity, int value) 
{
  RefEntity *temp_ref_entity = (RefEntity *) this_ref_entity;
  temp_ref_entity->marked(value);
}

  int  RefEntity_marked_2(void *this_ref_entity) 
{
  RefEntity *temp_ref_entity = (RefEntity *) this_ref_entity;
  return temp_ref_entity->marked ();
}

  void RefEntity_is_mergeable_1(void *this_ref_entity, int val) 
{
  RefEntity *temp_ref_entity = (RefEntity *) this_ref_entity;
  temp_ref_entity->is_mergeable(val);
}

  int  RefEntity_is_mergeable_2(void *this_ref_entity) 
{
  RefEntity *temp_ref_entity = (RefEntity *) this_ref_entity;
  return temp_ref_entity->is_mergeable();
}

  int RefEntity_children_mergeable(void *this_ref_entity) 
{
  RefEntity *temp_ref_entity = (RefEntity *) this_ref_entity;
  return temp_ref_entity->children_mergeable();
}
  
  int RefEntity_dimension(void *this_ref_entity) 
{
  RefEntity *temp_ref_entity = (RefEntity *) this_ref_entity;
  return temp_ref_entity->dimension();
}
  
  enum EntityType RefEntity_get_child_ref_entity_type(void *this_ref_entity) 
{
  RefEntity *temp_ref_entity = (RefEntity *) this_ref_entity;
  return temp_ref_entity->get_child_ref_entity_type();
}
  
  enum EntityType RefEntity_get_parent_ref_entity_type(void *this_ref_entity) 
{
  RefEntity *temp_ref_entity = (RefEntity *) this_ref_entity;
  return temp_ref_entity->get_parent_ref_entity_type();
}
  
  void RefEntity_get_child_ref_entities(void *this_ref_entity,
                                          /* DLRefEntityList & */ void ***entity_list, int *entity_list_size) 
{
  RefEntity *temp_ref_entity = (RefEntity *) this_ref_entity;
  DLRefEntityList temp_entity_list;
  
  temp_ref_entity->get_child_ref_entities(temp_entity_list);

  COPY_LIST_TO_ARRAY(temp_entity_list, *entity_list, *entity_list_size);
  
}
  
  void RefEntity_get_all_child_ref_entities(void *this_ref_entity,
                                              /* DLRefEntityList & */ void ***entity_list, int *entity_list_size) 
{
  RefEntity *temp_ref_entity = (RefEntity *) this_ref_entity;
  DLRefEntityList temp_entity_list;
  
  temp_ref_entity->get_all_child_ref_entities(temp_entity_list);

  COPY_LIST_TO_ARRAY(temp_entity_list, *entity_list, *entity_list_size);
  
}
  
  void RefEntity_gather_bdy_entities(void *,
                                              /* DLRefEntityList & */ void ***entity_list, int *entity_list_size, 
                                              /* DLRefEntityList & */ void ***bdy_list, int *bdy_list_size) 
{
  DLRefEntityList temp_bdy_list;
  DLRefEntityList temp_entity_list;
  COPY_ARRAY_TO_LIST(*entity_list, *entity_list_size, temp_entity_list);
  
  RefEntity::gather_bdy_entities(temp_entity_list, temp_bdy_list);

  COPY_LIST_TO_ARRAY(temp_bdy_list, *bdy_list, *bdy_list_size);
}

  void RefEntity_get_parent_ref_entities(void *this_ref_entity,
                                           /* DLRefEntityList & */ void ***entity_list, int *entity_list_size) 
{
  RefEntity *temp_ref_entity = (RefEntity *) this_ref_entity;
  DLRefEntityList temp_entity_list;
  
  temp_ref_entity->get_parent_ref_entities(temp_entity_list);

  COPY_LIST_TO_ARRAY(temp_entity_list, *entity_list, *entity_list_size);
}
  
  void RefEntity_change_to_parent_ref_entities(void *,
                                                        /* DLRefEntityList & */ void ***ancestors, int *ancestors_size) 
{
  DLRefEntityList temp_ancestors;
  
  COPY_ARRAY_TO_LIST(*ancestors, *ancestors_size, temp_ancestors);

  RefEntity::change_to_parent_ref_entities(temp_ancestors);

  COPY_LIST_TO_ARRAY(temp_ancestors, *ancestors, *ancestors_size);
}
  
    /* RefEntity * */ void *RefEntity_join_1(void *this_ref_entity,  /* RefEntity * */ void *ref_entity_2,
                                             /* DLRefEntityList & */ void ***join_set, int *join_set_size) 
{
  RefEntity *temp_ref_entity = (RefEntity *) this_ref_entity;
  DLRefEntityList temp_join_set;
  RefEntity *temp_ref_entity_2 = (RefEntity *) ref_entity_2;
  
  RefEntity *result = temp_ref_entity->join(temp_ref_entity_2, temp_join_set);

  COPY_LIST_TO_ARRAY(temp_join_set, *join_set, *join_set_size);

  return result;
}
  /* RefEntity * */ void *RefEntity_join_2(void *,
                                                  /* DLRefEntityList & */ void ***ref_entities, int *ref_entities_size, 
                                                  /* DLRefEntityList & */ void ***join_set, int *join_set_size) 
{
  DLRefEntityList temp_join_set;
  DLRefEntityList temp_ref_entities;
  COPY_ARRAY_TO_LIST(*ref_entities, *ref_entities_size, temp_ref_entities);
  
  RefEntity *result = RefEntity::join(temp_ref_entities, temp_join_set);

  COPY_LIST_TO_ARRAY(temp_join_set, *join_set, *join_set_size);
  
  return result;
}
  
    /* RefEntity * */ void *RefEntity_meet_1(void *this_ref_entity,  /* RefEntity * */ void *ref_entity_2,
                                             /* DLRefEntityList & */ void ***join_set, int *join_set_size) 
{
  RefEntity *temp_ref_entity = (RefEntity *) this_ref_entity;
  DLRefEntityList temp_join_set;
  RefEntity *temp_ref_entity_2 = (RefEntity *) ref_entity_2;
  
  RefEntity *result = temp_ref_entity->meet(temp_ref_entity_2, temp_join_set);

  COPY_LIST_TO_ARRAY(temp_join_set, *join_set, *join_set_size);

  return result;
}
  /* RefEntity * */ void *RefEntity_meet_2(void *,
                                                  /* DLRefEntityList & */ void ***ref_entities, int *ref_entities_size, 
                                                  /* DLRefEntityList & */ void ***join_set, int *join_set_size) 
{
  DLRefEntityList temp_join_set;
  DLRefEntityList temp_ref_entities;
  COPY_ARRAY_TO_LIST(*ref_entities, *ref_entities_size, temp_ref_entities);
  
  RefEntity *result = RefEntity::meet(temp_ref_entities, temp_join_set);

  COPY_LIST_TO_ARRAY(temp_join_set, *join_set, *join_set_size);
  
  return result;
}

  int RefEntity_valence(void *this_ref_entity, /* RefEntity * */ void *parent) 
{
  RefEntity *temp_ref_entity = (RefEntity *) this_ref_entity;
  RefEntity *temp_parent = (RefEntity *) parent;
  
  return temp_ref_entity->valence(temp_parent);
}

  enum CubitBoolean RefEntity_is_child(void *this_ref_entity, /* RefEntity * */ void *entity) 
{
  RefEntity *temp_ref_entity = (RefEntity *) this_ref_entity;
  RefEntity *temp_entity = (RefEntity *) entity;

  return temp_ref_entity->is_child(temp_entity);
}

  enum CubitBoolean RefEntity_is_parent(void *this_ref_entity, /* RefEntity * */ void *entity) 
{
  RefEntity *temp_ref_entity = (RefEntity *) this_ref_entity;
  RefEntity *temp_entity = (RefEntity *) entity;

  return temp_ref_entity->is_parent(temp_entity);
}

  struct CubitVectorStruct RefEntity_center_point(void *this_ref_entity) 
{
  RefEntity *temp_ref_entity = (RefEntity *) this_ref_entity;
  return temp_ref_entity->center_point();
}
  
  double RefEntity_measure(void *this_ref_entity) 
{
  RefEntity *temp_ref_entity = (RefEntity *) this_ref_entity;
  return temp_ref_entity->measure();
}
  
    /* CubitString */ const char *RefEntity_measure_label(void *this_ref_entity) 
{
  RefEntity *temp_ref_entity = (RefEntity *) this_ref_entity;
  return temp_ref_entity->measure_label().c_str();
}
  
  int RefEntity_validate(void *this_ref_entity) 
{
  RefEntity *temp_ref_entity = (RefEntity *) this_ref_entity;
  return temp_ref_entity->validate();
}
  
  enum EntityType RefEntity_entity_type(void *this_ref_entity) 
{
  RefEntity *temp_ref_entity = (RefEntity *) this_ref_entity;
  return temp_ref_entity->entity_type();
}
  
  void* RefEntity_get_address(void *this_ref_entity, enum EntityType inputEntityType) 
{
  RefEntity *temp_ref_entity = (RefEntity *) this_ref_entity;
  return temp_ref_entity->get_address(inputEntityType);
}
  
void RefEntity_notify(void *this_ref_entity, /* RefEntity * */ void *partner, enum EventType event) 
{
  RefEntity *temp_ref_entity = (RefEntity *) this_ref_entity;
  RefEntity *temp_partner = (RefEntity *) partner;

  temp_ref_entity->notify(temp_partner, event);
}
  
void RefEntity_add_compare_data(void *this_ref_entity, /* RefEntity * */ void *partner) 
{
  RefEntity *temp_ref_entity = (RefEntity *) this_ref_entity;
  RefEntity *temp_partner = (RefEntity *) partner;

  temp_ref_entity->add_compare_data(temp_partner);
}
  
void RefEntity_remove_compare_data(void *this_ref_entity) 
{
  RefEntity *temp_ref_entity = (RefEntity *) this_ref_entity;
  temp_ref_entity->remove_compare_data();
}

  void RefEntity_get_related_entity_list(void *this_ref_entity, enum EntityType related_entity_type,
                                           /* DLRefEntityList & */ void ***entity_list, int *entity_list_size) 
{
  RefEntity *temp_ref_entity = (RefEntity *) this_ref_entity;
  DLRefEntityList temp_entity_list;
  
  temp_ref_entity->get_related_entity_list(related_entity_type, temp_entity_list);

  COPY_LIST_TO_ARRAY(temp_entity_list, *entity_list, *entity_list_size);
  
}

void RefEntity_set_id(void *this_ref_entity, int i) 
{
  RefEntity *temp_ref_entity = (RefEntity *) this_ref_entity;
  temp_ref_entity->set_id(i);
}

void RefEntity_draw(void *this_ref_entity, int color) 
{
  RefEntity *temp_ref_entity = (RefEntity *) this_ref_entity;
  temp_ref_entity->draw(color);
}
