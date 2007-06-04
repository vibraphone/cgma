#ifndef __CGM_GEOM_HPP__
#define __CGM_GEOM_HPP__

#include "TSTTB_SNL.h"

#ifdef __cplusplus

extern "C" 
{
#else
  typedef int bool;
#endif

  typedef void* TSTTG_GentityHandle;
  typedef const void* TSTTG_CGentityHandle;
  typedef void* TSTTG_GentitySetHandle;
  typedef const void* TSTTG_CGentitySetHandle;
  typedef void* TSTTG_TagHandle;
  typedef void* TSTTG_Instance;
  typedef void* TSTTG_GentityIterator;
  typedef const void* TSTTG_CGentityIterator;

#define ARRAY_IN_DECL(a, b) \
  a *b, const int b ## _size  

#define CONST_ARRAY_IN_DECL(a, b) \
  const a *b, const int b ## _size  

#define ARRAY_INOUT_DECL(a, b) \
  a **b, int *b ## _allocated, int *b ## _size  

#define ARRAY_IN(b) \
  b, b ## _size  

#define ARRAY_INOUT(b) \
  b, b ## _allocated, b ## _size  

#define RETURN(a) {TSTTG_LAST_ERROR.error_type = a; return a;}
#define TSTTG_processError(a, b) {sprintf(TSTTG_LAST_ERROR.description, b); TSTTG_LAST_ERROR.error_type = a;}

  enum TSTTG_GentityType {       
    TSTTG_VERTEX, 
    TSTTG_EDGE, 
    TSTTG_FACE, 
    TSTTG_REGION,
    TSTTG_ALL_TYPES
  };

  enum TSTTG_TagValueType {
    TSTTG_INTEGER = 0,
    TSTTG_DOUBLE,
    TSTTG_ENTITY_HANDLE,
    TSTTG_BYTES
  };

  extern struct TSTTB_Error TSTTG_LAST_ERROR;

  enum TSTTB_ErrorType
  TSTTG_ctor(TSTTG_Instance *instance);

  enum TSTTB_ErrorType
  TSTTG_dtor(TSTTG_Instance instance);

  enum TSTTB_ErrorType
  TSTTG_createEntSet (TSTTG_Instance instance,
                      /*in*/ bool isList,
                      /*out*/ TSTTG_GentitySetHandle *entity_set);

  enum TSTTB_ErrorType
  TSTTG_destroyEntSet (TSTTG_Instance instance,
                       /*in*/ TSTTG_GentitySetHandle entity_set);

  bool
  TSTTG_isList (TSTTG_Instance instance,
                /*in*/ TSTTG_CGentitySetHandle entity_set);

  int
  TSTTG_getNumEntSets (TSTTG_Instance instance,
                       /*in*/ TSTTG_CGentitySetHandle entity_set,
                       /*in*/ int num_hops);

  enum TSTTB_ErrorType
  TSTTG_getEntSets (TSTTG_Instance instance,
                    /*in*/ TSTTG_CGentitySetHandle entity_set,
                    /*in*/ const int num_hops,
                    /*inout*/ ARRAY_INOUT_DECL(TSTTG_GentitySetHandle, 
                                               contained_entity_set_handles));

  enum TSTTB_ErrorType
  TSTTG_addEntToSet (TSTTG_Instance instance,
                     /*in*/ TSTTG_CGentityHandle entity_handle,
                     /*inout*/ TSTTG_GentitySetHandle *entity_set);

  enum TSTTB_ErrorType
  TSTTG_rmvEntFromSet (TSTTG_Instance instance,
                       /*in*/ TSTTG_CGentityHandle entity_handle,
                       /*inout*/ TSTTG_GentitySetHandle *entity_set);

  enum TSTTB_ErrorType
  TSTTG_addEntArrToSet (TSTTG_Instance instance,
                        /*in*/ ARRAY_IN_DECL(TSTTG_GentityHandle, entity_handles),
                        /*inout*/ TSTTG_GentitySetHandle *entity_set);

  enum TSTTB_ErrorType
  TSTTG_rmvEntArrFromSet (TSTTG_Instance instance,
                          /*in*/ ARRAY_IN_DECL(TSTTG_GentityHandle, entity_handles),
                          /*inout*/ TSTTG_GentitySetHandle *entity_set);

  enum TSTTB_ErrorType
  TSTTG_addEntSet (TSTTG_Instance instance,
                   /*in*/ TSTTG_CGentitySetHandle entity_set_to_add,
                   /*inout*/ TSTTG_GentitySetHandle *entity_set_handle);

  enum TSTTB_ErrorType
  TSTTG_rmvEntSet (TSTTG_Instance instance,
                   /*in*/ TSTTG_CGentitySetHandle entity_set_to_remove,
                   /*inout*/ TSTTG_GentitySetHandle *entity_set_handle);

  bool
  TSTTG_isEntContained (TSTTG_Instance instance,
                        /*in*/ TSTTG_CGentitySetHandle containing_entity_set,
                        /*in*/ TSTTG_CGentityHandle contained_entity);

  bool
  TSTTG_isEntSetContained (TSTTG_Instance instance,
                           /*in*/ TSTTG_CGentitySetHandle containing_entity_set,
                           /*in*/ TSTTG_CGentitySetHandle contained_entity_set);

  enum TSTTB_ErrorType
  TSTTG_addPrntChld (TSTTG_Instance instance,
                     /*inout*/ TSTTG_GentitySetHandle *parent_entity_set,
                     /*inout*/ TSTTG_GentitySetHandle *child_entity_set
                     );

  enum TSTTB_ErrorType
  TSTTG_rmvPrntChld (TSTTG_Instance instance,
                     /*inout*/ TSTTG_GentitySetHandle *parent_entity_set,
                     /*inout*/ TSTTG_GentitySetHandle *child_entity_set
                     );

  bool
  TSTTG_isChildOf (TSTTG_Instance instance,
                   /*in*/ TSTTG_CGentitySetHandle parent_entity_set,
                   /*in*/ TSTTG_CGentitySetHandle child_entity_set);

  int
  TSTTG_getNumChld (TSTTG_Instance instance,
                    /*in*/ TSTTG_CGentitySetHandle entity_set,
                    /*in*/ const int num_hops);

  int
  TSTTG_getNumPrnt (TSTTG_Instance instance,
                    /*in*/ TSTTG_CGentitySetHandle entity_set,
                    /*in*/ const int num_hops);

  enum TSTTB_ErrorType
  TSTTG_getChldn (TSTTG_Instance instance,
                  /*in*/ TSTTG_CGentitySetHandle from_entity_set,
                  /*in*/ const int num_hops,
                  /*inout*/ ARRAY_INOUT_DECL(TSTTG_GentitySetHandle, entity_set_handles));

  enum TSTTB_ErrorType
  TSTTG_getPrnts (TSTTG_Instance instance,
                  /*in*/ TSTTG_CGentitySetHandle from_entity_set,
                  /*in*/ const int num_hops,
                  /*inout*/ ARRAY_INOUT_DECL(TSTTG_GentitySetHandle, entity_set_handles));

  enum TSTTB_ErrorType
  TSTTG_gentitysetGetGentitiesOfType (TSTTG_Instance instance,
                                      /*in*/ TSTTG_CGentitySetHandle set_handle,
                                      /*in*/ const TSTTG_GentityType gentity_type,
                                      /*out*/ ARRAY_INOUT_DECL(TSTTG_GentityHandle, gentity_handles));

  int
  TSTTG_gentitysetGetNumberGentitiesOfType (TSTTG_Instance instance,
                                            /*in*/ TSTTG_CGentitySetHandle set_handle,
                                            /*in*/ const TSTTG_GentityType gentity_type);

  enum TSTTB_ErrorType
  TSTTG_gentityGetType (TSTTG_Instance instance,
                        /*in*/ ARRAY_IN_DECL(TSTTG_CGentityHandle, gentity_handles),
                        /*inout*/ ARRAY_INOUT_DECL(TSTTG_GentityType, gtype));

  enum TSTTB_ErrorType
  TSTTG_gentityGetAdjacencies (TSTTG_Instance instance,
                               /*in*/ TSTTG_CGentityHandle gentity_handle,
                               /*in*/ const int to_dimension,
                               /*inout*/ ARRAY_INOUT_DECL(TSTTG_GentityHandle, adj_gentities));

  enum TSTTB_ErrorType
  TSTTG_gentityGet2OAdjacencies (TSTTG_Instance instance,
                                 /*in*/ TSTTG_CGentityHandle gentity_handle,
                                 /*in*/ const int bridge_dimension,
                                 /*in*/ const int to_dimension,
                                 /*out*/ ARRAY_INOUT_DECL(TSTTG_GentityHandle, adjacent_gentities));

  enum TSTTB_ErrorType
  TSTTG_gentityIsAdjacent (TSTTG_Instance instance,
                           /*in*/ TSTTG_CGentityHandle gentity_handle1,
                           /*in*/ TSTTG_CGentityHandle gentity_handle2,
                           /*out*/ bool *are_adjacent);

  enum TSTTB_ErrorType
  TSTTG_Move (TSTTG_Instance instance,
              /*inout*/ TSTTG_GentityHandle *geom_entity,
              /*in*/ const double x,
              /*in*/ const double y,
              /*in*/ const double z);

  enum TSTTB_ErrorType
  TSTTG_Rotate (TSTTG_Instance instance,
                /*inout*/ TSTTG_GentityHandle *geom_entity,
                /*in*/ const double angle,
                /*in*/ const double axis_normal_x,
                /*in*/ const double axis_normal_y,
                /*in*/ const double axis_normal_z);

  enum TSTTB_ErrorType
  TSTTG_Reflect (TSTTG_Instance instance,
                 /*inout*/ TSTTG_GentityHandle *geom_entity,
                 /*in*/ const double plane_normal_x,
                 /*in*/ const double plane_normal_y,
                 /*in*/ const double plane_normal_z);

  enum TSTTB_ErrorType
  TSTTG_Scale (TSTTG_Instance instance,
               /*inout*/ TSTTG_GentityHandle *geom_entity,
               /*in*/ const double scale_x,
               /*in*/ const double scale_y,
               /*in*/ const double scale_z);

  enum TSTTB_ErrorType
  TSTTG_Copy (TSTTG_Instance instance,
              /*in*/ TSTTG_GentityHandle geom_entity,
              /*out*/ TSTTG_GentityHandle *geom_entity2);

  enum TSTTB_ErrorType
  TSTTG_SweepAboutAxis (TSTTG_Instance instance,
                        /*in*/ TSTTG_GentityHandle geom_entity,
                        /*in*/ const double angle,
                        /*in*/ const double axis_normal_x,
                        /*in*/ const double axis_normal_y,
                        /*in*/ const double axis_normal_z,
                        /*out*/ TSTTG_GentityHandle *geom_entity2);

  enum TSTTB_ErrorType
  TSTTG_Delete (TSTTG_Instance instance,
                /*in*/ TSTTG_GentityHandle geom_entity);

  enum TSTTB_ErrorType
  TSTTG_gLoad (TSTTG_Instance instance,
               /*in*/ const char *name,
               /*in*/ CONST_ARRAY_IN_DECL(char*, options));

  enum TSTTB_ErrorType
  TSTTG_gSave (TSTTG_Instance instance,
               /*in*/ const char *name,
               /*in*/ CONST_ARRAY_IN_DECL(char*, options));

  enum TSTTB_ErrorType
  TSTTG_gentityClosestPoint (TSTTG_Instance instance,
    /*in*/ ARRAY_IN_DECL(TSTTG_CGentityHandle, gentity_handles),
    /*in*/ CONST_ARRAY_IN_DECL(double, near_coordinates),
    /*out*/ ARRAY_INOUT_DECL(double, on_coordinates));

  enum TSTTB_ErrorType
  TSTTG_gentityNormal (TSTTG_Instance instance,
                       /*in*/ ARRAY_IN_DECL(TSTTG_CGentityHandle, gentity_handles),
                       /*in*/ CONST_ARRAY_IN_DECL(double, coordinates),
                       /*out*/ ARRAY_INOUT_DECL(double, normals));

  enum TSTTB_ErrorType
  TSTTG_gentityClosestPointAndNormal (TSTTG_Instance instance,
    /*in*/ ARRAY_IN_DECL(TSTTG_CGentityHandle, gentity_handles),
    /*in*/ CONST_ARRAY_IN_DECL(double, near_coordinates),
    /*out*/ ARRAY_INOUT_DECL(double, on_coordinates),
    /*out*/ ARRAY_INOUT_DECL(double, normals));

  enum TSTTB_ErrorType
  TSTTG_gentityTangent (TSTTG_Instance instance,
                        /*in*/ ARRAY_IN_DECL(TSTTG_CGentityHandle, gentity_handles),
                        /*in*/ CONST_ARRAY_IN_DECL(double, coordinates),
                        /*out*/ ARRAY_INOUT_DECL(double, tangents));

  enum TSTTB_ErrorType
  TSTTG_gentityBoundingBox (TSTTG_Instance instance,
                            /*in*/ ARRAY_IN_DECL(TSTTG_GentityHandle, gentity_handles),
                            /*out*/ ARRAY_INOUT_DECL(double, min_corner),
                            /*out*/ ARRAY_INOUT_DECL(double, max_corner));

  enum TSTTB_ErrorType
  TSTTG_getGvertexCoordinates (TSTTG_Instance instance,
                               /*in*/ ARRAY_IN_DECL(TSTTG_CGentityHandle, gentity_handles),
                               /*out*/ ARRAY_INOUT_DECL(double, coordinates));

  int
  TSTTG_getGnormalSense (TSTTG_Instance instance,
                         /*in*/ TSTTG_CGentityHandle gface,
                         /*in*/ TSTTG_CGentityHandle gregion);

  int
  TSTTG_getGtangentSense (TSTTG_Instance instance,
                          /*in*/ TSTTG_CGentityHandle gedge,
                          /*in*/ TSTTG_CGentityHandle gface);

  int
  TSTTG_getGvertexTangentSense (TSTTG_Instance instance,
                                /*in*/ TSTTG_CGentityHandle gedge,
                                /*in*/ TSTTG_CGentityHandle gvertex1,
                                /*in*/ TSTTG_CGentityHandle gvertex2);

  enum TSTTB_ErrorType
  TSTTG_createTag (TSTTG_Instance instance,
                   /*in*/ const char *tag_name,
                   /*in*/ const int tag_size,
                   /*in*/ const TSTTG_TagValueType tag_type,
                   /*out*/ TSTTG_TagHandle *tag_handle);

  enum TSTTB_ErrorType
  TSTTG_destroyTag (TSTTG_Instance instance,
                    /*in*/ const TSTTG_TagHandle tag_handle,
                    /*in*/ const bool forced);

  const char *
  TSTTG_getTagName (TSTTG_Instance instance,
                    /*in*/ const TSTTG_TagHandle tag_handle);

  int
  TSTTG_getTagSizeValues (TSTTG_Instance instance,
                          /*in*/ const TSTTG_TagHandle tag_handle);

  int
  TSTTG_getTagSizeBytes (TSTTG_Instance instance,
                         /*in*/ const TSTTG_TagHandle tag_handle);

  TSTTG_TagHandle
  TSTTG_getTagHandle (TSTTG_Instance instance,
                      /*in*/ const char *tag_name);

  TSTTG_TagValueType
  TSTTG_getTagType (TSTTG_Instance instance,
                    /*in*/ const TSTTG_TagHandle tag_handle);

  enum TSTTB_ErrorType
  TSTTG_getArrData (TSTTG_Instance instance,
                    /*in*/ ARRAY_IN_DECL(TSTTG_CGentityHandle, entity_handles),
                    /*in*/ const TSTTG_TagHandle tag_handle,
                    /*inout*/ ARRAY_INOUT_DECL(char, tag_value));

  enum TSTTB_ErrorType
  TSTTG_getIntArrData (TSTTG_Instance instance,
                       /*in*/ ARRAY_IN_DECL(TSTTG_CGentityHandle, entity_handles),
                       /*in*/ const TSTTG_TagHandle tag_handle,
                       /*inout*/ ARRAY_INOUT_DECL(int, tag_value));

  enum TSTTB_ErrorType
  TSTTG_getDblArrData (TSTTG_Instance instance,
                       /*in*/ ARRAY_IN_DECL(TSTTG_CGentityHandle, entity_handles),
                       /*in*/ const TSTTG_TagHandle tag_handle,
                       /*inout*/ ARRAY_INOUT_DECL(double, tag_value));

  enum TSTTB_ErrorType
  TSTTG_getEHArrData (TSTTG_Instance instance,
                      /*in*/ ARRAY_IN_DECL(TSTTG_CGentityHandle, entity_handles),
                      /*in*/ const TSTTG_TagHandle tag_handle,
                      /*inout*/ ARRAY_INOUT_DECL(TSTTG_GentityHandle, tag_value));

  enum TSTTB_ErrorType
  TSTTG_setArrData (TSTTG_Instance instance,
                    /*in*/ ARRAY_IN_DECL(TSTTG_GentityHandle, entity_handles),
                    /*in*/ const TSTTG_TagHandle tag_handle,
                    /*in*/ CONST_ARRAY_IN_DECL(char, tag_values));

  enum TSTTB_ErrorType
  TSTTG_setIntArrData (TSTTG_Instance instance,
                       /*in*/ ARRAY_IN_DECL(TSTTG_GentityHandle, entity_handles),
                       /*in*/ const TSTTG_TagHandle tag_handle,
                       /*in*/ CONST_ARRAY_IN_DECL(int, tag_values));

  enum TSTTB_ErrorType
  TSTTG_setDblArrData (TSTTG_Instance instance,
                       /*in*/ ARRAY_IN_DECL(TSTTG_GentityHandle, entity_handles),
                       /*in*/ const TSTTG_TagHandle tag_handle,
                       /*in*/ CONST_ARRAY_IN_DECL(double, tag_values));

  enum TSTTB_ErrorType
  TSTTG_setEHArrData (TSTTG_Instance instance,
                      /*in*/ ARRAY_IN_DECL(TSTTG_GentityHandle, entity_handles),
                      /*in*/ const TSTTG_TagHandle tag_handle,
                      /*in*/ ARRAY_IN_DECL(TSTTG_CGentityHandle, tag_values));

  enum TSTTB_ErrorType
  TSTTG_rmvArrTag (TSTTG_Instance instance,
                   /*in*/ ARRAY_IN_DECL(TSTTG_GentityHandle, entity_handles),
                   /*in*/ const TSTTG_TagHandle tag_handle);

  enum TSTTB_ErrorType
  TSTTG_getData (TSTTG_Instance instance,
                 /*in*/ TSTTG_CGentityHandle entity_handle,
                 /*in*/ const TSTTG_TagHandle tag_handle,
                 /*inout*/ ARRAY_INOUT_DECL(char, tag_value));

  int
  TSTTG_getIntData (TSTTG_Instance instance,
                    /*in*/ TSTTG_CGentityHandle entity_handle,
                    /*in*/ const TSTTG_TagHandle tag_handle);

  double
  TSTTG_getDblData (TSTTG_Instance instance,
                    /*in*/ TSTTG_CGentityHandle entity_handle,
                    /*in*/ const TSTTG_TagHandle tag_handle);

  TSTTG_GentityHandle
  TSTTG_getEHData (TSTTG_Instance instance,
                   /*in*/ TSTTG_CGentityHandle entity_handle,
                   /*in*/ const TSTTG_TagHandle tag_handle);

  enum TSTTB_ErrorType
  TSTTG_setData (TSTTG_Instance instance,
                 /*in*/ TSTTG_GentityHandle entity_handle,
                 /*in*/ const TSTTG_TagHandle tag_handle,
                 /*in*/ CONST_ARRAY_IN_DECL(char, tag_value));

  enum TSTTB_ErrorType
  TSTTG_setIntData (TSTTG_Instance instance,
                    /*in*/ TSTTG_GentityHandle entity_handle,
                    /*in*/ const TSTTG_TagHandle tag_handle,
                    /*in*/ const int tag_value);

  enum TSTTB_ErrorType
  TSTTG_setDblData (TSTTG_Instance instance,
                    /*in*/ TSTTG_GentityHandle entity_handle,
                    /*in*/ const TSTTG_TagHandle tag_handle,
                    /*in*/ const double tag_value);

  enum TSTTB_ErrorType
  TSTTG_setEHData (TSTTG_Instance instance,
                   /*in*/ TSTTG_GentityHandle entity_handle,
                   /*in*/ const TSTTG_TagHandle tag_handle,
                   /*in*/ TSTTG_CGentityHandle tag_value);

  enum TSTTB_ErrorType
  TSTTG_getAllTags (TSTTG_Instance instance,
                    /*in*/ TSTTG_CGentityHandle entity_handle,
                    /*inout*/ ARRAY_INOUT_DECL(void*, tag_handles));

  enum TSTTB_ErrorType
  TSTTG_rmvTag (TSTTG_Instance instance,
                /*in*/ TSTTG_GentityHandle entity_handle,
                /*in*/ const TSTTG_TagHandle tag_handle);

  enum TSTTB_ErrorType
  TSTTG_Brick (TSTTG_Instance instance,
               /*in*/ const double x,
               /*in*/ const double y,
               /*in*/ const double z,
               /*out*/ TSTTG_GentityHandle *geom_entity);

  enum TSTTB_ErrorType
  TSTTG_Cylinder (TSTTG_Instance instance,
                  /*in*/ const double height,
                  /*in*/ const double major_rad,
                  /*in*/ const double minor_rad,
                  /*out*/ TSTTG_GentityHandle *geom_entity);

  enum TSTTB_ErrorType
  TSTTG_Torus (TSTTG_Instance instance,
               /*in*/ const double major_rad,
               /*in*/ const double minor_rad,
               /*out*/ TSTTG_GentityHandle *geom_entity);

  enum TSTTB_ErrorType
  TSTTG_Unite (TSTTG_Instance instance,
               /*in*/ ARRAY_IN_DECL(TSTTG_GentityHandle, geom_entities),
               /*out*/ TSTTG_GentityHandle *geom_entity);

  enum TSTTB_ErrorType
  TSTTG_Subtract (TSTTG_Instance instance,
                  /*in*/ TSTTG_GentityHandle blank,
                  /*in*/ TSTTG_GentityHandle tool,
                  /*out*/ TSTTG_GentityHandle *geom_entity);

  enum TSTTB_ErrorType
  TSTTG_Section (TSTTG_Instance instance,
                 /*inout*/ TSTTG_GentityHandle *geom_entity,
                 /*in*/ const double plane_normal_x,
                 /*in*/ const double plane_normal_y,
                 /*in*/ const double plane_normal_z,
                 /*in*/ const double offset,
                 /*in*/ const bool reverse,
                 /*out*/ TSTTG_GentityHandle *geom_entity2);

  enum TSTTB_ErrorType
  TSTTG_Imprint (TSTTG_Instance instance,
                 /*inout*/ ARRAY_IN_DECL(TSTTG_GentityHandle, gentity_handles));
  
  enum TSTTB_ErrorType
  TSTTG_Merge (TSTTG_Instance instance,
               /*inout*/ ARRAY_IN_DECL(TSTTG_GentityHandle, gentity_handles),
               const double tolerance);
  
  enum TSTTB_ErrorType
  TSTTG_getGtolerance (TSTTG_Instance instance,
                       /*out*/ double *relative_tolerance,
                       /*out*/ double *absolute_tolerance
                       );

  enum TSTTB_ErrorType
  TSTTG_getGentityTolerance (TSTTG_Instance instance,
                             /*in*/ ARRAY_IN_DECL(TSTTG_CGentityHandle, gentity_handles),
                             /*out*/ ARRAY_INOUT_DECL(double, relative_tolerances),
                             /*out*/ ARRAY_INOUT_DECL(double, absolute_tolerances));

  enum TSTTB_ErrorType
  TSTTG_subtract (TSTTG_Instance instance,
                  /*in*/ TSTTG_CGentitySetHandle entity_set_1,
                  /*in*/ TSTTG_CGentitySetHandle entity_set_2,
                  /*out*/ TSTTG_GentitySetHandle *result_entity_set);

  enum TSTTB_ErrorType
  TSTTG_intersect (TSTTG_Instance instance,
                   /*in*/ TSTTG_CGentitySetHandle entity_set_1,
                   /*in*/ TSTTG_CGentitySetHandle entity_set_2,
                   /*out*/ TSTTG_GentitySetHandle *result_entity_set);

  enum TSTTB_ErrorType
  TSTTG_unite (TSTTG_Instance instance,
               /*in*/ TSTTG_CGentitySetHandle entity_set_1,
               /*in*/ TSTTG_CGentitySetHandle entity_set_2,
               /*out*/ TSTTG_GentitySetHandle *result_entity_set);

  enum TSTTB_ErrorType
  TSTTG_gentityIteratorInit (TSTTG_Instance instance,
                             /*in*/ const int gentity_dimension,
                             /*out*/ TSTTG_GentityIterator *gentity_iterator);

  bool
  TSTTG_gentityIteratorNext (TSTTG_Instance instance,
                             /*inout*/ TSTTG_GentityIterator *gentity_iterator,
                             /*out*/ TSTTG_GentityHandle *gentity_handle
                             );

  enum TSTTB_ErrorType
  TSTTG_gentityIteratorReset (TSTTG_Instance instance,
                              /*inout*/ TSTTG_GentityIterator *gentity_iterator
                              );

  enum TSTTB_ErrorType
  TSTTG_gentityIteratorDelete (TSTTG_Instance instance,
                               /*in*/ TSTTG_CGentityIterator Gentity_dim_iterator);

  int
  TSTTG_gentityIsParametric (TSTTG_Instance instance,
                             /*in*/ TSTTG_CGentityHandle gentity_handle);

  enum TSTTB_ErrorType
  TSTTG_gentityUvToXyz (TSTTG_Instance instance,
                        /*in*/ ARRAY_IN_DECL(TSTTG_CGentityHandle, gentity_handles),
                        /*in*/ CONST_ARRAY_IN_DECL(double, uv),
                        /*out*/ ARRAY_INOUT_DECL(double, coordinates));

  enum TSTTB_ErrorType
  TSTTG_gentityXyzToUv (TSTTG_Instance instance,
                        /*in*/ ARRAY_IN_DECL(TSTTG_CGentityHandle, gentity_handles),
                        /*in*/ CONST_ARRAY_IN_DECL(double, coordinates),
                        /*out*/ ARRAY_INOUT_DECL(double, uv));

  enum TSTTB_ErrorType
  TSTTG_gentityUvRange (TSTTG_Instance instance,
                        /*in*/ ARRAY_IN_DECL(TSTTG_CGentityHandle, gentity_handles),
                        /*out*/ ARRAY_INOUT_DECL(double, uv_min),
                        /*out*/ ARRAY_INOUT_DECL(double, uv_max));

  enum TSTTB_ErrorType
  TSTTG_Greparam_edge_face (TSTTG_Instance instance,
                            /*in*/ ARRAY_IN_DECL(TSTTG_GentityHandle, src_gentity_handles),
                            /*in*/ ARRAY_IN_DECL(double, src_uv),
                            /*in*/ ARRAY_IN_DECL(TSTTG_GentityHandle, trg_gentity_handles),
                            /*in*/ ARRAY_IN_DECL(double, trg_uv));

  enum TSTTB_ErrorType
  TSTTG_gentityNormalUv (TSTTG_Instance instance,
                         /*in*/ ARRAY_IN_DECL(TSTTG_CGentityHandle, gface_handles),
                         /*in*/ CONST_ARRAY_IN_DECL(double, parameters),
                         /*out*/ ARRAY_INOUT_DECL(double, normals));

  enum TSTTB_ErrorType
  TSTTG_gentityTangentU (TSTTG_Instance instance,
                         /*in*/ ARRAY_IN_DECL(TSTTG_CGentityHandle, gedge_handles),
                         /*in*/ CONST_ARRAY_IN_DECL(double, parameters),
                         /*out*/ ARRAY_INOUT_DECL(double, tangents));

  enum TSTTB_ErrorType
  TSTTG_setEntSetData (TSTTG_Instance instance,
                       /*in*/ TSTTG_GentitySetHandle entity_set,
                       /*in*/ const TSTTG_TagHandle tag_handle,
                       /*in*/ CONST_ARRAY_IN_DECL(char, tag_value));

  enum TSTTB_ErrorType
  TSTTG_setEntSetIntData (TSTTG_Instance instance,
                          /*in*/ TSTTG_GentitySetHandle entity_set,
                          /*in*/ const TSTTG_TagHandle tag_handle,
                          /*in*/ const int tag_value);

  enum TSTTB_ErrorType
  TSTTG_setEntSetDblData (TSTTG_Instance instance,
                          /*in*/ TSTTG_GentitySetHandle entity_set,
                          /*in*/ const TSTTG_TagHandle tag_handle,
                          /*in*/ const double tag_value);

  enum TSTTB_ErrorType
  TSTTG_setEntSetEHData (TSTTG_Instance instance,
                         /*in*/ TSTTG_GentitySetHandle entity_set,
                         /*in*/ const TSTTG_TagHandle tag_handle,
                         /*in*/ TSTTG_CGentityHandle tag_value);

  enum TSTTB_ErrorType
  TSTTG_getEntSetData (TSTTG_Instance instance,
                       /*in*/ TSTTG_CGentitySetHandle entity_set,
                       /*in*/ const TSTTG_TagHandle tag_handle,
                       /*inout*/ ARRAY_INOUT_DECL(char, tag_value));

  int
  TSTTG_getEntSetIntData (TSTTG_Instance instance,
                          /*in*/ TSTTG_CGentitySetHandle entity_set,
                          /*in*/ const TSTTG_TagHandle tag_handle);

  double
  TSTTG_getEntSetDblData (TSTTG_Instance instance,
                          /*in*/ TSTTG_CGentitySetHandle entity_set,
                          /*in*/ const TSTTG_TagHandle tag_handle);

  TSTTG_GentityHandle
  TSTTG_getEntSetEHData (TSTTG_Instance instance,
                         /*in*/ TSTTG_CGentitySetHandle entity_set,
                         /*in*/ const TSTTG_TagHandle tag_handle);

  enum TSTTB_ErrorType
  TSTTG_getAllEntSetTags (TSTTG_Instance instance,
                          /*in*/ TSTTG_CGentitySetHandle entity_set,
                          /*inout*/ ARRAY_INOUT_DECL(void*, tag_handles));

  enum TSTTB_ErrorType
  TSTTG_rmvEntSetTag (TSTTG_Instance instance,
                      /*in*/ TSTTG_GentitySetHandle entity_set,
                      /*in*/ const TSTTG_TagHandle tag_handle);


  enum TSTTB_ErrorType
  TSTTG_load_cub_geometry(const char *name);

  
#ifdef __cplusplus
} // extern "C"
#endif


#endif // #ifndef __CGM_GEOM_HPP__


  
  
