#ifndef GEOM_CCAPI_HPP
#define GEOM_CCAPI_HPP

/**
 * \file geom_ccapi.hpp
 *
 * \brief C api to classes in CUBIT geom subdirectory
 *
 * This file contains api functions for public members of the following classes:
 *
 * GeometryTool
 * MergeTool
 *
 * Classes should not be added to this list unless they are located in this
 * subdirectory (avoids illegal dependencies between directories)
 *
 * These functions are compiled with the C++ compiler (so that they
 * can call C++ member functions and understand C++ classes), but with
 * linkage for the C language (so they can be called from C functions)
 *
 * NAMING CONVENTIONS
 * In general, given a CUBIT class ClassName, which uses the CUBIT
 * convention for naming classes, the api functions that access member
 * functions will have the same name, but with a CN_ prefix.  For example,
 * ClassName::function_name() will have a corresponding api function
 * CN_function_name().  
 *
 * ARGUMENT AND RETURN TYPE CONVERSION CONVENTIONS
 * Some types of arguments and return types are not usable in C, namely
 * classes and their pointers.  Therefore, the following types of arguments
 * and return types are converted to or from (void *) before being input to or
 * output from these api functions:
 *
 * CubitClass *, CubitClass &, CubitClass *&, DLCubitEntList &, DLCubitEntList *&
 *
 * AUXILLIARY DATASTRUCTURES
 * Certain CUBIT objects are used extensively in the geometry interface; it
 * makes sense to define auxilliary structures which hold the data in these
 * classes, and which can be used directly in the application code.  Auxilliary
 * structures are named by appending 'Struct' to the CUBIT class name,
 * e.g. CubitVectorStruct.  The CUBIT classes for which auxilliary 
 * structures are defined are:
 *
 * CubitVector, CubitPlane, CubitBox
 *
 * In addition, (char *) is used in place of CubitString.
 *
 * The argument and return type conversions between these classes and 
 * their structures are defined as follows, where Xxx denotes the CUBIT
 * class name:
 *
 * Xxx *  => XxxStruct *
 * Xxx &  => XxxStruct &
 * DLXxxList &  => XxxStruct **
 * DLXxxList *& => XxxStruct ***
 *
 * \author Tim Tautges
 *
 * \date 1/2000
 *
 */

#include "CubitDefines.h"
#include "GeometryDefines.h"
#include "EntityType.h"
#include "CubitVectorStruct.h"
#include "CubitBoxStruct.h"
#include "CubitPlaneStruct.h"

#ifdef __cplusplus
extern "C" {
#endif

/* GeometryTool * */ 
extern void *GeometryTool_instance(/* SolidModelingEngine * */ void *SMEPtr);
  
enum CubitStatus GeometryTool_cubit_entity_list(const char *keyword, 
                                   /* DLRefEntityList & */ void ***entity_list, 
                                 int *entity_list_size,
                                                const enum CubitBoolean print_errors);
/* API function for GeometryTool::cubit_entity_list */


enum CubitStatus GeometryTool_ref_entity_list( char const* keyword,
                                  /* DLRefEntityList &*/ void ***entity_list,
                                int *entity_list_size,
                                               const enum CubitBoolean print_errors);

void GeometryTool_bodies (/*DLBodyList &*/ void ***bodies, int *bodies_size);

void GeometryTool_ref_volumes (/*DLRefVolumeList &*/ void ***ref_volumes, int *ref_volumes_size);

void GeometryTool_ref_groups (/*DLRefGroupList &*/ void ***ref_groups, int *ref_groups_size);

void GeometryTool_ref_parts (/*DLRefPartList &*/ void ***ref_parts, int *ref_parts_size);

void GeometryTool_ref_assemblies (/*DLRefAssemblyList &*/ void ***ref_assemblies, int *ref_assemblies_size);

void GeometryTool_ref_faces (/*DLRefFaceList &*/ void ***ref_faces, int *ref_faces_size);

void GeometryTool_ref_edges (/*DLRefEdgeList &*/ void ***ref_edges, int *ref_edges_size);

void GeometryTool_ref_vertices (/*DLRefVertexList &*/ void ***ref_vertices, int *ref_vertices_size);

void GeometryTool_ref_coordsys (/*DLRefCoordSysList &*/ void ***ref_coordsys, int *ref_coordsys_size);


int GeometryTool_num_bodies();
/* API function for GeometryTool::num_bodies */

int GeometryTool_num_ref_volumes();
/* API function for GeometryTool::num_ref_volumes */

int GeometryTool_num_ref_groups();
/* API function for GeometryTool::num_ref_groups */

int GeometryTool_num_ref_parts();
/* API function for GeometryTool::num_ref_parts */

int GeometryTool_num_ref_assemblies();
/* API function for GeometryTool::num_ref_assemblies */

int GeometryTool_num_ref_faces();
/* API function for GeometryTool::num_ref_faces */

int GeometryTool_num_ref_edges();
/* API function for GeometryTool::num_ref_edges */

int GeometryTool_num_ref_vertices();
/* API function for GeometryTool::num_ref_vertices */

int GeometryTool_num_ref_coordsys();
/* API function for GeometryTool::num_ref_coordsys */


/* RefEntity* */ void* GeometryTool_get_ref_entity_1 (const char *type, int id);
/* API function for GeometryTool::*get_ref_entity  */

/* RefEntity* */ void* GeometryTool_get_ref_entity_2 (const enum EntityType type, int id);
/* API function for GeometryTool::*get_ref_entity  */

/* Body* */ void* GeometryTool_get_body ( int id );
/* API function for GeometryTool::*get_body  */

/* RefVolume* */ void* GeometryTool_get_ref_volume ( int id );
/* API function for GeometryTool::*get_ref_volume  */

/* RefGroup* */ void* GeometryTool_get_ref_group ( int id );
/* API function for GeometryTool::*get_ref_group  */

/* RefPart* */ void* GeometryTool_get_ref_part ( int id );
/* API function for GeometryTool::*get_ref_part  */

/* RefAssembly* */ void* GeometryTool_get_ref_assembly ( int id );
/* API function for GeometryTool::*get_ref_assembly  */

/* RefFace* */ void* GeometryTool_get_ref_face ( int id );
/* API function for GeometryTool::*get_ref_face  */

/* RefEdge* */ void* GeometryTool_get_ref_edge ( int id );
/* API function for GeometryTool::*get_ref_edge  */

/* RefVertex* */ void* GeometryTool_get_ref_vertex ( int id );
/* API function for GeometryTool::*get_ref_vertex  */


/* Body* */ void* GeometryTool_get_first_body ();
/* API function for GeometryTool::*get_first_body  */

/* RefVolume* */ void* GeometryTool_get_first_ref_volume ();
/* API function for GeometryTool::*get_first_ref_volume  */

/* RefGroup* */ void* GeometryTool_get_first_ref_group ();
/* API function for GeometryTool::*get_first_ref_group  */

/* RefPart* */ void* GeometryTool_get_first_ref_part ();
/* API function for GeometryTool::*get_first_ref_part  */

/* RefAssembly* */ void* GeometryTool_get_first_ref_assembly ();
/* API function for GeometryTool::*get_first_ref_assembly  */

/* RefFace* */ void* GeometryTool_get_first_ref_face ();
/* API function for GeometryTool::*get_first_ref_face  */

/* RefEdge* */ void* GeometryTool_get_first_ref_edge ();
/* API function for GeometryTool::*get_first_ref_edge  */

/* RefVertex* */ void* GeometryTool_get_first_ref_vertex ();
/* API function for GeometryTool::*get_first_ref_vertex  */
  

/* Body* */ void* GeometryTool_get_next_body ();
/* API function for GeometryTool::*get_next_body  */

/* RefVolume* */ void* GeometryTool_get_next_ref_volume ();
/* API function for GeometryTool::*get_next_ref_volume  */

/* RefGroup* */ void* GeometryTool_get_next_ref_group ();
/* API function for GeometryTool::*get_next_ref_group  */

/* RefPart* */ void* GeometryTool_get_next_ref_part ();
/* API function for GeometryTool::*get_next_ref_part  */

/* RefAssembly* */ void* GeometryTool_get_next_ref_assembly ();
/* API function for GeometryTool::*get_next_ref_assembly  */

/* RefFace* */ void* GeometryTool_get_next_ref_face ();
/* API function for GeometryTool::*get_next_ref_face  */

/* RefEdge* */ void* GeometryTool_get_next_ref_edge ();
/* API function for GeometryTool::*get_next_ref_edge  */

/* RefVertex* */ void* GeometryTool_get_next_ref_vertex ();
/* API function for GeometryTool::*get_next_ref_vertex */
  

/* Body* */ void* GeometryTool_get_last_body ();
/* API function for GeometryTool::*get_last_body  */

/* RefVolume* */ void* GeometryTool_get_last_ref_volume ();
/* API function for GeometryTool::*get_last_ref_volume  */

/* RefGroup* */ void* GeometryTool_get_last_ref_group ();
/* API function for GeometryTool::*get_last_ref_group  */

/* RefPart* */ void* GeometryTool_get_last_ref_part ();
/* API function for GeometryTool::*get_last_ref_part  */

/* RefAssembly* */ void* GeometryTool_get_last_ref_assembly ();
/* API function for GeometryTool::*get_last_ref_assembly  */

/* RefFace* */ void* GeometryTool_get_last_ref_face ();
/* API function for GeometryTool::*get_last_ref_face  */

/* RefEdge* */ void* GeometryTool_get_last_ref_edge ();
/* API function for GeometryTool::*get_last_ref_edge  */

/* RefVertex* */ void* GeometryTool_get_last_ref_vertex ();
/* API function for GeometryTool::*get_last_ref_vertex  */

  /**< 
    <HR><H3>File operations (import, export)<H3>
  */

enum CubitStatus GeometryTool_export_solid_model( /* DLRefEntityList& */ void ***ref_entity_list,
                                   int *ref_entity_list_size,
                                   char const* file_name,
                                   char const* file_type,
                                   const /* CubitString & */ char *cubit_version,
                                   char const* logfile_name);

enum CubitStatus GeometryTool_import_solid_model_1(FILE* file_ptr, 
                                  const char* file_name,
                                  const char* file_type,
                                  const char* logfile_name,
                                                 enum CubitBoolean heal_step,
                                                 enum CubitBoolean import_bodies,
                                                 enum CubitBoolean import_surfaces,
                                                 enum CubitBoolean import_curves,
                                                 enum CubitBoolean import_vertices,
                                                 enum CubitBoolean free_surfaces);
/* API function for GeometryTool::import_solid_model */

enum CubitStatus GeometryTool_import_solid_model_2(enum SolidModelerType model_type);
/* API function for GeometryTool::import_solid_model */

enum CubitStatus GeometryTool_import_datum_curves(FILE* inputFile, 
                                   char const* filetype);
/* API function for GeometryTool::import_datum_curves */
   
enum CubitStatus GeometryTool_read_geometry_file(char const* fileName, 
                                  char const* includePath,
                                  char const* type);
/* API function for GeometryTool::read_geometry_file */
   
int GeometryTool_fire_ray(/* Body* */ void *body, 
                  /* const CubitVector */ struct CubitVectorStruct *ray_point,
                  /* const CubitVector */ struct CubitVectorStruct *unit,
                  /* int & */ int *num_hit, /* double *& */ double **ray_params,
                  /* DLRefEntityList * */ void ***entity_list,
                int *entity_list_size);
/* API function for GeometryTool::fire_ray */

  /**
   <HR><H3>Geometric primitive creation functions<H3>
  */
  int GeometryTool_sphere(double radius);
/* API function for GeometryTool::sphere */
   
  int GeometryTool_brick( double wid, double dep, double hi );
/* API function for GeometryTool::brick     */
   
  int GeometryTool_prism( double height, int sides, double major, double minor);
/* API function for GeometryTool::prism     */
   
int GeometryTool_pyramid( double height, int sides, double major, double minor,
                double top);
/* API function for GeometryTool::pyramid   */
   
int GeometryTool_cylinder( double hi, double r1, double r2, double r3 );
/* API function for GeometryTool::cylinder  */
   
int GeometryTool_torus( double r1, double r2 );
/* API function for GeometryTool::torus     */

int GeometryTool_planar_sheet ( struct CubitPlaneStruct plane, 
                                struct CubitBoxStruct bounding_box,
                      int extension_type, double extension);
/* API function for GeometryTool::planar_sheet  */
   
    /* <HR><H3>Geometry transformation (move, scale, etc.)<H3> */


/* Body* */ void* GeometryTool_copy_body ( /* Body* */ void *body_ptr);
/* API function for GeometryTool::copy_body  */

/* RefEntity* */ void* GeometryTool_copy_refentity( /* RefEntity * */ void *old_entity );
/* API function for GeometryTool::copy_refentity */

enum CubitStatus GeometryTool_reflect_body( /* Body * */ void *body_ptr,
                               /* Body *& */ void **new_body,
                             double x, double y,
                             double z,
                                            enum CubitBoolean keep_old);
/* API function for GeometryTool::reflect_body */

enum CubitStatus GeometryTool_align_body( /* Body * */ void *body_ptr, 
                             /* RefFace * */ void *my_face,
                             /* RefFace * */ void *target_face );
/* API function for GeometryTool::align_body */

   
    /* <HR><H3>Geometry modification (booleans, decomposition, etc.)<H3> */

enum CubitStatus GeometryTool_unite_1(/* DLBodyList & */ void ***bodies, int *bodies_size, 
                       /* Body *& */ void **newBody,
                     int keep_old);
/* API function for GeometryTool::unite */

enum CubitStatus GeometryTool_unite_2 ( /* Body * */ void *body1, /* Body * */ void *body2, 
                         /* Body *& */ void **newBody,
                       int keep_old);
/* API function for GeometryTool::unite  */

enum CubitStatus GeometryTool_subtract ( /* Body * */ void *tool_body_ptr, 
                            /* DLBodyList & */ void ***from_bodies,
                          int *from_bodies_size,
                            /* DLBodyList & */ void ***new_bodies, int *new_bodies_size,
                          int keep_old);
/* API function for GeometryTool::subtract  */

enum CubitStatus GeometryTool_intersect ( /* Body * */ void *tool_body_ptr, 
                             /* DLBodyList & */ void ***from_bodies,
                           int *from_bodies_size,
                             /* DLBodyList & */ void ***new_bodies,
                           int *new_bodies_size,
                           int keep_old);
/* API function for GeometryTool::intersect  */

enum CubitStatus GeometryTool_section( /* DLBodyList & */ void ***section_body_list,
                        int *section_body_list_size,
                          /* CubitVector & */ struct CubitVectorStruct *point_1,
                          /* CubitVector & */ struct CubitVectorStruct *point_2,
                          /* CubitVector & */ struct CubitVectorStruct *point_3,
                          /* DLBodyList & */ void ***new_body_list,
                        int *new_body_list_size,
                        enum CubitBoolean keep_normal_side,
                        enum CubitBoolean keep_old);
/* API function for GeometryTool::section */

enum CubitStatus GeometryTool_split_periodic( /* Body * */ void *body_ptr,
                                 /* Body *& */ void **new_body_ptr );
/* API function for GeometryTool::split_periodic */

  int GeometryTool_webcut_with_cylinder( /* DLBodyList& */ void ***webcut_body_list, int *webcut_body_list_size,
                                         double radius, /* CubitVector & */ struct CubitVectorStruct *axis,
                                           /* CubitVector& */ struct CubitVectorStruct *center, 
                                           /* DLBodyList & */ void ***results_list, int *results_list_size,
                                         enum CubitBoolean imprint,
                                         enum CubitBoolean merge);
  /**<  Webcuts the bodies in the list with a cutting cylinder.
    *  The cylinder is created by the given parameters.  This
    *  is done in the solid modeling engine to reduce the impact
    *  on body ids.
    */

  int  GeometryTool_webcut_across_translate( /* DLBodyList& */ void ***body_list, int *body_list_size,
                                               /* RefFace* */ void *plane_surf1, /* RefFace* */ void *plane_surf2,
                                               /* DLBodyList& */ void ***results_list, int *results_list_size,
                                             enum CubitBoolean imprint, enum CubitBoolean merge);
  /**<  Webcuts with a flat plate to make a body suitable for single-single
    *   sweeping.  Only experimental and available in the GUI.
    */
   
int GeometryTool_webcut_with_plane(/* DLBodyList& */ void ***webcut_body_list,
                         int *webcut_body_list_size,
                           /* CubitVector & */ struct CubitVectorStruct *vector1,
                           /* CubitVector & */ struct CubitVectorStruct *vector2,
                           /* CubitVector & */ struct CubitVectorStruct *vector3, 
                           /* DLBodyList & */ void ***results_list,
                         int *results_list_size,
                         enum CubitBoolean imprint,
                         enum CubitBoolean merge);
/* API function for GeometryTool::webcut_with_plane */

int GeometryTool_webcut_with_vertices( /* DLBodyList& */ void ***webcut_body_list,
                             int *webcut_body_list_size,
                               /* RefVertex* */ void *refVertex1,
                               /* RefVertex* */ void *refVertex2,
                               /* RefVertex* */ void *refVertex3, 
                               /* DLBodyList & */ void ***results_list,
                             int *results_list_size,
                                       enum CubitBoolean imprint,
                                       enum CubitBoolean merge);
/* < API function for GeometryTool::webcut_with_vertices */
   
  int GeometryTool_webcut_with_curve_loop(/* DLBodyList& */ void ***webcut_body_list, int *webcut_body_list_size,
                                            /* DLRefEdgeList& */ void ***refedge_list, int *refedge_list_size,
                                            /* DLBodyList& */ void ***results_list, int *results_list_size,
                                          enum CubitBoolean imprint);
   
int GeometryTool_webcut_with_surface(/* DLBodyList& */ void ***webcut_body_list,
                           int *webcut_body_list_size,
                             /* RefFace* */ void *refFace,
                             /* DLBodyList& */ void *** results_list,
                           int * results_list_size,
                           enum CubitBoolean imprint,
                                     enum CubitBoolean merge);
  
/* < API function for GeometryTool::webcut_with_surface */
   
enum CubitStatus GeometryTool_webcut_with_sheet(/* Body * */ void *webcut_body,
                                   /* Body * */ void *sheet_body,
                              /* DLBodyList & */ void ***new_bodies,
                                 int *new_bodies_size,
                                           enum CubitBoolean imprint );
  
/* < API function for GeometryTool::webcut_with_sheet */
   
int GeometryTool_webcut_with_body(/* DLBodyList& */ void ***webcut_body_list,
                        int *webcut_body_list_size,
                          /* Body* */ void *body,
                          /* DLBodyList& */ void ***results_list,
                        int * results_list_size,
                        enum CubitBoolean imprint,
                                  enum CubitBoolean merge) ;

/* < API function for GeometryTool::webcut_with_body */
   
enum CubitStatus GeometryTool_webcut_with_extended_surf( /* DLBodyList & */ void ***webcut_body_list,
                                          int *webcut_body_list_size,
                                            /* RefFace * */ void *face_to_extend,
                                            /* DLBodyList & */ void ***new_bodies,
                                          int *new_bodies_size,
                                            /* int & */ int *num_cut,
                                                    enum CubitBoolean imprint );

/* < API function for GeometryTool::webcut_with_extended_surf */

enum CubitStatus GeometryTool_remove_surfaces( /* DLRefFaceList & */ void ***ref_face_list,
                                int *ref_face_list_size,
                                  /* DLBodyList & */ void ***new_body_list,
                                int *new_body_list_size,
                                enum CubitBoolean extend_adjoining,
                                enum CubitBoolean keep_surface,
                                enum CubitBoolean keep_old_body );
/* API function for GeometryTool::remove_surfaces */

enum CubitStatus GeometryTool_offset_surfaces( /* DLRefFaceList& */ void ***ref_face_list,
                                int * ref_face_list_size,
                                  /* DLBodyList& */ void ***new_body_list,
                                int * new_body_list_size, double offset_distance,
                                enum CubitBoolean keep_old_body );
/* API function for GeometryTool::offset_surfaces */

enum CubitStatus GeometryTool_offset_curves( /* DLRefEdgeList& */ void ***ref_edge_list,
                              int *ref_edge_list_size, double offset_distance, 
                              /* CubitVector& */ struct CubitVectorStruct* offset_direction, 
                              int gap_type);
/* API function for GeometryTool::offset_curves */

enum CubitStatus GeometryTool_translate_surfaces( /* DLRefFaceList& */ void ***ref_face_list,
                                   int *ref_face_list_size,
                                     /* DLBodyList& */ void ***new_body_list,
                                   int *new_body_list_size,
                                     /* CubitVector & */ struct CubitVectorStruct *delta,
                                   enum CubitBoolean keep_old_body );
/* API function for GeometryTool::translate_surfaces */

enum CubitStatus GeometryTool_replace_surfaces( /* DLRefFaceList & */ void ***ref_face_list,
                                 int *ref_face_list_size,
                                   /* RefFace* */ void *tool_face_ptr,
                              /* DLBodyList & */ void ***new_body_list,
                                 int *new_body_list_size,
                                 enum CubitBoolean reverse_flg,
                                           enum CubitBoolean keep_old_body );
  
/* < API function for GeometryTool::replace_surfaces */
   
enum CubitStatus GeometryTool_trim_curve( /* RefEdge* */ void *trim_curve, 
                             /* CubitVector& */ struct CubitVectorStruct* trim_vector, 
                                       /* CubitVector& */ struct CubitVectorStruct* keep_vector );
  
/* < API function for GeometryTool::trim_curve */
   
  /* < <HR><H3>Topology modification (imprint, regularize, etc.)<H3> */

enum CubitStatus GeometryTool_imprint_1 (/* DLBodyList & */ void ***from_body_list,
                        int *from_body_list_size,
                          /* DLBodyList & */ void ***new_body_list,
                        int *new_body_list_size,
                          int keep_old);
  

enum CubitStatus GeometryTool_imprint_2( /* DLBodyList & */ void ***body_list,
                        int *body_list_size,
                          /* DLRefEdgeList & */ void ***ref_edge_list,
                        int *ref_edge_list_size,
                          /* DLBodyList& */ void ***new_body_list,
                        int *new_body_list_size,
                        int keep_old_body );
/* API function for GeometryTool::imprint */

enum CubitStatus GeometryTool_imprint_3( /* DLRefFaceList & */ void ***ref_face_list,
                        int *ref_face_list_size,
                          /* DLRefEdgeList & */ void ***ref_edge_list,
                        int *ref_edge_list_size,
                          /* DLBodyList& */ void ***new_body_list,
                        int *new_body_list_size,
                        int keep_old_body );
/* API function for GeometryTool::imprint */

enum CubitStatus GeometryTool_imprint_4( /* DLBodyList & */ void ***body_list,
                        int *body_list_size,
                          /* DLRefVertexList & */ void ***ref_vertex_list,
                        int *ref_vertex_list_size,
                          /* DLBodyList& */ void ***new_body_list,
                        int *new_body_list_size,
                        int keep_old_body );
/* API function for GeometryTool::imprint */

enum CubitStatus GeometryTool_regularize_body( /* Body * */ void *body_ptr, 
                                  /* Body *& */ void **new_body );
/* API function for GeometryTool::regularize_body */
   
enum CubitStatus GeometryTool_split_body( /* Body * */ void *body_ptr,
                             /* DLBodyList & */ void ***new_bodies,
                           int *new_bodies_size);
/* API function for GeometryTool::split_body */
   
    /* <HR><H3>GeometryTool options and settings<H3> */

void GeometryTool_group_imprint_1(enum CubitBoolean flag);
/* API function for GeometryTool::group_imprint */

enum CubitBoolean GeometryTool_group_imprint_2();
/* API function for GeometryTool::group_imprint */
   
void GeometryTool_set_all_edges_imprint( enum CubitBoolean flag);
/* API function for GeometryTool::void set_all_edges_imprint */

enum CubitBoolean GeometryTool_get_all_edges_imprint();
/* API function for GeometryTool::CubitBoolean get_all_edges_imprint */
   
    /* static */ void GeometryTool_booleans_after_merge_1( enum CubitBoolean flag );
    /* static */ enum CubitBoolean GeometryTool_booleans_after_merge_2();
  /**<  sets/gets the booleansAfterMerge flag.
    */
   
void GeometryTool_new_ids_1(enum CubitBoolean flag);

/* API function for GeometryTool::new_ids */

enum CubitBoolean GeometryTool_new_ids_2();
/* API function for GeometryTool::new_ids */
   
  void GeometryTool_geom_debug( /* DLTopologyEntityList */ void ***arg, int *arg_size);


/* API function for GeometryTool::void geom_debug */
   
  void GeometryTool_use_facet_bbox_1( enum CubitBoolean pass_flag );
  enum CubitBoolean GeometryTool_use_facet_bbox_2();
  /**< get and set the useFacetBBox flag. (compute our own bbox
    * instead of relying on bad ACIS geometry)
    */

   
#ifdef CUBIT_GUI
void GeometryTool_set_validate_file_ptr( FILE* file_ptr);
/* API function for GeometryTool::set_validate_file_ptr */

void GeometryTool_reset_validate_file_ptr();
/* API function for GeometryTool::reset_validate_file_ptr */

FILE* GeometryTool_get_validate_file_ptr();
/* API function for GeometryTool::get_validate_file_ptr */
#endif
   
    /* <HR><H3>SolidModelingEngine information<H3> */

char *GeometryTool_identify_modeling_engine();
/* API function for GeometryTool::identify_modeling_engine */
   
void GeometryTool_register_solid_modeling_engine(/* SolidModelingEngine* */ void *SMEPtr);
/* API function for GeometryTool::register_solid_modeling_engine */
     
void GeometryTool_register_datum_curve_sm_engine(/* SolidModelingEngine* */ void *SMEPtr);
/* API function for GeometryTool::register_datum_curve_sm_engine */

enum CubitStatus GeometryTool_set_engine_version(int version);
/* API function for GeometryTool::set_engine_version */
   
enum CubitStatus GeometryTool_get_engine_version(int *version);
/* API function for GeometryTool::get_engine_version */
   
enum CubitStatus GeometryTool_list_engine_versions(char **versions);
/* API function for GeometryTool::list_engine_versions */

double GeometryTool_get_sme_resabs_tolerance();
/* API function for GeometryTool::get_sme_resabs_tolerance */

double GeometryTool_set_sme_resabs_tolerance( double new_resabs);
/* API function for GeometryTool::set_sme_resabs_tolerance */

enum CubitStatus GeometryTool_set_steptools_path( const char* path);
/* API function for GeometryTool::set_steptools_path */

enum CubitStatus GeometryTool_set_igestools_path( const char* path);
/* API function for GeometryTool::set_igestools_path */

enum CubitStatus GeometryTool_set_sme_int_option( const char* opt_name, int val);
/* API function for GeometryTool::set_sme_int_option */

enum CubitStatus GeometryTool_set_sme_dbl_option( const char* opt_name, double val);
/* API function for GeometryTool::set_sme_dbl_option */

enum CubitStatus GeometryTool_set_sme_str_option( const char* opt_name, const char* val);
/* API function for GeometryTool::set_sme_str_option */

    /* <HR><H3>Topology/geometry creation functions<H3> */

/* RefVertex* */ void* GeometryTool_make_RefVertex(enum GeometryType ref_vertex_type,
                                         const struct CubitVectorStruct *point, int color);
/* API function for GeometryTool::make_RefVertex */
   
/* Body* */ void* GeometryTool_make_Body_1(/* Surface * */ void *surface);
/* API function for GeometryTool::*make_Body */

/* Body* */ void* GeometryTool_make_Body_2(/* BodySM * */ void *bodysm_ptr);
/* API function for GeometryTool::make_Body */

/* Shell* */ void* GeometryTool_make_Shell(/* ShellSM * */ void *shellsm_ptr );
/* API function for GeometryTool::make_Shell */

/* Loop* */ void* GeometryTool_make_Loop(/* LoopSM* */ void *loopsm_ptr );
/* API function for GeometryTool::make_Loop */

    /* RefEntity * */ void *GeometryTool_check_mergeable_refentity(/* TopologyBridge * */ void *bridge);
    /* check for mergeable ref entity, indicated by a merge attribute on the */
    /* bridge */

/* Chain* */ void* GeometryTool_make_Chain_1(/* DLPointList& */ void ***points, 
                                 int *points_size);
/* API function for GeometryTool::make_Chain */

/* Chain* */ void* GeometryTool_make_Chain_2(/* Curve * */ void *curve);
/* API function for GeometryTool::make_Chain */

/* RefEdge* */ void* GeometryTool_make_RefEdge_1( enum GeometryType ref_edge_type,
                                        /* RefVertex const* */ void *ref_vertex_1,
                                        /* RefVertex const* */ void *ref_vertex_2,
                                        /* DLCubitVectorList& */ struct CubitVectorStruct *vector_list,
                                      int vector_list_size,
                                        /* RefFace* */ void *reffaca_ptr );
/* API function for GeometryTool::make_RefEdge */

/* RefEdge* */ void* GeometryTool_make_RefEdge_2(  /* RefVertex const* */ void *ref_vertex_1,
                                         /* RefVertex const* */ void *ref_vertex_2,
                                         /* RefFace* */ void *ref_face_ptr,
                                         /* RefVertex const* */ void *ref_vertex_3 );
/* API function for GeometryTool::make_RefEdge */

    /* RefEdge* */ void *GeometryTool_make_RefEdge_3( /* RefEdge * */ void *ref_edge);

/* RefEdge* */ void* GeometryTool_make_RefEdge_4(enum GeometryType ref_edge_type,
                                       /* RefVertex const* */ void *ref_vertex_1,
                                       /* RefVertex const* */ void *ref_vertex_2,
                                     struct CubitVectorStruct const* intermediate_point,
                                                 enum CubitSense sense);
/* API function for GeometryTool::make_RefEdge */

/* RefFace* */ void* GeometryTool_make_RefFace_1(/* RefFace * */ void *from_ref_face,
                                     enum CubitBoolean extended_from);
/* API function for GeometryTool::make_RefFace */

/* RefFace* */ void* GeometryTool_make_RefFace_2(enum GeometryType ref_face_type,
                                       /* DLRefEdgeList& */ void ***ref_edge_list,
                                     int *ref_edge_list_size,
                                       /* RefFace * */ void *ref_face_ptr);
/* API function for GeometryTool::make_RefFace */

/* RefVolume* */ void* GeometryTool_make_RefVolume(enum GeometryType ref_volume_type,
                                           /* DLRefFaceList& */ void ***ref_face_list,
                                         int *ref_face_list_size);
/* API function for GeometryTool::make_RefVolume */

/* Body* */ void* GeometryTool_make_Body_3(/* DLRefVolumeList& */ void ***ref_volume_list, 
                               int *ref_volume_list_size);
/* API function for GeometryTool::make_Body */

/* Body* */ void* GeometryTool_make_Body_4(/* RefFace * */ void *from_ref_face,
                               enum CubitBoolean extended_from);
/* API function for GeometryTool::make_Body */

/* Body* */ void* GeometryTool_make_Body_5(enum GeometryType ref_face_type,
                                 /* DLRefEdgeList& */ void ***ref_edge_list,
                               int *ref_edge_list_size,
                                 /* RefFace * */ void *ref_face_ptr);
/* API function for GeometryTool::make_Body */

/* RefFace* */ void *GeometryTool_make_free_RefFace(/* SurfaceSM * */ void *surfacesm_ptr);

/* RefEdge* */ void *GeometryTool_make_free_RefEdge(/* CurveSM * */ void *curvesm_ptr);

/* RefVertex* */ void *GeometryTool_make_free_RefVertex(/* PointSM * */ void *pointsm_ptr);

enum CubitStatus GeometryTool_sweep_translational(/* DLRefFaceList& */ void ***faces_to_be_swept_list,
                                   int *faces_to_be_swept_list_size,
                                   struct CubitVectorStruct sweep_vector,
                                   double draft_angle,
                                   int draft_type,
                                   int switchside );
/* API function for GeometryTool::sweep_translational */

enum CubitStatus GeometryTool_sweep_rotational(/* DLRefFaceList& */ void ***faces_to_be_swept_list,
                                int *faces_to_be_swept_list_size,
                                struct CubitVectorStruct point,
                                struct CubitVectorStruct direction,
                                double angle,
                                int steps,
                                double draft_angle,
                                int draft_type,
                                int switchside );
/* API function for GeometryTool::sweep_rotational */

enum CubitStatus GeometryTool_sweep_along_curve(/* DLRefFaceList& */ void ***reffaces_to_be_swept_list,
                                 int *reffaces_to_be_swept_list_size,
                              /* DLRefEdgeList& */ void ***ref_edge_list,
                                 int *ref_edge_list_size,
                                 double draft_angle,
                                 int draft_type);
/* API function for GeometryTool::sweep_along_curve */

   
enum CubitStatus GeometryTool_create_body_from_surfs( /* DLRefFaceList & */ void ***ref_face_list,
                                       int *ref_face_list_size,
                                         /* Body *& */ void **new_body,
                                       int keep_old,
                                       int heal);
/* API function for GeometryTool::create_body_from_surfs */

enum CubitStatus GeometryTool_create_net_surface_1( /* DLRefFaceList& */ void ***ref_face_list,
                                   int *ref_face_list_size,
                                     /* Body *& */ void **new_body, 
                                     /* DLIList<DLCubitVectorList*> */ struct CubitVectorStruct *vec_lists_u, 
                                   int vec_lists_u_size,
                                     /* DLIList<DLCubitVectorList*> */ struct CubitVectorStruct *vec_lists_v, 
                                   int vec_lists_v_size,
                                   double net_tol,
                                   enum CubitBoolean heal );
/* API function for GeometryTool::create_net_surface */

enum CubitStatus GeometryTool_create_net_surface_2( /* DLRefEdgeList& */ void ***u_curves,
                                   int *u_curves_size,
                                     /* DLRefEdgeList& */ void ***v_curves,
                                   int *v_curves_size,
                                     /* Body *& */ void ** new_body,
                                   double net_tol, 
                                   enum CubitBoolean heal );
/* API function for GeometryTool::create_net_surface */

enum CubitStatus GeometryTool_create_offset_surface( /* RefFace* */ void *ref_face_ptr, 
                                        /* Body*& */ void ** new_body, 
                                      double offset_distance );
/* API function for GeometryTool::create_offset_surface */

enum CubitStatus GeometryTool_create_offset_body( /* Body * */ void *body_ptr, 
                                     /* Body*& */ void ** new_body, 
                                   double offset_distance );
/* API function for GeometryTool::create_offset_body */

enum CubitStatus GeometryTool_create_skin_surface( /* DLRefEdgeList& */ void ***curves,
                                    int *curves_size,
                                      /* Body*& */ void ** new_body );
/* API function for GeometryTool::create_skin_surface */

  enum CubitStatus GeometryTool_loft_surfaces( /* RefFace * */ void *face1, /* const double & */ double *takeoff1,
                                                 /* RefFace * */ void *face2, /* const double & */ double *takeoff2,
                                                 /* Body*& */ void **new_body,
                                               enum CubitBoolean arc_length_option,
                                               enum CubitBoolean twist_option,
                                               enum CubitBoolean align_direction,
                                               enum CubitBoolean perpendicular,
                                               enum CubitBoolean simplify_option);
    /* loft face1 to face2 */

enum CubitStatus GeometryTool_create_arc_three_1( /* RefVertex * */ void *ref_vertex1, 
                                   /* RefVertex * */ void *ref_vertex2,
                                   /* RefVertex * */ void *ref_vertex3, 
                                 enum CubitBoolean full );
/* API function for GeometryTool::create_arc_three */

enum CubitStatus GeometryTool_create_arc_three_2( /* RefEdge * */ void *ref_edge1, 
                                   /* RefEdge * */ void *ref_edge2,
                              /* RefEdge * */ void *ref_edge3, 
                                 enum CubitBoolean full );
/* API function for GeometryTool::create_arc_three */

enum CubitStatus GeometryTool_create_arc_center_edge( /* RefVertex* */ void *ref_vertex1, 
                                         /* RefVertex* */ void *ref_vertex2,
                                         /* RefVertex* */ void *ref_vertex3, 
                                       double radius,
                                       enum CubitBoolean full );
/* API function for GeometryTool::create_arc_center_edge */


    /* <HR><H3>Topology and geometry deletion<H3> */
void GeometryTool_delete_Body_1(/* DLBodyList& */ void ***body_list,
                    int *body_list_size,
                    enum CubitBoolean remove_solid_model_entities);
/* API function for GeometryTool::delete_Body */

enum CubitStatus GeometryTool_delete_Body_2(
    /* Body * */ void *body_ptr,
    enum CubitBoolean remove_solid_model_entities);
/* API function for GeometryTool::delete_Body */

enum CubitStatus GeometryTool_delete_RefAssembly(
    /* RefAssembly *& */ void *ref_assembly_ptr,
    enum CubitBoolean remove_solid_model_entities);
/* API function for GeometryTool::delete_RefAssembly */

enum CubitStatus GeometryTool_delete_RefPart(
    /* RefPart *& */ void *ref_part_ptr,
    enum CubitBoolean remove_solid_model_entities);
/* API function for GeometryTool::delete_RefPart */

enum CubitStatus GeometryTool_delete_RefEntity(
    /* RefEntity *& */ void *ref_entity_ptr,
    enum CubitBoolean remove_solid_model_entities,
    enum CubitBoolean remove_lower_entities );
/* API function for GeometryTool::delete_RefEntity */

void GeometryTool_cleanout_deactivated_geometry();
/* API function for GeometryTool::cleanout_deactivated_geometry */

void GeometryTool_cleanout_temporary_geometry ();
/* API function for GeometryTool::cleanout_temporary_geometry  */

void GeometryTool_delete_geometry();
/* API function for GeometryTool::delete_geometry */

    /* <HR><H3>Miscellaneous geometry evaluation functions<H3> */
enum CubitStatus GeometryTool_interpolate_along_surface( struct CubitVectorStruct *vector_1,
                                          struct CubitVectorStruct *vector_2,
                                            /* DLCubitVectorList & */ struct CubitVectorStruct *vector_list,
                                          int vector_list_size,
                                            /* RefFace * */ void *ref_face_ptr,
                                          int number_points );
/* API function for GeometryTool::interpolate_along_surface */

enum CubitBoolean GeometryTool_about_spatially_equal_1 (const /* CubitVector& */ struct CubitVectorStruct* Vec1, 
                                       const /* CubitVector& */ struct CubitVectorStruct* Vec2,
                                       double tolerance_factor);
/* API function for GeometryTool::about_spatially_equal  */

enum CubitBoolean GeometryTool_about_spatially_equal_2 (/* RefVertex * */ void *refVertex1, 
                                         /* RefVertex * */ void *refVertex2, 
                                       double tolerance_factor);
/* API function for GeometryTool::about_spatially_equal  */

double GeometryTool_geometric_angle_1(/* RefEdge * */ void *ref_edge_1, 
                                      /* RefEdge * */ void *ref_edge_2,
                            /* RefFace * */ void *ref_face);
/* API function for GeometryTool::geometric_angle */

  double GeometryTool_geometric_angle(/* CoEdge* */ void *co_edge_1, 
                                        /* CoEdge* */ void *co_edge_2 );
/* API function for GeometryTool::geometric_angle */

  double GeometryTool_surface_angle( /* RefEdge * */ void *ref_edge,
                                       /* RefFace * */ void *ref_face_1, 
                                       /* RefFace * */ void *ref_face_2,
                                       /* RefVolume * */ void *ref_volume );
    /**< Calculate dihedral angle at ref_edge between two faces of the
     *   volume. 
     */

enum CubitStatus GeometryTool_get_intersections( /* RefEdge * */ void *ref_edge1, 
                                    /* RefEdge * */ void *ref_edge2,
                                    /* DLCubitVectorList& */ struct CubitVectorStruct **intersection_list,
                                  int *intersection_list_size,
                                  enum CubitBoolean bounded,
                                  enum CubitBoolean closest );
/* API function for GeometryTool::get_intersections */

struct CubitBoxStruct GeometryTool_bounding_box_of_bodies();
/* API function for GeometryTool::bounding_box_of_bodies */


    /* <HR><H3>Healing functions<H3> */
enum CubitStatus GeometryTool_autoheal_bodies(/* DLBodyList & */ void ***body_list,
                               int *body_list_size,
                                 /* DLBodyList & */ void ***new_body_list,
                               int *new_body_list_size,
                               enum CubitBoolean rebuild, 
                               enum CubitBoolean keep_old, enum CubitBoolean make_tolerant,
                               FILE* logfile_ptr );
/* API function for GeometryTool::autoheal_bodies */

enum CubitStatus GeometryTool_healer_analyze_badgeom( /* DLBodyList & */ void ***body_list,
                                       int *body_list_size,
                                       FILE* logfile );
/* API function for GeometryTool::healer_analyze_badgeom */

enum CubitStatus GeometryTool_healer_show_badgeom( /* DLBodyList & */ void ***body_list,
                                    int *body_list_size);
/* API function for GeometryTool::healer_show_badgeom */

enum CubitStatus GeometryTool_healer_show_tcurves( /* DLBodyList & */ void ***body_list, int *body_list_size);
/* API function for GeometryTool::healer_show_tcurves */

enum CubitStatus GeometryTool_heal_incremental( /* DLBodyList & */ void ***body_list,
                                 int *body_list_size,
                              /* DLBodyList & */ void ***new_bodies,
                                 int *new_bodies_size,
                                 double simplify_tol,
                                 double stitch_min_tol,
                                 double stitch_max_tol,
                                 double geombuild_tol,
                                 double analytic_tol,
                                 double isospline_tol,
                                 double reblend_classify_tol,
                                 double reblend_tol,
                                 enum CubitBoolean keep_old,
                                 enum CubitBoolean make_tolerant,
                                 FILE* logfile_ptr);
/* API function for GeometryTool::heal_incremental */

void GeometryTool_healer_list_incremental();
/* API function for GeometryTool::healer_list_incremental */

void GeometryTool_healer_list_tolerances( /* DLBodyList & */ void ***body_list,
                                int *body_list_size);
/* API function for GeometryTool::healer_list_tolerances */

double GeometryTool_healer_get_default_simplify_tol();
/* API function for GeometryTool::healer_get_default_simplify_tol */

void GeometryTool_healer_set_default_simplify_tol( double tol);
/* API function for GeometryTool::healer_set_default_simplify_tol */

double GeometryTool_healer_get_default_stitch_min_tol();
/* API function for GeometryTool::healer_get_default_stitch_min_tol */

void GeometryTool_healer_set_default_stitch_min_tol( double tol);
/* API function for GeometryTool::healer_set_default_stitch_min_tol */

double GeometryTool_healer_get_default_stitch_max_tol();
/* API function for GeometryTool::healer_get_default_stitch_max_tol */

void GeometryTool_healer_set_default_stitch_max_tol( double tol);
/* API function for GeometryTool::healer_set_default_stitch_max_tol */

double GeometryTool_healer_get_default_geombuild_tol();
/* API function for GeometryTool::healer_get_default_geombuild_tol */

void GeometryTool_healer_set_default_geombuild_tol( double tol);
/* API function for GeometryTool::healer_set_default_geombuild_tol */

double GeometryTool_healer_get_default_analytic_tol();
/* API function for GeometryTool::healer_get_default_analytic_tol */

void GeometryTool_healer_set_default_analytic_tol( double tol);
/* API function for GeometryTool::healer_set_default_analytic_tol */

double GeometryTool_healer_get_default_isospline_tol();
/* API function for GeometryTool::healer_get_default_isospline_tol */

void GeometryTool_healer_set_default_isospline_tol( double tol);
/* API function for GeometryTool::healer_set_default_isospline_tol */

double GeometryTool_healer_get_default_reblend_classify_tol();
/* API function for GeometryTool::healer_get_default_reblend_classify_tol */

void GeometryTool_healer_set_default_reblend_classify_tol( double tol);
/* API function for GeometryTool::healer_set_default_reblend_classify_tol */

double GeometryTool_healer_get_default_reblend_tol();
/* API function for GeometryTool::healer_get_default_reblend_tol */

void GeometryTool_healer_set_default_reblend_tol( double tol);
/* API function for GeometryTool::healer_set_default_reblend_tol */

void GeometryTool_healer_reset_default_tolerances();
/* API function for GeometryTool::healer_reset_default_tolerances */

void GeometryTool_healer_list_default_tolerances();
/* API function for GeometryTool::healer_list_default_tolerances */

void GeometryTool_healer_clean_attributes( /* DLBodyList& */ void *body_list);
/* API function for GeometryTool::healer_clean_attributes */

enum CubitBoolean GeometryTool_healer_get_cleanatt_flg();
/* API function for GeometryTool::healer_get_cleanatt_flg */

void GeometryTool_healer_set_cleanatt_flg( enum CubitBoolean flg);
/* API function for GeometryTool::healer_set_cleanatt_flg */

int GeometryTool_healer_get_show_method();
/* API function for GeometryTool::healer_get_show_method */

void GeometryTool_healer_set_show_method( int method);
/* API function for GeometryTool::healer_set_show_method */

void GeometryTool_healer_set_show_summary_flg( enum CubitBoolean flg);
/* API function for GeometryTool::healer_set_show_summary_flg */

void GeometryTool_healer_set_show_details_flg( enum CubitBoolean flg);
/* API function for GeometryTool::healer_set_show_details_flg */

void GeometryTool_healer_set_show_on_heal_flg( enum CubitBoolean flg);
/* API function for GeometryTool::healer_set_show_on_heal_flg */

enum CubitBoolean GeometryTool_healer_get_check_vol_on_heal_flg();
/* API function for GeometryTool::healer_get_check_vol_on_heal_flg */

void GeometryTool_healer_set_check_vol_on_heal_flg( enum CubitBoolean flg);
/* API function for GeometryTool::healer_set_check_vol_on_heal_flg */

double GeometryTool_healer_get_vol_on_heal_limit();
/* API function for GeometryTool::healer_get_vol_on_heal_limit */

void GeometryTool_healer_set_vol_on_heal_limit( double limit);
/* API function for GeometryTool::healer_set_vol_on_heal_limit */

void GeometryTool_healer_set_check_surf_on_heal_flg( enum CubitBoolean flg);
/* API function for GeometryTool::healer_set_check_surf_on_heal_flg */

double GeometryTool_healer_get_surf_on_heal_limit();
/* API function for GeometryTool::healer_get_surf_on_heal_limit */

void GeometryTool_healer_set_surf_on_heal_limit( double limit);
/* API function for GeometryTool::healer_set_surf_on_heal_limit */

void GeometryTool_healer_set_check_curve_on_heal_flg( enum CubitBoolean flg);
/* API function for GeometryTool::healer_set_check_curve_on_heal_flg */

double GeometryTool_healer_get_curve_on_heal_limit();
/* API function for GeometryTool::healer_get_curve_on_heal_limit */

void GeometryTool_healer_set_curve_on_heal_limit( double limit);
/* API function for GeometryTool::healer_set_curve_on_heal_limit */

void GeometryTool_healer_set_show_bad_vertices_flg( enum CubitBoolean flg);
/* API function for GeometryTool::healer_set_show_bad_vertices_flg */

enum CubitBoolean GeometryTool_healer_get_show_bad_curves_flg();
/* API function for GeometryTool::healer_get_show_bad_curves_flg */

void GeometryTool_healer_set_show_bad_curves_flg( enum CubitBoolean flg);
/* API function for GeometryTool::healer_set_show_bad_curves_flg */

enum CubitBoolean GeometryTool_healer_get_show_bad_coedges_flg();
/* API function for GeometryTool::healer_get_show_bad_coedges_flg */

void GeometryTool_healer_set_show_bad_coedges_flg( enum CubitBoolean flg);
/* API function for GeometryTool::healer_set_show_bad_coedges_flg */

enum CubitBoolean GeometryTool_healer_get_show_bad_loops_flg();
/* API function for GeometryTool::healer_get_show_bad_loops_flg */

void GeometryTool_healer_set_show_bad_loops_flg( enum CubitBoolean flg);
/* API function for GeometryTool::healer_set_show_bad_loops_flg */

enum CubitBoolean GeometryTool_healer_get_show_bad_surfaces_flg();
/* API function for GeometryTool::healer_get_show_bad_surfaces_flg */

void GeometryTool_healer_set_show_bad_surfaces_flg( enum CubitBoolean flg);
/* API function for GeometryTool::healer_set_show_bad_surfaces_flg */

enum CubitBoolean GeometryTool_healer_get_show_bad_shells_flg();
/* API function for GeometryTool::healer_get_show_bad_shells_flg */

void GeometryTool_healer_set_show_bad_shells_flg( enum CubitBoolean flg);
/* API function for GeometryTool::healer_set_show_bad_shells_flg */

enum CubitBoolean GeometryTool_healer_get_show_bad_volumes_flg();
/* API function for GeometryTool::healer_get_show_bad_volumes_flg */

void GeometryTool_healer_set_show_bad_volumes_flg( enum CubitBoolean flg);
/* API function for GeometryTool::healer_set_show_bad_volumes_flg */

enum CubitBoolean GeometryTool_healer_get_show_bad_bodies_flg();
/* API function for GeometryTool::healer_get_show_bad_bodies_flg */

void GeometryTool_healer_set_show_bad_bodies_flg( enum CubitBoolean flg);
/* API function for GeometryTool::healer_set_show_bad_bodies_flg */

void GeometryTool_healer_list_onshow_flgs();
/* API function for GeometryTool::healer_list_onshow_flgs */

enum CubitBoolean GeometryTool_healer_get_inc_preprocess_flg();
/* API function for GeometryTool::healer_get_inc_preprocess_flg */

void GeometryTool_healer_set_inc_preprocess_flg( enum CubitBoolean flg);
/* API function for GeometryTool::healer_set_inc_preprocess_flg */

enum CubitBoolean GeometryTool_healer_get_inc_simplify_flg();
/* API function for GeometryTool::healer_get_inc_simplify_flg */

void GeometryTool_healer_set_inc_simplify_flg( enum CubitBoolean flg);
/* API function for GeometryTool::healer_set_inc_simplify_flg */

enum CubitBoolean GeometryTool_healer_get_inc_stitch_flg();
/* API function for GeometryTool::healer_get_inc_stitch_flg */

void GeometryTool_healer_set_inc_stitch_flg( enum CubitBoolean flg);
/* API function for GeometryTool::healer_set_inc_stitch_flg */

enum CubitBoolean GeometryTool_healer_get_inc_geombuild_flg();
/* API function for GeometryTool::healer_get_inc_geombuild_flg */

void GeometryTool_healer_set_inc_geombuild_flg( enum CubitBoolean flg);
/* API function for GeometryTool::healer_set_inc_geombuild_flg */

enum CubitBoolean GeometryTool_healer_get_inc_analytic_flg();
/* API function for GeometryTool::healer_get_inc_analytic_flg */

void GeometryTool_healer_set_inc_analytic_flg( enum CubitBoolean flg);
/* API function for GeometryTool::healer_set_inc_analytic_flg */

enum CubitBoolean GeometryTool_healer_get_inc_isospline_flg();
/* API function for GeometryTool::healer_get_inc_isospline_flg */

void GeometryTool_healer_set_inc_isospline_flg( enum CubitBoolean flg);
/* API function for GeometryTool::healer_set_inc_isospline_flg */

enum CubitBoolean GeometryTool_healer_get_inc_reblend_flg();
/* API function for GeometryTool::healer_get_inc_reblend_flg */

void GeometryTool_healer_set_inc_reblend_flg( enum CubitBoolean flg);
/* API function for GeometryTool::healer_set_inc_reblend_flg */

enum CubitBoolean GeometryTool_healer_get_inc_sharpedge_flg();
/* API function for GeometryTool::healer_get_inc_sharpedge_flg */

void GeometryTool_healer_set_inc_sharpedge_flg( enum CubitBoolean flg);
/* API function for GeometryTool::healer_set_inc_sharpedge_flg */

enum CubitBoolean GeometryTool_healer_get_inc_genericspline_flg();
/* API function for GeometryTool::healer_get_inc_genericspline_flg */

void GeometryTool_healer_set_inc_genericspline_flg( enum CubitBoolean flg);
/* API function for GeometryTool::healer_set_inc_genericspline_flg */

enum CubitBoolean GeometryTool_healer_get_inc_wrapup_flg();
/* API function for GeometryTool::healer_get_inc_wrapup_flg */

void GeometryTool_healer_set_inc_wrapup_flg( enum CubitBoolean flg);
/* API function for GeometryTool::healer_set_inc_wrapup_flg */

enum CubitBoolean GeometryTool_healer_get_inc_postprocess_flg();
/* API function for GeometryTool::healer_get_inc_postprocess_flg */

void GeometryTool_healer_set_inc_postprocess_flg( enum CubitBoolean flg);
/* API function for GeometryTool::healer_set_inc_postprocess_flg */

enum CubitStatus GeometryTool_force_simplify_to_plane( /* DLRefFaceList & */ void ***ref_face_list,
                                        int *ref_face_list_size,
                                          /* DLBodyList & */ void ***new_body_list,
                                        int *new_body_list_size, 
                                        enum CubitBoolean keep_old_body);
/* API function for GeometryTool::force_simplify_to_plane */

enum CubitStatus GeometryTool_force_simplify_to_cylinder( /* DLRefFaceList & */ void ***ref_face_list,
                                           int *ref_face_list_size,
                                             /* DLBodyList & */ void ***new_body_list,
                                           int *new_body_list_size, 
                                           enum CubitBoolean keep_old_body);
/* API function for GeometryTool::force_simplify_to_cylinder */

enum CubitStatus GeometryTool_force_simplify_to_cone( /* DLRefFaceList & */ void ***ref_face_list,
                                       int *ref_face_list_size,
                                         /* DLBodyList & */ void ***new_body_list,
                                       int *new_body_list_size, 
                                       enum CubitBoolean keep_old_body);
/* API function for GeometryTool::force_simplify_to_cone */

enum CubitStatus GeometryTool_force_simplify_to_sphere( /* DLRefFaceList & */ void ***ref_face_list,
                                         int *ref_face_list_size,
                                           /* DLBodyList & */ void ***new_body_list,
                                         int *new_body_list_size, 
                                         enum CubitBoolean keep_old_body);
/* API function for GeometryTool::force_simplify_to_sphere */

enum CubitStatus GeometryTool_force_simplify_to_torus( /* DLRefFaceList & */ void ***ref_face_list,
                                        int *ref_face_list_size,
                                          /* DLBodyList & */ void ***new_body_list,
                                        int *new_body_list_size, 
                                        enum CubitBoolean keep_old_body);
/* API function for GeometryTool::force_simplify_to_torus */

    /* <HR><H3>Merging functions.<H3> */
   
void GeometryTool_set_geometry_factor( double fac);
/* API function for GeometryTool::set_geometry_factor */

double GeometryTool_get_geometry_factor();
/* API function for GeometryTool::get_geometry_factor */

void GeometryTool_set_merge_test_bbox(enum CubitBoolean tof);
/* API function for GeometryTool::set_merge_test_bbox */

enum CubitBoolean GeometryTool_get_merge_test_bbox();
/* API function for GeometryTool::get_merge_test_bbox */

void GeometryTool_set_merge_test_internal(int tof);
/* API function for GeometryTool::set_merge_test_internal */

int GeometryTool_get_merge_test_internal();
/* API function for GeometryTool::get_merge_test_internal */

  enum CubitBoolean GeometryTool_same_engine(/* DLTopologyEntityList & */ void ***topo_list, int *topo_list_size);
  /**<  Returns CUBIT_TRUE if all the entities have the same geometry engine and
    *  if that is the same one as the default.
    */
    
    /* GeometricModelingEngine* */ void *GeometryTool_common_engine( 
        /* DLTopologyEntityList& */ void ***topology_list, int *topology_list_size,
          /* DLTopologyBridgeList& */ void ***engine_bridges, int *engine_bridges_size,
        enum CubitBoolean allow_virtual_engine);

  void GeometryTool_set_sep_after_webcut(enum CubitBoolean val);
  enum CubitBoolean GeometryTool_get_sep_after_webcut();
  /**< Gets/Sets the separate after webcut static flag.
    */ 
#ifdef __cplusplus
}
#endif

#endif
