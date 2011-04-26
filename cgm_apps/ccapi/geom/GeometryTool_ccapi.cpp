/**
 * \file geom_ccapi.cpp
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

#include "GeometryTool_ccapi.h"
#include "GeometryTool.hpp"
#include "CubitDefines.h"
#include "GeometryDefines.h"
#include "SolidModelingEngine.hpp"
#include "RefEntity.hpp"
#include "RefVertex.hpp"
#include "RefEdge.hpp"
#include "RefFace.hpp"
#include "RefVolume.hpp"
#include "RefAssembly.hpp"
#include "RefPart.hpp"
#include "RefCoordSys.hpp"
#include "Chain.hpp"
#include "Loop.hpp"
#include "Shell.hpp"
#include "Body.hpp"
#include "Lump.hpp"
#include "SurfaceSM.hpp"
#include "CurveSM.hpp"
#include "PointSM.hpp"
#include "BodySM.hpp"
#include "ShellSM.hpp"
#include "LoopSM.hpp"
#include "CoEdgeSM.hpp"
#include "CubitVector.hpp"
#include "CubitPlane.hpp"
#include "DLIList.hpp"
#include "DLPointList.hpp"
#include "DLBodyList.hpp"
#include "DLRefFaceList.hpp"
#include "DLRefEdgeList.hpp"
#include "DLRefVertexList.hpp"
#include "DLCubitEntityList.hpp"
#include "DLTopologyEntityList.hpp"
#include "CubitCompat.hpp"

#include "copy_defines.h"

/* GeometryTool * */ void *GeometryTool_instance(/* SolidModelingEngine * */ void *SMEPtr)
{
  SolidModelingEngine *temp_SMEPtr = (SolidModelingEngine *) SMEPtr;
  
  return GeometryTool::instance(temp_SMEPtr);
}
  
CubitStatus GeometryTool_cubit_entity_list(const char *keyword, 
                                   /* DLRefEntityList & */ void ***entity_list, 
                                 int *entity_list_size,
                                 const CubitBoolean print_errors)
/* < API function for GeometryTool::cubit_entity_list */
{

  DLCubitEntityList temp_list;
  CubitStatus status = GTI->cubit_entity_list(keyword, temp_list, 
                                              print_errors);
  
  COPY_LIST_TO_ARRAY(temp_list, *entity_list, *entity_list_size);
  
  return status;
}


CubitStatus GeometryTool_ref_entity_list( char const* keyword,
                                  /* DLRefEntityList &*/ void ***entity_list,
                                int *entity_list_size,
                                const CubitBoolean print_errors) 
{
/* < API function for GeometryTool::ref_entity_list */
  DLRefEntityList temp_list;

  CubitStatus status = GTI->ref_entity_list(keyword, temp_list, 
                                            print_errors);
  
  COPY_LIST_TO_ARRAY(temp_list, *entity_list, *entity_list_size);
  
  return status;
}

void GeometryTool_bodies (/*DLBodyList &*/ void ***bodies, int *bodies_size)
{
  DLBodyList temp_list;
  GTI->bodies(temp_list);
  COPY_LIST_TO_ARRAY(temp_list, *bodies, *bodies_size);
}

void GeometryTool_ref_volumes (/*DLRefVolumeList &*/ void ***ref_volumes, int *ref_volumes_size)
{
  DLRefVolumeList temp_list;
  GTI->ref_volumes(temp_list);
  COPY_LIST_TO_ARRAY(temp_list, *ref_volumes, *ref_volumes_size);
}

void GeometryTool_ref_groups (/*DLRefGroupList &*/ void ***ref_groups, int *ref_groups_size)
{
  DLRefGroupList temp_list;
  GTI->ref_groups(temp_list);
  COPY_LIST_TO_ARRAY(temp_list, *ref_groups, *ref_groups_size);
}

void GeometryTool_ref_parts (/*DLRefPartList &*/ void ***ref_parts, int *ref_parts_size)
{
  DLRefPartList temp_list;
  GTI->ref_parts(temp_list);
  COPY_LIST_TO_ARRAY(temp_list, *ref_parts, *ref_parts_size);
}

void GeometryTool_ref_assemblies (/*DLRefAssemblyList &*/ void ***ref_assemblies, int *ref_assemblies_size)
{
  DLRefAssemblyList temp_list;
  GTI->ref_assemblies(temp_list);
  COPY_LIST_TO_ARRAY(temp_list, *ref_assemblies, *ref_assemblies_size);
}

void GeometryTool_ref_faces (/*DLRefFaceList &*/ void ***ref_faces, int *ref_faces_size)
{
  DLRefFaceList temp_list;
  GTI->ref_faces(temp_list);
  COPY_LIST_TO_ARRAY(temp_list, *ref_faces, *ref_faces_size);
}

void GeometryTool_ref_edges (/*DLRefEdgeList &*/ void ***ref_edges, int *ref_edges_size)
{
  DLRefEdgeList temp_list;
  GTI->ref_edges(temp_list);
  COPY_LIST_TO_ARRAY(temp_list, *ref_edges, *ref_edges_size);
}

void GeometryTool_ref_vertices (/*DLRefVertexList &*/ void ***ref_vertices, int *ref_vertices_size)
{
  DLRefVertexList temp_list;
  GTI->ref_vertices(temp_list);
  COPY_LIST_TO_ARRAY(temp_list, *ref_vertices, *ref_vertices_size);
}

void GeometryTool_ref_coordsys (/*DLRefCoordSysList &*/ void ***ref_coordsys, int *ref_coordsys_size)
{
  DLRefCoordSysList temp_list;
  GTI->ref_coordsys(temp_list);
  COPY_LIST_TO_ARRAY(temp_list, *ref_coordsys, *ref_coordsys_size);
}

int GeometryTool_num_bodies() {return GTI->num_bodies();};
/* < API function for GeometryTool::num_bodies */

int GeometryTool_num_ref_volumes() {return GTI->num_ref_volumes();};
/* < API function for GeometryTool::num_ref_volumes */

int GeometryTool_num_ref_groups() {return GTI->num_ref_groups();};
/* < API function for GeometryTool::num_ref_groups */

int GeometryTool_num_ref_parts() {return GTI->num_ref_parts();};
/* < API function for GeometryTool::num_ref_parts */

int GeometryTool_num_ref_assemblies() {return GTI->num_ref_assemblies();};
/* < API function for GeometryTool::num_ref_assemblies */

int GeometryTool_num_ref_faces() {return GTI->num_ref_faces();};
/* < API function for GeometryTool::num_ref_faces */

int GeometryTool_num_ref_edges() {return GTI->num_ref_edges();};
/* < API function for GeometryTool::num_ref_edges */

int GeometryTool_num_ref_vertices() {return GTI->num_ref_vertices();};
/* < API function for GeometryTool::num_ref_vertices */

int GeometryTool_num_ref_coordsys() {return GTI->num_ref_coordsys();};
/* < API function for GeometryTool::num_ref_coordsys */


/* RefEntity* */ void* GeometryTool_get_ref_entity_1 (const char *type, int id) 
/* < API function for GeometryTool::*get_ref_entity  */
{return GTI->get_ref_entity (type, id);}

/* RefEntity* */ void* GeometryTool_get_ref_entity_2 (const EntityType type, int id) 
/* < API function for GeometryTool::*get_ref_entity  */
{return GTI->get_ref_entity (type, id);};


/* Body* */ void* GeometryTool_get_body ( int id ) {return GTI->get_body(id);};
/* < API function for GeometryTool::*get_body  */

/* RefVolume* */ void* GeometryTool_get_ref_volume ( int id ) {return GTI->get_ref_volume(id);};
/* < API function for GeometryTool::*get_ref_volume  */

/* RefGroup* */ void* GeometryTool_get_ref_group ( int id ) {return GTI->get_ref_group(id);};
/* < API function for GeometryTool::*get_ref_group  */

/* RefPart* */ void* GeometryTool_get_ref_part ( int id ) {return GTI->get_ref_part(id);};
/* < API function for GeometryTool::*get_ref_part  */

/* RefAssembly* */ void* GeometryTool_get_ref_assembly ( int id ) {return GTI->get_ref_assembly(id);};
/* < API function for GeometryTool::*get_ref_assembly  */

/* RefFace* */ void* GeometryTool_get_ref_face ( int id ) {return GTI->get_ref_face(id);};
/* < API function for GeometryTool::*get_ref_face  */

/* RefEdge* */ void* GeometryTool_get_ref_edge ( int id ) {return GTI->get_ref_edge(id);};
/* < API function for GeometryTool::*get_ref_edge  */

/* RefVertex* */ void* GeometryTool_get_ref_vertex ( int id ) {return GTI->get_ref_vertex(id);};
/* < API function for GeometryTool::*get_ref_vertex  */


/* Body* */ void* GeometryTool_get_first_body () {return GTI->get_first_body();}
/* < API function for GeometryTool::*get_first_body  */

/* RefVolume* */ void* GeometryTool_get_first_ref_volume () {return GTI->get_first_ref_volume();}
/* < API function for GeometryTool::*get_first_ref_volume  */

/* RefGroup* */ void* GeometryTool_get_first_ref_group () {return GTI->get_first_ref_group();}
/* < API function for GeometryTool::*get_first_ref_group  */

/* RefPart* */ void* GeometryTool_get_first_ref_part () {return GTI->get_first_ref_part();}
/* < API function for GeometryTool::*get_first_ref_part  */

/* RefAssembly* */ void* GeometryTool_get_first_ref_assembly () {return GTI->get_first_ref_assembly();}
/* < API function for GeometryTool::*get_first_ref_assembly  */

/* RefFace* */ void* GeometryTool_get_first_ref_face () {return GTI->get_first_ref_face();}
/* < API function for GeometryTool::*get_first_ref_face  */

/* RefEdge* */ void* GeometryTool_get_first_ref_edge () {return GTI->get_first_ref_edge();}
/* < API function for GeometryTool::*get_first_ref_edge  */

/* RefVertex* */ void* GeometryTool_get_first_ref_vertex () {return GTI->get_first_ref_vertex();}
/* < API function for GeometryTool::*get_first_ref_vertex  */
  

/* Body* */ void* GeometryTool_get_next_body () {return GTI->get_next_body();}
/* < API function for GeometryTool::*get_next_body  */

/* RefVolume* */ void* GeometryTool_get_next_ref_volume () {return GTI->get_next_ref_volume();}
/* < API function for GeometryTool::*get_next_ref_volume  */

/* RefGroup* */ void* GeometryTool_get_next_ref_group () {return GTI->get_next_ref_group();}
/* < API function for GeometryTool::*get_next_ref_group  */

/* RefPart* */ void* GeometryTool_get_next_ref_part () {return GTI->get_next_ref_part();}
/* < API function for GeometryTool::*get_next_ref_part  */

/* RefAssembly* */ void* GeometryTool_get_next_ref_assembly () {return GTI->get_next_ref_assembly();}
/* < API function for GeometryTool::*get_next_ref_assembly  */

/* RefFace* */ void* GeometryTool_get_next_ref_face () {return GTI->get_next_ref_face();}
/* < API function for GeometryTool::*get_next_ref_face  */

/* RefEdge* */ void* GeometryTool_get_next_ref_edge () {return GTI->get_next_ref_edge();}
/* < API function for GeometryTool::*get_next_ref_edge  */

/* RefVertex* */ void* GeometryTool_get_next_ref_vertex () {return GTI->get_next_ref_vertex();}
/* < API function for GeometryTool::*get_next_ref_vertex */
  

/* Body* */ void* GeometryTool_get_last_body () {return GTI->get_last_body();}
/* < API function for GeometryTool::*get_last_body  */

/* RefVolume* */ void* GeometryTool_get_last_ref_volume () {return GTI->get_last_ref_volume();}
/* < API function for GeometryTool::*get_last_ref_volume  */

/* RefGroup* */ void* GeometryTool_get_last_ref_group () {return GTI->get_last_ref_group();}
/* < API function for GeometryTool::*get_last_ref_group  */

/* RefPart* */ void* GeometryTool_get_last_ref_part () {return GTI->get_last_ref_part();}
/* < API function for GeometryTool::*get_last_ref_part  */

/* RefAssembly* */ void* GeometryTool_get_last_ref_assembly () {return GTI->get_last_ref_assembly();}
/* < API function for GeometryTool::*get_last_ref_assembly  */

/* RefFace* */ void* GeometryTool_get_last_ref_face () {return GTI->get_last_ref_face();}
/* < API function for GeometryTool::*get_last_ref_face  */

/* RefEdge* */ void* GeometryTool_get_last_ref_edge () {return GTI->get_last_ref_edge();}
/* < API function for GeometryTool::*get_last_ref_edge  */

/* RefVertex* */ void* GeometryTool_get_last_ref_vertex () {return GTI->get_last_ref_vertex();}
/* < API function for GeometryTool::*get_last_ref_vertex  */

  /**< 
    <HR><H3>File operations (import, export)</H3>
  */

CubitStatus GeometryTool_export_solid_model( /* DLRefEntityList& */ void ***ref_entity_list,
                                   int *ref_entity_list_size,
                                   char const* file_name,
                                   char const* file_type,
                                   const /* CubitString & */ char *cubit_version,
                                   char const* logfile_name) 
{
/* < API function for GeometryTool::export_solid_model */

  DLRefEntityList temp_list;
  COPY_ARRAY_TO_LIST(*ref_entity_list, *ref_entity_list_size, temp_list);

  return CubitCompat_export_solid_model(temp_list, file_name, file_type, 
                                 CubitString(cubit_version), logfile_name);
}

CubitStatus GeometryTool_import_solid_model_1(FILE* file_ptr, 
                                  const char* file_name,
                                  const char* file_type,
                                  const char* logfile_name,
                                  CubitBoolean heal_step,
                                  CubitBoolean import_bodies,
                                  CubitBoolean import_surfaces,
                                  CubitBoolean import_curves,
                                  CubitBoolean import_vertices,
                                  CubitBoolean free_surfaces)
{
  assert(file_ptr);
  CubitCompat_import_solid_model( file_name, file_type,
                                  logfile_name,
                                  heal_step,
                                  import_bodies, import_surfaces,
                                  import_curves, import_vertices,
                                  free_surfaces);
}

CubitStatus GeometryTool_import_solid_model_2(SolidModelerType model_type) 
/* < API function for GeometryTool::import_solid_model */
{
  return GTI->import_solid_model(model_type);
}
     
CubitStatus GeometryTool_import_datum_curves(FILE* inputFile, 
                                   char const* filetype)
/* < API function for GeometryTool::import_datum_curves */
{
  return GTI->import_datum_curves(inputFile, filetype);
}
   
CubitStatus GeometryTool_read_geometry_file(char const* fileName, 
                                  char const* includePath,
                                  char const* type)
/* < API function for GeometryTool::read_geometry_file */
{
  return GTI->read_geometry_file(fileName, includePath, type);
}
   
int GeometryTool_fire_ray(/* Body* */ void *body, 
                  /* const CubitVector */ CubitVectorStruct *ray_point,
                  /* const CubitVector */ CubitVectorStruct *unit,
                  /* int & */ int *num_hit, /* double *& */ double **ray_params,
                  /* DLRefEntityList * */ void ***entity_list,
                int *entity_list_size)
/* < API function for GeometryTool::fire_ray */
{
  CubitVector temp_ray_point(*ray_point);
  CubitVector temp_unit = *unit;
  DLRefEntityList temp_entity_list;
  Body *temp_body = (Body *) body;
  
  int result = GTI->fire_ray(temp_body, temp_ray_point, temp_unit, *num_hit,
                             *ray_params, &temp_entity_list);
  
  COPY_LIST_TO_ARRAY(temp_entity_list, *entity_list, *entity_list_size);
  
  return result;
}

  /**
   <HR><H3>Geometric primitive creation functions</H3>
  */
int GeometryTool_sphere(double radius) {return GTI->sphere(radius);}
/* < API function for GeometryTool::sphere */
   
int GeometryTool_brick( double wid, double dep, double hi ) {return GTI->brick(wid, dep, hi);}
/* < API function for GeometryTool::brick     */
   
int GeometryTool_prism( double height, int sides, double major, double minor) 
/* < API function for GeometryTool::prism     */
{
  return GTI->prism(height, sides, major, minor);
}
   
int GeometryTool_pyramid( double height, int sides, double major, double minor,
                double top)
/* < API function for GeometryTool::pyramid   */
{
  return GTI->pyramid(height, sides, major, minor, top);
}
   
int GeometryTool_cylinder( double hi, double r1, double r2, double r3 ) 
/* < API function for GeometryTool::cylinder  */
{
  return GTI->cylinder(hi, r1, r2, r3);
}
   
int GeometryTool_torus( double r1, double r2 ) {return GTI->torus(r1, r2);}
/* < API function for GeometryTool::torus     */

int GeometryTool_planar_sheet ( CubitPlaneStruct plane, CubitBoxStruct bounding_box,
                      int extension_type, double extension)
/* < API function for GeometryTool::planar_sheet  */
{
  CubitPlane temp_plane = plane;
  CubitBox temp_bounding_box = bounding_box;
  
  return GTI->planar_sheet(temp_plane, temp_bounding_box,
                           extension_type, extension);
}
   
  /* / <HR><H3>Geometry transformation (move, scale, etc.)</H3> */

/* Body* */ void* GeometryTool_copy_body ( /* Body* */ void *body_ptr)
/* < API function for GeometryTool::copy_body  */
{
  Body *temp_body_ptr = (Body *) body_ptr;
  return GTI->copy_body(temp_body_ptr);
}

   
/* RefEntity* */ void* GeometryTool_copy_refentity( /* RefEntity * */ void *old_entity )
/* < API function for GeometryTool::copy_refentity */
{
  RefEntity *temp_old_entity = (RefEntity *) old_entity;
  return GTI->copy_refentity(temp_old_entity);
}
   
CubitStatus GeometryTool_reflect_body( /* Body * */ void *body_ptr,
                               /* Body *& */ void **new_body,
                             double x, double y,
                             double z,
                             CubitBoolean keep_old)
/* < API function for GeometryTool::reflect_body */
{
  Body *temp_body_ptr = (Body *) body_ptr;
  Body *temp_new_body;
  CubitStatus result = GTI->reflect_body(temp_body_ptr, temp_new_body, x, y, z, keep_old);
  *new_body = temp_new_body;
  return result;
}

CubitStatus GeometryTool_align_body( /* Body * */ void *body_ptr, 
                             /* RefFace * */ void *my_face,
                             /* RefFace * */ void *target_face )
/* < API function for GeometryTool::align_body */
{
  Body *temp_body_ptr = (Body *) body_ptr;
  RefFace *temp_my_face = (RefFace *) my_face;
  RefFace *temp_target_face = (RefFace *) target_face;
  return GTI->align_body(temp_body_ptr, temp_my_face, temp_target_face);
}

   
  /* /< <HR><H3>Geometry modification (booleans, decomposition, etc.)</H3> */

CubitStatus GeometryTool_unite_1(/* DLBodyList & */ void ***bodies, int *bodies_size, 
                       /* Body *& */ void **newBody,
                     int keep_old)
/* < API function for GeometryTool::unite */
{
  DLBodyList temp_bodies;
  COPY_ARRAY_TO_LIST(*bodies, *bodies_size, temp_bodies);
  Body *temp_newBody;
  
  CubitStatus result = GTI->unite(temp_bodies, temp_newBody, keep_old);

  *newBody = temp_newBody;
  return result;
}

CubitStatus GeometryTool_unite_2 ( /* Body * */ void *body1, /* Body * */ void *body2, 
                         /* Body *& */ void **newBody,
                       int keep_old)
/* < API function for GeometryTool::unite  */
{
  Body *temp_body1 = (Body *) body1;
  Body *temp_body2 = (Body *) body2;
  Body *temp_newBody;
  
  CubitStatus result = GTI->unite(temp_body1, temp_body2, temp_newBody, keep_old);
  *newBody = temp_newBody;
  
  return result;
}

CubitStatus GeometryTool_subtract ( /* Body * */ void *tool_body_ptr, 
                            /* DLBodyList & */ void ***from_bodies,
                          int *from_bodies_size,
                            /* DLBodyList & */ void ***new_bodies, int *new_bodies_size,
                          int keep_old)
/* < API function for GeometryTool::subtract  */
{
  Body *temp_tool_body_ptr = (Body *) tool_body_ptr;
  DLBodyList temp_from_bodies;
  COPY_ARRAY_TO_LIST(*from_bodies, *from_bodies_size, temp_from_bodies);

  DLBodyList temp_new_bodies;

  CubitStatus result = GTI->subtract(temp_tool_body_ptr, temp_from_bodies, 
                                     temp_new_bodies, keep_old);
  
  COPY_LIST_TO_ARRAY(temp_new_bodies, *new_bodies, *new_bodies_size);
  
  return result;
}

CubitStatus GeometryTool_intersect ( /* Body * */ void *tool_body_ptr, 
                             /* DLBodyList & */ void ***from_bodies,
                           int *from_bodies_size,
                             /* DLBodyList & */ void ***new_bodies,
                           int *new_bodies_size,
                           int keep_old)
/* < API function for GeometryTool::intersect  */
{
  Body *temp_tool_body_ptr = (Body *) tool_body_ptr;
  DLBodyList temp_from_bodies;
  COPY_ARRAY_TO_LIST(*from_bodies, *from_bodies_size, temp_from_bodies);

  DLBodyList temp_new_bodies;
  
  CubitStatus result = GTI->intersect(temp_tool_body_ptr, temp_from_bodies, 
                                      temp_new_bodies, keep_old);
  
  COPY_LIST_TO_ARRAY(temp_new_bodies, *new_bodies, *new_bodies_size);
  
  return result;
}

CubitStatus GeometryTool_section( /* DLBodyList & */ void ***section_body_list,
                        int *section_body_list_size,
                          /* CubitVector & */ CubitVectorStruct *point_1,
                          /* CubitVector & */ CubitVectorStruct *point_2,
                          /* CubitVector & */ CubitVectorStruct *point_3,
                          /* DLBodyList & */ void ***new_body_list,
                        int *new_body_list_size,
                        CubitBoolean keep_normal_side,
                        CubitBoolean keep_old)
/* < API function for GeometryTool::section */
{
  DLBodyList temp_section_body_list;
  COPY_ARRAY_TO_LIST(*section_body_list, *section_body_list_size, temp_section_body_list);

  CubitVector temp_point_1 = *point_1;
  CubitVector temp_point_2 = *point_2;
  CubitVector temp_point_3 = *point_3;

  DLBodyList temp_new_body_list;

  CubitStatus result = GTI->section(temp_section_body_list, 
                                    temp_point_1, temp_point_2, 
                                    temp_point_3, temp_new_body_list, 
                                    keep_normal_side, keep_old);
  
  COPY_LIST_TO_ARRAY(temp_new_body_list, *new_body_list, *new_body_list_size);
  
  return result;
}

CubitStatus GeometryTool_split_periodic( /* Body * */ void *body_ptr,
                                 /* Body *& */ void **new_body_ptr )
/* < API function for GeometryTool::split_periodic */
{
  Body *temp_body_ptr = (Body *) body_ptr;
  Body *temp_new_body_ptr;

  CubitStatus result = GTI->split_periodic(temp_body_ptr, temp_new_body_ptr);
  
  *new_body_ptr = temp_new_body_ptr;
  
  return result;
}
   
  int  GeometryTool_webcut_across_translate( /* DLBodyList& */ void ***body_list, int *body_list_size,
                                               /* RefFace* */ void *plane_surf1, /* RefFace* */ void *plane_surf2,
                                               /* DLBodyList& */ void ***results_list, int *results_list_size,
                                             enum CubitBoolean imprint, enum CubitBoolean merge) 
{

  DLBodyList temp_body_list;
  RefFace *temp_plane_surf1 = (RefFace *) plane_surf1;
  RefFace *temp_plane_surf2 = (RefFace *) plane_surf2;
  DLBodyList temp_results_list;

  COPY_ARRAY_TO_LIST(*body_list, *body_list_size, temp_body_list);
  
  int result = GTI->webcut_across_translate(temp_body_list, temp_plane_surf1,
                                            temp_plane_surf2, temp_results_list,
                                            imprint, merge);
  
  COPY_LIST_TO_ARRAY(temp_results_list, *results_list, *results_list_size);
  
  return result;
}
   
int GeometryTool_webcut_with_cylinder( /* DLBodyList& */ void ***webcut_body_list,
                             int *webcut_body_list_size,
                             double radius, 
                               /* CubitVector & */ CubitVectorStruct *axis,
                               /* CubitVector& */ CubitVectorStruct* center, 
                               /* DLBodyList & */ void ***results_list,
                             int *results_list_size,
                             CubitBoolean imprint,
                             CubitBoolean merge)
/* < API function for GeometryTool::webcut_with_cylinder */
{
  DLBodyList temp_webcut_body_list;
  COPY_ARRAY_TO_LIST(*webcut_body_list, *webcut_body_list_size, temp_webcut_body_list);

  CubitVector temp_axis = *axis;
  CubitVector temp_center = *center;
  
  DLBodyList temp_results_list;
  int result = GTI->webcut_with_cylinder(temp_webcut_body_list, radius, 
                                         temp_axis, temp_center, 
                                         temp_results_list, 
                                         imprint, merge);

  COPY_LIST_TO_ARRAY(temp_results_list, *results_list, *results_list_size);

  return result;
}

int GeometryTool_webcut_with_plane(/* DLBodyList& */ void ***webcut_body_list,
                         int *webcut_body_list_size,
                           /* CubitVector & */ CubitVectorStruct *vector1,
                           /* CubitVector & */ CubitVectorStruct *vector2,
                           /* CubitVector & */ CubitVectorStruct *vector3, 
                           /* DLBodyList & */ void ***results_list,
                         int *results_list_size,
                         CubitBoolean imprint,
                         CubitBoolean merge) 
/* < API function for GeometryTool::webcut_with_plane */
{
  CubitVector temp_vector1 = *vector1;
  CubitVector temp_vector2 = *vector2;
  CubitVector temp_vector3 = *vector3;
  
  DLBodyList temp_webcut_body_list;
  COPY_ARRAY_TO_LIST(*webcut_body_list, *webcut_body_list_size, temp_webcut_body_list);
  
  DLBodyList temp_results_list;
  
  int result = GTI->webcut_with_plane(temp_webcut_body_list, 
                                      temp_vector1, temp_vector2, temp_vector3, 
                                      temp_results_list, imprint, merge);

  COPY_LIST_TO_ARRAY(temp_results_list, *results_list, *results_list_size);
  
  return result;
}
   
int GeometryTool_webcut_with_vertices( /* DLBodyList& */ void ***webcut_body_list,
                             int *webcut_body_list_size,
                               /* RefVertex* */ void *refVertex1,
                               /* RefVertex* */ void *refVertex2,
                               /* RefVertex* */ void *refVertex3, 
                               /* DLBodyList & */ void ***results_list,
                             int *results_list_size,
                             CubitBoolean imprint,
                             CubitBoolean merge)
/* < API function for GeometryTool::webcut_with_vertices */
{
  DLBodyList temp_webcut_body_list;
  COPY_ARRAY_TO_LIST(*webcut_body_list, *webcut_body_list_size, temp_webcut_body_list);

  RefVertex *temp_refVertex1 = (RefVertex *) refVertex1;
  RefVertex *temp_refVertex2 = (RefVertex *) refVertex2;
  RefVertex *temp_refVertex3 = (RefVertex *) refVertex3;
  
  DLBodyList temp_results_list;
  
  int result = GTI->webcut_with_vertices(temp_webcut_body_list, 
                                         temp_refVertex1, temp_refVertex2, temp_refVertex3, 
                                         temp_results_list, imprint, merge);

  COPY_LIST_TO_ARRAY(temp_results_list, *results_list, *results_list_size);
  
  return result;
}
   
int GeometryTool_webcut_with_curve_loop(/* DLBodyList& */ void ***webcut_body_list, int *webcut_body_list_size,
                                          /* DLRefEdgeList& */ void ***refedge_list, int *refedge_list_size,
                                            /* DLBodyList& */ void ***results_list, int *results_list_size,
                                          enum CubitBoolean imprint)
{
  DLBodyList temp_webcut_body_list, temp_results_list;
  DLRefEdgeList temp_refedge_list;
  
  COPY_ARRAY_TO_LIST(*webcut_body_list, *webcut_body_list_size, temp_webcut_body_list);
  COPY_ARRAY_TO_LIST(*refedge_list, *refedge_list_size, temp_refedge_list);

  int result = GTI->webcut_with_curve_loop(temp_webcut_body_list, temp_refedge_list, 
                                           temp_results_list, imprint);

  COPY_LIST_TO_ARRAY(temp_results_list, *results_list, *results_list_size);

  return result;
}

int GeometryTool_webcut_with_surface(/* DLBodyList& */ void ***webcut_body_list,
                           int *webcut_body_list_size,
                             /* RefFace* */ void *refFace,
                             /* DLBodyList& */ void *** results_list,
                           int * results_list_size,
                           CubitBoolean imprint,
                           CubitBoolean merge)
/* < API function for GeometryTool::webcut_with_surface */
{
  DLBodyList temp_webcut_body_list;
  COPY_ARRAY_TO_LIST(*webcut_body_list, *webcut_body_list_size, temp_webcut_body_list);
  
  RefFace *temp_refFace = (RefFace *) refFace;
  
  DLBodyList temp_results_list;
  
  int result = GTI->webcut_with_surface(temp_webcut_body_list, temp_refFace, 
                                        temp_results_list, 
                                        imprint, merge);

  COPY_LIST_TO_ARRAY(temp_results_list, *results_list, *results_list_size);
  
  return result;
}
   
CubitStatus GeometryTool_webcut_with_sheet(/* Body * */ void *webcut_body,
                                   /* Body * */ void *sheet_body,
                              /* DLBodyList & */ void ***new_bodies,
                                 int *new_bodies_size,
                                 CubitBoolean imprint )
/* < API function for GeometryTool::webcut_with_sheet */
{
  DLBodyList temp_new_bodies;

  Body *temp_webcut_body = (Body *) webcut_body;
  Body *temp_sheet_body = (Body *) sheet_body;
  
  CubitStatus result = GTI->webcut_with_sheet(temp_webcut_body, temp_sheet_body, 
                                              temp_new_bodies, imprint);

  COPY_LIST_TO_ARRAY(temp_new_bodies, *new_bodies, *new_bodies_size);
  
  return result;
}
   
int GeometryTool_webcut_with_body(/* DLBodyList& */ void ***webcut_body_list,
                        int *webcut_body_list_size,
                          /* Body* */ void *body,
                          /* DLBodyList& */ void ***results_list,
                        int * results_list_size,
                        CubitBoolean imprint,
                        CubitBoolean merge) 
/* < API function for GeometryTool::webcut_with_body */
{
  DLBodyList temp_webcut_body_list;
  COPY_ARRAY_TO_LIST(*webcut_body_list, *webcut_body_list_size, temp_webcut_body_list);

  Body *temp_body = (Body *) body;
  
  DLBodyList temp_results_list;
  
  int result = GTI->webcut_with_body(temp_webcut_body_list, temp_body, 
                                     temp_results_list, 
                                     imprint, merge);

  COPY_LIST_TO_ARRAY(temp_results_list, *results_list, *results_list_size);
  
  return result;
}
   
CubitStatus GeometryTool_webcut_with_extended_surf( /* DLBodyList & */ void ***webcut_body_list,
                                          int *webcut_body_list_size,
                                            /* RefFace * */ void *face_to_extend,
                                            /* DLBodyList & */ void ***new_bodies,
                                          int *new_bodies_size,
                                            /* int & */ int *num_cut,
                                          CubitBoolean imprint )
/* < API function for GeometryTool::webcut_with_extended_surf */
{
  DLBodyList temp_webcut_body_list;
  COPY_ARRAY_TO_LIST(*webcut_body_list, *webcut_body_list_size, temp_webcut_body_list);

  RefFace *temp_face_to_extend = (RefFace *) face_to_extend;
  
  DLBodyList temp_new_bodies;
  
  CubitStatus result = GTI->webcut_with_extended_surf(temp_webcut_body_list, 
                                                      temp_face_to_extend, 
                                                      temp_new_bodies, *num_cut, 
                                                      imprint);

  COPY_LIST_TO_ARRAY(temp_new_bodies, *new_bodies, *new_bodies_size);
  
  return result;
}

CubitStatus GeometryTool_remove_surfaces( /* DLRefFaceList & */ void ***ref_face_list,
                                int *ref_face_list_size,
                                  /* DLBodyList & */ void ***new_body_list,
                                int *new_body_list_size,
                                CubitBoolean extend_adjoining,
                                CubitBoolean keep_surface,
                                CubitBoolean keep_old_body )
/* < API function for GeometryTool::remove_surfaces */
{
  DLRefFaceList temp_ref_face_list;
  COPY_ARRAY_TO_LIST(*ref_face_list, *ref_face_list_size, temp_ref_face_list);
  
  DLBodyList temp_new_body_list;
  
  CubitStatus result = GTI->remove_surfaces(temp_ref_face_list, 
                                            temp_new_body_list, 
                                            extend_adjoining, keep_surface, keep_old_body);

  COPY_LIST_TO_ARRAY(temp_new_body_list, *new_body_list, *new_body_list_size);
  
  return result;
}

CubitStatus GeometryTool_offset_surfaces( /* DLRefFaceList& */ void ***ref_face_list,
                                int * ref_face_list_size,
                                  /* DLBodyList& */ void ***new_body_list,
                                int * new_body_list_size, double offset_distance,
                                CubitBoolean keep_old_body )
/* < API function for GeometryTool::offset_surfaces */
{
  DLRefFaceList temp_ref_face_list;
  COPY_ARRAY_TO_LIST(*ref_face_list, *ref_face_list_size, temp_ref_face_list);
  
  DLBodyList temp_new_body_list;
  
  CubitStatus result = GTI->offset_surfaces(temp_ref_face_list, 
                                            temp_new_body_list, 
                                            offset_distance, keep_old_body);

  COPY_LIST_TO_ARRAY(temp_new_body_list, *new_body_list, *new_body_list_size);
  
  return result;
}

CubitStatus GeometryTool_offset_curves( /* DLRefEdgeList& */ void ***ref_edge_list,
                              int *ref_edge_list_size, double offset_distance, 
                              /* CubitVector& */ CubitVectorStruct* offset_direction, 
                              int gap_type)
/* < API function for GeometryTool::offset_curves */
{
  DLRefEdgeList temp_ref_edge_list;
  COPY_ARRAY_TO_LIST(*ref_edge_list, *ref_edge_list_size, temp_ref_edge_list);

  CubitVector temp_offset_direction = *offset_direction;
  
  CubitStatus result = GTI->offset_curves(temp_ref_edge_list, 
                                          offset_distance, temp_offset_direction, 
                                          gap_type);
  return result;
}

CubitStatus GeometryTool_translate_surfaces( /* DLRefFaceList& */ void ***ref_face_list,
                                   int *ref_face_list_size,
                                     /* DLBodyList& */ void ***new_body_list,
                                   int *new_body_list_size,
                                     /* CubitVector & */ CubitVectorStruct *delta,
                                   CubitBoolean keep_old_body )
/* < API function for GeometryTool::translate_surfaces */
{
  DLRefFaceList temp_ref_face_list;
  COPY_ARRAY_TO_LIST(* ref_face_list, *ref_face_list_size, temp_ref_face_list);

  CubitVector temp_delta = *delta;
  
  DLBodyList temp_new_body_list;
  
  CubitStatus result = GTI->translate_surfaces(temp_ref_face_list, 
                                               temp_new_body_list, 
                                               temp_delta, keep_old_body);

  COPY_LIST_TO_ARRAY(temp_new_body_list, *new_body_list, *new_body_list_size);
  
  return result;
}
  
CubitStatus GeometryTool_replace_surfaces( /* DLRefFaceList & */ void ***ref_face_list,
                                 int *ref_face_list_size,
                                   /* RefFace* */ void *tool_face_ptr,
                              /* DLBodyList & */ void ***new_body_list,
                                 int *new_body_list_size,
                                 CubitBoolean reverse_flg,
                                 CubitBoolean keep_old_body )
/* < API function for GeometryTool::replace_surfaces */
{
  DLRefFaceList temp_ref_face_list;
  COPY_ARRAY_TO_LIST(*ref_face_list, *ref_face_list_size, temp_ref_face_list);
  
  RefFace *temp_tool_face_ptr = (RefFace *) tool_face_ptr;
  
  DLBodyList temp_new_body_list;
  
  CubitStatus result = GTI->replace_surfaces(temp_ref_face_list, 
                                             temp_tool_face_ptr, 
                                             temp_new_body_list, 
                                             reverse_flg, keep_old_body);

  COPY_LIST_TO_ARRAY(temp_new_body_list, *new_body_list, *new_body_list_size);
  
  return result;
}
   
CubitStatus GeometryTool_trim_curve( /* RefEdge* */ void *trim_curve, 
                             /* CubitVector& */ CubitVectorStruct* trim_vector, 
                             /* CubitVector& */ CubitVectorStruct* keep_vector )
/* < API function for GeometryTool::trim_curve */
{
  RefEdge *temp_trim_curve = (RefEdge *) trim_curve;
  
  CubitVector temp_trim_vector = *trim_vector;
  CubitVector temp_keep_vector = *keep_vector;
  
  CubitStatus result = GTI->trim_curve(temp_trim_curve, 
                                       temp_trim_vector, temp_keep_vector);
  return result;
}
   
  /* /< <HR><H3>Topology modification (imprint, regularize, etc.)</H3> */

CubitStatus GeometryTool_imprint_1 (/* DLBodyList & */ void ***from_body_list,
                        int *from_body_list_size,
                          /* DLBodyList & */ void ***new_body_list,
                        int *new_body_list_size,
                        int keep_old)
/* < API function for GeometryTool::imprint  */
{
  DLBodyList temp_from_body_list;
  COPY_ARRAY_TO_LIST(*from_body_list, *from_body_list_size, temp_from_body_list);
  
  DLBodyList temp_new_body_list;
  
  CubitStatus result = GTI->imprint (temp_from_body_list, 
                                     temp_new_body_list, keep_old);

  COPY_LIST_TO_ARRAY(temp_new_body_list, *new_body_list, *new_body_list_size);
  
  return result;
}

CubitStatus GeometryTool_imprint_2( /* DLBodyList & */ void ***body_list,
                        int *body_list_size,
                          /* DLRefEdgeList & */ void ***ref_edge_list,
                        int *ref_edge_list_size,
                          /* DLBodyList& */ void ***new_body_list,
                        int *new_body_list_size,
                        int keep_old_body )
/* < API function for GeometryTool::imprint */
{
  DLBodyList temp_body_list;
  COPY_ARRAY_TO_LIST(*body_list, *body_list_size, temp_body_list);
  
  DLRefEdgeList temp_ref_edge_list;
  COPY_ARRAY_TO_LIST(*ref_edge_list, *ref_edge_list_size, temp_ref_edge_list);
  
  DLBodyList temp_new_body_list;
  
  CubitStatus result = GTI->imprint(temp_body_list, temp_ref_edge_list, 
                                    temp_new_body_list, keep_old_body);

  COPY_LIST_TO_ARRAY(temp_new_body_list, *new_body_list, *new_body_list_size);
  
  return result;
}

CubitStatus GeometryTool_imprint_3( /* DLRefFaceList & */ void ***ref_face_list,
                        int *ref_face_list_size,
                          /* DLRefEdgeList & */ void ***ref_edge_list,
                        int *ref_edge_list_size,
                          /* DLBodyList& */ void ***new_body_list,
                        int *new_body_list_size,
                        int keep_old_body )
/* < API function for GeometryTool::imprint */
{
  DLRefFaceList temp_ref_face_list;
  COPY_ARRAY_TO_LIST(*ref_face_list, *ref_face_list_size, temp_ref_face_list);
  
  DLRefEdgeList temp_ref_edge_list;
  COPY_ARRAY_TO_LIST(*ref_edge_list, *ref_edge_list_size, temp_ref_edge_list);
  
  DLBodyList temp_new_body_list;
  
  CubitStatus result = GTI->imprint(temp_ref_face_list, temp_ref_edge_list, 
                                    temp_new_body_list, keep_old_body);

  COPY_LIST_TO_ARRAY(temp_new_body_list, *new_body_list, *new_body_list_size);
  
  return result;
}

CubitStatus GeometryTool_imprint_4( /* DLBodyList & */ void ***body_list,
                        int *body_list_size,
                          /* DLRefVertexList & */ void ***ref_vertex_list,
                        int *ref_vertex_list_size,
                          /* DLBodyList& */ void ***new_body_list,
                        int *new_body_list_size,
                        int keep_old_body )
/* < API function for GeometryTool::imprint */
{
  DLBodyList temp_body_list;
  COPY_ARRAY_TO_LIST(*body_list, *body_list_size, temp_body_list);
  
  DLRefVertexList temp_ref_vertex_list;
  COPY_ARRAY_TO_LIST(*ref_vertex_list, *ref_vertex_list_size, temp_ref_vertex_list);
  
  DLBodyList temp_new_body_list;
  
  CubitStatus result = GTI->imprint(temp_body_list, temp_ref_vertex_list, 
                                    temp_new_body_list, keep_old_body);

  COPY_LIST_TO_ARRAY(temp_new_body_list, *new_body_list, *new_body_list_size);
  
  return result;
}

CubitStatus GeometryTool_regularize_body( /* Body * */ void *body_ptr, 
                                  /* Body *& */ void **new_body )
/* < API function for GeometryTool::regularize_body */
{
  Body *temp_body_ptr = (Body *) body_ptr;
  
  Body *temp_new_body;
  
  CubitStatus result = GTI->regularize_body( temp_body_ptr, temp_new_body);

  *new_body = temp_new_body;
  
  return result;
}
   
CubitStatus GeometryTool_split_body( /* Body * */ void *body_ptr,
                             /* DLBodyList & */ void ***new_bodies,
                           int *new_bodies_size)
/* < API function for GeometryTool::split_body */
{
  Body *temp_body_ptr = (Body *) body_ptr;
  
  DLBodyList temp_new_bodies;
  
  CubitStatus result = GTI->split_body(temp_body_ptr,
                                       temp_new_bodies);

  COPY_LIST_TO_ARRAY(temp_new_bodies, *new_bodies, *new_bodies_size);
  
  return result;
}
   
  /* /< <HR><H3>GeometryTool options and settings</H3> */

void GeometryTool_group_imprint_1(CubitBoolean flag) {GTI->group_imprint(flag);}
/* < API function for GeometryTool::group_imprint */

CubitBoolean GeometryTool_group_imprint_2() {return GTI->group_imprint();}
/* < API function for GeometryTool::group_imprint */
   
void GeometryTool_set_all_edges_imprint( CubitBoolean flag) {GTI->set_all_edges_imprint(flag);}
/* < API function for GeometryTool::void set_all_edges_imprint */

CubitBoolean GeometryTool_get_all_edges_imprint() {return GTI->get_all_edges_imprint();}
/* < API function for GeometryTool::CubitBoolean get_all_edges_imprint */
   
    /* static */ void GeometryTool_booleans_after_merge( enum CubitBoolean flag ) {GTI->booleans_after_merge(flag);}
    /* static */ enum CubitBoolean GeometryTool_booleans_after_merge() {return GTI->booleans_after_merge();}

void GeometryTool_new_ids_1(CubitBoolean flag) {GTI->new_ids(flag);}

/* < API function for GeometryTool::new_ids */

CubitBoolean GeometryTool_new_ids_2() {return GTI->new_ids();}
/* < API function for GeometryTool::new_ids */
   
void GeometryTool_geom_debug( /* DLTopologyEntityList */ void ***arg, int *arg_size) 
{
  DLTopologyEntityList temp_arg;
  COPY_ARRAY_TO_LIST(*arg, *arg_size, temp_arg);
  
  GTI->geom_debug(temp_arg);

}

void GeometryTool_use_facet_bbox_1( enum CubitBoolean pass_flag ) {GTI->use_facet_bbox(pass_flag);}
enum CubitBoolean GeometryTool_use_facet_bbox_2() {return GTI->use_facet_bbox();}
/* < API function for GeometryTool::void geom_debug */
   
#ifdef CUBIT_GUI
void GeometryTool_set_validate_file_ptr( FILE* file_ptr) {GTI->set_validate_file_ptr(file_ptr);}
/* < API function for GeometryTool::set_validate_file_ptr */

void GeometryTool_reset_validate_file_ptr() {GTI->reset_validate_file_ptr()}
/* < API function for GeometryTool::reset_validate_file_ptr */

FILE* GeometryTool_get_validate_file_ptr() {GTI->get_validate_file_ptr()}
/* < API function for GeometryTool::get_validate_file_ptr */
#endif
   
  /* /< <HR><H3>SolidModelingEngine information</H3> */

char *GeometryTool_identify_modeling_engine()
/* < API function for GeometryTool::identify_modeling_engine */
{
  CubitString temp_string = GTI->identify_modeling_engine();
  return (char *) temp_string.c_str();
}
   
void GeometryTool_register_solid_modeling_engine(/* SolidModelingEngine* */ void *SMEPtr)  
/* < API function for GeometryTool::register_solid_modeling_engine */
{
  SolidModelingEngine *temp_SMEPtr = (SolidModelingEngine *) SMEPtr;
  GTI->register_solid_modeling_engine(temp_SMEPtr);
}
     
void GeometryTool_register_datum_curve_sm_engine(/* SolidModelingEngine* */ void *SMEPtr) 
/* < API function for GeometryTool::register_datum_curve_sm_engine */
{
  SolidModelingEngine *temp_SMEPtr = (SolidModelingEngine *) SMEPtr;
  GTI->register_datum_curve_sm_engine(temp_SMEPtr);
}

CubitStatus GeometryTool_set_engine_version(int version) {return GTI->set_engine_version(version);}
/* < API function for GeometryTool::set_engine_version */
   
CubitStatus GeometryTool_get_engine_version(int *version) {return GTI->get_engine_version(*version);}
/* < API function for GeometryTool::get_engine_version */
   
CubitStatus GeometryTool_list_engine_versions(char **versions) 
/* < API function for GeometryTool::list_engine_versions */
{
  CubitString temp_versions;
  
  CubitStatus result = GTI->list_engine_versions(temp_versions);

  *versions = (char *) temp_versions.c_str();

  return result;
}

double GeometryTool_get_sme_resabs_tolerance() {return GTI->get_sme_resabs_tolerance();}
/* < API function for GeometryTool::get_sme_resabs_tolerance */

double GeometryTool_set_sme_resabs_tolerance( double new_resabs) {return GTI->set_sme_resabs_tolerance( new_resabs);}
/* < API function for GeometryTool::set_sme_resabs_tolerance */

CubitStatus GeometryTool_set_steptools_path( const char* path) {return GTI->set_steptools_path(path);}
/* < API function for GeometryTool::set_steptools_path */

CubitStatus GeometryTool_set_igestools_path( const char* path) {return GTI->set_igestools_path(path);}
/* < API function for GeometryTool::set_igestools_path */

CubitStatus GeometryTool_set_sme_int_option( const char* opt_name, int val) {return GTI->set_sme_int_option(opt_name, val);}
/* < API function for GeometryTool::set_sme_int_option */

CubitStatus GeometryTool_set_sme_dbl_option( const char* opt_name, double val) {return GTI->set_sme_dbl_option(opt_name, val);}
/* < API function for GeometryTool::set_sme_dbl_option */

CubitStatus GeometryTool_set_sme_str_option( const char* opt_name, const char* val) {return GTI->set_sme_str_option(opt_name, val);}
/* < API function for GeometryTool::set_sme_str_option */

  /* /< <HR><H3>Topology/geometry creation functions</H3> */

/* RefVertex* */ void* GeometryTool_make_RefVertex(GeometryType ref_vertex_type,
                                         const CubitVectorStruct *point, int color)
/* < API function for GeometryTool::make_RefVertex */
{
  CubitVector temp_point = *point;
  
  return GTI->make_RefVertex(ref_vertex_type, temp_point, color);
}
   
/* Body* */ void* GeometryTool_make_Body_1(/* Surface * */ void *surface)
/* < API function for GeometryTool::*make_Body */
{
  Surface *temp_surface = (Surface *) surface;
  
  return GTI->make_Body(temp_surface);
}
/* Body* */ void* GeometryTool_make_Body_2(/* BodySM * */ void *bodysm_ptr)
/* < API function for GeometryTool::make_Body */
{
  BodySM *temp_bodysm_ptr = (BodySM *) bodysm_ptr;
  
  return GTI->make_Body(temp_bodysm_ptr);
}

/* Shell* */ void* GeometryTool_make_Shell(/* ShellSM * */ void *shellsm_ptr )
/* < API function for GeometryTool::make_Shell */
{
  ShellSM *temp_shellsm_ptr = (ShellSM *) shellsm_ptr;
  
  return GTI->make_Shell(temp_shellsm_ptr);
}

/* Loop* */ void* GeometryTool_make_Loop(/* LoopSM* */ void *loopsm_ptr )
/* < API function for GeometryTool::make_Loop */
{
  LoopSM *temp_loopsm_ptr = (LoopSM *) loopsm_ptr ;
  
  return GTI->make_Loop(temp_loopsm_ptr);
}

  /* RefEntity * */ void *GeometryTool_check_mergeable_refentity(/* TopologyBridge * */ void *bridge) 
{
  TopologyBridge *temp_bridge = (TopologyBridge *) bridge;
  return GTI->check_mergeable_refentity(temp_bridge);
}

/* Chain* */ void* GeometryTool_make_Chain_1(/* DLPointList& */ void ***points, 
                                 int *points_size)
/* < API function for GeometryTool::make_Chain */
{
  DLPointList temp_points;
  COPY_ARRAY_TO_LIST(*points, *points_size, temp_points);
  
  return GTI->make_Chain(temp_points);
}

/* Chain* */ void* GeometryTool_make_Chain_2(/* Curve * */ void *curve)
/* < API function for GeometryTool::make_Chain */
{
  Curve *temp_curve = (Curve *) curve;
  
  return GTI->make_Chain(temp_curve);
}

/* RefEdge* */ void* GeometryTool_make_RefEdge_1( GeometryType ref_edge_type,
                                        /* RefVertex const* */ void *ref_vertex_1,
                                        /* RefVertex const* */ void *ref_vertex_2,
                                        /* DLCubitVectorList& */ CubitVectorStruct *vector_list,
                                      int vector_list_size,
                                        /* RefFace* */ void *refface_ptr ) 
/* < API function for GeometryTool::make_RefEdge */
{
  RefVertex *temp_ref_vertex_1 = (RefVertex *) ref_vertex_1;
  RefVertex *temp_ref_vertex_2 = (RefVertex *) ref_vertex_2;

  DLCubitVectorList temp_vector_list;
  COPY_STRUCTARRAY_TO_LIST(vector_list, vector_list_size, 
                           temp_vector_list, CubitVector);
  
  RefFace *temp_refface_ptr = (RefFace *) refface_ptr;
  
  void *result =  GTI->make_RefEdge(ref_edge_type, temp_ref_vertex_1, temp_ref_vertex_2, 
                                    temp_vector_list, temp_refface_ptr);

  DELETE_STRUCTLIST(temp_vector_list);
  
  return result;
}

/* RefEdge* */ void* GeometryTool_make_RefEdge_2(  /* RefVertex const* */ void *ref_vertex_1,
                                         /* RefVertex const* */ void *ref_vertex_2,
                                         /* RefFace* */ void *ref_face_ptr,
                                         /* RefVertex const* */ void *ref_vertex_3 ) 
/* < API function for GeometryTool::make_RefEdge */
{
  RefVertex *temp_ref_vertex_1 = (RefVertex *) ref_vertex_1;
  RefVertex *temp_ref_vertex_2 = (RefVertex *) ref_vertex_2;

  RefFace *temp_ref_face_ptr = (RefFace *) ref_face_ptr;
  
  RefVertex *temp_ref_vertex_3 = (RefVertex *) ref_vertex_3;
  
  return GTI->make_RefEdge(temp_ref_vertex_1, temp_ref_vertex_2, temp_ref_face_ptr, 
                           temp_ref_vertex_3);
}

/* RefEdge* */ void *GeometryTool_make_RefEdge_3( /* RefEdge * */ void *ref_edge) 
{
  RefEdge *temp_ref_edge = (RefEdge *) ref_edge;
  return GTI->make_RefEdge(temp_ref_edge);
}

/* RefEdge* */ void* GeometryTool_make_RefEdge_4(GeometryType ref_edge_type,
                                       /* RefVertex const* */ void *ref_vertex_1,
                                       /* RefVertex const* */ void *ref_vertex_2,
                                       /* CubitVector * */ CubitVectorStruct const* intermediate_point,
                                     CubitSense sense) 
/* < API function for GeometryTool::make_RefEdge */
{
  RefVertex *temp_ref_vertex_1 = (RefVertex *) ref_vertex_1;
  RefVertex *temp_ref_vertex_2 = (RefVertex *) ref_vertex_2;

  CubitVector temp_intermediate_point = *intermediate_point;

  return GTI->make_RefEdge(ref_edge_type, temp_ref_vertex_1, temp_ref_vertex_2, 
                           &temp_intermediate_point, sense);
}

/* RefFace* */ void* GeometryTool_make_RefFace_1(/* RefFace * */ void *from_ref_face,
                                     CubitBoolean extended_from)
/* < API function for GeometryTool::make_RefFace */
{
  RefFace *temp_from_ref_face = (RefFace *) from_ref_face;
  
  return GTI->make_RefFace(temp_from_ref_face, extended_from);
}

/* RefFace* */ void* GeometryTool_make_RefFace_2(GeometryType ref_face_type,
                                       /* DLRefEdgeList& */ void ***ref_edge_list,
                                     int *ref_edge_list_size,
                                       /* RefFace * */ void *ref_face_ptr) 
/* < API function for GeometryTool::make_RefFace */
{
  DLRefEdgeList temp_ref_edge_list;
  COPY_ARRAY_TO_LIST(*ref_edge_list, *ref_edge_list_size, temp_ref_edge_list);
  
  RefFace *temp_ref_face_ptr = (RefFace *) ref_face_ptr;
  
  return GTI->make_RefFace(ref_face_type, temp_ref_edge_list, temp_ref_face_ptr);
}
/* RefVolume* */ void* GeometryTool_make_RefVolume(GeometryType ref_volume_type,
                                           /* DLRefFaceList& */ void ***ref_face_list,
                                         int *ref_face_list_size) 
/* < API function for GeometryTool::make_RefVolume */
{
  DLRefFaceList temp_ref_face_list;
  COPY_ARRAY_TO_LIST(*ref_face_list, *ref_face_list_size, temp_ref_face_list);
  
  return GTI->make_RefVolume(ref_volume_type, temp_ref_face_list);
}

/* Body* */ void* GeometryTool_make_Body_3(/* DLRefVolumeList& */ void ***ref_volume_list, 
                               int *ref_volume_list_size) 
/* < API function for GeometryTool::make_Body */
{
  DLRefVolumeList temp_ref_volume_list;
  COPY_ARRAY_TO_LIST(*ref_volume_list, *ref_volume_list_size, temp_ref_volume_list);
  
  return GTI->make_Body(temp_ref_volume_list);
}

/* Body* */ void* GeometryTool_make_Body_4(/* RefFace * */ void *from_ref_face,
                               CubitBoolean extended_from)
/* < API function for GeometryTool::make_Body */
{
  RefFace *temp_from_ref_face = (RefFace *) from_ref_face;
    
  return GTI->make_Body(temp_from_ref_face, extended_from);
}

/* Body* */ void* GeometryTool_make_Body_5(GeometryType ref_face_type,
                                 /* DLRefEdgeList& */ void ***ref_edge_list,
                               int *ref_edge_list_size,
                                 /* RefFace * */ void *ref_face_ptr) 
/* < API function for GeometryTool::make_Body */
{
  DLRefEdgeList temp_ref_edge_list;
  COPY_ARRAY_TO_LIST(*ref_edge_list, *ref_edge_list_size, temp_ref_edge_list);
  
  RefFace *temp_ref_face_ptr = (RefFace *) ref_face_ptr;
  
  return GTI->make_Body(ref_face_type, temp_ref_edge_list, temp_ref_face_ptr);
}

/* RefFace* */ void *GeometryTool_make_free_RefFace(/* SurfaceSM * */ void *surfacesm_ptr)
{
  SurfaceSM *temp_surfacesm_ptr = (SurfaceSM *) surfacesm_ptr;
    
  return GTI->make_free_RefFace(temp_surfacesm_ptr );
  
}

/* RefEdge* */ void *GeometryTool_make_free_RefEdge(/* CurveSM * */ void *curvesm_ptr)
{
  CurveSM *temp_curvesm_ptr = (CurveSM *) curvesm_ptr;
    
  return GTI->make_free_RefEdge(temp_curvesm_ptr );
  
}

/* RefVertex* */ void *GeometryTool_make_free_RefVertex(/* PointSM * */ void *pointsm_ptr)
{
  PointSM *temp_pointsm_ptr = (PointSM *) pointsm_ptr;
    
  return GTI->make_free_RefVertex(temp_pointsm_ptr);
  
}

CubitStatus GeometryTool_sweep_translational(/* DLRefFaceList& */ void ***faces_to_be_swept_list,
                                   int *faces_to_be_swept_list_size,
                                   CubitVectorStruct sweep_vector,
                                   double draft_angle,
                                   int draft_type,
                                   int switchside )
/* < API function for GeometryTool::sweep_translational */
{
  DLRefFaceList temp_faces_to_be_swept_list;
  COPY_ARRAY_TO_LIST(*faces_to_be_swept_list, *faces_to_be_swept_list_size, temp_faces_to_be_swept_list);
  
  CubitVector temp_sweep_vector = sweep_vector;
  
  return GTI->sweep_translational(temp_faces_to_be_swept_list, 
                                  &temp_sweep_vector, 
                                  draft_angle, draft_type, switchside);
}

CubitStatus GeometryTool_sweep_rotational(/* DLRefFaceList& */ void ***faces_to_be_swept_list,
                                int *faces_to_be_swept_list_size,
                                  /* CubitVector * */ CubitVectorStruct *point,
                                  /* CubitVector * */ CubitVectorStruct *direction,
                                double angle,
                                int steps,
                                double draft_angle,
                                int draft_type,
                                int switchside )  
/* < API function for GeometryTool::sweep_rotational */
{
  DLRefFaceList temp_faces_to_be_swept_list;
  COPY_ARRAY_TO_LIST(*faces_to_be_swept_list, *faces_to_be_swept_list_size, temp_faces_to_be_swept_list);

  CubitVector *temp_point = NULL;
  if (point) temp_point = new CubitVector(*point);
  
  CubitVector *temp_direction = NULL;
  if (direction) temp_direction = new CubitVector(*direction);
  
  CubitStatus result =  GTI->sweep_rotational(temp_faces_to_be_swept_list, 
                                              temp_point, temp_direction, 
                                              angle, steps, draft_angle, draft_type, switchside);
  if (temp_point) delete temp_point;
  if (temp_direction) delete temp_direction;

  return result;
}

CubitStatus GeometryTool_sweep_along_curve(/* DLRefFaceList& */ void ***reffaces_to_be_swept_list,
                                 int *reffaces_to_be_swept_list_size,
                              /* DLRefEdgeList& */ void ***ref_edge_list,
                                 int *ref_edge_list_size,
                                 double draft_angle,
                                 int draft_type) 
/* < API function for GeometryTool::sweep_along_curve */
{
  DLRefFaceList temp_reffaces_to_be_swept_list;
  COPY_ARRAY_TO_LIST(*reffaces_to_be_swept_list, *reffaces_to_be_swept_list_size, temp_reffaces_to_be_swept_list);
  
  DLRefEdgeList temp_ref_edge_list;
  COPY_ARRAY_TO_LIST(*ref_edge_list, *ref_edge_list_size, temp_ref_edge_list);
  
  return GTI->sweep_along_curve(temp_reffaces_to_be_swept_list, 
                                temp_ref_edge_list, draft_angle, draft_type);
}

   
CubitStatus GeometryTool_create_body_from_surfs( /* DLRefFaceList & */ void ***ref_face_list,
                                       int *ref_face_list_size,
                                         /* Body *& */ void **new_body,
                                       int keep_old,
                                       int heal)
/* < API function for GeometryTool::create_body_from_surfs */
{
  DLRefFaceList temp_ref_face_list;
  COPY_ARRAY_TO_LIST(*ref_face_list, *ref_face_list_size, temp_ref_face_list);
  
  Body *temp_new_body;
  
  CubitStatus result = GTI->create_body_from_surfs(temp_ref_face_list, 
                                                   temp_new_body, keep_old, heal);

  *new_body = temp_new_body;

  return result;
}

CubitStatus GeometryTool_create_net_surface_1( /* DLRefFaceList& */ void ***ref_face_list,
                                   int *ref_face_list_size,
                                     /* Body *& */ void **new_body, 
                                     /* DLIList<DLCubitVectorList*> */ CubitVectorStruct **vec_lists_u, 
                                   int *vec_lists_u_sizes, int num_vec_lists_u,
                                     /* DLIList<DLCubitVectorList*> */ CubitVectorStruct **vec_lists_v, 
                                   int *vec_lists_v_sizes, int num_vec_lists_v,
                                   double net_tol,
                                   CubitBoolean heal )
/* < API function for GeometryTool::create_net_surface */
{
  DLRefFaceList temp_ref_face_list;
  COPY_ARRAY_TO_LIST(*ref_face_list, *ref_face_list_size, temp_ref_face_list);

  DLIList<DLCubitVectorList *> temp_vec_lists_u;
  int i;
  for (i = 0; i < num_vec_lists_u; i++) {
    DLCubitVectorList *temp_list = new DLCubitVectorList(vec_lists_u_sizes[i]);
    temp_vec_lists_u.append(temp_list);
    COPY_STRUCTARRAY_TO_LIST(vec_lists_u[i], vec_lists_u_sizes[i], 
                             *temp_list, CubitVector);
  }
  
  DLIList<DLCubitVectorList *> temp_vec_lists_v;
  for (i = 0; i < num_vec_lists_v; i++) {
    DLCubitVectorList *temp_list = new DLCubitVectorList(vec_lists_v_sizes[i]);
    temp_vec_lists_v.append(temp_list);
    COPY_STRUCTARRAY_TO_LIST(vec_lists_v[i], vec_lists_v_sizes[i], 
                             *temp_list, CubitVector);
  }
  
  Body *temp_new_body;
  
  CubitStatus result = GTI->create_net_surface(temp_ref_face_list, temp_new_body, 
                                               temp_vec_lists_u, temp_vec_lists_v, 
                                               net_tol, heal);

  *new_body = temp_new_body;

  for (i = 0; i < num_vec_lists_u; i++) {
    DLCubitVectorList *temp_list = temp_vec_lists_u.get_and_step();
    DELETE_STRUCTLIST(*temp_list);
  }
  
  for (i = 0; i < num_vec_lists_v; i++) {
    DLCubitVectorList *temp_list = temp_vec_lists_v.get_and_step();
    DELETE_STRUCTLIST(*temp_list);
  }
  
  return result;
}

CubitStatus GeometryTool_create_net_surface_2( /* DLRefEdgeList& */ void ***u_curves,
                                   int *u_curves_size,
                                     /* DLRefEdgeList& */ void ***v_curves,
                                   int *v_curves_size,
                                     /* Body *& */ void ** new_body,
                                   double net_tol, 
                                   CubitBoolean heal )
/* < API function for GeometryTool::create_net_surface */
{
  DLRefEdgeList temp_u_curves;
  COPY_ARRAY_TO_LIST(*u_curves, *u_curves_size, temp_u_curves);
  
  DLRefEdgeList temp_v_curves;
  COPY_ARRAY_TO_LIST(*v_curves, *v_curves_size, temp_v_curves);
  
  Body *temp_new_body;
  
  CubitStatus result = GTI->create_net_surface(temp_u_curves, temp_v_curves, 
                                               temp_new_body, net_tol, heal);

  *new_body = temp_new_body;

  return result;
}

CubitStatus GeometryTool_create_offset_surface( /* RefFace* */ void *ref_face_ptr, 
                                        /* Body*& */ void ** new_body, 
                                      double offset_distance )
/* < API function for GeometryTool::create_offset_surface */
{
  Body *temp_new_body;
  
  RefFace *temp_ref_face_ptr = (RefFace *) ref_face_ptr;
  
  CubitStatus result = GTI->create_offset_surface(temp_ref_face_ptr, temp_new_body, 
                                                  offset_distance);

  *new_body = temp_new_body;

  return result;
}

CubitStatus GeometryTool_create_offset_body( /* Body * */ void *body_ptr, 
                                     /* Body*& */ void ** new_body, 
                                   double offset_distance )
/* < API function for GeometryTool::create_offset_body */
{
  Body *temp_new_body;
  
  Body *temp_body_ptr = (Body *) body_ptr;
  
  CubitStatus result = GTI->create_offset_body(temp_body_ptr, 
                                               temp_new_body, offset_distance);

  *new_body = temp_new_body;

  return result;
}

CubitStatus GeometryTool_create_skin_surface( /* DLRefEdgeList& */ void ***curves,
                                    int *curves_size,
                                      /* Body*& */ void ** new_body )
/* < API function for GeometryTool::create_skin_surface */
{
  DLRefEdgeList temp_curves;
  COPY_ARRAY_TO_LIST(*curves, *curves_size, temp_curves);
  
  Body *temp_new_body;
  
  CubitStatus result = GTI->create_skin_surface( temp_curves, temp_new_body);

  *new_body = temp_new_body;

  return result;
}

enum CubitStatus GeometryTool_loft_surfaces( /* RefFace * */ void *face1, /* const double & */ double *takeoff1,
                                               /* RefFace * */ void *face2, /* const double & */ double *takeoff2,
                                               /* Body*& */ void **new_body,
                                             enum CubitBoolean arc_length_option,
                                             enum CubitBoolean twist_option,
                                             enum CubitBoolean align_direction,
                                             enum CubitBoolean perpendicular,
                                             enum CubitBoolean simplify_option) 
{
  RefFace *temp_face1 = (RefFace *) face1;
  RefFace *temp_face2 = (RefFace *) face2;

  Body *temp_new_body = NULL;
  
  CubitStatus result = GTI->loft_surfaces(temp_face1, *takeoff1, temp_face2, *takeoff2,
                                          temp_new_body,
                                          arc_length_option, twist_option, align_direction, 
                                          perpendicular, simplify_option);
  
  *new_body = temp_new_body;
  
  return result;
}

CubitStatus GeometryTool_create_arc_three_1( /* RefVertex * */ void *ref_vertex1, 
                                   /* RefVertex * */ void *ref_vertex2,
                                   /* RefVertex * */ void *ref_vertex3, 
                                 CubitBoolean full )
/* < API function for GeometryTool::create_arc_three */
{
  RefVertex *temp_ref_vertex1 = (RefVertex *) ref_vertex1;
  RefVertex *temp_ref_vertex2 = (RefVertex *) ref_vertex2;
  RefVertex *temp_ref_vertex3 = (RefVertex *) ref_vertex3;
    
  return GTI->create_arc_three(temp_ref_vertex1, temp_ref_vertex2, 
                               temp_ref_vertex3, full);
}

CubitStatus GeometryTool_create_arc_three_2( /* RefEdge * */ void *ref_edge1, 
                                   /* RefEdge * */ void *ref_edge2,
                              /* RefEdge * */ void *ref_edge3, 
                                 CubitBoolean full )
/* < API function for GeometryTool::create_arc_three */
{
  RefEdge *temp_ref_edge1 = (RefEdge *) ref_edge1;
  RefEdge *temp_ref_edge2 = (RefEdge *) ref_edge2;
  RefEdge *temp_ref_edge3 = (RefEdge *) ref_edge3;
    
  return GTI->create_arc_three(temp_ref_edge1, temp_ref_edge2, 
                               temp_ref_edge3, full);
}

CubitStatus GeometryTool_create_arc_center_edge( /* RefVertex* */ void *ref_vertex1, 
                                         /* RefVertex* */ void *ref_vertex2,
                                         /* RefVertex* */ void *ref_vertex3, 
                                       double radius,
                                       CubitBoolean full )
/* < API function for GeometryTool::create_arc_center_edge */
{
  RefVertex *temp_ref_vertex1 = (RefVertex *) ref_vertex1;
  RefVertex *temp_ref_vertex2 = (RefVertex *) ref_vertex2;
  RefVertex *temp_ref_vertex3 = (RefVertex *) ref_vertex3;
    
  return GTI->create_arc_center_edge(temp_ref_vertex1, temp_ref_vertex2, 
                                     temp_ref_vertex3, radius, full);
}


  /* /<HR><H3>Topology and geometry deletion</H3> */
void GeometryTool_delete_Body_1(/* DLBodyList& */ void ***body_list,
                    int *body_list_size,
                    CubitBoolean remove_solid_model_entities)
/* < API function for GeometryTool::delete_Body */
{
  DLBodyList temp_body_list;
  COPY_ARRAY_TO_LIST(*body_list, *body_list_size, temp_body_list);
  
  GTI->delete_Body(temp_body_list,
                   remove_solid_model_entities);
}

CubitStatus GeometryTool_delete_Body_2(
    /* Body * */ void *body_ptr,
    CubitBoolean remove_solid_model_entities)
/* < API function for GeometryTool::delete_Body */
{
  Body *temp_body_ptr = (Body *) body_ptr;
    
  return GTI->delete_Body(temp_body_ptr, remove_solid_model_entities);
}

CubitStatus GeometryTool_delete_RefAssembly(
    /* RefAssembly *& */ void *ref_assembly_ptr,
    CubitBoolean remove_solid_model_entities)
/* < API function for GeometryTool::delete_RefAssembly */
{
  RefAssembly *temp_ref_assembly_ptr = (RefAssembly *) ref_assembly_ptr;
    
  return GTI->delete_RefAssembly(temp_ref_assembly_ptr, remove_solid_model_entities);
}

CubitStatus GeometryTool_delete_RefPart(
    /* RefPart *& */ void *ref_part_ptr,
    CubitBoolean remove_solid_model_entities)
/* < API function for GeometryTool::delete_RefPart */
{
  RefPart *temp_ref_part_ptr = (RefPart *) ref_part_ptr;
    
  return GTI->delete_RefPart(temp_ref_part_ptr, remove_solid_model_entities);
}

CubitStatus GeometryTool_delete_RefEntity(
    /* RefEntity *& */ void *ref_entity_ptr,
    CubitBoolean remove_solid_model_entities,
    CubitBoolean remove_lower_entities )
/* < API function for GeometryTool::delete_RefEntity */
{
  RefEntity *temp_ref_entity_ptr = (RefEntity *) ref_entity_ptr;
    
  return GTI->delete_RefEntity(temp_ref_entity_ptr, 
                               remove_solid_model_entities, remove_lower_entities);
}

void GeometryTool_cleanout_deactivated_geometry()
/* < API function for GeometryTool::cleanout_deactivated_geometry */
{
  GTI->cleanout_deactivated_geometry();
}

void GeometryTool_cleanout_temporary_geometry ()
/* < API function for GeometryTool::cleanout_temporary_geometry  */
{
  GTI->cleanout_temporary_geometry();
}

void GeometryTool_delete_geometry()
/* < API function for GeometryTool::delete_geometry */
{
  GTI->delete_geometry();
}

  /* /< <HR><H3>Miscellaneous geometry evaluation functions</H3> */
CubitStatus GeometryTool_interpolate_along_surface( CubitVectorStruct *vector_1,
                                          CubitVectorStruct *vector_2,
                                            /* DLCubitVectorList & */ CubitVectorStruct *vector_list,
                                          int vector_list_size,
                                            /* RefFace * */ void *ref_face_ptr,
                                          int number_points )
/* < API function for GeometryTool::interpolate_along_surface */
{
  CubitVector temp_vector_1 = *vector_1;
  CubitVector temp_vector_2 = *vector_2;
  
  DLCubitVectorList temp_vector_list;
  COPY_STRUCTARRAY_TO_LIST(vector_list, vector_list_size, 
                           temp_vector_list, CubitVector);
  
  RefFace *temp_ref_face_ptr = (RefFace *) ref_face_ptr;
  
  CubitStatus result = GTI->interpolate_along_surface(&temp_vector_1, &temp_vector_2, 
                                                      temp_vector_list, temp_ref_face_ptr, 
                                                      number_points);

  DELETE_STRUCTLIST(temp_vector_list);
  
  return result;
}

CubitBoolean GeometryTool_about_spatially_equal_1 (const /* CubitVector& */ CubitVectorStruct* Vec1, 
                                       const /* CubitVector& */ CubitVectorStruct* Vec2,
                                       double tolerance_factor)
/* < API function for GeometryTool::about_spatially_equal  */
{
  CubitVector temp_Vec1 = *Vec1;
  CubitVector temp_Vec2 = *Vec2;

  return GTI->about_spatially_equal (temp_Vec1, temp_Vec2, tolerance_factor);
}

CubitBoolean GeometryTool_about_spatially_equal_2 (/* RefVertex * */ void *refVertex1, 
                                         /* RefVertex * */ void *refVertex2, 
                                       double tolerance_factor)
/* < API function for GeometryTool::about_spatially_equal  */
{
  RefVertex *temp_refVertex1 = (RefVertex *) refVertex1;
  RefVertex *temp_refVertex2 = (RefVertex *) refVertex2;
  
  return GTI->about_spatially_equal (temp_refVertex1, temp_refVertex2, 
                                     tolerance_factor);
}

double GeometryTool_geometric_angle(/* RefEdge * */ void *ref_edge_1, 
                            /* RefEdge * */ void *ref_edge_2,
                            /* RefFace * */ void *ref_face)
/* < API function for GeometryTool::geometric_angle */
{
  RefEdge *temp_ref_edge_1 = (RefEdge *) ref_edge_1;
  RefEdge *temp_ref_edge_2 = (RefEdge *) ref_edge_2;
    
  RefFace *temp_ref_face = (RefFace *) ref_face;
  
  return GTI->geometric_angle(temp_ref_edge_1, temp_ref_edge_2, temp_ref_face);
}

  double GeometryTool_geometric_angle(/* CoEdge* */ void *co_edge_1, 
                                        /* CoEdge* */ void *co_edge_2 ) 
{

  CoEdge *temp_co_edge_1 = (CoEdge *) co_edge_1;
  CoEdge *temp_co_edge_2 = (CoEdge *) co_edge_2;

  return GTI->geometric_angle(temp_co_edge_1, temp_co_edge_2);
}

double GeometryTool_surface_angle( /* RefEdge * */ void *ref_edge,
                                     /* RefFace * */ void *ref_face_1, 
                                     /* RefFace * */ void *ref_face_2,
                                     /* RefVolume * */ void *ref_volume )
{
  RefEdge *temp_ref_edge = (RefEdge *) ref_edge;
  RefFace *temp_ref_face_1 = (RefFace *) ref_face_1;
  RefFace *temp_ref_face_2 = (RefFace *) ref_face_2;
  RefVolume *temp_ref_volume = (RefVolume *) ref_volume;

  return GTI->surface_angle(temp_ref_edge, temp_ref_face_1, temp_ref_face_2, 
                            temp_ref_volume);
}

CubitStatus GeometryTool_get_intersections( /* RefEdge * */ void *ref_edge1, 
                                    /* RefEdge * */ void *ref_edge2,
                                    /* DLCubitVectorList& */ CubitVectorStruct **intersection_list,
                                  int *intersection_list_size,
                                  CubitBoolean bounded,
                                  CubitBoolean closest )
/* < API function for GeometryTool::get_intersections */
{
  RefEdge *temp_ref_edge1 = (RefEdge *) ref_edge1;
  RefEdge *temp_ref_edge2 = (RefEdge *) ref_edge2;

  DLCubitVectorList temp_intersection_list;
  CubitStatus result = GTI->get_intersections(temp_ref_edge1, temp_ref_edge2, 
                                              temp_intersection_list, 
                                              bounded, closest);

  COPY_LIST_TO_STRUCTARRAY(temp_intersection_list, *intersection_list, 
                           *intersection_list_size, CubitVectorStruct);
  
  DELETE_STRUCTLIST(temp_intersection_list);
  
  return result;
}

CubitBoxStruct GeometryTool_bounding_box_of_bodies()
/* < API function for GeometryTool::bounding_box_of_bodies */
{
  return GTI->bounding_box_of_bodies();
}


  /* /< <HR><H3>Healing functions</H3> */
CubitStatus GeometryTool_autoheal_bodies(/* DLBodyList & */ void ***body_list,
                               int *body_list_size,
                                 /* DLBodyList & */ void ***new_body_list,
                               int *new_body_list_size,
                               CubitBoolean rebuild, 
                               CubitBoolean keep_old, CubitBoolean make_tolerant,
                               FILE* logfile_ptr )
/* < API function for GeometryTool::autoheal_bodies */
{
  DLBodyList temp_body_list;
  COPY_ARRAY_TO_LIST(*body_list, *body_list_size, temp_body_list);
  
  DLBodyList temp_new_body_list;
  
  CubitStatus result = GTI->autoheal_bodies(temp_body_list, temp_new_body_list, 
                                            rebuild, keep_old, make_tolerant, 
                                            logfile_ptr);

  COPY_LIST_TO_ARRAY(temp_new_body_list, *new_body_list, *new_body_list_size);

  return result;
}

CubitStatus GeometryTool_healer_analyze_badgeom( /* DLBodyList & */ void ***body_list,
                                       int *body_list_size,
                                       FILE* logfile )
/* < API function for GeometryTool::healer_analyze_badgeom */
{
  DLBodyList temp_body_list;
  COPY_ARRAY_TO_LIST(*body_list, *body_list_size, temp_body_list);
  
  return GTI->healer_analyze_badgeom(temp_body_list, logfile);
}

CubitStatus GeometryTool_healer_show_badgeom( /* DLBodyList & */ void ***body_list,
                                    int *body_list_size)
/* < API function for GeometryTool::healer_show_badgeom */
{
  DLBodyList temp_body_list;
  COPY_ARRAY_TO_LIST(*body_list, *body_list_size, temp_body_list);
  
  return GTI->healer_show_badgeom(temp_body_list);
}

CubitStatus GeometryTool_healer_show_tcurves( /* DLBodyList & */ void ***body_list, int *body_list_size)
/* < API function for GeometryTool::healer_show_tcurves */
{
  DLBodyList temp_body_list;
  COPY_ARRAY_TO_LIST(*body_list, *body_list_size, temp_body_list);
  
  return GTI->healer_show_tcurves(temp_body_list);
}

CubitStatus GeometryTool_heal_incremental( /* DLBodyList & */ void ***body_list,
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
                                 CubitBoolean keep_old,
                                 CubitBoolean make_tolerant,
                                 FILE* logfile_ptr)
/* < API function for GeometryTool::heal_incremental */
{
  DLBodyList temp_body_list;
  COPY_ARRAY_TO_LIST(*body_list, *body_list_size, temp_body_list);
  
  DLBodyList temp_new_bodies;
  
  CubitStatus result = GTI->heal_incremental(temp_body_list, temp_new_bodies, simplify_tol, 
                                             stitch_min_tol, stitch_max_tol, geombuild_tol, 
                                             analytic_tol, isospline_tol, reblend_classify_tol, 
                                             reblend_tol, keep_old, make_tolerant, logfile_ptr);

  COPY_LIST_TO_ARRAY(temp_new_bodies, *new_bodies, *new_bodies_size);
  

  return result;
}

void GeometryTool_healer_list_incremental()
/* < API function for GeometryTool::healer_list_incremental */
{
  GTI->healer_list_incremental();
}

void GeometryTool_healer_list_tolerances( /* DLBodyList & */ void ***body_list,
                                                   int *body_list_size)
/* < API function for GeometryTool::healer_list_tolerances */
{
  DLBodyList temp_body_list;
  COPY_ARRAY_TO_LIST(*body_list, *body_list_size, temp_body_list);
  
  GTI->healer_list_tolerances(temp_body_list);
}

double GeometryTool_healer_get_default_simplify_tol()
/* < API function for GeometryTool::healer_get_default_simplify_tol */
{
  return GTI->healer_get_default_simplify_tol();
}

void GeometryTool_healer_set_default_simplify_tol( double tol)
/* < API function for GeometryTool::healer_set_default_simplify_tol */
{
  GTI->healer_set_default_simplify_tol(tol);
}

double GeometryTool_healer_get_default_stitch_min_tol()
/* < API function for GeometryTool::healer_get_default_stitch_min_tol */
{
  return GTI->healer_get_default_stitch_min_tol();
}

void GeometryTool_healer_set_default_stitch_min_tol( double tol)
/* < API function for GeometryTool::healer_set_default_stitch_min_tol */
{
  GTI->healer_set_default_stitch_min_tol(tol);
}

double GeometryTool_healer_get_default_stitch_max_tol()
/* < API function for GeometryTool::healer_get_default_stitch_max_tol */
{
  return GTI->healer_get_default_stitch_max_tol();
}

void GeometryTool_healer_set_default_stitch_max_tol( double tol)
/* < API function for GeometryTool::healer_set_default_stitch_max_tol */
{
  GTI->healer_set_default_stitch_max_tol(tol);
}

double GeometryTool_healer_get_default_geombuild_tol()
/* < API function for GeometryTool::healer_get_default_geombuild_tol */
{
  return GTI->healer_get_default_geombuild_tol();
}

void GeometryTool_healer_set_default_geombuild_tol( double tol)
/* < API function for GeometryTool::healer_set_default_geombuild_tol */
{
  GTI->healer_set_default_geombuild_tol(tol);
}

double GeometryTool_healer_get_default_analytic_tol()
/* < API function for GeometryTool::healer_get_default_analytic_tol */
{
  return GTI->healer_get_default_analytic_tol();
}

void GeometryTool_healer_set_default_analytic_tol( double tol)
/* < API function for GeometryTool::healer_set_default_analytic_tol */
{
  GTI->healer_set_default_analytic_tol(tol);
}

double GeometryTool_healer_get_default_isospline_tol()
/* < API function for GeometryTool::healer_get_default_isospline_tol */
{
  return GTI->healer_get_default_isospline_tol();
}

void GeometryTool_healer_set_default_isospline_tol( double tol)
/* < API function for GeometryTool::healer_set_default_isospline_tol */
{
  GTI->healer_set_default_isospline_tol(tol);
}

double GeometryTool_healer_get_default_reblend_classify_tol()
/* < API function for GeometryTool::healer_get_default_reblend_classify_tol */
{
  return GTI->healer_get_default_reblend_classify_tol();
}

void GeometryTool_healer_set_default_reblend_classify_tol( double tol)
/* < API function for GeometryTool::healer_set_default_reblend_classify_tol */
{
  GTI->healer_set_default_reblend_classify_tol(tol);
}

double GeometryTool_healer_get_default_reblend_tol()
/* < API function for GeometryTool::healer_get_default_reblend_tol */
{
  return GTI->healer_get_default_reblend_tol();
}

void GeometryTool_healer_set_default_reblend_tol( double tol)
/* < API function for GeometryTool::healer_set_default_reblend_tol */
{
  GTI->healer_set_default_reblend_tol(tol);
}

void GeometryTool_healer_reset_default_tolerances()
/* < API function for GeometryTool::healer_reset_default_tolerances */
{
  GTI->healer_reset_default_tolerances();
}

void GeometryTool_healer_list_default_tolerances()
/* < API function for GeometryTool::healer_list_default_tolerances */
{
  GTI->healer_list_default_tolerances();
}

void GeometryTool_healer_clean_attributes( /* DLBodyList& */ void ***body_list,
                                                   int *body_list_size)
/* < API function for GeometryTool::healer_clean_attributes */
{
  DLBodyList temp_body_list;
  COPY_ARRAY_TO_LIST(*body_list, *body_list_size, temp_body_list);

  GTI->healer_clean_attributes(temp_body_list);
}

CubitBoolean GeometryTool_healer_get_cleanatt_flg()
/* < API function for GeometryTool::healer_get_cleanatt_flg */
{
  return GTI->healer_get_cleanatt_flg();
}

void GeometryTool_healer_set_cleanatt_flg( CubitBoolean flg)
/* < API function for GeometryTool::healer_set_cleanatt_flg */
{
  GTI->healer_set_cleanatt_flg(flg);
}

int GeometryTool_healer_get_show_method() 
/* < API function for GeometryTool::healer_get_show_method */
{
  return GTI->healer_get_show_method();
}

void GeometryTool_healer_set_show_method( int method)
/* < API function for GeometryTool::healer_set_show_method */
{
  GTI->healer_set_show_method(method);
}
   
CubitBoolean GeometryTool_healer_get_show_summary_flg()
/* < API function for GeometryTool::healer_get_show_summary_flg */
{
  return GTI->healer_get_show_summary_flg();
}

void GeometryTool_healer_set_show_summary_flg( CubitBoolean flg)
/* < API function for GeometryTool::healer_set_show_summary_flg */
{
  GTI->healer_set_show_summary_flg(flg);
}
   
CubitBoolean GeometryTool_healer_get_show_details_flg()
/* < API function for GeometryTool::healer_get_show_details_flg */
{
  return GTI->healer_get_show_details_flg();
}

void GeometryTool_healer_set_show_details_flg( CubitBoolean flg)
/* < API function for GeometryTool::healer_set_show_details_flg */
{
  GTI->healer_set_show_details_flg(flg);
}
   
CubitBoolean GeometryTool_healer_get_show_on_heal_flg()
/* < API function for GeometryTool::healer_get_show_on_heal_flg */
{
  return GTI->healer_get_show_on_heal_flg();
}

void GeometryTool_healer_set_show_on_heal_flg( CubitBoolean flg)
/* < API function for GeometryTool::healer_set_show_on_heal_flg */
{
  GTI->healer_set_show_on_heal_flg(flg);
}

CubitBoolean GeometryTool_healer_get_check_vol_on_heal_flg()
/* < API function for GeometryTool::healer_get_check_vol_on_heal_flg */
{
  return GTI->healer_get_check_vol_on_heal_flg();
}

void GeometryTool_healer_set_check_vol_on_heal_flg( CubitBoolean flg)
/* < API function for GeometryTool::healer_set_check_vol_on_heal_flg */
{
  GTI->healer_set_check_vol_on_heal_flg(flg);
}

double GeometryTool_healer_get_vol_on_heal_limit()
/* < API function for GeometryTool::healer_get_vol_on_heal_limit */
{
  return GTI->healer_get_vol_on_heal_limit();
}

void GeometryTool_healer_set_vol_on_heal_limit( double limit)
/* < API function for GeometryTool::healer_set_vol_on_heal_limit */
{
  GTI->healer_set_vol_on_heal_limit(limit);
}
   
CubitBoolean GeometryTool_healer_get_check_surf_on_heal_flg()
/* < API function for GeometryTool::healer_get_check_surf_on_heal_flg */
{
  return GTI->healer_get_check_surf_on_heal_flg();
}

void GeometryTool_healer_set_check_surf_on_heal_flg( CubitBoolean flg)
/* < API function for GeometryTool::healer_set_check_surf_on_heal_flg */
{
  GTI->healer_set_check_surf_on_heal_flg(flg);
}

double GeometryTool_healer_get_surf_on_heal_limit()
/* < API function for GeometryTool::healer_get_surf_on_heal_limit */
{
  return GTI->healer_get_surf_on_heal_limit();
}

void GeometryTool_healer_set_surf_on_heal_limit( double limit)
/* < API function for GeometryTool::healer_set_surf_on_heal_limit */
{
  GTI->healer_set_surf_on_heal_limit(limit);
}
   
CubitBoolean GeometryTool_healer_get_check_curve_on_heal_flg()
/* < API function for GeometryTool::healer_get_check_curve_on_heal_flg */
{
  return GTI->healer_get_check_curve_on_heal_flg();
}

void GeometryTool_healer_set_check_curve_on_heal_flg( CubitBoolean flg)
/* < API function for GeometryTool::healer_set_check_curve_on_heal_flg */
{
  GTI->healer_set_check_curve_on_heal_flg(flg);
}

double GeometryTool_healer_get_curve_on_heal_limit()
/* < API function for GeometryTool::healer_get_curve_on_heal_limit */
{
  return GTI->healer_get_curve_on_heal_limit();
}

void GeometryTool_healer_set_curve_on_heal_limit( double limit)
/* < API function for GeometryTool::healer_set_curve_on_heal_limit */
{
  GTI->healer_set_curve_on_heal_limit(limit);
}
   
CubitBoolean GeometryTool_healer_get_show_bad_vertices_flg()
/* < API function for GeometryTool::healer_get_show_bad_vertices_flg */
{
  return GTI->healer_get_show_bad_vertices_flg();
}

void GeometryTool_healer_set_show_bad_vertices_flg( CubitBoolean flg)
/* < API function for GeometryTool::healer_set_show_bad_vertices_flg */
{
  GTI->healer_set_show_bad_vertices_flg(flg);
}

CubitBoolean GeometryTool_healer_get_show_bad_curves_flg()
/* < API function for GeometryTool::healer_get_show_bad_curves_flg */
{
  return GTI->healer_get_show_bad_curves_flg();
}

void GeometryTool_healer_set_show_bad_curves_flg( CubitBoolean flg)
/* < API function for GeometryTool::healer_set_show_bad_curves_flg */
{
  GTI->healer_set_show_bad_curves_flg(flg);
}

CubitBoolean GeometryTool_healer_get_show_bad_coedges_flg()
/* < API function for GeometryTool::healer_get_show_bad_coedges_flg */
{
  return GTI->healer_get_show_bad_coedges_flg();
}

void GeometryTool_healer_set_show_bad_coedges_flg( CubitBoolean flg)
/* < API function for GeometryTool::healer_set_show_bad_coedges_flg */
{
  GTI->healer_set_show_bad_coedges_flg(flg);
}

CubitBoolean GeometryTool_healer_get_show_bad_loops_flg()
/* < API function for GeometryTool::healer_get_show_bad_loops_flg */
{
  return GTI->healer_get_show_bad_loops_flg();
}

void GeometryTool_healer_set_show_bad_loops_flg( CubitBoolean flg)
/* < API function for GeometryTool::healer_set_show_bad_loops_flg */
{
  GTI->healer_set_show_bad_loops_flg(flg);
}

CubitBoolean GeometryTool_healer_get_show_bad_surfaces_flg()
/* < API function for GeometryTool::healer_get_show_bad_surfaces_flg */
{
  return GTI->healer_get_show_bad_surfaces_flg();
}

void GeometryTool_healer_set_show_bad_surfaces_flg( CubitBoolean flg)
/* < API function for GeometryTool::healer_set_show_bad_surfaces_flg */
{
  GTI->healer_set_show_bad_surfaces_flg(flg);
}

CubitBoolean GeometryTool_healer_get_show_bad_shells_flg()
/* < API function for GeometryTool::healer_get_show_bad_shells_flg */
{
  return GTI->healer_get_show_bad_shells_flg();
}

void GeometryTool_healer_set_show_bad_shells_flg( CubitBoolean flg)
/* < API function for GeometryTool::healer_set_show_bad_shells_flg */
{
  GTI->healer_set_show_bad_shells_flg(flg);
}

CubitBoolean GeometryTool_healer_get_show_bad_volumes_flg()
/* < API function for GeometryTool::healer_get_show_bad_volumes_flg */
{
  return GTI->healer_get_show_bad_volumes_flg();
}

void GeometryTool_healer_set_show_bad_volumes_flg( CubitBoolean flg)
/* < API function for GeometryTool::healer_set_show_bad_volumes_flg */
{
  GTI->healer_set_show_bad_volumes_flg(flg);
}

CubitBoolean GeometryTool_healer_get_show_bad_bodies_flg()
/* < API function for GeometryTool::healer_get_show_bad_bodies_flg */
{
  return GTI->healer_get_show_bad_bodies_flg();
}

void GeometryTool_healer_set_show_bad_bodies_flg( CubitBoolean flg)
/* < API function for GeometryTool::healer_set_show_bad_bodies_flg */
{
  GTI->healer_set_show_bad_bodies_flg(flg);
}

void GeometryTool_healer_list_onshow_flgs()
/* < API function for GeometryTool::healer_list_onshow_flgs */
{
  GTI->healer_list_onshow_flgs();
}

CubitBoolean GeometryTool_healer_get_inc_preprocess_flg()
/* < API function for GeometryTool::healer_get_inc_preprocess_flg */
{
  return GTI->healer_get_inc_preprocess_flg();
}

void GeometryTool_healer_set_inc_preprocess_flg( CubitBoolean flg)
/* < API function for GeometryTool::healer_set_inc_preprocess_flg */
{
  GTI->healer_set_inc_preprocess_flg(flg);
}

CubitBoolean GeometryTool_healer_get_inc_simplify_flg()
/* < API function for GeometryTool::healer_get_inc_simplify_flg */
{
  return GTI->healer_get_inc_simplify_flg();
}

void GeometryTool_healer_set_inc_simplify_flg( CubitBoolean flg)
/* < API function for GeometryTool::healer_set_inc_simplify_flg */
{
  GTI->healer_set_inc_simplify_flg(flg);
}

CubitBoolean GeometryTool_healer_get_inc_stitch_flg()
/* < API function for GeometryTool::healer_get_inc_stitch_flg */
{
  return GTI->healer_get_inc_stitch_flg();
}

void GeometryTool_healer_set_inc_stitch_flg( CubitBoolean flg)
/* < API function for GeometryTool::healer_set_inc_stitch_flg */
{
  GTI->healer_set_inc_stitch_flg(flg);
}

CubitBoolean GeometryTool_healer_get_inc_geombuild_flg()
/* < API function for GeometryTool::healer_get_inc_geombuild_flg */
{
  return GTI->healer_get_inc_geombuild_flg();
}

void GeometryTool_healer_set_inc_geombuild_flg( CubitBoolean flg)
/* < API function for GeometryTool::healer_set_inc_geombuild_flg */
{
  GTI->healer_set_inc_geombuild_flg(flg);
}

CubitBoolean GeometryTool_healer_get_inc_analytic_flg()
/* < API function for GeometryTool::healer_get_inc_analytic_flg */
{
  return GTI->healer_get_inc_analytic_flg();
}

void GeometryTool_healer_set_inc_analytic_flg( CubitBoolean flg)
/* < API function for GeometryTool::healer_set_inc_analytic_flg */
{
  GTI->healer_set_inc_analytic_flg(flg);
}

CubitBoolean GeometryTool_healer_get_inc_isospline_flg()
/* < API function for GeometryTool::healer_get_inc_isospline_flg */
{
  return GTI->healer_get_inc_isospline_flg();
}

void GeometryTool_healer_set_inc_isospline_flg( CubitBoolean flg)
/* < API function for GeometryTool::healer_set_inc_isospline_flg */
{
  GTI->healer_set_inc_isospline_flg(flg);
}

CubitBoolean GeometryTool_healer_get_inc_reblend_flg()
/* < API function for GeometryTool::healer_get_inc_reblend_flg */
{
  return GTI->healer_get_inc_reblend_flg();
}

void GeometryTool_healer_set_inc_reblend_flg( CubitBoolean flg)
/* < API function for GeometryTool::healer_set_inc_reblend_flg */
{
  GTI->healer_set_inc_reblend_flg(flg);
}

CubitBoolean GeometryTool_healer_get_inc_sharpedge_flg()
/* < API function for GeometryTool::healer_get_inc_sharpedge_flg */
{
  return GTI->healer_get_inc_sharpedge_flg();
}

void GeometryTool_healer_set_inc_sharpedge_flg( CubitBoolean flg)
/* < API function for GeometryTool::healer_set_inc_sharpedge_flg */
{
  GTI->healer_set_inc_sharpedge_flg(flg);
}

CubitBoolean GeometryTool_healer_get_inc_genericspline_flg()
/* < API function for GeometryTool::healer_get_inc_genericspline_flg */
{
  return GTI->healer_get_inc_genericspline_flg();
}

void GeometryTool_healer_set_inc_genericspline_flg( CubitBoolean flg)
/* < API function for GeometryTool::healer_set_inc_genericspline_flg */
{
  GTI->healer_set_inc_genericspline_flg(flg);
}

CubitBoolean GeometryTool_healer_get_inc_wrapup_flg()
/* < API function for GeometryTool::healer_get_inc_wrapup_flg */
{
  return GTI->healer_get_inc_wrapup_flg();
}

void GeometryTool_healer_set_inc_wrapup_flg( CubitBoolean flg)
/* < API function for GeometryTool::healer_set_inc_wrapup_flg */
{
  GTI->healer_set_inc_wrapup_flg(flg);
}

CubitBoolean GeometryTool_healer_get_inc_postprocess_flg()
/* < API function for GeometryTool::healer_get_inc_postprocess_flg */
{
  return GTI->healer_get_inc_postprocess_flg();
}

void GeometryTool_healer_set_inc_postprocess_flg( CubitBoolean flg)
/* < API function for GeometryTool::healer_set_inc_postprocess_flg */
{
  GTI->healer_set_inc_postprocess_flg(flg);
}

CubitStatus GeometryTool_force_simplify_to_plane( /* DLRefFaceList & */ void ***ref_face_list,
                                        int *ref_face_list_size,
                                          /* DLBodyList & */ void ***new_body_list,
                                        int *new_body_list_size, 
                                        CubitBoolean keep_old_body)
/* < API function for GeometryTool::force_simplify_to_plane */
{
  DLRefFaceList temp_ref_face_list;
  COPY_ARRAY_TO_LIST(*ref_face_list, *ref_face_list_size, temp_ref_face_list);
  
  DLBodyList temp_new_body_list;
  
  CubitStatus result = GTI->force_simplify_to_plane(temp_ref_face_list, temp_new_body_list, keep_old_body);

  COPY_LIST_TO_ARRAY(temp_new_body_list, *new_body_list, *new_body_list_size);
  
  return result;
}

CubitStatus GeometryTool_force_simplify_to_cylinder( /* DLRefFaceList & */ void ***ref_face_list,
                                           int *ref_face_list_size,
                                             /* DLBodyList & */ void ***new_body_list,
                                           int *new_body_list_size, 
                                           CubitBoolean keep_old_body)
/* < API function for GeometryTool::force_simplify_to_cylinder */
{
  DLRefFaceList temp_ref_face_list;
  COPY_ARRAY_TO_LIST(*ref_face_list, *ref_face_list_size, temp_ref_face_list);
  
  DLBodyList temp_new_body_list;
  
  CubitStatus result = GTI->force_simplify_to_cylinder(temp_ref_face_list, temp_new_body_list, keep_old_body);

  COPY_LIST_TO_ARRAY(temp_new_body_list, *new_body_list, *new_body_list_size);

  return result;
}

CubitStatus GeometryTool_force_simplify_to_cone( /* DLRefFaceList & */ void ***ref_face_list,
                                       int *ref_face_list_size,
                                         /* DLBodyList & */ void ***new_body_list,
                                       int *new_body_list_size, 
                                       CubitBoolean keep_old_body)
/* < API function for GeometryTool::force_simplify_to_cone */
{
  DLRefFaceList temp_ref_face_list;
  COPY_ARRAY_TO_LIST(*ref_face_list, *ref_face_list_size, temp_ref_face_list);
  
  DLBodyList temp_new_body_list;
  
  CubitStatus result = GTI->force_simplify_to_cone(temp_ref_face_list, temp_new_body_list, keep_old_body);

  COPY_LIST_TO_ARRAY(temp_new_body_list, *new_body_list, *new_body_list_size);
  
  return result;
}

CubitStatus GeometryTool_force_simplify_to_sphere( /* DLRefFaceList & */ void ***ref_face_list,
                                         int *ref_face_list_size,
                                           /* DLBodyList & */ void ***new_body_list,
                                         int *new_body_list_size, 
                                         CubitBoolean keep_old_body)
/* < API function for GeometryTool::force_simplify_to_sphere */
{
  DLRefFaceList temp_ref_face_list;
  COPY_ARRAY_TO_LIST(*ref_face_list, *ref_face_list_size, temp_ref_face_list);
  
  DLBodyList temp_new_body_list;
  
  CubitStatus result = GTI->force_simplify_to_sphere(temp_ref_face_list, temp_new_body_list, keep_old_body);

  COPY_LIST_TO_ARRAY(temp_new_body_list, *new_body_list, *new_body_list_size);
  
  return result;
}

CubitStatus GeometryTool_force_simplify_to_torus( /* DLRefFaceList & */ void ***ref_face_list,
                                        int *ref_face_list_size,
                                          /* DLBodyList & */ void ***new_body_list,
                                        int *new_body_list_size, 
                                        CubitBoolean keep_old_body)
/* < API function for GeometryTool::force_simplify_to_torus */
{
  DLRefFaceList temp_ref_face_list;
  COPY_ARRAY_TO_LIST(*ref_face_list, *ref_face_list_size, temp_ref_face_list);
  
  DLBodyList temp_new_body_list;
  
  CubitStatus result = GTI->force_simplify_to_torus(temp_ref_face_list, temp_new_body_list, 
                                                    keep_old_body);

  COPY_LIST_TO_ARRAY(temp_new_body_list, *new_body_list, *new_body_list_size);
  
  return result;
}

  /* /< <HR><H3>Merging functions.</H3> */
   
void GeometryTool_set_geometry_factor( double fac)
/* < API function for GeometryTool::set_geometry_factor */
{
  GTI->set_geometry_factor(fac);
}

double GeometryTool_get_geometry_factor()
/* < API function for GeometryTool::get_geometry_factor */
{
  return GTI->get_geometry_factor();
}

void GeometryTool_set_merge_test_bbox(CubitBoolean tof)
/* < API function for GeometryTool::set_merge_test_bbox */
{
  GTI->set_merge_test_bbox(tof);
}

CubitBoolean GeometryTool_get_merge_test_bbox()
/* < API function for GeometryTool::get_merge_test_bbox */
{
  return GTI->get_merge_test_bbox();
}

void GeometryTool_set_merge_test_internal(int tof)
/* < API function for GeometryTool::set_merge_test_internal */
{
  GTI->set_merge_test_internal(tof);
}

int GeometryTool_get_merge_test_internal()
/* < API function for GeometryTool::get_merge_test_internal */
{
  return GTI->get_merge_test_internal();
}

enum CubitBoolean GeometryTool_same_engine(/* DLTopologyEntityList & */ void ***topo_list, int *topo_list_size) 
{
  DLTopologyEntityList temp_topo_list;
  COPY_ARRAY_TO_LIST(*topo_list, *topo_list_size, temp_topo_list);
  
  return GTI->same_engine(temp_topo_list);
}
    
  /* GeometricModelingEngine* */ void *GeometryTool_common_engine( 
        /* DLTopologyEntityList& */ void ***topology_list, int *topology_list_size,
          /* DLTopologyBridgeList& */ void ***engine_bridges, int *engine_bridges_size,
        enum CubitBoolean allow_virtual_engine) 
{
  DLTopologyEntityList temp_topology_list;
  COPY_ARRAY_TO_LIST(*topology_list, *topology_list_size, temp_topology_list);
  
  DLTopologyBridgeList temp_engine_bridges;
  void *dummy = GTI->common_engine(temp_topology_list, temp_engine_bridges, allow_virtual_engine);
  
  COPY_LIST_TO_ARRAY(temp_engine_bridges, *engine_bridges, *engine_bridges_size);
  
  return dummy;
}

  void GeometryTool_set_sep_after_webcut(enum CubitBoolean val) {GTI->set_sep_after_webcut(val);}
  enum CubitBoolean GeometryTool_get_sep_after_webcut() {return GTI->get_sep_after_webcut();}
