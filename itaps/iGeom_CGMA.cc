/**
 * Copyright 2006 Sandia Corporation.  Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Coroporation, the U.S. Government
 * retains certain rights in this software.
 * 
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * 
 */
 
/**\file iGeom_CGMA.cc
 *
 *\author Tim Tautges
 *\author Jason Kraftcheck
 *
 * Original file from SNL TSTT repository was named TSTTG_CGM.cpp.
 *
 * Renamed iGeom_CGMA.cc and added to ANL ITAPS repository by J.Kraftcheck,
 * 2007-6-15
 */
#include "iGeom.h"
#include "InitCGMA.hpp"

#include <iostream>
#include <math.h>
#include "GeometryQueryTool.hpp"

#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

const bool debug = false;

#define RETURN(a) do {iGeom_setLastError((*err = a)); return; } while(false)
#define ERROR(a, b) do {iGeom_setLastError((*err = a), b); return; } while(false)

#define ARRAY_IN(b) \
  b, b ## _size  

#define ARRAY_INOUT(b) \
  b, b ## _allocated, b ## _size  

// Check the array size, and allocate the array if necessary.
// Free the array upon leaving scope unless KEEP_ARRAY
// is invoked.
#define ALLOC_CHECK_ARRAY(array, this_size) \
  iGeomArrayManager array ## _manager ( reinterpret_cast<void**>(array), *(array ## _allocated), *(array ## _size), this_size, sizeof(**array), err ); \
  if (iBase_SUCCESS != *err) return

#define KEEP_ARRAY(array) \
  array ## _manager .keep_array()

// Check the array size, and allocate the array if necessary.
// Do NOT free the array upon leaving scope.
#define ALLOC_CHECK_ARRAY_NOFAIL(array, this_size) \
  ALLOC_CHECK_ARRAY(array, this_size); KEEP_ARRAY(array)

#define TAG_HANDLE(handle) reinterpret_cast<long>(handle)

#define ENTITY_HANDLE(handle) reinterpret_cast<RefEntity*>(handle)
#define SET_HANDLE(handle) reinterpret_cast<RefGroup*>(handle)

#define ENTITY_HANDLE_ARRAY(array) reinterpret_cast<RefEntity**>(array)
#define SET_HANDLE_ARRAY(array) reinterpret_cast<RefGroup**>(array)

#define ENTITY_HANDLE_CONST_ARRAY(array) reinterpret_cast<RefEntity* const*>(array)
#define SET_HANDLE_CONST_ARRAY(array) reinterpret_cast<RefGroup* const*>(array)

#define ENTITY_HANDLE_ARRAY_PTR(array) reinterpret_cast<RefEntity***>(array)
#define SET_HANDLE_ARRAY_PTR(array) reinterpret_cast<RefGroup***>(array)

#define CAST_TO_VOID(ptr) reinterpret_cast<void*>(ptr)
#define TAG_HANDLE_ARRAY_INOUT(a) reinterpret_cast<long**>(a), a ## _allocated, a ## _size

#define TM reinterpret_cast<CGMTagManager*>(instance)

const double RAD_TO_DEG = 180.0 / acos(-1.0);
const double DEG_TO_RAD = 1.0 / RAD_TO_DEG;

#include "CATag.hpp"
#include "CGMAIterator.hpp"
#include "iGeomError.h"

#include "RefEntityFactory.hpp"
#include "BasicTopologyEntity.hpp"
#include "RefGroup.hpp"
#include "Body.hpp"
#include "RefVertex.hpp"
#include "RefEdge.hpp"
#include "RefFace.hpp"
#include "RefVolume.hpp"
#include "CubitVector.hpp"

#ifdef USE_MPI
#ifdef HAVE_OCC
#include "CGMReadParallel.hpp"
#endif
#endif

#include "CGMApp.hpp"
#include "GeometryModifyTool.hpp"
#include "Surface.hpp"
#include "BasicTopologyEntity.hpp"
#include "CubitFile.hpp"
#include "CubitDefines.h"
#include "MergeTool.hpp"
#include "GeometryDefines.h"

#define gqt GeometryQueryTool::instance()
#define gmt GeometryModifyTool::instance()

const char *iGeom_entity_type_names[] = {"vertex", "curve", "surface", "body"};


// Implement RAII pattern for allocated arrays
class iGeomArrayManager
{
  void** arrayPtr;

public:


  iGeomArrayManager( void** array_ptr,
                     int& array_allocated_space,
                     int& array_size,
                     int count,
                     int val_size,
                     int* err ) : arrayPtr(0)
  {
    if (!*array_ptr) {
      *array_ptr = malloc(val_size * count);
      array_allocated_space = array_size = count;
      if (!*array_ptr) {
        ERROR(iBase_MEMORY_ALLOCATION_FAILED, "Couldn't allocate array.");
      }
      arrayPtr = array_ptr;
    }
    else {
      array_size = count;
      if (array_allocated_space < count) {
        ERROR(iBase_BAD_ARRAY_DIMENSION, 
          "Allocated array not large enough to hold returned contents.");
      }
    }
    RETURN(iBase_SUCCESS);
  }
  
  ~iGeomArrayManager() 
  {
    if (arrayPtr) {
      free(*arrayPtr);
      *arrayPtr = 0;
    }
  }
  
  void keep_array()
    { arrayPtr = 0; }
};


// declare private-type functions here, so they aren't visible outside
// this implementation file

/* Helper functions 
 * Declare these outside the 'extern "C"' block because they are not
 * extern.
 */
static void
iGeom_get_adjacent_entities( const RefEntity *from, 
                             const int to_dim,
                             DLIList<RefEntity*> &adj_ents,
                             int* err);

static void tokenize( const std::string& str, 
                      std::vector<std::string>& tokens );

// Expect option of the form "NAME=VALUE".
// If NAME portion matches, pass back VALUE and return true.
// Otherwise, leave 'value' unchanged and return false.
static bool match_option( const std::string& opt,
                          const char* name,
                          std::string& value );
                          
static void
iGeom_load_cub_geometry(const char *name, int* err) ;

static iBase_ErrorType 
process_attribs(iGeom_Instance instance) ;

static void
iGeom_get_adjacent_entities( const RefEntity *from, 
                             const int to_dim,
                             DLIList<RefEntity*> &adj_ents,
                             int* err ) ;
/* Common implementation for both single-entity and array functions. */
static CubitStatus iGeom_closest_point( RefEntity* this_entity,
                                        const CubitVector& near,
                                        CubitVector& on );

static CubitStatus
iGeom_closest_point_and_normal( RefEntity* this_entity, 
                                const CubitVector& near,
                                CubitVector& on,
                                CubitVector& normal );
static CubitStatus
iGeom_bounding_box( RefEntity* entity, CubitVector& minc, CubitVector& maxc );

/** Smits' algorithm */
static inline void box_min_max( double dir,
                                double min,
                                double max,
                                double pt,
                                double& tmin,
                                double& tmax );
static bool
iBase_intersect_ray_box( const CubitBox& box,
                         const CubitVector& point,
                         const CubitVector& direction );

static CubitStatus
iGeom_fire_ray( const CubitVector& point,
                const CubitVector& direction,
                DLIList<RefEntity*>& entities,
                DLIList<double>& ray_params );

static RefEntity*
iGeom_get_point_containment( const CubitVector& pt );

static int iGeom_get_nonmanifold_sense( const BasicTopologyEntity* child,
                                        const BasicTopologyEntity* parent,
                                        int* err );

static
int iGeom_edge_vertex_sense( const RefEdge* cedge,
                             const RefVertex* vtx1,
                             const RefVertex* vtx2,
                             int* err );

static int iGeom_is_parametric( RefEntity* entity );

static iBase_ErrorType
iGeom_get_vtx_to_u(RefVertex* vertex, RefEdge* edge, double& u);

static iBase_ErrorType
iGeom_get_vtx_to_uv(RefVertex* vertex, RefFace* face, double& u, double& uv);


static CubitStatus
iGeom_normal_from_uv( RefFace* face, double u, double v, CubitVector& normal );

static CubitStatus iGeom_is_periodic( RefEntity* entity, int& u, int& v );

static bool iGeom_is_face_degenerate( RefFace* face );

static
int count_ibase_type( int ibase_type, 
                      const DLIList<CubitEntity*>& list, 
                      int* err );

static 
void copy_ibase_type( int ibase_type, 
                      const DLIList<CubitEntity*>& list,
                      iBase_EntityHandle** entity_handles,
                      int* entity_handles_alloc,
                      int* entity_handles_size,
                      int* err );

static 
void append_ibase_type( int ibase_type, 
                        const DLIList<CubitEntity*>& source_list,
                        DLIList<RefEntity*>& target_list,
                        int* err );
static 
void append_all_ibase_type( int ibase_type, 
                            DLIList<RefEntity*>& target_list,
                            int* err );


static CubitStatus init_cgm( const std::string& engine )
{
  CubitStatus status;
  if (engine.empty()) 
    status = InitCGMA::initialize_cgma();
  else
    status = InitCGMA::initialize_cgma( engine.c_str() );
 
// sometimes can't have following, depending on CGM version
  // CGMApp::instance()->attrib_manager()->silent_flag(true);

  CGMApp::instance()->attrib_manager()->auto_flag(true);
  
  return status;
}


extern "C" {

void iGeom_getDescription( iGeom_Instance geom,
                           char* descr,
                           int* err,
                           int descr_len )
{
  iGeom_getLastError( *err, descr, descr_len );
}

void iGeom_newGeom( const char* options,
                    iGeom_Instance* instance_out,
                    int* err,
                    const int options_size) 
{
    // scan options for default engine option
  std::string engine;
  if (options && options_size) {
    std::string tmp(options, options_size);
    char f[] = ";engine="; f[0]=tmp[0]; // correct delimiter
    size_t p = tmp.find( f );
    if (p != std::string::npos) { // if we found engine option
      p += strlen(f); // advance to value (past '=')
      size_t e = tmp.find( tmp[0], p ); // find end delimiter
      if (e == std::string::npos) // if no end delim, then must be last option
        engine = tmp.substr( p, std::string::npos );
      else
        engine = tmp.substr( p, e-p );
    }
  }
  
    // initialize static var with result so that call happens only once
  static const CubitStatus status = init_cgm( engine );
    // but check the result for every call
  if (CUBIT_SUCCESS != status)
    RETURN (iBase_FAILURE);

    // return the tagmanager as the instance
  *instance_out = reinterpret_cast<iGeom_Instance>(&CGMTagManager::instance());
  RETURN(iBase_SUCCESS);
}


void
iGeom_dtor(iGeom_Instance instance, int* err) 
{
  if (NULL == TM) 
    ERROR(iBase_INVALID_ARGUMENT, "NUll Instance");
  
    // delete TM;
  
    // shut down CGM
  CGMApp::instance()->shutdown();

  RETURN(iBase_SUCCESS);
}
  
// user defined non-static methods:
/**
 * Load a model specified by name. Which formats are supported and the
 * specific meaning of this name string (e.g. file name, model name,
 * etc.) are implementation-dependent.  Options are also implementation-
 * dependent.
 * @param name Name of the model
 * @param options String options 
 */
void iGeom_load( iGeom_Instance instance,
                 /*in*/ const char *name,
                 /*in*/ const char *options,
                 int* err,
                 int name_len,
                 int options_size )
{
   // make sure we have a null-terminated string for file name
  std::string file_name( name, name_len );
  name = file_name.c_str();

  // check if work in parallel
  bool parallel = false;
  std::string parallel_opt;
  std::vector<std::string> opts;
  std::string tmp_options(options, options_size);
  tmp_options.insert(0, ";");
  tokenize(tmp_options, opts );
  for (std::vector<std::string>::iterator i = opts.begin(); i != opts.end(); ++i)
  {
    if (match_option( *i, "PARALLEL", parallel_opt )) parallel = true;
  }
  
  // parallel
  if (parallel) {
#ifdef USE_MPI
#ifdef HAVE_OCC
    CGMParallelComm* p_comm = new CGMParallelComm();
    CGMReadParallel* p_reader = new CGMReadParallel(GeometryQueryTool::instance(), p_comm);
    CubitStatus status = p_reader->load_file(name, options);
    if (CUBIT_SUCCESS != status) {
      ERROR(iBase_FAILURE, "Trouble loading geometry file in parallel.");
    }
#else
    ERROR(iBase_NOT_SUPPORTED, "Parallel not enabled in this version.");
#endif
#endif
  }
  else {
    CubitStatus status = CUBIT_SUCCESS;
    
    if (strstr(name, ".cub") != NULL) {
      iGeom_load_cub_geometry(name, err);
      if (iBase_SUCCESS != *err) {
	return;
      }
    }
    else {
      std::string file_type;
      if (strstr(name, ".brep") != NULL ||
	  strstr(name, ".BREP") != NULL ||
	  strstr(name, ".occ") != NULL ||
	  strstr(name, ".OCC") != NULL)
	file_type = "OCC";
      else if (strstr(name, ".stp") != NULL ||
	       strstr(name, ".STP") != NULL ||
	       strstr(name, ".step") != NULL ||
	       strstr(name, ".STEP") != NULL)
	file_type = "STEP";
      else if (strstr(name, ".igs") != NULL ||
	       strstr(name, ".IGS") != NULL ||
	       strstr(name, ".iges") != NULL ||
	       strstr(name, ".IGES") != NULL)
	file_type = "IGES";
      
      if (file_type.empty()) status = gqt->read_geometry_file(name);
      else status = gqt->read_geometry_file(name, NULL, file_type.c_str());

      if (CUBIT_SUCCESS != status) {
	ERROR(iBase_FILE_NOT_FOUND, "Trouble loading geometry file.");
      }
    }
  }

    // now process uncaught/unregistered attributes
  iBase_ErrorType result;
  result = process_attribs(instance);

  RETURN(result);
}

// user defined non-static methods:
/**
 * Save a model to file specified by name. Which formats are supported and the
 * specific meaning of this name string (e.g. file name, model name,
 * etc.) are implementation-dependent.  Options are also implementation-
 * dependent.
 * @param name Name of the model
 * @param options String options 
 */
void iGeom_save (iGeom_Instance instance,
                 /*in*/ const char *name,
                 /*in*/ const char* options,
                 int* err,
                 int name_len,
                 int options_len)
{
    // make sure strings are null terminated
  std::string name_buf(name, name_len);
  name = name_buf.c_str();
  
    // parse options
  std::string file_type;
  std::vector<std::string> opts;
  tokenize( std::string(options, options_len), opts );
  for (std::vector<std::string>::iterator i = opts.begin(); i != opts.end(); ++i)
  {
    if (!match_option( *i, "TYPE", file_type )) // e.g TYPE=ACIS_SAT
      ERROR( iBase_INVALID_ARGUMENT, i->c_str() );
  }

    // if no type, check file type for known extensions
  if (file_type.empty()) {
    size_t name_len = name_buf.length();
    if (name_buf.find(".igs") < name_len ||
        name_buf.find(".IGS") < name_len ||
        name_buf.find(".iges") < name_len ||
        name_buf.find(".IGES") < name_len)
      file_type = "IGES";
    else if (name_buf.find(".stp") < name_len ||
             name_buf.find(".STP") < name_len ||
             name_buf.find(".step") < name_len ||
             name_buf.find(".STEP") < name_len)
      file_type = "STEP";
    else if (name_buf.find(".sat") < name_len ||
             name_buf.find(".SAT") < name_len)
      file_type = "ACIS_SAT";
    else if (name_buf.find(".occ") < name_len ||
             name_buf.find(".brep") < name_len ||
             name_buf.find(".OCC") < name_len ||
             name_buf.find(".BREP") < name_len)
      file_type = "OCC";
    else
      ERROR(iBase_FAILURE, "Unknown geometry file extension and no file type.");
  }
  
  /** 
   * Work around AcisQueryEngine log file handling bug:
   * If the export format is iges or step, be sure the log file name is not null
   */
  const char* logfile_name = NULL;
  if( file_type == "IGES" ){
    logfile_name = "igeom_iges_export.log";
  }
  else if( file_type == "STEP" ){
    logfile_name = "igeom_step_export.log";
  }

    // process options (none right now...)
  DLIList<RefEntity*> bodies;
  int num_ents_exported;
  CubitString cubit_version(" (iGeom)");
  CubitStatus status = gqt->export_solid_model(bodies, name, file_type.c_str(),
                                               num_ents_exported, cubit_version, logfile_name );
  if (CUBIT_SUCCESS != status) 
    ERROR(iBase_FAILURE, "Trouble saving geometry file.");

  RETURN(iBase_SUCCESS);
}

void iGeom_getRootSet( iGeom_Instance, iBase_EntitySetHandle* root, int* err )
{
  *root = NULL;
  RETURN(iBase_SUCCESS);
}


void iGeom_getBoundBox( iGeom_Instance,
                        double* min_x,
                        double* min_y,
                        double* min_z,
                        double* max_x,
                        double* max_y,
                        double* max_z,
                        int* err ) 
{
  CubitBox box = GeometryQueryTool::instance()->model_bounding_box();
  *min_x = box.min_x();
  *min_y = box.min_y();
  *min_z = box.min_z();
  *max_x = box.max_x();
  *max_y = box.max_y();
  *max_z = box.max_z();
  RETURN(iBase_SUCCESS);
}


/**
 * Initialize an iterator over gentities of a specified dimension.
 * @param gentity_dimension Dimension of gentities to be iterated over
 * @param gentity_iterator Iterator initialized by this function
 */

void iGeom_initEntIter ( iGeom_Instance instance,
                         iBase_EntitySetHandle entity_set_handle,
                         const int dimension,
                         iGeom_EntityIterator *iterator,
                         int* err)
{
  iGeom_initEntArrIter( instance, entity_set_handle, dimension, 1, 
                        reinterpret_cast<iGeom_EntityArrIterator*>(iterator), 
                        err );
}

void
iGeom_initEntArrIter (iGeom_Instance instance,
                      /*in*/ iBase_EntitySetHandle entity_set_handle,
                      /*in*/ int type,
                      /*in*/ int requested_array_size,
                      /*out*/ iGeom_EntityArrIterator* entArr_iterator,
                      int* err)
{
  DLIList<RefEntity*> entities;
  if (RefGroup* group = reinterpret_cast<RefGroup*>(entity_set_handle)) {
    DLIList<CubitEntity*> centities;
    group->get_child_entities( centities );
    append_ibase_type( type, centities, entities, err );
  }
  else {
    append_all_ibase_type( type, entities, err );
  }
  
  if (*err == iBase_SUCCESS) {
    CGMAIterator* iter = new CGMAIterator(entities, requested_array_size);
    *entArr_iterator = reinterpret_cast<iGeom_EntityArrIterator>(iter);
  }
}



/**
 * Get the next entity for this iterator.
 * @param gentity_iterator Iterator being iterated over
 * @param gentity_handle Next gentity
 * @return If true, there were more gentities, if false, the iterator was
 *         already at the end.
 */
void
iGeom_getNextEntIter (iGeom_Instance instance,
                      /*in*/ iGeom_EntityIterator gentity_iterator,
                      /*out*/ iBase_EntityHandle *gentity_handle,
                      int* has_data,
                      int* err
                      )
{
  CGMAIterator* iterator = reinterpret_cast<CGMAIterator*>(gentity_iterator);
  RefEntity** out_handle = reinterpret_cast<RefEntity**>(gentity_handle);
  *has_data = !iterator->at_end();
  if (*has_data) {
    int count = 1;
    iterator->next( out_handle, count );
  }
  RETURN(iBase_SUCCESS);
}


void
iGeom_getNextEntArrIter(iGeom_Instance instance,
                        /*in*/ iGeom_EntityArrIterator entArr_iterator,
                        /*inout*/ iBase_EntityHandle **entity_handles,
                        int *entity_handles_allocated,
                        int *entity_handles_size,
                        int* has_data,
                        int* err
                        )
{
  CGMAIterator* iterator = reinterpret_cast<CGMAIterator*>(entArr_iterator);
  *has_data = !iterator->at_end();
  if (has_data) {
    ALLOC_CHECK_ARRAY_NOFAIL(entity_handles, iterator->size());
    iterator->next( (RefEntity**)*entity_handles, *entity_handles_size );
  }
  RETURN(iBase_SUCCESS);
}


/**
 * Reset an iterator back to the first gentity
 * @param gentity_iterator Iterator reset by this function
 */
void
iGeom_resetEntIter(iGeom_Instance instance,
                   /*in*/ iGeom_EntityIterator gentity_iterator,
                   int* err
                   )
{
  CGMAIterator* iterator = reinterpret_cast<CGMAIterator*>(gentity_iterator);
  iterator->reset();
  RETURN(iBase_SUCCESS);
}

void
iGeom_resetEntArrIter(iGeom_Instance instance,
                      /*in*/ iGeom_EntityArrIterator gentity_iterator,
                      int* err
                      )
{
  CGMAIterator* iterator = reinterpret_cast<CGMAIterator*>(gentity_iterator);
  iterator->reset();
  RETURN(iBase_SUCCESS);
}

/**
 * Delete an iterator
 * @param gentity_iterator Iterator deleted by this function
 */
void
iGeom_endEntIter (iGeom_Instance instance,
                  /*in*/ iGeom_EntityIterator gentity_iterator,
                  int* err)
{
  CGMAIterator* iterator = reinterpret_cast<CGMAIterator*>(gentity_iterator);
  delete iterator;
  RETURN(iBase_SUCCESS);
}

void
iGeom_endEntArrIter (iGeom_Instance instance,
                     /*in*/ iGeom_EntityArrIterator gentity_iterator,
                     int* err)
{
  CGMAIterator* iterator = reinterpret_cast<CGMAIterator*>(gentity_iterator);
  delete iterator;
  RETURN(iBase_SUCCESS);
}

/**
 *   Returns true if the gentity_sets are related through a parent/child
 *   relationship.
 */
void
iGeom_isChildOf (iGeom_Instance instance,
                 /*in*/ iBase_EntitySetHandle parent_entity_set,
                 /*in*/ iBase_EntitySetHandle child_entity_set,
                 int* is_child,
                 int* err)
{
  std::vector<RefGroup*> *par1, *ch1, *par2, *ch2;
  TM->pc_list(const_cast<RefGroup*>(SET_HANDLE(parent_entity_set)), par1, ch1, false);
  if (NULL == ch1 || ch1->empty()) {
    *is_child = false;
    RETURN (iBase_SUCCESS);
  }
  TM->pc_list(const_cast<RefGroup*>(SET_HANDLE(child_entity_set)), par2, ch2, false);
  if (NULL == par2 || par2->empty()) {
    *is_child = false;
    RETURN (iBase_SUCCESS);
  }
  
  const RefGroup *group1 = SET_HANDLE(parent_entity_set);
  const RefGroup *group2 = SET_HANDLE(child_entity_set);
  if ((std::find(ch1->begin(), ch1->end(), group2) != ch1->end())
      || (std::find(par2->begin(), par2->end(), group1) != par2->end()))
    *is_child = true;
  
  else *is_child = false;
  RETURN(iBase_SUCCESS);
}

/**
 *   Recursively gets the children of this gentity_set up to num_hops
 *   levels; if num_hops is set to -1 all children are returned
 */
void
iGeom_getChldn (iGeom_Instance instance,
                /*in*/ iBase_EntitySetHandle from_entity_set,
                /*in*/ const int num_hops,
                /*inout*/ iBase_EntitySetHandle** entity_set_handles,
                int* entity_set_handles_allocated,
                int* entity_set_handles_size,
                int* err)
{
  std::vector<RefGroup*> group_ptrs;
  const RefGroup *this_grp = SET_HANDLE(from_entity_set);
  TM->get_pc_groups(const_cast<RefGroup*>(this_grp), 1, num_hops, group_ptrs);
  ALLOC_CHECK_ARRAY_NOFAIL(entity_set_handles, group_ptrs.size());

  iBase_EntitySetHandle* ent_arr = reinterpret_cast<iBase_EntitySetHandle*>(&group_ptrs[0]);
  std::copy( ent_arr, ent_arr + group_ptrs.size(), *entity_set_handles);
  
  RETURN(iBase_SUCCESS);
}

/**
 *   Recursively gets the parents of this gentity_set up to num_hops
 *   levels; if num_hops is set to -1 all parents are returned
 */
void
iGeom_getPrnts (iGeom_Instance instance,
                /*in*/ iBase_EntitySetHandle from_entity_set,
                /*in*/ const int num_hops,
                /*inout*/ iBase_EntitySetHandle **entity_set_handles,
                int *entity_set_handles_allocated,
                int *entity_set_handles_size,
                int* err)
{
  std::vector<RefGroup*> group_ptrs;
  const RefGroup *this_grp = SET_HANDLE(from_entity_set);
  TM->get_pc_groups(const_cast<RefGroup*>(this_grp), 0, num_hops, group_ptrs);
  ALLOC_CHECK_ARRAY_NOFAIL(entity_set_handles, group_ptrs.size());

  iBase_EntitySetHandle* ent_arr = reinterpret_cast<iBase_EntitySetHandle*>(&group_ptrs[0]);
  std::copy( ent_arr, ent_arr + group_ptrs.size(), *entity_set_handles);
  RETURN(iBase_SUCCESS);
}

/**
 *   Returns the number of immediate children in the gentity_set (one
 *   level down only)
 */
void
iGeom_getNumChld (iGeom_Instance instance,
                  /*in*/ iBase_EntitySetHandle entity_set,
                  /*in*/ const int num_hops,
                  int* num_child, 
                  int* err)
{
    // HJK: num_hops has to be handled
  if (1 < num_hops) {
    ERROR(iBase_NOT_SUPPORTED, "Num_hops argument not yet supported.");
  }

  std::vector<RefGroup*> *my_children = TM->pc_list(const_cast<RefGroup*>(SET_HANDLE(entity_set)), 1, false);
  *num_child = my_children == NULL ? 0 : my_children->size();
  RETURN (iBase_SUCCESS);
}

/**
 *   Returns the number of immediate parents to the gentity_set (one
 *   level up only)
 */
void
iGeom_getNumPrnt (iGeom_Instance instance,
                  /*in*/ iBase_EntitySetHandle entity_set,
                  /*in*/ const int num_hops,
                  int* num_parent,
                  int* err)
{
    // HJK: num_hops has to be handled
  if (1 < num_hops) {
    ERROR(iBase_NOT_SUPPORTED, "Num_hops argument not yet supported.");
  }

  std::vector<RefGroup*> *my_parents = TM->pc_list(const_cast<RefGroup*>(SET_HANDLE(entity_set)), 0, false);
  *num_parent = my_parents == NULL ? 0 : my_parents->size();
  RETURN (iBase_SUCCESS);
}

/**
 *   Add a parent to the gentity_set 
 */
void
iGeom_addPrntChld (iGeom_Instance instance,
                   /*inout*/ iBase_EntitySetHandle parent_entity_set,
                   /*inout*/ iBase_EntitySetHandle child_entity_set,
                   int* err)
{
  std::vector<RefGroup*> *my_parents = 
    TM->pc_list(SET_HANDLE(child_entity_set), 0, true);
  std::vector<RefGroup*> *my_children = 
    TM->pc_list(SET_HANDLE(parent_entity_set), 1, true);
  RefGroup *par_group = SET_HANDLE(parent_entity_set);
  RefGroup *child_group = SET_HANDLE(child_entity_set);
  my_parents->push_back(par_group);
  my_children->push_back(child_group);
  RETURN(iBase_SUCCESS);
}

/**
 *   Remove a parent/child link between gentity_sets
 */
void
iGeom_rmvPrntChld (iGeom_Instance instance,
                   /*inout*/ iBase_EntitySetHandle parent_entity_set,
                   /*inout*/ iBase_EntitySetHandle child_entity_set,
                   int* err)
{
  RefGroup *parent = reinterpret_cast<RefGroup *>(parent_entity_set);
  RefGroup *child = reinterpret_cast<RefGroup *>(child_entity_set);
  std::vector<RefGroup*> *children = TM->pc_list(parent, 1, false),
    *parents = TM->pc_list(child, 0, false);
  if (NULL == children || NULL == parents) {
    RETURN(iBase_INVALID_ARGUMENT);
  }
  
  children->erase(std::remove(children->begin(), children->end(), child), children->end());
  parents->erase(std::remove(parents->begin(), parents->end(), parent), parents->end());
  RETURN(iBase_SUCCESS);
}

/**
 * Return gentities of specified dimension in this set, or in whole model.
 * @param set_handle Entity set being queried (if 0, whole model)
 * @param gentity_dimension Dimension of entities being queried
 * @param gentity_handles Entity handles
 */
void
iGeom_getEntities (iGeom_Instance instance,
                   /*in*/ iBase_EntitySetHandle set_handle,
                   /*in*/ int gentity_type,
                   /*out*/ iBase_EntityHandle **gentity_handles,
                   int *gentity_handles_allocated,
                   int *gentity_handles_size,
                   int* err)
{
  if (RefGroup *this_set = SET_HANDLE(set_handle)) {
    static DLIList<CubitEntity*> centities;
    centities.clean_out();
    this_set->get_child_entities(centities);
    copy_ibase_type( gentity_type, centities, 
                     gentity_handles,
                     gentity_handles_allocated,
                     gentity_handles_size,
                     err );
  }
  else {
    static DLIList<RefEntity*> dim_entities;
    dim_entities.clean_out();
    append_all_ibase_type( gentity_type, dim_entities, err );
    if (iBase_SUCCESS != *err)
      return;
    
    ALLOC_CHECK_ARRAY_NOFAIL(gentity_handles, dim_entities.size());
    dim_entities.copy_to((RefEntity**)*gentity_handles);
    RETURN(iBase_SUCCESS);
  }  
}

/**
 * Return number of gentities of specified dimension in this set, or in
 * whole model.
 * @param set_handle Entity set being queried (if 0, whole model)
 * @param gentity_dimension Dimension of entities being queried
 * @return Number of entities
 */
void
iGeom_getNumOfType (iGeom_Instance instance,
                    /*in*/ iBase_EntitySetHandle set_handle,
                    /*in*/ int gentity_type,
                    int* count,
                    int* err)
{
  const RefGroup *this_set = SET_HANDLE(set_handle);
  if (0 == this_set) {
    switch (gentity_type) {
      case iBase_ALL_TYPES:
        *count  = GeometryQueryTool::instance()->num_bodies();
        *count += GeometryQueryTool::instance()->num_ref_faces();
        *count += GeometryQueryTool::instance()->num_ref_edges();
        *count += GeometryQueryTool::instance()->num_ref_vertices();
        break;
      case iBase_REGION:
        *count = GeometryQueryTool::instance()->num_bodies();
        break;
      case iBase_FACE:
        *count = GeometryQueryTool::instance()->num_ref_faces();
        break;
      case iBase_EDGE:
        *count = GeometryQueryTool::instance()->num_ref_edges();
        break;
      case iBase_VERTEX:
        *count = GeometryQueryTool::instance()->num_ref_vertices();
        break;
      default:
        RETURN(iBase_BAD_TYPE_AND_TOPO);
        break;
    }
    RETURN (iBase_SUCCESS);
  }
  else {
    static DLIList<CubitEntity*> centities;
    centities.clean_out();
    const_cast<RefGroup*>(this_set)->get_child_entities(centities);
    *count = count_ibase_type( gentity_type, centities, err );
  }
}

void
iGeom_getEntType (iGeom_Instance instance, 
                  iBase_EntityHandle handle,
                  int* gtype,
                  int* err)
{
  RefEntity* entity = reinterpret_cast<RefEntity*>(handle);
  if (dynamic_cast<Body*>(entity))
    *gtype = iBase_REGION;
  else
    *gtype = static_cast<int>(entity->dimension());
  RETURN(iBase_SUCCESS);
}

/**
 *    Returns an integer array of topological dimensions for an input
 *    array of entity handles.
 */
void
iGeom_getArrType (iGeom_Instance instance,
                  /*in*/ iBase_EntityHandle const *gentity_handles,
                  int gentity_handles_size,
                  /*inout*/ int **gtype,
                  int *gtype_allocated,
                  int *gtype_size, 
                  int* err)
{
    // go through each entity and look up its dimension
  ALLOC_CHECK_ARRAY_NOFAIL(gtype, gentity_handles_size);

  const RefEntity **tmp_handles = (const RefEntity**)(gentity_handles);
  
  for (int i = 0; i < gentity_handles_size; i++) {
    if (dynamic_cast<const Body*>(tmp_handles[i]))
      (*gtype)[i] = iBase_REGION;
    else
      (*gtype)[i] = tmp_handles[i]->dimension();
  }

  RETURN(iBase_SUCCESS);
}

/**
 * Get the adjacent entities of a given dimension.
 * @param gentity_handle Entity for which adjacencies are requested
 * @param to_dimension Target dimension of adjacent entities
 * @param adj_gentities List returned with adjacent entities
 */
void
iGeom_getEntAdj (iGeom_Instance instance,
                 /*in*/ iBase_EntityHandle gentity_handle,
                 /*in*/int to_dimension,
                 /*inout*/ iBase_EntityHandle **adj_gentities,
                 int *adj_gentities_allocated,
                 int *adj_gentities_size,
                 int* err)
{
  const RefEntity *tmp_hndl = ENTITY_HANDLE(gentity_handle);

  if (tmp_hndl->dimension() == to_dimension) {
    ERROR(iBase_INVALID_ARGUMENT, "Can't get adjacencies to entities of same dimension.");
  }

  static DLIList<RefEntity*> tmp_ents;
  iGeom_get_adjacent_entities(tmp_hndl, to_dimension, tmp_ents, err);
  if (iBase_SUCCESS != *err) return;

  ALLOC_CHECK_ARRAY_NOFAIL(adj_gentities, tmp_ents.size());
  tmp_ents.copy_to((RefEntity**)*adj_gentities);
  RETURN(iBase_SUCCESS);
}

void
iGeom_getArrAdj(iGeom_Instance instance,
                /*in*/ iBase_EntityHandle const *entity_handles,
                const int entity_handles_size,
                /*in*/ int requested_entity_type,
                /*inout*/ iBase_EntityHandle **adj_entity_handles,
                int *adj_entity_handles_allocated,
                int *adj_entity_handles_size,
                /*inout*/ int **offset,
                int *offset_allocated,
                int *offset_size,
                int* err)
{
  ALLOC_CHECK_ARRAY(offset, entity_handles_size+1);
  DLIList<RefEntity*> temp_list, total_list;
  for (int i = 0; i < entity_handles_size; ++i) {
    (*offset)[i] = total_list.size();
    temp_list.clean_out();
    iGeom_get_adjacent_entities( (RefEntity*)(entity_handles[i]), requested_entity_type, temp_list, err );
    if (iBase_SUCCESS != *err) return;
    total_list += temp_list;
  }
  (*offset)[entity_handles_size] = total_list.size();

  ALLOC_CHECK_ARRAY_NOFAIL(adj_entity_handles, total_list.size());
  total_list.copy_to((RefEntity**)*adj_entity_handles);
  KEEP_ARRAY(offset);
  RETURN(iBase_SUCCESS);
}
  
/**
 * Get the "2nd order" adjacent entities, through a specified "bridge"
 * dimension, of a target dimension.  For example, given a region, return
 * the regions (to_dimension=3) sharing an edge (bridge_dimension=1)
 * with that region.  bridge_dimension must be less than dimension of 
 * gentity_handle, and to_dimension must be greater than bridge dimension.
 * 
 * @param gentity_handle Entity for which 2nd order adjacencies are requested
 * @param to_dimension Target dimension of 2nd order adjacent entities
 * @param bridge_dimension Dimension of "bridge" entities
 * @param adj_gentities List returned with 2nd order adjacent entities
 */
void
iGeom_getEnt2ndAdj (iGeom_Instance instance,
                    /*in*/ iBase_EntityHandle gentity_handle,
                    /*in*/ int bridge_dimension,
                    /*in*/ int to_dimension,
                    /*out*/ iBase_EntityHandle **adjacent_gentities,
                    int *adjacent_gentities_allocated,
                    int *adjacent_gentities_size,
                    int* err)
{
    // for better or worse, go to cgm as quickly as possible, to avoid working with
    // sidl arrays
  const RefEntity *gentity = ENTITY_HANDLE(gentity_handle);
  DLIList<RefEntity*> to_ents, bridge_ents, tmp_ents;

  iGeom_get_adjacent_entities(gentity, bridge_dimension, bridge_ents, err);
  if (iBase_SUCCESS != *err) return;
  
  for (int i = bridge_ents.size(); i > 0; i--) {
    iGeom_get_adjacent_entities(bridge_ents.get_and_step(), to_dimension, tmp_ents, err);
    if (iBase_SUCCESS != *err) return;

    to_ents += tmp_ents;
    tmp_ents.clean_out();
  }

  to_ents.uniquify_unordered();

  ALLOC_CHECK_ARRAY_NOFAIL(adjacent_gentities, to_ents.size());
  to_ents.copy_to((RefEntity**)*adjacent_gentities);

  RETURN(iBase_SUCCESS);
}

void
iGeom_getArr2ndAdj(iGeom_Instance instance,
                   /*in*/ iBase_EntityHandle const *entity_handles,
                   int entity_handles_size,
                   /*in*/ int order_adjacent_key,
                   /*in*/ int requested_entity_type,
                   /*inout*/ iBase_EntityHandle **adj_entity_handles,
                   int *adj_entity_handles_allocated,
                   int *adj_entity_handles_size,
                   /*inout*/ int **offset,
                   int *offset_allocated,
                   int *offset_size,
                   int *err)
{
  ALLOC_CHECK_ARRAY(offset, entity_handles_size+1);
  DLIList<RefEntity*> bridge_list, temp_list, entity_list, total_list;
   
  for (int i = 0; i < entity_handles_size; ++i) {
    bridge_list.clean_out();
    entity_list.clean_out();
    iGeom_get_adjacent_entities( (RefEntity*)(entity_handles[i]), order_adjacent_key, bridge_list, err );
    if (iBase_SUCCESS != *err) return;
    bridge_list.reset();
    for (int j = bridge_list.size(); j > 0; --j) {
      temp_list.clean_out();
      iGeom_get_adjacent_entities( bridge_list.get_and_step(), requested_entity_type, temp_list, err );
      if (iBase_SUCCESS != *err) return;
      entity_list += temp_list;
    }
    entity_list.uniquify_unordered();
    (*offset)[i] = total_list.size();
    total_list += entity_list;
  }
  (*offset)[entity_handles_size] = total_list.size();

  ALLOC_CHECK_ARRAY_NOFAIL(adj_entity_handles, total_list.size());
  total_list.copy_to((RefEntity**)*adj_entity_handles);
  KEEP_ARRAY(offset);
  RETURN(iBase_SUCCESS);
}

/**
 * Return whether or not entities are adjacent.
 * @param gentity_handle1 1st entity
 * @param gentity_handle2 2nd entity
 * @param are_adjacent If true, entities are adjacent
 */
void
iGeom_isEntAdj (iGeom_Instance instance,
                /*in*/ iBase_EntityHandle gentity_handle1,
                /*in*/ iBase_EntityHandle gentity_handle2,
                /*out*/ int *are_adjacent,
                int* err)
{
  const TopologyEntity *ent1 = dynamic_cast<const TopologyEntity*>(ENTITY_HANDLE(gentity_handle1));
  const TopologyEntity *ent2 = dynamic_cast<const TopologyEntity*>(ENTITY_HANDLE(gentity_handle2));
  if (ent1 != NULL && ent2 != NULL)
    *are_adjacent = const_cast<TopologyEntity*>(ent1)->is_directly_related(const_cast<TopologyEntity*>(ent2));
  else *are_adjacent = false;
  RETURN(iBase_SUCCESS);
}

void
iGeom_isArrAdj (iGeom_Instance instance,
                /*in*/ iBase_EntityHandle const *entity_handles_1,
                int entity_handles_1_size,
                /*in*/ iBase_EntityHandle const *entity_handles_2,
                int entity_handles_2_size,
                /*inout*/ int **is_adjacent_info,
                int *is_adjacent_info_allocated,
                int *is_adjacent_info_size,
                int* err)
{
  RefEntity **list_1_iter = (RefEntity**)entity_handles_1, 
    **list_2_iter = (RefEntity**)entity_handles_2;
  size_t list_1_step, list_2_step;
  int count;
    // If either list contains only 1 entry, compare that entry with
    // every entry in the other list.
  if (entity_handles_1_size == entity_handles_2_size) {
    list_1_step = list_2_step = 1;
    count = entity_handles_1_size;
  }
  else if (entity_handles_1_size == 1) {
    list_1_step = 0;
    list_2_step = 1;
    count = entity_handles_2_size;
  }
  else if (entity_handles_2_size == 1) {
    list_1_step = 1;
    list_2_step = 0;
    count = entity_handles_1_size;
  }
  else {
    RETURN(iBase_INVALID_ENTITY_COUNT);
  }
  
  ALLOC_CHECK_ARRAY_NOFAIL(is_adjacent_info, count);
  for (int i = 0; i < count; ++i)
  {
    TopologyEntity* ent1 = dynamic_cast<TopologyEntity*>(*list_1_iter);
    TopologyEntity* ent2 = dynamic_cast<TopologyEntity*>(*list_2_iter);
    (*is_adjacent_info)[i] = ent1->is_directly_related(ent2);
    list_1_iter += list_1_step;
    list_2_iter += list_2_step;
  }
  RETURN(iBase_SUCCESS);
}

void
iGeom_getTopoLevel (iGeom_Instance instance,
                     /*out*/ int* level,
                     int* err)
{
  *level = 2; // 0->basic entities only, 1->manifold, 2->non-manifold
  RETURN(iBase_SUCCESS);
}

/**
 * 
 *   Create a tag handle with a given name, size (in bytes), and
 *   default value.  The tag name is a unique string; if it duplicates
 *   an existing tag name, an error is returned.  The tag_handle is
 *   returned as an opaque value which is not associated with any mesh
 *   gentities until explicitly done so through one of the 'AddTag'
 *   functions defined later.  The implementation is assumed to
 *   allocate memory as needed to store the tag data.
 */
void
iGeom_createTag (iGeom_Instance instance,
                 /*in*/ const char *tag_name,
                 /*in*/ int tag_size,
                 /*in*/ int tag_type,
                 /*out*/ iBase_TagHandle *tag_handle,
                 int* err,
                 int tag_name_len)
{
  long new_tag;
  int this_size = tag_size;
  switch (tag_type) {
    case iBase_INTEGER:
      this_size *= sizeof(int);
      break;
    case iBase_DOUBLE:
      this_size *= sizeof(double);
      break;
    case iBase_ENTITY_HANDLE:
      this_size *= sizeof(iBase_EntityHandle);
      break;
    case iBase_BYTES:
      break;
  }
  
    // make sure string is null terminated
  std::string name_buf( tag_name, tag_name_len );
  tag_name = name_buf.c_str();
  
  iBase_ErrorType retval = TM->createTag(tag_name, this_size, tag_type, NULL, &new_tag);
  *tag_handle = (iBase_TagHandle) new_tag;
  RETURN(retval);
}

/**
 *   Delete a tag handle and the data associated with that tag.  The
 *   deletion can be forced or not forced.  If the deletion is forced,
 *   the tag and all of its associated data are deleted from the
 *   implementation even if the tag is still associated with mesh
 *   gentities.  If the deletion is not forced, the tag will not be
 *   deleted if it is still associated with one or more mesh gentities.
 *   In this case an error is returned asking the user to remove the
 *   tag from that gentity before deleting it.  If the underlying
 *   implementation does not support the requested deletion mechanism,
 *   an error will be returned.
 */
void
iGeom_destroyTag (iGeom_Instance instance,
                  /*in*/ iBase_TagHandle tag_handle,
                  /*in*/ int forced,
                  int* err)
{
  iBase_ErrorType retval = TM->destroyTag(TAG_HANDLE(tag_handle), forced);
  RETURN(retval);
}

/**
 *   Get the tag name associated with a given tag handle.
 */
void
iGeom_getTagName (iGeom_Instance instance,
                  /*in*/ iBase_TagHandle tag_handle,
                  char* name,
                  int* err,
                  int name_len)
{
  const char* name_ptr = TM->getTagName(TAG_HANDLE(tag_handle));
  int len = strlen(name_ptr) + 1;
  if (len > name_len)
    len = name_len;
  memcpy( name, name_ptr, len );
  RETURN(iBase_SUCCESS);
}

/**
 *   Get the size of the data associated with a given tag handle.
 */
void
iGeom_getTagSizeValues (iGeom_Instance instance,
                        /*in*/ iBase_TagHandle tag_handle,
                        int* tag_size,
                        int* err)
{
  *tag_size = TM->getTagSize(TAG_HANDLE(tag_handle));
  int type;
  iGeom_getTagType( instance, tag_handle, &type, err );
  if (iBase_SUCCESS != *err) return;
  
  switch (type) {
    case iBase_INTEGER:
      *tag_size /= sizeof(int);
      break;
    case iBase_DOUBLE:
      *tag_size /= sizeof(double);
      break;
    case iBase_ENTITY_HANDLE:
      *tag_size /= sizeof(iBase_EntityHandle);
      break;
    case iBase_BYTES:
      break;
  }
  RETURN(iBase_SUCCESS);
}

/**
 *   Get the size of the data associated with a given tag handle.
 */
void
iGeom_getTagSizeBytes (iGeom_Instance instance,
                       /*in*/ iBase_TagHandle tag_handle,
                       int* tag_size,
                       int* err)
{
  *tag_size = TM->getTagSize(TAG_HANDLE(tag_handle));
  RETURN(iBase_SUCCESS);
}

void
iGeom_getTagType (iGeom_Instance instance,
                  /*in*/ iBase_TagHandle tag_handle,
                  int* type,
                  int* err)
{
  *type = TM->getTagType(TAG_HANDLE(tag_handle));
  RETURN(iBase_SUCCESS);
}

/**
 *     Get the tag handle associated with a given string name.
 */
void
iGeom_getTagHandle (iGeom_Instance instance,
                    /*in*/ const char *tag_name,
                    iBase_TagHandle* tag_handle,
                    int* err,
                    int tag_name_len)
{
    // make sure string is null-terminated
  std::string tag_name_buf( tag_name, tag_name_len );
  tag_name = tag_name_buf.c_str();
  *tag_handle = reinterpret_cast<iBase_TagHandle>(static_cast<size_t>(TM->getTagHandle( tag_name )));

  *err = iGeom_getLastErrorType();
}

void
iGeom_rmvArrTag (iGeom_Instance instance,
                 const iBase_EntityHandle *entity_handles,
                 int entity_handles_size,
                 iBase_TagHandle tag_handle,
                 int* err) 
{
  RefEntity *const *tmp_entity_handles = reinterpret_cast<RefEntity*const *>(entity_handles);  
  iBase_ErrorType retval = TM->rmvArrTag(tmp_entity_handles, entity_handles_size, TAG_HANDLE(tag_handle));
  RETURN(retval);
}

/**
 *   Allows the user to disassociate the tag referenced by the tag
 *   handle from the specified gentities. The tag data is not deleted in
 *   this call, but can be deleted later using the deleteTag function
 *   defined above.
 */
void
iGeom_rmvTag (iGeom_Instance instance,
              /*in*/ iBase_EntityHandle entity_handle,
              /*in*/ iBase_TagHandle tag_handle,
              int* err)
{
  RefEntity *tmp_entity = ENTITY_HANDLE(entity_handle);
  iBase_ErrorType retval = TM->rmvArrTag(&tmp_entity, 1, TAG_HANDLE(tag_handle));
  RETURN(retval);
}

/**
 *   Get all tag handles associated with a given gentity.
 */
void
iGeom_getAllTags (iGeom_Instance instance,
                  /*in*/ iBase_EntityHandle entity_handle,
                  /*inout*/ iBase_TagHandle **tag_handles,
                  int *tag_handles_allocated,
                  int *tag_handles_size,
                  int* err)
{
  iBase_ErrorType retval = TM->getAllTags(SET_HANDLE(entity_handle), 
                                          TAG_HANDLE_ARRAY_INOUT(tag_handles));
  RETURN(retval);
}

void
iGeom_getArrData (iGeom_Instance instance,
                  /*in*/ iBase_EntityHandle const *entity_handles,
                  int entity_handles_size,
                  /*in*/ iBase_TagHandle tag_handle,
                  /*inout*/ char **tag_value,
                  int *tag_value_allocated,
                  int *tag_value_size,
                  int* err)
{
  iBase_ErrorType retval = TM->getArrData(reinterpret_cast<RefEntity*const*>(entity_handles),
                                          entity_handles_size,
                                          TAG_HANDLE(tag_handle), 
                                          ARRAY_INOUT(tag_value));
  RETURN(retval);
}


void
iGeom_getIntArrData (iGeom_Instance instance,
                     /*in*/ iBase_EntityHandle const *entity_handles,
                     int entity_handles_size,
                     /*in*/ iBase_TagHandle tag_handle,
                     /*inout*/ int **tag_value,
                     int *tag_value_allocated,
                     int *tag_value_size,
                     int* err)
{
  int tag_value_allocated_tmp = *tag_value_allocated * sizeof(int);
  int tag_value_size_tmp = *tag_value_size * sizeof(int);
  iGeom_getArrData(instance, entity_handles, 
                   entity_handles_size, tag_handle,
                   reinterpret_cast<char**>(tag_value), 
                   &tag_value_allocated_tmp, 
                   &tag_value_size_tmp,
                   err);
  *tag_value_allocated = tag_value_allocated_tmp / sizeof(int);
  *tag_value_size = tag_value_size_tmp / sizeof(int);
}


void
iGeom_getDblArrData (iGeom_Instance instance,
                     /*in*/ iBase_EntityHandle const *entity_handles,
                     int entity_handles_size,
                     /*in*/ iBase_TagHandle tag_handle,
                     /*inout*/ double **tag_value,
                     int *tag_value_allocated,
                     int *tag_value_size,
                     int* err)
{
  int tag_value_allocated_tmp = *tag_value_allocated * sizeof(double);
  int tag_value_size_tmp = *tag_value_size * sizeof(double);
  iGeom_getArrData(instance, entity_handles, 
                   entity_handles_size, tag_handle,
                   reinterpret_cast<char**>(tag_value), 
                   &tag_value_allocated_tmp, 
                   &tag_value_size_tmp,
                   err);
  *tag_value_allocated = tag_value_allocated_tmp / sizeof(double);
  *tag_value_size = tag_value_size_tmp / sizeof(double);
}

void
iGeom_getEHArrData (iGeom_Instance instance,
                    /*in*/ iBase_EntityHandle const *entity_handles,
                    int entity_handles_size,
                    /*in*/ iBase_TagHandle tag_handle,
                    /*inout*/ iBase_EntityHandle **tag_value,
                    int *tag_value_allocated,
                    int *tag_value_size,
                    int* err)
{
  int tag_value_allocated_tmp = *tag_value_allocated * sizeof(iBase_EntityHandle);
  int tag_value_size_tmp = *tag_value_size * sizeof(iBase_EntityHandle);
  iGeom_getArrData(instance, entity_handles, 
                   entity_handles_size, tag_handle,
                   reinterpret_cast<char**>(tag_value), 
                   &tag_value_allocated_tmp, 
                   &tag_value_size_tmp,
                   err);
  *tag_value_allocated = tag_value_allocated_tmp / sizeof(iBase_EntityHandle);
  *tag_value_size = tag_value_size_tmp / sizeof(iBase_EntityHandle);
}

void
iGeom_setArrData (iGeom_Instance instance,
                  const iBase_EntityHandle *entity_handles,
                  int entity_handles_size,
                  /*in*/ iBase_TagHandle tag_handle,
                  /*in*/ const char *tag_value,
                  const int tag_value_size,
                  int* err)
{
  iBase_ErrorType retval = TM->setArrData(ENTITY_HANDLE_CONST_ARRAY(entity_handles),
                                          entity_handles_size,
                                          TAG_HANDLE(tag_handle), 
                                          ARRAY_IN(tag_value));
  RETURN(retval);
}

void
iGeom_setIntArrData (iGeom_Instance instance,
                     const iBase_EntityHandle *entity_handles,
                     int entity_handles_size,
                     iBase_TagHandle tag_handle,
                     const int *tag_value,
                     int tag_value_size,
                     int* err)
{
  iGeom_setArrData(instance, entity_handles, 
                   entity_handles_size, tag_handle, 
                   reinterpret_cast<const char*>(tag_value), 
                   sizeof(int)*tag_value_size,
                   err);
}

void
iGeom_setDblArrData (iGeom_Instance instance,
                     const iBase_EntityHandle *entity_handles,
                     int entity_handles_size,
                     iBase_TagHandle tag_handle,
                     const double *tag_value,
                     int tag_value_size,
                     int* err)
{
  iGeom_setArrData(instance, entity_handles, 
                   entity_handles_size, tag_handle, 
                   reinterpret_cast<const char*>(tag_value), 
                   sizeof(double)*tag_value_size,
                   err);
}

void
iGeom_setEHArrData (iGeom_Instance instance,
                    const iBase_EntityHandle *entity_handles,
                    int entity_handles_size,
                    iBase_TagHandle tag_handle,
                    iBase_EntityHandle const *tag_values,
                    int tag_values_size,
                    int* err)
{
  iGeom_setArrData(instance, entity_handles, 
                   entity_handles_size, tag_handle, 
                   reinterpret_cast<const char*>(tag_values), 
                   sizeof(iBase_EntityHandle)*tag_values_size,
                   err);
}

/**
 *   Allows the user to retrieve an array of tag values associated with
 *   a tag handle from an input array of gentity handles
 */
void
iGeom_getData (iGeom_Instance instance,
               /*in*/ iBase_EntityHandle entity_handle,
               /*in*/ iBase_TagHandle tag_handle,
               /*inout*/ char **tag_value,
               int *tag_value_allocated,
               int *tag_value_size,
               int* err)
{
  RefEntity *tmp_entity = ENTITY_HANDLE(entity_handle);
  iBase_ErrorType retval = TM->getArrData(&tmp_entity, 1, TAG_HANDLE(tag_handle), 
                                          ARRAY_INOUT(tag_value));
  RETURN(retval);
}

void
iGeom_getIntData (iGeom_Instance instance,
                  iBase_EntityHandle entity_handle,
                  iBase_TagHandle tag_handle,
                  int* data_out,
                  int* err ) 
{
  char *val_ptr = reinterpret_cast<char*>(data_out);
  int val_size = sizeof(int);
  iGeom_getArrData(instance, &entity_handle, 1, 
                   tag_handle, &val_ptr, &val_size, &val_size, err);
}

void
iGeom_getDblData (iGeom_Instance instance,
                  iBase_EntityHandle entity_handle,
                  iBase_TagHandle tag_handle,
                  double* data_out,
                  int* err ) 
{
  char *val_ptr = reinterpret_cast<char*>(data_out);
  int val_size = sizeof(double);
  iGeom_getArrData(instance, &entity_handle, 1, 
                   tag_handle, &val_ptr, &val_size, &val_size, err);
}

void
iGeom_getEHData (iGeom_Instance instance,
                 iBase_EntityHandle entity_handle,
                 iBase_TagHandle tag_handle,
                 iBase_EntityHandle* data_out,
                 int* err ) 
{
  char *val_ptr = reinterpret_cast<char*>(data_out);
  int val_size = sizeof(iBase_EntityHandle);
  iGeom_getArrData(instance, &entity_handle, 1, 
                   tag_handle, &val_ptr, &val_size, &val_size, err);
}

/**
 *   Allows the user to set the tag data values on an array of gentity
 *   handles
 */
void
iGeom_setData (iGeom_Instance instance,
               /*in*/ iBase_EntityHandle entity_handle,
               /*in*/ iBase_TagHandle tag_handle,
               /*in*/ const char *tag_value,
               int tag_value_size,
               int* err)
{
  RefEntity *tmp_entity = ENTITY_HANDLE(entity_handle);
  iBase_ErrorType retval = TM->setArrData(&tmp_entity, 1, 
                                          TAG_HANDLE(tag_handle), 
                                          ARRAY_IN(tag_value));
  RETURN(retval);
}

void
iGeom_setIntData (iGeom_Instance instance,
                  /*in*/ iBase_EntityHandle entity_handle,
                  /*in*/ iBase_TagHandle tag_handle,
                  /*in*/ int tag_value,
                  int* err ) 
{
  iGeom_setArrData(instance, &entity_handle, 1, 
                   tag_handle, 
                   reinterpret_cast<const char*>(&tag_value), 
                   sizeof(int),
                   err);
}

void
iGeom_setDblData (iGeom_Instance instance,
                  /*in*/ iBase_EntityHandle entity_handle,
                  /*in*/ iBase_TagHandle tag_handle,
                  /*in*/ double tag_value,
                  int* err ) 
{
  iGeom_setArrData(instance, &entity_handle, 1, 
                   tag_handle, 
                   reinterpret_cast<const char*>(&tag_value), 
                   sizeof(double),
                   err);
}

void
iGeom_setEHData (iGeom_Instance instance,
                 /*in*/ iBase_EntityHandle entity_handle,
                 /*in*/ iBase_TagHandle tag_handle,
                 /*in*/ iBase_EntityHandle tag_value,
                 int* err ) 
{
  iGeom_setArrData(instance, &entity_handle, 1, 
                   tag_handle, reinterpret_cast<const char*>(&tag_value), 
                   sizeof(iBase_EntityHandle), err);
}

/**
 *   Remove the tag associated with the tag_handle from the mesh or
 *   gentity_set.  The tag data is not destroyed in this function, but
 *   can be destroyed using the deleteTag function.
 */
void
iGeom_rmvEntSetTag (iGeom_Instance instance,
                    /*in*/ iBase_EntitySetHandle entity_set,
                    /*in*/ iBase_TagHandle tag_handle,
                    int* err)
{
    // have to go through RefEntity* so that RefEntity** gets set right
  RefEntity *tmp_entity = SET_HANDLE(entity_set);
  iBase_ErrorType retval = TM->rmvArrTag(&tmp_entity, 1, TAG_HANDLE(tag_handle));
  RETURN(retval);
}

/**
 *   Get all tag handles associated with a given mesh or gentity_set.
 */
void
iGeom_getAllEntSetTags (iGeom_Instance instance,
                        /*in*/ iBase_EntitySetHandle entity_set,
                        /*inout*/ iBase_TagHandle **tag_handles,
                        int *tag_handles_allocated,
                        int *tag_handles_size,
                        int* err)
{
    // have to go through RefEntity* so that RefEntity** gets set right
  const RefEntity *tmp_entity = SET_HANDLE(entity_set);
  iBase_ErrorType retval = TM->getAllTags(tmp_entity,TAG_HANDLE_ARRAY_INOUT(tag_handles));
  RETURN(retval);
}

/**
 *   Get the tag data associated with a tag handle from the mesh or
 *   gentity_set.  It is assumed that the tag_value argument is
 *   allocated by the application before being passed into the getTag
 *   function.
 */
void
iGeom_getEntSetData (iGeom_Instance instance,
                     /*in*/ iBase_EntitySetHandle entity_set,
                     /*in*/ iBase_TagHandle tag_handle,
                     /*inout*/ char **tag_value,
                     int *tag_value_allocated,
                     int *tag_value_size,
                     int* err)
{
    // have to go through RefEntity* so that RefEntity** gets set right
  RefEntity *tmp_entity = SET_HANDLE(entity_set);
  iBase_ErrorType retval = TM->getArrData(&tmp_entity, 1, TAG_HANDLE(tag_handle), 
                                          ARRAY_INOUT(tag_value));
  RETURN(retval);
}

void
iGeom_getEntSetIntData (iGeom_Instance instance,
                        /*in*/ iBase_EntitySetHandle entity_set,
                        /*in*/ iBase_TagHandle tag_handle,
                        int* tag_ptr,
                        int* err ) 
{
  int tag_size = sizeof(int);
  char* data_ptr = reinterpret_cast<char*>(tag_ptr);
  iGeom_getEntSetData(instance, entity_set, tag_handle, &data_ptr, 
                      &tag_size, &tag_size, err);
}

void
iGeom_getEntSetDblData (iGeom_Instance instance,
                        /*in*/ iBase_EntitySetHandle entity_set,
                        /*in*/ iBase_TagHandle tag_handle,
                        double* tag_ptr,
                        int* err ) 
{
  int tag_size = sizeof(double);
  char* data_ptr = reinterpret_cast<char*>(tag_ptr);
  iGeom_getEntSetData(instance, entity_set, tag_handle, &data_ptr, 
                      &tag_size, &tag_size, err);
}
void
iGeom_getEntSetEHData (iGeom_Instance instance,
                       /*in*/ iBase_EntitySetHandle entity_set,
                       /*in*/ iBase_TagHandle tag_handle,
                       iBase_EntityHandle* tag_ptr,
                       int* err ) 
{
  int tag_size = sizeof(iBase_EntityHandle);
  char* data_ptr = reinterpret_cast<char*>(tag_ptr);
  iGeom_getEntSetData(instance, entity_set, tag_handle, &data_ptr, 
                      &tag_size, &tag_size, err);
}

/**
 *   Set the tag data associated with a given tag handle on the mesh or
 *   gentity_set
 */
void
iGeom_setEntSetData (iGeom_Instance instance,
                     /*in*/ iBase_EntitySetHandle entity_set,
                     /*in*/ iBase_TagHandle tag_handle,
                     /*in*/ const char *tag_value,
                     int tag_value_size,
                     int* err)
{
    // have to go through RefEntity* so that RefEntity** gets set right
  RefEntity *tmp_entity = SET_HANDLE(entity_set);
  iBase_ErrorType retval = TM->setArrData(&tmp_entity, 1, TAG_HANDLE(tag_handle), 
                                          ARRAY_IN(tag_value));
  RETURN(retval);
}

void
iGeom_setEntSetIntData (iGeom_Instance instance,
                        /*in*/ iBase_EntitySetHandle entity_set,
                        /*in*/ iBase_TagHandle tag_handle,
                        /*in*/ int tag_value,
                        int* err ) 
{
  iGeom_setEntSetData(instance, entity_set, tag_handle, 
                      reinterpret_cast<const char*>(&tag_value), 
                      sizeof(int),
                      err);
}

void
iGeom_setEntSetDblData (iGeom_Instance instance,
                        /*in*/ iBase_EntitySetHandle entity_set,
                        /*in*/ iBase_TagHandle tag_handle,
                        /*in*/ double tag_value,
                        int* err ) 
{
  iGeom_setEntSetData(instance, entity_set, tag_handle, 
                      reinterpret_cast<const char*>(&tag_value),
                      sizeof(double), err);
}

void
iGeom_setEntSetEHData (iGeom_Instance instance,
                       /*in*/ iBase_EntitySetHandle entity_set,
                       /*in*/ iBase_TagHandle tag_handle,
                       /*in*/ iBase_EntityHandle tag_value,
                       int* err ) 
{
  iGeom_setEntSetData(instance, entity_set, tag_handle, 
                      reinterpret_cast<const char*>(&tag_value), 
                      sizeof(iBase_EntityHandle), err);
}

void
iGeom_getEntClosestPt(iGeom_Instance instance,
                      /*in*/ iBase_EntityHandle entity_handle,
                      /*in*/ double near_x,
                      /*in*/ double near_y,
                      /*in*/ double near_z,
                      /*out*/ double* on_x,
                      /*out*/ double* on_y, 
                      /*out*/ double* on_z,
                      int* err)
{
  RefEntity* entity = (RefEntity*)entity_handle;
  CubitVector on, near(near_x, near_y, near_z);
  CubitStatus status = iGeom_closest_point( entity, near, on );
  if (status == CUBIT_FAILURE) {
    ERROR(iBase_FAILURE, "Problems getting closest point for some entity.");
  }
  on.get_xyz( *on_x, *on_y, *on_z );
  RETURN(iBase_SUCCESS);
}

/**
 * Return a points on specified entities closest to specified points
 * in space.  Input coordinates and output points are interleaved in 
 * the arrays.
 * @param gentity_handles The gentities being queried
 * @param near_coordinates Input coordinates
 * @param on_coordinates Closest point on gentity
 */
void
iGeom_getArrClosestPt (iGeom_Instance instance,
                       /*in*/ iBase_EntityHandle const *gentity_handles,
                       int gentity_handles_size,
                       /*in*/ int storage_order,
                       /*in*/ const double *near_coordinates,
                       int near_coordinates_size,
                       /*out*/ double **on_coordinates,
                       int *on_coordinates_allocated,
                       int *on_coordinates_size,
                       int* err)
{
    /* set up iteration according to storage order.
       allow either gentity_handles or near_coordinates to contain
       only one value, where that single value is applied for every
       entry in the other list.
    */
  size_t near_step, on_step = 1, ent_step;
  int count;
  if (3*gentity_handles_size == near_coordinates_size) {
    near_step = ent_step = 1;
    count = gentity_handles_size;
  }
  else if (near_coordinates_size == 3) {
    near_step = 0;
    ent_step = 1;
    count = gentity_handles_size;
  }
  else if (gentity_handles_size == 1) {
    near_step = 1;
    ent_step = 0;
    count = near_coordinates_size / 3;
  }
  else {
    ERROR( iBase_INVALID_ENTITY_COUNT, "Mismatched array sizes" );
  }
  ALLOC_CHECK_ARRAY( on_coordinates, 3*count );
  
  const double *near_x, *near_y, *near_z;
  double *on_x, *on_y, *on_z;
  if (storage_order == iBase_BLOCKED) {
    near_x = near_coordinates;
    near_y = near_x + near_coordinates_size/3;
    near_z = near_y + near_coordinates_size/3;
    on_x = *on_coordinates;
    on_y = on_x + count;
    on_z = on_y + count;
    on_step = 1;
  }
  else {
    storage_order = iBase_INTERLEAVED; /* set if unspecified */
    near_x = near_coordinates;
    near_y = near_x+1;
    near_z = near_x+2;
    on_x = *on_coordinates;
    on_y = on_x+1;
    on_z = on_x+2;
    near_step *= 3;
    on_step = 3;
  }
  
  RefEntity** ent_iter = (RefEntity**)(gentity_handles);
  CubitStatus final_result = CUBIT_SUCCESS;
  for (int i = 0; i < count; ++i) {
    CubitVector on, near( *near_x, *near_y, *near_z );
    CubitStatus status = iGeom_closest_point( *ent_iter, near, on );
    if (status != CUBIT_SUCCESS)
      final_result = status;
    on.get_xyz( *on_x, *on_y, *on_z );
    
    ent_iter += ent_step;
    near_x += near_step;
    near_y += near_step;
    near_z += near_step;
    on_x += on_step;
    on_y += on_step;
    on_z += on_step;
  }
  
  if (final_result == CUBIT_FAILURE) {
    ERROR(iBase_FAILURE, "Problems getting closest point for some entity.");
  }

  KEEP_ARRAY(on_coordinates);
  RETURN(iBase_SUCCESS);
}


void
iGeom_getEntNrmlPlXYZ(iGeom_Instance instance,
                      /*in*/ iBase_EntityHandle entity_handle,
                      /*in*/ double x,
                      /*in*/ double y,
                      /*in*/ double z,
                      /*out*/ double* pt_x,
                      /*out*/ double* pt_y, 
                      /*out*/ double* pt_z,
                      /*out*/ double* nmrl_i,
                      /*out*/ double* nmrl_j, 
                      /*out*/ double* nmrl_k,
                      int* err)
{
  RefEntity* entity = (RefEntity*)entity_handle;
  CubitVector pt, nmrl, near(x, y, z);
  CubitStatus status = iGeom_closest_point_and_normal( entity, near, pt, nmrl );
  if (status == CUBIT_FAILURE) {
    ERROR(iBase_FAILURE, "Problems getting closest point for some entity.");
  }
  pt.get_xyz( *pt_x, *pt_y, *pt_z );
  nmrl.get_xyz( *nmrl_i, *nmrl_j, *nmrl_k );
  RETURN(iBase_SUCCESS);
}

/**
 * Return a points on specified entities closest to specified points
 * in space, and normals at that point.  Input coordinates and output points are 
 * interleaved in 
 * the arrays.
 * @param gentity_handles The gentities being queried
 * @param near_coordinates Input coordinates
 * @param on_coordinates Closest point on gentity
 */
void
iGeom_getArrNrmlPlXYZ (iGeom_Instance instance,
                       /*in*/ iBase_EntityHandle const *gentity_handles,
                       const int gentity_handles_size,
                       /*in*/ int storage_order,
                       /*in*/ const double *near_coordinates,
                       const int near_coordinates_size,
                       /*out*/ double **on_coordinates,
                       int *on_coordinates_allocated,
                       int *on_coordinates_size,
                       /*out*/ double **normals,
                       int *normals_allocated,
                       int *normals_size,
                       int* err)
{
    /* set up iteration according to storage order.
       allow either gentity_handles or near_coordinates to contain
       only one value, where that single value is applied for every
       entry in the other list.
    */
  size_t near_step, on_step = 1, ent_step;
  int count;
  if (3*gentity_handles_size == near_coordinates_size) {
    near_step = ent_step = 1;
    count = gentity_handles_size;
  }
  else if (near_coordinates_size == 3) {
    near_step = 0;
    ent_step = 1;
    count = gentity_handles_size;
  }
  else if (gentity_handles_size == 1) {
    near_step = 1;
    ent_step = 0;
    count = near_coordinates_size / 3;
  }
  else {
    ERROR( iBase_INVALID_ENTITY_COUNT, "Mismatched array sizes" );
  }
  ALLOC_CHECK_ARRAY( on_coordinates, 3*count );
  ALLOC_CHECK_ARRAY( normals, 3*count );
  
  const double *near_x, *near_y, *near_z;
  double *on_x, *on_y, *on_z;
  double *norm_x, *norm_y, *norm_z;
  if (storage_order == iBase_BLOCKED) {
    near_x = near_coordinates;
    near_y = near_x + near_coordinates_size/3;
    near_z = near_y + near_coordinates_size/3;
    on_x = *on_coordinates;
    on_y = on_x + count;
    on_z = on_y + count;
    norm_x = *normals;
    norm_y = norm_x + count;
    norm_z = norm_y + count;
    on_step = 1;
  }
  else {
    storage_order = iBase_INTERLEAVED; /* set if unspecified */
    near_x = near_coordinates;
    near_y = near_x+1;
    near_z = near_x+2;
    on_x = *on_coordinates;
    on_y = on_x+1;
    on_z = on_x+2;
    norm_x = *normals;
    norm_y = norm_x+1;
    norm_z = norm_x+2;
    near_step *= 3;
    on_step = 3;
  }
  
  RefEntity** ent_iter = (RefEntity**)(gentity_handles);
  CubitStatus final_result = CUBIT_SUCCESS;
  for (int i = 0; i < count; ++i) {
    CubitVector on, norm, near( *near_x, *near_y, *near_z );
    CubitStatus status = iGeom_closest_point_and_normal( *ent_iter, near, on, norm );
    if (status != CUBIT_SUCCESS)
      final_result = status;
    on.get_xyz( *on_x, *on_y, *on_z );
    norm.get_xyz( *norm_x, *norm_y, *norm_z );
    
    ent_iter += ent_step;
    near_x += near_step;
    near_y += near_step;
    near_z += near_step;
    on_x += on_step;
    on_y += on_step;
    on_z += on_step;
    norm_x += on_step;
    norm_y += on_step;
    norm_z += on_step;
  }
  
  if (final_result == CUBIT_FAILURE) {
    ERROR(iBase_FAILURE, "Problems getting closest point for some entity.");
  }

  KEEP_ARRAY(on_coordinates);
  KEEP_ARRAY(normals);
  RETURN(iBase_SUCCESS);
}

void
iGeom_getEntNrmlXYZ (iGeom_Instance instance,
                     /*in*/ iBase_EntityHandle entity_handle,
                     /*in*/ double near_x,
                     /*in*/ double near_y,
                     /*in*/ double near_z,
                     /*out*/ double* nmrl_i,
                     /*out*/ double* nmrl_j, 
                     /*out*/ double* nmrl_k,
                     int* err)
{
  RefEntity* entity = (RefEntity*)entity_handle;
  RefFace* face = dynamic_cast<RefFace*>(entity);
  if (NULL == face) {
    ERROR(iBase_INVALID_ENTITY_TYPE, "Entities passed into gentityNormal must be faces.");
  }
  
  CubitVector normal, near( near_x, near_y, near_z );
  normal = face->normal_at( near );
  normal.get_xyz( *nmrl_i, *nmrl_j, *nmrl_k );
  RETURN(iBase_SUCCESS);
} 

/**
 * Return the normals at point on specified entities.  Returns error
 * if any input entity is not a gface.  Input coordinates and normals
 * are interleaved in the arrays.
 * @param gentity_handles The gentities being queried
 * @param coordinates Input coordinates, interleaved
 * @param normals The normals at the specified points, interleaved
 */
void
iGeom_getArrNrmlXYZ (iGeom_Instance instance,
                     /*in*/ iBase_EntityHandle const *gentity_handles,
                     int gentity_handles_size,
                     int storage_order,
                     /*in*/ const double *coordinates,
                     int coordinates_size,
                     /*out*/ double **normals,
                     int *normals_allocated,
                     int *normals_size,
                     int* err)
{
    /* set up iteration according to storage order.
       allow either gentity_handles or near_coordinates to contain
       only one value, where that single value is applied for every
       entry in the other list.
    */
  size_t coord_step, norm_step = 1, ent_step;
  int count;
  if (3*gentity_handles_size == coordinates_size) {
    coord_step = ent_step = 1;
    count = gentity_handles_size;
  }
  else if (coordinates_size == 3) {
    coord_step = 0;
    ent_step = 1;
    count = gentity_handles_size;
  }
  else if (gentity_handles_size == 1) {
    coord_step = 1;
    ent_step = 0;
    count = coordinates_size / 3;
  }
  else {
    ERROR( iBase_INVALID_ENTITY_COUNT, "Mismatched array sizes" );
  }

    // check or pre-allocate the coordinate arrays
  ALLOC_CHECK_ARRAY( normals, 3*count );
  
  const double *coord_x, *coord_y, *coord_z;
  double *norm_x, *norm_y, *norm_z;
  if (storage_order == iBase_BLOCKED) {
    coord_x = coordinates;
    coord_y = coord_x + coordinates_size/3;
    coord_z = coord_y + coordinates_size/3;
    norm_x = *normals;
    norm_y = norm_x + count;
    norm_z = norm_y + count;
    norm_step = 1;
  }
  else {
    storage_order = iBase_INTERLEAVED; /* set if unspecified */
    coord_x = coordinates;
    coord_y = coord_x+1;
    coord_z = coord_x+2;
    norm_x = *normals;
    norm_y = norm_x+1;
    norm_z = norm_x+2;
    coord_step *= 3;
    norm_step = 3;
  }
  
  RefEntity** entities = (RefEntity**)(gentity_handles);
  for (int i = 0; i < count; ++i)
  {
    RefFace* face = dynamic_cast<RefFace*>(*entities);
    if (NULL == face) {
      ERROR(iBase_INVALID_ENTITY_TYPE, "Entities passed into gentityNormal must be faces.");
    }
    else {
      CubitVector normal, coords( *coord_x, *coord_y, *coord_z );
      normal = face->normal_at( coords );
      normal.get_xyz( *norm_x, *norm_y, *norm_z );
    }
    entities += ent_step;
    coord_x += coord_step;
    coord_y += coord_step;
    coord_z += coord_step;
    norm_x += norm_step;
    norm_y += norm_step;
    norm_z += norm_step;
  }

  KEEP_ARRAY(normals);
  RETURN(iBase_SUCCESS);
}


void
iGeom_getEntTgntXYZ (iGeom_Instance instance,
                     /*in*/ iBase_EntityHandle entity_handle,
                     /*in*/ double x,
                     /*in*/ double y,
                     /*in*/ double z,
                     /*out*/ double* tgnt_i,
                     /*out*/ double* tgnt_j, 
                     /*out*/ double* tgnt_k,
                     int* err)
{
  RefEntity* entity = (RefEntity*)entity_handle;
  RefEdge* edge = dynamic_cast<RefEdge*>(entity);
  if (NULL == edge) {
    ERROR(iBase_INVALID_ENTITY_TYPE, "Entities passed into gentityTangent must be edges.");
  }
  
  CubitVector tangent, near( x, y, x );
  edge->tangent( near, tangent );
  tangent.get_xyz( *tgnt_i, *tgnt_j, *tgnt_k );
  RETURN(iBase_SUCCESS);
} 

/**
 * Return the tangent at point on specified entities.  Returns error
 * if any input entity is not a gedge.  Input coordinates and tangents
 * are interleaved in the arrays.
 * @param gentity_handles The gentities being queried
 * @param coordinates Input coordinates, interleaved
 * @param tangents The tangents at the specified points, interleaved
 */
void
iGeom_getArrTgntXYZ (iGeom_Instance instance,
                     /*in*/ iBase_EntityHandle const *gentity_handles,
                     int gentity_handles_size,
                     /*in*/ int storage_order,
                     /*in*/ const double *coordinates,
                     int coordinates_size,
                     /*out*/ double **tangents,
                     int *tangents_allocated,
                     int *tangents_size,
                     int* err)
{
    /* set up iteration according to storage order.
       allow either gentity_handles or near_coordinates to contain
       only one value, where that single value is applied for every
       entry in the other list.
    */
  size_t coord_step, tan_step = 1, ent_step;
  int count;
  if (3*gentity_handles_size == coordinates_size) {
    coord_step = ent_step = 1;
    count = gentity_handles_size;
  }
  else if (coordinates_size == 3) {
    coord_step = 0;
    ent_step = 1;
    count = gentity_handles_size;
  }
  else if (gentity_handles_size == 1) {
    coord_step = 1;
    ent_step = 0;
    count = coordinates_size / 3;
  }
  else {
    ERROR( iBase_INVALID_ENTITY_COUNT, "Mismatched array sizes" );
  }

    // check or pre-allocate the coordinate arrays
  ALLOC_CHECK_ARRAY( tangents, 3*count );
  
  const double *coord_x, *coord_y, *coord_z;
  double *tan_x, *tan_y, *tan_z;
  if (storage_order == iBase_BLOCKED) {
    coord_x = coordinates;
    coord_y = coord_x + coordinates_size/3;
    coord_z = coord_y + coordinates_size/3;
    tan_x = *tangents;
    tan_y = tan_x + count;
    tan_z = tan_y + count;
    tan_step = 1;
  }
  else {
    storage_order = iBase_INTERLEAVED; /* set if unspecified */
    coord_x = coordinates;
    coord_y = coord_x+1;
    coord_z = coord_x+2;
    tan_x = *tangents;
    tan_y = tan_x+1;
    tan_z = tan_x+2;
    coord_step *= 3;
    tan_step = 3;
  }
  
  RefEntity** entities = (RefEntity**)(gentity_handles);
  for (int i = 0; i < count; ++i)
  {
    RefEdge* edge = dynamic_cast<RefEdge*>(*entities);
    if (NULL == edge) {
      ERROR(iBase_INVALID_ENTITY_TYPE, "Entities passed into gentityTangent must be edges.");
    }
    else {
      CubitVector tangent, coords( *coord_x, *coord_y, *coord_z );
      edge->tangent( coords, tangent );
      tangent.get_xyz( *tan_x, *tan_y, *tan_z );
    }
    
    entities += ent_step;
    coord_x += coord_step;
    coord_y += coord_step;
    coord_z += coord_step;
    tan_x += tan_step;
    tan_y += tan_step;
    tan_z += tan_step;
  }

  KEEP_ARRAY(tangents);
  RETURN(iBase_SUCCESS);
}


void
iGeom_getFcCvtrXYZ(iGeom_Instance instance,
                   /*in*/ iBase_EntityHandle face_handle,
                   /*in*/ double x,
                   /*in*/ double y,
                   /*in*/ double z,
                   /*out*/ double* cvtr1_i,
                   /*out*/ double* cvtr1_j,
                   /*out*/ double* cvtr1_k,
                   /*out*/ double* cvtr2_i,
                   /*out*/ double* cvtr2_j,
                   /*out*/ double* cvtr2_k,
                   int* err)
{
  RefEntity* entity = (RefEntity*)face_handle;
  RefFace* face = dynamic_cast<RefFace*>(entity);
  if (!face)
    RETURN(iBase_INVALID_ENTITY_TYPE);
  
  Surface* surf = face->get_surface_ptr();
  CubitVector loc(x,y,z), curv1, curv2;
  CubitStatus status = surf->closest_point( loc, 0, 0, &curv1, &curv2 );
  if (surf->bridge_sense() == CUBIT_REVERSED) {
    curv1 = -curv1;
    curv2 = -curv2;
  }
  curv1.get_xyz( *cvtr1_i, *cvtr1_j, *cvtr1_k );
  curv2.get_xyz( *cvtr2_i, *cvtr2_j, *cvtr2_k );
  RETURN( (status == CUBIT_SUCCESS ? iBase_SUCCESS : iBase_FAILURE) );
}

void
iGeom_getEgCvtrXYZ(iGeom_Instance instance,
                   /*in*/ iBase_EntityHandle edge_handle,
                   /*in*/ double x,
                   /*in*/ double y,
                   /*in*/ double z,
                   /*out*/ double* cvtr_i,
                   /*out*/ double* cvtr_j,
                   /*out*/ double* cvtr_k,
                   int* err)
{
  RefEntity* entity = (RefEntity*)edge_handle;
  RefEdge* edge = dynamic_cast<RefEdge*>(entity);
  if (!edge)
    RETURN(iBase_INVALID_ENTITY_TYPE);
  
  CubitVector loc(x,y,z), junk, curv;
  CubitStatus status = edge->closest_point( loc, junk, 0, &curv );
  curv.get_xyz( *cvtr_i, *cvtr_j, *cvtr_k );
  RETURN( (status == CUBIT_SUCCESS ? iBase_SUCCESS : iBase_FAILURE) );
}

void
iGeom_getEntArrCvtrXYZ(iGeom_Instance instance,
                       /*in*/ const iBase_EntityHandle *entity_handles,
                       int entity_handles_size,
                       /*in*/ int storage_order,
                       /*in*/ const double *coordinates,
                       int coordinates_size,
                       /*inout*/ double **cvtr_1,
                       int *cvtr_1_allocated,
                       int *cvtr_1_size,
                       /*inout*/ double **cvtr_2,
                       int *cvtr_2_allocated,
                       int *cvtr_2_size,
                       int* err)
{
  RefEntity** entities = (RefEntity**)(entity_handles);
  
    /* set up iteration according to storage order.
       allow either gentity_handles or near_coordinates to contain
       only one value, where that single value is applied for every
       entry in the other list.
    */
  size_t coord_step, cvtr_step = 1, ent_step;
  int count;
  if (3*entity_handles_size == coordinates_size) {
    coord_step = ent_step = 1;
    count = entity_handles_size;
  }
  else if (coordinates_size == 3) {
    coord_step = 0;
    ent_step = 1;
    count = entity_handles_size;
  }
  else if (entity_handles_size == 1) {
    coord_step = 1;
    ent_step = 0;
    count = coordinates_size / 3;
  }
  else {
    ERROR( iBase_INVALID_ENTITY_COUNT, "Mismatched array sizes" );
  }
  
    /* check if input list contains any surfaces. */
  bool have_surfs = false;
  for (int s = 0; s < entity_handles_size; ++s)
    if (dynamic_cast<RefFace*>(entities[s])) {
      have_surfs = true;
      break;
    }

    // check or pre-allocate the coordinate arrays
  ALLOC_CHECK_ARRAY( cvtr_1, 3*count );
  ALLOC_CHECK_ARRAY( cvtr_2, have_surfs ? 3*count : 0 );
  
  const double *coord_x, *coord_y, *coord_z;
  double *c1x, *c1y, *c1z, *c2x, *c2y, *c2z;
  if (storage_order == iBase_BLOCKED) {
    coord_x = coordinates;
    coord_y = coord_x + coordinates_size/3;
    coord_z = coord_y + coordinates_size/3;
    c1x = *cvtr_1;
    c1y = c1x + count;
    c1z = c1y + count;
    c2x = *cvtr_2;
    c2y = c2x + count;
    c2z = c2y + count;
    cvtr_step = 1;
  }
  else {
    storage_order = iBase_INTERLEAVED; /* set if unspecified */
    coord_x = coordinates;
    coord_y = coord_x+1;
    coord_z = coord_x+2;
    c1x = *cvtr_1;
    c1y = c1x+1;
    c1z = c1x+2;
    c2x = *cvtr_2;
    c2y = c2x+1;
    c2z = c2x+2;
    coord_step *= 3;
    cvtr_step = 3;
  }
  
  RefFace *face;
  RefEdge *edge;
  CubitStatus result = CUBIT_SUCCESS;
  for (int i = 0; i < count; ++i) {
    CubitStatus status;
    const CubitVector coords( *coord_x, *coord_y, *coord_z );
    if ((face = dynamic_cast<RefFace*>(*entities))) {
      Surface* surf = face->get_surface_ptr();
      CubitVector curv1, curv2;
      status = surf->closest_point( coords, 0, 0, &curv1, &curv2 );
      curv1.get_xyz( *c1x, *c1y, *c1z );
      curv2.get_xyz( *c2x, *c2y, *c2z );
    }
    else if ((edge = dynamic_cast<RefEdge*>(*entities))) {
      CubitVector junk, curv;
      status = edge->closest_point( coords, junk, 0, &curv );
      curv.get_xyz( *c1x, *c1y, *c1z );
    }
    else {
      status = CUBIT_FAILURE;
    }

    if (CUBIT_SUCCESS != status)
      result = status;
    
    entities += ent_step;
    coord_x += coord_step;
    coord_y += coord_step;
    coord_z += coord_step;
    c1x += cvtr_step;
    c1y += cvtr_step;
    c1z += cvtr_step;
    c2x += cvtr_step;
    c2y += cvtr_step;
    c2z += cvtr_step;
  }
  
  if (result == CUBIT_FAILURE)
    RETURN(iBase_FAILURE);
  
  KEEP_ARRAY( cvtr_1 );
  KEEP_ARRAY( cvtr_2 );
  RETURN( iBase_SUCCESS );
}


void
iGeom_getEgEvalXYZ(iGeom_Instance instance,
                   /*in*/ iBase_EntityHandle edge_handle,
                   /*in*/ double x,
                   /*in*/ double y,
                   /*in*/ double z,
                   /*out*/ double* on_x,
                   /*out*/ double* on_y,
                   /*out*/ double* on_z,
                   /*out*/ double* tan_i,
                   /*out*/ double* tan_j,
                   /*out*/ double* tan_k,
                   /*out*/ double* cvtr_i,
                   /*out*/ double* cvtr_j,
                   /*out*/ double* cvtr_k,
                   int* err)
{
  RefEntity* entity = (RefEntity*)edge_handle;
  RefEdge* edge = dynamic_cast<RefEdge*>(entity);
  if (!edge)
    RETURN(iBase_INVALID_ENTITY_TYPE);
  
  CubitVector loc(x,y,z), on, tan, curv;
  CubitStatus status = edge->closest_point( loc, on, &tan, &curv );
  on.get_xyz( *on_x, *on_y, *on_z );
  tan.get_xyz( *tan_i, *tan_j, *tan_k );
  curv.get_xyz( *cvtr_i, *cvtr_j, *cvtr_k );
  RETURN( (status == CUBIT_SUCCESS ? iBase_SUCCESS : iBase_FAILURE) );
}

void
iGeom_getFcEvalXYZ(iGeom_Instance instance,
                   /*in*/ iBase_EntityHandle face_handle,
                   /*in*/ double x,
                   /*in*/ double y,
                   /*in*/ double z,
                   /*out*/ double* on_x,
                   /*out*/ double* on_y,
                   /*out*/ double* on_z,
                   /*out*/ double* norm_i,
                   /*out*/ double* norm_j,
                   /*out*/ double* norm_k,
                   /*out*/ double* cvtr1_i,
                   /*out*/ double* cvtr1_j,
                   /*out*/ double* cvtr1_k,
                   /*out*/ double* cvtr2_i,
                   /*out*/ double* cvtr2_j,
                   /*out*/ double* cvtr2_k,
                   int* err)
{
  RefEntity* entity = (RefEntity*)face_handle;
  RefFace* face = dynamic_cast<RefFace*>(entity);
  if (!face)
    RETURN(iBase_INVALID_ENTITY_TYPE);
  
  Surface* surf = face->get_surface_ptr();
  CubitVector loc(x,y,z), on, norm, curv1, curv2;
  CubitStatus status = surf->closest_point( loc, &on, &norm, &curv1, &curv2 );
  if (surf->bridge_sense() == CUBIT_REVERSED) {
    norm = -norm;
    curv1 = -curv1;
    curv2 = -curv2;
  }
  on.get_xyz( *on_x, *on_y, *on_z );
  norm.get_xyz( *norm_i, *norm_j, *norm_k );
  curv1.get_xyz( *cvtr1_i, *cvtr1_j, *cvtr1_k );
  curv2.get_xyz( *cvtr2_i, *cvtr2_j, *cvtr2_k );
  RETURN( (status == CUBIT_SUCCESS ? iBase_SUCCESS : iBase_FAILURE) );
}
                       
void
iGeom_getArrEgEvalXYZ (iGeom_Instance instance,
                       /*in*/ const iBase_EntityHandle *edge_handles,
                       int edge_handles_size,
                       /*in*/ int storage_order,
                       /*in*/ const double *coordinates,
                       int coordinates_size,
                       /*inout*/ double **on_coords,
                       int *on_coords_allocated,
                       int *on_coords_size,
                       /*inout*/ double **tangents,
                       int *tangents_allocated,
                       int *tangents_size,
                       /*inout*/ double **curvatures,
                       int *curvatures_allocated,
                       int *curvatures_size,
                       int* err)
{
    /* set up iteration according to storage order.
       allow either gentity_handles or near_coordinates to contain
       only one value, where that single value is applied for every
       entry in the other list.
    */
  size_t coord_step, on_step = 1, ent_step;
  int count;
  if (3*edge_handles_size == coordinates_size) {
    coord_step = ent_step = 1;
    count = edge_handles_size;
  }
  else if (coordinates_size == 3) {
    coord_step = 0;
    ent_step = 1;
    count = edge_handles_size;
  }
  else if (edge_handles_size == 1) {
    coord_step = 1;
    ent_step = 0;
    count = coordinates_size / 3;
  }
  else {
    ERROR( iBase_INVALID_ENTITY_COUNT, "Mismatched array sizes" );
  }

    // check or pre-allocate the coordinate arrays
  ALLOC_CHECK_ARRAY( on_coords, 3*count );
  ALLOC_CHECK_ARRAY( tangents, 3*count );
  ALLOC_CHECK_ARRAY( curvatures, 3*count );
  
  const double *coord_x, *coord_y, *coord_z;
  double *on_x, *on_y, *on_z;
  double *tan_x, *tan_y, *tan_z;
  double *curv_x, *curv_y, *curv_z;
  if (storage_order == iBase_BLOCKED) {
    coord_x = coordinates;
    coord_y = coord_x + coordinates_size/3;
    coord_z = coord_y + coordinates_size/3;
    on_x = *on_coords;
    on_y = on_x + count;
    on_z = on_y + count;
    tan_x = *tangents;
    tan_y = tan_x + count;
    tan_z = tan_y + count;
    curv_x = *curvatures;
    curv_y = curv_x + count;
    curv_z = curv_y + count;
    on_step = 1;
  }
  else {
    storage_order = iBase_INTERLEAVED; /* set if unspecified */
    coord_x = coordinates;
    coord_y = coord_x+1;
    coord_z = coord_x+2;
    on_x = *on_coords;
    on_y = on_x + 1;
    on_z = on_x + 2;
    tan_x = *tangents;
    tan_y = tan_x+1;
    tan_z = tan_x+2;
    curv_x = *curvatures;
    curv_y = curv_x + 1;
    curv_z = curv_x + 2;
    coord_step *= 3;
    on_step = 3;
  }
  
  RefEntity** ent_iter = (RefEntity**)edge_handles;
  iBase_ErrorType result = iBase_SUCCESS;
  for (int i = 0; i < count; ++i) {
    RefEdge* edge = dynamic_cast<RefEdge*>(*ent_iter);
    if (!edge) {
      ERROR(iBase_INVALID_ENTITY_TYPE, "Non-edge input handle.");
    }
    else {
      const CubitVector coords( *coord_x, *coord_y, *coord_z );
      CubitVector on, tan, curv;
      CubitStatus s = edge->closest_point( coords, on, &tan, &curv );
      if (s == CUBIT_FAILURE)
        result = iBase_FAILURE;
      on.get_xyz( *on_x, *on_y, *on_z );
      tan.get_xyz( *tan_x, *tan_y, *tan_z );
      curv.get_xyz( *curv_x, *curv_y, *curv_z );
    }
    
    ent_iter += ent_step;
    coord_x += coord_step;
    coord_y += coord_step;
    coord_z += coord_step;
    on_x += on_step;
    on_y += on_step;
    on_z += on_step;
    tan_x += on_step;
    tan_y += on_step;
    tan_z += on_step;
    curv_x += on_step;
    curv_y += on_step;
    curv_z += on_step;
  }
  
  KEEP_ARRAY( on_coords );
  KEEP_ARRAY( tangents );
  KEEP_ARRAY( curvatures );
  RETURN(result);
}
                       
void
iGeom_getArrFcEvalXYZ (iGeom_Instance instance,
                       /*in*/ const iBase_EntityHandle *face_handles,
                       int face_handles_size,
                       /*in*/ int storage_order,
                       /*in*/ const double *coordinates,
                       int coordinates_size,
                       /*inout*/ double **on_coords,
                       int *on_coords_allocated,
                       int *on_coords_size,
                       /*inout*/ double **normals,
                       int *normals_allocated,
                       int *normals_size,
                       /*inout*/ double **curvatures_1,
                       int *curvatures_1_allocated,
                       int *curvatures_1_size,
                       /*inout*/ double **curvatures_2,
                       int *curvatures_2_allocated,
                       int *curvatures_2_size,
                       int* err)
{
    /* set up iteration according to storage order.
       allow either gentity_handles or near_coordinates to contain
       only one value, where that single value is applied for every
       entry in the other list.
    */
  size_t coord_step, on_step = 1, ent_step;
  int count;
  if (3*face_handles_size == coordinates_size) {
    coord_step = ent_step = 1;
    count = face_handles_size;
  }
  else if (coordinates_size == 3) {
    coord_step = 0;
    ent_step = 1;
    count = face_handles_size;
  }
  else if (face_handles_size == 1) {
    coord_step = 1;
    ent_step = 0;
    count = coordinates_size / 3;
  }
  else {
    ERROR( iBase_INVALID_ENTITY_COUNT, "Mismatched array sizes" );
  }

    // check or pre-allocate the coordinate arrays
  ALLOC_CHECK_ARRAY( on_coords, 3*count );
  ALLOC_CHECK_ARRAY( normals, 3*count );
  ALLOC_CHECK_ARRAY( curvatures_1, 3*count );
  ALLOC_CHECK_ARRAY( curvatures_2, 3*count );
  
  const double *coord_x, *coord_y, *coord_z;
  double *on_x, *on_y, *on_z;
  double *norm_x, *norm_y, *norm_z;
  double *curv1_x, *curv1_y, *curv1_z;
  double *curv2_x, *curv2_y, *curv2_z;
  if (storage_order == iBase_BLOCKED) {
    coord_x = coordinates;
    coord_y = coord_x + coordinates_size/3;
    coord_z = coord_y + coordinates_size/3;
    on_x = *on_coords;
    on_y = on_x + count;
    on_z = on_y + count;
    norm_x = *normals;
    norm_y = norm_x + count;
    norm_z = norm_y + count;
    curv1_x = *curvatures_1;
    curv1_y = curv1_x + count;
    curv1_z = curv1_y + count;
    curv2_x = *curvatures_2;
    curv2_y = curv2_x + count;
    curv2_z = curv2_y + count;
    on_step = 1;
  }
  else {
    storage_order = iBase_INTERLEAVED; /* set if unspecified */
    coord_x = coordinates;
    coord_y = coord_x+1;
    coord_z = coord_x+2;
    on_x = *on_coords;
    on_y = on_x + 1;
    on_z = on_x + 2;
    norm_x = *normals;
    norm_y = norm_x+1;
    norm_z = norm_x+2;
    curv1_x = *curvatures_1;
    curv1_y = curv1_x + 1;
    curv1_z = curv1_x + 2;
    curv2_x = *curvatures_2;
    curv2_y = curv2_x + 1;
    curv2_z = curv2_x + 2;
    coord_step *= 3;
    on_step = 3;
  }
  
  RefEntity** ent_iter = (RefEntity**)face_handles;
  iBase_ErrorType result = iBase_SUCCESS;
  for (int i = 0; i < count; ++i) {
    RefFace* face = dynamic_cast<RefFace*>(*ent_iter);
    if (!face) {
      ERROR(iBase_INVALID_ENTITY_TYPE, "Non-face input handle.");
      result = iBase_INVALID_ENTITY_TYPE;
    }
    else {
      const CubitVector coords( *coord_x, *coord_y, *coord_z );
      CubitVector on, norm, curv1, curv2;
      Surface* surf = face->get_surface_ptr();
      CubitStatus s = surf->closest_point( coords, &on, &norm, &curv1, &curv2 );
      if (s == CUBIT_FAILURE)
        result = iBase_FAILURE;
      on.get_xyz( *on_x, *on_y, *on_z );
      norm.get_xyz( *norm_x, *norm_y, *norm_z );
      curv1.get_xyz( *curv1_x, *curv1_y, *curv1_z );
      curv2.get_xyz( *curv2_x, *curv2_y, *curv2_z );
    }
    
    ent_iter += ent_step;
    coord_x += coord_step;
    coord_y += coord_step;
    coord_z += coord_step;
    on_x += on_step;
    on_y += on_step;
    on_z += on_step;
    norm_x += on_step;
    norm_y += on_step;
    norm_z += on_step;
    curv1_x += on_step;
    curv1_y += on_step;
    curv1_z += on_step;
    curv2_x += on_step;
    curv2_y += on_step;
    curv2_z += on_step;
  }
  
  KEEP_ARRAY( on_coords );
  KEEP_ARRAY( normals );
  KEEP_ARRAY( curvatures_1 );
  KEEP_ARRAY( curvatures_2 );
  RETURN(result);
}

/**
 * Return the arc length / area / volume of the entities
 * @param gentity_handles Entities for which measure is requested
 * @param gentity_handles_size Number of gentities
 * @param measures Arc length / area / volume of the entities
 * @param measures_length Number of entries in measures
 */
void
iGeom_measure( iGeom_Instance instance,
               /*in*/ const iBase_EntityHandle *gentity_handles,
               int gentity_handles_size,
               /*out*/ double **measures,
               int *measures_allocated,
               int *measures_size,
               int* err) 
{
  RefEntity **handle_array = (RefEntity**)(gentity_handles);

    // check or pre-allocate the measure arrays
  ALLOC_CHECK_ARRAY_NOFAIL(measures, gentity_handles_size);
  for (int i = 0; i < gentity_handles_size; i++)
    (*measures)[i] = handle_array[i]->measure();
  
  RETURN(iBase_SUCCESS);
}

void
iGeom_getEntBoundBox(iGeom_Instance instance,
                     /*in*/ iBase_EntityHandle entity_handle,
                     /*out*/ double* min_x,
                     /*out*/ double* min_y,
                     /*out*/ double* min_z,
                     /*out*/ double* max_x,
                     /*out*/ double* max_y,
                     /*out*/ double* max_z,
                     int* err)
{
  RefEntity* entity = (RefEntity*)entity_handle;
  CubitVector minc, maxc;
  CubitStatus status = iGeom_bounding_box( entity, minc, maxc );
  minc.get_xyz( *min_x, *min_y, *min_z );
  maxc.get_xyz( *max_x, *max_y, *max_z );
  RETURN( (status == CUBIT_SUCCESS ? iBase_SUCCESS : iBase_FAILURE) );
}

/**
 * Return the type of surface as a string; if not a surface, an error is returned
 * @param face_handle Face for which the type is requested
 * @param face_type Type of face, returned as a string
 */
void
iGeom_getFaceType( iGeom_Instance instance,
                   /*in*/ iBase_EntityHandle gentity_handle,
                   /*out*/ char *face_type,
                   int* err,
                   int *face_type_length) 
{
  static const char *surf_types[] = {"cone", "plane", "sphere", "spline", 
                                     "torus", "best_fit", "facet", "undefined"};
  
  RefEntity *this_ent = ENTITY_HANDLE(gentity_handle);
  RefFace *this_face = dynamic_cast<RefFace*>(this_ent);
  if (!this_face)
    RETURN(iBase_INVALID_ENTITY_TYPE);
  GeometryType this_type = this_face->get_surface_ptr()->geometry_type();
  if (this_type < CONE_SURFACE_TYPE || this_type > UNDEFINED_SURFACE_TYPE) {
    RETURN(iBase_FAILURE);
  }
  
  const char* result = surf_types[this_type - CONE_SURFACE_TYPE];
  const int len = strlen(result);
  if (len < *face_type_length) {
    strcpy(face_type, result);
    *face_type_length = len;
  }
  else {
    strncpy(face_type, result, *face_type_length-1);
    face_type[*face_type_length] = '\0';
  }

  RETURN(iBase_SUCCESS);
}

/**
 * Return the bounding boxex of given entities; coordinates returned
 * interleaved.
 * @param gentity_handles The gentities being queried
 * @param min_corners Minimum corner coordinates of the boxes, interleaved
 * @param max_corners Maximum corner coordinates of the boxes, interleaved
 */
void
iGeom_getArrBoundBox (iGeom_Instance instance,
                      /*in*/ const iBase_EntityHandle *gentity_handles,
                      int gentity_handles_size,
                      /*in*/ int storage_order,
                      /*out*/ double **min_corner,
                      int *min_corner_allocated,
                      int *min_corner_size,
                      /*out*/ double **max_corner,
                      int *max_corner_allocated,
                      int *max_corner_size,
                      int* err)
{
    // check or pre-allocate the coordinate arrays
  ALLOC_CHECK_ARRAY(min_corner, 3*gentity_handles_size);
  ALLOC_CHECK_ARRAY(max_corner, 3*gentity_handles_size);
  
  size_t step, init;
  if (storage_order == iBase_BLOCKED) {
    step = 1;
    init = gentity_handles_size;
  }
  else {
    step = 3;
    init = 1;
  }
  double *min_x, *min_y, *min_z, *max_x, *max_y, *max_z;
  min_x = *min_corner;
  max_x = *max_corner;
  min_y = min_x + init;
  max_y = max_x + init;
  min_z = min_y + init;
  max_z = max_y + init;
  
  iBase_ErrorType result = iBase_SUCCESS;
  RefEntity** entities = (RefEntity**)gentity_handles;
  for (int i = 0; i < gentity_handles_size; ++i)
  {
    CubitVector min_c, max_c;
    CubitStatus s = iGeom_bounding_box( entities[i], min_c, max_c );
    if (s != CUBIT_SUCCESS)
      result = iBase_FAILURE;
    min_c.get_xyz( *min_x, *min_y, *min_z );
    max_c.get_xyz( *max_x, *max_y, *max_z );
    
    min_x += step;
    max_x += step;
    min_y += step;
    max_y += step;
    min_z += step;
    max_z += step;
  }

  KEEP_ARRAY(min_corner);
  KEEP_ARRAY(max_corner);
  RETURN(result);
}

void
iGeom_getVtxCoord (iGeom_Instance instance,
                   /*in*/ iBase_EntityHandle vertex_handle,
                   /*out*/ double* x,
                   /*out*/ double* y,
                   /*out*/ double* z,
                   int* err)
{
  RefEntity* entity = (RefEntity*)vertex_handle;
  RefVertex* vertex = dynamic_cast<RefVertex*>(entity);
  if (!vertex) {
    ERROR(iBase_INVALID_ENTITY_TYPE, "Entity not a vertex");
  }
  vertex->coordinates().get_xyz(*x, *y, *z);
  RETURN(iBase_SUCCESS);
}

/**
 * Return the coordinates of the specified vertices; returns error if any
 * of the entities are not gvertices.  Coordinates returned interleaved.
 * @param gentity_handles The gentities being queried
 * @param coordinates The coordinates of the gvertices, interleaved.
 */
void
iGeom_getVtxArrCoords (iGeom_Instance instance,
                       /*in*/ const iBase_EntityHandle *gentity_handles,
                       int gentity_handles_size,
                       /*in*/ int storage_order,
                       /*out*/ double **coordinates,
                       int *coordinates_allocated,
                       int *coordinates_size,
                       int* err)
{
  const RefEntity **handle_array = (const RefEntity**)(gentity_handles);

    // check or pre-allocate the coordinate arrays
  ALLOC_CHECK_ARRAY(coordinates, 3*gentity_handles_size);
  
  CubitVector dumvec;
  const RefVertex *this_vertex;

  double *x, *y, *z;
  size_t step;
  if (storage_order == iBase_BLOCKED) {
    x = *coordinates;
    y = x + gentity_handles_size;
    z = y + gentity_handles_size;
    step = 1;
  }
  else {
    x = *coordinates;
    y = x + 1;
    z = x + 2;
    step = 3;
  }

  for (int i = 0; i < gentity_handles_size; i++) {
    this_vertex = dynamic_cast<const RefVertex*>(handle_array[i]);
    if (NULL == this_vertex) {
      ERROR(iBase_INVALID_ENTITY_TYPE, "Entities passed into getGvertexCoordinates must be vertices.");
    }
    else {
      dumvec = this_vertex->coordinates();
      dumvec.get_xyz(*x,*y,*z);
      x += step;
      y += step;
      z += step;
    }
  }

  KEEP_ARRAY(coordinates);
  RETURN(iBase_SUCCESS);
}

void
iGeom_getPntRayIntsct( iGeom_Instance instance,
                       /*in*/ double x,
                       /*in*/ double y,
                       /*in*/ double z,
                       /*in*/ double dx,
                       /*in*/ double dy,
                       /*in*/ double dz,
                       /*inout*/ iBase_EntityHandle **intersect_entity_handles,
                       int *intersect_entity_handles_allocated,
                       int *intersect_entity_handles_size,
                       /*in*/ int storage_order,
                       /*inout*/ double **intersect_coords,
                       int *intersect_coords_allocated,
                       int *intersect_coords_size,
                       /*inout*/ double **param_coords,
                       int *param_coords_allocated,
                       int *param_coords_size,
                       int* err)
{
  DLIList<double> ray_params;
  DLIList<RefEntity*> entities;
  CubitVector point(x,y,z), dir(dx,dy,dz);
  CubitStatus status = iGeom_fire_ray( point, dir, entities, ray_params );
  if (status != CUBIT_SUCCESS)
    RETURN(iBase_FAILURE);
  
  ALLOC_CHECK_ARRAY_NOFAIL( intersect_entity_handles, entities.size() );
  ALLOC_CHECK_ARRAY_NOFAIL( intersect_coords, 3*ray_params.size() );
  ALLOC_CHECK_ARRAY_NOFAIL( param_coords, ray_params.size() );
  
  size_t init, step;
  if (storage_order == iBase_BLOCKED) {
    init = ray_params.size();
    step = 1;
  }
  else {
    init = 1;
    step = 3;
  }
  
  double *x_iter = *intersect_coords;
  double *y_iter = x_iter + init;
  double *z_iter = y_iter + init;
  for (size_t i = ray_params.size(); i > 0; --i) {
    double t = ray_params.get_and_step();
    CubitVector pos = t * dir + point;
    pos.get_xyz( *x_iter, *y_iter, *z_iter );
    x_iter += step;
    y_iter += step;
    z_iter += step;
  }
  ray_params.copy_to( *param_coords );
  entities.copy_to( (RefEntity**)*intersect_entity_handles );
         
  RETURN(iBase_SUCCESS);  
}

void
iGeom_getPntArrRayIntsct( iGeom_Instance instance,
                          /*in*/ int storage_order,
                          /*in*/ const double *points,
                          int points_size,
                          /*in*/ const double *directions,
                          int directions_size,
                          /*inout*/ iBase_EntityHandle **intersect_entity_handles,
                          int *intersect_entity_handles_allocated,
                          int *intersect_entity_handles_size,
                          /*inout*/ int **offset,
                          int *offset_allocated,
                          int *offset_size,
                          /*inout*/ double **intersect_coords,
                          int *intersect_coords_allocated,
                          int *intersect_coords_size,
                          /*inout*/ double **param_coords,
                          int *param_coords_allocated,
                          int *param_coords_size,
                          int* err)
{
  if (points_size != directions_size || points_size % 3) {
    ERROR(iBase_INVALID_ARGUMENT, "Mismatched or invalid input array size");
  }
  
  const int count = points_size / 3;
  ALLOC_CHECK_ARRAY( offset, count );
  
  const double *px, *py, *pz, *dx, *dy, *dz;
  size_t init, step;
  if (storage_order == iBase_BLOCKED) {
    init = count;
    step = 1;
  }
  else {
    storage_order = iBase_INTERLEAVED;
    init = 1;
    step = 3;
  }
  px = points;
  py = px + init;
  pz = py + init;
  dx = directions;
  dy = dx + init;
  dz = dy + init;
  
  DLIList<RefEntity*> entities, tmp_entities;
  DLIList<double> params, tmp_params;
  std::vector<CubitVector> coords;
  for (int i = 0; i < count; ++i)
  {
    tmp_entities.clean_out();
    tmp_params.clean_out();
    (*offset)[i] = params.size();
    const CubitVector point(*px, *py, *pz), dir(*dx, *dy, *dz);
    CubitStatus s = iGeom_fire_ray( point, dir, tmp_entities, tmp_params );
    if (CUBIT_SUCCESS != s) {
      RETURN(iBase_FAILURE);
    }

    entities += tmp_entities;
    params += tmp_params;
    tmp_params.reset();
    for (int j = tmp_params.size(); j > 0; --j) 
      coords.push_back( tmp_params.get_and_step() * dir + point );
    
    px += step;
    py += step;
    pz += step;
    dx += step;
    dy += step;
    dz += step;
  }
  
  ALLOC_CHECK_ARRAY_NOFAIL( intersect_entity_handles, entities.size() );
  ALLOC_CHECK_ARRAY_NOFAIL( intersect_coords, coords.size() );
  ALLOC_CHECK_ARRAY_NOFAIL( param_coords, params.size() );
  entities.copy_to( (RefEntity**)*intersect_entity_handles );
  params.copy_to( *param_coords );
  
  double *x = *intersect_coords;
  double *y = x + init;
  double *z = y + init;
  for (std::vector<CubitVector>::const_iterator k = coords.begin(); k != coords.end(); ++k)
  {
    k->get_xyz( *x, *y, *z );
    x += step;
    y += step;
    z += step;
  }
  
  KEEP_ARRAY(offset);
  RETURN(iBase_SUCCESS);
}

void
iGeom_getPntClsf (iGeom_Instance instance,
                  /*in*/ double x,
                  /*in*/ double y,
                  /*in*/ double z,
                  /*out*/ iBase_EntityHandle* entity_handle,
                  int* err)
{
  RefEntity** ptr = (RefEntity**)entity_handle;
  const CubitVector pt(x,y,z);
  *ptr = iGeom_get_point_containment( pt );
  RETURN( *ptr ? iBase_SUCCESS : iBase_FAILURE );
}

void
iGeom_getPntArrClsf (iGeom_Instance instance,
                     /*in*/ int storage_order,
                     /*in*/ const double *coords,
                     int coords_size,
                     /*inout*/ iBase_EntityHandle **entity_handles,
                     int *entity_handles_allocated,
                     int *entity_handles_size,
                     int* err)
{
  size_t init, step;
  int count = coords_size / 3;
  if (storage_order == iBase_BLOCKED) {
    init = count;
    step = 1;
  }
  else {
    storage_order = iBase_INTERLEAVED;
    init = 1;
    step = 3;
  }
  const double *x = coords;
  const double *y = x + init;
  const double *z = y + init;
  
  ALLOC_CHECK_ARRAY( entity_handles, count );
  
  RefEntity** array = (RefEntity**)*entity_handles;
  for (int i = 0; i < count; ++i)
  {
    const CubitVector pt( *x, *y, *z );
    array[i] = iGeom_get_point_containment( pt );
    if (!array[i]) {
      RETURN(iBase_FAILURE);
    }
    x += step;
    y += step;
    z += step;
  }
  
  KEEP_ARRAY(entity_handles);
  RETURN(iBase_SUCCESS);
}

/**
 * Return the sense of a face with respect to a region.  Sense is either
 * forward (=1), reverse (=-1), both (=0).  Error is returned
 * if first entity is not a gface or second entity is not a gregion.
 * @param gface face whose sense is being queried.
 * @param gregion region gface is being queried with respect to
 */
void
iGeom_getEntNrmlSense (iGeom_Instance instance,
                       /*in*/ iBase_EntityHandle gface,
                       /*in*/ iBase_EntityHandle gregion,
                       int* rel_sense,
                       int* err)
{
  const RefFace *face_ent = dynamic_cast<const RefFace*>(ENTITY_HANDLE(gface));
  if (NULL == face_ent) {
    ERROR(iBase_INVALID_ENTITY_TYPE, "1st argument to getEntNrmlSense must be a face.");
  }

  // XXX: workaround; remove this when we switch iBase_REGIONs to RefVolume
  const RefVolume *volume_ent;
  Body *body_ent = dynamic_cast<Body*>(ENTITY_HANDLE(gregion));
  if (NULL != body_ent) {
    DLIList<RefEntity*> children;
    body_ent->get_child_ref_entities(children);
    if (children.size() != 1) {
      ERROR(iBase_FAILURE, "Can only support bodies with one volume");
    }
    volume_ent = dynamic_cast<const RefVolume*>(children[0]);
  }
  else {
    volume_ent = dynamic_cast<const RefVolume*>(ENTITY_HANDLE(gregion));
  }

  if (NULL == volume_ent) {
    ERROR(iBase_INVALID_ENTITY_TYPE, "2nd argument to getEntNrmlSense must be a region.");
  }

  *rel_sense = iGeom_get_nonmanifold_sense( face_ent, volume_ent, err );
}

void
iGeom_getArrNrmlSense(iGeom_Instance instance,
                      /*in*/ iBase_EntityHandle const *faces,
                      int faces_size,
                      /*in*/ iBase_EntityHandle const *regions,
                      int regions_size,
                      /*inout*/ int **senses,
                      int *senses_allocated,
                      int *senses_size,
                      int* err)
{
  size_t faces_step, regions_step;
  int count;
  if (faces_size == regions_size) {
    faces_step = regions_step = 1;
    count = faces_size;
  }
  else if (faces_size == 1) {
    faces_step = 0;
    regions_step = 1;
    count = regions_size;
  }
  else if (regions_size == 1) {
    faces_step = 1;
    regions_step = 0;
    count = faces_size;
  }
  else {
    ERROR( iBase_INVALID_ENTITY_COUNT, "Mismatched input array sizes" );
    RETURN (iBase_INVALID_ENTITY_COUNT);
  }
  
  RefEntity** face_iter = (RefEntity**)faces;
  RefEntity** region_iter = (RefEntity**)regions;
  ALLOC_CHECK_ARRAY( senses, count );
  for (int i = 0; i < count; ++i) {
    RefFace *face_ent = dynamic_cast<RefFace*>(*face_iter);
    if (NULL == face_ent) {
      ERROR(iBase_INVALID_ENTITY_TYPE, "1st argument to getArrNrmlSense must be a face.");
    }

    // XXX: workaround; remove this when we switch iBase_REGIONs to RefVolume
    const RefVolume *volume_ent;
    Body *body_ent = dynamic_cast<Body*>(*region_iter);
    if (NULL != body_ent) {
      DLIList<RefEntity*> children;
      body_ent->get_child_ref_entities(children);
      if (children.size() != 1) {
        ERROR(iBase_FAILURE, "Can only support bodies with one volume");
      }
      volume_ent = dynamic_cast<const RefVolume*>(children[0]);
    }
    else {
      volume_ent = dynamic_cast<const RefVolume*>(*region_iter);
    }

    if (NULL == volume_ent) {
      ERROR(iBase_INVALID_ENTITY_TYPE, "2nd argument to getArrNrmlSense must be a region.");
    }

    (*senses)[i] = iGeom_get_nonmanifold_sense( face_ent, volume_ent, err );
    if (iBase_SUCCESS != *err)
      return;

    face_iter += faces_step;
    region_iter += regions_step;
  }
  
  KEEP_ARRAY(senses);
  RETURN(iBase_SUCCESS);
}



/**
 * Return the sense of a gedge with respect to a gface.  Sense is either
 * forward (=1), reverse (=-1), both (=0), or unknown (=2).  Error is returned
 * if first entity is not a gedge or second entity is not a gface.
 * @param gedge Gedge whose sense is being queried.
 * @param gface Gface gedge is being queried with respect to
 */
void
iGeom_getEgFcSense (iGeom_Instance instance,
                    /*in*/ iBase_EntityHandle gedge,
                    /*in*/ iBase_EntityHandle gface,
                    int* rel_sense,
                    int* err)
{
  const RefEdge *edge_ent = dynamic_cast<const RefEdge*>(ENTITY_HANDLE(gedge));
  if (NULL == edge_ent) {
    ERROR(iBase_INVALID_ENTITY_TYPE, "1st argument to getGtangentSense must be an edge.");
  }
  const RefFace *face_ent = dynamic_cast<const RefFace*>(ENTITY_HANDLE(gface));
  if (NULL == face_ent) {
    ERROR(iBase_INVALID_ENTITY_TYPE, "2nd argument to getGtangentSense must be a face.");
  }
  *rel_sense = iGeom_get_nonmanifold_sense( edge_ent, face_ent, err );
}

void
iGeom_getEgFcArrSense(iGeom_Instance instance,
                      /*in*/ iBase_EntityHandle const *edges,
                      int edges_size,
                      /*in*/ iBase_EntityHandle const *faces,
                      int faces_size,
                      /*inout*/ int **senses,
                      int *senses_allocated,
                      int *senses_size,
                      int* err)
{
  size_t faces_step, edges_step;
  int count;
  if (faces_size == edges_size) {
    faces_step = edges_step = 1;
    count = faces_size;
  }
  else if (faces_size == 1) {
    faces_step = 0;
    edges_step = 1;
    count = edges_size;
  }
  else if (edges_size == 1) {
    faces_step = 1;
    edges_step = 0;
    count = faces_size;
  }
  else {
    ERROR( iBase_INVALID_ENTITY_COUNT, "Mismatched input array sizes" );
  }
  
  RefEntity** face_iter = (RefEntity**)faces;
  RefEntity** edge_iter = (RefEntity**)edges;
  ALLOC_CHECK_ARRAY( senses, count );
  for (int i = 0; i < count; ++i) {
    RefEdge *edge_ent = dynamic_cast<RefEdge*>(*edge_iter);
    if (NULL == edge_ent) {
      ERROR(iBase_INVALID_ENTITY_TYPE, "1st argument to getGnormalSense must be a face.");
    }
    RefFace *face_ent = dynamic_cast<RefFace*>(*face_iter);
    if (NULL == face_ent) {
      ERROR(iBase_INVALID_ENTITY_TYPE, "2nd argument to getGnormalSense must be a region.");
    }
    (*senses)[i] = iGeom_get_nonmanifold_sense( edge_ent, face_ent, err );
    if (iBase_SUCCESS != *err)
      return;

    face_iter += faces_step;
    edge_iter += edges_step;
  }
  
  KEEP_ARRAY(senses);
  RETURN(iBase_SUCCESS);
}

/**
 * Return the sense of a gedge with respect to a specified order of
 * vertices bounding the gedge.  Sense is either forward (=1), reverse (=-1), 
 * or closed (=0).  Error is returned if any gentities are not the expected
 * type or if the gedge is bounded by only one gvertex (in this case, use
 * getGtangentSense).
 * @param gedge Gedge whose sense is being queried.
 * @param gvertex1 First gvertex
 * @param gvertex2 Second gvertex
 */
void
iGeom_getEgVtxSense (iGeom_Instance instance,
                     /*in*/ iBase_EntityHandle gedge,
                     /*in*/ iBase_EntityHandle gvertex1,
                     /*in*/ iBase_EntityHandle gvertex2,
                     int* rel_sense,
                     int* err)
{
  const RefEdge *this_edge = dynamic_cast<const RefEdge*>(ENTITY_HANDLE(gedge));
  const RefVertex *vertex1 = dynamic_cast<const RefVertex*>(ENTITY_HANDLE(gvertex1));
  const RefVertex *vertex2 = dynamic_cast<const RefVertex*>(ENTITY_HANDLE(gvertex2));
  if (NULL == this_edge || NULL == vertex1 || NULL == vertex2) {
    ERROR(iBase_INVALID_ENTITY_TYPE, "Bad entity argument to getGvertexTangentSense.");
  }
  *rel_sense = iGeom_edge_vertex_sense( this_edge, vertex1, vertex2, err );
}

void
iGeom_getEgVtxArrSense (iGeom_Instance instance,
                        /*in*/ iBase_EntityHandle const *edges ,
                        int edges_size,
                        /*in*/ iBase_EntityHandle const *start_vertices ,
                        int start_vertices_size,
                        /*in*/ iBase_EntityHandle const *end_vertices ,
                        int end_vertices_size,
                        /*inout*/  int **senses ,
                        int *senses_allocated,
                        int *senses_size,
                        int* err)
{
  int count;
  size_t edge_step, start_step, end_step;
  edge_step = edges_size > 1;
  start_step = start_vertices_size > 1;
  end_step = end_vertices_size > 1;
  count = edges_size;
  if (count != 1) {
    if (start_vertices_size != 1 && start_vertices_size != count) {
      ERROR( iBase_INVALID_ENTITY_COUNT, "Mismatched input array sizes" );
    }
  }
  else
    count = start_vertices_size;
  if (count != 1) {
    if (end_vertices_size != 1 && end_vertices_size != count) {
      ERROR( iBase_INVALID_ENTITY_COUNT, "Mismatched input array sizes" );
    }
  }
  else
    count = end_vertices_size;
    
  ALLOC_CHECK_ARRAY( senses, count );
  
  RefEntity** edge_iter = (RefEntity**)edges;
  RefEntity** start_iter = (RefEntity**)start_vertices;
  RefEntity** end_iter = (RefEntity**)end_vertices;
  for (int i = 0; i < count; ++i) {
    RefEdge* edge = dynamic_cast<RefEdge*>(*edge_iter);
    RefVertex* vtx1 = dynamic_cast<RefVertex*>(*start_iter);
    RefVertex* vtx2 = dynamic_cast<RefVertex*>(*end_iter);
    if (!edge || !vtx1 || !vtx2) {
      ERROR(iBase_INVALID_ENTITY_TYPE, "Bad entity argument to getGvertexTangentSense.");
    }
    
    (*senses)[i] = iGeom_edge_vertex_sense( edge, vtx1, vtx2, err );
    if (iBase_SUCCESS != *err)
      return;
    edge_iter += edge_step;
    start_iter += start_step;
    end_iter += end_step;
  }
  
  KEEP_ARRAY(senses);
  RETURN(iBase_SUCCESS);
}
  
  
    
      

/**
 * Returns the number of gentity_sets contained in a given model or
 * gentity_set one level deep
 * @param gentity_set_handle Entity set being queried
 * @return Number of entity sets in gentity_set_handle
 */
void
iGeom_getNumEntSets (iGeom_Instance instance,
                     /*in*/ iBase_EntitySetHandle entity_set,
                     /*in*/ int num_hops,
                     int* num_ent_sets,
                     int* err)
{
    // HJK: num_hops has to be handled
  if (1 < num_hops) {
    ERROR(iBase_NOT_SUPPORTED, "Num_hops argument not yet supported.");
  }
  
  const RefGroup *this_set = SET_HANDLE(entity_set);
  if (NULL == this_set)
    *num_ent_sets = RefEntityFactory::instance()->num_ref_groups();
  else {
    DLIList<RefEntity*> tmp_ents;
    DLIList<RefGroup*> groups;
    const_cast<RefGroup*>(this_set)->get_child_ref_entities(tmp_ents);
    CAST_LIST(tmp_ents, groups, RefGroup);
    *num_ent_sets = groups.size();
  }

  RETURN(iBase_SUCCESS);
}

/**
 * Returns the gentity_set handles contained in a given model or
 * gentity_set one level deep
 * @param gentity_set_handle Entity set being queried
 * @param contained_gentity_set_handles Number of entity sets in 
 *  gentity_set_handle
 */
void
iGeom_getEntSets (iGeom_Instance instance,
                  /*in*/ iBase_EntitySetHandle entity_set,
                  /*in*/ int num_hops,
                  /*inout*/ iBase_EntitySetHandle **contained_entity_set_handles,
                  int *contained_entity_set_handles_allocated,
                  int *contained_entity_set_handles_size,
                  int *err)
{
  if (1 < num_hops) {
    ERROR(iBase_NOT_SUPPORTED, "Num_hops argument not yet supported.");
  }

  const RefGroup *this_set = SET_HANDLE(entity_set);
  DLIList<RefEntity*> tmp_ents;
  DLIList<RefGroup*> groups;
  if (NULL == this_set)
    RefEntityFactory::instance()->ref_groups(groups);
  else {
    const_cast<RefGroup*>(this_set)->get_child_ref_entities(tmp_ents);
    CAST_LIST(tmp_ents, groups, RefGroup);
  }
  
  ALLOC_CHECK_ARRAY_NOFAIL(contained_entity_set_handles, groups.size());

  groups.copy_to(*SET_HANDLE_ARRAY_PTR(contained_entity_set_handles));
  
  RETURN(iBase_SUCCESS);
}

/**
 * This function allows a new gentity_set to be created.  The user may 
 * set the multiset, ordered, isMesh, flags as needed; otherwise default values 
 * (all false) will be used.  On creation, Entitysets are empty of
 * entities and contained in the parent geometry interface.  They must be
 * explicitly filled with entities using the addGentities call and
 * relationships with other Entitysets must be done through the
 * addEntityset and parent/child relationship calls.
 * @param multiset If true, gentities can appear more than once in this gentity_set
 * @param ordered If true, order of addition and removal is maintained for
 *   this gentity_set
 * @param gentity_set_created Entity_set created by this function
 */
void
iGeom_createEntSet (iGeom_Instance instance,
                    /*in*/ int isList,
                    /*out*/ iBase_EntitySetHandle *entity_set,
                    int* err)
{
  RefGroup* grp = RefEntityFactory::instance()->construct_RefGroup();
  *entity_set = reinterpret_cast<iBase_EntitySetHandle>(grp);
    // need to set a tag denoting multiset or not...
  if (*entity_set == NULL) {
    RETURN(iBase_FAILURE);
  }
  
  else {
    RETURN(iBase_SUCCESS);
  }
}

/**
 *   Destroy the gentity set.  This method only destroys the grouping of
 *   gentities, not the gentities themselves.
 * @param gentity_set Entity_set to be destroyed
 */
void
iGeom_destroyEntSet (iGeom_Instance instance,
                     /*in*/ iBase_EntitySetHandle entity_set,
                     int* err)
{
  if (NULL == entity_set) {
    ERROR(iBase_INVALID_ARGUMENT, "Can't destroy interface set.");
  }
      
  CubitStatus result = RefGroup::delete_group(SET_HANDLE(entity_set));
  if (CUBIT_SUCCESS == result) {
    RETURN(iBase_SUCCESS);
  }
  
  else {
    RETURN(iBase_FAILURE);
  }
  
}

void
iGeom_isList (iGeom_Instance instance,
              /*in*/ iBase_EntitySetHandle entity_set,
              int* result,
              int* err) 
{
  *result = true;
  RETURN (iBase_SUCCESS);
}

/**
 *   Allows the user to explicitly add one or more gentity_sets to
 *   another.  This automatically sets the contained in relationship,
 *   but not the parent/child relationships.  All gentity_set handles
 *   are automatically contained in the parent mesh interface, so
 *   passing in NULL as the first argument results in no action.
 * @param gentity_set Set to which other sets are being added
 * @param gentity_set_handles Sets added to gentity_set
 */
void
iGeom_addEntSet (iGeom_Instance instance,
                 /*in*/ iBase_EntitySetHandle entity_set_to_add,
                 /*inout*/ iBase_EntitySetHandle entity_set_handle,
                 int* err)
{
  iGeom_addEntToSet(instance, reinterpret_cast<iBase_EntityHandle>(entity_set_to_add), entity_set_handle, err);
}

/**
 *   Allows the user to remove one or more gentity_sets from another
 *   gentity_set.  Users cannot delete a contained in relationship of an
 *   gentity set with the parent mesh interface so passing in a NULL
 *   value for the first argument results in no action.
 * @param gentity_set Set from which other sets are being removed
 * @param gentity_set_handles Sets being removed from gentity_set
 */
void
iGeom_rmvEntSet (iGeom_Instance instance,
                 /*in*/ iBase_EntitySetHandle entity_set_to_remove,
                 /*inout*/ iBase_EntitySetHandle entity_set_handle,
                 int* err)
{
  iGeom_rmvEntFromSet(instance, reinterpret_cast<iBase_EntityHandle>(entity_set_to_remove), entity_set_handle, err);
}

/**
 *   Allows the user to explicitly add one or more gentity_sets to
 *   another.  This automatically sets the contained in relationship,
 *   but not the parent/child relationships.  All gentity_set handles
 *   are automatically contained in the parent mesh interface, so
 *   passing in NULL as the first argument results in no action.
 * @param gentity_set Set to which other sets are being added
 * @param gentity_set_handles Sets added to gentity_set
 */
void
iGeom_addEntToSet (iGeom_Instance instance,
                   /*in*/ iBase_EntityHandle entity_to_add,
                   /*inout*/ iBase_EntitySetHandle entity_set_handle,
                   int* err)
{
  if (NULL == entity_to_add) RETURN(iBase_INVALID_ARGUMENT);
  
  CubitStatus status = SET_HANDLE(entity_set_handle)->
    add_ref_entity(const_cast<RefEntity*>(ENTITY_HANDLE(entity_to_add)));
  
  if (CUBIT_SUCCESS != status) {
    ERROR(iBase_FAILURE, "Problem adding entity to another set.");
  }

  RETURN(iBase_SUCCESS);
}

/**
 *   Allows the user to remove one or more gentity_sets from another
 *   gentity_set.  Users cannot delete a contained in relationship of an
 *   gentity set with the parent mesh interface so passing in a NULL
 *   value for the first argument results in no action.
 * @param gentity_set Set from which other sets are being removed
 * @param gentity_set_handles Sets being removed from gentity_set
 */
void
iGeom_rmvEntFromSet (iGeom_Instance instance,
                     /*in*/ iBase_EntityHandle entity_to_remove,
                     /*inout*/ iBase_EntitySetHandle entity_set_handle,
                     int* err)
{
  if (NULL == entity_set_handle) RETURN(iBase_INVALID_ARGUMENT);
  
  CubitStatus status = SET_HANDLE(entity_set_handle)->
    remove_ref_entity(const_cast<RefEntity*>(ENTITY_HANDLE(entity_to_remove)));
  if (CUBIT_SUCCESS != status) {
    ERROR(iBase_FAILURE, "Problem removing entity from a set.");
  }
  
  RETURN(iBase_SUCCESS);
}

void
iGeom_isEntContained (iGeom_Instance instance,
                      /*in*/ iBase_EntitySetHandle containing_entity_set,
                      /*in*/ iBase_EntityHandle contained_entity,
                      int* is_contained,
                      int* err) 
{
    // everything is contained in root set
  if (NULL == containing_entity_set && NULL != contained_entity) {
    *is_contained = true;
    RETURN( iBase_SUCCESS );
  }

  DLIList<RefGroup*> containing_groups;
  RefEntity *contained_entity_temp = const_cast<RefEntity*>(ENTITY_HANDLE(contained_entity));
  RefGroup::get_groups_within(contained_entity_temp, containing_groups, CUBIT_FALSE);
  *is_contained = containing_groups.size() > 0 &&
      containing_groups.move_to(const_cast<RefGroup*>(SET_HANDLE(containing_entity_set))) ;

  RETURN(iBase_SUCCESS);
}

void iGeom_isEntArrContained( iGeom_Instance instance,
                       /*in*/ iBase_EntitySetHandle containing_set,
                       /*in*/ const iBase_EntityHandle* entity_handles,
                       /*in*/ int num_entity_handles,
                    /*inout*/ int** is_contained,
                    /*inout*/ int* is_contained_allocated,
                      /*out*/ int* is_contained_size,
                      /*out*/ int* err )
{
    // go through each entity and look up its dimension
  ALLOC_CHECK_ARRAY(is_contained, num_entity_handles);

  for (int i = 0; i < num_entity_handles; ++i) {
    iGeom_isEntContained( instance, containing_set, entity_handles[i], (*is_contained)+i, err );
    if (iBase_SUCCESS != *err)
      return;
  }

  KEEP_ARRAY(is_contained);
  RETURN(iBase_SUCCESS);
}

void
iGeom_isEntSetContained (iGeom_Instance instance,
                         /*in*/ iBase_EntitySetHandle containing_entity_set,
                         /*in*/ iBase_EntitySetHandle contained_entity_set,
                         int* is_contained,
                         int* err) 
{
  iGeom_isEntContained( instance, containing_entity_set, 
                        reinterpret_cast<iBase_EntityHandle>(contained_entity_set),
                        is_contained, err);
}

/**
 *   Add existing gentities to the gentity_set (do not create them).
 *   Note that if a gentity of dimension d>0 is added to the gentityset,
 *   the lower-dimensional gentities that bound it are not
 *   automatically associated with the gentityset.
 * @param gentity_set Set being added to
 * @param gentity_handles Gentities being added to gentity_set
 */
void
iGeom_addEntArrToSet (iGeom_Instance instance,
                      /*in*/ iBase_EntityHandle const* entity_handles,
                      int entity_handles_size,
                      /*inout*/ iBase_EntitySetHandle entity_set,
                      int* err)
{
  if (NULL == entity_set) RETURN(iBase_INVALID_ARGUMENT);

  RefGroup *this_set = SET_HANDLE(entity_set);
  RefEntity **ent_array = (RefEntity**)(entity_handles);
  CubitStatus status = CUBIT_SUCCESS, tmp_status;
  for (int i = 0; i < entity_handles_size; i++) {
    tmp_status = this_set->add_ref_entity(ent_array[i]);
    if (CUBIT_SUCCESS != tmp_status) status = tmp_status;
  }

  if (CUBIT_SUCCESS != status) {
    ERROR(iBase_FAILURE, "Problem adding entities to a set.");
  }
  
  RETURN(iBase_SUCCESS);
}

/**
 *   Remove existing gentities from the gentity_set (do not delete them)
 * @param gentity_set Set being removed from
 * @param gentity_handles Gentities being removed from gentity_set
 */
void
iGeom_rmvEntArrFromSet (iGeom_Instance instance,
                        /*in*/ iBase_EntityHandle const* entity_handles,
                        int entity_handles_size,
                        /*inout*/ iBase_EntitySetHandle entity_set,
                        int* err)
{
  if (NULL == entity_set) RETURN(iBase_INVALID_ARGUMENT);

  RefGroup *this_set = SET_HANDLE(entity_set);
  RefEntity **ent_array = (RefEntity**)(entity_handles);
  CubitStatus status = CUBIT_SUCCESS, tmp_status;
  for (int i = 0; i < entity_handles_size; i++) {
    tmp_status = this_set->remove_ref_entity(ent_array[i]);
    if (CUBIT_SUCCESS != tmp_status) status = tmp_status;
  }

  if (CUBIT_SUCCESS != status) {
    ERROR(iBase_FAILURE, "Problem removing entities from a set.");
  }
  
  RETURN(iBase_SUCCESS);
}

/**
 * Return the relative and absolute tolerances at the modeler level.  If
 * model does not have a modeler-wide tolerance, zero is returned for both
 * values.
 * @param relative_tolerance Relative tolerance for model as a whole
 * @param absolute_tolerance Absolute tolerance for model as a whole
 */
void
iGeom_getTolerance (iGeom_Instance instance,
                    /*out*/ int* type,
                    /*out*/ double* tolerance,
                    int* err)
{
  *type = 1;
  *tolerance = gqt->get_sme_resabs_tolerance();
  RETURN(iBase_SUCCESS);
}

void
iGeom_getEntTolerance (iGeom_Instance instance,
                       /*in*/ iBase_EntityHandle entity_handle,
                       double* tolerance,
                       int* err)
{
  *tolerance = gqt->get_sme_resabs_tolerance();
  RETURN(iBase_SUCCESS);
}

/**
 * Return the relative and absolute tolerances for specified gentities.  If
 * a gentity does not have a specific tolerance, zero is returned for both
 * values.
 * @param gentity_handles Gentities being queried
 * @param relative_tolerances Relative tolerances
 * @param absolute_tolerances Absolute tolerances
 */
void
iGeom_getArrTolerance (iGeom_Instance instance,
                       /*in*/ iBase_EntityHandle const *gentity_handles,
                       const int gentity_handles_size,
                       /*out*/ double **tolerances,
                       int *tolerances_allocated,
                       int *tolerances_size,
                       int* err)
{
  ALLOC_CHECK_ARRAY_NOFAIL(tolerances, gentity_handles_size);
  double dum_abs_tol = gqt->get_sme_resabs_tolerance();
  for (int i = 0; i < gentity_handles_size; i++) {
    (*tolerances)[i] = dum_abs_tol;
  }
  RETURN(iBase_SUCCESS);
}

void
iGeom_getParametric(iGeom_Instance instance, int* is_parametric, int* err)
{
  *is_parametric = true;
  RETURN(iBase_SUCCESS);
}

/**
 * Return whether a given gentity is parametric or not.  If a gentity
 * is not parametric, all of the following functions will return an error
 * when called on that entity.
 * @param gentity_handle Entity being queried.
 */
void
iGeom_isEntParametric (iGeom_Instance instance,
                       /*in*/ iBase_EntityHandle gentity_handle,
                       int* is_parametric,
                       int* err)
{
  *is_parametric = iGeom_is_parametric( (RefEntity*)gentity_handle );
  RETURN (iBase_SUCCESS);
}

void
iGeom_isArrParametric (iGeom_Instance instance,
                       iBase_EntityHandle const *entity_handles,
                       int entity_handles_size,
                       int **is_parametric,
                       int *is_parametric_allocated,
                       int *is_parametric_size,
                       int* err)
{
  ALLOC_CHECK_ARRAY_NOFAIL( is_parametric, entity_handles_size );
  RefEntity** const ent_array = (RefEntity**)entity_handles;
  for (int i = 0; i < entity_handles_size; ++i)
    (*is_parametric)[i] = iGeom_is_parametric( ent_array[i] );
  RETURN(iBase_SUCCESS);
}

void
iGeom_getEntUVtoXYZ (iGeom_Instance instance,
                     /*in*/ iBase_EntityHandle entity_handle,
                     /*in*/ double u,
                     /*in*/ double v,
                     /*out*/ double* x,
                     /*out*/ double* y,
                     /*out*/ double* z,
                     int* err)
{
  RefEntity* entity = (RefEntity*)entity_handle;
  RefFace* face = dynamic_cast<RefFace*>(entity);
  if (!face) {
    ERROR(iBase_INVALID_ENTITY_TYPE, "Expected face for UV method.");
  }
  
  CubitVector xyz = face->position_from_u_v( u, v );
  xyz.get_xyz( *x, *y, *z );
  RETURN(iBase_SUCCESS);
}

/**
 * Given sets of parametric coordinates, return the corresponding real
 * space coordinates on the gentities.  Input and output coordinates are
 * interleaved.
 * @param gentity_handles Gentities being queried.
 * @param uv Input parametric coordinates
 * @param xyz Output real space coordinates
 */
void
iGeom_getArrUVtoXYZ (iGeom_Instance instance,
                     /*in*/ iBase_EntityHandle const *gentity_handles,
                     int gentity_handles_size,
                     /*in*/ int storage_order,
                     /*in*/ const double *uv,
                     int uv_size,
                     /*out*/ double **coordinates,
                     int *coordinates_allocated,
                     int *coordinates_size,
                     int* err)
{
  size_t ent_step, coord_step, uv_step;
  int count;
  if (2*gentity_handles_size == uv_size) {
    ent_step = uv_step = 1;
    count = gentity_handles_size;
  }
  else if (uv_size == 2) {
    ent_step = 1;
    uv_step = 0;
    count = gentity_handles_size;
  }
  else if (gentity_handles_size == 1) {
    ent_step = 0;
    uv_step = 1;
    count = uv_size/2;
  }
  else {
    ERROR(iBase_INVALID_ENTITY_COUNT, "Mismatched input array sizes.");
  }

  ALLOC_CHECK_ARRAY( coordinates, 3*count );
    
  const double *u, *v;
  double *x, *y, *z;
  u = uv;
  x = *coordinates;
  if (storage_order == iBase_BLOCKED) {
    v = u + (uv_step ? count : 1);
    y = x + count;
    z = y + count;
    coord_step = 1;
  } 
  else {
    storage_order = iBase_INTERLEAVED;
    v = u + 1;
    y = x + 1;
    z = x + 2;
    coord_step = 3;
    uv_step *= 2;
  }
  
  RefEntity** ent = (RefEntity**)gentity_handles;
  for (int i = 0; i < count; ++i) {
    RefFace* face = dynamic_cast<RefFace*>(*ent);
    if (!face) {
      ERROR(iBase_INVALID_ENTITY_TYPE, "Expected face for UV method.");
    }
    
    CubitVector xyz = face->position_from_u_v( *u, *v );
    xyz.get_xyz( *x, *y, *z );
    
    ent += ent_step;
    u += uv_step;
    v += uv_step;
    x += coord_step;
    y += coord_step;
    z += coord_step;
  }
  
  KEEP_ARRAY(coordinates);
  RETURN(iBase_SUCCESS);
}

void
iGeom_getEntUtoXYZ (iGeom_Instance instance,
                    /*in*/ iBase_EntityHandle entity_handle,
                    /*in*/ double u,
                    /*out*/ double* x,
                    /*out*/ double* y,
                    /*out*/ double* z,
                    int* err)
{
  RefEntity* entity = (RefEntity*)entity_handle;
  RefEdge* edge = dynamic_cast<RefEdge*>(entity);
  if (!edge) {
    ERROR(iBase_INVALID_ENTITY_TYPE, "Expected edge for 1-param method.");
  }
  
  CubitVector xyz;
  CubitStatus s = edge->position_from_u( u, xyz );
  xyz.get_xyz( *x, *y, *z );
  RETURN( (s == CUBIT_SUCCESS ? iBase_SUCCESS : iBase_FAILURE) );
}

void
iGeom_getArrUtoXYZ (iGeom_Instance instance,
                    /*in*/ iBase_EntityHandle const *gentity_handles,
                    int gentity_handles_size,
                    /*in*/ const double *u,
                    int u_size,
                    /*in*/ int storage_order,
                    /*out*/ double **coordinates,
                    int *coordinates_allocated,
                    int *coordinates_size,
                    int* err)
{
  int count;
  size_t ent_step, coord_step, u_step;
  if (gentity_handles_size == u_size) {
    ent_step = u_step = 1;
    count = gentity_handles_size;
  }
  else if (u_size == 1) {
    ent_step = 1;
    u_step = 0;
    count = gentity_handles_size;
  }
  else if (gentity_handles_size == 1) {
    ent_step = 0;
    u_step = 1;
    count = u_size;
  }
  else {
    ERROR(iBase_INVALID_ENTITY_COUNT, "Mismatched input array sizes.");
  }
    
  ALLOC_CHECK_ARRAY( coordinates, 3*count );
  
  const double *u_iter;
  double *x, *y, *z;
  u_iter = u;
  x = *coordinates;
  if (storage_order == iBase_BLOCKED) {
    y = x + count;
    z = y + count;
    coord_step = 1;
  } 
  else {
    y = x + 1;
    z = x + 2;
    coord_step = 3;
  }
  
  iBase_ErrorType result = iBase_SUCCESS;
  RefEntity** ent = (RefEntity**)gentity_handles;
  for (int i = 0; i < count; ++i) {
    RefEdge* edge = dynamic_cast<RefEdge*>(*ent);
    if (!edge) {
      ERROR(iBase_INVALID_ENTITY_TYPE, "Expected edge for 1-param method.");
    }
    
    CubitVector xyz;
    CubitStatus s = edge->position_from_u( *u_iter, xyz );
    if (CUBIT_SUCCESS != s)
      result = iBase_FAILURE;
    xyz.get_xyz( *x, *y, *z );
    
    ent += ent_step;
    u_iter += u_step;
    x += coord_step;
    y += coord_step;
    z += coord_step;
  }
  
  KEEP_ARRAY(coordinates);
  RETURN(result);
}

void
iGeom_getEntXYZtoUV (iGeom_Instance instance,
                     /*in*/ iBase_EntityHandle entity_handle,
                     /*in*/ double x,
                     /*in*/ double y,
                     /*in*/ double z,
                     /*out*/ double* u,
                     /*out*/ double* v,
                     int* err)
{
  RefEntity* entity = (RefEntity*)entity_handle;
  RefFace* face = dynamic_cast<RefFace*>(entity);
  if (!face) {
    ERROR(iBase_INVALID_ENTITY_TYPE, "Expected face for UV method.");
  }
  
  const CubitVector xyz( x, y, z );
  CubitStatus s = face->u_v_from_position( xyz, *u, *v );
  RETURN( (s == CUBIT_SUCCESS ? iBase_SUCCESS : iBase_FAILURE) );
}


/**
 * Given sets of real space coordinates, return the corresponding 
 * parametric coordinates on the gentities.  Input and output coordinates 
 * are interleaved.
 * @param gentity_handles Gentities being queried.
 * @param xyz Input real space coordinates
 * @param uv Output parametric coordinates
 */
void
iGeom_getArrXYZtoUV (iGeom_Instance instance,
                     /*in*/ iBase_EntityHandle const *gentity_handles,
                     int gentity_handles_size,
                     /*in*/ int storage_order,
                     /*in*/ const double *coordinates,
                     int coordinates_size,
                     /*out*/ double **uv,
                     int *uv_allocated,
                     int *uv_size,
                     int* err)
{
  int count;
  size_t ent_step, coord_step, uv_step;
  if (3*gentity_handles_size == coordinates_size) {
    ent_step = coord_step = 1;
    count = gentity_handles_size;
  }
  else if (coordinates_size == 3) {
    ent_step = 1;
    coord_step = 0;
    count = gentity_handles_size;
  }
  else if (gentity_handles_size == 1) {
    ent_step = 0;
    coord_step = 1;
    count = coordinates_size/3;
  }
  else {
    ERROR(iBase_INVALID_ENTITY_COUNT, "Mismatched input array sizes.");
  }

  ALLOC_CHECK_ARRAY( uv, 2*count );

  double *u, *v;
  const double *x, *y, *z;
  u = *uv;
  x = coordinates;
  if (storage_order == iBase_BLOCKED) {
    v = u + count;
    y = x + (coord_step ? count : 1);
    z = y + (coord_step ? count : 1);
    uv_step = 1;
  } 
  else {
    storage_order = iBase_INTERLEAVED;
    v = u + 1;
    y = x + 1;
    z = x + 2;
    coord_step *= 3;
    uv_step = 2;
  }
  
  RefEntity** ent = (RefEntity**)gentity_handles;
  for (int i = 0; i < count; ++i) {
    RefFace* face = dynamic_cast<RefFace*>(*ent);
    if (!face) {
      ERROR(iBase_INVALID_ENTITY_TYPE, "Expected face for UV method.");
    }
    
    CubitVector xyz( *x, *y, *z );
    CubitStatus s = face->u_v_from_position( xyz, *u, *v );
    if (CUBIT_SUCCESS != s)
      RETURN(iBase_FAILURE);
    
    ent += ent_step;
    u += uv_step;
    v += uv_step;
    x += coord_step;
    y += coord_step;
    z += coord_step;
  }
  
  KEEP_ARRAY(uv);
  RETURN(iBase_SUCCESS);
}


void
iGeom_getEntXYZtoU (iGeom_Instance instance,
                    /*in*/ iBase_EntityHandle entity_handle,
                    /*in*/ double x,
                    /*in*/ double y,
                    /*in*/ double z,
                    /*out*/ double* u,
                    int* err)
{
  RefEntity* entity = (RefEntity*)entity_handle;
  RefEdge* edge = dynamic_cast<RefEdge*>(entity);
  if (!edge) {
    ERROR(iBase_INVALID_ENTITY_TYPE, "Expected edge for 1-param method.");
  }
  
  const CubitVector xyz( x, y, z );
  *u = edge->u_from_position( xyz );
  RETURN(iBase_SUCCESS);
}


void
iGeom_getArrXYZtoU (iGeom_Instance instance,
                    /*in*/ iBase_EntityHandle const *gentity_handles,
                    int gentity_handles_size,
                    /*in*/ int storage_order,
                    /*in*/ const double *coordinates,
                    int coordinates_size,
                    /*out*/ double **u,
                    int *u_allocated,
                    int *u_size,
                    int* err)
{
  int count;
  size_t ent_step, coord_step;
  if (3*gentity_handles_size == coordinates_size) {
    ent_step = coord_step = 1;
    count = gentity_handles_size;
  }
  else if (coordinates_size == 3) {
    ent_step = 1;
    coord_step = 0;
    count = gentity_handles_size;
  }
  else if (gentity_handles_size == 1) {
    ent_step = 0;
    coord_step = 1;
    count = coordinates_size/3;
  }
  else {
    ERROR(iBase_INVALID_ENTITY_COUNT, "Mismatched input array sizes.");
  }
    
  const double *x, *y, *z;
  x = coordinates;
  if (storage_order == iBase_BLOCKED) {
    y = x + (coord_step ? count : 1);
    z = y + (coord_step ? count : 1);
  } 
  else {
    storage_order = iBase_INTERLEAVED;
    y = x + 1;
    z = x + 2;
    coord_step *= 3;
  }
  
  ALLOC_CHECK_ARRAY( u, count );
  
  RefEntity** ent = (RefEntity**)gentity_handles;
  for (int i = 0; i < count; ++i) {
    RefEdge* edge = dynamic_cast<RefEdge*>(*ent);
    if (!edge) {
      ERROR(iBase_INVALID_ENTITY_TYPE, "Expected edge for 1-param method.");
    }
    
    CubitVector xyz( *x, *y, *z );
    (*u)[i] = edge->u_from_position( xyz );
    
    ent += ent_step;
    x += coord_step;
    y += coord_step;
    z += coord_step;
  }
  
  KEEP_ARRAY(u);
  RETURN(iBase_SUCCESS);
}

void
iGeom_getEntXYZtoUVHint (iGeom_Instance instance,
                         /*in*/ iBase_EntityHandle entity_handle,
                         /*in*/ double x,
                         /*in*/ double y,
                         /*in*/ double z,
                         /*inout*/ double* u,
                         /*inout*/ double* v,
                         int* err)
{
  RefEntity* entity = (RefEntity*)entity_handle;
  RefFace* face = dynamic_cast<RefFace*>(entity);
  if (!face) {
    ERROR(iBase_INVALID_ENTITY_TYPE, "Expected face for UV method.");
  }
  
  CubitVector xyz( x, y, z );
  face->move_to_surface( xyz, u, v );
  RETURN(iBase_SUCCESS);
}


/**
 * Given sets of real space coordinates, return the corresponding 
 * parametric coordinates on the gentities.  Input and output coordinates 
 * are interleaved.
 * @param gentity_handles Gentities being queried.
 * @param xyz Input real space coordinates
 * @param uv Output parametric coordinates
 */
void
iGeom_getArrXYZtoUVHint (iGeom_Instance instance,
                         /*in*/ iBase_EntityHandle const *gentity_handles,
                         int gentity_handles_size,
                         /*in*/ int storage_order,
                         /*in*/ const double *coordinates,
                         int coordinates_size,
                         /*inout*/ double **uv,
                         int *uv_allocated,
                         int *uv_size,
                         int* err)
{
  int count;
  size_t ent_step, coord_step, uv_step;
  if (3*gentity_handles_size == coordinates_size) {
    ent_step = coord_step = 1;
    count = gentity_handles_size;
  }
  else if (coordinates_size == 3) {
    ent_step = 1;
    coord_step = 0;
    count = gentity_handles_size;
  }
  else if (gentity_handles_size == 1) {
    ent_step = 0;
    coord_step = 1;
    count = coordinates_size/3;
  }
  else {
    ERROR(iBase_INVALID_ENTITY_COUNT, "Mismatched input array sizes.");
  }
    
  double *u, *v;
  const double *x, *y, *z;
  u = *uv;
  x = coordinates;
  if (storage_order == iBase_BLOCKED) {
    v = u + count;
    y = x + (coord_step ? count : 1);
    z = y + (coord_step ? count : 1);
    uv_step = 1;
  } 
  else {
    storage_order = iBase_INTERLEAVED;
    v = u + 1;
    y = x + 1;
    z = x + 2;
    coord_step *= 3;
    uv_step = 2;
  }
  
  ALLOC_CHECK_ARRAY( uv, 2*count );
  
  RefEntity** ent = (RefEntity**)gentity_handles;
  for (int i = 0; i < count; ++i) {
    RefFace* face = dynamic_cast<RefFace*>(*ent);
    if (!face) {
      ERROR(iBase_INVALID_ENTITY_TYPE, "Expected face for UV method.");
    }
    
    CubitVector xyz( *x, *y, *z );
    face->move_to_surface( xyz, u, v );
    
    ent += ent_step;
    u += uv_step;
    v += uv_step;
    x += coord_step;
    y += coord_step;
    z += coord_step;
  }
  
  KEEP_ARRAY(uv);
  RETURN(iBase_SUCCESS);
}

void
iGeom_getEntUVRange (iGeom_Instance instance,
                     /*in*/ iBase_EntityHandle entity_handle,
                     /*out*/ double *u_min,
                     /*out*/ double *v_min,
                     /*out*/ double *u_max,
                     /*out*/ double *v_max,
                     int* err)
{
  RefEntity* entity = (RefEntity*)entity_handle;
  RefFace* face = dynamic_cast<RefFace*>(entity);
  if (!face) {
    ERROR(iBase_INVALID_ENTITY_TYPE, "Expected face for UV method.");
  }
  
  CubitBoolean r1 = face->get_param_range_U(*u_min, *u_max);
  CubitBoolean r2 = face->get_param_range_V(*v_min, *v_max);
  RETURN( ((r1 && r2) ? iBase_SUCCESS : iBase_FAILURE) );
}

void
iGeom_getEntURange (iGeom_Instance instance,
                    /*in*/ iBase_EntityHandle entity_handle,
                    /*out*/ double *u_min,
                    /*out*/ double *u_max,
                    int* err)
{
  RefEntity* entity = (RefEntity*)entity_handle;
  RefEdge* edge = dynamic_cast<RefEdge*>(entity);
  if (!edge) {
    ERROR(iBase_INVALID_ENTITY_TYPE, "Expected edge for 1-param method.");
  }
  
  CubitBoolean r1 = edge->get_param_range(*u_min, *u_max);
  RETURN( (r1 ? iBase_SUCCESS : iBase_FAILURE) );
}
    

/**
 * Return the uv range of the specified gentities.  Parameters are interleaved.
 * @param gentity_handles Gentities being queried.
 * @param uv_min Minimum parameters of gentities, interleaved
 * @param uv_max Maximum parameters of gentities, interleaved
 */
void
iGeom_getArrUVRange (iGeom_Instance instance,
                     /*in*/ iBase_EntityHandle const *gentity_handles,
                     int gentity_handles_size,
                     /*in*/ int storage_order,
                     /*out*/ double **uv_min,
                     int *uv_min_allocated,
                     int *uv_min_size,
                     /*out*/ double **uv_max,
                     int *uv_max_allocated,
                     int *uv_max_size,
                     int* err)
{
  ALLOC_CHECK_ARRAY(uv_min, 2*gentity_handles_size);
  ALLOC_CHECK_ARRAY(uv_max, 2*gentity_handles_size);
  
  size_t init, step;
  if (storage_order == iBase_BLOCKED) {
    init = gentity_handles_size;
    step = 1;
  }
  else {
    init = 1;
    step = 2;
  }
  
  double *u_min = *uv_min;
  double *v_min = u_min + init;
  double *u_max = *uv_max;
  double *v_max = u_max + init;
  
  RefEntity** entities = (RefEntity**)gentity_handles;
  for (int i = 0; i < gentity_handles_size; ++i) {
    RefFace* face = dynamic_cast<RefFace*>(entities[i]);
    if (!face) {
      ERROR(iBase_INVALID_ENTITY_TYPE, "Expected face for UV method.");
    }
    
    CubitBoolean r1 = face->get_param_range_U(*u_min, *u_max);
    CubitBoolean r2 = face->get_param_range_V(*v_min, *v_max);
    if (!(r1 && r2)) {
      RETURN(iBase_FAILURE);
    }
    
    u_min += step;
    v_min += step;
    u_max += step;
    v_max += step;
  }
  
  KEEP_ARRAY(uv_min);
  KEEP_ARRAY(uv_max);
  RETURN(iBase_SUCCESS);
}

void
iGeom_getArrURange (iGeom_Instance instance,
                    /*in*/ iBase_EntityHandle const *gentity_handles,
                    const int gentity_handles_size,
                    /*out*/ double **u_min,
                    int *u_min_allocated,
                    int *u_min_size,
                    /*out*/ double **u_max,
                    int *u_max_allocated,
                    int *u_max_size,
                    int* err)
{
  ALLOC_CHECK_ARRAY(u_min, gentity_handles_size);
  ALLOC_CHECK_ARRAY(u_max, gentity_handles_size);
  
  RefEntity** entities = (RefEntity**)gentity_handles;
  for (int i = 0; i < gentity_handles_size; ++i) {
    RefEdge* edge = dynamic_cast<RefEdge*>(entities[i]);
    if (!edge) {
      ERROR(iBase_INVALID_ENTITY_TYPE, "Expected edge for 1-param method.");
    }
    
    CubitBoolean r1 = edge->get_param_range((*u_min)[i], (*u_max)[i]);
    if (!r1) {
      RETURN(iBase_FAILURE);
    }
  }
  
  KEEP_ARRAY(u_min);
  KEEP_ARRAY(u_max);
  RETURN(iBase_SUCCESS);
}

void
iGeom_getEntUtoUV (iGeom_Instance instance,
                   /*in*/ iBase_EntityHandle edge_handle,
                   /*in*/ iBase_EntityHandle face_handle,
                   /*in*/ double in_u,
                   /*out*/ double* u,
                   /*out*/ double* v,
                   int* err)
{
  RefEdge* edge = dynamic_cast<RefEdge*>((RefEntity*)edge_handle);
  RefFace* face = dynamic_cast<RefFace*>((RefEntity*)face_handle);
  if (!edge || !face) {
    RETURN(iBase_INVALID_ENTITY_TYPE);
  }
  
  CubitVector xyz;
  CubitStatus s;
  s = edge->position_from_u( in_u, xyz );
  if (s != CUBIT_SUCCESS)
    RETURN(iBase_FAILURE);

  s = face->u_v_from_position( xyz, *u, *v );
  RETURN( (s == CUBIT_SUCCESS ? iBase_SUCCESS : iBase_FAILURE) );
}

static bool
iGeom_check_array_size(int size1, int size2)
{
    return size1 == 1 || size2 == 1 || size1 == size2;        
}

void
iGeom_getArrUtoUV (iGeom_Instance instance,
                   /*in*/ iBase_EntityHandle const *edge_handles,
                   int edge_handles_size,
                   /*in*/ iBase_EntityHandle const *face_handles,
                   int face_handles_size,
                   /*in*/ const double *in_u,
                   int in_u_size,
                   /*in*/ int storage_order,
                   /*inout*/ double **uv,
                   int *uv_allocated,
                   int *uv_size,
                   int* err)
{
  int count;
  size_t edge_step, face_step, coord_step, in_u_step;

  if (!(iGeom_check_array_size(edge_handles_size, face_handles_size) &&
        iGeom_check_array_size(edge_handles_size, in_u_size) &&
        iGeom_check_array_size(face_handles_size, in_u_size))) {
    ERROR(iBase_INVALID_ENTITY_COUNT, "Mismatched input array sizes.");
  }

  edge_step = (edge_handles_size == 1) ? 0:1;
  face_step = (face_handles_size == 1) ? 0:1;
  in_u_step = (in_u_size == 1) ? 0:1;

  count = std::max(edge_handles_size, std::max(face_handles_size, in_u_size));

  ALLOC_CHECK_ARRAY( uv, 2*count );
  
  const double *in_u_iter;
  double *u, *v;
  in_u_iter = in_u;
  u = *uv;
  if (storage_order == iBase_BLOCKED) {
    v = u + count;
    coord_step = 1;
  } 
  else {
    v = u + 1;
    coord_step = 2;
  }
  
  RefEntity** edge_ent = (RefEntity**)edge_handles;
  RefEntity** face_ent = (RefEntity**)face_handles;
  for (int i = 0; i < count; ++i) {
    RefEdge* edge = dynamic_cast<RefEdge*>(*edge_ent);
    RefFace* face = dynamic_cast<RefFace*>(*face_ent);
    if (!edge || !face) {
      RETURN(iBase_INVALID_ENTITY_TYPE);
    }
    
    CubitVector xyz;
    CubitStatus s;
    s = edge->position_from_u( *in_u_iter, xyz );
    if (CUBIT_SUCCESS != s)
      RETURN(iBase_FAILURE);
    s = face->u_v_from_position( xyz, *u, *v );
    if (CUBIT_SUCCESS != s)
      RETURN(iBase_FAILURE);

    edge_ent += edge_step;
    face_ent += face_step;
    in_u_iter += in_u_step;
    u += coord_step;
    v += coord_step;
  }
  
  KEEP_ARRAY(uv);
  RETURN(iBase_SUCCESS);
}

void
iGeom_getVtxToUV (iGeom_Instance instance,
                  /*in*/ iBase_EntityHandle vertex_handle,
                  /*in*/ iBase_EntityHandle face_handle,
                  /*out*/ double* u,
                  /*out*/ double* v,
                  int* err)
{
  RefVertex* vtx = dynamic_cast<RefVertex*>((RefEntity*)vertex_handle);
  RefFace* face = dynamic_cast<RefFace*>((RefEntity*)face_handle);
  if (!vtx || !face) {
    RETURN(iBase_INVALID_ENTITY_TYPE);
  }

  iBase_ErrorType result = iGeom_get_vtx_to_uv( vtx, face, *u, *v );
  RETURN(result);
}  

void
iGeom_getVtxArrToUV (iGeom_Instance instance,
                     /*in*/ iBase_EntityHandle const *vertex_handles,
                     int vertex_handles_size,
                     /*in*/ iBase_EntityHandle const *face_handles,
                     int face_handles_size,
                     /*in*/ int storage_order,
                     /*inout*/ double **uv,
                     int *uv_allocated,
                     int *uv_size,
                     int* err)
{
  int count;
  size_t vtx_step, face_step, uv_step;
  if (vertex_handles_size == face_handles_size) {
    count = vertex_handles_size;
    vtx_step = face_step = 1;
  }
  else if (face_handles_size == 1) {
    count = vertex_handles_size;
    vtx_step = 1;
    face_step = 0;
  }
  else if (vertex_handles_size == 1) {
    count = face_handles_size;
    vtx_step = 0;
    face_step = 1;
  }
  else {
    ERROR(iBase_INVALID_ENTITY_COUNT, "Mismatched input array sizes.");
  }
  
  ALLOC_CHECK_ARRAY( uv, 2*count );

  double *u, *v;
  u = *uv;
  if (storage_order == iBase_BLOCKED) {
    v = u + count;
    uv_step = 1;
  } 
  else {
    storage_order = iBase_INTERLEAVED;
    v = u + 1;
    uv_step = 2;
  }
  
  RefEntity** vtx_iter = (RefEntity**)vertex_handles;
  RefEntity** face_iter = (RefEntity**)face_handles;
  for (int i = 0; i < count; ++i) {
    RefVertex* vtx = dynamic_cast<RefVertex*>(*vtx_iter);
    RefFace* face = dynamic_cast<RefFace*>(*face_iter);
    if (!vtx || !face) {
      RETURN(iBase_INVALID_ENTITY_TYPE);
    }
    
    iBase_ErrorType rval = iGeom_get_vtx_to_uv( vtx, face, *u, *v );
    if (iBase_SUCCESS != rval) {
      RETURN(rval);
    }

    vtx_iter += vtx_step;
    face_iter += face_step;
    u += uv_step;
    v += uv_step;
  }
  KEEP_ARRAY(uv);
  RETURN(iBase_SUCCESS);
}

void
iGeom_getVtxToU (iGeom_Instance instance,
                 /*in*/ iBase_EntityHandle vertex_handle,
                 /*in*/ iBase_EntityHandle edge_handle,
                 /*out*/ double* u,
                 int* err)
{
  RefVertex* vtx = dynamic_cast<RefVertex*>((RefEntity*)vertex_handle);
  RefEdge* edge = dynamic_cast<RefEdge*>((RefEntity*)edge_handle);
  if (!vtx || !edge) {
    RETURN(iBase_INVALID_ENTITY_TYPE);
  }
  
  iBase_ErrorType result = iGeom_get_vtx_to_u( vtx, edge, *u );
  RETURN(result);
}

void
iGeom_getVtxArrToU (iGeom_Instance instance,
                    /*in*/ iBase_EntityHandle const *vertex_handles ,
                    int vertex_handles_size,
                    /*in*/ iBase_EntityHandle const *edge_handles ,
                    int edge_handles_size,
                    /*inout*/  double **u ,
                    int *u_allocated,
                    int *u_size,
                    int* err)
{
  int count;
  size_t vtx_step, edge_step;
  if (vertex_handles_size == edge_handles_size) {
    count = vertex_handles_size;
    vtx_step = edge_step = 1;
  }
  else if (edge_handles_size == 1) {
    count = vertex_handles_size;
    vtx_step = 1;
    edge_step = 0;
  }
  else if (vertex_handles_size == 1) {
    count = edge_handles_size;
    vtx_step = 0;
    edge_step = 1;
  }
  else {
    ERROR(iBase_INVALID_ENTITY_COUNT, "Mismatched input array sizes.");
  }
  
  ALLOC_CHECK_ARRAY( u, count );
  
  RefEntity** vtx_iter = (RefEntity**)vertex_handles;
  RefEntity** edge_iter = (RefEntity**)edge_handles;
  for (int i = 0; i < count; ++i) {
    RefVertex* vtx = dynamic_cast<RefVertex*>(*vtx_iter);
    RefEdge* edge = dynamic_cast<RefEdge*>(*edge_iter);
    if (!vtx || !edge) {
      RETURN(iBase_INVALID_ENTITY_TYPE);
    }
    
    iBase_ErrorType rval = iGeom_get_vtx_to_u( vtx, edge, (*u)[i] );
    if (iBase_SUCCESS != rval) {
      RETURN(rval);
    }
    
    vtx_iter += vtx_step;
    edge_iter += edge_step;
  }
  KEEP_ARRAY(u);
  RETURN(iBase_SUCCESS);
}

void
iGeom_getEntNrmlUV (iGeom_Instance instance,
                    /*in*/ iBase_EntityHandle entity_handle,
                    /*in*/ double u,
                    /*in*/ double v,
                    /*out*/ double* nrml_i,
                    /*out*/ double* nrml_j,
                    /*out*/ double* nrml_k,
                    int* err)
{
  RefFace* face = dynamic_cast<RefFace*>((RefEntity*)entity_handle);
  if (!face) {
    RETURN (iBase_INVALID_ENTITY_TYPE);
  }
  
  CubitVector normal;
  CubitStatus s = iGeom_normal_from_uv( face, u, v, normal );
  normal.get_xyz( *nrml_i, *nrml_j, *nrml_k );
  RETURN ((s == CUBIT_SUCCESS ? iBase_SUCCESS : iBase_FAILURE));
}

/**
 * Return the normals at specified uv positions on gfaces.  If any
 * gentity input is not a face, returns error.  Input parameters and 
 * output normals are interleaved.
 * @param gface_handles The entities being queried
 * @param parameters The uv parameters of points being queried, interleaved
 * @param normals Normals at specified points, interleaved
 */
void
iGeom_getArrNrmlUV (iGeom_Instance instance,
                    /*in*/ iBase_EntityHandle const *gface_handles,
                    int gface_handles_size,
                    /*in*/ int storage_order,
                    /*in*/ const double *parameters,
                    int parameters_size,
                    /*out*/ double **normals,
                    int *normals_allocated,
                    int *normals_size,
                    int* err)
{
  int count;
  size_t face_step, param_step, norm_step;
  if (2*gface_handles_size == parameters_size) {
    count = gface_handles_size;
    face_step = param_step = 1;
  }
  else if (parameters_size == 2) {
    count = gface_handles_size;
    face_step = 1;
    param_step = 0;
  }
  else if (gface_handles_size == 1) {
    count = parameters_size/2;
    face_step = 0;
    param_step = 1;
  }
  else {
    ERROR(iBase_INVALID_ENTITY_COUNT, "Mismatched input array sizes.");
  }
  
    // check or pre-allocate the coordinate arrays
  ALLOC_CHECK_ARRAY( normals, 3*count );

  const double *u, *v;
  double *x, *y, *z;
  u = parameters;
  x = *normals;
  if (storage_order == iBase_BLOCKED) {
    v = u + (param_step ? count : 1);
    y = x + count;
    z = y + count;
    norm_step = 1;
  }
  else {
    storage_order = iBase_INTERLEAVED;
    v = u + 1;
    y = x + 1;
    z = x + 2;
    norm_step = 3;
    param_step *= 2;
  }
  
  RefEntity** face_iter = (RefEntity**)gface_handles;
  for (int i = 0; i < count; ++i) {
    RefFace* face = dynamic_cast<RefFace*>(*face_iter);
    if (!face) { RETURN (iBase_INVALID_ENTITY_TYPE); }
    
    CubitVector normal;
    CubitStatus s = iGeom_normal_from_uv( face, *u, *v, normal );
    normal.get_xyz( *x, *y, *z );
    if (CUBIT_SUCCESS != s) {
      RETURN(iBase_FAILURE);
    }
    
    face_iter += face_step;
    x += norm_step;
    y += norm_step;
    z += norm_step;
    u += param_step;
    v += param_step;
  }
  
  KEEP_ARRAY(normals);
  RETURN(iBase_SUCCESS);
}


void
iGeom_getEntTgntU (iGeom_Instance instance,
                   /*in*/ iBase_EntityHandle entity_handle,
                   /*in*/ double param_coord,
                   /*out*/ double* tngt_i,
                   /*out*/ double* tngt_j,
                   /*out*/ double* tngt_k,
                   int* err)
{
  double x, y, z;
  iGeom_getEntUtoXYZ( instance, entity_handle, param_coord, &x, &y, &z, err );
  if (iBase_SUCCESS == *err)
    iGeom_getEntTgntXYZ( instance, entity_handle, x, y, z, tngt_i, tngt_j, tngt_k, err );
}

/**
 * Return the tangents at specified u positions on gedges.  If any
 * gentity input is not a face, returns error.  Output normals are 
 * interleaved.
 * @param gentity_handles The gedges being queried
 * @param parameters The u parameters of points being queried
 * @param tangents Tangents at specified points, interleaved
 */
void
iGeom_getArrTgntU (iGeom_Instance instance,
                   /*in*/ iBase_EntityHandle const *gedge_handles,
                   int gedge_handles_size,
                   /*in*/ int storage_order,
                   /*in*/ const double *parameters,
                   int parameters_size,
                   /*out*/ double **tangents,
                   int *tangents_allocated,
                   int *tangents_size,
                   int* err)
{
  iGeom_getArrUtoXYZ(instance, ARRAY_IN(gedge_handles),
                     ARRAY_IN(parameters), storage_order, 
                     ARRAY_INOUT(tangents), err); 
  if (iBase_SUCCESS == *err)
    iGeom_getArrTgntXYZ(instance, ARRAY_IN(gedge_handles), storage_order,
                        *tangents, *tangents_size, ARRAY_INOUT(tangents), err);
}

void
iGeom_getEnt1stDrvt (iGeom_Instance instance,
                     /*in*/ iBase_EntityHandle entity_handle,
                     /*in*/ double u,
                     /*in*/ double v,
                     /*inout*/  double **dvrt_u ,
                     int *dvrt_u_allocated,
                     int *dvrt_u_size,
                     /*inout*/  double **dvrt_v ,
                     int *dvrt_v_allocated,
                     int *dvrt_v_size,
                     int* err)
{
  RefFace* face = dynamic_cast<RefFace*>((RefEntity*)entity_handle);
  if (!face) { 
    ERROR( iBase_INVALID_ENTITY_TYPE, "Derivatives only for faces." );
  }
  
  CubitVector du, dv;
  CubitStatus s = face->get_surface_ptr()->uv_derivitives( u, v, du, dv );
  if (CUBIT_SUCCESS != s) {
    RETURN(iBase_FAILURE);
  }
  
  ALLOC_CHECK_ARRAY_NOFAIL( dvrt_u, 3 );
  ALLOC_CHECK_ARRAY_NOFAIL( dvrt_v, 3 );
  du.get_xyz( *dvrt_u );
  dv.get_xyz( *dvrt_v );
  RETURN (iBase_SUCCESS);
}

void
iGeom_getArr1stDrvt (iGeom_Instance instance,
                     /*in*/ iBase_EntityHandle const *entity_handles ,
                     int entity_handles_size,
                     /*in*/ int storage_order,
                     /*in*/ const  double *uv ,
                     int uv_size,
                     /*inout*/  double **drvt_u ,
                     int *drvt_u_allocated,
                     int *drvt_u_size,
                     /*inout*/  int **u_offset ,
                     int *u_offset_allocated,
                     int *u_offset_size,  
                     /*inout*/  double **drvt_v ,
                     int *drvt_v_allocated,
                     int *drvt_v_size,
                     /*inout*/  int **v_offset ,
                     int *v_offset_allocated,
                     int *v_offset_size,
                     int* err)
{
  if (2*entity_handles_size != uv_size) {
    ERROR(iBase_INVALID_ENTITY_COUNT, "Mismatched input array sizes.");
    RETURN(iBase_INVALID_ENTITY_COUNT);
  }
  
  ALLOC_CHECK_ARRAY( drvt_u, 3*entity_handles_size );
  ALLOC_CHECK_ARRAY( drvt_v, 3*entity_handles_size );
  ALLOC_CHECK_ARRAY( u_offset, entity_handles_size+1 );
  ALLOC_CHECK_ARRAY( v_offset, entity_handles_size+1 );
  
  size_t u_step, du_step, init;
  if (storage_order == iBase_BLOCKED) {
    u_step = du_step = 1;
    init = entity_handles_size;
  }
  else {
    storage_order = iBase_INTERLEAVED;
    u_step = 2;
    du_step = 3;
    init = 1;
  }
  
  const double *u = uv;
  const double *v = u + init;
  double *du_x = *drvt_u;
  double *du_y = du_x + init;
  double *du_z = du_y + init;
  double *dv_x = *drvt_v;
  double *dv_y = dv_x + init;
  double *dv_z = dv_y + init;
  
  int off = 0;
  RefEntity** entities = (RefEntity**)entity_handles;
  for (int i = 0; i < entity_handles_size; ++i) {
    RefFace* face = dynamic_cast<RefFace*>(entities[i]);
    if (!face) { 
      ERROR( iBase_INVALID_ENTITY_TYPE, "Derivatives only for faces." );
    }
    
    CubitVector du, dv;
    CubitStatus s = face->get_surface_ptr()->uv_derivitives( *u, *v, du, dv );
    if (s != CUBIT_SUCCESS) { RETURN(iBase_FAILURE); }
    du.get_xyz( *du_x, *du_y, *du_z );
    dv.get_xyz( *dv_x, *dv_y, *dv_z );
  
    u += u_step;
    v += u_step;
    du_x += du_step;
    du_y += du_step;
    du_z += du_step;
    dv_x += du_step;
    dv_y += du_step;
    dv_z += du_step;

    (*u_offset)[i] = off;
    (*v_offset)[i] = off;
    off += du_step;
  }
  (*u_offset)[entity_handles_size] = off;
  (*v_offset)[entity_handles_size] = off; 

  KEEP_ARRAY( drvt_u );
  KEEP_ARRAY( drvt_v );
  KEEP_ARRAY( u_offset );
  KEEP_ARRAY( v_offset );
  RETURN(iBase_SUCCESS);
}

void
iGeom_getEnt2ndDrvt (iGeom_Instance instance,
                     /*in*/ iBase_EntityHandle entity_handle,
                     /*in*/ double u,
                     /*in*/ double v,
                     /*inout*/  double **dvrt_uu ,
                     int *dvrt_uu_allocated,
                     int *dvrt_uu_size,
                     /*inout*/  double **dvrt_uv ,
                     int *dvrt_uv_allocated,
                     int *dvrt_uv_size,
                     /*inout*/  double **dvrt_vv ,
                     int *dvrt_vv_allocated,
                     int *dvrt_vv_size,
                     int* err)
{
  RETURN(iBase_NOT_SUPPORTED);
}

void
iGeom_getArr2ndDrvt (iGeom_Instance instance,
                     /*in*/ iBase_EntityHandle const *entity_handles ,
                     int entity_handles_size,
                     /*in*/ int storage_order,
                     /*in*/ const  double *uv ,
                     int uv_size,
                     /*inout*/  double **drvt_uu ,
                     int *drvt_uu_allocated,
                     int *drvt_uu_size,
                     /*inout*/  int **uu_offset ,
                     int *uu_offset_allocated,
                     int *uu_offset_size,  
                     /*inout*/  double **drvt_uv ,
                     int *drvt_uv_allocated,
                     int *drvt_uv_size,
                     /*inout*/  int **uv_offset ,
                     int *uv_offset_allocated,
                     int *uv_offset_size,  
                     /*inout*/  double **drvt_vv ,
                     int *drvt_vv_allocated,
                     int *drvt_vv_size,
                     /*inout*/  int **vv_offset ,
                     int *vv_offset_allocated,
                     int *vv_offset_size,
                     int* err)
{
  RETURN(iBase_NOT_SUPPORTED);
}

void
iGeom_getFcCvtrUV (iGeom_Instance instance,
                   /*in*/ iBase_EntityHandle entity_handle,
                   /*in*/ double u,
                   /*in*/ double v,
                   /*out*/ double* cvtr1_i,
                   /*out*/ double* cvtr1_j,
                   /*out*/ double* cvtr1_k,
                   /*out*/ double* cvtr2_i,
                   /*out*/ double* cvtr2_j,
                   /*out*/ double* cvtr2_k,
                   int* err)
{
  double x, y, z;
  iGeom_getEntUVtoXYZ(instance, entity_handle, u, v, &x, &y, &z, err );
  if (*err == iBase_SUCCESS)
    iGeom_getFcCvtrXYZ(instance, entity_handle, x, y, z, 
                       cvtr1_i, cvtr1_j, cvtr1_k,
                       cvtr2_i, cvtr2_j, cvtr2_k,
                       err );
}


void
iGeom_getFcArrCvtrUV (iGeom_Instance instance,
                      /*in*/ const  iBase_EntityHandle *face_handles ,
                      int face_handles_size,
                      /*in*/ int storage_order,
                      /*in*/ const  double *uv ,
                      int uv_size,
                      /*inout*/  double **cvtr_1 ,
                      int *cvtr_1_allocated,
                      int *cvtr_1_size,
                      /*inout*/  double **cvtr_2 ,
                      int *cvtr_2_allocated,
                      int *cvtr_2_size,
                      int* err)
{
  iGeom_getArrUVtoXYZ(instance, ARRAY_IN(face_handles), storage_order, 
                      ARRAY_IN(uv), ARRAY_INOUT(cvtr_1), err);
  if (*err == iBase_SUCCESS)
    iGeom_getEntArrCvtrXYZ(instance, ARRAY_IN(face_handles), storage_order,
                           *cvtr_1, *cvtr_1_size, 
                           ARRAY_INOUT(cvtr_1), ARRAY_INOUT(cvtr_2), err);
}

void
iGeom_isEntPeriodic (iGeom_Instance instance, 
                     /*in*/ iBase_EntityHandle entity_handle,
                     /*out*/ int* in_u,
                     /*out*/ int* in_v,
                     int* err)
{
  CubitStatus s = iGeom_is_periodic( (RefEntity*)entity_handle, *in_u, *in_v );
  RETURN( (s == CUBIT_SUCCESS ? iBase_SUCCESS : iBase_FAILURE) );
}

void
iGeom_isArrPeriodic (iGeom_Instance instance,
                     /*in*/ iBase_EntityHandle const *entity_handles ,
                     const int entity_handles_size,
                     /*inout*/  int **in_uv ,
                     int *in_uv_allocated,
                     int *in_uv_size,
                     int* err)
{
  ALLOC_CHECK_ARRAY( in_uv, 2*entity_handles_size );
  RefEntity** ents = (RefEntity**)entity_handles;
  for (int i = 0; i < entity_handles_size; ++i) {
    CubitStatus s = iGeom_is_periodic( ents[i], (*in_uv)[2*i], (*in_uv)[2*i+1] );
    if (s != CUBIT_SUCCESS) {
      RETURN(iBase_FAILURE);
    }
  }
  KEEP_ARRAY(in_uv);
  RETURN(iBase_SUCCESS);
}

void
iGeom_isFcDegenerate (iGeom_Instance instance,
                      /*in*/ iBase_EntityHandle entity_handle,
                      int* is_degenerate,
                      int* err)
{
  RefFace* face = dynamic_cast<RefFace*>((RefEntity*)entity_handle);
  if (!face) { RETURN (iBase_INVALID_ENTITY_TYPE); }
  *is_degenerate = iGeom_is_face_degenerate(face);
  RETURN (iBase_SUCCESS);
}

void
iGeom_isFcArrDegenerate (iGeom_Instance instance,
                         /*in*/ iBase_EntityHandle const *face_handles ,
                         const int face_handles_size,
                         /*inout*/  int **degenerate ,
                         int *degenerate_allocated,
                         int *degenerate_size,
                         int* err)
{
  ALLOC_CHECK_ARRAY( degenerate, face_handles_size );
  RefEntity** faces = (RefEntity**)face_handles;
  for (int i = 0; i < face_handles_size; ++i) {
    RefFace* face = dynamic_cast<RefFace*>(faces[i]);
    if (!face) { RETURN (iBase_INVALID_ENTITY_TYPE); }
    (*degenerate)[i] = iGeom_is_face_degenerate( face );
  }
  KEEP_ARRAY(degenerate);
  RETURN(iBase_SUCCESS);
}
               

/**
 *   Subtract the gentities in gentity_set_2 from the gentities in
 *   gentity_set_1; the result is returned in gentity_set_1.  If
 *   gentity_set 1 is a multiset and contains m entries of a particular
 *   gentity, say gentity k, and gentity_set 2 contains n entry of gentity
 *   k, the result contains m-n entries of gentity k.  (What about
 *   ordered sets? What if gentity_set 2 is not a subset of gentity set
 *   1?  Should be consistent with STL... )
 */
void
iGeom_subtract (iGeom_Instance instance,
                /*in*/ iBase_EntitySetHandle entity_set_1,
                /*in*/ iBase_EntitySetHandle entity_set_2,
                /*out*/ iBase_EntitySetHandle *result_entity_set,
                int* err)
{
  const RefGroup *set1 = SET_HANDLE(entity_set_1);
  const RefGroup *set2 = SET_HANDLE(entity_set_2);
  RefGroup *set3 = SET_HANDLE(*result_entity_set);
  if (NULL == set3) {
    set3 = RefEntityFactory::instance()->construct_RefGroup();
    *result_entity_set = reinterpret_cast<iBase_EntitySetHandle>(set3);
  }
    
  const_cast<RefGroup*>(set1)->subtract(const_cast<RefGroup*>(set2), set3);
  RETURN(iBase_SUCCESS);
}

/**
 *   Intersect gentity_set_1 and gentity_set_2; the result is returned in
 *   gentity_set_1.  If both are multisets and both contain 2 entries of
 *   gentity k, the intersection contains 2 entries as well.
 */
void
iGeom_intersect (iGeom_Instance instance,
                 /*in*/ iBase_EntitySetHandle entity_set_1,
                 /*in*/ iBase_EntitySetHandle entity_set_2,
                 /*out*/ iBase_EntitySetHandle *result_entity_set,
                 int* err)
{
  const RefGroup *set1 = SET_HANDLE(entity_set_1);
  const RefGroup *set2 = SET_HANDLE(entity_set_2);
  RefGroup *set3 = SET_HANDLE(*result_entity_set);
  if (NULL == set3) {
    set3 = RefEntityFactory::instance()->construct_RefGroup();
    *result_entity_set = reinterpret_cast<iBase_EntitySetHandle>(set3);
  }

  const_cast<RefGroup*>(set1)->intersect(const_cast<RefGroup*>(set2), set3);
  RETURN(iBase_SUCCESS);
}

/**
 *   Union the gentities in gentity_set_1 and gentity_set_2; the result is
 *   returned in gentity_set_1.  The multiset flag in gentityset 1
 *   determines if duplicate entries are allowed.
 */
void
iGeom_unite (iGeom_Instance instance,
             /*in*/ iBase_EntitySetHandle entity_set_1,
             /*in*/ iBase_EntitySetHandle entity_set_2,
             /*out*/ iBase_EntitySetHandle *result_entity_set,
             int* err)
{
  const RefGroup *set1 = SET_HANDLE(entity_set_1);
  const RefGroup *set2 = SET_HANDLE(entity_set_2);
  RefGroup *set3 = SET_HANDLE(*result_entity_set);
  if (NULL == set3) {
    set3 = RefEntityFactory::instance()->construct_RefGroup();
    *result_entity_set = reinterpret_cast<iBase_EntitySetHandle>(set3);
  }

  const_cast<RefGroup*>(set1)->unite(const_cast<RefGroup*>(set2), set3);
  RETURN(iBase_SUCCESS);
}

void
iGeom_copyEnt (iGeom_Instance instance,
               /*in*/ iBase_EntityHandle geom_entity,
               /*out*/ iBase_EntityHandle *geom_entity2,
               int* err)
{
  Body *this_body = dynamic_cast<Body*>(ENTITY_HANDLE(geom_entity));
  RefVolume *this_vol = dynamic_cast<RefVolume*>(ENTITY_HANDLE(geom_entity));
  if (NULL != this_vol || NULL != this_body) {
      // need to get the associated body, since cgm only supports copying bodies,
      // not volumes
    if (NULL == this_body) {
      this_body = this_vol->get_body_ptr();
      if (NULL == this_body) {
        ERROR(iBase_FAILURE, "Can't get body from volume.");
      }
    }

    RefEntity *temp_entity = gmt->copy_body(this_body);
    *geom_entity2 = reinterpret_cast<iBase_EntityHandle>(temp_entity);
  }
  else {
    RefEntity *this_ent = ENTITY_HANDLE(geom_entity);
    RefEntity *temp_entity = gmt->copy_refentity(this_ent);
    *geom_entity2 = reinterpret_cast<iBase_EntityHandle>(temp_entity);
  }

  if (NULL == *geom_entity2) {
    ERROR(iBase_FAILURE, "NULL returned from CGM copy.");
  }
  
  RETURN(iBase_SUCCESS);
}
      
void
iGeom_sweepEntAboutAxis (iGeom_Instance instance,
                         /*in*/ iBase_EntityHandle geom_entity,
                         /*in*/ const double angle,
                         /*in*/ const double axis_normal_x,
                         /*in*/ const double axis_normal_y,
                         /*in*/ const double axis_normal_z,
                         /*out*/ iBase_EntityHandle *geom_entity2,
                         int* err)
{
  DLIList<RefEntity*> ents;
  RefEntity *this_ent = ENTITY_HANDLE(geom_entity);
  if (NULL == this_ent) {
    ERROR(iBase_FAILURE, "Can't get ref entity from specified entity.");
  }
  else if (this_ent->dimension() > 2) {
    ERROR(iBase_INVALID_ENTITY_TYPE, "Can't sweep 3-d entities.");
  }

  RefEntity *new_ent = gmt->copy_refentity(this_ent);
  ents.append(new_ent);
  
  CubitVector origin(0.0, 0.0, 0.0);
  CubitVector direction(axis_normal_x, axis_normal_y, axis_normal_z);

  double radian_angle = angle * DEG_TO_RAD;
  
  CubitStatus result = gmt->sweep_rotational(ents, origin, direction,
                                             radian_angle);
  if (CUBIT_FAILURE == result) {
      // destroy entity we copied before
    gqt->delete_RefEntity(new_ent);
    ERROR(iBase_FAILURE, "Sweep returned failure.");
  }
  
  else {
      // HACK to get last entity created, because cgm doesn't return
      // body created by sweep
    RefEntity *new_body = RefEntityFactory::instance()->get_last_body();
    *geom_entity2 = reinterpret_cast<iBase_EntityHandle>(new_body);

      // now we know it's succeeded, delete original body
    Body *this_body = dynamic_cast<TopologyEntity*>(this_ent)->body();
    if (NULL != this_body)
      gqt->delete_RefEntity(this_body);
    else
      gqt->delete_RefEntity(this_ent);
  }


  if (NULL == *geom_entity2) {
    RETURN(iBase_FAILURE);
  }
  
  RETURN(iBase_SUCCESS);
}

void
iGeom_deleteAll( iGeom_Instance , int* err )
{
  GeometryQueryTool::instance()->delete_geometry();
  *err = iBase_SUCCESS;
}


void
iGeom_deleteEnt (iGeom_Instance instance,
                 /*in*/ iBase_EntityHandle geom_entity,
                 int* err)
{
  RefEntity *this_ent = ENTITY_HANDLE(geom_entity);

    // special case: if this is a volume, delete the body instead
  CubitStatus result;
  RefVolume *this_vol = CAST_TO(this_ent, RefVolume);
  if (NULL != this_vol)
    result = gqt->delete_Body(this_vol->body());
  else    
    result = gqt->delete_RefEntity(this_ent);

  if (CUBIT_FAILURE == result) {
    ERROR(iBase_FAILURE, "Problems deleting entity.");
  }

    // check to see if this was last thing deleted; if so, reset ids
  RefEntityFactory *rfi = RefEntityFactory::instance();
  if (rfi->num_bodies() == 0 && 
      rfi->num_ref_volumes() == 0 &&
      rfi->num_ref_faces() == 0 &&
      rfi->num_ref_edges() == 0 &&
      rfi->num_ref_vertices() == 0)
    rfi->reset_ids();
  
  RETURN(iBase_SUCCESS);
}

void
iGeom_createSphere( iGeom_Instance instance,
                    double radius,
                    iBase_EntityHandle *geom_entity,
                    int* err )
{
  if (radius <= 0.0) {
    ERROR(iBase_INVALID_ARGUMENT, "Sphere radius must be must be positive.");
  }
  
  RefEntity* tmp_body = gmt->sphere( radius );
  *geom_entity = reinterpret_cast<iBase_EntityHandle>(tmp_body);
  RETURN ((tmp_body ? iBase_SUCCESS : iBase_FAILURE));
}


void
iGeom_createPrism( iGeom_Instance instance,
                   /*in*/ double height,
		   /*in*/ int n_sides,
		   /*in*/ double major_rad,
		   /*in*/ double minor_rad,
		   /*out*/ iBase_EntityHandle *geom_entity,
		   int* err )
{
  if ( 0.0>=height ) {
    ERROR(iBase_INVALID_ARGUMENT, "Prism height must be positive.");
  } else if ( 3>n_sides ) {
    ERROR(iBase_INVALID_ARGUMENT, "Prism must have at least three sides.");
  }
  
  RefEntity* tmp_body = gmt->prism( height, n_sides, major_rad, minor_rad );
  *geom_entity = reinterpret_cast<iBase_EntityHandle>(tmp_body);
  RETURN ((tmp_body ? iBase_SUCCESS :iBase_FAILURE));
}


void
iGeom_createBrick (iGeom_Instance instance,
                   /*in*/ double x,
                   /*in*/ double y,
                   /*in*/ double z,
                   /*out*/ iBase_EntityHandle *geom_entity,
                   int* err)
{
  double tmp_x = x;
  double tmp_y = y;
  double tmp_z = z;
  
  if (0.0 == y && 0.0 == z) {
    tmp_y = x;
    tmp_z = x;
  }
  
  if (0.0 >= tmp_x || 0.0 >= tmp_y || 0.0 >= tmp_z) {
    ERROR(iBase_INVALID_ARGUMENT, "Dimensions must be >= 0, or y & z must both be zero.");
  }
    
  RefEntity *temp_body = gmt->brick(tmp_x, tmp_y, tmp_z);
  *geom_entity = reinterpret_cast<iBase_EntityHandle>(temp_body);

  if (NULL == *geom_entity) {
    RETURN(iBase_FAILURE);
  }

  RETURN(iBase_SUCCESS);
}
          
void
iGeom_createCylinder (iGeom_Instance instance,
                      /*in*/ double height,
                      /*in*/ double major_rad,
                      /*in*/ double minor_rad,
                      /*out*/ iBase_EntityHandle *geom_entity,
                      int* err)
{
  double tmp_minor = (0.0 == minor_rad ? major_rad : minor_rad);
  RefEntity *temp_body = 
    gmt->cylinder(height, major_rad, tmp_minor, major_rad);
  *geom_entity = reinterpret_cast<iBase_EntityHandle>(temp_body);


  if (NULL == *geom_entity) {
    RETURN(iBase_FAILURE);
  }

  RETURN(iBase_SUCCESS);
}
      
void
iGeom_createCone (iGeom_Instance instance,
		  /*in*/ double height,
		  /*in*/ double major_rad_base,
		  /*in*/ double minor_rad_base,
		  /*in*/ double rad_top,
		  /*out*/ iBase_EntityHandle *geom_entity,
		  int* err)
{
  double tmp_minor = (0.0 == minor_rad_base ? major_rad_base : minor_rad_base);
  RefEntity *temp_body = 
    gmt->cylinder(height, major_rad_base, tmp_minor, rad_top);
  *geom_entity = reinterpret_cast<iBase_EntityHandle>(temp_body);


  if (NULL == *geom_entity) {
    RETURN(iBase_FAILURE);
  }

  RETURN(iBase_SUCCESS);
}
    
void
iGeom_createTorus (iGeom_Instance instance,
                   /*in*/ double major_rad,
                   /*in*/ double minor_rad,
                   /*out*/ iBase_EntityHandle *geom_entity,
                   int* err)
{
  if (minor_rad >= major_rad) {
    ERROR(iBase_INVALID_ARGUMENT, "Major radius must be greater than minor radius for tori.");
  }
  
  RefEntity *temp_body = gmt->torus(major_rad, minor_rad);
  *geom_entity = reinterpret_cast<iBase_EntityHandle>(temp_body);
   
  if (NULL == *geom_entity) {
    RETURN(iBase_FAILURE);
  }

  RETURN(iBase_SUCCESS);
}

void
iGeom_moveEnt (iGeom_Instance instance,
               /*inout*/ iBase_EntityHandle geom_entity,
               /*in*/ double x,
               /*in*/ double y,
               /*in*/ double z,
               int* err)
{
  CubitVector vec(x, y, z);
  Body *this_bod = dynamic_cast<Body*>(ENTITY_HANDLE(geom_entity));
  CubitStatus result;
  if (NULL != this_bod) {
    result = gqt->translate(this_bod, vec);
    if (CUBIT_SUCCESS != result) {
      ERROR(iBase_FAILURE, "Failed to move body.");
      RETURN(iBase_FAILURE);
    }
    
    RETURN(iBase_SUCCESS);
  }
  
  BasicTopologyEntity *this_bte = dynamic_cast<BasicTopologyEntity*>(ENTITY_HANDLE(geom_entity));
  if (NULL != this_bte) {
      // non-body move; check to see if there are any siblings to this entity in the
      // same body; if so, we can't move it; if not, get the body and move that; if
      // there is no body, it's a free entity and we can move it anyway
    Body *this_body = this_bte->body();
    if (NULL == this_body) {
      result = gqt->translate(this_bte, vec);
      if (CUBIT_SUCCESS != result) {
        ERROR(iBase_FAILURE, "Failed to move entity.");
      }
    }
    else {
      int num_sibs = -1;
      switch (this_bte->dimension()) {
        case 0: num_sibs = this_body->num_ref_vertices(); break;
        case 1: num_sibs = this_body->num_ref_edges(); break;
        case 2: num_sibs = this_body->num_ref_faces(); break;
      }
      if (num_sibs == 1) {
          // ok to move the body instead
        result = gqt->translate(this_body, vec);
        if (CUBIT_SUCCESS != result) {
          ERROR(iBase_FAILURE, "Failed to move body even only one entity of that"
                             " dimension in the body.");
        }
      }
      else {
        ERROR(iBase_FAILURE, "Too many siblings for an entity to move it.");
      }
    }

    RETURN(iBase_SUCCESS);
  }
  
  ERROR(iBase_INVALID_ENTITY_TYPE, "Wrong type of entity specified for move.");
}
      
void
iGeom_rotateEnt (iGeom_Instance instance,
                 /*inout*/ iBase_EntityHandle geom_entity,
                 /*in*/ double angle,
                 /*in*/ double axis_normal_x,
                 /*in*/ double axis_normal_y,
                 /*in*/ double axis_normal_z,
                 int* err)
{
  CubitVector this_axis(axis_normal_x, axis_normal_y, axis_normal_z);
  Body *this_bod = dynamic_cast<Body*>(ENTITY_HANDLE(geom_entity));
  CubitStatus result;
  if (NULL != this_bod) {
    result = gqt->rotate(this_bod, this_axis, angle);
    if (CUBIT_SUCCESS != result) {
      ERROR(iBase_FAILURE, "Failed to rotate body.");
    }
    
    RETURN(iBase_SUCCESS);
  }
  
  BasicTopologyEntity *this_bte = dynamic_cast<BasicTopologyEntity*>(ENTITY_HANDLE(geom_entity));
  if (NULL != this_bte) {
    result = gqt->rotate(this_bte, this_axis, angle);
    if (CUBIT_SUCCESS != result) {
      ERROR(iBase_FAILURE, "Failed to rotate entity.");
    }
    
    RETURN(iBase_SUCCESS);
  }
  
  ERROR(iBase_INVALID_ENTITY_TYPE, "Wrong type of entity specified for move.");
}
      
void
iGeom_reflectEnt (iGeom_Instance instance,
                  /*inout*/ iBase_EntityHandle geom_entity,
                  /*in*/ double plane_normal_x,
                  /*in*/ double plane_normal_y,
                  /*in*/ double plane_normal_z,
                  int* err)
{
  CubitVector this_plane(plane_normal_x, plane_normal_y, plane_normal_z);
  Body *this_bod = dynamic_cast<Body*>(ENTITY_HANDLE(geom_entity));
  DLIList<Body*> bods;
  bods.append(this_bod);
  CubitStatus result;
  if (NULL != this_bod) {
    result = gqt->reflect(bods, this_plane);
    if (CUBIT_SUCCESS != result) {
      ERROR(iBase_FAILURE, "Failed to reflect body.");
    }
    
    RETURN(iBase_SUCCESS);
  }
  
  BasicTopologyEntity *this_bte = dynamic_cast<BasicTopologyEntity*>(ENTITY_HANDLE(geom_entity));
  if (NULL != this_bte) {
    result = gqt->reflect(this_bte, this_plane);
    if (CUBIT_SUCCESS != result) {
      ERROR(iBase_FAILURE, "Failed to reflect entity.");
    }
    
    RETURN(iBase_SUCCESS);
  }
  
  ERROR(iBase_INVALID_ENTITY_TYPE, "Wrong type of entity specified for reflect.");
}

void
iGeom_scaleEnt (iGeom_Instance instance,
                /*inout*/ iBase_EntityHandle geom_entity,
                /*in*/ double scale_x,
                /*in*/ double scale_y,
                /*in*/ double scale_z,
                int* err) 
{
  CubitVector factor(scale_x, scale_y, scale_z);
  Body *this_bod = dynamic_cast<Body*>(ENTITY_HANDLE(geom_entity));
  CubitStatus result;
  if (NULL != this_bod) {
    result = gqt->scale(this_bod, factor);
    if (CUBIT_SUCCESS != result) {
      ERROR(iBase_FAILURE, "Failed to scale body.");
    }
    
    RETURN(iBase_SUCCESS);
  }
  
  BasicTopologyEntity *this_bte = dynamic_cast<BasicTopologyEntity*>(ENTITY_HANDLE(geom_entity));
    // non-body move; check to see if there are any siblings to this entity in the
    // same body; if so, we can't move it; if not, get the body and move that; if
    // there is no body, it's a free entity and we can move it anyway
  Body *this_body = this_bte->body();
  if (NULL == this_body && NULL != this_bte) {
    result = gqt->scale(this_bte, factor);
    if (CUBIT_SUCCESS != result) {
      ERROR(iBase_FAILURE, "Failed to scale entity.");
    }
  }
  else if (NULL != this_body && NULL != this_bte) {
    int num_sibs = -1;
    switch (this_bte->dimension()) {
      case 0: num_sibs = this_body->num_ref_vertices(); break;
      case 1: num_sibs = this_body->num_ref_edges(); break;
      case 2: num_sibs = this_body->num_ref_faces(); break;
          // for TSTT, always scale volumes
      case 3: num_sibs = 1; break;
    }
    if (num_sibs == 1) {
        // ok to scale the body instead
      result = gqt->scale(this_body, factor);
      if (CUBIT_SUCCESS != result) {
        ERROR(iBase_FAILURE, "Failed to scale body even only one entity of that"
                           " dimension in the body.");
      }
    }
    else {
      ERROR(iBase_FAILURE, "Too many siblings for an entity to scale it.");
    }
  }

  RETURN(iBase_SUCCESS);
}

void
iGeom_uniteEnts (iGeom_Instance instance,
                 /*in*/ iBase_EntityHandle const* geom_entities,
                 int geom_entities_size,
                 /*out*/ iBase_EntityHandle *geom_entity,
                 int* err)
{
  DLIList<Body*> bods, orig_bods;
  RefEntity* const* handle_array = ENTITY_HANDLE_CONST_ARRAY(geom_entities);
  for (int i = 0; i < geom_entities_size; i++) {
    Body *this_body = dynamic_cast<Body*>(handle_array[i]);
    if (NULL != this_body) {
      Body *new_body = gmt->copy_body(this_body);
      if (NULL != new_body) {
        bods.append(this_body);
        orig_bods.append(new_body);
      }
    }
  }
  if (bods.size() < geom_entities_size) {
    ERROR(iBase_INVALID_ARGUMENT, "Not all entities input were regions.");
    for (int i = bods.size(); i > 0; i--)
      gqt->delete_RefEntity(bods.get_and_step());
    
    RETURN(iBase_SUCCESS);
  }
  
  DLIList<Body*> new_bods;
  CubitStatus result = gmt->unite(bods, new_bods, false);
  if (CUBIT_SUCCESS != result || 1 != new_bods.size()) {
    ERROR(iBase_FAILURE, "Unite failed.");
  }
    
  else {
    *geom_entity = reinterpret_cast<iBase_EntityHandle>(dynamic_cast<RefEntity*>(new_bods.get()));
    for (int i = orig_bods.size(); i > 0; i--)
      gqt->delete_RefEntity(orig_bods.get_and_step());
  }

  RETURN(iBase_SUCCESS);
}
      
void
iGeom_subtractEnts (iGeom_Instance instance,
                    /*in*/ iBase_EntityHandle blank,
                    /*in*/ iBase_EntityHandle tool,
                    /*out*/ iBase_EntityHandle *geom_entity,
                    int* err)
{
  Body *this_blank = dynamic_cast<Body*>(ENTITY_HANDLE(blank));
  Body *blank_copy = gmt->copy_body(this_blank);
  if (NULL == blank_copy) {
    ERROR(iBase_FAILURE, "Trouble copying blank.");
  }
  Body *this_tool = dynamic_cast<Body*>(ENTITY_HANDLE(tool));
  Body *tool_copy = gmt->copy_body(this_tool);
  if (NULL == tool_copy) {
    ERROR(iBase_FAILURE, "Trouble copying tool.");
    gqt->delete_RefEntity(blank_copy);
    RETURN(iBase_FAILURE);
  }

  DLIList<Body*> blank_list, new_body_list;
  blank_list.append(blank_copy);
  
  RefEntity *new_body = NULL;
  CubitStatus result = gmt->subtract(tool_copy, blank_list, new_body_list);
  if (CUBIT_SUCCESS != result || 0 == new_body_list.size()) {
    ERROR(iBase_FAILURE, "Subtract failed.");
  }
  else {
    new_body = new_body_list.get();
    *geom_entity = reinterpret_cast<iBase_EntityHandle>(new_body);
    gqt->delete_RefEntity(this_blank);
    gqt->delete_RefEntity(this_tool);
  }

  RETURN(iBase_SUCCESS);
}

void
iGeom_intersectEnts ( iGeom_Instance instance,
                     /*in*/ iBase_EntityHandle ent1,
		     /*in*/ iBase_EntityHandle ent2,
		     /*out*/ iBase_EntityHandle *geom_entity,
		     int* err )
{
  Body *this_ent1 = dynamic_cast<Body*>(ENTITY_HANDLE(ent1));
  Body *ent1_copy = gmt->copy_body(this_ent1);
  if (NULL == ent1_copy) {
    ERROR(iBase_FAILURE, "Trouble copying blank.");
  }
  Body *this_ent2 = dynamic_cast<Body*>(ENTITY_HANDLE(ent2));
  Body *ent2_copy = gmt->copy_body(this_ent2);
  if (NULL == ent2_copy) {
    ERROR(iBase_FAILURE, "Trouble copying tool.");
    gqt->delete_RefEntity(ent1_copy);
    RETURN(iBase_FAILURE);
  }

  DLIList<Body*> ent1_list, new_body_list;
  ent1_list.append(ent1_copy);
  
  RefEntity *new_body = NULL;
  CubitStatus result = gmt->intersect(ent2_copy, ent1_list, new_body_list);
  if (CUBIT_SUCCESS != result || 0 == new_body_list.size()) {
    ERROR(iBase_FAILURE, "Intersect failed.");
  }
  else {
    new_body = new_body_list.get();
    *geom_entity = reinterpret_cast<iBase_EntityHandle>(new_body);
    gqt->delete_RefEntity(this_ent2);
    gqt->delete_RefEntity(this_ent1);
  }

  RETURN(iBase_SUCCESS);
}

void
iGeom_sectionEnt (iGeom_Instance instance,
                  /*inout*/ iBase_EntityHandle geom_entity,
                  /*in*/ double plane_normal_x,
                  /*in*/ double plane_normal_y,
                  /*in*/ double plane_normal_z,
                  /*in*/ double offset,
                  /*in*/ int reverse,
                  /*out*/ iBase_EntityHandle *geom_entity2,
                  int* err)
{
  Body *this_body = dynamic_cast<Body*>(ENTITY_HANDLE(geom_entity));
  if (NULL == this_body) {
    RefVolume *this_vol = dynamic_cast<RefVolume*>(ENTITY_HANDLE(geom_entity));
    if (NULL != this_vol)
      this_body = this_vol->get_body_ptr();
  }
  if (NULL == this_body) {
    ERROR(iBase_INVALID_ARGUMENT, "Can only section bodies.");
  }

  CubitVector normal(plane_normal_x, plane_normal_y, plane_normal_z);
  if (normal.length_squared() == 0.0) {
    ERROR(iBase_INVALID_ARGUMENT, "Zero-length vector input.");
  }
  
  CubitVector point1 = normal * CubitVector(1.0, 0.0, 0.0);
  if (point1.length_squared() == 0.0)
    point1 = normal * CubitVector(0.0, 1.0, 0.0);
    
  CubitVector point2 = normal * point1;
  CubitVector point3(0.0, 0.0, 0.0);

  if (0.0 != offset) {
    normal.normalize();
    normal *= offset;
    point1 += normal;
    point2 += normal;
    point3 += normal;
  }

  
  DLIList<Body*> blank_list, new_body_list;
  blank_list.append(gmt->copy_body(this_body));
  CubitStatus result = gmt->section(blank_list, point1, point2, point3,
                                    new_body_list, !reverse,
                                    false);
  if (CUBIT_SUCCESS != result || 0 == new_body_list.size()) {
    gqt->delete_RefEntity(blank_list.get());
    ERROR(iBase_FAILURE, "Section failed.");
  }
    
  else {
      // need to assign it to a RE* first so the void cast gets done right
    RefEntity *new_body = new_body_list.get();
    *geom_entity2 = reinterpret_cast<iBase_EntityHandle>(new_body);
      // also, delete the original body, now that the section worked
    gqt->delete_RefEntity(this_body);
  }

  RETURN(iBase_SUCCESS);
}

void
iGeom_imprintEnts (iGeom_Instance instance,
                   /*in*/ iBase_EntityHandle const* gentity_handles,
                   int gentity_handles_size,
                   int* err) 
{
  if (gentity_handles_size < 1) // GMT::imprint segfaults if passed an empty list
    RETURN(iBase_SUCCESS);

  DLIList<Body*> bods;
  DLIList<RefVolume*> vols, temp_vols;
  RefEntity* const* handle_array = ENTITY_HANDLE_CONST_ARRAY(gentity_handles);
  CubitStatus status = CUBIT_SUCCESS;
  for (int i = 0; i < gentity_handles_size; i++) {
    Body *temp_bod = dynamic_cast<Body*>(handle_array[i]);
    if (NULL != temp_bod) {
      bods.append_unique(temp_bod);
      continue;
    }
    
    RefVolume *temp_vol = dynamic_cast<RefVolume*>(handle_array[i]);
    if (NULL != temp_vol) {
      TopologyEntity *topo_ent = dynamic_cast<TopologyEntity*>(handle_array[i]);
      if (NULL == topo_ent) {
        status = CUBIT_FAILURE;
        continue;
      }
      temp_bod = topo_ent->body();
      if (NULL == temp_bod) {
        status = CUBIT_FAILURE;
        continue;
      }
      bods.append_unique(temp_bod);
      continue;
    }

      // if we've gotten here, it's an error
    status = CUBIT_FAILURE;
  }
  
  if (CUBIT_SUCCESS != status) RETURN(iBase_FAILURE);

  DLIList<Body*> temp_bods;
  status = GeometryModifyTool::instance()->imprint(bods, temp_bods, false);
  
  RETURN(iBase_SUCCESS);
}
  
void
iGeom_mergeEnts (iGeom_Instance instance,
                 /*in*/ iBase_EntityHandle const* gentity_handles,
                 int gentity_handles_size,
                 double tolerance,
                 int* err) 
{
  double old_factor = GeometryQueryTool::instance()->get_geometry_factor();
  if (tolerance != old_factor) 
    GeometryQueryTool::instance()->set_geometry_factor(tolerance*1.0e6);
  else old_factor = 0.0;
  
  DLIList<Body*> bods;
  DLIList<RefVolume*> vols, temp_vols;
  DLIList<RefFace*> faces;
  DLIList<RefEdge*> edges;
  DLIList<RefVertex*> verts;
  RefEntity* const* handle_array = ENTITY_HANDLE_CONST_ARRAY(gentity_handles);
  for (int i = 0; i < gentity_handles_size; i++) {
    TopologyEntity *topo_ent = dynamic_cast<TopologyEntity*>(handle_array[i]);
    if (NULL == topo_ent) continue;
    Body *temp_bod;
    RefVolume *temp_vol;
    RefFace *temp_face;
    RefEdge *temp_edge;
    RefVertex *temp_vert;
    switch (handle_array[i]->dimension()) {
      case -1:
          // it should be a body
        temp_bod = dynamic_cast<Body*>(handle_array[i]);
        if (NULL == temp_bod) RETURN(iBase_FAILURE);
        temp_vols.clean_out();
        topo_ent->ref_volumes(temp_vols);
        vols += temp_vols;
        break;
      case 0:
        temp_vert = dynamic_cast<RefVertex*>(handle_array[i]);
        if (NULL == temp_vert) RETURN(iBase_FAILURE);
        verts.append(temp_vert);
        break;
      case 1:
        temp_edge = dynamic_cast<RefEdge*>(handle_array[i]);
        if (NULL == temp_edge) RETURN(iBase_FAILURE);
        edges.append(temp_edge);
        break;
      case 2:
        temp_face = dynamic_cast<RefFace*>(handle_array[i]);
        if (NULL == temp_face) RETURN(iBase_FAILURE);
        faces.append(temp_face);
        break;
      case 3:
        temp_vol = dynamic_cast<RefVolume*>(handle_array[i]);
        if (NULL == temp_vol) RETURN(iBase_FAILURE);
        vols.append(temp_vol);
        break;
    }
  }
  
  CubitStatus status = CUBIT_SUCCESS, temp_status;
    
  if (verts.size() != 0) {
    temp_status = MergeTool::instance()->merge_refvertices(verts, false);
    if (CUBIT_SUCCESS != temp_status) status = temp_status;
  }
    
  if (edges.size() != 0) {
    temp_status = MergeTool::instance()->merge_refedges(edges, true, false);
    if (CUBIT_SUCCESS != temp_status) status = temp_status;
  }
    
  if (faces.size() != 0) {
    temp_status = MergeTool::instance()->merge_reffaces(faces, false);
    if (CUBIT_SUCCESS != temp_status) status = temp_status;
  }
    
  if (vols.size() != 0) {
    temp_status = MergeTool::instance()->merge_volumes(vols, false);
    if (CUBIT_SUCCESS != temp_status) status = temp_status;
  }
    
  if (bods.size() != 0) {
    temp_status = MergeTool::instance()->merge_bodies(bods);
    if (CUBIT_SUCCESS != temp_status) status = temp_status;
  }

  if (0 != old_factor)
    GeometryQueryTool::instance()->set_geometry_factor(old_factor);
    
  if (CUBIT_SUCCESS != status) {
    RETURN(iBase_FAILURE);
  }
  
  else {
    RETURN(iBase_SUCCESS);
  }
}

} // extern "C"

/********************* HELPER FUNCTION IMPLEMENTATIONS ***************************/

static void tokenize( const std::string& str, 
                      std::vector<std::string>& tokens )
{
  char delim =  str[0];
  std::string::size_type last = str.find_first_not_of( delim, 1 );
  std::string::size_type pos  = str.find_first_of( delim, last );
  while (std::string::npos != pos && std::string::npos != last) {
    tokens.push_back( str.substr( last, pos - last ) );
    last = str.find_first_not_of( delim, pos );
    pos  = str.find_first_of( delim, last ); 
  }
}

// Expect option of the form "NAME=VALUE".
// If NAME portion matches, pass back VALUE and return true.
// Otherwise, leave 'value' unchanged and return false.
static bool match_option( const std::string& opt,
                          const char* name,
                          std::string& value )
{
  std::string::size_type len = strlen( name );
  if (opt[len] != '=')
    return false;
  if (opt.compare( 0, len, name, len ))
    return false;
  value = opt.substr( len + 1 );
  return true;
}

void
iGeom_load_cub_geometry(const char *name, int* err) 
{
  FILE *cubfile = fopen(name, "rb");
  if (NULL == cubfile) RETURN(iBase_FILE_NOT_FOUND);

    // first get the length
  //int result = fseek(cubfile, 0, SEEK_END);
  int result = fseek(cubfile, 0, 2);
  if (result) {
    ERROR(iBase_FAILURE, "Couldn't seek to end of .cub file.");
  }
  
  long endpos = ftell(cubfile);
  
  char magic_str[4] = {'\0', '\0', '\0', '\0'};
  fread(magic_str, 1, 4, cubfile);
  if (!strcmp(magic_str, "CUBE")) {
    ERROR(iBase_NOT_SUPPORTED, "Wrong magic string in .cub file.");
  }
  
    // get the model header
  //result = fseek(cubfile, 4, SEEK_SET);
  result = fseek(cubfile, 4, 0);
  if (result) {
    ERROR(iBase_FAILURE, "Seek failed reading cub file header.");
  }
  int header[6];
  fread(header, 4, 6, cubfile);
  int num_models = header[2];
  int model_table_offset = header[3];

    // get the model table
  int model_entries[36], model_offset[6], 
    model_length[6], model_type[6];
  if (num_models > 6) {
    ERROR(iBase_INVALID_ARGUMENT, "Too many models in .cub file.");
  }

  if (model_table_offset+24*num_models-1 > endpos) {
    ERROR(iBase_FAILURE, "Reading model table will go past end of file.");
  }
    
  //result = fseek(cubfile, model_table_offset, SEEK_SET);
  result = fseek(cubfile, model_table_offset, 0);
  if (result) {
    ERROR(iBase_FAILURE, "Seek failed seeking to model table.");
  }
  fread(model_entries, 4, 6*num_models, cubfile);
  for (int i = 0; i < num_models; i++) {
    model_type[i] = model_entries[6*i];
    model_offset[i] = model_entries[6*i+1];
    model_length[i] = model_entries[6*i+2];
  }
  
    // now read the actual models; but, need to fclose first
  fclose(cubfile);
  cubfile = NULL;
  
  const char *model_type_str[6] = {"", "ACIS_SAT", "", "FACET", "", "GRANITE"};
  for (int i = 0; i < num_models; i++) {
    if ((1 != model_type[i] && 3 != model_type[i] && 5 != model_type[i]) ||
        model_length[i] == 0) continue;

    if (model_offset[i] > endpos) {
      ERROR(iBase_FAILURE, "Reading actual model will go past end of file.");
    }

      // extract the geom file from the cub file
    char *tmp_name = tmpnam(NULL);
    FILE *tmp_file = fopen(tmp_name, "w+b");
    if(NULL == tmp_file){
      ERROR(iBase_FAILURE, "Couldn't create temporary file.");
    }
    NCubitFile::CCubitFile ccf;
    unsigned int retval = ccf.Open(name, NULL, NULL);
    if (0 != retval) {
      std::cerr << "Failed to open." << std::endl;
      ERROR(iBase_FAILURE, "Failed to begin read of temporary file.");
    }
    retval = ccf.BeginReadGeomModel(model_type[i], tmp_file);
    if (0 != retval) {
      std::cerr << "Failed to begin read." << std::endl;
      ERROR(iBase_FAILURE, "Failed to begin read of temporary file.");
    }
      
    retval = ccf.EndReadGeomModel();
    if (0 != retval) {
      std::cerr << "Failed to end read." << std::endl;
      ERROR(iBase_FAILURE, "Failed to end read of temporary file.");
    }
    fclose(tmp_file);

    CubitStatus status = gqt->
      import_solid_model(tmp_name, model_type_str[model_type[i]], 
                         NULL, false);
  
    if (CUBIT_SUCCESS != status) {
      std::string msg = std::string("Trouble loading geometry file of type ") + 
        model_type_str[model_type[i]] + std::string(".");
      ERROR(iBase_FAILURE, msg.c_str());
    }
  }

  if (NULL != cubfile) fclose(cubfile);

  RETURN(iBase_SUCCESS);
}

static void
iGeom_get_adjacent_entities( const RefEntity *from, 
                             const int to_dim,
                             DLIList<RefEntity*> &adj_ents,
                             int* err ) 
{
  TopologyEntity *topo_ent = const_cast<TopologyEntity*>(dynamic_cast<const TopologyEntity*>(from));
  if (NULL == topo_ent) {
    ERROR(iBase_INVALID_ARGUMENT, "Bad entity input.");
  }
  
  adj_ents.clean_out();
  static DLIList<RefVertex*> tmp_verts;
  static DLIList<RefEdge*> tmp_edges;
  static DLIList<RefFace*> tmp_faces;
  static DLIList<Body*> tmp_bodies;
  
  switch (to_dim) {
    case 0:
      tmp_verts.clean_out();
      topo_ent->ref_vertices(tmp_verts);
      CAST_LIST_TO_PARENT(tmp_verts, adj_ents);
      break;
    case 1:
      tmp_edges.clean_out();
      topo_ent->ref_edges(tmp_edges);
      CAST_LIST_TO_PARENT(tmp_edges, adj_ents);
      break;
    case 2:
      tmp_faces.clean_out();
      topo_ent->ref_faces(tmp_faces);
      CAST_LIST_TO_PARENT(tmp_faces, adj_ents);
      break;
    case 3:
      tmp_bodies.clean_out();
      topo_ent->bodies(tmp_bodies);
      CAST_LIST_TO_PARENT(tmp_bodies, adj_ents);
      break;
    default:
      ERROR(iBase_INVALID_ARGUMENT, "Bad input dimension getting adjacent entities.");
  }
  
  RETURN(iBase_SUCCESS);
}

/* Common implementation for both single-entity and array functions. */
static CubitStatus iGeom_closest_point( RefEntity* this_entity,
                                        const CubitVector& near,
                                        CubitVector& on )
{
  RefEdge *this_edge;
  RefFace *this_face;
  Surface *this_surf;
  CubitStatus status;
 
  switch (this_entity->dimension()) {
    case 0:
      on = dynamic_cast<RefVertex*>(this_entity)->coordinates();
      status = CUBIT_SUCCESS;
      break;
    case 1:
      this_edge = dynamic_cast<RefEdge*>(this_entity);
      if (NULL == this_edge) return CUBIT_FAILURE;

      status = this_edge->closest_point(near, on);

      if (debug) {
        std::cout << "Edge " << this_edge->id() << " closest point to (" 
                  << near.x() << ", " << near.y() << ", " << near.z() 
                  << ") is "
                  << on.x() << ", " << on.y() << ", " << on.z() 
                  << ")" << std::endl;
      }
      break;
    case 2:
      this_face = dynamic_cast<RefFace*>(this_entity);
      if (NULL == this_face) return CUBIT_FAILURE;
      this_surf = this_face->get_surface_ptr();
      if (NULL == this_surf) return CUBIT_FAILURE;
      status = this_surf->closest_point( near, &on );
      break;
    default:
        // just copy over the coordinates
      on = near;
      status = CUBIT_SUCCESS;
      break;
  }
  
  return status;
}

static CubitStatus
iGeom_closest_point_and_normal( RefEntity* this_entity, 
                                const CubitVector& near,
                                CubitVector& on,
                                CubitVector& normal )
{
  RefEdge* this_edge;
  RefFace* this_face;
  Surface* this_surf;
  CubitStatus status;
  switch (this_entity->dimension()) {
    case 0:
      on = dynamic_cast<RefVertex*>(this_entity)->coordinates();
      normal.set( 0, 0, 0 );
      status = CUBIT_SUCCESS;
      break;
    case 1:
      this_edge = dynamic_cast<RefEdge*>(this_entity);
      if (NULL == this_edge) return CUBIT_FAILURE;
      
      status = this_edge->closest_point( near, on );
      if (debug) {
        std::cout << "Edge " << this_edge->id() << " closest point to (" 
                  << near.x() << ", " << near.y() << ", " << near.z() 
                  << ") is "
                  << on.x() << ", " << on.y() << ", " << on.x()
                  << ")" << std::endl;
      }
      break;
    case 2:
      this_face = dynamic_cast<RefFace*>(this_entity);
      if (NULL == this_face) return CUBIT_FAILURE;
      this_surf = this_face->get_surface_ptr();
      if (NULL == this_surf) return CUBIT_FAILURE;
      status = this_surf->closest_point( near, &on, &normal );
      if (this_surf->bridge_sense() == CUBIT_REVERSED)
        normal = -normal;
      break;
    default:
        // just copy over the coordinates
      on = near;
      normal.set( 0, 0, 0 );
      status = CUBIT_SUCCESS;
      break;
  }
  
  return status;
}

static CubitStatus
iGeom_bounding_box( RefEntity* entity, CubitVector& minc, CubitVector& maxc )
{
  CubitBox box;
  if (BasicTopologyEntity* bte = dynamic_cast<BasicTopologyEntity*>(entity))
    box = bte->bounding_box();
  else if(Body* body = dynamic_cast<Body*>(entity))
    box = body->bounding_box();
  else {
    iGeom_setLastError(iBase_INVALID_ENTITY_HANDLE, "Entities passed into gentityBoundingBox must be vertex, edge, face, or region."); 
    return CUBIT_FAILURE;
  }
  
  minc = box.minimum();
  maxc = box.maximum();
  return CUBIT_SUCCESS;
}

/** We could use Smits' algorithm here, but only if we turned of
    floating-point exceptions */
static inline void box_min_max( double dir,
                                double min,
                                double max,
                                double pt,
                                double& tmin,
                                double& tmax )
{
  if (dir == 0) {
    tmin = -INFINITY;
    tmax = INFINITY;
  }
  else if (dir > 0.0) {
    tmin = (min - pt) / dir;
    tmax = (max - pt) / dir;
  }
  else {
    tmin = (max - pt) / dir;
    tmax = (min - pt) / dir;
  }
}

static bool
iBase_intersect_ray_box( const CubitBox& box,
                         const CubitVector& point,
                         const CubitVector& direction )
{
  double txmin, txmax, tymin, tymax, tzmin, tzmax;
  box_min_max( direction.x(), box.minimum().x(), box.maximum().x(), point.x(), txmin, txmax );
  box_min_max( direction.y(), box.minimum().y(), box.maximum().y(), point.y(), tymin, tymax );
  if (txmin > tymax || tymin > txmax)
    return false;
  
  if (tymin > txmin)
    txmin = tymin;
  if (tymax < txmax)
    txmax = tymax;
  
  box_min_max( direction.z(), box.minimum().z(), box.maximum().z(), point.z(), tzmin, tzmax );
  if (txmin > txmax || tzmin > txmax)
    return false;
  
  return true;
}


static CubitStatus
iGeom_fire_ray( const CubitVector& point,
                const CubitVector& direction,
                DLIList<RefEntity*>& entities,
                DLIList<double>& ray_params )
{
  const double EPSILON = 0.0;
  CubitStatus s;
  
    // get all free entities in model
  DLIList<RefEntity*> target_entities;
  s = GeometryQueryTool::instance()->get_free_ref_entities( target_entities );
  if (CUBIT_SUCCESS != s) return s;
  DLIList<Body*> bodies;
  GeometryQueryTool::instance()->bodies( bodies );
  CAST_LIST_TO_PARENT( bodies, target_entities );
    
    // do ray fire at list of free entities
  return GeometryQueryTool::instance()->
    fire_ray( point, direction, target_entities, ray_params, 0, EPSILON, &entities );
}

static RefEntity* point_classification( const CubitVector& pt, RefVertex* vtx )
{
  return (pt - vtx->coordinates()).length_squared() > GEOMETRY_RESABS*GEOMETRY_RESABS ? 0 : vtx;
}

static RefEntity* point_classification( const CubitVector& pt, RefEdge* edge )
{
  CubitVector closest;
  edge->closest_point_trimmed( pt, closest );
  if ((pt - closest).length_squared() > GEOMETRY_RESABS*GEOMETRY_RESABS)
    return 0;
  
  if (RefEntity* vtx = point_classification( pt, edge->start_vertex() ))
    return vtx;
  else if (RefEntity* vtx = point_classification( pt, edge->end_vertex() ))
    return vtx;
  else
    return edge;
}

static RefEntity* point_classification( const CubitVector& pt, RefFace* face )
{
  CubitBox extents = face->bounding_box();
  if (extents.distance_squared(pt) > GEOMETRY_RESABS*GEOMETRY_RESABS)
    return 0;
  
  CubitVector closest;
  face->find_closest_point_trimmed( pt, closest );
  if ((pt - closest).length_squared() > GEOMETRY_RESABS*GEOMETRY_RESABS)
    return 0;
  
  DLIList<RefEdge*> edges;
  face->ref_edges( edges );
  edges.last();
  for (int i = 0; i < edges.size(); ++i)
    if (RefEntity* ent = point_classification( pt, edges.step_and_get() ))
      return ent;
  
  return face;
}

static RefEntity* point_classification( const CubitVector& pt, Body* body )
{
  CubitVector nonconst_pt(pt);
  CubitPointContainment pc = body->point_containment( nonconst_pt );
  if (CUBIT_PNT_INSIDE == pc) 
    return body;
  else if (CUBIT_PNT_BOUNDARY != pc)
    return 0;
  
    // If we're here, then we're on the boundary.  
    // Find which boundary entity we're on.
  DLIList<RefFace*> faces;
  body->ref_faces( faces );
  faces.last();
  for (int i = 0; i < faces.size(); ++i)
    if (RefEntity* ent = point_classification( pt, faces.step_and_get() ))
      return ent;
  
    // We don't appear to be on any face.  Is the tolerance for
    // Body::point_classification something other than GEOMETRY_RESABS??
  return body;
}



static RefEntity*
iGeom_get_point_containment( const CubitVector& pt )
{
  DLIList<RefEntity*> ents;
  gqt->get_free_ref_entities( ents );
  
  ents.reset();
  for (int i = 0; i < ents.size(); ++i)
    if (RefVertex* vtx = dynamic_cast<RefVertex*>(ents.get_and_step()))
      if (RefEntity* ent = point_classification( pt, vtx ))
        return ent;
  for (int i = 0; i < ents.size(); ++i)
    if (RefEdge* edge = dynamic_cast<RefEdge*>(ents.get_and_step()))
      if (RefEntity* ent = point_classification( pt, edge ))
        return ent;
  for (int i = 0; i < ents.size(); ++i)
    if (RefFace* face = dynamic_cast<RefFace*>(ents.get_and_step()))
      if (RefEntity* ent = point_classification( pt, face ))
        return ent;
        
  DLIList<Body*> bodies;
  gqt->bodies( bodies );
  bodies.reset();
  for (int i = 0; i < bodies.size(); ++i)
   if (RefEntity* ent = point_classification( pt, bodies.get_and_step() ))
     return ent;
  
  return 0;
}

static int iGeom_get_nonmanifold_sense( const BasicTopologyEntity* child,
                                        const BasicTopologyEntity* parent,
                                        int* err )
{
  DLIList<SenseEntity*> se_list(1);
  const_cast<BasicTopologyEntity*>(child)
    ->get_sense_entities( se_list, const_cast<BasicTopologyEntity*>(parent) );
  if (se_list.size() == 0) {
    iGeom_setLastError((*err = iBase_INVALID_ENTITY_HANDLE), "Relative senes of unrelated entities"); 
    return 2;
  }
  se_list.reset();
  CubitSense sense = se_list.get_and_step()->get_sense();
  for (int i = se_list.size() - 1; i > 0; --i)
    if (se_list.get_and_step()->get_sense() != sense)
      return 0;
  return sense == CUBIT_FORWARD ? 1 : -1;
}

static
int iGeom_edge_vertex_sense( const RefEdge* cedge,
                             const RefVertex* vtx1,
                             const RefVertex* vtx2,
                             int* err )
{
  RefEdge* edge = const_cast<RefEdge*>(cedge);
  if (edge->start_vertex() == vtx1 && edge->end_vertex() == vtx2) 
    return vtx1 == vtx2 ? 0 : 1;
  else if (edge->start_vertex() == vtx2 && edge->end_vertex() == vtx1)
    return -1;
  iGeom_setLastError((*err = iBase_INVALID_ENTITY_HANDLE), "Relative senes of unrelated entities"); 
  return 2;
}

static int iGeom_is_parametric( RefEntity* entity )
{
  if (RefFace* face = dynamic_cast<RefFace*>(entity))
    return face->is_parametric();
  else if (dynamic_cast<RefEdge*>(entity))
    return true;
  else
    return false;
}

static iBase_ErrorType
iGeom_get_vtx_to_u(RefVertex* vertex,
                   RefEdge* edge,
                   double& u)
{
  if (edge->start_vertex() == vertex)
    u = edge->start_param();
  else if (edge->end_vertex() == vertex)
    u = edge->end_param();
  else
    return iBase_INVALID_ARGUMENT;
  return iBase_SUCCESS;
}

static iBase_ErrorType
iGeom_get_vtx_to_uv(RefVertex* vertex, RefFace* face, double& u, double& v)
{
  DLIList<RefVertex*> pts;
  TopologyEntity *topo_ent = dynamic_cast<const TopologyEntity*>(face);
  topo_ent->ref_vertices(pts);

  pts.reset();
  for(int i=0; i<pts.size(); i++) {
    if (pts.get_and_step() == vertex) {
      face->u_v_from_position( vertex->coordinates(), u, v );         
      return iBase_SUCCESS;
    }
  }
  return iBase_INVALID_ARGUMENT;
}

static CubitStatus
iGeom_normal_from_uv( RefFace* face, double u, double v, CubitVector& normal )
{
  
  CubitVector coords = face->position_from_u_v( u, v );
  Surface* surf_ptr = face->get_surface_ptr();
  CubitStatus rval = surf_ptr->closest_point_uv_guess( coords, u, v, 0, &normal );
  if (surf_ptr->bridge_sense() == CUBIT_REVERSED)
    normal = -normal;
  return rval;
}

static CubitStatus
iGeom_is_periodic( RefEntity* entity, int& u, int& v )
{
  if (RefEdge* edge = dynamic_cast<RefEdge*>(entity)) {
    u = edge->is_periodic();
    v = CUBIT_FALSE;
  }
  else if (RefFace* face = dynamic_cast<RefFace*>(entity)) {
    double period;
    u = face->is_periodic_in_U(period);
    v = face->is_periodic_in_V(period);
  }
  else {
    u = v = CUBIT_FALSE;
  }
  
  return CUBIT_SUCCESS;
}

static bool
iGeom_is_face_degenerate( RefFace* face )
{
  double param = 0;
  CubitBoolean b1 = face->is_singular_in_U(param);
  CubitBoolean b2 = face->is_singular_in_V(param);
  return b1 || b2;
}

static iBase_ErrorType 
process_attribs(iGeom_Instance instance) 
{
    // go through all entities, checking for remaining simple attribs
  DLIList<RefEntity*> ref_list;
  DLIList<CubitSimpleAttrib*> csa_list;
  const char *type_names[] = {"vertex", "curve", "surface", "volume", "body"};
  iBase_ErrorType result = iBase_SUCCESS;

  for (int dim = 0; dim < 5; dim++) {
    ref_list.clean_out();
    CubitStatus status = gqt->ref_entity_list(type_names[dim], ref_list, false);
    if (CUBIT_SUCCESS != status) return iBase_FAILURE;
    
    for (int i = ref_list.size(); i > 0; i--) {
      csa_list.clean_out();

        // get all the csa's still on this entity
      RefEntity *this_ent = ref_list.get_and_step();
      TopologyEntity *topo_ent = dynamic_cast<TopologyEntity*>(this_ent);
      topo_ent->bridge_manager()->topology_bridge()->get_simple_attribute(csa_list);

      for (int csa = csa_list.size(); csa > 0; csa--) {
        // see if there's a tag for this csa already
        CubitSimpleAttrib *csa_ptr = csa_list.get_and_step();
        const char *tag_name = csa_ptr->character_type().c_str();
        long tag = TM->getTagHandle(tag_name);
        if (tag < 0) continue;
        else if (0 == tag) {
            // don't have a tag handle for this yet - make one
          result = TM->create_csa_tag(tag_name, &tag);
          if (iBase_SUCCESS != result || 0 == tag) return iBase_FAILURE;
        }
        
          // now set the tag
        iBase_ErrorType tmp_result = TM->set_csa_tag(this_ent, tag, csa_ptr);
        if (iBase_SUCCESS != tmp_result) result = tmp_result;
    
          // now destroy this csa
        delete csa_ptr;
      }
    }
  }

  return iBase_SUCCESS;
}


template <typename T> static inline 
int count_type( const DLIList<CubitEntity*>& list )
{
  int count = 0, size = list.size();
  for (int i = 0; i < size; ++i)
    if (dynamic_cast<T*>(list[i]))
      ++count;
  return count;
}

static inline
int count_ibase_type( int ibase_type, const DLIList<CubitEntity*>& list, int* err )
{
  *err = iBase_SUCCESS;
  switch (ibase_type) {
    case iBase_ALL_TYPES: return list.size() - count_type<RefGroup>(list);
    case iBase_REGION:    return count_type<Body>(list);
    case iBase_FACE:      return count_type<RefFace>(list);
    case iBase_EDGE:      return count_type<RefEdge>(list);
    case iBase_VERTEX:    return count_type<RefVertex>(list);
    default:
      *err = iBase_INVALID_ENTITY_TYPE;
      iGeom_setLastError( *err );
      return -1;
  }
}

template <typename TARGET_TYPE, typename LIST_TYPE> static inline 
void append_type( const DLIList<CubitEntity*>& source_list,
                  DLIList<LIST_TYPE*>& target_list )
{
  int size = source_list.size();
  for (int i = 0; i < size; ++i) 
    if (TARGET_TYPE* ent = dynamic_cast<TARGET_TYPE*>(source_list[i]))
      target_list.append(ent);
}

template <typename SKIP_TYPE, typename LIST_TYPE> static inline 
void append_not_type( const DLIList<CubitEntity*>& source_list,
                      DLIList<LIST_TYPE*>& target_list )
{
  int size = source_list.size();
  for (int i = 0; i < size; ++i) 
    if (!dynamic_cast<SKIP_TYPE*>(source_list[i]))
      if (RefEntity* ent = dynamic_cast<LIST_TYPE*>(source_list[i]))
        target_list.append(ent);
}

// Returns count of entities appended.
template <typename TARGET_TYPE> static inline 
int append_type( const DLIList<CubitEntity*>& source_list,
                 iBase_EntityHandle* array, int array_size )
{
  RefEntity* re_ptr;
  int len = source_list.size();
  int count = 0;
  for (int i = 0; i < len; ++i) {
    if (TARGET_TYPE* ent = dynamic_cast<TARGET_TYPE*>(source_list[i])) {
      if (count < array_size)
        array[count] = reinterpret_cast<iBase_EntityHandle>(re_ptr = ent);
      ++count;
    }
  }
  return count;
}

template <typename SKIP_TYPE> static inline 
int append_not_type( const DLIList<CubitEntity*>& source_list,
                     iBase_EntityHandle* array, int array_size )
{
  int len = source_list.size();
  int count = 0;
  for (int i = 0; i < len; ++i) {
    if (!dynamic_cast<SKIP_TYPE*>(source_list[i])) {
      if (count == array_size) 
        return -1;
      else if (RefEntity* ent = dynamic_cast<RefEntity*>(source_list[i]))
        array[count++] = reinterpret_cast<iBase_EntityHandle>(ent);
    }
  }
  return count;
}

static 
void copy_ibase_type( int ibase_type, 
                      const DLIList<CubitEntity*>& list,
                      iBase_EntityHandle** entity_handles,
                      int* entity_handles_alloc,
                      int* entity_handles_size,
                      int* err )
{
  int count;
  if (*entity_handles_alloc == 0) {
    count = count_ibase_type( ibase_type, list, err );
    if (count < 0)
      return;
    *entity_handles = (iBase_EntityHandle*)malloc( count * sizeof(iBase_EntityHandle) );
    if (!*entity_handles) 
      RETURN(iBase_MEMORY_ALLOCATION_FAILED);
    *entity_handles_alloc = count;
  }
  
  switch (ibase_type) {
    case iBase_ALL_TYPES:
      count = append_not_type<RefGroup>(list,*entity_handles, *entity_handles_alloc);
      break;
    case iBase_REGION:
      count = append_type<Body>(list,*entity_handles, *entity_handles_alloc);
      break;
    case iBase_FACE:
      count = append_type<RefFace>(list,*entity_handles, *entity_handles_alloc);
      break;
    case iBase_EDGE:
      count = append_type<RefEdge>(list,*entity_handles, *entity_handles_alloc);
      break;
    case iBase_VERTEX:
      count = append_type<RefVertex>(list,*entity_handles, *entity_handles_alloc);
      break;
    default:
      RETURN(iBase_INVALID_ENTITY_TYPE);
      break;
  }
  
  *entity_handles_size = count;
  if (count > *entity_handles_alloc)
    RETURN(iBase_BAD_ARRAY_DIMENSION);
  
  RETURN(iBase_SUCCESS);
}

static 
void append_ibase_type( int ibase_type, 
                        const DLIList<CubitEntity*>& source_list,
                        DLIList<RefEntity*>& target_list,
                        int* err )
{
  switch (ibase_type) {
    case iBase_ALL_TYPES:
      append_not_type<RefGroup>(source_list, target_list);
      break;
    case iBase_REGION:
      append_type<Body>(source_list, target_list);
      break;
    case iBase_FACE:
      append_type<RefFace>(source_list, target_list);
      break;
    case iBase_EDGE:
      append_type<RefEdge>(source_list, target_list);
      break;
    case iBase_VERTEX:
      append_type<RefVertex>(source_list, target_list);
      break;
    default:
      RETURN(iBase_INVALID_ENTITY_TYPE);
      break;
  }
  
  RETURN(iBase_SUCCESS);
}

static 
void append_all_ibase_type( int ibase_type, 
                            DLIList<RefEntity*>& target_list,
                            int* err )
{
  RefEntityFactory *const ref = RefEntityFactory::instance();
  if (ibase_type == iBase_ALL_TYPES) {
    for (int i = 0; i < 4; ++i) {
      DLIList<RefEntity*> tmp;
      ref->ref_entity_list( iGeom_entity_type_names[i], tmp );
      target_list += tmp;
    }
  }
  else if (abs(ibase_type) < iBase_ALL_TYPES) {
    ref->ref_entity_list( iGeom_entity_type_names[ibase_type], target_list );
  }
  else {
    RETURN(iBase_INVALID_ENTITY_TYPE);
  }
  
  RETURN(iBase_SUCCESS);
}
