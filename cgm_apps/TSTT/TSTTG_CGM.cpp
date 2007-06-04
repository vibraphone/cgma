#include "TSTTG_CGM.h"

#include <iostream>
#include <math.h>

const bool debug = false;

#define CHECK_SIZE(array, type, this_size)  \
  if (0 == array ## _allocated || array ## _allocated < this_size) {\
    if (NULL != array) free(array); \
    array = (type*)malloc(this_size*sizeof(type));\
    array ## _allocated=this_size;\
    if (NULL == array) {TSTTG_processError(TSTTB_MEMORY_ALLOCATION_FAILED, \
          "Couldn't allocate array.");RETURN(TSTTB_MEMORY_ALLOCATION_FAILED); }\
  }; \
  array ## _size = this_size;

// TAG_CHECK_SIZE is like CHECK_SIZE except it checks for and makes the allocated memory
// size a multiple of sizeof(void*), and the pointer is assumed to be type char*
#define TAG_CHECK_SIZE(array, allocated, size)  \
  if (allocated < size) {\
    if (NULL != array) free(array); \
    allocated=size; \
    if (allocated%sizeof(void*) != 0) allocated=(size/sizeof(void*)+1)*sizeof(void*);\
    array = (char*)malloc(allocated); \
    if (NULL == array) {TSTTG_processError(TSTTB_MEMORY_ALLOCATION_FAILED, \
          "Couldn't allocate array.");RETURN(TSTTB_MEMORY_ALLOCATION_FAILED); }\
  }
#define TAG_HANDLE(handle) reinterpret_cast<int>(handle)
#define CONST_TAG_HANDLE(handle) reinterpret_cast<const int>(const_cast<TSTTG_TagHandle>(handle))

#define ENTITY_HANDLE(handle) reinterpret_cast<RefEntity*>(handle)
#define SET_HANDLE(handle) reinterpret_cast<RefGroup*>(handle)

#define CONST_ENTITY_HANDLE(handle) reinterpret_cast<const RefEntity*>(handle)
#define CONST_SET_HANDLE(handle) reinterpret_cast<const RefGroup*>(handle)

#define ENTITY_HANDLE_ARRAY(array) reinterpret_cast<RefEntity**>(array)
#define SET_HANDLE_ARRAY(array) reinterpret_cast<RefGroup**>(array)

#define CONST_ENTITY_HANDLE_ARRAY(array) reinterpret_cast<const RefEntity**>(array)
#define CONST_SET_HANDLE_ARRAY(array) reinterpret_cast<const RefGroup**>(array)

#define ENTITY_HANDLE_ARRAY_PTR(array) reinterpret_cast<RefEntity***>(array)
#define SET_HANDLE_ARRAY_PTR(array) reinterpret_cast<RefGroup***>(array)

#define CAST_TO_VOID(ptr) reinterpret_cast<void*>(ptr)
#define TAG_HANDLE_ARRAY_INOUT(a) reinterpret_cast<int**>(a), a ## _allocated, a ## _size

#define TM reinterpret_cast<CGMTagManager*>(instance)

// map from CGM's dimension to TSTT's entity type
const TSTTG_GentityType tsttg_type_table[] =
{
  TSTTG_VERTEX,         // RefVertex
  TSTTG_EDGE,           // RefEdge
  TSTTG_FACE,           // RefFace
  TSTTG_REGION          // RefVolume
};

const TSTTG_TagValueType tstt_tag_value_table[] = 
{
  TSTTG_BYTES,
  TSTTG_INTEGER,
  TSTTG_DOUBLE,
  TSTTG_ENTITY_HANDLE
};

const double RAD_TO_DEG = 180.0 / acos(-1.0);
const double DEG_TO_RAD = 1.0 / RAD_TO_DEG;

#include "CATSTT.hpp"
#include "RefEntityFactory.hpp"
#include "BasicTopologyEntity.hpp"
#include "RefGroup.hpp"
#include "Body.hpp"
#include "RefVertex.hpp"
#include "RefEdge.hpp"
#include "RefFace.hpp"
#include "RefVolume.hpp"
#include "CubitVector.hpp"

#include "FacetQueryEngine.hpp"
#include "AcisQueryEngine.hpp"
#include "AcisModifyEngine.hpp"
#include "CGMApp.hpp"
#include "GeometryQueryTool.hpp"
#include "GeometryModifyTool.hpp"
#include "Surface.hpp"
#include "BasicTopologyEntity.hpp"
#include "CubitFile.hpp"
#include "MergeTool.hpp"

// declare private-type functions here, so they aren't visible outside
// this implementation file
void TSTTG_get_adjacent_entities(const RefEntity *from, const int to_dim,
                                 DLIList<RefEntity*> &adj_ents);

#define gqt GeometryQueryTool::instance()
#define gmt GeometryModifyTool::instance()

const char *TSTTG_entity_type_names[] = {"vertex", "curve", "surface", "volume", "body"};

TSTTB_Error TSTTG_LAST_ERROR;

// user defined constructor
enum TSTTB_ErrorType
TSTTG_ctor(TSTTG_Instance *instance) 
{
    // initialize CGM
  CGMApp::instance()->startup(0, NULL);

    // assume for now we want acis
  AcisQueryEngine::instance();
  AcisModifyEngine::instance();

  CGMTagManager *tm = new CGMTagManager();

  CGMApp::instance()->attrib_manager()->silent_flag(true);
  CGMApp::instance()->attrib_manager()->auto_flag(true);

    // return the tagmanager as the instance
  *instance = tm;

  RETURN(TSTTB_SUCCESS);
}

enum TSTTB_ErrorType
TSTTG_dtor(TSTTG_Instance instance) 
{
  if (NULL != TM) delete TM;
  
    // shut down CGM
  CGMApp::instance()->shutdown();

  RETURN(TSTTB_SUCCESS);
}

// user defined static methods: (none)

// user defined non-static methods:
/**
 * Load a model specified by name. Which formats are supported and the
 * specific meaning of this name string (e.g. file name, model name,
 * etc.) are implementation-dependent.  Options are also implementation-
 * dependent.
 * @param name Name of the model
 * @param options String options 
 */
enum TSTTB_ErrorType
TSTTG_gLoad (TSTTG_Instance instance,
             /*in*/ const char *name,
             /*in*/ CONST_ARRAY_IN_DECL(char*, options))
{
    // process options
  bool parallel = false;
  int par_load_option = -1;
  for (int i = 0; i < options_size; i++) {
    if (!strcmp(options[i], "PARALLEL") || !strcmp(options[i], "parallel")) {
      parallel = true;
    }
    
      // look for parallel load option
    if (!strcmp(options[i], "GEOM_BCAST") || !strcmp(options[i], "geom_bcast"))
      par_load_option = 0;
    
    if (!strcmp(options[i], "GEOM_BCAST_AND_DELETE") || !strcmp(options[i], "geom_bcast_and_delete"))
      par_load_option = 1;
    
    if (!strcmp(options[i], "GEOM_SCATTER") || !strcmp(options[i], "geom_scatter"))
      par_load_option = 2;
  }
    
  if (strstr(name, ".cub") != NULL) {
    return TSTTG_load_cub_geometry(name);
  }
  
  CubitStatus status = gqt->read_geometry_file(name);
  if (CUBIT_SUCCESS != status) 
    TSTTG_processError(TSTTB_FILE_NOT_FOUND, "Trouble loading geometry file.");

  RETURN(TSTTB_SUCCESS);
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
enum TSTTB_ErrorType
TSTTG_gSave (TSTTG_Instance instance,
             /*in*/ const char *name,
             /*in*/ CONST_ARRAY_IN_DECL(char*, options))
{
    // process options (none right now...)
  DLIList<RefEntity*> bodies;
  int num_ents_exported;
  CubitString cubit_version(" (TSTTG)");
  CubitStatus status = gqt->export_solid_model(bodies, name, "ACIS_SAT",
                                               num_ents_exported, cubit_version);
  if (CUBIT_SUCCESS != status) 
    TSTTG_processError(TSTTB_FAILURE, "Trouble saving geometry file.");

  RETURN(TSTTB_SUCCESS);
}

enum TSTTB_ErrorType
TSTTG_load_cub_geometry(const char *name) 
{
  FILE *cubfile = fopen(name, "rb");
  if (NULL == cubfile) RETURN(TSTTB_FILE_NOT_FOUND);

    // first get the length
  int result = fseek(cubfile, 0, SEEK_END);
  if (result) {
    TSTTG_processError(TSTTB_FAILURE, "Couldn't seek to end of .cub file.");
    RETURN(TSTTB_FAILURE);
  }
  
  long endpos = ftell(cubfile);
  
  char magic_str[4] = {'\0', '\0', '\0', '\0'};
  fread(magic_str, 1, 4, cubfile);
  if (!strcmp(magic_str, "CUBE")) {
    TSTTG_processError(TSTTB_NOT_SUPPORTED, "Wrong magic string in .cub file.");
    RETURN(TSTTB_NOT_SUPPORTED);
  }
  
    // get the model header
  result = fseek(cubfile, 4, SEEK_SET);
  if (result) {
    TSTTG_processError(TSTTB_FAILURE, "Seek failed reading cub file header.");
    RETURN(TSTTB_FAILURE);
  }
  int header[6];
  fread(header, 4, 6, cubfile);
  int num_models = header[2];
  int model_table_offset = header[3];

    // get the model table
  int model_entries[36], model_offset[6], 
    model_length[6], model_type[6];
  if (num_models > 6) {
    TSTTG_processError(TSTTB_INVALID_ARGUMENT, "Too many models in .cub file.");
    RETURN(TSTTB_INVALID_ARGUMENT);
  }

  if (model_table_offset+24*num_models-1 > endpos) {
    TSTTG_processError(TSTTB_FAILURE, "Reading model table will go past end of file.");
    RETURN(TSTTB_FAILURE);
  }
    
  result = fseek(cubfile, model_table_offset, SEEK_SET);
  if (result) {
    TSTTG_processError(TSTTB_FAILURE, "Seek failed seeking to model table.");
    RETURN(TSTTB_FAILURE);
  }
  fread(model_entries, 4, 6*num_models, cubfile);
  for (int i = 0; i < num_models; i++) {
    model_type[i] = model_entries[6*i];
    model_offset[i] = model_entries[6*i+1];
    model_length[i] = model_entries[6*i+2];
  }
  
    // now read the actual models
  const char *model_type_str[6] = {"", "ACIS_SAT", "", "FACET", "", "GRANITE"};
  for (int i = 0; i < num_models; i++) {
    if ((1 != model_type[i] && 3 != model_type[i] && 5 != model_type[i]) ||
        model_length[i] == 0) continue;

    if (model_offset[i] > endpos) {
      TSTTG_processError(TSTTB_FAILURE, "Reading actual model will go past end of file.");
      RETURN(TSTTB_FAILURE);
    }
    result = fseek(cubfile, model_offset[i], SEEK_SET);
    if (result) {
      TSTTG_processError(TSTTB_FAILURE, "Seek failed seeking to a model.");
      RETURN(TSTTB_FAILURE);
    }
      // if the type is facet, need special handling
    if (3 == model_type[i]) {
        // first create an instance of the facet engine
      FacetQueryEngine::instance();
        // now copy the file contents to a scratch file
      fclose(cubfile);
      cubfile = tmpfile();
      NCubitFile::CCubitFile ccf;
      ccf.Open(name, NULL, NULL);
      ccf.BeginReadGeomModel(3, cubfile);
      ccf.EndReadGeomModel();
    }

    CubitStatus status = gqt->import_temp_geom_file(cubfile, model_type_str[model_type[i]]);
    if (CUBIT_SUCCESS != status) {
      std::string msg = std::string("Trouble loading geometry file of type ") + 
        model_type_str[model_type[i]] + std::string(".");
      TSTTG_processError(TSTTB_FAILURE, msg.c_str());
    }

      // finish up if reading a facet file
    if (3 == model_type[i]) {
      fclose(cubfile);
      cubfile = fopen(name, "rb");
      if (NULL == cubfile) 
        RETURN(TSTTB_SUCCESS);
    }
  }

  if (NULL != cubfile) fclose(cubfile);

  RETURN(TSTTB_SUCCESS);
}
                               
/**
 * Initialize an iterator over gentities of a specified dimension.
 * @param gentity_dimension Dimension of gentities to be iterated over
 * @param gentity_iterator Iterator initialized by this function
 */

enum TSTTB_ErrorType
TSTTG_gentityIteratorInit (TSTTG_Instance instance,
                           /*in*/ const int gentity_dimension,
                           /*out*/ TSTTG_GentityIterator *gentity_iterator)
{
  RETURN(TSTTB_NOT_SUPPORTED);
}

/**
 * Get the next entity for this iterator.
 * @param gentity_iterator Iterator being iterated over
 * @param gentity_handle Next gentity
 * @return If true, there are more gentities, if false, this is the last one
 */
bool
TSTTG_gentityIteratorNext (TSTTG_Instance instance,
                           /*inout*/ TSTTG_GentityIterator *gentity_iterator,
                           /*out*/ TSTTG_GentityHandle *gentity_handle
                           )
{
  TSTTG_processError(TSTTB_NOT_SUPPORTED, "Operation not supported.");
  return false;
}

/**
 * Reset an iterator back to the first gentity
 * @param gentity_iterator Iterator reset by this function
 */
enum TSTTB_ErrorType
TSTTG_gentityIteratorReset(TSTTG_Instance instance,
                           /*inout*/ TSTTG_GentityIterator *gentity_iterator
                           )
{
  RETURN(TSTTB_NOT_SUPPORTED);
}

/**
 * Delete an iterator
 * @param gentity_iterator Iterator deleted by this function
 */
enum TSTTB_ErrorType
TSTTG_gentityIteratorDelete (TSTTG_Instance instance,
                             /*in*/ TSTTG_CGentityIterator Gentity_dim_iterator)
{
  RETURN(TSTTB_NOT_SUPPORTED);
}

/**
 *   Returns true if the gentity_sets are related through a parent/child
 *   relationship.
 */
bool
TSTTG_isChildOf (TSTTG_Instance instance,
                 /*in*/ TSTTG_CGentitySetHandle parent_entity_set,
                 /*in*/ TSTTG_CGentitySetHandle child_entity_set)
{
  std::vector<RefGroup*> *par1, *ch1, *par2, *ch2;
  TM->pc_list(const_cast<RefGroup*>(CONST_SET_HANDLE(parent_entity_set)), par1, ch1, false);
  if (NULL == ch1 || ch1->empty()) return false;
  TM->pc_list(const_cast<RefGroup*>(CONST_SET_HANDLE(child_entity_set)), par2, ch2, false);
  if (NULL == par2 || par2->empty()) return false;
  
  const RefGroup *group1 = CONST_SET_HANDLE(parent_entity_set);
  const RefGroup *group2 = CONST_SET_HANDLE(child_entity_set);
  if ((std::find(ch1->begin(), ch1->end(), group2) != ch1->end())
      || (std::find(par2->begin(), par2->end(), group1) != par2->end()))
    return true;
  
  else return false;
}

/**
 *   Recursively gets the children of this gentity_set up to num_hops
 *   levels; if num_hops is set to -1 all children are returned
 */
enum TSTTB_ErrorType
TSTTG_getChldn (TSTTG_Instance instance,
                /*in*/ TSTTG_CGentitySetHandle from_entity_set,
                /*in*/ const int num_hops,
                /*inout*/ ARRAY_INOUT_DECL(TSTTG_GentitySetHandle, 
                                           entity_set_handles))
{
  std::vector<RefGroup*> group_ptrs;
  const RefGroup *this_grp = CONST_SET_HANDLE(from_entity_set);
  TM->get_pc_groups(const_cast<RefGroup*>(this_grp), 1, num_hops, group_ptrs);
  CHECK_SIZE(*entity_set_handles, TSTTG_GentityHandle, 
             (int)group_ptrs.size());

  std::copy(group_ptrs.begin(), group_ptrs.end(), *entity_set_handles);
  
  RETURN(TSTTB_SUCCESS);
}

/**
 *   Recursively gets the parents of this gentity_set up to num_hops
 *   levels; if num_hops is set to -1 all parents are returned
 */
enum TSTTB_ErrorType
TSTTG_getPrnts (TSTTG_Instance instance,
                /*in*/ TSTTG_CGentitySetHandle from_entity_set,
                /*in*/ const int num_hops,
                /*inout*/ ARRAY_INOUT_DECL(TSTTG_GentitySetHandle, entity_set_handles))
{
  std::vector<RefGroup*> group_ptrs;
  const RefGroup *this_grp = CONST_SET_HANDLE(from_entity_set);
  TM->get_pc_groups(const_cast<RefGroup*>(this_grp), 0, num_hops, group_ptrs);
  CHECK_SIZE(*entity_set_handles, TSTTG_GentityHandle, 
             (int)group_ptrs.size());

  std::copy(group_ptrs.begin(), group_ptrs.end(), *entity_set_handles);
  RETURN(TSTTB_SUCCESS);
}

/**
 *   Returns the number of immediate children in the gentity_set (one
 *   level down only)
 */
int
TSTTG_getNumChld (TSTTG_Instance instance,
                  /*in*/ TSTTG_CGentitySetHandle entity_set,
                  /*in*/ const int num_hops)
{
    // HJK: num_hops has to be handled
  if (1 < num_hops) {
    TSTTG_processError(TSTTB_NOT_SUPPORTED, "Num_hops argument not yet supported.");
    return -1;
  }

  std::vector<RefGroup*> *my_children = TM->pc_list(const_cast<RefGroup*>(CONST_SET_HANDLE(entity_set)), 1, false);
  if (NULL == my_children || my_children->empty()) return 0;
  else return my_children->size();
}

/**
 *   Returns the number of immediate parents to the gentity_set (one
 *   level up only)
 */
int
TSTTG_getNumPrnt (TSTTG_Instance instance,
                  /*in*/ TSTTG_CGentitySetHandle entity_set,
                  /*in*/ const int num_hops)
{
    // HJK: num_hops has to be handled
  if (1 < num_hops) {
    TSTTG_processError(TSTTB_NOT_SUPPORTED, "Num_hops argument not yet supported.");
    return -1;
  }

  std::vector<RefGroup*> *my_parents = TM->pc_list(const_cast<RefGroup*>(CONST_SET_HANDLE(entity_set)), 0, false);
  if (NULL == my_parents || my_parents->empty()) return 0;
  else return my_parents->size();
}

/**
 *   Add a parent to the gentity_set 
 */
enum TSTTB_ErrorType
TSTTG_addPrntChld (TSTTG_Instance instance,
                   /*inout*/ TSTTG_GentitySetHandle *parent_entity_set,
                   /*inout*/ TSTTG_GentitySetHandle *child_entity_set)
{
  std::vector<RefGroup*> *my_parents = 
    TM->pc_list(SET_HANDLE(*child_entity_set), 0, true);
  std::vector<RefGroup*> *my_children = 
    TM->pc_list(SET_HANDLE(*parent_entity_set), 1, true);
  RefGroup *par_group = SET_HANDLE(*parent_entity_set);
  RefGroup *child_group = SET_HANDLE(*child_entity_set);
  my_parents->push_back(par_group);
  my_children->push_back(child_group);
  RETURN(TSTTB_SUCCESS);
}

/**
 *   Remove a parent/child link between gentity_sets
 */
enum TSTTB_ErrorType
TSTTG_rmvPrntChld (TSTTG_Instance instance,
                   /*inout*/ TSTTG_GentitySetHandle *parent_entity_set,
                   /*inout*/ TSTTG_GentitySetHandle *child_entity_set)
{
  RefGroup *parent = reinterpret_cast<RefGroup *>(*parent_entity_set);
  RefGroup *child = reinterpret_cast<RefGroup *>(*child_entity_set);
  std::vector<RefGroup*> *children = TM->pc_list(parent, 1, false),
    *parents = TM->pc_list(child, 0, false);
  if (NULL == children || NULL == parents) {
    RETURN(TSTTB_INVALID_ARGUMENT);
  }
  
  children->erase(std::remove(children->begin(), children->end(), child), children->end());
  parents->erase(std::remove(parents->begin(), parents->end(), parent), parents->end());
  RETURN(TSTTB_SUCCESS);
}

/**
 * Return gentities of specified dimension in this set, or in whole model.
 * @param set_handle Entity set being queried (if 0, whole model)
 * @param gentity_dimension Dimension of entities being queried
 * @param gentity_handles Gentity handles
 */
enum TSTTB_ErrorType
TSTTG_gentitysetGetGentitiesOfType (TSTTG_Instance instance,
                                    /*in*/ TSTTG_CGentitySetHandle set_handle,
                                    /*in*/ const TSTTG_GentityType gentity_type,
                                    /*out*/ ARRAY_INOUT_DECL(TSTTG_GentityHandle, gentity_handles))
{
  const RefGroup *this_set = CONST_SET_HANDLE(set_handle);
  static DLIList<RefEntity*> dim_entities;
  dim_entities.clean_out();
  if (0 == this_set) {
    if (gentity_type == TSTTG_ALL_TYPES) {
      for (int i = 0; i < 4; i++) {
        static DLIList<RefEntity*> temp_entities;
        RefEntityFactory::instance()->ref_entity_list(TSTTG_entity_type_names[i], 
                                                      temp_entities, CUBIT_FALSE);
        dim_entities += temp_entities;
      }
    }
    else {
      RefEntityFactory::instance()->ref_entity_list(TSTTG_entity_type_names[gentity_type], 
                                                    dim_entities, CUBIT_FALSE);
    }
  }
  else {
    static DLIList<CubitEntity*> centities;
    static DLIList<RefEntity*> entities;
    centities.clean_out();
    entities.clean_out();
    const_cast<RefGroup*>(this_set)->get_child_entities(centities);
    CAST_LIST(centities, entities, RefEntity);
    if (gentity_type == 4) {
      dim_entities += entities;
    }
    else {
      entities.reset();
      int num_ents = entities.size();
      for (int i = 0; i < num_ents; i++) {
        if (entities.get()->dimension() == gentity_type)
          dim_entities.append(entities.get_and_step());
        else
          entities.step();
      }
    }
  }
    
  CHECK_SIZE(*gentity_handles, TSTTG_GentityHandle, dim_entities.size());
  dim_entities.copy_to((RefEntity**)*gentity_handles);
  RETURN(TSTTB_SUCCESS);
}

/**
 * Return number of gentities of specified dimension in this set, or in
 * whole model.
 * @param set_handle Entity set being queried (if 0, whole model)
 * @param gentity_dimension Dimension of entities being queried
 * @return Number of entities
 */
int
TSTTG_gentitysetGetNumberGentitiesOfType (TSTTG_Instance instance,
                                          /*in*/ TSTTG_CGentitySetHandle set_handle,
                                          /*in*/ const TSTTG_GentityType gentity_type)
{
  const RefGroup *this_set = CONST_SET_HANDLE(set_handle);
  static DLIList<RefEntity*> *dim_entities = NULL;
  if (NULL == dim_entities) dim_entities = new DLIList<RefEntity*>();
  dim_entities->clean_out();
  if (0 == this_set) {
    RefEntityFactory::instance()->ref_entity_list(TSTTG_entity_type_names[gentity_type], 
                                                  *dim_entities, CUBIT_FALSE);
  }
  else {
    static DLIList<CubitEntity*> centities;
    static DLIList<RefEntity*> entities;
    centities.clean_out();
    entities.clean_out();
    const_cast<RefGroup*>(this_set)->get_child_entities(centities);
    CAST_LIST(centities, entities, RefEntity);
    entities.reset();
    int num_ents = entities.size();
    for (int i = 0; i < num_ents; i++) {
      if (entities.get()->dimension() == gentity_type)
        dim_entities->append(entities.get_and_step());
      else
        entities.step();
    }
  }
    
  return dim_entities->size();
}

/**
 *    Returns an integer array of topological dimensions for an input
 *    array of entity handles.
 */
enum TSTTB_ErrorType
TSTTG_gentityGetType (TSTTG_Instance instance,
                      /*in*/ ARRAY_IN_DECL(TSTTG_CGentityHandle, gentity_handles),
                      /*inout*/ ARRAY_INOUT_DECL(TSTTG_GentityType, gtype))
{
    // go through each entity and look up its dimension
  CHECK_SIZE(*gtype, TSTTG_GentityType, gentity_handles_size);

  const RefEntity **tmp_handles = CONST_ENTITY_HANDLE_ARRAY(gentity_handles);
  
  for (int i = 0; i < gentity_handles_size; i++) {
    (*gtype)[i] = (TSTTG_GentityType) tmp_handles[i]->dimension();
  }

  RETURN(TSTTB_SUCCESS);
}

/**
 * Get the adjacent entities of a given dimension.
 * @param gentity_handle Entity for which adjacencies are requested
 * @param to_dimension Target dimension of adjacent entities
 * @param adj_gentities List returned with adjacent entities
 */
enum TSTTB_ErrorType
TSTTG_gentityGetAdjacencies (TSTTG_Instance instance,
                             /*in*/ TSTTG_CGentityHandle gentity_handle,
                             /*in*/ const int to_dimension,
                             /*inout*/ ARRAY_INOUT_DECL(TSTTG_GentityHandle, adj_gentities))
{
  const RefEntity *tmp_hndl = CONST_ENTITY_HANDLE(gentity_handle);

  if (tmp_hndl->dimension() == to_dimension) {
    TSTTG_processError(TSTTB_INVALID_ARGUMENT, "Can't get adjacencies to entities of same dimension.");
    return TSTTB_INVALID_ARGUMENT;
  }

  static DLIList<RefEntity*> tmp_ents;
  TSTTG_get_adjacent_entities(tmp_hndl, to_dimension, tmp_ents);

  CHECK_SIZE(*adj_gentities, TSTTG_GentityHandle, (int) tmp_ents.size());
  tmp_ents.copy_to((RefEntity**)*adj_gentities);
  RETURN(TSTTB_SUCCESS);
}

void TSTTG_get_adjacent_entities(const RefEntity *from, const int to_dim,
                                 DLIList<RefEntity*> &adj_ents) 
{
  TopologyEntity *topo_ent = const_cast<TopologyEntity*>(dynamic_cast<const TopologyEntity*>(from));
  if (NULL == topo_ent) return;
  
  adj_ents.clean_out();
  static DLIList<RefVertex*> tmp_verts;
  static DLIList<RefEdge*> tmp_edges;
  static DLIList<RefFace*> tmp_faces;
  static DLIList<RefVolume*> tmp_volumes;
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
      tmp_volumes.clean_out();
      topo_ent->ref_volumes(tmp_volumes);
      CAST_LIST_TO_PARENT(tmp_volumes, adj_ents);
      break;
  }
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
enum TSTTB_ErrorType
TSTTG_gentityGet2OAdjacencies (TSTTG_Instance instance,
                               /*in*/ TSTTG_CGentityHandle gentity_handle,
                               /*in*/ const int bridge_dimension,
                               /*in*/ const int to_dimension,
                               /*out*/ ARRAY_INOUT_DECL(TSTTG_GentityHandle, adjacent_gentities))
{
    // for better or worse, go to cgm as quickly as possible, to avoid working with
    // sidl arrays
  const RefEntity *gentity = CONST_ENTITY_HANDLE(gentity_handle);
  DLIList<RefEntity*> to_ents, bridge_ents, tmp_ents;

  TSTTG_get_adjacent_entities(gentity, bridge_dimension, bridge_ents);
  
  for (int i = bridge_ents.size(); i > 0; i--) {
    TSTTG_get_adjacent_entities(bridge_ents.get_and_step(), to_dimension, tmp_ents);
    to_ents += tmp_ents;
    tmp_ents.clean_out();
  }
  
  to_ents.uniquify_unordered();

  CHECK_SIZE(*adjacent_gentities, TSTTG_GentityHandle, (int) to_ents.size());
  to_ents.copy_to((RefEntity**)*adjacent_gentities);

  RETURN(TSTTB_SUCCESS);
}

/**
 * Return whether or not entities are adjacent.
 * @param gentity_handle1 1st entity
 * @param gentity_handle2 2nd entity
 * @param are_adjacent If true, entities are adjacent
 */
enum TSTTB_ErrorType
TSTTG_gentityIsAdjacent (TSTTG_Instance instance,
                         /*in*/ TSTTG_CGentityHandle gentity_handle1,
                         /*in*/ TSTTG_CGentityHandle gentity_handle2,
                         /*out*/ bool *are_adjacent)
{
  const TopologyEntity *ent1 = dynamic_cast<const TopologyEntity*>(CONST_ENTITY_HANDLE(gentity_handle1));
  const TopologyEntity *ent2 = dynamic_cast<const TopologyEntity*>(CONST_ENTITY_HANDLE(gentity_handle2));
  if (ent1 != NULL && ent2 != NULL)
    *are_adjacent = const_cast<TopologyEntity*>(ent1)->is_directly_related(const_cast<TopologyEntity*>(ent2));
  else *are_adjacent = false;
  RETURN(TSTTB_SUCCESS);
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
enum TSTTB_ErrorType
TSTTG_createTag (TSTTG_Instance instance,
                 /*in*/ const char *tag_name,
                 /*in*/ const int tag_size,
                 /*in*/ const TSTTG_TagValueType tag_type,
                 /*out*/ TSTTG_TagHandle *tag_handle)
{
  int new_tag;
  int this_size = tag_size;
  switch (tag_type) {
    case TSTTG_INTEGER:
      this_size *= sizeof(int);
      break;
    case TSTTG_DOUBLE:
      this_size *= sizeof(double);
      break;
    case TSTTG_ENTITY_HANDLE:
      this_size *= sizeof(TSTTG_GentityHandle);
      break;
    case TSTTG_BYTES:
      break;
  }
  TM->createTag(tag_name, this_size, tag_type, NULL, &new_tag);
  *tag_handle = (TSTTG_TagHandle) new_tag;
  RETURN(TSTTB_SUCCESS);
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
enum TSTTB_ErrorType
TSTTG_destroyTag (TSTTG_Instance instance,
                  /*in*/ const TSTTG_TagHandle tag_handle,
                  /*in*/ const bool forced)
{
  TM->destroyTag(CONST_TAG_HANDLE(tag_handle), forced);
  RETURN(TSTTB_SUCCESS);
}

/**
 *   Get the tag name associated with a given tag handle.
 */
const char *
TSTTG_getTagName (TSTTG_Instance instance,
                  /*in*/ const TSTTG_TagHandle tag_handle)
{
  return TM->getTagName(CONST_TAG_HANDLE(tag_handle));
}

/**
 *   Get the size of the data associated with a given tag handle.
 */
int
TSTTG_getTagSizeValues (TSTTG_Instance instance,
                        /*in*/ const TSTTG_TagHandle tag_handle)
{
  int this_size = TM->getTagSize(CONST_TAG_HANDLE(tag_handle));
  switch (TSTTG_getTagType(instance, tag_handle)) {
    case TSTTG_INTEGER:
      this_size /= sizeof(int);
      break;
    case TSTTG_DOUBLE:
      this_size /= sizeof(double);
      break;
    case TSTTG_ENTITY_HANDLE:
      this_size /= sizeof(TSTTG_GentityHandle);
      break;
    case TSTTG_BYTES:
      break;
  }
  return this_size;
}

/**
 *   Get the size of the data associated with a given tag handle.
 */
int
TSTTG_getTagSizeBytes (TSTTG_Instance instance,
                       /*in*/ const TSTTG_TagHandle tag_handle)
{
  return TM->getTagSize(CONST_TAG_HANDLE(tag_handle));
}

TSTTG_TagValueType
TSTTG_getTagType (TSTTG_Instance instance,
                    /*in*/ const TSTTG_TagHandle tag_handle)
{
  return (TSTTG_TagValueType) TM->getTagType(CONST_TAG_HANDLE(tag_handle));
}

/**
 *     Get the tag handle associated with a given string name.
 */
TSTTG_TagHandle
TSTTG_getTagHandle (TSTTG_Instance instance,
                    /*in*/ const char *tag_name)
{
  return (TSTTG_TagHandle) TM->getTagHandle(tag_name);
}

enum TSTTB_ErrorType
TSTTG_rmvArrTag (TSTTG_Instance instance,
                 /*in*/ ARRAY_IN_DECL(TSTTG_GentityHandle, entity_handles),
                 /*in*/ const TSTTG_TagHandle tag_handle) 
{
  RefEntity **tmp_entity_handles = ENTITY_HANDLE_ARRAY(entity_handles);  
  TM->rmvArrTag(tmp_entity_handles, entity_handles_size, CONST_TAG_HANDLE(tag_handle));
  RETURN(TSTTB_SUCCESS);
}

/**
 *   Allows the user to disassociate the tag referenced by the tag
 *   handle from the specified gentities. The tag data is not deleted in
 *   this call, but can be deleted later using the deleteTag function
 *   defined above.
 */
enum TSTTB_ErrorType
TSTTG_rmvTag (TSTTG_Instance instance,
              /*in*/ TSTTG_GentityHandle entity_handle,
              /*in*/ const TSTTG_TagHandle tag_handle)
{
  RefEntity *tmp_entity = ENTITY_HANDLE(entity_handle);
  TM->rmvArrTag(&tmp_entity, 1, CONST_TAG_HANDLE(tag_handle));
  RETURN(TSTTB_SUCCESS);
}

/**
 *   Get all tag handles associated with a given gentity.
 */
enum TSTTB_ErrorType
TSTTG_getAllTags (TSTTG_Instance instance,
                  /*in*/ TSTTG_CGentityHandle entity_handle,
                  /*inout*/ ARRAY_INOUT_DECL(TSTTG_TagHandle, tag_handles))
{
  TM->getAllTags(CONST_SET_HANDLE(entity_handle), 
                 TAG_HANDLE_ARRAY_INOUT(tag_handles));
  RETURN(TSTTB_SUCCESS);
}

enum TSTTB_ErrorType
TSTTG_getArrData (TSTTG_Instance instance,
                    /*in*/ ARRAY_IN_DECL(TSTTG_CGentityHandle, entity_handles),
                    /*in*/ const TSTTG_TagHandle tag_handle,
                    /*inout*/ ARRAY_INOUT_DECL(char, tag_value))
{
  const RefEntity **tmp_entity_handles = CONST_ENTITY_HANDLE_ARRAY(entity_handles);  
  TM->getArrData(tmp_entity_handles, entity_handles_size,
                 CONST_TAG_HANDLE(tag_handle), 
                 ARRAY_INOUT(tag_value));
  RETURN(TSTTB_SUCCESS);
}


enum TSTTB_ErrorType
TSTTG_getIntArrData (TSTTG_Instance instance,
                       /*in*/ ARRAY_IN_DECL(TSTTG_CGentityHandle, entity_handles),
                       /*in*/ const TSTTG_TagHandle tag_handle,
                       /*inout*/ ARRAY_INOUT_DECL(int, tag_value))
{
  *tag_value_allocated *= sizeof(int);
  *tag_value_size *= sizeof(int);
  TSTTB_ErrorType retval = TSTTG_getArrData(instance, entity_handles, 
                                           entity_handles_size, tag_handle,
                                           reinterpret_cast<char**>(tag_value), 
                                           tag_value_allocated, 
                                           tag_value_size);
  *tag_value_allocated /= sizeof(int);
  *tag_value_size /= sizeof(int);
  return retval;
}


enum TSTTB_ErrorType
TSTTG_getDblArrData (TSTTG_Instance instance,
                       /*in*/ ARRAY_IN_DECL(TSTTG_CGentityHandle, entity_handles),
                       /*in*/ const TSTTG_TagHandle tag_handle,
                       /*inout*/ ARRAY_INOUT_DECL(double, tag_value))
{
  *tag_value_allocated *= sizeof(double);
  *tag_value_size *= sizeof(double);
  TSTTB_ErrorType retval = TSTTG_getArrData(instance, entity_handles, 
                                           entity_handles_size, tag_handle,
                                           reinterpret_cast<char**>(tag_value), 
                                           tag_value_allocated, tag_value_size);
  *tag_value_allocated /= sizeof(double);
  *tag_value_size /= sizeof(double);
  return retval;
}

enum TSTTB_ErrorType
TSTTG_getEHArrData (TSTTG_Instance instance,
                      /*in*/ ARRAY_IN_DECL(TSTTG_CGentityHandle, entity_handles),
                      /*in*/ const TSTTG_TagHandle tag_handle,
                      /*inout*/ ARRAY_INOUT_DECL(void*, tag_value))
{
  *tag_value_allocated *= sizeof(TSTTG_GentityHandle);
  *tag_value_size *= sizeof(TSTTG_GentityHandle);
  TSTTB_ErrorType retval = TSTTG_getArrData(instance, entity_handles, 
                                           entity_handles_size, tag_handle,
                                           reinterpret_cast<char**>(tag_value), 
                                           tag_value_allocated, 
                                           tag_value_size);
  *tag_value_allocated /= sizeof(TSTTG_GentityHandle);
  *tag_value_size /= sizeof(TSTTG_GentityHandle);
  return retval;
}

enum TSTTB_ErrorType
TSTTG_setArrData (TSTTG_Instance instance,
                    /*in*/ ARRAY_IN_DECL(TSTTG_GentityHandle, entity_handles),
                    /*in*/ const TSTTG_TagHandle tag_handle,
                    /*in*/ CONST_ARRAY_IN_DECL(char, tag_value))
{
  CubitStatus status = TM->setArrData(ENTITY_HANDLE_ARRAY(entity_handles),
                                      entity_handles_size,
                                      TAG_HANDLE(tag_handle), 
                                      ARRAY_IN(tag_value));

  if (CUBIT_SUCCESS == status) {
    RETURN(TSTTB_SUCCESS);
  }
  else {
    RETURN(TSTTB_FAILURE);
  }
}

enum TSTTB_ErrorType
TSTTG_setIntArrData (TSTTG_Instance instance,
                       /*in*/ ARRAY_IN_DECL(TSTTG_GentityHandle, entity_handles),
                       /*in*/ const TSTTG_TagHandle tag_handle,
                       /*in*/ CONST_ARRAY_IN_DECL(int, tag_value))
{
  return TSTTG_setArrData(instance, entity_handles, 
                          entity_handles_size, tag_handle, 
                          reinterpret_cast<const char*>(tag_value), 
                          sizeof(int)*tag_value_size);
}

enum TSTTB_ErrorType
TSTTG_setDblArrData (TSTTG_Instance instance,
                       /*in*/ ARRAY_IN_DECL(TSTTG_GentityHandle, entity_handles),
                       /*in*/ const TSTTG_TagHandle tag_handle,
                       /*in*/ CONST_ARRAY_IN_DECL(double, tag_value))
{
  return TSTTG_setArrData(instance, entity_handles, 
                   entity_handles_size, tag_handle, 
                          reinterpret_cast<const char*>(tag_value), 
                          sizeof(double)*tag_value_size);
}

enum TSTTB_ErrorType
TSTTG_setEHArrData (TSTTG_Instance instance,
                      /*in*/ ARRAY_IN_DECL(TSTTG_GentityHandle, entity_handles),
                      /*in*/ const TSTTG_TagHandle tag_handle,
                      /*in*/ ARRAY_IN_DECL(TSTTG_CGentityHandle, tag_values))
{
  return TSTTG_setArrData(instance, entity_handles, 
                   entity_handles_size, tag_handle, 
                          reinterpret_cast<const char*>(tag_values), 
                          sizeof(TSTTG_GentityHandle)*tag_values_size);
}

/**
 *   Allows the user to retrieve an array of tag values associated with
 *   a tag handle from an input array of gentity handles
 */
enum TSTTB_ErrorType
TSTTG_getData (TSTTG_Instance instance,
               /*in*/ TSTTG_CGentityHandle entity_handle,
               /*in*/ const TSTTG_TagHandle tag_handle,
               /*inout*/ ARRAY_INOUT_DECL(char, tag_value))
{
  const RefEntity *tmp_entity = CONST_ENTITY_HANDLE(entity_handle);
  TM->getArrData(&tmp_entity, 1, CONST_TAG_HANDLE(tag_handle), 
                 ARRAY_INOUT(tag_value));
  RETURN(TSTTB_SUCCESS);
}

int
TSTTG_getIntData (TSTTG_Instance instance,
                  /*in*/ TSTTG_CGentityHandle entity_handle,
                  /*in*/ const TSTTG_TagHandle tag_handle ) 
{
  int val;
  char *val_ptr = reinterpret_cast<char*>(&val);
  int val_size = sizeof(int);
  TSTTG_getArrData(instance, &entity_handle, 1, 
                   tag_handle, &val_ptr, &val_size, &val_size);
  return val;
}

double
TSTTG_getDblData (TSTTG_Instance instance,
                  /*in*/ TSTTG_CGentityHandle entity_handle,
                  /*in*/ const TSTTG_TagHandle tag_handle ) 
{
  double val;
  char *val_ptr = reinterpret_cast<char*>(&val);
  int val_size = sizeof(double);
  TSTTG_getArrData(instance, &entity_handle, 1, 
                   tag_handle, &val_ptr, &val_size, &val_size);
  return val;
}

TSTTG_GentityHandle
TSTTG_getEHData (TSTTG_Instance instance,
                 /*in*/ TSTTG_CGentityHandle entity_handle,
                 /*in*/ const TSTTG_TagHandle tag_handle ) 
{
  TSTTG_GentityHandle val;
  char *val_ptr = reinterpret_cast<char*>(&val);
  int dum = sizeof(TSTTG_GentityHandle);
  TSTTG_getArrData(instance, &entity_handle, 1, 
                   tag_handle, &val_ptr, &dum, &dum);
  return val;
}

/**
 *   Allows the user to set the tag data values on an array of gentity
 *   handles
 */
enum TSTTB_ErrorType
TSTTG_setData (TSTTG_Instance instance,
               /*in*/ TSTTG_GentityHandle entity_handle,
               /*in*/ const TSTTG_TagHandle tag_handle,
               /*in*/ CONST_ARRAY_IN_DECL(char, tag_value))
{
  RefEntity *tmp_entity = ENTITY_HANDLE(entity_handle);
  TM->setArrData(&tmp_entity, 1, 
                 CONST_TAG_HANDLE(tag_handle), 
                 ARRAY_IN(tag_value));
  RETURN(TSTTB_SUCCESS);
}

TSTTB_ErrorType
TSTTG_setIntData (TSTTG_Instance instance,
                  /*in*/ TSTTG_GentityHandle entity_handle,
                  /*in*/ const TSTTG_TagHandle tag_handle,
                  /*in*/ const int tag_value ) 
{
  return TSTTG_setArrData(instance, &entity_handle, 1, 
                          tag_handle, 
                          reinterpret_cast<const char*>(&tag_value), 
                          sizeof(int));
}

TSTTB_ErrorType
TSTTG_setDblData (TSTTG_Instance instance,
                   
                  /*in*/ TSTTG_GentityHandle entity_handle,
                  /*in*/ const TSTTG_TagHandle tag_handle,
                  /*in*/ const double tag_value ) 
{
  return TSTTG_setArrData(instance, &entity_handle, 1, 
                          tag_handle, 
                          reinterpret_cast<const char*>(&tag_value), 
                          sizeof(double));
}

TSTTB_ErrorType
TSTTG_setEHData (TSTTG_Instance instance,
                 /*in*/ TSTTG_GentityHandle entity_handle,
                 /*in*/ const TSTTG_TagHandle tag_handle,
                 /*in*/ TSTTG_CGentityHandle tag_value ) 
{
  return TSTTG_setArrData(instance, &entity_handle, 1, 
                          tag_handle, reinterpret_cast<const char*>(&tag_value), 
                          sizeof(TSTTG_GentityHandle));
}

/**
 *   Remove the tag associated with the tag_handle from the mesh or
 *   gentity_set.  The tag data is not destroyed in this function, but
 *   can be destroyed using the deleteTag function.
 */
enum TSTTB_ErrorType
TSTTG_rmvEntSetTag (TSTTG_Instance instance,
                    /*in*/ TSTTG_GentitySetHandle entity_set,
                    /*in*/ const TSTTG_TagHandle tag_handle)
{
    // have to go through RefEntity* so that RefEntity** gets set right
  RefEntity *tmp_entity = SET_HANDLE(entity_set);
  TM->rmvArrTag(&tmp_entity, 1, CONST_TAG_HANDLE(tag_handle));
  RETURN(TSTTB_SUCCESS);
}

/**
 *   Get all tag handles associated with a given mesh or gentity_set.
 */
enum TSTTB_ErrorType
TSTTG_getAllEntSetTags (TSTTG_Instance instance,
                        /*in*/ TSTTG_CGentitySetHandle entity_set,
                        /*inout*/ ARRAY_INOUT_DECL(TSTTG_TagHandle, tag_handles))
{
    // have to go through RefEntity* so that RefEntity** gets set right
  const RefEntity *tmp_entity = CONST_SET_HANDLE(entity_set);
  TM->getAllTags(tmp_entity,TAG_HANDLE_ARRAY_INOUT(tag_handles));
  RETURN(TSTTB_SUCCESS);
}

/**
 *   Get the tag data associated with a tag handle from the mesh or
 *   gentity_set.  It is assumed that the tag_value argument is
 *   allocated by the application before being passed into the getTag
 *   function.
 */
enum TSTTB_ErrorType
TSTTG_getEntSetData (TSTTG_Instance instance,
                     /*in*/ TSTTG_CGentitySetHandle entity_set,
                     /*in*/ const TSTTG_TagHandle tag_handle,
                     /*inout*/ ARRAY_INOUT_DECL(char, tag_value))
{
    // have to go through RefEntity* so that RefEntity** gets set right
  const RefEntity *tmp_entity = CONST_SET_HANDLE(entity_set);
  TM->getArrData(&tmp_entity, 1, CONST_TAG_HANDLE(tag_handle), 
                 ARRAY_INOUT(tag_value));
  RETURN(TSTTB_SUCCESS);
}

int
TSTTG_getEntSetIntData (TSTTG_Instance instance,
                        /*in*/ TSTTG_CGentityHandle entity_set,
                        /*in*/ const TSTTG_TagHandle tag_handle ) 
{
  int tag_value;
  char *tag_ptr = reinterpret_cast<char*>(&tag_value);
  int tag_size = 4;
  TSTTG_getEntSetData(instance, entity_set, tag_handle, &tag_ptr, 
                      &tag_size, &tag_size);
  return tag_value;
}

double
TSTTG_getEntSetDblData (TSTTG_Instance instance,
                        /*in*/ TSTTG_CGentityHandle entity_set,
                        /*in*/ const TSTTG_TagHandle tag_handle ) 
{
  double tag_value;
  char *tag_ptr = reinterpret_cast<char*>(&tag_value);
  int tag_size = 8;
  TSTTG_getEntSetData(instance, entity_set, tag_handle, &tag_ptr, 
                      &tag_size, &tag_size);
  return tag_value;
}

TSTTG_GentityHandle
TSTTG_getEntSetEHData (TSTTG_Instance instance,
                       /*in*/ TSTTG_CGentityHandle entity_set,
                       /*in*/ const TSTTG_TagHandle tag_handle ) 
{
  TSTTG_GentityHandle tag_value;
  char* tag_ptr = reinterpret_cast<char*>(&tag_value);
  int tag_size;
  TSTTG_getEntSetData(instance, entity_set, tag_handle, &tag_ptr, 
                      &tag_size, &tag_size);
  return tag_value;
}

/**
 *   Set the tag data associated with a given tag handle on the mesh or
 *   gentity_set
 */
enum TSTTB_ErrorType
TSTTG_setEntSetData (TSTTG_Instance instance,
                     /*in*/ TSTTG_GentitySetHandle entity_set,
                     /*in*/ const TSTTG_TagHandle tag_handle,
                     /*in*/ CONST_ARRAY_IN_DECL(char, tag_value))
{
    // have to go through RefEntity* so that RefEntity** gets set right
  RefEntity *tmp_entity = SET_HANDLE(entity_set);
  TM->setArrData(&tmp_entity, 1, CONST_TAG_HANDLE(tag_handle), 
                 ARRAY_IN(tag_value));
  RETURN(TSTTB_SUCCESS);
}

TSTTB_ErrorType
TSTTG_setEntSetIntData (TSTTG_Instance instance,
                          /*in*/ TSTTG_GentityHandle entity_set,
                          /*in*/ const TSTTG_TagHandle tag_handle,
                          /*in*/ const int tag_value ) 
{
  return TSTTG_setEntSetData(instance, entity_set, tag_handle, 
                             reinterpret_cast<const char*>(&tag_value), 
                             sizeof(int));
}

TSTTB_ErrorType
TSTTG_setEntSetDblData (TSTTG_Instance instance,
                          /*in*/ TSTTG_GentityHandle entity_set,
                          /*in*/ const TSTTG_TagHandle tag_handle,
                          /*in*/ const double tag_value ) 
{
  return TSTTG_setEntSetData(instance, entity_set, tag_handle, 
                             reinterpret_cast<const char*>(&tag_value),
                             sizeof(double));
}

TSTTB_ErrorType
TSTTG_setEntSetEHData (TSTTG_Instance instance,
                         /*in*/ TSTTG_GentitySetHandle entity_set,
                         /*in*/ const TSTTG_TagHandle tag_handle,
                         /*in*/ TSTTG_CGentityHandle tag_value ) 
{
  return TSTTG_setEntSetData(instance, entity_set, tag_handle, 
                             reinterpret_cast<const char*>(&tag_value), 
                             sizeof(TSTTG_GentityHandle));
}

/**
 * Return a points on specified entities closest to specified points
 * in space.  Input coordinates and output points are interleaved in 
 * the arrays.
 * @param gentity_handles The gentities being queried
 * @param near_coordinates Input coordinates
 * @param on_coordinates Closest point on gentity
 */
enum TSTTB_ErrorType
TSTTG_gentityClosestPoint (TSTTG_Instance instance,
                           /*in*/ ARRAY_IN_DECL(TSTTG_CGentityHandle, gentity_handles),
                           /*in*/ CONST_ARRAY_IN_DECL(double, near_coordinates),
                           /*out*/ ARRAY_INOUT_DECL(double, on_coordinates))
{
  const RefEntity **handle_array = CONST_ENTITY_HANDLE_ARRAY(gentity_handles);

    // check or pre-allocate the coordinate arrays
  CHECK_SIZE(*on_coordinates, double, 3*gentity_handles_size);
  
  CubitVector dumvec;
  const RefEdge *this_edge;
  const RefFace *this_face;
  Surface *this_surf;
  CubitStatus status = CUBIT_SUCCESS, temp_status;

  const double *near_vec = near_coordinates;
  double *on_vec = *on_coordinates;

  for (int i = 0; i < gentity_handles_size; i++) {
    const RefEntity *this_entity = handle_array[i];
  
    switch (this_entity->dimension()) {
      case 0:
        dynamic_cast<const RefVertex*>(this_entity)->coordinates().get_xyz(&on_vec[3*i]);
        break;
      case 1:
        this_edge = dynamic_cast<const RefEdge*>(this_entity);
        if (NULL == this_edge) status = CUBIT_FAILURE;
        else {
          if (debug) {
            std::cout << "Edge " << this_edge->id() << " closest point to (" 
                      << near_vec[3*i] << ", " << near_vec[3*i+1] << ", " << near_vec[3*i+2] 
                      << ") is ";
          }
        
          temp_status = 
            const_cast<RefEdge*>(this_edge)->closest_point(CubitVector(&near_vec[3*i]), dumvec);
          if (temp_status != CUBIT_SUCCESS) status = temp_status;
          else dumvec.get_xyz(&on_vec[3*i]);
        }
        if (debug) {
          std::cout << on_vec[3*i] << ", " << on_vec[3*i+1] << ", " << on_vec[3*i+2] 
                    << ")" << std::endl;
        }
        break;
      case 2:
        this_face = dynamic_cast<const RefFace*>(this_entity);
        if (NULL != this_face) this_surf = const_cast<RefFace*>(this_face)->get_surface_ptr();
      
        if (NULL == this_face) status = CUBIT_FAILURE;
        else {
          temp_status = 
            this_surf->closest_point(CubitVector(&near_vec[3*i]), &dumvec);
          if (temp_status != CUBIT_SUCCESS) status = temp_status;
          else dumvec.get_xyz(&on_vec[3*i]);
        }
        break;
      case 3:
          // just copy over the coordinates
        on_vec[3*i] = near_vec[3*i];
        on_vec[3*i+1] = near_vec[3*i+1];
        on_vec[3*i+2] = near_vec[3*i+2];
        break;
    }
  }
  
  if (status == CUBIT_FAILURE) {
    TSTTG_processError(TSTTB_FAILURE, "Problems getting closest point for some entity.");
    RETURN(TSTTB_FAILURE);
  }

  RETURN(TSTTB_SUCCESS);
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
enum TSTTB_ErrorType
TSTTG_gentityClosestPointAndNormal (TSTTG_Instance instance,
                                    /*in*/ ARRAY_IN_DECL(TSTTG_CGentityHandle, gentity_handles),
                                    /*in*/ CONST_ARRAY_IN_DECL(double, near_coordinates),
                                    /*out*/ ARRAY_INOUT_DECL(double, on_coordinates),
                                    /*out*/ ARRAY_INOUT_DECL(double, normals))
{
  const RefEntity **handle_array = CONST_ENTITY_HANDLE_ARRAY(gentity_handles);

    // check or pre-allocate the coordinate arrays
  CHECK_SIZE(*on_coordinates, double, 3*gentity_handles_size);
  CHECK_SIZE(*normals, double, 3*gentity_handles_size);
  
  CubitVector dumvec, dumnormal;
  const RefEdge *this_edge;
  const RefFace *this_face;
  Surface *this_surf;
  CubitStatus status = CUBIT_SUCCESS, temp_status;

  const double *near_vec = near_coordinates;
  double *on_vec = *on_coordinates;
  double *normals_vec = *normals;

  for (int i = 0; i < gentity_handles_size; i++) {
    const RefEntity *this_entity = handle_array[i];
  
    switch (this_entity->dimension()) {
      case 0:
        dynamic_cast<const RefVertex*>(this_entity)->coordinates().get_xyz(&on_vec[3*i]);
        normals_vec[3*i] = normals_vec[3*i+1] = normals_vec[3*i+2] = 0.0;
        break;
      case 1:
        this_edge = dynamic_cast<const RefEdge*>(this_entity);
        if (NULL == this_edge) status = CUBIT_FAILURE;
        else {
          if (debug) {
            std::cout << "Edge " << this_edge->id() << " closest point to (" 
                      << near_vec[3*i] << ", " << near_vec[3*i+1] << ", " << near_vec[3*i+2] 
                      << ") is ";
          }
        
          temp_status = 
            const_cast<RefEdge*>(this_edge)->closest_point(CubitVector(&near_vec[3*i]), dumvec);
          if (temp_status != CUBIT_SUCCESS) status = temp_status;
          else dumvec.get_xyz(&on_vec[3*i]);
          normals_vec[3*i] = normals_vec[3*i+1] = normals_vec[3*i+2] = 0.0;
        }
        if (debug) {
          std::cout << on_vec[3*i] << ", " << on_vec[3*i+1] << ", " << on_vec[3*i+2] 
                    << ")" << std::endl;
        }
        break;
      case 2:
        this_face = dynamic_cast<const RefFace*>(this_entity);
      
        if (NULL == this_face) 
          status = CUBIT_FAILURE;
        else {
          this_surf = const_cast<RefFace*>(this_face)->get_surface_ptr();
          temp_status = 
            this_surf->closest_point(CubitVector(&near_vec[3*i]), &dumvec, &dumnormal);
          if (temp_status != CUBIT_SUCCESS) 
            status = temp_status;
          else {
            if (this_surf->bridge_sense() == CUBIT_REVERSED)
              dumnormal = -dumnormal;
            dumvec.get_xyz(&on_vec[3*i]);
            dumnormal.get_xyz(&normals_vec[3*i]);
          }
        }
        break;
      case 3:
          // just copy over the coordinates
        on_vec[3*i] = near_vec[3*i];
        on_vec[3*i+1] = near_vec[3*i+1];
        on_vec[3*i+2] = near_vec[3*i+2];
        normals_vec[3*i] = normals_vec[3*i+1] = normals_vec[3*i+2] = 0.0;
        break;
    }
  }
  
  if (status == CUBIT_FAILURE) {
    TSTTG_processError(TSTTB_FAILURE, "Problems getting closest point for some entity.");
    RETURN(TSTTB_FAILURE);
  }

  RETURN(TSTTB_SUCCESS);
}

/**
 * Return the normals at point on specified entities.  Returns error
 * if any input entity is not a gface.  Input coordinates and normals
 * are interleaved in the arrays.
 * @param gentity_handles The gentities being queried
 * @param coordinates Input coordinates, interleaved
 * @param normals The normals at the specified points, interleaved
 */
enum TSTTB_ErrorType
TSTTG_gentityNormal (TSTTG_Instance instance,
                     /*in*/ ARRAY_IN_DECL(TSTTG_CGentityHandle, gentity_handles),
                     /*in*/ CONST_ARRAY_IN_DECL(double, coordinates),
                     /*out*/ ARRAY_INOUT_DECL(double, normals))
{
  const RefEntity **handle_array = CONST_ENTITY_HANDLE_ARRAY(gentity_handles);

    // check or pre-allocate the coordinate arrays
  CHECK_SIZE(*normals, double, 3*gentity_handles_size);
  
  CubitVector dumvec;
  const RefFace *this_face;

  const double *coord_vec = coordinates;
  double *normal_vec = *normals;

  for (int i = 0; i < gentity_handles_size; i++) {
    this_face = dynamic_cast<const RefFace*>(handle_array[i]);
    if (NULL == this_face) {
      TSTTG_processError(TSTTB_INVALID_ENTITY_TYPE, "Entities passed into gentityNormal must be faces.");
      continue;
    }
    else {
      dumvec = const_cast<RefFace*>(this_face)->normal_at(CubitVector(&coord_vec[3*i]));
      dumvec.get_xyz(&normal_vec[3*i]);
    }
  }

  RETURN(TSTTB_SUCCESS);
}

/**
 * Return the tangent at point on specified entities.  Returns error
 * if any input entity is not a gedge.  Input coordinates and tangents
 * are interleaved in the arrays.
 * @param gentity_handles The gentities being queried
 * @param coordinates Input coordinates, interleaved
 * @param tangents The tangents at the specified points, interleaved
 */
enum TSTTB_ErrorType
TSTTG_gentityTangent (TSTTG_Instance instance,
                      /*in*/ ARRAY_IN_DECL(TSTTG_CGentityHandle, gentity_handles),
                      /*in*/ CONST_ARRAY_IN_DECL(double, coordinates),
                      /*out*/ ARRAY_INOUT_DECL(double, tangents))
{
  const RefEntity **handle_array = CONST_ENTITY_HANDLE_ARRAY(gentity_handles);

    // check or pre-allocate the coordinate arrays
  CHECK_SIZE(*tangents, double, 3*gentity_handles_size);
  
  CubitVector dumvec;
  const RefEdge *this_edge;

  const double *coord_vec = coordinates;
  double *tangent_vec = *tangents;

  for (int i = 0; i < gentity_handles_size; i++) {
    this_edge = dynamic_cast<const RefEdge*>(handle_array[i]);
    if (NULL == this_edge) {
      TSTTG_processError(TSTTB_INVALID_ENTITY_TYPE, "Entities passed into gentityTangent must be edges.");
      continue;
    }
    else {
      const_cast<RefEdge*>(this_edge)->tangent(CubitVector(&coord_vec[3*i]), dumvec);
      dumvec.get_xyz(&tangent_vec[3*i]);
    }
  }

  RETURN(TSTTB_SUCCESS);
}

/**
 * Return the bounding boxex of given entities; coordinates returned
 * interleaved.
 * @param gentity_handles The gentities being queried
 * @param min_corners Minimum corner coordinates of the boxes, interleaved
 * @param max_corners Maximum corner coordinates of the boxes, interleaved
 */
enum TSTTB_ErrorType
TSTTG_gentityBoundingBox (TSTTG_Instance instance,
                          /*in*/ ARRAY_IN_DECL(TSTTG_GentityHandle, gentity_handles),
                          /*out*/ ARRAY_INOUT_DECL(double, min_corner),
                          /*out*/ ARRAY_INOUT_DECL(double, max_corner))
{
  CubitBox temp_box;
  RefEntity **handle_array = ENTITY_HANDLE_ARRAY(gentity_handles);

    // check or pre-allocate the coordinate arrays
  CHECK_SIZE(*min_corner, double, 3*gentity_handles_size);
  CHECK_SIZE(*max_corner, double, 3*gentity_handles_size);
  
  double *min_vec = *min_corner;
  double *max_vec = *max_corner;

  for (int i = 0; i < gentity_handles_size; i++) {
    BasicTopologyEntity *this_bte = dynamic_cast<BasicTopologyEntity*>(handle_array[i]);
    if (NULL != this_bte)
      temp_box = this_bte->bounding_box();
    else {
      Body *this_body = dynamic_cast<Body*>(handle_array[i]);
      if (NULL != this_body)
        temp_box = this_body->bounding_box();
      else {
        TSTTG_processError(TSTTB_INVALID_ENTITY_TYPE, "Entities passed into gentityBoundingBox must be vertex, edge, face, or region.");
        continue;
      }
    }

    temp_box.minimum().get_xyz(&min_vec[3*i]);
    temp_box.maximum().get_xyz(&max_vec[3*i]);
  }

  RETURN(TSTTB_SUCCESS);
}

/**
 * Return the coordinates of the specified vertices; returns error if any
 * of the entities are not gvertices.  Coordinates returned interleaved.
 * @param gentity_handles The gentities being queried
 * @param coordinates The coordinates of the gvertices, interleaved.
 */
enum TSTTB_ErrorType
TSTTG_getGvertexCoordinates (TSTTG_Instance instance,
                             /*in*/ ARRAY_IN_DECL(TSTTG_CGentityHandle, gentity_handles),
                             /*out*/ ARRAY_INOUT_DECL(double, coordinates))
{
  const RefEntity **handle_array = CONST_ENTITY_HANDLE_ARRAY(gentity_handles);

    // check or pre-allocate the coordinate arrays
  CHECK_SIZE(*coordinates, double, 3*gentity_handles_size);
  
  CubitVector dumvec;
  const RefVertex *this_vertex;

  double *coord_vec = *coordinates;

  for (int i = 0; i < gentity_handles_size; i++) {
    this_vertex = dynamic_cast<const RefVertex*>(handle_array[i]);
    if (NULL == this_vertex) {
      TSTTG_processError(TSTTB_INVALID_ENTITY_TYPE, "Entities passed into getGvertexCoordinates must be vertices.");
      continue;
    }
    else {
      dumvec = this_vertex->coordinates();
      dumvec.get_xyz(&coord_vec[3*i]);
    }
  }

  RETURN(TSTTB_SUCCESS);
}

/**
 * Return the sense of a gface with respect to a gregion.  Sense is either
 * forward (=1), reverse (=-1), both (=2), or unknown (=0).  Error is returned
 * if first entity is not a gface or second entity is not a gregion.
 * @param gface Gface whose sense is being queried.
 * @param gregion Gregion gface is being queried with respect to
 */
int
TSTTG_getGnormalSense (TSTTG_Instance instance,
                       /*in*/ TSTTG_CGentityHandle gface,
                       /*in*/ TSTTG_CGentityHandle gregion)
{
  const RefFace *face_ent = dynamic_cast<const RefFace*>(CONST_ENTITY_HANDLE(gface));
  if (NULL == face_ent) {
    TSTTG_processError(TSTTB_INVALID_ENTITY_TYPE, "1st argument to getGnormalSense must be a face.");
    return 0;
  }
  const RefVolume *volume_ent = dynamic_cast<const RefVolume*>(CONST_ENTITY_HANDLE(gregion));
  if (NULL == face_ent) {
    TSTTG_processError(TSTTB_INVALID_ENTITY_TYPE, "2nd argument to getGnormalSense must be a region.");
    return 0;
  }
  CubitSense this_sense = const_cast<RefFace*>(face_ent)->sense(const_cast<RefVolume*>(volume_ent));
  int tstt_sense = 0;
  if (CUBIT_FORWARD == this_sense) tstt_sense = 1;
  else if (CUBIT_REVERSED == this_sense) tstt_sense = -1;
  
  return tstt_sense;
}

/**
 * Return the sense of a gedge with respect to a gface.  Sense is either
 * forward (=1), reverse (=-1), both (=2), or unknown (=0).  Error is returned
 * if first entity is not a gedge or second entity is not a gface.
 * @param gedge Gedge whose sense is being queried.
 * @param gface Gface gedge is being queried with respect to
 */
int
TSTTG_getGtangentSense (TSTTG_Instance instance,
                        /*in*/ TSTTG_CGentityHandle gedge,
                        /*in*/ TSTTG_CGentityHandle gface)
{
  const RefEdge *edge_ent = dynamic_cast<const RefEdge*>(CONST_ENTITY_HANDLE(gedge));
  if (NULL == edge_ent) {
    TSTTG_processError(TSTTB_INVALID_ENTITY_TYPE, "1st argument to getGtangentSense must be an edge.");
    return 0;
  }
  const RefFace *face_ent = dynamic_cast<const RefFace*>(CONST_ENTITY_HANDLE(gface));
  if (NULL == face_ent) {
    TSTTG_processError(TSTTB_INVALID_ENTITY_TYPE, "2nd argument to getGtangentSense must be a face.");
    return 0;
  }
  CubitSense this_sense = const_cast<RefEdge*>(edge_ent)->sense(const_cast<RefFace*>(face_ent));
  int tstt_sense = 0;
  if (CUBIT_FORWARD == this_sense) tstt_sense = 1;
  else if (CUBIT_REVERSED == this_sense) tstt_sense = -1;
  
  return tstt_sense;
}

/**
 * Return the sense of a gedge with respect to a specified order of
 * vertices bounding the gedge.  Sense is either forward (=1), reverse (=-1), 
 * or unknown (=0).  Error is returned if any gentities are not the expected
 * type or if the gedge is bounded by only one gvertex (in this case, use
 * getGtangentSense).
 * @param gedge Gedge whose sense is being queried.
 * @param gvertex1 First gvertex
 * @param gvertex2 Second gvertex
 */
int
TSTTG_getGvertexTangentSense (TSTTG_Instance instance,
                              /*in*/ TSTTG_CGentityHandle gedge,
                              /*in*/ TSTTG_CGentityHandle gvertex1,
                              /*in*/ TSTTG_CGentityHandle gvertex2)
{
  const RefEdge *this_edge = dynamic_cast<const RefEdge*>(CONST_ENTITY_HANDLE(gedge));
  const RefVertex *vertex1 = dynamic_cast<const RefVertex*>(CONST_ENTITY_HANDLE(gvertex1));
  const RefVertex *vertex2 = dynamic_cast<const RefVertex*>(CONST_ENTITY_HANDLE(gvertex2));
  if (NULL == this_edge || NULL == vertex1 || NULL == vertex2) {
    TSTTG_processError(TSTTB_INVALID_ENTITY_TYPE, "Bad entity argument to getGvertexTangentSense.");
    return 0;
  }
  else if (const_cast<RefEdge*>(this_edge)->start_vertex() == vertex1 &&
           const_cast<RefEdge*>(this_edge)->end_vertex() == vertex2)
      // forward sense
    return 1;
  else if (const_cast<RefEdge*>(this_edge)->start_vertex() == vertex2 &&
           const_cast<RefEdge*>(this_edge)->end_vertex() == vertex1)
    return -1;
  else
    return 0;
}

/**
 * Returns the number of gentity_sets contained in a given model or
 * gentity_set one level deep
 * @param gentity_set_handle Gentity set being queried
 * @return Number of entity sets in gentity_set_handle
 */
int
TSTTG_getNumEntSets (TSTTG_Instance instance,
                     /*in*/ TSTTG_CGentitySetHandle entity_set,
                     /*in*/ int num_hops)
{
    // HJK: num_hops has to be handled
  if (1 < num_hops) {
    TSTTG_processError(TSTTB_NOT_SUPPORTED, "Num_hops argument not yet supported.");
    return -1;
  }
  
  const RefGroup *this_set = CONST_SET_HANDLE(entity_set);
  if (NULL == this_set)
    return RefEntityFactory::instance()->num_ref_groups();
  else {
    DLIList<RefEntity*> tmp_ents;
    DLIList<RefGroup*> groups;
    const_cast<RefGroup*>(this_set)->get_child_ref_entities(tmp_ents);
    CAST_LIST(tmp_ents, groups, RefGroup);
    return groups.size();
  }

  return -1;
}

/**
 * Returns the gentity_set handles contained in a given model or
 * gentity_set one level deep
 * @param gentity_set_handle Gentity set being queried
 * @param contained_gentity_set_handles Number of entity sets in 
 *  gentity_set_handle
 */
enum TSTTB_ErrorType
TSTTG_getEntSets (TSTTG_Instance instance,
                  /*in*/ TSTTG_CGentitySetHandle entity_set,
                  /*in*/ const int num_hops,
                  /*inout*/ ARRAY_INOUT_DECL(TSTTG_GentitySetHandle, 
                                             contained_entity_set_handles))
{
  if (1 < num_hops) {
    TSTTG_processError(TSTTB_NOT_SUPPORTED, "Num_hops argument not yet supported.");
    RETURN(TSTTB_NOT_SUPPORTED);
  }

  const RefGroup *this_set = CONST_SET_HANDLE(entity_set);
  DLIList<RefEntity*> tmp_ents;
  DLIList<RefGroup*> groups;
  if (NULL == this_set)
    RefEntityFactory::instance()->ref_groups(groups);
  else {
    const_cast<RefGroup*>(this_set)->get_child_ref_entities(tmp_ents);
    CAST_LIST(tmp_ents, groups, RefGroup);
  }
  
  CHECK_SIZE(*contained_entity_set_handles, TSTTG_GentitySetHandle, groups.size());

  groups.copy_to(*SET_HANDLE_ARRAY_PTR(contained_entity_set_handles));
  
  RETURN(TSTTB_SUCCESS);
}

/**
 * This function allows a new gentity_set to be created.  The user may 
 * set the multiset, ordered, isMesh, flags as needed; otherwise default values 
 * (all false) will be used.  On creation, Gentitysets are empty of
 * entities and contained in the parent geometry interface.  They must be
 * explicitly filled with entities using the addGentities call and
 * relationships with other Gentitysets must be done through the
 * addGentityset and parent/child relationship calls.
 * @param multiset If true, gentities can appear more than once in this gentity_set
 * @param ordered If true, order of addition and removal is maintained for
 *   this gentity_set
 * @param gentity_set_created Gentity_set created by this function
 */
enum TSTTB_ErrorType
TSTTG_createEntSet (TSTTG_Instance instance,
                    /*in*/ bool isList,
                    /*out*/ TSTTG_GentitySetHandle *entity_set)
{
  *entity_set = RefEntityFactory::instance()->construct_RefGroup();
    // need to set a tag denoting multiset or not...
  RETURN(TSTTB_SUCCESS);
}

/**
 *   Destroy the gentity set.  This method only destroys the grouping of
 *   gentities, not the gentities themselves.
 * @param gentity_set Gentity_set to be destroyed
 */
enum TSTTB_ErrorType
TSTTG_destroyEntSet (TSTTG_Instance instance,
                     /*in*/ TSTTG_GentitySetHandle entity_set)
{
  if (NULL == entity_set) {
    TSTTG_processError(TSTTB_INVALID_ARGUMENT, "Can't destroy interface set.");
  }
      
  RefGroup::delete_group(SET_HANDLE(entity_set));
  RETURN(TSTTB_SUCCESS);
}

bool
TSTTG_isList (TSTTG_Instance instance,
                /*in*/ TSTTG_CGentitySetHandle entity_set) 
{
  return true;
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
enum TSTTB_ErrorType
TSTTG_addEntSet (TSTTG_Instance instance,
                 /*in*/ TSTTG_CGentitySetHandle entity_set_to_add,
                 /*inout*/ TSTTG_GentitySetHandle *entity_set_handle)
{
  const RefEntity *tmp_ent = CONST_SET_HANDLE(entity_set_to_add);
  return TSTTG_addEntToSet(instance, tmp_ent, entity_set_handle);
}

/**
 *   Allows the user to remove one or more gentity_sets from another
 *   gentity_set.  Users cannot delete a contained in relationship of an
 *   gentity set with the parent mesh interface so passing in a NULL
 *   value for the first argument results in no action.
 * @param gentity_set Set from which other sets are being removed
 * @param gentity_set_handles Sets being removed from gentity_set
 */
enum TSTTB_ErrorType
TSTTG_rmvEntSet (TSTTG_Instance instance,
                 /*in*/ TSTTG_CGentitySetHandle entity_set_to_remove,
                 /*inout*/ TSTTG_GentitySetHandle *entity_set_handle)
{
  const RefEntity *tmp_ent = CONST_SET_HANDLE(entity_set_to_remove);
  return TSTTG_rmvEntFromSet(instance, tmp_ent, entity_set_handle);
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
enum TSTTB_ErrorType
TSTTG_addEntToSet (TSTTG_Instance instance,
                 /*in*/ TSTTG_CGentityHandle entity_to_add,
                 /*inout*/ TSTTG_GentitySetHandle *entity_set_handle)
{
  if (NULL == entity_to_add) RETURN(TSTTB_INVALID_ARGUMENT);
  
  CubitStatus status = SET_HANDLE(*entity_set_handle)->
    add_ref_entity(const_cast<RefEntity*>(CONST_ENTITY_HANDLE(entity_to_add)));
  
  if (CUBIT_SUCCESS != status)
    TSTTG_processError(TSTTB_FAILURE, "Problem adding entity to another set.");

  RETURN(TSTTB_SUCCESS);
}

/**
 *   Allows the user to remove one or more gentity_sets from another
 *   gentity_set.  Users cannot delete a contained in relationship of an
 *   gentity set with the parent mesh interface so passing in a NULL
 *   value for the first argument results in no action.
 * @param gentity_set Set from which other sets are being removed
 * @param gentity_set_handles Sets being removed from gentity_set
 */
enum TSTTB_ErrorType
TSTTG_rmvEntFromSet (TSTTG_Instance instance,
                     /*in*/ TSTTG_CGentityHandle entity_to_remove,
                     /*inout*/ TSTTG_GentitySetHandle *entity_set_handle)
{
  if (NULL == *entity_set_handle) RETURN(TSTTB_INVALID_ARGUMENT);
  
  CubitStatus status = SET_HANDLE(*entity_set_handle)->
    remove_ref_entity(const_cast<RefEntity*>(CONST_ENTITY_HANDLE(entity_to_remove)));
  if (CUBIT_SUCCESS != status) {
    TSTTG_processError(TSTTB_FAILURE, "Problem removing entity from a set.");
    RETURN(TSTTB_FAILURE);
  }
  
  RETURN(TSTTB_SUCCESS);
}

bool
TSTTG_isEntContained (TSTTG_Instance instance,
                      /*in*/ TSTTG_CGentitySetHandle containing_entity_set,
                      /*in*/ TSTTG_CGentityHandle contained_entity) 
{
  if (NULL == containing_entity_set && NULL != contained_entity) return true;

  DLIList<RefGroup*> containing_groups;
  RefEntity *contained_entity_temp = const_cast<RefEntity*>(CONST_ENTITY_HANDLE(contained_entity));
  RefGroup::get_groups_within(contained_entity_temp, containing_groups, CUBIT_FALSE);
  if (containing_groups.size() > 0 &&
      containing_groups.move_to(const_cast<RefGroup*>(CONST_SET_HANDLE(containing_entity_set)))) return true;
  
  return false;
  
  RETURN(TSTTB_SUCCESS);
}

bool
TSTTG_isEntSetContained (TSTTG_Instance instance,
                         /*in*/ TSTTG_CGentitySetHandle containing_entity_set,
                         /*in*/ TSTTG_CGentitySetHandle contained_entity_set) 
{
  if (NULL == containing_entity_set && NULL != contained_entity_set) return true;

  DLIList<RefGroup*> containing_groups;
  RefGroup *contained_set = const_cast<RefGroup*>(CONST_SET_HANDLE(contained_entity_set));
  RefGroup::get_groups_within(const_cast<RefGroup*>(contained_set), containing_groups, 
                              CUBIT_FALSE);
  if (containing_groups.size() > 0 &&
      containing_groups.move_to(const_cast<RefGroup*>(CONST_SET_HANDLE(containing_entity_set)))) return true;
  
  return false;
  
  RETURN(TSTTB_SUCCESS);
}

/**
 *   Add existing gentities to the gentity_set (do not create them).
 *   Note that if a gentity of dimension d>0 is added to the gentityset,
 *   the lower-dimensional gentities that bound it are not
 *   automatically associated with the gentityset.
 * @param gentity_set Set being added to
 * @param gentity_handles Gentities being added to gentity_set
 */
enum TSTTB_ErrorType
TSTTG_addEntArrToSet (TSTTG_Instance instance,
                      /*in*/ ARRAY_IN_DECL(TSTTG_GentityHandle, entity_handles),
                      /*inout*/ TSTTG_GentitySetHandle *entity_set)
{
  if (NULL == entity_set) RETURN(TSTTB_INVALID_ARGUMENT);

  RefGroup *this_set = SET_HANDLE(*entity_set);
  RefEntity **ent_array = ENTITY_HANDLE_ARRAY(entity_handles);
  CubitStatus status = CUBIT_SUCCESS, tmp_status;
  for (int i = 0; i < entity_handles_size; i++) {
    tmp_status = this_set->add_ref_entity(ent_array[i]);
    if (CUBIT_SUCCESS != tmp_status) status = tmp_status;
  }

  if (CUBIT_SUCCESS != status) {
    TSTTG_processError(TSTTB_FAILURE, "Problem adding entities to a set.");
    RETURN(TSTTB_FAILURE);
  }
  
  RETURN(TSTTB_SUCCESS);
}

/**
 *   Remove existing gentities from the gentity_set (do not delete them)
 * @param gentity_set Set being removed from
 * @param gentity_handles Gentities being removed from gentity_set
 */
enum TSTTB_ErrorType
TSTTG_rmvEntArrFromSet (TSTTG_Instance instance,
                        /*in*/ ARRAY_IN_DECL(TSTTG_GentityHandle, entity_handles),
                        /*inout*/ TSTTG_GentitySetHandle *entity_set)
{
  if (NULL == entity_set) RETURN(TSTTB_INVALID_ARGUMENT);

  RefGroup *this_set = SET_HANDLE(*entity_set);
  RefEntity **ent_array = ENTITY_HANDLE_ARRAY(entity_handles);
  CubitStatus status = CUBIT_SUCCESS, tmp_status;
  for (int i = 0; i < entity_handles_size; i++) {
    tmp_status = this_set->remove_ref_entity(ent_array[i]);
    if (CUBIT_SUCCESS != tmp_status) status = tmp_status;
  }

  if (CUBIT_SUCCESS != status) {
    TSTTG_processError(TSTTB_FAILURE, "Problem removing entities from a set.");
    RETURN(TSTTB_FAILURE);
  }
  
  RETURN(TSTTB_SUCCESS);
}

/**
 * Return the relative and absolute tolerances at the modeler level.  If
 * model does not have a modeler-wide tolerance, zero is returned for both
 * values.
 * @param relative_tolerance Relative tolerance for model as a whole
 * @param absolute_tolerance Absolute tolerance for model as a whole
 */
enum TSTTB_ErrorType
TSTTG_getGtolerance (TSTTG_Instance instance,
                     /*out*/ double *relative_tolerance,
                     /*out*/ double *absolute_tolerance)
{
  *absolute_tolerance = gqt->get_sme_resabs_tolerance();
  *relative_tolerance = 0.0;
  RETURN(TSTTB_SUCCESS);
}

/**
 * Return the relative and absolute tolerances for specified gentities.  If
 * a gentity does not have a specific tolerance, zero is returned for both
 * values.
 * @param gentity_handles Gentities being queried
 * @param relative_tolerances Relative tolerances
 * @param absolute_tolerances Absolute tolerances
 */
enum TSTTB_ErrorType
TSTTG_getGentityTolerance (TSTTG_Instance instance,
                           /*in*/ ARRAY_IN_DECL(TSTTG_CGentityHandle, gentity_handles),
                           /*out*/ ARRAY_INOUT_DECL(double, relative_tolerances),
                           /*out*/ ARRAY_INOUT_DECL(double, absolute_tolerances))
{
  CHECK_SIZE(*relative_tolerances, double, gentity_handles_size);
  CHECK_SIZE(*absolute_tolerances, double, gentity_handles_size);
  double dum_abs_tol = gqt->get_sme_resabs_tolerance();
  for (int i = 0; i < gentity_handles_size; i++) {
    (*relative_tolerances[i]) = 0.0;
    (*absolute_tolerances[i]) = dum_abs_tol;
  }
  RETURN(TSTTB_SUCCESS);
}

/**
 * Return whether a given gentity is parametric or not.  If a gentity
 * is not parametric, all of the following functions will return an error
 * when called on that entity.
 * @param gentity_handle Gentity being queried.
 */
int
TSTTG_gentityIsParametric (TSTTG_Instance instance,
                           /*in*/ TSTTG_CGentityHandle gentity_handle)
{
  const RefFace *this_face = dynamic_cast<const RefFace*>(CONST_ENTITY_HANDLE(gentity_handle));
  if (NULL != this_face)
    return const_cast<RefFace*>(this_face)->is_parametric();

  const RefEdge *this_edge = dynamic_cast<const RefEdge*>(CONST_ENTITY_HANDLE(gentity_handle));
  if (NULL != this_edge)
    return true;

    // must not be either; return false
  return false;
}

/**
 * Given sets of parametric coordinates, return the corresponding real
 * space coordinates on the gentities.  Input and output coordinates are
 * interleaved.
 * @param gentity_handles Gentities being queried.
 * @param uv Input parametric coordinates
 * @param xyz Output real space coordinates
 */
enum TSTTB_ErrorType
TSTTG_gentityUvToXyz (TSTTG_Instance instance,
                      /*in*/ ARRAY_IN_DECL(TSTTG_CGentityHandle, gentity_handles),
                      /*in*/ CONST_ARRAY_IN_DECL(double, uv),
                      /*out*/ ARRAY_INOUT_DECL(double, coordinates))
{
  const RefEntity **handle_array = CONST_ENTITY_HANDLE_ARRAY(gentity_handles);

    // check or pre-allocate the coordinate arrays
  CHECK_SIZE(*coordinates, double, 3*gentity_handles_size);
  
  CubitVector dumvec;
  CubitStatus status = CUBIT_SUCCESS, temp_status;

  const double *uv_vec = uv;
  double *coord_vec = *coordinates;

  for (int i = 0; i < gentity_handles_size; i++) {
    const RefEdge *this_edge = dynamic_cast<const RefEdge*>(handle_array[i]);
    if (NULL != this_edge) {
      temp_status = const_cast<RefEdge*>(this_edge)->position_from_u(uv_vec[2*i], dumvec); 
      if (CUBIT_SUCCESS == temp_status)
        dumvec.get_xyz(&coord_vec[3*i]);
      else status = temp_status;
      continue;
    }
    const RefFace *this_face = dynamic_cast<const RefFace*>(handle_array[i]);
    if (NULL != this_face) {
      dumvec = const_cast<RefFace*>(this_face)->position_from_u_v(uv_vec[2*i], uv_vec[2*i+1]);
      dumvec.get_xyz(&coord_vec[3*i]);
      continue;
    }
      // if we got here, zero out the points
    coord_vec[3*i] = coord_vec[3*i+1] = coord_vec[3*i+2] = 0.0;
    status = CUBIT_FAILURE;
  }

  if (CUBIT_FAILURE == status) {
    TSTTG_processError(TSTTB_FAILURE, "Problem getting xyz from uv.");
    RETURN(TSTTB_FAILURE);
  }
  
  RETURN(TSTTB_SUCCESS);
}

/**
 * Given sets of real space coordinates, return the corresponding 
 * parametric coordinates on the gentities.  Input and output coordinates 
 * are interleaved.
 * @param gentity_handles Gentities being queried.
 * @param xyz Input real space coordinates
 * @param uv Output parametric coordinates
 */
enum TSTTB_ErrorType
TSTTG_gentityXyzToUv (TSTTG_Instance instance,
                      /*in*/ ARRAY_IN_DECL(TSTTG_CGentityHandle, gentity_handles),
                      /*in*/ CONST_ARRAY_IN_DECL(double, coordinates),
                      /*out*/ ARRAY_INOUT_DECL(double, uv))
{
  const RefEntity **handle_array = CONST_ENTITY_HANDLE_ARRAY(gentity_handles);

    // check or pre-allocate the coordinate arrays
  CHECK_SIZE(*uv, double, 2*gentity_handles_size);
  
  CubitStatus status = CUBIT_SUCCESS, temp_status;

  double *uv_vec = *uv;
  const double *coord_vec = coordinates;

  for (int i = 0; i < gentity_handles_size; i++) {
    const RefEdge *this_edge = dynamic_cast<const RefEdge*>(handle_array[i]);
    if (NULL != this_edge) {
      uv_vec[2*i] = const_cast<RefEdge*>(this_edge)->u_from_position(CubitVector(&coord_vec[3*i])); 
      continue;
    }
    const RefFace *this_face = dynamic_cast<const RefFace*>(handle_array[i]);
    if (NULL != this_face) {
      temp_status = const_cast<RefFace*>(this_face)->u_v_from_position(CubitVector(&coord_vec[3*i]),
                                                 uv_vec[2*i], uv_vec[2*i+1]);
      continue;
    }
      // if we got here, zero out the points
    uv_vec[2*i] = uv_vec[2*i+1] = 0.0;
    status = CUBIT_FAILURE;
  }

  if (CUBIT_FAILURE == status) {
    TSTTG_processError(TSTTB_FAILURE, "Problem getting uv from xyz.");
    RETURN(TSTTB_FAILURE);
  }
  
  RETURN(TSTTB_SUCCESS);
}

/**
 * Return the uv range of the specified gentities.  Parameters are interleaved.
 * @param gentity_handles Gentities being queried.
 * @param uv_min Minimum parameters of gentities, interleaved
 * @param uv_max Maximum parameters of gentities, interleaved
 */
enum TSTTB_ErrorType
TSTTG_gentityUvRange (TSTTG_Instance instance,
                      /*in*/ ARRAY_IN_DECL(TSTTG_CGentityHandle, gentity_handles),
                      /*out*/ ARRAY_INOUT_DECL(double, uv_min),
                      /*out*/ ARRAY_INOUT_DECL(double, uv_max))
{
  const RefEntity **handle_array = CONST_ENTITY_HANDLE_ARRAY(gentity_handles);

    // check or pre-allocate the coordinate arrays
  CHECK_SIZE(*uv_min, double, 2*gentity_handles_size);
  CHECK_SIZE(*uv_max, double, 2*gentity_handles_size);
  
  CubitStatus status = CUBIT_SUCCESS;

  double *uvmin_vec = *uv_min;
  double *uvmax_vec = *uv_max;

  for (int i = 0; i < gentity_handles_size; i++) {
    const RefEdge *this_edge = dynamic_cast<const RefEdge*>(handle_array[i]);
    if (NULL != this_edge) {
      CubitBoolean result = const_cast<RefEdge*>(this_edge)->get_param_range(uvmin_vec[2*i], uvmax_vec[2*i]);
      if (!result) status = CUBIT_FAILURE;
      continue;
    }
    const RefFace *this_face = dynamic_cast<const RefFace*>(handle_array[i]);
    if (NULL != this_face) {
      CubitBoolean result = const_cast<RefFace*>(this_face)->get_param_range_U(uvmin_vec[2*i], uvmax_vec[2*i]);
      if (!result) status = CUBIT_FAILURE;
      result = const_cast<RefFace*>(this_face)->get_param_range_V(uvmin_vec[2*i+1], uvmax_vec[2*i+1]);
      if (!result) status = CUBIT_FAILURE;
      continue;
    }
      // if we got here, zero out the points
    uvmin_vec[2*i] = uvmin_vec[2*i+1] = uvmax_vec[2*i] = uvmax_vec[2*i+1] = 0.0;
    status = CUBIT_FAILURE;
  }

  if (CUBIT_FAILURE == status) {
    TSTTG_processError(TSTTB_FAILURE, "Problem getting uv range.");
    RETURN(TSTTB_FAILURE);
  }
  
  RETURN(TSTTB_SUCCESS);
}

/**
 * Given source gentities, parametric positions on those gentities, and 
 * bounding gentities, return the parametric positions on the bounding
 * gentities.  If a source gentity is a gvertex, parametric positions for
 * that entry are ignored.  In cases where multiple source entities are
 * input, two input and output parameters per gentity is assumed, even if
 * input consists only of gedges or gvertices.
 * @param src_gentity_handles Source gentities
 * @param src_uv Uv positions on source gentities
 * @param trg_gentity_handles Target gentities
 * @param trg_uv Uv positions on target gentities
 */
enum TSTTB_ErrorType
TSTTG_Greparam_edge_face (TSTTG_Instance instance,
                          /*in*/ ARRAY_IN_DECL(TSTTG_GentityHandle, src_gentity_handles),
                          /*in*/ ARRAY_IN_DECL(double, src_uv),
                          /*in*/ ARRAY_IN_DECL(TSTTG_GentityHandle, trg_gentity_handles),
                          /*in*/ ARRAY_IN_DECL(double, trg_uv))
{
  TSTTG_processError(TSTTB_NOT_SUPPORTED, "Not implemented.");
  RETURN(TSTTB_NOT_SUPPORTED);
}

/**
 * Return the normals at specified uv positions on gfaces.  If any
 * gentity input is not a face, returns error.  Input parameters and 
 * output normals are interleaved.
 * @param gface_handles The entities being queried
 * @param parameters The uv parameters of points being queried, interleaved
 * @param normals Normals at specified points, interleaved
 */
enum TSTTB_ErrorType
TSTTG_gentityNormalUv (TSTTG_Instance instance,
                       /*in*/ ARRAY_IN_DECL(TSTTG_CGentityHandle, gface_handles),
                       /*in*/ CONST_ARRAY_IN_DECL(double, parameters),
                       /*out*/ ARRAY_INOUT_DECL(double, normals))
{
  enum TSTTB_ErrorType result = 
    TSTTG_gentityUvToXyz(instance, ARRAY_IN(gface_handles), ARRAY_IN(parameters), 
                         ARRAY_INOUT(normals));
  if (TSTTB_SUCCESS == result)
    result = TSTTG_gentityNormal(instance, ARRAY_IN(gface_handles), *normals, *normals_size, 
                           ARRAY_INOUT(normals));
  RETURN(result);
}

/**
 * Return the tangents at specified u positions on gedges.  If any
 * gentity input is not a face, returns error.  Output normals are 
 * interleaved.
 * @param gentity_handles The gedges being queried
 * @param parameters The u parameters of points being queried
 * @param tangents Tangents at specified points, interleaved
 */
enum TSTTB_ErrorType
TSTTG_gentityTangentU (TSTTG_Instance instance,
                       /*in*/ ARRAY_IN_DECL(TSTTG_CGentityHandle, gedge_handles),
                       /*in*/ CONST_ARRAY_IN_DECL(double, parameters),
                       /*out*/ ARRAY_INOUT_DECL(double, tangents))
{
  enum TSTTB_ErrorType result =
    TSTTG_gentityUvToXyz(instance, ARRAY_IN(gedge_handles), ARRAY_IN(parameters), 
                         ARRAY_INOUT(tangents));
  if (TSTTB_SUCCESS == result)
    result = TSTTG_gentityNormal(instance, ARRAY_IN(gedge_handles), *tangents, *tangents_size,
                                 ARRAY_INOUT(tangents));
  RETURN(result);
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
enum TSTTB_ErrorType
TSTTG_subtract (TSTTG_Instance instance,
                /*in*/ TSTTG_CGentitySetHandle entity_set_1,
                /*in*/ TSTTG_CGentitySetHandle entity_set_2,
                /*out*/ TSTTG_GentitySetHandle *result_entity_set)
{
  const RefGroup *set1 = CONST_SET_HANDLE(entity_set_1);
  const RefGroup *set2 = CONST_SET_HANDLE(entity_set_2);
  RefGroup *set3 = SET_HANDLE(*result_entity_set);
  if (NULL == set3) {
    set3 = RefEntityFactory::instance()->construct_RefGroup();
    *result_entity_set = set3;
  }
    
  set3->subtract(const_cast<RefGroup*>(set2), const_cast<RefGroup*>(set1));
  RETURN(TSTTB_SUCCESS);
}

/**
 *   Intersect gentity_set_1 and gentity_set_2; the result is returned in
 *   gentity_set_1.  If both are multisets and both contain 2 entries of
 *   gentity k, the intersection contains 2 entries as well.
 */
enum TSTTB_ErrorType
TSTTG_intersect (TSTTG_Instance instance,
                 /*in*/ TSTTG_CGentitySetHandle entity_set_1,
                 /*in*/ TSTTG_CGentitySetHandle entity_set_2,
                 /*out*/ TSTTG_GentitySetHandle *result_entity_set)
{
  const RefGroup *set1 = CONST_SET_HANDLE(entity_set_1);
  const RefGroup *set2 = CONST_SET_HANDLE(entity_set_2);
  RefGroup *set3 = SET_HANDLE(*result_entity_set);
  if (NULL == set3) {
    set3 = RefEntityFactory::instance()->construct_RefGroup();
    *result_entity_set = set3;
  }

  set3->intersect(const_cast<RefGroup*>(set2), const_cast<RefGroup*>(set1));
  RETURN(TSTTB_SUCCESS);
}

/**
 *   Union the gentities in gentity_set_1 and gentity_set_2; the result is
 *   returned in gentity_set_1.  The multiset flag in gentityset 1
 *   determines if duplicate entries are allowed.
 */
enum TSTTB_ErrorType
TSTTG_unite (TSTTG_Instance instance,
             /*in*/ TSTTG_CGentitySetHandle entity_set_1,
             /*in*/ TSTTG_CGentitySetHandle entity_set_2,
             /*out*/ TSTTG_GentitySetHandle *result_entity_set)
{
  const RefGroup *set1 = CONST_SET_HANDLE(entity_set_1);
  const RefGroup *set2 = CONST_SET_HANDLE(entity_set_2);
  RefGroup *set3 = SET_HANDLE(*result_entity_set);
  if (NULL == set3) {
    set3 = RefEntityFactory::instance()->construct_RefGroup();
    *result_entity_set = set3;
  }

  set3->unite(const_cast<RefGroup*>(set2), const_cast<RefGroup*>(set1));
  RETURN(TSTTB_SUCCESS);
}

enum TSTTB_ErrorType
TSTTG_Copy (TSTTG_Instance instance,
            /*in*/ TSTTG_GentityHandle geom_entity,
            /*out*/ TSTTG_GentityHandle *geom_entity2)
{
  Body *this_body = dynamic_cast<Body*>(ENTITY_HANDLE(geom_entity));
  RefVolume *this_vol = dynamic_cast<RefVolume*>(ENTITY_HANDLE(geom_entity));
  if (NULL != this_vol || NULL != this_body) {
      // need to get the associated body, since cgm only supports copying bodies,
      // not volumes
    if (NULL == this_body) {
      this_body = this_vol->get_body_ptr();
      if (NULL == this_body) {
        TSTTG_processError(TSTTB_FAILURE, "Can't get body from volume.");
        RETURN(TSTTB_INVALID_ARGUMENT);
      }
    }

    RefEntity *temp_entity = gmt->copy_body(this_body);
    *geom_entity2 = temp_entity;
  }
  else {
    RefEntity *this_ent = ENTITY_HANDLE(geom_entity);
    RefEntity *temp_entity = gmt->copy_refentity(this_ent);
    *geom_entity2 = temp_entity;
  }

  if (NULL == *geom_entity2) {
    TSTTG_processError(TSTTB_FAILURE, "NULL returned from CGM copy.");
    RETURN(TSTTB_FAILURE);
  }
  
  RETURN(TSTTB_SUCCESS);
}
      
enum TSTTB_ErrorType
TSTTG_SweepAboutAxis (TSTTG_Instance instance,
                      /*in*/ TSTTG_GentityHandle geom_entity,
                      /*in*/ const double angle,
                      /*in*/ const double axis_normal_x,
                      /*in*/ const double axis_normal_y,
                      /*in*/ const double axis_normal_z,
                      /*out*/ TSTTG_GentityHandle *geom_entity2)
{
  DLIList<RefEntity*> ents;
  RefEntity *this_ent = ENTITY_HANDLE(geom_entity);
  if (NULL == this_ent) {
    TSTTG_processError(TSTTB_FAILURE, "Can't get ref entity from specified entity.");
    RETURN(TSTTB_INVALID_ARGUMENT);
  }
  else if (this_ent->dimension() > 2) {
    TSTTG_processError(TSTTB_INVALID_ENTITY_TYPE, "Can't sweep 3-d entities.");
    RETURN(TSTTB_INVALID_ARGUMENT);
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
    TSTTG_processError(TSTTB_FAILURE, "Sweep returned failure.");
    RETURN(TSTTB_FAILURE);
  }
  
  else {
      // HACK to get last entity created, because cgm doesn't return
      // body created by sweep
    RefEntity *new_body = RefEntityFactory::instance()->get_last_body();
    *geom_entity2 = new_body;

      // now we know it's succeeded, delete original body
    Body *this_body = dynamic_cast<TopologyEntity*>(this_ent)->body();
    if (NULL != this_body)
      gqt->delete_RefEntity(this_body);
    else
      gqt->delete_RefEntity(this_ent);
  }


  if (NULL == *geom_entity2) {
    RETURN(TSTTB_FAILURE);
  }
  
  RETURN(TSTTB_SUCCESS);
}
      
enum TSTTB_ErrorType
TSTTG_Delete (TSTTG_Instance instance,
              /*in*/ TSTTG_GentityHandle geom_entity)
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
    TSTTG_processError(TSTTB_FAILURE, "Problems deleting entity.");
    RETURN(TSTTB_FAILURE);
  }

    // check to see if this was last thing deleted; if so, reset ids
  RefEntityFactory *rfi = RefEntityFactory::instance();
  if (rfi->num_bodies() == 0 && 
      rfi->num_ref_volumes() == 0 &&
      rfi->num_ref_faces() == 0 &&
      rfi->num_ref_edges() == 0 &&
      rfi->num_ref_vertices() == 0)
    rfi->reset_ids();
  
  RETURN(TSTTB_SUCCESS);
}

enum TSTTB_ErrorType
TSTTG_Brick (TSTTG_Instance instance,
             /*in*/ const double x,
             /*in*/ const double y,
             /*in*/ const double z,
             /*out*/ TSTTG_GentityHandle *geom_entity)
{
  double tmp_x = x;
  double tmp_y = y;
  double tmp_z = z;
  
  if (0.0 == y && 0.0 == z) {
    tmp_y = x;
    tmp_z = x;
  }
  
  if (0.0 >= tmp_x || 0.0 >= tmp_y || 0.0 >= tmp_z) {
    TSTTG_processError(TSTTB_INVALID_ARGUMENT, "Dimensions must be >= 0, or y & z must both be zero.");
    RETURN(TSTTB_INVALID_ARGUMENT);
  }
    
  RefEntity *temp_body = gmt->brick(tmp_x, tmp_y, tmp_z);
  *geom_entity = temp_body;

  if (NULL == *geom_entity) {
    RETURN(TSTTB_FAILURE);
  }
  
  RETURN(TSTTB_SUCCESS);
}
          
enum TSTTB_ErrorType
TSTTG_Cylinder (TSTTG_Instance instance,
                /*in*/ const double height,
                /*in*/ const double major_rad,
                /*in*/ const double minor_rad,
                /*out*/ TSTTG_GentityHandle *geom_entity)
{
  double tmp_minor = (0.0 == minor_rad ? major_rad : minor_rad);
  RefEntity *temp_body = 
    gmt->cylinder(height, major_rad, tmp_minor, major_rad);
  *geom_entity = temp_body;


  if (NULL == *geom_entity) {
    RETURN(TSTTB_FAILURE);
  }
  
  RETURN(TSTTB_SUCCESS);
}
          
enum TSTTB_ErrorType
TSTTG_Torus (TSTTG_Instance instance,
             /*in*/ const double major_rad,
             /*in*/ const double minor_rad,
             /*out*/ TSTTG_GentityHandle *geom_entity)
{
  if (minor_rad >= major_rad) {
    TSTTG_processError(TSTTB_INVALID_ARGUMENT, "Major radius must be greater than minor radius for tori.");
    RETURN(TSTTB_INVALID_ARGUMENT);
  }
  
  RefEntity *temp_body = gmt->torus(major_rad, minor_rad);
  *geom_entity = temp_body;

  if (NULL == *geom_entity) {
    RETURN(TSTTB_FAILURE);
  }
  
  RETURN(TSTTB_SUCCESS);
}

enum TSTTB_ErrorType
TSTTG_Move (TSTTG_Instance instance,
            /*inout*/ TSTTG_GentityHandle *geom_entity,
            /*in*/ const double x,
            /*in*/ const double y,
            /*in*/ const double z)
{
  CubitVector vec(x, y, z);
  Body *this_bod = dynamic_cast<Body*>(ENTITY_HANDLE(*geom_entity));
  CubitStatus result;
  if (NULL != this_bod) {
    result = gqt->translate(this_bod, vec);
    if (CUBIT_SUCCESS != result) {
      TSTTG_processError(TSTTB_FAILURE, "Failed to move body.");
      RETURN(TSTTB_FAILURE);
    }
    
    RETURN(TSTTB_SUCCESS);
  }
  
  BasicTopologyEntity *this_bte = dynamic_cast<BasicTopologyEntity*>(ENTITY_HANDLE(*geom_entity));
  if (NULL != this_bte) {
      // non-body move; check to see if there are any siblings to this entity in the
      // same body; if so, we can't move it; if not, get the body and move that; if
      // there is no body, it's a free entity and we can move it anyway
    Body *this_body = this_bte->body();
    if (NULL == this_body) {
      result = gqt->translate(this_bte, vec);
      if (CUBIT_SUCCESS != result) {
        TSTTG_processError(TSTTB_FAILURE, "Failed to move entity.");
        RETURN(TSTTB_FAILURE);
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
          TSTTG_processError(TSTTB_FAILURE, "Failed to move body even only one entity of that"
                             " dimension in the body.");
          RETURN(TSTTB_FAILURE);
        }
      }
      else {
        TSTTG_processError(TSTTB_FAILURE, "Too many siblings for an entity to move it.");
        RETURN(TSTTB_FAILURE);
      }
    }

    RETURN(TSTTB_SUCCESS);
  }
  
  TSTTG_processError(TSTTB_INVALID_ENTITY_TYPE, "Wrong type of entity specified for move.");
  RETURN(TSTTB_INVALID_ENTITY_TYPE);
}
      
enum TSTTB_ErrorType
TSTTG_Rotate (TSTTG_Instance instance,
              /*inout*/ TSTTG_GentityHandle *geom_entity,
              /*in*/ const double angle,
              /*in*/ const double axis_normal_x,
              /*in*/ const double axis_normal_y,
              /*in*/ const double axis_normal_z)
{
  CubitVector this_axis(axis_normal_x, axis_normal_y, axis_normal_z);
  Body *this_bod = dynamic_cast<Body*>(ENTITY_HANDLE(*geom_entity));
  CubitStatus result;
  if (NULL != this_bod) {
    result = gqt->rotate(this_bod, this_axis, angle);
    if (CUBIT_SUCCESS != result) {
      TSTTG_processError(TSTTB_FAILURE, "Failed to rotate body.");
      RETURN(TSTTB_FAILURE);
    }
    
    RETURN(TSTTB_SUCCESS);
  }
  
  BasicTopologyEntity *this_bte = dynamic_cast<BasicTopologyEntity*>(ENTITY_HANDLE(*geom_entity));
  if (NULL != this_bte) {
    result = gqt->rotate(this_bte, this_axis, angle);
    if (CUBIT_SUCCESS != result) {
      TSTTG_processError(TSTTB_FAILURE, "Failed to rotate entity.");
      RETURN(TSTTB_FAILURE);
    }
    
    RETURN(TSTTB_SUCCESS);
  }
  
  TSTTG_processError(TSTTB_INVALID_ENTITY_TYPE, "Wrong type of entity specified for move.");
  RETURN(TSTTB_INVALID_ENTITY_TYPE);
}
      
enum TSTTB_ErrorType
TSTTG_Reflect (TSTTG_Instance instance,
               /*inout*/ TSTTG_GentityHandle *geom_entity,
               /*in*/ const double plane_normal_x,
               /*in*/ const double plane_normal_y,
               /*in*/ const double plane_normal_z)
{
  CubitVector this_plane(plane_normal_x, plane_normal_y, plane_normal_z);
  Body *this_bod = dynamic_cast<Body*>(ENTITY_HANDLE(*geom_entity));
  DLIList<Body*> bods;
  bods.append(this_bod);
  CubitStatus result;
  if (NULL != this_bod) {
    result = gqt->reflect(bods, this_plane);
    if (CUBIT_SUCCESS != result) {
      TSTTG_processError(TSTTB_FAILURE, "Failed to reflect body.");
      RETURN(TSTTB_FAILURE);
    }
    
    RETURN(TSTTB_SUCCESS);
  }
  
  BasicTopologyEntity *this_bte = dynamic_cast<BasicTopologyEntity*>(ENTITY_HANDLE(*geom_entity));
  if (NULL != this_bte) {
    result = gqt->reflect(this_bte, this_plane);
    if (CUBIT_SUCCESS != result) {
      TSTTG_processError(TSTTB_FAILURE, "Failed to reflect entity.");
      RETURN(TSTTB_FAILURE);
    }
    
    RETURN(TSTTB_SUCCESS);
  }
  
  TSTTG_processError(TSTTB_INVALID_ENTITY_TYPE, "Wrong type of entity specified for reflect.");
  RETURN(TSTTB_INVALID_ENTITY_TYPE);
}

enum TSTTB_ErrorType
TSTTG_Scale (TSTTG_Instance instance,
             /*inout*/ TSTTG_GentityHandle *geom_entity,
             /*in*/ const double scale_x,
             /*in*/ const double scale_y,
             /*in*/ const double scale_z) 
{
  CubitVector factor(scale_x, scale_y, scale_z);
  Body *this_bod = dynamic_cast<Body*>(ENTITY_HANDLE(*geom_entity));
  CubitStatus result;
  if (NULL != this_bod) {
    result = gqt->scale(this_bod, factor);
    if (CUBIT_SUCCESS != result) {
      TSTTG_processError(TSTTB_FAILURE, "Failed to scale body.");
      RETURN(TSTTB_FAILURE);
    }
    
    RETURN(TSTTB_SUCCESS);
  }
  
  BasicTopologyEntity *this_bte = dynamic_cast<BasicTopologyEntity*>(ENTITY_HANDLE(*geom_entity));
    // non-body move; check to see if there are any siblings to this entity in the
    // same body; if so, we can't move it; if not, get the body and move that; if
    // there is no body, it's a free entity and we can move it anyway
  Body *this_body = this_bte->body();
  if (NULL == this_body && NULL != this_bte) {
    result = gqt->scale(this_bte, factor);
    if (CUBIT_SUCCESS != result) {
      TSTTG_processError(TSTTB_FAILURE, "Failed to scale entity.");
      RETURN(TSTTB_FAILURE);
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
        TSTTG_processError(TSTTB_FAILURE, "Failed to scale body even only one entity of that"
                           " dimension in the body.");
        RETURN(TSTTB_FAILURE);
      }
    }
    else {
      TSTTG_processError(TSTTB_FAILURE, "Too many siblings for an entity to scale it.");
      RETURN(TSTTB_FAILURE);
    }
  }

  RETURN(TSTTB_SUCCESS);
}

enum TSTTB_ErrorType
TSTTG_Unite (TSTTG_Instance instance,
             /*in*/ ARRAY_IN_DECL(TSTTG_GentityHandle, geom_entities),
             /*out*/ TSTTG_GentityHandle *geom_entity)
{
  DLIList<Body*> bods, orig_bods;
  RefEntity **handle_array = ENTITY_HANDLE_ARRAY(geom_entities);
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
    TSTTG_processError(TSTTB_INVALID_ARGUMENT, "Not all entities input were regions.");
    for (int i = bods.size(); i > 0; i--)
      gqt->delete_RefEntity(bods.get_and_step());
    
    RETURN(TSTTB_SUCCESS);
  }
  
  DLIList<Body*> new_bods;
  CubitStatus result = gmt->unite(bods, new_bods, false);
  if (CUBIT_SUCCESS != result || 1 != new_bods.size()) {
    TSTTG_processError(TSTTB_FAILURE, "Unite failed.");
    RETURN(TSTTB_FAILURE);
  }
    
  else {
    *geom_entity = dynamic_cast<RefEntity*>(new_bods.get());
    for (int i = orig_bods.size(); i > 0; i--)
      gqt->delete_RefEntity(orig_bods.get_and_step());
  }

  RETURN(TSTTB_SUCCESS);
}
      
enum TSTTB_ErrorType
TSTTG_Subtract (TSTTG_Instance instance,
                /*in*/ TSTTG_GentityHandle blank,
                /*in*/ TSTTG_GentityHandle tool,
                /*out*/ TSTTG_GentityHandle *geom_entity)
{
  Body *this_blank = dynamic_cast<Body*>(ENTITY_HANDLE(blank));
  Body *blank_copy = gmt->copy_body(this_blank);
  if (NULL == blank_copy) {
    TSTTG_processError(TSTTB_FAILURE, "Trouble copying blank.");
    RETURN(TSTTB_FAILURE);
  }
  Body *this_tool = dynamic_cast<Body*>(ENTITY_HANDLE(tool));
  Body *tool_copy = gmt->copy_body(this_tool);
  if (NULL == tool_copy) {
    TSTTG_processError(TSTTB_FAILURE, "Trouble copying tool.");
    gqt->delete_RefEntity(blank_copy);
    RETURN(TSTTB_FAILURE);
  }

  DLIList<Body*> blank_list, new_body_list;
  blank_list.append(blank_copy);
  
  RefEntity *new_body = NULL;
  CubitStatus result = gmt->subtract(tool_copy, blank_list, new_body_list);
  if (CUBIT_SUCCESS != result || 0 == new_body_list.size()) {
    TSTTG_processError(TSTTB_FAILURE, "Subtract failed.");
    RETURN(TSTTB_FAILURE);
  }
  else {
    new_body = new_body_list.get();
    *geom_entity = new_body;
    gqt->delete_RefEntity(this_blank);
    gqt->delete_RefEntity(this_tool);
  }

  RETURN(TSTTB_SUCCESS);
}

enum TSTTB_ErrorType
TSTTG_Section (TSTTG_Instance instance,
               /*inout*/ TSTTG_GentityHandle *geom_entity,
               /*in*/ const double plane_normal_x,
               /*in*/ const double plane_normal_y,
               /*in*/ const double plane_normal_z,
               /*in*/ const double offset,
               /*in*/ const bool reverse,
               /*out*/ TSTTG_GentityHandle *geom_entity2)
{
  Body *this_body = dynamic_cast<Body*>(ENTITY_HANDLE(*geom_entity));
  if (NULL == this_body) {
    RefVolume *this_vol = dynamic_cast<RefVolume*>(ENTITY_HANDLE(*geom_entity));
    if (NULL != this_vol)
      this_body = this_vol->get_body_ptr();
  }
  if (NULL == this_body) {
    TSTTG_processError(TSTTB_INVALID_ARGUMENT, "Can only section bodies.");
    RETURN(TSTTB_INVALID_ARGUMENT);
  }

  CubitVector normal(plane_normal_x, plane_normal_y, plane_normal_z);
  if (normal.length_squared() == 0.0) {
    TSTTG_processError(TSTTB_INVALID_ARGUMENT, "Zero-length vector input.");
    RETURN(TSTTB_INVALID_ARGUMENT);
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
    TSTTG_processError(TSTTB_FAILURE, "Section failed.");
    gqt->delete_RefEntity(blank_list.get());
    RETURN(TSTTB_FAILURE);
  }
    
  else {
      // need to assign it to a RE* first so the void cast gets done right
    RefEntity *new_body = new_body_list.get();
    *geom_entity2 = new_body;
      // also, delete the original body, now that the section worked
    gqt->delete_RefEntity(this_body);
  }

  RETURN(TSTTB_SUCCESS);
}

enum TSTTB_ErrorType
TSTTG_Imprint (TSTTG_Instance instance,
               /*in*/ ARRAY_IN_DECL(TSTTG_GentityHandle, gentity_handles)) 
{
  DLIList<Body*> bods;
  DLIList<RefVolume*> vols, temp_vols;
  RefEntity **handle_array = ENTITY_HANDLE_ARRAY(gentity_handles);
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
  
  if (CUBIT_SUCCESS != status) RETURN(TSTTB_FAILURE);

  DLIList<Body*> temp_bods;
  status = GeometryModifyTool::instance()->imprint(bods, temp_bods, false);
  
  RETURN(TSTTB_SUCCESS);
}
  
enum TSTTB_ErrorType
TSTTG_Merge (TSTTG_Instance instance,
             /*in*/ ARRAY_IN_DECL(TSTTG_GentityHandle, gentity_handles),
             const double tolerance) 
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
  RefEntity **handle_array = ENTITY_HANDLE_ARRAY(gentity_handles);
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
        if (NULL == temp_bod) RETURN(TSTTB_FAILURE);
        temp_vols.clean_out();
        topo_ent->ref_volumes(temp_vols);
        vols += temp_vols;
        break;
      case 0:
        temp_vert = dynamic_cast<RefVertex*>(handle_array[i]);
        if (NULL == temp_vert) RETURN(TSTTB_FAILURE);
        verts.append(temp_vert);
        break;
      case 1:
        temp_edge = dynamic_cast<RefEdge*>(handle_array[i]);
        if (NULL == temp_edge) RETURN(TSTTB_FAILURE);
        edges.append(temp_edge);
        break;
      case 2:
        temp_face = dynamic_cast<RefFace*>(handle_array[i]);
        if (NULL == temp_face) RETURN(TSTTB_FAILURE);
        faces.append(temp_face);
        break;
      case 3:
        temp_vol = dynamic_cast<RefVolume*>(handle_array[i]);
        if (NULL == temp_vol) RETURN(TSTTB_FAILURE);
        vols.append(temp_vol);
        break;
    }
  }
  
  CubitStatus status = CUBIT_SUCCESS, temp_status;
    
  if (bods.size() != 0) {
    temp_status = MergeTool::instance()->merge_bodies(bods);
    if (CUBIT_SUCCESS != temp_status) status = temp_status;
  }
    
  if (vols.size() != 0) {
    temp_status = MergeTool::instance()->merge_volumes(vols, false);
    if (CUBIT_SUCCESS != temp_status) status = temp_status;
  }
    
  if (faces.size() != 0) {
    temp_status = MergeTool::instance()->merge_reffaces(faces, false);
    if (CUBIT_SUCCESS != temp_status) status = temp_status;
  }
    
  if (edges.size() != 0) {
    temp_status = MergeTool::instance()->merge_refedges(edges, true, false);
    if (CUBIT_SUCCESS != temp_status) status = temp_status;
  }
    
  if (verts.size() != 0) {
    temp_status = MergeTool::instance()->merge_refvertices(verts, false);
    if (CUBIT_SUCCESS != temp_status) status = temp_status;
  }

  if (0 != old_factor)
    GeometryQueryTool::instance()->set_geometry_factor(old_factor);
    
  if (CUBIT_SUCCESS != status) {
    RETURN(TSTTB_FAILURE);
  }
  
  else {
    RETURN(TSTTB_SUCCESS);
  }
}
  
