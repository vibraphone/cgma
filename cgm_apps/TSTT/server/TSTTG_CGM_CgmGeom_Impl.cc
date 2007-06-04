// 
// File:          TSTTG_CGM_CgmGeom_Impl.cc
// Symbol:        TSTTG_CGM.CgmGeom-v0.1
// Symbol Type:   class
// Babel Version: 0.9.8
// sidl Created:  20051011 16:26:00 GMT-06:00
// Generated:     20051011 16:26:03 GMT-06:00
// Description:   Server-side implementation for TSTTG_CGM.CgmGeom
// 
// WARNING: Automatically generated; only changes within splicers preserved
// 
// babel-version = 0.9.8
// source-line   = 5
// source-url    = file:/home/tjtautg/tstt/cvs/TSTTG/TSTTG_CGM.sidl
// 
#include "TSTTG_CGM_CgmGeom_Impl.hh"

// DO-NOT-DELETE splicer.begin(TSTTG_CGM.CgmGeom._includes)
// Put additional includes or other arbitrary code here...
#include "sidlArray.h"
#include "TSTTB.hh"
#include "TSTTG_CGM.h"
#include "TSTTG_GentityType_IOR.h"
#include "TSTTB_Error.hh"
#include <iostream>

TSTTG_Instance tsttgInstance = NULL;
#define PROCESS_ERROR if (TSTTG_LAST_ERROR.error_type != TSTTB_SUCCESS) this->processError()
#define GENERATE_ERROR(a) {TSTTG_LAST_ERROR.error_type = a; this->processError();}

// need this for definitions in TSTTB_SNL_SIDL_defs.h
#define LOCAL_TSTTB_ERROR TSTTG_LAST_ERROR
#include "TSTTB_SNL_SIDL_defs.h"
#include <string>
// DO-NOT-DELETE splicer.end(TSTTG_CGM.CgmGeom._includes)

// user defined constructor
void TSTTG_CGM::CgmGeom_impl::_ctor() {
    // DO-NOT-DELETE splicer.begin(TSTTG_CGM.CgmGeom._ctor)
    // add construction details here
  TSTTG_ctor(&tsttgInstance);
  PROCESS_ERROR;
    // DO-NOT-DELETE splicer.end(TSTTG_CGM.CgmGeom._ctor)
}

// user defined destructor
void TSTTG_CGM::CgmGeom_impl::_dtor() {
    // DO-NOT-DELETE splicer.begin(TSTTG_CGM.CgmGeom._dtor)
    // add destruction details here
  TSTTG_dtor(tsttgInstance);
  PROCESS_ERROR;
    // DO-NOT-DELETE splicer.end(TSTTG_CGM.CgmGeom._dtor)
}

// user defined static methods: (none)

// user defined non-static methods:
/**
 * Return gentities of specified dimension in this set, or in whole model.
 * @param set_handle Entity set being queried (if 0, whole model)
 * @param gentity_dimension Dimension of entities being queried
 * @param gentity_handles Gentity handles
 */
void
TSTTG_CGM::CgmGeom_impl::gentitysetGetGentitiesOfType (
  /*in*/ void* set_handle,
  /*in*/ ::TSTTG::GentityType gentity_type,
  /*inout*/ ::sidl::array<void*>& gentity_handles,
  /*out*/ int32_t& gentity_handles_size ) 
throw ( 
  ::TSTTB::Error
){
    // DO-NOT-DELETE splicer.begin(TSTTG_CGM.CgmGeom.gentitysetGetGentitiesOfType)
    // insert implementation here
  if (gentity_type > ::TSTTG::GentityType_GALL_TYPES || 
      gentity_type < ::TSTTG::GentityType_GVERTEX) {
    GENERATE_ERROR(TSTTB_INVALID_ENTITY_TYPE);
    return;
  }
    
  CREATE_TEMP_ARRAY(void*, gentity_handles);

  TSTTG_gentitysetGetGentitiesOfType (tsttgInstance, 
                                      set_handle,
                                      (TSTTG_GentityType)gentity_type,
                                      TEMP_ARRAY_INOUT(gentity_handles));
  PROCESS_ERROR;

  ASSIGN_ARRAY(gentity_handles);
    // DO-NOT-DELETE splicer.end(TSTTG_CGM.CgmGeom.gentitysetGetGentitiesOfType)
}

/**
 * Return number of gentities of specified dimension in this set, or in
 * whole model.
 * @param set_handle Entity set being queried (if 0, whole model)
 * @param gentity_dimension Dimension of entities being queried
 * @return Number of entities
 */
int32_t
TSTTG_CGM::CgmGeom_impl::gentitysetGetNumberGentitiesOfType (
  /*in*/ void* set_handle,
  /*in*/ ::TSTTG::GentityType gentity_type ) 
throw ( 
  ::TSTTB::Error
){
    // DO-NOT-DELETE splicer.begin(TSTTG_CGM.CgmGeom.gentitysetGetNumberGentitiesOfType)
    // insert implementation here
  if (gentity_type > ::TSTTG::GentityType_GALL_TYPES || 
      gentity_type < ::TSTTG::GentityType_GVERTEX) {
    GENERATE_ERROR(TSTTB_INVALID_ENTITY_TYPE);
    return -1;
  }
  int32_t retval = TSTTG_gentitysetGetNumberGentitiesOfType (tsttgInstance, 
                                                             set_handle,
                                                             (TSTTG_GentityType)gentity_type);
  PROCESS_ERROR;
  return retval;
    // DO-NOT-DELETE splicer.end(TSTTG_CGM.CgmGeom.gentitysetGetNumberGentitiesOfType)
}

/**
 *    Returns an integer array of topological dimensions for an input
 *    array of entity handles.
 */
void
TSTTG_CGM::CgmGeom_impl::gentityGetType (
  /*in*/ ::sidl::array<void*> gentity_handles,
  /*in*/ int32_t gentity_handles_size,
  /*inout*/ ::sidl::array< ::TSTTG::GentityType>& gtype,
  /*inout*/ int32_t& gtype_size ) 
throw ( 
  ::TSTTB::Error
){
    // DO-NOT-DELETE splicer.begin(TSTTG_CGM.CgmGeom.gentityGetType)
    // insert implementation here
  CREATE_TEMP_ENUM_ARRAY(TSTTG_GentityType, gtype);
  
  TSTTG_gentityGetType (tsttgInstance, 
                        (const void**) TEMP_ARRAY_IN(gentity_handles), 
                        TEMP_ARRAY_INOUT(gtype));
  PROCESS_ERROR;

  ASSIGN_ENUM_ARRAY(gtype);
    // DO-NOT-DELETE splicer.end(TSTTG_CGM.CgmGeom.gentityGetType)
}

/**
 * Get the adjacent entities of a given dimension.
 * @param gentity_handle Entity for which adjacencies are requested
 * @param to_dimension Target dimension of adjacent entities
 * @param adj_gentities List returned with adjacent entities
 */
void
TSTTG_CGM::CgmGeom_impl::gentityGetAdjacencies (
  /*in*/ void* gentity_handle,
  /*in*/ int32_t to_dimension,
  /*inout*/ ::sidl::array<void*>& adj_gentities,
  /*inout*/ int32_t& adj_gentities_size ) 
throw ( 
  ::TSTTB::Error
){
    // DO-NOT-DELETE splicer.begin(TSTTG_CGM.CgmGeom.gentityGetAdjacencies)
    // insert implementation here
  CREATE_TEMP_ARRAY(void*, adj_gentities);
  
  TSTTG_gentityGetAdjacencies (tsttgInstance, 
                               gentity_handle,
                               to_dimension,
                               TEMP_ARRAY_INOUT(adj_gentities));
  PROCESS_ERROR;

  ASSIGN_ARRAY(adj_gentities);
    // DO-NOT-DELETE splicer.end(TSTTG_CGM.CgmGeom.gentityGetAdjacencies)
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
TSTTG_CGM::CgmGeom_impl::gentityGet2OAdjacencies (
  /*in*/ void* gentity_handle,
  /*in*/ int32_t bridge_dimension,
  /*in*/ int32_t to_dimension,
  /*inout*/ ::sidl::array<void*>& adjacent_gentities,
  /*out*/ int32_t& adjacent_gentities_size ) 
throw ( 
  ::TSTTB::Error
){
    // DO-NOT-DELETE splicer.begin(TSTTG_CGM.CgmGeom.gentityGet2OAdjacencies)
    // insert implementation here
  CREATE_TEMP_ARRAY(void*, adjacent_gentities);
  
  TSTTG_gentityGet2OAdjacencies (tsttgInstance, 
                                 gentity_handle,
                                 bridge_dimension,
                                 to_dimension,
                                 TEMP_ARRAY_INOUT(adjacent_gentities));
  PROCESS_ERROR;
  
  ASSIGN_ARRAY(adjacent_gentities);
  
    // DO-NOT-DELETE splicer.end(TSTTG_CGM.CgmGeom.gentityGet2OAdjacencies)
}

/**
 * Return whether or not entities are adjacent.
 * @param gentity_handle1 1st entity
 * @param gentity_handle2 2nd entity
 * @param are_adjacent If true, entities are adjacent
 */
void
TSTTG_CGM::CgmGeom_impl::gentityIsAdjacent (
  /*in*/ void* gentity_handle1,
  /*in*/ void* gentity_handle2,
  /*out*/ bool& are_adjacent ) 
throw ( 
  ::TSTTB::Error
){
    // DO-NOT-DELETE splicer.begin(TSTTG_CGM.CgmGeom.gentityIsAdjacent)
    // insert implementation here
  TSTTG_gentityIsAdjacent (tsttgInstance, 
                           gentity_handle1,
                           gentity_handle2,
                           &are_adjacent);
  PROCESS_ERROR;
    // DO-NOT-DELETE splicer.end(TSTTG_CGM.CgmGeom.gentityIsAdjacent)
}

/**
 * Load a model specified by name. Which formats are supported and the
 * specific meaning of this name string (e.g. file name, model name,
 * etc.) are implementation-dependent.  Options are also implementation-
 * dependent.
 * @param name Name of the model
 * @param options String options 
 */
void
TSTTG_CGM::CgmGeom_impl::gLoad (
  /*in*/ const ::std::string& name,
  /*in*/ ::sidl::array< ::std::string> options,
  /*in*/ int32_t options_size ) 
throw ( 
  ::TSTTB::Error
){
    // DO-NOT-DELETE splicer.begin(TSTTG_CGM.CgmGeom.gLoad)
    // insert implementation here
  const char **tmp_options = NULL;
  if (options_size > 0) {
    tmp_options = new const char*[options_size];
    for (int i = 0; i < options_size; i++)
      tmp_options[i] = (options.get(i)).c_str();
  }
  
  TSTTG_gLoad(tsttgInstance, name.c_str(), tmp_options, 
              options_size);
  
  delete [] tmp_options;

  PROCESS_ERROR;
    // DO-NOT-DELETE splicer.end(TSTTG_CGM.CgmGeom.gLoad)
}

/**
 * Save a model to a file specified by name. Which formats are supported and the
 * specific meaning of this name string (e.g. file name, model name,
 * etc.) are implementation-dependent.  Options are also implementation-
 * dependent.
 * @param name Name of the file to save to
 * @param options String options 
 */
void
TSTTG_CGM::CgmGeom_impl::gSave (
  /*in*/ const ::std::string& name,
  /*in*/ ::sidl::array< ::std::string> options,
  /*in*/ int32_t options_size ) 
throw ( 
  ::TSTTB::Error
){
  // DO-NOT-DELETE splicer.begin(TSTTG_CGM.CgmGeom.gSave)
  // insert implementation here
  const char **tmp_options = NULL;
  if (options_size > 0) {
    tmp_options = new const char*[options_size];
    for (int i = 0; i < options_size; i++)
      tmp_options[i] = (options.get(i)).c_str();
  }
  
  TSTTG_gSave(tsttgInstance, name.c_str(), tmp_options, 
              options_size);
  
  delete [] tmp_options;

  PROCESS_ERROR;
  // DO-NOT-DELETE splicer.end(TSTTG_CGM.CgmGeom.gSave)
}

/**
 * Method:  createTag[]
 */
void
TSTTG_CGM::CgmGeom_impl::createTag (
  /*in*/ const ::std::string& tag_name,
  /*in*/ int32_t number_of_values,
  /*in*/ ::TSTTB::TagValueType tag_type,
  /*out*/ void*& tag_handle ) 
throw ( 
  ::TSTTB::Error
){
    // DO-NOT-DELETE splicer.begin(TSTTG_CGM.CgmGeom.createTag)
    // insert implementation here
  if (tag_type > ::TSTTB::TagValueType_BYTES || 
      tag_type < ::TSTTB::TagValueType_INTEGER) {
    GENERATE_ERROR(TSTTB_INVALID_ARGUMENT);
    return;
  }
  TSTTG_createTag(tsttgInstance, tag_name.c_str(), number_of_values, (TSTTG_TagValueType)tag_type, 
                  &tag_handle);
  PROCESS_ERROR;
    // DO-NOT-DELETE splicer.end(TSTTG_CGM.CgmGeom.createTag)
}

/**
 * Method:  destroyTag[]
 */
void
TSTTG_CGM::CgmGeom_impl::destroyTag (
  /*in*/ void* tag_handle,
  /*in*/ bool forced ) 
throw ( 
  ::TSTTB::Error
){
    // DO-NOT-DELETE splicer.begin(TSTTG_CGM.CgmGeom.destroyTag)
    // insert implementation here
  TSTTG_destroyTag(tsttgInstance, tag_handle, forced);
  PROCESS_ERROR;
    // DO-NOT-DELETE splicer.end(TSTTG_CGM.CgmGeom.destroyTag)
}

/**
 * Method:  getTagName[]
 */
::std::string
TSTTG_CGM::CgmGeom_impl::getTagName (
  /*in*/ void* tag_handle ) 
throw ( 
  ::TSTTB::Error
){
    // DO-NOT-DELETE splicer.begin(TSTTG_CGM.CgmGeom.getTagName)
    // insert implementation here
  std::string retval = TSTTG_getTagName(tsttgInstance, tag_handle);
  PROCESS_ERROR;
  return retval;
    // DO-NOT-DELETE splicer.end(TSTTG_CGM.CgmGeom.getTagName)
}

/**
 * Method:  getTagSizeValues[]
 */
int32_t
TSTTG_CGM::CgmGeom_impl::getTagSizeValues (
  /*in*/ void* tag_handle ) 
throw ( 
  ::TSTTB::Error
){
  // DO-NOT-DELETE splicer.begin(TSTTG_CGM.CgmGeom.getTagSizeValues)
  // insert implementation here
  int32_t retval = TSTTG_getTagSizeValues(tsttgInstance, tag_handle);
  PROCESS_ERROR;
  return retval;
  // DO-NOT-DELETE splicer.end(TSTTG_CGM.CgmGeom.getTagSizeValues)
}

/**
 * Method:  getTagSizeBytes[]
 */
int32_t
TSTTG_CGM::CgmGeom_impl::getTagSizeBytes (
  /*in*/ void* tag_handle ) 
throw ( 
  ::TSTTB::Error
){
  // DO-NOT-DELETE splicer.begin(TSTTG_CGM.CgmGeom.getTagSizeBytes)
  // insert implementation here
  int32_t retval = TSTTG_getTagSizeBytes(tsttgInstance, tag_handle);
  PROCESS_ERROR;
  return retval;
  // DO-NOT-DELETE splicer.end(TSTTG_CGM.CgmGeom.getTagSizeBytes)
}

/**
 * Method:  getTagHandle[]
 */
void*
TSTTG_CGM::CgmGeom_impl::getTagHandle (
  /*in*/ const ::std::string& tag_name ) 
throw ( 
  ::TSTTB::Error
){
    // DO-NOT-DELETE splicer.begin(TSTTG_CGM.CgmGeom.getTagHandle)
    // insert implementation here
  void *retval = TSTTG_getTagHandle(tsttgInstance, tag_name.c_str());
  PROCESS_ERROR;
  return retval;
    // DO-NOT-DELETE splicer.end(TSTTG_CGM.CgmGeom.getTagHandle)
}

/**
 * Method:  getTagType[]
 */
::TSTTB::TagValueType
TSTTG_CGM::CgmGeom_impl::getTagType (
  /*in*/ void* tag_handle ) 
throw ( 
  ::TSTTB::Error
){
    // DO-NOT-DELETE splicer.begin(TSTTG_CGM.CgmGeom.getTagType)
    // insert implementation here
  ::TSTTB::TagValueType retval = (::TSTTB::TagValueType) TSTTG_getTagType(tsttgInstance,
                                                                        tag_handle );
  PROCESS_ERROR;
  return retval;
    // DO-NOT-DELETE splicer.end(TSTTG_CGM.CgmGeom.getTagType)
}

/**
 * Method:  getArrData[]
 */
void
TSTTG_CGM::CgmGeom_impl::getArrData (
  /*in*/ ::sidl::array<void*> entity_handles,
  /*in*/ int32_t entity_handles_size,
  /*in*/ void* tag_handle,
  /*inout*/ ::sidl::array<char>& value_array,
  /*out*/ int32_t& value_array_size ) 
throw ( 
  ::TSTTB::Error
){
    // DO-NOT-DELETE splicer.begin(TSTTG_CGM.CgmGeom.getArrData)
    // insert implementation here
  CREATE_TEMP_TAG_ARRAY(value_array);
  
  TSTTG_getArrData(tsttgInstance, 
                   (const void **) TEMP_ARRAY_IN(entity_handles), 
                   tag_handle, 
                   TEMP_ARRAY_INOUT(value_array));
  PROCESS_ERROR;

  ASSIGN_TAG_ARRAY(value_array);
    // DO-NOT-DELETE splicer.end(TSTTG_CGM.CgmGeom.getArrData)
}

/**
 * Method:  getIntArrData[]
 */
void
TSTTG_CGM::CgmGeom_impl::getIntArrData (
  /*in*/ ::sidl::array<void*> entity_handles,
  /*in*/ int32_t entity_handles_size,
  /*in*/ void* tag_handle,
  /*inout*/ ::sidl::array<int32_t>& value_array,
  /*out*/ int32_t& value_array_size ) 
throw ( 
  ::TSTTB::Error
){
    // DO-NOT-DELETE splicer.begin(TSTTG_CGM.CgmGeom.getIntArrData)
    // insert implementation here
  CREATE_TEMP_ARRAY(int32_t, value_array);
  
  TSTTG_getIntArrData (tsttgInstance,
                       (const void **) TEMP_ARRAY_IN(entity_handles),
                       tag_handle,
                       TEMP_ARRAY_INOUT(value_array));
  PROCESS_ERROR;

  ASSIGN_ARRAY(value_array);
    // DO-NOT-DELETE splicer.end(TSTTG_CGM.CgmGeom.getIntArrData)
}

/**
 * Method:  getDblArrData[]
 */
void
TSTTG_CGM::CgmGeom_impl::getDblArrData (
  /*in*/ ::sidl::array<void*> entity_handles,
  /*in*/ int32_t entity_handles_size,
  /*in*/ void* tag_handle,
  /*inout*/ ::sidl::array<double>& value_array,
  /*out*/ int32_t& value_array_size ) 
throw ( 
  ::TSTTB::Error
){
    // DO-NOT-DELETE splicer.begin(TSTTG_CGM.CgmGeom.getDblArrData)
    // insert implementation here
  CREATE_TEMP_ARRAY(double, value_array);
  
  TSTTG_getDblArrData (tsttgInstance,
                       (const void **) TEMP_ARRAY_IN(entity_handles),
                       tag_handle,
                       TEMP_ARRAY_INOUT(value_array));
  PROCESS_ERROR;

  ASSIGN_ARRAY(value_array);
    // DO-NOT-DELETE splicer.end(TSTTG_CGM.CgmGeom.getDblArrData)
}

/**
 * Method:  getEHArrData[]
 */
void
TSTTG_CGM::CgmGeom_impl::getEHArrData (
  /*in*/ ::sidl::array<void*> entity_handles,
  /*in*/ int32_t entity_handles_size,
  /*in*/ void* tag_handle,
  /*inout*/ ::sidl::array<void*>& value_array,
  /*out*/ int32_t& value_array_size ) 
throw ( 
  ::TSTTB::Error
){
    // DO-NOT-DELETE splicer.begin(TSTTG_CGM.CgmGeom.getEHArrData)
    // insert implementation here

  CREATE_TEMP_ARRAY(void*, value_array);
  
  TSTTG_getEHArrData (tsttgInstance,
                      (const void**) TEMP_ARRAY_IN(entity_handles),
                      tag_handle,
                      TEMP_ARRAY_INOUT(value_array));
  PROCESS_ERROR;

  ASSIGN_ARRAY(value_array);
    // DO-NOT-DELETE splicer.end(TSTTG_CGM.CgmGeom.getEHArrData)
}

/**
 * Method:  setArrData[]
 */
void
TSTTG_CGM::CgmGeom_impl::setArrData (
  /*in*/ ::sidl::array<void*> entity_handles,
  /*in*/ int32_t entity_handles_size,
  /*in*/ void* tag_handle,
  /*in*/ ::sidl::array<char> value_array,
  /*in*/ int32_t value_array_size ) 
throw ( 
  ::TSTTB::Error
){
    // DO-NOT-DELETE splicer.begin(TSTTG_CGM.CgmGeom.setArrData)
    // insert implementation here
  TSTTG_setArrData(tsttgInstance, 
                   TEMP_ARRAY_IN(entity_handles),
                   tag_handle, 
                   TEMP_TAG_ARRAY_IN(value_array));
  PROCESS_ERROR;
    // DO-NOT-DELETE splicer.end(TSTTG_CGM.CgmGeom.setArrData)
}

/**
 * Method:  setIntArrData[]
 */
void
TSTTG_CGM::CgmGeom_impl::setIntArrData (
  /*in*/ ::sidl::array<void*> entity_handles,
  /*in*/ int32_t entity_handles_size,
  /*in*/ void* tag_handle,
  /*in*/ ::sidl::array<int32_t> value_array,
  /*in*/ int32_t value_array_size ) 
throw ( 
  ::TSTTB::Error
){
    // DO-NOT-DELETE splicer.begin(TSTTG_CGM.CgmGeom.setIntArrData)
    // insert implementation here
  TSTTG_setIntArrData (tsttgInstance,
                       TEMP_ARRAY_IN(entity_handles),
                       tag_handle,
                       TEMP_ARRAY_IN(value_array));
  PROCESS_ERROR;
    // DO-NOT-DELETE splicer.end(TSTTG_CGM.CgmGeom.setIntArrData)
}

/**
 * Method:  setDblArrData[]
 */
void
TSTTG_CGM::CgmGeom_impl::setDblArrData (
  /*in*/ ::sidl::array<void*> entity_handles,
  /*in*/ int32_t entity_handles_size,
  /*in*/ void* tag_handle,
  /*in*/ ::sidl::array<double> value_array,
  /*in*/ int32_t value_array_size ) 
throw ( 
  ::TSTTB::Error
){
    // DO-NOT-DELETE splicer.begin(TSTTG_CGM.CgmGeom.setDblArrData)
    // insert implementation here
  TSTTG_setDblArrData (tsttgInstance,
                       TEMP_ARRAY_IN(entity_handles),
                       tag_handle,
                       TEMP_ARRAY_IN(value_array));
  PROCESS_ERROR;
    // DO-NOT-DELETE splicer.end(TSTTG_CGM.CgmGeom.setDblArrData)
}

/**
 * Method:  setEHArrData[]
 */
void
TSTTG_CGM::CgmGeom_impl::setEHArrData (
  /*in*/ ::sidl::array<void*> entity_handles,
  /*in*/ int32_t entity_handles_size,
  /*in*/ void* tag_handle,
  /*in*/ ::sidl::array<void*> value_array,
  /*in*/ int32_t value_array_size ) 
throw ( 
  ::TSTTB::Error
){
    // DO-NOT-DELETE splicer.begin(TSTTG_CGM.CgmGeom.setEHArrData)
    // insert implementation here
  TSTTG_setEHArrData (tsttgInstance,
                      TEMP_ARRAY_IN(entity_handles),
                      tag_handle,
                      (const void **) TEMP_ARRAY_IN(value_array));
  PROCESS_ERROR;
    // DO-NOT-DELETE splicer.end(TSTTG_CGM.CgmGeom.setEHArrData)
}

/**
 * Method:  rmvArrTag[]
 */
void
TSTTG_CGM::CgmGeom_impl::rmvArrTag (
  /*in*/ ::sidl::array<void*> entity_handles,
  /*in*/ int32_t entity_handles_size,
  /*in*/ void* tag_handle ) 
throw ( 
  ::TSTTB::Error
){
    // DO-NOT-DELETE splicer.begin(TSTTG_CGM.CgmGeom.rmvArrTag)
    // insert implementation here
  TSTTG_rmvArrTag (tsttgInstance,
                   TEMP_ARRAY_IN(entity_handles),
                   tag_handle );
  PROCESS_ERROR;
    // DO-NOT-DELETE splicer.end(TSTTG_CGM.CgmGeom.rmvArrTag)
}

/**
 * Initialize an iterator over gentities of a specified dimension.
 * @param gentity_dimension Dimension of gentities to be iterated over
 * @param gentity_iterator Iterator initialized by this function
 */
void
TSTTG_CGM::CgmGeom_impl::gentityIteratorInit (
  /*in*/ int32_t gentity_dimension,
  /*out*/ void*& gentity_iterator ) 
throw ( 
  ::TSTTB::Error
){
    // DO-NOT-DELETE splicer.begin(TSTTG_CGM.CgmGeom.gentityIteratorInit)
    // insert implementation here
  TSTTG_gentityIteratorInit (tsttgInstance, gentity_dimension, &gentity_iterator);
  PROCESS_ERROR;
    // DO-NOT-DELETE splicer.end(TSTTG_CGM.CgmGeom.gentityIteratorInit)
}

/**
 * Get the next entity for this iterator.
 * @param gentity_iterator Iterator being iterated over
 * @param gentity_handle Next gentity
 * @return If true, there are more gentities, if false, this is the last one
 */
bool
TSTTG_CGM::CgmGeom_impl::gentityIteratorNext (
  /*inout*/ void*& gentity_iterator,
  /*out*/ void*& gentity_handle ) 
throw ( 
  ::TSTTB::Error
){
    // DO-NOT-DELETE splicer.begin(TSTTG_CGM.CgmGeom.gentityIteratorNext)
    // insert implementation here
  bool retval = TSTTG_gentityIteratorNext (tsttgInstance, 
                                           &gentity_iterator,
                                           &gentity_handle);
  PROCESS_ERROR;
  return retval;
    // DO-NOT-DELETE splicer.end(TSTTG_CGM.CgmGeom.gentityIteratorNext)
}

/**
 * Reset an iterator back to the first gentity
 * @param gentity_iterator Iterator reset by this function
 */
void
TSTTG_CGM::CgmGeom_impl::gentityIteratorReset (
  /*inout*/ void*& gentity_iterator ) 
throw ( 
  ::TSTTB::Error
){
    // DO-NOT-DELETE splicer.begin(TSTTG_CGM.CgmGeom.gentityIteratorReset)
    // insert implementation here
  TSTTG_gentityIteratorReset (tsttgInstance, 
                              &gentity_iterator);
  PROCESS_ERROR;
    // DO-NOT-DELETE splicer.end(TSTTG_CGM.CgmGeom.gentityIteratorReset)
}

/**
 * Delete an iterator
 * @param gentity_iterator Iterator deleted by this function
 */
void
TSTTG_CGM::CgmGeom_impl::gentityIteratorDelete (
  /*in*/ void* Gentity_dim_iterator ) 
throw ( 
  ::TSTTB::Error
){
    // DO-NOT-DELETE splicer.begin(TSTTG_CGM.CgmGeom.gentityIteratorDelete)
    // insert implementation here
  TSTTG_gentityIteratorDelete (tsttgInstance, 
                               Gentity_dim_iterator);
  PROCESS_ERROR;
    // DO-NOT-DELETE splicer.end(TSTTG_CGM.CgmGeom.gentityIteratorDelete)
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
TSTTG_CGM::CgmGeom_impl::gentityClosestPoint (
  /*in*/ ::sidl::array<void*> gentity_handles,
  /*in*/ int32_t gentity_handles_size,
  /*in*/ ::sidl::array<double> near_coordinates,
  /*in*/ int32_t near_coordinates_size,
  /*inout*/ ::sidl::array<double>& on_coordinates,
  /*out*/ int32_t& on_coordinates_size ) 
throw ( 
  ::TSTTB::Error
){
    // DO-NOT-DELETE splicer.begin(TSTTG_CGM.CgmGeom.gentityClosestPoint)
    // insert implementation here
  CREATE_TEMP_ARRAY(double, on_coordinates);
  
  TSTTG_gentityClosestPoint (tsttgInstance, 
                             (const void **) TEMP_ARRAY_IN(gentity_handles),
                             TEMP_ARRAY_IN(near_coordinates),
                             TEMP_ARRAY_INOUT(on_coordinates));
  PROCESS_ERROR;

  ASSIGN_ARRAY(on_coordinates);
    // DO-NOT-DELETE splicer.end(TSTTG_CGM.CgmGeom.gentityClosestPoint)
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
TSTTG_CGM::CgmGeom_impl::gentityNormal (
  /*in*/ ::sidl::array<void*> gentity_handles,
  /*in*/ int32_t gentity_handles_size,
  /*in*/ ::sidl::array<double> coordinates,
  /*in*/ int32_t coordinates_size,
  /*inout*/ ::sidl::array<double>& normals,
  /*out*/ int32_t& normals_size ) 
throw ( 
  ::TSTTB::Error
){
    // DO-NOT-DELETE splicer.begin(TSTTG_CGM.CgmGeom.gentityNormal)
    // insert implementation here
  CREATE_TEMP_ARRAY(double, normals);
  
  TSTTG_gentityNormal (tsttgInstance, 
                       (const void **) TEMP_ARRAY_IN(gentity_handles),
                       TEMP_ARRAY_IN(coordinates),
                       TEMP_ARRAY_INOUT(normals));
  PROCESS_ERROR;

  ASSIGN_ARRAY(normals);
    // DO-NOT-DELETE splicer.end(TSTTG_CGM.CgmGeom.gentityNormal)
}

/**
 * Return points and normals on specified entities closest to specified points
 * in space.  Input coordinates and output points are interleaved in 
 * the arrays.
 * @param gentity_handles The gentities being queried
 * @param near_coordinates Input coordinates
 * @param on_coordinates Closest point on gentity
 * @param normals Normals on gentity
 */
void
TSTTG_CGM::CgmGeom_impl::gentityClosestPointAndNormal (
  /*in*/ ::sidl::array<void*> gentity_handles,
  /*in*/ int32_t gentity_handles_size,
  /*in*/ ::sidl::array<double> near_coordinates,
  /*in*/ int32_t near_coordinates_size,
  /*inout*/ ::sidl::array<double>& on_coordinates,
  /*out*/ int32_t& on_coordinates_size,
  /*inout*/ ::sidl::array<double>& normals,
  /*out*/ int32_t& normals_size ) 
throw ( 
  ::TSTTB::Error
){
  // DO-NOT-DELETE splicer.begin(TSTTG_CGM.CgmGeom.gentityClosestPointAndNormal)
  // insert implementation here
  CREATE_TEMP_ARRAY(double, on_coordinates);
  CREATE_TEMP_ARRAY(double, normals);
  
  TSTTG_gentityClosestPointAndNormal (tsttgInstance, 
                             (const void **) TEMP_ARRAY_IN(gentity_handles),
                             TEMP_ARRAY_IN(near_coordinates),
                             TEMP_ARRAY_INOUT(on_coordinates),
                             TEMP_ARRAY_INOUT(normals));
  PROCESS_ERROR;

  ASSIGN_ARRAY(on_coordinates);
  ASSIGN_ARRAY(normals);
  // DO-NOT-DELETE splicer.end(TSTTG_CGM.CgmGeom.gentityClosestPointAndNormal)
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
TSTTG_CGM::CgmGeom_impl::gentityTangent (
  /*in*/ ::sidl::array<void*> gentity_handles,
  /*in*/ int32_t gentity_handles_size,
  /*in*/ ::sidl::array<double> coordinates,
  /*in*/ int32_t coordinates_size,
  /*inout*/ ::sidl::array<double>& tangents,
  /*out*/ int32_t& tangents_size ) 
throw ( 
  ::TSTTB::Error
){
    // DO-NOT-DELETE splicer.begin(TSTTG_CGM.CgmGeom.gentityTangent)
    // insert implementation here
  CREATE_TEMP_ARRAY(double, tangents);
  
  TSTTG_gentityTangent (tsttgInstance, 
                        (const void **) TEMP_ARRAY_IN(gentity_handles),
                        TEMP_ARRAY_IN(coordinates),
                        TEMP_ARRAY_INOUT(tangents));
  PROCESS_ERROR;

  ASSIGN_ARRAY(tangents);
    // DO-NOT-DELETE splicer.end(TSTTG_CGM.CgmGeom.gentityTangent)
}

/**
 * Return the bounding boxex of given entities; coordinates returned
 * interleaved.
 * @param gentity_handles The gentities being queried
 * @param min_corners Minimum corner coordinates of the boxes, interleaved
 * @param max_corners Maximum corner coordinates of the boxes, interleaved
 */
void
TSTTG_CGM::CgmGeom_impl::gentityBoundingBox (
  /*in*/ ::sidl::array<void*> gentity_handles,
  /*in*/ int32_t gentity_handles_size,
  /*inout*/ ::sidl::array<double>& min_corner,
  /*out*/ int32_t& min_corner_size,
  /*inout*/ ::sidl::array<double>& max_corner,
  /*out*/ int32_t& max_corner_size ) 
throw ( 
  ::TSTTB::Error
){
    // DO-NOT-DELETE splicer.begin(TSTTG_CGM.CgmGeom.gentityBoundingBox)
    // insert implementation here
  CREATE_TEMP_ARRAY(double, max_corner);
  CREATE_TEMP_ARRAY(double, min_corner);
  
  TSTTG_gentityBoundingBox (tsttgInstance, 
                            TEMP_ARRAY_IN(gentity_handles),
                            TEMP_ARRAY_INOUT(min_corner),
                            TEMP_ARRAY_INOUT(max_corner));
  PROCESS_ERROR;

  ASSIGN_ARRAY(min_corner);
  ASSIGN_ARRAY(max_corner);
    // DO-NOT-DELETE splicer.end(TSTTG_CGM.CgmGeom.gentityBoundingBox)
}

/**
 * Return the coordinates of the specified vertices; returns error if any
 * of the entities are not gvertices.  Coordinates returned interleaved.
 * @param gentity_handles The gentities being queried
 * @param coordinates The coordinates of the gvertices, interleaved.
 */
void
TSTTG_CGM::CgmGeom_impl::getGvertexCoordinates (
  /*in*/ ::sidl::array<void*> gentity_handles,
  /*in*/ int32_t gentity_handles_size,
  /*inout*/ ::sidl::array<double>& coordinates,
  /*out*/ int32_t& coordinates_size ) 
throw ( 
  ::TSTTB::Error
){
    // DO-NOT-DELETE splicer.begin(TSTTG_CGM.CgmGeom.getGvertexCoordinates)
    // insert implementation here
  CREATE_TEMP_ARRAY(double, coordinates);
  
  TSTTG_getGvertexCoordinates (tsttgInstance, 
                               (const void **) TEMP_ARRAY_IN(gentity_handles),
                               TEMP_ARRAY_INOUT(coordinates));
  PROCESS_ERROR;

  ASSIGN_ARRAY(coordinates);
    // DO-NOT-DELETE splicer.end(TSTTG_CGM.CgmGeom.getGvertexCoordinates)
}

/**
 * Return the sense of a gface with respect to a gregion.  Sense is either
 * forward (=1), reverse (=-1), both (=2), or unknown (=0).  Error is returned
 * if first entity is not a gface or second entity is not a gregion.
 * @param gface Gface whose sense is being queried.
 * @param gregion Gregion gface is being queried with respect to
 */
int32_t
TSTTG_CGM::CgmGeom_impl::getGnormalSense (
  /*in*/ void* gface,
  /*in*/ void* gregion ) 
throw ( 
  ::TSTTB::Error
){
    // DO-NOT-DELETE splicer.begin(TSTTG_CGM.CgmGeom.getGnormalSense)
    // insert implementation here
  int32_t retval = TSTTG_getGnormalSense (tsttgInstance, 
                                          gface,
                                          gregion);
  PROCESS_ERROR;
  return retval;
    // DO-NOT-DELETE splicer.end(TSTTG_CGM.CgmGeom.getGnormalSense)
}

/**
 * Return the sense of a gedge with respect to a gface.  Sense is either
 * forward (=1), reverse (=-1), both (=2), or unknown (=0).  Error is returned
 * if first entity is not a gedge or second entity is not a gface.
 * @param gedge Gedge whose sense is being queried.
 * @param gface Gface gedge is being queried with respect to
 */
int32_t
TSTTG_CGM::CgmGeom_impl::getGtangentSense (
  /*in*/ void* gedge,
  /*in*/ void* gface ) 
throw ( 
  ::TSTTB::Error
){
    // DO-NOT-DELETE splicer.begin(TSTTG_CGM.CgmGeom.getGtangentSense)
    // insert implementation here
  int32_t retval = TSTTG_getGtangentSense (tsttgInstance, 
                                           gedge,
                                           gface);
  PROCESS_ERROR;
  return retval;
    // DO-NOT-DELETE splicer.end(TSTTG_CGM.CgmGeom.getGtangentSense)
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
int32_t
TSTTG_CGM::CgmGeom_impl::getGvertexTangentSense (
  /*in*/ void* gedge,
  /*in*/ void* gvertex1,
  /*in*/ void* gvertex2 ) 
throw ( 
  ::TSTTB::Error
){
    // DO-NOT-DELETE splicer.begin(TSTTG_CGM.CgmGeom.getGvertexTangentSense)
    // insert implementation here
  int32_t retval = TSTTG_getGvertexTangentSense (tsttgInstance, 
                                                 gedge,
                                                 gvertex1,
                                                 gvertex2);
  PROCESS_ERROR;
  return retval;
    // DO-NOT-DELETE splicer.end(TSTTG_CGM.CgmGeom.getGvertexTangentSense)
}

/**
 * Return whether a given gentity is parametric or not.  If a gentity
 * is not parametric, all of the following functions will return an error
 * when called on that entity.
 * @param gentity_handle Gentity being queried.
 */
int32_t
TSTTG_CGM::CgmGeom_impl::gentityIsParametric (
  /*in*/ void* gentity_handle ) 
throw ( 
  ::TSTTB::Error
){
    // DO-NOT-DELETE splicer.begin(TSTTG_CGM.CgmGeom.gentityIsParametric)
    // insert implementation here
  int32_t retval = TSTTG_gentityIsParametric (tsttgInstance, 
                                              gentity_handle);
  PROCESS_ERROR;
  return retval;
    // DO-NOT-DELETE splicer.end(TSTTG_CGM.CgmGeom.gentityIsParametric)
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
TSTTG_CGM::CgmGeom_impl::gentityUvToXyz (
  /*in*/ ::sidl::array<void*> gentity_handles,
  /*in*/ int32_t gentity_handles_size,
  /*in*/ ::sidl::array<double> uv,
  /*in*/ int32_t uv_size,
  /*inout*/ ::sidl::array<double>& coordinates,
  /*out*/ int32_t& coordinates_size ) 
throw ( 
  ::TSTTB::Error
){
    // DO-NOT-DELETE splicer.begin(TSTTG_CGM.CgmGeom.gentityUvToXyz)
    // insert implementation here
  CREATE_TEMP_ARRAY(double, coordinates);
  
  TSTTG_gentityUvToXyz (tsttgInstance, 
                        (const void **) TEMP_ARRAY_IN(gentity_handles),
                        TEMP_ARRAY_IN(uv),
                        TEMP_ARRAY_INOUT(coordinates));
  PROCESS_ERROR;

  ASSIGN_ARRAY(coordinates);
    // DO-NOT-DELETE splicer.end(TSTTG_CGM.CgmGeom.gentityUvToXyz)
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
TSTTG_CGM::CgmGeom_impl::gentityXyzToUv (
  /*in*/ ::sidl::array<void*> gentity_handles,
  /*in*/ int32_t gentity_handles_size,
  /*in*/ ::sidl::array<double> coordinates,
  /*in*/ int32_t coordinates_size,
  /*inout*/ ::sidl::array<double>& uv,
  /*out*/ int32_t& uv_size ) 
throw ( 
  ::TSTTB::Error
){
    // DO-NOT-DELETE splicer.begin(TSTTG_CGM.CgmGeom.gentityXyzToUv)
    // insert implementation here
  CREATE_TEMP_ARRAY(double, uv);
  
  TSTTG_gentityXyzToUv (tsttgInstance, 
                        (const void **) TEMP_ARRAY_IN(gentity_handles),
                        TEMP_ARRAY_IN(coordinates),
                        TEMP_ARRAY_INOUT(uv));
  PROCESS_ERROR;

  ASSIGN_ARRAY(uv);
    // DO-NOT-DELETE splicer.end(TSTTG_CGM.CgmGeom.gentityXyzToUv)
}

/**
 * Return the uv range of the specified gentities.  Parameters are interleaved.
 * @param gentity_handles Gentities being queried.
 * @param uv_min Minimum parameters of gentities, interleaved
 * @param uv_max Maximum parameters of gentities, interleaved
 */
void
TSTTG_CGM::CgmGeom_impl::gentityUvRange (
  /*in*/ ::sidl::array<void*> gentity_handles,
  /*in*/ int32_t gentity_handles_size,
  /*inout*/ ::sidl::array<double>& uv_min,
  /*out*/ int32_t& uv_min_size,
  /*inout*/ ::sidl::array<double>& uv_max,
  /*out*/ int32_t& uv_max_size ) 
throw ( 
  ::TSTTB::Error
){
    // DO-NOT-DELETE splicer.begin(TSTTG_CGM.CgmGeom.gentityUvRange)
    // insert implementation here
  CREATE_TEMP_ARRAY(double, uv_max);
  CREATE_TEMP_ARRAY(double, uv_min);
  
  TSTTG_gentityUvRange (tsttgInstance, 
                        (const void **) TEMP_ARRAY_IN(gentity_handles),
                        TEMP_ARRAY_INOUT(uv_min),
                        TEMP_ARRAY_INOUT(uv_max));
  PROCESS_ERROR;

  ASSIGN_ARRAY(uv_max);
  ASSIGN_ARRAY(uv_min);
    // DO-NOT-DELETE splicer.end(TSTTG_CGM.CgmGeom.gentityUvRange)
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
void
TSTTG_CGM::CgmGeom_impl::Greparam_edge_face (
  /*in*/ ::sidl::array<void*> src_gentity_handles,
  /*in*/ int32_t src_gentity_handles_size,
  /*in*/ ::sidl::array<double> src_uv,
  /*in*/ int32_t src_uv_size,
  /*in*/ ::sidl::array<void*> trg_gentity_handles,
  /*in*/ int32_t trg_gentity_handles_size,
  /*in*/ ::sidl::array<double> trg_uv,
  /*in*/ int32_t trg_uv_size ) 
throw ( 
  ::TSTTB::Error
){
    // DO-NOT-DELETE splicer.begin(TSTTG_CGM.CgmGeom.Greparam_edge_face)
    // insert implementation here
  TSTTG_Greparam_edge_face (tsttgInstance, 
                            TEMP_ARRAY_IN(src_gentity_handles),
                            TEMP_ARRAY_IN(src_uv),
                            TEMP_ARRAY_IN(trg_gentity_handles),
                            TEMP_ARRAY_IN(trg_uv));
  PROCESS_ERROR;
    // DO-NOT-DELETE splicer.end(TSTTG_CGM.CgmGeom.Greparam_edge_face)
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
TSTTG_CGM::CgmGeom_impl::gentityNormalUv (
  /*in*/ ::sidl::array<void*> gface_handles,
  /*in*/ int32_t gface_handles_size,
  /*in*/ ::sidl::array<double> parameters,
  /*in*/ int32_t parameters_size,
  /*inout*/ ::sidl::array<double>& normals,
  /*in*/ int32_t normals_size ) 
throw ( 
  ::TSTTB::Error
){
    // DO-NOT-DELETE splicer.begin(TSTTG_CGM.CgmGeom.gentityNormalUv)
    // insert implementation here
  CREATE_TEMP_ARRAY(double, normals);
  
  TSTTG_gentityNormalUv (tsttgInstance, 
                         (const void **) TEMP_ARRAY_IN(gface_handles),
                         TEMP_ARRAY_IN(parameters),
                         TEMP_ARRAY_INOUT(normals));
  PROCESS_ERROR;

  ASSIGN_ARRAY(normals);
    // DO-NOT-DELETE splicer.end(TSTTG_CGM.CgmGeom.gentityNormalUv)
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
TSTTG_CGM::CgmGeom_impl::gentityTangentU (
  /*in*/ ::sidl::array<void*> gedge_handles,
  /*in*/ int32_t gedge_handles_size,
  /*in*/ ::sidl::array<double> parameters,
  /*in*/ int32_t parameters_size,
  /*inout*/ ::sidl::array<double>& tangents,
  /*out*/ int32_t& tangents_size ) 
throw ( 
  ::TSTTB::Error
){
    // DO-NOT-DELETE splicer.begin(TSTTG_CGM.CgmGeom.gentityTangentU)
    // insert implementation here
  CREATE_TEMP_ARRAY(double, tangents);
  
  TSTTG_gentityTangentU (tsttgInstance, 
                         (const void **) TEMP_ARRAY_IN(gedge_handles),
                         TEMP_ARRAY_IN(parameters),
                         TEMP_ARRAY_INOUT(tangents));
  PROCESS_ERROR;

  ASSIGN_ARRAY(tangents);
    // DO-NOT-DELETE splicer.end(TSTTG_CGM.CgmGeom.gentityTangentU)
}

/**
 * Method:  getData[]
 */
void
TSTTG_CGM::CgmGeom_impl::getData (
  /*in*/ void* entity_handle,
  /*in*/ void* tag_handle,
  /*inout*/ ::sidl::array<char>& tag_value,
  /*out*/ int32_t& tag_value_size ) 
throw ( 
  ::TSTTB::Error
){
    // DO-NOT-DELETE splicer.begin(TSTTG_CGM.CgmGeom.getData)
    // insert implementation here
  CREATE_TEMP_TAG_ARRAY(tag_value);
  
  TSTTG_getData (tsttgInstance,
                 entity_handle,
                 tag_handle,
                 TEMP_ARRAY_INOUT(tag_value));
  PROCESS_ERROR;

  ASSIGN_TAG_ARRAY(tag_value);
    // DO-NOT-DELETE splicer.end(TSTTG_CGM.CgmGeom.getData)
}

/**
 * Method:  getIntData[]
 */
int32_t
TSTTG_CGM::CgmGeom_impl::getIntData (
  /*in*/ void* entity_handle,
  /*in*/ void* tag_handle ) 
throw ( 
  ::TSTTB::Error
){
    // DO-NOT-DELETE splicer.begin(TSTTG_CGM.CgmGeom.getIntData)
    // insert implementation here
  int32_t retval = TSTTG_getIntData(tsttgInstance,
                                    entity_handle,
                                    tag_handle );
  PROCESS_ERROR;
  return retval;
    // DO-NOT-DELETE splicer.end(TSTTG_CGM.CgmGeom.getIntData)
}

/**
 * Method:  getDblData[]
 */
double
TSTTG_CGM::CgmGeom_impl::getDblData (
  /*in*/ void* entity_handle,
  /*in*/ void* tag_handle ) 
throw ( 
  ::TSTTB::Error
){
    // DO-NOT-DELETE splicer.begin(TSTTG_CGM.CgmGeom.getDblData)
    // insert implementation here
  double retval = TSTTG_getDblData(tsttgInstance,
                                   entity_handle,
                                   tag_handle );  
  PROCESS_ERROR;
  return retval;
    // DO-NOT-DELETE splicer.end(TSTTG_CGM.CgmGeom.getDblData)
}

/**
 * Method:  getEHData[]
 */
void*
TSTTG_CGM::CgmGeom_impl::getEHData (
  /*in*/ void* entity_handle,
  /*in*/ void* tag_handle ) 
throw ( 
  ::TSTTB::Error
){
    // DO-NOT-DELETE splicer.begin(TSTTG_CGM.CgmGeom.getEHData)
    // insert implementation here
  void *retval = TSTTG_getEHData(tsttgInstance,
                                 entity_handle,
                                 tag_handle );
  PROCESS_ERROR;
  return retval;
    // DO-NOT-DELETE splicer.end(TSTTG_CGM.CgmGeom.getEHData)
}

/**
 * Method:  setData[]
 */
void
TSTTG_CGM::CgmGeom_impl::setData (
  /*in*/ void* entity_handle,
  /*in*/ void* tag_handle,
  /*in*/ ::sidl::array<char> tag_value,
  /*in*/ int32_t tag_value_size ) 
throw ( 
  ::TSTTB::Error
){
    // DO-NOT-DELETE splicer.begin(TSTTG_CGM.CgmGeom.setData)
    // insert implementation here
  TSTTG_setData(tsttgInstance, 
                entity_handle,
                tag_handle, 
                TEMP_TAG_ARRAY_IN(tag_value));
  PROCESS_ERROR;
    // DO-NOT-DELETE splicer.end(TSTTG_CGM.CgmGeom.setData)
}

/**
 * Method:  setIntData[]
 */
void
TSTTG_CGM::CgmGeom_impl::setIntData (
  /*in*/ void* entity_handle,
  /*in*/ void* tag_handle,
  /*in*/ int32_t tag_value ) 
throw ( 
  ::TSTTB::Error
){
    // DO-NOT-DELETE splicer.begin(TSTTG_CGM.CgmGeom.setIntData)
    // insert implementation here
  TSTTG_setIntData (tsttgInstance,
                    entity_handle,
                    tag_handle,
                    tag_value );
  PROCESS_ERROR;
    // DO-NOT-DELETE splicer.end(TSTTG_CGM.CgmGeom.setIntData)
}

/**
 * Method:  setDblData[]
 */
void
TSTTG_CGM::CgmGeom_impl::setDblData (
  /*in*/ void* entity_handle,
  /*in*/ void* tag_handle,
  /*in*/ double tag_value ) 
throw ( 
  ::TSTTB::Error
){
    // DO-NOT-DELETE splicer.begin(TSTTG_CGM.CgmGeom.setDblData)
    // insert implementation here
  TSTTG_setDblData (tsttgInstance,
                    entity_handle,
                    tag_handle,
                    tag_value );
  PROCESS_ERROR;
    // DO-NOT-DELETE splicer.end(TSTTG_CGM.CgmGeom.setDblData)
}

/**
 * Method:  setEHData[]
 */
void
TSTTG_CGM::CgmGeom_impl::setEHData (
  /*in*/ void* entity_handle,
  /*in*/ void* tag_handle,
  /*in*/ void* tag_value ) 
throw ( 
  ::TSTTB::Error
){
    // DO-NOT-DELETE splicer.begin(TSTTG_CGM.CgmGeom.setEHData)
    // insert implementation here
  TSTTG_setEHData (tsttgInstance,
                   entity_handle,
                   tag_handle,
                   tag_value );
  PROCESS_ERROR;
    // DO-NOT-DELETE splicer.end(TSTTG_CGM.CgmGeom.setEHData)
}

/**
 * Method:  getAllTags[]
 */
void
TSTTG_CGM::CgmGeom_impl::getAllTags (
  /*in*/ void* entity_handle,
  /*inout*/ ::sidl::array<void*>& tag_handles,
  /*out*/ int32_t& tag_handles_size ) 
throw ( 
  ::TSTTB::Error
){
    // DO-NOT-DELETE splicer.begin(TSTTG_CGM.CgmGeom.getAllTags)
    // insert implementation here
  CREATE_TEMP_ARRAY(void*, tag_handles);
  
  TSTTG_getAllTags(tsttgInstance, entity_handle, 
                   TEMP_ARRAY_INOUT(tag_handles));
  PROCESS_ERROR;

  ASSIGN_ARRAY(tag_handles);
    // DO-NOT-DELETE splicer.end(TSTTG_CGM.CgmGeom.getAllTags)
}

/**
 * Method:  rmvTag[]
 */
void
TSTTG_CGM::CgmGeom_impl::rmvTag (
  /*in*/ void* entity_handle,
  /*in*/ void* tag_handle ) 
throw ( 
  ::TSTTB::Error
){
    // DO-NOT-DELETE splicer.begin(TSTTG_CGM.CgmGeom.rmvTag)
    // insert implementation here
  TSTTG_rmvTag (tsttgInstance,
                entity_handle,
                tag_handle );
  PROCESS_ERROR;
    // DO-NOT-DELETE splicer.end(TSTTG_CGM.CgmGeom.rmvTag)
}

/**
 * Method:  createEntSet[]
 */
void
TSTTG_CGM::CgmGeom_impl::createEntSet (
  /*in*/ bool isList,
  /*out*/ void*& entity_set ) 
throw ( 
  ::TSTTB::Error
){
    // DO-NOT-DELETE splicer.begin(TSTTG_CGM.CgmGeom.createEntSet)
    // insert implementation here
  TSTTG_createEntSet(tsttgInstance, isList, &entity_set);
  PROCESS_ERROR;
    // DO-NOT-DELETE splicer.end(TSTTG_CGM.CgmGeom.createEntSet)
}

/**
 * Method:  destroyEntSet[]
 */
void
TSTTG_CGM::CgmGeom_impl::destroyEntSet (
  /*in*/ void* entity_set ) 
throw ( 
  ::TSTTB::Error
){
    // DO-NOT-DELETE splicer.begin(TSTTG_CGM.CgmGeom.destroyEntSet)
    // insert implementation here
  TSTTG_destroyEntSet(tsttgInstance, entity_set);
  PROCESS_ERROR;
    // DO-NOT-DELETE splicer.end(TSTTG_CGM.CgmGeom.destroyEntSet)
}

/**
 * Method:  isList[]
 */
bool
TSTTG_CGM::CgmGeom_impl::isList (
  /*in*/ void* entity_set ) 
throw ( 
  ::TSTTB::Error
){
    // DO-NOT-DELETE splicer.begin(TSTTG_CGM.CgmGeom.isList)
    // insert implementation here
  bool retval = TSTTG_isList(tsttgInstance, entity_set);
  PROCESS_ERROR;
  return retval;
    // DO-NOT-DELETE splicer.end(TSTTG_CGM.CgmGeom.isList)
}

/**
 * Method:  getNumEntSets[]
 */
int32_t
TSTTG_CGM::CgmGeom_impl::getNumEntSets (
  /*in*/ void* entity_set,
  /*in*/ int32_t num_hops ) 
throw ( 
  ::TSTTB::Error
){
    // DO-NOT-DELETE splicer.begin(TSTTG_CGM.CgmGeom.getNumEntSets)
    // insert implementation here
  int32_t retval = TSTTG_getNumEntSets(tsttgInstance, entity_set, num_hops);
  PROCESS_ERROR;
  return retval;
    // DO-NOT-DELETE splicer.end(TSTTG_CGM.CgmGeom.getNumEntSets)
}

/**
 * Method:  getEntSets[]
 */
void
TSTTG_CGM::CgmGeom_impl::getEntSets (
  /*in*/ void* entity_set,
  /*in*/ int32_t num_hops,
  /*inout*/ ::sidl::array<void*>& contained_entset_handles,
  /*out*/ int32_t& contained_entset_handles_size ) 
throw ( 
  ::TSTTB::Error
){
    // DO-NOT-DELETE splicer.begin(TSTTG_CGM.CgmGeom.getEntSets)
    // insert implementation here
  CREATE_TEMP_ARRAY(void*, contained_entset_handles);

  TSTTG_getEntSets(tsttgInstance, entity_set, num_hops, 
                   TEMP_ARRAY_INOUT(contained_entset_handles));
  PROCESS_ERROR;
  ASSIGN_ARRAY(contained_entset_handles);
    // DO-NOT-DELETE splicer.end(TSTTG_CGM.CgmGeom.getEntSets)
}

/**
 * Method:  addEntToSet[]
 */
void
TSTTG_CGM::CgmGeom_impl::addEntToSet (
  /*in*/ void* entity_handle,
  /*inout*/ void*& entity_set ) 
throw ( 
  ::TSTTB::Error
){
    // DO-NOT-DELETE splicer.begin(TSTTG_CGM.CgmGeom.addEntToSet)
    // insert implementation here
  TSTTG_addEntToSet (tsttgInstance,
                     entity_handle,
                     &entity_set);
  PROCESS_ERROR;
    // DO-NOT-DELETE splicer.end(TSTTG_CGM.CgmGeom.addEntToSet)
}

/**
 * Method:  rmvEntFromSet[]
 */
void
TSTTG_CGM::CgmGeom_impl::rmvEntFromSet (
  /*in*/ void* entity_handle,
  /*inout*/ void*& entity_set ) 
throw ( 
  ::TSTTB::Error
){
    // DO-NOT-DELETE splicer.begin(TSTTG_CGM.CgmGeom.rmvEntFromSet)
    // insert implementation here
  TSTTG_rmvEntFromSet (tsttgInstance,
                       entity_handle,
                       &entity_set );
  PROCESS_ERROR;
    // DO-NOT-DELETE splicer.end(TSTTG_CGM.CgmGeom.rmvEntFromSet)
}

/**
 * Method:  addEntArrToSet[]
 */
void
TSTTG_CGM::CgmGeom_impl::addEntArrToSet (
  /*in*/ ::sidl::array<void*> entity_handles,
  /*in*/ int32_t entity_handles_size,
  /*inout*/ void*& entity_set ) 
throw ( 
  ::TSTTB::Error
){
    // DO-NOT-DELETE splicer.begin(TSTTG_CGM.CgmGeom.addEntArrToSet)
    // insert implementation here
  entity_handles.addRef();
  void **tmp_ptr1 = entity_handles._get_ior()->d_firstElement;
    //void **tmp_ptr3 = entity_handles._get_ior()->d_firstElement;
  TSTTG_addEntArrToSet(tsttgInstance, tmp_ptr1,
                       entity_handles_size,
                       &entity_set);
  PROCESS_ERROR;
    // DO-NOT-DELETE splicer.end(TSTTG_CGM.CgmGeom.addEntArrToSet)
}

/**
 * Method:  rmvEntArrFromSet[]
 */
void
TSTTG_CGM::CgmGeom_impl::rmvEntArrFromSet (
  /*in*/ ::sidl::array<void*> entity_handles,
  /*in*/ int32_t entity_handles_size,
  /*inout*/ void*& entity_set ) 
throw ( 
  ::TSTTB::Error
){
    // DO-NOT-DELETE splicer.begin(TSTTG_CGM.CgmGeom.rmvEntArrFromSet)
    // insert implementation here
  TSTTG_rmvEntArrFromSet(tsttgInstance, TEMP_ARRAY_IN(entity_handles), 
                         &entity_set);
  PROCESS_ERROR;
    // DO-NOT-DELETE splicer.end(TSTTG_CGM.CgmGeom.rmvEntArrFromSet)
}

/**
 * Method:  addEntSet[]
 */
void
TSTTG_CGM::CgmGeom_impl::addEntSet (
  /*in*/ void* entity_set_to_add,
  /*inout*/ void*& entity_set_handle ) 
throw ( 
  ::TSTTB::Error
){
    // DO-NOT-DELETE splicer.begin(TSTTG_CGM.CgmGeom.addEntSet)
    // insert implementation here
  TSTTG_addEntSet(tsttgInstance, entity_set_to_add, 
                  &entity_set_handle);
  PROCESS_ERROR;
    // DO-NOT-DELETE splicer.end(TSTTG_CGM.CgmGeom.addEntSet)
}

/**
 * Method:  rmvEntSet[]
 */
void
TSTTG_CGM::CgmGeom_impl::rmvEntSet (
  /*in*/ void* entity_set_to_remove,
  /*inout*/ void*& entity_set_handle ) 
throw ( 
  ::TSTTB::Error
){
    // DO-NOT-DELETE splicer.begin(TSTTG_CGM.CgmGeom.rmvEntSet)
    // insert implementation here
  TSTTG_rmvEntSet (tsttgInstance,
                   entity_set_to_remove,
                   &entity_set_handle );
  PROCESS_ERROR;
    // DO-NOT-DELETE splicer.end(TSTTG_CGM.CgmGeom.rmvEntSet)
}

/**
 * Method:  isEntContained[]
 */
bool
TSTTG_CGM::CgmGeom_impl::isEntContained (
  /*in*/ void* containing_entity_set,
  /*in*/ void* entity_handle ) 
throw ( 
  ::TSTTB::Error
){
  // DO-NOT-DELETE splicer.begin(TSTTG_CGM.CgmGeom.isEntContained)
  // insert implementation here
  bool retval = TSTTG_isEntContained(tsttgInstance,
                                     containing_entity_set,
                                     entity_handle);
  PROCESS_ERROR;
  return retval;
  // DO-NOT-DELETE splicer.end(TSTTG_CGM.CgmGeom.isEntContained)
}

/**
 * Method:  isEntSetContained[]
 */
bool
TSTTG_CGM::CgmGeom_impl::isEntSetContained (
  /*in*/ void* containing_entity_set,
  /*in*/ void* contained_entity_set ) 
throw ( 
  ::TSTTB::Error
){
  // DO-NOT-DELETE splicer.begin(TSTTG_CGM.CgmGeom.isEntSetContained)
  // insert implementation here
  bool retval = TSTTG_isEntSetContained(tsttgInstance,
                                        containing_entity_set,
                                        contained_entity_set);
  PROCESS_ERROR;
  return retval;
  // DO-NOT-DELETE splicer.end(TSTTG_CGM.CgmGeom.isEntSetContained)
}

/**
 * Method:  subtract[]
 */
void
TSTTG_CGM::CgmGeom_impl::subtract (
  /*in*/ void* entity_set_1,
  /*in*/ void* entity_set_2,
  /*out*/ void*& result_entity_set ) 
throw ( 
  ::TSTTB::Error
){
    // DO-NOT-DELETE splicer.begin(TSTTG_CGM.CgmGeom.subtract)
    // insert implementation here
  TSTTG_subtract(tsttgInstance, entity_set_1, entity_set_2, &result_entity_set);
  PROCESS_ERROR;
    // DO-NOT-DELETE splicer.end(TSTTG_CGM.CgmGeom.subtract)
}

/**
 * Method:  intersect[]
 */
void
TSTTG_CGM::CgmGeom_impl::intersect (
  /*in*/ void* entity_set_1,
  /*in*/ void* entity_set_2,
  /*out*/ void*& result_entity_set ) 
throw ( 
  ::TSTTB::Error
){
    // DO-NOT-DELETE splicer.begin(TSTTG_CGM.CgmGeom.intersect)
    // insert implementation here
  TSTTG_intersect(tsttgInstance, entity_set_1, entity_set_2, &result_entity_set);
  PROCESS_ERROR;
    // DO-NOT-DELETE splicer.end(TSTTG_CGM.CgmGeom.intersect)
}

/**
 * Method:  unite[]
 */
void
TSTTG_CGM::CgmGeom_impl::unite (
  /*in*/ void* entity_set_1,
  /*in*/ void* entity_set_2,
  /*out*/ void*& result_entity_set ) 
throw ( 
  ::TSTTB::Error
){
    // DO-NOT-DELETE splicer.begin(TSTTG_CGM.CgmGeom.unite)
    // insert implementation here
  TSTTG_unite(tsttgInstance, entity_set_1, entity_set_2, &result_entity_set);
  PROCESS_ERROR;
    // DO-NOT-DELETE splicer.end(TSTTG_CGM.CgmGeom.unite)
}

/**
 * Return the relative and absolute tolerances at the modeler level.  If
 * model does not have a modeler-wide tolerance, zero is returned for both
 * values.
 * @param relative_tolerance Relative tolerance for model as a whole
 * @param absolute_tolerance Absolute tolerance for model as a whole
 */
void
TSTTG_CGM::CgmGeom_impl::getGtolerance (
  /*out*/ double& relative_tolerance,
  /*out*/ double& absolute_tolerance ) 
throw ( 
  ::TSTTB::Error
){
    // DO-NOT-DELETE splicer.begin(TSTTG_CGM.CgmGeom.getGtolerance)
    // insert implementation here
  TSTTG_getGtolerance (tsttgInstance, 
                       &relative_tolerance,
                       &absolute_tolerance);
  PROCESS_ERROR;
    // DO-NOT-DELETE splicer.end(TSTTG_CGM.CgmGeom.getGtolerance)
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
TSTTG_CGM::CgmGeom_impl::getGentityTolerance (
  /*in*/ ::sidl::array<void*> gentity_handles,
  /*in*/ int32_t gentity_handles_size,
  /*inout*/ ::sidl::array<double>& relative_tolerances,
  /*out*/ int32_t& relative_tolerances_size,
  /*inout*/ ::sidl::array<double>& absolute_tolerances,
  /*out*/ int32_t& absolute_tolerances_size ) 
throw ( 
  ::TSTTB::Error
){
    // DO-NOT-DELETE splicer.begin(TSTTG_CGM.CgmGeom.getGentityTolerance)
    // insert implementation here
  CREATE_TEMP_ARRAY(double, absolute_tolerances);
  CREATE_TEMP_ARRAY(double, relative_tolerances);
  
  TSTTG_getGentityTolerance (tsttgInstance, 
                             (const void **) TEMP_ARRAY_IN(gentity_handles),
                             TEMP_ARRAY_INOUT(relative_tolerances),
                             TEMP_ARRAY_INOUT(absolute_tolerances));
  PROCESS_ERROR;

  ASSIGN_ARRAY(absolute_tolerances);
  ASSIGN_ARRAY(relative_tolerances);
    // DO-NOT-DELETE splicer.end(TSTTG_CGM.CgmGeom.getGentityTolerance)
}

/**
 * Method:  setEntSetData[]
 */
void
TSTTG_CGM::CgmGeom_impl::setEntSetData (
  /*in*/ void* entity_set,
  /*in*/ void* tag_handle,
  /*inout*/ ::sidl::array<char>& tag_value,
  /*in*/ int32_t tag_value_size ) 
throw ( 
  ::TSTTB::Error
){
    // DO-NOT-DELETE splicer.begin(TSTTG_CGM.CgmGeom.setEntSetData)
    // insert implementation here
  TSTTG_setEntSetData(tsttgInstance, entity_set, tag_handle,
                      TEMP_TAG_ARRAY_IN(tag_value));
  PROCESS_ERROR;
    // DO-NOT-DELETE splicer.end(TSTTG_CGM.CgmGeom.setEntSetData)
}

/**
 * Method:  setEntSetIntData[]
 */
void
TSTTG_CGM::CgmGeom_impl::setEntSetIntData (
  /*in*/ void* entity_set,
  /*in*/ void* tag_handle,
  /*in*/ int32_t tag_value ) 
throw ( 
  ::TSTTB::Error
){
    // DO-NOT-DELETE splicer.begin(TSTTG_CGM.CgmGeom.setEntSetIntData)
    // insert implementation here
  TSTTG_setEntSetIntData (tsttgInstance,
                          entity_set,
                          tag_handle,
                          tag_value );
  PROCESS_ERROR;
    // DO-NOT-DELETE splicer.end(TSTTG_CGM.CgmGeom.setEntSetIntData)
}

/**
 * Method:  setEntSetDblData[]
 */
void
TSTTG_CGM::CgmGeom_impl::setEntSetDblData (
  /*in*/ void* entity_set,
  /*in*/ void* tag_handle,
  /*in*/ double tag_value ) 
throw ( 
  ::TSTTB::Error
){
    // DO-NOT-DELETE splicer.begin(TSTTG_CGM.CgmGeom.setEntSetDblData)
    // insert implementation here
  TSTTG_setEntSetDblData (tsttgInstance,
                          entity_set,
                          tag_handle,
                          tag_value );
  PROCESS_ERROR;
    // DO-NOT-DELETE splicer.end(TSTTG_CGM.CgmGeom.setEntSetDblData)
}

/**
 * Method:  setEntSetEHData[]
 */
void
TSTTG_CGM::CgmGeom_impl::setEntSetEHData (
  /*in*/ void* entity_set,
  /*in*/ void* tag_handle,
  /*in*/ void* tag_value ) 
throw ( 
  ::TSTTB::Error
){
    // DO-NOT-DELETE splicer.begin(TSTTG_CGM.CgmGeom.setEntSetEHData)
    // insert implementation here
  TSTTG_setEntSetEHData (tsttgInstance,
                         entity_set,
                         tag_handle,
                         tag_value );
  PROCESS_ERROR;
    // DO-NOT-DELETE splicer.end(TSTTG_CGM.CgmGeom.setEntSetEHData)
}

/**
 * Method:  getEntSetData[]
 */
void
TSTTG_CGM::CgmGeom_impl::getEntSetData (
  /*in*/ void* entity_set,
  /*in*/ void* tag_handle,
  /*inout*/ ::sidl::array<char>& tag_value,
  /*out*/ int32_t& tag_value_size ) 
throw ( 
  ::TSTTB::Error
){
    // DO-NOT-DELETE splicer.begin(TSTTG_CGM.CgmGeom.getEntSetData)
    // insert implementation here
  CREATE_TEMP_TAG_ARRAY(tag_value);
  
  TSTTG_getEntSetData (tsttgInstance,
                       entity_set,
                       tag_handle,
                       TEMP_ARRAY_INOUT(tag_value));
  PROCESS_ERROR;

  ASSIGN_TAG_ARRAY(tag_value);
    // DO-NOT-DELETE splicer.end(TSTTG_CGM.CgmGeom.getEntSetData)
}

/**
 * Method:  getEntSetIntData[]
 */
int32_t
TSTTG_CGM::CgmGeom_impl::getEntSetIntData (
  /*in*/ void* entity_set,
  /*in*/ void* tag_handle ) 
throw ( 
  ::TSTTB::Error
){
    // DO-NOT-DELETE splicer.begin(TSTTG_CGM.CgmGeom.getEntSetIntData)
    // insert implementation here
  int32_t retval = TSTTG_getEntSetIntData(tsttgInstance,
                                          entity_set,
                                          tag_handle );
  PROCESS_ERROR;
  return retval;
    // DO-NOT-DELETE splicer.end(TSTTG_CGM.CgmGeom.getEntSetIntData)
}

/**
 * Method:  getEntSetDblData[]
 */
double
TSTTG_CGM::CgmGeom_impl::getEntSetDblData (
  /*in*/ void* entity_set,
  /*in*/ void* tag_handle ) 
throw ( 
  ::TSTTB::Error
){
    // DO-NOT-DELETE splicer.begin(TSTTG_CGM.CgmGeom.getEntSetDblData)
    // insert implementation here
  double retval = TSTTG_getEntSetDblData(tsttgInstance,
                                         entity_set,
                                         tag_handle );
  PROCESS_ERROR;
  return retval;
    // DO-NOT-DELETE splicer.end(TSTTG_CGM.CgmGeom.getEntSetDblData)
}

/**
 * Method:  getEntSetEHData[]
 */
void*
TSTTG_CGM::CgmGeom_impl::getEntSetEHData (
  /*in*/ void* entity_set,
  /*in*/ void* tag_handle ) 
throw ( 
  ::TSTTB::Error
){
    // DO-NOT-DELETE splicer.begin(TSTTG_CGM.CgmGeom.getEntSetEHData)
    // insert implementation here
  void *retval = TSTTG_getEntSetEHData(tsttgInstance,
                                       entity_set,
                                       tag_handle );
  PROCESS_ERROR;
  return retval;
    // DO-NOT-DELETE splicer.end(TSTTG_CGM.CgmGeom.getEntSetEHData)
}

/**
 * Method:  getAllEntSetTags[]
 */
void
TSTTG_CGM::CgmGeom_impl::getAllEntSetTags (
  /*in*/ void* entity_set,
  /*inout*/ ::sidl::array<void*>& tag_handles,
  /*out*/ int32_t& tag_handles_size ) 
throw ( 
  ::TSTTB::Error
){
    // DO-NOT-DELETE splicer.begin(TSTTG_CGM.CgmGeom.getAllEntSetTags)
    // insert implementation here
  CREATE_TEMP_ARRAY(void*, tag_handles);
  
  TSTTG_getAllEntSetTags (tsttgInstance,
                          entity_set,
                          TEMP_ARRAY_INOUT(tag_handles));
  PROCESS_ERROR;

  ASSIGN_ARRAY(tag_handles);
    // DO-NOT-DELETE splicer.end(TSTTG_CGM.CgmGeom.getAllEntSetTags)
}

/**
 * Method:  rmvEntSetTag[]
 */
void
TSTTG_CGM::CgmGeom_impl::rmvEntSetTag (
  /*in*/ void* entity_set,
  /*in*/ void* tag_handle ) 
throw ( 
  ::TSTTB::Error
){
    // DO-NOT-DELETE splicer.begin(TSTTG_CGM.CgmGeom.rmvEntSetTag)
    // insert implementation here
  TSTTG_rmvEntSetTag (tsttgInstance,
                      entity_set,
                      tag_handle );
  PROCESS_ERROR;
    // DO-NOT-DELETE splicer.end(TSTTG_CGM.CgmGeom.rmvEntSetTag)
}

/**
 * Method:  Brick[]
 */
void
TSTTG_CGM::CgmGeom_impl::Brick (
  /*in*/ double x,
  /*in*/ double y,
  /*in*/ double z,
  /*out*/ void*& geom_entity ) 
throw ( 
  ::TSTTB::Error
){
    // DO-NOT-DELETE splicer.begin(TSTTG_CGM.CgmGeom.Brick)
    // insert implementation here
  TSTTG_Brick(tsttgInstance, x, y, z, &geom_entity);
  PROCESS_ERROR;
    // DO-NOT-DELETE splicer.end(TSTTG_CGM.CgmGeom.Brick)
}

/**
 * Method:  Cylinder[]
 */
void
TSTTG_CGM::CgmGeom_impl::Cylinder (
  /*in*/ double height,
  /*in*/ double major_rad,
  /*in*/ double minor_rad,
  /*out*/ void*& geom_entity ) 
throw ( 
  ::TSTTB::Error
){
    // DO-NOT-DELETE splicer.begin(TSTTG_CGM.CgmGeom.Cylinder)
    // insert implementation here
  TSTTG_Cylinder(tsttgInstance, height, major_rad, minor_rad, &geom_entity);
  PROCESS_ERROR;
    // DO-NOT-DELETE splicer.end(TSTTG_CGM.CgmGeom.Cylinder)
}

/**
 * Method:  Torus[]
 */
void
TSTTG_CGM::CgmGeom_impl::Torus (
  /*in*/ double major_rad,
  /*in*/ double minor_rad,
  /*out*/ void*& geom_entity ) 
throw ( 
  ::TSTTB::Error
){
    // DO-NOT-DELETE splicer.begin(TSTTG_CGM.CgmGeom.Torus)
    // insert implementation here
  TSTTG_Torus(tsttgInstance, major_rad, minor_rad, &geom_entity);
  PROCESS_ERROR;
    // DO-NOT-DELETE splicer.end(TSTTG_CGM.CgmGeom.Torus)
}

/**
 * Method:  addPrntChld[]
 */
void
TSTTG_CGM::CgmGeom_impl::addPrntChld (
  /*inout*/ void*& parent_entity_set,
  /*inout*/ void*& child_entity_set ) 
throw ( 
  ::TSTTB::Error
){
    // DO-NOT-DELETE splicer.begin(TSTTG_CGM.CgmGeom.addPrntChld)
    // insert implementation here
  TSTTG_addPrntChld(tsttgInstance, &parent_entity_set, &child_entity_set);
  PROCESS_ERROR;
    // DO-NOT-DELETE splicer.end(TSTTG_CGM.CgmGeom.addPrntChld)
}

/**
 * Method:  rmvPrntChld[]
 */
void
TSTTG_CGM::CgmGeom_impl::rmvPrntChld (
  /*inout*/ void*& parent_entity_set,
  /*inout*/ void*& child_entity_set ) 
throw ( 
  ::TSTTB::Error
){
    // DO-NOT-DELETE splicer.begin(TSTTG_CGM.CgmGeom.rmvPrntChld)
    // insert implementation here
  TSTTG_rmvPrntChld(tsttgInstance, &parent_entity_set, &child_entity_set);
  PROCESS_ERROR;
    // DO-NOT-DELETE splicer.end(TSTTG_CGM.CgmGeom.rmvPrntChld)
}

/**
 * Method:  isChildOf[]
 */
bool
TSTTG_CGM::CgmGeom_impl::isChildOf (
  /*in*/ void* parent_entity_set,
  /*in*/ void* child_entity_set ) 
throw ( 
  ::TSTTB::Error
){
    // DO-NOT-DELETE splicer.begin(TSTTG_CGM.CgmGeom.isChildOf)
    // insert implementation here
  bool retval = TSTTG_isChildOf(tsttgInstance, parent_entity_set, 
                                child_entity_set);
  PROCESS_ERROR;
  return retval;
    // DO-NOT-DELETE splicer.end(TSTTG_CGM.CgmGeom.isChildOf)
}

/**
 * Method:  getNumChld[]
 */
int32_t
TSTTG_CGM::CgmGeom_impl::getNumChld (
  /*in*/ void* entity_set,
  /*in*/ int32_t num_hops ) 
throw ( 
  ::TSTTB::Error
){
    // DO-NOT-DELETE splicer.begin(TSTTG_CGM.CgmGeom.getNumChld)
    // insert implementation here
  int32_t retval = TSTTG_getNumChld(tsttgInstance, entity_set, num_hops);
  PROCESS_ERROR;
  return retval;
    // DO-NOT-DELETE splicer.end(TSTTG_CGM.CgmGeom.getNumChld)
}

/**
 * Method:  getNumPrnt[]
 */
int32_t
TSTTG_CGM::CgmGeom_impl::getNumPrnt (
  /*in*/ void* entity_set,
  /*in*/ int32_t num_hops ) 
throw ( 
  ::TSTTB::Error
){
    // DO-NOT-DELETE splicer.begin(TSTTG_CGM.CgmGeom.getNumPrnt)
    // insert implementation here
  int32_t retval = TSTTG_getNumPrnt(tsttgInstance, entity_set, num_hops);
  PROCESS_ERROR;
  return retval;
    // DO-NOT-DELETE splicer.end(TSTTG_CGM.CgmGeom.getNumPrnt)
}

/**
 * Method:  getChldn[]
 */
void
TSTTG_CGM::CgmGeom_impl::getChldn (
  /*in*/ void* from_entity_set,
  /*in*/ int32_t num_hops,
  /*inout*/ ::sidl::array<void*>& child_handles,
  /*out*/ int32_t& child_handles_size ) 
throw ( 
  ::TSTTB::Error
){
    // DO-NOT-DELETE splicer.begin(TSTTG_CGM.CgmGeom.getChldn)
    // insert implementation here
  CREATE_TEMP_ARRAY(void*, child_handles);
  
  TSTTG_getChldn(tsttgInstance, from_entity_set, num_hops, 
                 TEMP_ARRAY_INOUT(child_handles));
  PROCESS_ERROR;

  ASSIGN_ARRAY(child_handles);
    // DO-NOT-DELETE splicer.end(TSTTG_CGM.CgmGeom.getChldn)
}

/**
 * Method:  getPrnts[]
 */
void
TSTTG_CGM::CgmGeom_impl::getPrnts (
  /*in*/ void* from_entity_set,
  /*in*/ int32_t num_hops,
  /*inout*/ ::sidl::array<void*>& parent_handles,
  /*out*/ int32_t& parent_handles_size ) 
throw ( 
  ::TSTTB::Error
){
    // DO-NOT-DELETE splicer.begin(TSTTG_CGM.CgmGeom.getPrnts)
    // insert implementation here
  CREATE_TEMP_ARRAY(void*, parent_handles);
  
  TSTTG_getPrnts(tsttgInstance, from_entity_set, num_hops, 
                 TEMP_ARRAY_INOUT(parent_handles));
  PROCESS_ERROR;
  
  ASSIGN_ARRAY(parent_handles);
    // DO-NOT-DELETE splicer.end(TSTTG_CGM.CgmGeom.getPrnts)
}

/**
 * Method:  Move[]
 */
void
TSTTG_CGM::CgmGeom_impl::Move (
  /*inout*/ void*& geom_entity,
  /*in*/ double x,
  /*in*/ double y,
  /*in*/ double z ) 
throw ( 
  ::TSTTB::Error
){
    // DO-NOT-DELETE splicer.begin(TSTTG_CGM.CgmGeom.Move)
    // insert implementation here
  TSTTG_Move(tsttgInstance, &geom_entity, x, y, z);
  PROCESS_ERROR;
    // DO-NOT-DELETE splicer.end(TSTTG_CGM.CgmGeom.Move)
}

/**
 * Method:  Rotate[]
 */
void
TSTTG_CGM::CgmGeom_impl::Rotate (
  /*inout*/ void*& geom_entity,
  /*in*/ double angle,
  /*in*/ double axis_normal_x,
  /*in*/ double axis_normal_y,
  /*in*/ double axis_normal_z ) 
throw ( 
  ::TSTTB::Error
){
    // DO-NOT-DELETE splicer.begin(TSTTG_CGM.CgmGeom.Rotate)
    // insert implementation here
  TSTTG_Rotate(tsttgInstance, &geom_entity, angle, axis_normal_x, 
               axis_normal_y, axis_normal_z);
  PROCESS_ERROR;
    // DO-NOT-DELETE splicer.end(TSTTG_CGM.CgmGeom.Rotate)
}

/**
 * Method:  Reflect[]
 */
void
TSTTG_CGM::CgmGeom_impl::Reflect (
  /*inout*/ void*& geom_entity,
  /*in*/ double plane_normal_x,
  /*in*/ double plane_normal_y,
  /*in*/ double plane_normal_z ) 
throw ( 
  ::TSTTB::Error
){
    // DO-NOT-DELETE splicer.begin(TSTTG_CGM.CgmGeom.Reflect)
    // insert implementation here
  TSTTG_Reflect(tsttgInstance, &geom_entity, plane_normal_x, 
                plane_normal_y, plane_normal_z);
  PROCESS_ERROR;
    // DO-NOT-DELETE splicer.end(TSTTG_CGM.CgmGeom.Reflect)
}

/**
 * Method:  Scale[]
 */
void
TSTTG_CGM::CgmGeom_impl::Scale (
  /*inout*/ void*& geom_entity,
  /*in*/ double scale_x,
  /*in*/ double scale_y,
  /*in*/ double scale_z ) 
throw ( 
  ::TSTTB::Error
){
  // DO-NOT-DELETE splicer.begin(TSTTG_CGM.CgmGeom.Scale)
  // insert implementation here
  TSTTG_Scale(tsttgInstance, &geom_entity, scale_x, scale_y, scale_z);
  PROCESS_ERROR;
  // DO-NOT-DELETE splicer.end(TSTTG_CGM.CgmGeom.Scale)
}

/**
 * Method:  Copy[]
 */
void
TSTTG_CGM::CgmGeom_impl::Copy (
  /*in*/ void* geom_entity,
  /*out*/ void*& geom_entity2 ) 
throw ( 
  ::TSTTB::Error
){
    // DO-NOT-DELETE splicer.begin(TSTTG_CGM.CgmGeom.Copy)
    // insert implementation here
  TSTTG_Copy(tsttgInstance, geom_entity, &geom_entity2);
  PROCESS_ERROR;
    // DO-NOT-DELETE splicer.end(TSTTG_CGM.CgmGeom.Copy)
}

/**
 * Method:  SweepAboutAxis[]
 */
void
TSTTG_CGM::CgmGeom_impl::SweepAboutAxis (
  /*in*/ void* geom_entity,
  /*in*/ double angle,
  /*in*/ double axis_normal_x,
  /*in*/ double axis_normal_y,
  /*in*/ double axis_normal_z,
  /*out*/ void*& geom_entity2 ) 
throw ( 
  ::TSTTB::Error
){
    // DO-NOT-DELETE splicer.begin(TSTTG_CGM.CgmGeom.SweepAboutAxis)
    // insert implementation here
  TSTTG_SweepAboutAxis(tsttgInstance, geom_entity, angle,
                       axis_normal_x, axis_normal_y, axis_normal_z, 
                       &geom_entity2);
  PROCESS_ERROR;
    // DO-NOT-DELETE splicer.end(TSTTG_CGM.CgmGeom.SweepAboutAxis)
}

/**
 * Method:  Delete[]
 */
void
TSTTG_CGM::CgmGeom_impl::Delete (
  /*in*/ void* geom_entity ) 
throw ( 
  ::TSTTB::Error
){
    // DO-NOT-DELETE splicer.begin(TSTTG_CGM.CgmGeom.Delete)
    // insert implementation here
  TSTTG_Delete(tsttgInstance, geom_entity);
  PROCESS_ERROR;
    // DO-NOT-DELETE splicer.end(TSTTG_CGM.CgmGeom.Delete)
}

/**
 * Method:  Unite[]
 */
void
TSTTG_CGM::CgmGeom_impl::Unite (
  /*in*/ ::sidl::array<void*> geom_entities,
  /*in*/ int32_t geom_entities_size,
  /*out*/ void*& geom_entity ) 
throw ( 
  ::TSTTB::Error
){
    // DO-NOT-DELETE splicer.begin(TSTTG_CGM.CgmGeom.Unite)
    // insert implementation here
  TSTTG_Unite(tsttgInstance, 
              TEMP_ARRAY_IN(geom_entities), &geom_entity);
  PROCESS_ERROR;
    // DO-NOT-DELETE splicer.end(TSTTG_CGM.CgmGeom.Unite)
}

/**
 * Method:  Subtract[]
 */
void
TSTTG_CGM::CgmGeom_impl::Subtract (
  /*in*/ void* blank,
  /*in*/ void* tool,
  /*out*/ void*& geom_entity ) 
throw ( 
  ::TSTTB::Error
){
    // DO-NOT-DELETE splicer.begin(TSTTG_CGM.CgmGeom.Subtract)
    // insert implementation here
  TSTTG_Subtract(tsttgInstance, blank, tool, &geom_entity);
  PROCESS_ERROR;
    // DO-NOT-DELETE splicer.end(TSTTG_CGM.CgmGeom.Subtract)
}

/**
 * Method:  Section[]
 */
void
TSTTG_CGM::CgmGeom_impl::Section (
  /*inout*/ void*& geom_entity,
  /*in*/ double plane_normal_x,
  /*in*/ double plane_normal_y,
  /*in*/ double plane_normal_z,
  /*in*/ double offset,
  /*in*/ bool reverse,
  /*out*/ void*& geom_entity2 ) 
throw ( 
  ::TSTTB::Error
){
    // DO-NOT-DELETE splicer.begin(TSTTG_CGM.CgmGeom.Section)
    // insert implementation here
  TSTTG_Section(tsttgInstance, &geom_entity, plane_normal_x, 
                plane_normal_y, plane_normal_z, 
                offset, reverse, &geom_entity2);
  PROCESS_ERROR;
    // DO-NOT-DELETE splicer.end(TSTTG_CGM.CgmGeom.Section)
}

/**
 * Method:  Imprint[]
 */
void
TSTTG_CGM::CgmGeom_impl::Imprint (
  /*in*/ ::sidl::array<void*> geom_entities,
  /*in*/ int32_t geom_entities_size ) 
throw ( 
  ::TSTTB::Error
){
  // DO-NOT-DELETE splicer.begin(TSTTG_CGM.CgmGeom.Imprint)
  // insert implementation here
  TSTTG_Imprint (tsttgInstance,
                 TEMP_ARRAY_IN(geom_entities));
  PROCESS_ERROR;

  // DO-NOT-DELETE splicer.end(TSTTG_CGM.CgmGeom.Imprint)
}

/**
 * Method:  Merge[]
 */
void
TSTTG_CGM::CgmGeom_impl::Merge (
  /*in*/ ::sidl::array<void*> geom_entities,
  /*in*/ int32_t geom_entities_size,
  /*in*/ double tolerance ) 
throw ( 
  ::TSTTB::Error
){
  // DO-NOT-DELETE splicer.begin(TSTTG_CGM.CgmGeom.Merge)
  // insert implementation here
  TSTTG_Merge (tsttgInstance,
               TEMP_ARRAY_IN(geom_entities),
               tolerance);
  PROCESS_ERROR;

  // DO-NOT-DELETE splicer.end(TSTTG_CGM.CgmGeom.Merge)
}


// DO-NOT-DELETE splicer.begin(TSTTG_CGM.CgmGeom._misc)
// Put miscellaneous code here
// call this function when you want to throw an error
void
TSTTG_CGM::CgmGeom_impl::processError() throw(::TSTTB::Error)
{
  static void *behavior_tag = NULL;
  static ::TSTTB::ErrorActions action = ::TSTTB::ErrorActions_THROW_ERROR;

    // save this info before calling tag get function to get behavior
  std::string this_desc(TSTTG_LAST_ERROR.description);
  TSTTB::ErrorType this_err = (TSTTB::ErrorType)TSTTG_LAST_ERROR.error_type;

  if (0 == behavior_tag)
    behavior_tag = TSTTG_getTagHandle(tsttgInstance, "Error_Behavior");

  if (0 != behavior_tag) {
    ::TSTTB::ErrorActions temp = (::TSTTB::ErrorActions) TSTTG_getIntData (tsttgInstance, NULL,
                                                                         behavior_tag);
    if (TSTTB_SUCCESS == TSTTG_LAST_ERROR.error_type) action = temp;
  }

  ::TSTTB::Error this_error = ::TSTTB::Error::_create();
  switch (action) {
    case ::TSTTB::ErrorActions_THROW_ERROR:
      this_error.set(this_err, this_desc);
      throw this_error;
      break;
    case ::TSTTB::ErrorActions_WARN_ONLY:
      std::cerr << this_desc << std::endl;
      break;
    case ::TSTTB::ErrorActions_SILENT:
      break;
  }
}

// DO-NOT-DELETE splicer.end(TSTTG_CGM.CgmGeom._misc)
