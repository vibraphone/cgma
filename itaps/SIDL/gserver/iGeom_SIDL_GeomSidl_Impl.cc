// 
// File:          iGeom_SIDL_GeomSidl_Impl.cc
// Symbol:        iGeom_SIDL.GeomSidl-v0.1
// Symbol Type:   class
// Babel Version: 0.10.10
// sidl Created:  20090126 13:13:24 CST
// Generated:     20090126 13:13:26 CST
// Description:   Server-side implementation for iGeom_SIDL.GeomSidl
// 
// WARNING: Automatically generated; only changes within splicers preserved
// 
// babel-version = 0.10.10
// source-line   = 5
// source-url    = file:/home/jason/meshkit/cgm/itaps/SIDL/iGeom_SIDL.sidl
// 
#include "iGeom_SIDL_GeomSidl_Impl.hh"

// DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl._includes)
// Insert-Code-Here {iGeom_SIDL.GeomSidl._includes} (additional includes or code)
#include "sidlArray.h"
#include "iBase.hh"
#include "iGeom.h"
#include <iostream>

#define PROCESS_ERROR if (iGeom_LAST_ERROR.error_type != iBase_SUCCESS) this->processError()
#define PROCESS_ERROR_MSG(a,b) \
   iGeom_LAST_ERROR.error_type = a; \
   sprintf(iGeom_LAST_ERROR.description, "%s", b);\
   this->processError()
#define GENERATE_ERROR(a) {iGeom_LAST_ERROR.error_type = a; this->processError();}

// need this for definitions in TSTTB_SNL_SIDL_defs.h
#define LOCAL_iBase_ERROR iGeom_LAST_ERROR
#include "iBase_SIDL_defs.h"
#include <string>
// DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl._includes)

// user-defined constructor.
void iGeom_SIDL::GeomSidl_impl::_ctor() {
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl._ctor)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl._ctor} (constructor)
  iGeom_ctor(&igeomInstance);
  PROCESS_ERROR;
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl._ctor)
}

// user-defined destructor.
void iGeom_SIDL::GeomSidl_impl::_dtor() {
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl._dtor)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl._dtor} (destructor)
  iGeom_dtor(igeomInstance);
  PROCESS_ERROR;
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl._dtor)
}

// static class initializer.
void iGeom_SIDL::GeomSidl_impl::_load() {
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl._load)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl._load} (class initialization)
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl._load)
}

// user-defined static methods: (none)

// user-defined non-static methods:
/**
 * Method:  createTag[]
 */
void
iGeom_SIDL::GeomSidl_impl::createTag (
  /* in */ const ::std::string& tag_name,
  /* in */ int32_t number_of_values,
  /* in */ ::iBase::TagValueType tag_type,
  /* out */ void*& tag_handle ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.createTag)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.createTag} (createTag method)
  if (tag_type > ::iBase::TagValueType_BYTES || 
      tag_type < ::iBase::TagValueType_INTEGER) {
    PROCESS_ERROR_MSG(iBase_INVALID_ARGUMENT, "Wrong tag value type");
    return;
  }
  iGeom_createTag(igeomInstance, tag_name.c_str(), number_of_values, (iBase_TagValueType)tag_type, 
                  &tag_handle);
  PROCESS_ERROR;
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.createTag)
}

/**
 * Method:  destroyTag[]
 */
void
iGeom_SIDL::GeomSidl_impl::destroyTag (
  /* in */ void* tag_handle,
  /* in */ bool forced ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.destroyTag)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.destroyTag} (destroyTag method)
  iGeom_destroyTag(igeomInstance, tag_handle, forced);
  PROCESS_ERROR;
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.destroyTag)
}

/**
 * Method:  getTagName[]
 */
void
iGeom_SIDL::GeomSidl_impl::getTagName (
  /* in */ void* tag_handle,
  /* out */ ::std::string& tag_name ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.getTagName)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.getTagName} (getTagName method)
  std::string retval = iGeom_getTagName(igeomInstance, tag_handle);
  PROCESS_ERROR;
  return retval;
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.getTagName)
}

/**
 * Method:  getTagSizeValues[]
 */
void
iGeom_SIDL::GeomSidl_impl::getTagSizeValues (
  /* in */ void* tag_handle,
  /* out */ int32_t& size_values ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.getTagSizeValues)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.getTagSizeValues} (getTagSizeValues method)
  int32_t retval = iGeom_getTagSizeValues(igeomInstance, tag_handle);
  PROCESS_ERROR;
  return retval;
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.getTagSizeValues)
}

/**
 * Method:  getTagSizeBytes[]
 */
void
iGeom_SIDL::GeomSidl_impl::getTagSizeBytes (
  /* in */ void* tag_handle,
  /* out */ int32_t& size_bytes ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.getTagSizeBytes)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.getTagSizeBytes} (getTagSizeBytes method)
  int32_t retval = iGeom_getTagSizeBytes(igeomInstance, tag_handle);
  PROCESS_ERROR;
  return retval;
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.getTagSizeBytes)
}

/**
 * Method:  getTagHandle[]
 */
void
iGeom_SIDL::GeomSidl_impl::getTagHandle (
  /* in */ const ::std::string& tag_name,
  /* out */ void*& tag_handle ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.getTagHandle)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.getTagHandle} (getTagHandle method)
  void *retval = iGeom_getTagHandle(igeomInstance, tag_name.c_str());
  PROCESS_ERROR;
  return retval;
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.getTagHandle)
}

/**
 * Method:  getTagType[]
 */
void
iGeom_SIDL::GeomSidl_impl::getTagType (
  /* in */ void* tag_handle,
  /* out */ ::iBase::TagValueType& tag_data_type ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.getTagType)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.getTagType} (getTagType method)
  ::iBase::TagValueType retval = (::iBase::TagValueType) iGeom_getTagType(igeomInstance,
                                                                        tag_handle );
  PROCESS_ERROR;
  return retval;
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.getTagType)
}

/**
 * Method:  getData[]
 */
void
iGeom_SIDL::GeomSidl_impl::getData (
  /* in */ void* entity_handle,
  /* in */ void* tag_handle,
  /* inout */ ::sidl::array<char>& tag_value,
  /* out */ int32_t& tag_value_size ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.getData)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.getData} (getData method)
  CREATE_TEMP_TAG_ARRAY(tag_value);
  
  iGeom_getData (igeomInstance,
                 entity_handle,
                 tag_handle,
                 TEMP_ARRAY_INOUT(tag_value));
  PROCESS_ERROR;

  ASSIGN_TAG_ARRAY(tag_value);
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.getData)
}

/**
 * Method:  getIntData[]
 */
void
iGeom_SIDL::GeomSidl_impl::getIntData (
  /* in */ void* entity_handle,
  /* in */ void* tag_handle,
  /* out */ int32_t& int_data ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.getIntData)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.getIntData} (getIntData method)
  int32_t retval = iGeom_getIntData(igeomInstance,
                                    entity_handle,
                                    tag_handle );
  PROCESS_ERROR;
  return retval;
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.getIntData)
}

/**
 * Method:  getDblData[]
 */
void
iGeom_SIDL::GeomSidl_impl::getDblData (
  /* in */ void* entity_handle,
  /* in */ void* tag_handle,
  /* out */ double& dbl_data ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.getDblData)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.getDblData} (getDblData method)
  double retval = iGeom_getDblData(igeomInstance,
                                   entity_handle,
                                   tag_handle );  
  PROCESS_ERROR;
  return retval;
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.getDblData)
}

/**
 * Method:  getEHData[]
 */
void
iGeom_SIDL::GeomSidl_impl::getEHData (
  /* in */ void* entity_handle,
  /* in */ void* tag_handle,
  /* out */ void*& eh_data ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.getEHData)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.getEHData} (getEHData method)
  void *retval = iGeom_getEHData(igeomInstance,
                                 entity_handle,
                                 tag_handle );
  PROCESS_ERROR;
  return retval;
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.getEHData)
}

/**
 * Method:  setData[]
 */
void
iGeom_SIDL::GeomSidl_impl::setData (
  /* in */ void* entity_handle,
  /* in */ void* tag_handle,
  /* in */ ::sidl::array<char> tag_value,
  /* in */ int32_t tag_value_size ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.setData)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.setData} (setData method)
  iGeom_setData(igeomInstance, 
                entity_handle,
                tag_handle, 
                TEMP_TAG_ARRAY_IN(tag_value));
  PROCESS_ERROR;
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.setData)
}

/**
 * Method:  setIntData[]
 */
void
iGeom_SIDL::GeomSidl_impl::setIntData (
  /* in */ void* entity_handle,
  /* in */ void* tag_handle,
  /* in */ int32_t tag_value ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.setIntData)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.setIntData} (setIntData method)
  iGeom_setIntData (igeomInstance,
                    entity_handle,
                    tag_handle,
                    tag_value );
  PROCESS_ERROR;
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.setIntData)
}

/**
 * Method:  setDblData[]
 */
void
iGeom_SIDL::GeomSidl_impl::setDblData (
  /* in */ void* entity_handle,
  /* in */ void* tag_handle,
  /* in */ double tag_value ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.setDblData)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.setDblData} (setDblData method)
  iGeom_setDblData (igeomInstance,
                    entity_handle,
                    tag_handle,
                    tag_value );
  PROCESS_ERROR;
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.setDblData)
}

/**
 * Method:  setEHData[]
 */
void
iGeom_SIDL::GeomSidl_impl::setEHData (
  /* in */ void* entity_handle,
  /* in */ void* tag_handle,
  /* in */ void* tag_value ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.setEHData)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.setEHData} (setEHData method)
  iGeom_setEHData (igeomInstance,
                   entity_handle,
                   tag_handle,
                   tag_value );
  PROCESS_ERROR;
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.setEHData)
}

/**
 * Method:  getAllTags[]
 */
void
iGeom_SIDL::GeomSidl_impl::getAllTags (
  /* in */ void* entity_handle,
  /* inout */ ::sidl::array<void*>& tag_handles,
  /* out */ int32_t& tag_handles_size ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.getAllTags)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.getAllTags} (getAllTags method)
  CREATE_TEMP_ARRAY(void*, tag_handles);
  
  iGeom_getAllTags(igeomInstance, entity_handle, 
                   TEMP_ARRAY_INOUT(tag_handles));
  PROCESS_ERROR;

  ASSIGN_ARRAY(tag_handles);
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.getAllTags)
}

/**
 * Method:  rmvTag[]
 */
void
iGeom_SIDL::GeomSidl_impl::rmvTag (
  /* in */ void* entity_handle,
  /* in */ void* tag_handle ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.rmvTag)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.rmvTag} (rmvTag method)
  iGeom_rmvTag (igeomInstance,
                entity_handle,
                tag_handle );
  PROCESS_ERROR;
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.rmvTag)
}

/**
 * Method:  getArrData[]
 */
void
iGeom_SIDL::GeomSidl_impl::getArrData (
  /* in */ ::sidl::array<void*> entity_handles,
  /* in */ int32_t entity_handles_size,
  /* in */ void* tag_handle,
  /* inout */ ::sidl::array<char>& value_array,
  /* out */ int32_t& value_array_size ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.getArrData)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.getArrData} (getArrData method)
  CREATE_TEMP_TAG_ARRAY(value_array);
  
  iGeom_getArrData(igeomInstance, 
                   TEMP_ARRAY_IN(entity_handles), 
                   tag_handle, 
                   TEMP_ARRAY_INOUT(value_array));
  PROCESS_ERROR;

  ASSIGN_TAG_ARRAY(value_array);
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.getArrData)
}

/**
 * Method:  getIntArrData[]
 */
void
iGeom_SIDL::GeomSidl_impl::getIntArrData (
  /* in */ ::sidl::array<void*> entity_handles,
  /* in */ int32_t entity_handles_size,
  /* in */ void* tag_handle,
  /* inout */ ::sidl::array<int32_t>& value_array,
  /* out */ int32_t& value_array_size ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.getIntArrData)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.getIntArrData} (getIntArrData method)
  CREATE_TEMP_ARRAY(int32_t, value_array);
  
  iGeom_getIntArrData (igeomInstance,
                       TEMP_ARRAY_IN(entity_handles),
                       tag_handle,
                       TEMP_ARRAY_INOUT(value_array));
  PROCESS_ERROR;

  ASSIGN_ARRAY(value_array);
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.getIntArrData)
}

/**
 * Method:  getDblArrData[]
 */
void
iGeom_SIDL::GeomSidl_impl::getDblArrData (
  /* in */ ::sidl::array<void*> entity_handles,
  /* in */ int32_t entity_handles_size,
  /* in */ void* tag_handle,
  /* inout */ ::sidl::array<double>& value_array,
  /* out */ int32_t& value_array_size ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.getDblArrData)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.getDblArrData} (getDblArrData method)
  CREATE_TEMP_ARRAY(double, value_array);
  
  iGeom_getDblArrData (igeomInstance,
                       TEMP_ARRAY_IN(entity_handles),
                       tag_handle,
                       TEMP_ARRAY_INOUT(value_array));
  PROCESS_ERROR;

  ASSIGN_ARRAY(value_array);
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.getDblArrData)
}

/**
 * Method:  getEHArrData[]
 */
void
iGeom_SIDL::GeomSidl_impl::getEHArrData (
  /* in */ ::sidl::array<void*> entity_handles,
  /* in */ int32_t entity_handles_size,
  /* in */ void* tag_handle,
  /* inout */ ::sidl::array<void*>& value_array,
  /* out */ int32_t& value_array_size ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.getEHArrData)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.getEHArrData} (getEHArrData method)
  CREATE_TEMP_ARRAY(void*, value_array);
  
  iGeom_getEHArrData (igeomInstance,
                      TEMP_ARRAY_IN(entity_handles),
                      tag_handle,
                      TEMP_ARRAY_INOUT(value_array));
  PROCESS_ERROR;

  ASSIGN_ARRAY(value_array);
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.getEHArrData)
}

/**
 * Method:  setArrData[]
 */
void
iGeom_SIDL::GeomSidl_impl::setArrData (
  /* in */ ::sidl::array<void*> entity_handles,
  /* in */ int32_t entity_handles_size,
  /* in */ void* tag_handle,
  /* in */ ::sidl::array<char> value_array,
  /* in */ int32_t value_array_size ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.setArrData)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.setArrData} (setArrData method)
  iGeom_setArrData(igeomInstance, 
                   TEMP_ARRAY_IN(entity_handles),
                   tag_handle, 
                   TEMP_TAG_ARRAY_IN(value_array));
  PROCESS_ERROR;
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.setArrData)
}

/**
 * Method:  setIntArrData[]
 */
void
iGeom_SIDL::GeomSidl_impl::setIntArrData (
  /* in */ ::sidl::array<void*> entity_handles,
  /* in */ int32_t entity_handles_size,
  /* in */ void* tag_handle,
  /* in */ ::sidl::array<int32_t> value_array,
  /* in */ int32_t value_array_size ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.setIntArrData)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.setIntArrData} (setIntArrData method)
  iGeom_setIntArrData (igeomInstance,
                       TEMP_ARRAY_IN(entity_handles),
                       tag_handle,
                       TEMP_ARRAY_IN(value_array));
  PROCESS_ERROR;
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.setIntArrData)
}

/**
 * Method:  setDblArrData[]
 */
void
iGeom_SIDL::GeomSidl_impl::setDblArrData (
  /* in */ ::sidl::array<void*> entity_handles,
  /* in */ int32_t entity_handles_size,
  /* in */ void* tag_handle,
  /* in */ ::sidl::array<double> value_array,
  /* in */ int32_t value_array_size ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.setDblArrData)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.setDblArrData} (setDblArrData method)
  iGeom_setDblArrData (igeomInstance,
                       TEMP_ARRAY_IN(entity_handles),
                       tag_handle,
                       TEMP_ARRAY_IN(value_array));
  PROCESS_ERROR;
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.setDblArrData)
}

/**
 * Method:  setEHArrData[]
 */
void
iGeom_SIDL::GeomSidl_impl::setEHArrData (
  /* in */ ::sidl::array<void*> entity_handles,
  /* in */ int32_t entity_handles_size,
  /* in */ void* tag_handle,
  /* in */ ::sidl::array<void*> value_array,
  /* in */ int32_t value_array_size ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.setEHArrData)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.setEHArrData} (setEHArrData method)
  iGeom_setEHArrData (igeomInstance,
                      TEMP_ARRAY_IN(entity_handles),
                      tag_handle,
                      TEMP_ARRAY_IN(value_array));
  PROCESS_ERROR;
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.setEHArrData)
}

/**
 * Method:  rmvArrTag[]
 */
void
iGeom_SIDL::GeomSidl_impl::rmvArrTag (
  /* in */ ::sidl::array<void*> entity_handles,
  /* in */ int32_t entity_handles_size,
  /* in */ void* tag_handle ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.rmvArrTag)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.rmvArrTag} (rmvArrTag method)
  iGeom_rmvArrTag (igeomInstance,
                   TEMP_ARRAY_IN(entity_handles),
                   tag_handle );
  PROCESS_ERROR;
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.rmvArrTag)
}

/**
 * Method:  setEntSetData[]
 */
void
iGeom_SIDL::GeomSidl_impl::setEntSetData (
  /* in */ void* entity_set,
  /* in */ void* tag_handle,
  /* inout */ ::sidl::array<char>& tag_value,
  /* in */ int32_t tag_value_size ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.setEntSetData)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.setEntSetData} (setEntSetData method)
  iGeom_setEntSetData(igeomInstance, entity_set, tag_handle,
                      TEMP_TAG_ARRAY_IN(tag_value));
  PROCESS_ERROR;
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.setEntSetData)
}

/**
 * Method:  setEntSetIntData[]
 */
void
iGeom_SIDL::GeomSidl_impl::setEntSetIntData (
  /* in */ void* entity_set,
  /* in */ void* tag_handle,
  /* in */ int32_t tag_value ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.setEntSetIntData)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.setEntSetIntData} (setEntSetIntData method)
  iGeom_setEntSetIntData (igeomInstance,
                          entity_set,
                          tag_handle,
                          tag_value );
  PROCESS_ERROR;
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.setEntSetIntData)
}

/**
 * Method:  setEntSetDblData[]
 */
void
iGeom_SIDL::GeomSidl_impl::setEntSetDblData (
  /* in */ void* entity_set,
  /* in */ void* tag_handle,
  /* in */ double tag_value ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.setEntSetDblData)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.setEntSetDblData} (setEntSetDblData method)
  iGeom_setEntSetDblData (igeomInstance,
                          entity_set,
                          tag_handle,
                          tag_value );
  PROCESS_ERROR;
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.setEntSetDblData)
}

/**
 * Method:  setEntSetEHData[]
 */
void
iGeom_SIDL::GeomSidl_impl::setEntSetEHData (
  /* in */ void* entity_set,
  /* in */ void* tag_handle,
  /* in */ void* tag_value ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.setEntSetEHData)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.setEntSetEHData} (setEntSetEHData method)
  iGeom_setEntSetEHData (igeomInstance,
                         entity_set,
                         tag_handle,
                         tag_value );
  PROCESS_ERROR;
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.setEntSetEHData)
}

/**
 * Method:  getEntSetData[]
 */
void
iGeom_SIDL::GeomSidl_impl::getEntSetData (
  /* in */ void* entity_set,
  /* in */ void* tag_handle,
  /* inout */ ::sidl::array<char>& tag_value,
  /* out */ int32_t& tag_value_size ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.getEntSetData)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.getEntSetData} (getEntSetData method)
  CREATE_TEMP_TAG_ARRAY(tag_value);
  
  iGeom_getEntSetData (igeomInstance,
                       entity_set,
                       tag_handle,
                       TEMP_ARRAY_INOUT(tag_value));
  PROCESS_ERROR;

  ASSIGN_TAG_ARRAY(tag_value);
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.getEntSetData)
}

/**
 * Method:  getEntSetIntData[]
 */
void
iGeom_SIDL::GeomSidl_impl::getEntSetIntData (
  /* in */ void* entity_set,
  /* in */ void* tag_handle,
  /* out */ int32_t& int_data ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.getEntSetIntData)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.getEntSetIntData} (getEntSetIntData method)
  int32_t retval = iGeom_getEntSetIntData(igeomInstance,
                                          entity_set,
                                          tag_handle );
  PROCESS_ERROR;
  return retval;
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.getEntSetIntData)
}

/**
 * Method:  getEntSetDblData[]
 */
void
iGeom_SIDL::GeomSidl_impl::getEntSetDblData (
  /* in */ void* entity_set,
  /* in */ void* tag_handle,
  /* out */ double& dbl_data ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.getEntSetDblData)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.getEntSetDblData} (getEntSetDblData method)
  double retval = iGeom_getEntSetDblData(igeomInstance,
                                         entity_set,
                                         tag_handle );
  PROCESS_ERROR;
  return retval;
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.getEntSetDblData)
}

/**
 * Method:  getEntSetEHData[]
 */
void
iGeom_SIDL::GeomSidl_impl::getEntSetEHData (
  /* in */ void* entity_set,
  /* in */ void* tag_handle,
  /* out */ void*& eh_data ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.getEntSetEHData)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.getEntSetEHData} (getEntSetEHData method)
  void *retval = iGeom_getEntSetEHData(igeomInstance,
                                       entity_set,
                                       tag_handle );
  PROCESS_ERROR;
  return retval;
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.getEntSetEHData)
}

/**
 * Method:  getAllEntSetTags[]
 */
void
iGeom_SIDL::GeomSidl_impl::getAllEntSetTags (
  /* in */ void* entity_set,
  /* inout */ ::sidl::array<void*>& tag_handles,
  /* out */ int32_t& tag_handles_size ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.getAllEntSetTags)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.getAllEntSetTags} (getAllEntSetTags method)
  CREATE_TEMP_ARRAY(void*, tag_handles);
  
  iGeom_getAllEntSetTags (igeomInstance,
                          entity_set,
                          TEMP_ARRAY_INOUT(tag_handles));
  PROCESS_ERROR;

  ASSIGN_ARRAY(tag_handles);
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.getAllEntSetTags)
}

/**
 * Method:  rmvEntSetTag[]
 */
void
iGeom_SIDL::GeomSidl_impl::rmvEntSetTag (
  /* in */ void* entity_set,
  /* in */ void* tag_handle ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.rmvEntSetTag)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.rmvEntSetTag} (rmvEntSetTag method)
  iGeom_rmvEntSetTag (igeomInstance,
                      entity_set,
                      tag_handle );
  PROCESS_ERROR;
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.rmvEntSetTag)
}

/**
 * Method:  createEntSet[]
 */
void
iGeom_SIDL::GeomSidl_impl::createEntSet (
  /* in */ bool isList,
  /* out */ void*& entity_set ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.createEntSet)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.createEntSet} (createEntSet method)
  iGeom_createEntSet(igeomInstance, isList, &entity_set);
  PROCESS_ERROR;
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.createEntSet)
}

/**
 * Method:  destroyEntSet[]
 */
void
iGeom_SIDL::GeomSidl_impl::destroyEntSet (
  /* in */ void* entity_set ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.destroyEntSet)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.destroyEntSet} (destroyEntSet method)
  iGeom_destroyEntSet(igeomInstance, entity_set);
  PROCESS_ERROR;
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.destroyEntSet)
}

/**
 * Method:  isList[]
 */
void
iGeom_SIDL::GeomSidl_impl::isList (
  /* in */ void* entity_set,
  /* out */ int32_t& is_list ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.isList)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.isList} (isList method)
  bool retval = iGeom_isList(igeomInstance, entity_set);
  PROCESS_ERROR;
  return retval;
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.isList)
}

/**
 * Method:  getNumEntSets[]
 */
void
iGeom_SIDL::GeomSidl_impl::getNumEntSets (
  /* in */ void* entity_set,
  /* in */ int32_t num_hops,
  /* out */ int32_t& num_sets ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.getNumEntSets)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.getNumEntSets} (getNumEntSets method)
  int32_t retval = iGeom_getNumEntSets(igeomInstance, entity_set, num_hops);
  PROCESS_ERROR;
  return retval;
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.getNumEntSets)
}

/**
 * Method:  getEntSets[]
 */
void
iGeom_SIDL::GeomSidl_impl::getEntSets (
  /* in */ void* entity_set,
  /* in */ int32_t num_hops,
  /* inout */ ::sidl::array<void*>& contained_entset_handles,
  /* out */ int32_t& contained_entset_handles_size ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.getEntSets)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.getEntSets} (getEntSets method)
  CREATE_TEMP_ARRAY(void*, contained_entset_handles);

  iGeom_getEntSets(igeomInstance, entity_set, num_hops, 
                   TEMP_ARRAY_INOUT(contained_entset_handles));
  PROCESS_ERROR;
  ASSIGN_ARRAY(contained_entset_handles);
  CREATE_TEMP_ARRAY(void*, contained_entset_handles);

  iGeom_getEntSets(igeomInstance, entity_set, num_hops, 
                   TEMP_ARRAY_INOUT(contained_entset_handles));
  PROCESS_ERROR;
  ASSIGN_ARRAY(contained_entset_handles);
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.getEntSets)
}

/**
 * Method:  addEntToSet[]
 */
void
iGeom_SIDL::GeomSidl_impl::addEntToSet (
  /* in */ void* entity_handle,
  /* in */ void* entity_set ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.addEntToSet)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.addEntToSet} (addEntToSet method)
  iGeom_addEntToSet (igeomInstance,
                     entity_handle,
                     &entity_set);
  PROCESS_ERROR;
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.addEntToSet)
}

/**
 * Method:  rmvEntFromSet[]
 */
void
iGeom_SIDL::GeomSidl_impl::rmvEntFromSet (
  /* in */ void* entity_handle,
  /* in */ void* entity_set ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.rmvEntFromSet)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.rmvEntFromSet} (rmvEntFromSet method)
  iGeom_rmvEntFromSet (igeomInstance,
                       entity_handle,
                       &entity_set );
  PROCESS_ERROR;
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.rmvEntFromSet)
}

/**
 * Method:  addEntArrToSet[]
 */
void
iGeom_SIDL::GeomSidl_impl::addEntArrToSet (
  /* in */ ::sidl::array<void*> entity_handles,
  /* in */ int32_t entity_handles_size,
  /* in */ void* entity_set ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.addEntArrToSet)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.addEntArrToSet} (addEntArrToSet method)
  entity_handles.addRef();
  void **tmp_ptr1 = entity_handles._get_ior()->d_firstElement;
    //void **tmp_ptr3 = entity_handles._get_ior()->d_firstElement;
  iGeom_addEntArrToSet(igeomInstance, tmp_ptr1,
                       entity_handles_size,
                       &entity_set);
  PROCESS_ERROR;
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.addEntArrToSet)
}

/**
 * Method:  rmvEntArrFromSet[]
 */
void
iGeom_SIDL::GeomSidl_impl::rmvEntArrFromSet (
  /* in */ ::sidl::array<void*> entity_handles,
  /* in */ int32_t entity_handles_size,
  /* in */ void* entity_set ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.rmvEntArrFromSet)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.rmvEntArrFromSet} (rmvEntArrFromSet method)
  iGeom_rmvEntArrFromSet(igeomInstance, TEMP_ARRAY_IN(entity_handles), 
                         &entity_set);
  PROCESS_ERROR;
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.rmvEntArrFromSet)
}

/**
 * Method:  addEntSet[]
 */
void
iGeom_SIDL::GeomSidl_impl::addEntSet (
  /* in */ void* entity_set_to_add,
  /* in */ void* entity_set_handle ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.addEntSet)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.addEntSet} (addEntSet method)
  iGeom_addEntSet(igeomInstance, entity_set_to_add, 
                  &entity_set_handle);
  PROCESS_ERROR;
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.addEntSet)
}

/**
 * Method:  rmvEntSet[]
 */
void
iGeom_SIDL::GeomSidl_impl::rmvEntSet (
  /* in */ void* entity_set_to_remove,
  /* in */ void* entity_set_handle ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.rmvEntSet)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.rmvEntSet} (rmvEntSet method)
  iGeom_rmvEntSet (igeomInstance,
                   entity_set_to_remove,
                   &entity_set_handle );
  PROCESS_ERROR;
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.rmvEntSet)
}

/**
 * Method:  isEntContained[]
 */
void
iGeom_SIDL::GeomSidl_impl::isEntContained (
  /* in */ void* containing_entity_set,
  /* in */ void* entity_handle,
  /* out */ int32_t& is_contained ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.isEntContained)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.isEntContained} (isEntContained method)
  bool retval = iGeom_isEntContained(igeomInstance,
                                     containing_entity_set,
                                     entity_handle);
  PROCESS_ERROR;
  return retval;
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.isEntContained)
}

/**
 * Method:  isEntArrContained[]
 */
void
iGeom_SIDL::GeomSidl_impl::isEntArrContained (
  /* in */ void* containing_set,
  /* in */ ::sidl::array<void*> entity_handles,
  /* in */ int32_t entity_handles_size,
  /* inout */ ::sidl::array<int32_t>& is_contained,
  /* out */ int32_t& is_contained_size ) 
throw () 
{
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.isEntArrContained)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.isEntArrContained} (isEntArrContained method)
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.isEntArrContained)
}

/**
 * Method:  isEntSetContained[]
 */
void
iGeom_SIDL::GeomSidl_impl::isEntSetContained (
  /* in */ void* containing_entity_set,
  /* in */ void* contained_entity_set,
  /* out */ int32_t& is_contained ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.isEntSetContained)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.isEntSetContained} (isEntSetContained method)
  bool retval = iGeom_isEntSetContained(igeomInstance,
                                        containing_entity_set,
                                        contained_entity_set);
  PROCESS_ERROR;
  return retval;
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.isEntSetContained)
}

/**
 * Method:  addPrntChld[]
 */
void
iGeom_SIDL::GeomSidl_impl::addPrntChld (
  /* in */ void* parent_entity_set,
  /* in */ void* child_entity_set ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.addPrntChld)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.addPrntChld} (addPrntChld method)
  iGeom_addPrntChld(igeomInstance, &parent_entity_set, &child_entity_set);
  PROCESS_ERROR;
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.addPrntChld)
}

/**
 * Method:  rmvPrntChld[]
 */
void
iGeom_SIDL::GeomSidl_impl::rmvPrntChld (
  /* in */ void* parent_entity_set,
  /* in */ void* child_entity_set ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.rmvPrntChld)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.rmvPrntChld} (rmvPrntChld method)
  iGeom_rmvPrntChld(igeomInstance, &parent_entity_set, &child_entity_set);
  PROCESS_ERROR;
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.rmvPrntChld)
}

/**
 * Method:  isChildOf[]
 */
void
iGeom_SIDL::GeomSidl_impl::isChildOf (
  /* in */ void* parent_entity_set,
  /* in */ void* child_entity_set,
  /* out */ int32_t& is_child ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.isChildOf)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.isChildOf} (isChildOf method)
  bool retval = iGeom_isChildOf(igeomInstance, parent_entity_set, 
                                child_entity_set);
  PROCESS_ERROR;
  return retval;
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.isChildOf)
}

/**
 * Method:  getNumChld[]
 */
void
iGeom_SIDL::GeomSidl_impl::getNumChld (
  /* in */ void* entity_set,
  /* in */ int32_t num_hops,
  /* out */ int32_t& num_child ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.getNumChld)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.getNumChld} (getNumChld method)
  int32_t retval = iGeom_getNumChld(igeomInstance, entity_set, num_hops);
  PROCESS_ERROR;
  return retval;
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.getNumChld)
}

/**
 * Method:  getNumPrnt[]
 */
void
iGeom_SIDL::GeomSidl_impl::getNumPrnt (
  /* in */ void* entity_set,
  /* in */ int32_t num_hops,
  /* out */ int32_t& num_parent ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.getNumPrnt)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.getNumPrnt} (getNumPrnt method)
  int32_t retval = iGeom_getNumPrnt(igeomInstance, entity_set, num_hops);
  PROCESS_ERROR;
  return retval;
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.getNumPrnt)
}

/**
 * Method:  getChldn[]
 */
void
iGeom_SIDL::GeomSidl_impl::getChldn (
  /* in */ void* from_entity_set,
  /* in */ int32_t num_hops,
  /* inout */ ::sidl::array<void*>& child_handles,
  /* out */ int32_t& child_handles_size ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.getChldn)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.getChldn} (getChldn method)
  CREATE_TEMP_ARRAY(void*, child_handles);
  
  iGeom_getChldn(igeomInstance, from_entity_set, num_hops, 
                 TEMP_ARRAY_INOUT(child_handles));
  PROCESS_ERROR;

  ASSIGN_ARRAY(child_handles);
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.getChldn)
}

/**
 * Method:  getPrnts[]
 */
void
iGeom_SIDL::GeomSidl_impl::getPrnts (
  /* in */ void* from_entity_set,
  /* in */ int32_t num_hops,
  /* inout */ ::sidl::array<void*>& parent_handles,
  /* out */ int32_t& parent_handles_size ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.getPrnts)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.getPrnts} (getPrnts method)
  CREATE_TEMP_ARRAY(void*, parent_handles);
  
  iGeom_getPrnts(igeomInstance, from_entity_set, num_hops, 
                 TEMP_ARRAY_INOUT(parent_handles));
  PROCESS_ERROR;
  
  ASSIGN_ARRAY(parent_handles);
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.getPrnts)
}

/**
 * Method:  subtract[]
 */
void
iGeom_SIDL::GeomSidl_impl::subtract (
  /* in */ void* entity_set_1,
  /* in */ void* entity_set_2,
  /* out */ void*& result_entity_set ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.subtract)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.subtract} (subtract method)
  iGeom_subtract(igeomInstance, entity_set_1, entity_set_2, &result_entity_set);
  PROCESS_ERROR;
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.subtract)
}

/**
 * Method:  intersect[]
 */
void
iGeom_SIDL::GeomSidl_impl::intersect (
  /* in */ void* entity_set_1,
  /* in */ void* entity_set_2,
  /* out */ void*& result_entity_set ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.intersect)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.intersect} (intersect method)
  iGeom_intersect(igeomInstance, entity_set_1, entity_set_2, &result_entity_set);
  PROCESS_ERROR;
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.intersect)
}

/**
 * Method:  unite[]
 */
void
iGeom_SIDL::GeomSidl_impl::unite (
  /* in */ void* entity_set_1,
  /* in */ void* entity_set_2,
  /* out */ void*& result_entity_set ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.unite)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.unite} (unite method)
  iGeom_unite(igeomInstance, entity_set_1, entity_set_2, &result_entity_set);
  PROCESS_ERROR;
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.unite)
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
iGeom_SIDL::GeomSidl_impl::load (
  /* in */ const ::std::string& name,
  /* in */ const ::std::string& options ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.load)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.load} (load method)
  const char **tmp_options = NULL;
  if (options_size > 0) {
    tmp_options = new const char*[options_size];
    for (int i = 0; i < options_size; i++)
      tmp_options[i] = (options.get(i)).c_str();
  }
  
  iGeom_load(igeomInstance, name.c_str(), tmp_options, 
             options_size);
  
  delete [] tmp_options;

  PROCESS_ERROR;
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.load)
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
iGeom_SIDL::GeomSidl_impl::save (
  /* in */ const ::std::string& name,
  /* in */ const ::std::string& options ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.save)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.save} (save method)
  const char **tmp_options = NULL;
  if (options_size > 0) {
    tmp_options = new const char*[options_size];
    for (int i = 0; i < options_size; i++)
      tmp_options[i] = (options.get(i)).c_str();
  }
  
  iGeom_save(igeomInstance, name.c_str(), tmp_options, 
              options_size);
  
  delete [] tmp_options;

  PROCESS_ERROR;
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.save)
}

/**
 * Return gentities of specified dimension in this set, or in whole model.
 * @param set_handle Entity set being queried (if 0, whole model)
 * @param gentity_dimension Dimension of entities being queried
 * @param gentity_handles Gentity handles
 */
void
iGeom_SIDL::GeomSidl_impl::getEntities (
  /* in */ void* set_handle,
  /* in */ ::iBase::EntityType gentity_type,
  /* inout */ ::sidl::array<void*>& gentity_handles,
  /* inout */ int32_t& gentity_handles_size ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.getEntities)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.getEntities} (getEntities method)
  if (gentity_type > ::iGeom::EntityType_ALL_TYPES || 
      gentity_type < ::iGeom::EntityType_VERTEX) {
    PROCESS_ERROR_MSG(iBase_INVALID_ARGUMENT, "Wrong entity type");
    return;
  }
    
  CREATE_TEMP_ARRAY(void*, gentity_handles);

  iGeom_getEntities (igeomInstance, 
                     set_handle,
                     (iGeom_EntityType)gentity_type,
                     TEMP_ARRAY_INOUT(gentity_handles));
  PROCESS_ERROR;

  ASSIGN_ARRAY(gentity_handles);
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.getEntities)
}

/**
 * Return number of gentities of specified dimension in this set, or in
 * whole model.
 * @param set_handle Entity set being queried (if 0, whole model)
 * @param gentity_dimension Dimension of entities being queried
 * @return Number of entities
 */
void
iGeom_SIDL::GeomSidl_impl::getNumOfType (
  /* in */ void* set_handle,
  /* in */ ::iBase::EntityType gentity_type,
  /* out */ int32_t& num_type ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.getNumOfType)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.getNumOfType} (getNumOfType method)
  if (gentity_type > ::iGeom::EntityType_ALL_TYPES || 
      gentity_type < ::iGeom::EntityType_VERTEX) {
    PROCESS_ERROR_MSG(iBase_INVALID_ARGUMENT, "Wrong entity type");
    return -1;
  }
  int32_t retval = iGeom_getNumOfType(igeomInstance, 
                                      set_handle,
                                      (iGeom_EntityType)gentity_type);
  PROCESS_ERROR;
  return retval;
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.getNumOfType)
}

/**
 * Method:  getEntType[]
 */
::iBase::EntityType
iGeom_SIDL::GeomSidl_impl::getEntType (
  /* in */ void* handle ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.getEntType)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.getEntType} (getEntType method)
  iGeom_EntityType gtype;
  iGeom_getEntType (igeomInstance, handle, &gtype );

  PROCESS_ERROR;

  return static_cast< ::iGeom::EntityType >(gtype);
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.getEntType)
}

/**
 *    Returns an integer array of topological dimensions for an input
 *    array of entity handles.
 */
void
iGeom_SIDL::GeomSidl_impl::getArrType (
  /* in */ ::sidl::array<void*> gentity_handles,
  /* in */ int32_t gentity_handles_size,
  /* inout */ ::sidl::array< ::iBase::EntityType>& gtype,
  /* inout */ int32_t& gtype_size ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.getArrType)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.getArrType} (getArrType method)
  CREATE_TEMP_ENUM_ARRAY(iGeom_EntityType, gtype);
  
  iGeom_getArrType (igeomInstance, 
                    TEMP_ARRAY_IN(gentity_handles), 
                    TEMP_ARRAY_INOUT(gtype));
  PROCESS_ERROR;

  ASSIGN_ENUM_ARRAY(gtype);
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.getArrType)
}

/**
 * Get the adjacent entities of a given dimension.
 * @param gentity_handle Entity for which adjacencies are requested
 * @param to_dimension Target dimension of adjacent entities
 * @param adj_gentities List returned with adjacent entities
 */
void
iGeom_SIDL::GeomSidl_impl::getEntAdj (
  /* in */ void* entity_handle,
  /* in */ ::iBase::EntityType to_dimension,
  /* inout */ ::sidl::array<void*>& adj_gentities,
  /* inout */ int32_t& adj_gentities_size ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.getEntAdj)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.getEntAdj} (getEntAdj method)
  CREATE_TEMP_ARRAY(void*, adj_gentities);
  
  iGeom_getEntAdj (igeomInstance, 
                   gentity_handle,
                   static_cast<iGeom_EntityType>(to_dimension),
                   TEMP_ARRAY_INOUT(adj_gentities));
  PROCESS_ERROR;

  ASSIGN_ARRAY(adj_gentities);
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.getEntAdj)
}

/**
 * Get the adjacent entities of a given dimension.
 * @param gentity_handle Entity for which adjacencies are requested
 * @param to_dimension Target dimension of adjacent entities
 * @param adj_gentities List returned with adjacent entities
 */
void
iGeom_SIDL::GeomSidl_impl::getArrAdj (
  /* in */ ::sidl::array<void*> entity_handles,
  /* in */ int32_t entity_handles_size,
  /* in */ ::iBase::EntityType requested_entity_type,
  /* inout */ ::sidl::array<void*>& adj_entity_handles,
  /* inout */ int32_t& adj_entity_handles_size,
  /* inout */ ::sidl::array<int32_t>& offset,
  /* inout */ int32_t& offset_size ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.getArrAdj)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.getArrAdj} (getArrAdj method)
  CREATE_TEMP_ARRAY(iGeom_EntityHandle, adj_entity_handles);
  CREATE_TEMP_ARRAY(int32_t, offset);
  
  iGeom_getArrAdj (igeomInstance, 

//                   TEMP_ARRAY_IN(entity_handles)
                   (entity_handles._get_ior() == NULL ? NULL : entity_handles._get_ior()->d_firstElement),
                   entity_handles_size,

                   static_cast<iGeom_EntityType>(requested_entity_type),

//                   TEMP_ARRAY_INOUT(adj_entity_handles),
                   &adj_entity_handles_temp, 
                   &adj_entity_handles_allocated_size, 
                   &adj_entity_handles_size, 

//                   TEMP_ARRAY_INOUT(offset)
                   &offset_temp, 
                   &offset_allocated_size, 
                   &offset_size 

                   );
  PROCESS_ERROR;

  ASSIGN_ARRAY(adj_entity_handles);
  ASSIGN_ARRAY(offset);
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.getArrAdj)
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
iGeom_SIDL::GeomSidl_impl::getEnt2ndAdj (
  /* in */ void* gentity_handle,
  /* in */ ::iBase::EntityType bridge_dimension,
  /* in */ ::iBase::EntityType to_dimension,
  /* inout */ ::sidl::array<void*>& adjacent_gentities,
  /* inout */ int32_t& adjacent_gentities_size ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.getEnt2ndAdj)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.getEnt2ndAdj} (getEnt2ndAdj method)
  CREATE_TEMP_ARRAY(void*, adjacent_gentities);
  
  iGeom_getEnt2ndAdj (igeomInstance, 
                      gentity_handle,
                      bridge_dimension,
                      to_dimension,
                      TEMP_ARRAY_INOUT(adjacent_gentities));
  PROCESS_ERROR;
  
  ASSIGN_ARRAY(adjacent_gentities);
  
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.getEnt2ndAdj)
}

/**
 * Get the "2nd order" adjacent entities, through a specified "bridge"
 * dimension, of a target dimension.  For example, given a region, return
 * the regions (to_dimension=3) sharing an edge (bridge_dimension=1)
 * with that region.  bridge_dimension must be less than dimension of 
 * gentity_handle, and to_dimension must be greater than bridge dimension.
 * 
 * @param gentity_handle Entity for which 2nd order adjacencies are requested
 * @param order_adjacent_key Dimension of "bridge" entities
 * @param requested_entity_type Target dimension of 2nd order adjacent entities
 * @param adj_entity_handles List returned with 2nd order adjacent entities
 */
void
iGeom_SIDL::GeomSidl_impl::getArr2ndAdj (
  /* in */ ::sidl::array<void*> entity_handles,
  /* in */ int32_t entity_handles_size,
  /* in */ ::iBase::EntityType order_adjacent_key,
  /* in */ ::iBase::EntityType requested_entity_type,
  /* inout */ ::sidl::array<void*>& adj_entity_handles,
  /* inout */ int32_t& adj_entity_handles_size,
  /* inout */ ::sidl::array<int32_t>& offset,
  /* inout */ int32_t& offset_size ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.getArr2ndAdj)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.getArr2ndAdj} (getArr2ndAdj method)
  CREATE_TEMP_ARRAY(void*, adj_entity_handles);
  CREATE_TEMP_ARRAY(int32_t, offset);
  
  iGeom_getArr2ndAdj (igeomInstance, 
                   TEMP_ARRAY_IN(entity_handles),
                   static_cast<iGeom_EntityType>(order_adjacent_key),
                   static_cast<iGeom_EntityType>(requested_entity_type),
                   TEMP_ARRAY_INOUT(adj_entity_handles),
                   TEMP_ARRAY_INOUT(offset));
  PROCESS_ERROR;

  ASSIGN_ARRAY(adj_entity_handles);
  ASSIGN_ARRAY(offset);
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.getArr2ndAdj)
}

/**
 * Return whether or not entities are adjacent.
 * @param gentity_handle1 1st entity
 * @param gentity_handle2 2nd entity
 * @param are_adjacent If true, entities are adjacent
 */
int32_t
iGeom_SIDL::GeomSidl_impl::isEntAdj (
  /* in */ void* gentity_handle1,
  /* in */ void* gentity_handle2,
  /* out */ bool& are_adjacent ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.isEntAdj)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.isEntAdj} (isEntAdj method)
  iGeom_isEntAdj (igeomInstance, 
                  gentity_handle1,
                  gentity_handle2,
                  &are_adjacent);
  PROCESS_ERROR;
  return are_adjacent;
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.isEntAdj)
}

/**
 * Return whether or not arrays of entities are adjacent.
 * @param gentity_handle1 1st entity
 * @param gentity_handle2 2nd entity
 * @param are_adjacent If true, entities are adjacent
 */
void
iGeom_SIDL::GeomSidl_impl::isArrAdj (
  /* in */ ::sidl::array<void*> entity_handles_1,
  /* in */ int32_t entity_handles_1_size,
  /* in */ ::sidl::array<void*> entity_handles_2,
  /* in */ int32_t entity_handles_2_size,
  /* inout */ ::sidl::array<int32_t>& is_adjacent_info,
  /* inout */ int32_t& is_adjacent_info_size ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.isArrAdj)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.isArrAdj} (isArrAdj method)
  CREATE_TEMP_ARRAY(int32_t, is_adjacent_info);
  iGeom_isArrAdj (igeomInstance, 
                  TEMP_ARRAY_IN(entity_handles_1),
                  TEMP_ARRAY_IN(entity_handles_2),
                  TEMP_ARRAY_INOUT(is_adjacent_info));
  PROCESS_ERROR;
  ASSIGN_ARRAY(is_adjacent_info);
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.isArrAdj)
}

/**
 * Method:  getTopoLevel[]
 */
int32_t
iGeom_SIDL::GeomSidl_impl::getTopoLevel (
  /* in */ const ::std::string& model_name ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.getTopoLevel)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.getTopoLevel} (getTopoLevel method)
  int level;
  iGeom_getTopoLevel (igeomInstance, model_name.c_str(), &level);
  PROCESS_ERROR;
  return level;
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.getTopoLevel)
}

/**
 * Method:  getEntClosestPt[]
 */
void
iGeom_SIDL::GeomSidl_impl::getEntClosestPt (
  /* in */ void* entity_handle,
  /* in */ double near_x,
  /* in */ double near_y,
  /* in */ double near_z,
  /* out */ double& on_x,
  /* out */ double& on_y,
  /* out */ double& on_z ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.getEntClosestPt)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.getEntClosestPt} (getEntClosestPt method)
  iGeom_getEntClosestPt (igeomInstance, 
                         entity_handle,
                         near_x, near_y, near_z, 
                         &on_x, &on_y, &on_z);
  PROCESS_ERROR;
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.getEntClosestPt)
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
iGeom_SIDL::GeomSidl_impl::getArrClosestPt (
  /* in */ ::sidl::array<void*> gentity_handles,
  /* in */ int32_t gentity_handles_size,
  /* in */ ::iBase::StorageOrder storage_order,
  /* in */ ::sidl::array<double> near_coordinates,
  /* in */ int32_t near_coordinates_size,
  /* inout */ ::sidl::array<double>& on_coordinates,
  /* inout */ int32_t& on_coordinates_size ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.getArrClosestPt)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.getArrClosestPt} (getArrClosestPt method)
  CREATE_TEMP_ARRAY(double, on_coordinates);
  
  iGeom_getArrClosestPt (igeomInstance, 
                             TEMP_ARRAY_IN(gentity_handles),
                             static_cast<iBase_StorageOrder>(storage_order),
                             TEMP_ARRAY_IN(near_coordinates),
                             TEMP_ARRAY_INOUT(on_coordinates));
  PROCESS_ERROR;

  ASSIGN_ARRAY(on_coordinates);
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.getArrClosestPt)
}

/**
 * Method:  getEntNrmlXYZ[]
 */
void
iGeom_SIDL::GeomSidl_impl::getEntNrmlXYZ (
  /* in */ void* entity_handle,
  /* in */ double x,
  /* in */ double y,
  /* in */ double z,
  /* out */ double& nrml_i,
  /* out */ double& nrml_j,
  /* out */ double& nrml_k ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.getEntNrmlXYZ)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.getEntNrmlXYZ} (getEntNrmlXYZ method)
  iGeom_getEntNrmlXYZ (igeomInstance, 
                       entity_handle,
                       x, y, z, 
                       &nrml_i, &nrml_j, &nrml_k);
  PROCESS_ERROR;
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.getEntNrmlXYZ)
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
iGeom_SIDL::GeomSidl_impl::getArrNrmlXYZ (
  /* in */ ::sidl::array<void*> gentity_handles,
  /* in */ int32_t gentity_handles_size,
  /* in */ ::iBase::StorageOrder storage_order,
  /* in */ ::sidl::array<double> coordinates,
  /* in */ int32_t coordinates_size,
  /* inout */ ::sidl::array<double>& normals,
  /* inout */ int32_t& normals_size ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.getArrNrmlXYZ)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.getArrNrmlXYZ} (getArrNrmlXYZ method)
  CREATE_TEMP_ARRAY(double, normals);
  
  iGeom_getArrNrmlXYZ (igeomInstance, 
                       TEMP_ARRAY_IN(gentity_handles),
                       static_cast<iBase_StorageOrder>(storage_order),
                       TEMP_ARRAY_IN(coordinates),
                       TEMP_ARRAY_INOUT(normals));
  PROCESS_ERROR;

  ASSIGN_ARRAY(normals);
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.getArrNrmlXYZ)
}

/**
 * Method:  getEntNrmlPlXYZ[]
 */
void
iGeom_SIDL::GeomSidl_impl::getEntNrmlPlXYZ (
  /* in */ void* entity_handle,
  /* in */ double x,
  /* in */ double y,
  /* in */ double z,
  /* out */ double& pt_x,
  /* out */ double& pt_y,
  /* out */ double& pt_z,
  /* out */ double& nrml_i,
  /* out */ double& nrml_j,
  /* out */ double& nrml_k ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.getEntNrmlPlXYZ)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.getEntNrmlPlXYZ} (getEntNrmlPlXYZ method)
  iGeom_getEntNrmlPlXYZ (igeomInstance, 
                       entity_handle,
                       x, y, z,
                       &pt_x, &pt_y, &pt_z, 
                       &nrml_i, &nrml_j, &nrml_k);
  PROCESS_ERROR;
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.getEntNrmlPlXYZ)
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
iGeom_SIDL::GeomSidl_impl::getArrNrmlPlXYZ (
  /* in */ ::sidl::array<void*> gentity_handles,
  /* in */ int32_t gentity_handles_size,
  /* in */ ::iBase::StorageOrder storage_order,
  /* in */ ::sidl::array<double> near_coordinates,
  /* in */ int32_t near_coordinates_size,
  /* inout */ ::sidl::array<double>& on_coordinates,
  /* inout */ int32_t& on_coordinates_size,
  /* inout */ ::sidl::array<double>& normals,
  /* inout */ int32_t& normals_size ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.getArrNrmlPlXYZ)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.getArrNrmlPlXYZ} (getArrNrmlPlXYZ method)
  CREATE_TEMP_ARRAY(double, on_coordinates);
  CREATE_TEMP_ARRAY(double, normals);
  
  iGeom_getArrNrmlPlXYZ (igeomInstance, 
                             TEMP_ARRAY_IN(gentity_handles),
                             static_cast<iBase_StorageOrder>(storage_order),
                             TEMP_ARRAY_IN(near_coordinates),
                             TEMP_ARRAY_INOUT(on_coordinates),
                             TEMP_ARRAY_INOUT(normals));
  PROCESS_ERROR;

  ASSIGN_ARRAY(on_coordinates);
  ASSIGN_ARRAY(normals);
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.getArrNrmlPlXYZ)
}

/**
 * Method:  getEntTgntXYZ[]
 */
void
iGeom_SIDL::GeomSidl_impl::getEntTgntXYZ (
  /* in */ void* entity_handle,
  /* in */ double x,
  /* in */ double y,
  /* in */ double z,
  /* out */ double& tgnt_i,
  /* out */ double& tgnt_j,
  /* out */ double& tgnt_k ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.getEntTgntXYZ)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.getEntTgntXYZ} (getEntTgntXYZ method)
  iGeom_getEntTgntXYZ (igeomInstance, 
                       entity_handle,
                       x, y, z, 
                       &tgnt_i, &tgnt_j, &tgnt_k);
  PROCESS_ERROR;
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.getEntTgntXYZ)
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
iGeom_SIDL::GeomSidl_impl::getArrTgntXYZ (
  /* in */ ::sidl::array<void*> gentity_handles,
  /* in */ int32_t gentity_handles_size,
  /* in */ ::iBase::StorageOrder storage_order,
  /* in */ ::sidl::array<double> coordinates,
  /* in */ int32_t coordinates_size,
  /* inout */ ::sidl::array<double>& tangents,
  /* inout */ int32_t& tangents_size ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.getArrTgntXYZ)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.getArrTgntXYZ} (getArrTgntXYZ method)
  CREATE_TEMP_ARRAY(double, tangents);
  
  iGeom_getArrTgntXYZ (igeomInstance, 
                        TEMP_ARRAY_IN(gentity_handles),
                        static_cast<iBase_StorageOrder>(storage_order),
                        TEMP_ARRAY_IN(coordinates),
                        TEMP_ARRAY_INOUT(tangents));
  PROCESS_ERROR;

  ASSIGN_ARRAY(tangents);
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.getArrTgntXYZ)
}

/**
 * Method:  getFcCvtrXYZ[]
 */
void
iGeom_SIDL::GeomSidl_impl::getFcCvtrXYZ (
  /* in */ void* face_handle,
  /* in */ double x,
  /* in */ double y,
  /* in */ double z,
  /* out */ double& cvtr1_i,
  /* out */ double& cvtr1_j,
  /* out */ double& cvtr1_k,
  /* out */ double& cvtr2_i,
  /* out */ double& cvtr2_j,
  /* out */ double& cvtr2_k ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.getFcCvtrXYZ)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.getFcCvtrXYZ} (getFcCvtrXYZ method)
  iGeom_getFcCvtrXYZ (igeomInstance, 
                       face_handle,
                       x, y, z, 
                       &cvtr1_i, &cvtr1_j, &cvtr1_k,
                       &cvtr2_i, &cvtr2_j, &cvtr2_k);
  PROCESS_ERROR;
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.getFcCvtrXYZ)
}

/**
 * Method:  getEgCvtrXYZ[]
 */
void
iGeom_SIDL::GeomSidl_impl::getEgCvtrXYZ (
  /* in */ void* edge_handle,
  /* in */ double x,
  /* in */ double y,
  /* in */ double z,
  /* out */ double& cvtr_i,
  /* out */ double& cvtr_j,
  /* out */ double& cvtr_k ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.getEgCvtrXYZ)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.getEgCvtrXYZ} (getEgCvtrXYZ method)
  iGeom_getEgCvtrXYZ ( igeomInstance, 
                       edge_handle,
                       x, y, z, 
                       &cvtr_i, &cvtr_j, &cvtr_k);
  PROCESS_ERROR;
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.getEgCvtrXYZ)
}

/**
 * Method:  getEntArrCvtrXYZ[]
 */
void
iGeom_SIDL::GeomSidl_impl::getEntArrCvtrXYZ (
  /* in */ ::sidl::array<void*> entity_handles,
  /* in */ int32_t entity_handles_size,
  /* in */ ::iBase::StorageOrder storage_order,
  /* in */ ::sidl::array<double> coords,
  /* in */ int32_t coords_size,
  /* inout */ ::sidl::array<double>& cvtr_1,
  /* inout */ int32_t& cvtr_1_size,
  /* inout */ ::sidl::array<double>& cvtr_2,
  /* inout */ int32_t& cvtr_2_size ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.getEntArrCvtrXYZ)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.getEntArrCvtrXYZ} (getEntArrCvtrXYZ method)
  CREATE_TEMP_ARRAY(double, cvtr_1);
  CREATE_TEMP_ARRAY(double, cvtr_2);
  
  iGeom_getEntArrCvtrXYZ (igeomInstance, 
                        TEMP_ARRAY_IN(entity_handles),
                        static_cast<iBase_StorageOrder>(storage_order),
                        TEMP_ARRAY_IN(coords),
                        TEMP_ARRAY_INOUT(cvtr_1),
                        TEMP_ARRAY_INOUT(cvtr_2));
  PROCESS_ERROR;

  ASSIGN_ARRAY(cvtr_1);
  ASSIGN_ARRAY(cvtr_2);
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.getEntArrCvtrXYZ)
}

/**
 * Method:  getEgEvalXYZ[]
 */
void
iGeom_SIDL::GeomSidl_impl::getEgEvalXYZ (
  /* in */ void* edge_handle,
  /* in */ double x,
  /* in */ double y,
  /* in */ double z,
  /* out */ double& on_x,
  /* out */ double& on_y,
  /* out */ double& on_z,
  /* out */ double& tgnt_i,
  /* out */ double& tgnt_j,
  /* out */ double& tgnt_k,
  /* out */ double& cvtr_i,
  /* out */ double& cvtr_j,
  /* out */ double& cvtr_k ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.getEgEvalXYZ)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.getEgEvalXYZ} (getEgEvalXYZ method)
  iGeom_getEgEvalXYZ (igeomInstance, 
                       edge_handle,
                       x, y, z, 
                       &on_x, &on_y, &on_z,
                       &tgnt_i, &tgnt_j, &tgnt_k,
                       &cvtr_i, &cvtr_j, &cvtr_k );
  PROCESS_ERROR;
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.getEgEvalXYZ)
}

/**
 * Method:  getFcEvalXYZ[]
 */
void
iGeom_SIDL::GeomSidl_impl::getFcEvalXYZ (
  /* in */ void* face_handle,
  /* in */ double x,
  /* in */ double y,
  /* in */ double z,
  /* out */ double& on_x,
  /* out */ double& on_y,
  /* out */ double& on_z,
  /* out */ double& nrml_i,
  /* out */ double& nrml_j,
  /* out */ double& nrml_k,
  /* out */ double& cvtr1_i,
  /* out */ double& cvtr1_j,
  /* out */ double& cvtr1_k,
  /* out */ double& cvtr2_i,
  /* out */ double& cvtr2_j,
  /* out */ double& cvtr2_k ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.getFcEvalXYZ)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.getFcEvalXYZ} (getFcEvalXYZ method)
  iGeom_getFcEvalXYZ (igeomInstance, 
                       face_handle,
                       x, y, z, 
                       &on_x, &on_y, &on_z,
                       &nrml_i, &nrml_j, &nrml_k,
                       &cvtr1_i, &cvtr1_j, &cvtr1_k,
                       &cvtr2_i, &cvtr2_j, &cvtr2_k);
  PROCESS_ERROR;
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.getFcEvalXYZ)
}

/**
 * Method:  getArrEgEvalXYZ[]
 */
void
iGeom_SIDL::GeomSidl_impl::getArrEgEvalXYZ (
  /* in */ ::sidl::array<void*> edge_handles,
  /* in */ int32_t edge_handles_size,
  /* in */ ::iBase::StorageOrder storage_order,
  /* in */ ::sidl::array<double> coords,
  /* in */ int32_t coords_size,
  /* inout */ ::sidl::array<double>& on_coords,
  /* inout */ int32_t& on_coords_size,
  /* inout */ ::sidl::array<double>& tangent,
  /* inout */ int32_t& tangent_size,
  /* inout */ ::sidl::array<double>& cvtr,
  /* inout */ int32_t& cvtr_size ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.getArrEgEvalXYZ)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.getArrEgEvalXYZ} (getArrEgEvalXYZ method)
  CREATE_TEMP_ARRAY(double, on_coords);
  CREATE_TEMP_ARRAY(double, tangent);
  CREATE_TEMP_ARRAY(double, cvtr);
  
  iGeom_getArrEgEvalXYZ (igeomInstance, 
                        TEMP_ARRAY_IN(edge_handles),
                        static_cast<iBase_StorageOrder>(storage_order),
                        TEMP_ARRAY_IN(coords),
                        TEMP_ARRAY_INOUT(on_coords),
                        TEMP_ARRAY_INOUT(tangent),
                        TEMP_ARRAY_INOUT(cvtr));
  PROCESS_ERROR;

  ASSIGN_ARRAY(on_coords);
  ASSIGN_ARRAY(tangent);
  ASSIGN_ARRAY(cvtr);
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.getArrEgEvalXYZ)
}

/**
 * Method:  getArrFcEvalXYZ[]
 */
void
iGeom_SIDL::GeomSidl_impl::getArrFcEvalXYZ (
  /* in */ ::sidl::array<void*> face_handles,
  /* in */ int32_t face_handles_size,
  /* in */ ::iBase::StorageOrder storage_order,
  /* in */ ::sidl::array<double> coords,
  /* in */ int32_t coords_size,
  /* inout */ ::sidl::array<double>& on_coords,
  /* inout */ int32_t& on_coords_size,
  /* inout */ ::sidl::array<double>& normal,
  /* inout */ int32_t& normal_size,
  /* inout */ ::sidl::array<double>& cvtr_1,
  /* inout */ int32_t& cvtr_1_size,
  /* inout */ ::sidl::array<double>& cvtr_2,
  /* inout */ int32_t& cvtr_2_size ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.getArrFcEvalXYZ)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.getArrFcEvalXYZ} (getArrFcEvalXYZ method)
  CREATE_TEMP_ARRAY(double, on_coords);
  CREATE_TEMP_ARRAY(double, normal);
  CREATE_TEMP_ARRAY(double, cvtr_1);
  CREATE_TEMP_ARRAY(double, cvtr_2);
  
  iGeom_getArrFcEvalXYZ (igeomInstance, 
                        TEMP_ARRAY_IN(face_handles),
                        static_cast<iBase_StorageOrder>(storage_order),
                        TEMP_ARRAY_IN(coords),
                        TEMP_ARRAY_INOUT(on_coords),
                        TEMP_ARRAY_INOUT(normal),
                        TEMP_ARRAY_INOUT(cvtr_1),
                        TEMP_ARRAY_INOUT(cvtr_2));
  PROCESS_ERROR;

  ASSIGN_ARRAY(on_coords);
  ASSIGN_ARRAY(normal);
  ASSIGN_ARRAY(cvtr_1);
  ASSIGN_ARRAY(cvtr_2);
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.getArrFcEvalXYZ)
}

/**
 * Method:  getEntBoundBox[]
 */
void
iGeom_SIDL::GeomSidl_impl::getEntBoundBox (
  /* in */ void* entity_handle,
  /* out */ double& min_x,
  /* out */ double& min_y,
  /* out */ double& min_z,
  /* out */ double& max_x,
  /* out */ double& max_y,
  /* out */ double& max_z ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.getEntBoundBox)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.getEntBoundBox} (getEntBoundBox method)
  iGeom_getEntBoundBox (igeomInstance, 
                       entity_handle,
                       &min_x, &min_y, &min_z,
                       &max_x, &max_y, &max_z);
  PROCESS_ERROR;
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.getEntBoundBox)
}

/**
 * Return the bounding boxex of given entities; coordinates returned
 * interleaved.
 * @param gentity_handles The gentities being queried
 * @param min_corners Minimum corner coordinates of the boxes, interleaved
 * @param max_corners Maximum corner coordinates of the boxes, interleaved
 */
void
iGeom_SIDL::GeomSidl_impl::getArrBoundBox (
  /* in */ ::sidl::array<void*> gentity_handles,
  /* in */ int32_t gentity_handles_size,
  /* inout */ ::iBase::StorageOrder& storage_order,
  /* inout */ ::sidl::array<double>& min_corner,
  /* inout */ int32_t& min_corner_size,
  /* inout */ ::sidl::array<double>& max_corner,
  /* inout */ int32_t& max_corner_size ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.getArrBoundBox)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.getArrBoundBox} (getArrBoundBox method)
  CREATE_TEMP_ARRAY(double, max_corner);
  CREATE_TEMP_ARRAY(double, min_corner);
  
  iGeom_getArrBoundBox (igeomInstance, 
                            TEMP_ARRAY_IN(gentity_handles),
                            static_cast<iBase_StorageOrder>(storage_order),
                            TEMP_ARRAY_INOUT(min_corner),
                            TEMP_ARRAY_INOUT(max_corner));
  PROCESS_ERROR;

  ASSIGN_ARRAY(min_corner);
  ASSIGN_ARRAY(max_corner);
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.getArrBoundBox)
}

/**
 * Method:  getVtxCoord[]
 */
void
iGeom_SIDL::GeomSidl_impl::getVtxCoord (
  /* in */ void* vertex_handle,
  /* out */ double& x,
  /* out */ double& y,
  /* out */ double& z ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.getVtxCoord)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.getVtxCoord} (getVtxCoord method)
  
  iGeom_getVtxCoord (igeomInstance, vertex_handle, &x, &y, &z);
  PROCESS_ERROR;
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.getVtxCoord)
}

/**
 * Return the coordinates of the specified vertices; returns error if any
 * of the entities are not gvertices.  Coordinates returned interleaved.
 * @param gentity_handles The gentities being queried
 * @param coordinates The coordinates of the gvertices, interleaved.
 */
void
iGeom_SIDL::GeomSidl_impl::getVtxArrCoords (
  /* in */ ::sidl::array<void*> gentity_handles,
  /* in */ int32_t gentity_handles_size,
  /* inout */ ::iBase::StorageOrder& storage_order,
  /* inout */ ::sidl::array<double>& coordinates,
  /* inout */ int32_t& coordinates_size ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.getVtxArrCoords)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.getVtxArrCoords} (getVtxArrCoords method)
  CREATE_TEMP_ARRAY(double, coordinates);
  
  iGeom_getVtxArrCoords (igeomInstance, 
                               TEMP_ARRAY_IN(gentity_handles),
                               static_cast<iBase_StorageOrder>(storage_order),
                               TEMP_ARRAY_INOUT(coordinates));
  PROCESS_ERROR;

  ASSIGN_ARRAY(coordinates);
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.getVtxArrCoords)
}

/**
 * Method:  getPntIntsct[]
 */
void
iGeom_SIDL::GeomSidl_impl::getPntIntsct (
  /* in */ double x,
  /* in */ double y,
  /* in */ double z,
  /* in */ double dir_x,
  /* in */ double dir_y,
  /* in */ double dir_z,
  /* inout */ ::sidl::array<void*>& intersect_entity_handles,
  /* inout */ int32_t& intersect_entity_handles_size,
  /* inout */ ::iBase::StorageOrder& storage_order,
  /* inout */ ::sidl::array<double>& intersect_coords,
  /* inout */ int32_t& intersect_coords_size,
  /* inout */ ::sidl::array<double>& param_coords,
  /* inout */ int32_t& param_coords_size ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.getPntIntsct)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.getPntIntsct} (getPntIntsct method)
  CREATE_TEMP_ARRAY(void*, intersect_entity_handles);
  CREATE_TEMP_ARRAY(double, intersect_coords);
  CREATE_TEMP_ARRAY(double, param_coords);
  
  iBase_StorageOrder storage_order_tmp;
  iGeom_getPntIntsct (igeomInstance, x, y, z, dir_x, dir_y, dir_z,
                      TEMP_ARRAY_INOUT( intersect_entity_handles ),
                      &storage_order_tmp,
                      TEMP_ARRAY_INOUT( intersect_coords ),
                      TEMP_ARRAY_INOUT( param_coords ) );
  PROCESS_ERROR;

  storage_order = static_cast< ::iBase::StorageOrder >(storage_order_tmp);
  ASSIGN_ARRAY(intersect_entity_handles);
  ASSIGN_ARRAY(intersect_coords);
  ASSIGN_ARRAY(param_coords);
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.getPntIntsct)
}

/**
 * Method:  getPntArrRayIntsct[]
 */
void
iGeom_SIDL::GeomSidl_impl::getPntArrRayIntsct (
  /* in */ ::iBase::StorageOrder storage_order,
  /* in */ ::sidl::array<double> coords,
  /* in */ int32_t coords_size,
  /* in */ ::sidl::array<double> directions,
  /* in */ int32_t directions_size,
  /* inout */ ::sidl::array<void*>& intersect_entity_handles,
  /* inout */ int32_t& intersect_entity_handles_size,
  /* inout */ ::sidl::array<int32_t>& offset,
  /* inout */ int32_t& offset_size,
  /* inout */ ::sidl::array<double>& intersect_coords,
  /* inout */ int32_t& intersect_coords_size,
  /* inout */ ::sidl::array<double>& param_coords,
  /* inout */ int32_t& param_coords_size ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.getPntArrRayIntsct)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.getPntArrRayIntsct} (getPntArrRayIntsct method)
  CREATE_TEMP_ARRAY(void*, intersect_entity_handles);
  CREATE_TEMP_ARRAY(int32_t, offset);
  CREATE_TEMP_ARRAY(double, intersect_coords);
  CREATE_TEMP_ARRAY(double, param_coords);
  
  iGeom_getPntArrRayIntsct (igeomInstance, 
                      static_cast<iBase_StorageOrder>(storage_order),
                      TEMP_ARRAY_IN( coords ),
                      TEMP_ARRAY_IN( directions ),
                      TEMP_ARRAY_INOUT( intersect_entity_handles ),
                      TEMP_ARRAY_INOUT( offset ),
                      TEMP_ARRAY_INOUT( intersect_coords ),
                      TEMP_ARRAY_INOUT( param_coords ) );
  PROCESS_ERROR;

  ASSIGN_ARRAY(intersect_entity_handles);
  ASSIGN_ARRAY(offset);
  ASSIGN_ARRAY(intersect_coords);
  ASSIGN_ARRAY(param_coords);
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.getPntArrRayIntsct)
}

/**
 * Return the sense of a gface with respect to a gregion.  Sense is either
 * forward (=1), reverse (=-1), both (=2), or unknown (=0).  Error is returned
 * if first entity is not a gface or second entity is not a gregion.
 * @param gface Gface whose sense is being queried.
 * @param gregion Gregion gface is being queried with respect to
 */
int32_t
iGeom_SIDL::GeomSidl_impl::getEntNrmlSense (
  /* in */ void* gface,
  /* in */ void* gregion ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.getEntNrmlSense)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.getEntNrmlSense} (getEntNrmlSense method)
  int32_t retval = iGeom_getEntNrmlSense (igeomInstance, 
                                          gface,
                                          gregion);
  PROCESS_ERROR;
  return retval;
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.getEntNrmlSense)
}

/**
 * Method:  getArrNrmlSense[]
 */
void
iGeom_SIDL::GeomSidl_impl::getArrNrmlSense (
  /* in */ ::sidl::array<void*> face_handles,
  /* in */ int32_t face_handles_size,
  /* in */ ::sidl::array<void*> region_handles,
  /* in */ int32_t region_handles_size,
  /* inout */ ::sidl::array<int32_t>& sense,
  /* inout */ int32_t& sense_size ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.getArrNrmlSense)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.getArrNrmlSense} (getArrNrmlSense method)
  CREATE_TEMP_ARRAY(int32_t, sense);
  iGeom_getArrNrmlSense (igeomInstance,
                       TEMP_ARRAY_IN(face_handles), 
                       TEMP_ARRAY_IN(region_handles), 
                       TEMP_ARRAY_INOUT(sense));
  PROCESS_ERROR;
  ASSIGN_ARRAY(sense);
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.getArrNrmlSense)
}

/**
 * Return the sense of a gedge with respect to a gface.  Sense is either
 * forward (=1), reverse (=-1), both (=2), or unknown (=0).  Error is returned
 * if first entity is not a gedge or second entity is not a gface.
 * @param gedge Gedge whose sense is being queried.
 * @param gface Gface gedge is being queried with respect to
 */
int32_t
iGeom_SIDL::GeomSidl_impl::getEgFcSense (
  /* in */ void* gedge,
  /* in */ void* gface ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.getEgFcSense)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.getEgFcSense} (getEgFcSense method)
  int32_t retval = iGeom_getEgFcSense (igeomInstance, 
                                           gedge,
                                           gface);
  PROCESS_ERROR;
  return retval;
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.getEgFcSense)
}

/**
 * Method:  getEgFcArrSense[]
 */
void
iGeom_SIDL::GeomSidl_impl::getEgFcArrSense (
  /* in */ ::sidl::array<void*> edge_handles,
  /* in */ int32_t edge_handles_size,
  /* in */ ::sidl::array<void*> face_handles,
  /* in */ int32_t face_handles_size,
  /* inout */ ::sidl::array<int32_t>& sense,
  /* inout */ int32_t& sense_size ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.getEgFcArrSense)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.getEgFcArrSense} (getEgFcArrSense method)
  CREATE_TEMP_ARRAY(int32_t, sense);
  iGeom_getEgFcArrSense (igeomInstance,
                       TEMP_ARRAY_IN(edge_handles), 
                       TEMP_ARRAY_IN(face_handles), 
                       TEMP_ARRAY_INOUT(sense));
  PROCESS_ERROR;
  ASSIGN_ARRAY(sense);
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.getEgFcArrSense)
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
iGeom_SIDL::GeomSidl_impl::getEgVtxSense (
  /* in */ void* gedge,
  /* in */ void* gvertex1,
  /* in */ void* gvertex2 ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.getEgVtxSense)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.getEgVtxSense} (getEgVtxSense method)
  int32_t retval = iGeom_getEgVtxSense (igeomInstance, 
                                                 gedge,
                                                 gvertex1,
                                                 gvertex2);
  PROCESS_ERROR;
  return retval;
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.getEgVtxSense)
}

/**
 * Method:  getEgVtxArrSense[]
 */
void
iGeom_SIDL::GeomSidl_impl::getEgVtxArrSense (
  /* in */ ::sidl::array<void*> edge_handles,
  /* in */ int32_t edge_handles_size,
  /* in */ ::sidl::array<void*> vertex_handles_1,
  /* in */ int32_t vertex_handles_1_size,
  /* in */ ::sidl::array<void*> vertex_handles_2,
  /* in */ int32_t vertex_handles_2_size,
  /* inout */ ::sidl::array<int32_t>& sense,
  /* inout */ int32_t& sense_size ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.getEgVtxArrSense)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.getEgVtxArrSense} (getEgVtxArrSense method)
  CREATE_TEMP_ARRAY(int32_t, sense);
  iGeom_getEgVtxArrSense (igeomInstance,
                       TEMP_ARRAY_IN(edge_handles), 
                       TEMP_ARRAY_IN(vertex_handles_1), 
                       TEMP_ARRAY_IN(vertex_handles_2), 
                       TEMP_ARRAY_INOUT(sense));
  PROCESS_ERROR;
  ASSIGN_ARRAY(sense);
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.getEgVtxArrSense)
}

/**
 * Return the arc length / area / volume of the entities
 * @param gentity_handles Entities for which measure is requested
 * @param gentity_handles_size Number of gentities
 * @param measures Arc length / area / volume of the entities
 * @param measures_length Number of entries in measures
 */
void
iGeom_SIDL::GeomSidl_impl::measure (
  /* in */ ::sidl::array<void*> gentity_handles,
  /* in */ int32_t gentity_handles_size,
  /* inout */ ::sidl::array<double>& measures,
  /* inout */ int32_t& measures_size ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.measure)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.measure} (measure method)
  CREATE_TEMP_ARRAY(int32_t, sense);
  iGeom_getEgVtxArrSense (igeomInstance,
                       TEMP_ARRAY_IN(edge_handles), 
                       TEMP_ARRAY_IN(vertex_handles_1), 
                       TEMP_ARRAY_IN(vertex_handles_2), 
                       TEMP_ARRAY_INOUT(sense));
  PROCESS_ERROR;
  ASSIGN_ARRAY(sense);
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.measure)
}

/**
 * Return the type of surface as a string; if not a surface, an error is returned
 * @param face_handle Face for which the type is requested
 * @param face_type Type of face, returned as a string
 */
void
iGeom_SIDL::GeomSidl_impl::getFaceType (
  /* in */ void* gface_handle,
  /* inout */ ::std::string& face_type ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.getFaceType)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.getFaceType} (getFaceType method)
  
  iGeom_measure(igeomInstance, 
                       TEMP_ARRAY_IN(gentity_handles),
                       TEMP_ARRAY_INOUT(measures));
  PROCESS_ERROR;

  ASSIGN_ARRAY(measures);
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.getFaceType)
}

/**
 * Method:  getParametric[]
 */
int32_t
iGeom_SIDL::GeomSidl_impl::getParametric ()
throw ( 
  ::iBase::Error
)
{
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.getParametric)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.getParametric} (getParametric method)
  CREATE_TEMP_ARRAY(int32_t, sense);
  iGeom_getEgVtxArrSense (igeomInstance,
                       TEMP_ARRAY_IN(edge_handles), 
                       TEMP_ARRAY_IN(vertex_handles_1), 
                       TEMP_ARRAY_IN(vertex_handles_2), 
                       TEMP_ARRAY_INOUT(sense));
  PROCESS_ERROR;
  ASSIGN_ARRAY(sense);
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.getParametric)
}

/**
 * Return whether a given gentity is parametric or not.  If a gentity
 * is not parametric, all of the following functions will return an error
 * when called on that entity.
 * @param gentity_handle Gentity being queried.
 */
int32_t
iGeom_SIDL::GeomSidl_impl::isEntParametric (
  /* in */ void* gentity_handle ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.isEntParametric)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.isEntParametric} (isEntParametric method)
  iGeom_getFaceType(igeomInstance, gface_handle, &this_type);
  if (NULL != this_type) face_type = this_type;
  PROCESS_ERROR;
    // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.isEntParametric)
}

/**
 * Method:  isArrParametric[]
 */
void
iGeom_SIDL::GeomSidl_impl::isArrParametric (
  /* in */ ::sidl::array<void*> entity_handles,
  /* in */ int32_t entity_handles_size,
  /* inout */ ::sidl::array<int32_t>& is_parametric,
  /* inout */ int32_t& is_parametric_size ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.isArrParametric)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.isArrParametric} (isArrParametric method)
  CREATE_TEMP_ARRAY(int32_t, is_parametric);
  
  iGeom_isArrParametric (igeomInstance, 
                         TEMP_ARRAY_IN(entity_handles),
                         TEMP_ARRAY_INOUT(is_parametric));
  PROCESS_ERROR;
  ASSIGN_ARRAY(is_parametric);
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.isArrParametric)
}

/**
 * Method:  getEntUVtoXYZ[]
 */
void
iGeom_SIDL::GeomSidl_impl::getEntUVtoXYZ (
  /* in */ void* entity_handle,
  /* in */ double u,
  /* in */ double v,
  /* out */ double& x,
  /* out */ double& y,
  /* out */ double& z ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.getEntUVtoXYZ)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.getEntUVtoXYZ} (getEntUVtoXYZ method)
  iGeom_getEntUVtoXYZ (igeomInstance, entity_handle, u, v, &x, &y, &z );
  PROCESS_ERROR;
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.getEntUVtoXYZ)
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
iGeom_SIDL::GeomSidl_impl::getArrUVtoXYZ (
  /* in */ ::sidl::array<void*> gentity_handles,
  /* in */ int32_t gentity_handles_size,
  /* in */ ::iBase::StorageOrder storage_order,
  /* in */ ::sidl::array<double> uv,
  /* in */ int32_t uv_size,
  /* inout */ ::sidl::array<double>& coordinates,
  /* inout */ int32_t& coordinates_size ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.getArrUVtoXYZ)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.getArrUVtoXYZ} (getArrUVtoXYZ method)
  CREATE_TEMP_ARRAY(double, coordinates);
  
  iGeom_getArrUVtoXYZ (igeomInstance, 
                        TEMP_ARRAY_IN(gentity_handles),
                        static_cast<iBase_StorageOrder>(storage_order),
                        TEMP_ARRAY_IN(uv),
                        TEMP_ARRAY_INOUT(coordinates));
  PROCESS_ERROR;

  ASSIGN_ARRAY(coordinates);
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.getArrUVtoXYZ)
}

/**
 * Method:  getEntUtoXYZ[]
 */
void
iGeom_SIDL::GeomSidl_impl::getEntUtoXYZ (
  /* in */ void* entity_handle,
  /* in */ double u,
  /* out */ double& x,
  /* out */ double& y,
  /* out */ double& z ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.getEntUtoXYZ)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.getEntUtoXYZ} (getEntUtoXYZ method)
  iGeom_getEntUtoXYZ (igeomInstance, entity_handle, u, &x, &y, &z );
  PROCESS_ERROR;
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.getEntUtoXYZ)
}

/**
 * Method:  getArrUtoXYZ[]
 */
void
iGeom_SIDL::GeomSidl_impl::getArrUtoXYZ (
  /* in */ ::sidl::array<void*> entity_handles,
  /* in */ int32_t entity_handles_size,
  /* in */ ::sidl::array<double> u,
  /* in */ int32_t u_size,
  /* inout */ ::iBase::StorageOrder& storage_order,
  /* inout */ ::sidl::array<double>& on_coords,
  /* inout */ int32_t& on_coords_size ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.getArrUtoXYZ)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.getArrUtoXYZ} (getArrUtoXYZ method)
  CREATE_TEMP_ARRAY(double, on_coords);
  iGeom_getArrUtoXYZ (igeomInstance, 
                        TEMP_ARRAY_IN(entity_handles),
                        static_cast<iBase_StorageOrder>(storage_order),
                        TEMP_ARRAY_IN(u),
                        TEMP_ARRAY_INOUT(on_coords));
  PROCESS_ERROR;
  ASSIGN_ARRAY(on_coords);
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.getArrUtoXYZ)
}

/**
 * Method:  getEntXYZtoUV[]
 */
void
iGeom_SIDL::GeomSidl_impl::getEntXYZtoUV (
  /* in */ void* entity_handle,
  /* in */ double x,
  /* in */ double y,
  /* in */ double z,
  /* out */ double& u,
  /* out */ double& v ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.getEntXYZtoUV)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.getEntXYZtoUV} (getEntXYZtoUV method)
  iGeom_getEntXYZtoUV (igeomInstance, entity_handle, x, y, z, &u, &v );
  PROCESS_ERROR;
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.getEntXYZtoUV)
}

/**
 * Method:  getEntXYZtoU[]
 */
void
iGeom_SIDL::GeomSidl_impl::getEntXYZtoU (
  /* in */ void* entity_handle,
  /* in */ double x,
  /* in */ double y,
  /* in */ double z,
  /* out */ double& u ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.getEntXYZtoU)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.getEntXYZtoU} (getEntXYZtoU method)
  iGeom_getEntXYZtoU (igeomInstance, entity_handle, x, y, z, &u);
  PROCESS_ERROR;
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.getEntXYZtoU)
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
iGeom_SIDL::GeomSidl_impl::getArrXYZtoUV (
  /* in */ ::sidl::array<void*> gentity_handles,
  /* in */ int32_t gentity_handles_size,
  /* in */ ::iBase::StorageOrder storage_order,
  /* in */ ::sidl::array<double> coordinates,
  /* in */ int32_t coordinates_size,
  /* inout */ ::sidl::array<double>& uv,
  /* inout */ int32_t& uv_size ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.getArrXYZtoUV)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.getArrXYZtoUV} (getArrXYZtoUV method)
  CREATE_TEMP_ARRAY(double, uv);
  
  iGeom_getArrXYZtoUV (igeomInstance, 
                        TEMP_ARRAY_IN(gentity_handles),
                        static_cast<iBase_StorageOrder>(storage_order),
                        TEMP_ARRAY_IN(coordinates),
                        TEMP_ARRAY_INOUT(uv));
  PROCESS_ERROR;

  ASSIGN_ARRAY(uv);
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.getArrXYZtoUV)
}

/**
 * Method:  getArrXYZtoU[]
 */
void
iGeom_SIDL::GeomSidl_impl::getArrXYZtoU (
  /* in */ ::sidl::array<void*> gentity_handles,
  /* in */ int32_t gentity_handles_size,
  /* in */ ::iBase::StorageOrder storage_order,
  /* in */ ::sidl::array<double> coordinates,
  /* in */ int32_t coordinates_size,
  /* inout */ ::sidl::array<double>& u,
  /* inout */ int32_t& u_size ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.getArrXYZtoU)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.getArrXYZtoU} (getArrXYZtoU method)
  CREATE_TEMP_ARRAY(double, u);
  
  iGeom_getArrXYZtoU (igeomInstance, 
                        TEMP_ARRAY_IN(gentity_handles),
                        static_cast<iBase_StorageOrder>(storage_order),
                        TEMP_ARRAY_IN(coordinates),
                        TEMP_ARRAY_INOUT(u));
  PROCESS_ERROR;

  ASSIGN_ARRAY(u);
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.getArrXYZtoU)
}

/**
 * Method:  getEntXYZtoUVHint[]
 */
void
iGeom_SIDL::GeomSidl_impl::getEntXYZtoUVHint (
  /* in */ void* entity_handle,
  /* in */ double x,
  /* in */ double y,
  /* in */ double z,
  /* inout */ double& u,
  /* inout */ double& v ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.getEntXYZtoUVHint)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.getEntXYZtoUVHint} (getEntXYZtoUVHint method)
  iGeom_getEntXYZtoUVHint (igeomInstance, entity_handle, x, y, z, &u, &v );
  PROCESS_ERROR;
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.getEntXYZtoUVHint)
}

/**
 * Method:  getArrXYZtoUVHint[]
 */
void
iGeom_SIDL::GeomSidl_impl::getArrXYZtoUVHint (
  /* in */ ::sidl::array<void*> entity_handles,
  /* in */ int32_t entity_handles_size,
  /* in */ ::iBase::StorageOrder storage_order,
  /* in */ ::sidl::array<double> coords,
  /* in */ int32_t coords_size,
  /* inout */ ::sidl::array<double>& uv,
  /* inout */ int32_t& uv_size ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.getArrXYZtoUVHint)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.getArrXYZtoUVHint} (getArrXYZtoUVHint method)
  CREATE_TEMP_ARRAY(double, uv);
  
  iGeom_getArrXYZtoUVHint (igeomInstance, 
                        TEMP_ARRAY_IN(entity_handles),
                        static_cast<iBase_StorageOrder>(storage_order),
                        TEMP_ARRAY_IN(coords),
                        TEMP_ARRAY_INOUT(uv));
  PROCESS_ERROR;

  ASSIGN_ARRAY(uv);
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.getArrXYZtoUVHint)
}

/**
 * Method:  getEntUVRange[]
 */
void
iGeom_SIDL::GeomSidl_impl::getEntUVRange (
  /* in */ void* entity_handle,
  /* out */ double& u_min,
  /* out */ double& v_min,
  /* out */ double& u_max,
  /* out */ double& v_max ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.getEntUVRange)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.getEntUVRange} (getEntUVRange method)
  iGeom_getEntUVRange (igeomInstance, 
                       entity_handle,
                       &u_min, &v_min, &u_max, &v_max);
  PROCESS_ERROR;
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.getEntUVRange)
}

/**
 * Method:  getEntURange[]
 */
void
iGeom_SIDL::GeomSidl_impl::getEntURange (
  /* in */ void* entity_handle,
  /* out */ double& u_min,
  /* out */ double& u_max ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.getEntURange)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.getEntURange} (getEntURange method)
  iGeom_getEntURange (igeomInstance, 
                       entity_handle,
                       &u_min, &u_max);
  PROCESS_ERROR;
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.getEntURange)
}

/**
 * Return the uv range of the specified gentities.  Parameters are interleaved.
 * @param gentity_handles Gentities being queried.
 * @param uv_min Minimum parameters of gentities, interleaved
 * @param uv_max Maximum parameters of gentities, interleaved
 */
void
iGeom_SIDL::GeomSidl_impl::getArrUVRange (
  /* in */ ::sidl::array<void*> gentity_handles,
  /* in */ int32_t gentity_handles_size,
  /* inout */ ::iBase::StorageOrder& storage_order,
  /* inout */ ::sidl::array<double>& uv_min,
  /* inout */ int32_t& uv_min_size,
  /* inout */ ::sidl::array<double>& uv_max,
  /* inout */ int32_t& uv_max_size ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.getArrUVRange)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.getArrUVRange} (getArrUVRange method)
  CREATE_TEMP_ARRAY(double, uv_max);
  CREATE_TEMP_ARRAY(double, uv_min);
  
  iGeom_getArrUVRange (igeomInstance, 
                        TEMP_ARRAY_IN(gentity_handles),
                        static_cast<iBase_StorageOrder>(storage_order),
                        TEMP_ARRAY_INOUT(uv_min),
                        TEMP_ARRAY_INOUT(uv_max));
  PROCESS_ERROR;

  ASSIGN_ARRAY(uv_max);
  ASSIGN_ARRAY(uv_min);
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.getArrUVRange)
}

/**
 * Method:  getArrURange[]
 */
void
iGeom_SIDL::GeomSidl_impl::getArrURange (
  /* in */ ::sidl::array<void*> entity_handles,
  /* in */ int32_t entity_handles_size,
  /* inout */ ::sidl::array<double>& u_min,
  /* inout */ int32_t& u_min_size,
  /* inout */ ::sidl::array<double>& u_max,
  /* inout */ int32_t& u_max_size ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.getArrURange)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.getArrURange} (getArrURange method)
  CREATE_TEMP_ARRAY(double, u_max);
  CREATE_TEMP_ARRAY(double, u_min);
  
  iGeom_getArrURange (igeomInstance, 
                        TEMP_ARRAY_IN(entity_handles),
                         TEMP_ARRAY_INOUT(u_min),
                        TEMP_ARRAY_INOUT(u_max));
  PROCESS_ERROR;

  ASSIGN_ARRAY(u_max);
  ASSIGN_ARRAY(u_min);
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.getArrURange)
}

/**
 * Method:  getEntUtoUV[]
 */
void
iGeom_SIDL::GeomSidl_impl::getEntUtoUV (
  /* in */ void* edge_handle,
  /* in */ void* face_handle,
  /* in */ double in_u,
  /* out */ double& u,
  /* out */ double& v ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.getEntUtoUV)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.getEntUtoUV} (getEntUtoUV method)
  iGeom_getEntUtoUV (igeomInstance, 
                     edge_handle,
                     face_handle,
                     in_u, &u, &v);
  PROCESS_ERROR;
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.getEntUtoUV)
}

/**
 * Method:  getVtxToUV[]
 */
void
iGeom_SIDL::GeomSidl_impl::getVtxToUV (
  /* in */ void* vertex_handle,
  /* in */ void* face_handle,
  /* out */ double& u,
  /* out */ double& v ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.getVtxToUV)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.getVtxToUV} (getVtxToUV method)
  iGeom_getVtxToUV (igeomInstance, 
                     vertex_handle,
                     face_handle,
                     &u, &v);
  PROCESS_ERROR;
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.getVtxToUV)
}

/**
 * Method:  getVtxToU[]
 */
void
iGeom_SIDL::GeomSidl_impl::getVtxToU (
  /* in */ void* vertex_handle,
  /* in */ void* edge_handle,
  /* out */ double& u ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.getVtxToU)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.getVtxToU} (getVtxToU method)
  iGeom_getVtxToU (igeomInstance, 
                     vertex_handle,
                     edge_handle,
                     &u);
  PROCESS_ERROR;
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.getVtxToU)
}

/**
 * Method:  getArrUtoUV[]
 */
void
iGeom_SIDL::GeomSidl_impl::getArrUtoUV (
  /* in */ ::sidl::array<void*> edge_handles,
  /* in */ int32_t edge_handles_size,
  /* in */ ::sidl::array<void*> face_handles,
  /* in */ int32_t face_handles_size,
  /* in */ ::sidl::array<double> u_in,
  /* in */ int32_t u_in_size,
  /* inout */ ::iBase::StorageOrder& storage_order,
  /* inout */ ::sidl::array<double>& uv,
  /* inout */ int32_t& uv_size ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.getArrUtoUV)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.getArrUtoUV} (getArrUtoUV method)
  CREATE_TEMP_ARRAY(double, uv);
  iBase_StorageOrder storage_order_tmp;
  iGeom_getArrUtoUV (igeomInstance, 
                     TEMP_ARRAY_IN(edge_handles),
                     TEMP_ARRAY_IN(face_handles),
                     TEMP_ARRAY_IN(u_in),
                     &storage_order_tmp,
                     TEMP_ARRAY_INOUT(uv));
  PROCESS_ERROR;
  storage_order = static_cast< ::iBase::StorageOrder >(storage_order_tmp);
  ASSIGN_ARRAY(uv);
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.getArrUtoUV)
}

/**
 * Method:  getVtxArrToUV[]
 */
void
iGeom_SIDL::GeomSidl_impl::getVtxArrToUV (
  /* in */ ::sidl::array<void*> vertex_handles,
  /* in */ int32_t vertex_handles_size,
  /* in */ ::sidl::array<void*> face_handles,
  /* in */ int32_t face_handles_size,
  /* inout */ ::iBase::StorageOrder& storage_order,
  /* inout */ ::sidl::array<double>& uv,
  /* inout */ int32_t& uv_size ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.getVtxArrToUV)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.getVtxArrToUV} (getVtxArrToUV method)
  CREATE_TEMP_ARRAY(double, uv);
  iBase_StorageOrder storage_order_tmp;
  iGeom_getVtxArrToUV (igeomInstance, 
                       TEMP_ARRAY_IN(vertex_handles),
                       TEMP_ARRAY_IN(face_handles),
                       &storage_order_tmp,
                       TEMP_ARRAY_INOUT(uv));
  PROCESS_ERROR;
  storage_order = static_cast< ::iBase::StorageOrder >(storage_order_tmp);
  ASSIGN_ARRAY(uv);
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.getVtxArrToUV)
}

/**
 * Method:  getVtxArrToU[]
 */
void
iGeom_SIDL::GeomSidl_impl::getVtxArrToU (
  /* in */ ::sidl::array<void*> vertex_handles,
  /* in */ int32_t vertex_handles_size,
  /* in */ ::sidl::array<void*> edge_handles,
  /* in */ int32_t edge_handles_size,
  /* inout */ ::sidl::array<double>& u,
  /* inout */ int32_t& u_size ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.getVtxArrToU)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.getVtxArrToU} (getVtxArrToU method)
  CREATE_TEMP_ARRAY(double, u);
  iGeom_getVtxArrToU (igeomInstance, 
                       TEMP_ARRAY_IN(vertex_handles),
                       TEMP_ARRAY_IN(edge_handles),
                       TEMP_ARRAY_INOUT(u));
  PROCESS_ERROR;
  ASSIGN_ARRAY(u);
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.getVtxArrToU)
}

/**
 * Method:  getEntNrmlUV[]
 */
void
iGeom_SIDL::GeomSidl_impl::getEntNrmlUV (
  /* in */ void* entity_handle,
  /* in */ double u,
  /* in */ double v,
  /* out */ double& nrml_i,
  /* out */ double& nrml_j,
  /* out */ double& nrml_k ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.getEntNrmlUV)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.getEntNrmlUV} (getEntNrmlUV method)
  iGeom_getEntNrmlUV (igeomInstance, entity_handle, u, v, &nrml_i, &nrml_j, &nrml_k );
  PROCESS_ERROR;
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.getEntNrmlUV)
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
iGeom_SIDL::GeomSidl_impl::getArrNrmlUV (
  /* in */ ::sidl::array<void*> gface_handles,
  /* in */ int32_t gface_handles_size,
  /* in */ ::iBase::StorageOrder storage_order,
  /* in */ ::sidl::array<double> parameters,
  /* in */ int32_t parameters_size,
  /* inout */ ::sidl::array<double>& normals,
  /* in */ int32_t normals_size ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.getArrNrmlUV)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.getArrNrmlUV} (getArrNrmlUV method)
  CREATE_TEMP_ARRAY(double, normals);
  
  iGeom_getArrNrmlUV (igeomInstance, 
                         TEMP_ARRAY_IN(gface_handles),
                         static_cast<iBase_StorageOrder>(storage_order),
                         TEMP_ARRAY_IN(parameters),
                         TEMP_ARRAY_INOUT(normals));
  PROCESS_ERROR;

  ASSIGN_ARRAY(normals);
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.getArrNrmlUV)
}

/**
 * Method:  getEntTgntU[]
 */
void
iGeom_SIDL::GeomSidl_impl::getEntTgntU (
  /* in */ void* entity_handle,
  /* in */ double param_coord,
  /* out */ double& tngt_i,
  /* out */ double& tngt_j,
  /* out */ double& tngt_k ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.getEntTgntU)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.getEntTgntU} (getEntTgntU method)
  iGeom_getEntTgntU (igeomInstance, entity_handle, param_coord, &tngt_i, &tngt_j, &tngt_k );
  PROCESS_ERROR;
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.getEntTgntU)
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
iGeom_SIDL::GeomSidl_impl::getArrTgntU (
  /* in */ ::sidl::array<void*> gedge_handles,
  /* in */ int32_t gedge_handles_size,
  /* in */ ::iBase::StorageOrder storage_order,
  /* in */ ::sidl::array<double> parameters,
  /* in */ int32_t parameters_size,
  /* inout */ ::sidl::array<double>& tangents,
  /* inout */ int32_t& tangents_size ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.getArrTgntU)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.getArrTgntU} (getArrTgntU method)
  CREATE_TEMP_ARRAY(double, tangents);
  
  iGeom_getArrTgntU (igeomInstance, 
                         TEMP_ARRAY_IN(gedge_handles),
                         static_cast<iBase_StorageOrder>(storage_order),
                         TEMP_ARRAY_IN(parameters),
                         TEMP_ARRAY_INOUT(tangents));
  PROCESS_ERROR;

  ASSIGN_ARRAY(tangents);
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.getArrTgntU)
}

/**
 * Method:  getEnt1stDrvt[]
 */
void
iGeom_SIDL::GeomSidl_impl::getEnt1stDrvt (
  /* in */ void* entity_handle,
  /* in */ double u,
  /* in */ double v,
  /* inout */ ::sidl::array<double>& drvt_u,
  /* inout */ int32_t& drvt_u_size,
  /* inout */ ::sidl::array<double>& drvt_v,
  /* inout */ int32_t& drvt_v_size ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.getEnt1stDrvt)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.getEnt1stDrvt} (getEnt1stDrvt method)
  CREATE_TEMP_ARRAY(double, drvt_u);
  CREATE_TEMP_ARRAY(double, drvt_v);
  
  iGeom_getEnt1stDrvt (igeomInstance, 
                       entity_handle, u, v,
                       TEMP_ARRAY_INOUT(drvt_u),
                       TEMP_ARRAY_INOUT(drvt_v));
  PROCESS_ERROR;

  ASSIGN_ARRAY(drvt_u);
  ASSIGN_ARRAY(drvt_v);
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.getEnt1stDrvt)
}

/**
 * Method:  getArr1stDrvt[]
 */
void
iGeom_SIDL::GeomSidl_impl::getArr1stDrvt (
  /* in */ ::sidl::array<void*> entity_handles,
  /* in */ int32_t entity_handles_size,
  /* in */ ::iBase::StorageOrder storage_order,
  /* in */ ::sidl::array<double> uv,
  /* in */ int32_t uv_size,
  /* inout */ ::sidl::array<double>& drvt_u,
  /* inout */ int32_t& drvt_u_size,
  /* inout */ ::sidl::array<int32_t>& u_offset,
  /* inout */ int32_t& u_offset_size,
  /* inout */ ::sidl::array<double>& drvt_v,
  /* inout */ int32_t& drvt_v_size,
  /* inout */ ::sidl::array<int32_t>& v_offset,
  /* inout */ int32_t& v_offset_size ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.getArr1stDrvt)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.getArr1stDrvt} (getArr1stDrvt method)
  CREATE_TEMP_ARRAY(double, drvt_u);
  CREATE_TEMP_ARRAY(int32_t, u_offset);
  CREATE_TEMP_ARRAY(double, drvt_v);
  CREATE_TEMP_ARRAY(int32_t, v_offset);
  
  iGeom_getArr1stDrvt (igeomInstance, 
                       TEMP_ARRAY_IN(entity_handles),
                       static_cast<iBase_StorageOrder>(storage_order),
                       TEMP_ARRAY_IN(uv),
                       TEMP_ARRAY_INOUT(drvt_u),
                       TEMP_ARRAY_INOUT(u_offset),
                       TEMP_ARRAY_INOUT(drvt_v),
                       TEMP_ARRAY_INOUT(v_offset));
  PROCESS_ERROR;

  ASSIGN_ARRAY(drvt_u);
  ASSIGN_ARRAY(u_offset);
  ASSIGN_ARRAY(drvt_v);
  ASSIGN_ARRAY(v_offset);
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.getArr1stDrvt)
}

/**
 * Method:  getEnt2ndDrvt[]
 */
void
iGeom_SIDL::GeomSidl_impl::getEnt2ndDrvt (
  /* in */ void* entity_handle,
  /* in */ double u,
  /* in */ double v,
  /* inout */ ::sidl::array<double>& drvt_uu,
  /* inout */ int32_t& drvt_uu_size,
  /* inout */ ::sidl::array<double>& drvt_vv,
  /* inout */ int32_t& drvt_vv_size,
  /* inout */ ::sidl::array<double>& drvt_uv,
  /* inout */ int32_t& drvt_uv_size ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.getEnt2ndDrvt)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.getEnt2ndDrvt} (getEnt2ndDrvt method)
  CREATE_TEMP_ARRAY(double, drvt_uu);
  CREATE_TEMP_ARRAY(double, drvt_uv);
  CREATE_TEMP_ARRAY(double, drvt_vv);
 
  iGeom_getEnt2ndDrvt (igeomInstance, 
                       entity_handle, u, v,
                       TEMP_ARRAY_INOUT(drvt_uu),
                       TEMP_ARRAY_INOUT(drvt_uv),
                       TEMP_ARRAY_INOUT(drvt_vv));
  PROCESS_ERROR;

  ASSIGN_ARRAY(drvt_uu);
  ASSIGN_ARRAY(drvt_uv);
  ASSIGN_ARRAY(drvt_vv);
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.getEnt2ndDrvt)
}

/**
 * Method:  getArr2ndDrvt[]
 */
void
iGeom_SIDL::GeomSidl_impl::getArr2ndDrvt (
  /* in */ ::sidl::array<void*> entity_handles,
  /* in */ int32_t entity_handles_size,
  /* in */ ::iBase::StorageOrder storage_order,
  /* in */ ::sidl::array<double> uv,
  /* in */ int32_t uv_size,
  /* inout */ ::sidl::array<double>& drvt_uu,
  /* inout */ int32_t& drvt_uu_size,
  /* inout */ ::sidl::array<int32_t>& uu_offset,
  /* inout */ int32_t& uu_offset_size,
  /* inout */ ::sidl::array<double>& drvt_vv,
  /* inout */ int32_t& drvt_vv_size,
  /* inout */ ::sidl::array<int32_t>& vv_offset,
  /* inout */ int32_t& vv_offset_size,
  /* inout */ ::sidl::array<double>& drvt_uv,
  /* inout */ int32_t& drvt_uv_size,
  /* inout */ ::sidl::array<int32_t>& uv_offset,
  /* inout */ int32_t& uv_offset_size ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.getArr2ndDrvt)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.getArr2ndDrvt} (getArr2ndDrvt method)
  CREATE_TEMP_ARRAY(double, drvt_uu);
  CREATE_TEMP_ARRAY(int32_t, uu_offset);
  CREATE_TEMP_ARRAY(double, drvt_uv);
  CREATE_TEMP_ARRAY(int32_t, uv_offset);
  CREATE_TEMP_ARRAY(double, drvt_vv);
  CREATE_TEMP_ARRAY(int32_t, vv_offset);
  
  iGeom_getArr2ndDrvt (igeomInstance, 
                       TEMP_ARRAY_IN(entity_handles),
                       static_cast<iBase_StorageOrder>(storage_order),
                       TEMP_ARRAY_IN(uv),
                       TEMP_ARRAY_INOUT(drvt_uu),
                       TEMP_ARRAY_INOUT(uu_offset),
                       TEMP_ARRAY_INOUT(drvt_uv),
                       TEMP_ARRAY_INOUT(uv_offset),
                       TEMP_ARRAY_INOUT(drvt_vv),
                       TEMP_ARRAY_INOUT(vv_offset));
  PROCESS_ERROR;

  ASSIGN_ARRAY(drvt_uu);
  ASSIGN_ARRAY(uu_offset);
  ASSIGN_ARRAY(drvt_uv);
  ASSIGN_ARRAY(uv_offset);
  ASSIGN_ARRAY(drvt_vv);
  ASSIGN_ARRAY(vv_offset);
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.getArr2ndDrvt)
}

/**
 * Method:  getFcCvtrUV[]
 */
void
iGeom_SIDL::GeomSidl_impl::getFcCvtrUV (
  /* in */ void* entity_handle,
  /* in */ double u,
  /* in */ double v,
  /* out */ double& cvtr1_i,
  /* out */ double& cvtr1_j,
  /* out */ double& cvtr1_k,
  /* out */ double& cvtr2_i,
  /* out */ double& cvtr2_j,
  /* out */ double& cvtr2_k ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.getFcCvtrUV)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.getFcCvtrUV} (getFcCvtrUV method)
  iGeom_getFcCvtrUV (igeomInstance, entity_handle, u, v, 
                     &cvtr1_i, &cvtr1_j, &cvtr1_k,
                     &cvtr2_i, &cvtr2_j, &cvtr2_k );
  PROCESS_ERROR;
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.getFcCvtrUV)
}

/**
 * Method:  getFcArrCvtrUV[]
 */
void
iGeom_SIDL::GeomSidl_impl::getFcArrCvtrUV (
  /* in */ ::sidl::array<void*> face_handles,
  /* in */ int32_t face_handles_size,
  /* in */ ::iBase::StorageOrder storage_order,
  /* in */ ::sidl::array<double> uv,
  /* in */ int32_t uv_size,
  /* inout */ ::sidl::array<double>& cvtr_1,
  /* inout */ int32_t& cvtr_1_size,
  /* inout */ ::sidl::array<double>& cvtr_2,
  /* inout */ int32_t& cvtr_2_size ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.getFcArrCvtrUV)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.getFcArrCvtrUV} (getFcArrCvtrUV method)
  CREATE_TEMP_ARRAY(double, cvtr_1);
  CREATE_TEMP_ARRAY(double, cvtr_2);
  
  iGeom_getFcArrCvtrUV (igeomInstance, 
                         TEMP_ARRAY_IN(face_handles),
                         static_cast<iBase_StorageOrder>(storage_order),
                         TEMP_ARRAY_IN(uv),
                         TEMP_ARRAY_INOUT(cvtr_1),
                         TEMP_ARRAY_INOUT(cvtr_2));
  PROCESS_ERROR;

  ASSIGN_ARRAY(cvtr_1);
  ASSIGN_ARRAY(cvtr_2);
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.getFcArrCvtrUV)
}

/**
 * Method:  isEntPeriodic[]
 */
void
iGeom_SIDL::GeomSidl_impl::isEntPeriodic (
  /* in */ void* entity_handle,
  /* out */ int32_t& in_u,
  /* out */ int32_t& in_v ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.isEntPeriodic)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.isEntPeriodic} (isEntPeriodic method)
  iGeom_isEntPeriodic (igeomInstance, 
                       entity_handle,
                       &in_u, &in_v);
  PROCESS_ERROR;
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.isEntPeriodic)
}

/**
 * Method:  isArrPeriodic[]
 */
void
iGeom_SIDL::GeomSidl_impl::isArrPeriodic (
  /* in */ ::sidl::array<void*> entity_handles,
  /* in */ int32_t entity_handles_size,
  /* inout */ ::sidl::array<int32_t>& in_uv,
  /* inout */ int32_t& in_uv_size ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.isArrPeriodic)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.isArrPeriodic} (isArrPeriodic method)
  CREATE_TEMP_ARRAY(int32_t, in_uv);
  iGeom_isArrPeriodic (igeomInstance, 
                       TEMP_ARRAY_IN(entity_handles),
                       TEMP_ARRAY_INOUT(in_uv));
  PROCESS_ERROR;
  ASSIGN_ARRAY(in_uv);
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.isArrPeriodic)
}

/**
 * Method:  isFcDegenerate[]
 */
int32_t
iGeom_SIDL::GeomSidl_impl::isFcDegenerate (
  /* in */ void* face_handle ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.isFcDegenerate)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.isFcDegenerate} (isFcDegenerate method)
  int32_t retval = iGeom_isFcDegenerate (igeomInstance, face_handle);
  PROCESS_ERROR;
  return retval;
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.isFcDegenerate)
}

/**
 * Method:  isFcArrDegenerate[]
 */
void
iGeom_SIDL::GeomSidl_impl::isFcArrDegenerate (
  /* in */ ::sidl::array<void*> face_handles,
  /* in */ int32_t face_handles_size,
  /* inout */ ::sidl::array<int32_t>& degenerate,
  /* inout */ int32_t& degenerate_size ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.isFcArrDegenerate)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.isFcArrDegenerate} (isFcArrDegenerate method)
  CREATE_TEMP_ARRAY(int32_t, degenerate);
  iGeom_isFcArrDegenerate (igeomInstance, 
                       TEMP_ARRAY_IN(face_handles),
                       TEMP_ARRAY_INOUT(degenerate));
  PROCESS_ERROR;
  ASSIGN_ARRAY(degenerate);
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.isFcArrDegenerate)
}

/**
 * Return the relative and absolute tolerances at the modeler level.  If
 * model does not have a modeler-wide tolerance, zero is returned for both
 * values.
 * @param relative_tolerance Relative tolerance for model as a whole
 * @param absolute_tolerance Absolute tolerance for model as a whole
 */
void
iGeom_SIDL::GeomSidl_impl::getTolerance (
  /* out */ int32_t& type,
  /* out */ double& tolerance ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.getTolerance)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.getTolerance} (getTolerance method)
  iGeom_getTolerance (igeomInstance, 
                       &type,
                       &tolerance);
  PROCESS_ERROR;
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.getTolerance)
}

/**
 * Method:  getEntTolerance[]
 */
double
iGeom_SIDL::GeomSidl_impl::getEntTolerance (
  /* in */ void* entity_handle ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.getEntTolerance)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.getEntTolerance} (getEntTolerance method)
  double retval = iGeom_getEntTolerance (igeomInstance, entity_handle);
  PROCESS_ERROR;
  return retval;
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.getEntTolerance)
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
iGeom_SIDL::GeomSidl_impl::getArrTolerance (
  /* in */ ::sidl::array<void*> entity_handles,
  /* in */ int32_t entity_handles_size,
  /* inout */ ::sidl::array<double>& tolerances,
  /* inout */ int32_t& tolerances_size ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.getArrTolerance)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.getArrTolerance} (getArrTolerance method)
  CREATE_TEMP_ARRAY(double, tolerances);
  
  iGeom_getArrTolerance (igeomInstance, 
                         TEMP_ARRAY_IN(entity_handles),
                         TEMP_ARRAY_INOUT(tolerances));
  PROCESS_ERROR;

  ASSIGN_ARRAY(tolerances);
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.getArrTolerance)
}

/**
 * Initialize an iterator over gentities of a specified dimension.
 * @param gentity_dimension Dimension of gentities to be iterated over
 * @param gentity_iterator Iterator initialized by this function
 */
void
iGeom_SIDL::GeomSidl_impl::initEntIter (
  /* in */ int32_t gentity_dimension,
  /* out */ void*& gentity_iterator ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.initEntIter)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.initEntIter} (initEntIter method)
  iGeom_initEntIter (igeomInstance, gentity_dimension, &gentity_iterator);
  PROCESS_ERROR;
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.initEntIter)
}

/**
 * Method:  initEntArrIter[]
 */
bool
iGeom_SIDL::GeomSidl_impl::initEntArrIter (
  /* in */ void* entity_set_handle,
  /* in */ ::iBase::EntityType requested_entity_type,
  /* in */ int32_t requested_array_size,
  /* out */ void*& entArr_iterator ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.initEntArrIter)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.initEntArrIter} (initEntArrIter method)
  bool result = iGeom_initEntArrIter (igeomInstance, 
                        entity_set_handle,
                        static_cast<iGeom_EntityType>(requested_entity_type),
                        requested_array_size,
                        &entArr_iterator);
  PROCESS_ERROR;
  return static_cast<bool>(result);
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.initEntArrIter)
}

/**
 * Get the next entity for this iterator.
 * @param gentity_iterator Iterator being iterated over
 * @param gentity_handle Next gentity
 * @return If true, there are more gentities, if false, this is the last one
 */
bool
iGeom_SIDL::GeomSidl_impl::getNextEntIter (
  /* inout */ void*& gentity_iterator,
  /* out */ void*& gentity_handle ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.getNextEntIter)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.getNextEntIter} (getNextEntIter method)
  bool retval = iGeom_getNextEntIter (igeomInstance, 
                                           &gentity_iterator,
                                           &gentity_handle);
  PROCESS_ERROR;
  return retval;
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.getNextEntIter)
}

/**
 * Method:  getNextEntArrIter[]
 */
bool
iGeom_SIDL::GeomSidl_impl::getNextEntArrIter (
  /* in */ void* entArr_iterator,
  /* inout */ ::sidl::array<void*>& entity_handles,
  /* inout */ int32_t& entity_handles_size ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.getNextEntArrIter)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.getNextEntArrIter} (getNextEntArrIter method)
  CREATE_TEMP_ARRAY(void*, entity_handles);
  bool result = iGeom_getNextEntArrIter (igeomInstance, 
                           entArr_iterator,
                           TEMP_ARRAY_INOUT(entity_handles));
  PROCESS_ERROR;
  ASSIGN_ARRAY(entity_handles);
  return result;
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.getNextEntArrIter)
}

/**
 * Reset an iterator back to the first gentity
 * @param gentity_iterator Iterator reset by this function
 */
void
iGeom_SIDL::GeomSidl_impl::resetEntIter (
  /* in */ void* gentity_iterator ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.resetEntIter)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.resetEntIter} (resetEntIter method)
  iGeom_resetEntIter (igeomInstance, gentity_iterator);
  PROCESS_ERROR;
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.resetEntIter)
}

/**
 * Method:  resetEntArrIter[]
 */
void
iGeom_SIDL::GeomSidl_impl::resetEntArrIter (
  /* in */ void* entArr_iterator ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.resetEntArrIter)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.resetEntArrIter} (resetEntArrIter method)
  iGeom_resetEntArrIter (igeomInstance, entArr_iterator);
  PROCESS_ERROR;
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.resetEntArrIter)
}

/**
 * Delete an iterator
 * @param gentity_iterator Iterator deleted by this function
 */
void
iGeom_SIDL::GeomSidl_impl::endEntIter (
  /* in */ void* Gentity_dim_iterator ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.endEntIter)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.endEntIter} (endEntIter method)
  iGeom_endEntIter (igeomInstance, 
                               Gentity_dim_iterator);
  PROCESS_ERROR;
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.endEntIter)
}

/**
 * Method:  endEntArrIter[]
 */
void
iGeom_SIDL::GeomSidl_impl::endEntArrIter (
  /* in */ void* entArr_iterator ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.endEntArrIter)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.endEntArrIter} (endEntArrIter method)
  iGeom_endEntArrIter (igeomInstance, entArr_iterator);
  PROCESS_ERROR;
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.endEntArrIter)
}

/**
 * Method:  Copy[]
 */
void
iGeom_SIDL::GeomSidl_impl::Copy (
  /* in */ void* geom_entity,
  /* out */ void*& geom_entity2 ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.Copy)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.Copy} (Copy method)
  iGeom_Copy(igeomInstance, geom_entity, &geom_entity2);
  PROCESS_ERROR;
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.Copy)
}

/**
 * Method:  SweepAboutAxis[]
 */
void
iGeom_SIDL::GeomSidl_impl::SweepAboutAxis (
  /* in */ void* geom_entity,
  /* in */ double angle,
  /* in */ double axis_normal_x,
  /* in */ double axis_normal_y,
  /* in */ double axis_normal_z,
  /* out */ void*& geom_entity2 ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.SweepAboutAxis)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.SweepAboutAxis} (SweepAboutAxis method)
  iGeom_SweepAboutAxis(igeomInstance, geom_entity, angle,
                       axis_normal_x, axis_normal_y, axis_normal_z, 
                       &geom_entity2);
  PROCESS_ERROR;
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.SweepAboutAxis)
}

/**
 * Method:  Delete[]
 */
void
iGeom_SIDL::GeomSidl_impl::Delete (
  /* in */ void* geom_entity ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.Delete)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.Delete} (Delete method)
  iGeom_Delete(igeomInstance, geom_entity);
  PROCESS_ERROR;
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.Delete)
}

/**
 * Method:  Sphere[]
 */
void
iGeom_SIDL::GeomSidl_impl::Sphere (
  /* in */ double radius,
  /* out */ void*& geom_entity ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.Sphere)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.Sphere} (Sphere method)
  iGeom_Sphere(igeomInstance, radius, &geom_entity);
  PROCESS_ERROR;
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.Sphere)
}

/**
 * Method:  Brick[]
 */
void
iGeom_SIDL::GeomSidl_impl::Brick (
  /* in */ double x,
  /* in */ double y,
  /* in */ double z,
  /* out */ void*& geom_entity ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.Brick)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.Brick} (Brick method)
  iGeom_Brick(igeomInstance, x, y, z, &geom_entity);
  PROCESS_ERROR;
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.Brick)
}

/**
 * Method:  Cylinder[]
 */
void
iGeom_SIDL::GeomSidl_impl::Cylinder (
  /* in */ double height,
  /* in */ double major_rad,
  /* in */ double minor_rad,
  /* out */ void*& geom_entity ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.Cylinder)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.Cylinder} (Cylinder method)
  iGeom_Cylinder(igeomInstance, height, major_rad, minor_rad, &geom_entity);
  PROCESS_ERROR;
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.Cylinder)
}

/**
 * Method:  Torus[]
 */
void
iGeom_SIDL::GeomSidl_impl::Torus (
  /* in */ double major_rad,
  /* in */ double minor_rad,
  /* out */ void*& geom_entity ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.Torus)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.Torus} (Torus method)
  iGeom_Torus(igeomInstance, major_rad, minor_rad, &geom_entity);
  PROCESS_ERROR;
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.Torus)
}

/**
 * Method:  Move[]
 */
void
iGeom_SIDL::GeomSidl_impl::Move (
  /* inout */ void*& geom_entity,
  /* in */ double x,
  /* in */ double y,
  /* in */ double z ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.Move)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.Move} (Move method)
  iGeom_Move(igeomInstance, &geom_entity, x, y, z);
  PROCESS_ERROR;
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.Move)
}

/**
 * Method:  Rotate[]
 */
void
iGeom_SIDL::GeomSidl_impl::Rotate (
  /* inout */ void*& geom_entity,
  /* in */ double angle,
  /* in */ double axis_normal_x,
  /* in */ double axis_normal_y,
  /* in */ double axis_normal_z ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.Rotate)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.Rotate} (Rotate method)
  iGeom_Rotate(igeomInstance, &geom_entity, angle, axis_normal_x, 
               axis_normal_y, axis_normal_z);
  PROCESS_ERROR;
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.Rotate)
}

/**
 * Method:  Reflect[]
 */
void
iGeom_SIDL::GeomSidl_impl::Reflect (
  /* inout */ void*& geom_entity,
  /* in */ double plane_normal_x,
  /* in */ double plane_normal_y,
  /* in */ double plane_normal_z ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.Reflect)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.Reflect} (Reflect method)
  iGeom_Reflect(igeomInstance, &geom_entity, plane_normal_x, 
                plane_normal_y, plane_normal_z);
  PROCESS_ERROR;
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.Reflect)
}

/**
 * Method:  Scale[]
 */
void
iGeom_SIDL::GeomSidl_impl::Scale (
  /* inout */ void*& geom_entity,
  /* in */ double scale_x,
  /* in */ double scale_y,
  /* in */ double scale_z ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.Scale)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.Scale} (Scale method)
  iGeom_Scale(igeomInstance, &geom_entity, scale_x, scale_y, scale_z);
  PROCESS_ERROR;
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.Scale)
}

/**
 * Method:  Unite[]
 */
void
iGeom_SIDL::GeomSidl_impl::Unite (
  /* in */ ::sidl::array<void*> geom_entities,
  /* in */ int32_t geom_entities_size,
  /* out */ void*& geom_entity ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.Unite)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.Unite} (Unite method)
  iGeom_Unite(igeomInstance, 
              TEMP_ARRAY_IN(geom_entities), &geom_entity);
  PROCESS_ERROR;
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.Unite)
}

/**
 * Method:  Subtract[]
 */
void
iGeom_SIDL::GeomSidl_impl::Subtract (
  /* in */ void* blank,
  /* in */ void* tool,
  /* out */ void*& geom_entity ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.Subtract)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.Subtract} (Subtract method)
  iGeom_Subtract(igeomInstance, blank, tool, &geom_entity);
  PROCESS_ERROR;
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.Subtract)
}

/**
 * Method:  Section[]
 */
void
iGeom_SIDL::GeomSidl_impl::Section (
  /* inout */ void*& geom_entity,
  /* in */ double plane_normal_x,
  /* in */ double plane_normal_y,
  /* in */ double plane_normal_z,
  /* in */ double offset,
  /* in */ bool reverse,
  /* out */ void*& geom_entity2 ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.Section)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.Section} (Section method)
  iGeom_Section(igeomInstance, &geom_entity, plane_normal_x, 
                plane_normal_y, plane_normal_z, 
                offset, reverse, &geom_entity2);
  PROCESS_ERROR;
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.Section)
}

/**
 * Method:  Imprint[]
 */
void
iGeom_SIDL::GeomSidl_impl::Imprint (
  /* in */ ::sidl::array<void*> geom_entities,
  /* in */ int32_t geom_entities_size ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.Imprint)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.Imprint} (Imprint method)
  iGeom_Imprint (igeomInstance,
                 TEMP_ARRAY_IN(geom_entities));
  PROCESS_ERROR;
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.Imprint)
}

/**
 * Method:  Merge[]
 */
void
iGeom_SIDL::GeomSidl_impl::Merge (
  /* in */ ::sidl::array<void*> geom_entities,
  /* in */ int32_t geom_entities_size,
  /* in */ double tolerance ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl.Merge)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl.Merge} (Merge method)
  iGeom_Merge (igeomInstance,
               TEMP_ARRAY_IN(geom_entities),
               tolerance);
  PROCESS_ERROR;
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl.Merge)
}


// DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl._misc)
// Insert-Code-Here {iGeom_SIDL.GeomSidl._misc} (miscellaneous code)
// call this function when you want to throw an error
void
iGeom_SIDL::GeomSidl_impl::processError() throw(::iBase::Error)
{
  static void *behavior_tag = NULL;
  static ::iBase::ErrorActions action = ::iBase::ErrorActions_THROW_ERROR;

    // save this info before calling tag get function to get behavior
  std::string this_desc(iGeom_LAST_ERROR.description);
  iBase::ErrorType this_err = (iBase::ErrorType)iGeom_LAST_ERROR.error_type;

  if (0 == behavior_tag)
    behavior_tag = iGeom_getTagHandle(igeomInstance, "Error_Behavior");

  if (0 != behavior_tag) {
    ::iBase::ErrorActions temp = (::iBase::ErrorActions) iGeom_getIntData (igeomInstance, NULL,
                                                                         behavior_tag);
    if (iBase_SUCCESS == iGeom_LAST_ERROR.error_type) action = temp;
  }

  ::iBase::Error this_error = ::iBase::Error::_create();
  switch (action) {
    case ::iBase::ErrorActions_THROW_ERROR:
      this_error.set(this_err, this_desc);
      throw this_error;
      break;
    case ::iBase::ErrorActions_WARN_ONLY:
      std::cerr << this_desc << std::endl;
      break;
    case ::iBase::ErrorActions_SILENT:
      break;
  }
}
// DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl._misc)

#error File has unused splicer blocks.
/**
 * ================= BEGIN UNREFERENCED METHOD(S) ================
 * The following code segment(s) belong to unreferenced method(s).
 * This can result from a method rename/removal in the sidl file.
 * Move or remove the code in order to compile cleanly.
 */
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl_impl.getFaceType)
  const char *this_type = NULL;
  int32_t retval = iGeom_isEntParametric (igeomInstance, 
                                              gentity_handle);
  PROCESS_ERROR;
  return retval;
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl_impl.getFaceType)
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl_impl.measure)
  CREATE_TEMP_ARRAY(double, measures);
  int32_t retval = iGeom_getParametric();
  PROCESS_ERROR;
  return retval;
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl_impl.measure)
// ================== END UNREFERENCED METHOD(S) =================
