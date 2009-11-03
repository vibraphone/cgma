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
 
/**\file CATag.cpp
 *
 *\author Tim Tautges
 *\author Jason Kraftcheck
 *
 * Original file from SNL TSTT repository was named CATag.
 *
 * Renamed CATag and added to ANL ITAPS repository by J.Kraftcheck,
 * 2007-6-15
 */

#include <algorithm>
 
#include "CATag.hpp"
#include "RefEntity.hpp"
#include "RefEntityName.hpp"
#include "RefEntityFactory.hpp"
#include "CubitAttribManager.hpp"
#include "RefGroup.hpp"
#include "TDUniqueId.hpp"
#include "CGMApp.hpp"
#include "iGeomError.h"
#include "TopologyEntity.hpp"

#define CHECK_SIZE(array, type, this_size, retval)  \
  if (0 == array ## _allocated || array ## _allocated < this_size) {\
    if (NULL != array) free(array); \
    array = (type*)malloc(this_size*sizeof(type));\
    array ## _allocated=this_size;\
    if (NULL == array) { \
      iGeom_setLastError(iBase_MEMORY_ALLOCATION_FAILED);\
      return retval; \
    }\
  }; \
  array ## _size = this_size;

static inline iBase_ErrorType tag_check_size( char*& array,
                                              int& allocated,
                                              int size )
{
  if (!array || !allocated) {
    allocated = size;
    array = (char*)malloc(allocated);
    if (!array) 
      return iBase_MEMORY_ALLOCATION_FAILED;
  }
  else if (allocated < size) {
    return iBase_BAD_ARRAY_SIZE;
  }
  
  return iBase_SUCCESS;
}

#define TAG_CHECK_SIZE(array, allocated, size)  \
  if (iBase_ErrorType tag_check_size_tmp = tag_check_size( array, allocated, size )) \
    return tag_check_size_tmp;

#define RETURN(a) {iGeom_setLastError(a); return a;}

const char *CGMTagManager::CATag_NAME = "ITAPS_Tag";
const char *CGMTagManager::CATag_NAME_INTERNAL = "ITAPS_TAG";


class CSATagData 
{
public:
  CSATagData() 
      : dblData(0), dblDataSize(0), intData(0), intDataSize(0),
        stringData(0), stringDataSize(0), indivStringSize(0) 
    {}

  ~CSATagData() 
    {
    }

  void reset();

  double *dblData;
  int dblDataSize;
  
  int *intData;
  int intDataSize;
  
  char *stringData;
  int stringDataSize;
  int indivStringSize;
};


CubitAttrib *CGMTagManager::CATag_creator(RefEntity* entity, CubitSimpleAttrib *p_csa)
{
  CATag *this_ca = new CATag(&instance(), entity, p_csa);
  return this_ca;
}
 

static CGMTagManager::TagInfo preset_tag_list[] = {
   // tag size      tag name           tag data type  default active
 { 0,              "",                 iBase_BYTES,   NULL,  false },
 { 32,             "NAME",             iBase_BYTES,   NULL,   true },
 { sizeof(int),    "GLOBAL_ID",        iBase_INTEGER, NULL,   true },
 { sizeof(int),    "UNIQUE_ID",        iBase_INTEGER, NULL,   true },
 { sizeof(int),    "MESH_INTERVAL",    iBase_INTEGER, NULL,   true },
 { sizeof(double), "MESH_SIZE",        iBase_DOUBLE,  NULL,   true },
 { 4,              "SIZE_FIRMNESS",    iBase_BYTES,   NULL,   true } };
 

CGMTagManager::TagInfo* const CGMTagManager::presetTagInfo = preset_tag_list;
const int CGMTagManager::numPresetTag = sizeof(preset_tag_list)/sizeof(preset_tag_list[0]);

CGMTagManager::CGMTagManager() 
    : interfaceGroup(NULL)
{
  pcTag = 0;
  tagInfo.push_back(preset_tag_list[0]);
  
    // get the tag number for CATag
  DLIList<int> tag_types;
  int max_type = 0;
  CubitAttribManager *cam = CGMApp::instance()->attrib_manager();
  cam->get_registered_types(tag_types);
  for (int i = 0; i < tag_types.size(); i++) {
    int this_type = tag_types.get_and_step();
    max_type = (max_type < this_type ? this_type : max_type);
  }
  
  max_type++;
  CubitStatus status = cam->register_attrib_type(max_type, CATag_NAME, CATag_NAME_INTERNAL,
                                                 &CGMTagManager::CATag_creator, false, true, 
                                                 true, true, true, false);
  if (CUBIT_FAILURE == status) {
    iGeom_setLastError( iBase_FAILURE, "Couldn't create cgm attribute for tags." );
  }
  else
    CATag_att_type = max_type;

    // create preset tags, for CGM attributes we want to be visible as tags
    // name - make same as in MBTagConventions
  for (int i = 1; i < numPresetTag; ++i)
    tagNameMap[presetTagInfo[i].tagName] = -i; // neg handles beginning with -1
}

CGMTagManager::~CGMTagManager() 
{

}

long CGMTagManager::pc_tag(const bool create_if_missing) 
{
  if (0 == pcTag && create_if_missing) {
    pcTag = getTagHandle("__cgm_parent_child_tag");
    if (0 == pcTag) {
        // ok, one doesn't exist, create one
      void *def_val[] = {NULL, NULL};
      
      createTag("__cgm_parent_child_tag", 2*sizeof(std::vector<RefGroup*>*),
                iBase_BYTES, (char*)def_val, &pcTag);
    }
  }
  
  return pcTag;
}

void CSATagData::reset()  
{
  if (0 != dblData) free(dblData);
  if (0 != intData) free(intData);
  if (0 != stringData) free(stringData);
}
  
iBase_ErrorType CGMTagManager::create_csa_tag(const char *tag_name, long *tag_handle) 
{
    // make a special tag type to hold csa data from CGM; this tag holds, in this order:
    // a) dbl_data (double*), dbl_data_size (int)
    // b) int_data (int*), int_data_size (int)
    // c) string_data (char*), string_data_size (int), indiv_string_size (int)
    //
    // For string data, multiple strings are packed into a single string, with each string
    // occupying indiv_string_size bytes and string_size indicating number of such strings
    // in packed string
    // 
    // So, size of this tag is 3*sizeof(void*) + 4*sizeof(int), assuming pointers occupy
    // same space regardless of what they point to

  CSATagData dum_data;
  return createTag(tag_name, sizeof(CSATagData), iBase_BYTES, (char*)&dum_data, tag_handle);
}

iBase_ErrorType CGMTagManager::set_csa_tag(RefEntity *this_ent,
                                           long tag_handle,
                                           CubitSimpleAttrib *csa_ptr) 
{
  CSATagData this_data;
  CSATagData *this_data_ptr = &this_data;
  int this_data_size = sizeof(CSATagData), this_data_alloc = this_data_size;

    // if data was set on this tag, reset in an attempt to not leak memory
  iBase_ErrorType result = getArrData(&this_ent, 1, tag_handle, 
                                      (char**)&this_data_ptr, &this_data_alloc, &this_data_size);
  if (iBase_SUCCESS != result && iBase_TAG_NOT_FOUND != result) return result;
  
  if (iBase_TAG_NOT_FOUND != result) this_data.reset();
  
    // set the data
  this_data.dblDataSize = csa_ptr->double_data_list()->size();
  this_data.dblData = (double*) malloc(this_data.dblDataSize*sizeof(double));
  csa_ptr->double_data_list()->reset();
  for (int i = 0; i < this_data.dblDataSize; i++)
    this_data.dblData[i] = *(csa_ptr->double_data_list()->get_and_step());
  
  this_data.intDataSize = csa_ptr->int_data_list()->size();
  this_data.intData = (int*) malloc(this_data.intDataSize*sizeof(int));
  csa_ptr->int_data_list()->reset();
  for (int i = 0; i < this_data.intDataSize; i++)
    this_data.intData[i] = *(csa_ptr->int_data_list()->get_and_step());
  
    // find longest string, then allocate
  DLIList<CubitString*> *sd = csa_ptr->string_data_list();
  this_data.indivStringSize = 0;
  for (int i = 0; i < sd->size(); i++) {
    if (((int) sd->get()->length()) > this_data.indivStringSize) 
      this_data.indivStringSize = sd->get()->length();
    sd->step();
  }
    // round to next highest multiple of sizeof(int)
  if (this_data.indivStringSize%sizeof(int) != 0) 
    this_data.indivStringSize = ((this_data.indivStringSize/sizeof(int))+1)*sizeof(int);
    // now allocate
  this_data.stringDataSize = sd->size()*this_data.indivStringSize;
  this_data.stringData = (char*) malloc(this_data.stringDataSize);
  sd->reset();
  for (int i = 0; i < sd->size(); i++)
    strncpy(this_data.stringData+i*this_data.indivStringSize, sd->get_and_step()->c_str(),
            this_data.indivStringSize);

  result = setArrData(&this_ent, 1, tag_handle, (const char*)&this_data_ptr, sizeof(CSATagData));
  
  return result;
}
  
iBase_ErrorType CGMTagManager::createTag (/*in*/ const char *tag_name,
                                          /*in*/ const int tag_length,
                                          /*in*/ const int tag_type,
                                          /*in*/ char* default_value,
                                          /*out*/ long *tag_handle)
{
  std::string tmp_name(tag_name);
  TagInfo tmp_info = {tag_length, tmp_name, tag_type, NULL, true};

  std::map<std::string,long>::iterator mit = tagNameMap.find(tmp_name);
  if (mit != tagNameMap.end()) {
    // we found a tag with this name; is it still active?
    bool active = (mit->second > 0 ? tagInfo[mit->second] :
                   presetTagInfo[-mit->second]).isActive;
    *tag_handle = mit->second;
    if (active) {
      iGeom_setLastError( iBase_TAG_ALREADY_EXISTS );
      return iBase_TAG_ALREADY_EXISTS;
    }

    tagInfo[*tag_handle] = tmp_info;
  }
  else {
    // create a new tag entirely
    tagInfo.push_back(tmp_info);
    *tag_handle = tagInfo.size() - 1;

    // put the name and handle into the map too
    tagNameMap[std::string(tag_name)] = *tag_handle;
  }

  if (default_value != NULL) {
    tagInfo[*tag_handle].defaultValue = (char *) malloc(tag_length);
    memcpy(tagInfo[*tag_handle].defaultValue, default_value, tag_length);
  }

  RETURN(iBase_SUCCESS);
}


iBase_ErrorType CGMTagManager::destroyTag (/*in*/ const long tag_handle,
                                           /*in*/ const bool forced)
{
  if (!forced) {
      // see whether this tag is still assigned anywhere
      // not implemented yet
    RETURN(iBase_NOT_SUPPORTED);
  }
  
    // if we got here, assume we can delete it
  TagInfo *tinfo = (tag_handle > 0 ? &tagInfo[tag_handle] : &presetTagInfo[-tag_handle]);
  tinfo->isActive = false;

  RETURN(iBase_SUCCESS);
}

const char *CGMTagManager::getTagName (/*in*/ const long tag_handle)
{
  iGeom_clearLastError();
  return (tag_handle > 0 ? tagInfo[tag_handle].tagName.c_str() : 
          presetTagInfo[-tag_handle].tagName.c_str());
}

int CGMTagManager::getTagSize (/*in*/ const long tag_handle)
{
  iGeom_clearLastError();
  return (tag_handle > 0 ? 
          tagInfo[tag_handle].tagLength : 
          presetTagInfo[-tag_handle].tagLength);
}

long CGMTagManager::getTagHandle (/*in*/ const char *tag_name)
{
  std::map<std::string,long>::iterator it =
    tagNameMap.find(std::string(tag_name));
  if (it != tagNameMap.end()) {
    bool active = (it->second > 0 ? tagInfo[it->second] :
                   presetTagInfo[-it->second]).isActive;
    if (active) {
      iGeom_clearLastError();
      return it->second;
    }
  }

  iGeom_setLastError( iBase_TAG_NOT_FOUND );
  return 0;
}

int CGMTagManager::getTagType (/*in*/ const long tag_handle) 
{
  iGeom_clearLastError();
  return (tag_handle > 0 ? 
          tagInfo[tag_handle].tagType : 
          presetTagInfo[-tag_handle].tagType);
}

iBase_ErrorType CGMTagManager::getArrData (ARRAY_IN_DECL(RefEntity*, entity_handles),
                                           /*in*/ const long tag_handle,
                                           /*inout*/ ARRAY_INOUT_DECL(char, tag_value))
{
  TagInfo *tinfo = (tag_handle > 0 ? &tagInfo[tag_handle] : &presetTagInfo[-tag_handle]);
  int tag_size = tinfo->tagLength;
    // either way, we have to have that many bytes when we leave this function
  const bool allocated_data_arr = (*tag_value_allocated == 0);
  TAG_CHECK_SIZE(*tag_value, *tag_value_allocated, entity_handles_size*tinfo->tagLength);
  char *val_ptr = *tag_value;
  if (tag_handle < 0) {
    for (int i = 0; i < entity_handles_size; i++) {
      bool result;
      if (NULL == entity_handles[i])
        result = getPresetTagData(interface_group(), tag_handle, val_ptr, tinfo->tagLength);
      else
        result = getPresetTagData(entity_handles[i], tag_handle, val_ptr, tinfo->tagLength);
      if (!result) {
        if (allocated_data_arr) {
          free(*tag_value);
          *tag_value = 0;
          *tag_value_allocated = 0;
        }
        RETURN(iBase_TAG_NOT_FOUND);
      }
      val_ptr += tinfo->tagLength;
    }
    *tag_value_size = entity_handles_size*tinfo->tagLength;
    RETURN(iBase_SUCCESS);
  }

  iBase_ErrorType result = iBase_SUCCESS, tmp_result;

  for (int i = 0; i < entity_handles_size; i++) {
      // ok to cast away const-ness because "false" passed in for create_if_missing
    RefEntity *this_ent = (NULL == entity_handles[i] ? interface_group() :
                           const_cast<RefEntity*>(entity_handles[i]));
    CATag *catag = get_catag(this_ent);
    if (NULL != catag) {
      tmp_result = catag->get_tag_data(tag_handle, val_ptr);
      if (iBase_SUCCESS != tmp_result)
        iGeom_setLastError( tmp_result, "Problem getting tag data." );
    }
    else if (NULL != tinfo->defaultValue) {
      memcpy(val_ptr, tinfo->defaultValue, tinfo->tagLength);
      tmp_result = iBase_SUCCESS;
    }
    else {
      tmp_result = iBase_TAG_NOT_FOUND;
      iGeom_setLastError( iBase_TAG_NOT_FOUND );
    }

    if (iBase_SUCCESS != tmp_result) result = tmp_result;
    
    val_ptr += tag_size;
  }

  if (iBase_SUCCESS != result)
    *tag_value_size = 0;
  else
    *tag_value_size = entity_handles_size*tinfo->tagLength;

  RETURN(result);
}

iBase_ErrorType CGMTagManager::setArrData (/*in*/ ARRAY_IN_DECL(RefEntity*, entity_handles),
                                           /*in*/ const long tag_handle,
                                           /*in*/ const char *tag_values, const int tag_values_size)
{
  TagInfo *tinfo = (tag_handle > 0 ? &tagInfo[tag_handle] : &presetTagInfo[-tag_handle]);
  int tag_size = tinfo->tagLength;
  
  const char *val_ptr = tag_values;

  iBase_ErrorType result = iBase_SUCCESS, tmp_result;
  
  if (tag_handle < 0) {
    for (int i = 0; i < entity_handles_size; i++) {
      if (NULL == entity_handles[i])
        tmp_result = setPresetTagData(interface_group(), tag_handle, val_ptr, tag_size);
      else
        tmp_result = setPresetTagData(entity_handles[i], tag_handle, val_ptr, tag_size);

      val_ptr += tag_size;
      if (iBase_SUCCESS != tmp_result) result = tmp_result;
    }
    RETURN(result);
  }

  for (int i = 0; i < entity_handles_size; i++) {
    RefEntity *this_ent = (NULL == entity_handles[i] ? interface_group() : 
                           entity_handles[i]);
    CATag *catag = get_catag(this_ent, true);
    assert(NULL != catag);
    catag->set_tag_data(tag_handle, val_ptr);
    val_ptr += tag_size;
  }

  RETURN(iBase_SUCCESS);
}

iBase_ErrorType CGMTagManager::rmvArrTag (/*in*/ ARRAY_IN_DECL(RefEntity*, entity_handles),
                                          /*in*/ const long tag_handle)
{
  for (int i = 0; i < entity_handles_size; i++) {
    CATag *catag = get_catag((entity_handles[i] == NULL ? 
                                 interface_group() : entity_handles[i]));
    if (NULL != catag) catag->remove_tag(tag_handle);
  }

  RETURN(iBase_SUCCESS);
}

iBase_ErrorType CGMTagManager::getAllTags (/*in*/ const RefEntity* entity_handle,
                                           /*inout*/ ARRAY_INOUT_DECL(long, tag_handles))
{
  int i = 0, uid = 0, tag_size;
  char *uid_ptr = (char*) &uid;
  bool has_uid = getPresetTagData(entity_handle, -3, uid_ptr, tag_size);
  int num_tags = (has_uid ? 3 : 2);

  RefEntity *this_ent = (NULL == entity_handle ? interface_group() : 
                         const_cast<RefEntity*>(entity_handle));
  
    // const-cast because we're passing in false for create_if_missing
  CATag *catag = get_catag(this_ent);
  if (NULL != catag) {
      // need to check whether entity has a uid
    num_tags += catag->tagData.size();
    CHECK_SIZE(*tag_handles, long, num_tags, iBase_FAILURE);
    for (std::map<int,void*>::iterator tag_it = catag->tagData.begin();
         tag_it != catag->tagData.end(); tag_it++)
      (*tag_handles)[i++] = (*tag_it).first;
  }
  else {
    CHECK_SIZE(*tag_handles, long, num_tags, iBase_FAILURE);
  }
  (*tag_handles)[i++] = -1;
  (*tag_handles)[i++] = -2;
  if (has_uid) (*tag_handles)[i++] = -3;

  RETURN(iBase_SUCCESS);
}

CATag *CGMTagManager::get_catag(RefEntity *ent, 
                                  const bool create_if_missing) 
{
  CubitAttrib *this_attrib = ent->get_cubit_attrib(CATag_att_type, create_if_missing);
  if (NULL != this_attrib)
    return dynamic_cast<CATag*>(this_attrib);
  else
    return NULL;
}

RefGroup *CGMTagManager::interface_group(const bool create_if_missing) 
{
  if (NULL == interfaceGroup) 
    interfaceGroup = 
      dynamic_cast<RefGroup*>(RefEntityName::instance()->get_refentity("interface_group"));
  
  if (NULL == interfaceGroup && create_if_missing)
    interfaceGroup = RefEntityFactory::instance()->construct_RefGroup("interface_group");

  return interfaceGroup;
}
  
iBase_ErrorType CGMTagManager::setPresetTagData(RefEntity *entity, 
                                                const long tag_handle, 
                                                const char *tag_value, 
                                                const int tag_size) 
{
  switch (-tag_handle) {
    case 1:
        // entity name
      if (presetTagInfo[-tag_handle].tagLength != tag_size) {
        std::string tmp_str = "Tag of type '";
        tmp_str += presetTagInfo[-tag_handle].tagName + "' is the wrong size.";
        iGeom_setLastError(iBase_INVALID_ARGUMENT, tmp_str.c_str());
        return iBase_INVALID_ARGUMENT;
      }
      entity->entity_name(CubitString(tag_value));
      RETURN(iBase_SUCCESS);
    case 2:
        // entity id
      iGeom_setLastError( iBase_NOT_SUPPORTED, "Can't set id of entities with this implementation." );
      return iBase_NOT_SUPPORTED;
    case 3:
        // unique id
      iGeom_setLastError( iBase_NOT_SUPPORTED, "Can't set unique id of entities with this implementation." );
      return iBase_NOT_SUPPORTED;
    case 4: // mesh interval
    case 5: // mesh size
    case 6: // mesh interval firmness
    default:
      iGeom_setLastError( iBase_NOT_SUPPORTED, "Can't set this tag on entities with this implementation." );
      return iBase_NOT_SUPPORTED;
  }

  iGeom_setLastError( iBase_TAG_NOT_FOUND );
  return iBase_TAG_NOT_FOUND;
}

CubitSimpleAttrib* CGMTagManager::get_simple_attrib(RefEntity* entity,
                                                    const char* name )
{
  TopologyEntity* te_ptr = dynamic_cast<TopologyEntity*>(entity);
  TopologyBridge* tb_ptr = te_ptr ? te_ptr->bridge_manager()->topology_bridge() : 0;
  if (!tb_ptr) {
    iGeom_setLastError( iBase_INVALID_ENTITY_HANDLE, "Entity not topology" );
    return 0;
  }
  DLIList<CubitSimpleAttrib*> attr_list;
  tb_ptr->get_simple_attribute("MESH_INTERVAL", attr_list);
  if (attr_list.size() == 0) {
    iGeom_setLastError( iBase_TAG_NOT_FOUND, "No MESH_INTERVAL attribute" );
    return 0;
  }
  CubitSimpleAttrib* result = attr_list.pop();
  while (attr_list.size() != 0)
    delete attr_list.pop();
  return result;
}


bool CGMTagManager::getPresetTagData(const RefEntity *entity, 
                                     const long tag_handle, 
                                     char *tag_value, 
                                     int &tag_size) 
{
  const char *this_name;
  int name_len, val;
  double dval;
  int *this_id;
  int *this_uid;
  
  if (-tag_handle >= numPresetTag || tag_handle >= 0) {
    iGeom_setLastError( iBase_INVALID_TAG_HANDLE, "Invalid tag handle" );
    return false;
  }
  
  const TagInfo& info = presetTagInfo[-tag_handle];
  tag_size = info.tagLength;
  CubitSimpleAttrib* csa;
  CubitString* str;
  
  switch (-tag_handle) {
    case 1:
        // entity name
      this_name = entity->entity_name().c_str();
      name_len = strlen(this_name);
        // if name is too long, truncate
      if (name_len > info.tagLength)
        name_len = info.tagLength;
      strncpy( tag_value, this_name, name_len );
        // if name is too short, pad with zero bytes
      memset( tag_value + name_len, 0, info.tagLength );
      return true;
    case 2:
        // entity id
      tag_size = sizeof(int);
      this_id = reinterpret_cast<int*>(tag_value);
      *this_id = entity->id();
      return true;
    case 3:
        // unique id
      tag_size = sizeof(int);
      this_uid = reinterpret_cast<int*>(tag_value);
        // const_cast because we're passing false for create_if_missing
      *this_uid = TDUniqueId::get_unique_id(const_cast<RefEntity*>(entity), false);
      return (*this_uid == 0 ? false : true);
    case 4: // mesh interval
      csa = get_simple_attrib( const_cast<RefEntity*>(entity), "MESH_INTERVAL" );
      if (!csa)
        return false;
      csa->int_data_list()->reset();
      val = *csa->int_data_list()->get_and_step();
        // check if interval is set
      csa->string_data_list()->reset();
        // If a) the size is set and b) the firmness is LIMP, then
        // the interval count has not been set.
      if ( csa->string_data_list()->size() && 
          *csa->string_data_list()->get() == "LIMP" && 
           csa->int_data_list()->size() > 1 && 
         !*csa->int_data_list()->get())
        val = 0;
      delete csa;
      if (val == 0 || val == CUBIT_INT_MIN) {
        if (info.defaultValue)
          val = *(int*)info.defaultValue;
        else {
          iGeom_setLastError( iBase_TAG_NOT_FOUND, "Interval not set" );
          return 0;
        }
      }
      *(int*)tag_value = val;
      return true;
    case 5: // mesh size
      csa = get_simple_attrib( const_cast<RefEntity*>(entity), "MESH_INTERVAL" );
      if (!csa)
        return false;
      csa->int_data_list()->reset();
      csa->int_data_list()->step();
      csa->double_data_list()->reset();
        // if size value is invalid or flag indicates size has not been set...
      if (csa->double_data_list()->size() == 0 ||
         *csa->double_data_list()->get() == CUBIT_DBL_MIN ||
         (csa->int_data_list()->size() > 1 && !*csa->int_data_list()->get())) {
        if (info.defaultValue)
          dval = *(double*)info.defaultValue;
        else {
          delete csa;
          iGeom_setLastError( iBase_TAG_NOT_FOUND, "Mesh size not set" );
          return false;
        }
      }
      else
        dval = *csa->double_data_list()->get();
      delete csa;
      *(double*)tag_value = dval;
      return true;
    case 6: // interval firmness
      csa = get_simple_attrib( const_cast<RefEntity*>(entity), "MESH_INTERVAL" );
      if (!csa)
        return false;
      if (csa->string_data_list()->size() < 2) {
        delete csa;
        iGeom_setLastError( iBase_TAG_NOT_FOUND, "Interval firmness not set" );
        return false;
      }
      csa->string_data_list()->reset();
      str = csa->string_data_list()->step_and_get();
      memcpy( tag_value, str->c_str(), 4 );
      delete csa;
      return true;
  }

  return false;
}
        
std::vector<RefGroup*> *CGMTagManager::pc_list(RefEntity *gentity, int list_no, 
                                               const bool create_if_missing) 
{
  if (NULL == gentity) return NULL;
  
  int dum_tag_size = sizeof(std::vector<RefGroup*>*);
  int dum = 2*dum_tag_size;
  std::vector<RefGroup*> *pc_lists[2];
  char *dum_ptr = (char*) pc_lists;
  iBase_ErrorType result =
    getArrData(&gentity, 1, pc_tag(create_if_missing), &dum_ptr, 
               &dum, &dum);
  assert(iBase_SUCCESS == result);
  
  if (0 > list_no || 1 < list_no) return NULL;
  
  if (NULL == pc_lists[list_no] && create_if_missing) {
    pc_lists[list_no] = new std::vector<RefGroup*>();
    result =
      setArrData(&gentity, 1, pc_tag(), (char*)pc_lists, 2*dum_tag_size);
  }

  return pc_lists[list_no];
}

void CGMTagManager::pc_list(RefEntity *gentity, std::vector<RefGroup*> *&parents,
                            std::vector<RefGroup*> *&children,
                            const bool create_if_missing)
{
  if (NULL == gentity) return;

  int dum_tag_size = sizeof(std::vector<RefGroup*>*);
  int dum = 2*dum_tag_size;
  std::vector<RefGroup*> *pc_lists[2];
  char *dum_ptr = (char*) pc_lists;
  iBase_ErrorType result =
    getArrData(&gentity, 1, pc_tag(), 
               &dum_ptr, &dum, &dum);
  assert(iBase_SUCCESS == result);

  if ((NULL == pc_lists[0] || NULL == pc_lists[1]) && create_if_missing) {
    bool must_set = false;
    for (int i = 0; i < 2; i++) {
      if (NULL == pc_lists[i]) {
        pc_lists[i] = new std::vector<RefGroup*>();
        must_set = true;
      }
    }
    if (must_set) {
      result = setArrData(&gentity, 1, pc_tag(), (char*)pc_lists, 
                          2*dum_tag_size);
      assert(iBase_SUCCESS == result);
    }
  }
  
  parents = pc_lists[0];
  children = pc_lists[1];
}

void CGMTagManager::get_pc_groups(RefGroup *this_grp, const int p_or_c, const int num_hops,
                                  std::vector<RefGroup *> &group_ptrs) 
{
  if (NULL == this_grp) return;
  
  int next_hop = num_hops - 1;
  std::vector<RefGroup*> tmp_groups;
  
    // get my children
  std::vector<RefGroup*> *my_children = pc_list(this_grp, p_or_c, false);
  if (NULL != my_children)
    tmp_groups = *my_children;
  
    // get their children if we're not out of hops
  std::vector<RefGroup*>::iterator vit;
  if (0 < next_hop) {
    for (vit = tmp_groups.begin(); vit != tmp_groups.end(); vit++)
      get_pc_groups(*vit, p_or_c, next_hop, group_ptrs);
  }
  
    // now add mine to the list; make sure it isn't there already
  for (vit = tmp_groups.begin(); vit != tmp_groups.end(); vit++) {
    if (std::find(group_ptrs.begin(), group_ptrs.end(), *vit) == group_ptrs.end())
      group_ptrs.push_back(*vit); 
  }
}

CATag::~CATag() 
{
  for (std::map<int, void*>::iterator 
         mit = tagData.begin(); mit != tagData.end(); mit++)
    if (NULL != (*mit).second) free ((*mit).second);
}


CATag::CATag(CGMTagManager *manager, RefEntity *entity) 
    : CubitAttrib(entity), myManager(manager)
{
}

CATag::CATag(CGMTagManager *manager, RefEntity *owner, CubitSimpleAttrib *csa_ptr) 
    : CubitAttrib(owner), myManager(manager)
{
  if (NULL != csa_ptr) add_csa_data(csa_ptr);
}

CubitStatus CATag::reset()
{
  for (std::map<int, void*>::iterator 
         mit = tagData.begin(); mit != tagData.end(); mit++)
    if (NULL != (*mit).second) free ((*mit).second);

  tagData.clear();

  return CUBIT_SUCCESS;
}

CubitSimpleAttrib* CATag::cubit_simple_attrib() 
{
    //if (tagData.size() == 0) return NULL;
  
  DLIList<int> int_data;
  DLIList<CubitString*> str_data;
  DLIList<double> dbl_data;

  str_data.append(new CubitString(myManager->CATag_NAME_INTERNAL));

    // int data first gets the # tags on this entity
  int_data.append(tagData.size());

    // for each tag:
  for (std::map<int, void*>::iterator 
         mit = tagData.begin(); mit != tagData.end(); mit++) {
    long tag_handle = (*mit).first;
    CGMTagManager::TagInfo *tinfo = (tag_handle > 0 ? 
                                     &(myManager->tagInfo[tag_handle]) :
                                     &(myManager->presetTagInfo[-tag_handle]));

      // store the name
    str_data.append(new CubitString(tinfo->tagName.c_str()));
    
      // store the length in bytes
    int_data.append(tinfo->tagLength);
    
      // now the data
      // store the raw memory interpreted as an array of ints, padded to a full int
    int tag_ints = tinfo->tagLength/4;
    if (tinfo->tagLength % 4 != 0) tag_ints++;
    
    int *tag_data = reinterpret_cast<int*>((*mit).second);
    for (int i = 0; i < tag_ints; i++)
      int_data.append(tag_data[i]);
  }

    // store the data on the csa
  CubitSimpleAttrib *csa = new CubitSimpleAttrib(&str_data, &dbl_data, &int_data);

  for (int i = 0; i < str_data.size(); i++)
    delete str_data.get_and_step();
  
  return csa;
}

void CATag::add_csa_data(CubitSimpleAttrib *csa_ptr) 
{
    // make sure it's a CATag
  static CubitString my_type("CA_TAG");
  if (csa_ptr->character_type() != my_type) 
    return;

  csa_ptr->int_data_list()->reset();
  csa_ptr->string_data_list()->reset();
  int num_attribs = *(csa_ptr->int_data_list()->get_and_step());

  int *tmp_data;
  
  for (int i = 0; i < num_attribs; i++) {

      // map the attrib name to a tag
    std::map<std::string,long>::iterator pos =
      myManager->tagNameMap.find(std::string(csa_ptr->string_data_list()->get()->c_str()));

    long thandle = 0;
    
    if (pos == myManager->tagNameMap.end()) {
        // tag doesn't exist - create one
      myManager->createTag(csa_ptr->string_data_list()->get()->c_str(),
                           *(csa_ptr->int_data_list()->get()), iBase_BYTES,
                           NULL, &thandle);
    }
    else thandle = (*pos).second;

    
    long tag_handle = thandle;

      // copy the ints to a temporary space we can get a ptr to...
    int int_length = *(csa_ptr->int_data_list()->get_and_step());
    if (int_length % 4 != 0) int_length++;
    tmp_data = (int*) malloc(int_length*sizeof(int));
    for (int j = 0; j < int_length; j++) 
      tmp_data[j] = *(csa_ptr->int_data_list()->get_and_step());

      // now actually set the data
    this->set_tag_data(tag_handle, tmp_data, true);
  }
}
    
void CATag::print() 
{
  std::cout << "This entity has " << tagData.size() << " tags.  Types are: " << std::endl;
  for (std::map<int,void*>::iterator mit = tagData.begin(); mit != tagData.end(); mit++) 
  {
    if ((*mit).first > 0)
      std::cout << myManager->tagInfo[(*mit).first].tagName << std::endl;
    else
      std::cout << myManager->presetTagInfo[-(*mit).first].tagName << std::endl;
  }
}

iBase_ErrorType CATag::get_tag_data(long tag_handle, void *tag_data) 
{
  assert(NULL != tag_data);

  CGMTagManager::TagInfo *tinfo = (tag_handle > 0 ? 
                                   &(myManager->tagInfo[tag_handle]) : 
                                   &(myManager->presetTagInfo[-tag_handle]));
  
    // check if this attribute has this tag
  std::map<int, void*>::iterator tdpos = tagData.find(tag_handle);
  if (tdpos == tagData.end()) {
    if (NULL != tinfo->defaultValue)
      memcpy(tag_data, tinfo->defaultValue, tinfo->tagLength);
    else {
      iGeom_setLastError( iBase_TAG_NOT_FOUND );
      return iBase_TAG_NOT_FOUND;
    }
    
  }
  
  else
    memcpy(tag_data, (*tdpos).second, tinfo->tagLength);

  RETURN(iBase_SUCCESS);
}
  
iBase_ErrorType CATag::set_tag_data(long tag_handle, const void *tag_data, 
                                     const bool can_shallow_copy)
{
  CGMTagManager::TagInfo *tinfo = (tag_handle > 0 ? 
                                   &(myManager->tagInfo[tag_handle]) : 
                                   &(myManager->presetTagInfo[-tag_handle]));
  
    // check if this attribute has this tag
  std::map<int, void*>::iterator tdpos = tagData.find(tag_handle);
  if (tdpos == tagData.end())
    tdpos = tagData.insert(tagData.end(),
                           std::pair<int,void*>(tag_handle, NULL));
    
  if (!can_shallow_copy) {
      // need to copy the data; might need to allocate first
    if ((*tdpos).second == NULL)
      (*tdpos).second = malloc(tinfo->tagLength);

    memcpy((*tdpos).second, tag_data, tinfo->tagLength);
  }
  else {
      // should shallow copy; might have to delete what's there already
    if ((*tdpos).second != NULL) free((*tdpos).second);
  
      // if shallow copying, caller is saying we can copy, so cast away const
    (*tdpos).second = const_cast<void*>(tag_data);
  }

  RETURN(iBase_SUCCESS);
}

void CATag::remove_tag(long tag_handle)
{
  tagData.erase(tag_handle);
}

CubitStatus CATag::update() 
{
  if (tagData.empty())
    this->delete_attrib(true);

  return CUBIT_SUCCESS;
}
