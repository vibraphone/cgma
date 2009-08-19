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
 
/**\file CATag.hpp
 *
 *\author Tim Tautges
 *\author Jason Kraftcheck
 *
 * Original file from SNL TSTT repository was named CATSTT.
 *
 * Renamed CATag and added to ANL ITAPS repository by J.Kraftcheck,
 * 2007-6-15
 */

#ifndef CA_TAG_HPP
#define CA_TAG_HPP

#include "iBase.h"
#include "CubitAttrib.hpp"
#include "DLIList.hpp"
#include "CubitDefines.h"

#include <string>
#include <map>
#include <vector>

#define ARRAY_IN_DECL(a, b) \
  a const* b, const int b ## _size  

#define ARRAY_INOUT_DECL(a, b) \
  a **b, int *b ## _allocated, int *b ## _size  

class RefEntity;
class RefGroup;
class CATag;
class CGMParallelComm;

class CGMTagManager 
{
public:

  friend class CATag;
  
  ~CGMTagManager();
  
  struct TagInfo
  {
    int tagLength;
    std::string tagName;
    int tagType;
    char *defaultValue;
    bool isActive;
  };
  
  static CubitAttrib* CATag_creator(RefEntity* entity, CubitSimpleAttrib *p_csa);

  iBase_ErrorType createTag (/*in*/ const char *tag_name,
                             /*in*/ const int tag_size,
                             /*in*/ const int tag_type,
                             /*in*/ char* default_value,
                             /*out*/ long *tag_handle);

  iBase_ErrorType destroyTag (/*in*/ const long tag_handle,
                              /*in*/ const bool forced);

  const char *getTagName (/*in*/ const long tag_handle);

  int getTagSize (/*in*/ const long tag_handle);

  long getTagHandle (/*in*/ const char *tag_name);

  int getTagType (/*in*/ const long tag_handle);

  iBase_ErrorType getArrData ( ARRAY_IN_DECL(RefEntity*, entity_handles),
                              /*in*/ const long tag_handle,
                              /*inout*/ ARRAY_INOUT_DECL(char, tag_value));

  iBase_ErrorType setArrData (/*in*/ ARRAY_IN_DECL(RefEntity*, entity_handles),
                              /*in*/ const long tag_handle,
                              /*in*/ const char *tag_values,
                              const int tag_values_size);

  iBase_ErrorType rmvArrTag (/*in*/ ARRAY_IN_DECL(RefEntity*, entity_handles),
                             /*in*/ const long tag_handle);

  iBase_ErrorType getAllTags (/*in*/ const RefEntity* entity_handle,
                              /*inout*/ ARRAY_INOUT_DECL(long, tag_handles));

  std::vector<RefGroup*> *pc_list(RefEntity *gentity, int list_no, 
                                  const bool create_if_missing);

  void pc_list(RefEntity *gentity, std::vector<RefGroup*> *&parents,
               std::vector<RefGroup*> *&children,
               const bool create_if_missing);

  void get_pc_groups(RefGroup *this_grp, const int p_or_c, const int num_hops,
                     std::vector<RefGroup *> &group_ptrs);

  iBase_ErrorType create_csa_tag(const char *tag_name, long *tag_handle);
  
  iBase_ErrorType set_csa_tag(RefEntity *this_ent,
                              long tag_handle,
                              CubitSimpleAttrib *csa_ptr);
  
  static inline CGMTagManager& instance()
  {
    static CGMTagManager static_instance;
    return static_instance;
  }

  inline std::vector<CGMParallelComm*>& get_pc_array() { return pc_array; }

private:
  CGMTagManager();

  int CATag_att_type;
  long pcTag;
  std::vector<TagInfo> tagInfo;
  static TagInfo* const presetTagInfo;
  static const int numPresetTag;
  std::map<std::string, long> tagNameMap;
  static const char *CATag_NAME;
  static const char *CATag_NAME_INTERNAL;
  RefGroup *interfaceGroup;

  // list of parallel comm
  std::vector<CGMParallelComm*> pc_array;

  bool getPresetTagData(const RefEntity *entity, const long tag_num, 
                        char *tag_value, int &tag_size);

  iBase_ErrorType setPresetTagData(RefEntity *entity, const long tag_num, 
                                   const char *tag_value, const int tag_size);
  
  CATag *get_catag(RefEntity *ent, 
                   const bool create_if_missing = false);

  long pc_tag(const bool create_if_missing = false);
  
  RefGroup *interface_group(const bool create_if_missing = true);
  
  CubitSimpleAttrib* get_simple_attrib( RefEntity* entity, const char* name );  
};

class CATag: public CubitAttrib
{
private:
  friend class CGMTagManager;
 
  std::map<int, void*> tagData;

  CGMTagManager *myManager;

  CATag(CGMTagManager *manager, RefEntity *owner);

  CATag(CGMTagManager *manager, RefEntity *owner, CubitSimpleAttrib *csa_ptr);
    //- create a CATag from a simple attribute

public:

  virtual ~CATag();

    //HEADER- RTTI and safe casting functions.
  virtual const type_info& entity_type_info() const
    { return typeid(CATag);}
    //R- The geometric modeler type
    //- This function returns the type of the geometric modeler.

  CubitStatus actuate() {return CUBIT_SUCCESS;}

  CubitStatus update();

  CubitStatus reset();

  CubitSimpleAttrib *split_owner() {return NULL;}

  void merge_owner(CubitAttrib *deletable_attrib) {}

  CubitSimpleAttrib* cubit_simple_attrib();

  int int_attrib_type() {return myManager->CATag_att_type;}

  void add_csa_data(CubitSimpleAttrib *csa_ptr);

  void print();

  iBase_ErrorType get_tag_data(long tag_num, void *tag_data);
  
  iBase_ErrorType set_tag_data(long tag_num, const void *tag_data,
                               const bool can_shallow_copy = false);

  void remove_tag(long tag_num);
};

#endif

