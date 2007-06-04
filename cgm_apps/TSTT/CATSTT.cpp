#include "CATSTT.hpp"
#include "RefEntity.hpp"
#include "RefEntityName.hpp"
#include "RefEntityFactory.hpp"
#include "CubitAttribManager.hpp"
#include "RefGroup.hpp"
#include "TDUniqueId.hpp"
#include "CGMApp.hpp"
#include "TSTTB_SNL.h"
#include "TSTTG_CGM.h"

#define CHECK_SIZE(array, type, this_size, retval)  \
  if (0 == array ## _allocated || array ## _allocated < this_size) {\
    if (NULL != array) free(array); \
    array = (type*)malloc(this_size*sizeof(type));\
    array ## _allocated=this_size;\
    if (NULL == array) {TSTTG_processError(TSTTB_MEMORY_ALLOCATION_FAILED, \
          "Couldn't allocate array.");return retval; }\
  }; \
  array ## _size = this_size;

#define TAG_CHECK_SIZE(array, allocated, size)  \
  if (allocated < size) {\
    if (NULL != array) free(array); \
    allocated=size; \
    array = (char*)malloc(allocated); \
    if (NULL == array) {TSTTG_processError(TSTTB_MEMORY_ALLOCATION_FAILED, \
          "Couldn't allocate array.");return CUBIT_FAILURE; }\
  }

const char *CGMTagManager::CATSTT_NAME = "TSTTB_Tag";
const char *CGMTagManager::CATSTT_NAME_INTERNAL = "TSTTB_TAG";

CGMTagManager *CGMTagManager::staticCGMTagManager = NULL;

CubitAttrib *CGMTagManager::CATSTT_creator(RefEntity* entity, CubitSimpleAttrib *p_csa)
{
  CATSTT *this_ca = new CATSTT(staticCGMTagManager, entity, p_csa);
  return this_ca;
}

CGMTagManager::CGMTagManager() 
    : interfaceGroup(NULL)
{
  staticCGMTagManager = this;

  pcTag = 0;
  
    // get the tag number for CATSTT
  DLIList<int> tag_types;
  int max_type = 0;
  CubitAttribManager *cam = CGMApp::instance()->attrib_manager();
  cam->get_registered_types(tag_types);
  for (int i = 0; i < tag_types.size(); i++) {
    int this_type = tag_types.get_and_step();
    max_type = (max_type < this_type ? this_type : max_type);
  }
  
  max_type++;
  CubitStatus status = cam->register_attrib_type(max_type, CATSTT_NAME, CATSTT_NAME_INTERNAL,
                                                 &CGMTagManager::CATSTT_creator, false, true, 
                                                 true, true, true, false);
  if (CUBIT_FAILURE == status) {
    TSTTG_processError(TSTTB_FAILURE, "Couldn't create cgm attribute for tags.");
  }
  else
    CATSTT_att_type = max_type;

    // create the 0th tag, since this index represents no tag
  TagInfo tmp_info = {0, std::string("NULL tag"), TSTTG_BYTES, NULL, false};
  tagInfo.push_back(tmp_info);
  presetTagInfo.push_back(tmp_info);

    // create preset tags, for CGM attributes we want to be visible as tags
    // name
  tmp_info.tagLength = sizeof(std::string);
  tmp_info.tagName = std::string("ENTITY_NAME");
  tmp_info.defaultValue = NULL;
  tmp_info.tagType = TSTTG_BYTES;
  tmp_info.isActive = true;
  presetTagInfo.push_back(tmp_info);
  tagNameMap[std::string("ENTITY_NAME")] = -1;
    // id
  tmp_info.tagLength = sizeof(int);
  tmp_info.tagName = std::string("GLOBAL_ID");
  tmp_info.tagType = TSTTG_INTEGER;
  presetTagInfo.push_back(tmp_info);
  tagNameMap[std::string("GLOBAL_ID")] = -2;
    // uid
  tmp_info.tagName = std::string("UNIQUE_ID");
  tmp_info.tagType = TSTTG_INTEGER;
  presetTagInfo.push_back(tmp_info);
  tagNameMap[std::string("UNIQUE_ID")] = -3;
}

CGMTagManager::~CGMTagManager() 
{
  staticCGMTagManager = NULL;
}

int CGMTagManager::pc_tag(const bool create_if_missing) 
{
  if (0 == pcTag && create_if_missing) {
    pcTag = getTagHandle("__cgm_parent_child_tag");
    if (0 == pcTag) {
        // ok, one doesn't exist, create one
      void *def_val[] = {NULL, NULL};
      
      createTag("__cgm_parent_child_tag", 2*sizeof(std::vector<RefGroup*>*),
                TSTTG_BYTES, (char*)def_val, &pcTag);
    }
  }
  
  return pcTag;
}

CubitStatus CGMTagManager::createTag (/*in*/ const char *tag_name,
                                      /*in*/ const int tag_length,
                                      /*in*/ const int tag_type,
                                      /*in*/ char* default_value,
                                      /*out*/ int *tag_handle)
{
  std::string tmp_name(tag_name);
  std::map<std::string,int>::iterator mit = tagNameMap.find(tmp_name);
  if (mit != tagNameMap.end()) {
    *tag_handle = mit->second;
    TSTTG_processError(TSTTB_TAG_ALREADY_EXISTS, "Tag already created.");
    return CUBIT_SUCCESS;
  }
    
  TagInfo tmp_info = {tag_length, std::string(tag_name), tag_type, NULL, true};
  tagInfo.push_back(tmp_info);
  *tag_handle = tagInfo.size()-1;
  if (default_value != NULL) {
    tagInfo[*tag_handle].defaultValue = (char *) malloc(tag_length);
    memcpy(tagInfo[*tag_handle].defaultValue, default_value, tag_length);
  }

    // put the name and handle into the map too
  tagNameMap[std::string(tag_name)] = *tag_handle;

  return CUBIT_SUCCESS;
}


CubitStatus CGMTagManager::destroyTag (/*in*/ const int tag_handle,
                                       /*in*/ const bool forced)
{
  if (!forced) {
      // see whether this tag is still assigned anywhere
      // not implemented yet
    return CUBIT_FAILURE;
  }
  
    // if we got here, assume we can delete it
  TagInfo *tinfo = (tag_handle >= 0 ? &tagInfo[tag_handle] : &presetTagInfo[-tag_handle]);
  tinfo->isActive = false;

  return CUBIT_SUCCESS;
}

const char *CGMTagManager::getTagName (/*in*/ const int tag_handle)
{
  return (tag_handle >= 0 ? tagInfo[tag_handle].tagName.c_str() : 
          presetTagInfo[-tag_handle].tagName.c_str());
}

int CGMTagManager::getTagSize (/*in*/ const int tag_handle)
{
  return (tag_handle >= 0 ? 
          tagInfo[tag_handle].tagLength : 
          presetTagInfo[-tag_handle].tagLength);
}

int CGMTagManager::getTagHandle (/*in*/ const char *tag_name)
{
  std::map<std::string,int>::iterator it =
    tagNameMap.find(std::string(tag_name));
  if (it != tagNameMap.end())
    return (*it).second;
  else return 0;
}

int CGMTagManager::getTagType (/*in*/ const int tag_handle) 
{
  return (tag_handle >= 0 ? 
          tagInfo[tag_handle].tagType : 
          presetTagInfo[-tag_handle].tagType);
}

CubitStatus CGMTagManager::getArrData (/*in*/ const RefEntity**entity_handles, 
                                       const int entity_handles_size,
                        /*in*/ const int tag_handle,
                        /*inout*/ ARRAY_INOUT_DECL(char, tag_value))
{
  TagInfo *tinfo = (tag_handle >= 0 ? &tagInfo[tag_handle] : &presetTagInfo[-tag_handle]);
  int tag_size = tinfo->tagLength;
    // either way, we have to have that many bytes when we leave this function
  *tag_value_size = entity_handles_size*tinfo->tagLength;
  TAG_CHECK_SIZE(*tag_value, *tag_value_allocated, *tag_value_size);
  char *val_ptr = *tag_value;
  if (tag_handle < 0) {
    for (int i = 0; i < entity_handles_size; i++) {
      if (NULL == entity_handles[i])
        getPresetTagData(interface_group(), tag_handle, val_ptr, tinfo->tagLength);
      else
        getPresetTagData(entity_handles[i], tag_handle, val_ptr, tinfo->tagLength);
      val_ptr += tinfo->tagLength;
    }
    return CUBIT_SUCCESS;
  }

  CubitStatus tmp_result = CUBIT_SUCCESS;

  for (int i = 0; i < entity_handles_size; i++) {
      // ok to cast away const-ness because "false" passed in for create_if_missing
    RefEntity *this_ent = (NULL == entity_handles[i] ? interface_group() :
                           const_cast<RefEntity*>(entity_handles[i]));
    CATSTT *catstt = get_catstt(this_ent);
    if (NULL != catstt) tmp_result = catstt->get_tag_data(tag_handle, val_ptr);
    else if (NULL != tinfo->defaultValue)
      memcpy(val_ptr, tinfo->defaultValue, tinfo->tagLength);
    else
      tmp_result = CUBIT_FAILURE;

    val_ptr += tag_size;
  }

  return CUBIT_SUCCESS;
}

CubitStatus CGMTagManager::setArrData (/*in*/ ARRAY_IN_DECL(RefEntity*, entity_handles),
                        /*in*/ const int tag_handle,
                        /*in*/ const char *tag_values, const int tag_values_size)
{
  TagInfo *tinfo = (tag_handle >= 0 ? &tagInfo[tag_handle] : &presetTagInfo[-tag_handle]);
  int tag_size = tinfo->tagLength;
  
  const char *val_ptr = tag_values;

  if (tag_handle < 0) {
    for (int i = 0; i < entity_handles_size; i++) {
      if (NULL == entity_handles[i])
        setPresetTagData(interface_group(), tag_handle, val_ptr, tag_size);
      else
        setPresetTagData(entity_handles[i], tag_handle, val_ptr, tag_size);
      val_ptr += tag_size;
    }
    return CUBIT_SUCCESS;
  }

  for (int i = 0; i < entity_handles_size; i++) {
    RefEntity *this_ent = (NULL == entity_handles[i] ? interface_group() : 
                           entity_handles[i]);
    CATSTT *catstt = get_catstt(this_ent, true);
    catstt->set_tag_data(tag_handle, val_ptr);
    val_ptr += tag_size;
  }

  return CUBIT_SUCCESS;
}

CubitStatus CGMTagManager::rmvArrTag (/*in*/ ARRAY_IN_DECL(RefEntity*, entity_handles),
                       /*in*/ const int tag_handle)
{
  for (int i = 0; i < entity_handles_size; i++) {
    CATSTT *catstt = get_catstt((entity_handles[i] == NULL ? 
                                 interface_group() : entity_handles[i]));
    if (NULL != catstt) catstt->remove_tag(tag_handle);
  }

  return CUBIT_SUCCESS;
}

CubitStatus CGMTagManager::getAllTags (/*in*/ const RefEntity* entity_handle,
                        /*inout*/ ARRAY_INOUT_DECL(int, tag_handles))
{
  int i = 0, uid = 0, tag_size;
  char *uid_ptr = (char*) &uid;
  bool has_uid = getPresetTagData(entity_handle, -3, uid_ptr, tag_size);
  int num_tags = (has_uid ? 3 : 2);

  RefEntity *this_ent = (NULL == entity_handle ? interface_group() : 
                         const_cast<RefEntity*>(entity_handle));
  
    // const-cast because we're passing in false for create_if_missing
  CATSTT *catstt = get_catstt(this_ent);
  if (NULL != catstt) {
      // need to check whether entity has a uid
    num_tags += catstt->tagData.size();
    CHECK_SIZE(*tag_handles, int, num_tags, CUBIT_FAILURE);
    for (std::map<int,void*>::iterator tag_it = catstt->tagData.begin();
         tag_it != catstt->tagData.end(); tag_it++)
      (*tag_handles)[i++] = (*tag_it).first;
  }
  else {
    CHECK_SIZE(*tag_handles, int, num_tags, CUBIT_FAILURE);
  }
  (*tag_handles)[i++] = -1;
  (*tag_handles)[i++] = -2;
  if (has_uid) (*tag_handles)[i++] = -3;

  return CUBIT_SUCCESS;
}

CATSTT *CGMTagManager::get_catstt(RefEntity *ent, 
                                  const bool create_if_missing) 
{
  CubitAttrib *this_attrib = ent->get_cubit_attrib(CATSTT_att_type, create_if_missing);
  if (NULL != this_attrib)
    return dynamic_cast<CATSTT*>(this_attrib);
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
  
void CGMTagManager::setPresetTagData(RefEntity *entity, 
                                     const int tag_handle, 
                                     const char *tag_value, 
                                     const int tag_size) 
{
  const std::string *this_name;
  switch (-tag_handle) {
    case 1:
        // entity name
      if (sizeof(std::string) != tag_size) {
        TSTTG_processError(TSTTB_INVALID_ARGUMENT, "Tag of type 'name' is the wrong size.");
        return;
      }
      this_name = reinterpret_cast<const std::string*>(tag_value);
      entity->entity_name(CubitString(this_name->c_str()));
      return;
    case 2:
        // entity id
      TSTTG_processError(TSTTB_NOT_SUPPORTED, "Can't set id of entities with this implementation.");
      return;
    case 3:
        // unique id
      TSTTG_processError(TSTTB_NOT_SUPPORTED, "Can't set unique id of entities with this implementation.");
      return;
  }
}
        
bool CGMTagManager::getPresetTagData(const RefEntity *entity, 
                                     const int tag_handle, 
                                     char *tag_value, 
                                     int &tag_size) 
{
  std::string *this_name;
  int *this_id;
  int *this_uid;
  
  switch (-tag_handle) {
    case 1:
        // entity name
      tag_size = sizeof(std::string);
      this_name = reinterpret_cast<std::string*>(tag_value);
      *this_name = std::string(entity->entity_name().c_str());
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
  CubitStatus result =
    getArrData(const_cast<const RefEntity**>(&gentity), 1, pc_tag(create_if_missing), &dum_ptr, 
               &dum, &dum);
  assert(CUBIT_SUCCESS == result);
  
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
  CubitStatus result =
    getArrData(const_cast<const RefEntity**>(&gentity), 1, pc_tag(), 
               &dum_ptr, &dum, &dum);
  assert(CUBIT_SUCCESS == result);

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
      assert(CUBIT_SUCCESS == result);
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

CATSTT::~CATSTT() 
{
  for (std::map<int, void*>::iterator 
         mit = tagData.begin(); mit != tagData.end(); mit++)
    if (NULL != (*mit).second) free ((*mit).second);
}


CATSTT::CATSTT(CGMTagManager *manager, RefEntity *entity) 
    : CubitAttrib(entity), myManager(manager)
{
}

CATSTT::CATSTT(CGMTagManager *manager, RefEntity *owner, CubitSimpleAttrib *csa_ptr) 
    : CubitAttrib(owner), myManager(manager)
{
  if (NULL != csa_ptr) add_csa_data(csa_ptr);
}

CubitStatus CATSTT::reset()
{
  for (std::map<int, void*>::iterator 
         mit = tagData.begin(); mit != tagData.end(); mit++)
    if (NULL != (*mit).second) free ((*mit).second);

  tagData.clear();

  return CUBIT_SUCCESS;
}

CubitSimpleAttrib* CATSTT::cubit_simple_attrib() 
{
    //if (tagData.size() == 0) return NULL;
  
  DLIList<int> int_data;
  DLIList<CubitString*> str_data;
  DLIList<double> dbl_data;

  str_data.append(new CubitString(myManager->CATSTT_NAME_INTERNAL));

    // int data first gets the # tags on this entity
  int_data.append(tagData.size());

    // for each tag:
  for (std::map<int, void*>::iterator 
         mit = tagData.begin(); mit != tagData.end(); mit++) {
    int tag_handle = (*mit).first;
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

void CATSTT::add_csa_data(CubitSimpleAttrib *csa_ptr) 
{
    // make sure it's a CATSTT
  static CubitString my_type("CA_TSTT");
  if (csa_ptr->character_type() != my_type) 
    return;

  csa_ptr->int_data_list()->reset();
  csa_ptr->string_data_list()->reset();
  int num_attribs = *(csa_ptr->int_data_list()->get_and_step());

  int *tmp_data;
  
  for (int i = 0; i < num_attribs; i++) {

      // map the attrib name to a tag
    std::map<std::string,int>::iterator pos =
      myManager->tagNameMap.find(std::string(csa_ptr->string_data_list()->get()->c_str()));

    int thandle = 0;
    
    if (pos == myManager->tagNameMap.end()) {
        // tag doesn't exist - create one
      myManager->createTag(csa_ptr->string_data_list()->get()->c_str(),
                           *(csa_ptr->int_data_list()->get()), TSTTG_BYTES,
                           NULL, &thandle);
    }
    else thandle = (*pos).second;

    
    int tag_handle = reinterpret_cast<int>(thandle);

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
    
void CATSTT::print() 
{
  std::cout << "This entity has " << tagData.size() << " tags.  Types are: " << std::endl;
  for (std::map<int,void*>::iterator mit = tagData.begin(); mit != tagData.end(); mit++) 
  {
    if ((*mit).first < 0)
      std::cout << myManager->tagInfo[-(*mit).first].tagName << std::endl;
    else
      std::cout << myManager->tagInfo[(*mit).first].tagName << std::endl;
  }
}

CubitStatus CATSTT::get_tag_data(int tag_handle, void *tag_data) 
{
  assert(NULL != tag_data);

  CGMTagManager::TagInfo *tinfo = (tag_handle >= 0 ? 
                                   &(myManager->tagInfo[tag_handle]) : 
                                   &(myManager->presetTagInfo[-tag_handle]));
  
    // check if this attribute has this tag
  std::map<int, void*>::iterator tdpos = tagData.find(tag_handle);
  if (tdpos == tagData.end()) {
    if (NULL != tinfo->defaultValue)
      memcpy(tag_data, tinfo->defaultValue, tinfo->tagLength);
    else return CUBIT_FAILURE;
  }
  
  else
    memcpy(tag_data, (*tdpos).second, tinfo->tagLength);

  return CUBIT_SUCCESS;
}
  
CubitStatus CATSTT::set_tag_data(int tag_handle, const void *tag_data, 
                                 const bool can_shallow_copy)
{
  CGMTagManager::TagInfo *tinfo = (tag_handle >= 0 ? 
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

  return CUBIT_SUCCESS;
}

void CATSTT::remove_tag(int tag_handle)
{
  tagData.erase(tag_handle);
}

CubitStatus CATSTT::update() 
{
  if (tagData.empty())
    this->delete_attrib(true);

  return CUBIT_SUCCESS;
}
