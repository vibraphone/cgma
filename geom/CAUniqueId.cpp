//- Class:          CAUniqueId
//- Owner:          Tim Tautges
//- Description:    Cubit Attribute for unique ids
//- Checked By:
//- Version:

#include "CAUniqueId.hpp"
#include "CubitSimpleAttrib.hpp"
#include "RefEntity.hpp"
#include "TDUniqueId.hpp"
#include "GSaveOpen.hpp"
#include "GeometryQueryTool.hpp"

DLIList<CAUniqueId *> CAUniqueId::allCAUniqueIds;
bool CAUniqueId::autoUniqueId = false;
UIDMap CAUniqueId::oldUIDToNewUID;

CubitAttrib* CAUniqueId_creator(RefEntity* entity, const CubitSimpleAttrib &p_csa)
{
  return new CAUniqueId(entity, p_csa);
}

CAUniqueId::~CAUniqueId()
{
  if (allCAUniqueIds.move_to(this))
    allCAUniqueIds.extract();
}

CAUniqueId::CAUniqueId(RefEntity* new_attrib_owner,
                               const CubitSimpleAttrib &csa_ptr)
        : CubitAttrib(new_attrib_owner)
{
  uniqueId = -1;

  if(!csa_ptr.isEmpty())
  {
    uniqueId = csa_ptr.int_data_list()[0];
  }
  allCAUniqueIds.append(this);
}

CubitStatus CAUniqueId::actuate()
{

  if (hasActuated == CUBIT_TRUE) return CUBIT_SUCCESS;

    // create a TDUniqueId for the entity, if it doesn't already
    // exist
  TDUniqueId *uid = (TDUniqueId *) attrib_owner()->get_TD(&TDUniqueId::is_unique_id);

  if (uid != NULL) {
      // check to make sure it's the same unique id
    if (uid->unique_id() != uniqueId) {
      PRINT_ERROR("Different unique id found for %s %d.\n",
                  attrib_owner()->class_name(), attrib_owner()->id());
      return CUBIT_FAILURE;
    }
  }
  else 
  {
    if( !GSaveOpen::performingUndo && 
         GeometryQueryTool::importingSolidModel &&
        !GeometryQueryTool::mergeGloballyOnImport)
    {
      //Is there an entity that already has this id?
      ToolDataUser *tdu = TDUniqueId::find_td_unique_id(uniqueId);
      if( tdu )
      {
        //is it already in the map
        UIDMap::iterator iter;
        iter = oldUIDToNewUID.find( uniqueId );
        
        if( iter != oldUIDToNewUID.end() )
          uniqueId = (*iter).second;
        else
        {
          int new_unique_id = TDUniqueId::generate_unique_id();
          UIDMap::value_type this_pair(uniqueId, new_unique_id);
          oldUIDToNewUID.insert(this_pair);
          uniqueId = new_unique_id;
        }
      }
    }

      // else make a new one
    uid = new TDUniqueId(attrib_owner(), uniqueId);
  }
  
  delete_attrib(CUBIT_TRUE);
  hasActuated = CUBIT_TRUE;
  
  return CUBIT_SUCCESS;
}

CubitStatus CAUniqueId::update()
{
  if (hasUpdated) return CUBIT_SUCCESS;
  
    // set the updated flag
  hasUpdated = CUBIT_TRUE;

    // if the owner has a unique id, save it, otherwise delete this one
  TDUniqueId *td_uid = (TDUniqueId *) attrib_owner()->get_TD(&TDUniqueId::is_unique_id);

  if (NULL == td_uid && autoUniqueId) 
    td_uid = new TDUniqueId(attrib_owner());

  if (td_uid == NULL) {
    delete_attrib(CUBIT_TRUE);
  }

  else {
    uniqueId = td_uid->unique_id();
    if (delete_attrib() == CUBIT_TRUE) delete_attrib(CUBIT_FALSE);
  }
  
  return CUBIT_SUCCESS;
}

CubitSimpleAttrib CAUniqueId::cubit_simple_attrib()
{
  return CubitSimpleAttrib(att_internal_name(), "", "",
                               &uniqueId);
}

CubitStatus CAUniqueId::actuate_all()
{
    //- actuate all the CAUI's on the list, then empty the list
  for (int i = allCAUniqueIds.size(); i > 0; i--) {
    CAUniqueId *cauid = allCAUniqueIds.get();
    if (cauid->actuate() == CUBIT_SUCCESS) allCAUniqueIds.extract();
    else allCAUniqueIds.step();
  }
  
  return CUBIT_SUCCESS;
}

void CAUniqueId::print()
{
    // print info on this attribute
  
  PRINT_INFO("CAUniqueId: owner = %s %d: uid=%d\n",
             attribOwnerEntity->class_name(), attribOwnerEntity->id(),
             uniqueId);
}

void CAUniqueId::clear_out_old_to_new_map()
{
  oldUIDToNewUID.clear(); 
}

