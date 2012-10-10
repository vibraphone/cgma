#include <stdlib.h>
#include "TDUniqueId.hpp"
#include "RandomMersenne.hpp"
#include "DLIList.hpp"
#include "ToolDataUser.hpp"
#include "CubitAttribUser.hpp"
#include "CubitAttrib.hpp"
#include "CADefines.hpp"
#include "RefEntity.hpp"
#include "TopologyEntity.hpp"
#include "CastTo.hpp"
#include "GeometryModifyTool.hpp"
#include "GeometryQueryTool.hpp"
#include "GSaveOpen.hpp"
#include "CAUniqueId.hpp"

#include <stdio.h>
#include <time.h>

extern int init_genrand(unsigned long s);
extern long genrand_int31(void);

TDUIDList TDUniqueId::uniqueIdList;
COPYUIDMap TDUniqueId::mapForCopying;

int TDUniqueId::generate_unique_id()
{
    int result;

// initialize the random number generator, using time as a seed; do it here
// instead of as a static variable of the class to add a little more randomness
// to the seed (which uses the time)
    static RandomMersenne rng(time(NULL));

      // generate random no. and check for duplicate; keep generating until
      // we find a non-duplicate; initialize found to true so that we can
      // test for duplicates in the loop below
    bool found = false;
    TDUIDList::iterator found_value;
    do {
      result = (int)(rng.genrand_int31());
      found_value = uniqueIdList.lower_bound(result);
      if(found_value != uniqueIdList.end())
      {
        found = ((*found_value).first == result);
        if (found)
          DIAGNOSTIC("Duplicate UID found - very unusual!\n");
      }
    }
    while (found);

    return result;
}


TDUniqueId::TDUniqueId(ToolDataUser *owner, const int id)
{
  ownerEntity = owner;

  uniqueId = id == 0 ? generate_unique_id() : id;
  TDUIDList::value_type this_pair(uniqueId, this);
  unique_id_list().insert(this_pair);

  ownerEntity->add_TD(this);
  

    // update the attribute if this is a CAU
  CubitAttribUser *cau = CAST_TO(owner, CubitAttribUser);
  if (cau) {
    CubitAttrib *attrib = cau->get_cubit_attrib(CA_UNIQUE_ID);
    attrib->has_updated(CUBIT_FALSE);
    attrib->update();
  }
}

int TDUniqueId::is_unique_id(const ToolData* td)
{
  return (CAST_TO(const_cast<ToolData*>(td), TDUniqueId) != NULL);
}

TDUniqueId::~TDUniqueId()
{
    // remove this from the list
	//- With multimaps, keys do not have to be unique.  With the erase(Key) function,
	//- the multimap code will remove all entries with that key from the multimap.  
	//- The erase(iterator) function only removes the information at that iteration.
	//- Here, we search for the proper entry if there is more than one entry for a 
	//- particular key.  Increment the iterator and compare the element to 'this'.  If 
	//- it is a match, erase it by iterator.  This will leave the other entry with
	//- the similar key in the list, which is necessary for save/restore, and remove the 
	//- deactivated geometry due to the merge.

  std::pair<TDUIDList::iterator, TDUIDList::iterator>
    bounds_pair = unique_id_list().equal_range(uniqueId);

  TDUIDList::iterator it;
  for(it = bounds_pair.first; it != bounds_pair.second; ++it)
  {
    if(it->second == this)
    {
      unique_id_list().erase( it );
      break;
    }
  }
}
    
int TDUniqueId::get_unique_id(ToolDataUser *owner,
                              const CubitBoolean create_new)
{
  assert(owner != 0);
  TDUniqueId *uid = (TDUniqueId *) owner->get_TD(&TDUniqueId::is_unique_id);
  if (!uid && create_new == CUBIT_TRUE) {
    uid = new TDUniqueId(owner);
  }

  return (uid ? uid->unique_id() : 0);
}

int TDUniqueId::find_td_unique_id(const int temp_id,
                                  DLIList<ToolDataUser*> &td_list,
                                  const RefEntity *related_entity)
{
  td_list.clean_out();

  int unique_id = temp_id;

  //if we are not doing an undo and importing and merging within a file...
  if( !GSaveOpen::performingUndo && 
       GeometryQueryTool::importingSolidModel &&
      !GeometryQueryTool::mergeGloballyOnImport)
  {
    //see if the old id maps to a new id...if so, use the new id
    UIDMap old_uid_to_new_uid_map = CAUniqueId::get_old_to_new_uid_map();
    
    UIDMap::iterator iter;
    iter = old_uid_to_new_uid_map.find( unique_id );
    
    if( iter != old_uid_to_new_uid_map.end() )
      unique_id = (*iter).second;
  }
  
  std::pair<TDUIDList::iterator, TDUIDList::iterator> 
    bounds_pair = unique_id_list().equal_range(unique_id);

  TDUIDList::iterator
    it = bounds_pair.first, upper = bounds_pair.second;
 
  if(it == unique_id_list().end())
    return 0;

  if ((*it).first == unique_id ) {
      // the lower bound key is equal to unique_id, so this id is in the list

      // look for duplicate id's, return one that's directly related
      // get all td's with that id
    for (; it != upper; it++) 
    {
      bool related = true;
      ToolDataUser *temp_tdu = (*it).second->owner_entity();
      if (NULL != related_entity) {
        TopologyEntity *topo_ent = CAST_TO(temp_tdu, TopologyEntity);
        RefEntity* temp_entity = const_cast<RefEntity*>(related_entity);
        if (!topo_ent ||
            !topo_ent->is_directly_related(CAST_TO(temp_entity, TopologyEntity))) 
          related = false;
      }
      if (related) td_list.append(temp_tdu);
    }
  }

  return td_list.size();
}

ToolDataUser *TDUniqueId::find_td_unique_id(const int temp_id,
                                            const RefEntity *related_entity)
{
  //ToolDataUser *tdu = NULL;
  DLIList<ToolDataUser*> td_list;
  find_td_unique_id(temp_id, td_list, related_entity);
  if (0 == td_list.size()) return NULL;
  else {
    td_list.reset();
    return td_list.get();
  }
}

int TDUniqueId::get_unique_id_for_copy( int original_id )
{
  //check in the map if the original is already there
  //if so, get the corresponding id in the map (this is the new 
  //unique id on the copy)
  COPYUIDMap::iterator iter;
  iter = mapForCopying.find( original_id );
  if( iter != mapForCopying.end() ) 
    return (*iter).second;
   
  //else, a) generate another unique id b) put the original id in the map with
  //the new id directly across from it 
  int new_uid = generate_unique_id();
  COPYUIDMap::value_type tmp_pair( original_id, new_uid );
  mapForCopying.insert( tmp_pair );
  return new_uid;
}
          
  
void TDUniqueId::clear_copy_map()
{
  //remove TDs off of entities
  COPYUIDMap::iterator iter = mapForCopying.begin();
  for(; iter != mapForCopying.end(); iter++ )
  {
    DLIList<ToolDataUser*> tool_data_users;
    int num_tool_datas = find_td_unique_id( (*iter).second, tool_data_users );
    
    for( ; num_tool_datas--; ) 
    {
      ToolDataUser *td_user = tool_data_users.get_and_step();

      //make sure the pointer isn't null
      ToolData *tool_data = NULL;
      if(td_user)
        tool_data = td_user->remove_TD(TDUniqueId::is_unique_id);
      
      //delete the TDUniqueId as well
      if( tool_data )
        delete tool_data;

      td_user = find_td_unique_id( (*iter).second );
    } 
  } 

  //clear out the map 
  if (mapForCopying.empty()) return;
    mapForCopying.clear();

}
 
int TDUniqueId::unique_id()
{
  //if we're doing a copy, don't use the Unique id on the original.
  //generate another one
  if( GeometryModifyTool::instance()->get_copy_entity() )
  {
    int unique_id_for_copy = get_unique_id_for_copy( uniqueId );
    return unique_id_for_copy;
  }
  else 
    return uniqueId;
}
