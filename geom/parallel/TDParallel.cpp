#include <stdlib.h>
#include "TDParallel.hpp"
#include "DLIList.hpp"
#include "ToolDataUser.hpp"
#include "CubitAttribUser.hpp"
#include "CubitAttrib.hpp"
#include "RefEntity.hpp"
#include "TopologyEntity.hpp"
#include "CastTo.hpp"
#include "Body.hpp"
#include "TDUniqueId.hpp"
#include "CADefines.hpp"
#ifdef USE_MPI
#include "ProcData.hpp"
#include "CGMMemFile.hpp"
#endif

#include <stdio.h>
#include <time.h>
  
TDParallel::TDParallel(ToolDataUser *owner, DLIList<int> *id_list,
		       DLIList<int> *shared_list, int id,
		       int counter_id, int merge_id, int unique_id)
{
  CubitStatus status;
  ownerEntity = owner;

  if (id_list == NULL) {
    status = set_body_list();
    procId = id;
    counterProcId = counter_id;
    mergeId = merge_id;
    uniqueId = unique_id;
  }
  else {
    int size = id_list->size();
    id_list->reset();
    bodyUniqueIdList.clean_out();
    for (int i = 0; i < size; i++) {
      bodyUniqueIdList.append(id_list->get_and_step());
    }
    procId = id;
    counterProcId = counter_id;
    mergeId = merge_id;
    uniqueId = unique_id;
  }

  //cout << "shared_list=" << shared_list << endl;
  if (shared_list != NULL) {
    int shared_size = shared_list->size();
    //cout << "shared_size=" << shared_size << endl;
    shared_list->reset();
    sharedIdList.clean_out();
    for (int i = 0; i < shared_size; i++) {
      sharedIdList.append(shared_list->get_and_step());
    }
  }

  set_local_non_local_list();
  ownerEntity->add_TD(this);
  
    // update the attribute if this is a CAU
  CubitAttribUser *cau = CAST_TO(owner, CubitAttribUser);
  
  if (cau) {
    CubitAttrib *attrib = cau->get_cubit_attrib(CA_BODIES);
    attrib->has_updated(CUBIT_FALSE);
    attrib->update();
  }
}

int TDParallel::is_parallel(const ToolData* td)
{
  //cout << "td=" << td << endl;
  //cout << "const_cast<ToolData*>(td)=" << temp << endl;
  //cout << "CAST_TO(const_cast<ToolData*>(td), TDParallel)="
  //   << dynamic_cast<TDParallel*>(temp) << endl;
  return (CAST_TO(const_cast<ToolData*>(td), TDParallel) != NULL);
}

TDParallel::~TDParallel()
{
}

void TDParallel::set_local_non_local_list()
{
  DLIList<int> *body_to_proc = NULL;
#ifdef USE_MPI
  body_to_proc = CGMMemFile::get_body_to_proc();
#endif
  if (body_to_proc != NULL) {
    int size = body_to_proc->size();
    if (size > 0) {
      int id;
      int myId = -1;
#ifdef USE_MPI
      myId = ProcData::instance()->myRank;
#endif
      int size1 = bodyUniqueIdList.size();
      localList.clean_out();
      nonLocalList.clean_out();
      bodyUniqueIdList.reset();
      
      for (int i = 0; i < size1; i++) {
	id = bodyUniqueIdList.get_and_step();
	body_to_proc->reset();
	for (int j = 0; j < size; j += 2) {
	  if (body_to_proc->get_and_step() == id) {
	    if (body_to_proc->get_and_step() == myId) {
	      localList.append(id);
	    }
	    else {
	      nonLocalList.append(id);
	    }
	  }
	  else {
	    body_to_proc->get_and_step();
	  }
	}
      }
    }
  }
}

DLIList<int>* TDParallel::body_unique_id_list()
{
  return &bodyUniqueIdList;
}

CubitStatus TDParallel::set_body_list()
{  
  // get the bodies containing ownerEntity
  Body *body;
  CubitStatus status;
  DLIList<Body*> body_list;
  int i;

  status = (dynamic_cast<TopologyEntity*> (ownerEntity))->bodies(body_list);

  // write the data for the bodies to this TD
  for(i = body_list.size(); i > 0; i--)
  {
    // get the body which gets assigned to this TDParallel
    body = body_list.get_and_step();

    // append to this TDParallel the id of the body;
    //bodyIdList.append(body->id());
    bodyUniqueIdList.append(TDUniqueId::get_unique_id(body));
  }

  return status;
}

void TDParallel::set_shared_list(DLIList<int> list) {
  sharedIdList.clean_out();
  list.reset();
  for (int i = list.size(); i > 0; i--)
    sharedIdList.append(list.get_and_step());
}

bool TDParallel::is_shared(int id) {
  return sharedIdList.is_in_list(id);
}
