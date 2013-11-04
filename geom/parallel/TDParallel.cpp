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

#include <stdio.h>
#include <time.h>

TDParallel::TDParallel(ToolDataUser *owner, DLIList<int> *shared_bodies,
		       DLIList<int> *shared_procs, DLIList<int> *ghost_procs,
                       int unique_id, int interface)
  : m_uniqueId(unique_id), m_interface(interface)
{
  CubitStatus status;
  ownerEntity = owner;

  // shared bodies
  if (shared_bodies == NULL) {
    status = set_shared_body_list();
  }
  else {
    int size = shared_bodies->size();
    shared_bodies->reset();
    m_sharedBodyList.clean_out();
    for (int i = 0; i < size; i++) {
      m_sharedBodyList.append(shared_bodies->get_and_step());
    }
  }

  if (shared_procs != NULL) {
    int shared_size = shared_procs->size();
    shared_procs->reset();
    m_sharedProcList.clean_out();
    for (int i = 0; i < shared_size; i++) {
      m_sharedProcList.append(shared_procs->get_and_step());
    }
  }

  if (ghost_procs != NULL) {
    int ghost_size = ghost_procs->size();
    ghost_procs->reset();
    m_ghostProcList.clean_out();
    for (int i = 0; i < ghost_size; i++) {
      m_ghostProcList.append(ghost_procs->get_and_step());
    }
  }

  //set_local_non_local_body_list();
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
  return (CAST_TO(const_cast<ToolData*>(td), TDParallel) != NULL);
}

TDParallel::~TDParallel()
{
}

DLIList<int>* TDParallel::get_shared_body_list()
{
  return &m_sharedBodyList;
}

CubitStatus TDParallel::set_shared_body_list()
{  
  // get the bodies containing ownerEntity
  Body *body;
  DLIList<Body*> body_list;
  int i;

  CubitStatus result = (dynamic_cast<TopologyEntity*> (ownerEntity))->bodies(body_list);
  if (CUBIT_SUCCESS != result) return result;

  // write the data for the bodies to this TD
  for(i = body_list.size(); i > 0; i--)
  {
    // get the body which gets assigned to this TDParallel
    body = body_list.get_and_step();

    // append to this TDParallel the id of the body;
    m_sharedBodyList.append(TDUniqueId::get_unique_id(body));
  }

  return CUBIT_SUCCESS;
}

void TDParallel::set_shared_proc_list(DLIList<int> list)
{
  m_sharedProcList.clean_out();
  list.reset();
  for (int i = list.size(); i > 0; i--) {
    m_sharedProcList.append(list.get_and_step());
  }
}
