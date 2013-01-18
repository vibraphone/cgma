//- Class: TDParallel
//- Description: This tool data generates local and non-local bodies 
//-              which the entity is part of
//- Author: Hong-Jun Kim
//- Checked By: Tim Tautges
//- Version:

#ifndef TD_PARALLEL
#define TD_PARALLEL

#include "ToolData.hpp"
#include "MemoryManager.hpp"
#include "CubitDefines.h"
#include "DLIList.hpp"
#ifdef PARALLEL
#include "CGMMemFile.hpp"
#endif

template <class X> class DLIList;
class ToolDataUser;
class RefEntity;

class TDParallel : public ToolData
{
private:

  DLIList<int> bodyUniqueIdList;
  DLIList<int> bodyToProcList;
  DLIList<int> localList;
  DLIList<int> nonLocalList;
  DLIList<int> sharedIdList;

  ToolDataUser *ownerEntity;
    //- back pointer to the owning entity

  int procId;
  //- processor id which takes in charge of meshing this entity

  int counterProcId;
  //- processor id by which this entity is shared

  int mergeId;
  //- unique id merged

  int uniqueId;
  //- unique id

public:
  
  TDParallel(ToolDataUser *owner = NULL, DLIList<int> *id_list = NULL,
	     DLIList<int> *shared_list = NULL, int id = -1,
	     int counter_id = -1, int merge_id = -1, int unique_id = -1);
  
  virtual ~TDParallel();
    //-constructor and destructor

  static int is_parallel(const ToolData* td);
  
  DLIList<int>* body_unique_id_list();
    //- return bodyUniqueIdList;

  DLIList<int>* shared_id_list();

  DLIList<int>* local_body_list();

  DLIList<int>* non_local_body_list();

  void set_local_non_local_list();

  ToolDataUser *owner_entity() {return ownerEntity;};

  void owner_entity(ToolDataUser *owner) {ownerEntity = owner;};
    //- get/set functions for ownerEntity
  
  CubitStatus set_body_list();

  int get_proc_id();

  void set_proc_id(int id);
  
  int get_counter_proc_id();

  void set_counter_proc_id(int id);

  int get_merge_id();

  void set_merge_id(int id);

  int get_unique_id();

  void set_unique_id(int id);

  void set_shared_list(DLIList<int> list);

  int get_n_shared_list();

  void add_shared_id(int id);

  bool is_shared(int id);
};

inline DLIList<int>* TDParallel::local_body_list() {
  return &localList;
}

inline DLIList<int>* TDParallel::non_local_body_list() {
  return &nonLocalList;
}

inline DLIList<int>* TDParallel::shared_id_list() {
  return &sharedIdList;
}

inline int TDParallel::get_proc_id() {
  return procId;
}

inline void TDParallel::set_proc_id(int id) {
  procId = id;
}

inline int TDParallel::get_counter_proc_id() {
  return counterProcId;
}

inline void TDParallel::set_counter_proc_id(int id) {
  counterProcId = id;
}

inline int TDParallel::get_merge_id() {
  return mergeId;
}

inline void TDParallel::set_merge_id(int id) {
  mergeId = id;
}

inline int TDParallel::get_unique_id() {
  return uniqueId;
}

inline void TDParallel::set_unique_id(int id) {
  uniqueId = id;
}

inline int TDParallel::get_n_shared_list() {
  return sharedIdList.size();
}

inline void TDParallel::add_shared_id(int id) {
  sharedIdList.append_unique(id);
}

#endif 














