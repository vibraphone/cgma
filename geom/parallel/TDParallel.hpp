//- Class: TDParallel
//- Description: This tool data generates local and non-local bodies 
//-              which the entity is part of

#ifndef TD_PARALLEL
#define TD_PARALLEL

#include "ToolData.hpp"
#include "MemoryManager.hpp"
#include "CubitDefines.h"
#include "DLIList.hpp"

template <class X> class DLIList;
class ToolDataUser;
class RefEntity;

class TDParallel : public ToolData
{
private:

  // bodies sharing this entity
  DLIList<int> m_sharedBodyList;

  // for body to processor mapping
  DLIList<int> m_bodyToProcList;

  // bodies sharing this entity in the same processor
  DLIList<int> m_localBodyList;

  // bodies sharing this entity in the remote processor
  DLIList<int> m_nonLocalBodyList;

  // processors sharing this entity
  // if this entity is surface, the first processor will mesh it
  DLIList<int> m_sharedProcList;

  DLIList<int> m_ghostProcList;

  // back pointer to the owning entity
  ToolDataUser *ownerEntity;

  // unique id
  int m_uniqueId;

  // if it is interface entity shared by processors
  int m_interface;

public:
  
  TDParallel(ToolDataUser *owner = NULL, DLIList<int> *shared_bodies = NULL,
	     DLIList<int> *shared_procs = NULL, DLIList<int> *ghost_procs = NULL,
             int unique_id = -1, int interface = 0);
  
  virtual ~TDParallel();
    //-constructor and destructor

  static int is_parallel(const ToolData* td);
  
  DLIList<int>* get_shared_body_list();
  CubitStatus set_shared_body_list();

  DLIList<int>* get_shared_proc_list();
  void set_shared_proc_list(DLIList<int> list);

  DLIList<int>* get_ghost_proc_list();

  DLIList<int>* get_local_body_list();
  DLIList<int>* get_non_local_body_list();

  ToolDataUser *owner_entity() {return ownerEntity;};

  void owner_entity(ToolDataUser *owner) {ownerEntity = owner;};
  
  //- get/set functions for ownerEntity
  unsigned int get_charge_proc();
  void set_charge_proc(unsigned int proc);
  
  unsigned int get_counter_proc();
  void set_counter_proc(unsigned int proc);

  unsigned int get_unique_id();
  void set_unique_id(unsigned int id);

  int is_interface();

  unsigned int get_n_shared_procs();

  void add_shared_proc(unsigned int proc);

  void add_ghost_proc(unsigned int proc);

  bool is_shared_proc(unsigned int proc);
};

inline int TDParallel::is_interface() {
  return m_interface;
}

inline DLIList<int>* TDParallel::get_local_body_list() {
  return &m_localBodyList;
}

inline DLIList<int>* TDParallel::get_non_local_body_list() {
  return &m_nonLocalBodyList;
}

inline DLIList<int>* TDParallel::get_shared_proc_list() {
  return &m_sharedProcList;
}

inline DLIList<int>* TDParallel::get_ghost_proc_list() {
  return &m_ghostProcList;
}

inline unsigned int TDParallel::get_charge_proc()
{
  return m_sharedProcList[0];
}

inline void TDParallel::set_charge_proc(unsigned int proc)
{
  m_sharedProcList.insert_first(proc);
}

inline unsigned int TDParallel::get_counter_proc()
{
  return m_sharedProcList[1];
}

inline void TDParallel::set_counter_proc(unsigned int proc)
{
  m_sharedProcList.append(proc);
}

inline unsigned int TDParallel::get_unique_id()
{
  return m_uniqueId;
}

inline void TDParallel::set_unique_id(unsigned int id)
{
   m_uniqueId = id;
}

inline unsigned int TDParallel::get_n_shared_procs()
{
  return m_sharedProcList.size();
}

inline void TDParallel::add_shared_proc(unsigned int proc)
{
  m_sharedProcList.append_unique(proc);
}

inline void TDParallel::add_ghost_proc(unsigned int proc)
{
  if (proc != get_charge_proc()) {
    m_ghostProcList.append_unique(proc);
  }
}

inline bool TDParallel::is_shared_proc(unsigned int p)
{
  return m_sharedProcList.is_in_list(p);
}

#endif 














