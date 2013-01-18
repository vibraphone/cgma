//
// class CGMMemFile
//
// This class mimics the Acis MemFile construct, but for general CGM
// entities.  This class supports the saving and restoring of entities
// to/from memory buffers, and the passing of these memory buffers to
// processors using mpi.
//

#ifndef CGMMEMFILE
#define CGMMEMFILE
//#include <vector>
#include "CubitDefines.h"
#include "DLIList.hpp"

class RefEntity;

class CGMMemFile 
{
public:
  CGMMemFile(void* buffer = NULL,
             size_t bufSize = 0);
    //- constructor, with optional buffer passed in; if buffer is passed in, it
    //- will be used for the writing operation; if NULL, the next writing operation
    //- will create a buffer of the right size.  NOTE: buffer is deleted in the
    //- destructor, even if it was passed in

  virtual ~CGMMemFile();
    //- destructor; deletes m_pBuffer
  
  CubitStatus bcast_entity_list(DLIList<RefEntity*> &ref_entity_list);
    //- broadcasts the entity list; root takes list as input, 
    //- others write as output; uses read_refentity_list and write_refentity_list
    //- functions
  
  CubitStatus bcast_and_delete_entity_list(DLIList<RefEntity*> &ref_entity_list,
					   DLIList<RefEntity*> &ref_entity_list_master);
  //- broadcasts the entity list; delete the rest that this processor doesn't need

  //CubitStatus send_body_to_procs(DLIList<RefEntity*> &ref_entity_list);

  //CubitStatus send_body_to_procs_balanced(DLIList<RefEntity*> &ref_entity_list);

  //CubitStatus vol_balace_send_body_processors(DLIList<RefEntity*> &ref_entity_list);
  
  CubitStatus scatter_entity_list(DLIList<RefEntity*> &ref_entity_list);
  //- scatter the exact amount of entity that each processors need;
  
  CubitStatus scatter_balanced_entity_list(DLIList<RefEntity*> &ref_entity_list,
					   DLIList<RefEntity*> &ref_entity_list_master,
					   DLIList<RefEntity*> **balanced_lists);
  
  CubitStatus scatter_test(DLIList<RefEntity*> &ref_entity_list);
  
  //static int *get_body_to_proc();
  static DLIList <int> *get_body_to_proc();
  //static std::vector<int> *get_body_to_proc();

  //static std::vector<int> get_num_body_to_proc();


protected:
  //int balanceSurf;
  //- load balance surfaces

  //double *totalLoads;

  //DLIList<RefEntity*> **balancedLists;

  unsigned char *m_pBuffer;
    //- the memory buffer

  size_t m_sizeBuffer;
    //- the size of the buffer in bytes
  
  unsigned long m_currentPosition;
    //- the current position in the buffer

  static DLIList <int> bodyToProc;
    //static std::vector<int> bodyToProc;
    //- body to processor assignment vector

  //static int *bodyToProc;

  //static int nBodyToProc;

  //static DLIList <int> *bodyToProcPtr;
    //- body to processor assignment vector pointer

  virtual CubitStatus read_refentity_list(DLIList<RefEntity*> &ref_entity_list) = 0;
    //- read a RefEntity list from the buffer into ref_entity_list; relies on engine-
    //- specific implementation

  virtual CubitStatus write_refentity_list(DLIList<RefEntity*> &ref_entity_list) = 0;
    //- write a RefEntity list to the buffer from ref_entity_list; relies on engine-
    //- specific implementation
  
  virtual CubitStatus append_refentity_list(DLIList<RefEntity*> &ref_entity_list,
                                            int &buffer_size) = 0;

  virtual CubitStatus get_refentity_list(DLIList<RefEntity*> &ref_entity_list) = 0;

  virtual int get_refentity_list_size(DLIList<RefEntity*> ref_entity_list) = 0;
  //- get the size of the requested entity list in the buffer

  CubitStatus bcast_buffer();
    //- broadcasts the buffer contained in this object; uses ProcData singleton class

  CubitStatus check_size(int &target_size, CubitBoolean keep = CUBIT_FALSE);
    //- check the size of the buffer and, if necessary, expand it
};

#endif
