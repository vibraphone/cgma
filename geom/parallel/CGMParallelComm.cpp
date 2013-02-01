#include "CGMParallelComm.hpp"
#include "TopologyEntity.hpp"
#include "GeometryQueryEngine.hpp"
#include "RefEntity.hpp"
#include "GeometryQueryTool.hpp"
#include "TDParallel.hpp"

#include <algorithm>

#define INITIAL_BUFF_SIZE 1024
#define RRA(a) if (CUBIT_SUCCESS != result) {	\
    std::string tmp_str;			\
    tmp_str.append("\n"); tmp_str.append(a);	\
    PRINT_ERROR(tmp_str.c_str());		\
    return result;}

std::vector<CGMParallelComm*> CGMParallelComm::instanceList;

CGMParallelComm::CGMParallelComm(MPI_Comm comm, int* id )
  : procConfig(comm)
{
  myBuffer.resize(INITIAL_BUFF_SIZE);
  
  int flag = 1;
  int retval = MPI_Initialized(&flag);
  if (MPI_SUCCESS != retval || !flag) {
    int argc = 0;
    char **argv = NULL;
    
    // mpi not initialized yet - initialize here
    retval = MPI_Init(&argc, &argv);
  }
  
  pcommID = add_pcomm(this);
  if (id)
    *id = pcommID;

  m_pBuffer = NULL;
  m_nBufferSize = 0;
  m_currentPosition = 0;
  set_master(0);
}

CGMParallelComm::CGMParallelComm(std::vector<unsigned char> &tmp_buff, 
				 MPI_Comm comm,
				 int* id) 
  : procConfig(comm)
{
  myBuffer.swap(tmp_buff);
  int flag = 1;
  int retval = MPI_Initialized(&flag);
  if (MPI_SUCCESS != retval || !flag) {
    int argc = 0;
    char **argv = NULL;
    
    // mpi not initialized yet - initialize here
    retval = MPI_Init(&argc, &argv);
  }

  pcommID = add_pcomm(this);
  if (id)
    *id = pcommID;

  m_pBuffer = NULL;
  m_nBufferSize = 0;
  m_currentPosition = 0;
  set_master(0);
}

CGMParallelComm::~CGMParallelComm() 
{
  remove_pcomm(this);
  delete m_pBuffer;
}

int CGMParallelComm::add_pcomm(CGMParallelComm *pc) 
{
  instanceList.push_back(pc);
  return instanceList.size() - 1;
}

void CGMParallelComm::remove_pcomm(CGMParallelComm *pc) 
{
  std::vector<CGMParallelComm*>::iterator pc_it = 
    std::find(instanceList.begin(), instanceList.end(), pc);
  assert(pc_it != instanceList.end());
  instanceList.erase(pc_it);
}

//! get the indexed pcomm object from the interface
CGMParallelComm *CGMParallelComm::get_pcomm(const int index) 
{  
  if (instanceList.size() < (unsigned int) (index + 1)) return NULL;
  else return instanceList[index];
}

CubitStatus CGMParallelComm::get_all_pcomm(std::vector<CGMParallelComm*>& list )
{
  list = instanceList;
  return CUBIT_SUCCESS;
}
  

//! get the indexed pcomm object from the interface
CGMParallelComm *CGMParallelComm::get_pcomm(//MBEntityHandle prtn,
					    const MPI_Comm* comm ) 
{
  CGMParallelComm* result = 0;
  int pcomm_id;
  result= get_pcomm( pcomm_id );

  return result;
}

CubitStatus CGMParallelComm::bcast_buffer(const unsigned int from_proc) 
{
  //- broadcasts the buffer contained in this object
  if (procConfig.proc_rank() == from_proc) {
    printf("Broadcasting buffer size from %d.\n", from_proc);
    MPI_Bcast(&m_nBufferSize, 1, MPI_INT, from_proc, MPI_COMM_WORLD);
    printf("Broadcasting buffer from %d, %d bytes.\n", from_proc,
	   m_nBufferSize);
    MPI_Bcast(m_pBuffer, m_nBufferSize, MPI_BYTE, from_proc, 
              MPI_COMM_WORLD);
  }
  else {
    int this_size;
    printf("Broadcasting buffer size from proc %d.\n",
	   procConfig.proc_rank());
    MPI_Bcast(&this_size, 1, MPI_INT, from_proc, 
              MPI_COMM_WORLD);
    printf("Processor %d: received size of %d.\n", procConfig.proc_rank(), this_size);
    check_size(this_size);
    printf("Broadcasting buffer from proc %d, %d bytes.\n", 
	   procConfig.proc_rank(), this_size);
    MPI_Bcast(m_pBuffer, this_size, MPI_BYTE, from_proc, 
              MPI_COMM_WORLD);
  }
 
  return CUBIT_SUCCESS;
}

CubitStatus CGMParallelComm::broadcast_entities(const unsigned int from_proc,
						DLIList<RefEntity*> &ref_entity_list)
{
#ifndef USE_MPI
  return CUBIT_FAILURE;
#else
  CubitStatus result = CUBIT_SUCCESS;

  if (procConfig.proc_rank() == from_proc) {
    int nBufferSize = 0;
    result = write_buffer(ref_entity_list, m_pBuffer, nBufferSize, false);
    RRA("Failed to write ref entity list to buffer.");

    result = check_size(nBufferSize);
    RRA("Failed to write ref entity list to buffer.");

    result = write_buffer(ref_entity_list, m_pBuffer, nBufferSize, true);
    RRA("Failed to write ref entity list to buffer.");
  }
  
  result = bcast_buffer(from_proc);
  RRA("Failed to broadcast buffer to processors.");
  
  if ((int)procConfig.proc_rank() != from_proc) {
    result = read_buffer(ref_entity_list, m_pBuffer, m_nBufferSize);
    RRA("Failed to read ref entity list from buffer.");
  }
  
  return CUBIT_SUCCESS;
#endif
}

// scatter exact amount of geometry information to each processors
CubitStatus CGMParallelComm::scatter_entities(const unsigned int from_proc,
					      DLIList<RefEntity*> &ref_entity_list)
{
#ifndef USE_MPI
  return CUBIT_FAILURE;
#else
  CubitStatus result = CUBIT_SUCCESS;
  int i, j, mySendCount, nEntity;
  int nProcs = procConfig.proc_size();
  int *sendCounts = new int[nProcs];
  int *displacements = new int[nProcs];
  displacements[0] = 0;

  if (procConfig.proc_rank() == from_proc) {
    // make a balanced entity lists
    int sum = 0;
    DLIList<RefEntity*> **balancedLists = new DLIList<RefEntity*>*[nProcs];
    
    for (i = 0; i < nProcs; i++) {
      balancedLists[i] = new DLIList<RefEntity*>;
    }
    
    nEntity = ref_entity_list.size();
    ref_entity_list.reset();
    for (i = 0; i < nEntity; i++) {
      RefEntity* entity = ref_entity_list.get_and_step();
      TDParallel *td_par = (TDParallel *) entity->get_TD(&TDParallel::is_parallel);
      
      if (td_par == NULL) {
	PRINT_ERROR("Partitioned entities should have TDParallel data.");
	return CUBIT_FAILURE;
      }
      int charge_p = td_par->get_charge_proc();
      if (charge_p != from_proc) { // only to compute processors
        balancedLists[charge_p]->append(entity); // add charge processor
      }
      
      DLIList<int>* ghost_procs = td_par->get_ghost_proc_list();
      int n_ghost = ghost_procs->size();
      ghost_procs->reset();
      for (j = 0; j < n_ghost; j++) { // add ghost processors
        int ghost_p = ghost_procs->get_and_step();
        if (ghost_p != from_proc) balancedLists[ghost_p]->append(entity);
      }
    }
    
    // add buffer size for each processors
    for (i = 0; i < nProcs; i++) {
      result = write_buffer(*(balancedLists[i]), m_pBuffer, sendCounts[i], false);
      RRA("Failed to write ref entity list to buffer.");
      sum += sendCounts[i];
    }
  
    // check the size of the buffer and resize if necessary
    check_size(sum);
    
    // now append the real information
    ref_entity_list.reset();
    for (i = 0; i < nProcs; i++) {
      append_to_buffer(*(balancedLists[i]), sendCounts[i]);
    }

    delete [] balancedLists;
  }

  // broadcast buffer size array
  printf("Broadcasting buffer size array from master.\n");
  MPI_Bcast(sendCounts, nProcs, MPI_INT, from_proc, MPI_COMM_WORLD);
  
  for (i = 1; i < nProcs; i++) {
    displacements[i] = displacements[i-1] + sendCounts[i-1];
  }
  
  mySendCount = sendCounts[procConfig.proc_rank()];

  if (procConfig.proc_rank() != from_proc) check_size(mySendCount);

  printf("Scattering buffer from master.\n");

  // scatter geometry
  MPI_Scatterv(m_pBuffer, sendCounts, displacements, MPI_BYTE, m_pBuffer, 
	       mySendCount, MPI_BYTE, from_proc, MPI_COMM_WORLD);

  if (procConfig.proc_rank() != from_proc) {
    result = read_buffer(ref_entity_list, m_pBuffer, mySendCount);
    RRA("Failed to read ref entity list from buffer.");
  }

  return CUBIT_SUCCESS;
#endif
}

CubitStatus CGMParallelComm::write_buffer(DLIList<RefEntity*> &ref_entity_list,
					  char* pBuffer,
					  int& n_buffer_size,
					  bool b_write_buffer)
{
#ifndef USE_MPI
  return CUBIT_FAILURE;
#else
  if (ref_entity_list.size() == 0) {
    n_buffer_size = 0;
    return CUBIT_SUCCESS;
  }

#ifdef HAVE_OCC
  CubitStatus result = GeometryQueryTool::instance()->export_solid_model(ref_entity_list, pBuffer,
									 n_buffer_size, b_write_buffer);
  RRA("Failed to write ref entities to buffer.");
#endif

  if (b_write_buffer) m_currentPosition += n_buffer_size;
  return CUBIT_SUCCESS;
#endif
}

CubitStatus CGMParallelComm::read_buffer(DLIList<RefEntity*> &ref_entity_list,
					 const char* pBuffer,
					 const int n_buffer_size)
{
#ifndef USE_MPI
  return CUBIT_FAILURE;
#else
  if (n_buffer_size == 0) return CUBIT_SUCCESS;

#ifdef HAVE_OCC
  CubitStatus result = GeometryQueryTool::instance()->import_solid_model(&ref_entity_list, pBuffer,
									 n_buffer_size);
  RRA("Failed to read ref entities from buffer.");
#endif
  
  return CUBIT_SUCCESS;
#endif
}

CubitStatus CGMParallelComm::check_size(int& target_size, const CubitBoolean keep) 
{
  printf("Checking buffer size on proc %d, target size %d.\n", 
                  procConfig.proc_rank(), target_size);

  if (m_nBufferSize < target_size) {
    printf("Increasing buffer size on proc %d.\n", procConfig.proc_rank());
    void *temp_buffer = malloc(target_size);
    if (keep && 0 != m_currentPosition) memcpy(temp_buffer, m_pBuffer, m_currentPosition);
    delete m_pBuffer;
    m_pBuffer = (char *) temp_buffer;
    m_nBufferSize = target_size;
  }

  return CUBIT_SUCCESS;
}

CubitStatus CGMParallelComm::append_to_buffer(DLIList<RefEntity*> &ref_entity_list,
					      int add_size) 
{
  if (m_currentPosition + add_size > m_nBufferSize) return CUBIT_FAILURE;
  CubitStatus result = write_buffer(ref_entity_list, m_pBuffer + m_currentPosition, add_size, true);
  RRA("Failed to append ref entity list to buffer.");
  
  return CUBIT_SUCCESS;
}
