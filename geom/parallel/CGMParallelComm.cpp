#include "CGMParallelComm.hpp"
//#include "CubitAttrib.hpp"
//#include "OCCQueryEngine.hpp"
#include "TopologyEntity.hpp"
#include "GeometryQueryEngine.hpp"
#include "RefEntity.hpp"

//#include <iostream>
#include <sstream>
#include <algoritm>

#ifdef USE_MPI
#include "CGMmpi.h"
#endif

#define INITIAL_BUFF_SIZE 1024
#define RRA(a) if (CUBIT_SUCCESS != result) {	\
    std::string tmp_str;			\
    tmp_str.append("\n"); tmp_str.append(a);	\
    PRINT_ERROR(tmp_str.c_str());		\
    return result;}

CGMParallelComm::CGMParallelComm(CGMTagManager *impl,
				 MPI_Comm comm, int* id )
  : cgmImpl(impl), procConfig(comm)
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
}

CGMParallelComm::CGMParallelComm(CGMTagManager *impl,
				 std::vector<unsigned char> &tmp_buff, 
				 MPI_Comm comm,
				 int* id) 
  : cgmImpl(impl), procConfig(comm)
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
}

CGMParallelComm::~CGMParallelComm() 
{
  remove_pcomm(this);
  delete m_pBuffer;
}

int CGMParallelComm::add_pcomm(CGMParallelComm *pc) 
{
  // add this pcomm to instance tag
  std::vector<CGMParallelComm *> pc_array;
  pc_array = cgmImpl->get_pc_array();
  pc_array.push_back(pc);

  return pc_array.size() - 1;
}

void CGMParallelComm::remove_pcomm(CGMParallelComm *pc) 
{
  // remove this pcomm from instance tag
  std::vector<CGMParallelComm *> pc_array;
  pc_array = cgmImpl->get_pc_array();

  std::vector<CGMParallelComm*>::iterator pc_it = 
    std::find(pc_array.begin(), pc_array.end(), pc);
  assert(pc_it != pc_array.end());
  pc_array.erase(pc_it);
}

//! get the indexed pcomm object from the interface
CGMParallelComm *CGMParallelComm::get_pcomm(CGMTagManager *impl,
					    const int index) 
{  
  std::vector<CGMParallelComm *> pc_array;
  pc_array = impl->get_pc_array();

  if (pc_array.size() < (unsigned int) (index + 1)) return NULL;
  else return pc_array[index];
}

CubitStatus CGMParallelComm::get_all_pcomm(CGMTagManager *impl,
					   std::vector<CGMParallelComm*>& list )
{
  list = impl->get_pc_array();
  return CUBIT_SUCCESS;
}
  

//! get the indexed pcomm object from the interface
CGMParallelComm *CGMParallelComm::get_pcomm(CGMTagManager *impl,
					    //MBEntityHandle prtn,
					    const MPI_Comm* comm ) 
{
  //MBErrorCode rval;
  CGMParallelComm* result = 0;
  
  //MBTag prtn_tag;
  //rval = impl->tag_create( PARTITIONING_PCOMM_TAG_NAME, 
  //                       sizeof(int),
  //                       MB_TAG_SPARSE,
  //                       MB_TYPE_INTEGER,
  //                       prtn_tag,
  //                       0, true );
  //if (MB_SUCCESS != rval)
  //return 0;
  
  int pcomm_id;
  //rval = impl->tag_get_data( prtn_tag, &prtn, 1, &pcomm_id );
  //if (MB_SUCCESS == rval) {
    result= get_pcomm(impl, 
		      pcomm_id );
    //}
    /*
  else if (MB_TAG_NOT_FOUND == rval && comm) {
    result = new MBParallelComm( impl, *comm, &pcomm_id );
    if (!result)
      return 0;
    result->set_partitioning( prtn );
    
    rval = impl->tag_set_data( prtn_tag, &prtn, 1, &pcomm_id );
    if (MB_SUCCESS != rval) {
      delete result;
      result = 0;
    }
    }*/
  
  return result;
}

CubitStatus CGMParallelComm::bcast_buffer(const unsigned int from_proc) 
{
  //- broadcasts the buffer contained in this object

  if ((int)procConfig.proc_rank() == from_proc) {
    PRINT_DEBUG_100("Broadcasting buffer size from %d.\n", from_proc);
    MPI_Bcast(&m_nBufferSize, 1, MPI_INT, from_proc, MPI_COMM_WORLD);
    PRINT_DEBUG_100("Broadcasting buffer from %d, %d bytes.\n", from_proc,
		    m_nBufferSize);
    MPI_Bcast(m_pBuffer, m_nBufferSize, MPI_BYTE, from_proc, 
              MPI_COMM_WORLD);
  }
  else {
    int this_size;
    PRINT_DEBUG_100("Broadcasting buffer size from proc %d.\n",
		    procConfig.proc_rank());
    MPI_Bcast(&this_size, 1, MPI_INT, from_proc, 
              MPI_COMM_WORLD);
    PRINT_DEBUG_100("Processor %d: received size of %d.\n", procConfig.proc_rank(), this_size);
    check_size(this_size);
    PRINT_DEBUG_100("Broadcasting buffer from proc %d, %d bytes.\n", 
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

  if ((int)procConfig.proc_rank() == from_proc) {
    int nBufferSize = 0;
    result = write_buffer(ref_entity_list, m_pBuffer, nBufferSize, false);
    RRA("Failed to write ref entity list to buffer.");

    result = check_size(nBufferSize);
    RRA("Failed to write ref entity list to buffer.");

    result = write_buffer(ref_entity_list, m_pBuffer, nBufferSize, true);
    //m_currentPosition = m_nBufferSize;
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
  CubitStatus status;
  double tScatter;
  int mySendCount;
  int nEntity;
  DLIList<RefEntity*> temp_list, temp_ref_list;

  int nProcs = procConfig.proc_size();
  int *sendCounts = new int[nProcs];
  int *displacements = new int[nProcs];
  int nBarEntity;
  int restEntity;
  int nEndEntity;
  displacements[0] = 0;

  if (procConfig.proc_rank() == from_proc) {
    int curPosition = 0;
    nEntity = ref_entity_list.size();
    nBarEntity = nEntity/nProcs;
    restEntity = nEntity%nProcs;
    nEndEntity = nBarEntity + restEntity;
    
    ref_entity_list.reset();
    int sum = 0;

    // make temporary lists to contain geometry information for each processors
    for (int i = 0; i < nProcs; i++) {
      if (i == from_proc) {
	ref_entity_list.step(nBarEntity);
	if (i < restEntity) ref_entity_list.step();
	sendCounts[i] = 0;
      }
      else {
	for ( int j = 0; j < nBarEntity; j++) {
	  RefEntity* body_ptr = ref_entity_list.get_and_step();
	  temp_list.append(body_ptr);
	}
	
	if (i < restEntity) {
	  RefEntity* body_ptr = ref_entity_list.get_and_step();
	  temp_list.append(body_ptr);
	}
	
	//sendCounts[i] = get_ref_list_buffer_size(temp_list);
	result = write_buffer(temp_list, m_pBuffer, sendCounts[i], false);
	RRA("Failed to write ref entity list to buffer.");
	
	sum += sendCounts[i];
	temp_list.clean_out();
      }
    }

    // check the size of the buffer and resize if necessary
    check_size(sum);
    
    // now append the real information
    ref_entity_list.reset();
    for (int i = 0; i < nProcs; i++) {
      if (i == from_proc) {
	ref_entity_list.step(nBarEntity);
	if (i < restEntity) ref_entity_list.step();
      }
      else {
	for ( int j = 0; j < nBarEntity; j++) {
	  temp_list.append(ref_entity_list.get_and_step());
	}
	
	if (i < restEntity) {
	  temp_list.append(ref_entity_list.get_and_step());
	}
	
	append_to_buffer(temp_list, sendCounts[i]);
	temp_list.clean_out();
      }
    }
  }

  // broadcast buffer size array
  PRINT_DEBUG_100("Broadcasting buffer size array from master.\n");
  MPI_Bcast(sendCounts, nProcs, MPI_INT, from_proc, MPI_COMM_WORLD);
  
  for (int i = 1; i < nProcs; i++) {
    displacements[i] = displacements[i-1] + sendCounts[i-1];
  }
  
  mySendCount = sendCounts[procConfig.proc_rank()];

  if (procConfig.proc_rank() != from_proc) check_size(mySendCount);

  PRINT_DEBUG_100("Scattering buffer from master.\n");

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
/*
CubitStatus CGMParallelComm::write_buffer(DLIList<RefEntity*> &ref_entity_list,
					 std::ofstream& os)
{
#ifndef USE_MPI
  return CUBIT_FAILURE;
#else
  //ofstream os;
  //CubitStatus result = GeometryQueryTool::instance()->export_solid_model(ref_entity_list, p_buffer);
  CubitStatus result = GeometryQueryTool::instance()->export_solid_model(ref_entity_list, os);
  RRA("Failed to compute buffer size in broadcast_entities.");
  
  return CUBIT_SUCCESS;
#endif
}

CubitStatus CGMParallelComm::write_buffer(DLIList<RefEntity*> &ref_entity_list,
					 std::ostringstream& os)
{
#ifndef USE_MPI
  return CUBIT_FAILURE;
#else
  CubitStatus result = GeometryQueryTool::instance()->export_solid_model(ref_entity_list, os);
  RRA("Failed to compute buffer size in broadcast_entities.");
  
  return CUBIT_SUCCESS;
#endif
}
*/
CubitStatus CGMParallelComm::write_buffer(DLIList<RefEntity*> &ref_entity_list,
					  char* pBuffer,
					  int& n_buffer_size,
					  bool b_write_buffer)
{
#ifndef USE_MPI
  return CUBIT_FAILURE;
#else
  CubitStatus result = GeometryQueryTool::instance()->export_solid_model(ref_entity_list, pBuffer,
									 n_buffer_size, b_write_buffer);
  RRA("Failed to write ref entities to buffer.");
  
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
  CubitStatus result = GeometryQueryTool::instance()->import_solid_model(&ref_entity_list, pBuffer,
									 n_buffer_size);
  RRA("Failed to read ref entities from buffer.");
  
  return CUBIT_SUCCESS;
#endif
}

CubitStatus CGMParallelComm::check_size(int& target_size, const CubitBoolean keep) 
{
  PRINT_DEBUG_100("Checking buffer size on proc %d, target size %d.\n", 
                  procConfig.proc_rank(), target_size);

  if (m_nBufferSize < target_size) {
    PRINT_DEBUG_100("Increasing buffer size on proc %d.\n", procConfig.proc_rank());
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
  
  //m_currentPosition += add_size;

  return CUBIT_SUCCESS;
}
