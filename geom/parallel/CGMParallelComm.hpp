/**
 * \class CGMParallelComm
 * \brief Parallel communications in CGM
 * \author Hong-Jun Kim, copied from MBParallelComm.hpp
 *
 *  This class implements methods to communicate geometry between processors
 *
 */

#ifndef CGM_PARALLEL_COMM_HPP
#define CGM_PARALLEL_COMM_HPP

//#include "CGMForward.hpp"
#include "GeometryQueryTool.hpp"
//#include "CGMRange.hpp"
#include "CGMProcConfig.hpp"
#include "CATag.hpp"
#include <map>
#include <set>
#include <vector>
#include "math.h"
#include "CGMmpi.h"


extern "C" {
  struct tuple_list;
}

//class TagServer;
//class SequenceManager;
//template <typename KeyType, typename ValType, ValType NullVal> class RangeMap;
//typedef RangeMap<CGMEntityHandle, CGMEntityHandle, 0> HandleMap;

#define MAX_SHARING_PROCS 64

class CGMParallelComm 
{
public:

    //! constructor
  CGMParallelComm(CGMTagManager *impl,
		  MPI_Comm comm = MPI_COMM_WORLD,
		  int* pcomm_id_out = 0);
  
  //! constructor taking buffer, for testing
  CGMParallelComm(CGMTagManager *impl,
		  std::vector<unsigned char> &tmp_buff,
		  MPI_Comm comm = MPI_COMM_WORLD,
		  int* pcomm_id_out = 0);
  
    //! Get ID used to reference this PCOMM instance
  int get_id() const { return pcommID; }

    //! get the indexed pcomm object from the interface
  static CGMParallelComm *get_pcomm(CGMTagManager *impl,
				    const int index);
  
    //! Get CGMParallelComm instance associated with partition handle
    //! Will create CGMParallelComm instance if a) one does not already
    //! exist and b) a valid value for MPI_Comm is passed.
  static CGMParallelComm *get_pcomm(CGMTagManager *impl,
                                    //CGMEntityHandle partitioning,
                                    const MPI_Comm* comm = 0 );

  static CubitStatus get_all_pcomm(CGMTagManager *impl,
				   std::vector<CGMParallelComm*>& list );
  
  //! destructor
  ~CGMParallelComm();
  
  //static unsigned char PROC_SHARED, PROC_OWNER;
  /*
    //! assign a global id space, for largest-dimension or all entities (and
    //! in either case for vertices too)
  CGMErrorCode assign_global_ids(CGMEntityHandle this_set,
                                const int dimension,
                                const int start_id = 1,
                                const bool largest_dim_only = true,
                                const bool parallel = true);

    //! check for global ids; based only on tag handle being there or not;
    //! if it's not there, create them for the specified dimensions
  CGMErrorCode check_global_ids(CGMEntityHandle this_set,
                               const int dimension, 
                               const int start_id = 1,
                               const bool largest_dim_only = true,
                               const bool parallel = true);
  */
	
  //! return partition ref entity list
  DLIList<RefEntity*> &partition_surf_list() {return partitioningSurfList;}
  const DLIList<RefEntity*> &partition_surf_list() const {return partitioningSurfList;}
  DLIList<RefEntity*> &partition_body_list() {return partitioningBodyList;}
  const DLIList<RefEntity*> &partition_body_list() const {return partitioningBodyList;}

  //! Get proc config for this communication object
  const CGMProcConfig &proc_config() const {return procConfig;}
  
  //! Get proc config for this communication object
  CGMProcConfig &proc_config() {return procConfig;}

  CubitStatus broadcast_entities(const unsigned int from_proc,
				 DLIList<RefEntity*> &ref_entity_list);

  CubitStatus scatter_entities(const unsigned int from_proc,
			       DLIList<RefEntity*> &ref_entity_list);
  /*
  CubitStatus write_buffer(DLIList<RefEntity*> &ref_entity_list,
			  //const unsigned char* p_buffer);
			  //char* p_buffer);
			  std::ofstream& os);

  CubitStatus write_buffer(DLIList<RefEntity*> &ref_entity_list,
			  std::ostringstream& os);
  */
  CubitStatus write_buffer(DLIList<RefEntity*> &ref_entity_list,
			  char* pBuffer,
			  int& n_buffer_size,
			  bool b_export_buffer);

  CubitStatus read_buffer(DLIList<RefEntity*> &ref_entity_list,
			    const char* pBuffer,
			    const int n_buffer_size);
    
  CubitStatus bcast_buffer(const unsigned int from_proc);

  CubitStatus append_to_buffer(DLIList<RefEntity*> &ref_entity_list,
			       int add_size);
  
private:  

    //! add a pc to the iface instance tag PARALLEL_COMM
  int add_pcomm(CGMParallelComm *pc);
  
    //! remove a pc from the iface instance tag PARALLEL_COMM
  void remove_pcomm(CGMParallelComm *pc);

  CubitStatus check_size(int& target_size, const CubitBoolean keep = CUBIT_FALSE);
  
    //! CGM query tool interface associated with this writer
  GeometryQueryTool *gqt;

    //! CGM tag manager interface associated with this writer
  CGMTagManager *cgmImpl;

    //! Proc config object, keeps info on parallel stuff
  CGMProcConfig procConfig;

    //! Tag server, so we can get more info about tags
  //TagServer *tagServer;

    //! Sequence manager, to get more efficient access to entities
  //SequenceManager *sequenceManager;
  
  
    //! data buffer used to communicate
  std::vector<unsigned char> myBuffer;

  char* m_pBuffer;
  
  int m_nBufferSize;
  
  int m_currentPosition;

  //std::ofstream mOfstream;

    //! more data buffers, proc-specific
  //std::vector<unsigned char> ownerRBuffs[MAX_SHARING_PROCS],
  //ownerSBuffs[MAX_SHARING_PROCS], ghostRBuffs[MAX_SHARING_PROCS],
  // ghostSBuffs[MAX_SHARING_PROCS];

    //! request objects, may be used if store_remote_handles is used
  //MPI_Request sendReqs[2*MAX_SHARING_PROCS];

    //! processor rank for each buffer index
  //std::vector<int> buffProcs;

    //! the partition, interface sets for this comm'n instance
  //CGMRange partitionSets, interfaceSets;
  
    //! local entities ghosted to other procs
  //std::map<unsigned int, CGMRange> ghostedEnts;
  
    //! tags used to save sharing procs and handles
  //CGMTag sharedpTag, sharedpsTag, sharedhTag, sharedhsTag, pstatusTag, 
  //ifaceSetsTag, partitionTag;
    
  //int globalPartCount; //!< Cache of global part count
  
  //CGMEntityHandle partitioningSet; //!< entity set containing all parts
  DLIList<RefEntity*> partitioningSurfList; // ref entity list containing all parts
  DLIList<RefEntity*> partitioningBodyList; // ref entity list containing all parts
  
  int pcommID;
};

#endif
