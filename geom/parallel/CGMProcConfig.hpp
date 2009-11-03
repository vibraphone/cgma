#ifndef CGM_PROC_CONFIG_HPP
#define CGM_PROC_CONFIG_HPP

//#include "MBTypes.h"
//#include "MBRange.hpp"

//class MBInterface;


#ifdef USE_MPI
#  include "CGMmpi.h"
#else
typedef int MPI_Comm;
#define MPI_COMM_WORLD 0
//typedef void* crystal_data;
#endif

/**\brief Multi-CPU information for parallel CGM */
class CGMProcConfig {
public:

  CGMProcConfig(MPI_Comm proc_comm = MPI_COMM_WORLD);
  
  ~CGMProcConfig();
  
    //! Get the current processor number
  unsigned proc_rank() const 
    { return procRank; }
      
    //! Get the number of processors
  unsigned proc_size() const 
    { return procSize; }
      
    //! get a crystal router for this parallel job
  //crystal_data *crystal_router(bool construct_if_missing = true);

    //! get/set the communicator for this proc config
  const MPI_Comm proc_comm() const {return procComm;}
  void proc_comm(MPI_Comm this_comm) {procComm = this_comm;}
  
private:

    //! MPI communicator set for this instance
  MPI_Comm procComm;

    //! rank of this processor
  unsigned procRank;
  
    //! number of processors
  unsigned procSize;

    //! whether the crystal router's been initialized or not
  bool crystalInit;
  
    //! crystal router for this parallel job
  //crystal_data crystalData;
  
};

#endif
