#ifndef CGM_PROC_CONFIG_HPP
#define CGM_PROC_CONFIG_HPP

#ifdef USE_MPI
#  include "CGMmpi.h"
#else
typedef int MPI_Comm;
#define MPI_COMM_WORLD 0
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
      
    //! get/set the communicator for this proc config
  const MPI_Comm proc_comm() const {return procComm;}
  void proc_comm(MPI_Comm this_comm) {procComm = this_comm;}

  void set_master(unsigned int proc);
  unsigned int get_master();
  
private:

    //! MPI communicator set for this instance
  MPI_Comm procComm;

    //! rank of this processor
  unsigned procRank;
  
    //! number of processors
  unsigned procSize;

    //! whether the crystal router's been initialized or not
  bool crystalInit;
  
  unsigned int master;
  
};

inline void CGMProcConfig::set_master(unsigned int proc) {
  master = proc;
}

inline unsigned int CGMProcConfig::get_master() {
  return master;
}

#endif
