#include "CGMProcConfig.hpp"

//! Constructor
CGMProcConfig::CGMProcConfig(MPI_Comm proc_comm) 
  : procComm(proc_comm),
    crystalInit(false)
{
#ifdef USE_MPI
  int rank, size;
  MPI_Comm_rank(procComm, &rank); 
  procRank = (unsigned int) rank;
  MPI_Comm_size(procComm, &size); 
  procSize = (unsigned int) size;
#else
  procRank = 0;
  procSize = 1;
#endif
}

CGMProcConfig::~CGMProcConfig() 
{
}

