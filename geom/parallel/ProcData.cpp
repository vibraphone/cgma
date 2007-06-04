//
//
#include "ProcData.hpp"

ProcData *ProcData::instance_ = 0;
int ProcData::isInitialized = 0;

int ProcData::initialize(int &argc, char **&argv) 
{
  if (isInitialized == 1) return 1;

#ifdef USE_MPI  
  int ierror;
  ierror = MPI_Init(&argc, &argv);
  ierror = MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
  ierror = MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
#else
  myRank = 0;
  numProcs = 1;
#endif
  
/* initialize designated rank of master */  
  masterRank = 0;

  isInitialized = 1;

  return ierror;
}

ProcData::~ProcData() 
{
#ifdef USE_MPI
  MPI_Finalize();
#endif

  isInitialized = 0;
}

