#ifndef PROBDATA_HPP
#define PROBDATA_HPP

#ifdef USE_MPI
#include "CGMmpi.h"
#endif

class ProcData 
{
public:

  ~ProcData();
    // destructor; calls MPI_Finalize

  static ProcData *instance();
    // singleton instance
  
  int initialize(int &argc, char **&argv);
    // initialize some communication parameters for this run

  int is_master() const {return myRank == masterRank;};
    // return whether this processor is the master or not
  
  int numProcs;
    // number of processors in this run

  int myRank;
    // rank of this processor
  
  int masterRank;
    // rank of the master process

private:

  static ProcData *instance_;
    // static singleton instance

  static int isInitialized;
    // is this procdata initialized?

  ProcData() {};
    // private constructor, since this is a singleton

};

inline ProcData *ProcData::instance() 
{
   if (instance_ == 0 ) {
     instance_ = new ProcData;
   }
   return instance_;
}

#endif
