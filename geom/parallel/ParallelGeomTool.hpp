#ifndef PARALLELGEOMTOOL_HPP
#define PARALLELGEOMTOOL_HPP

#include <vector>
#include "CubitDefines.h"
#include "CGMParallelComm.hpp"
#include "FileOptions.hpp"

class ParallelGeomTool 
{
public:

  ParallelGeomTool(CGMParallelComm *pc = NULL);

  //static ParallelGeomTool *instance();

  enum ParallelLoadOption {BCAST, BCAST_AND_DELETE, SCATTER};
  
  // load a file
  CubitStatus load_file(const char *file_names, const char *options,
			const char* set_tag_name = 0,
			const int* set_tag_values = 0,
			int num_set_tag_values = 0);
  
  CubitStatus load_file(const char *file_name,
                        int parallel_mode, 
                        std::string &partition_tag_name, 
                        std::vector<int> &partition_tag_vals, 
                        bool distrib,
                        std::vector<int> &pa_vec,
                        const FileOptions &opts,
                        const char* set_tag_name,
                        const int* set_tag_values,
                        const int num_set_tag_values,
                        const int reader_rank,
			const bool surf_partition,
			const bool body_partition);
  
  static const char *CGMparallelOptsNames[];

  static const char *CGMpartitionOptsNames[];
  
  enum CGMParallelOpts {POPT_NONE=0, POPT_READ, POPT_READ_DELETE, POPT_BCAST,
			POPT_BCAST_DELETE, PORT_SCATTER,
			PORT_SCATTER_DELETE, POPT_READ_PARALLEL,
			POPT_FORMAT, POPT_DEFAULT};

private:

  // surf ref entity list to be partitioned
  DLIList<RefEntity*> msurf_entity_list;

  // body ref entity list to be partitioned
  DLIList<RefEntity*> mbody_entity_list;

  CubitStatus delete_nonlocal_entities(std::string &ptag_name,
                                       std::vector<int> &ptag_vals,
				       DLIList<RefEntity*> surf_entity_list,
				       DLIList<RefEntity*> body_entity_list,
				       bool distribute);
  
  //CubitStatus delete_nonlocal_entities(MBEntityHandle file_set);


  //static ParallelGeomTool *instance_;

  // each reader can keep track of its own pcomm
  CGMParallelComm *myPcomm;

  //GeometryQueryTool *gqt;
};

/*
inline ParallelGeomTool *ParallelGeomTool::instance()
{
  if (!instance_) instance_ = new ParallelGeomTool();
  return instance_;
  }*/

#endif
