#ifndef CGM_READ_PARALLEL_HPP
#define CGM_READ_PARALLEL_HPP

#include <vector>
#include "CubitDefines.h"
#include "GeometryQueryTool.hpp"
#include "CGMParallelComm.hpp"
#include "CGMFileOptions.hpp"

enum BALANCE_METHOD {
  ROUND_ROBIN = 0,
  PARTITION_STATIC,
  PARTITION_DYNAMIC
};

class CGMReadParallel 
{
public:

  CGMReadParallel(GeometryQueryTool* gqt = NULL, CGMParallelComm *pc = NULL);

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
                        std::vector<int> &pa_vec,
                        const CGMFileOptions &opts,
                        const char* set_tag_name,
                        const int* set_tag_values,
                        const int num_set_tag_values,
                        const int reader_rank,
			const bool surf_partition
			);
  
  static const char *CGMparallelOptsNames[];

  static const char *CGMpartitionOptsNames[];
  
  enum CGMParallelOpts {POPT_NONE=0, POPT_READ, POPT_READ_DELETE, POPT_BCAST,
			POPT_BCAST_DELETE, PORT_SCATTER,
			PORT_SCATTER_DELETE, POPT_READ_PARALLEL,
			POPT_FORMAT, POPT_DEFAULT};

  void set_reader(unsigned int reader);

private:

  GeometryQueryTool* m_gqt;

  CGMParallelComm *m_pcomm;

  bool m_scatter, m_reader;

  unsigned int m_rank, m_proc_size;

  BALANCE_METHOD m_bal_method;

  // surf ref entity list to be partitioned
  DLIList<RefEntity*> m_surf_entity_list;

  CubitStatus read_entities(const char* file_name);

  // balance and save the information as attribute
  CubitStatus balance_round_robin();

  CubitStatus delete_nonlocal_entities(int reader,
				       std::string &ptag_name,
                                       std::vector<int> &ptag_vals);

  CubitStatus check_partition_info();
};
#endif
