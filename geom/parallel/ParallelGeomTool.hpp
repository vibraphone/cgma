#ifndef PARALLELGEOMTOOL_HPP
#define PARALLELGEOMTOOL_HPP

#include "TSTTG_CGM.h"

class ParallelGeomTool 
{
public:
  static ParallelGeomTool *instance();

  enum ParallelLoadOption {BCAST, BCAST_AND_DELETE, SCATTER};
  
    // load geometry to parallel processors, using ParallelLoadOption
    // to determine how 
  int load_parallel(TSTTG_Instance geom,
                    const char *name, 
                    int par_load_option);
  
private:
  ParallelGeomTool() {}

  static ParallelGeomTool *instance_;
};

  
inline ParallelGeomTool *ParallelGeomTool::instance()
{
  if (!instance_) instance_ = new ParallelGeomTool();
  return instance_;
}

#endif
