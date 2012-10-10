//-------------------------------------------------------------------------
// Filename		 : SyncTool.hpp
//
// Purpose		 : Defines the SyncTool class
//		   
//
// Special Notes :
//
// Creator       : Derek Quam
//
// Creation Date : 07/03/2007
//
// Owner         : ???
//-------------------------------------------------------------------------
#ifndef SYNC_TOOL_HPP
#define SYNC_TOOL_HPP

#include "CubitDefines.h"

class SyncTool
{
public:
  virtual CubitStatus sync_geom() = 0;
  virtual CubitStatus sync_fem() = 0;
  CubitStatus notify_updated_surfaces();
  CubitStatus notify_updated_curves();
  CubitStatus notify_virtual_curves();
  CubitStatus notify_virtual_surfaces();
};

#endif //SYNC_TOOL_HPP