//-------------------------------------------------------------------------
// Filename      : LocalToleranceTool.hpp
//
// Purpose       : Determines smart local tolerances at ref entities
//
// Creator       : William Roshan Quadros
//
// Creation Date :  11/16/2010
//
// Owner         : 
//-------------------------------------------------------------------------

#ifndef LOCAL_TOLERANCE_TOOL_HPP
#define LOCAL_TOLERANCE_TOOL_HPP

#include "CubitDefines.h"
#include "CubitGeomConfigure.h" 
#include "DLIList.hpp"
class BodySM;

class CUBIT_GEOM_EXPORT LocalToleranceTool
{
private:
   //Static pointer to the unique instance of this class.
  static LocalToleranceTool* instance_;

  // Smartly set local toleraceces on all the ref entities recursively 
  bool calculate_local_tolerances_automatically( DLIList<BodySM*> body_sm_list );

  
public:
  LocalToleranceTool();
  ~LocalToleranceTool();
  
  // Return LocalToleranceTool* - Pointer to the singleton LocalToleranceTool object
  static LocalToleranceTool* instance( void );
  static void delete_instance();

  // Main interface for calculating local tolerances
  bool calculate_local_tolerances( DLIList<BodySM*> body_sm_list );

  // Debugging function to print local tolerances at all ref_entities of the input bodies
  bool print_local_tolerances( DLIList<BodySM*> body_sm_list );

};

#endif

