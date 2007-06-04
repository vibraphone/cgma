//-------------------------------------------------------------------------
// Filename      : ShellACIS.cpp
//
// Purpose       : 
//
// Special Notes :
//
// Creator       : Xuechen Liu
//
// Creation Date : 08/06/96
//
// Owner         : Malcolm J. Panthaki
//-------------------------------------------------------------------------

// ********** BEGIN STANDARD INCLUDES      **********
#include <stddef.h>
// ********** END STANDARD INCLUDES        **********

// ********** BEGIN ACIS INCLUDES          **********
#if CUBIT_ACIS_VERSION < 1100
#include "kernel/kernapi/api/kernapi.hxx"
#include "kernel/kerndata/data/datamsc.hxx"
#include "kernel/kerndata/top/shell.hxx"
#include "kernel/kerndata/top/lump.hxx"
#else
#include "kernapi.hxx"
#include "datamsc.hxx"
#include "shell.hxx"
#include "lump.hxx"
#endif
// ********** END ACIS INCLUDES            **********

// ********** BEGIN CUBIT INCLUDES         **********
#include "attrib_cubit_owner.hpp"
#include "attrib_snl_simple.hpp"

#include "CastTo.hpp"

#include "AcisQueryEngine.hpp"
#include "ShellACIS.hpp"
#include "ShellSM.hpp"

#include "Lump.hpp"
#include "Surface.hpp"

// ********** END CUBIT INCLUDES           **********

// ********** BEGIN STATIC DECLARATIONS    **********
// ********** END STATIC DECLARATIONS      **********

// ********** BEGIN PUBLIC FUNCTIONS       **********

//-------------------------------------------------------------------------
// Purpose       : A constructor with a pointer to a ACIS SHELL.
//
// Special Notes :
//
// Creator       : Xuechen Liu
//
// Creation Date : 08/06/96
//-------------------------------------------------------------------------
ShellACIS::ShellACIS(SHELL* SHELL_ptr)
    : AcisBridge(SHELL_ptr)
{}

//-------------------------------------------------------------------------
// Purpose       : The destructor.
//
// Special Notes :
//
// Creator       : Raikanta Sahu
//
// Creation Date : 09/06/96
//-------------------------------------------------------------------------
ShellACIS::~ShellACIS()
{}

void ShellACIS::set_SHELL_ptr(SHELL* SHELL_ptr)
{
  ENTITY_ptr(SHELL_ptr);
}

//-------------------------------------------------------------------------
// Purpose       : Get geometry modeling engine: AcisQueryEngine
//
// Special Notes :
//
// Creator       : Stephen J. Verzi
//
// Creation Date : 02/26/97
//-------------------------------------------------------------------------
GeometryQueryEngine* 
                 ShellACIS::get_geometry_query_engine() const
{
   return get_acis_query_engine();   
}                 

//-------------------------------------------------------------------------
// Purpose       : The purpose of this function is to append a
//                 attribute to the OSME. The name is attached to the 
//                 underlying solid model entity this one points to.
//
//
// Special Notes : 
//
// Creator       : Malcolm J. Panthaki
//
// Creation Date : 11/21/96
//-------------------------------------------------------------------------
void ShellACIS::append_simple_attribute_virt(CubitSimpleAttrib* csattrib_ptr)
{
  AcisBridge::append_simple_attribute_virt(csattrib_ptr);
}

//-------------------------------------------------------------------------
// Purpose       : The purpose of this function is to remove a simple 
//                 attribute attached to this geometry entity. The name is 
//                 removed from the underlying BODY this points to.
//
// Special Notes : 
//
// Creator       : David R. White
//
// Creation Date : 03/18/97
//-------------------------------------------------------------------------
void ShellACIS::remove_simple_attribute_virt(CubitSimpleAttrib* csattrib_ptr)
{
  AcisBridge::remove_simple_attribute_virt(csattrib_ptr);
}

//-------------------------------------------------------------------------
// Purpose       : The purpose of this function is to remove all simple 
//                 attributes attached to this geometry entity.  Also
//                 removes lingering GTC attributes.
//
//
// Special Notes : 
//
// Creator       : Greg Nielson
//
// Creation Date : 07/10/98
//-------------------------------------------------------------------------
void ShellACIS::remove_all_simple_attribute_virt()
{
  AcisBridge::remove_all_simple_attribute_virt();
}

//-------------------------------------------------------------------------
// Purpose       : The purpose of this function is to get the  
//                 attributes attached to this geometry entity. The name is 
//                 attached to the underlying BODY this points to.
//
// Special Notes : 
//
// Creator       : Malcolm J. Panthaki
//
// Creation Date : 01/23/97
//-------------------------------------------------------------------------
CubitStatus ShellACIS::get_simple_attribute(DLIList<CubitSimpleAttrib*>&
                                           cubit_simple_attrib_list)
{
  return AcisBridge::get_simple_attribute(cubit_simple_attrib_list);
}
CubitStatus ShellACIS::get_simple_attribute(const CubitString& name,
                                       DLIList<CubitSimpleAttrib*>& list)
  { return AcisBridge::get_simple_attribute(name,list); }

/*
void ShellACIS::bodysms(DLIList<BodySM*> &bodies) 
{
  AcisBridge::bodysms(bodies);
}

void ShellACIS::lumps(DLIList<Lump*> &lumps)
{
  AcisBridge::lumps(lumps);
}

void ShellACIS::shellsms(DLIList<ShellSM*> &shellsms)
{
  AcisBridge::shellsms(shellsms);
}

void ShellACIS::surfaces(DLIList<Surface*> &surfaces)
{
  AcisBridge::surfaces(surfaces);
}

void ShellACIS::loopsms(DLIList<LoopSM*> &loopsms)
{
  AcisBridge::loopsms(loopsms);
}

void ShellACIS::curves(DLIList<Curve*> &curves)
{
  AcisBridge::curves(curves);
}

void ShellACIS::coedgesms(DLIList<CoEdgeSM*> &coedgesms)
{
  AcisBridge::coedgesms(coedgesms);
}

void ShellACIS::points(DLIList<Point*> &points)
{
  AcisBridge::points(points);
}
*/


void ShellACIS::get_parents_virt( DLIList<TopologyBridge*>& parents )
{
  LUMP* lump = get_SHELL_ptr()->lump();
  if (lump)
  { 
    AcisBridge* acis_bridge = ATTRIB_CUBIT_OWNER::cubit_owner( lump );
    TopologyBridge* bridge = dynamic_cast<TopologyBridge*>(acis_bridge);
    if (bridge)
      parents.append( bridge );
  }
}

void ShellACIS::get_children_virt( DLIList<TopologyBridge*>& children )
{
  ENTITY_LIST entities;
  api_get_faces( get_SHELL_ptr(), entities );
  ATTRIB_CUBIT_OWNER::cubit_owner( entities, children );
}

// ********** END PUBLIC FUNCTIONS         **********

// ********** BEGIN PROTECTED FUNCTIONS    **********
// ********** END PROTECTED FUNCTIONS      **********

// ********** BEGIN PRIVATE FUNCTIONS      **********
// ********** END PRIVATE FUNCTIONS        **********

// ********** BEGIN HELPER CLASSES         **********
// ********** END HELPER CLASSES           **********

// ********** BEGIN EXTERN FUNCTIONS       **********
// ********** END EXTERN FUNCTIONS         **********

// ********** BEGIN STATIC FUNCTIONS       **********
// ********** END STATIC FUNCTIONS         **********

