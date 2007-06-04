//-------------------------------------------------------------------------
// Filename      : LoopACIS.cpp
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
#include "kernel/kerndata/top/loop.hxx"
#include "kernel/kerndata/top/coedge.hxx"
#else
#include "kernapi.hxx"
#include "datamsc.hxx"
#include "loop.hxx"
#include "coedge.hxx"
#endif

#include "attrib_cubit_owner.hpp"
#include "attrib_snl_simple.hpp"
// ********** END ACIS INCLUDES            **********

// ********** BEGIN CUBIT INCLUDES         **********
#include "CastTo.hpp"
#include "RefEntity.hpp"
#include "AcisQueryEngine.hpp"
#include "GeometryQueryTool.hpp"

#include "LoopACIS.hpp"
#include "LoopSM.hpp"
#include "Loop.hpp"
#include "CoEdge.hpp"
#include "CoEdgeACIS.hpp"
#include "DLIList.hpp"

#include "CoEdgeSM.hpp"
#include "Surface.hpp"

// ********** END CUBIT INCLUDES           **********

// ********** BEGIN FORWARD DECLARATIONS   **********
// ********** END FORWARD DECLARATIONS     **********

// ********** BEGIN STATIC DECLARATIONS    **********
// ********** END STATIC DECLARATIONS      **********

// ********** BEGIN PUBLIC FUNCTIONS       **********

//-------------------------------------------------------------------------
// Purpose       : A constructor with a pointer to a ACIS LOOP.
//
// Special Notes :
//
// Creator       : Xuechen Liu
//
// Creation Date : 08/06/96
//-------------------------------------------------------------------------
LoopACIS::LoopACIS(LOOP* LOOP_ptr)
    : AcisBridge(LOOP_ptr)
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
LoopACIS::~LoopACIS()
{}


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
LoopACIS::get_geometry_query_engine() const
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
void LoopACIS::append_simple_attribute_virt(CubitSimpleAttrib* csattrib_ptr)
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
void LoopACIS::remove_simple_attribute_virt(CubitSimpleAttrib* csattrib_ptr)
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
void LoopACIS::remove_all_simple_attribute_virt()
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
CubitStatus LoopACIS::get_simple_attribute(DLIList<CubitSimpleAttrib*>&
                                           cubit_simple_attrib_list)
{
  return AcisBridge::get_simple_attribute(cubit_simple_attrib_list);
}
CubitStatus LoopACIS::get_simple_attribute(const CubitString& name,
                                       DLIList<CubitSimpleAttrib*>& list)
  { return AcisBridge::get_simple_attribute(name, list); }

LOOP* LoopACIS::get_LOOP_ptr() const
{
  return (LOOP*)ENTITY_ptr();
}

void LoopACIS::set_LOOP_ptr(LOOP* LOOP_ptr)
{
  ENTITY_ptr(LOOP_ptr);
}

void LoopACIS::get_parents_virt( DLIList<TopologyBridge*>& parents )
{
  FACE* face = get_LOOP_ptr()->face();
  if( face )
  {
    AcisBridge* acis_bridge = ATTRIB_CUBIT_OWNER::cubit_owner( face );
    TopologyBridge* bridge = dynamic_cast<TopologyBridge*>(acis_bridge);
    if (bridge)
      parents.append( bridge );
  }
}

void LoopACIS::get_children_virt( DLIList<TopologyBridge*>& children )
{
  ENTITY_LIST entities;
  api_get_coedges( get_LOOP_ptr(), entities );
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

