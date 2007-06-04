//-------------------------------------------------------------------------
// Filename      : CoEdgeACIS.cpp
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
// ********** END STANDARD INCLUDES        **********

// ********** BEGIN ACIS INCLUDES          **********
#if CUBIT_ACIS_VERSION < 1100
#include "kernel/kernapi/api/kernapi.hxx"
#include "kernel/kerndata/data/datamsc.hxx"
#include "kernel/kerndata/top/coedge.hxx"
#include "kernel/kerndata/top/loop.hxx"
#include "kernel/kerndata/top/edge.hxx"
#else
#include "kernapi.hxx"
#include "datamsc.hxx"
#include "coedge.hxx"
#include "loop.hxx"
#include "edge.hxx"
#endif
// ********** END ACIS INCLUDES            **********

// ********** BEGIN CUBIT INCLUDES         **********
#include "attrib_cubit_owner.hpp"
#include "CoEdgeACIS.hpp"
#include "AcisQueryEngine.hpp"
#include "CastTo.hpp"
#include "CoEdge.hpp"
#include "LoopSM.hpp"
#include "Curve.hpp"
// ********** END CUBIT INCLUDES           **********

// ********** BEGIN FORWARD DECLARATIONS   **********
// ********** END FORWARD DECLARATIONS     **********

// ********** BEGIN STATIC DECLARATIONS    **********
// ********** END STATIC DECLARATIONS      **********

// ********** BEGIN PUBLIC FUNCTIONS       **********

//-------------------------------------------------------------------------
// Purpose       : A constructor with a pointer to a ACIS COEDGE.
//
// Special Notes :
//
// Creator       : Xuechen Liu
//
// Creation Date : 08/06/96
//-------------------------------------------------------------------------
CoEdgeACIS::CoEdgeACIS(COEDGE* COEDGEPtr)
    : AcisBridge(COEDGEPtr)
{}

//-------------------------------------------------------------------------
// Purpose       : The destructor
//
// Special Notes :
//
// Creator       : Raikanta Sahu
//
// Creation Date : 09/06/96
//-------------------------------------------------------------------------
CoEdgeACIS::~CoEdgeACIS()
{}

COEDGE* CoEdgeACIS::get_COEDGE_ptr() const
{
  return (COEDGE*)ENTITY_ptr();
}

void CoEdgeACIS::set_COEDGE_ptr(COEDGE* COEDGE_ptr)
{
  ENTITY_ptr(COEDGE_ptr);
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
CoEdgeACIS::get_geometry_query_engine() const
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
void CoEdgeACIS::append_simple_attribute_virt(CubitSimpleAttrib* csattrib_ptr)
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
void CoEdgeACIS::remove_simple_attribute_virt(CubitSimpleAttrib* csattrib_ptr)
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
void CoEdgeACIS::remove_all_simple_attribute_virt()
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
CubitStatus CoEdgeACIS::get_simple_attribute(DLIList<CubitSimpleAttrib*>&
                                             cubit_simple_attrib_list)
{
  return AcisBridge::get_simple_attribute(cubit_simple_attrib_list);
}
CubitStatus CoEdgeACIS::get_simple_attribute(const CubitString& name,
                                       DLIList<CubitSimpleAttrib*>& list)
  { return AcisBridge::get_simple_attribute(name, list); }


CubitSense CoEdgeACIS::sense()
{
  COEDGE *co_edge = get_COEDGE_ptr();
  if (co_edge && (co_edge->sense() == FORWARD))
    return CUBIT_FORWARD;
  else if (co_edge && (co_edge->sense() == REVERSED))
    return CUBIT_REVERSED;
  else return CUBIT_UNKNOWN;
}

void CoEdgeACIS::get_parents_virt( DLIList<TopologyBridge*>& parents )
{
  LOOP* loop = get_COEDGE_ptr()->loop();
  if( loop )
  {
    AcisBridge* acis_bridge = ATTRIB_CUBIT_OWNER::cubit_owner( loop );
    TopologyBridge* bridge = dynamic_cast<TopologyBridge*>(acis_bridge);
    if (bridge)
      parents.append( bridge );
  }
}

void CoEdgeACIS::get_children_virt( DLIList<TopologyBridge*>& children )
{
  EDGE* edge = get_COEDGE_ptr()->edge();
  if( edge )
  {
    AcisBridge* acis_bridge = ATTRIB_CUBIT_OWNER::cubit_owner( edge );
    TopologyBridge* bridge = dynamic_cast<TopologyBridge*>(acis_bridge);
    if (bridge)
      children.append( bridge );
  }
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

