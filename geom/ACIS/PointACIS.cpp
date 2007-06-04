//-------------------------------------------------------------------------
// Filename      : PointACIS.cpp
//
// Purpose       : 
//
// Special Notes :
//
// Creator       : Xuechen Liu
//
// Creation Date : 08/02/96
//
// Owner         : Malcolm J. Panthaki
//-------------------------------------------------------------------------

// ********** BEGIN STANDARD INCLUDES      **********
#include <assert.h>
// ********** END STANDARD INCLUDES        **********

// ********** BEGIN ACIS INCLUDES          **********
#if CUBIT_ACIS_VERSION < 1100
#include "kernel/kernapi/api/kernapi.hxx"
#include "kernel/kerndata/data/datamsc.hxx"
#include "kernel/kerndata/top/vertex.hxx"
#include "kernel/kerndata/geom/point.hxx"
#else
#include "kernapi.hxx"
#include "datamsc.hxx"
#include "vertex.hxx"
#include "point.hxx"
#endif

#include "attrib_cubit_owner.hpp"
#include "attrib_snl_simple.hpp"
// ********** END ACIS INCLUDES            **********

// ********** BEGIN CUBIT INCLUDES         **********
#include "CubitMessage.hpp"

#include "PointACIS.hpp"
#include "RefVertex.hpp"
#include "Curve.hpp"

#include "AcisQueryEngine.hpp"
#include "CastTo.hpp"
// ********** END CUBIT INCLUDES           **********

// ********** BEGIN STATIC DECLARATIONS    **********
// ********** END STATIC DECLARATIONS      **********

// ********** BEGIN PUBLIC FUNCTIONS       **********

//-------------------------------------------------------------------------
// Purpose       : The constructor with a pointer to the associated VERTEX . 
//
// Special Notes :
//
// Creator       : Xuechen Liu
//
// Creation Date : 08/02/96
//-------------------------------------------------------------------------
PointACIS::PointACIS(VERTEX* VERTEX_ptr)
    : AcisBridge(VERTEX_ptr)
{}

//-------------------------------------------------------------------------
// Purpose       : The destructor. 
//
// Special Notes :
//
// Creator       : Raikanta Sahu
//
// Creation Date : 10/22/96
//-------------------------------------------------------------------------
PointACIS::~PointACIS() 
{}

//-------------------------------------------------------------------------
// Purpose       : This function adds a VERTEX pointer to the list of the
//                 VERTEX pointers associated with this object.  
//
// Special Notes :
//
// Creator       : Xuechen Liu
//
// Creation Date : 08/02/96
//-------------------------------------------------------------------------
void PointACIS::set_VERTEX_ptr(VERTEX* VERTEX_ptr)
{
  ENTITY_ptr(VERTEX_ptr);
}


VERTEX *PointACIS::get_VERTEX_ptr() const
{
  return (VERTEX*)ENTITY_ptr();
}

//-------------------------------------------------------------------------
// Purpose       : The purpose of this function is to append a
//                 attribute to the GE. The name is attached to the 
//                 underlying solid model entity this one points to.
//
//
// Special Notes : 
//
// Creator       : Malcolm J. Panthaki
//
// Creation Date : 11/21/96
//-------------------------------------------------------------------------
void PointACIS::append_simple_attribute_virt(CubitSimpleAttrib* csattrib_ptr)
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
void PointACIS::remove_simple_attribute_virt(CubitSimpleAttrib* csattrib_ptr)
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
void PointACIS::remove_all_simple_attribute_virt()
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
CubitStatus PointACIS::get_simple_attribute(DLIList<CubitSimpleAttrib*>&
                                            cubit_simple_attrib_list)
{
  return AcisBridge::get_simple_attribute(cubit_simple_attrib_list);
}
CubitStatus PointACIS::get_simple_attribute(const CubitString& name,
                                       DLIList<CubitSimpleAttrib*>& list)
  { return AcisBridge::get_simple_attribute(name, list); }

//-------------------------------------------------------------------------
// Purpose       : Returns the coordinates of this Point. 
//
// Special Notes :
//
// Creator       : Malcolm J. Panthaki
//
// Creation Date : 10/12/96
//-------------------------------------------------------------------------
CubitVector PointACIS::coordinates() const
{
    // Get the first ACIS VERTEX associated with this PointACIS object.
  VERTEX* VERTEX_ptr = get_VERTEX_ptr();
  
    // Get the coordinates and return them
  if (VERTEX_ptr != NULL)
  {
      // Get the coordinates of the VERTEX
    APOINT* point = VERTEX_ptr->geometry() ;
    SPAposition pos = point->coords() ;
    
    CubitVector coords( pos.x(), pos.y(), pos.z() ) ;
    return coords;
  }
  else
  {
    PRINT_ERROR("In PointACIS::coordinates\n"
                "       No ACIS VERTEX associated with this PointACIS object.\n"
                "       THIS IS A BUG - PLEASE REPORT IT.\n");
    assert(VERTEX_ptr != NULL);
    return CubitVector();
  }
}

//-------------------------------------------------------------------------
// Purpose       : Get geometry modeling engine: AcisQueryEngine
//
// Special Notes :
//
// Creator       : jihong Ma
//
// Creation Date : 10/22/96
//-------------------------------------------------------------------------
GeometryQueryEngine* PointACIS::get_geometry_query_engine() const
{
  return get_acis_query_engine();   
}                 

//-------------------------------------------------------------------------
// Purpose       : Get the bounding box of the object.
//
// Special Notes :
//
// Creator       : Raikanta Sahu
//
// Creation Date : 10/23/96
//-------------------------------------------------------------------------
CubitBox PointACIS::bounding_box() const 
{
  CubitVector temp_vector = this->coordinates();
  CubitBox temp_box(temp_vector);
  return temp_box;
}

void PointACIS::get_parents_virt( DLIList<TopologyBridge*>& parents )
{
  ENTITY_LIST entities;
  api_get_edges( get_VERTEX_ptr(), entities );
  ATTRIB_CUBIT_OWNER::cubit_owner( entities, parents );
}

void PointACIS::get_children_virt( DLIList<TopologyBridge*>& )
{
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
