//-------------------------------------------------------------------------
// Filename      : OCCCoEdge.cpp
//
// Purpose       : 
//
// Special Notes :
//
// Creator       : Steven J. Owen
//
// Creation Date : 07/18/00
//
// Owner         : Steven J. Owen
//-------------------------------------------------------------------------

// ********** BEGIN STANDARD INCLUDES      **********
// ********** END STANDARD INCLUDES        **********
#include "config.h"
// ********** BEGIN CUBIT INCLUDES         **********
#include "CubitDefines.h"
#include "CastTo.hpp"
#include "OCCCoEdge.hpp"
#include "OCCLoop.hpp"
#include "OCCQueryEngine.hpp"
#include "CubitUtil.hpp"

#include "OCCBody.hpp"
#include "OCCLump.hpp"
#include "OCCShell.hpp"
#include "OCCSurface.hpp"
#include "OCCCurve.hpp"
// ********** END CUBIT INCLUDES           **********

// ********** BEGIN FORWARD DECLARATIONS   **********
// ********** END FORWARD DECLARATIONS     **********

// ********** BEGIN STATIC DECLARATIONS    **********
// ********** END STATIC DECLARATIONS      **********

// ********** BEGIN PUBLIC FUNCTIONS       **********


OCCCoEdge::OCCCoEdge( TopoDS_Edge* edge, Curve *curv_ptr, 
	  	      LoopSM *loop_ptr, CubitSense sense )
	 	    : myTopoDSEdge(edge), myCurve(curv_ptr), 
 	  	      myLoop(loop_ptr),edgeSense(sense)
{
}

//-------------------------------------------------------------------------
// Purpose       : The destructor
//
// Special Notes :
//
// Creator       : Steve Owen
//
// Creation Date : 07/18/00
//-------------------------------------------------------------------------
OCCCoEdge::~OCCCoEdge()
{
  myTopoDSEdge = (TopoDS_Edge*) NULL;
  myCurve = (Curve *) NULL;
  myLoop = (LoopSM *)NULL;
}

//-------------------------------------------------------------------------
// Purpose       : Get geometry modeling engine: OCCQueryEngine
//
// Special Notes :
//
// Creator       : Steve Owen
//
// Creation Date : 07/18/00
//-------------------------------------------------------------------------
GeometryQueryEngine* OCCCoEdge::get_geometry_query_engine() const
{
  return OCCQueryEngine::instance();
}  

//-------------------------------------------------------------------------
// Purpose       : The purpose of this function is to append a
//                 attribute to the OSME. The name is attached to the 
//                 underlying solid model entity this one points to.
//
//
// Special Notes : 
//
// Creator       : Steve Owen
//
// Creation Date : 07/18/00
//-------------------------------------------------------------------------
void OCCCoEdge::append_simple_attribute_virt(CubitSimpleAttrib* /*csattrib_ptr*/)
{
  //PRINT_ERROR("OCCCoEdge::append_simple_attribute_virt not implemented\n");
  return;
}


//-------------------------------------------------------------------------
// Purpose       : The purpose of this function is to remove a simple 
//                 attribute attached to this geometry entity. The name is 
//                 removed from the underlying BODY this points to.
//
// Special Notes : 
//
// Creator       : Steve Owen
//
// Creation Date : 07/18/00
//-------------------------------------------------------------------------
void OCCCoEdge::remove_simple_attribute_virt(CubitSimpleAttrib* /*csattrib_ptr*/)
{
  //PRINT_ERROR("OCCCoEdge::remove_simple_attribute_virt not implemented\n");
  return;
}

//-------------------------------------------------------------------------
// Purpose       : The purpose of this function is to remove all simple 
//                 attributes attached to this geometry entity.  Also
//                 removes lingering GTC attributes.
//
//
// Special Notes : 
//
// Creator       : Steve Owen
//
// Creation Date : 07/18/00
//-------------------------------------------------------------------------
void OCCCoEdge::remove_all_simple_attribute_virt()
{
  //PRINT_ERROR("OCCCoEdge::remove_all_simple_attribute_virt not implemented\n");
  return;
}

//-------------------------------------------------------------------------
// Purpose       : The purpose of this function is to get the  
//                 attributes attached to this geometry entity. The name is 
//                 attached to the underlying BODY this points to.
//
// Special Notes : 
//
// Creator       : Steve Owen
//
// Creation Date : 07/18/00
//-------------------------------------------------------------------------
CubitStatus OCCCoEdge::get_simple_attribute(DLIList<CubitSimpleAttrib*>&
                                              /*cubit_simple_attrib_list*/)
{
  //PRINT_ERROR("OCCCoEdge::get_simple_attribute not implemented\n");
  return CUBIT_FAILURE;
}
CubitStatus OCCCoEdge::get_simple_attribute(const CubitString&,
                                              DLIList<CubitSimpleAttrib*>&)
  { return CUBIT_FAILURE; }


void OCCCoEdge::reverse_sense()
{
  edgeSense = CubitUtil::opposite_sense( edgeSense );
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

