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

//-------------------------------------------------------------------------
// Purpose       : A constructor with a pointer to a Facet CoEdge.
//
// Special Notes :
//
// Creator       : Steve Owen
//
// Creation Date : 07/18/00
//-------------------------------------------------------------------------
OCCCoEdge::OCCCoEdge( TopoDS_Edge *theEdge, Curve *curv_ptr )
{
  myTopoDSEdge = theEdge;
  myCurve = curv_ptr;
}

OCCCoEdge::OCCCoEdge( Curve *curv_ptr, LoopSM *loop_ptr, CubitSense sense )
{
  assert(0);
  myCurve = curv_ptr;
  myLoop = loop_ptr;
  edgeSense = sense; 
}

//-------------------------------------------------------------------------
// Purpose       : A constructor with a pointer to a Facet CoEdge.
//
// Special Notes : Constructor used for save/restore
//
// Creator       : Corey Ernst 
//
// Creation Date : 02/03/03
//-------------------------------------------------------------------------
OCCCoEdge::OCCCoEdge( Curve *curv_ptr, CubitSense sense )
{
  assert(0);
  myCurve = curv_ptr;
  myLoop = NULL; 
  edgeSense = sense; 
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


CubitSense OCCCoEdge::sense()
{
  TopAbs_Orientation d = myTopoDSEdge->Orientation();
  if (d == TopAbs_FORWARD) return CUBIT_FORWARD;
  else if (d == TopAbs_REVERSED) return CUBIT_REVERSED;
  else {
	  printf("Check Orientation");
	  return CUBIT_UNKNOWN;
  }
}
/*
void OCCCoEdge::bodysms(DLIList<BodySM*> &bodies)
{
  myLoop->bodysms(bodies);
}

void OCCCoEdge::lumps(DLIList<Lump*> &lumps)
{
  myLoop->lumps(lumps);
}

void OCCCoEdge::shellsms(DLIList<ShellSM*> &shellsms)
{
  myLoop->shellsms(shellsms);
}

void OCCCoEdge::surfaces(DLIList<Surface*> &surfaces)
{
  myLoop->surfaces( surfaces );
}

void OCCCoEdge::loopsms(DLIList<LoopSM*> &loopsms)
{
  loopsms.append_unique( myLoop );
}

void OCCCoEdge::coedgesms(DLIList<CoEdgeSM*> &coedgesms)
{
  coedgesms.append_unique( this );
}

void OCCCoEdge::curves(DLIList<Curve*> &curves)
{
  curves.append_unique( myCurve );
}

void OCCCoEdge::points(DLIList<Point*> &points)
{
  myCurve->points( points );
}
*/

void OCCCoEdge::get_parents_virt( DLIList<TopologyBridge*>& parents )
  { parents.append( myLoop ); }

void OCCCoEdge::get_children_virt( DLIList<TopologyBridge*>& children )
  { children.append( myCurve ); }

void OCCCoEdge::reverse_sense()
{
  edgeSense = CubitUtil::opposite_sense( edgeSense );
}


void OCCCoEdge::get_lumps( DLIList<OCCLump*>& result_list )
{
  DLIList<OCCShell*> shell_list;
  get_shells( shell_list );
  shell_list.reset();
  for ( int i = shell_list.size(); i--; )
  {
    OCCShell* shell = shell_list.get_and_step();
    shell->get_lumps( result_list );
    OCCLump* lump = dynamic_cast<OCCLump*>(shell->get_lump());
    if (lump)
      result_list.append_unique(lump);
  }
}

void OCCCoEdge::get_shells( DLIList<OCCShell*>& result_list )
{
  if ( OCCLoop* loop = dynamic_cast<OCCLoop*>(myLoop) )
    if ( OCCSurface* surf = dynamic_cast<OCCSurface*>(loop->get_surface()) )
      surf->get_shells( result_list );
}

void OCCCoEdge::get_curves( DLIList<OCCCurve*>& result_list )
{
  if (OCCCurve* curve = dynamic_cast<OCCCurve*>(myCurve))
    result_list.append(curve);
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

