//-------------------------------------------------------------------------
// Filename      : FacetCoEdge.cpp
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

// ********** BEGIN CUBIT INCLUDES         **********
#include "CubitDefines.h"
#include "CastTo.hpp"
#include "FacetCoEdge.hpp"
#include "FacetLoop.hpp"
#include "FacetQueryEngine.hpp"
#include "CubitUtil.hpp"

#include "FacetBody.hpp"
#include "FacetLump.hpp"
#include "FacetShell.hpp"
#include "FacetSurface.hpp"
#include "FacetCurve.hpp"
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
FacetCoEdge::FacetCoEdge( Curve *curv_ptr, LoopSM *loop_ptr, CubitSense sense )
{
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
FacetCoEdge::FacetCoEdge( Curve *curv_ptr, CubitSense sense )
{
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
FacetCoEdge::~FacetCoEdge()
{
}

//-------------------------------------------------------------------------
// Purpose       : Get geometry modeling engine: FacetQueryEngine
//
// Special Notes :
//
// Creator       : Steve Owen
//
// Creation Date : 07/18/00
//-------------------------------------------------------------------------
GeometryQueryEngine* FacetCoEdge::get_geometry_query_engine() const
{
  return FacetQueryEngine::instance();
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
void FacetCoEdge::append_simple_attribute_virt(const CubitSimpleAttrib& /*csattrib_ptr*/)
{
  //PRINT_ERROR("FacetCoEdge::append_simple_attribute_virt not implemented\n");
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
void FacetCoEdge::remove_simple_attribute_virt(const CubitSimpleAttrib& /*csattrib_ptr*/)
{
  //PRINT_ERROR("FacetCoEdge::remove_simple_attribute_virt not implemented\n");
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
void FacetCoEdge::remove_all_simple_attribute_virt()
{
  //PRINT_ERROR("FacetCoEdge::remove_all_simple_attribute_virt not implemented\n");
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
CubitStatus FacetCoEdge::get_simple_attribute(DLIList<CubitSimpleAttrib>&
                                              /*cubit_simple_attrib_list*/)
{
  //PRINT_ERROR("FacetCoEdge::get_simple_attribute not implemented\n");
  return CUBIT_FAILURE;
}
CubitStatus FacetCoEdge::get_simple_attribute(const CubitString&,
                                              DLIList<CubitSimpleAttrib>&)
  { return CUBIT_FAILURE; }


CubitSense FacetCoEdge::sense()
{
  return edgeSense;
}
/*
void FacetCoEdge::bodysms(DLIList<BodySM*> &bodies)
{
  myLoop->bodysms(bodies);
}

void FacetCoEdge::lumps(DLIList<Lump*> &lumps)
{
  myLoop->lumps(lumps);
}

void FacetCoEdge::shellsms(DLIList<ShellSM*> &shellsms)
{
  myLoop->shellsms(shellsms);
}

void FacetCoEdge::surfaces(DLIList<Surface*> &surfaces)
{
  myLoop->surfaces( surfaces );
}

void FacetCoEdge::loopsms(DLIList<LoopSM*> &loopsms)
{
  loopsms.append_unique( myLoop );
}

void FacetCoEdge::coedgesms(DLIList<CoEdgeSM*> &coedgesms)
{
  coedgesms.append_unique( this );
}

void FacetCoEdge::curves(DLIList<Curve*> &curves)
{
  curves.append_unique( myCurve );
}

void FacetCoEdge::points(DLIList<Point*> &points)
{
  myCurve->points( points );
}
*/

void FacetCoEdge::get_parents_virt( DLIList<TopologyBridge*>& parents )
  { parents.append( myLoop ); }

void FacetCoEdge::get_children_virt( DLIList<TopologyBridge*>& children )
  { children.append( myCurve ); }

void FacetCoEdge::reverse_sense()
{
  edgeSense = CubitUtil::opposite_sense( edgeSense );
}


void FacetCoEdge::get_lumps( DLIList<FacetLump*>& result_list )
{
  DLIList<FacetShell*> shell_list;
  get_shells( shell_list );
  shell_list.reset();
  for ( int i = shell_list.size(); i--; )
  {
    FacetShell* shell = shell_list.get_and_step();
    shell->get_lumps( result_list );
    FacetLump* lump = dynamic_cast<FacetLump*>(shell->get_lump());
    if (lump)
      result_list.append_unique(lump);
  }
}

void FacetCoEdge::get_shells( DLIList<FacetShell*>& result_list )
{
  if ( FacetLoop* loop = dynamic_cast<FacetLoop*>(myLoop) )
    if ( FacetSurface* surf = dynamic_cast<FacetSurface*>(loop->get_surface()) )
      surf->get_shells( result_list );
}

void FacetCoEdge::get_curves( DLIList<FacetCurve*>& result_list )
{
  if (FacetCurve* curve = dynamic_cast<FacetCurve*>(myCurve))
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

