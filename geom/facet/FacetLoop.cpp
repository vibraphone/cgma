//-------------------------------------------------------------------------
// Filename      : FacetLoop.cpp
//
// Purpose       : Loops for Facet-based geometry module
//
// Special Notes :
//
// Creator       : Steven J. Owen
//
// Creation Date : 12/06/00
//
// Owner         : Steven J. Owen
//-------------------------------------------------------------------------

#include "CastTo.hpp"
#include "FacetLoop.hpp"
#include "FacetQueryEngine.hpp"
#include "CoEdgeSM.hpp"

#include "FacetBody.hpp"
#include "FacetLump.hpp"
#include "FacetShell.hpp"
#include "FacetSurface.hpp"
#include "FacetCoEdge.hpp"
#include "FacetCurve.hpp"
#include "FacetPoint.hpp"
// ********** END CUBIT INCLUDES           **********

// ********** BEGIN STATIC DECLARATIONS    **********
// ********** END STATIC DECLARATIONS      **********

// ********** BEGIN PUBLIC FUNCTIONS       **********

//-------------------------------------------------------------------------
// Purpose       : The constructor with a pointer to the FacetEvalTool.
//
// Special Notes :
//
//-------------------------------------------------------------------------
FacetLoop::FacetLoop( Surface *surf_ptr,
                      DLIList<CoEdgeSM*> &coedge_list )
{
  mySurface = surf_ptr;
  myCoEdges += coedge_list;
}

//-------------------------------------------------------------------------
// Purpose       : The constructor with coedges.
//
// Special Notes :
//
//-------------------------------------------------------------------------
FacetLoop::FacetLoop( DLIList<CoEdgeSM*> &coedge_list )
{
  myCoEdges += coedge_list;
}

//-------------------------------------------------------------------------
// Purpose       : The default destructor.
//
// Special Notes :
//
//-------------------------------------------------------------------------
FacetLoop::~FacetLoop()
{
}

//-------------------------------------------------------------------------
// Purpose       : The purpose of this function is to see if a loop is an external
//                  or internal loop of a surface.
//
// Special Notes : 
//
// Creator       : Jonathan Bugman
//
// Creation Date : 9/9/2008
//-------------------------------------------------------------------------
CubitBoolean FacetLoop::is_external()
{
		  PRINT_ERROR( "This command is not supported with this engine.\n");
          return CUBIT_FAILURE;
}

LoopType FacetLoop::loop_type()
{
  return LOOP_TYPE_UNKNOWN;
}

//-------------------------------------------------------------------------
// Purpose       : The purpose of this function is to append a
//                 attribute to the GE. The name is attached to the
//                 underlying solid model entity this one points to.
//
//
// Special Notes :
//
//-------------------------------------------------------------------------
void FacetLoop::append_simple_attribute_virt(const CubitSimpleAttrib& /*csattrib_ptr*/)
{
  //PRINT_ERROR("FacetLoop::append_simple_attribute_virt not defined\n");
  return;

}

//-------------------------------------------------------------------------
// Purpose       : The purpose of this function is to remove a simple
//                 attribute attached to this geometry entity. The name is
//                 removed from the underlying BODY this points to.
//
// Special Notes :
//
//-------------------------------------------------------------------------
void FacetLoop::remove_simple_attribute_virt(const CubitSimpleAttrib& /*csattrib_ptr*/)
{
  //PRINT_ERROR("FacetLoop::remove_simple_attribute_virt not defined\n");
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
//-------------------------------------------------------------------------
void FacetLoop::remove_all_simple_attribute_virt()
{
  //PRINT_ERROR(" FacetLoop::remove_all_simple_attribute_virt not defined\n");
  return;
}

//-------------------------------------------------------------------------
// Purpose       : The purpose of this function is to get the
//                 attributes attached to this geometry entity. The name is
//                 attached to the underlying BODY this points to.
//
// Special Notes :
//
//-------------------------------------------------------------------------
CubitStatus FacetLoop::get_simple_attribute(DLIList<CubitSimpleAttrib>&
                                                 /*cubit_simple_attrib_list*/)
{
  //PRINT_ERROR("FacetLoop::get_simple_attribute not defined\n");
  return CUBIT_FAILURE;
}
CubitStatus FacetLoop::get_simple_attribute(const CubitString&,
                                              DLIList<CubitSimpleAttrib>&)
  { return CUBIT_FAILURE; }

//-------------------------------------------------------------------------
// Purpose       : compute bounding box of loop
//
// Special Notes :
//
//-------------------------------------------------------------------------
CubitBox FacetLoop::bounding_box() const
{
   CubitBox box;
   //PRINT_ERROR("FacetLoop::bounding_box not implemented\n");
   return box;
}

//-------------------------------------------------------------------------
// Purpose       : Get geometry modeling engine: FacetQueryEngine
//
// Special Notes :
//
//-------------------------------------------------------------------------
GeometryQueryEngine* FacetLoop::get_geometry_query_engine() const
{
   return FacetQueryEngine::instance();
}                

/*
void FacetLoop::bodysms(DLIList<BodySM*> &bodies)
{
  mySurface->bodysms(bodies);
}

void FacetLoop::lumps(DLIList<Lump*> &lumps)
{
  mySurface->lumps(lumps);
}

void FacetLoop::shellsms(DLIList<ShellSM*> &shellsms)
{
  mySurface->shellsms(shellsms);
}

void FacetLoop::surfaces(DLIList<Surface*> &surfaces)
{
  surfaces.append_unique( mySurface );
}

void FacetLoop::loopsms(DLIList<LoopSM*> &loopsms)
{
  loopsms.append_unique( this );
}
void FacetLoop::coedgesms(DLIList<CoEdgeSM*> &coedgesms)
{
  int ii;
  for ( ii = myCoEdges.size(); ii > 0; ii-- )
  {
    coedgesms.append_unique(myCoEdges.get_and_step());
  }
}

void FacetLoop::curves(DLIList<Curve*> &curves)
{
  int ii;
  for ( ii = myCoEdges.size(); ii > 0; ii-- )
  {
    myCoEdges.get_and_step()->curves(curves);
  }
}
void FacetLoop::points(DLIList<Point*> &points)
{
  int ii;
  for ( ii = myCoEdges.size(); ii > 0; ii-- )
  {
    myCoEdges.get_and_step()->points(points);
  }
}
*/


void FacetLoop::get_parents_virt( DLIList<TopologyBridge*>& parents )
  { parents.append( mySurface ); }
void FacetLoop::get_children_virt( DLIList<TopologyBridge*>& children )
  { CAST_LIST_TO_PARENT( myCoEdges, children ); }


void FacetLoop::get_lumps( DLIList<FacetLump*>& result_list )
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

void FacetLoop::get_shells( DLIList<FacetShell*>& result_list )
{
  if ( FacetSurface* surf = dynamic_cast<FacetSurface*>(mySurface) )
    surf->get_shells( result_list );
}

void FacetLoop::get_coedges( DLIList<FacetCoEdge*>& result_list )
{
  myCoEdges.reset();
  for ( int i = 0; i < myCoEdges.size(); i++ )
    if (FacetCoEdge* coedge = dynamic_cast<FacetCoEdge*>(myCoEdges.next(i)))
      result_list.append( coedge );
}

void FacetLoop::get_curves( DLIList<FacetCurve*>& result_list )
{
  DLIList<FacetCoEdge*> coedge_list;
  get_coedges( coedge_list );
  coedge_list.reset();
  for ( int i = coedge_list.size(); i--; )
  {
    FacetCoEdge* coedge = coedge_list.get_and_step();
    FacetCurve* curve = dynamic_cast<FacetCurve*>(coedge->curve());
    if (curve)
      result_list.append_unique(curve);
  }
}


//-------------------------------------------------------------------------
// Purpose       : Tear down topology
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 09/29/03
//-------------------------------------------------------------------------
void FacetLoop::disconnect_all_coedges()
{
  myCoEdges.reset();
  for (int i = myCoEdges.size(); i--; )
  {
    CoEdgeSM* sm_ptr = myCoEdges.get_and_step();
    FacetCoEdge* coedge = dynamic_cast<FacetCoEdge*>(sm_ptr);
    if (coedge)
    {
      assert(coedge->get_loop() == this);
      coedge->remove_loop();
    }
  }
  myCoEdges.clean_out();
}


