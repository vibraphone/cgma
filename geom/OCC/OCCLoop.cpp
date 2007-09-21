//-------------------------------------------------------------------------
// Filename      : OCCLoop.cpp
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
#include "OCCLoop.hpp"
#include "OCCQueryEngine.hpp"
#include "CoEdgeSM.hpp"

#include "OCCBody.hpp"
#include "OCCLump.hpp"
#include "OCCShell.hpp"
#include "OCCSurface.hpp"
#include "OCCCoEdge.hpp"
#include "OCCCurve.hpp"
#include "OCCPoint.hpp"

#include <TopExp.hxx>
#include <TopTools_IndexedMapOfShape.hxx>
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
OCCLoop::OCCLoop( TopoDS_Wire *theWire )
{
  myTopoDSWire = theWire;
}
OCCLoop::OCCLoop( Surface *surf_ptr,
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
OCCLoop::OCCLoop( DLIList<CoEdgeSM*> &coedge_list )
{
  myCoEdges += coedge_list;
}

//-------------------------------------------------------------------------
// Purpose       : The default destructor.
//
// Special Notes :
//
//-------------------------------------------------------------------------
OCCLoop::~OCCLoop()
{
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
void OCCLoop::append_simple_attribute_virt(CubitSimpleAttrib* /*csattrib_ptr*/)
{
  //PRINT_ERROR("OCCLoop::append_simple_attribute_virt not defined\n");
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
void OCCLoop::remove_simple_attribute_virt(CubitSimpleAttrib* /*csattrib_ptr*/)
{
  //PRINT_ERROR("OCCLoop::remove_simple_attribute_virt not defined\n");
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
void OCCLoop::remove_all_simple_attribute_virt()
{
  //PRINT_ERROR(" OCCLoop::remove_all_simple_attribute_virt not defined\n");
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
CubitStatus OCCLoop::get_simple_attribute(DLIList<CubitSimpleAttrib*>&
                                                 /*cubit_simple_attrib_list*/)
{
  //PRINT_ERROR("OCCLoop::get_simple_attribute not defined\n");
  return CUBIT_FAILURE;
}
CubitStatus OCCLoop::get_simple_attribute(const CubitString&,
                                              DLIList<CubitSimpleAttrib*>&)
  { return CUBIT_FAILURE; }

//-------------------------------------------------------------------------
// Purpose       : compute bounding box of loop
//
// Special Notes :
//
//-------------------------------------------------------------------------
CubitBox OCCLoop::bounding_box() const
{
   CubitBox box;
   PRINT_ERROR("OCCLoop::bounding_box not implemented\n");
   return box;
}

//-------------------------------------------------------------------------
// Purpose       : Get geometry modeling engine: OCCQueryEngine
//
// Special Notes :
//
//-------------------------------------------------------------------------
GeometryQueryEngine* OCCLoop::get_geometry_query_engine() const
{
   return OCCQueryEngine::instance();
}                

/*
void OCCLoop::bodysms(DLIList<BodySM*> &bodies)
{
  mySurface->bodysms(bodies);
}

void OCCLoop::lumps(DLIList<Lump*> &lumps)
{
  mySurface->lumps(lumps);
}

void OCCLoop::shellsms(DLIList<ShellSM*> &shellsms)
{
  mySurface->shellsms(shellsms);
}

void OCCLoop::surfaces(DLIList<Surface*> &surfaces)
{
  surfaces.append_unique( mySurface );
}

void OCCLoop::loopsms(DLIList<LoopSM*> &loopsms)
{
  loopsms.append_unique( this );
}
void OCCLoop::coedgesms(DLIList<CoEdgeSM*> &coedgesms)
{
  int ii;
  for ( ii = myCoEdges.size(); ii > 0; ii-- )
  {
    coedgesms.append_unique(myCoEdges.get_and_step());
  }
}

void OCCLoop::curves(DLIList<Curve*> &curves)
{
  int ii;
  for ( ii = myCoEdges.size(); ii > 0; ii-- )
  {
    myCoEdges.get_and_step()->curves(curves);
  }
}
void OCCLoop::points(DLIList<Point*> &points)
{
  int ii;
  for ( ii = myCoEdges.size(); ii > 0; ii-- )
  {
    myCoEdges.get_and_step()->points(points);
  }
}
*/


void OCCLoop::get_parents_virt( DLIList<TopologyBridge*>& parents )
  { parents.append( mySurface ); }
void OCCLoop::get_children_virt( DLIList<TopologyBridge*>& children )
{
  TopTools_IndexedMapOfShape M;
  TopExp::MapShapes(*myTopoDSWire, TopAbs_EDGE, M);
  int ii;
  for (ii=1; ii<=M.Extent(); ii++) {
	  TopologyBridge *curve = OCCQueryEngine::occ_to_cgm(M(ii));
	  children.append_unique(curve);
  }
}


void OCCLoop::get_lumps( DLIList<OCCLump*>& result_list )
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

void OCCLoop::get_shells( DLIList<OCCShell*>& result_list )
{
  if ( OCCSurface* surf = dynamic_cast<OCCSurface*>(mySurface) )
    surf->get_shells( result_list );
}

void OCCLoop::get_coedges( DLIList<OCCCoEdge*>& result_list )
{
  TopTools_IndexedMapOfShape M;
  TopExp::MapShapes(*myTopoDSWire, TopAbs_EDGE, M);
  int ii;
  for (ii=1; ii<=M.Extent(); ii++) {
	  TopologyBridge *curve = OCCQueryEngine::occ_to_cgm(M(ii));
	  result_list.append_unique(dynamic_cast<OCCCoEdge*>(curve));
  }
}

void OCCLoop::get_curves( DLIList<OCCCurve*>& result_list )
{
  DLIList<OCCCoEdge*> coedge_list;
  get_coedges( coedge_list );
  coedge_list.reset();
  for ( int i = coedge_list.size(); i--; )
  {
    OCCCoEdge* coedge = coedge_list.get_and_step();
    OCCCurve* curve = dynamic_cast<OCCCurve*>(coedge->curve());
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
void OCCLoop::disconnect_all_coedges()
{
  myCoEdges.reset();
  for (int i = myCoEdges.size(); i--; )
  {
    CoEdgeSM* sm_ptr = myCoEdges.get_and_step();
    OCCCoEdge* coedge = dynamic_cast<OCCCoEdge*>(sm_ptr);
    if (coedge)
    {
      assert(coedge->get_loop() == this);
      coedge->remove_loop();
    }
  }
  myCoEdges.clean_out();
}


