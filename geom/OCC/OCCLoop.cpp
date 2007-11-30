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
#include "config.h"

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
#include "TopExp_Explorer.hxx"
#include "TopoDS.hxx"
#include "TopTools_ListIteratorOfListOfShape.hxx"
#include "TopTools_DataMapOfShapeInteger.hxx"
#include "TopTools_IndexedDataMapOfShapeListOfShape.hxx"
// ********** END CUBIT INCLUDES           **********

// ********** BEGIN STATIC DECLARATIONS    **********
// ********** END STATIC DECLARATIONS      **********

// ********** BEGIN PUBLIC FUNCTIONS       **********

//-------------------------------------------------------------------------
// Purpose       : The constructor with a pointer to the TopoDS_Wire.
//
// Special Notes :
//
//-------------------------------------------------------------------------
OCCLoop::OCCLoop( TopoDS_Wire *theWire )
{
  myTopoDSWire = theWire;
}

//-------------------------------------------------------------------------
// Purpose       : The default destructor.
//
// Special Notes :
//
//-------------------------------------------------------------------------
OCCLoop::~OCCLoop()
{
  disconnect_all_curves();
}

//-------------------------------------------------------------------------
// Purpose       : Tear down topology
//
// Special Notes :
//
// Creator       : Jane Hu
//
// Creation Date : 11/29/07
//-------------------------------------------------------------------------
void OCCLoop::disconnect_all_curves()
{
  myCoEdgeList.reset();
  for (int i = myCoEdgeList.size(); i--; )
  {
    OCCCoEdge* coedge = myCoEdgeList.get_and_step();
    assert(coedge->loop() == this);
    OCCCurve* curve = CAST_TO(coedge->curve(), OCCCurve);
    curve->remove_loop(this);
  }
  myCoEdgeList.clean_out();
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

void OCCLoop::get_parents_virt( DLIList<TopologyBridge*>& parents )
{
  OCCQueryEngine* oqe = (OCCQueryEngine*) get_geometry_query_engine();
  OCCSurface * surf = NULL;
  DLIList <OCCSurface* > *surfs = oqe->SurfaceList;
  TopTools_IndexedDataMapOfShapeListOfShape M;
  for(int i = 0; i <  surfs->size(); i++)
  {
     surf = surfs->get_and_step();
     TopExp::MapShapesAndAncestors(*(surf->get_TopoDS_Face()),
                                   TopAbs_WIRE, TopAbs_FACE, M);
     const TopTools_ListOfShape& ListOfShapes =
                                M.FindFromKey(*(get_TopoDS_Wire()));
     if (!ListOfShapes.IsEmpty())
     {
         TopTools_ListIteratorOfListOfShape it(ListOfShapes) ;
         for (;it.More(); it.Next())
         {
           TopoDS_Face Face = TopoDS::Face(it.Value());
           int k = oqe->OCCMap->Find(Face);
           parents.append((OCCSurface*)(oqe->OccToCGM->find(k))->second);
         }
     }
  }
}

void OCCLoop::get_children_virt( DLIList<TopologyBridge*>& children )
{
  TopTools_IndexedMapOfShape M;
  TopExp::MapShapes(*myTopoDSWire, TopAbs_EDGE, M);
  int ii;
  for (ii=1; ii<=M.Extent(); ii++) {
          TopologyBridge *curve = OCCQueryEngine::instance()->occ_to_cgm(M(ii));
          children.append_unique(curve);
  }
}

//-------------------------------------------------------------------------
// Purpose       : compute bounding box of loop
//
// Special Notes :
//
//-------------------------------------------------------------------------
CubitBox OCCLoop::bounding_box() const
{
   CubitBox box;
   for (int i = myCoEdgeList.size(); i > 0; i--)
   {
      DLIList<OCCCoEdge*> coedges = myCoEdgeList;
      OCCCoEdge* coedge = coedges.get_and_step();
      box |= coedge->curve()->bounding_box();
   }
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



