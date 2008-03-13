//-------------------------------------------------------------------------
// Filename      : OCCShell.cpp
//
// Purpose       : 
//
// Special Notes :
//
// Creator       : David White
//
// Creation Date : 7/18/2000
//
//-------------------------------------------------------------------------

// ********** BEGIN STANDARD INCLUDES      **********
#include <stddef.h>
// ********** END STANDARD INCLUDES        **********
#include "config.h"
// ********** BEGIN CUBIT INCLUDES         **********
#include "CastTo.hpp"
#include "CubitUtil.hpp"

#include "OCCQueryEngine.hpp"
#include "OCCShell.hpp"
#include "ShellSM.hpp"
#include "Lump.hpp"
#include "Surface.hpp"

#include "OCCBody.hpp"
#include "OCCLump.hpp"
#include "OCCSurface.hpp"
#include "GfxDebug.hpp"

#include <TopExp.hxx>
#include <TopTools_IndexedMapOfShape.hxx>
#include "TopTools_IndexedDataMapOfShapeListOfShape.hxx"
#include "TopoDS.hxx"
#include "TopTools_ListIteratorOfListOfShape.hxx"
#include "TopTools_DataMapOfShapeInteger.hxx"
#include "TopTools_ListOfShape.hxx"
#include "BRepBuilderAPI_Transform.hxx"
// ********** END CUBIT INCLUDES           **********

// ********** BEGIN STATIC DECLARATIONS    **********
// ********** END STATIC DECLARATIONS      **********

// ********** BEGIN PUBLIC FUNCTIONS       **********


//-------------------------------------------------------------------------
// Purpose       : A constructor with a pointer to an OCC SHELL.
//
// Special Notes :
//
//-------------------------------------------------------------------------
OCCShell::OCCShell(TopoDS_Shell *theShell, OCCSurface* surface)
{
  myTopoDSShell = theShell;
  mySheetSurface = surface;
}

//-------------------------------------------------------------------------
// Purpose       : The destructor.
//
// Special Notes :
//
//-------------------------------------------------------------------------
OCCShell::~OCCShell()
{
}

//-------------------------------------------------------------------------
// Purpose       : Get geometry modeling engine: OCCQueryEngine
//
// Special Notes :
//
//-------------------------------------------------------------------------
GeometryQueryEngine* 
                 OCCShell::get_geometry_query_engine() const
{
   return OCCQueryEngine::instance();
}                 

void OCCShell::append_simple_attribute_virt(CubitSimpleAttrib*)
{
}
void OCCShell::remove_simple_attribute_virt(CubitSimpleAttrib* )
{
}
void OCCShell::remove_all_simple_attribute_virt()
{
}
CubitStatus OCCShell::get_simple_attribute(DLIList<CubitSimpleAttrib*>&)
{
  return CUBIT_FAILURE;
}
CubitStatus OCCShell::get_simple_attribute(const CubitString&,
                                              DLIList<CubitSimpleAttrib*>&)
  { return CUBIT_FAILURE; }

//-------------------------------------------------------------------------
// Purpose       : Query solid modeler topology
//
// Special Notes : 
//
// Author        : Jane Hu 
//
// Creation Date : 11/28/07
//-------------------------------------------------------------------------
void OCCShell::get_parents_virt( DLIList<TopologyBridge*>& parents ) 
{ 
  if(mySheetSurface)
  {
    parents.append(mySheetSurface->my_body());
    return;
  }

  OCCQueryEngine* oqe = (OCCQueryEngine*) get_geometry_query_engine();
  OCCBody * body = NULL;
  DLIList <OCCBody* > *bodies = oqe->BodyList;
  TopTools_IndexedDataMapOfShapeListOfShape M;
  for(int i = 0; i <  bodies->size(); i++)
  {
     body = bodies->get_and_step();
     TopExp::MapShapesAndAncestors(*(body->get_TopoDS_Shape()),
				   TopAbs_SHELL, TopAbs_SOLID, M);
     if (!M.Contains(*(get_TopoDS_Shell())))
	continue;

     const TopTools_ListOfShape& ListOfShapes = 
				M.FindFromKey(*(get_TopoDS_Shell()));
     if (!ListOfShapes.IsEmpty()) 
     {
         TopTools_ListIteratorOfListOfShape it(ListOfShapes) ;
         for (;it.More(); it.Next())
         {
	   TopoDS_Solid Solid = TopoDS::Solid(it.Value());
           int k = oqe->OCCMap->Find(Solid);
	   parents.append((OCCLump*)(oqe->OccToCGM->find(k))->second);
	 }
     } 
  }
}


void OCCShell::get_children_virt( DLIList<TopologyBridge*>& children )
{
  if(mySheetSurface)
  {
    children.append(mySheetSurface);
    return;
  }
  TopTools_IndexedMapOfShape M;
  TopExp::MapShapes(*myTopoDSShell, TopAbs_FACE, M);
  int ii;
  for (ii=1; ii<=M.Extent(); ii++) {
	  TopologyBridge *surface = OCCQueryEngine::instance()->occ_to_cgm(M(ii));
	  children.append_unique(surface);
  }
}

//----------------------------------------------------------------
// Function: to update the core Shell
//           for any movement of the body.
// Author: Jane Hu
//----------------------------------------------------------------
CubitStatus OCCShell::update_OCC_entity( BRepBuilderAPI_Transform &aBRepTrsf)
{
  if(mySheetSurface)
    return CUBIT_FAILURE;

  TopoDS_Shape shape = aBRepTrsf.ModifiedShape(*get_TopoDS_Shell());
  TopoDS_Shell shell = TopoDS::Shell(shape);

  int k = OCCQueryEngine::instance()->OCCMap->Find(*myTopoDSShell);
  assert (k > 0 && k <= OCCQueryEngine::instance()->iTotalTBCreated);
  OCCQueryEngine::instance()->OCCMap->UnBind(*myTopoDSShell);
  OCCQueryEngine::instance()->OCCMap->Bind(shell, k);

  //set the lumps
  DLIList<TopologyBridge *> surfaces;
  this->get_children_virt(surfaces);
  for (int i = 1; i <= surfaces.size(); i++)
  {
     OCCSurface *surface = CAST_TO(surfaces.get_and_step(), OCCSurface);
     surface->update_OCC_entity(aBRepTrsf);
  }
  set_TopoDS_Shell(shell);
  return CUBIT_SUCCESS;
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

