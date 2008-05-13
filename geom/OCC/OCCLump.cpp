//-------------------------------------------------------------------------
// Filename      : OCCLump.cpp
//
// Purpose       : 
//
// Creator       : David White
//
// Creation Date : 7/18/2000
//
//-------------------------------------------------------------------------

// ********** BEGIN STANDARD INCLUDES      **********
#include <assert.h>
// ********** END STANDARD INCLUDES        **********
#include "config.h"
// ********** BEGIN CUBIT INCLUDES         **********
#include "OCCQueryEngine.hpp"
#include "OCCLump.hpp"
#include "CastTo.hpp"
#include "Surface.hpp"
#include "DLIList.hpp"
#include "CubitFacet.hpp"
#include "CubitPoint.hpp"
#include "CubitVector.hpp"
#include "CubitString.hpp"
#include "ShellSM.hpp"
#include "BodySM.hpp"
#include "Body.hpp"

#include "OCCBody.hpp"
#include "OCCShell.hpp"
#include "OCCSurface.hpp"
#include "OCCLoop.hpp"
#include "OCCCoEdge.hpp"
#include "OCCCurve.hpp"
#include "OCCPoint.hpp"

#include <TopExp.hxx>
#include <TopTools_IndexedMapOfShape.hxx>
#include "BRepBuilderAPI_Transform.hxx"
#include "TopTools_DataMapOfShapeInteger.hxx"
#include "BRepAlgoAPI_BooleanOperation.hxx"
#include "Bnd_Box.hxx"
#include "BRepBndLib.hxx"
#include "GProp_GProps.hxx"
#include "BRepGProp.hxx"
#include "BOP_SolidClassifier.hxx"
#include "TopExp_Explorer.hxx"
#include "TopoDS.hxx"
#include "BRep_Tool.hxx"
// ********** END CUBIT INCLUDES           **********

// ********** BEGIN FORWARD DECLARATIONS   **********
class RefVolume;
// ********** END FORWARD DECLARATIONS     **********

// ********** BEGIN STATIC DECLARATIONS    **********
// ********** END STATIC DECLARATIONS      **********

// ********** BEGIN PUBLIC FUNCTIONS       **********

//-------------------------------------------------------------------------
// Purpose       : The constructor with a pointer to TopoDS_Solid
//
// Special Notes :
//
//-------------------------------------------------------------------------
OCCLump::OCCLump(TopoDS_Solid *theSolid, OCCSurface* surface, OCCShell* shell)
{
  myTopoDSSolid = theSolid;
  mySheetSurface = surface;
  myShell = shell;
}

OCCLump::~OCCLump()
{ 
  if (myTopoDSSolid)
    delete myTopoDSSolid;
}

void OCCLump::set_TopoDS_Solid(TopoDS_Solid solid)
{
  if(myTopoDSSolid)
    *myTopoDSSolid = solid;
  else
    myTopoDSSolid = new TopoDS_Solid(solid);
}
//-------------------------------------------------------------------------
// Purpose       : Find centroid
//
// Special Notes :
//
// Author       : Jane Hu
//
// Creation Date : 12/06/07
//-------------------------------------------------------------------------
CubitStatus OCCLump::mass_properties( CubitVector& centroid,
                                      double& volume )
{
  if (mySheetSurface || myShell)
    return CUBIT_FAILURE;

  GProp_GProps myProps;
  BRepGProp::VolumeProperties(*myTopoDSSolid, myProps);
  volume = myProps.Mass();
  gp_Pnt pt = myProps.CentreOfMass();
  centroid.set(pt.X(), pt.Y(), pt.Z());

  return CUBIT_SUCCESS;
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
void OCCLump::append_simple_attribute_virt(CubitSimpleAttrib *csa)
  { attribSet.append_attribute(csa); }

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
void OCCLump::remove_simple_attribute_virt(CubitSimpleAttrib *csa )
  { attribSet.remove_attribute(csa); }

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
void OCCLump::remove_all_simple_attribute_virt()
{ attribSet.remove_all_attributes(); }

//-------------------------------------------------------------------------
// Purpose       : The purpose of this function is to get the  
//                 attributes attached to this geometry entity. The name is 
//                 attached to the underlying BODY this points to.
//
// Special Notes : 
//
//-------------------------------------------------------------------------
CubitStatus OCCLump::get_simple_attribute(DLIList<CubitSimpleAttrib*>& csa_list)
  { return attribSet.get_attributes( csa_list ); }

CubitStatus OCCLump::get_simple_attribute( const CubitString& name,
                                        DLIList<CubitSimpleAttrib*>& csa_list )
  { return attribSet.get_attributes( name, csa_list ); }

CubitStatus OCCLump::save_attribs( FILE *file_ptr )
  { return attribSet.save_attributes( file_ptr ); }

CubitStatus OCCLump::restore_attribs( FILE *file_ptr, unsigned int endian )
  { return attribSet.restore_attributes( file_ptr, endian ); }




//-------------------------------------------------------------------------
// Purpose       : Get the bounding box of the object.
//
// Special Notes :
//
//-------------------------------------------------------------------------
CubitBox OCCLump::bounding_box() const 
{
  Bnd_Box box;
  TopoDS_Shape shape;
  if(mySheetSurface)
    shape = *(mySheetSurface->get_TopoDS_Face());
  else if(myShell)
    shape = *(myShell->get_TopoDS_Shell());
  else
    shape =*myTopoDSSolid;

  //calculate the bounding box
  BRepBndLib::Add(shape, box);
  double min[3], max[3];

  //get values
  box.Get(min[0], min[1], min[2], max[0], max[1], max[2]);

  //update boundingbox.
  CubitBox cBox(min, max);
  return cBox;
}

//-------------------------------------------------------------------------
// Purpose       : Get geometry modeling engine: AcisGeometryEngine
//
// Special Notes :
//
//-------------------------------------------------------------------------
GeometryQueryEngine* 
                 OCCLump::get_geometry_query_engine() const
{
   return OCCQueryEngine::instance();
}                 

//-------------------------------------------------------------------------
// Purpose       : Returns the volume of the Lump
//
// Special Notes :
//
// Creator       : 
//
// Creation Date : 
//-------------------------------------------------------------------------
double OCCLump::measure()
{
  if(mySheetSurface) 
    return mySheetSurface->measure();
  else if(myShell)
    return myShell->measure();

  GProp_GProps myProps;
  BRepGProp::VolumeProperties(*myTopoDSSolid, myProps);
  return myProps.Mass();
}

void OCCLump::get_parents_virt(DLIList<TopologyBridge*> &bodies) 
{
  if(mySheetSurface)
    bodies.append(mySheetSurface->my_body());
  else if (myShell)
    bodies.append(myShell->my_body());
  else
    bodies.append(myBodyPtr);
}

void OCCLump::get_children_virt(DLIList<TopologyBridge*> &shellsms)
{
  if (mySheetSurface)
  {
    shellsms.append(mySheetSurface->my_shell());
    return;
  }
  else if(myShell)
  {
    shellsms.append(myShell);
    return;
  }
  TopTools_IndexedMapOfShape M;
  TopExp::MapShapes(*myTopoDSSolid, TopAbs_SHELL, M);
  int ii;
  for (ii=1; ii<=M.Extent(); ii++) {
	  TopologyBridge *shell = OCCQueryEngine::instance()->occ_to_cgm(M(ii));
          if (shell)
	    shellsms.append_unique(shell);
  }
}



CubitPointContainment OCCLump::point_containment( const CubitVector &point )
{
  if (mySheetSurface || myShell)
    return CUBIT_PNT_UNKNOWN;

  BOP_SolidClassifier ps;
  TopoDS_Solid * solid = get_TopoDS_Solid();
  gp_Pnt pnt(point.x(), point.y(), point.z());
  
  //use face tolerance at the tolerence to see if the point is on.
  TopExp_Explorer Ex;
  Ex.Init(*solid, TopAbs_FACE, TopAbs_SOLID);
  TopoDS_Face face = TopoDS::Face(Ex.Current());
  
  double dtol = BRep_Tool::Tolerance(face);
  TopAbs_State state = ps.Classify(*solid, pnt, dtol); 
  
  if (state == TopAbs_IN)
     return CUBIT_PNT_INSIDE;
  else if (state == TopAbs_OUT)
     return CUBIT_PNT_OUTSIDE;
  else if (state == TopAbs_ON)
     return CUBIT_PNT_BOUNDARY;

  return CUBIT_PNT_UNKNOWN;
  
}

//----------------------------------------------------------------
// Function: to update the core Solid
//           for any movement of the lump.
// Author: Jane Hu
//----------------------------------------------------------------
CubitStatus OCCLump::update_OCC_entity( BRepBuilderAPI_Transform *aBRepTrsf,
                                        BRepAlgoAPI_BooleanOperation *op)
{
  if(mySheetSurface || myShell)
    return CUBIT_FAILURE;

  assert(aBRepTrsf != NULL || op != NULL);
 
  TopoDS_Shape shape;
  if (aBRepTrsf)
    shape = aBRepTrsf->ModifiedShape(*get_TopoDS_Solid());
  else
  {
    TopTools_ListOfShape shapes;
    shapes.Assign(op->Modified(*get_TopoDS_Solid()));
    if(shapes.Extent() > 0)
      shape = shapes.First();
    else if(op->IsDeleted(*get_TopoDS_Solid()))
      ;
    else
      return CUBIT_SUCCESS;
  }
  
  TopoDS_Solid solid;
  if(aBRepTrsf || !op->IsDeleted(*get_TopoDS_Solid()))
    solid = TopoDS::Solid(shape);

  //set the lumps
  DLIList<TopologyBridge *> shells;
  this->get_children_virt(shells);
  for (int i = 1; i <= shells.size(); i++)
  {
     OCCShell *shell = CAST_TO(shells.get_and_step(), OCCShell);
     shell->update_OCC_entity(aBRepTrsf, op);
  }
  OCCQueryEngine::instance()->update_OCC_map(*myTopoDSSolid, solid);

}

//----------------------------------------------------------------
// Function: TopoDS_Shape level function to update the core Lump
//           for any movement  or Boolean operation of the body.
// Author: Jane Hu
//----------------------------------------------------------------
CubitStatus OCCLump::update_OCC_entity(TopoDS_Solid& old_solid,
                                       TopoDS_Shape& new_shape,
                                       BRepAlgoAPI_BooleanOperation *op)
{
  //set the Shells
  TopTools_IndexedMapOfShape M;
  TopoDS_Shape shape;
  TopExp::MapShapes(new_shape, TopAbs_SOLID,M);
  CubitBoolean is_null_new_shape = CUBIT_FALSE;
  if(M.Extent() > 1)
    is_null_new_shape = CUBIT_TRUE;

  M.Clear();
  TopExp::MapShapes(old_solid, TopAbs_SHELL, M);
  TopTools_ListOfShape shapes;
 
  for(int ii=1; ii<=M.Extent(); ii++)
  {
    TopoDS_Shell shell = TopoDS::Shell(M(ii));

    TopTools_ListOfShape shapes;
    shapes.Assign(op->Modified(shell));
    if (shapes.Extent() == 1)
      shape = shapes.First();

    else if(shapes.Extent() > 1)
      shape.Nullify();

    else if(op->IsDeleted(shell))
    {
       TopTools_IndexedMapOfShape M_new;
       TopExp::MapShapes(new_shape, TopAbs_SHELL, M_new);
       if (M_new.Extent()== 1)
         shape = M_new(1);
       else
         shape.Nullify();
    }
    else
    {
       shape = shell;
       continue;
    }
 
    if(shapes.Extent() > 0 || op->IsDeleted(shell))
      OCCShell::update_OCC_entity(shell, shape, op);
  }
  TopoDS_Solid new_solid;
  if(!is_null_new_shape && !op->IsDeleted(old_solid))
    new_solid = TopoDS::Solid(new_shape);
  OCCQueryEngine::instance()->update_OCC_map(old_solid, new_solid);
}
