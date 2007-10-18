//-------------------------------------------------------------------------
// Filename      : OCCBody.cpp
//
// Purpose       : 
//
// Special Notes :
//
// Creator       : David White
//
// Creation Date : 7/18/00
//
//-------------------------------------------------------------------------

// ********** BEGIN STANDARD INCLUDES      **********
#include <assert.h>
// ********** END STANDARD INCLUDES        **********
#include "config.h"

// ********** BEGIN CUBIT INCLUDES         **********
#include "CubitDefines.h"
#include "CubitString.hpp"
#include "CubitPoint.hpp"
#include "CastTo.hpp"
#include "BodySM.hpp"
#include "Body.hpp"
#include "OCCBody.hpp"
#include "CubitSimpleAttrib.hpp"
#include "OCCQueryEngine.hpp"
#include "DLIList.hpp"
#include "FacetEvalTool.hpp"
#include "Surface.hpp"
#include "OCCSurface.hpp"
#include "CubitTransformMatrix.hpp"
#include "OCCPoint.hpp"
#include "OCCCurve.hpp"
#include "OCCCoEdge.hpp"
#include "OCCLoop.hpp"
#include "OCCShell.hpp"
#include "OCCLump.hpp"
#include "CubitFacetEdge.hpp"
#include "OCCModifyEngine.hpp"
#include "OCCAttrib.hpp"

#include <TopExp.hxx>
#include <TopTools_IndexedMapOfShape.hxx>
#include "BRepBuilderAPI_Transform.hxx"
#include "gp_Ax1.hxx"
#include "Bnd_Box.hxx"
#include "BRepBndLib.hxx"

//-------------------------------------------------------------------------
// Purpose       : A constructor with a list of lumps that are attached.
//
// Special Notes :
//
//-------------------------------------------------------------------------
OCCBody::OCCBody(TopoDS_Shape *theShape)
{
  myTopoDSShape = theShape;
  update_bounding_box();
}

OCCBody::~OCCBody() 
{
    //Not sure what to do..
}

GeometryQueryEngine* OCCBody::get_geometry_query_engine() const
{
  return OCCQueryEngine::instance();
}

void OCCBody::append_simple_attribute_virt(CubitSimpleAttrib *csa)
  { attribSet.append_attribute(csa); }
  
void OCCBody::remove_simple_attribute_virt(CubitSimpleAttrib *csa)
  { attribSet.remove_attribute(csa); }

void OCCBody::remove_all_simple_attribute_virt()
  { attribSet.remove_all_attributes(); }
  
CubitStatus OCCBody::get_simple_attribute(DLIList<CubitSimpleAttrib*>& csa_list)
  { return attribSet.get_attributes(csa_list); }

CubitStatus OCCBody::get_simple_attribute( const CubitString& name,
                                          DLIList<CubitSimpleAttrib*>& csa_list )
  { return attribSet.get_attributes( name, csa_list ); }

CubitStatus OCCBody::save_attribs( FILE *file_ptr )
  { return attribSet.save_attributes( file_ptr); }

CubitStatus OCCBody::restore_attribs( FILE *file_ptr, unsigned int endian )
  { return attribSet.restore_attributes( file_ptr, endian ); }



//----------------------------------------------------------------
// Function: copy
// Description: create a new copy of the body.
// Author: sjowen
//----------------------------------------------------------------
BodySM* OCCBody::copy()
{
  return (BodySM*)NULL;
}
//---------------------------------------------------------------- 
// Function: can_be_deleted 
// Description: determine if the body can be deleted 
// 
// Author: sjowen 
//---------------------------------------------------------------- 
CubitBoolean OCCBody::can_be_deleted( DLIList <Body*> &body_list ) 
{ 
  CubitBoolean delete_ok = CUBIT_TRUE; 
  DLIList<OCCSurface *>surf_list; 
  //get_surfaces(surf_list); 
  int ii; 
  for (ii=0; ii<surf_list.size() && delete_ok; ii++) 
  { 
    OCCSurface *surf_ptr = surf_list.get_and_step(); 
    DLIList<OCCBody*>my_body_list; 
    surf_ptr->get_bodies(my_body_list); 
    int jj; 
    if (my_body_list.size() >= 2) 
    { 
      for (jj=0; jj<my_body_list.size() && delete_ok; jj++) 
      { 
        BodySM *my_body_ptr = my_body_list.get_and_step(); 
        if (my_body_ptr != this) 
        { 
          int kk; 
          int found = 0; 
          for (kk=0; kk<body_list.size() && !found; kk++) 
          { 
            Body *body_ptr = body_list.get_and_step(); 
            OCCBody* fbody_ptr = CAST_TO(body_ptr->get_body_sm_ptr(), OCCBody); 
            if (fbody_ptr) 
            { 
              if (my_body_ptr == fbody_ptr) 
                found = 1; 
            } 
          } 
          if (!found) 
          { 
            delete_ok = CUBIT_FALSE; 
            PRINT_ERROR("Body cannot be deleted because it is merged with adjacent Body\n"); 
            PRINT_INFO("    Mesh Based Geometry entities cannot be unmerged.\n" 
              "    Try using the no_merge option when importing the mesh\n"); 
          } 
        } 
      } 
    } 
  } 
  return delete_ok; 
} 
    
//----------------------------------------------------------------
// Function: move
// Description: translate the body and its child entities
//
// Author: sjowen
//----------------------------------------------------------------
CubitStatus OCCBody::move(double dx, double dy, double dz)
{
  CubitTransformMatrix tfmat;
  tfmat.translate( dx, dy, dz );

  CubitStatus stat = transform( tfmat, CUBIT_FALSE );

  if (stat == CUBIT_SUCCESS)
    myTransforms.translate( dx, dy, dz );

  return stat;
}


//----------------------------------------------------------------
// Function: rotate
// Description: rotate the body and its child entities
//
// Author: sjowen
//----------------------------------------------------------------
CubitStatus OCCBody::rotate( double x, double y, double z, 
                               double angle_in_degrees )
{

  CubitTransformMatrix rotmat;
  CubitVector axis( x, y, z );
  rotmat.rotate( angle_in_degrees, axis );

  CubitStatus stat = transform( rotmat, CUBIT_TRUE );

  if (stat == CUBIT_SUCCESS)
    myTransforms.rotate( angle_in_degrees, axis );

  return stat;
}

//----------------------------------------------------------------
// Function: scale
// Description: scale the body and its child entities
//              use a constant scale factor
//
// Author: sjowen
//----------------------------------------------------------------
CubitStatus OCCBody::scale(double scale_factor )
{
  return scale(scale_factor,scale_factor,scale_factor);
}

//----------------------------------------------------------------
// Function: scale
// Description: scale the body and its child entities
//
// Author: sjowen
//----------------------------------------------------------------
CubitStatus OCCBody::scale(double scale_factor_x,
                             double scale_factor_y,
                             double scale_factor_z )
{
  CubitTransformMatrix scalemat;
  scalemat.scale_about_origin( scale_factor_x, 
                               scale_factor_y, 
                               scale_factor_z );

  CubitStatus stat = transform( scalemat, CUBIT_FALSE );

  if (stat == CUBIT_SUCCESS)
    myTransforms.scale_about_origin( scale_factor_x, 
                                     scale_factor_y, 
                                     scale_factor_z );

  // scale the facetcurve

  DLIList<OCCCurve *> curve_list;
  //get_curves(curve_list); 
  Curve *curv_ptr;
  for (int ii=0; ii<curve_list.size(); ii++)
  {
    curv_ptr = curve_list.get_and_step();
    OCCCurve *fcurve = CAST_TO( curv_ptr, OCCCurve );
    if (fcurve)
    {
      fcurve->reset_length();
    }
  }

  return stat;
}

//----------------------------------------------------------------
// Function: restore
// Description: restore the body and its child entities
//              to its original coordinates using the inverse
//              transformation matrix
//
// Author: sjowen
//----------------------------------------------------------------
CubitStatus OCCBody::restore()
{
  // invert the transformation matrix and apply to entities 
  // (assumes an orthogonal matrix (ie. no shear or non-uniform scaling)

  CubitTransformMatrix inverse_mat;
  inverse_mat = myTransforms.inverse();

  CubitStatus stat = transform( inverse_mat, CUBIT_TRUE );

  if (stat == CUBIT_SUCCESS)
    myTransforms.set_to_identity();

  return stat;
}

//----------------------------------------------------------------
// Function: reflect
// Description: reflect the body about a exis
//
// Author: sjowen
//----------------------------------------------------------------
CubitStatus OCCBody::reflect( double reflect_axis_x,
                                double reflect_axis_y,
                                double reflect_axis_z )
{
  gp_Pnt aOrigin(0,0,0);
  gp_Dir aDir(reflect_axis_x, reflect_axis_y,reflect_axis_z); 
  gp_Ax1 anAxis(aOrigin, aDir);

  gp_Trsf aTrsf;
  aTrsf.SetMirror(anAxis);

  BRepBuilderAPI_Transform aBRepTrsf(*myTopoDSShape, aTrsf); 
  update_bounding_box();
  return CUBIT_SUCCESS;
}

//----------------------------------------------------------------
// Function: update_bounding_box
// Description: calculate for bounding box of this OCCBody
//
// Author: janehu
//----------------------------------------------------------------
void OCCBody::update_bounding_box() 
{
  Bnd_Box box;
  const TopoDS_Shape shape=*myTopoDSShape;
  //calculate the bounding box
  BRepBndLib::Add(shape, box);
  double min[3], max[3];
  
  //get values
  box.Get(min[0], min[1], min[2], max[0], max[1], max[2]);

  //update boundingbox.
  CubitBox cBox(min, max);
  boundingbox = cBox;
}

//----------------------------------------------------------------
// Function: get_bounding_box
// Description: get the  bounding box of this OCCBody
//
// Author: janehu
//----------------------------------------------------------------
CubitBox OCCBody::get_bounding_box()
{
  return boundingbox ;
}

//----------------------------------------------------------------
// Function: transform
// Description: transform the body based on a transformation matrix
//              main function for applying transformations to 
//              facet-based bodies
//
// Author: sjowen
//----------------------------------------------------------------
CubitStatus OCCBody::transform( CubitTransformMatrix &tfmat, 
                                  CubitBoolean is_rotation )
{
  return CUBIT_SUCCESS;
}

CubitStatus OCCBody::get_transforms( CubitTransformMatrix &tfm ) 
{
  tfm = myTransforms;
  return CUBIT_SUCCESS;
}

CubitStatus OCCBody::set_transforms( CubitTransformMatrix tfm ) 
{
  myTransforms = tfm;
  return CUBIT_SUCCESS;
}

int OCCBody::validate(const CubitString &, DLIList <TopologyEntity*>&)
{
  PRINT_ERROR("This option is not available for mesh defined geometry.\n");
  return 0;
}

void OCCBody::get_parents_virt( DLIList<TopologyBridge*>& ) 
  {}
  
void OCCBody::get_children_virt( DLIList<TopologyBridge*>& lumps )
{
  TopTools_IndexedMapOfShape M;
  TopExp::MapShapes(*myTopoDSShape, TopAbs_SOLID, M);
  int ii;
  for (ii=1; ii<=M.Extent(); ii++) {
	  TopologyBridge *lump = OCCQueryEngine::occ_to_cgm(M(ii));
	  lumps.append_unique(lump);
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
void OCCBody::disconnect_all_lumps()
{
  myLumps.reset();
  for (int i = myLumps.size(); i--; )
  {
    Lump* sm_ptr = myLumps.get_and_step();
    OCCLump* lump = dynamic_cast<OCCLump*>(sm_ptr);
    if (lump)
    {
      assert(lump->get_body() == this);
      lump->remove_body();
    }
  }
  myLumps.clean_out();
}

//-------------------------------------------------------------------------
// Purpose       : Find centroid
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 05/10/04
//-------------------------------------------------------------------------
CubitStatus OCCBody::mass_properties( CubitVector& centroid, 
                                        double& volume )
{
  centroid.set( 0.0, 0.0, 0.0 );
  volume = 0.0;
  
  DLIList<OCCLump*> lumps (myLumps.size());
  CAST_LIST( myLumps, lumps, OCCLump );
  assert( myLumps.size() == lumps.size() );
  for (int i = lumps.size(); i--; )
  {
    CubitVector cent;
    double vol;
    if (CUBIT_SUCCESS != lumps.get_and_step()->mass_properties(cent,vol))
      return CUBIT_FAILURE;
    centroid += vol*cent;
    volume += vol;
  }
  if (volume > CUBIT_RESABS)
  {
    centroid /= volume;
  }
  else
  {
    centroid.set( 0.0, 0.0, 0.0 );
    volume = 0.0;
  }

  return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Purpose       : Used to be OCCQueryEngine::is_point_in_body
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 05/10/04
//-------------------------------------------------------------------------
CubitPointContainment OCCBody::point_containment( const CubitVector &point )
{
  CubitPointContainment pc_value; 
  OCCLump *facet_lump;

  int i;
  for(i=myLumps.size(); i--;)
  {
    facet_lump = dynamic_cast<OCCLump*>(myLumps.get_and_step()); 
    pc_value = facet_lump->point_containment( point );
    if( pc_value == CUBIT_PNT_INSIDE )
      return CUBIT_PNT_INSIDE;
    else if( pc_value == CUBIT_PNT_BOUNDARY )
      return CUBIT_PNT_BOUNDARY;
  }

  return CUBIT_PNT_OUTSIDE;
}


