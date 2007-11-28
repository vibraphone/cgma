//-------------------------------------------------------------------------
// Filename      : OCCSurface.cpp
//
// Purpose       : 
//
// Special Notes :
//
// Creator       : Alexander Danilov
//
// Creation Date : 
//
// Owner         : 
//-------------------------------------------------------------------------

// ********** BEGIN OCC INCLUDES           **********
#include "config.h"
#include "OCCSurface.hpp"
#include "RefFace.hpp"
#include "OCCQueryEngine.hpp"
#include "OCCAttrib.hpp"

#include "OCCBody.hpp"
#include "OCCLump.hpp"
#include "OCCShell.hpp"
#include "OCCLoop.hpp"
#include "OCCCoEdge.hpp"
#include "OCCCurve.hpp"
#include "OCCPoint.hpp"

// ********** END OCC INCLUDES           **********

// ********** BEGIN CUBIT INCLUDES       **********

#include "CubitSimpleAttrib.hpp"
#include "CubitVector.hpp"
#include "GeometryDefines.h"
#include "CubitUtil.hpp"
#include "CastTo.hpp"
#include "RefVolume.hpp"
#include "GeometryQueryEngine.hpp"
#include "DLIList.hpp"
#include "ShellSM.hpp"
#include "Lump.hpp"
#include "LoopSM.hpp"
#include "CubitPointData.hpp"


// ********** END CUBIT INCLUDES           **********


// ********** BEGIN OpenCascade INCLUDES   **********

#include <BRepAdaptor_Surface.hxx>
#include <TopExp.hxx>
#include <TopTools_IndexedMapOfShape.hxx>
#include <Bnd_Box.hxx>
#include <BndLib_AddSurface.hxx>
#include <Precision.hxx>
#include <TopoDS.hxx>
#include <Extrema_ExtPS.hxx>
#include <BRepLProp_SLProps.hxx>
#include <BRepGProp.hxx>
#include <GProp_GProps.hxx>

// ********** END OpenCascade INCLUDES      **********


// ********** BEGIN STATIC DECLARATIONS    **********
// ********** END STATIC DECLARATIONS      **********


// ********** BEGIN PUBLIC FUNCTIONS       **********
//-------------------------------------------------------------------------
// Purpose       : The constructor with a pointer to the FacetEvalTool. 
//
// Special Notes :
//
//-------------------------------------------------------------------------
OCCSurface::OCCSurface(TopoDS_Face *theFace)
{
  myTopoDSFace = theFace;

/*  printf("Yeah!\n");
  TopoDS_Face face = *myTopoDSFace;
  BRepAdaptor_Surface asurface(face);
  Bnd_Box aBox;
  BndLib_AddSurface::Add(asurface, Precision::Approximation(), aBox);
  double min[3], max[3];
  aBox.Get( min[0], min[1], min[2], max[0], max[1], max[2]);
  printf(".  .  .  box: %lf %lf %lf, %lf %lf %lf\n", min[0], min[1], min[2], max[0], max[1], max[2]);*/

}


OCCSurface::~OCCSurface() 
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
void OCCSurface::append_simple_attribute_virt(CubitSimpleAttrib *csa)
  { attribSet.append_attribute(csa); }


//-------------------------------------------------------------------------
// Purpose       : The purpose of this function is to remove a simple 
//                 attribute attached to this geometry entity. The name is 
//                 removed from the underlying BODY this points to.
//
// Special Notes : 
//
//-------------------------------------------------------------------------
void OCCSurface::remove_simple_attribute_virt(CubitSimpleAttrib *csa)
  { attribSet.remove_attribute( csa ); }


//-------------------------------------------------------------------------
// Purpose       : The purpose of this function is to remove all simple 
//                 attributes attached to this geometry entity.  Also
//                 removes lingering GTC attributes.
//
//
// Special Notes : 
//
//-------------------------------------------------------------------------
void OCCSurface::remove_all_simple_attribute_virt()
  { attribSet.remove_all_attributes(); }


//-------------------------------------------------------------------------
// Purpose       : The purpose of this function is to get the  
//                 attributes attached to this geometry entity. The name is 
//                 attached to the underlying BODY this points to.
//
// Special Notes : 
//
//-------------------------------------------------------------------------
CubitStatus OCCSurface::get_simple_attribute(DLIList<CubitSimpleAttrib*>&
                                                 csa_list)
  { return attribSet.get_attributes(csa_list); }

CubitStatus OCCSurface::get_simple_attribute(const CubitString& name,
                                        DLIList<CubitSimpleAttrib*>& csa_list )
  { return attribSet.get_attributes( name, csa_list ); }


CubitStatus OCCSurface::save_attribs( FILE *file_ptr )
  { return attribSet.save_attributes(file_ptr); }
  
CubitStatus OCCSurface::restore_attribs( FILE *file_ptr, unsigned int endian )
  { return attribSet.restore_attributes( file_ptr, endian ); }



//-------------------------------------------------------------------------
// Purpose       : Get geometry modeling engine: OCCQueryEngine
//
// Special Notes :
//
//-------------------------------------------------------------------------
GeometryQueryEngine* 
                 OCCSurface::get_geometry_query_engine() const
{
   return OCCQueryEngine::instance();
}                 

//-------------------------------------------------------------------------
// Purpose       : Get the bounding box of the object.
//
// Special Notes :
//
//-------------------------------------------------------------------------
CubitBox OCCSurface::bounding_box() const 
{
  TopoDS_Face face = *myTopoDSFace;
  BRepAdaptor_Surface asurface(face);
  Bnd_Box aBox;
  BndLib_AddSurface::Add(asurface, Precision::Approximation(), aBox);
  double min[3], max[3];
  aBox.Get( min[0], min[1], min[2], max[0], max[1], max[2]);
  return CubitBox(min, max);
}


CubitStatus OCCSurface::get_point_normal( CubitVector& location,
                                            CubitVector& normal )
{
  assert(0);
  return closest_point( location, NULL, &normal );
}   

CubitStatus OCCSurface::closest_point_uv_guess(  
          CubitVector const& location,
          double& , double& ,
          CubitVector* closest_location,
          CubitVector* unit_normal )
{
  // don't use u and v guesses
 return closest_point(location, closest_location, unit_normal);
}


//-------------------------------------------------------------------------
// Purpose       : Computes the closest_point on the surface to the input 
//                 location.  Optionally, it also computes and returns
//                 the normal to the surface and the principal curvatures
//                 at closest_location.
//
//-------------------------------------------------------------------------
CubitStatus OCCSurface::closest_point( CubitVector const& location, 
                                         CubitVector* closest_location,
                                         CubitVector* unit_normal_ptr,
                                         CubitVector* curvature1_ptr,
                                         CubitVector* curvature2_ptr)
{
  BRepAdaptor_Surface asurface(*myTopoDSFace);
  gp_Pnt p(location.x(), location.y(), location.z()), newP(0.0, 0.0, 0.0);
  double minDist=0.0, u, v;
  int i;
  BRepLProp_SLProps SLP(asurface, 2, Precision::PConfusion());
  Extrema_ExtPS ext(p, asurface, Precision::Approximation(), Precision::Approximation());
  if (ext.IsDone() && (ext.NbExt() > 0)) {
	  for ( i = 1 ; i < ext.NbExt() ; i++ ) {
		  if ( (i==1) || (p.Distance(ext.Point(i).Value()) < minDist) ) {
			  minDist = p.Distance(ext.Point(i).Value());
			  newP = ext.Point(i).Value();
			  ext.Point(i).Parameter(u, v);
			  SLP.SetParameters(u, v);
		  }
	  }
  }
  *closest_location = CubitVector(newP.X(), newP.Y(), newP.Z());
  if (unit_normal_ptr != NULL) {
	  gp_Dir normal;
	  if (SLP.IsNormalDefined()) {
		  normal = SLP.Normal();
		  *unit_normal_ptr = CubitVector(normal.X(), normal.Y(), normal.Z()); 
	  }
  }
  
  return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Purpose       : Computes the closest_point on the trimmed surface to the 
//                 input location. 
//
// Special Notes : 
//-------------------------------------------------------------------------
void OCCSurface::closest_point_trimmed( CubitVector from_point, 
                                         CubitVector& point_on_surface)
{
  BRepAdaptor_Surface asurface(*myTopoDSFace);
  gp_Pnt p(from_point.x(), from_point.y(), from_point.z()), newP(0.0, 0.0, 0.0);
  double minDist=0.0;
  int i;
  Extrema_ExtPS ext(p, asurface, Precision::Approximation(), Precision::Approximation());
  if (ext.IsDone() && (ext.NbExt() > 0)) {
	  for ( i = 1 ; i < ext.NbExt() ; i++ ) {
		  if ( (i==1) || (p.Distance(ext.Point(i).Value()) < minDist) ) {
			  minDist = p.Distance(ext.Point(i).Value());
			  newP = ext.Point(i).Value();
		  }
	  }
  }
  point_on_surface = CubitVector(newP.X(), newP.Y(), newP.Z());
}

//-------------------------------------------------------------------------
// Purpose       : This functions computes the point on the surface that is 
//                 closest to the input location and then calculates the 
//                 magnitudes of the principal curvatures at this (possibly, 
//                 new) point on the surface. Specifying the RefVolume for 
//                 reference is optional.
//
// Special Notes :
//
//-------------------------------------------------------------------------

CubitStatus OCCSurface::principal_curvatures(
  CubitVector const& location, 
  double& curvature_1,
  double& curvature_2,
  CubitVector* closest_location )
{
  BRepAdaptor_Surface asurface(*myTopoDSFace);
  gp_Pnt p(location.x(), location.y(), location.z()), newP(0.0, 0.0, 0.0);
  double minDist=0.0, u, v;
  int i;
  BRepLProp_SLProps SLP(asurface, 2, Precision::PConfusion());
  Extrema_ExtPS ext(p, asurface, Precision::Approximation(), Precision::Approximation());
  if (ext.IsDone() && (ext.NbExt() > 0)) {
	  for ( i = 1 ; i < ext.NbExt() ; i++ ) {
		  if ( (i==1) || (p.Distance(ext.Point(i).Value()) < minDist) ) {
			  minDist = p.Distance(ext.Point(i).Value());
			  newP = ext.Point(i).Value();
			  ext.Point(i).Parameter(u, v);
			  SLP.SetParameters(u, v);
		  }
	  }
  }
  *closest_location = CubitVector(newP.X(), newP.Y(), newP.Z());
  curvature_1 = SLP.MinCurvature();
  curvature_2 = SLP.MaxCurvature();
  return CUBIT_SUCCESS;
}


//-------------------------------------------------------------------------
// Purpose       : Given values of the two parameters, get the position.
//
// Special Notes :
//
//-------------------------------------------------------------------------
CubitVector OCCSurface::position_from_u_v (double u, double v)
{
  BRepAdaptor_Surface asurface(*myTopoDSFace);
  gp_Pnt p = asurface.Value(u, v);
  return CubitVector (p.X(), p.Y(), p.Z());
}

//-------------------------------------------------------------------------
// Purpose       : This function returns the {u, v} coordinates of the point 
//                 on the Surface closest to the input point (specified in 
//                 global space). The closest_location is also returned.
//
// Special Notes :
//
//-------------------------------------------------------------------------
CubitStatus OCCSurface::u_v_from_position( CubitVector const& location,
                                             double& u,
                                             double& v,
                                             CubitVector* closest_location )
{
  BRepAdaptor_Surface asurface(*myTopoDSFace);
  gp_Pnt p(location.x(), location.y(), location.z()), newP(0.0, 0.0, 0.0);
  double minDist=0.0;
  int i;
  Extrema_ExtPS ext(p, asurface, Precision::Approximation(), Precision::Approximation());
  if (ext.IsDone() && (ext.NbExt() > 0)) {
	  for ( i = 1 ; i <= ext.NbExt() ; i++ ) {
		  if ( (i==1) || (p.Distance(ext.Point(i).Value()) < minDist) ) {
			  minDist = p.Distance(ext.Point(i).Value());
			  newP = ext.Point(i).Value();
			  ext.Point(i).Parameter(u, v);
		  }
	  }
  }
  if (closest_location != NULL) *closest_location = CubitVector(newP.X(), newP.Y(), newP.Z());
  return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Purpose       : Determines whether the Facet surface is periodic. Not
//                 available yet.
//                 
//
//-------------------------------------------------------------------------
CubitBoolean OCCSurface::is_periodic() 
{
  BRepAdaptor_Surface asurface(*myTopoDSFace);
  return (asurface.IsUPeriodic() || asurface.IsVPeriodic())?CUBIT_TRUE:CUBIT_FALSE;
}

//-------------------------------------------------------------------------
// Purpose       : Determines if the face is periodic in the given parameter
//                 direction.  Not available yet.
//
//-------------------------------------------------------------------------
CubitBoolean OCCSurface::is_periodic_in_U( double& period ) 
{
  BRepAdaptor_Surface asurface(*myTopoDSFace);
  return (asurface.IsUPeriodic())?CUBIT_TRUE:CUBIT_FALSE;
}

//-------------------------------------------------------------------------
// Purpose       : Determines if the face is periodic in the given parameter
//                 direction.  Not available yet.
//
//
//-------------------------------------------------------------------------
CubitBoolean OCCSurface::is_periodic_in_V( double& period ) 
{
  BRepAdaptor_Surface asurface(*myTopoDSFace);
  return (asurface.IsVPeriodic())?CUBIT_TRUE:CUBIT_FALSE;
}

//-------------------------------------------------------------------------
// Purpose       : Determines if the face is singular in the given parameter
//                 direction.  Not available yet.
//
//-------------------------------------------------------------------------
CubitBoolean OCCSurface::is_singular_in_U( double )
{
  assert(0);

  //PRINT_ERROR("OCCSurface::is_singular_in_U not implemented yet\n");
  return CUBIT_FALSE;
}  

//-------------------------------------------------------------------------
// Purpose       : Determines if the face is singular in the given parameter
//                 direction.  Not available yet.
//-------------------------------------------------------------------------
CubitBoolean OCCSurface::is_singular_in_V( double )
{
  assert(0);

  //PRINT_ERROR("OCCSurface::is_singular_in_V not implemented yet\n");
  return CUBIT_FALSE;
}

//-------------------------------------------------------------------------
// Purpose       : Determines if the face is closed in the U parameter.
//
//-------------------------------------------------------------------------
CubitBoolean OCCSurface::is_closed_in_U()
{
  BRepAdaptor_Surface asurface(*myTopoDSFace);
  return (asurface.IsUClosed())?CUBIT_TRUE:CUBIT_FALSE;
}

//-------------------------------------------------------------------------
// Purpose       : Determines if the face is closed in the V parameter.
//-------------------------------------------------------------------------
CubitBoolean OCCSurface::is_closed_in_V()
{
  BRepAdaptor_Surface asurface(*myTopoDSFace);
  return (asurface.IsUClosed())?CUBIT_TRUE:CUBIT_FALSE;
}

//-------------------------------------------------------------------------
// Purpose       : Calculates the derivitives at a given parameter location.
//
//-------------------------------------------------------------------------
CubitStatus OCCSurface::uv_derivitives( double u,
                                          double v,
                                          CubitVector &du,
                                          CubitVector &dv )
{
  BRepAdaptor_Surface asurface(*myTopoDSFace);
  gp_Pnt p;
  gp_Vec d1u, d1v;
  asurface.D1(u, v, p, d1u, d1v);
  du = CubitVector(d1u.X(), d1u.Y(), d1u.Z());
  dv = CubitVector(d1v.X(), d1v.Y(), d1v.Z());
  return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Purpose       : Determines whether the surface is parametrically defined.
//                 Hopefully later this will be available.
//
//-------------------------------------------------------------------------
CubitBoolean OCCSurface::is_parametric() 
{
  return CUBIT_TRUE;
}

//-------------------------------------------------------------------------
// Purpose       : Returns the lower and upper parametric bounds of the 
//                 surface in U, if it is parametric.  Otherwise, it returns
//                 CUBIT_FALSE and zeroes for the upper and lower parametric
//                 bounds.
//
// Creator       : 
//
// Creation Date : 
//-------------------------------------------------------------------------
CubitBoolean OCCSurface::get_param_range_U( double& lower_bound,
                                             double& upper_bound )
{
  BRepAdaptor_Surface asurface(*myTopoDSFace);
  lower_bound = asurface.FirstUParameter();
  upper_bound = asurface.LastUParameter();
  return CUBIT_TRUE;
}

//-------------------------------------------------------------------------
// Purpose       : Returns the lower and upper parametric bounds of the 
//                 surface in V, if it is parametric.  Otherwise, it returns
//                 CUBIT_FALSE and zeroes for the upper and lower parametric
//                 bounds.
//
//-------------------------------------------------------------------------
CubitBoolean OCCSurface::get_param_range_V( double& lower_bound,
                                             double& upper_bound )
{
  BRepAdaptor_Surface asurface(*myTopoDSFace);
  lower_bound = asurface.FirstVParameter();
  upper_bound = asurface.LastVParameter();
  return CUBIT_TRUE;
}

//-------------------------------------------------------------------------
// Purpose       : Returns the area of the Surface
//
//-------------------------------------------------------------------------
double OCCSurface::measure() 
{
  GProp_GProps myProps;
  BRepGProp::SurfaceProperties(*myTopoDSFace, myProps);
  return myProps.Mass();
}

//-------------------------------------------------------------------------
// Purpose       : This function tests the passed in position to see if
//                 is on the underlying surface.
//
//-------------------------------------------------------------------------
CubitBoolean OCCSurface::is_position_on( CubitVector &test_position )
{
  CubitVector new_point;
  CubitStatus stat = closest_point(test_position, &new_point, NULL,NULL,NULL);
  if ( !stat )
    return CUBIT_FALSE;
  CubitVector result_vec = test_position - new_point;
  if ( result_vec.length_squared() < GEOMETRY_RESABS )
    return CUBIT_TRUE;
  return CUBIT_FALSE;
}

CubitPointContainment OCCSurface::point_containment( const CubitVector &/*point*/ )
{
  assert(0);
   return CUBIT_PNT_UNKNOWN;
}

CubitPointContainment OCCSurface::point_containment( double /*u_param*/, 
                                                       double /*v_param*/ )
{
  assert(0);
  return CUBIT_PNT_UNKNOWN; 
}

//CubitPointContainment OCCSurface::point_containment( const CubitVector &/*point*/, 
//                                                       double /*u_param*/,
//                                                       double /*v_param*/ )
//{
//   return CUBIT_PNT_UNKNOWN;
//}

CubitSense OCCSurface::get_geometry_sense()
{
  TopAbs_Orientation d = myTopoDSFace->Orientation();
  if (d == TopAbs_FORWARD) return CUBIT_FORWARD;
  else if (d == TopAbs_REVERSED) return CUBIT_REVERSED;
  else {
          printf("Check Orientation (surface)");
          return CUBIT_UNKNOWN;
  }
}

void OCCSurface::get_parents_virt( DLIList<TopologyBridge*>& parents )
  { /*CAST_LIST_TO_PARENT( myShells, parents );*/ }

void OCCSurface::get_children_virt( DLIList<TopologyBridge*>& children )
{
  TopTools_IndexedMapOfShape M;
  TopExp::MapShapes(*myTopoDSFace, TopAbs_WIRE, M);
  int ii;
  for (ii=1; ii<=M.Extent(); ii++) {
	  TopologyBridge *loop = OCCQueryEngine::occ_to_cgm(M(ii));
	  children.append_unique(loop);
  }
}

// return the sense with respect to the given shell
CubitSense OCCSurface::get_shell_sense( ShellSM* shell_ptr ) const
{
  return CUBIT_FORWARD;
}


void OCCSurface::get_loops( DLIList<OCCLoop*>& result_list )
{
  TopTools_IndexedMapOfShape M;
  TopExp::MapShapes(*myTopoDSFace, TopAbs_WIRE, M);
  int ii;
  for (ii=1; ii<=M.Extent(); ii++) {
	  TopologyBridge *loop = OCCQueryEngine::occ_to_cgm(M(ii));
	  result_list.append_unique(dynamic_cast<OCCLoop*>(loop));
  }
}

void OCCSurface::get_coedges( DLIList<OCCCoEdge*>& result_list )
{
  DLIList<OCCLoop*> loop_list;
  get_loops( loop_list );
  loop_list.reset();
  for ( int i = 0; i < loop_list.size(); i++ )
    result_list += loop_list.next(i)->coedges( );
}

void OCCSurface::get_curves( DLIList<OCCCurve*>& result_list )
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