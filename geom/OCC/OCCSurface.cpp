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
#include "CubitEvaluator.hpp"
#include "CubitUtil.hpp"
#include "CastTo.hpp"
#include "RefVolume.hpp"
#include "GeometryQueryEngine.hpp"
#include "DLIList.hpp"
#include "ShellSM.hpp"
#include "Lump.hpp"
#include "LoopSM.hpp"
#include "CubitPointData.hpp"

////// #include "FacetEvalTool.hpp"
////// #include "CubitFacetData.hpp"
////// #include "CubitFacetEdge.hpp"

// ********** END CUBIT INCLUDES           **********

//// #include "CubitEvaluator.hpp"
//// #include "SphereEvaluator.hpp"
//// #include "CylinderEvaluator.hpp"

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

//// OCCSurface::OCCSurface(FacetEvalTool *facet_tool,
////                           DLIList<ShellSM*> &shellsms,
////                           DLIList<LoopSM*> &loopsms)
//// {
////  assert(0);
////    // Calculate a bounding box if there isn't one already
////  facetEvalTool = facet_tool;
////    //sense_ = CUBIT_FORWARD;
////  myShells += shellsms;
////  myLoops += loopsms;
////  myShellSense = CUBIT_UNKNOWN;
////  myEvaluator = NULL;
//// }

//// //-------------------------------------------------------------------------
//// // Purpose       : The constructor with a pointer to the FacetEvalTool. 
//// //
//// // Special Notes : Used for save/restore
//// //
//// //-------------------------------------------------------------------------
//// OCCSurface::OCCSurface( const SphereEvaluatorData *sphere_data,
////                             FacetEvalTool *facet_tool,
////                             DLIList<ShellSM*> &shellsms,
////                             DLIList<LoopSM*> &loopsms )
//// {
////   assert(0);
////   facetEvalTool = facet_tool;
////     //sense_ = CUBIT_FORWARD;
////   myShells += shellsms;
////   myLoops += loopsms;
////   myShellSense = CUBIT_UNKNOWN;
//// 
////   myEvaluator = new SphereEvaluator( sphere_data );
//// }
//// //-------------------------------------------------------------------------
//// // Purpose       : The constructor with a pointer to the FacetEvalTool. 
//// //
//// // Special Notes : Used for save/restore
//// //
//// //-------------------------------------------------------------------------
//// OCCSurface::OCCSurface( const CylinderEvaluatorData *cylinder_data,
////                             FacetEvalTool *facet_tool,
////                             DLIList<ShellSM*> &shellsms,
////                             DLIList<LoopSM*> &loopsms )
//// {
////   assert(0);
////   facetEvalTool = facet_tool;
////     //sense_ = CUBIT_FORWARD;
////   myShells += shellsms;
////   myLoops += loopsms;
////   myShellSense = CUBIT_UNKNOWN;
//// 
////   myEvaluator = new CylinderEvaluator( cylinder_data );
//// }

//// //-------------------------------------------------------------------------
//// // Purpose       : The constructor with a pointer to the FacetEvalTool. 
//// //
//// // Special Notes : Used for save/restore
//// //
//// //-------------------------------------------------------------------------
//// OCCSurface::OCCSurface(FacetEvalTool *facet_tool,
////                            CubitSense sense,
////                            CubitSense shell_sense0,
////                            CubitBoolean use_facets,
////                            DLIList<LoopSM*> &loopsms)
//// {
////   assert(0);
////     // Calculate a bounding box if there isn't one already
////   facetEvalTool = facet_tool;
////     //sense_ = CUBIT_FORWARD;
////   myLoops += loopsms;
////   myShellSense = shell_sense0;
////   myEvaluator = NULL;
//// }
//// //-------------------------------------------------------------------------
//// // Purpose       : The default destructor. 
//// //
//// // Special Notes :
//// //
//// //-------------------------------------------------------------------------

OCCSurface::~OCCSurface() 
{
  ////  if ( facetEvalTool )
  //// {
  ////  delete facetEvalTool;
  //// }
}


//-------------------------------------------------------------------------
// Purpose       : get the interpolation order of the FacetEvalTool 
//
// Special Notes :
//
//-------------------------------------------------------------------------
//// int OCCSurface::interp_order() 
//// {
////   assert(0);
////   assert(facetEvalTool != NULL);
////  
////  return facetEvalTool->interp_order();
//// }

//-------------------------------------------------------------------------
// Purpose       : get the min dot of the FacetEvalTool 
//
// Special Notes :
//
//-------------------------------------------------------------------------
//// double OCCSurface::min_dot() 
//// {
////  assert(0);
////  assert(facetEvalTool != NULL);
////  
////  return facetEvalTool->get_min_dot();
////}


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
  if ( myEvaluator )
      return myEvaluator->is_singular_in_U();

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
  if ( myEvaluator )
      return myEvaluator->is_singular_in_V();

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
// Purpose       : Returns a surface type ID
//
//-------------------------------------------------------------------------
GeometryType OCCSurface::geometry_type()
{
  return UNDEFINED_SURFACE_TYPE;
/*    if ( is_flat() )
    {
        return PLANE_SURFACE_TYPE;
    }
    else if ( is_spherical() )
    {
        return SPHERE_SURFACE_TYPE;
    }
    else if ( is_conical() )
    {
        return CONE_SURFACE_TYPE;
    }
    else
    {
        return UNDEFINED_SURFACE_TYPE;
    }*/
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

//// //-------------------------------------------------------------------------
//// // Purpose       : Updates the (cached) area of the Surface so that if 
//// //                 measure() is called, the correct value is returned.
//// //-------------------------------------------------------------------------
//// void OCCSurface::update_measurement() 
//// {
////  assert(0);
////  facetEvalTool->calculate_area();
//// }

//// //-------------------------------------------------------------------------
//// // Purpose       : Returns whether the facet surface is completely flat
//// //
//// //-------------------------------------------------------------------------
//// CubitBoolean OCCSurface::is_flat() 
//// {
////   assert(0);
////   // Danilov: try to use BRepAdaptor_Surface::GetType
////   return (facetEvalTool->is_flat() == 1) ? CUBIT_TRUE : CUBIT_FALSE; 
//// }

//// //-------------------------------------------------------------------------
//// // Purpose       : Returns whether the facet surface is spherical
//// //
//// //-------------------------------------------------------------------------
//// CubitBoolean OCCSurface::is_spherical()
//// {
////   assert(0);
////   // Danilov: try to use BRepAdaptor_Surface::GetType
////     if ( myEvaluator && myEvaluator->ask_type() == SPHERE_SURFACE_TYPE )
////     {
////         return CUBIT_TRUE;
////     }
////     return CUBIT_FALSE;
//// }
//// 

//// //-------------------------------------------------------------------------
//// // Purpose       : Returns whether the facet surface is spherical
//// //
//// //-------------------------------------------------------------------------
//// CubitBoolean OCCSurface::is_conical()
//// {
////   assert(0);
////   // Danilov: try to use BRepAdaptor_Surface::GetType
////     if ( myEvaluator && myEvaluator->ask_type() == CONE_SURFACE_TYPE )
////     {
////         return CUBIT_TRUE;
////     }
////     return CUBIT_FALSE;
//// }

//// const CubitEvaluatorData *OCCSurface::evaluator_data( void )
//// {
////   assert(0);
////     if ( myEvaluator )
////     {
////         return myEvaluator->evaluator_data();
////     }
////     return NULL;
//// }

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
  assert(0);
  CubitSense sense = CUBIT_FORWARD; //get_relative_surface_sense();
  return sense;
}

// CubitSense OCCSurface::get_relative_surface_sense()
// {
//     //not sure if this is right for the facet surface...
//   return sense_;
// }

/*

void OCCSurface::bodysms(DLIList<BodySM*> &bodies) 
{
  int ii;
  for ( ii = myShells.size(); ii > 0; ii-- )
  {
    myShells.get_and_step()->bodysms(bodies);
  }
}

void OCCSurface::lumps(DLIList<Lump*> &lumps)
{
  int ii;
  for ( ii = myShells.size(); ii > 0; ii-- )
  {
    myShells.get_and_step()->lumps(lumps);
  }
}

void OCCSurface::shellsms(DLIList<ShellSM*> &shellsms)
{
  int ii;
  for ( ii = myShells.size(); ii > 0; ii-- )
  {
    shellsms.append_unique(myShells.get_and_step());
  }
}

void OCCSurface::surfaces(DLIList<Surface*> &surfaces)
{
  surfaces.append_unique(this);
}

void OCCSurface::loopsms(DLIList<LoopSM*> &loopsms)
{
  int ii;
  for ( ii = myLoops.size(); ii > 0; ii-- )
  {
    loopsms.append_unique(myLoops.get_and_step());
  }
}
void OCCSurface::coedgesms(DLIList<CoEdgeSM*> &coedgesms)
{
  int ii;
  for ( ii = myLoops.size(); ii > 0; ii-- )
  {
    myLoops.get_and_step()->coedgesms(coedgesms);
  }
}

void OCCSurface::curves(DLIList<Curve*> &curves)
{
  int ii;
  for ( ii = myLoops.size(); ii > 0; ii-- )
  {
    myLoops.get_and_step()->curves(curves);
  }
}
void OCCSurface::points(DLIList<Point*> &points)
{
  int ii;
  for ( ii = myLoops.size(); ii > 0; ii-- )
  {
    myLoops.get_and_step()->points(points);
  }
}
*/


void OCCSurface::get_parents_virt( DLIList<TopologyBridge*>& parents )
  { CAST_LIST_TO_PARENT( myShells, parents ); }

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


//// CubitStatus OCCSurface::get_my_facets(DLIList<CubitFacet*> &facet_list,
////                                         DLIList<CubitPoint*>& point_list)
//// {
////   facetEvalTool->get_facets(facet_list);
////   facetEvalTool->get_points(point_list);
////   return CUBIT_SUCCESS;
//// }

//// void OCCSurface::tris(DLIList<CubitFacet *> &facet_list)
//// {
//// 	facetEvalTool->get_facets(facet_list);
//// }

//// void OCCSurface::get_my_points(DLIList<CubitPoint*> &point_list)
//// {
////   facetEvalTool->get_points(point_list);
//// }

//// void OCCSurface::get_my_facetedges(DLIList<CubitFacetEdge*> &edge_list)
//// {
////  facetEvalTool->get_edges(edge_list);
//// }

void OCCSurface::get_shell_sense( CubitSense &sense0)
{
  TopAbs_Orientation d = myTopoDSFace->Orientation();
  if (d == TopAbs_FORWARD) sense0 = CUBIT_FORWARD;
  else if (d == TopAbs_REVERSED) sense0 = CUBIT_REVERSED;
  else {
	  printf("Check Orientation (surface)");
	  sense0 = CUBIT_UNKNOWN;
  }
}


// return the sense with respect to the given shell
CubitSense OCCSurface::get_shell_sense( ShellSM* shell_ptr ) const
{
  TopAbs_Orientation d = myTopoDSFace->Orientation();
  if (d == TopAbs_FORWARD) return CUBIT_FORWARD;
  else if (d == TopAbs_REVERSED) return CUBIT_REVERSED;
  else {
	  printf("Check Orientation (surface)");
	  return CUBIT_UNKNOWN;
  }
}

// set the sense of the surface with respect to the shell
void OCCSurface::set_shell_sense( OCCShell *facet_shell, 
                                    CubitSense thesense )
{
  assert(0);
//    if(thesense == CUBIT_REVERSED){
//        PRINT_INFO("should not do this.");
//    }
  int idx = myShells.get_index();
  if(idx > 0){
    PRINT_ERROR("Multiple shells attached to a single surface.\n");
    return;
  }
  ShellSM *shell_ptr = (ShellSM *)facet_shell;
  ShellSM *ashell = myShells.get();
  if (ashell == shell_ptr)
  {
    myShellSense = thesense;
  }
//   else
//   {
//     myShells.step();
//     idx = myShells.get_index();
//     ashell = myShells.get();
//     if (ashell == shell_ptr)
//     {
//       myShellSense[idx] = thesense;
//     }
//   }
}

//// //----------------------------------------------------------------
//// // Function: copy_facets
//// // Description: copy the points and facets
//// //
//// // Author: sjowen
//// //----------------------------------------------------------------
//// CubitStatus OCCSurface::copy_facets(DLIList<CubitFacet*>&copy_facet_list,
////                                       DLIList<CubitPoint*>&copy_point_list)
//// {
////   assert(0);
////   if (!facetEvalTool)
////   {
////     PRINT_ERROR("Couldn't copy facets.");
////     return CUBIT_FAILURE;
////   }
////   int ii;
////   DLIList<CubitFacet*>facet_list;
////   DLIList<CubitPoint*>point_list;
////   facetEvalTool->get_facets( facet_list );
////   facetEvalTool->get_points( point_list );
////   CubitPoint **point_array = new CubitPoint* [point_list.size()];
//// 
////   //- copy the points
//// 
////   point_list.reset();
////   CubitPoint *new_point, *the_point;
////   for(ii=0; ii<point_list.size(); ii++)
////   {
////     the_point = point_list.get_and_step();
////     new_point = new CubitPointData( the_point->coordinates() );
////     the_point->marked( ii );
////     copy_point_list.append( new_point );
////     point_array[ii] = new_point;
////   }
//// 
////   //- copy the facets
//// 
////   int jj, idx;
////   CubitFacet *new_facet, *the_facet;
////   CubitPoint *points[3];
////   for (ii=0; ii<facet_list.size(); ii++)
////   {
////     the_facet = facet_list.get_and_step();
////     for (jj=0; jj<3; jj++)
////     {
////       idx = the_facet->point(jj)->marked();
////       points[jj] = point_array[idx];
////     }
////     new_facet = new CubitFacetData( points[0], points[1], points[2] );
////     copy_facet_list.append( new_facet );
////   }
//// 
////   delete [] point_array;
//// 
////   return CUBIT_SUCCESS;
//// }

void OCCSurface::get_bodies( DLIList<OCCBody*>& result_list )
{
  DLIList<OCCLump*> lump_list;
  get_lumps( lump_list );
  lump_list.reset();
  for ( int i = lump_list.size(); i--; )
  {
    OCCLump* lump = lump_list.get_and_step();
    OCCBody* body = dynamic_cast<OCCBody*>(lump->get_body());
    if (body)
      result_list.append_unique(body);
  }
}

void OCCSurface::get_lumps( DLIList<OCCLump*>& result_list )
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

void OCCSurface::get_shells( DLIList<OCCShell*>& result_list )
{
  myShells.reset();
  for ( int i = 0; i < myShells.size(); i++ )
    if ( OCCShell* shell = dynamic_cast<OCCShell*>(myShells.next(i)) )
      result_list.append(shell);
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
    loop_list.next(i)->get_coedges( result_list );
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


//-------------------------------------------------------------------------
// Purpose       : Remove Shell from shell list
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 09/29/03
//-------------------------------------------------------------------------
CubitStatus OCCSurface::remove_shell(OCCShell* shell)
{
  assert(0);
    // Something strange here -- A DLIList of Shells and a 
    // two-element array for senses?  Try to keep the senses
    // intact anyway...
  myShells.reset();
  if (myShells.get() == shell)
    myShellSense = CUBIT_UNKNOWN;
  
  if (!myShells.move_to(shell))
    return CUBIT_FAILURE;
  
  myShells.remove();
  return CUBIT_SUCCESS;
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
void OCCSurface::disconnect_all_loops()
{
  assert(0);
  myLoops.reset();
  for (int i = myLoops.size(); i--; )
  {
    LoopSM* sm_ptr = myLoops.get_and_step();
    OCCLoop* loop = dynamic_cast<OCCLoop*>(sm_ptr);
    if (loop)
    {
      assert(loop->get_surface() == this);
      loop->remove_surface();
    }
  }
  myLoops.clean_out();
}
/*
void OCCSurface::add_transformation( CubitTransformMatrix &tfmat )
{
  assert(0);
    if ( myEvaluator )
        myEvaluator->add_transformation( tfmat );
}
*/
// ********** END PUBLIC FUNCTIONS         **********

// ********** BEGIN PROTECTED FUNCTIONS    **********
// ********** END PROTECTED FUNCTIONS      **********

// ********** BEGIN PRIVATE FUNCTIONS      **********

//// void  OCCSurface::reverse_sense()
//// {
////   assert(0);
////   facetEvalTool->reverse_facets();
////   myLoops.reset();
////   int i,j;
////   OCCLoop* this_loop;
////   LoopSM* this_loop_sm;
////   DLIList<OCCCoEdge *> this_coedge_list;
////   for(i=0;i<myLoops.size();i++){
////     this_loop_sm= myLoops.get_and_step();
////     this_loop = dynamic_cast<OCCLoop*>(this_loop_sm);
////     if(!this_loop){
////       PRINT_ERROR("Unexpected null pointer for loop.\n");
////       return;
////     }
////    this_loop->reverse();
////    this_coedge_list.clean_out();
////    this_loop->get_coedges(this_coedge_list);
////    for(j=0; j<this_coedge_list.size(); j++){
////      this_coedge_list.get_and_step()->reverse_sense();
////    }
////   }
////   
////   //sense_ = CubitUtil::opposite_sense( sense_ );
////   myShellSense = CubitUtil::opposite_sense( myShellSense );
////   //myShellSense[1] = CubitUtil::opposite_sense( myShellSense[1] );
//// }

// ********** END PRIVATE FUNCTIONS        **********

// ********** BEGIN HELPER CLASSES         **********
// ********** END HELPER CLASSES           **********

// ********** BEGIN EXTERN FUNCTIONS       **********
// ********** END EXTERN FUNCTIONS         **********

// ********** BEGIN STATIC FUNCTIONS       **********
// ********** END STATIC FUNCTIONS         **********
