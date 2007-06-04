//-------------------------------------------------------------------------
// Filename      : SurfaceACIS.cpp
//
// Purpose       : 
//
// Special Notes :
//
// Creator       : Malcolm J. Panthaki
//
// Creation Date : 08/02/96
//
// Owner         : Malcolm J. Panthaki
//-------------------------------------------------------------------------

// ********** BEGIN STANDARD INCLUDES      **********
// ********** END STANDARD INCLUDES        **********

// ********** BEGIN ACIS INCLUDES          **********
#if CUBIT_ACIS_VERSION < 1100
#include "kernel/acis.hxx"
#include "kernel/kernapi/api/kernapi.hxx"
#include "kernel/kerndata/data/datamsc.hxx"
#include "constrct/kernapi/api/cstrapi.hxx"
#include "intersct/kernapi/api/intrapi.hxx"
#include "kernel/kernapi/api/api.hxx"
#include "kernel/kerndata/top/face.hxx"
#include "kernel/kerndata/top/alltop.hxx"
#include "kernel/kerndata/geom/allsurf.hxx"
#include "kernel/kerndata/lists/lists.hxx"
#include "kernel/kerngeom/surface/surdef.hxx"
#include "kernel/sg_husk/query/q_wire.hxx"
#include "kernel/sg_husk/face/faceutil.hxx"
#include "kernel/kernint/d3_chk/chk_stat.hxx"
#include "intersct/kerndata/ptinface/ptinface.hxx"
#include "baseutil/vector/position.hxx"
#include "baseutil/vector/param.hxx"
#include "baseutil/vector/vector.hxx"
#include "baseutil/vector/unitvec.hxx"
#include "baseutil/vector/interval.hxx"

#else
#include "acis.hxx"
#include "kernapi.hxx"
#include "datamsc.hxx"
#include "cstrapi.hxx"
#include "intrapi.hxx"
#include "api.hxx"
#include "face.hxx"
#include "alltop.hxx"
#include "allsurf.hxx"
#include "lists.hxx"
#include "surdef.hxx"
#include "q_wire.hxx"
#include "chk_stat.hxx"
#include "ptinface.hxx"
#include "position.hxx"
#include "param.hxx"
#include "vector.hxx"
#include "unitvec.hxx"
#include "interval.hxx"
#include "getbox.hxx"
#include "faceutil.hxx"
#include "raytest.hxx"
#endif
// ********** END ACIS INCLUDES            **********

// ********** BEGIN CUBIT INCLUDES         **********
#include "attrib_cubit_owner.hpp"
#include "attrib_snl_simple.hpp"

#include "CubitVector.hpp"
#include "GeometryDefines.h"
#include "SurfaceACIS.hpp"
#include "RefFace.hpp"

#include "ShellSM.hpp"
#include "LoopSM.hpp"

#include "CubitUtil.hpp"
#include "CastTo.hpp"
#include "AcisQueryEngine.hpp"
#include "RefVolume.hpp"
#include "DLIList.hpp"
// ********** END CUBIT INCLUDES           **********

// ********** BEGIN STATIC DECLARATIONS    **********
// ********** END STATIC DECLARATIONS      **********

// ********** BEGIN PUBLIC FUNCTIONS       **********
//-------------------------------------------------------------------------
// Purpose       : The constructor with a pointer to the first FACE. 
//
// Special Notes :
//
// Creator       : Xuechen Liu
//
// Creation Date : 08/02/96
//-------------------------------------------------------------------------
SurfaceACIS::SurfaceACIS(FACE* FACE_ptr)
    : AcisBridge(FACE_ptr)
{
    // Calculate a bounding box if there isn't one already
  get_acis_query_engine()->bounding_box(FACE_ptr);
}

//-------------------------------------------------------------------------
// Purpose       : The default destructor. 
//
// Special Notes :
//
// Creator       : Raikanta Sahu
//
// Creation Date : 09/06/96
//-------------------------------------------------------------------------
SurfaceACIS::~SurfaceACIS() 
{
}

FACE* SurfaceACIS::get_FACE_ptr() const
{
  return (FACE*)ENTITY_ptr();
}

void SurfaceACIS::set_FACE_ptr(FACE* FACE_ptr)
{
  ENTITY_ptr(FACE_ptr);
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
void SurfaceACIS::append_simple_attribute_virt(CubitSimpleAttrib* csattrib_ptr)
{
  AcisBridge::append_simple_attribute_virt(csattrib_ptr);
}

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
void SurfaceACIS::remove_simple_attribute_virt(CubitSimpleAttrib* csattrib_ptr)
{
  AcisBridge::remove_simple_attribute_virt(csattrib_ptr);
}

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
void SurfaceACIS::remove_all_simple_attribute_virt()
{
  AcisBridge::remove_all_simple_attribute_virt();
}

//-------------------------------------------------------------------------
// Purpose       : The purpose of this function is to get the  
//                 attributes attached to this geometry entity. The name is 
//                 attached to the underlying BODY this points to.
//
// Special Notes : 
//
// Creator       : Malcolm J. Panthaki
//
// Creation Date : 01/23/97
//-------------------------------------------------------------------------
CubitStatus SurfaceACIS::get_simple_attribute(DLIList<CubitSimpleAttrib*>&
                                           cubit_simple_attrib_list)
{
  return AcisBridge::get_simple_attribute(cubit_simple_attrib_list);
}
CubitStatus SurfaceACIS::get_simple_attribute(const CubitString& name,
                                         DLIList<CubitSimpleAttrib*>& list )
  { return AcisBridge::get_simple_attribute(name,list); }

//-------------------------------------------------------------------------
// Purpose       : Get geometry modeling engine: AcisQueryEngine
//
// Special Notes :
//
// Creator       : jihong Ma
//
// Creation Date : 10/22/96
//-------------------------------------------------------------------------

GeometryQueryEngine* 
                 SurfaceACIS::get_geometry_query_engine() const
{
   return get_acis_query_engine();   
}                 

//-------------------------------------------------------------------------
// Purpose       : Get the bounding box of the object.
//
// Special Notes :
//
// Creator       : Raikanta Sahu
//
// Creation Date : 10/23/96
//-------------------------------------------------------------------------
CubitBox SurfaceACIS::bounding_box() const 
{
    // Calculate a bounding box if there isn't one already.
  SPAbox ACIS_box = get_acis_query_engine()->bounding_box(get_FACE_ptr());
  
    // Convert to a CubitBox and return it
  return get_acis_query_engine()->bounding_box(ACIS_box);
}

CubitStatus SurfaceACIS::get_point_normal( CubitVector& point, CubitVector& normal )
{
   if( geometry_type() != PLANE_SURFACE_TYPE )
      return CUBIT_FAILURE;
   
   FACE* FACE_ptr = get_FACE_ptr();
   
   surface const* acis_surface = &(FACE_ptr->geometry()->equation());
   
   //plane const* plane_surface = (plane *)acis_surface;
   
   SPAposition pos_surf_orig;
   SPAunit_vector uv_surf_norm;
   acis_surface->planar( pos_surf_orig, uv_surf_norm );
   
   point.set( pos_surf_orig.x(), pos_surf_orig.y(), pos_surf_orig.z() );
   normal.set( uv_surf_norm.x(), uv_surf_norm.y(), uv_surf_norm.z() );

   return CUBIT_SUCCESS;
}   


//-------------------------------------------------------------------------
// Creator       : Will Frost
//
// Creation Date : 07/15/2003
//-------------------------------------------------------------------------
CubitStatus SurfaceACIS::closest_point_uv_guess(  
    CubitVector const& location, 
    double& u, double& v,
    CubitVector* closest_location,
    CubitVector* unit_normal_ptr )
{


    // If there are no underlying ACIS FACE, yell bloody murder
    // and die a horrible death! :-)
  if (ENTITY_ptr() == NULL)
  {
    PRINT_ERROR("In SurfaceACIS::closest_point\n"
                "       No underlying ACIS FACEs.\n");
    assert (ENTITY_ptr() != NULL);
    return CUBIT_FAILURE;
  }
  
    // Get the FACE to be used for computing the closest location, normal, etc.
  FACE* FACE_ptr = NULL;   
  
    // Use the first FACE in the FACE list.
  FACE_ptr = get_FACE_ptr();
    // If no valid FACE was found, crash and burn :)
  if (FACE_ptr == NULL) 
  {
    PRINT_ERROR("In SurfaceACIS::closest_point\n"
                "  No first FACE associated with a SurfaceACIS.\n"
                "  THIS IS A BUG - PLEASE REPORT IT!\n");
    return CUBIT_FAILURE;
  }
  
    // Now compute the point on the FACE closest to the input location,
    // and any other information that is required (normal and/or curvatures)
  SPAposition point ( location.x(), location.y(), location.z() );
  SPAposition closest_position;
  SPAunit_vector normal;
  
  SPApar_pos uv_guess(u, v);
  SPApar_pos uv_actual(0.0, 0.0);

  if(unit_normal_ptr != NULL)
  {
    (FACE_ptr->geometry()->equation()).point_perp ( point, 
                                                    closest_position,
                                                    normal,
                                                    uv_guess,
                                                    uv_actual );
    if ( get_FACE_sense() == CUBIT_REVERSED )
    {
      normal = -normal;
    }
    unit_normal_ptr->set( normal.x(), normal.y(), normal.z() );
  }
  else
  {
    (FACE_ptr->geometry()->equation()).point_perp ( point, 
                                                    closest_position, 
                                                    uv_guess,
                                                    uv_actual );
  }

  u = uv_actual.u;
  v = uv_actual.v;

    // Fill in the closest_location object if it was passed in
  if (closest_location)
    closest_location->set( closest_position.x(), 
                           closest_position.y(), 
                           closest_position.z() );

  return CUBIT_SUCCESS;
}




//-------------------------------------------------------------------------
// Purpose       : Computes the closest_point on the surface to the input 
//                 location.  Optionally, it also computes and returns
//                 the normal to the surface and the principal curvatures
//                 at closest_location.
//
// Special Notes : The querying is done on the *first* FACE.
//
//                 If the normal and/or the principal curvatures are 
//                 needed by the calling code, it must allocate space
//                 for these CubitVectors and pass the relevant non_NULL
//                 pointers in.  These are optional.
//
// Creator       : Malcolm J. Panthaki
//
// Creation Date : 12/11/96
//-------------------------------------------------------------------------
CubitStatus SurfaceACIS::closest_point( CubitVector const& location, 
                                        CubitVector* closest_location,
                                        CubitVector* unit_normal_ptr,
                                        CubitVector* curvature1_ptr,
                                        CubitVector* curvature2_ptr)
{
    // If there are no underlying ACIS FACE, yell bloody murder
    // and die a horrible death! :-)
  if (ENTITY_ptr() == NULL)
  {
    PRINT_ERROR("In SurfaceACIS::closest_point\n"
                "       No underlying ACIS FACEs.\n");
    assert (ENTITY_ptr() != NULL);
    return CUBIT_FAILURE;
  }
  
    // Get the FACE to be used for computing the closest location, normal, etc.
  FACE* FACE_ptr = NULL;   
  
    // Use the first FACE in the FACE list.
  FACE_ptr = get_FACE_ptr();
    // If no valid FACE was found, crash and burn :)
  if (FACE_ptr == NULL) 
  {
    PRINT_ERROR("In SurfaceACIS::closest_point\n"
                "  No first FACE associated with a SurfaceACIS.\n"
                "  THIS IS A BUG - PLEASE REPORT IT!\n");
    return CUBIT_FAILURE;
  }
  
    // Now compute the point on the FACE closest to the input location,
    // and any other information that is required (normal and/or curvatures)
  SPAposition point ( location.x(), location.y(), location.z() );
  SPAposition closest_position;
  SPAunit_vector normal;
  surf_princurv principal_curvatures;
  
    // Only need to compute the closest location
  if (unit_normal_ptr == NULL && curvature1_ptr == NULL && curvature2_ptr == NULL)
  {
    (FACE_ptr->geometry()->equation()).point_perp ( point, 
                                                   closest_position );
  }
  
    // Need to compute the closest location and the normal, but not the
    // curvatures
  else if ( (unit_normal_ptr != NULL) && 
            (curvature1_ptr == NULL && curvature2_ptr == NULL) )
  {
    (FACE_ptr->geometry()->equation()).point_perp ( 
      point, closest_position, normal );
  }
  
    // The greedy caller wants *everything*! We aim to please, of course :-)
  else
  {
    (FACE_ptr->geometry()->equation()).point_perp ( 
      point, 
      closest_position, 
      normal,
      principal_curvatures);
  }   
  
  // Fill in the closest_location object if it was passed in
  if(closest_location)
    closest_location->set( closest_position.x(), 
                           closest_position.y(), 
                           closest_position.z() ); 
  
    // Fill in the normal, if necessary.
  if (unit_normal_ptr != NULL)
  {
      // Reverse the normal, if the sense of the FACE is reversed wrt the 
      // underlying ACIS surface object.
    if ( get_FACE_sense() == CUBIT_REVERSED )
    {
      normal = -normal;
    }
    unit_normal_ptr->set( normal.x(), normal.y(), normal.z() );
  }
  
    // Fill in the first principal curvature, if necessary
  if (curvature1_ptr != NULL)
  {
      // First set curvature1 using the direction of this curvature
    curvature1_ptr->set( principal_curvatures.dir1.x(),
                         principal_curvatures.dir1.y(), 
                         principal_curvatures.dir1.z() );
    if( !(fabs(principal_curvatures.dir1.x()) < CUBIT_RESABS &&
          fabs(principal_curvatures.dir1.y()) < CUBIT_RESABS && 
          fabs(principal_curvatures.dir1.z()) < CUBIT_RESABS ) )
    {
        // Now multiply it with the magnitude of the curvature
      *(curvature1_ptr) = (principal_curvatures.curv1) * 
        (*(curvature1_ptr));
        // Reverse the normal, if the sense of the FACE is reversed wrt the 
        // underlying ACIS surface object.
    }
    if ( get_FACE_sense() == CUBIT_REVERSED )
    {
      *curvature1_ptr = (*curvature1_ptr) * -1.0;
    }
    
  }
  
    // Fill in the second principal curvature, if necessary
  if (curvature2_ptr != NULL)
  {
      // First set curvature2 using the direction of this curvature
    curvature2_ptr->set( principal_curvatures.dir2.x(),
                         principal_curvatures.dir2.y(), 
                         principal_curvatures.dir2.z() );
    
      // Now multiply it with the magnitude of the curvature
    if( !(fabs(principal_curvatures.dir2.x()) < CUBIT_RESABS &&
          fabs(principal_curvatures.dir2.y()) < CUBIT_RESABS && 
          fabs(principal_curvatures.dir2.z()) < CUBIT_RESABS ) )
    {
        // Now multiply it with the magnitude of the curvature
      *(curvature2_ptr) = (principal_curvatures.curv2) * 
        (*(curvature2_ptr));
        // Reverse the normal, if the sense of the FACE is reversed wrt the 
        // underlying ACIS surface object.
    }
    if ( get_FACE_sense() == CUBIT_REVERSED )
    {
      *(curvature2_ptr) = (*curvature2_ptr) * -1.0;
    }
  }
  
  return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Purpose       : Computes the closest_point on the trimmed surface to the 
//                 input location. 
//
// Special Notes : 
// Creator       : Brett W. Clark
//
// Creation Date : 3/26/97
//-------------------------------------------------------------------------
void SurfaceACIS::closest_point_trimmed( CubitVector from_point, 
                                         CubitVector& point_on_surface)
{
    // If there are no underlying ACIS FACEs, yell bloody murder
    // and die a horrible death! :-)
   if (ENTITY_ptr() == NULL)
   {
      PRINT_ERROR("In SurfaceACIS::closest_point_trimmed\n"
                  "       No underlying ACIS FACEs.\n");
      assert (ENTITY_ptr() != NULL);
   }
   
     // Get the FACE to be used for computing the closest location.
   FACE* FACE_ptr = get_FACE_ptr();
   const SPAposition ACIS_from_point (from_point.x(), from_point.y(), 
                                   from_point.z() );
   SPAposition closest_point;
   outcome result =
     api_find_cls_ptto_face(ACIS_from_point, FACE_ptr, closest_point );
   if ( !result.ok() )
   {
     PRINT_ERROR("Problems with ACIS closest point trimmed.\n");
   }
   point_on_surface.x(closest_point.x());
   point_on_surface.y(closest_point.y());
   point_on_surface.z(closest_point.z());
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
// Creator       : Malcolm J. Panthaki
//
// Creation Date : 03/04/97
//-------------------------------------------------------------------------
CubitStatus SurfaceACIS::principal_curvatures(
  CubitVector const& location, 
  double& curvature_1,
  double& curvature_2,
  CubitVector* closest_location )
{
    // Get the principal curvature vectors
  CubitVector curvature1_vector;
  CubitVector curvature2_vector;
  CubitStatus result = this->closest_point( location,
                                            closest_location,
                                            (CubitVector *)NULL,
                                            &curvature1_vector,
                                            &curvature2_vector );
  
  if ( result == CUBIT_FAILURE )
  {
    PRINT_ERROR("In SurfaceACIS::principal_curvatures\n"
                "       Could not compute the principal curvature vectors.\n");
    return CUBIT_FAILURE;
  }
  
    // Extract the magnitudes of the principal curvatures
  curvature_1 = curvature1_vector.length();
  curvature_2 = curvature2_vector.length();
  
  return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Purpose       : Given values of the two parameters, get the position.
//
// Special Notes :
//
// Creator       : Raikanta Sahu
//
// Creation Date : 01/05/97
//-------------------------------------------------------------------------
CubitVector SurfaceACIS::position_from_u_v (double u, double v)
{
    // Get the first FACE
  FACE* FACE_ptr = get_FACE_ptr();
  
    // Construct a parametric position
  SPApar_pos node_par_pos (u, v);
  
    // Use the parametric position to get the SPAposition of the 
    // parametric point
  SPAposition node_pos;
  
  node_pos = (FACE_ptr->geometry()->equation()).
    eval_position ( node_par_pos );
  
  return CubitVector(node_pos.x(), node_pos.y(), node_pos.z());
}

//-------------------------------------------------------------------------
// Purpose       : This function returns the {u, v} coordinates of the point 
//                 on the Surface closest to the input point (specified in 
//                 global space). The closest_location is also returned.
//
// Special Notes :
//
// Creator       : Malcolm J. Panthaki
//
// Creation Date : 03/04/97
//-------------------------------------------------------------------------
CubitStatus SurfaceACIS::u_v_from_position (CubitVector const& location,
                                            double& u, 
                                            double& v,
                                            CubitVector* closest_location )
{
    // A RefVolume was provided for reference and the calling code needs 
    // at least one of the following to be generated: normal or the
    // principal curvatures (i.e., at least one of those pointers is
    // non-NULL)
  FACE* FACE_ptr = get_FACE_ptr();
    // If no valid FACE was found, crash and burn :)
  if (FACE_ptr == NULL)
  {
    PRINT_ERROR("In SurfaceACIS::u_v_from_position\n"
                "  THIS IS A BUG - PLEASE REPORT IT!\n");
    return CUBIT_FAILURE;
  }
  
    // Get the parameter values of the closest point, based on the FACE 
    // just extracted, above
  surface const *surface_ptr = &(FACE_ptr->geometry()->equation());
  
  if ( closest_location )
  {
      // Get the closest point
    CubitStatus result = closest_point( location,
                                        closest_location );
    if ( result == CUBIT_FAILURE )
    {
      PRINT_ERROR("In SurfaceACIS::u_v_from_position\n"
                  "       Could not compute the closest point to the input "
                  "point {%f %f %f}.\n",
                  closest_location->x(), 
                  closest_location->y(), 
                  closest_location->z());
      assert ( result == CUBIT_SUCCESS );
      return CUBIT_FAILURE;
    }
    SPAposition point( closest_location->x(), 
                    closest_location->y(), 
                    closest_location->z() );
    SPApar_pos params;
    params = surface_ptr->param( point );
    u = params.u;
    v = params.v;
  }
  else
  {
    SPAposition point( location.x(), 
                    location.y(), 
                    location.z() );
    SPApar_pos params;
    
    params = surface_ptr->param( point );
    u = params.u;
    v = params.v;
  }
  return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Purpose       : Determines whether the ACIS surface object associated
//                 with one of the FACEs of this SurfaceACIS object is 
//                 periodic or not.
//
// Special Notes : MJP Note:
//                 This code is very ACIS-specific and could change with 
//                 new versions of ACIS.
//
//                 The querying is done on the *first* FACE.
//
// Creator       : Malcolm J. Panthaki
//
// Creation Date : 12/12/96
//-------------------------------------------------------------------------
CubitBoolean SurfaceACIS::is_periodic() 
{
    // Get the number of LOOPs that bound a FACE of this object.
  int num_loops = number_of_LOOPs();
  
  if ( num_loops == 2 &&
       ( this->geometry_type() == CONE_SURFACE_TYPE   ||
         this->geometry_type() == SPHERE_SURFACE_TYPE ||
         this->geometry_type() == TORUS_SURFACE_TYPE ) )
  {
    return CUBIT_TRUE;
  }
  
  else
  {
    return CUBIT_FALSE;
  }
}

//-------------------------------------------------------------------------
// Purpose       : Determines whether the ACIS surface object associated
//                 with one of the FACEs of this SurfaceACIS object is 
//                 periodic in the U direction or not.  If it is, it
//                 returns CUBIT_TRUE and the value of the period. Otherwise,
//                 it returns CUBIT_FALSE and a value of 0.0 or the period.
//
// Special Notes : MJP Note:
//                 This code is very ACIS-specific and could change with 
//                 new versions of ACIS.
//
//                 The querying is done on the *first* FACE.
//
// Creator       : Malcolm J. Panthaki
//
// Creation Date : 12/17/96
//-------------------------------------------------------------------------
CubitBoolean SurfaceACIS::is_periodic_in_U( double& period_U ) 
{
    // Get the first FACE in the list of FACEs associated with this
    // SurfaceACIS object
  FACE* FACE_ptr = get_FACE_ptr();
  
  surface const *surface_ptr = &(FACE_ptr->geometry()->equation());
  
    // Now ask ACIS whether its underlying surface is periodic in the
    // U direction
  logical is_periodic;
  is_periodic = surface_ptr->periodic_u();
  
  if (is_periodic)
  {
      // Now get the value of the period and return it
    period_U = surface_ptr->param_period_u();
    return CUBIT_TRUE;
  }
  else
  {
    period_U = 0.0;
    return CUBIT_FALSE;
  }
}

//-------------------------------------------------------------------------
// Purpose       : Determines whether the ACIS surface object associated
//                 with one of the FACEs of this SurfaceACIS object is 
//                 periodic in the V direction or not. If it is, it
//                 returns CUBIT_TRUE and the value of the period. Otherwise,
//                 it returns CUBIT_FALSE and a value of 0.0 or the period.
//
// Special Notes : MJP Note:
//                 This code is very ACIS-specific and could change with 
//                 new versions of ACIS.
//
//                 The querying is done on the *first* FACE.
//
// Creator       : Malcolm J. Panthaki
//
// Creation Date : 12/17/96
//-------------------------------------------------------------------------
CubitBoolean SurfaceACIS::is_periodic_in_V( double& period_V ) 
{
    // Get the first FACE in the list of FACEs associated with this
    // SurfaceACIS object
  FACE* FACE_ptr = get_FACE_ptr();
  
  surface const *surface_ptr = &(FACE_ptr->geometry()->equation());
    // Now ask ACIS whether its underlying surface is periodic in the
    // V direction
  int is_periodic;
  is_periodic = surface_ptr->periodic_v();
  
  if (is_periodic)
  {
      // Now get the value of the period and return it
    period_V = surface_ptr->param_period_v();
    return CUBIT_TRUE;
  }
  else
  {
    period_V = 0.0;
    return CUBIT_FALSE;
  }
}

//-------------------------------------------------------------------------
// Purpose       : Determines whether the ACIS surface object associated
//                 with one of the FACEs of this SurfaceACIS object is 
//                 singular in the U direction or not at a given u parameter.
//
// Note          : The assumption is made that the u_param is in the
//                 bounds of the surface.
//
// Creator       : David White
//
// Creation Date : 10/08/97
//-------------------------------------------------------------------------
CubitBoolean SurfaceACIS::is_singular_in_U( double u_param )
{
    // Get the first FACE in the list of FACEs associated with this
    // SurfaceACIS object
  FACE* FACE_ptr = get_FACE_ptr();
  
  double u_lower, u_upper;
  get_param_range_U( u_lower, u_upper );
  if ( u_param < u_lower - CUBIT_RESABS ||
       u_param > u_upper + CUBIT_RESABS )
  {
    PRINT_ERROR("v parameter is outside parameter bounds.\n");
    return CUBIT_FALSE;
  }
  
    // Now ask ACIS whether its underlying surface is periodic in the
    // V direction
  surface const *surface_ptr = &(FACE_ptr->geometry()->equation());
  
  int is_singular = surface_ptr->singular_u( u_param );
  
  if (is_singular)
  {
    return CUBIT_TRUE;
  }
  
  else
  {
    return CUBIT_FALSE;
  }
}  

//-------------------------------------------------------------------------
// Purpose       : Determines whether the ACIS surface object associated
//                 with one of the FACEs of this SurfaceACIS object is 
//                 singular in the V direction or not at a given v parameter.
//
// Creator       : David White
//
// Creation Date : 10/08/97
//-------------------------------------------------------------------------
CubitBoolean SurfaceACIS::is_singular_in_V( double v_param )
{
    // Get the first FACE in the list of FACEs associated with this
    // SurfaceACIS object
  FACE* FACE_ptr = get_FACE_ptr();
  double v_lower, v_upper;
  get_param_range_V( v_lower, v_upper );
  if ( v_param < v_lower - CUBIT_RESABS ||
       v_param > v_upper + CUBIT_RESABS )
  {
    PRINT_ERROR("v parameter is outside parameter bounds.\n");
    return CUBIT_FALSE;
  }
  surface const *surface_ptr = &(FACE_ptr->geometry()->equation());
  
    // Now ask ACIS whether its underlying surface is periodic in the
    // V direction
  int is_singular = surface_ptr->singular_v( v_param );
  
  if (is_singular)
  {
    return CUBIT_TRUE;
  }
  else
  {
    return CUBIT_FALSE;
  }
}

//-------------------------------------------------------------------------
// Purpose       : Determines whether the ACIS surface object associated
//                 with one of the FACEs of this SurfaceACIS object is 
//                 closed in the U direction or not.
//
// Creator       : David White
//
// Creation Date : 10/08/97
//-------------------------------------------------------------------------
CubitBoolean SurfaceACIS::is_closed_in_U()
{
    // Get the first FACE in the list of FACEs associated with this
    // SurfaceACIS object
  FACE* FACE_ptr = get_FACE_ptr();
  
    // Now ask ACIS whether its underlying surface is periodic in the
    // U direction
  int is_closed = (&(FACE_ptr->geometry()->equation()))->closed_u();
  
  if (is_closed)
  {
    return CUBIT_TRUE;
  }
  
  else
  {
    return CUBIT_FALSE;
  }
}

//-------------------------------------------------------------------------
// Purpose       : Determines whether the ACIS surface object associated
//                 with one of the FACEs of this SurfaceACIS object is 
//                 closed in the V direction or not.
//
// Creator       : David White
//
// Creation Date : 10/08/97
//-------------------------------------------------------------------------
CubitBoolean SurfaceACIS::is_closed_in_V()
{
    // Get the first FACE in the list of FACEs associated with this
    // SurfaceACIS object
  FACE* FACE_ptr = get_FACE_ptr();
  
    // Now ask ACIS whether its underlying surface is periodic in the
    // V direction
  surface const *surface_ptr = &(FACE_ptr->geometry()->equation());
  
  int is_closed = surface_ptr->closed_v();
  
  if (is_closed)
  {
    return CUBIT_TRUE;
  }
  
  else
  {
    return CUBIT_FALSE;
  }
}

//-------------------------------------------------------------------------
// Purpose       : Calculates the derivitives at a given parameter location.
//
// Creator       : David White
//
// Creation Date : 10/08/97
//-------------------------------------------------------------------------
CubitStatus SurfaceACIS::uv_derivitives( double u_param,
                                         double v_param,
                                         CubitVector &du,
                                         CubitVector &dv )
{
    //Get the first face and work of it...
  FACE* FACE_ptr = get_FACE_ptr();
  
  surface const *surface_ptr = &(FACE_ptr->geometry()->equation());
  
  SPApar_pos params( u_param, v_param );
  SPAposition xyz_position;
  SPAvector first_derives[2];
  surface_ptr->eval( params, xyz_position, first_derives );
  
  du.set( first_derives[0].x(), first_derives[0].y(),
          first_derives[0].z() );
  dv.set( first_derives[1].x(), first_derives[1].y(),
          first_derives[1].z() );
  return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Purpose       : Determines whether the ACIS surface object associated
//                 with one of the FACEs of this SurfaceACIS object is 
//                 parametrically defined or not.
//
// Special Notes : MJP Note:
//                 This code is very ACIS-specific and could change with 
//                 new versions of ACIS.
//
//                 The querying is done on the *first* FACE.
//
// Creator       : Malcolm J. Panthaki
//
// Creation Date : 12/17/96
//-------------------------------------------------------------------------
CubitBoolean SurfaceACIS::is_parametric() 
{
  return CUBIT_TRUE;
}

//-------------------------------------------------------------------------
// Purpose       : Returns the lower and upper parametric bounds of the 
//                 surface in U, if it is parametric.  Otherwise, it returns
//                 CUBIT_FALSE and zeroes for the upper and lower parametric
//                 bounds.
//
// Special Notes : MJP Note:
//                 This code is very ACIS-specific and could change with 
//                 new versions of ACIS.
//
//                 The querying is done on the *first* FACE.
//
// Creator       : Malcolm J. Panthaki
//
// Creation Date : 12/17/96
//-------------------------------------------------------------------------
CubitBoolean SurfaceACIS::get_param_range_U( double& lower_bound,
                                             double& upper_bound )
{
    // Get the first FACE in the list of FACEs associated with this
    // SurfaceACIS object
  FACE* FACE_ptr = get_FACE_ptr();
  
    // Now get the parametric range in U from ACIS, if this is a parametric
    // surface
  CubitBoolean alls_well = CUBIT_TRUE;

//  SPApar_box face_range;

//   API_BEGIN
//     sg_get_face_par_box( FACE_ptr, face_range );   
//   API_END

//     // Make sure it is bounded and return the appropriate values
//   if ( face_range.u_range().bounded() )
//   {
//     lower_bound = face_range.u_range().start_pt();
//     upper_bound = face_range.u_range().end_pt();
//   }
  SPAbox *bbox = FACE_ptr->bound();
  SPAinterval u_face_range;
  u_face_range = FACE_ptr->geometry()->equation().param_range_u(*bbox);
// Make sure it is bounded and return the appropriate values
  if ( u_face_range.bounded() )
  {
    lower_bound = u_face_range.start_pt();
    upper_bound = u_face_range.end_pt();
  }
  
    // The parameter range is unbounded in one or both directions.
  else
  {
    alls_well = CUBIT_FALSE;
  }
  
    // The surface is not parametric or the range is unbounded
  if (alls_well == CUBIT_FALSE)
  {
    lower_bound = 0.0;
    upper_bound = 0.0;
    return CUBIT_FALSE;
  }
  
  return CUBIT_TRUE;
}

//-------------------------------------------------------------------------
// Purpose       : Returns the lower and upper parametric bounds of the 
//                 surface in V, if it is parametric.  Otherwise, it returns
//                 CUBIT_FALSE and zeroes for the upper and lower parametric
//                 bounds.
//
// Special Notes : MJP Note:
//                 This code is very ACIS-specific and could change with 
//                 new versions of ACIS.
//
//                 The querying is done on the *first* FACE.
//
// Creator       : Malcolm J. Panthaki
//
// Creation Date : 12/17/96
//-------------------------------------------------------------------------
CubitBoolean SurfaceACIS::get_param_range_V( double& lower_bound,
                                             double& upper_bound )
{
    // Get the first FACE in the list of FACEs associated with this
    // SurfaceACIS object
  FACE* FACE_ptr = get_FACE_ptr();
  
    // Now get the parametric range in V from ACIS, if this is a parametric
    // surface
  CubitBoolean alls_well = CUBIT_TRUE;

//   SPApar_box face_range;

//   API_BEGIN
//     sg_get_face_par_box( FACE_ptr, face_range );   
//   API_END

//     // Make sure it is bounded and return the appropriate values
//   if ( face_range.v_range().bounded() )
//   {
//     lower_bound = face_range.v_range().start_pt();
//     upper_bound = face_range.v_range().end_pt();
//   }
  SPAbox *bbox = FACE_ptr->bound();
  SPAinterval v_face_range;
  v_face_range = FACE_ptr->geometry()->equation().param_range_v(*bbox);
// Make sure it is bounded and return the appropriate values
  if ( v_face_range.bounded() )
  {
    lower_bound = v_face_range.start_pt();
    upper_bound = v_face_range.end_pt();
  }
  
  else
  {
    alls_well = CUBIT_FALSE;
  }
  
    // The surface is not parametric or the range is unbounded
  if (alls_well == CUBIT_FALSE)
  {
    lower_bound = 0.0;
    upper_bound = 0.0;
    return CUBIT_FALSE;
  }
  
  return CUBIT_TRUE;
}

//-------------------------------------------------------------------------
// Purpose       : Returns a surface type ID -- the values of these are
//                 determined by ACIS.
//
// Special Notes : MJP Note:
//                 This code is very ACIS-specific and could change with 
//                 new versions of ACIS.  There are #defines for the
//                 various surface type ID's.  These are defined in the 
//                 header files of each of the specific surface classes
//                 in ACIS.
//
//                 The querying is done on the *first* FACE.
//
// Creator       : Malcolm J. Panthaki
//
// Creation Date : 12/12/96
//-------------------------------------------------------------------------
GeometryType SurfaceACIS::geometry_type()
{
    // Get the first FACE associated with this SurfaceACIS object
  FACE* FACE_ptr = get_FACE_ptr();
  
    // Ask ACIS for the type of the underlying surface of that FACE
  int surface_type = (&(FACE_ptr->geometry()->equation()))->type();
  
  GeometryType local_type;
  
  switch (surface_type)
  {
    case cone_type:
      local_type = CONE_SURFACE_TYPE;
      break;
    case plane_type:
      local_type = PLANE_SURFACE_TYPE;
      break;
    case sphere_type:
      local_type = SPHERE_SURFACE_TYPE;
      break;
    case spline_type:
      local_type = SPLINE_SURFACE_TYPE;
      break;
    case torus_type:
      local_type = TORUS_SURFACE_TYPE;
      break;
    default:
      local_type = UNDEFINED_SURFACE_TYPE;
      break;
  }
  
  return local_type;
}

//-------------------------------------------------------------------------
// Purpose       : Returns the area of the Surface
//
// Special Notes :
//
// Creator       : Raikanta Sahu
//
// Creation Date : 12/19/96
//-------------------------------------------------------------------------
double SurfaceACIS::measure() 
{
    // Get a FACE from the list of FACES
  FACE* FACE_ptr = get_FACE_ptr();
  
  double result_accuracy = 0.0;
  double area = 0.0;
  outcome result = api_ent_area ( FACE_ptr,
                                  1.0e-4,
                                  area,
                                  result_accuracy );
  if (!result.ok())
  {
    PRINT_ERROR("In SurfaceACIS::area.\n"
                "       ACIS api_ent_area function failed\n");
    get_acis_query_engine()->ACIS_API_error (result);
    return 0.0;
  }
  else
    return area;
}

//-------------------------------------------------------------------------
// Purpose       : This function tests the passed in position to see if
//                 is on the underlying surface.
//
// Special Notes :
//
// Creator       : David R. White
//
// Creation Date : 04/08/97
//-------------------------------------------------------------------------
CubitBoolean SurfaceACIS::is_position_on( CubitVector &test_position )
{
  SPAposition pos ( test_position.x(), test_position.y(), test_position.z() );

  surface const *acis_surface;
  acis_surface = &(get_FACE_ptr()->geometry()->equation());
  
  if ( acis_surface->test_point( pos ) )
  {
    return CUBIT_TRUE;
  }
  else
  {
    return CUBIT_FALSE;
  }
}

//-------------------------------------------------------------------------
// Purpose       : Determine whether a point, known to lie on the surface 
//                 lies INSIDE, OUTSIDE, or on the BOUNDARY curves of the 
//                 surface. 
//                  
//
// Special Notes :
//
// Creator       : Steve Storm
//
// Creation Date : 12/01/00
//-------------------------------------------------------------------------
CubitPointContainment SurfaceACIS::point_containment( const CubitVector &point )
{
   FACE *FACE_ptr = get_FACE_ptr();
   
   SPAposition test_point( point.x(), point.y(), point.z() );
   
   SPAtransf ftrans;

   //SPApar_pos params;
   //params = acis_surface->param( test_point );
   //point_face_containment pf_rel = point_in_face( test_point, FACE_ptr, ftrans,
   //                                               params);

   // Let the function figure-out uv coordinates
   point_face_containment pf_rel;
#if CUBIT_ACIS_VERSION < 1100
   pf_rel = point_in_face( test_point, FACE_ptr, ftrans );
#else
   outcome r = api_point_in_face( test_point, FACE_ptr, ftrans, pf_rel );
   if (!r.ok())
   {
     PRINT_ERROR("ACIS api_point_in_face failed.\n");
     return CUBIT_PNT_UNKNOWN;
   }
#endif
   
   switch( pf_rel )
   {
   case point_outside_face:
      return CUBIT_PNT_OUTSIDE;
   case point_inside_face:
      return CUBIT_PNT_INSIDE;
   case point_boundary_face:
      return CUBIT_PNT_BOUNDARY;
   case point_unknown_face:
      return CUBIT_PNT_UNKNOWN;
   }
   return CUBIT_PNT_UNKNOWN;
}

CubitPointContainment SurfaceACIS::point_containment( double u_param, double v_param )
{
   FACE *FACE_ptr = get_FACE_ptr();
   surface const *acis_surface = &(FACE_ptr->geometry()->equation());
   
   SPApar_pos test_uv( u_param, v_param );
   
   SPAposition test_point = acis_surface->eval_position(test_uv);
   
   SPAtransf ftrans;
   
   point_face_containment pf_rel;
   
#if CUBIT_ACIS_VERSION < 1100
   pf_rel = point_in_face( test_point, FACE_ptr, ftrans, test_uv);
#else
   outcome r = api_point_in_face( test_point, FACE_ptr, ftrans, pf_rel, test_uv );
   if (!r.ok())
   {
     PRINT_ERROR("ACIS api_point_in_face failed.\n");
     return CUBIT_PNT_UNKNOWN;
   }
#endif
   
   switch( pf_rel )
   {
   case point_outside_face:
      return CUBIT_PNT_OUTSIDE;
   case point_inside_face:
      return CUBIT_PNT_INSIDE;
   case point_boundary_face:
      return CUBIT_PNT_BOUNDARY;
   case point_unknown_face:
      return CUBIT_PNT_UNKNOWN;
   }
   return CUBIT_PNT_UNKNOWN;
}

CubitPointContainment SurfaceACIS::point_containment( const CubitVector &point, 
                                                      double u_param, 
                                                      double v_param )
{
   FACE *FACE_ptr = get_FACE_ptr();
   
   SPAposition test_point( point.x(), point.y(), point.z() );
   
   SPApar_pos test_uv( u_param, v_param );
   
   SPAtransf ftrans;
   
   point_face_containment pf_rel;
#if CUBIT_ACIS_VERSION < 1100
   pf_rel = point_in_face( test_point, FACE_ptr, ftrans, test_uv);
#else
   outcome r = api_point_in_face( test_point, FACE_ptr, ftrans, pf_rel, test_uv );
   if (!r.ok())
   {
     PRINT_ERROR("ACIS api_point_in_face failed.\n");
     return CUBIT_PNT_UNKNOWN;
   }
#endif

   switch( pf_rel )
   {
   case point_outside_face:
      return CUBIT_PNT_OUTSIDE;
   case point_inside_face:
      return CUBIT_PNT_INSIDE;
   case point_boundary_face:
      return CUBIT_PNT_BOUNDARY;
   case point_unknown_face:
      return CUBIT_PNT_UNKNOWN;
   }
   return CUBIT_PNT_UNKNOWN;
}

int SurfaceACIS::validate(const CubitString &user_name,
                          DLIList <TopologyEntity*> &bad_entities)
{
  int error = check_FACE((FACE*)ENTITY_ptr(), user_name);
  if ( error > 0 )
  {
    bad_entities.append(this->topology_entity());
  }
    
  return error;
  
}

int SurfaceACIS::check_FACE(FACE *face,
                            const CubitString &refentity_name)
{
  int error = 0;
  assert(face != NULL);
  
  check_status_list *list = NULL;
  api_check_face(face, list);
  while(list != NULL){
    error++;
    switch(list->status()){
      case check_irregular:
        PRINT_ERROR("\ttwisted or scrunched up surface for %s\n",
                    refentity_name.c_str());
        break;
      case check_self_intersects:
        PRINT_ERROR("\tself-intersecting surface for %s\n",
                    refentity_name.c_str());
        break;
      case check_bad_closure:
        PRINT_ERROR("\tsurface closure is wrong for %s\n",
                    refentity_name.c_str());
        break;
      case check_bs3_null:
        PRINT_ERROR("\tno bs3 surface for %s\n",
                    refentity_name.c_str());
        break;
      case check_bs3_coi_verts:
        PRINT_ERROR("\terror in control point coincidence on surface for %s\n",
                    refentity_name.c_str());
        break;
      case check_bad_degeneracies:
        PRINT_ERROR("\tdegenerate edges on %s\n",
                    refentity_name.c_str());
        break;
      case check_untreatable_singularity:
        PRINT_ERROR("\tuntreatable singularities in surface for %s\n",
                    refentity_name.c_str());
        break;
      case check_non_G0:
        PRINT_ERROR("\tsurface for %s is not G0\n",
                    refentity_name.c_str());
        break;
      case check_non_G1:
        PRINT_ERROR("\tsurface for %s is not G1\n",
                    refentity_name.c_str());
        break;
      case check_non_G2:
        PRINT_ERROR("\tsurface for %s is not G2\n",
                    refentity_name.c_str());
        break;
      case check_non_C1:
        PRINT_ERROR("\tsurface for %s is not C1\n",
                    refentity_name.c_str());
        break;
      case check_unknown:
        break;
      default:
        break;
    }
    list = list->next();
  }
  return error;
}

CubitSense SurfaceACIS::get_geometry_sense()
{
  CubitSense sense = get_FACE_sense();
  if (this->get_FACE_ptr()->geometry()->equation().left_handed_uv())
  {
    if (sense == CUBIT_FORWARD)
      sense = CUBIT_REVERSED;
    else
      sense = CUBIT_FORWARD;
  }
  return sense;
}


void SurfaceACIS::get_parents_virt( DLIList<TopologyBridge*>& parents )
{
  ENTITY_LIST entities;
  api_get_shells( get_FACE_ptr(), entities );
  ATTRIB_CUBIT_OWNER::cubit_owner( entities, parents );
}

void SurfaceACIS::get_children_virt( DLIList<TopologyBridge*>& children )
{
  ENTITY_LIST entities;
  api_get_loops( get_FACE_ptr(), entities );
  ATTRIB_CUBIT_OWNER::cubit_owner( entities, children );
}


//-------------------------------------------------------------------------
// Purpose       : Get CoFace sense
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 08/27/02
//-------------------------------------------------------------------------
CubitSense SurfaceACIS::get_shell_sense( ShellSM* ) const
{
  FACE* face = get_FACE_ptr();
  if( !face || face->sides() == DOUBLE_SIDED )
    return CUBIT_UNKNOWN;
  
  return CUBIT_FORWARD;
}  

// ********** END PUBLIC FUNCTIONS         **********

// ********** BEGIN PROTECTED FUNCTIONS    **********
// ********** END PROTECTED FUNCTIONS      **********

// ********** BEGIN PRIVATE FUNCTIONS      **********
//-------------------------------------------------------------------------
// Purpose       : Return the sense of the first FACE (in the list of FACEs
//                 in this Surface) wrt its underlying ACIS surface.
//
// Special Notes : 
//
// Creator       : Malcolm J. Panthaki
//
// Creation Date : 04/03/97
//-------------------------------------------------------------------------
CubitSense SurfaceACIS::get_FACE_sense()
{
    // Get the first FACE
  FACE* FACEPtr = this->get_FACE_ptr();
  
    // Return the sense value
  REVBIT FACEsense = FACEPtr->sense();
  if (FACEsense == FORWARD)
  {
    return CUBIT_FORWARD;
  }
  else if (FACEsense == REVERSED)
  {
    return CUBIT_REVERSED;
  }
  else
  {
    PRINT_ERROR("In SurfaceACIS::get_FACE_sense\n"
                "       Could not get the sense of the FACE "
                "wrt its surface.\n");
    return CUBIT_FORWARD;
  }
}

int SurfaceACIS::number_of_LOOPs()
{
  FACE *FACE_ptr = get_FACE_ptr();
  
  int number_of_LOOPs = 0;
  LOOP* LOOP_ptr = FACE_ptr->loop();
  while ( LOOP_ptr != NULL )
  {
    number_of_LOOPs++;
    LOOP_ptr = LOOP_ptr->next();
  }
  
  return number_of_LOOPs;
}

CubitStatus SurfaceACIS::fire_ray(const CubitVector &ray_point,
                                  const CubitVector &unit,
                                  DLIList<double>& ray_params) const
{
   ray_params.clean_out();
     
     // fire a ray at the specified surface, returning the hits and
     // the parameters along the ray; return non-zero if error
   
   SPAposition pos(ray_point[0], ray_point[1], ray_point[2]);
   
  SPAunit_vector vector(unit[0], unit[1], unit[2]);
  hit *temp_hit;
  ray this_ray(pos, vector, GEOMETRY_RESABS, 2);
  hit *hit_list = raytest_face(this_ray, get_FACE_ptr());
  hit *orig_hit_list = hit_list;
  if ( hit_list != NULL ){
    while (hit_list != NULL) {
      ray_params.append(hit_list->ray_param);
      temp_hit = hit_list->next;
      delete hit_list;
      hit_list = temp_hit;
    }
  }

  return CUBIT_SUCCESS;
}

// ********** END PRIVATE FUNCTIONS        **********

// ********** BEGIN HELPER CLASSES         **********
// ********** END HELPER CLASSES           **********

// ********** BEGIN EXTERN FUNCTIONS       **********
// ********** END EXTERN FUNCTIONS         **********

// ********** BEGIN STATIC FUNCTIONS       **********
// ********** END STATIC FUNCTIONS         **********
