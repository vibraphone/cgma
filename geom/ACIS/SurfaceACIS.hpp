//-------------------------------------------------------------------------
// Filename      : SurfaceACIS.hpp
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

#ifndef SURFACE_ACIS_HPP
#define SURFACE_ACIS_HPP

// ********** BEGIN STANDARD INCLUDES      **********
// ********** END STANDARD INCLUDES        **********

// ********** BEGIN ACIS INCLUDES          **********
#if CUBIT_ACIS_VERSION < 1100
#include "kernel/kerndata/top/face.hxx"
#include "kernel/kerndata/geom/allsurf.hxx"
#else
#include "face.hxx"
#include "allsurf.hxx"
#endif
// ********** END ACIS INCLUDES            **********

// ********** BEGIN CUBIT INCLUDES         **********
#include "CubitDefines.h"
#include "Surface.hpp"
#include "AcisBridge.hpp"
// ********** END CUBIT INCLUDES           **********

// ********** BEGIN FORWARD DECLARATIONS   **********
class TopologyEntity;
class RefVolume;
class RefFace;
class RefVolume;
template <class X> class DLIList;

class FACE;
// ********** END FORWARD DECLARATIONS     **********

class SurfaceACIS : public Surface, public AcisBridge
{
public :
  
  SurfaceACIS(FACE* FACE_ptr = NULL);
    //I- FACEPtr*
    //I- A pointer to the first FACE which the object will be associated with.
    //- This constructor takes a pointer to the first FACE to which it will be
    //- connected.
  
  virtual ~SurfaceACIS() ;
    //- The destructor
  
#ifdef BOYD14
  static CubitStatus get_FACEs_of_RefFace( RefFace *ref_face,
                                           DLIList<FACE*>& output_FACE_list);
    //R CubitStatus
    //R- CUBIT_SUCCESS/CUBIT_FAILURE
    //I ref_face
    //I- a RefFace to find the associated ACIS FACEs
    //O output_FACE_list
    //O- a list of ACIS FACE, each item is associated with the corrsponding
    //O- item of the input RefFace list.
    // Given the input RefFace, find its parent ACIS FACEs and
    // return them in a list.
#endif
  
  FACE* get_FACE_ptr() const;
#ifdef BOYD14
  static FACE* get_first_FACE(const RefFace *ref_face);
    //R FACE*
    //R- Returned FACE pointer
    //- This function returns the FACE associated with this object,
    //- or the first FACE associated with ref_face, if it is a merged surface.
#endif

  void set_FACE_ptr(FACE* FACE_ptr);
    //- I FACE_ptr - The FACE to be associated with this object.
    //- It replaces any previous FACE.
  
  virtual void append_simple_attribute_virt(CubitSimpleAttrib*);
    //R void
    //I 
    //I- 
    //I- that is to be appended to this OSME object.
    //- The purpose of this function is to append a 
    //- attribute to the OSME. The  is attached to each of the 
    //- underlying solid model entities this one points to.
  
  virtual void remove_simple_attribute_virt(CubitSimpleAttrib*);
    //R void
    //I CubitSimpleAttrib*
    //I- A reference to a CubitSimpleAttrib object which is the object
    //I- that is to be removed to this OSME object.
    //- The purpose of this function is to remove a simple
    //- attribute from the OSME. The attribute is attached to each of the
    //- underlying solid model entities this one points to.
  
  virtual void remove_all_simple_attribute_virt();
    //R void
    //I-
    //- The purpose of this function is to remove all simple
    //- attributes from the OSME. 
  
  virtual CubitStatus get_simple_attribute(DLIList<CubitSimpleAttrib*>&);
  virtual CubitStatus get_simple_attribute(const CubitString& name,
                                           DLIList<CubitSimpleAttrib*>&);
    //R CubitSimpleAttrib*
    //R- the returned cubit simple attribute.
    //- The purpose of this function is to get the attributes
    //- of the geometry entity. The name is attached to the underlying solid
    //- model entity(ies) this one points to.
    //- MJP Note:
    //- This is the code that implements the requirement that names
    //- of VGI Entities propagate across solid model boolean
    //- operations.  The success of this relies, of course, on the underlying
    //- solid modeler being able to propagate attributes across
    //- such operations on its entities. If it cannot, then "names"
    //- of VGI entities will not propagate.
  
  virtual CubitBox bounding_box() const ;
    // see comments in GeometryEntity.hpp
  
  virtual GeometryQueryEngine* 
  get_geometry_query_engine() const;
    //R GeometryQueryEngine*
    //R- A pointer to the geometric modeling engine associated with
    //R- the object.
    //- This function returns a pointer to the geometric modeling engine
    //- associated with the object.
  
  // Added by CAT
  virtual CubitStatus get_point_normal( CubitVector& point,
                                        CubitVector& normal );
    //- Only valid for planar surfaces
    //- Finds the underlying plane's origin and normal vector
    //- Returns CubitFailure if not a plane.  The origin and normal 
    //- are returned directly from the underlying format for a plane.
  
  virtual void closest_point_trimmed(CubitVector from_point,
                                     CubitVector& point_on_surface);
    //R void
    //I CubitVector
    //I- point from which to find closest point on trimmed surface
    //O CubitVector
    //O- point on trimmed surface closest to passed-in point
    //- This function finds the closest point on a TRIMMED surface to the
    //- passed-in point.

  virtual CubitStatus closest_point_uv_guess(  
    CubitVector const& location,
    double& u, double& v,
    CubitVector* closest_location = NULL,
    CubitVector* unit_normal_ptr = NULL );


  
  virtual CubitStatus closest_point(  
    CubitVector const& location, 
    CubitVector* closest_location = NULL,
    CubitVector* unit_normal_ptr = NULL,
    CubitVector* curvature1_ptr = NULL,
    CubitVector* curvature2_ptr = NULL);
    //R CubitStatus
    //R- CUBIT_SUCCESS/FAILURE
    //I location
    //I- The point to which the closest point on the surface is desired.
    //O closest_location
    //O- The point on the Surface, closest to the 
    //O- input location (which might not be on the Surface).  This is
    //O- input as a reference so that the function can modify its
    //O- contents.
    //O unit_normal_ptr
    //O- The normal (represented as a unit vector) at the closest_location.
    //O- If this pointer is NULL, the normal is not returned.
    //O curvature1_ptr
    //O- The first principal curvature of the surface at closest_location.
    //O- If this pointer is NULL, this curvature is not returned.
    //O curvature2_ptr
    //O- The second principal curvature of the surface at closest_location.
    //O- If this pointer is NULL, this curvature is not returned.
    //- This function computes the point on the surface closest to the input 
    //- location -- i.e., closest_location. 
    //- The first ACIS FACE in the list
    //- is queried.
    //-
    //- If the input pointer values of unit_normal, curvature1 and
    //- curvature2
    //- are non-NULL, the normal and principal curvatures, too, are
    //- returned.  These are computed at closest_location, not at the
    //- input location.
    //-
    //- NOTE:
    //- It is assumed that if the calling code needs the normal or the 
    //- principal curvatures, it will *allocate* space for the CubitVectors
    //- before sending in the pointers.
  
  virtual CubitStatus principal_curvatures(
    CubitVector const& location, 
    double& curvature_1,
    double& curvature_2,
    CubitVector* closest_location = NULL );
    //R CubitStatus
    //R- CUBIT_SUCCESS/FAILURE
    //I location
    //I- The point at which the curvatures are being requested -- it is also
    //I- the point to which the closest point on the surface is returned.
    //I- curvatures.
    //O closest_location
    //O- The point on the surface, closest to the input location (this
    //O- might not be on the surface).  This is input as a reference 
    //O- so that the function can modify its contents.
    //O curvature_1/2
    //O- Returned principal curvature magnitudes.
    //- This functions computes the point on the surface that is closest
    //- to the input location and then calculates the magnitudes of the
    //- principal curvatures at this (possibly, new) point on the surface. 
  
  virtual CubitVector position_from_u_v (double u, double v);
    //R CubitVector
    //R- Returned position vector.
    //I u, v
    //I- Input point in {u.v} space
    //- This function returns the coordinates in world space of a point
    //- in the parameter space of this Surface object.
  
  virtual CubitStatus u_v_from_position (CubitVector const& location,
                                         double& u, 
                                         double& v,
                                         CubitVector*
                                         closest_location = NULL );
    //R CubitStatus
    //R- CUBIT_SUCCESS/FAILURE
    //I location
    //I- The input point in global space
    //O closest_point
    //O- The point on the Surface closest to the input location
    //O u, v
    //O- The returned u, v coordinate values (in local parametric space)
    //O- of the closest_point
    //I refvolume_ptr
    //- This function returns the {u, v} coordinates of the point 
    //- on the Surface closest to the input point (specified in global
    //- space). The closest_location is also returned.
  
  virtual CubitBoolean is_periodic();
    //R CubitBoolean
    //R- CUBIT_TRUE/CUBIT_FALSE
    //- This function determines whether the underlying geometry of the
    //- SurfaceACIS is periodic or not.  Returns CUBIT_TRUE if it is and 
    //- CUBIT_FALSE if it is not.
    //- MJP NOTE: 
    //- The first ACIS FACE in the list is queried.  It is assumed
    //- that all the FACEs have the same underlying surface.
  
    virtual CubitBoolean is_periodic_in_U( double& period );
    //R CubitBoolean
    //R- CUBIT_TRUE/CUBIT_FALSE
    //O period
    //O- The value of the period in the U direction.
    //- Determines whether the ACIS surface object associated
    //- with one of the FACEs of this SurfaceACIS object is 
    //- periodic in the U direction or not.  If it is, it
    //- returns CUBIT_TRUE and the value of the period. Otherwise,
    //- it returns CUBIT_FALSE and a value of 0.0 or the period.
    //- MJP NOTE: 
    //- The first ACIS FACE in the list is queried.  It is assumed
    //- that all the FACEs have the same underlying surface.
  
  virtual CubitBoolean is_periodic_in_V( double& period );
    //R CubitBoolean
    //R- CUBIT_TRUE/CUBIT_FALSE
    //O period
    //O- The value of the period in the V direction.
    //- Determines whether the ACIS surface object associated
    //- with one of the FACEs of this SurfaceACIS object is 
    //- periodic in the V direction or not.  If it is, it
    //- returns CUBIT_TRUE and the value of the period. Otherwise,
    //- it returns CUBIT_FALSE and a value of 0.0 or the period.
    //- MJP NOTE: 
    //- The first ACIS FACE in the list is queried.  It is assumed
    //- that all the FACEs have the same underlying surface.
  
  virtual CubitBoolean is_singular_in_U( double u_param );
  virtual CubitBoolean is_singular_in_V( double v_param );
    //R CubitBoolean
    //R- CUBIT_TRUE/CUBIT_FALSE
    //I double u/v parameter value.
    //- Determines if the surface is singular in a given direction
    //- at a given parameter value.
  
  virtual CubitBoolean is_closed_in_U();  
  virtual CubitBoolean is_closed_in_V();
    //R CubitBoolean
    //R- CUBIT_TRUE/CUBIT_FALSE
    //- Determines if the surface is closed, smoothly or not in the
    //- given parameter direction.
    //- A periodic surface is always closed but a closed surface is
    //- is not always periodic.
  
  virtual CubitStatus uv_derivitives( double u_param,
                                      double v_param,
                                      CubitVector &du,
                                      CubitVector &dv );
    //R CubitStatus
    //R- CUBIT_SUCCESS/CUBIT_FAILURE
    //O- du, dv
    //- Determines the u and v derivitives from the given parameter
    //- values.
  
  virtual CubitBoolean is_parametric();
    //R CubitBoolean
    //R- CUBIT_TRUE/CUBIT_FALSE
 			//- This method returns CUBIT_TRUE if a parametric representation
			//- is available for the surface
  
  virtual CubitBoolean get_param_range_U( double& lower_bound,
                                          double& upper_bound );
    //R CubitBoolean
    //R- CUBIT_TRUE/CUBIT_FALSE
    //O lower_bound
    //O- The lower bound of the parametric range in the U direction.
    //O- This is set to 0.0 if the surface is not parametric.
    //O upper_bound
    //O- The upper bound of the parametric range in the U direction.
    //O- This is set to 0.0 if the surface is not parametric.
    //- Returns the lower and upper parametric bounds of the 
    //- surface in U, if it is parametric.  Otherwise, it returns
    //- CUBIT_FALSE and zeroes for the upper and lower parametric
    //- bounds.
    //- MJP NOTE: 
    //- The first ACIS FACE in the list is queried.  It is assumed
    //- that all the FACEs have the same underlying surface.
  
  virtual CubitBoolean get_param_range_V( double& lower_bound,
                                          double& upper_bound );
    //R CubitBoolean
    //R- CUBIT_TRUE/CUBIT_FALSE
    //O lower_bound
    //O- The lower bound of the parametric range in the V direction.
    //O- This is set to 0.0 if the surface is not parametric.
    //O upper_bound
    //O- The upper bound of the parametric range in the V direction.
    //O- This is set to 0.0 if the surface is not parametric.
    //- Returns the lower and upper parametric bounds of the 
    //- surface in V, if it is parametric.  Otherwise, it returns
    //- CUBIT_FALSE and zeroes for the upper and lower parametric
    //- bounds.
    //- MJP NOTE: 
    //- The first ACIS FACE in the list is queried.  It is assumed
    //- that all the FACEs have the same underlying surface.
  
  virtual CubitBoolean is_position_on( CubitVector &test_position );
    //R CubitBoolean
    //R- CUBIT_TRUE/CUBIT_FALSE
    //I CubitVector
    //I- position, point where we want to test, whether or not it
    //- is on the surface.

  virtual CubitPointContainment point_containment( const CubitVector &point );
  virtual CubitPointContainment point_containment( double u, double v );
  virtual CubitPointContainment point_containment( const CubitVector &point, 
                                                   double u, double v );
    //R CubitPointContainment - is the point outside, inside or on the boundary?
    //R- CUBIT_PNT_OUTSIDE, CUBIT_PNT_INSIDE, CUBIT_PNT_BOUNDARY, 
    //   CUBIT_PNT_UNKNOWN
    //I CubitVector
    //I- position to check, known to be on the Surface
    //I double
    //I- u coordinate, if known (significantly faster, if this is known - however
    //                           if not known let the function figure it out)
    //I double
    //I- v coordinate, if known (significantly faster, if this is known - however
    //                           if not known let the function figure it out)
    // NOTE: POINT MUST LIE ON THE SURFACE FOR THIS FUNCTION TO WORK PROPERLY.
  
  GeometryType geometry_type();
    //R GeometryType (enum)
    //R- The enumerated type of the geometric representation
  
  virtual double measure();
    //R double
    //R- The numeric value of the measure (its units depend on the dimension
    //R- of the RefEntity being "measured")
    //- A generic geometric extent function.
    //- Returns volume for Lump, area for Surface, length for Curve and 
    //- 1.0 for Point
  
  virtual int validate(const CubitString &user_name,
                       DLIList <TopologyEntity*> &bad_entities);
    //- Check that entity is valid. Returns number of problems detected.
  
  int check_FACE(FACE *face, const CubitString &refentity_name);
    //- perform check on an individual FACE
 
  virtual CubitSense get_geometry_sense();
    //- Return the relative surface sense. (see below)

  void get_parents_virt( DLIList<TopologyBridge*>& parents );
  void get_children_virt( DLIList<TopologyBridge*>& children );

  virtual CubitSense get_shell_sense( ShellSM* shell_ptr ) const;
  
  virtual CubitStatus fire_ray(const CubitVector &ray_point,
                               const CubitVector &unit,
                               DLIList<double>& ray_params) const;
  
protected: 
  
private:
  
#ifdef BOYD14
  surface const* get_ACIS_surface();
    //R surface const&
    //R- Returned ACIS surface object 
    //- This function returns the ACIS surface object associated with 
    //- the first FACE in FACEPtrList_
#endif
  
  CubitSense get_FACE_sense();
    //R CubitSense
    //R- Returned sense value
    //- This function returns the sense of the first ACIS FACE in FACEPtrList_
    //- wrt its underlying ACIS surface.
    //- If there is an error getting the sense value, then CUBIT_FORWARD
    //- is returned.

  int number_of_LOOPs();
    //- returns the number of acis loops in the first acis face
};


// ********** BEGIN INLINE FUNCTIONS       **********
// ********** END INLINE FUNCTIONS         **********

// ********** BEGIN FRIEND FUNCTIONS       **********
// ********** END FRIEND FUNCTIONS         **********

// ********** BEGIN EXTERN FUNCTIONS       **********
// ********** END EXTERN FUNCTIONS         **********

#endif

