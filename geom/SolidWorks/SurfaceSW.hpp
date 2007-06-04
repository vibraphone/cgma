//-------------------------------------------------------------------------
// Filename      : SurfaceSW.hpp
//
// Purpose       : 
//
// Special Notes :
//
// Creator       : Joel Kopp, Malcolm J. Panthaki
//
// Creation Date : 8/22/00
//
// Owner         : Joel Kopp, Malcolm J. Panthaki
//-------------------------------------------------------------------------

#ifndef SURFACE_SW_HPP
#define SURFACE_SW_HPP

#include "CubitDefines.h"
#include "CubitEntity.hpp"
#include "SurfaceSM.hpp"
#include "SWPart.hpp"

class TopologyEntity;
class RefVolume;
class RefFace;

class GMem;

struct IFace2;
template <class X> class DLIList;


class SurfaceSW : public SurfaceSM
{
public :
  
  SurfaceSW(SWPart *pPart);
  
  virtual ~SurfaceSW() ;
    //- The destructor
  
  IFace2 *get_FACE_ptr() const;

  void set_FACE_ptr(IFace2 *face);
    //- I faceDisp - The FACE to be associated with this object.
    //- It replaces any previous FACE.
  
  virtual GeometryQueryEngine* get_geometry_query_engine() const;
    //R GeometryQueryEngine*
    //R- A pointer to the geometric modeling engine associated with
    //R- the object.
    //- This function returns a pointer to the geometric modeling engine
    //- associated with the object.
  
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
  
  virtual CubitStatus get_simple_attribute(const CubitString& name,
                                           DLIList<CubitSimpleAttrib*>&);
  virtual CubitStatus get_simple_attribute(DLIList<CubitSimpleAttrib*>&);

  
  virtual void get_parents_virt(DLIList<TopologyBridge*> &parents );
  virtual void get_children_virt(DLIList<TopologyBridge*> &children );



  virtual CubitBox bounding_box() const ;
    // see comments in GeometryEntity.hpp
  
//  virtual GeometricModelingEngine* 
//  get_geometric_modeling_engine() const;
    //R GeometricModelingEngine*
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
    //- The first SW FACE in the list
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
    //- SurfaceSW is periodic or not.  Returns CUBIT_TRUE if it is and 
    //- CUBIT_FALSE if it is not.
    //- MJP NOTE: 
    //- The first SW FACE in the list is queried.  It is assumed
    //- that all the FACEs have the same underlying surface.
  
    virtual CubitBoolean is_periodic_in_U( double& period_U );
    //R CubitBoolean
    //R- CUBIT_TRUE/CUBIT_FALSE
    //O period
    //O- The value of the period in the U direction.
    //- Determines whether the SW surface object associated
    //- with one of the FACEs of this SurfaceSW object is 
    //- periodic in the U direction or not.  If it is, it
    //- returns CUBIT_TRUE and the value of the period. Otherwise,
    //- it returns CUBIT_FALSE and a value of 0.0 or the period.
    //- MJP NOTE: 
    //- The first SW FACE in the list is queried.  It is assumed
    //- that all the FACEs have the same underlying surface.
  
  virtual CubitBoolean is_periodic_in_V( double& period_V );
    //R CubitBoolean
    //R- CUBIT_TRUE/CUBIT_FALSE
    //O period
    //O- The value of the period in the V direction.
    //- Determines whether the SW surface object associated
    //- with one of the FACEs of this SurfaceSW object is 
    //- periodic in the V direction or not.  If it is, it
    //- returns CUBIT_TRUE and the value of the period. Otherwise,
    //- it returns CUBIT_FALSE and a value of 0.0 or the period.
    //- MJP NOTE: 
    //- The first SW FACE in the list is queried.  It is assumed
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
    //- The first SW FACE in the list is queried.  It is assumed
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
    //- The first SW FACE in the list is queried.  It is assumed
    //- that all the FACEs have the same underlying surface.
  
  virtual CubitBoolean is_position_on( CubitVector &test_position );
    //R CubitBoolean
    //R- CUBIT_TRUE/CUBIT_FALSE
    //I CubitVector
    //I- position, point where we want to test, whether or not it
    //- is on the surface.

    virtual CubitPointContainment point_containment( const CubitVector &point );
    virtual CubitPointContainment point_containment( double u, double v );
    virtual CubitPointContainment point_containment( CubitVector &point,
                                                     double u, double v );


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
  
  CubitStatus compare_alignment( IFace2 *first_face,
                                 IFace2 *second_face,
                                 CubitSense &sense );
    //R CubitStaus
    //R- CUBIT_SUCCESS/CUBIT_FAILURE.
    //I- pointer to first SW face.
    //I- pointer to next SW face.
    //O- sense of the two faces.
    //- This function will compare the alignment of the first and
    //- second SW FACES.  This is done by computing the normals
    //- for the faces at a point determined by the bounding box of the
    //- first FACE.  This point is used because of two assumptions:
    //- 1) the two faces are spacially equivalent, so the center of
    //- one will be near the center of the other.
    //- 2) the point_normal function takes the input point and
    //- finds the closes point on the actual surface to do the normal
    //- calculation.
    //- Once the normals are found, they are put into the correct sense
    //- orientation. They are then dotted to determine the sense.
    //- This function is currently used for determining if the sense_
    //- should change when the first_face is removed.
  
  virtual CubitSense get_geometry_sense();
    //- Return the relative surface sense. (see below)
  
  virtual void reverse_sense();
    //- Switch the sense of this Surface wrt the RefFace that owns it:
    //- For SW, this means switch the sense of the RefFace that
    //- owns this Surface with respect to the positive sense of the
    //- first FACE in FACEPtrList_.
  
  void bodysms(DLIList<BodySM*> &bodies);
  void lumps(DLIList<Lump*> &lumps);
  void shellsms(DLIList<ShellSM*> &shellsms);
  void surfaces(DLIList<Surface*> &surfaces);
  void loopsms(DLIList<LoopSM*> &loopsms);
  void curves(DLIList<Curve*> &curves);
  void coedgesms(DLIList<CoEdgeSM*> &coedgesms);
  void points(DLIList<Point*> &points);
    //- topology traversal of TB's; need to implement at this level 'cuz
    //- don't know how many SW entities per TB entity


  virtual CubitSense get_shell_sense( ShellSM* shell_ptr ) const;
  virtual CubitStatus closest_point_uv_guess( CubitVector const& location,
                                              double& u_guess, double& v_guess,
                                              CubitVector* closest_location = NULL,
                                              CubitVector* unit_normal = NULL );


    CubitStatus facet_face(int &number_triangles, int &number_points,
                           int &facet_list_size, GMem *gMem,
                           unsigned short normal_tolerance, float distance_tolerance);

protected: 
  
private:
  ISurface *get_SW_surface();
    //R surface const&
    //R- Returned SW surface object 
    //- This function returns the SW surface object associated with 
    //- the first FACE in FACEPtrList_
  
  CubitSense get_FACE_sense();
    //R CubitSense
    //R- Returned sense value
    //- This function returns the sense of the first SW FACE in FACEPtrList_
    //- wrt its underlying SW surface.
    //- If there is an error getting the sense value, then CUBIT_FORWARD
    //- is returned.
  
  CubitSense get_relative_surface_sense();
    //R CubitSense
    //R- Returned sense value
    //- Returns the sense of the RefFace with respect to the underlying
    //- SW surface.
  
  CubitSense sense_;
    //- The sense of the RefFace that owns this Surface with respect
    //- to the positive sense of the first FACE in FACEPtrList_.
    //- When a Surface is first constructed, this value is arbitrarily
    //- set to CUBIT_FORWARD.
    //- In the case of Surfaces, the normal is used in determining
    //- the relative sense value.
    //- MJP NOTE:
    //- Not only does the RefFace have a sense wrt its Surface, but each
    //- SW FACE has a sense wrt its underlying "surface" object.


//
// utility routine to determine periodic properties of a face
//
void SurfaceSW::periodicInfo(CubitBoolean &bPeriodicU, double &dPeriodU,
                             CubitBoolean &bPeriodicV, double &dPeriodV);


    IFace2 *m_pSWFace;
    SWPart *m_pSWPart;
};


#endif

