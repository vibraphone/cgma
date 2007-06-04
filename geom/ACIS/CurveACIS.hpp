//-------------------------------------------------------------------------
// Filename      : CurveACIS.hpp
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

#ifndef CURVE_ACIS_HPP
#define CURVE_ACIS_HPP

// ********** BEGIN STANDARD INCLUDES      **********
// ********** END STANDARD INCLUDES        **********

// ********** BEGIN ACIS INCLUDES          **********
#if CUBIT_ACIS_VERSION < 1100
#include "kernel/kerndata/top/edge.hxx"
#include "kernel/kerndata/geom/allcurve.hxx"
#else
#include "edge.hxx"
#include "allcurve.hxx"
#endif
// ********** END ACIS INCLUDES            **********

// ********** BEGIN CUBIT INCLUDES         **********
#include "CubitDefines.h"
#include "Curve.hpp"
#include "AcisBridge.hpp"
// ********** END CUBIT INCLUDES           **********

// ********** BEGIN FORWARD DECLARATIONS   **********
class TopologyEntity;
// ********** END FORWARD DECLARATIONS     **********

class CurveACIS : public Curve, public AcisBridge
{
public :
  
  CurveACIS(EDGE* EDGEPtr = NULL);
    //I- EDGE*
    //I- A pointer to the first EDGE which the object will be associated with.
    //- This constructor takes a pointer to the first EDGE to which it will be
    //- connected.
  
  virtual ~CurveACIS() ;
    //- The destructor
    
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
    //- see comments in GeometryEntity.hpp
    //- 
    //- The geometric computations are done using the first underlying
    //- EDGE in the list of EDGEs associated with this CurveACIS.
  
  virtual GeometryQueryEngine* 
  get_geometry_query_engine() const;
    //R GeometryQueryEngine*
    //R- A pointer to the geometric modeling engine associated with
    //R- the object.
    //- This function returns a pointer to the geometric modeling engine
    //- associated with the object.
  
  virtual double measure();
    //R double
    //R- The numeric value of the measure (its units depend on the dimension
    //R- of the RefEntity being "measured")
    //- A generic geometric extent function.
    //- Returns volume for Lump, area for Surface, length for Curve and 
    //- 1.0 for Point
    //-
    //- If there is an error computing the length a value of -1.0 is
    //- returned.
    //- 
    //- The geometric computations are done using the first underlying
    //- EDGE in the list of EDGEs associated with this CurveACIS.
  
  virtual double length_from_u( double parameter1,
                                double parameter2 );
    //R double
    //R- Returned length value
    //I parameter1
    //I- The first parameter value
    //I parameter2
    //I- The second parameter value
    //- This function returns the arc length along the Curve starting from
    //- the point represented by the parameter1 going to the point represented
    //- by parameter2.
    //-
    //- The sign of the returned length value is always positive.
    //- 
    //- The geometric computations are done using the first underlying
    //- EDGE in the list of EDGEs associated with this CurveACIS.
  
  virtual CubitBoolean is_periodic( double& period);
    //R CubitBoolean
    //R- CUBIT_TRUE/CUBIT_FALSE
    //O period
    //O- Returned period value
    //- This function determines whether the underlying geometry of the
    //- Curve is periodic or not.  Returns CUBIT_TRUE if it is and 
    //- CUBIT_FALSE if it is not.
    //-
    //- If it is periodic, then it returns the period in the input
    //- reference variable, "period". This value is set to 0.0 if
    //- the Curve is not periodic.
    //- 
    //- The geometric query is done on the first underlying
    //- EDGE in the list of EDGEs associated with this CurveACIS.
  
  virtual CubitBoolean get_param_range( double& lower_bound,
                                        double& upper_bound );
    //R CubitBoolean
    //R- CUBIT_TRUE/CUBIT_FALSE
    //O lower_bound
    //O- The lower bound of the parametric range.
    //O upper_bound
    //O- The upper bound of the parametric range.
    //- Returns the lower and upper parametric bounds of the 
    //- Curve, based on the extent of the its first underlying EDGE.
    //- All ACIS curves are parametric, so the returned value, which
    //- specifies whether the curve is parametrically defined or not,  
    //- is always CUBIT_TRUE.
    //-
    //- IMPORTANT NOTE:
    //- Note that the lower bound is the parameter value of the start
    //- location of the RefEdge that uses this Curve and the upper
    //- bound is the parameter value of the end location of the RefEdge
    //- that uses this Curve.  This takes into account the sense of the
    //- RefEdge with respect to the Curve (which could be REVERSED).
    //- Hence, the numerical value of the lower parameter bound could be
    //- greater than that of the upper parameter bound.
    //- 
    //- The geometric computations are done using the first underlying
    //- EDGE in the list of EDGEs associated with this CurveACIS.
  
  virtual CubitStatus get_interior_extrema(DLIList<CubitVector*>& interior_points,
                                           CubitSense& return_sense);
    //- Finds the extrema along this RefEdge.  An extremum is defined as
    //- a local min or max in the direction of one of the primary axial directions.
    //- O-interior_points: list of coordinates where the extrema occur.
    //- O-return_sense: Whether the interior extrema are ordered in the
    //-                 FORWARD or REVERSED direction of this RefEdge.
    //-
    //- ***IMPORTANT!!!***
    //-    This function dynamically allocates the CubitVectors appended to
    //-    interior_points.  It is the responsibility of the calling code to
    //-    delete these CubitVectors (or in the case of RefEdge, to make sure
    //-    that *it's* calling code knows that it should delete the CubitVectors)!
  
  virtual CubitStatus closest_point( CubitVector const& location, 
                                     CubitVector& closest_location,
                                     CubitVector* tangent_ptr = NULL,
                                     CubitVector* curvature_ptr = NULL,
                                     double *param = NULL);
    //R void
    //I location
    //I- The point to which the closest point on the Curve is desired.
    //O closest_location
    //O- The point on the Curve, closest to the input location which
    //O- might not be on the Curve.  This is input as a reference 
    //O- so that the function can modify its contents.
    //O tangent_ptr
    //O- The tangent to the Curve (output as a unit vector) at the 
    //O- closest_location.
    //O curvature_ptr
    //O- The curvature of the Curve at the closest_location.
    //- This function computes the point on the Curve closest to the input 
    //- location.
    //-
    //- If the tangent and/or curvature is required, then the calling code
    //- is responsible for allocating space for the CubitVector(s) and
    //- sending in the relevant non-NULL pointers.  If either of these
    //- pointers is NULL, the related quantity is not computed.
    //-
    //- Notes:
    //- The tangent direction is always in the positive direction of the 
    //- *owning RefEdge*, regardless of the positive direction of the
    //- underlying solid model entities.
    //- 
    //- The geometric computations are done using the first underlying
    //- EDGE in the list of EDGEs associated with this CurveACIS.
  
    //**** Added by Jason Kraftcheck, 07/17/98 ****
//   virtual CubitStatus closest_point_trimmed( CubitVector const& from_pt,
//                                              CubitVector & result_pt );
    //R CubitStatus
    //R- CUBIT_SUCCESS (always)
    //I from_pt
    //I- The position for which to find the closest point on the curve.
    //O result_pt
    //O- The closest point on the curve.
    //- This method calls the private method 
    //- AcisGeometryEngine::closest_point( Edge*, CubitVector*) 
    //- to calculate the closest point on a bounded curve.

  virtual CubitPointContainment point_containment( const CubitVector &point );
    //R CubitPointContainment - is the point on bounds of the curve?
    //R- CUBIT_PNT_OFF, CUBIT_PNT_ON, CUBIT_PNT_UNKNOWN
    //I CubitVector
    //I- position to check, known to be on the Curve
    // NOTE: POINT MUST LIE ON THE CURVE FOR THIS FUNCTION TO WORK PROPERLY.
  
#ifdef BOYD14
  void get_tangent( CubitVector const& location, 
                    CubitVector& tangent);
    //- this function returns the tangent vector at the given location
  
  void get_curvature( CubitVector const& location, 
                      CubitVector& curvature);
    //- this function returns the curvature vector at the given location
#endif
  
  EDGE* get_EDGE_ptr() const;
  void set_EDGE_ptr(EDGE* EDGE_ptr);
    // set/get the EDGE associated with this object.
  
#ifdef BOYD14
  static EDGE* get_first_EDGE(const RefEdge *ref_edge);
    //R EDGE*
    //R- Returned EDGE pointer
    //- This function returns a pointer to the first EDGE in the list of
    //- ACIS EDGEs associated with this Curve.
#endif
  
  virtual CubitStatus position_from_u (double u_value,
                                       CubitVector& output_position);
    //R CubitStatus
    //R- CUBIT_SUCCESS/FAILURE
    //I u_value
    //I- The input u parameter value
    //O output_position
    //O- The output position
    //- This function returns the coordinates of a point in the global
    //- (world) space that corresponds to the input parametric position 
    //- in the local space.
    //-
    //- If the input parameter value is not defined for the Curve, then 
    //- the input CubitVector is not modified and CUBIT_FAILURE is
    //- returned. Otherwise, position is appropriately modified and
    //- CUBIT_SUCCESS is returned.
    //-
    //- If the curve is periodic, the input u_value is first "normalized"
    //- to the fundamental period of the Curve before its position
    //- in global space is determined.
    //- 
    //- The geometric computations are done using the first underlying
    //- EDGE in the list of EDGEs associated with this CurveACIS.
  
  double u_from_position (const CubitVector& input_position);
  
  virtual double u_from_arc_length ( double root_param,
                                     double arc_length );
    //R double
    //R- Returned parameter value
    //I root_param
    //I- The parameter value of the "root point"
    //I arc_length
    //I- A distance along the Curve
    //- This function returns the parameter value of the point that is
    //- "arc_length" away from the root point in the
    //- positive sense direction of the owning RefEdge.
    //-
    //- If arc_length is negative, the new point (whose parameter value is
    //- being computed) is in the negative sense direction (along the 
    //- RefEdge) from the root point (whose parameter value is root_param).
    //- 
    //- If the curve is not periodic and the new point, "arc_length" away
    //- from the root point in the appropriate direction, goes beyond
    //- the end point of the first EDGE, that end point is used to generate
    //- the returned parameter value.
    //-
    //- If the curve is periodic and the new point, "arc_length" away
    //- from the root point in the appropriate direction, goes beyond
    //- the end point of the first EDGE, wrap around is done.
    //-
    //- NOTE:
    //- The important assumption that is made in this routine is that
    //- the end points of the RefEdge that owns this CurveACIS are the same
    //- as the end points of the first ACIS EDGE in the list of EDGEs
    //- associated with this CurveACIS.
  
  virtual CubitBoolean is_position_on( const CubitVector &test_position );
    //R CubitBoolean
    //R- CUBIT_TRUE/CUBIT_FALSE
    //I CubitVector
    //I- position, point where we want to test, whether or not it
    //- is on the curve.
  
  GeometryType geometry_type();
    //R GeometryType (enum)
    //R- The enumerated type of the geometric representation
  
  CubitStatus get_point_direction( CubitVector& origin, CubitVector& direction );
  //- Only valid for straight lines
  //- Finds the underlying line's origin and direction unit vector
  //- Returns CUBIT_FAILURE if curve is not a line

  CubitStatus get_center_radius( CubitVector& center, double& radius );
  //- Only valid for arcs
  //- Finds the underlying arc's center point and radius
  //- Returns CUBIT_FAILURE if curve is not an arc

  virtual double start_param();
    //R double parameter
    //R- start parameter of curve with respect to refEdge.
  
  virtual double end_param();
    //R double parameter
    //R- start parameter of curve with respect to refEdge.
    
  virtual CubitBoolean G1_discontinuous( double param,
                                         CubitVector* minus_tangent = NULL,
                                         CubitVector* plus_tangent = NULL );
  
  virtual int validate(const CubitString &user_name,
                       DLIList <TopologyEntity*> &bad_entities);
    //- Check that entity is valid. Returns number of problems detected.

//=========  Add Code by SRS of Cat,  3/3/99 2:37:45 PM  =========
  CubitBoolean is_tolerant();

//=========  Code End by SRS of Cat,  3/3/99 2:37:45 PM  =========
  
  int check_EDGE(EDGE *edge, const CubitString &refentity_name);
    //- check an individual EDGE

  void get_parents_virt( DLIList<TopologyBridge*>& parents );
  void get_children_virt( DLIList<TopologyBridge*>& children );

protected: 
  
private:
  
  const curve* get_ACIS_curve();
    //R curve const&
    //R- Returned ACIS curve object 
    //- This function returns the ACIS curve object associated with 
    //- the first EDGE in EDGEPtrList_
  
  CubitSense get_EDGE_sense();
    //R CubitSense
    //R- Returned sense value
    //- This function returns the sense of the first ACIS EDGE in EDGEPtrList_
    //- wrt its underlying ACIS curve.
    //- If there is an error getting the sense value, then CUBIT_FORWARD
    //- is returned.
  
  void adjust_periodic_parameter(double& param);
  
  friend void run_test_function();
};

// ********** BEGIN INLINE FUNCTIONS       **********
// ********** END INLINE FUNCTIONS         **********

// ********** BEGIN FRIEND FUNCTIONS       **********
// ********** END FRIEND FUNCTIONS         **********

// ********** BEGIN EXTERN FUNCTIONS       **********
// ********** END EXTERN FUNCTIONS         **********

#endif

