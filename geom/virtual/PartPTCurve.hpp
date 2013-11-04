//-------------------------------------------------------------------------
// Filename      : PartPtCurve.hpp
//
// Purpose       : Point Curve for Partition Geometry
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 04/27/02
//-------------------------------------------------------------------------

#ifndef PART_PT_CURVE
#define PART_PT_CURVE

#include "PartitionCurve.hpp"

class PartitionSurface;

class PartPTCurve : public PartitionCurve
{

public:

  PartPTCurve( PartitionSurface* owner );
  
  static PartPTCurve* construct( const CubitSimpleAttrib& attrib,
                                 PartitionSurface* owner );
  
  virtual ~PartPTCurve();
  
  PartitionCurve* split( double );
  
  CubitStatus combine( PartitionCurve* curve );
  
  CubitStatus get_graphics( GMem& result,
                            double angle_tolerance=0,
                            double distance_tolerance=0,
                            double max_edge_length=0);
  
  CubitVector coordinates() const;

  virtual CubitStatus save( CubitSimpleAttrib& attrib );

    /*** Methods from TopologyBridge ***/

  void append_simple_attribute_virt( const CubitSimpleAttrib& );
  void remove_simple_attribute_virt( const CubitSimpleAttrib& );
  void remove_all_simple_attribute_virt();
  CubitStatus get_simple_attribute( DLIList<CubitSimpleAttrib>& );
  CubitStatus get_simple_attribute( const CubitString& name, 
                                    DLIList<CubitSimpleAttrib>& );

    /*** Methods from GeometryEntity ***/

  CubitBox bounding_box() const;
  double measure();
  GeometryType geometry_type();

    /*** Methods from Curve ***/
  
  CubitStatus closest_point( CubitVector const& from, CubitVector& closest,
			     CubitVector* tangent = 0,  CubitVector* curvature = 0,
			     double* param = 0 );
  CubitPointContainment point_containment( const CubitVector& point );
  CubitBoolean is_position_on( const CubitVector& position );
  CubitBoolean G1_discontinuous( double u, CubitVector* minus_tangent = 0,
                                           CubitVector* plus_tangent = 0 );
  CubitStatus get_interior_extrema( DLIList<CubitVector*>& points, CubitSense& sense );
  CubitStatus get_center_radius( CubitVector& center, double& radius );

  CubitBoolean get_param_range( double& lower, double& upper );
  double start_param();
  double end_param();
  CubitBoolean is_periodic( double& period );
  
  double length_from_u( double u1, double u2 );
  double u_from_position( const CubitVector& position );
  CubitStatus position_from_u( double u, CubitVector& result );
  double u_from_arc_length( double root, double length );
  
  void reverse_sense();
  
  //R CubitStatus
  //O- true or false if spline is rational or not.
  //O- the degree of this spline
  //O- the control points
  //O- If rational, weight for each control point
  //O- whether underlying spline is reversed
  virtual CubitStatus get_spline_params( bool &rational,
                                         int &degree,
                                         DLIList<CubitVector> &cntrl_pts,
                                         DLIList<double> &cntrl_pt_weights,
                                         DLIList<double> &knots,
                                         bool &spline_is_reversed
                                       ) const;
  //R CubitStatus
  //O- center - ellipse center point
  //O- normal - normal of the plane of the ellipse
  //O- major_axis - major axis of the ellipse
  //O- radius_ratio - ratio of the length of the major to minor axis.
  virtual CubitStatus get_ellipse_params( CubitVector &center,
                                          CubitVector &normal,
                                          CubitVector &major_axis,
                                          double &radius_ratio ) const;

  private:
  
  PartPTCurve( PartitionSurface* owner, 
               const CubitSimpleAttrib& attrib,
               DLIList<int>& vertex_connectivity );
  
};


#endif
  
  
  
