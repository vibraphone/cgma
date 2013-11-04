//-------------------------------------------------------------------------
// Filename      : SegmentedCurve.hpp
//
// Purpose       : 
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 04/18/02
//-------------------------------------------------------------------------

#ifndef SEGMENTED_CURVE_HPP
#define SEGMENTED_CURVE_HPP

#include "PartitionCurve.hpp"
#include "PartitionLoop.hpp"
#include "VGArray.hpp"

class PartitionSurface;
class PartitionLump;
class CubitSimpleAttrib;
class CubitTransformMatrix;

class SegmentedCurve : public PartitionCurve
{
public:

  PartitionCurve* split( double param );
  CubitStatus combine( PartitionCurve* dead_curve );
  
  static SegmentedCurve* construct( const CubitSimpleAttrib& attrib,
                                    PartitionEntity* parent );
  
  SegmentedCurve( PartitionSurface* owning_surf, DLIList<CubitVector*>& segments );
  SegmentedCurve( PartitionLump* owning_vol, DLIList<CubitVector*>& segments );
  
  virtual ~SegmentedCurve( );
			
  CubitBoolean is_linear( );
  
  virtual CubitStatus save( CubitSimpleAttrib& attrib );
			
  virtual CubitStatus get_point_direction( CubitVector& origin,
                                           CubitVector& direction );
  
  virtual CubitStatus get_center_radius( CubitVector& center, double& radius );
		
  virtual double length_from_u( double param1, double param2 );
  
  virtual CubitBoolean is_periodic( double& period);
  CubitBoolean is_periodic() const;
  
  virtual CubitBoolean get_param_range( double& lower_bound,
                                        double& upper_bound );
  
  virtual CubitStatus closest_point( 
    CubitVector const& location, 
    CubitVector& closest_location,
    CubitVector* tangent_ptr = NULL,
    CubitVector* curvature_ptr = NULL,
    double* param = NULL);
  
  virtual CubitStatus closest_point_trimmed( 
    CubitVector const& from_pt, CubitVector& result_pt );
  
  
  virtual double u_from_position (const CubitVector& position);
  
  virtual CubitStatus position_from_u (double u_value,
                                       CubitVector& position);
  
  virtual double u_from_arc_length ( double param, double length );
  
  virtual CubitBoolean is_position_on( const CubitVector &position );
  
  virtual double start_param(){ return 0.0; }
  virtual double end_param(){ return point_count() - 1; }
  
  virtual CubitBox bounding_box() const;
  
  virtual double measure();
  
  virtual CubitStatus get_segments( DLIList<CubitVector*>& points );
  
  virtual CubitStatus get_interior_extrema( DLIList<CubitVector*>& points,
                                            CubitSense& return_sense );
  
  virtual int point_count() const;
    
  virtual CubitBoolean G1_discontinuous( double param,
                        CubitVector* minus_tangent = NULL,
                        CubitVector* plus_tangent = NULL );
                        
  CubitPointContainment point_containment( const CubitVector& vect );
  
  CubitStatus get_graphics( GMem& result,
                            double angle_tolerance=0,
                            double distance_tolerance=0,
                            double max_edge_length=0); 
  
  void reverse_sense();
  
  virtual void print_debug_info( const char* prefix, bool pss = true ) const;
  

  virtual void transform( const CubitTransformMatrix& );
  
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
protected:
  
  int segment_from_u( double& u );
    // Return the index of the segment corresponding to the
    // passed parameter value, and change the passed double
    // to be a parameter on the segment.  Returns -1 to
    // indicate an invalid parameter value.
  
  int closest_segment( CubitVector const& location,
                       double* fraction_on_segment = NULL,
                       CubitVector* closest_pt = NULL );
    //Find the segment of the curve closest to the passed location.
    //The return value is the index of the start point of the 
    //segment.  Optional values passed back are the fraction along
    //the segment at which the passed location was closest, and the
    //corresponding position.
  
  static double closest_point_on_segment( 
                                  const CubitVector& segment_base,
                                  const CubitVector& segment_direction,
                                  const CubitVector& from_position );
    //Find the closest location on a bounded segment defined by
    //segment_base and segment_direction to from_position.  The
    //return value is the fraction along segment_direction from 
    //segment_base.  Or a parameter value in the range [0,1], where
    //zero is segment_base, and 1 is segment_base+segment_direction.
  
  double segment_length( int segment_no ) const;
  
  CubitVector position( int index ) const;
  
  void normalize_periodic_parameter( double& u );
  
private:
  SegmentedCurve( PartitionEntity* owning_vol, 
                  const CubitSimpleAttrib& attrib,
                  DLIList<int>& pass_back_vertex_conn );

  SegmentedCurve( SegmentedCurve* split_from, int num_points );
  
  VGArray<CubitVector> point_list;
};

#endif
