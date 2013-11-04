//-------------------------------------------------------------------------
// Filename      : SubSurface.hpp
//
// Purpose       : Geoemtry for the split of a face.
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 06/11/98
//-------------------------------------------------------------------------

#ifndef SUB_SURFACE_HPP
#define SUB_SURFACE_HPP

#include "PartitionSurface.hpp"

class SubSurface : public PartitionSurface
{
public:
		
  SubSurface( Surface* real_surf );
    //- Constructor
    
  CubitStatus combine( PartitionSurface* dead_surface );
		 
  virtual ~SubSurface();
    //- destructor.
		
  inline Surface* partitioned_surface() const;
  
  virtual CubitStatus save(CubitSimpleAttrib&);
  
  //CubitBox bounding_box() const;
  //Use facet-based version inherited from PartitionSurface

  double measure();
  GeometryType geometry_type();


  CubitStatus get_point_normal( CubitVector& point, CubitVector& normal );

  virtual CubitStatus closest_point_uv_guess(  
      CubitVector const& location,
      double& u, double& v,
      CubitVector* closest_location = NULL,
      CubitVector* unit_normal = NULL );

  virtual CubitStatus evaluate( double u, double v,
                                CubitVector *position,                                   
                                CubitVector *unit_normal,
                                CubitVector *curvature1,
                                CubitVector *curvature2 );

  CubitStatus closest_point( CubitVector const& pos, 
                             CubitVector* close = 0, CubitVector* norm = 0,
                             CubitVector* curv1 = 0, CubitVector* curv2 = 0);
  void closest_point_trimmed( CubitVector from, CubitVector& result );

  CubitStatus closest_point_along_vector( CubitVector& from_point, 
                                          CubitVector& along_vector,
                                          CubitVector& point_on_surface);
                            
  CubitStatus principal_curvatures( CubitVector const& loc, 
                                    double& curv1, double& curv2,
                                    CubitVector* closest_location=0 );

  CubitVector position_from_u_v( double u, double v );
  CubitStatus u_v_from_position( CubitVector const& location,
                                 double& u, double& v,
                                 CubitVector* closest_location = NULL );

  CubitBoolean is_periodic();
  CubitBoolean is_periodic_in_U( double& period );
  CubitBoolean is_periodic_in_V( double& period );
  CubitBoolean is_singular_in_U( double u_param );
  CubitBoolean is_singular_in_V( double v_param );
  CubitBoolean is_closed_in_U();
  CubitBoolean is_closed_in_V();

  CubitStatus uv_derivitives( double u, double v,
                              CubitVector &du, CubitVector &dv );
  CubitBoolean is_parametric();
  CubitBoolean get_param_range_U( double& lo, double& hi );
  CubitBoolean get_param_range_V( double& lo, double& hi );
 
  CubitBoolean is_position_on( CubitVector &test_position );

  CubitSense get_geometry_sense();

  CubitStatus get_sphere_params( CubitVector &center,
                                 double &radius ) const;
    //- Only valid for spherical surfaces
    //O center
    //O- The center of the sphere
    //O radius
    //O- The radius of the sphere

  CubitStatus get_cone_params( CubitVector &center,
                               CubitVector &normal,
                               CubitVector &major_axis,
                               double &radius_ratio,
                               double &sine_angle,
                               double &cos_angle ) const;
    //- Only valid for conical surfaces.  Cylinders are a special case of conicals.
    //O center
    //O- 
    //O normal
    //O- 
    //O major_axis
    //O- 
    //O radius_ratio
    //O- 
    //O sine_angle
    //O- 
    //O cos_angle
    //O- 

  virtual CubitStatus get_torus_params( CubitVector &center,
                                        CubitVector &normal,
                                        double &major_radius,
                                        double &minor_radius ) const;
    //- Only valid for torus surfaces.
    //O center
    //O- 
    //O normal
    //O- 
    //O major_radius
    //O- 
    //O minor_radius
    //O- 

  virtual CubitStatus get_nurb_params( bool &rational,
                                       int &degree_u,
                                       int &degree_v,
                                       int &num_cntrl_pts_u,
                                       int &num_cntrl_pts_v,
                                       DLIList<CubitVector> &cntrl_pts,
                                       DLIList<double> &weights,
                                       DLIList<double> &u_knots,
                                       DLIList<double> &v_knots ) const;
  //- Only valid for nurbs surfaces
  //O rational
  //O-   True if the nurb is rational
  //O degree_u
  //O-   The degree of the nurb in the u direction
  //O degree_v
  //O-   The degree of the nurb in the v direction
  //O num_cntrl_pts_u
  //O-   Number of control points in the u direction
  //O num_cntrl_pts_v
  //O-   Number of control points in the v direction
  //O cntrl_pts
  //O-   The control points stored as
  //O-           cntrl_pts[0                ] = pt[u=0][v=0]
  //O-           cntrl_pts[1                ] = pt[u=1][v=0]
  //O-               ...
  //O-           cntrl_pts[num_cntrl_pts_u-1] = pt[u=?][v=0]
  //O-           cntrl_pts[num_cntrl_pts_u  ] = pt[u=0][v=1]
  //O-               ...
  //O weights
  //O-   If rational, weights for each control point, stored in the same
  //O-   order as the control points.  No weights are returned if
  //O-   rational == false
  //O u_knots
  //O-   knot vector in the u direction
  //O v_knots
  //O-   knot vector in the v direction

  void reverse_sense();

protected:

  virtual PartitionSurface* copy();
  
private:

  SubSurface( SubSurface* split_from );
};

inline Surface* SubSurface::partitioned_surface() const 
  { return dynamic_cast<Surface*>(partitioned_entity()); }


#endif

