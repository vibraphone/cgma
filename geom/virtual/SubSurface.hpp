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

  CubitStatus closest_point( CubitVector const& pos, 
                             CubitVector* close = 0, CubitVector* norm = 0,
                             CubitVector* curv1 = 0, CubitVector* curv2 = 0);
  void closest_point_trimmed( CubitVector from, CubitVector& result );
                            
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
  void reverse_sense();

protected:

  virtual PartitionSurface* copy();
  
private:

  SubSurface( SubSurface* split_from );
};

inline Surface* SubSurface::partitioned_surface() const 
  { return dynamic_cast<Surface*>(partitioned_entity()); }


#endif

