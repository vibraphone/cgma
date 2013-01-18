//-------------------------------------------------------------------------
// Filename      : SubSurface.cpp
//
// Purpose       : 
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 06/11/98
//-------------------------------------------------------------------------

#include "SubSurface.hpp"
#include "PartitionEngine.hpp"
#include "PartitionLoop.hpp"
#include "PartitionCurve.hpp"
#include "VGLoopTool.hpp"
#include "CubitFacetData.hpp"

#include "GMem.hpp"
#include "GeometryQueryTool.hpp"

//-------------------------------------------------------------------------
// Purpose       : Constructor
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck 
//
// Creation Date : 04/23/02
//-------------------------------------------------------------------------
SubSurface::SubSurface( Surface* real_surf ) 
{
  assert( dynamic_cast<SubEntitySet*>(real_surf->owner()) == 0 );
  new SubEntitySet( real_surf, this );

  // ensure that the sense is the same on both entities
  //if( real_surf->bridge_sense() != this->bridge_sense() )
  //  reverse_bridge_sense();

  geometry_sense = CUBIT_FORWARD;
}

SubSurface::SubSurface( SubSurface* split_from ) 
{
  split_from->sub_entity_set().add_partition( this );
  geometry_sense = CUBIT_FORWARD;
}

PartitionSurface* SubSurface::copy()
{
  return new SubSurface(this);
}

//-------------------------------------------------------------------------
// Purpose       : Destructor
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck 
//
// Creation Date : 06/15/98
//-------------------------------------------------------------------------
SubSurface::~SubSurface()
  { }



//-------------------------------------------------------------------------
// Purpose       : get bounding box
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck 
//
// Creation Date : 06/11/98
//-------------------------------------------------------------------------
//CubitBox SubSurface::bounding_box() const
//  { return partitioned_surface()->bounding_box(); }


//-------------------------------------------------------------------------
// Purpose       : Calculate surface area
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck 
//
// Creation Date : 06/11/98
//-------------------------------------------------------------------------
double SubSurface::measure()
{
  if( ! sub_entity_set().has_multiple_sub_entities() )
    return partitioned_surface()->measure();
  else
    return PartitionSurface::measure();
}
  

GeometryType SubSurface::geometry_type()
  { return partitioned_surface()->geometry_type(); }


CubitStatus SubSurface::get_point_normal( CubitVector& origin, CubitVector& normal )
  { return partitioned_surface()->get_point_normal( origin, normal ); }


CubitStatus SubSurface::closest_point_uv_guess(  
    CubitVector const& location,
    double &u, double &v,
    CubitVector* closest_location,
    CubitVector* unit_normal )
{
  CubitStatus ret = partitioned_surface()->
    closest_point_uv_guess(location, u, v, closest_location, unit_normal);

  if(unit_normal && geometry_sense == CUBIT_REVERSED)
    *unit_normal = -(*unit_normal);

  return ret;
}


//-------------------------------------------------------------------------
// Purpose       : Return the closest point on the surface.
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck 
//
// Creation Date : 06/11/98
//-------------------------------------------------------------------------
CubitStatus SubSurface::closest_point( CubitVector const& a, 
                                       CubitVector* b,
                                       CubitVector* c,
                                       CubitVector* d,
                                       CubitVector* e )
{ 
  CubitStatus ret = partitioned_surface()->closest_point( a, b, c, d, e ); 

  if(c && geometry_sense == CUBIT_REVERSED)
    *c = -(*c);

  return ret;
}

//-------------------------------------------------------------------------
// Purpose       : closest point trimmed
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 03/28/03
//-------------------------------------------------------------------------
void SubSurface::closest_point_trimmed( CubitVector from, CubitVector& result )
{
  PartitionCurve* curve = 0;
  if ( point_containment( from, curve ) == CUBIT_PNT_INSIDE )
    partitioned_surface()->closest_point( from, &result );
  else
    curve->closest_point_trimmed( from, result );
}

//-------------------------------------------------------------------------
// Purpose       : closest point trimmed
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 03/28/03
//-------------------------------------------------------------------------
CubitStatus SubSurface::closest_point_along_vector( CubitVector& from_point, 
                                             CubitVector& along_vector,
                                             CubitVector& point_on_surface)
{
  return partitioned_surface()->closest_point_along_vector( from_point, along_vector, point_on_surface );
}


//-------------------------------------------------------------------------
// Purpose       : get the principal curvatues at a point
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck 
//
// Creation Date : 06/11/98
//-------------------------------------------------------------------------
CubitStatus SubSurface::principal_curvatures( CubitVector const& a,
                                              double& b, double& c,
                                              CubitVector* d )
  { return partitioned_surface()->principal_curvatures(a,b,c,d); }


//-------------------------------------------------------------------------
// Purpose       : return the position for a set of parameter values.
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck 
//
// Creation Date : 06/11/98
//-------------------------------------------------------------------------
CubitVector SubSurface::position_from_u_v( double u, double v )
  { return partitioned_surface()->position_from_u_v( u, v ); }

//-------------------------------------------------------------------------
// Purpose       : return the u and v values at a position
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck 
//
// Creation Date : 06/11/98
//-------------------------------------------------------------------------
CubitStatus SubSurface::u_v_from_position( CubitVector const& location,
                                                 double& u, double& v,
																								 CubitVector* closest )
  { return partitioned_surface()->u_v_from_position(location,u,v,closest); }

//-------------------------------------------------------------------------
// Purpose       : Check if surface is periodic.
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck 
//
// Creation Date : 03/31/03
//-------------------------------------------------------------------------
CubitBoolean SubSurface::is_periodic()
{ 
  double u_period = 0.0, v_period = 0.0;
  bool u_periodic = is_periodic_in_U(u_period);
  bool v_periodic = is_periodic_in_V(v_period);
  if( !u_periodic && !v_periodic )
    return CUBIT_FALSE;
    
  if ( !next_loop(0) ) // sphere
    return CUBIT_TRUE; 
    
  if( u_period < 0.0 )
    u_period = -u_period;
  if( v_period < 0.0 )
    v_period = -v_period; 
    
    // If a single step exceeds this delta of the parameter
    // range, assume the step crossed the "international date
    // line" of the parameterization.  Could probably safely
    // use a smaller value.
  double u_wrap = 0.75 * u_period;
  double v_wrap = 0.75 * v_period;
  
  int num_periodic_loops = 0;
  std::vector<CubitVector> loop_polyline;
  PartitionLoop* loop = 0;
  while( (loop = next_loop(loop) ) != NULL )
  {
      // Get a polyline representation of the loop
    VGLoopTool<PartitionSurface,PartitionLoop,PartitionCoEdge,PartitionCurve,PartitionPoint>
      ::get_loop_polyline( loop->first_coedge(), loop_polyline );
    
    if ( loop_polyline.size() < 2 )
      continue;
    
      // Sum steps (in parameter space) between points in polyline
    std::vector<CubitVector>::iterator itor = loop_polyline.begin(),
                                       end = loop_polyline.end();
    double u_sum = 0, v_sum = 0;
    double u_prev, v_prev;
    u_v_from_position( loop_polyline[loop_polyline.size()-1], u_prev, v_prev );
    for ( ; itor != end; ++itor )
    {
      double u, v;
      u_v_from_position( *itor, u, v );
      
      double u_step = u - u_prev;
      if ( u_step > u_wrap )
        u_step = u_period - u_step;
      else if ( u_step < -u_wrap )
        u_step = -u_period - u_step;
      
      double v_step = v - v_prev;
      if ( v_step > v_wrap )
        v_step = v_period - v_step;
      else if ( v_step < -v_wrap )
        v_step = -v_period - v_step;
        
      u_sum += u_step;
      v_sum += v_step;
      u_prev = u;
      v_prev = v;
    }

      // Sum should be either zero or +/-period if loop
      // is non-periodic or periodic respectively.
    if ( (u_periodic && fabs(u_sum) > 0.5*u_period) ||
         (v_periodic && fabs(v_sum) > 0.5*v_period) )
      num_periodic_loops++;
  }
  
  return CubitBoolean(num_periodic_loops > 0);
}
      
      
    

CubitBoolean SubSurface::is_periodic_in_U( double& period )
  { return partitioned_surface()->is_periodic_in_U(period); }
CubitBoolean SubSurface::is_periodic_in_V( double& period )
  { return partitioned_surface()->is_periodic_in_V(period); }

//-------------------------------------------------------------------------
// Purpose       : check if the surface is singular in either parameter
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck 
//
// Creation Date : 06/11/98
//-------------------------------------------------------------------------
CubitBoolean SubSurface::is_singular_in_U( double u_param )
  {  return partitioned_surface()->is_singular_in_U(u_param); }
CubitBoolean SubSurface::is_singular_in_V( double v_param )
  {  return partitioned_surface()->is_singular_in_V(v_param); }

//-------------------------------------------------------------------------
// Purpose       : Check if the surface is closed along either parameter.
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck 
//
// Creation Date : 06/11/98
//-------------------------------------------------------------------------
CubitBoolean SubSurface::is_closed_in_U()
  { return CUBIT_FALSE; }
CubitBoolean SubSurface::is_closed_in_V()
  { return CUBIT_FALSE; }

//-------------------------------------------------------------------------
// Purpose       : Get uv derivitives
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck 
//
// Creation Date : 06/11/98
//-------------------------------------------------------------------------
CubitStatus SubSurface::uv_derivitives( double u, double v,
                                CubitVector& du, CubitVector& dv )
  { return partitioned_surface()->uv_derivitives(u,v,du,dv); }

//-------------------------------------------------------------------------
// Purpose       : Check if surface is parameterized
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck 
//
// Creation Date : 06/11/98
//-------------------------------------------------------------------------
CubitBoolean SubSurface::is_parametric()
  { return partitioned_surface()->is_parametric(); }

//-------------------------------------------------------------------------
// Purpose       : return the parameter ranges
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck 
//
// Creation Date : 06/11/98
//-------------------------------------------------------------------------
CubitBoolean SubSurface::get_param_range_U(  double& lower, double& upper )
  { return partitioned_surface()->get_param_range_U(lower,upper); }
CubitBoolean SubSurface::get_param_range_V(  double& lower, double& upper )
  { return partitioned_surface()->get_param_range_V(lower,upper); }


//-------------------------------------------------------------------------
// Purpose       : Check if a position lies on this surface.
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck 
//
// Creation Date : 06/11/98
//-------------------------------------------------------------------------
CubitBoolean SubSurface::is_position_on( CubitVector& position )
  { return partitioned_surface()->is_position_on( position ); }



CubitStatus SubSurface::combine( PartitionSurface* dead_surface )
{
  SubSurface* dead = dynamic_cast<SubSurface*>(dead_surface );
  if( !dead )
  {
    assert(0);
    return CUBIT_FAILURE;
  }
  
  return PartitionSurface::combine(dead_surface);
}

CubitSense SubSurface::get_geometry_sense() 
{
  /*
#ifdef ALPHA_TREADSWEEP
  return geometry_sense;
#else
  */
  return CUBIT_FORWARD;
  /*
#endif
  */
}

void SubSurface::reverse_sense()
{
  reverse_loops();
  if( owner() )
    owner()->notify_reversed(this);

  if(geometry_sense == CUBIT_FORWARD)
    geometry_sense = CUBIT_REVERSED;
  else
    geometry_sense = CUBIT_FORWARD;

  int j;
  DLIList<CubitFacetData*> surf_facets;
  this->get_facet_data( surf_facets );
  for(j=surf_facets.size(); j--;)
    surf_facets.get_and_step()->flip();
}

CubitStatus SubSurface::save( CubitSimpleAttrib& attrib )
{
  DLIList<int> curves;
  get_save_topology(curves);
  int id = sub_entity_set().get_id(this);
  return sub_entity_set().save_geometry( id, 2, 0, 0, &curves, 0, attrib );
}


CubitStatus SubSurface::evaluate( double u, double v,
                                  CubitVector *position,                                   
                                  CubitVector *unit_normal,
                                  CubitVector *curvature1,
                                  CubitVector *curvature2 )
{
  CubitStatus ret = partitioned_surface()->evaluate(u, v, position, unit_normal, curvature1, curvature2 ); 

  if(unit_normal && geometry_sense == CUBIT_REVERSED)
    *unit_normal = -(*unit_normal);

  return ret;
}

CubitStatus SubSurface::get_sphere_params
(
  CubitVector &center,
  double &radius
) const
{
  PRINT_ERROR("Currently, Cubit is unable to determine sphere parameters for a SubSurface.\n");
  return CUBIT_FAILURE;
}

CubitStatus SubSurface::get_cone_params
(
   CubitVector &center,
   CubitVector &normal,
   CubitVector &major_axis,
   double &radius_ratio,
   double &sine_angle,
   double &cos_angle
) const
{
  PRINT_ERROR("Currently, Cubit is unable to determine cone parameters for SubSurfaces.\n");
  return CUBIT_FAILURE;
}

CubitStatus SubSurface::get_torus_params
(
  CubitVector &center,
  CubitVector &normal,
  double &major_radius,
  double &minor_radius
) const
{
  PRINT_ERROR("Currently, Cubit is unable to determine torus parameters for SubSurface.\n");
  return CUBIT_FAILURE;
}

CubitStatus SubSurface::get_nurb_params
(
  bool &rational,
  int &degree_u,
  int &degree_v,
  int &num_cntrl_pts_u,
  int &num_cntrl_pts_v,
  DLIList<CubitVector> &cntrl_pts,
  DLIList<double> &cntrl_pt_weights,
  DLIList<double> &u_knots,
  DLIList<double> &v_knots
) const
{
  PRINT_ERROR("Currently, Cubit is unable to determine nurbs parameters for SubSurface.\n");
  return CUBIT_FAILURE;
}
