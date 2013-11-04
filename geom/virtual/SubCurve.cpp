#include "SubCurve.hpp"
#include "GMem.hpp"
#include "GeometryQueryEngine.hpp"

PartitionCurve* SubCurve::split( double param )
{
  if( startParam <= endParam )
  {
    if( (param <= startParam) || (param >= endParam) )
      return 0;
  }
  else if( (param >= startParam) || (param <= endParam) )
    return 0;
    
  return new SubCurve( this, param );
}

CubitStatus SubCurve::combine( PartitionCurve* dead_curve )
{
  SubCurve* curve = dynamic_cast<SubCurve*>(dead_curve);
  if( !curve || curve->real_curve() != this->real_curve() )
    return CUBIT_FAILURE;
  
  double dse = fabs(start_param() - dead_curve->end_param());
  double des = fabs(end_param() - dead_curve->start_param());
  
  if( des <= dse )
    endParam = curve->endParam;
  else 
    startParam = curve->startParam;
  
  return CUBIT_SUCCESS;
}
  
  
       

SubCurve::SubCurve( Curve* curve_ptr ) 
{
  assert( dynamic_cast<SubEntitySet*>(curve_ptr->owner()) == 0 );
  new SubEntitySet( curve_ptr, this );
  
  
  startParam = curve_ptr->start_param();
  endParam = curve_ptr->end_param();
  
  double period;
  if( fabs( startParam - endParam ) < CUBIT_RESABS )
  {
    if( ! curve_ptr->is_periodic( period ) )
    {
      assert( 0 );
    }
    endParam = startParam + period ;
  }
  
  if( curve_ptr->bridge_sense() != this->bridge_sense() )
    reverse_bridge_sense();
}

SubCurve::SubCurve( SubCurve* split_from, double start )
{
  split_from->sub_entity_set().add_partition( this, split_from );
  
  endParam = split_from->endParam;
  startParam = split_from->endParam = start;

  if( CUBIT_REVERSED == split_from->bridge_sense() )
    this->reverse_bridge_sense();
}

SubCurve::~SubCurve()
{
}

int SubCurve::num_partitions() const
{
  DLIList<TopologyBridge*> tmp;
  get_all_partitions( tmp );
  return tmp.size();
}

void SubCurve::get_all_partitions( DLIList<TopologyBridge*>& curves ) const
{
  sub_entity_set().get_owners(curves);
}

CubitBox SubCurve::bounding_box() const
{
/*
  CubitVector start, end;
  SubCurve* me = (SubCurve*)this;
  me->position_from_u( me->start_param(), start );
  me->position_from_u( me->  end_param(),   end );
  CubitBox box( start, end );
  
  DLIList<CubitVector*> extrema;
  CubitSense sense;
  me->get_interior_extrema( extrema, sense );
  if( sense == CUBIT_REVERSED )
    extrema.reverse();
  while( extrema.size() )
  {
    CubitVector* v = extrema.pop();
    box |= *v;
    delete v;
  }
*/  
  return dynamic_cast<Curve*>(partitioned_entity())->bounding_box();
}

double SubCurve::measure()
{
  return fabs(real_curve()->length_from_u( startParam, endParam ));
}

GeometryType SubCurve::geometry_type()
{
  return real_curve()->geometry_type();
}

CubitStatus SubCurve::closest_point( const CubitVector& from, 
                                     CubitVector& closest,
				     CubitVector* tangent,
				     CubitVector* curvature,
				     double* param )
{
  double u;
  CubitStatus result =
     real_curve()->closest_point( from, closest, tangent, curvature, &u );
  if( param )
    *param = u;

  // Snap position to the ends of the SubCurve if
  // necessary.
  if(result)
  {
    double lower, upper;
    double period;

    // If curve is periodic we first need to see if our location
    // can be mapped to the subcurve range just by adding or
    // subtracting multiples of the period.  If no then we
    // will snap it to one of the ends.
    int is_per = is_periodic(period);

    // Get the param range of the SubCurve (will probably 
    // be different from that of the underlying real curve).
    get_ordered_param_range(lower, upper);

    // If we are past the low end...
    if(u < lower)
    {
      // Adjust for periodic curve.
      if(is_per)
      {
        double tmp_u = u;
        while(tmp_u < upper + GEOMETRY_RESABS)
        {
          tmp_u += period;
          if(tmp_u > lower - GEOMETRY_RESABS &&
            tmp_u < upper + GEOMETRY_RESABS)
          {
            u = tmp_u;
          }
        }
      }
      // If the value of u was not adjusted because
      // of a periodic curve go ahead and snap to
      // the lower bound.
      if(u < lower)
      {
        u = lower;
      }
      // Update the position we will return.
      position_from_u(u, closest);
    }
    // If we are past the high end...
    else if(u > upper)
    {
      // Adjust for periodic curve.
      if(is_per)
      {
        double tmp_u = u;
        while(tmp_u > upper - GEOMETRY_RESABS)
        {
          tmp_u -= period;
          if(tmp_u > lower - GEOMETRY_RESABS &&
            tmp_u < upper + GEOMETRY_RESABS)
          {
            u = tmp_u;
          }
        }
      }
      // If the value of u was not adjusted because
      // of a periodic curve go ahead and snap to
      // the upper bound.
      if(u > upper)
      {
        u = upper;
      }
      // Update the position we will return.
      position_from_u(u, closest);
    }
  }

  if( result && tangent )
  {
    CubitVector start, end;
    if( !position_from_u( start_param(), start ) ||
        !position_from_u( end_param(), end ) )
      return CUBIT_FAILURE;
      
    double ds = (start - from).length_squared();
    double de = (end   - from).length_squared();
    const double rsq = GEOMETRY_RESABS * GEOMETRY_RESABS;
    if( ds < rsq && G1_discontinuous( u, &end, &start ) )
      *tangent = start;
    else if( de < rsq && G1_discontinuous( u, &end, &start ) )
      *tangent = end;
  }
  
  return result;
}

CubitPointContainment SubCurve::point_containment( const CubitVector& point )
{
  CubitPointContainment cpc = real_curve()->point_containment( point );
  if( (cpc == CUBIT_PNT_INSIDE) || (cpc == CUBIT_PNT_ON) ||
      (cpc == CUBIT_PNT_BOUNDARY) )
  {
    double u = real_curve()->u_from_position( point );
    double lower, upper;
    get_ordered_param_range( lower, upper );
    if( u < lower || u > upper )
      cpc = CUBIT_PNT_OFF;
  }

  if( (cpc == CUBIT_PNT_INSIDE) || (cpc == CUBIT_PNT_ON) )
  {
    CubitVector start, end;
    position_from_u( start_param(), start );
    position_from_u( end_param(), end );
    if( (point-start).length_squared() < GEOMETRY_RESABS*GEOMETRY_RESABS ||
        (point-end  ).length_squared() < GEOMETRY_RESABS*GEOMETRY_RESABS )
    {
      cpc = CUBIT_PNT_BOUNDARY;
    }
  }
  
  return cpc;
}

CubitBoolean SubCurve::is_position_on( const CubitVector& position )
{
  return real_curve()->is_position_on( position );
}

CubitBoolean SubCurve::G1_discontinuous( double param, 
                                         CubitVector* minus_tan,
					 CubitVector* plus_tan )
{
  double lower, upper;
  get_ordered_param_range( lower, upper );
  if( param <= lower || param >= upper )
    return CUBIT_FALSE;
  
  CubitVector pos, start, end;
  if( !position_from_u( lower, start ) ||
      !position_from_u( upper, end   ) ||
      !position_from_u( param, pos   ) )
    assert(0);
  
  double ds = (pos - start).length_squared();
  double de = (pos - end  ).length_squared();
  const double grs = GEOMETRY_RESABS * GEOMETRY_RESABS;
  if( ds < grs || de < grs )
    return CUBIT_FALSE;
  
  return real_curve()->G1_discontinuous( param, minus_tan, plus_tan );
}

CubitStatus SubCurve::get_interior_extrema( DLIList<CubitVector*>& points, 
                                            CubitSense& sense )
{
  assert( !points.size() );
  
  double lower, upper;
  get_ordered_param_range( lower, upper );
  
  CubitVector start, end;
  if( ! real_curve()->position_from_u( lower, start ) ||
      ! real_curve()->position_from_u( upper, end ) )
    return CUBIT_FAILURE;
  
  DLIList<CubitVector*> list;
  if( ! real_curve()->get_interior_extrema( list, sense ) ) 
    return CUBIT_FAILURE;
  if( list.size() == 0 )
    return CUBIT_SUCCESS;
  
  list.last();
  double param = real_curve()->u_from_position( *list.get() );
  while( (param <= lower) || (param >= upper) )
  {
    delete list.pop();
    if( list.size() == 0 )
      return CUBIT_SUCCESS;
    list.last();
    param = real_curve()->u_from_position( *list.get() );
  }
  
  double ds = (*list.get() - start).length_squared();
  double de = (*list.get() = end  ).length_squared();
  const double grs = GEOMETRY_RESABS * GEOMETRY_RESABS;
  if( ds < grs || de < grs )
  {
    delete list.pop();
    if( list.size() == 0 )
      return CUBIT_SUCCESS;
    list.last();
    param = real_curve()->u_from_position( *list.get() );
  }
  
  while( (param > lower) && (param < upper) )
  {
    points.append( list.pop() );
    if( list.size() == 0 )
      break;
    list.last();
    param = real_curve()->u_from_position( *list.get() );
  }
  
  while( list.size() )
    delete list.pop();
  
  if( sense == CUBIT_FORWARD )
    sense = CUBIT_REVERSED;
  else if( sense == CUBIT_REVERSED )
    sense = CUBIT_FORWARD;
  
  return CUBIT_SUCCESS;
}

CubitStatus SubCurve::get_center_radius( CubitVector& center, double& radius )
  { return real_curve()->get_center_radius(center,radius); }


CubitBoolean SubCurve::get_param_range( double& lower, double& upper )
{
  lower = startParam;
  upper = endParam;
  return CUBIT_TRUE;
}

CubitBoolean SubCurve::get_ordered_param_range( double& lower, double& upper )
{
  if( startParam < endParam )
  {
    lower = startParam;
    upper = endParam;
  }
  else
  {
    lower = endParam;
    upper = startParam;
  }
  return CUBIT_TRUE;
}

double SubCurve::start_param()
  { return startParam; }
double SubCurve::end_param()
  { return endParam; }

CubitBoolean SubCurve::is_periodic( double& period )
{
  return real_curve()->is_periodic( period );
}

double SubCurve::length_from_u( double u1, double u2 )
{
  return real_curve()->length_from_u( u1, u2 );
}

double SubCurve::u_from_position( const CubitVector& p )
{
  double result = real_curve()->u_from_position( p );
  fixup_periodic_param(result);
  return result;
}

CubitStatus SubCurve::position_from_u( double u, CubitVector& p )
{
  return real_curve()->position_from_u( u, p );
}

double SubCurve::u_from_arc_length( double root, double length )
{
  return real_curve()->u_from_arc_length( root, length );
}

void SubCurve::fixup_periodic_param( double& param ) const
{
  double period;
  if( ! real_curve()->is_periodic(period) )
    return;
  
  assert(period > GEOMETRY_RESABS);
  double start, end;
  if ( startParam < endParam ) {
    start = startParam;
    end = endParam;
  } else {
    start = endParam;
    end = startParam;
  }
  
  if( fabs( start-param ) < GEOMETRY_RESABS )
    param = start;
  if( fabs( end-param ) < GEOMETRY_RESABS )
    param = end;

  while ( param > end )
    param -= period;
  while ( param < start )
    param += period;
}

CubitStatus SubCurve::get_graphics( GMem& result,
                            double angle_tolerance,
                            double distance_tolerance,
                            double max_edge_length )
{
  if (POINT_CURVE_TYPE == this->geometry_type())
  {
    result.pointListCount = 0;
    return CUBIT_SUCCESS;
  }

  int i;
  if( !real_curve()->get_geometry_query_engine()->
    get_graphics( real_curve(), &result, angle_tolerance, distance_tolerance, max_edge_length ) )
    return CUBIT_FAILURE;

  if (0 == result.pointListCount)
    return CUBIT_FAILURE;
  
  double lo, hi, start_param, end_param;
  get_ordered_param_range( lo, hi );
  get_param_range(start_param, end_param);
  DLIList<double> param_list( result.pointListCount );

  // don't use first and last points -- this avoids projection
  // close to a seam on a periodic curve.  The projection can
  // jump the seam and cause endpoints of the curve to incorrectly
  // fall inside the parameter range causing extra segments to be drawn
  for( i = 1; i < result.pointListCount-1; i++ )
  {
    CubitVector v( result.point_list()[i].x, 
                   result.point_list()[i].y,
                   result.point_list()[i].z );
    double u = real_curve()->u_from_position( v );
    if( (u >= lo) && (u <= hi) ) 
      param_list.append( u );
  }
  
  CubitVector start_pt, end_pt;
  position_from_u( start_param, start_pt );
  position_from_u( end_param, end_pt );
  
  bool append_start = true, append_end = true;
  if( param_list.size() > 0 )
  {
    CubitVector first_pt, last_pt;
    param_list.last();
    position_from_u( param_list.get(), last_pt );
    param_list.reset();
    position_from_u( param_list.get(), first_pt );
  
    double d1 = (start_pt - first_pt).length_squared();
    double d2 = (start_pt - last_pt ).length_squared();
    double d4 = (  end_pt - first_pt).length_squared();
    double d3 = (  end_pt - last_pt ).length_squared();
    double forward = d1 + d3;
    double reverse = d2 + d4;
    if( reverse < forward )
    {
      CubitVector tmp( start_pt );
      start_pt = end_pt;
      end_pt = tmp;
      d1 = d2;
      d3 = d4;
    }
    
    if( d1 < (GEOMETRY_RESABS*GEOMETRY_RESABS) )
      append_start = false;
    if( d3 < (GEOMETRY_RESABS*GEOMETRY_RESABS) )
      append_end = false;
  }
      
  int index = 0;
  result.allocate_polylines( param_list.size() + 1 );

  if( append_start )
  {
    result.point_list()[index].x = (float)(start_pt.x());
    result.point_list()[index].y = (float)(start_pt.y());
    result.point_list()[index].z = (float)(start_pt.z());
    index++;
  }
  
  param_list.reset();
  CubitVector position;
  for( i = param_list.size(); i--; )
  {
    position_from_u( param_list.get_and_step(), position );
    result.point_list()[index].x = (float)(position.x());
    result.point_list()[index].y = (float)(position.y());
    result.point_list()[index].z = (float)(position.z());
    index++;
  }

  if( append_end )
  {
    result.point_list()[index].x = (float)(end_pt.x());
    result.point_list()[index].y = (float)(end_pt.y());
    result.point_list()[index].z = (float)(end_pt.z());
    index++;
  }
  
  result.pointListCount = index;
  return CUBIT_SUCCESS;
}

  
void SubCurve::reverse_sense()
{
  reverse_point_order();
  if( owner() )
    owner()->notify_reversed(this);
}

CubitStatus SubCurve::save( CubitSimpleAttrib& attrib )
{
  int id = sub_entity_set().get_id(this);
  if( id <= 0 ) return CUBIT_FAILURE;
  
  DLIList<int> end_points(4);
  get_save_topology(end_points);

  return sub_entity_set().save_geometry( id, 1, 0, 0, &end_points, 0, attrib );
}

CubitStatus SubCurve::get_spline_params
(
  bool &rational,    // return true/false
  int &degree,       // the degree of this spline
  DLIList<CubitVector> &cntrl_pts,  // xyz position of controlpoints
  DLIList<double> &cntrl_pt_weights, // if rational, a weight for each cntrl point.
  DLIList<double> &knots,   // There should be order+cntrl_pts.size()-2 knots
  bool &spline_is_reversed
) const
{
  PRINT_ERROR("Currently, Cubit is unable to determine spline parameters for SubCurves.\n");
  return CUBIT_FAILURE;
}

CubitStatus SubCurve::get_ellipse_params
(
  CubitVector &center_vec,
  CubitVector &normal,
  CubitVector &major_axis,
  double &radius_ratio
) const
{
  PRINT_ERROR("Currently, Cubit is unable to determine ellipse parameters for SubCurves.\n");
  return CUBIT_FAILURE;
}
