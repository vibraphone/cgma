//-------------------------------------------------------------------------
// Filename      : SegmentedCurve.cc
//
// Purpose       : 
//
// Special Notes : The parameterization of this is defined such that the
//                 parameter range between two adjacent nodes is one.
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 06/10/98
//-------------------------------------------------------------------------

#include "SegmentedCurve.hpp"

#include "PartitionPoint.hpp"
#include "PartitionCoEdge.hpp"
#include "PartitionSurface.hpp"
#include "PartitionLump.hpp"
#include "GMem.hpp"
#include "PartitionEngine.hpp"
#include "CubitTransformMatrix.hpp"

int SegmentedCurve::segment_from_u( double& param )
{
  if( is_periodic() ) 
    normalize_periodic_parameter( param );
  
  int num_segments = point_count() - 1;
  if( param < -CUBIT_RESABS ||
      param > (num_segments + CUBIT_RESABS) )
    return -1;
  
    // Could call trunc() here, but why?  Need to cast to an
    // int anyway.  Also, the cast to int should be save 
    // because the magnitude of param is tested above.
  int segment_index = (int)param;
  if( segment_index == num_segments )
    --segment_index;
  
  param -= segment_index;
  return segment_index;
}

PartitionCurve* SegmentedCurve::split( double param )
{
  int seg = segment_from_u(param);
  if( seg < 0 )
    return 0;  
  
  CubitVector beg = position(seg);
  CubitVector end = position(seg+1);
  CubitVector pos = beg + param * (end - beg);
  
  if( (pos - end).length_squared() < GEOMETRY_RESABS*GEOMETRY_RESABS )
  {
    seg++;
    param = 0;
    beg = pos = end;
    
    if( seg == point_count() - 1 )
      return 0;
  }
  else if( (pos - beg).length_squared() < GEOMETRY_RESABS*GEOMETRY_RESABS )
  {
    if( seg == 0 )
      return 0;
      
    pos = beg;
    param = 0;
  }
  
  SegmentedCurve* new_curve = new SegmentedCurve( this, point_count() - seg );
  
  new_curve->point_list[0] = pos;
  for( int i = 1; i < new_curve->point_list.size(); i++ )
    new_curve->point_list[i] = point_list[seg + i];
  if( param == 0 )
  {
    point_list.size( seg + 1 );
  }
  else
  {
    point_list.size( seg + 2 );
    point_list[seg+1] = pos;
  }
    
  return new_curve;
}

CubitStatus SegmentedCurve::combine( PartitionCurve* dead_curve )
{
  SegmentedCurve* dead = dynamic_cast<SegmentedCurve*>(dead_curve);
  if( !dead || dead->partitioned_entity() != this->partitioned_entity() )
    return CUBIT_FAILURE;
  
  CubitVector this_start = this->position(0);
  CubitVector dead_start = dead->position(0);
  CubitVector this_end   = this->position( this->point_count() - 1 );
  CubitVector dead_end   = dead->position( dead->point_count() - 1 );
  
  int i;
  bool append, reverse;
  const double tol_sqr = GEOMETRY_RESABS * GEOMETRY_RESABS;
  if( (this_end - dead_start).length_squared() < tol_sqr )
  {
    append = true;
    reverse = false;
  }
  else if( (this_end - dead_end).length_squared() < tol_sqr )
  {
    append = true;
    reverse = true;
  }
  else if( (this_start - dead_end).length_squared()  < tol_sqr )
  {
    append = false;
    reverse = false;
  }
  else if( (this_start - dead_start).length_squared() < tol_sqr )
  {
    append = false;
    reverse = true;
  }
  else
  {
    assert(0);
    return CUBIT_FAILURE;
  }
  
  if( reverse )
    dead->point_list.reverse();
  
  int this_point_count = this->point_list.size();
  int dead_point_count = dead->point_list.size();
  int new_point_count = this_point_count + dead_point_count - 1;
  this->point_list.size(new_point_count);
  
  if( append )
  {
    for( i = 1; i < dead_point_count; i++ )
      point_list[this_point_count + i - 1] = dead->point_list[i];
  }
  else
  {
    for( i = 1; i <= this_point_count; i++ )
      point_list[new_point_count-i] = point_list[this_point_count-1];
    
    for( i = 0; i < dead_point_count; i++ )
      this->point_list[i] = dead->point_list[i];
  }
    
  return CUBIT_SUCCESS;
}  
      
  

CubitBoolean SegmentedCurve::is_linear()
{
  int i, count = point_count();
  if( count == 2 ) return CUBIT_TRUE;
  if( is_periodic() ) return CUBIT_FALSE;
  
  CubitVector mean(0.0,0.0,0.0);
  for( i = 1; i < count; i++ ) mean += (position(i)-position(i-1));
  mean.normalize();
  for( i = 1; i < count; i++ )
  {
    if( (1-(mean % ~(position(i)-position(i-1)))) > CUBIT_RESABS )
      return CUBIT_FALSE;
  }
  return CUBIT_TRUE;
}

CubitStatus SegmentedCurve::get_point_direction( 
  CubitVector& origin, CubitVector& direction )
{
  if( ! is_linear() ) return CUBIT_FAILURE;
  origin = position(0);
  direction = ~(position(point_count()-1) - origin);
  return CUBIT_SUCCESS;
}

double SegmentedCurve::length_from_u( double param1, double param2 )
{
  if( param1 > param2 )
  {
    double temp = param1;
    param1 = param2;
    param2 = temp;
  }
  
  if( is_periodic() )
  {
    normalize_periodic_parameter( param1 );
    normalize_periodic_parameter( param2 );
  }
  
  int first = (int)floor( param1 );
  int last = (int)floor( param2 );
  int seg_count = point_count() - 1;

  //check for errors, round a little
  if( first < 0 )
  {
    if( param1 > -CUBIT_RESABS ) first = 0;
    else return -1.0;
  }
  if( last < 0 )
  {
    if( param2 > -CUBIT_RESABS ) last = 0;
    else return -1.0;
  }
  if( last >= seg_count )
  {
    if( param2 < (CUBIT_RESABS + seg_count) )
      last--;
    else return -1.0;
  }
  if( first >= seg_count )
  {
    if( param1 < (CUBIT_RESABS + seg_count) )
      first--;
    else return -1.0;
  }
    
  //calculate the remaining length on the first segment
  
  double length_sum = segment_length( first ) * (1 - (param1 - first) );
  
  //add up the remaining segments upto and including last
  
  for( int i = first + 1; i <= last; i++ )
    length_sum += segment_length( i );
    
  // by including the last segment, we (may have) gone too far,
  // so subtract back off the extra on the last segment
  
  length_sum -= segment_length(last) * (1 - (param2 - last) );
  
  return length_sum;

}

CubitBoolean SegmentedCurve::get_param_range( double& lower, double& upper )
{
  lower = 0.0;
  upper = double(point_count() - 1);
  return CUBIT_TRUE;
}

CubitStatus SegmentedCurve::closest_point(
  CubitVector const& location,
  CubitVector& closest_location,
  CubitVector* tangent_ptr,
  CubitVector* curvature_ptr,
  double* param)
{
  double seg_param;
  int segment_no = closest_segment( location, &seg_param, &closest_location );
  
  if ( Surface* surf = dynamic_cast<Surface*>(partitioned_entity()) )
  {
    CubitVector seg_point(closest_location);
    surf->closest_point( seg_point, &closest_location );
  }
  
  if( tangent_ptr != NULL )
  {
    //First check if the closest point is at the join of
    //two segments.  If so, return the average of the two 
    //tangets.
    if( (seg_param < CUBIT_RESABS) && (segment_no > 0) )
    {
      CubitVector seg1_start = position( segment_no - 1 );
      CubitVector seg2_end   = position( segment_no + 1 );
      *tangent_ptr = seg2_end - seg1_start;
    }
    
    //ditto
    else if( ( (1.0 - seg_param) < CUBIT_RESABS ) && 
             ( segment_no < (point_count() - 2) ) )
    {
      CubitVector seg1_start = position( segment_no );
      CubitVector seg2_end   = position( segment_no + 2 );
      *tangent_ptr = seg2_end - seg1_start;
    }
    
    //The segment_no is the index of the point beginning the segment.
    //A segment_no equal to point_count() - 1, the last point is 
    //possible.  So, for any segment point but the last, return the
    //tangent of the segment beginning with the point at segment_no.
    else if( segment_no < point_count() - 1 )
    {
      *tangent_ptr = this->position( segment_no + 1 ) 
        - this->position( segment_no );
    }
    
    else if( segment_no == 0 )
    {
      *tangent_ptr = this->position( segment_no + 1 )
        - this->position( segment_no );
    }
    
    //If the point is the last one, return the direction of
    //the last segment (this point and the previous one.)
    else 
    {
      *tangent_ptr = this->position( segment_no ) 
        - this->position( segment_no - 1 );
    }
    
    tangent_ptr->normalize();
  }
  
  //There is never a valid curvature for a segmented curve.  
  //The curvature is either infinity in the interior of segments,
  //or zero at the points between segments.  Just return zero.
  if( curvature_ptr != NULL )
  {
    curvature_ptr->set( 0., 0., 0. );
  }

  if (param)
  {
    *param = double(segment_no) + seg_param;
  }
  
  return CUBIT_SUCCESS;
}

CubitStatus SegmentedCurve::closest_point_trimmed( 
  CubitVector const& from_pt, CubitVector& result_pt ) 
{
  return closest_point( from_pt, result_pt );
}

int SegmentedCurve::closest_segment( CubitVector const& location,
                                     double* seg_fraction,
                                     CubitVector* closest_pt )
{  
  int length_ = point_count();
  
  int segment_no = 0;
  CubitVector start = position( 0 );
  CubitVector end = position( 1 );
  CubitVector dir = end - start;
  double seg_param = closest_point_on_segment( start, dir, location );
  CubitVector point = start + seg_param * dir;
  if( closest_pt ) *closest_pt = point;
  double shortest = (point - location ).length_squared();
   
  for( int i = 2; i < length_; i++ )
  {
    start = end;
    end = position( i );
    dir = end - start;
    double param = closest_point_on_segment( start, dir, location );
    point = start + param * dir;
    double dist = (point - location).length_squared();
    if( dist < shortest )
    {
      shortest = dist;
      if( closest_pt ) *closest_pt = point;
      segment_no = i - 1;
      seg_param = param;
    }
  }
  
  if( seg_fraction ) *seg_fraction = seg_param;  
  return segment_no;
}

double SegmentedCurve::u_from_position( const CubitVector& position )
{
  double seg_param;
  int segment = closest_segment( position, &seg_param );
  return double(segment) + seg_param;
}

CubitStatus SegmentedCurve::position_from_u( double u, CubitVector& p)
{
  double fraction = u;
  int segment_no = segment_from_u(fraction);
  if( segment_no < 0 )
    return CUBIT_FAILURE;
    
  CubitVector s = position( segment_no );
  CubitVector e = position( segment_no+1 );
  p = s + fraction * ( e - s );
  
  if ( Surface* surf = dynamic_cast<Surface*>(partitioned_entity()) )
  {
    CubitVector seg_pos(p);
    surf->closest_point( seg_pos, &p );
  }

  return CUBIT_SUCCESS;
}

double SegmentedCurve::u_from_arc_length( double param, double length )
{
  if( fabs(length) < CUBIT_RESABS ) 
    return param;
  
  double fraction = param;  
  int segment = segment_from_u( fraction );
  if( segment < 0 )
    return param;
  
  //Find the length of the remaining portion of the first
  //segment (after the base param.)
  double seg_len = segment_length( segment );
// BWC  double len_sum = seg_len * ( 1.0 - param + segment );
  double len_sum;
  
  if(length < 0.0)
    len_sum = seg_len*(param-segment);
  else
    len_sum = seg_len * ( 1.0 - param + segment );
  
  //If the passed base param and length correspond to a
  //portion of the curve which lies entirely on a single
  //segment, return the passed base param plus the 
  //fraction of the segment (the parameter range for a 
  //segment is one.)
//BWC  if( len_sum >= length ) 
  if(length < 0.0)
  {
    if(-len_sum < length)
      return param + ((seg_len > CUBIT_DBL_MIN ) ? length / seg_len : 0.0);
  }
  else
  {
    if( len_sum >= length ) 
      return param + ((seg_len > CUBIT_DBL_MIN ) ? length / seg_len : 0.0);
  }
  
  //Increment the total length until we have passed up
  //the specified arc length.
  /*
  BWC
  while( length > len_sum)
  {
    segment++;
    if( segment >= point_count() - 1 ) return -1.0;
    seg_len = segment_length( segment );
    len_sum += seg_len;
    if(fabs(len_sum-length) < GEOMETRY_RESABS)
      len_sum = length;
 }
 */
  if(length < 0.0)
  {
    while( fabs(length) > len_sum )
    {
      segment--;
      if( segment < 0 ) return -1.0;
      seg_len = segment_length( segment );
      len_sum += seg_len;
      if(fabs(len_sum+length) < GEOMETRY_RESABS)
        len_sum = -length;
    }
  }
  else
  {
    while( length > len_sum )
    {
      segment++;
      if( segment >= point_count() - 1 ) return -1.0;
      seg_len = segment_length( segment );
      len_sum += seg_len;
      if(fabs(len_sum-length) < GEOMETRY_RESABS)
        len_sum = length;
    }
  }
  
  //Now subtract off the extra length on the last segment,
  //the amount passed the specified arc length.  The ratio
  //of this length to the total length of the last segment
  //is the about we overshot the desired parameter value.
  if(length < 0.0)
  {
    if( seg_len > CUBIT_DBL_MIN )
      return (double)(segment) + (len_sum + length)/seg_len;
    else 
      return (double)(segment);
  }
  else
  {
    if( seg_len > CUBIT_DBL_MIN )
      return (double)(segment + 1) - (len_sum - length)/seg_len;
    else 
      return (double)(segment + 1);
  }
}

CubitBoolean SegmentedCurve::is_position_on( const CubitVector& pos )
{
  //Get the segment to test on
  CubitVector v;
  closest_segment( pos, NULL, &v);
  
  return (pos - v).length() <= GEOMETRY_RESABS ? CUBIT_TRUE : CUBIT_FALSE;
}

CubitBox SegmentedCurve::bounding_box() const
{
  CubitBox b( position( 0 ) );
  for( int i = 1; i < point_count(); i++ )
    b |= position( i );

  return b;
}

double SegmentedCurve::measure()
{
  double len = 0;
  for( int i = 0; i < point_count() - 1; i++ )
  {
    len += segment_length( i );
  }
  return len;
}

CubitStatus SegmentedCurve::get_interior_extrema( 
                                DLIList<CubitVector*>& interior_points,
                                CubitSense& return_sense )
{
  return_sense = CUBIT_FORWARD;
  
  //For each point, check the segments to either side.  If any of deltaX, 
  //deltaY or deltaZ have opposite signs for the two segments, then the
  //point is an extrema in x, y, or z, respectively.
  int count = point_count() - 1; //stop with second-to-last point
  int prev_index = 0;

  for( int i = 1; i < count; i++ )
  {
    CubitVector curr = position(i+1) - position(i);
    CubitVector prev = position( i ) - position( prev_index );
    if( (prev.x()*curr.x() < 0.0) || (prev.y()*curr.y()<0.0) || (prev.z()*curr.z()<0.0) )
    {
      interior_points.append( new CubitVector( position(i) ) );
      prev_index = i;
    }
  }
  return CUBIT_SUCCESS;
}

CubitStatus SegmentedCurve::get_center_radius( CubitVector&, double& )
  { return CUBIT_FAILURE; }
  

double SegmentedCurve::closest_point_on_segment( 
                                              const CubitVector& base,
                                              const CubitVector& dir,
                                              const CubitVector& location )
{
  if( dir.length_squared() < CUBIT_DBL_MIN ) return 0.0;
  
  //Find closest location to infinite line
  double param = dir % ( location - base ) / dir.length_squared();
  
  //Trim to segment
  if( param < 0.0 ) param = 0.0;
  else if( param > 1.0 ) param = 1.0;
  
  return param;
}

//-------------------------------------------------------------------------
// Purpose       : Check for G1 discontinuities
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 05/03/00
//-------------------------------------------------------------------------
CubitBoolean SegmentedCurve::G1_discontinuous( double param,
                                            CubitVector* minus_tangent,
                                            CubitVector* plus_tangent )
{
  if( is_periodic() ) normalize_periodic_parameter( param );
  
  double lower = floor( param );
  double upper = ceil( param );
  int first_point = -1;
  int last = point_count() - 1;
  if( last < 2 ) return CUBIT_FALSE;

  //If the parameter value corresponds to the start of a segment,
  //and that segment is not the first one
  if( ((param - lower) < CUBIT_RESABS) && (int(lower) != 0) )
    first_point = int(lower) - 1;
  //If the parameter value corresponds to the end of a segment,
  //and that segment is not the last one
  else if( ((upper - param) < CUBIT_RESABS) && (int(upper) < last ) )
    first_point = int(upper);
  //If the parameter value corresponds to the start or end of the
  //curve, and the curve is a closed, periodic curve.
  else if( is_periodic() && (fabs(start_param() - end_param()) < CUBIT_RESABS) &&
           ((((param - lower) < CUBIT_RESABS) && (int(lower) == 0)) ||
            (((upper - param) < CUBIT_RESABS) && (int(upper) == last)) ) )
    first_point = last - 1;
  //Otherwise the point is in the interior of a segment
  else return CUBIT_FALSE;
  
  int third_point = first_point + 2;
  if( third_point > last )
  {
    if( is_periodic() )
    {
      third_point = (third_point + 1) % point_count();
    }
    else
    {
      return CUBIT_FALSE;
    }
  }

  CubitVector tan1, tan2;
  tan1 = position( first_point + 1 ) - position( first_point );
  tan2 = position( third_point ) - position( first_point + 1 );
  if( (tan1 * tan2).length_squared() < CUBIT_RESABS )
    return CUBIT_FALSE;
  
  if( minus_tangent ) *minus_tangent = tan1;
  if( plus_tangent )  *plus_tangent = tan2;
  return CUBIT_TRUE;
}


//-------------------------------------------------------------------------
// Purpose       : If the curve is periodic, adjust the param value to
//                 be in the base range.
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 05/15/00
//-------------------------------------------------------------------------
void SegmentedCurve::normalize_periodic_parameter( double& u ) 
{
  //The parameter range for a segmented curve is always (by my implementation)
  //from zero to some positive n.  And the period is equal to that n.  This 
  //simplifies the following a lot. -jk
  double period;
  if( is_periodic( period ) )
  {
    if( u > period )
      u = fmod( u, period );
    else if( u < -period ) 
      u = fmod( u, period ) + period;
    else if( u < 0.0 )
      u += period;  //fmod() seems to be expensive, so avoid it here.
  }
}

//-------------------------------------------------------------------------
// Purpose       : Constructor
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 05/15/00
//-------------------------------------------------------------------------
SegmentedCurve::SegmentedCurve( PartitionSurface* surf, 
                                DLIList<CubitVector*>& points )
  : point_list( points.size() )
{
  points.reset();
  for( int i = 0; i < points.size(); i++ )
    point_list[i] = *points.get_and_step();
  
  surf->sub_entity_set().add_lower_order( this );
}
SegmentedCurve::SegmentedCurve( PartitionLump* vol, 
                                DLIList<CubitVector*>& points )
  : point_list( points.size() )
{
  points.reset();
  for( int i = 0; i < points.size(); i++ )
    point_list[i] = *points.get_and_step();
  
  vol->sub_entity_set().add_lower_order( this );
}


SegmentedCurve::SegmentedCurve( SegmentedCurve* curve, int ptcount )
  : point_list(ptcount)
{
  curve->sub_entity_set().add_lower_order(this);
}
 

SegmentedCurve::~SegmentedCurve()
{
}

int SegmentedCurve::point_count() const
{
  return point_list.size();
}

double SegmentedCurve::segment_length( int segment_no ) const
{
  return (point_list[segment_no+1] - point_list[segment_no]).length();
}

CubitVector SegmentedCurve::position( int index ) const
{
  return point_list[index];
}

  
CubitStatus SegmentedCurve::get_graphics( GMem& result )
{
  result.allocate_polylines( point_count()-1 );
  for( int i = 0; i < point_count(); i++ )
  {
    CubitVector v = position(i);
    result.point_list()[i].x = (float)(v.x());
    result.point_list()[i].y = (float)(v.y());
    result.point_list()[i].z = (float)(v.z());
  }
  result.pointListCount = point_count();
  return CUBIT_SUCCESS;
}

CubitBoolean SegmentedCurve::is_periodic( double& period )
{
  if( is_periodic() )
  {
    period = end_param() - start_param();
    return CUBIT_TRUE;
  }
  return CUBIT_FALSE;
}

CubitBoolean SegmentedCurve::is_periodic() const
{
  return start_point() == end_point();
}

CubitPointContainment SegmentedCurve::point_containment( const CubitVector& pos )
{
  const double tol_sqr = GEOMETRY_RESABS * GEOMETRY_RESABS;
  if( (start_point()->coordinates() - pos).length_squared() < tol_sqr )
    return CUBIT_PNT_BOUNDARY;
  else if( (end_point()->coordinates() - pos).length_squared() < tol_sqr )
    return CUBIT_PNT_BOUNDARY;
  else if( is_position_on( pos ) )
    return CUBIT_PNT_ON;
  else 
    return CUBIT_PNT_OFF;
}


CubitStatus SegmentedCurve::get_segments( DLIList<CubitVector*>& list )
{
  for( int i = 0; i < point_count(); i++ )
    list.append( new CubitVector( position(i) ) );
  return CUBIT_SUCCESS;
}

void SegmentedCurve::reverse_sense()
{
  point_list.reverse();
  reverse_point_order();
  if( owner() )
    owner()->notify_reversed(this);
}

CubitStatus SegmentedCurve::save( CubitSimpleAttrib& attrib )
{
  DLIList<CubitVector*> segments;
  get_segments(segments);
  
  DLIList<int> topo;
  get_save_topology( topo );
  
  int id = sub_entity_set().get_id(this);
  
  return sub_entity_set().save_geometry( id, 1, &segments, 0, &topo, 0, attrib );
}

SegmentedCurve* SegmentedCurve::construct( CubitSimpleAttrib* attrib,
                                           PartitionEntity* parent )
{
  PartitionSurface* owning_surf = dynamic_cast<PartitionSurface*>(parent);
  PartitionLump* owning_lump = dynamic_cast<PartitionLump*>(parent);
  if( !owning_surf && !owning_lump )
    return 0;
  
  DLIList<int> vertex_conn;
  SegmentedCurve* result = new SegmentedCurve( parent, *attrib, vertex_conn );
  
  if( vertex_conn.size() != 4 )
  {
    delete result;
    return 0;
  }
  
  PartitionPoint *start = 0, *end = 0;
  PartitionEntity* ent;
  vertex_conn.reset();
  int set_id = vertex_conn.get_and_step();
  int ent_id = vertex_conn.get_and_step();
  ent = PartitionEngine::instance().entity_from_id(set_id,ent_id,parent->sub_entity_set());
  start = dynamic_cast<PartitionPoint*>(ent);
  set_id = vertex_conn.get_and_step();
  ent_id = vertex_conn.get_and_step();
  ent = PartitionEngine::instance().entity_from_id(set_id,ent_id,parent->sub_entity_set());
  end = dynamic_cast<PartitionPoint*>(ent);
  if( !start || !end )
  {
    delete result;
    return 0;
  }
  
  result->start_point(start);
  result->end_point(end);
  return result;
}
  
  
SegmentedCurve::SegmentedCurve( PartitionEntity* vol, 
                                CubitSimpleAttrib& attrib,
                                DLIList<int>& vertex_conn )
{
  DLIList<CubitVector*> points;
  DLIList<int> junk;
  vol->sub_entity_set().add_lower_order( this, attrib, 1, points, junk, vertex_conn, junk );
  
  points.reset();
  point_list.size(points.size());
  for( int i = 0; i < points.size(); i++ )
  {
    CubitVector* pt = points.get_and_step();
    point_list[i] = *pt;
    delete pt;
  }
}

void SegmentedCurve::transform( const CubitTransformMatrix& xform )
{
  PartitionCurve::transform(xform);
  for( int i = 0; i < point_list.size(); i++ )
    point_list[i] = xform * point_list[i];
}

void SegmentedCurve::print_debug_info( const char* prefix, bool pss ) const
{
  if( prefix == 0 ) prefix = "";
  PartitionCurve::print_debug_info( prefix, pss );
  PRINT_INFO("%sSegmentedCurve %p\n", prefix, static_cast<void*>(this) );
  PRINT_INFO("%s%d segment points:\n", prefix, point_list.size());
  int i;
  for( i = 0; i  < point_list.size() - 1; i+= 2 )
  {
    PRINT_INFO("%s  (%f, %f, %f), (%f, %f, %f),\n",
      prefix, point_list[i  ].x(), point_list[i  ].y(), point_list[i  ].z(),
              point_list[i+1].x(), point_list[i+1].y(), point_list[i+1].z());
  }
  if( i == point_list.size() - 1 )
  {
    PRINT_INFO("%s  (%f, %f, %f)\n", prefix, point_list[i].x(), 
       point_list[i].y(), point_list[i].z());
  }
}

