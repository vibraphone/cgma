//-------------------------------------------------------------------------
// Filename      : CompositeCurve.cpp
//
// Purpose       : Geometry defined as the joining of a chain of curves.
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 12/19/01
//-------------------------------------------------------------------------

#include <math.h>

#include "CompositeCurve.hpp"
#include "CompositeCoEdge.hpp"
#include "CompositeLoop.hpp"
#include "CompositePoint.hpp"
#include "VirtualQueryEngine.hpp"
#include "CompositeEngine.hpp"

//#include "GfxDebug.hpp"
#include "GMem.hpp"

//-------------------------------------------------------------------------
// Purpose       : Constructor
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 12/19/01
//-------------------------------------------------------------------------
CompositeCurve::CompositeCurve( Curve* curve )
  : hiddenSet(0), 
    firstCoEdge(0),
    startPoint(0), 
    endPoint(0), 
    startNext(0), 
    endNext(0),
    stitchNext(0),
    HadBridgeRemoved(0)
{
  compGeom = new CompositeGeom(1);
  compGeom->append( curve, CUBIT_FORWARD );
  if( curve->owner() )
    curve->owner()->swap_bridge( curve, this, false );
  curve->owner(this);
}

CompositeCurve::CompositeCurve( CompositeGeom* geometry )
  : compGeom( geometry ),
    hiddenSet( 0 ),
    firstCoEdge( 0 ),
    startPoint( 0 ),
    endPoint( 0 ),
    startNext( 0 ),
    endNext( 0 ),
    stitchNext(0),
    HadBridgeRemoved(0)
{
  for( int i = 0; i < compGeom->num_entities(); i++ )
  {
    GeometryEntity* entity = compGeom->entity(i);
    assert( !entity->owner() );
    entity->owner(this);
  }
}

//-------------------------------------------------------------------------
// Purpose       : Point-curve constructor
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 02/27/04
//-------------------------------------------------------------------------
CompositeCurve::CompositeCurve( CompositePoint* point )
  : hiddenSet(0), 
    firstCoEdge(0),
    startPoint(0), 
    endPoint(0), 
    startNext(0), 
    endNext(0),
    stitchNext(0),
    HadBridgeRemoved(0)
{
  compGeom = new CompositeGeom(1);
  //compGeom->append( point, CUBIT_FORWARD );
  start_point(point);
  end_point(point);
}


CompositeCurve::~CompositeCurve()
{
  while( firstCoEdge )
    remove( firstCoEdge );
  
  for( int i = 0; i < num_curves(); i++ )
    if( get_curve(i)->owner() == this )
      get_curve(i)->owner(0);
      
  start_point(0);
  end_point(0);
  
  unstitch_all();
  
  delete hiddenSet;
  delete compGeom;
  hiddenSet = (HiddenEntitySet*)0xbdbdbdbd;
  compGeom = (CompositeGeom*)0xbdbdbdbd;
  assert(!startNext);
  assert(!endNext);
}


CompositeCurve* CompositeCurve::next( const CompositePoint* around )
{
  return around == startPoint ? startNext :
         around == endPoint ? endNext : 0;
}

CubitStatus CompositeCurve::start_point( CompositePoint* pt )
{
  return set_point( true, pt );
}

CubitStatus CompositeCurve::end_point( CompositePoint* pt )
{
  return set_point( false, pt );
}

CubitStatus CompositeCurve::set_point( bool start, CompositePoint* pt )
{
  CompositePoint*& my_point = start ? startPoint : endPoint;
  CompositeCurve*& next_crv = start ? startNext : endNext;

  if( pt == my_point )
    return CUBIT_SUCCESS;
  
  if( my_point && startPoint != endPoint )
  {
    CompositeCurve* prev = my_point->firstCurve;
    assert( prev != NULL );

    if( prev == this )
    {
      my_point->firstCurve = next_crv;
    }
    else
    {
      CompositeCurve* curve = prev->next(my_point);
      while( curve != this )
      {
        prev = curve;
	      curve = prev->next(my_point);
      }
      
      if( prev->startPoint == my_point )
      {
        assert(prev->startNext == this);
      	prev->startNext = next_crv;
      }
      
      if( prev->endPoint == my_point )
      {
        assert(prev->endNext == this);
        prev->endNext = next_crv;
      }
    }
  }

  my_point = pt;
  if (!pt)
  {
    next_crv = 0;
  }
  else if (startPoint == endPoint)
  {
    next_crv = start ? endNext : startNext;
  }
  else
  {
    next_crv = pt->firstCurve;
    pt->firstCurve = this;
  }

  return CUBIT_SUCCESS;
}
    

CompositePoint* CompositeCurve::other_point( CompositePoint* pt )
{ 
  CompositePoint* sp = start_point();
  CompositePoint* ep = end_point();
  return (pt == sp) ? ep : (pt == ep) ? sp : 0;
}

CompositeCurve* CompositeCurve::split( Curve* curve )
{
  int index = index_of( curve );
  if( (index < 0) || (index == (num_curves()-1)) )
    return 0;
  
  CompositeGeom* new_geom = compGeom->split( index );
  if( !new_geom )
    return 0;
    
  for( int i = 0; i < new_geom->num_entities(); i++ )
    new_geom->entity(i)->owner( 0 );
  
  return new CompositeCurve( new_geom );
}


CubitStatus CompositeCurve::combine( CompositeCurve* curve, bool prepend )
{
  int start, stop;
  if ( prepend ) {
    start = 0;
    stop = curve->compGeom->num_entities();
  } else {
    start = compGeom->num_entities();
    stop = start + curve->compGeom->num_entities();
  }
  
  compGeom->merge( *(curve->compGeom), prepend );
  if( curve->hiddenSet )
    hidden_entities().merge( curve->hiddenSet );
    
  for( int i = start; i < stop; i++ )
  {
    TopologyBridge* bridge = compGeom->entity(i);
    assert( bridge->owner() == curve );
    bridge->owner(this);
  }
  
  return CUBIT_SUCCESS;
}

bool CompositeCurve::has_parent_composite_surface() const
{
  CompositeCoEdge* coedge;
  for( coedge = first_coedge(); coedge; coedge = next_coedge(coedge) )
    if( coedge->get_loop() )
      return true;
  return false;
}
  
void CompositeCurve::reverse()
{
  compGeom->reverse();
}

CompositePoint* CompositeCurve::common_point( CompositeCurve* curve )
{
  CompositePoint* result = 0;
  CompositePoint* this_sp = start_point();
  CompositePoint* this_ep = end_point();
  CompositePoint* othr_sp = curve->start_point();
  CompositePoint* othr_ep = curve->end_point();
  
  if( this_sp == othr_sp || this_sp == othr_ep )
    result = this_sp;
  else if( this_ep == othr_sp || this_ep == othr_ep )
    result = this_ep;
  
  return result;
}
/*
CompositeCoEdge* CompositeCurve::find_coedge( CompositeSurface* surface )
{
  CoEdge* coedge = 0;
  while( coedge = next_coedge( coedge ) )
    if( coedge->get_loop() && coedge->get_loop()->get_surface() == surface )
      break;
  return coedge;
}
*/
CubitStatus CompositeCurve::add( CompositeCoEdge* coedge )
{
  if( coedge->myCurve )  
  {
    assert(0);
    return CUBIT_FAILURE;
  }
  
  coedge->myCurve = this;
  coedge->nextOnCurve = firstCoEdge;
  firstCoEdge = coedge;
  
  return CUBIT_SUCCESS;
}

CubitStatus CompositeCurve::remove( CompositeCoEdge* coedge )
{
  if( coedge->myCurve != this )
    return CUBIT_FAILURE;
  
  if( coedge == firstCoEdge )
  {
    firstCoEdge = coedge->nextOnCurve;
  }
  else
  {
    CompositeCoEdge *prev = firstCoEdge,
                    *next = firstCoEdge->nextOnCurve;
    while( next != coedge )
    {
      assert(next != NULL);
      prev = next;
      next = next->nextOnCurve;
    }
  
    prev->nextOnCurve = next->nextOnCurve;
  }
  
  coedge->nextOnCurve = 0;
  coedge->myCurve = 0;
  return CUBIT_SUCCESS;                
}


CubitBox CompositeCurve::bounding_box() const
{ return compGeom->bounding_box(); }

double CompositeCurve::measure() 
{ return compGeom->measure(); }

GeometryType CompositeCurve::geometry_type()
{
  switch (num_curves()) {
    case 0  : return POINT_CURVE_TYPE;
    case 1  : return get_curve(0)->geometry_type();
    default : return UNDEFINED_CURVE_TYPE;
  }
}

GeometryQueryEngine* CompositeCurve::get_geometry_query_engine() const
{ return VirtualQueryEngine::instance(); }


double CompositeCurve::start_param()
{ return 0; }

double CompositeCurve::end_param()
{ return measure(); }

CubitBoolean CompositeCurve::get_param_range( double& lower, double& upper )
{
  lower = 0;
  upper = measure();
  return CUBIT_TRUE;
}

CubitBoolean CompositeCurve::is_periodic( double& period )
{
  if (geometry_type() == POINT_CURVE_TYPE)
    return CUBIT_FALSE;
  
  if( startPoint == endPoint && num_curves() > 0 )
  {
    period = measure();
    return CUBIT_TRUE;
  }
  else
  {
    period = 0;
    return CUBIT_FALSE;
  }
}

CubitStatus CompositeCurve::position_from_u( double u, CubitVector& position )
{
  int index;
  double ui;
  
  if (num_curves() == 0) // point curve
    return CUBIT_FAILURE;
  
  if( ! curve_param( u, ui, index ) || index < 0 )
    return CUBIT_FAILURE;
  
  Curve* curve_ptr = get_curve( index );
  return curve_ptr ? curve_ptr->position_from_u( ui, position ) : CUBIT_FAILURE;
}  

double CompositeCurve::u_from_position( const CubitVector& position )
{
  if (num_curves() == 0) // point curve
    return 0.0;
  
  int index = closest_curve(position );
  double param = get_curve( index )->u_from_position( position );
  return composite_param( index, param );
}

double CompositeCurve::u_from_arc_length( double param, double length )
{ return param + length; }
double CompositeCurve::length_from_u( double u1, double u2 )
{ return u2 - u1; }

CubitStatus CompositeCurve::closest_point( const CubitVector& location,
                                           CubitVector& closest,
                                           CubitVector* tangent_ptr,
                                           CubitVector* curvature_ptr,
                                           double* param )
{
  if (num_curves() == 0) // point curve
  {
    if (!startPoint)
      return CUBIT_FAILURE;
    
    closest = startPoint->coordinates();
    if (tangent_ptr)
      tangent_ptr->set(0.0,0.0,0.0);
    if (curvature_ptr)
      curvature_ptr->set(0.0,0.0,0.0);
    if (param)
      *param = 0.0;
  }
  
  if( compGeom->num_entities() == 1 )
  {
    CubitStatus result = get_curve(0)
      ->closest_point( location, closest, tangent_ptr, curvature_ptr, param );

    if( tangent_ptr && compGeom->sense(0) == CUBIT_REVERSED )
      *tangent_ptr *= -1.0;
    
    if( param )
      *param = composite_param( 0, *param );
      
    return result;
  }
  
  int index = closest_curve( location, &closest );
  Curve* curve = get_curve( index );
  if( tangent_ptr || curvature_ptr || param )
  {
    curve->closest_point( location, closest, tangent_ptr, 
                          curvature_ptr, param );
    if( param )
      *param = composite_param( index, *param );
    if( tangent_ptr && compGeom->sense(index) == CUBIT_REVERSED )
      *tangent_ptr *= -1.0;
  }
  
  return CUBIT_SUCCESS;
}


CubitStatus CompositeCurve::closest_point_trimmed( 
                                    const CubitVector& position,
                                    CubitVector& closest )
{
  if (num_curves() == 0) 
    return closest_point( position, closest );
    
  if( compGeom->num_entities() == 1 )
    return get_curve(0)->closest_point_trimmed( position, closest );
  
  closest_curve( position, &closest );
  return CUBIT_SUCCESS;
}

CubitBoolean CompositeCurve::is_position_on( const CubitVector& position )
{
    // point curve
  if (num_curves() == 0)
  {
    CompositePoint* point = start_point();
    double lensqr = (position - point->coordinates()).length_squared();
    return lensqr <= (GEOMETRY_RESABS*GEOMETRY_RESABS);
  }
  
  if( compGeom->num_entities() == 1 )
    return get_curve(0)->is_position_on( position );
  
  int index = closest_curve( position );
  return get_curve(index)->is_position_on( position );
}

CubitPointContainment CompositeCurve::point_containment( const CubitVector& point )
{
  int index;
  
  if (num_curves() == 0) 
    return is_position_on(point) ? CUBIT_PNT_ON : CUBIT_PNT_OFF;
  
  for( index = 0; index < compGeom->num_entities(); index++ )
  {
    if( get_curve(index)->point_containment( point ) == CUBIT_PNT_ON )
      return CUBIT_PNT_ON;
  }
  
  return CUBIT_PNT_OFF;  //not on any Curve
}

CubitStatus CompositeCurve::get_point_direction( CubitVector& origin,
                                                 CubitVector& direction )
{
  if (num_curves() == 0) // point curve
    return CUBIT_FAILURE;
  
  if( compGeom->num_entities() == 1 )
    return get_curve(0)->get_point_direction( origin, direction );
  
  int count = compGeom->num_entities();
  CubitVector* vect_list = new CubitVector[count];
  double RESABS_SQUARED = CUBIT_RESABS * CUBIT_RESABS;
  direction.set(0,0,0);
  for( int i = 0; i < count; i++ )
  {
    if( ! get_curve(i)->get_point_direction(origin,vect_list[i]) ||
        (vect_list[i].length_squared() < RESABS_SQUARED) )
    {
      delete [] vect_list;
      return CUBIT_FAILURE;
    }
    if( get_sense(i) == CUBIT_REVERSED )
      vect_list[i] *= -1;
    direction += vect_list[i];
  }
  //If we reach this point, then all of the underlying curves are linear.
  //Next check if they are colinear.
  if( direction.length_squared() < RESABS_SQUARED )
  {
    delete [] vect_list;
    return CUBIT_FAILURE;
  }
  CubitVector mean = ~direction;
  for( int j = 0; j < count; j++ )
  {
    if( fabs( 1 - (mean % ~vect_list[j]) ) > CUBIT_RESABS )
    {
      delete [] vect_list;
      return CUBIT_FAILURE;
    }
  }
  
  delete [] vect_list;
  get_curve(0)->get_point_direction(origin,direction);
  direction = mean;
  return CUBIT_SUCCESS;
}

CubitStatus CompositeCurve::get_interior_extrema( 
                                    DLIList<CubitVector*>& interior_points,
                                    CubitSense& return_sense )
{
  return_sense = CUBIT_FORWARD;
  if (num_curves() == 0) // point curve
    return CUBIT_SUCCESS;

    // Go through each curve in the composite
  DLIList<CubitVector*> curve_point_list;
  for (int i = compGeom->num_entities() - 1; i >= 0; i--)
  {
      // Get the next curve's extrema
    Curve* cur_curve = get_curve(i);
    cur_curve->get_interior_extrema(curve_point_list,return_sense);
      // See which order to put them into the return list
    if (return_sense == get_sense(i))
    {
      interior_points += curve_point_list;
    }
    else
    {
      curve_point_list.last();
      for (int j = curve_point_list.size(); j--; )
        interior_points.append(curve_point_list.get_and_back());
    }
      // Unless this is the last curve, put in the point
      // between this and the next curve
    if (i != 0)
    {
      CubitVector* endpoint = new CubitVector(0,0,0);
      if (get_sense(i) == CUBIT_FORWARD)
        cur_curve->position_from_u(cur_curve->end_param(), *endpoint);
      else
        cur_curve->position_from_u(cur_curve->end_param(), *endpoint);
      interior_points.append(endpoint);
    }
      // clean out the list for the next curve
    curve_point_list.clean_out();
  }

  return CUBIT_SUCCESS;
}

CubitStatus CompositeCurve::get_center_radius( CubitVector& c, double& r )
{
  if (num_curves() == 0) // point curve
    return CUBIT_FAILURE;
  
  if (!get_curve(0)->get_center_radius(c,r) )
    return CUBIT_FAILURE;
    
  for (int i = compGeom->num_entities() - 1; i > 0; i-- )
  {
    CubitVector c2;
    double r2;
    if (!get_curve(i)->get_center_radius(c2,r2) ||
        fabs(r - r2) > GEOMETRY_RESABS ||
        (c2 - c).length_squared() > GEOMETRY_RESABS*GEOMETRY_RESABS)
    {
      return CUBIT_FAILURE;
    }
  }
  return CUBIT_SUCCESS;
}
  
CubitBoolean CompositeCurve::G1_discontinuous( double param,
  CubitVector* minus_tangent, CubitVector* plus_tangent )
{
  if (num_curves() == 0) // point curve
    return CUBIT_FALSE;
  
  const double half_resabs = CUBIT_RESABS / 2.0;
  double curve_param;
  int curve_index;
  CubitStatus s = this->curve_param( param, curve_param, curve_index );
  if( ! s )
  {
    PRINT_ERROR("Parameter %f out of range [%f,%f] in "
     "CompositeCurve::G1_discontinous(..)\n",param,start_param(),end_param());
    return CUBIT_FALSE;
  }
  
  CubitVector tan1, tan2, location, cp;
  CubitBoolean result = CUBIT_FALSE;
  if( (curve_index > 0) && 
      (fabs(compGeom->measure(curve_index-1) - param) < half_resabs) )
  {
    get_curve( curve_index )->position_from_u( curve_param, location );
    get_curve( curve_index-1 )->closest_point( location, cp, &tan1 );
    if( get_sense( curve_index-1 ) == CUBIT_REVERSED ) tan1 *= -1.0;
    get_curve( curve_index   )->closest_point( location, cp, &tan2 );
    if( get_sense( curve_index   ) == CUBIT_REVERSED ) tan2 *= -1.0;
    double cross = (tan1 * tan2).length_squared();
    double sum   = (tan1 + tan2).length_squared();
    double diff  = (tan1 - tan2).length_squared();
    if( (cross > CUBIT_RESABS) || (sum < diff) )
    {
      if( minus_tangent ) *minus_tangent = tan1;
      if( plus_tangent ) *plus_tangent = tan2;
      result = CUBIT_TRUE;
    }
    else
      result = CUBIT_FALSE;
    
  }
  else if( (curve_index < (compGeom->num_entities()-1)) &&
           (fabs(compGeom->measure(curve_index) - param) < half_resabs) )
  {
    get_curve( curve_index )->position_from_u( curve_param, location );
    get_curve( curve_index   )->closest_point( location, cp, &tan1 );
    if( get_sense( curve_index   ) == CUBIT_REVERSED ) tan1 *= -1.0;
    get_curve( curve_index+1 )->closest_point( location, cp, &tan2 );
    if( get_sense( curve_index+1 ) == CUBIT_REVERSED ) tan2 *= -1.0;
    double cross = (tan1 * tan2).length_squared();
    double sum   = (tan1 + tan2).length_squared();
    double diff  = (tan1 - tan2).length_squared();
    if( (cross > CUBIT_RESABS) || (sum < diff) )
    {
      if( minus_tangent ) *minus_tangent = tan1;
      if( plus_tangent ) *plus_tangent = tan2;
      result = CUBIT_TRUE;
    }
    else
      result = CUBIT_FALSE;
  }
  else if( get_sense(curve_index) == CUBIT_FORWARD )
  {
    result = get_curve( curve_index )
      ->G1_discontinuous( curve_param, minus_tangent, plus_tangent );
  }
  else
  {
    result = get_curve( curve_index )
      ->G1_discontinuous( curve_param, plus_tangent, minus_tangent );
    if( result )
    {
      if( plus_tangent )  *plus_tangent  *= -1.0;
      if( minus_tangent) *minus_tangent *= -1.0;
    }
  }
  return result;  
}


int CompositeCurve::closest_curve( const CubitVector& location,
                                   CubitVector *point )
{  
  double shortest_distance_sqr, current_distance_sqr;
  CubitVector closest_point, current_point;
  int closest_curve, current_curve;
  
  closest_curve = compGeom->closest_box( location );
  
  get_curve( closest_curve )->closest_point_trimmed( location, closest_point );
  shortest_distance_sqr = (location - closest_point).length_squared();
  
  while( ( current_curve = compGeom->next_box_within_dist( shortest_distance_sqr ) ) >= 0 )
  {
    get_curve( current_curve )->closest_point_trimmed( location, current_point );
    current_distance_sqr = (location - current_point).length_squared();
    if( current_distance_sqr < shortest_distance_sqr )
    {
      closest_curve = current_curve;
      shortest_distance_sqr = current_distance_sqr;
      closest_point = current_point;
    }
  }
    
  if( point ) *point = closest_point;
  return closest_curve;
}



//-------------------------------------------------------------------------
// Purpose       : Parameter Conversion
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck 
//
// Creation Date : 06/08/98
//-------------------------------------------------------------------------
double CompositeCurve::composite_param( int index, double u ) const
{
  if( (index < 0) || (index >= compGeom->num_entities()) ) return -1.0;
  Curve* curve_ptr = get_curve( index );
  double sum = lengthUntilI( index );

  double u_min, u_max;
  curve_ptr->get_param_range( u_min, u_max );

  CubitSense sense = this->get_sense(index);
  double root_param = (sense == CUBIT_FORWARD) ? u_min : u_max;
  double lfu = (root_param < u) ?
    curve_ptr->length_from_u( root_param, u ) :
    curve_ptr->length_from_u( u, root_param );
  return sum + fabs( lfu );
}
  
  
CubitStatus CompositeCurve::curve_param( double uc, double& ui, int& index ) const
{
  double sum = 0.0, period;
  if( const_cast<CompositeCurve*>(this)->is_periodic(period) )
    fixup_periodic_param( uc );

  int max = compGeom->num_entities();
  int min = 0;
  index = max / 2;
  sum = lengthUntilI( max - 1);
  if( uc >= sum ) index = max - 1;
  else while( min != max )
  {
    sum = lengthUntilI(index);
    if( sum > uc )
      max = index;
    else if( lengthUntilI(index+1) > uc ) 
      break;
    else 
      min = index;
    
    index = (min + max) / 2;
  }
      
  if( index >= compGeom->num_entities() ) return CUBIT_FAILURE;

  double start_param, end_param;
  get_curve(index)->get_param_range( start_param, end_param );

  double arc_len = uc - sum;
  if( get_sense(index) == CUBIT_REVERSED ) 
    arc_len = get_curve(index)->measure() - arc_len;

  ui = get_curve(index)->u_from_arc_length( start_param, arc_len );
  return CUBIT_SUCCESS;
}




//-------------------------------------------------------------------------
// Purpose       : Calculate the sum of the lengths of curves 0 to i-1.
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck 
//
// Creation Date : 06/09/98
//-------------------------------------------------------------------------
double CompositeCurve::lengthUntilI( int i ) const
{
  assert( i >= 0 && i < compGeom->num_entities() );
  return ( i == 0 ) ? 0.0 : compGeom->measure(i - 1);
}

//-------------------------------------------------------------------------
// Purpose       : TopologyBridge queries
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 03/04/02
//-------------------------------------------------------------------------
void CompositeCurve::get_parents_virt( DLIList<TopologyBridge*>& list )
{ 
  for( CompositeCoEdge* coedge = firstCoEdge; 
       coedge != 0;
       coedge = coedge->nextOnCurve )
    if (!dynamic_cast<HiddenEntitySet*>(coedge->owner()))
      list.append( coedge );
  
  if (stitchNext)
    stitchNext->get_parents_virt(list);
}
void CompositeCurve::get_children_virt( DLIList<TopologyBridge*>& list )
{
  CompositePoint* sp = start_point();
  CompositePoint* ep = end_point();
  
  if( sp ) 
    list.append( sp );
  if( ep && ep != sp )
    list.append( ep );
}

//-------------------------------------------------------------------------
// Purpose       : TBOwner methods
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 03/04/02
//-------------------------------------------------------------------------
CubitStatus CompositeCurve::remove_bridge( TopologyBridge* bridge )
{
  int i = compGeom->index_of(bridge);
  if ( i >= 0 )
  {
    bridge->owner(0);
    if ( !compGeom->remove(i,true) )
      return CUBIT_FAILURE;

    if (compGeom->num_entities() == 0)
      CompositeEngine::instance().notify_deactivated(this);
    HadBridgeRemoved = 1;
    return CUBIT_SUCCESS;
  }
  else
  {
    for (CompositeCurve* ptr = this; ptr->stitchNext; ptr = ptr->stitchNext)
    {
      if (ptr->stitchNext == bridge)
      {
        ptr->stitchNext = ((CompositeCurve*)bridge)->stitchNext;
HadBridgeRemoved = 1;
        return CUBIT_SUCCESS;
      }
    }
  }
  
  return CUBIT_FAILURE;
}

Curve* CompositeCurve::remove_curve( int index )
{
  Curve* curve = get_curve(index);
  if (!curve || !compGeom->remove(index,false) )
    return 0;
  curve->owner(0);
  return curve;
}

  

CubitStatus CompositeCurve::swap_bridge( TopologyBridge* o,
                                         TopologyBridge* n,
                                         bool reversed )
{
  int i = compGeom->index_of(o);
  GeometryEntity* ge = dynamic_cast<GeometryEntity*>(n);
  if( i >= 0 && ge != 0 )
  {
    o->owner(0);
    n->owner(this);
    if ( ! compGeom->swap( i, ge ) )
      return CUBIT_FAILURE;
    
    if (reversed)
      compGeom->reverse_sense(i);
    return CUBIT_SUCCESS;
  }
  else
    return CUBIT_FAILURE;
}
CubitBoolean CompositeCurve::contains_bridge( TopologyBridge* bridge ) const
{
  if (!compGeom)
    return CUBIT_FALSE;
  
  return (CubitBoolean)(compGeom->index_of(bridge) >= 0);
}


//-------------------------------------------------------------------------
// Purpose       : Attach a CubitSimpleAttribute
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 03/10/98
//-------------------------------------------------------------------------
void CompositeCurve::append_simple_attribute_virt(
                const CubitSimpleAttrib& simple_attrib_ptr )
{
  compGeom->add_attribute( simple_attrib_ptr );
}

//-------------------------------------------------------------------------
// Purpose       : Remove an attached CubitSimpleAttrib
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 03/10/98
//-------------------------------------------------------------------------
void CompositeCurve::remove_simple_attribute_virt(
                const CubitSimpleAttrib& simple_attrib_ptr )
{
  compGeom->rem_attribute( simple_attrib_ptr );
}


//-------------------------------------------------------------------------
// Purpose       : Remove an all attached CubitSimpleAttrib
//
// Special Notes : 
//
// Creator       : Greg Nielson
//
// Creation Date : 07/10/98
//-------------------------------------------------------------------------
void CompositeCurve::remove_all_simple_attribute_virt()
{
  compGeom->rem_all_attributes( );
}


//-------------------------------------------------------------------------
// Purpose       : Return the attached CubitSimpleAttribs.
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 03/10/98
//-------------------------------------------------------------------------
CubitStatus CompositeCurve::get_simple_attribute(
               DLIList<CubitSimpleAttrib>& attrib_list )
{
  compGeom->get_attributes( attrib_list );
  return CUBIT_SUCCESS;
}



//-------------------------------------------------------------------------
// Purpose       : Get attribs by name
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 03/03/03
//-------------------------------------------------------------------------
CubitStatus CompositeCurve::get_simple_attribute(
          const CubitString& name, DLIList<CubitSimpleAttrib>& attrib_list )
{
  compGeom->get_attributes( name.c_str(), attrib_list );
  return CUBIT_SUCCESS;
}

void CompositeCurve::fixup_periodic_param( double& param ) const
{
  double period = compGeom->measure();
  if( param < 0.0 )
    param = period - param;
  if( param >= period )
    param = fmod( param, period );
}

//-------------------------------------------------------------------------
// Purpose       : Write attributes to underlying geometry entity
//
// Special Notes : Special case for point-curves
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 02/28/04
//-------------------------------------------------------------------------
void CompositeCurve::write_attributes()
{
  if (num_curves()) // not a point-curve
  {
    compGeom->write_attributes();
  }
  else
  {
    CompositePoint* point = start_point();
    assert(point == end_point());
    compGeom->write_attributes(point);
  }
}

//-------------------------------------------------------------------------
// Purpose       : Read attributes from underlying geometry
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 02/28/04
//-------------------------------------------------------------------------
void CompositeCurve::read_attributes()
{
  if (num_curves()) // not a point-curve
  {
    compGeom->read_attributes();
  }
  else
  {
    CompositePoint* point = start_point();
    assert(point == end_point());
    compGeom->read_attributes(point);
  }
}


void CompositeCurve::get_hidden_points( DLIList<TBPoint*>& points )
{
  if( hiddenSet )
  {
    hiddenSet->hidden_points( points );
  }
}



void CompositeCurve::print_debug_info( const char* prefix, 
                                       bool brief ) const
{
  if( prefix == 0 ) prefix = "";
  
  if ( brief )
  {
#ifdef TOPOLOGY_BRIDGE_IDS
    PRINT_INFO("%sCompositeCurve %d : points (%d,%d) ", prefix, get_id(),
      startPoint ? startPoint->get_id() : -1, endPoint ? endPoint->get_id() : -1 );
    if ( num_curves() == 1 )
      PRINT_INFO(": %s %d\n", fix_type_name(typeid(*get_curve(0)).name()), get_curve(0)->get_id());
    else
      PRINT_INFO(": %d curves.\n", num_curves() );
#else
    PRINT_INFO("%sCompositeCurve %p : points (%p,%p) ", prefix, this,
      startPoint, endPoint );
    if ( num_curves() == 1 )
      PRINT_INFO(": %s %p\n", fix_type_name(typeid(*get_curve(0)).name()), get_curve(0));
    else
      PRINT_INFO(": %d curves.\n", num_curves() );
#endif
    return;
  }
    
  char* new_prefix = new char[strlen(prefix)+3];
  strcpy( new_prefix, prefix );
  strcat( new_prefix, "  ");
  
#ifdef TOPOLOGY_BRIDGE_IDS
  PRINT_INFO("%sCompositeCurve %d\n", prefix, get_id() );
#else
  PRINT_INFO("%sCompositeCurve %p\n", prefix, this );
#endif
  compGeom->print_debug_info( new_prefix );
  if( hiddenSet ) hiddenSet->print_debug_info( new_prefix );
  else PRINT_INFO("%s  No Hidden Entities.\n", prefix );
  
  CompositeCoEdge* coedge = first_coedge();
  DLIList<TopologyBridge*> loops(1), surface_list(1);
  
  while( coedge )
  {
    TopologyBridge* surf_ptr = 0;
    loops.clean_out();
    surface_list.clean_out();
    coedge->get_parents_virt(loops);
    if ( loops.size() == 1 )
    {
      loops.get()->get_parents_virt(surface_list);
      if( surface_list.size() == 1 )
        surf_ptr = surface_list.get();
    }

#ifdef TOPOLOGY_BRIDGE_IDS
    PRINT_INFO("%s  %s on Surface %d\n", prefix, 
      coedge->sense() == CUBIT_FORWARD ? "FORWARD" :
      coedge->sense() == CUBIT_REVERSED ? "REVERSE" : "UNKNOWN",
      surf_ptr->get_id() );
#else
    PRINT_INFO("%s  %s on Surface %p\n", prefix, 
      coedge->sense() == CUBIT_FORWARD ? "FORWARD" :
      coedge->sense() == CUBIT_REVERSED ? "REVERSE" : "UNKNOWN",
      surf_ptr );
#endif

    if (!coedge->get_loop())
      coedge->print_debug_info( new_prefix, true );
    coedge = next_coedge(coedge);
  }
  
  if ( !start_point() )
    PRINT_INFO("%s  NULL START POINT\n", prefix );
  else
    start_point()->print_debug_info( new_prefix, true );
  
  if ( !end_point() )
    PRINT_INFO("%s  NULL END POINT\n", prefix );
  else if ( end_point() == start_point() )
    PRINT_INFO("%s  end point SAME as start point\n", prefix ); 
  else
    end_point()->print_debug_info( new_prefix, true );
  
  delete [] new_prefix;
}


//-------------------------------------------------------------------------
// Purpose       : Update composite for change in sense of underlying curve
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 
//-------------------------------------------------------------------------
void CompositeCurve::notify_reversed( TopologyBridge* bridge )
{
  int index = compGeom->index_of(bridge);
  if( index >= 0 )
    compGeom->reverse_sense(index);
}

//-------------------------------------------------------------------------
// Purpose       : Merge
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 12/02/02
//-------------------------------------------------------------------------
CubitStatus CompositeCurve::stitch( CompositeCurve* dead )
{
    // check to make sure end points are already stitched.
  bool reversed;
  
  if( dead->start_point() == this->start_point() )
    reversed = false;
  else if( dead->end_point() == this->start_point() )
    reversed = true;
  else
    return CUBIT_FAILURE;
  
  if( ( reversed && (dead->start_point() != this->end_point())) ||
      (!reversed && (dead->end_point()   != this->end_point())) )
    return CUBIT_FAILURE;

    // shouldn't be merging hidden curves
  if( dynamic_cast<CompositeCurve*>(this->owner()) ||
      dynamic_cast<CompositeCurve*>(dead->owner()) )
  {
    assert(0);
    return CUBIT_FAILURE;
  }
  
    // DAG / next level up should already be merged
  if( dead->owner() != this->owner() )
  {
    assert(0);
    return CUBIT_FAILURE;
  }
  
    // update owner 
  if( dead->owner() )
    dead->owner()->notify_merged( dead, this );
  assert( !dead->owner() );
  dead->owner( this );
  
    // merge lists
  CompositeCurve* end = dead;
  while (end->stitchNext) 
  {
    end->stitchNext->owner( this );
    end = end->stitchNext;
  }
  end->stitchNext = stitchNext;
  stitchNext = end;
  
  return CUBIT_SUCCESS;
}
      
//-------------------------------------------------------------------------
// Purpose       : Get the visible curve of a set of stitched curves
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 12/02/02
//-------------------------------------------------------------------------
CompositeCurve* CompositeCurve::primary_stitched_curve()
{
  CompositeCurve* result = dynamic_cast<CompositeCurve*>(owner());
  return result ? result : this;
}

//-------------------------------------------------------------------------
// Purpose       : Is this curve stitched with any others
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 12/02/02
//-------------------------------------------------------------------------
bool CompositeCurve::is_stitched()
{
  return stitchNext || dynamic_cast<CompositeCurve*>(owner());
}

//-------------------------------------------------------------------------
// Purpose       : Get list of stitched curves
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 12/16/02
//-------------------------------------------------------------------------
void CompositeCurve::get_stitched( DLIList<CompositeCurve*>& list )
{
  for (CompositeCurve* curve = primary_stitched_curve();
       curve; curve = curve->stitchNext)
    list.append( curve );
}

//-------------------------------------------------------------------------
// Purpose       : Remove all stitched curves
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 12/02/02
//-------------------------------------------------------------------------
void CompositeCurve::unstitch_all()
{
    // should only call on visible curve
  assert (this == primary_stitched_curve());
  
  while (stitchNext)
  {
    stitchNext->owner(0);
    if (owner())
      owner()->notify_copied( stitchNext, this );
    stitchNext = stitchNext->stitchNext;
  }
}

  
//-------------------------------------------------------------------------
// Purpose       : Update for split in underlying curve
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 12/19/02
//-------------------------------------------------------------------------
void CompositeCurve::notify_split( TopologyBridge* new_bridge,
                                   TopologyBridge* old_bridge )
{
  Curve* old_curve = dynamic_cast<Curve*>(old_bridge);
  Curve* new_curve = dynamic_cast<Curve*>(new_bridge);
  assert( old_curve && new_curve && index_of(old_curve) >= 0 );
  
    // get start and end points
  DLIList<TopologyBridge*> bridges;
  old_curve->get_children(bridges,true,old_curve->layer());
  bridges.reset();
  TBPoint* old_start = dynamic_cast<TBPoint*>(bridges.get());
  TBPoint* old_end   = dynamic_cast<TBPoint*>(bridges.next());
  bridges.clean_out();
  new_curve->get_children(bridges,true,new_curve->layer());
  TBPoint* new_start = dynamic_cast<TBPoint*>(bridges.get());
  TBPoint* new_end   = dynamic_cast<TBPoint*>(bridges.next());
  bridges.clean_out();

  
    // find new point
  bool sp = (new_start == old_start || new_start == old_end);
  bool ep = (new_end   == old_start || new_end   == old_end);
  assert( ep != sp ); // one must be true and one must be false
  TBPoint* new_pt = sp ? new_start : new_end;
  
    // find relative sense of curves
  int old_index = index_of(old_curve);
  bool is_start = (old_start == new_pt);
  bool is_forward = (get_sense(old_index) == CUBIT_FORWARD);
  bool prepend = (is_start == is_forward);
  is_start = (new_start == new_pt);
  assert( is_start || new_pt == new_end );
  bool reversed = (is_start == prepend);
  
  CompositePoint* new_cpt = new CompositePoint(new_pt);
  
  hidden_entities().hide(new_cpt);
  
  CubitSense new_sense = reversed ? CUBIT_REVERSED : CUBIT_FORWARD;
  int insert_index = old_index;
  if( !prepend ) insert_index++;
  if( ! compGeom->insert( insert_index, new_curve, new_sense ) )
  {
    assert(0);
    delete new_cpt;
    return;
  }
  new_curve->owner(this);
  
  DLIList<CoEdgeSM*> new_coedges, old_coedges;
  bridges.clean_out();
  new_curve->get_parents_virt( bridges );
  CAST_LIST( bridges, new_coedges, CoEdgeSM );
  assert(bridges.size() == new_coedges.size());
  bridges.clean_out();
  old_curve->get_parents_virt( bridges );
  CAST_LIST( bridges, old_coedges, CoEdgeSM );
  assert(bridges.size() == old_coedges.size());
  
  for( CompositeCoEdge* coedge = first_coedge();
       coedge;
       coedge = next_coedge(coedge) )
  {
      // for each coedge of the old curve
    int i;
    CoEdgeSM* old_coedge = 0;
    for( i = 0; i < coedge->num_coedges(); i++ )
    {
      old_coedge = coedge->get_coedge(i);
      if( old_coedges.is_in_list(old_coedge) )
        break;
    }
    assert(i < coedge->num_coedges() && old_coedge);
    
    bridges.clean_out();
    old_coedge->get_parents_virt(bridges);
    assert(bridges.size() == 1);
    LoopSM* loopsm = dynamic_cast<LoopSM*>(bridges.get());
    assert(0 != loopsm);
      
    bridges.clean_out();
    loopsm->get_children(bridges, true, old_curve->layer());
    bridges.move_to(old_coedge);
    assert(bridges.get() == old_coedge);
      
      // Determine if new_coedge (the one we are looking for)
      // should occur before or after old_coedge in the loop.
      // If new_curve is after old_curve (old_end == new_pt), 
      // then the we want the coedge after old_coedge iff the
      // sense of old_coedge is forward, otherwise the one before.
      // Invert that if new_curve is before old_curve (old_start
      // == new_pt).
    bool coe_reversed = (old_coedge->sense() == CUBIT_REVERSED);
    bool curve_prepend = (old_start == new_pt);
    bool previous = coe_reversed != curve_prepend;
    CoEdgeSM* new_coedge = 
      dynamic_cast<CoEdgeSM*>(previous ? bridges.prev() : bridges.next());
    assert(new_coedges.is_in_list(new_coedge));
    
    CubitStatus s = coedge->insert_coedge(insert_index, new_coedge);
    assert(s);
  }
}

CubitStatus CompositeCurve::get_spline_params
(
  bool &rational,    // return true/false
  int &degree,       // the degree of this spline
  DLIList<CubitVector> &cntrl_pts,  // xyz position of controlpoints
  DLIList<double> &cntrl_pt_weights, // if rational, a weight for each cntrl point.
  DLIList<double> &knots,   // There should be order+cntrl_pts.size()-2 knots
  bool &spline_is_reversed
) const
{
  PRINT_ERROR("Currently, Cubit is unable to determine spline parameters for CompositeCurves.\n");
  return CUBIT_FAILURE;
}

CubitStatus CompositeCurve::get_ellipse_params
(
  CubitVector &center_vec,
  CubitVector &normal,
  CubitVector &major_axis,
  double &radius_ratio
) const
{
  PRINT_ERROR("Currently, Cubit is unable to determine ellipse parameters for CompositeCurves.\n");
  return CUBIT_FAILURE;
}

/*
void CompositeCurve::draw( int color )
{
  int num_pts;
  GMem gmem;
  
  if (!VirtualQueryEngine::instance()->get_graphics(this, num_pts, &gmem))
    return;
  
    GfxDebug::draw_polyline( gmem.point_list(),
                             gmem.pointListCount,
                             color );
    GfxDebug::flush();
}
*/
