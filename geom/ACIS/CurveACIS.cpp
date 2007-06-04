//-------------------------------------------------------------------------
// Filename      : CurveACIS.cpp
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

// ********** BEGIN STANDARD INCLUDES      **********
#include <assert.h>
// ********** END STANDARD INCLUDES        **********

// ********** BEGIN ACIS INCLUDES          **********
#if CUBIT_ACIS_VERSION < 1100
#include "kernel/acis.hxx"
#include "kernel/kernapi/api/api.hxx"
#include "kernel/kernapi/api/kernapi.hxx"
#include "intersct/kernapi/api/intrapi.hxx"
#include "kernel/kerndata/data/entity.hxx"
#include "kernel/kerndata/top/edge.hxx"
#include "kernel/kerndata/top/vertex.hxx"
#include "kernel/kerndata/geom/point.hxx"
#include "kernel/kerndata/geom/allcurve.hxx"
#include "kernel/kerndata/lists/lists.hxx"
#include "kernel/kerndata/top/alltop.hxx" // Required for is_TEDGE
#include "kernel/kernint/d3_chk/chk_stat.hxx"
#include "intersct/sg_husk/query/sgquery.hxx"
#include "intersct/sg_husk/query/sgquertn.hxx"
#include "baseutil/vector/position.hxx"
#include "baseutil/vector/interval.hxx"
#else
#include "acis.hxx"
#include "api.hxx" 
#include "kernapi.hxx" 
#include "intrapi.hxx"
#include "entity.hxx" 
#include "edge.hxx"
#include "vertex.hxx"
#include "point.hxx"
#include "allcurve.hxx"
#include "lists.hxx"
#include "alltop.hxx" // Required for is_TEDGE
#include "chk_stat.hxx"
#include "sgquery.hxx"
#include "sgquertn.hxx"
#include "position.hxx"
#include "interval.hxx"
#endif

/*
//Added for G1_discontinous(..)
#include "acistype.hxx"
#include "sps3crtn.hxx"
//This declaration is incorrect in the above header for ACIS 5.3.  
//The following declaration agrees with the ACIS 5.3 documentation, 
//and agrees with libspline.
DECL_SPLINE void
bs3_curve_smoothness( bs3_curve, double, logical&, logical& );
*/

// ********** END ACIS INCLUDES            **********

// ********** BEGIN CUBIT INCLUDES         **********
#include "CubitDefines.h"
#include "GeometryDefines.h"
#include "GeometryQueryTool.hpp"
#include "CubitVector.hpp"

#include "AcisQueryEngine.hpp"

#include "RefEdge.hpp"
#include "RefVolume.hpp"
#include "CurveACIS.hpp"

#include "CoEdgeSM.hpp"
#include "Point.hpp"

#include "CastTo.hpp"
#include "DLIList.hpp"
#include "attrib_cubit_owner.hpp"

// ********** END CUBIT INCLUDES           **********

// ********** BEGIN FORWARD DECLARATIONS   **********
// ********** END FORWARD DECLARATIONS     **********

// ********** BEGIN STATIC DECLARATIONS    **********
// ********** END STATIC DECLARATIONS      **********

// ********** BEGIN PUBLIC FUNCTIONS       **********

//-------------------------------------------------------------------------
// Purpose       : The constructor with a pointer to the first EDGE.
//
// Special Notes :
//
// Creator       : Xuechen Liu
//
// Creation Date : 08/02/96
//-------------------------------------------------------------------------
CurveACIS::CurveACIS(EDGE* EDGE_ptr)
    : AcisBridge(EDGE_ptr)
{
    // Calculate a bounding box if there isn't one already
  get_acis_query_engine()->bounding_box(EDGE_ptr);
}

//-------------------------------------------------------------------------
// Purpose       : The destructor. 
//
// Special Notes :
//
// Creator       : Raikanta Sahu
//
// Creation Date : 09/06/96
//-------------------------------------------------------------------------
CurveACIS::~CurveACIS() 
{
}

//-------------------------------------------------------------------------
// Purpose       : The purpose of this function is to append a
//                 attribute to the GE. The name is attached to the 
//                 underlying solid model entity this one points to.
//
//
// Special Notes : 
//
// Creator       : Malcolm J. Panthaki
//
// Creation Date : 11/21/96
//-------------------------------------------------------------------------
void CurveACIS::append_simple_attribute_virt(CubitSimpleAttrib* csattrib_ptr)
{
  AcisBridge::append_simple_attribute_virt(csattrib_ptr);
}

//-------------------------------------------------------------------------
// Purpose       : The purpose of this function is to remove a simple 
//                 attribute attached to this geometry entity. The name is 
//                 removed from the underlying BODY this points to.
//
// Special Notes : 
//
// Creator       : David R. White
//
// Creation Date : 03/18/97
//-------------------------------------------------------------------------
void CurveACIS::remove_simple_attribute_virt(CubitSimpleAttrib* csattrib_ptr)
{
  AcisBridge::remove_simple_attribute_virt(csattrib_ptr);
}

//-------------------------------------------------------------------------
// Purpose       : The purpose of this function is to remove all simple 
//                 attributes attached to this geometry entity.  Also
//                 removes lingering GTC attributes.
//
//
// Special Notes : 
//
// Creator       : Greg Nielson
//
// Creation Date : 07/10/98
//-------------------------------------------------------------------------
void CurveACIS::remove_all_simple_attribute_virt()
{
  AcisBridge::remove_all_simple_attribute_virt();
}

//-------------------------------------------------------------------------
// Purpose       : The purpose of this function is to get the  
//                 attributes attached to this geometry entity. The name is 
//                 attached to the underlying BODY this points to.
//
// Special Notes : 
//
// Creator       : Malcolm J. Panthaki
//
// Creation Date : 01/23/97
//-------------------------------------------------------------------------
CubitStatus CurveACIS::get_simple_attribute(DLIList<CubitSimpleAttrib*>&
                                            cubit_simple_attrib_list)
{
  return AcisBridge::get_simple_attribute(cubit_simple_attrib_list);
}
CubitStatus CurveACIS::get_simple_attribute(const CubitString& name,
                                       DLIList<CubitSimpleAttrib*>& list)
  { return AcisBridge::get_simple_attribute(name, list); }

//-------------------------------------------------------------------------
// Purpose       : Get geometry modeling engine: AcisQueryEngine
//
// Special Notes :
//
// Creator       : jihong Ma
//
// Creation Date : 10/22/96
//-------------------------------------------------------------------------
GeometryQueryEngine* 
CurveACIS::get_geometry_query_engine() const
{
   return get_acis_query_engine();   
}                 

//-------------------------------------------------------------------------
// Purpose       : Get the bounding box of the object.
//
// Special Notes :
//
// Creator       : Raikanta Sahu
//
// Creation Date : 10/23/96
//-------------------------------------------------------------------------
CubitBox CurveACIS::bounding_box() const 
{
    // Calculate a bounding box if there isn't one already and return
    // the result.
  SPAbox ACIS_box = get_acis_query_engine()->bounding_box(get_EDGE_ptr());
  
    // Convert to a CubitBox and return it
  return get_acis_query_engine()->bounding_box(ACIS_box);
}

/*  
//-------------------------------------------------------------------------
// Purpose       : merges "this" with input GeometryEntity
//
// Special Notes :
//
// Creator       : Jihong Ma
//
// Creation Date : 11/07/96
//-------------------------------------------------------------------------
CubitStatus CurveACIS::merge(GeometryEntity* )
{ return CUBIT_FAILURE; }

TopologyEntity* CurveACIS::unmerge( DLIList<RefVolume*>)
{
  PRINT_ERROR("Unmerging of Acis-based entities currently disabled.\n");
  return NULL;
}
*/
//-------------------------------------------------------------------------
// Purpose       : Return the length of the curve.
//
// Special Notes :
//
// Creator       : Raikanta Sahu
//
// Creation Date : 12/19/96
//-------------------------------------------------------------------------
double CurveACIS::measure()
{
    // Get the first EDGE from the list of EDGES
  EDGE const* EDGE_ptr = this->get_EDGE_ptr();
   
    // First get the start and end parameter values for this ACIS EDGE
  if ( EDGE_ptr == NULL || EDGE_ptr->geometry() == NULL )
  {
    return 0.0;
  }
  double t_low, t_high;
  
  if ( EDGE_ptr -> sense() == REVERSED )
  {
    t_low  = - (EDGE_ptr->end_param());
    t_high = - (EDGE_ptr->start_param());
  }
  
  else
  {
    t_low  = EDGE_ptr->start_param();
    t_high = EDGE_ptr->end_param();
  }

    // Although t_low should automatically be less than t_high,
    // let's be sure, and attempt to fix it if it's periodic
  if ( t_high <= t_low )
  {
    t_low += EDGE_ptr->geometry()->equation().param_period();
    if ( t_high <= t_low )
    {
      PRINT_ERROR("In CurveACIS::length\n"
                  "       Couldn't resolve parameterization for "
                  "the input EDGE\n");
      return 0.0;
    }
  }
    // Get the ACIS curve associated with this EDGE
  return EDGE_ptr->geometry()->equation().length(t_low, t_high);
}

//-------------------------------------------------------------------------
// Purpose       : Return the arc length along the Curve starting from
//                 the point represented by the parameter1 going to the 
//                 point represented by parameter2.
//
// Special Notes : The sign of the returned length value is always positive.
//                 Parameter1 and parameter2 are with respect to the EDGE.
//
// Creator       : Malcolm J. Panthaki
//
// Creation Date : 2/26/97
//-------------------------------------------------------------------------
double CurveACIS::length_from_u( double parameter1, double parameter2 )
{
     // Get the ACIS curve associated with the first EDGE
   curve const* curve_ptr = this->get_ACIS_curve();
   if ( curve_ptr == NULL )
      return 0.0;
   
   EDGE const* EDGE_ptr = this->get_EDGE_ptr();
   if ( EDGE_ptr->sense() == REVERSED )
   {
      parameter1 = -(parameter1);
      parameter2 = -(parameter2);
   }
     // Get the arc length along this curve
   return fabs( curve_ptr->length ( parameter1, parameter2 ) );   
}

//-------------------------------------------------------------------------
// Purpose       : Returns CUBIT_TRUE and the associated period value, if 
//                 the ACIS curve associated with the first EDGE is periodic.
//                 Otherwise returns CUBIT_FALSE and a value of 0.0 for
//                 the period.
//
// Special Notes :  
//
// Creator       : Malcolm J. Panthaki
//
// Creation Date : 2/26/97
//-------------------------------------------------------------------------
CubitBoolean CurveACIS::is_periodic(double& period)
{
     // Get the ACIS curve associated with the first EDGE
   curve const* curve_ptr = this->get_ACIS_curve();
   if ( curve_ptr == NULL )
       return CUBIT_FALSE;
   
     // Find out if it is periodic and if so, return the period as well
   if ( (curve_ptr->periodic()) == TRUE )
   {
      period = curve_ptr->param_period();
      return CUBIT_TRUE;
   }
   
   else
   {
      period = 0.0;
      return CUBIT_FALSE;
   }
}

//------------------------------------------------------------------
// Purpose: Returns CUBIT_TRUE and the associated parametric values, 
//          if the ACIS curve associated with the first EDGE is 
//          parametric.
//          Otherwise returns CUBIT_FALSE and the values of 
//          the lower and upper parametric bounds are undetermined.
//
// Special Notes:
//   *All* ACIS curves are considered to be parametric,
//   even those that have an analytic definition (such as
//   a straight line).  Hence, this function always returns
//   CUBIT_TRUE and values for the upper and lower bounds.
//
//   Note that the lower bound is the parameter value of the start
//   location of the RefEdge that uses this Curve and the upper
//   bound is the parameter value of the end location of the RefEdge
//   that uses this Curve.  This takes into account the sense of the
//   RefEdge with respect to the Curve (which could be REVERSED).
//   Hence, the numerical value of the lower parameter bound could be
//   greater than that of the upper parameter bound.
//
// Creator       : Malcolm J. Panthaki
//
// Creation Date : 2/26/97
//-------------------------------------------------------------------
CubitBoolean CurveACIS::get_param_range( double& lower_bound,
                                         double& upper_bound )
{
    // Get the EDGE
  EDGE const* EDGE_ptr = this->get_EDGE_ptr();
  
    // Get the start and end parameters of this EDGE
  SPAparameter start_param = EDGE_ptr->start_param();
  SPAparameter end_param = EDGE_ptr->end_param();
  
    // Make sure that we're taking into account the sense of the
    // RefEdge with respect to this curve
  lower_bound = (double)start_param;
  upper_bound = (double)end_param;

    // All ACIS curves are parametrically defined so...
  return CUBIT_TRUE;
}

// Finds the extrema along this RefEdge.  It is the responsibility of the
// calling code to delete the CubitVectors added to interior_points!
CubitStatus CurveACIS::get_interior_extrema(DLIList<CubitVector*>& interior_points,
                                            CubitSense& return_sense)
{
    // The points go into the list in order of increasing curve parameter,
    // NOT increasing EDGE parameter. So, the return sense is the same as
    // the sense of the RefEdge wrt to the curve (not EDGE).
  EDGE* edge = get_EDGE_ptr();
  return_sense = get_EDGE_sense();
  
    // Hardpoints have NULL geometry, no internal extrema
  if (edge->geometry() == NULL)
    return CUBIT_SUCCESS;
  
  const curve* acis_curve = &edge->geometry()->equation();
  
    // If it's a straight curve, no internal extrema
  if (acis_curve->type() == straight_type)
    return CUBIT_SUCCESS;

    // Get the parameter interval and period, accounting for sense
  SPAinterval edge_range;
  API_BEGIN;
  edge_range = edge->param_range();
  API_END;
  
  if (return_sense == CUBIT_REVERSED)
    edge_range.negate();
  double period = acis_curve->param_period();
  DLIList<double> point_param_list;
 
    // Prepare to check for extrema in 3 principal directions
  SPAunit_vector n(1.0, 0.0, 0.0);
  for (int i = 0; i < 3; i++) 
  {
      // Set which direction we're looking in
    n.set_component((i+2)%3, 0.0);
    n.set_component(i, 1.0);
    
      // Look for the extrema
    curve_extremum* extremum = acis_curve->find_extrema(n);
    while (extremum) 
    {
      double extreme_param = extremum->param;
      
        // Map the parameter into the right range.  Points that don't lie
        // in the edge_range will end up lower than edge_range.
        // Acis should handle this correctly already but it doesn't!!
      if (period != 0.) 
      {
        assert (period > 0);
        while ( extreme_param < edge_range )
          extreme_param += period;
        while ( extreme_param > edge_range )
          extreme_param -= period;
      }
      
        // if extremum is within the edge, save it               
      if (edge_range >> extreme_param)
        point_param_list.append(extreme_param);
      
        // go to the next extremum
      curve_extremum* prev_extremum = extremum;
      extremum = extremum->next;
      delete prev_extremum;
    }
  }
  
    // Now take the parameters we've collected
    // and convert them into CubitVectors
  if (point_param_list.size())
  {
      // The points should be at least this far apart to be worth putting in the list
    const double epsilon = 30.* GEOMETRY_RESABS;
    const double epsilon_squared = epsilon*epsilon;
    
    SPAposition start = acis_curve->eval_position(edge->start_param());
    SPAposition end = acis_curve->eval_position(edge->end_param());
    int j;
    CubitVector* cubit_position = NULL;
    
    point_param_list.sort();
    point_param_list.reset();
    for (j = point_param_list.size(); j--; ) 
    {
      SPAposition extremum_position =
        acis_curve->eval_position(point_param_list.get_and_step());
      
        // save if not equal to an endpoint, or prior point
      if (!((extremum_position - start).is_zero(epsilon)) &&
          !((extremum_position - end).is_zero(epsilon)))
      {
        CubitVector temp_position( extremum_position.x(),
                                   extremum_position.y(),
                                   extremum_position.z() );
        if (!cubit_position ||
            ((temp_position - *cubit_position).length_squared() > epsilon_squared))
        {
          cubit_position = new CubitVector( temp_position );
          interior_points.append( cubit_position );
        } // If point isn't close to previous point
      } // If point isn't at an endpoint
    } // for each point
  } // End of conversion to CubitVectors
  
  return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Purpose       : This function computes the point on the curve closest 
//                 to the input location.  Optionally, it can also compute
//                 the tangent and curvature on the Curve at the point on
//                 on the Curve closest to the input location.
//
// Special Notes : The tangent direction is always in the positive direction of the 
//                 owning RefEdge, regardless of the positive direction of the
//                 underlying solid model entities.
//
//                 If the calling code needs the tangent and/or the curvature,
//                 it is responsible for allocating the memory for these
//                 CubitVector(s) and sending in the relevant non-NULL
//                 pointers to this routine.
//
// Creator       : Malcolm J. Panthaki
//
// Creation Date : 2/25/97
//-------------------------------------------------------------------------
CubitStatus CurveACIS::closest_point( CubitVector const& location, 
                                      CubitVector& closest_location,
                                      CubitVector* tangent_ptr,
                                      CubitVector* curvature_ptr,
                                      double* param)
{
    // Get the ACIS curve associated with the first EDGE in this Curve and
    // call the appropriate (ACIS) point_perp function that calculates the 
    // new location (on the curve) as well as any other requested data
    // such as the tangent and/or the curvature
  curve const* curve_ptr = this->get_ACIS_curve();
  
    // Special case : point-curve
  if (curve_ptr == NULL )
  {
    if ( tangent_ptr != NULL )
      tangent_ptr->set(0.0,0.0,0.0);
    if ( curvature_ptr != NULL )
      curvature_ptr->set(0.0,0.0,0.0);
    position_from_u( start_param(), closest_location);
    if(param != NULL)
       *param = start_param();
    
    return CUBIT_SUCCESS;
  }

    // Evaluate to ACIS geometry
  SPAposition acis_point ( location.x(), location.y(), location.z() );
  SPAposition result_pos;
  SPAparameter result_param;

    // Need curvature?  If so, must get tangent also.
  if (curvature_ptr)
  {
    SPAunit_vector tangent;
    SPAvector curvature;
    curve_ptr->point_perp( acis_point, 
                           result_pos, 
                           tangent, 
                           curvature,
                           *(SPAparameter*)NULL_REF,
                           result_param );
    curvature_ptr->set( curvature.x(), curvature.y(), curvature.z() );
      // Wantend tangent?
    if (tangent_ptr)
      tangent_ptr->set( tangent.x(), tangent.y(), tangent.z() );
  }
    // Tangent but not curvature
  else if(tangent_ptr)
  {
    SPAunit_vector tangent;
    curve_ptr->point_perp( acis_point,
                           result_pos,
                           tangent,
                           *(SPAparameter*)NULL_REF,
                           result_param );
    tangent_ptr->set( tangent.x(), tangent.y(), tangent.z() );
  }
    // Neither tangent nor curvature
  else
  {
    curve_ptr->point_perp( acis_point, 
                           result_pos,
                           *(SPAparameter*)NULL_REF,
                           result_param );
  }

    // Copy result position to output argument
  closest_location.set( result_pos.x(), result_pos.y(), result_pos.z() );

    // Pass back parameter of closest_location
  if (param) 
  {
      // Adjust for EDGE sense
    if (get_EDGE_sense() == CUBIT_FORWARD)
      *param = result_param;
    else
      *param = -result_param;

    adjust_periodic_parameter(*param);
  }

    // Update passed-back tangent and curvature pointers for
    // relative sense.
  if (get_EDGE_sense() == CUBIT_REVERSED)
  {
    if (tangent_ptr)
      *tangent_ptr = -*tangent_ptr;
    if (curvature_ptr)
      *curvature_ptr = -*curvature_ptr;
  }

  return CUBIT_SUCCESS;
}
/*
//-------------------------------------------------------------------------
// Purpose       : Evaluate closest point to bounded EDGE
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 01/27/04
//-------------------------------------------------------------------------
CubitStatus CurveACIS::closest_point_trimmed( CubitVector const& from_pt,
                                              CubitVector& result )
{
  curve const* curve_ptr = this->get_ACIS_curve();
  if (!curve_ptr)
    return position_from_u( start_param(), result );
  
  SPAposition input_pos ( from_pt.x(), from_pt.y(), from_pt.z() );
  SPAposition result_pos, pos;
  EDGE const* edge_ptr = get_EDGE_ptr();
  SPAinterval range = edge_ptr->param_range();
  SPAparameter param, param2, period = curve_ptr->param_period();
  
    // First try regular point-perp
  curve_ptr->point_perp( input_pos, result_pos, *(SPAparameter*)NULL_REF, param );
  
    // Adjust param if periodic curve
  if (period != 0.0)
  {
    while (param < range)
      param += period;
    while (param > range)
      param -= period;
  }
  
    // If result is on bounded edge, done
  if ( range >> param )
  {
    result.set( result_pos.x(), result_pos.y(), result_pos.z() );
    return CUBIT_SUCCESS;
  }
  
    // Unless curve is a line, look for other potential results
  double dist_sqr = CUBIT_DBL_MAX;
  if (curve_ptr->type() != straight_type)
  {
    for (int axis = 0; axis < 3; axis++ )
    {
      SPAunit_vector direction( 0, 0, 0 );
      direction.set_component( axis, 1 );
      curve_extremum* extrema = curve_ptr->find_extrema(direction);

      while( extrema )
      {
          // get paramter value of extrema, and advance pointer
        param2 = extrema->param;
        curve_extremum* dead = extrema;
        extrema = extrema->next;
        delete dead;

        curve_ptr->point_perp( input_pos, pos, param2, param );

          // Adjust param if periodic curve
        if (period != 0.0)
        {
          while (param < range)
            param += period;
          while (param > range)
            param -= period;
        }

          // If not within range, skip it
        if ( !(range >> param) )
          continue;

        SPAvector diff( pos.x() - input_pos.x(),
                        pos.y() - input_pos.y(),
                        pos.z() - input_pos.z() );

        if (diff.len_sq() < dist_sqr)
        {
          dist_sqr = diff.len_sq();
          result_pos = pos;
        }
      }
    }
  }
  
    // Check start and end positions
  pos = curve_ptr->eval_position( edge_ptr->start_param() );
  SPAvector start_diff( pos.x() - input_pos.x(),
                        pos.y() - input_pos.y(),
                        pos.z() - input_pos.z() );
  if (start_diff.len_sq() < dist_sqr)
  {
    dist_sqr = start_diff.len_sq();
    result_pos = pos;
  }
  
  pos = curve_ptr->eval_position( edge_ptr->end_param() );
  SPAvector end_diff( pos.x() - input_pos.x(),
                      pos.y() - input_pos.y(),
                      pos.z() - input_pos.z() );
  if (end_diff.len_sq() < dist_sqr)
  {
    dist_sqr = end_diff.len_sq();
    result_pos = pos;
  }
  
  if (dist_sqr == CUBIT_DBL_MAX)
    return CUBIT_FAILURE;
  
  result.set( result_pos.x(), result_pos.y(), result_pos.z() );
  return CUBIT_SUCCESS;
}
*/
//------------------------------------------------------------------
// Purpose       : Determine if given location is on bounds of curve
//                 
//
// Special Notes :
//
// Creator       : Steve Storm
//
// Creation Date : 12/12/00
//------------------------------------------------------------------
CubitPointContainment CurveACIS::point_containment( const CubitVector &point )
{
   EDGE *EDGE_ptr = get_EDGE_ptr();
   if (EDGE_ptr->geometry() == NULL)
     return CUBIT_PNT_UNKNOWN;
  
   //const curve* acis_curve = &EDGE_ptr->geometry()->equation();
   
   SPAposition test_point( point.x(), point.y(), point.z() );
   
   SPAtransf ftrans;

   point_edge_containment pe_rel;

   ENTITY *ENTITY_ptr = NULL;
   SPAparameter param;

   pe_rel = sg_point_in_edge( test_point, EDGE_ptr, ftrans, ENTITY_ptr, param );

   switch( pe_rel )
   {
   case point_off_edge:
      return CUBIT_PNT_OFF;
   case point_on_edge:
      return CUBIT_PNT_ON;
   }
   return CUBIT_PNT_UNKNOWN;
}

//------------------------------------------------------------------
// Purpose       : Return a pointer to the ACIS EDGE
//                 in this Curve.
//
// Special Notes :
//
// Creator       : Malcolm J. Panthaki
//
// Creation Date : 2/25/97
//------------------------------------------------------------------
EDGE* CurveACIS::get_EDGE_ptr() const
{
  return (EDGE*)ENTITY_ptr();
}

void CurveACIS::set_EDGE_ptr(EDGE* EDGE_ptr)
{
  ENTITY_ptr(EDGE_ptr);
}


//------------------------------------------------------------------
// Purpose: This function returns the coordinate of a point in the local
//          parametric (u) space that corresponds to the input position 
//          in global (world) space.  The input point is first moved to 
//          the closest point on the Curve and the parameter value of 
//          that point is determined. 
//
// Creator       : Malcolm J. Panthaki
//
// Creation Date : 2/25/97
//-------------------------------------------------------------------
CubitStatus CurveACIS::position_from_u (double u_value,
                                        CubitVector& output_position)
{
    // Get the first ACIS curve
  curve const* curve_ptr = this->get_ACIS_curve();
  EDGE const* EDGE_ptr = this->get_EDGE_ptr();
  
    // If there is no geometry associated with this curve,
    // just return the coords of an endpoint
  if ( curve_ptr == NULL )
  {
    SPAposition position_on_curve;
    position_on_curve = EDGE_ptr->start()->geometry()->coords();
    output_position.x( position_on_curve.x() );
    output_position.y( position_on_curve.y() );
    output_position.z( position_on_curve.z() );
    return CUBIT_SUCCESS;
  }
  
  CubitStatus status = CUBIT_SUCCESS;
  
    // Get the SPAparameter range of this Curve
  SPAinterval param_interval = EDGE_ptr->param_range();
  
    // Make sure the requested u_value is either within the range or, if
    // the Curve is periodic, then reduce the input value down to the
    // fundamental range
  if ((param_interval >> u_value) == FALSE )
  {
    adjust_periodic_parameter(u_value);
    if ((param_interval >> u_value) == FALSE)
    {
//      assert(param_interval >> u_value == TRUE);
        // find the closest endpoint and return failure
      double old_u = u_value;
      if (param_interval.start_pt() < param_interval.end_pt())
      {
        if (u_value < param_interval.start_pt())
          u_value = param_interval.start_pt();
        else
          u_value = param_interval.end_pt();
      }
      else
      {
        if (u_value < param_interval.end_pt())
          u_value = param_interval.end_pt();
        else
          u_value = param_interval.start_pt();
      }
      if (fabs(u_value - old_u) > GEOMETRY_RESABS*10.0) {
        status = CUBIT_FAILURE;
        PRINT_ERROR("In CurveACIS::position_from_u\n"
                    "       Input parameter value %f is not within the"
                    " parameter range of this Curve (%f to %f)\n",
                    old_u, param_interval.start_pt(), param_interval.end_pt());
      }
    }
  }
  
    // Now that we have a "valid" parameter value, get its global location
    // on the Curve
  SPAposition position_on_curve;
  
    //Now we assume that the u_value is with respect to the
    //EDGE, we must switch it according the the EDGE->curve sense.
  if (EDGE_ptr->sense() == REVERSED )
  {
    u_value = -(u_value);
  }
    //now the u_value is with respect to curve so we are safe in
    //passing this parameter for evaluation on the curve.
  curve_ptr->eval( u_value, position_on_curve );
  output_position.x( position_on_curve.x() );
  output_position.y( position_on_curve.y() );
  output_position.z( position_on_curve.z() );
  return status;
}

//-------------------------------------------------------------------------
// Purpose       : This function returns the coordinate of a point in the local
//                 parametric (u) space that corresponds to the input position 
//                 in global (world) space.  The input point is first moved to 
//                 the closest point on the Curve and the parameter value of 
//                 that point is determined. 
//
// Special Notes : 
//
// Creator       : Malcolm J. Panthaki
//
// Creation Date : 2/25/97
//-------------------------------------------------------------------------
double CurveACIS::u_from_position (const CubitVector& input_position)
{
    // Get the closest point on the Curve to the input position
  CubitVector closest_point;
  double u_val;
  this->closest_point(input_position, closest_point,
                      NULL, NULL, &u_val);
    // closest_point already makes adjustments for sense and periodicity
  
  return u_val;
}

//------------------------------------------------------------------
// Purpose: This function returns the parameter value of the point 
//          that is "arc_length" away from the root point, in the
//          positive sense direction of the owning RefEdge.
//
// Special Notes : 
//   If arc_length is negative, the new point (whose parameter value
//   is being computed) is in the negative sense direction (along
//   the RefEdge) from the root point (whose parameter value is
//   root_param).
//
//   If the curve is not periodic and the new point, "arc_length"
//   away from the root point in the appropriate direction, goes
//   beyond the end point of the first EDGE, that end point is used
//   to generate the returned parameter value.
//
// If the curve is periodic and the new point, "arc_length" away
// from the root point in the appropriate direction, goes beyond
// the end point of the first EDGE, wrap around is done.  After
// wrap around, the point is treated as with other curves
//
// NOTE:
// The important assumption that is made in this routine is that
// the end points of the RefEdge that owns this CurveACIS are the
// same as the end points of the first ACIS EDGE in the list of EDGEs
// associated with this CurveACIS.
//
// Assume that the parameter root_param is with respect to the
// RefEdge as well as arc_length.  Before calling the ACIS "curve",
// we need to get them with respect to the curve.   
//
// Creator       : Malcolm J. Panthaki
//
// Creation Date : 2/28/97
//------------------------------------------------------------------
double CurveACIS::u_from_arc_length ( double root_param,
                                      double arc_length )
{
    // Sanity check
  if (arc_length == 0.0)
  {
    return root_param;
  }
    // get the EDGE and curve
  EDGE const* EDGE_ptr = this->get_EDGE_ptr();
  const curve* curve_ptr = this->get_ACIS_curve();
    // Check for NULL geometry
  if ( curve_ptr == NULL )
    return start_param();
  
    // Get the lower and upper parameter values for the first EDGE.
    // Note that the lower and upper bounds are with respect to the start and
    // end locations of the EDGE, not the RefEdge.
  double EDGE_start_param = (double)(EDGE_ptr->start_param());
  double EDGE_end_param = (double)(EDGE_ptr->end_param());
  
    // If the Curve is periodic, then let ACIS do the work.  Keep in mind 
    // that ACIS EDGEs also have a sense with respect to their underlying 
    // ACIS curves.
  double period = 0.0;
  if ( this->is_periodic(period) == CUBIT_TRUE ) 
  {
      // If the sense of the EDGE with respect to its curve is REVERSED,
      // switch the sign of arc_length so we are searching in the right
      // direction
    if ( EDGE_ptr->sense() == REVERSED )
    {
      arc_length = -(arc_length);
      root_param = -(root_param);
    }
    
      // Do the needed work on the curve and return the parameter value.
      // NOTE: The following function does not work correctly if the
      //       curve is bounded and the required arc_length puts you over
      //       the end of the curve.  The result, in that case, is *NOT 
      //       DEFINED* and there is no error reported!! We can use it
      //       as the curve is periodic.
    double param = curve_ptr->length_param(root_param, arc_length);
    
      // Now this is the parameter with respect to the curve.  Reverse sign if
      // necessary to make it wrt the EDGE.
    if ( EDGE_ptr->sense() == REVERSED )
      param = -(param);
    adjust_periodic_parameter(param);
    return param;
  }
  
    // Now, let's deal with non-periodic Curves...
  
    // Find out how much "headroom" there is between the root_point and the 
    // "relevant" end of the EDGE.
    //
    // If arc_length is negative, then we have to find the 
    // distance between root point and the start of the EDGE. If arc_length
    // is positive, then we have to find the distance between root_point and
    // the end of the EDGE.
  double headroom = CUBIT_DBL_MAX;
  if(arc_length < 0.0)
    headroom = this->length_from_u( root_param, EDGE_start_param );
  else
    headroom = this->length_from_u( root_param, EDGE_end_param );
  
    // If the arc_length specifies a point past the "relevant" end of the 
    // EDGE, choose that end -- i.e., don't go beyond the end of the
    // EDGE.
  if (fabs(arc_length) >  headroom) 
  {
      // This EDGE is not "periodic", so stop at the relevant endpoint
      // and return its parameter value.
      // NOTE: We have already taken into account the sense of the EDGE
      //       with respect to its underlying curve...whew!!! :-)
    if (arc_length > 0.0)
    {
      return EDGE_end_param;
    }
    else
    {
      return EDGE_start_param;
    }
  }
  
    // The EDGE is non-periodic but arc_length doesn't put us past the end of 
    // the EDGE.
    // NOTE: The following function does not work correctly if the
    //       curve is bounded and the required arc_length puts you over
    //       the end of the curve.  The result, in that case, is *NOT 
    //       DEFINED* and there is no error reported!!  We can use it
    //       as we've already checked to make sure we don't go off the
    //       edge of the EDGE :-)
  else
  {
      // Assert that root_param lies on the bounded EDGE
    assert( EDGE_end_param > EDGE_start_param ?
            root_param >= (EDGE_start_param - GEOMETRY_RESABS*10) &&
            root_param <= (EDGE_end_param + GEOMETRY_RESABS*10) :
            root_param <= (EDGE_start_param + GEOMETRY_RESABS*10) &&
            root_param >= (EDGE_end_param - GEOMETRY_RESABS*10) );
    
      // Change sign for calculation, if needed
    if (EDGE_ptr->sense() == REVERSED)
    {
      arc_length = -(arc_length);
      root_param = -(root_param);
    }
    double param = curve_ptr->length_param( root_param, arc_length );
      // change sign of result, if needed
    if ( EDGE_ptr->sense() == REVERSED )
      param = -(param);
    else if (param < EDGE_start_param || param > EDGE_end_param) {
        // problem using parameterization - just estimate instead
      param = EDGE_start_param + (EDGE_end_param-EDGE_start_param)*
        (arc_length/this->length_from_u(EDGE_start_param, EDGE_end_param));
    }
    
    return param;
  }
}

//-------------------------------------------------------------------------
// Purpose       : This function tests the passed in position to see if
//                 is on the underlying curve. 
//
// Special Notes :
//
// Creator       : David R. White
//
// Creation Date : 04/08/97
//-------------------------------------------------------------------------
CubitBoolean CurveACIS::is_position_on( const CubitVector &test_position )
{
  if ( get_EDGE_ptr()->geometry() == NULL )
  {
    SPAposition some_point = get_EDGE_ptr()->start()->geometry()->coords();
    CubitVector test_two;
    test_two.set( some_point.x(), some_point.y(), some_point.z());
    return GeometryQueryTool::instance()->
      about_spatially_equal( test_two,
                             test_position,
                             GeometryQueryTool::get_geometry_factor() );
  }
  
  SPAposition pos ( test_position.x(), test_position.y(), test_position.z() );
  return
    get_EDGE_ptr()->geometry()->equation().test_point(pos) ? CUBIT_TRUE : CUBIT_FALSE;
}

//-------------------------------------------------------------------------
// Purpose       : This function returns the type of underlying curve. 
//
// Special Notes : It checks to see if *any* of the ACIS curves associated
//                 with the EDGEs in the list of EDGEs of this Curve is of
//                 a particular type and returns the appropriate value
//                 of the enum, CurveType.
//
// Creator       : Malcolm J. Panthaki
//
// Creation Date : 3/3/97
//-------------------------------------------------------------------------
GeometryType CurveACIS::geometry_type()
{
  EDGE* EDGE_ptr = get_EDGE_ptr();
  
    // Get its type (note that "straight_type", etc. are defined by
    // ACIS)
  if ( EDGE_ptr->geometry() == NULL )
    return POINT_CURVE_TYPE;
  
  int type = (&(EDGE_ptr->geometry()->equation()))->type();
  
  GeometryType local_type;
  
  switch (type)
  {
    case straight_type:
      local_type = STRAIGHT_CURVE_TYPE;
      break;
    case ellipse_type:
    {
      ELLIPSE *ellipse = static_cast<ELLIPSE*>( EDGE_ptr->geometry() );
      if( fabs( ellipse->radius_ratio() - 1) < GEOMETRY_RESABS )
        local_type = ARC_CURVE_TYPE;
      else
        local_type = ELLIPSE_CURVE_TYPE;
      break;
    }
    case intcurve_type:
      local_type = SPLINE_CURVE_TYPE;
      break;
    default:
      local_type = UNDEFINED_CURVE_TYPE;
      break;
  }
  
  return local_type;
}

CubitStatus CurveACIS::get_point_direction( CubitVector& point, 
                                            CubitVector& direction )
{
  if( geometry_type() != STRAIGHT_CURVE_TYPE )
    return CUBIT_FAILURE;
  
  EDGE* EDGE_ptr = get_EDGE_ptr();
  
  curve const* acis_curve = &(EDGE_ptr->geometry()->equation());
  
  straight const* straight_curve = (straight *)acis_curve;
  
  SPAposition pos_curve_orig;
  SPAunit_vector uv_curve_dir;
  pos_curve_orig = straight_curve->root_point;
  uv_curve_dir = straight_curve->direction;
  
  point.set( pos_curve_orig.x(), pos_curve_orig.y(), pos_curve_orig.z() );
  direction.set( uv_curve_dir.x(), uv_curve_dir.y(), uv_curve_dir.z() );
  
  return CUBIT_SUCCESS;
}

CubitStatus CurveACIS::get_center_radius( CubitVector& center, 
                                         double& radius )
{
  if( geometry_type() != ELLIPSE_CURVE_TYPE && 
      geometry_type() != ARC_CURVE_TYPE )
    return CUBIT_FAILURE;
  
  EDGE* EDGE_ptr = get_EDGE_ptr();
  
  curve const* acis_curve = &(EDGE_ptr->geometry()->equation());
  
  ellipse const* ellipse_curve = (ellipse *)acis_curve;
  
  SPAposition pos_curve_orig;
  pos_curve_orig = ellipse_curve->centre;

  SPAvector major_axis;
  major_axis = ellipse_curve->major_axis;

  radius = major_axis.len();
  
  center.set( pos_curve_orig.x(), pos_curve_orig.y(), pos_curve_orig.z() );
  
  return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Purpose       : This function returns the start parameter.
//
// Special Notes : The start param is with respect to the ref_edge.
//
// Creator       : David White
//
// Creation Date : 03/24/97
//-------------------------------------------------------------------------
double CurveACIS::start_param()
{
   double start = 0.0, end = 0.0;
   
   get_param_range( start, end );
   return start;
}
//-------------------------------------------------------------------------
// Purpose       : This function returns the end parameter.
//
// Special Notes : The end param is with respect to the ref_edge.
//
// Creator       : David White
//
// Creation Date : 03/24/97
//-------------------------------------------------------------------------
double CurveACIS::end_param()
{
   double start = 0.0, end = 0.0;
   
   get_param_range( start, end );
   return end;
}

int CurveACIS::validate(const CubitString &user_name,
                        DLIList <TopologyEntity*> &bad_entities)
{
  int error = check_EDGE(get_EDGE_ptr(), user_name);
  if ( error > 0 )
    bad_entities.append(this->topology_entity());
  return error;
}

//=========  Add Code by SRS of Cat,  3/3/99 2:38:26 PM  =========
CubitBoolean CurveACIS::is_tolerant()
{
   EDGE* EDGE_ptr = this->get_EDGE_ptr(); 
#if CUBIT_ACIS_VERSION < 500
   return CUBIT_FALSE;
#else
   if( is_TEDGE( EDGE_ptr ) )
      return CUBIT_TRUE;
   else
      return CUBIT_FALSE;
#endif
}
//=========  Code End by SRS of Cat,  3/3/99 2:38:26 PM  =========

int CurveACIS::check_EDGE(EDGE *edge,
                          const CubitString &refentity_name)
{
  int error = 0;
  assert(edge != NULL);
  
  check_status_list *list = NULL;
  api_check_edge(edge, list );
  while(list != NULL){
    error++;
    switch(list->status()){
      case check_irregular:
        PRINT_WARNING("twisted or scrunched up curve for %s\n",
                      refentity_name.c_str());
        break;
      case check_self_intersects:
        PRINT_WARNING("self-intersecting curve for %s\n",
                      refentity_name.c_str());
        break;
      case check_bad_closure:
        PRINT_WARNING("closure is wrong for %s\n",
                      refentity_name.c_str());
        break;
      case check_bs3_null:
        PRINT_WARNING("no bs3 curve for %s\n",
                      refentity_name.c_str());
        break;
      case check_bs3_coi_verts:
        PRINT_WARNING("error in control point coincidence for %s\n",
                      refentity_name.c_str());
        break;
      case check_bad_degeneracies:
          // This case is never returned for curves
        break;
      case check_untreatable_singularity:
          // This case is never returned for curves
        break;
      case check_non_G0:
        PRINT_WARNING("curve for %s is not G0\n",
                      refentity_name.c_str());
        break;
      case check_non_G1:
        PRINT_WARNING("curve for %s is not G1\n",
                      refentity_name.c_str());
        break;
      case check_non_G2:
        PRINT_WARNING("curve for %s is not G2\n",
                      refentity_name.c_str());
        break;
      case check_non_C1:
        PRINT_WARNING("curve for %s is not C1\n",
                      refentity_name.c_str());
        break;
      case check_unknown:
        break;
      default:
        break;
    }
    list = list->next();
  }
  return error;
}
/*
void CurveACIS::bodysms(DLIList<BodySM*> &bodies) 
{
  AcisBridge::bodysms(bodies);
}

void CurveACIS::lumps(DLIList<Lump*> &lumps)
{
  AcisBridge::lumps(lumps);
}

void CurveACIS::shellsms(DLIList<ShellSM*> &shellsms)
{
  AcisBridge::shellsms(shellsms);
}

void CurveACIS::surfaces(DLIList<Surface*> &surfaces)
{
  AcisBridge::surfaces(surfaces);
}

void CurveACIS::loopsms(DLIList<LoopSM*> &loopsms)
{
  AcisBridge::loopsms(loopsms);
}

void CurveACIS::curves(DLIList<Curve*> &curves)
{
  AcisBridge::curves(curves);
}

void CurveACIS::coedgesms(DLIList<CoEdgeSM*> &coedgesms)
{
  AcisBridge::coedgesms(coedgesms);
}

void CurveACIS::points(DLIList<Point*> &points)
{
  AcisBridge::points(points);
}
*/

void CurveACIS::get_parents_virt( DLIList<TopologyBridge*>& parents )
{
  ENTITY_LIST entities;
  api_get_coedges( get_EDGE_ptr(), entities );
  ATTRIB_CUBIT_OWNER::cubit_owner( entities, parents );
}

void CurveACIS::get_children_virt( DLIList<TopologyBridge*>& children )
{
  ENTITY_LIST entities;
  api_get_vertices( get_EDGE_ptr(), entities );
  ATTRIB_CUBIT_OWNER::cubit_owner( entities, children );
}


//-------------------------------------------------------------------------
// Purpose       : Check for G1 discontinuity
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 05/16/00
//-------------------------------------------------------------------------
CubitBoolean CurveACIS::G1_discontinuous( 
      double /*param*/, CubitVector* /*mtan*/, CubitVector* /*ptan*/ )
{ 
/*
  const EDGE* EDGE_ptr = get_EDGE_ptr();
  const CURVE* CURVE_ptr = EDGE_ptr->geometry();
  if( CURVE_ptr == NULL ) //point curve
    return CUBIT_FALSE;
  
  const curve* curve_ptr = &(CURVE_ptr->equation());
  if( !is_intcurve( curve_ptr ) )
    return CUBIT_FALSE;
    
  const intcurve* intcurve_ptr = (const intcurve*)curve_ptr;
  bs3_curve spline_curve_ptr = intcurve_ptr->cur();
  if( spline_curve_ptr == NULL )
    return CUBIT_FALSE;
  
  if (EDGE_ptr->sense() == REVERSED )
  {
    param = -(param);
  }

  logical G1, C1;
  bs3_curve_smoothness( spline_curve_ptr, param, C1, G1 );
  return G1 ? CUBIT_FALSE : Curve::G1_discontinuous( param, mtan, ptan );
*/
  return CUBIT_FALSE;
}

// ********** END PUBLIC FUNCTIONS         **********

// ********** BEGIN PROTECTED FUNCTIONS    **********
// ********** END PROTECTED FUNCTIONS      **********

// ********** BEGIN PRIVATE FUNCTIONS      **********

//-------------------------------------------------------------------------
// Purpose       : Return a pointer to the ACIS curve associated with
//                 the first EDGE in the list of EDGEs in this Curve.
//
// Special Notes :
//
// Creator       : Malcolm J. Panthaki
//
// Creation Date : 2/25/97
//-------------------------------------------------------------------------
curve const* CurveACIS::get_ACIS_curve()
{
     // Get the first EDGE
   EDGE* EDGE_ptr = this->get_EDGE_ptr();
   
     // Get its curve
   if ( EDGE_ptr->geometry() == NULL )
       return (curve const*) NULL;
   return &(EDGE_ptr->geometry()->equation());
}

//-------------------------------------------------------------------------
// Purpose       : Return the sense of the first EDGE (in the list of EDGEs
//                 in this Curve) wrt its underlying ACIS curve.
//
// Special Notes : 
//
// Creator       : Malcolm J. Panthaki
//
// Creation Date : 2/25/97
//-------------------------------------------------------------------------
CubitSense CurveACIS::get_EDGE_sense()
{
     // Get the first EDGE
   EDGE* EDGE_ptr = this->get_EDGE_ptr();
   
     // Return the sense value
   REVBIT EDGEsense = EDGE_ptr->sense();
   if (EDGEsense == FORWARD)
   {
      return CUBIT_FORWARD;
   }
   
   else if (EDGEsense == REVERSED)
   {
      return CUBIT_REVERSED;
   }
   
   else
   {
      PRINT_ERROR("In CurveACIS::get_EDGE_sense\n"
                  "       Could not get the sense of the EDGE wrt its curve.\n");
      return CUBIT_FORWARD;
   }
}

//----------------------------------------------------------------
// Adjusts the input parameter so that it falls within the
// parameter range of this Curve, if possible.  Necessary for
// periodic curves.
//----------------------------------------------------------------
void CurveACIS::adjust_periodic_parameter(double& param)
{
    // Adjustment only legal if this is a periodic curve.
  double period;
  if ( this->is_periodic(period) && (fabs(period) > CUBIT_RESABS))
  {
    EDGE const* EDGE_ptr = this->get_EDGE_ptr();

    SPAinterval edge_range = EDGE_ptr->param_range();
		assert( edge_range.length() > CUBIT_RESABS * 100 );

    SPAinterval range( edge_range.start_pt() - SPAresabs,
		                edge_range.  end_pt() + SPAresabs );
    
      // Make sure period is positive
    if (period < 0)
      period = -period;

      // Move the parameter above the low param
    while (param < range)
      param += period;
      // Move the parameter below the high param
    while (param > range)
      param -= period;
  }
}

// ********** END PRIVATE FUNCTIONS        **********

// ********** BEGIN HELPER CLASSES         **********
// ********** END HELPER CLASSES           **********

// ********** BEGIN EXTERN FUNCTIONS       **********
// ********** END EXTERN FUNCTIONS         **********

// ********** BEGIN STATIC FUNCTIONS       **********
// ********** END STATIC FUNCTIONS         **********
