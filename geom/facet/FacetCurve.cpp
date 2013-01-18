//-------------------------------------------------------------------------
// Filename      : FacetCurve.cpp
//
// Purpose       : 
//
// Special Notes :
//
// Creator       : Steven J. Owen
//
// Creation Date : 07/14/00
//
// Owner         : Steven J. Owen
//-------------------------------------------------------------------------

// ********** BEGIN STANDARD INCLUDES      **********

// ********** END STANDARD INCLUDES        **********

// ********** BEGIN CUBIT INCLUDES         **********

#include "CastTo.hpp"
#include "CubitVector.hpp"
#include "CubitBox.hpp"
#include "GeometryDefines.h"
#include "FacetCurve.hpp"
#include "FacetAttrib.hpp"
#include "FacetEvalTool.hpp"
#include "CurveFacetEvalTool.hpp"
#include "GeometryQueryEngine.hpp"
#include "FacetQueryEngine.hpp"
#include "CoEdgeSM.hpp"
#include "CubitPoint.hpp"

#include "FacetBody.hpp"
#include "FacetLump.hpp"
#include "FacetShell.hpp"
#include "FacetSurface.hpp"
#include "FacetLoop.hpp"
#include "FacetCoEdge.hpp"
#include "FacetPoint.hpp"

// ********** END CUBIT INCLUDES           **********

// ********** BEGIN FORWARD DECLARATIONS   **********
// ********** END FORWARD DECLARATIONS     **********

// ********** BEGIN STATIC DECLARATIONS    **********
// ********** END STATIC DECLARATIONS      **********

// ********** BEGIN PUBLIC FUNCTIONS       **********

//-------------------------------------------------------------------------
// Purpose       : The default constructor
//
// Special Notes :
//
// Creator:      : Steve Owen
//
// Creation Date : 07/14/00
//-------------------------------------------------------------------------
FacetCurve::FacetCurve(CurveFacetEvalTool *curve_facet_tool,
                       TBPoint *start_ptr, TBPoint *end_ptr,
                       DLIList<CoEdgeSM*> &coedge_list )
                       : sense_(CUBIT_FORWARD)
{
  static int counter = 0;
  myId = counter++;
    // Calculate a bounding box if there isn't one already
  curveFacetEvalTool = curve_facet_tool; 
  myStartPoint = start_ptr;
  myEndPoint = end_ptr;
  myCoEdges += coedge_list;
  periodic = start_ptr == end_ptr;
}

//-------------------------------------------------------------------------
// Purpose       : Another constructor
//
// Special Notes : Implemented for save/restore
//
// Creator:      : Corey Ernst 
//
// Creation Date : 02/03/03
//-------------------------------------------------------------------------
FacetCurve::FacetCurve(CurveFacetEvalTool *curve_facet_tool,
                       TBPoint *start_ptr, TBPoint *end_ptr,
                       CubitSense sense )
{
  static int counter = 0;
  myId = counter++;
    // Calculate a bounding box if there isn't one already
  curveFacetEvalTool = curve_facet_tool; 
  myStartPoint = start_ptr;
  myEndPoint = end_ptr;
  periodic = start_ptr == end_ptr;
  sense_ = sense;
}
//-------------------------------------------------------------------------
// Purpose       : The destructor. 
//
// Special Notes :
//
// Creator       : Steve Owen
//
// Creation Date : 07/14/00
//-------------------------------------------------------------------------
FacetCurve::~FacetCurve() 
{
    if(this->curveFacetEvalTool)
        delete this->curveFacetEvalTool;
  
}

//-------------------------------------------------------------------------
// Purpose       : The purpose of this function is to append a
//                 attribute to the GE. The name is attached to the 
//                 underlying solid model entity this one points to.
//
//
// Special Notes : 
//
// Creator       : Steve Owen
//
// Creation Date : 07/14/00
//-------------------------------------------------------------------------
void FacetCurve::append_simple_attribute_virt(const CubitSimpleAttrib &csa)
  { attribSet.append_attribute(csa); }

//-------------------------------------------------------------------------
// Purpose       : The purpose of this function is to remove a simple 
//                 attribute attached to this geometry entity. The name is 
//                 removed from the underlying BODY this points to.
//
// Special Notes : 
//
// Creator       : Steve Owen
//
// Creation Date : 07/14/00
//-------------------------------------------------------------------------
void FacetCurve::remove_simple_attribute_virt(const CubitSimpleAttrib &csa)
  { attribSet.remove_attribute(csa); }

//-------------------------------------------------------------------------
// Purpose       : The purpose of this function is to remove all simple 
//                 attributes attached to this geometry entity.  Also
//                 removes lingering GTC attributes.
//
//
// Special Notes : 
//
// Creator       : Steve Owen
//
// Creation Date : 07/14/00
//-------------------------------------------------------------------------
void FacetCurve::remove_all_simple_attribute_virt()
  { attribSet.remove_all_attributes(); }

//-------------------------------------------------------------------------
// Purpose       : The purpose of this function is to get the  
//                 attributes attached to this geometry entity. The name is 
//                 attached to the underlying BODY this points to.
//
// Special Notes : 
//
// Creator       : Steve Owen
//
// Creation Date : 07/14/00
//-------------------------------------------------------------------------
CubitStatus FacetCurve::get_simple_attribute(DLIList<CubitSimpleAttrib>&
                                               csa_list)
  { return attribSet.get_attributes(csa_list); }
  
CubitStatus FacetCurve::get_simple_attribute( const CubitString& name,
                                      DLIList<CubitSimpleAttrib>& csa_list)
  { return attribSet.get_attributes( name, csa_list ); }


CubitStatus FacetCurve::save_attribs( FILE *file_ptr )
  { return attribSet.save_attributes(file_ptr); }

CubitStatus FacetCurve::restore_attribs( FILE *file_ptr, unsigned int endian )
  { return attribSet.restore_attributes(file_ptr, endian); }



//-------------------------------------------------------------------------
// Purpose       : Get geometry modeling engine: FacetQueryEngine
//
// Special Notes :
//
// Creator       : Steve Owen
//
// Creation Date : 07/14/00
//-------------------------------------------------------------------------
GeometryQueryEngine* FacetCurve::get_geometry_query_engine() const
{
  return FacetQueryEngine::instance();
}                 

//-------------------------------------------------------------------------
// Purpose       : Get the bounding box of the object.
//
// Special Notes :
//
// Creator       : Steve Owen
//
// Creation Date : 10/23/96
//-------------------------------------------------------------------------
CubitBox FacetCurve::bounding_box() const 
{
  return curveFacetEvalTool->bounding_box();
}


//-------------------------------------------------------------------------
// Purpose       : Return the length of the curve.
//
// Special Notes :
//
// Creator       : Steve Owen
//
// Creation Date : 07/14/00
//-------------------------------------------------------------------------
double FacetCurve::measure()
{
  return curveFacetEvalTool->length();
}

//-------------------------------------------------------------------------
// Purpose       : Return the arc length along the Curve starting from
//                 the point represented by the parameter1 going to the 
//                 point represented by parameter2.
//
// Special Notes : Parameter1 and parameter2 are with respect to the EDGE.
//
// Creator       : Steve Owen
//
// Creation Date : 07/14/00
//-------------------------------------------------------------------------
double FacetCurve::length_from_u( double parameter1, double parameter2 )
{

    //Don't perform the periodic adjustment so that
    // this function will satisfy the properties of
    // the function as defined in RefEdge.hpp.
    // Also, take the fabs of the arc length for the same reason.
  //if( periodic )
  //{
  //  adjust_periodic_parameter( parameter1 );
  //  adjust_periodic_parameter( parameter2 );
  //}
  return curveFacetEvalTool->length_from_u( parameter1, parameter2 );
}

//-------------------------------------------------------------------------
// Purpose       : Returns CUBIT_TRUE and the associated period value. Not
//                 implemented yet
//
// Special Notes :  
//
// Creator       : Steve Owen
//
// Creation Date : 07/14/00
//-------------------------------------------------------------------------
CubitBoolean FacetCurve::is_periodic(double& period)
{
  if( periodic )
  {
    period = 1.0;
    return CUBIT_TRUE;
  }
  return CUBIT_FALSE;
}

//------------------------------------------------------------------
// Purpose: Returns CUBIT_TRUE and the associated parametric values, 
//          if the facet curve associated with the first EDGE is 
//          parametric.
//          Otherwise returns CUBIT_FALSE and the values of 
//          the lower and upper parametric bounds are undetermined.
//          NOT IMPLEMENTED YET
//
// Creator       : Steve Owen
//
// Creation Date : 07/14/00
//-------------------------------------------------------------------
CubitBoolean FacetCurve::get_param_range( double& lower_bound,
                                          double& upper_bound )
{
  lower_bound = 0.0;
  upper_bound = 1.0;
  return CUBIT_TRUE;
}

//------------------------------------------------------------------
// Purpose:        Finds the extrema along this Curve. 
//
// Special Notes : It is the responsibility of the
//                 calling code to delete the CubitVectors added to 
//                 interior_points!
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 05/29/01
//-------------------------------------------------------------------
CubitStatus FacetCurve::get_interior_extrema(
  DLIList<CubitVector*>& interior_points,
  CubitSense& return_sense )
{
  // Need curveFacetEvalTool to get point list
  if( ! curveFacetEvalTool ) 
    return CUBIT_FAILURE;
  
  // Get list of points defining this curve
  DLIList<CubitPoint*> point_list;
  get_points( point_list );
  
  // If there are only 2 points, then the curve is a line and there
  // are no interior extrema
  if( point_list.size() < 3 )
    return CUBIT_SUCCESS;
  
  // Return sense is whatever the sense of this curve is.
  return_sense = sense_;
  
  // Get a vector between the first two points
  point_list.reset();
  CubitVector prev_pt = point_list.get_and_step()->coordinates();
  CubitVector curr_pt = point_list.get_and_step()->coordinates();
  CubitVector prev_vct = curr_pt - prev_pt;
  CubitVector next_vct;
  
  for( int i = point_list.size(); i > 2; i-- )
  {
    // Get a vector between the next two points
    next_vct = point_list.get()->coordinates() - curr_pt;
    
    // In CurveACIS::get_interior_extrema, the extrema seem to
    // be evaluated with respect to the principle axes, so do
    // the same here.  The extrema are points at which the
    // derivitive in the specified direction (principle axis)
    // is zero.  So look for a sign change in the slope across
    // a point wrt each principle direction.
    if( (prev_vct.x() * next_vct.x() < 0.) ||  // x extrema
        (prev_vct.y() * next_vct.y() < 0.) ||  // y extrema
        (prev_vct.z() * next_vct.z() < 0.)  )  // z extrema
      interior_points.append( new CubitVector( curr_pt ) );
    
    // Advance to next point.
    prev_vct = next_vct;
    curr_pt = point_list.get_and_step()->coordinates();
  }
  
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
// Creator       : Steve Owen
//
// Creation Date : 07/14/00
//-------------------------------------------------------------------------
CubitStatus FacetCurve::closest_point( 
  CubitVector const& location, 
  CubitVector& closest_location,
  CubitVector* tangent_ptr,
  CubitVector* curvature_ptr,
  double* param)
{  
    // Only the closest point is required
  if (tangent_ptr == NULL && curvature_ptr == NULL)
  {
    CubitVector temp = location;
    curveFacetEvalTool->closest_point( temp,
                                       closest_location, NULL, NULL, param );
  }
  
    // The closest point and the tangent are required
  else if (tangent_ptr != NULL && curvature_ptr == NULL)
  {
    CubitVector temp = location;
    curveFacetEvalTool->closest_point( temp, closest_location, tangent_ptr, NULL, param );
  }
  
    // Everything is required
    // NOTE: If the curvature is required but not the tangent, the tangent,
    //       too, will be computed.
  else
  {
    CubitVector temp = location;
    curveFacetEvalTool->closest_point( temp, closest_location, tangent_ptr, curvature_ptr, param );
  }
  
  if (tangent_ptr != NULL)
  {
    if (sense_ == CUBIT_REVERSED)
      *tangent_ptr = -(*tangent_ptr);
  }
  
    // Set the curvature, if necessary
  if (curvature_ptr != NULL)
  {
      // get the sense wrt the acis curve, consists of 2 sense values
    if (sense_ == CUBIT_REVERSED)
      *curvature_ptr = -(*curvature_ptr);
  }
  
  if (param != NULL)
  {
    adjust_periodic_parameter(*param);
  }
  
  return CUBIT_SUCCESS;
}


//------------------------------------------------------------------
// Purpose: This function returns the coordinate of a point in the local
//          parametric (u) space that corresponds to the input position 
//          in global (world) space.  The input point is first moved to 
//          the closest point on the Curve and the parameter value of 
//          that point is determined. 
//
// Creator       : Steve Owen
//
// Creation Date : 07/14/00
//-------------------------------------------------------------------
CubitStatus FacetCurve::position_from_u (double u_value,
                                        CubitVector& output_position)
{
  
  CubitStatus status = CUBIT_SUCCESS;
  
    // Get the parameter range of this Curve
  double lower_bound, upper_bound;
  this->get_param_range( lower_bound, upper_bound );
  double param_interval = upper_bound - lower_bound;
  
    // Make sure the requested u_value is either within the range or, if
    // the Curve is periodic, then reduce the input value down to the
    // fundamental range
  if (u_value > upper_bound || u_value < lower_bound)
  {
    adjust_periodic_parameter(u_value);
  }
  
    // Now that we have a "valid" parameter value, get its global location
    // on the Curve
  
    // Now we assume that the u_value is with respect to the
    // EDGE, we must switch it according the the EDGE->curve sense.
  
  if (sense_ == CUBIT_REVERSED)
  {
    u_value = -(u_value);
  }
    //now the u_value is with respect to curve so we are safe in
    //passing this parameter for evaluation on the curve.

  double fraction = (u_value - lower_bound) / param_interval;
  status = curveFacetEvalTool->position_from_fraction( fraction, output_position );

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
double FacetCurve::u_from_position (const CubitVector& input_position)
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
double FacetCurve::u_from_arc_length ( double root_param,
                                       double arc_length )
{

  if (!curveFacetEvalTool)
  {
    PRINT_ERROR("curve facet evaluation tool not defined in FacetCurve::u_from_arc_length\n");
    return 0.0;
  }
  return curveFacetEvalTool->u_from_arc_length( root_param, arc_length );

}

//-------------------------------------------------------------------------
// Purpose       : This function tests the passed in position to see if
//                 is on the underlying curve. 
//
// Special Notes :
//
// Creator       : Steve Owen
//
// Creation Date : 07/14/00
//-------------------------------------------------------------------------
CubitBoolean FacetCurve::is_position_on( const CubitVector &test_position )
{
  CubitVector new_point;
  CubitStatus stat = closest_point(test_position, new_point, NULL,NULL,NULL);

  if ( !stat )
    return CUBIT_FALSE;
  CubitVector result_vec = test_position - new_point;
  if ( result_vec.length_squared() < GEOMETRY_RESABS )
    return CUBIT_TRUE;
  return CUBIT_FALSE;
}

//-------------------------------------------------------------------------
// Purpose       : This function returns the type of underlying curve. 
//
// Special Notes : It checks to see if *any* of the ACIS curves associated
//                 with the EDGEs in the list of EDGEs of this Curve is of
//                 a particular type and returns the appropriate value
//                 of the enum, CurveType.
//
// Creator       : Steve Owen
//
// Creation Date : 07/14/00
//-------------------------------------------------------------------------
GeometryType FacetCurve::geometry_type()
{
  return SEGMENTED_CURVE_TYPE;
}

//-------------------------------------------------------------------------
// Purpose       : Return direction of point on curve
//
// Special Notes : not currently implemented
//
// Creator       : Steve Owen
//
// Creation Date : 07/14/00
//-------------------------------------------------------------------------
CubitStatus FacetCurve::get_point_direction( CubitVector& point, 
                                             CubitVector& direction )
{
  point = point;
  direction = direction;
  PRINT_DEBUG_122("FacetCurve::get_point_direction currently not implemented.\n");
  return CUBIT_FAILURE;
}

//-------------------------------------------------------------------------
// Purpose       : Return the center and radius of an arc
//
// Special Notes : not currently implemented
//
// Creator       : Steve Owen
//
// Creation Date : 07/14/00
//-------------------------------------------------------------------------
CubitStatus FacetCurve::get_center_radius( CubitVector& center, 
                                           double& radius )
{
  center = center;
  radius = radius;
  PRINT_DEBUG_122("FacetCurve::get_center_radius currently not implemented.\n");
  return CUBIT_FAILURE;
}

//-------------------------------------------------------------------------
// Purpose       : This function returns the start parameter.
//
// Special Notes : The start param is with respect to the ref_edge.
//
// Creator       : Steve Owen
//
// Creation Date : 07/14/00
//-------------------------------------------------------------------------
double FacetCurve::start_param()
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
// Creator       : Steve Owen
//
// Creation Date : 07/14/00
//-------------------------------------------------------------------------
double FacetCurve::end_param()
{
   double start = 0.0, end = 0.0;
   
   get_param_range( start, end );
   return end;
}

/*
void FacetCurve::bodysms(DLIList<BodySM*> &bodies)
{
  int ii;
  for ( ii = myCoEdges.size(); ii > 0; ii-- )
  {
    myCoEdges.get_and_step()->bodysms(bodies);
  }
}

void FacetCurve::lumps(DLIList<Lump*> &lumps)
{
  int ii;
  for ( ii = myCoEdges.size(); ii > 0; ii-- )
  {
    myCoEdges.get_and_step()->lumps(lumps);
  }
}

void FacetCurve::shellsms(DLIList<ShellSM*> &shellsms)
{
  int ii;
  for ( ii = myCoEdges.size(); ii > 0; ii-- )
  {
    myCoEdges.get_and_step()->shellsms(shellsms);
  }
}

void FacetCurve::surfaces(DLIList<Surface*> &surfaces)
{
  int ii;
  for ( ii = myCoEdges.size(); ii > 0; ii-- )
  {
    myCoEdges.get_and_step()->surfaces(surfaces);
  }
}

void FacetCurve::loopsms(DLIList<LoopSM*> &loopsms)
{
  int ii; 
  for ( ii = myCoEdges.size(); ii > 0; ii-- )
  {
    myCoEdges.get_and_step()->loopsms(loopsms);
  }
}


void FacetCurve::coedgesms(DLIList<CoEdgeSM*> &coedgesms)
{
  int ii; 
  for ( ii = myCoEdges.size(); ii > 0; ii-- )
  {
    coedgesms.append_unique( myCoEdges.get_and_step() );
  } 
}

void FacetCurve::curves(DLIList<Curve*> &curves)
{
  curves.append_unique( this );
}

void FacetCurve::points(DLIList<TBPoint*> &points)
{
  points.append_unique( myStartPoint );
  points.append_unique( myEndPoint );
}
*/


void FacetCurve::get_parents_virt( DLIList<TopologyBridge*>& parents ) 
  { CAST_LIST_TO_PARENT( myCoEdges, parents ); }
void FacetCurve::get_children_virt( DLIList<TopologyBridge*>& children ) 
  { 
    children.append( myStartPoint ); 
    if( myStartPoint != myEndPoint )
      children.append( myEndPoint );
  }
  



//-------------------------------------------------------------------------
// Purpose       : Check for G1 discontinuity
//
// Special Notes : not implemented
//
// Creator       : Steve Owen
//
// Creation Date : 07/14/00
//-------------------------------------------------------------------------
CubitBoolean FacetCurve::G1_discontinuous( 
      double param, CubitVector* mtan, CubitVector* ptan )
{ 
  DLIList<CubitPoint*> point_list;
  curveFacetEvalTool->get_points( point_list );
  CubitVector position;
  position_from_u( param, position );
  point_list.reset();
  CubitPoint* prev = point_list.get_and_step();
  for( int i = point_list.size(); i > 2; i--)
  {
    CubitPoint* point = point_list.get_and_step();
    if( (point->coordinates() - position).length_squared() < 
        (GEOMETRY_RESABS*GEOMETRY_RESABS) )
    {
      if( mtan )
      {
        *mtan = point->coordinates() - prev->coordinates();
      }
      if( ptan )
      {
        *ptan = point_list.get()->coordinates() - point->coordinates();
      }
      return CUBIT_TRUE;
    }
  }
  return CUBIT_FALSE;
}

//-------------------------------------------------------------------------
// Purpose       : retreive the facet edge list for this curve
//
// Special Notes :
//
// Creator       : Steve J. Owen
//
// Creation Date : 12/10/00
//-------------------------------------------------------------------------
void FacetCurve::get_facets(DLIList<CubitFacetEdge*>& facet_list)
{
  if (curveFacetEvalTool)
  {
    curveFacetEvalTool->get_facets(facet_list);
  }
  else
  {
    PRINT_ERROR("curve facet evaluation tool not defined for FacetCurve\n");
  }
}

//-------------------------------------------------------------------------
// Purpose       : retreive the facet point list for this curve
//
// Special Notes :
//
// Creator       : Steve J. Owen
//
// Creation Date : 12/10/00
//-------------------------------------------------------------------------
void FacetCurve::get_points(DLIList<CubitPoint*>& point_list)
{
  if (curveFacetEvalTool)
  {
    curveFacetEvalTool->get_points(point_list);
  }
  else
  {
    PRINT_ERROR("curve facet evaluation tool not defined for FacetCurve\n");
  }
}

//-------------------------------------------------------------------------
// Purpose       : set the facetLength in the CurveFacetEvalTool
//
// Special Notes :
//
// Creator       : Steve J. Owen
//
// Creation Date : 03/19/02
//-------------------------------------------------------------------------
void FacetCurve::reset_length()
{
  if (curveFacetEvalTool)
  {
    curveFacetEvalTool->set_length();
  }
  else
  {
    PRINT_ERROR("curve facet evaluation tool not defined for FacetCurve\n");
  }
}


void FacetCurve::get_lumps( DLIList<FacetLump*>& result_list )
{
  DLIList<FacetShell*> shell_list;
  get_shells( shell_list );
  shell_list.reset();
  for ( int i = shell_list.size(); i--; )
  {
    FacetShell* shell = shell_list.get_and_step();
    shell->get_lumps( result_list );
    FacetLump* lump = dynamic_cast<FacetLump*>(shell->get_lump());
    if (lump)
      result_list.append_unique(lump);
  }
}

void FacetCurve::get_shells( DLIList<FacetShell*>& result_list )
{
  DLIList<FacetSurface*> surface_list;
  DLIList<FacetShell*> temp_list;
  get_surfaces( surface_list );
  surface_list.reset();
  for ( int i = surface_list.size(); i--; )
  {
    FacetSurface* surface = surface_list.get_and_step();
    temp_list.clean_out();
    surface->get_shells( temp_list );
    result_list.merge_unique( temp_list );
  }
}

void FacetCurve::get_surfaces( DLIList<FacetSurface*>& result_list )
{
  DLIList<FacetLoop*> loop_list;
  get_loops( loop_list );
  loop_list.reset();
  for ( int i = loop_list.size(); i--; )
  {
    FacetLoop* loop = loop_list.get_and_step();
    FacetSurface* surface = dynamic_cast<FacetSurface*>(loop->get_surface());
    if (surface)
      result_list.append_unique(surface);
  }
}

void FacetCurve::get_loops( DLIList<FacetLoop*>& result_list )
{
  DLIList<FacetCoEdge*> coedge_list;
  get_coedges( coedge_list );
  coedge_list.reset();
  for ( int i = coedge_list.size(); i--; )
  {
    FacetCoEdge* coedge = coedge_list.get_and_step();
    FacetLoop* loop = dynamic_cast<FacetLoop*>(coedge->get_loop());
    if (loop)
      result_list.append_unique(loop);
  }
}

void FacetCurve::get_coedges( DLIList<FacetCoEdge*>& result_list )
{
  myCoEdges.reset();
  for ( int i = 0; i < myCoEdges.size(); i++ )
    if ( FacetCoEdge* coedge = dynamic_cast<FacetCoEdge*>(myCoEdges.next(i)) )
      result_list.append(coedge);
}

void FacetCurve::get_points( DLIList<FacetPoint*>& result_list )
{
  FacetPoint* point;
  if ( (point = dynamic_cast<FacetPoint*>(myStartPoint)) )
    result_list.append( point );
  if ( (myStartPoint != myEndPoint) &&
       (point = dynamic_cast<FacetPoint*>(myEndPoint)) )
    result_list.append( point );
}


// ********** END PUBLIC FUNCTIONS         **********

// ********** BEGIN PROTECTED FUNCTIONS    **********
// ********** END PROTECTED FUNCTIONS      **********

// ********** BEGIN PRIVATE FUNCTIONS      **********


//----------------------------------------------------------------
// Adjusts the input parameter so that it falls within the
// parameter range of this Curve, if possible.  Necessary for
// periodic curves.
//----------------------------------------------------------------
void FacetCurve::adjust_periodic_parameter(double& param)
{
    // Adjustment only legal if this is a periodic curve.
  double period;
  if ( this->is_periodic(period) && (fabs(period) > CUBIT_RESABS))
  {
    double upper_bound, lower_bound;
    this->get_param_range( lower_bound, upper_bound );
    double edge_range = upper_bound - lower_bound;
		assert( edge_range > CUBIT_RESABS * 100 );

    lower_bound -= CUBIT_RESABS;
    upper_bound += CUBIT_RESABS;
    
      // Make sure period is positive
    if (period < 0.)
      period = -period;

      // Move the parameter above the low param
    while (param < lower_bound)
      param += period;
      // Move the parameter below the high param
    while (param > upper_bound)
      param -= period;
  }
}
CubitPointContainment FacetCurve::point_containment( const CubitVector &/*point*/ )
{
   return CUBIT_PNT_UNKNOWN;
}
CubitPointContainment FacetCurve::point_containment( double /*u_param*/, 
                                                       double /*v_param*/ )
{
  return CUBIT_PNT_UNKNOWN; 
}
CubitPointContainment FacetCurve::point_containment( CubitVector &/*point*/, 
                                                       double /*u_param*/,
                                                       double /*v_param*/ )
{
   return CUBIT_PNT_UNKNOWN;
}

//-------------------------------------------------------------------------
// Purpose       : Tear down topology
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 09/29/03
//-------------------------------------------------------------------------
CubitStatus FacetCurve::disconnect_coedge( FacetCoEdge* coedge )
{
  if (!myCoEdges.move_to(coedge))
    return CUBIT_FAILURE;
  myCoEdges.remove();

  assert(coedge->curve() == this);
  coedge->remove_curve();
  
  return CUBIT_SUCCESS;
}

CubitStatus FacetCurve::get_spline_params
(
  bool &rational,    // return true/false
  int &degree,       // the degree of this spline
  DLIList<CubitVector> &cntrl_pts,  // xyz position of controlpoints
  DLIList<double> &cntrl_pt_weights, // if rational, a weight for each cntrl point.
  DLIList<double> &knots,   // There should be order+cntrl_pts.size()-2 knots
  bool &spline_is_reversed
) const
{
  PRINT_ERROR("Currently, Cubit is unable to determine spline parameters for FacetCurves.\n");
  return CUBIT_FAILURE;
}

CubitStatus FacetCurve::get_ellipse_params
(
  CubitVector &center_vec,
  CubitVector &normal,
  CubitVector &major_axis,
  double &radius_ratio
) const
{
  PRINT_ERROR("Currently, Cubit is unable to determine ellipse parameters for FacetCurves.\n");
  return CUBIT_FAILURE;
}

// ********** END PRIVATE FUNCTIONS        **********

// ********** BEGIN HELPER CLASSES         **********
// ********** END HELPER CLASSES           **********

// ********** BEGIN EXTERN FUNCTIONS       **********
// ********** END EXTERN FUNCTIONS         **********

// ********** BEGIN STATIC FUNCTIONS       **********
// ********** END STATIC FUNCTIONS         **********
