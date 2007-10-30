//-------------------------------------------------------------------------
// Filename      : OCCCurve.cpp
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
#include "config.h"
// ********** BEGIN CUBIT INCLUDES         **********

#include "CastTo.hpp"
#include "CubitVector.hpp"
#include "CubitBox.hpp"
#include "GeometryDefines.h"
#include "OCCCurve.hpp"
#include "OCCAttrib.hpp"
#include "GeometryQueryEngine.hpp"
#include "OCCQueryEngine.hpp"
#include "CoEdgeSM.hpp"
#include "CubitPoint.hpp"

#include "OCCBody.hpp"
#include "OCCLump.hpp"
#include "OCCShell.hpp"
#include "OCCSurface.hpp"
#include "OCCLoop.hpp"
#include "OCCCoEdge.hpp"
#include "OCCPoint.hpp"

#include <BRepAdaptor_Curve.hxx>
#include <TopExp.hxx>
#include <TopTools_IndexedMapOfShape.hxx>
#include <GCPnts_AbscissaPoint.hxx>
#include <Bnd_Box.hxx>
#include <BndLib_Add3dCurve.hxx>
#include <Precision.hxx>
#include <Extrema_ExtPC.hxx>
#include <BRepLProp_CLProps.hxx>
#include <BRep_Tool.hxx>
#include <TopoDS.hxx>
#include "GeomLProp_CurveTool.hxx"
#include "GeomAPI_ExtremaCurveCurve.hxx"
#include "Geom_Line.hxx"
#include "GC_MakeLine.hxx"
//#include "TopOpeBRep_ShapeIntersector.hxx"
//#include "TopOpeBRep_Point2d.hxx"
//#include "TopOpeBRep_EdgesIntersector.hxx"
//#include "TopOpeBRepTool_ShapeTool.hxx"
//#include "BRepPrimAPI_MakePrism.hxx"
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
OCCCurve::OCCCurve( TopoDS_Edge *theEdge )
{
  myTopoDSEdge = theEdge;
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
OCCCurve::~OCCCurve() 
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
// Creator       : Steve Owen
//
// Creation Date : 07/14/00
//-------------------------------------------------------------------------
void OCCCurve::append_simple_attribute_virt(CubitSimpleAttrib *csa)
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
void OCCCurve::remove_simple_attribute_virt(CubitSimpleAttrib *csa)
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
void OCCCurve::remove_all_simple_attribute_virt()
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
CubitStatus OCCCurve::get_simple_attribute(DLIList<CubitSimpleAttrib*>&
                                               csa_list)
  { return attribSet.get_attributes(csa_list); }
  
CubitStatus OCCCurve::get_simple_attribute( const CubitString& name,
                                      DLIList<CubitSimpleAttrib*>& csa_list)
  { return attribSet.get_attributes( name, csa_list ); }


CubitStatus OCCCurve::save_attribs( FILE *file_ptr )
  { return attribSet.save_attributes(file_ptr); }

CubitStatus OCCCurve::restore_attribs( FILE *file_ptr, unsigned int endian )
  { return attribSet.restore_attributes(file_ptr, endian); }



//-------------------------------------------------------------------------
// Purpose       : Get geometry modeling engine: OCCQueryEngine
//
// Special Notes :
//
// Creator       : Steve Owen
//
// Creation Date : 07/14/00
//-------------------------------------------------------------------------
GeometryQueryEngine* OCCCurve::get_geometry_query_engine() const
{
  return OCCQueryEngine::instance();
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
CubitBox OCCCurve::bounding_box() const 
{
  BRepAdaptor_Curve acurve(*myTopoDSEdge);
  Bnd_Box aBox;
  BndLib_Add3dCurve::Add(acurve, Precision::Approximation(), aBox);
  double min[3], max[3];
  aBox.Get( min[0], min[1], min[2], max[0], max[1], max[2]);
  return CubitBox(min, max);
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
double OCCCurve::measure()
{
  BRepAdaptor_Curve acurve(*myTopoDSEdge);
  return GCPnts_AbscissaPoint::Length(acurve);
}

//-------------------------------------------------------------------------
// Purpose       : Return the arc length along the Curve starting from
//                 the point represented by the parameter1 going to the 
//                 point represented by parameter2.
//
// Special Notes : The sign of the returned length value is always positive.
//                 Parameter1 and parameter2 are with respect to the EDGE.
//
// Creator       : Steve Owen
//
// Creation Date : 07/14/00
//-------------------------------------------------------------------------
double OCCCurve::length_from_u( double parameter1, double parameter2 )
{
  BRepAdaptor_Curve acurve(*myTopoDSEdge);
  return GCPnts_AbscissaPoint::Length(acurve, parameter1, parameter2);
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
CubitBoolean OCCCurve::is_periodic(double& period)
{
  BRepAdaptor_Curve acurve(*myTopoDSEdge);
  if (acurve.IsPeriodic())
  {
    period = acurve.Period();
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
CubitBoolean OCCCurve::get_param_range( double& lower_bound,
                                          double& upper_bound )
{
  BRepAdaptor_Curve acurve(*myTopoDSEdge);
  lower_bound = acurve.FirstParameter();
  upper_bound = acurve.LastParameter();
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
CubitStatus OCCCurve::get_interior_extrema(
  DLIList<CubitVector*>& interior_points,
  CubitSense& return_sense )
{
  // Danilov: try to use GeomAPI_ExtremaCurveCurve
  // Will do 3 primary directions seperately. 
  assert(0);
  
  CubitVector x(1.0, 0.0, 0.0);
  get_interior_extrema_in_direction(interior_points, x);
  CubitVector y(0.0, 1.0, 0.0);
  get_interior_extrema_in_direction(interior_points, y);
  CubitVector z(0.0, 0.0, 1.0);
  get_interior_extrema_in_direction(interior_points, z );
 
  // Return sense is whatever the sense of this curve is.
  TopAbs_Orientation sense = myTopoDSEdge->Orientation();
  return_sense = (sense == TopAbs_FORWARD ? CUBIT_FORWARD : CUBIT_REVERSED);
 
  return CUBIT_SUCCESS;
}

CubitStatus OCCCurve::get_interior_extrema_in_direction(
		 	DLIList<CubitVector*>& interior_points,
			CubitVector dir)
{
  //Create a straight line.
  gp_Pnt origin(0.0, 0.0, 0.0);
  gp_Dir adir(dir.x(), dir.y(), dir.z());
  Handle(Geom_Line) line = GC_MakeLine(origin, adir);

  //get the Geom_Curve of the OCCCurve
  Standard_Real first;
  Standard_Real last;
  Handle(Geom_Curve) myCurve = BRep_Tool::Curve(*myTopoDSEdge, first, last);
  
  GeomAPI_ExtremaCurveCurve extrema(myCurve, line);
  int nPnt = extrema.NbExtrema();
  for (int i = 1; i <= nPnt ; i++)
  {
    gp_Pnt P1, P2;
    extrema.Points(i, P1, P2);
    CubitVector* v = new CubitVector(P1.X(), P1.Y(), P1.Z());
    interior_points.append(v);
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
CubitStatus OCCCurve::closest_point( 
  CubitVector const& location, 
  CubitVector& closest_location,
  CubitVector* tangent_ptr,
  CubitVector* curvature_ptr,
  double* param)
{  
  BRepAdaptor_Curve acurve(*myTopoDSEdge);
  gp_Pnt p(location.x(), location.y(), location.z()), newP(0.0, 0.0, 0.0);
  double newVal=0.0;
  int i;
  BRepLProp_CLProps CLP(acurve, 2, Precision::PConfusion());
  Extrema_ExtPC ext(p, acurve, Precision::Approximation());
  if (ext.IsDone() && (ext.NbExt() > 0)) {
	  for ( i = 1 ; i <= ext.NbExt() ; i++ ) {
		  if ( ext.IsMin(i) ) {
			  newVal = ext.Point(i).Parameter();
			  newP = ext.Point(i).Value();
			  CLP.SetParameter(newVal);
		  }
	  }
  }
  closest_location = CubitVector(newP.X(), newP.Y(), newP.Z());
  if (tangent_ptr != NULL) {
	  gp_Dir tangent;
	  if (CLP.IsTangentDefined()) {
		  CLP.Tangent(tangent);
		  *tangent_ptr = CubitVector(tangent.X(), tangent.Y(), tangent.Z()); 
	  }
  }
  if (curvature_ptr != NULL) {
	  double curvature = CLP.Curvature();
	  // Danilov: confused here
  }
  if (param != NULL) {
	  *param = newVal;
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
CubitStatus OCCCurve::position_from_u (double u_value,
                                        CubitVector& output_position)
{
  BRepAdaptor_Curve acurve(*myTopoDSEdge);
  gp_Pnt p = acurve.Value(u_value);
  output_position.x(p.X());
  output_position.y(p.Y());
  output_position.z(p.Z());
  return CUBIT_SUCCESS;
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
double OCCCurve::u_from_position (const CubitVector& input_position)
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
double OCCCurve::u_from_arc_length ( double root_param,
                                       double arc_length )
{
  BRepAdaptor_Curve acurve(*myTopoDSEdge);
  GCPnts_AbscissaPoint abs(acurve, arc_length, root_param);
  if (abs.IsDone()) return abs.Parameter();
  else return 0.0;
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
CubitBoolean OCCCurve::is_position_on( const CubitVector &test_position )
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
GeometryType OCCCurve::geometry_type()
{
  BRepAdaptor_Curve acurve(*myTopoDSEdge);
  if (acurve.GetType() == GeomAbs_BezierCurve)
     return BEZIER_CURVE_TYPE;
  if (acurve.GetType() == GeomAbs_BSplineCurve)
     return SPLINE_CURVE_TYPE;
  if (acurve.GetType() == GeomAbs_Line)
     return STRAIGHT_CURVE_TYPE;
  if (acurve.GetType() == GeomAbs_Parabola)
     return PARABOLA_CURVE_TYPE;
  if (acurve.GetType() == GeomAbs_Hyperbola)
     return HYPERBOLA_CURVE_TYPE;
  if (acurve.GetType() == GeomAbs_Circle)
      return ARC_CURVE_TYPE;
  if (acurve.GetType() == GeomAbs_Ellipse)
     return ELLIPSE_CURVE_TYPE;
  return UNDEFINED_CURVE_TYPE;
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
CubitStatus OCCCurve::get_point_direction( CubitVector& point, 
                                             CubitVector& direction )
{
  point = point;
  direction = direction;
  PRINT_DEBUG_122("OCCCurve::get_point_direction currently not implemented.\n");
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
CubitStatus OCCCurve::get_center_radius( CubitVector& center, 
                                           double& radius )
{
  if( geometry_type() != ELLIPSE_CURVE_TYPE &&
      geometry_type() != ARC_CURVE_TYPE )
    return CUBIT_FAILURE;
  
  center = center;
  radius = radius;
  PRINT_DEBUG_122("OCCCurve::get_center_radius currently not implemented.\n");
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
double OCCCurve::start_param()
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
double OCCCurve::end_param()
{
   double start = 0.0, end = 0.0;
   
   get_param_range( start, end );
   return end;
}

/*
void OCCCurve::bodysms(DLIList<BodySM*> &bodies)
{
  int ii;
  for ( ii = myCoEdges.size(); ii > 0; ii-- )
  {
    myCoEdges.get_and_step()->bodysms(bodies);
  }
}

void OCCCurve::lumps(DLIList<Lump*> &lumps)
{
  int ii;
  for ( ii = myCoEdges.size(); ii > 0; ii-- )
  {
    myCoEdges.get_and_step()->lumps(lumps);
  }
}

void OCCCurve::shellsms(DLIList<ShellSM*> &shellsms)
{
  int ii;
  for ( ii = myCoEdges.size(); ii > 0; ii-- )
  {
    myCoEdges.get_and_step()->shellsms(shellsms);
  }
}

void OCCCurve::surfaces(DLIList<Surface*> &surfaces)
{
  int ii;
  for ( ii = myCoEdges.size(); ii > 0; ii-- )
  {
    myCoEdges.get_and_step()->surfaces(surfaces);
  }
}

void OCCCurve::loopsms(DLIList<LoopSM*> &loopsms)
{
  int ii; 
  for ( ii = myCoEdges.size(); ii > 0; ii-- )
  {
    myCoEdges.get_and_step()->loopsms(loopsms);
  }
}


void OCCCurve::coedgesms(DLIList<CoEdgeSM*> &coedgesms)
{
  int ii; 
  for ( ii = myCoEdges.size(); ii > 0; ii-- )
  {
    coedgesms.append_unique( myCoEdges.get_and_step() );
  } 
}

void OCCCurve::curves(DLIList<Curve*> &curves)
{
  curves.append_unique( this );
}

void OCCCurve::points(DLIList<Point*> &points)
{
  points.append_unique( myStartPoint );
  points.append_unique( myEndPoint );
}
*/


void OCCCurve::get_parents_virt( DLIList<TopologyBridge*>& parents ) 
  { /*CAST_LIST_TO_PARENT( myCoEdges, parents ); */}
void OCCCurve::get_children_virt( DLIList<TopologyBridge*>& children ) 
{
	TopTools_IndexedMapOfShape M;
	TopExp::MapShapes(*myTopoDSEdge, TopAbs_VERTEX, M);
	TopologyBridge *point1, *point2;
	if (M.Extent()==1) {
		point1 = OCCQueryEngine::occ_to_cgm(M(1));
		children.append_unique(point1);
	} else if (M.Extent()==2) {
		if (  fabs(BRep_Tool::Parameter(TopoDS::Vertex(M(1)), *myTopoDSEdge)-start_param()) > 
				fabs(BRep_Tool::Parameter(TopoDS::Vertex(M(2)), *myTopoDSEdge)-start_param())  ) {
			point1 = OCCQueryEngine::occ_to_cgm(M(2));
			point2 = OCCQueryEngine::occ_to_cgm(M(1));
		} else {
			point1 = OCCQueryEngine::occ_to_cgm(M(1));
			point2 = OCCQueryEngine::occ_to_cgm(M(2));
		}
		if (point1 == point2) {
			children.append_unique(point1);
		} else {
			children.append_unique(point1);
			children.append_unique(point2);
		}
	}
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
CubitBoolean OCCCurve::G1_discontinuous( 
      double param, CubitVector* mtan, CubitVector* ptan )
{ 
  assert(0);
  return CUBIT_FALSE;
}

void OCCCurve::get_lumps( DLIList<OCCLump*>& result_list )
{
  DLIList<OCCShell*> shell_list;
  get_shells( shell_list );
  shell_list.reset();
  for ( int i = shell_list.size(); i--; )
  {
    OCCShell* shell = shell_list.get_and_step();
    shell->get_lumps( result_list );
    OCCLump* lump = dynamic_cast<OCCLump*>(shell->get_lump());
    if (lump)
      result_list.append_unique(lump);
  }
}

void OCCCurve::get_shells( DLIList<OCCShell*>& result_list )
{
  DLIList<OCCSurface*> surface_list;
  DLIList<OCCShell*> temp_list;
  get_surfaces( surface_list );
  surface_list.reset();
  for ( int i = surface_list.size(); i--; )
  {
    OCCSurface* surface = surface_list.get_and_step();
    temp_list.clean_out();
    surface->get_shells( temp_list );
    result_list.merge_unique( temp_list );
  }
}

void OCCCurve::get_surfaces( DLIList<OCCSurface*>& result_list )
{
  DLIList<OCCLoop*> loop_list;
  get_loops( loop_list );
  loop_list.reset();
  for ( int i = loop_list.size(); i--; )
  {
    OCCLoop* loop = loop_list.get_and_step();
    OCCSurface* surface = dynamic_cast<OCCSurface*>(loop->get_surface());
    if (surface)
      result_list.append_unique(surface);
  }
}

void OCCCurve::get_loops( DLIList<OCCLoop*>& result_list )
{
  DLIList<OCCCoEdge*> coedge_list;
  get_coedges( coedge_list );
  coedge_list.reset();
  for ( int i = coedge_list.size(); i--; )
  {
    OCCCoEdge* coedge = coedge_list.get_and_step();
    OCCLoop* loop = dynamic_cast<OCCLoop*>(coedge->get_loop());
    if (loop)
      result_list.append_unique(loop);
  }
}

void OCCCurve::get_points( DLIList<OCCPoint*>& result_list )
{
  TopTools_IndexedMapOfShape M;
  TopExp::MapShapes(*myTopoDSEdge, TopAbs_VERTEX, M);
  int ii;
  for (ii=M.Extent(); ii>0; ii--) {
	  TopologyBridge *point = OCCQueryEngine::occ_to_cgm(M(ii));
	  result_list.append_unique(dynamic_cast<OCCPoint*>(point));
  }
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
void OCCCurve::adjust_periodic_parameter(double& param)
{
  assert(0);
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
CubitPointContainment OCCCurve::point_containment( const CubitVector &/*point*/ )
{
   return CUBIT_PNT_UNKNOWN;
}


// ********** END PRIVATE FUNCTIONS        **********

// ********** BEGIN HELPER CLASSES         **********
// ********** END HELPER CLASSES           **********

// ********** BEGIN EXTERN FUNCTIONS       **********
// ********** END EXTERN FUNCTIONS         **********

// ********** BEGIN STATIC FUNCTIONS       **********
// ********** END STATIC FUNCTIONS         **********
