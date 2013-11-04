#include "PartPTCurve.hpp"
#include "PartitionSurface.hpp"
#include "PartitionPoint.hpp"
#include "GMem.hpp"
#include "PartitionEngine.hpp"

PartPTCurve::PartPTCurve( PartitionSurface* owner )
{
  owner->sub_entity_set().add_lower_order( this );
}

PartPTCurve* PartPTCurve::construct( const CubitSimpleAttrib& attrib,
                                     PartitionSurface* owner )
{
  DLIList<int> vertex_conn;
  PartPTCurve* result = new PartPTCurve( owner, attrib, vertex_conn );
  
  if ( vertex_conn.size() != 4 )
  {
    delete result;
    return 0;
  }
  
  PartitionPoint *start = 0, *end = 0;
  PartitionEntity* ent;
  vertex_conn.reset();
  int set_id = vertex_conn.get_and_step();
  int ent_id = vertex_conn.get_and_step();
  ent = PartitionEngine::instance().entity_from_id(set_id,ent_id,owner->sub_entity_set());
  start = dynamic_cast<PartitionPoint*>(ent);
  set_id = vertex_conn.get_and_step();
  ent_id = vertex_conn.get_and_step();
  ent = PartitionEngine::instance().entity_from_id(set_id,ent_id,owner->sub_entity_set());
  end = dynamic_cast<PartitionPoint*>(ent);
  
  if ( !start || !end || start != end )
  {
    delete result;
    return 0;
  }
  
  result->start_point( end );
  result->end_point( end );
  return result;
}

PartPTCurve::PartPTCurve( PartitionSurface* surface,
                          const CubitSimpleAttrib& attrib,
                          DLIList<int>& vertex_conn )
{
  DLIList<CubitVector*> points;
  DLIList<int> junk;
  surface->sub_entity_set().add_lower_order( this, attrib, 1, points, junk, vertex_conn, junk );
  assert( points.size() == 0 );
}

PartPTCurve::~PartPTCurve()
  {}

PartitionCurve* PartPTCurve::split( double )
  { assert(0); return 0; }
  
CubitStatus PartPTCurve::combine( PartitionCurve* curve )
{
  assert( dynamic_cast<PartPTCurve*>(curve) != NULL );
  assert( 0 );
  return CUBIT_FAILURE;
}

CubitStatus PartPTCurve::get_graphics( GMem& result, 
                                    double /*angle_tolerance*/,
                                    double /*distance_tolerance*/,
                                    double /*max_edge_length*/) 
{
  result.pointListCount = 0;
  return CUBIT_SUCCESS;
}


void PartPTCurve::append_simple_attribute_virt(const CubitSimpleAttrib& csa)
  { sub_entity_set().add_attribute( this, csa ); }
void PartPTCurve::remove_simple_attribute_virt(const CubitSimpleAttrib& csa)
  { sub_entity_set().rem_attribute( this, csa ); }
void PartPTCurve::remove_all_simple_attribute_virt()
  { sub_entity_set().rem_all_attrib( this ); }
CubitStatus PartPTCurve::get_simple_attribute(DLIList<CubitSimpleAttrib>& list)
{ 
  sub_entity_set().get_attributes( this, list ); 
  return CUBIT_SUCCESS;
}
CubitStatus PartPTCurve::get_simple_attribute(const CubitString& name,
                                       DLIList<CubitSimpleAttrib>& list)
{ 
  sub_entity_set().get_attributes( this, name.c_str(), list ); 
  return CUBIT_SUCCESS;
}

CubitBox PartPTCurve::bounding_box() const
  { return CubitBox( coordinates() ); }
double PartPTCurve::measure()
  { return 0.0; }
GeometryType PartPTCurve::geometry_type()
  { return POINT_CURVE_TYPE; }

CubitStatus PartPTCurve::closest_point( CubitVector const& ,
                                        CubitVector& closest,
                                        CubitVector* tangent,
                                        CubitVector* curvature,
                                        double* param )
{
  closest = coordinates();
  if( tangent ) tangent->set(0.,0.,0.);
  if( curvature) curvature->set(0.,0.,0.);
  if( param ) *param = 0.0;
  return CUBIT_SUCCESS;
}

CubitPointContainment PartPTCurve::point_containment( const CubitVector& position )
  { return is_position_on(position) ? CUBIT_PNT_BOUNDARY : CUBIT_PNT_OFF; }

CubitBoolean  PartPTCurve::is_position_on( const CubitVector& position )
{
  const double tolsqr = GEOMETRY_RESABS * GEOMETRY_RESABS;
  return (position - coordinates()).length_squared() < tolsqr;
}

CubitBoolean PartPTCurve::G1_discontinuous( double, CubitVector* min, CubitVector* plu )
{
  if( min ) min->set( 0., 0., 0. );
  if( plu ) plu->set( 0., 0., 0. );
  return CUBIT_TRUE;
}

CubitStatus PartPTCurve::get_interior_extrema( DLIList<CubitVector*>&, CubitSense& )
  { return CUBIT_SUCCESS; }

CubitStatus PartPTCurve::get_center_radius( CubitVector&, double& )
  { return CUBIT_FAILURE; }

CubitBoolean PartPTCurve::get_param_range( double& lower, double& upper )
  { lower = upper = 0.0; return CUBIT_TRUE; }

double PartPTCurve::start_param()
  { return 0.0; }
double PartPTCurve::end_param()
  { return 0.0; }

CubitBoolean PartPTCurve::is_periodic( double& period )
  { period = 0.0; return CUBIT_TRUE; }

double PartPTCurve::length_from_u( double, double )
  { return  0.0; }

double PartPTCurve::u_from_position( const CubitVector& )
  { return 0.0; }

CubitStatus PartPTCurve::position_from_u( double, CubitVector& result )
  { result = coordinates(); return CUBIT_SUCCESS; }

double PartPTCurve::u_from_arc_length( double, double )
  { return 0.0; }

CubitVector PartPTCurve::coordinates() const
{
  if( start_point() )
    return start_point()->coordinates();
  
  if( end_point() )
    return end_point()->coordinates();
  
  return CubitVector(0.,0.,0.);
}

void PartPTCurve::reverse_sense()
  { }
  
CubitStatus PartPTCurve::save( CubitSimpleAttrib& attrib )
{
  int id = sub_entity_set().get_id(this);
  if( id <= 0 ) return CUBIT_FAILURE;
  
  DLIList<int> end_points(4);
  get_save_topology(end_points);

  return sub_entity_set().save_geometry( id, 1, 0, 0, &end_points, 0, attrib );
}

CubitStatus PartPTCurve::get_spline_params
(
  bool &rational,    // return true/false
  int &degree,       // the degree of this spline
  DLIList<CubitVector> &cntrl_pts,  // xyz position of controlpoints
  DLIList<double> &cntrl_pt_weights, // if rational, a weight for each cntrl point.
  DLIList<double> &knots,   // There should be order+cntrl_pts.size()-2 knots
  bool &spline_is_reversed
) const
{
  PRINT_ERROR("Currently, Cubit is unable to determine spline parameters for PartPTCurves.\n");
  return CUBIT_FAILURE;
}

CubitStatus PartPTCurve::get_ellipse_params
(
  CubitVector &center_vec,
  CubitVector &normal,
  CubitVector &major_axis,
  double &radius_ratio
) const
{
  PRINT_ERROR("Currently, Cubit is unable to determine ellipse parameters for PartPTCurves.\n");
  return CUBIT_FAILURE;
}
