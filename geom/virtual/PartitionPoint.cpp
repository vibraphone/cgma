#include "PartitionPoint.hpp"
#include "PartitionCurve.hpp"
#include "VirtualQueryEngine.hpp"
#include "PartitionEngine.hpp"

#include "CubitPointData.hpp"
#include "TDVGFacetOwner.hpp"
#include "CubitTransformMatrix.hpp"
#include "Surface.hpp"

PartitionPoint::PartitionPoint( const CubitVector& position,
                                PartitionEntity* owner )
  : firstCurve(0), curveCount(0), myPosition( position ), facetPoint(0)
{ 
  owner->sub_entity_set().add_lower_order( this );
}
  
PartitionPoint::PartitionPoint( Point* real_point )
  : firstCurve(0), curveCount(0), 
    myPosition( real_point->coordinates() ), 
    facetPoint(0)
{
  assert( dynamic_cast<SubEntitySet*>(real_point->owner()) == 0 );
  new SubEntitySet( real_point, this );
}

PartitionPoint::PartitionPoint( CubitSimpleAttrib& attrib,
                                PartitionEntity* owner )
  : firstCurve(0), curveCount(0), facetPoint(0)
{ 
  DLIList<int> junk;
  DLIList<CubitVector*> pt;
  owner->sub_entity_set().add_lower_order( this, attrib, 0, pt, junk, junk, junk);
  assert(pt.size() == 1);
  myPosition = *pt.get();
  delete pt.get();
}


PartitionPoint::~PartitionPoint()
{
  while( firstCurve )
  {
    PartitionCurve* curve = firstCurve;
    if( curve->start_point() == this )
      curve->start_point(0);
    if( curve->end_point() == this )
      curve->end_point(0);
    assert( firstCurve != curve );
  }
  assert( curveCount == 0 );
  facet_point(0);
}

CubitVector PartitionPoint::coordinates() const
  { return myPosition; }
/*
int PartitionPoint::num_curves() const
{
  int count = 0;
  for( PartitionCurve* curve = firstCurve; curve; curve=curve->next_curve(this) )
    count++;
  return count;
}
*/
PartitionCurve* PartitionPoint::next_curve( PartitionCurve* prev ) const
  { return prev ? prev->next_curve(this) : firstCurve; }

PartitionCurve* PartitionPoint::common_curve( PartitionPoint* other ) const
{
  PartitionCurve* curve = 0;
  if( other == this )
  {
    while( (curve = next_curve(curve)) )
      if( curve->start_point() == curve->end_point() )
        return curve;
  }
  else
  {
    while( (curve = next_curve(curve)) )
      if( curve->other_point(other) )
        return curve;
  }
  return 0;
}

CubitStatus PartitionPoint::move( CubitVector& /*delta*/ )
{
  //if( curveSet.size() )
    return CUBIT_FAILURE;
  
  //myPosition += delta;
  //return CUBIT_SUCCESS;
}

void PartitionPoint::get_parents_virt( DLIList<TopologyBridge*>& parents )
{
  parents.clean_out();
  if( real_point() )
  {
      // Get all the real point's parent curves.
    sub_entity_set().get_entity()->get_parents_virt( parents );

      // Remove all real parents hidden by partitions.  We will add
      // the partitions hiding these curves later.  
    for( int i = parents.size(); i--; )
      if( dynamic_cast<SubEntitySet*>(parents.step_and_get()->owner() ) )
        parents.change_to(0);
    parents.remove_all_with_value(0);
  }
  
    // Now add any parent partition curves
  for( PartitionCurve* curve = firstCurve; curve; curve = curve->next_curve(this) )
    parents.append( curve );
    
}

void PartitionPoint::get_children_virt( DLIList<TopologyBridge*>& )
{
}

Point* PartitionPoint::real_point() const
{
  return dynamic_cast<Point*>(sub_entity_set().get_entity());
}

GeometryQueryEngine* PartitionPoint::get_geometry_query_engine() const
{
  return VirtualQueryEngine::instance();
}

void PartitionPoint::print_debug_info( const char* prefix,
                                       bool ent_set ) const
{
  if( !prefix ) prefix = "";
  char* new_prefix = new char[strlen(prefix)+3];
  strcpy( new_prefix, prefix );
  strcat( new_prefix, "  ");
  CubitVector p = coordinates();
  PRINT_INFO("%sPartitionPoint %p at (%f,%f,%f)\n", 
    prefix, static_cast<const void*>(this), p.x(), p.y(), p.z() );
  DLIList<Curve*> curve_list;
  const_cast<PartitionPoint*>(this)->TopologyBridge::curves( curve_list );
  PRINT_INFO("%s  %d Curves (%d PartitionCurves).\n", prefix, 
    curve_list.size(), num_curves() );
  
  if ( facet_point() ) {
    p = facet_point()->coordinates();
    PRINT_INFO("%s  CubitPoint %p at [%f,%f,%f] (%f)\n", prefix,
      static_cast<void*>(facet_point()), p.x(), p.y(), p.z(), 
      (coordinates() - facet_point()->coordinates()).length());
  }
  
  if( ent_set )
    sub_entity_set().print_debug_info( new_prefix );
  else 
    print_partitioned_entity(new_prefix);
  delete [] new_prefix;
}


void PartitionPoint::append_simple_attribute_virt(CubitSimpleAttrib* csa)
{ sub_entity_set().add_attribute( this, csa ); }
void PartitionPoint::remove_simple_attribute_virt(CubitSimpleAttrib* csa)
{ sub_entity_set().rem_attribute( this, csa ); }
void PartitionPoint::remove_all_simple_attribute_virt()
{ sub_entity_set().rem_all_attrib( this ); }
CubitStatus PartitionPoint::get_simple_attribute(DLIList<CubitSimpleAttrib*>& list)
{ 
  sub_entity_set().get_attributes( this, list ); 
  return CUBIT_SUCCESS;
}
CubitStatus PartitionPoint::get_simple_attribute(const CubitString& name,
                                       DLIList<CubitSimpleAttrib*>& list)
{ 
  sub_entity_set().get_attributes( this, name.c_str(), list ); 
  return CUBIT_SUCCESS;
}

void PartitionPoint::reverse_sense()
  { }

void PartitionPoint::notify_split( FacetEntity* , FacetEntity* )
  { assert(0); }

void PartitionPoint::facet_point( CubitPointData* set )
{
  if( facetPoint )
    TDVGFacetOwner::remove(facetPoint);
  facetPoint = set;
  if( set )
  {
    assert((set->coordinates() - coordinates()).length_squared() 
               < GEOMETRY_RESABS*GEOMETRY_RESABS);
    TDVGFacetOwner::set(facetPoint,this);
  }
}

CubitStatus PartitionPoint::save( CubitSimpleAttrib& attrib )
{
  int id = sub_entity_set().get_id(this);
  if( id <= 0 )
    return CUBIT_FAILURE;
  
  DLIList<CubitVector*> pt_list(1);
  pt_list.append( new CubitVector(coordinates()) );
  
  return sub_entity_set().save_geometry( id, 0, &pt_list, 0, 0, 0, attrib );
}

CubitBox PartitionPoint::bounding_box() const
  { return CubitBox(coordinates()); }

CubitStatus PartitionPoint::move_to_geometry( CubitVector& pos )
{
  pos = coordinates();
  return CUBIT_SUCCESS;
}

void PartitionPoint::transform( const CubitTransformMatrix& xform )
{
  if( Point* point = dynamic_cast<Point*>(partitioned_entity()) )
    myPosition = point->coordinates();
  else
    myPosition = xform * myPosition;
  
  if( facetPoint )
    facetPoint->set(myPosition);
}
