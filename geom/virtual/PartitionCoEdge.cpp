#include "PartitionCoEdge.hpp"
#include "PartitionSurface.hpp"
#include "PartitionCurve.hpp"
#include "PartitionLoop.hpp"
#include "PartitionEngine.hpp"
#include "VirtualQueryEngine.hpp"

PartitionCoEdge::PartitionCoEdge(PartitionSurface* surf, CubitSense sense)
  : mySense(sense),
    myLoop(0), 
    loopPrev(0),
    loopNext(0),
    myCurve(0), 
    curveNext(0)
{ 
  surf->sub_entity_set().add_lower_order( this );
}

PartitionCoEdge::PartitionCoEdge(CoEdgeSM* coedge)
  : myLoop(0), 
    loopPrev(0),
    loopNext(0),
    myCurve(0), 
    curveNext(0)
{
  mySense = coedge->sense();
  assert( dynamic_cast<SubEntitySet*>(coedge->owner()) == 0 );
  new SubEntitySet( coedge, this );
}

PartitionCoEdge::PartitionCoEdge(PartitionCoEdge* split_from)
  : myLoop(0), 
    loopPrev(0),
    loopNext(0),
    myCurve(0), 
    curveNext(0)
{
  mySense = split_from->mySense;
  if( split_from->real_coedge() )
  {
    split_from->sub_entity_set().add_partition( this, split_from );
  }
  else
  {
    split_from->sub_entity_set().add_lower_order( this );
  }
}

PartitionCoEdge::~PartitionCoEdge()
{
  if( myCurve )
    myCurve->remove(this);
  if( myLoop )
    myLoop->remove(this);
  
  assert( !myCurve && !myLoop );
}

void PartitionCoEdge::reverse_sense()
{
  mySense = mySense == CUBIT_FORWARD  ? CUBIT_REVERSED :
            mySense == CUBIT_REVERSED ? CUBIT_FORWARD  :
                                        CUBIT_UNKNOWN  ;
} 

TopologyBridge* PartitionCoEdge::find_parent_loop() const
{
  if( get_loop() )
    return get_loop();
  
  CoEdgeSM* coedge = real_coedge();
  if( !coedge )
    return 0;
  
  DLIList<TopologyBridge*> list;
  coedge->get_parents_virt( list );
  assert( list.size() == 1 );
  return list.get();
}

CoEdgeSM* PartitionCoEdge::real_coedge() const
{
  return dynamic_cast<CoEdgeSM*>(partitioned_entity());
    // this will be null for coedges of segmented (split) curves.
}

void PartitionCoEdge::get_children_virt( DLIList<TopologyBridge*>& list )
{
  assert( myCurve != NULL );
  list.append( myCurve );
}


void PartitionCoEdge::get_parents_virt( DLIList<TopologyBridge*>& list )
{
  TopologyBridge* result = find_parent_loop();
  assert( result != NULL );
  list.append( result );
}


    
void PartitionCoEdge::append_simple_attribute_virt( const CubitSimpleAttrib& )
  { }
void PartitionCoEdge::remove_simple_attribute_virt( const CubitSimpleAttrib& )
  { }
void PartitionCoEdge::remove_all_simple_attribute_virt()
  { }
CubitStatus PartitionCoEdge::get_simple_attribute( DLIList<CubitSimpleAttrib>& )
  { return CUBIT_FAILURE; }
CubitStatus PartitionCoEdge::get_simple_attribute(const CubitString& ,
                                       DLIList<CubitSimpleAttrib>& )
  { return CUBIT_FAILURE; }



GeometryQueryEngine* PartitionCoEdge::get_geometry_query_engine() const
{
  return VirtualQueryEngine::instance();
}

PartitionPoint* PartitionCoEdge::start_point() const
{
  PartitionPoint* result = 0;
  if( myCurve )
  {
    if( mySense == CUBIT_FORWARD )
      result = myCurve->start_point();
    else if( mySense == CUBIT_REVERSED )
      result = myCurve->end_point();
  }
  return result;
}

PartitionPoint* PartitionCoEdge::end_point() const
{
  PartitionPoint* result = 0;
  if( myCurve )
  {
    if( mySense == CUBIT_FORWARD )
      result = myCurve->end_point();
    else if( mySense == CUBIT_REVERSED )
      result = myCurve->start_point();
  }
  return result;
}

void PartitionCoEdge::notify_split( FacetEntity*, FacetEntity* )
  { assert(0); }

CubitBox PartitionCoEdge::bounding_box() const
  { return ((Curve*)get_curve())->bounding_box(); }

