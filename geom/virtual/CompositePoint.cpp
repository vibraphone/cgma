//-------------------------------------------------------------------------
// Filename      : CompositePoint.cpp
//
// Purpose       : Point decorator for composite topology
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 03/11/02
//-------------------------------------------------------------------------

#include "CompositePoint.hpp"
#include "CompositeEngine.hpp"
#include "VirtualQueryEngine.hpp"

CompositePoint::CompositePoint( TBPoint* real_pt )
  : firstCurve(0), realPoint(real_pt), stitchNext(0), HadBridgeRemoved(0)
{
  if( real_pt->owner() )
    real_pt->owner()->swap_bridge( real_pt, this, false );
  real_pt->owner(this);
}

CompositePoint::~CompositePoint()
{
    // remove from all curves
  while( firstCurve )
  {
    CompositeCurve* curve = firstCurve;
    if( curve->start_point() == this )
      curve->start_point(0);
    if( curve->end_point() == this )
      curve->end_point(0);
    assert(firstCurve != curve);
  }

  if( stitchNext )
  {
      // unmerge
    unstitch_all();
  }
  
  if (realPoint)
  {
      // update owner
    assert(realPoint->owner() == this);
    realPoint->owner(0);
    if( owner() )  
      owner()->swap_bridge( this, realPoint, false );
    
    realPoint = 0;
  }
}

void CompositePoint::unstitch_all()
{
  assert( !firstCurve );
  while (stitchNext)
  {
    stitchNext->owner(0);
    if (owner())
      owner()->notify_copied( stitchNext, this );
    stitchNext = stitchNext->stitchNext;
  }
}


void CompositePoint::get_parents_virt( DLIList<TopologyBridge*>& list )
{
  DLIList<TopologyBridge*> point_parents;
  realPoint->get_parents_virt( point_parents );
  for( int i = point_parents.size(); i--; )
  {
    TopologyBridge* tb = point_parents.get_and_step();
    if( CompositeCurve* curve = dynamic_cast<CompositeCurve*>(tb->owner()) )
    {
      if( ! dynamic_cast<HiddenEntitySet*>(curve->owner()) )
        list.append_unique( curve );
    }
    else if( ! dynamic_cast<HiddenEntitySet*>(tb->owner()) )
      list.append( tb );
  }

    // get point-curves also
  CompositeCurve* curve = 0;
  while ((curve = next_curve(curve)))
    if (curve->num_curves() == 0)
      list.append(curve);
  
  if (stitchNext)
  {
    point_parents.clean_out();
    stitchNext->get_parents_virt( point_parents );
    point_parents.reset();
    for (int j = point_parents.size(); j--; )
    {
      TopologyBridge* bridge = point_parents.get_and_step();
      CompositeCurve* curv = dynamic_cast<CompositeCurve*>(bridge);
      if (curv)
        list.append_unique( curv->primary_stitched_curve() );
      else
        list.append_unique( bridge );
    }
  }
}

void CompositePoint::get_children_virt( DLIList<TopologyBridge*>& )
{
}


CubitStatus CompositePoint::remove_bridge( TopologyBridge* bridge )
{
  if( bridge->owner() != this  )
  {
    assert(0);
    return CUBIT_FAILURE;
  }
  
  if( bridge == realPoint )
  {
    realPoint = 0;
  }
  else
  {
    return CUBIT_FAILURE;
  }
  
  bridge->owner(0);
  
  if (!realPoint)
    CompositeEngine::instance().notify_deactivated(this);

  HadBridgeRemoved = 1;

  return CUBIT_SUCCESS;
}

CubitStatus CompositePoint::swap_bridge( TopologyBridge* oldtb, 
                                         TopologyBridge* newtb, 
                                         bool )
{
  TBPoint* oldpt = dynamic_cast<TBPoint*>(oldtb);
  TBPoint* newpt = dynamic_cast<TBPoint*>(newtb);
  if( !(oldpt && newpt) || newpt->owner() )
    return CUBIT_FAILURE;

  assert(oldpt == realPoint );

  realPoint = newpt;

  newtb->owner(this);
  oldtb->owner(0);
  return CUBIT_SUCCESS;
}


void CompositePoint::print_debug_info( const char* prefix, bool brief ) const
{
  if( prefix == 0 ) prefix = "";

#ifdef TOPOLOGY_BRIDGE_IDS
  PRINT_INFO("%sCompositePoint %d : %s %d\n", prefix, get_id(),
    realPoint ? fix_type_name(typeid(*realPoint).name()) : "NO REAL POINT", 
    realPoint ? realPoint->get_id() : 0 );
#else  
  PRINT_INFO("%sCompositePoint %p : %s %p\n", prefix, this,
    realPoint ? fix_type_name(typeid(*realPoint).name()) : "NO REAL POINT", 
    realPoint);
#endif
  
  if ( !brief )
  {  
    CubitVector p = coordinates();
    PRINT_INFO("%s  (%f,%f,%f)\n", prefix, p.x(), p.y(), p.z() );
  }
}

void CompositePoint::notify_reversed( TopologyBridge* )
  { }

CubitStatus CompositePoint::stitch( CompositePoint* point )
{
  if( !point || point->owner() != this->owner() )
  {
    assert(0);
    return CUBIT_FAILURE;
  }
  
  if( point->owner() )
    point->owner()->notify_merged( point, this );
  
  if( point->owner() )
    point->owner()->remove_bridge( point );

  point->owner( this );
  CompositePoint* end = point;
  while (end->stitchNext)
  {
    end->stitchNext->owner( this );
    end = end->stitchNext;
  }
  
  end->stitchNext = stitchNext;
  stitchNext = end;
  
  return CUBIT_SUCCESS;
}

void CompositePoint::get_stitched( DLIList<CompositePoint*>& result )
{
  CompositePoint* pt = dynamic_cast<CompositePoint*>(owner());
  if (!pt)
    pt = this;
  while (pt)
  {
    result.append(pt);
    pt = pt->stitchNext;
  }
}

GeometryQueryEngine* CompositePoint::get_geometry_query_engine() const
{ return VirtualQueryEngine::instance(); }
