//-------------------------------------------------------------------------
// Filename      : HiddenEntitySet.cpp
//
// Purpose       : A class to hold a list of hidden entites and act as
//                 their owner.
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 01/14/02
//-------------------------------------------------------------------------

#include "VGDefines.h"
#include "HiddenEntitySet.hpp"

#include "BodySM.hpp"
#include "Lump.hpp"
#include "ShellSM.hpp"
#include "Surface.hpp"
#include "LoopSM.hpp"
#include "CoEdgeSM.hpp"
#include "Curve.hpp"
#include "Point.hpp"

// for debug output
#include "CompositePoint.hpp"
#include "CompositeCurve.hpp"
#include "CompositeCoEdge.hpp"
#include "CompositeSurface.hpp"


HiddenEntitySet::~HiddenEntitySet()
{
  std::vector<TopologyBridge*>::iterator iter;
  for (iter=hiddenList.begin(); iter!=hiddenList.end(); iter++)
    (*iter)->owner(0);

  hiddenList.clear();
}

CubitStatus HiddenEntitySet::hide( TopologyBridge* bridge )
{
  if( bridge->owner() && !bridge->owner()->remove_bridge(bridge) )
    return CUBIT_FAILURE;
    
  bridge->owner( this );

  // It would be simpler to just add new hidden entities at the back of the list
  // but this code mimics the previous behavior of inserting new entities 
  // after the first entity in the list.  Changing this changes the ids of
  // geometry entities when a composite is removed.  BWH - 11/09/05
  //
  // TODO - BWH - do we want to change the order now to make things simpler or come up
  // with a way to keep ids after removing composites
  //
  // TODO - BWH - the list functions in HiddenEntitySet are almost identical to the list functions in BridgeManager
  // the main difference is the order that entities are inserted.  If we want to make the behavior consistent, maybe a TBList class
  // could be implemented to use in both places.

  std::vector<TopologyBridge*>::iterator iter;
  iter = hiddenList.begin();
  if (hiddenList.size() > 0)
    iter++;
  hiddenList.insert(iter, bridge);
  
  return CUBIT_SUCCESS;
}
 
CubitStatus HiddenEntitySet::restore( TopologyBridge* bridge )
{
  if( bridge->owner() != this )
    return CUBIT_FAILURE;
  
  std::vector<TopologyBridge*>::iterator iter;
  iter = std::find(hiddenList.begin(), hiddenList.end(), bridge);
  if (iter != hiddenList.end())
    hiddenList.erase(iter);
  else
    return CUBIT_FAILURE;

  bridge->owner( 0 );
    
  return CUBIT_SUCCESS;
}

CubitStatus HiddenEntitySet::remove_bridge( TopologyBridge* bridge )
{
  return restore(bridge);
}

CubitStatus HiddenEntitySet::swap_bridge( TopologyBridge* old_tb,
                                          TopologyBridge* new_tb,
                                          bool )
{
  if ( new_tb->owner() || !restore(old_tb) )
  {
    assert(0);
    return CUBIT_FAILURE;
  }
  
  return hide(new_tb);
}

CubitStatus HiddenEntitySet::merge( HiddenEntitySet* other )
{
  std::vector<TopologyBridge*>::iterator iter;
  std::vector<TopologyBridge*> other_list;

  // copy the list being merged so we can operate on the real list
  // without invalidating the iterator
  other_list = other->hiddenList;
  for (iter=other_list.begin(); iter!=other_list.end(); iter++ )
  {
    other->restore( *iter );
    hide( *iter );
  }

  return CUBIT_SUCCESS;
}


void HiddenEntitySet::hidden_surfaces( DLIList<Surface*>& result )
{
  std::vector<TopologyBridge*>::iterator iter;
  for (iter=hiddenList.begin(); iter!=hiddenList.end(); iter++)
  {
    if (Surface* ptr = dynamic_cast<Surface*>( *iter ))
      result.append(ptr);
  }
}

void HiddenEntitySet::hidden_coedges( DLIList<CoEdgeSM*>& result )
{
  std::vector<TopologyBridge*>::iterator iter;
  for (iter=hiddenList.begin(); iter!=hiddenList.end(); iter++)
  {
    if (CoEdgeSM* ptr = dynamic_cast<CoEdgeSM*>( *iter ))
      result.append(ptr);
  }
}

void HiddenEntitySet::hidden_curves( DLIList<Curve*>& result )
{
  std::vector<TopologyBridge*>::iterator iter;
  for (iter=hiddenList.begin(); iter!=hiddenList.end(); iter++)
  {
    if (Curve* ptr = dynamic_cast<Curve*>( *iter ))
      result.append(ptr);
  }
}

void HiddenEntitySet::hidden_points( DLIList<Point*>& result )
{
  std::vector<TopologyBridge*>::iterator iter;
  for (iter=hiddenList.begin(); iter!=hiddenList.end(); iter++)
  {
    if (Point* ptr = dynamic_cast<Point*>( *iter ))
      result.append(ptr);
  }
}

void HiddenEntitySet::print_debug_info( const char* prefix ) const
{
  if (!prefix)
    prefix = "";
  
  PRINT_INFO("%sHiddenEntitySet %p owned by %s %p\n", 
    prefix, this, 
    myOwner ? fix_type_name(typeid(*myOwner).name()) : "(null)",
    myOwner );
    
  char* new_prefix = new char[strlen(prefix)+3];
  strcpy( new_prefix, prefix );
  strcat( new_prefix, "  ");

  std::vector<TopologyBridge*>::const_iterator iter;
  for (iter=hiddenList.begin(); iter!=hiddenList.end(); iter++)
  {
    if( CompositePoint* cp = dynamic_cast<CompositePoint*>(*iter) )
      cp->print_debug_info(new_prefix, true);
    else if( CompositeCurve* cc = dynamic_cast<CompositeCurve*>(*iter) )
      cc->print_debug_info(new_prefix, true);
    else if( CompositeCoEdge* ce = dynamic_cast<CompositeCoEdge*>(*iter) )
      ce->print_debug_info(new_prefix, true);
    else if( CompositeSurface* cs = dynamic_cast<CompositeSurface*>(*iter) )
      cs->print_debug_info(new_prefix, true);
    else
#ifdef TOPOLOGY_BRIDGE_IDS
      PRINT_INFO("%s%s %d\n", new_prefix, fix_type_name(typeid(*(*iter)).name()), (*iter)->get_id() );
#else
      PRINT_INFO("%s%s %p\n", new_prefix, fix_type_name(typeid(*(*iter)).name()), (*iter) );
#endif
  }
  delete [] new_prefix;
}

void HiddenEntitySet::notify_reversed( TopologyBridge* )
  {} 
