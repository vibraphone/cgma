//-------------------------------------------------------------------------
// Filename      : TopologyBridge.cpp
//
// Purpose       : 
//                 
// Creator       : Tim Tautges
//
// Creation Date : 03/22/99
//
// Owner         : Darryl Melander
//-------------------------------------------------------------------------
#include "TopologyBridge.hpp"
#include "TopologyEntity.hpp"
#include "DLIList.hpp"
#include "BodySM.hpp"
#include "Lump.hpp"
#include "ShellSM.hpp"
#include "Surface.hpp"
#include "LoopSM.hpp"
#include "Curve.hpp"
#include "CoEdgeSM.hpp"
#include "Point.hpp"
#include "CubitSimpleAttrib.hpp"
#include "BridgeManager.hpp"
#include "TBOwner.hpp"
#include "TBOwnerSet.hpp"
#include "CADefines.hpp"

TopologyBridge::~TopologyBridge()
{
  if (bridgeOwner)
    bridgeOwner->bridge_destroyed(this);
}

TopologyEntity* TopologyBridge::topology_entity() const
{
  BridgeManager* bridgeManager = bridge_manager();
  if (bridgeManager)
    return bridgeManager->topology_entity();
  else
    return NULL;
}

BodySM *TopologyBridge::bodysm()
{
  DLIList<BodySM*> bodies;
  bodysms(bodies);
  if (bodies.size() > 0) return bodies.get();
  else return NULL;
}

Lump *TopologyBridge::lump()
{
  DLIList<Lump*> lump_list;
  lumps(lump_list);
  if (lump_list.size() > 0) return lump_list.get();
  else return NULL;
}


LoopSM *TopologyBridge::loopsm()
{
  DLIList<LoopSM*> loopsm_list;
  loopsms(loopsm_list);
  if (loopsm_list.size() > 0) return loopsm_list.get();
  else return NULL;
}




void TopologyBridge::bodysms( DLIList<BodySM*>& bodies, bool unique )
{
  DLIList<TopologyBridge*> parents;
  BodySM* this_body = dynamic_cast<BodySM*>(this);
  if( this_body )
  {
    bodies.append( this_body );
  }
  else if( dynamic_cast<Lump*>(this) )
  {
    get_parents( parents );
    parents.reset();
    for ( int i = parents.size(); i--; )
      bodies.append (dynamic_cast<BodySM*>(parents.get_and_step()));
  }
  else
  {
    get_parents( parents );
    for( int i = parents.size(); i--; )
      parents.get_and_step()->bodysms( bodies, false );
    
    if (unique)
      bodies.uniquify_ordered();
  }
}
    
void TopologyBridge::lumps( DLIList<Lump*>& lumps, bool unique )
{
  DLIList<TopologyBridge*> related;
  Lump* this_lump = dynamic_cast<Lump*>(this);
  if( this_lump )
  {
    lumps.append( this_lump );
  }
  else if( dynamic_cast<BodySM*>(this) )
  {
    get_children( related );
    related.reset();
    for ( int i = related.size(); i--; )
      lumps.append( dynamic_cast<Lump*>(related.get_and_step()) );
  }
  else if( dynamic_cast<ShellSM*>(this) )
  {
    get_parents( related );
    related.reset();
    for ( int i = related.size(); i--; )
      lumps.append( dynamic_cast<Lump*>(related.get_and_step()) );
  }
  else
  {
    get_parents( related );
    for( int i = related.size(); i--; )
      related.get_and_step()->lumps( lumps, false );
      
    if (unique)
      lumps.uniquify_ordered();
  }
}

void TopologyBridge::shellsms( DLIList<ShellSM*>& shells, bool unique )
{
  DLIList<TopologyBridge*> related;
  ShellSM* this_shell = dynamic_cast<ShellSM*>(this);
  if( this_shell )
  {
    shells.append(this_shell);
  }
  else if( dynamic_cast<Lump*>(this) )
  {
    get_children( related );
    related.reset();
    for ( int i = related.size(); i--; )
      shells.append( dynamic_cast<ShellSM*>(related.get_and_step()) );
  }
  else if( dynamic_cast<Surface*>(this) )
  {
    get_parents( related );
    related.reset();
    for ( int i = related.size(); i--; )
      shells.append( dynamic_cast<ShellSM*>(related.get_and_step()) );
  }
  else if( dynamic_cast<BodySM*>(this) )
  {
    get_children( related );
    for ( int i = related.size(); i--; )
      related.get_and_step()->shellsms( shells );
  }
  else 
  {
    get_parents( related );
      
    for( int i = related.size(); i--; )
      related.get_and_step()->shellsms( shells, false );
    
    if ( unique )
      shells.uniquify_ordered();
  }
}  

void TopologyBridge::surfaces( DLIList<Surface*>& surfaces, bool unique )
{
  DLIList<TopologyBridge*> related;
  Surface* this_surf = dynamic_cast<Surface*>(this);
  if( this_surf )
  {
    surfaces.append(this_surf);
  }
  else if( dynamic_cast<LoopSM*>(this) )
  {
    get_parents( related );
    related.reset();
    for (int i = related.size(); i--; )
      surfaces.append( dynamic_cast<Surface*>(related.get_and_step()) );
  }
  else if( dynamic_cast<ShellSM*>(this) )
  {
    get_children( related );
    related.reset();
    for (int i = related.size(); i--; )
      surfaces.append( dynamic_cast<Surface*>(related.get_and_step()) );
  }
  else if ( dynamic_cast<BodySM*>(this) || dynamic_cast<Lump*>(this) )
  {
    get_children(related);
    related.reset();
    for (int i = related.size(); i--; )
      related.get_and_step()->surfaces( surfaces, false );
      
    if (unique)
      surfaces.uniquify_ordered();
  }
  else
  {
    get_parents(related);
    related.reset();
    for (int i = related.size(); i--; )
      related.get_and_step()->surfaces( surfaces, false );

    if (unique && !dynamic_cast<CoEdgeSM*>(this))
      related.uniquify_ordered();
  }
}  

void TopologyBridge::loopsms( DLIList<LoopSM*>& loops, bool unique )
{
  DLIList<TopologyBridge*> related;
  
  LoopSM* this_loop = dynamic_cast<LoopSM*>(this);
  if( this_loop )
  {
    loops.append(this_loop);
  }
  else if( dynamic_cast<CoEdgeSM*>(this ) )
  {
    get_parents( related );
    related.reset();
    for (int i = related.size(); i--; )
      loops.append( dynamic_cast<LoopSM*>(related.get_and_step()) );
  }
  else if( dynamic_cast<Surface*>(this) )
  {
    get_children( related );
    related.reset();
    for (int i = related.size(); i--; )
      loops.append( dynamic_cast<LoopSM*>(related.get_and_step()) );
  }
  else
  {
    if( dynamic_cast<Curve*>(this) || dynamic_cast<TBPoint*>(this) )
      get_parents( related );
    else
      get_children( related );
      
    related.reset();
    for( int i = related.size(); i--; )
      related.get_and_step()->loopsms( loops, false );
    
    if (unique)
      loops.uniquify_ordered();
  }
}  


void TopologyBridge::coedgesms( DLIList<CoEdgeSM*>& coedges, bool unique )
{
  DLIList<TopologyBridge*> related;
  
  CoEdgeSM* this_coedge = dynamic_cast<CoEdgeSM*>(this);
  if( this_coedge )
  {
    coedges.append(this_coedge);
  }
  else if( dynamic_cast<Curve*>(this) )
  {
    get_parents( related );
    related.reset();
    for (int i = related.size(); i--; )
      coedges.append( dynamic_cast<CoEdgeSM*>(related.get_and_step()) );
  }
  else if( dynamic_cast<LoopSM*>(this) )
  {
    get_children( related );
    related.reset();
    for (int i = related.size(); i--; )
      coedges.append( dynamic_cast<CoEdgeSM*>(related.get_and_step()) );
  }
  else
  {
    if( dynamic_cast<TBPoint*>(this) )
      get_parents( related );
    else
      get_children( related );
    
    related.reset();  
    for( int i = related.size(); i--; )
      related.get_and_step()->coedgesms(coedges, false);
    
    if (unique)
      coedges.uniquify_ordered();
  }
}  

void TopologyBridge::curves_ignore_virtual( DLIList<Curve*>& curves, bool unique )
{
  int i;
  if(dynamic_cast<BodySM*>(this))
  {
    DLIList<TopologyBridge*> lumps;
    this->get_children_virt(lumps);
    for(i=lumps.size(); i>0; i--)
      lumps.get_and_step()->curves_ignore_virtual(curves, unique);
  }
  else if(dynamic_cast<Lump*>(this))
  {
    DLIList<TopologyBridge*> shells;
    this->get_children_virt(shells);
    for(i=shells.size(); i>0; i--)
      shells.get_and_step()->curves_ignore_virtual(curves, unique);
  }
  else if(dynamic_cast<ShellSM*>(this))
  {
    DLIList<TopologyBridge*> surfs;
    this->get_children_virt(surfs);
    for(i=surfs.size(); i>0; i--)
      surfs.get_and_step()->curves_ignore_virtual(curves, unique);
  }
  else if(dynamic_cast<Surface*>(this))
  {
    DLIList<TopologyBridge*> loops;
    this->get_children_virt(loops);
    for(i=loops.size(); i>0; i--)
      loops.get_and_step()->curves_ignore_virtual(curves, unique);
  }
  else if(dynamic_cast<LoopSM*>(this))
  {
    DLIList<TopologyBridge*> coedges;
    this->get_children_virt(coedges);
    for(i=coedges.size(); i>0; i--)
      coedges.get_and_step()->curves_ignore_virtual(curves, unique);
  }
  else if(dynamic_cast<CoEdgeSM*>(this))
  {
    DLIList<TopologyBridge*> tmp_curves;
    this->get_children_virt(tmp_curves);
    for(i=tmp_curves.size(); i>0; i--)
      tmp_curves.get_and_step()->curves_ignore_virtual(curves, unique);
  }
  else if(dynamic_cast<Curve*>(this))
  {
    curves.append(dynamic_cast<Curve*>(this));
  }
  if(unique)
    curves.uniquify_ordered();
}

void TopologyBridge::surfaces_ignore_virtual( DLIList<Surface*>& surfaces, bool unique )
{
  int i;
  if(dynamic_cast<BodySM*>(this))
  {
    DLIList<TopologyBridge*> lumps;
    this->get_children_virt(lumps);
    for(i=lumps.size(); i>0; i--)
      lumps.get_and_step()->surfaces_ignore_virtual(surfaces, unique);
  }
  else if(dynamic_cast<Lump*>(this))
  {
    DLIList<TopologyBridge*> shells;
    this->get_children_virt(shells);
    for(i=shells.size(); i>0; i--)
      shells.get_and_step()->surfaces_ignore_virtual(surfaces, unique);
  }
  else if(dynamic_cast<ShellSM*>(this))
  {
    DLIList<TopologyBridge*> surfs;
    this->get_children_virt(surfs);
    for(i=surfs.size(); i>0; i--)
      surfs.get_and_step()->surfaces_ignore_virtual(surfaces, unique);
  }
  else if(dynamic_cast<Surface*>(this))
  {
    surfaces.append(dynamic_cast<Surface*>(this));
  }
  if(unique)
    surfaces.uniquify_ordered();
}

void TopologyBridge::curves( DLIList<Curve*>& curves, bool unique )
{
  DLIList<TopologyBridge*> related;
  
  Curve* this_curve = dynamic_cast<Curve*>(this);
  if( this_curve )
  {
    curves.append(this_curve);
  }
  else if( dynamic_cast<TBPoint*>(this) )
  {
    get_parents( related );
    related.reset();
    for (int i = related.size(); i--; )
      curves.append( dynamic_cast<Curve*>(related.get_and_step()) );
  }
  else if( dynamic_cast<CoEdgeSM*>(this) )
  {
    get_children( related );
    related.reset();
    for (int i = related.size(); i--; )
      curves.append( dynamic_cast<Curve*>(related.get_and_step()) );
  }
  else
  {
    get_children( related );
    
    related.reset();
    for (int i = related.size(); i--; )
      related.get_and_step()->curves (curves, false);
    
    if (unique)
      curves.uniquify_ordered();
  }
}  

void TopologyBridge::points( DLIList<TBPoint*>& points, bool unique )
{
  DLIList<TopologyBridge*> children;
  
  TBPoint* this_point = dynamic_cast<TBPoint*>(this);
  if( this_point )
  {
    points.append(this_point);
  }
  else if( dynamic_cast<Curve*>(this) )
  {
    get_children( children );
    children.reset();
    for (int i = children.size(); i--; )
      points.append( dynamic_cast<TBPoint*>(children.get_and_step()) );
  }
  else
  {
    get_children( children );
      
    for( int i = children.size(); i--; )
      children.get_and_step()->points(points, false);
    
    if (unique)
      points.uniquify_ordered();
  }
}  

BridgeManager* TopologyBridge::bridge_manager() const
{ return dynamic_cast<BridgeManager*>(bridgeOwner); }

void TopologyBridge::bridge_manager( BridgeManager* manager )
{ bridgeOwner = manager; }

void TopologyBridge::get_parents( DLIList<TopologyBridge*>& parents )
{
  parents.clean_out();
  get_parents_virt( parents );
  for ( int i = parents.size(); i--; )
  {
    TopologyBridge* tb_ptr = parents.step_and_get();
    TBOwnerSet* partition_body = dynamic_cast<TBOwnerSet*>(tb_ptr->owner());
    if ( partition_body )
    {
      DLIList<TopologyBridge*> owner_list;
      partition_body->get_owners(owner_list);
      assert(owner_list.size() == 1);
      parents.change_to(owner_list.get());
    }
  }
}

/*
void TopologyBridge::get_parents( DLIList<TopologyBridge*>& parents,
                                  int layer, bool return_hidden )
{
  assert( this->layer() <= layer );
  DLIList<TopologyBridge*> parents_real;
  get_parents_virt( parents_real );
  while( parents_real.size() )
  {
      // order of parent lists not important (i.e. okay if 
      // result list is reverse of parents_real.)
    TopologyBridge* parent_real = parents_real.pop();
    TBOwner* owner = parent_real->owner();

    TBOwnerSet* set = 0;
    TopologyBridge* composite = 0;

    if( !owner || dynamic_cast<BridgeManager*>(owner) )
    {
      parents.append_unique( parent_real );
    }
    else if( set = dynamic_cast<TBOwnerSet*>(owner) )
    {
      if( set->get_owner_layer() > layer )
      {
        parents.append_unique( parent_real );
      }
      else
      {
        DLIList<TopologyBridge*> bridge_set, children;
        set->get_owners( bridge_set );
        for( int i = bridge_set.size(); i--; )
        {
          TopologyBridge* partition = bridge_set.get_and_step();
          children.clean_out();
          partition->get_children_virt( children );
          if( children.is_in_list( this ) )
            parents_real.append_unique( partition );
        }
      }
    }
    else if( composite = dynamic_cast<TopologyBridge*>(owner) )
    {
      if( composite->layer() > layer )
        parents.append_unique(parent_real);
      else
        parents_real.append_unique( composite );
    }
    else if( return_hidden )
    {
      parents.append( parent_real );
    }
    else
    {
      ; // Do nothing.  This entity is hidden by some composite.
    }
  }
}
*/
            
void TopologyBridge::get_children( DLIList<TopologyBridge*>& children,
                                   bool return_hidden, int layer )
{
  assert(this->layer() <= layer);
  
  DLIList<TopologyBridge*> child_list[2], partitions;
  int current = 0;
  
  get_children_virt( child_list[current] );
  // Note: be careful to return result list in the same order
  // as the child list we got above.
  
  bool done = false;
  while( !done )
  {
    done = true;
    
    int next = 1 - current;
    child_list[current].reset();
    child_list[next].clean_out();
    for( int i = child_list[current].size(); i--; )
    {
      TopologyBridge* child = child_list[current].get_and_step();
      TBOwner* owner = child->owner();
      
      TopologyBridge* composite = 0;
      TBOwnerSet* partition = 0;
      
      if( !owner || dynamic_cast<BridgeManager*>(owner) )
      {
        child_list[next].append( child );
      }
      else if( (partition = dynamic_cast<TBOwnerSet*>(owner) ) != NULL )
      {
        if( partition->get_owner_layer() > layer )
        {
          child_list[next].append(child);
        }
        else
        {
          done = false;
          partitions.clean_out();
          partition->get_owners( partitions );
          partitions.reset();
          child_list[next] += partitions;
        }
      }
      else if( (composite = dynamic_cast<TopologyBridge*>(owner) ) != NULL )
      {
        if( composite->layer() > layer )
        {
           child_list[next].append(child);
        }
        else
        {
          done = false;
          if( !child_list[current].is_in_list(composite) )
            child_list[next].append_unique( composite );
        }
      }
      else if( return_hidden )
      {
        child_list[next].append( child );
      }
      else
      {
        ; // Do nothing. This entity is hidden by some composite.
      }
    }
      
    current = next;
  }
  
  children = child_list[current];
}
