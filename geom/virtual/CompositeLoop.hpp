//-------------------------------------------------------------------------
// Filename      : CompositeLoop.hpp
//
// Purpose       : LoopSM used by composite geometry
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 01/11/02
//-------------------------------------------------------------------------

#ifndef COMPOSITE_LOOP_HPP
#define COMPOSITE_LOOP_HPP

#include "LoopSM.hpp"
#include "CompositeCoEdge.hpp"

class CompositeSurface;

class CompositeLoop : public LoopSM
{
friend class CompositeSurface;
public:

  CompositeLoop();
  virtual ~CompositeLoop();
  
  CompositeSurface* get_surface() const;
  CompositeLoop* next_loop() const;
  
  CompositeCoEdge* first_coedge( ) const;
  CompositeCoEdge* next_coedge( CompositeCoEdge* after_this ) const;
  CompositeCoEdge* prev_coedge( CompositeCoEdge* before_this ) const;
  int num_coedges() const;
  
  CubitStatus insert_after( CompositeCoEdge* coedge, CompositeCoEdge* after );
  CubitStatus insert_before( CompositeCoEdge* coedge, CompositeCoEdge* before );
  CubitStatus remove( CompositeCoEdge* coedge );
  CubitStatus remove_all_coedges( DLIList<CompositeCoEdge*>* removed = 0);
  
  void get_parents_virt( DLIList<TopologyBridge*>& parents );
  void get_children_virt( DLIList<TopologyBridge*>& children );
  GeometryQueryEngine* get_geometry_query_engine() const;
  int layer() const { return COMPOSITE_LAYER; }
  
  virtual LoopType loop_type() ;
  virtual CubitBoolean is_external() ;
    //R CubitBoolean
    //R- CUBIT_TRUE/CUBIT_FALSE
    //- Returns CUBIT_TRUE if the Loop is an external Loop and CUBIT_FALSE
    //- otherwise.

  void append_simple_attribute_virt( const CubitSimpleAttrib& );
  void remove_simple_attribute_virt( const CubitSimpleAttrib& );
  void remove_all_simple_attribute_virt();
  CubitStatus get_simple_attribute( DLIList<CubitSimpleAttrib>& );
  CubitStatus get_simple_attribute( const CubitString& name,
                                    DLIList<CubitSimpleAttrib>& attrib_list );

  // reverse the direction of the loop
  // if b_reverse_coedges is true, the function will also reverse the coedges
  // but if b_reverse_coedges is false, the caller needs to make sure coedges
  // get reversed properly
  void reverse(bool b_reverse_coedges);

  void print_debug_info( const char* line_prefix = 0 );

private:

  CompositeSurface* mySurface;
  CompositeCoEdge* myCoedge;
  CompositeLoop* loopNext;
  int numCoedges;
};

inline CompositeSurface* CompositeLoop::get_surface() const
  { return mySurface; }

inline CompositeLoop* CompositeLoop::next_loop() const
  { return loopNext; }

inline CompositeCoEdge* CompositeLoop::first_coedge( ) const
  { return myCoedge; }

inline int CompositeLoop::num_coedges() const
  { return numCoedges; }

inline CompositeCoEdge* CompositeLoop::next_coedge( CompositeCoEdge* prev ) const
  { return prev->myLoop == this ? prev->nextCoedge : 0; }

inline CompositeCoEdge* CompositeLoop::prev_coedge( CompositeCoEdge* next ) const
  { return next->myLoop == this ? next->prevCoedge : 0; }
  

#endif
