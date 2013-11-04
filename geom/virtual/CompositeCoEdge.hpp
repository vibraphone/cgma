//-------------------------------------------------------------------------
// Filename      : CompositeCoEdge.hpp
//
// Purpose       : Combined set of CoEdgeSMs
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 01/11/02
//-------------------------------------------------------------------------

#ifndef COMPOSITE_COEDGE_HPP
#define COMPOSITE_COEDGE_HPP

#include "VGDefines.h"
#include "CoEdgeSM.hpp"
#include "TBOwner.hpp"
#include "VGArray.hpp"

class CompositeCurve;
class CompositeLoop;
class CompositePoint;

class CompositeCoEdge : public CoEdgeSM, public TBOwner
{
friend class CompositeLoop;
friend class CompositeCurve;

public:

  CompositeCoEdge( CoEdgeSM* coedge_ptr );
  CompositeCoEdge( CompositeCurve* point_curve );
  
  ~CompositeCoEdge();
  
  int num_coedges() const;
  
  int index_of( CoEdgeSM* coedge_ptr ) const;
  int index_of( Curve* curve_ptr ) const;
  
  CubitSense get_sense( int index ) const;
  CoEdgeSM* get_coedge( int index ) const;
  
  CompositeCoEdge* next() const;
  CompositeCoEdge* prev() const;
    // Return next or previous coedge in a CompositeLoop.
    // If there is not a parent CompositeSurface, then
    // there is not a parent CompositeLoop and these will
    // return NULL.
  
  CubitStatus combine( CompositeCoEdge* dead_coedge, bool prepend );
  CubitStatus insert_coedge( int index, CoEdgeSM* coedge );
  CubitStatus remove_coedge( int index );
  CompositeCoEdge* split( int index );
  
  CompositePoint* start_point();
  CompositePoint* end_point();
/*
  CoEdgeSM* remove_first();
  CoEdgeSM* remove_last();
*/  
  virtual CubitSense sense();
  
  void sense( CubitSense sense );
  
  CompositeLoop* get_loop() const;
//  void loop( CompositeLoop* );
  
  CompositeCurve* get_curve() const;
//  void curve( CompositeCurve* );
  
  LoopSM* get_parent_loop();
    //- Get the parent loop at the composite level of
    //- the topo bridge graph.  I.E. Get parent CompositeLoop
    //- if there is one.  Otherwise get the SM-level loop.
  
  void get_parents_virt( DLIList<TopologyBridge*>& parents );
  void get_children_virt( DLIList<TopologyBridge*>& children );
  GeometryQueryEngine* get_geometry_query_engine() const;
  int layer() const { return COMPOSITE_LAYER; }

  void append_simple_attribute_virt( const CubitSimpleAttrib& simple_attrib_ptr );
  void remove_simple_attribute_virt( const CubitSimpleAttrib& simple_attrib_ptr );
  void remove_all_simple_attribute_virt();
  CubitStatus get_simple_attribute( DLIList<CubitSimpleAttrib>& attrib_list );
  CubitStatus get_simple_attribute( const CubitString& name,
                                    DLIList<CubitSimpleAttrib>& attrib_list );
  
  CubitStatus remove_bridge( TopologyBridge* bridge );
  CubitStatus swap_bridge( TopologyBridge* old_tb, 
                           TopologyBridge* new_tb,
                           bool reversed );
  CubitBoolean contains_bridge( TopologyBridge* bridge ) const;
  void notify_reversed( TopologyBridge* bridge );
  
  void reverse();
  
  void print_debug_info( const char* line_prefix = 0, bool brief = false );
  
private:

  CompositeCoEdge();
  
  VGArray<CoEdgeSM*> coedgeSet;
  
  CubitSense mySense;
  
  CompositeLoop*  myLoop;
  CompositeCoEdge* nextCoedge;
  CompositeCoEdge* prevCoedge;

  CompositeCurve* myCurve;
  CompositeCoEdge* nextOnCurve;
};

inline int CompositeCoEdge::num_coedges() const
  { return coedgeSet.size(); }

inline CoEdgeSM* CompositeCoEdge::get_coedge( int index ) const
  { return coedgeSet[index]; }

inline int CompositeCoEdge::index_of( CoEdgeSM* coedge ) const
  { return coedgeSet.find( coedge ); }

inline void CompositeCoEdge::sense( CubitSense sense )
  { assert( sense != CUBIT_UNKNOWN ); mySense = sense; }

inline CubitSense CompositeCoEdge::sense()
  { return mySense; }

inline CompositeCoEdge* CompositeCoEdge::next() const
  { return nextCoedge; }

inline CompositeLoop* CompositeCoEdge::get_loop() const
  { return myLoop; }

inline CompositeCurve* CompositeCoEdge::get_curve() const
  { return myCurve; }
#endif
