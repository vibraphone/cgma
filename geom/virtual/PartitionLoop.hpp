//-------------------------------------------------------------------------
// Filename      : PartitionLoop.hpp
//
// Purpose       : LoopSM used by Partition geometry
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 01/11/02
//-------------------------------------------------------------------------

#ifndef PARTITION_LOOP_HPP
#define PARTITION_LOOP_HPP

#include "LoopSM.hpp"
#include "PartitionCoEdge.hpp"

class PartitionLoop : public LoopSM
{
friend class PartitionSurface;
public:

  PartitionLoop( );
  virtual ~PartitionLoop();
  
  PartitionSurface* get_surface() const;
  //void surface( PartitionSurface* );
  
  virtual CubitBoolean is_external() ;
    //R CubitBoolean
    //R- CUBIT_TRUE/CUBIT_FALSE
    //- Returns CUBIT_TRUE if the Loop is an external Loop and CUBIT_FALSE
    //- otherwise.
  
  virtual LoopType loop_type() ;

  PartitionCoEdge* first_coedge( );
  PartitionCoEdge* next_coedge( PartitionCoEdge* after_this );
  PartitionCoEdge* prev_coedge( PartitionCoEdge* before_this );
  int num_coedges();
  
  CubitStatus insert_after( PartitionCoEdge* coedge, PartitionCoEdge* after );
  CubitStatus insert_before( PartitionCoEdge* coedge, PartitionCoEdge* before );
  CubitStatus remove( PartitionCoEdge* coedge );
  CubitStatus remove_all_coedges( DLIList<PartitionCoEdge*>* removed = 0);
  
  void get_parents_virt( DLIList<TopologyBridge*>& parents );
  void get_children_virt( DLIList<TopologyBridge*>& children );
  int layer() const ;
  GeometryQueryEngine* get_geometry_query_engine() const;
  
  void append_simple_attribute_virt( const CubitSimpleAttrib& );
  void remove_simple_attribute_virt( const CubitSimpleAttrib& );
  void remove_all_simple_attribute_virt();
  CubitStatus get_simple_attribute( DLIList<CubitSimpleAttrib>& );
  CubitStatus get_simple_attribute( const CubitString&,
                                    DLIList<CubitSimpleAttrib>& );
 
  void reverse();

  void print_debug_info( const char* line_prefix = 0 );

private:

  PartitionLoop( const PartitionLoop& );

  PartitionSurface* mySurface;
  PartitionCoEdge* firstCoedge;
  PartitionLoop* nextInSurface;
  int numCoedges;
};

inline PartitionSurface* PartitionLoop::get_surface() const
  { return mySurface; }

inline PartitionCoEdge* PartitionLoop::first_coedge( )
  { return firstCoedge; }

inline int PartitionLoop::num_coedges()
  { return numCoedges; }

inline PartitionCoEdge* PartitionLoop::next_coedge( PartitionCoEdge* prev )
  { return prev->myLoop == this ? prev->loopNext : 0; }

inline PartitionCoEdge* PartitionLoop::prev_coedge( PartitionCoEdge* next )
  { return next->myLoop == this ? next->loopPrev : 0; }

#endif
