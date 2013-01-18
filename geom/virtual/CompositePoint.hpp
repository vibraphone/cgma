//-------------------------------------------------------------------------
// Filename      : CompositePoint.hpp
//
// Purpose       : Decorator for Points owned by higher-order Composites
//
// Special Notes : This object is used a) to complete the child topology
//                 graph of curves and b) to stitch points.
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 03/07/02
//-------------------------------------------------------------------------

#ifndef COMPOSITE_POINT_HPP
#define COMPOSITE_POINT_HPP

#include "Point.hpp"
#include "TBOwner.hpp"
#include "CompositeCurve.hpp"
#include "VGDefines.h"
template <class X> class DLIList;

class CompositePoint : public TBPoint, public TBOwner
{
friend class CompositeCurve;

public:
  int HadBridgeRemoved;

  CompositePoint( TBPoint* real_pt );
  virtual ~CompositePoint();
  
  CompositeCurve* next_curve( CompositeCurve* prev = 0 ) const
    { return prev ? prev->next( this ) : firstCurve; }
  
  TBPoint* get_point() const
    { return realPoint; }
  
  void append_simple_attribute_virt( const CubitSimpleAttrib& csa )
    { get_point()->append_simple_attribute_virt( csa ); }
  void remove_simple_attribute_virt( const CubitSimpleAttrib& csa )
    { get_point()->remove_simple_attribute_virt( csa ); }
  void remove_all_simple_attribute_virt()
    { get_point()->remove_all_simple_attribute_virt(); }
  CubitStatus get_simple_attribute( DLIList<CubitSimpleAttrib>& list )
    { return get_point()->get_simple_attribute( list ); }
  CubitStatus get_simple_attribute( const CubitString& name,
                                    DLIList<CubitSimpleAttrib>& attrib_list )
    { return get_point()->get_simple_attribute( name, attrib_list ); }
  
  GeometryQueryEngine* get_geometry_query_engine() const;

  void get_parents_virt( DLIList<TopologyBridge*>& );
  void get_children_virt( DLIList<TopologyBridge*>& );
  int layer() const { return COMPOSITE_LAYER; }
  
  CubitVector coordinates() const
    { return get_point()->coordinates(); }
  
  CubitStatus remove_bridge( TopologyBridge* bridge );
  CubitStatus swap_bridge( TopologyBridge* old_tb, TopologyBridge* new_tb, bool );
  void notify_reversed( TopologyBridge* bridge );

  void print_debug_info( const char* prefix = 0, bool brief = false ) const;

  CubitStatus stitch( CompositePoint* point );
  void unstitch_all();
  void get_stitched( DLIList<CompositePoint*>& result );

private:

  CompositeCurve* firstCurve;
  
  TBPoint* realPoint;
  
  CompositePoint* stitchNext;
};

#endif
