//-------------------------------------------------------------------------
// Filename      : PartitionCurve-new.hpp
//
// Purpose       : 
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 04/10/02
//-------------------------------------------------------------------------
#ifndef PARTITION_CURVE_HPP
#define PARTITION_CURVE_HPP

#include "Curve.hpp"
#include "PartitionPoint.hpp"
#include "PartitionCoEdge.hpp"

class PartitionPoint;
class PartitionSurface;
class GMem;
class CubitFacetEdgeData;

class PartitionCurve : public Curve, public PartitionEntity
{

public:

  virtual PartitionCurve* split( double param ) = 0;
  
  virtual CubitStatus combine( PartitionCurve* dead_curve ) = 0;

  virtual ~PartitionCurve();

  PartitionPoint* start_point() const
    { return startPoint; }
  
  PartitionPoint* end_point() const
    { return endPoint; }
  
  PartitionPoint* other_point( const PartitionPoint* pt ) const
    { return pt == startPoint ? endPoint : pt == endPoint ? startPoint : 0; }
    
  int num_coedges() const;
  
  PartitionCoEdge* next_coedge( const PartitionCoEdge* prev = 0 ) const
    { return !prev ? firstCoEdge : prev->myCurve == this ? prev->curveNext : 0; }
  
  CubitStatus add( PartitionCoEdge* coedge );
  CubitStatus remove( PartitionCoEdge* coedge );
  void remove_all_coedges();
  
  CubitStatus start_point( PartitionPoint* point );
  CubitStatus end_point  ( PartitionPoint* point );
  
  virtual CubitStatus get_graphics( GMem& result,
                                    double angle_tolerance=0,
                                    double distance_tolerance=0,
                                    double max_edge_length=0) = 0; 
  
  PartitionCurve* next_curve( const PartitionPoint* about_this ) const
    { return about_this == startPoint ? startNext :
             about_this ==   endPoint ?   endNext : 0; }
  
  bool is_nonmanifold( const PartitionSurface* in_this_surf ) const;
  bool is_in_surface( const PartitionSurface* surf, 
                      bool manifold_only = false ) const;
  
  CubitStatus move_to_geometry( CubitVector& position );
  CubitStatus move_to_geometry( CubitPoint* facetPoint );
  
  void append_simple_attribute_virt( const CubitSimpleAttrib& );
  void remove_simple_attribute_virt( const CubitSimpleAttrib& );
  void remove_all_simple_attribute_virt();
  CubitStatus get_simple_attribute( DLIList<CubitSimpleAttrib>& );
  CubitStatus get_simple_attribute( const CubitString& name,
                                    DLIList<CubitSimpleAttrib>& );

  void get_parents_virt( DLIList<TopologyBridge*>& );
  void get_children_virt( DLIList<TopologyBridge*>& );
  int layer() const { return sub_entity_set().get_owner_layer(); }
  GeometryQueryEngine* get_geometry_query_engine() const;

		
  virtual void print_debug_info( const char* prefix = 0,
                                 bool print_subentity_set = true ) const;
  

  void get_facet_data( DLIList<CubitFacetEdgeData*>& result_list ) const;
  void set_facet_data( const DLIList<CubitFacetEdgeData*>& new_list );
  bool has_facet_data() const { return facetEdges.size() > 0; }
  void remove_facet_data();
  void replace_facet(CubitFacetEdgeData* dead_facet,
                     DLIList<CubitFacetEdgeData*> &new_facets);
  virtual void notify_split( FacetEntity* old_entity, FacetEntity* new_entity );
  CubitStatus fix_facet_data( PartitionCurve* new_curve );
  void remove_dead_facet( CubitFacetEdgeData* edge );
    //- update facet data for curve partition operation
    //- must be called after end points have been updated.
    
  void do_facet_cleanup();
    //- remove small edges, etc.
 
  void draw_facets(int color);
  
  virtual void transform( const CubitTransformMatrix& );

protected:

  PartitionCurve( );
  
  CubitStatus get_save_topology( DLIList<int>& points );
  
  void reverse_point_order();

private:

  DLIList<CubitFacetEdgeData*> facetEdges;

  CubitStatus remove_start_point();
  CubitStatus remove_end_point();
  CubitStatus remove_from_point( PartitionPoint* point, PartitionCurve* next );

  PartitionCoEdge* firstCoEdge;
  
  PartitionPoint *startPoint, *endPoint;
  PartitionCurve *startNext, *endNext;
};

#endif

  
  
