//-------------------------------------------------------------------------
// Filename      : CompositeCurve.hpp
//
// Purpose       : Geometry defined as the joining of a chain of curves.
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 12/19/01
//-------------------------------------------------------------------------

#ifndef COMPOSITE_CURVE_HPP
#define COMPOSITE_CURVE_HPP

#include "VGDefines.h"
#include "Curve.hpp"
#include "CompositeGeom.hpp"
#include "HiddenEntitySet.hpp"
#include "CompositeCoEdge.hpp"

class CompositeSurface;
class CompositePoint;
class Point;

class CompositeCurve : public Curve, public TBOwner
{

  public:
  
    int HadBridgeRemoved;
    //CompositeCurve( DLIList<Curve*>& curve_list, bool periodic );
    //CompositeCurve( );
    CompositeCurve( Curve* curve );
    CompositeCurve( CompositeGeom* geometry );
    CompositeCurve( CompositePoint* point ); // create point-curve
    
    virtual ~CompositeCurve();
   
    inline Curve* get_curve( int index ) const;
    inline CubitSense get_sense( int index ) const;
    inline int num_curves( ) const;
    inline int index_of( Curve* curve_ptr ) const;
    inline void update();
    Curve* remove_curve( int index );
    
/*
    inline CubitStatus append_curve( Curve* curve_ptr );
    inline CubitStatus prepend_curve( Curve* curve_ptr );
    CubitStatus insert_curve( Curve* curve_ptr, int index );
    Curve* remove_curve( Curve* curve_ptr );
    Curve* remove_curve( int index );
*/
    CompositeCurve* split( Curve* curve_ptr );
    CubitStatus combine( CompositeCurve* curve_ptr, bool prepend );
    void reverse();
    
    CompositeCoEdge* first_coedge() const;
    CompositeCoEdge* next_coedge( CompositeCoEdge* after_this ) const;
#ifdef BOYD15
    CompositeCoEdge* find_coedge( CompositeSurface* surface ) const;
#endif
    CubitStatus add( CompositeCoEdge* coedge );
    CubitStatus remove( CompositeCoEdge* coedge );
    
    CompositePoint* start_point() const;
    CompositePoint* end_point() const;
    CubitStatus start_point( CompositePoint* pt );
    CubitStatus end_point( CompositePoint* pt );
    CompositePoint* other_point( CompositePoint* pt );
#ifdef BOYD15
    CompositePoint* closest_end_point( const CubitVector& pos );
#endif
    CompositePoint* common_point( CompositeCurve* curve );
    CompositeCurve* next( const CompositePoint* around_this );

    HiddenEntitySet& hidden_entities();
    void get_hidden_points( DLIList<Point*>& points );
    
    bool has_parent_composite_surface() const;
    
    /**************** Methods from TopologyBridge *******************/
    
    void append_simple_attribute_virt( CubitSimpleAttrib* csa );
    void remove_simple_attribute_virt(CubitSimpleAttrib* csa );
    void remove_all_simple_attribute_virt();
    CubitStatus get_simple_attribute(DLIList<CubitSimpleAttrib*>& csa_list );
    CubitStatus get_simple_attribute( const CubitString& name,
                                  DLIList<CubitSimpleAttrib*>& attrib_list );
    GeometryQueryEngine* get_geometry_query_engine() const;
    
    void get_parents_virt( DLIList<TopologyBridge*>& parents );
    void get_children_virt( DLIList<TopologyBridge*>& children );
    int layer() const { return COMPOSITE_LAYER; }

    /******************** Methods from TBOwner **********************/
   
    CubitStatus remove_bridge( TopologyBridge* bridge );
    CubitStatus swap_bridge( TopologyBridge* old, 
                             TopologyBridge* neww, 
                             bool reversed );
    CubitBoolean contains_bridge( TopologyBridge* bridge ) const;
    void notify_reversed( TopologyBridge* bridge );
    void notify_split( TopologyBridge* new_bridge, TopologyBridge* split_from );
//    void notify_joined( TopologyBridge* dead, TopologyBridge* keep );

    /**************** Methods from GeometryEntity *******************/

	  CubitBox bounding_box() const;
    double measure( ) ;
    GeometryType geometry_type();
    

    /******************** Methods from Curve ***********************/

    CubitBoolean get_param_range( double& lower, double& upper );
    CubitBoolean is_periodic( double& period);
    double start_param();
    double end_param();

    CubitStatus position_from_u( double u_value, CubitVector& pos );
    double u_from_position (const CubitVector& input_position);

    double length_from_u( double parameter1, double parameter2 );
    double u_from_arc_length ( double root_param, double arc_length );

    CubitStatus closest_point( CubitVector const& location, 
                               CubitVector& closest_location,
                               CubitVector* tangent_ptr = NULL,
                               CubitVector* curvature_ptr = NULL,
                               double *param = NULL);
#ifdef BOYD15
    void get_tangent( CubitVector const& location, CubitVector& tangent );
    void get_curvature( CubitVector const& location, CubitVector& curvature );
#endif

    CubitStatus closest_point_trimmed( CubitVector const& from_pt,
	                                   CubitVector& result_pt );

    CubitBoolean is_position_on( const CubitVector &test_position );
    CubitPointContainment point_containment( const CubitVector &point );

    CubitStatus get_point_direction( CubitVector& origin,
	                                   CubitVector& direction );

    CubitStatus get_interior_extrema( 
                            DLIList<CubitVector*>& interior_points,
                            CubitSense& return_sense);

    CubitStatus get_center_radius( CubitVector& c, double& r );

    CubitBoolean G1_discontinuous( double param,
                CubitVector* minus_tangent = NULL,
                CubitVector* plus_tangent = NULL );


    void print_debug_info( const char* prefix = 0, bool brief = false ) const;

    CubitStatus stitch( CompositeCurve* merge_with );
    void unstitch_all();
    CompositeCurve* primary_stitched_curve();
    bool is_stitched();
    void get_stitched( DLIList<CompositeCurve*>& list );

    CubitStatus curve_param( double u_composite, double& u_curve, int& index ) const;
 	  double composite_param( int index, double param ) const;
    
    //void draw( int color );
  
    void read_attributes() ; //{ compGeom->read_attributes(); }
    void write_attributes() ; //{ compGeom->write_attributes(); }

  protected: 
  
    CompositeCurve();

    int closest_curve( CubitVector const& location,
                       CubitVector *point = NULL );
    //R int
    //R- The index of the closest curve
    //I location
    //I- A position for which the closest curve is desired.
    //O point
    //O- The closest point on the curve
    //- This function finds the closest underlying Curve of this 
    //- CompositeCurve to a given point
  
 	  double lengthUntilI( int index ) const;
  
    void fixup_periodic_param( double& u ) const;
    
    CubitStatus set_point( bool set_start_point, CompositePoint* point );

private:

      // these have no implementation, just private delcarations
      // to prevent the compiler from generating default implementations
    CompositeCurve& operator=(const CompositeCurve&);
    CompositeCurve(const CompositeCurve&);

    CompositeGeom* compGeom;
    
    HiddenEntitySet* hiddenSet;
    
    CompositeCoEdge* firstCoEdge;
    
    CompositePoint* startPoint;
    CompositePoint* endPoint;
    
    CompositeCurve* startNext;
    CompositeCurve* endNext;
    
    CompositeCurve* stitchNext;
};

inline Curve* CompositeCurve::get_curve( int index ) const
  { return dynamic_cast<Curve*>(compGeom->entity(index)); }

inline CubitSense CompositeCurve::get_sense( int index ) const
  { return compGeom->sense( index ); }

inline int CompositeCurve::num_curves( ) const
  { return compGeom->num_entities(); }

inline int CompositeCurve::index_of( Curve* curve ) const
  { return compGeom->index_of( curve ); }

inline void CompositeCurve::update()
  { compGeom->update_cached_data();}

/*
inline CubitStatus CompositeCurve::prepend_curve( Curve* curve )
  { return insert_curve( curve, 0 ); }

inline CubitStatus CompositeCurve::append_curve( Curve* curve )
  { return insert_curve( curve, num_curves() ); }
*/
inline CompositeCoEdge* CompositeCurve::first_coedge( ) const
  { return firstCoEdge; }

inline CompositeCoEdge* CompositeCurve::next_coedge( CompositeCoEdge* coedge ) const
  { return !coedge ? firstCoEdge : coedge->myCurve == this ? coedge->nextOnCurve : 0; }

inline CompositePoint* CompositeCurve::start_point() const
  { return startPoint; }

inline CompositePoint* CompositeCurve::end_point() const
  { return endPoint; }



inline HiddenEntitySet& CompositeCurve::hidden_entities()
  { if( !hiddenSet ) 
      hiddenSet = new HiddenEntitySet(this);
    return *hiddenSet; 
  }


#endif
