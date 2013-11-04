//-------------------------------------------------------------------------
// Filename      : CompositeSurface.hpp
//
// Purpose       : A class for composited surface geometry
//
// Special Notes : 
// 
// Creator       : Jason Kraftcheck
//
// Creation Date : 03/10/98
//
// Owner         : Jason Kraftcheck
//-------------------------------------------------------------------------

#ifndef COMPOSITE_SURFACE_HPP
#define COMPOSITE_SURFACE_HPP

#include "VGDefines.h"
#include "CompositeGeom.hpp"
#include "Surface.hpp"
#include "TBOwner.hpp"
#include "HiddenEntitySet.hpp"
#include "CompositeLoop.hpp"
#include "CompositeCoSurf.hpp"

class CompositeCoSurf;
class CompositeEngine;
class CompositePoint;
class CompositeCurve;
class CompositeShell;
class CompositeLump;
class CompSurfFacets;
class GMem;

class CompositeSurface: public Surface, public TBOwner
{
  friend class CompositeEngine;

public:
  CompositeSurface( Surface* surf_ptr );
  CompositeSurface( CompositeGeom* geom_ptr );
  //- Constructor 
		
  virtual ~CompositeSurface();
  //- Destructor
	
  //CubitStatus add( Surface* surface, CubitSense relative_sense );
  //CubitStatus remove( Surface* surface );	
  void get_ignored_surfs(DLIList<Surface*> &surfs);
  int num_surfs() const;
  int index_of( Surface* surface ) const;
  void update();
  Surface* get_surface( int index ) const;
  CubitSense get_sense( int index ) const;
  Surface* remove_surface( int index );
  void ignore_surface(int surface_id);
  void ignore_surface(Surface *surf);
  void unignore_surface(int surface_id);
  CompositeLoop* first_loop() const;
  CompositeLoop* next_loop( CompositeLoop* after_this = 0 ) const;
  CubitStatus add( CompositeLoop* loop );
  CubitStatus remove( CompositeLoop* loop );
  
  CubitStatus add( CompositeCoSurf* cosurf );
  CubitStatus remove( CompositeCoSurf* cosurf );
  
  CompositeCoSurf* next_co_surface( CompositeCoSurf* prev = 0 ) const;
  
  CompositeCoSurf* find_first( CompositeShell* shell ) const;
    // Get the first CoSurface connecting this surface w/ the passed shell
  CompositeCoSurf* find_first( CompositeLump* lump ) const;
    // Get the first CoSurface connecting this surface w/ a
    // shell owned by the passed volume
  CompositeCoSurf* find_next( CompositeCoSurf* cosurf ) const;
    // Get the next CoSurface after the passed CoSurface and
    // with the same Shell as the passed CoSurface.
  
  void get_curves( DLIList<CompositeCurve*>& curves );
  
  HiddenEntitySet& hidden_entities();
  bool has_hidden_entities() const;

  CompositeSurface* split( VGArray<int>& indices_to_move );
  CubitStatus combine( CompositeSurface* dead_surface );

  CubitStatus stitch( CompositeSurface* stitch_partner );
  CompositeSurface* unstitch( );
  CompositeSurface* get_stitch_partner() const
    { return stitchPartner; }
  

  /************************************************ 
   **** Methods Inhereted from GeometryEntity ***** 
   ************************************************/
  CubitBox bounding_box() const;
  double measure();
  GeometryQueryEngine* get_geometry_query_engine() const;
  GeometryType geometry_type();


  /************************************************ 
   ******** Methods Inhereted from TBOwner ******** 
   ************************************************/
  CubitStatus remove_bridge( TopologyBridge* bridge );
  CubitStatus swap_bridge( TopologyBridge* old_tb, TopologyBridge* new_tb,
                           bool reversed );
  CubitBoolean contains_bridge( TopologyBridge* bridge ) const;
  void notify_reversed( TopologyBridge* bridge );
  void notify_split( TopologyBridge* new_bridge, TopologyBridge* old_bridge );
  void notify_topology_modified( TopologyBridge* bridge );

  /************************************************ 
   **** Methods Inhereted from TopologyBridge ***** 
   ************************************************/
  void append_simple_attribute_virt( const CubitSimpleAttrib& simple_attrib_ptr );
  void remove_simple_attribute_virt( const CubitSimpleAttrib& simple_attrib_ptr );
  void remove_all_simple_attribute_virt();
  CubitStatus get_simple_attribute( DLIList<CubitSimpleAttrib>& attrib_list );
  CubitStatus get_simple_attribute( const CubitString& name,
                                    DLIList<CubitSimpleAttrib>& attrib_list );
  

  void get_parents_virt( DLIList<TopologyBridge*>& parents );
  void get_children_virt( DLIList<TopologyBridge*>& children );
  int layer() const { return COMPOSITE_LAYER; }

  /********************************************************** 
   **** Methods Inhereted from Surface, and implemented ***** 
   **********************************************************/
  virtual CubitSense get_shell_sense( ShellSM* shell_ptr ) const;
   
  virtual void closest_point_trimmed( CubitVector from_point, 
                                      CubitVector& result );
  virtual CubitStatus get_point_normal( CubitVector& origin, 
                                        CubitVector& normal );

  virtual CubitStatus closest_point_uv_guess(  
      CubitVector const& location,
      double &u, double &v,
      CubitVector* closest_location = NULL,
      CubitVector* unit_normal = NULL );

  virtual CubitStatus evaluate( double u, double v,
                                CubitVector *position,                                   
                                CubitVector *normal,
                                CubitVector *curvature1,
                                CubitVector *curvature2 );

  virtual CubitStatus get_projected_distance_on_surface( CubitVector *pos1,
                                                         CubitVector *pos2, 
                                                         double &distance );

  virtual CubitStatus closest_point( CubitVector const& location,
                                     CubitVector* closest_location = NULL,
                                     CubitVector* unit_normal = NULL,
                                     CubitVector* curvature1 = NULL,
                                     CubitVector* curvature2 = NULL );

  virtual CubitStatus closest_point_along_vector(CubitVector& from_point, 
                                         CubitVector& along_vector,
                                         CubitVector& point_on_surface);

  virtual CubitStatus principal_curvatures( CubitVector const& location,
                                            double& curvature_1, 
                                            double& curvature_2,
                                            CubitVector* closest_location );
  virtual CubitBoolean is_parametric();
  virtual CubitBoolean is_position_on(CubitVector& position);
  virtual CubitPointContainment point_containment( const CubitVector &point );
  virtual CubitPointContainment point_containment( double u, double v );
//  virtual CubitPointContainment point_containment( CubitVector &point, 
//                                                   double u, double v );
  virtual CubitSense get_geometry_sense();
  virtual void reverse_sense();
		

  /**************************************************************
   **** Methods Inhereted from Surface, and not implemented ***** 
   **************************************************************/
  virtual CubitVector position_from_u_v( double u, double v );
  virtual CubitStatus u_v_from_position( CubitVector const& position, 
                                         double& u, double& v,
                                         CubitVector* closest_point );
  virtual CubitBoolean is_periodic();
  virtual CubitBoolean is_periodic_in_U( double& period );
  virtual CubitBoolean is_periodic_in_V( double& period );
  virtual CubitBoolean is_singular_in_U( double u_value );
  virtual CubitBoolean is_singular_in_V( double v_value );
  virtual CubitBoolean is_closed_in_U();
  virtual CubitBoolean is_closed_in_V();
  virtual CubitStatus uv_derivitives( double u, double v,
                                      CubitVector& du, 
                                      CubitVector& dv );
  virtual CubitBoolean get_param_range_U( double& lower, double& upper );
  virtual CubitBoolean get_param_range_V( double& lower, double& upper );
		
  virtual CubitStatus get_sphere_params( CubitVector &center,
                                         double &radius ) const;
    //- Only valid for spherical surfaces
    //O center
    //O- The center of the sphere
    //O radius
    //O- The radius of the sphere

  virtual CubitStatus get_cone_params( CubitVector &center,
                                       CubitVector &normal,
                                       CubitVector &major_axis,
                                       double &radius_ratio,
                                       double &sine_angle,
                                       double &cos_angle ) const;
    //- Only valid for conical surfaces.  Cylinders are a special case of conicals.
    //O center
    //O- 
    //O normal
    //O- 
    //O major_axis
    //O- 
    //O radius_ratio
    //O- 
    //O sine_angle
    //O- 
    //O cos_angle
    //O- 

    virtual CubitStatus get_torus_params( CubitVector &center,
                                        CubitVector &normal,
                                        double &major_radius,
                                        double &minor_radius ) const;
    //- Only valid for torus surfaces.
    //O center
    //O- 
    //O normal
    //O- 
    //O major_radius
    //O- 
    //O minor_radius
    //O- 

  virtual CubitStatus get_nurb_params( bool &rational,
                                       int &degree_u,
                                       int &degree_v,
                                       int &num_cntrl_pts_u,
                                       int &num_cntrl_pts_v,
                                       DLIList<CubitVector> &cntrl_pts,
                                       DLIList<double> &weights,
                                       DLIList<double> &u_knots,
                                       DLIList<double> &v_knots ) const;
  //- Only valid for nurbs surfaces
  //O rational
  //O-   True if the nurb is rational
  //O degree_u
  //O-   The degree of the nurb in the u direction
  //O degree_v
  //O-   The degree of the nurb in the v direction
  //O num_cntrl_pts_u
  //O-   Number of control points in the u direction
  //O num_cntrl_pts_v
  //O-   Number of control points in the v direction
  //O cntrl_pts
  //O-   The control points stored as
  //O-           cntrl_pts[0                ] = pt[u=0][v=0]
  //O-           cntrl_pts[1                ] = pt[u=1][v=0]
  //O-               ...
  //O-           cntrl_pts[num_cntrl_pts_u-1] = pt[u=?][v=0]
  //O-           cntrl_pts[num_cntrl_pts_u  ] = pt[u=0][v=1]
  //O-               ...
  //O weights
  //O-   If rational, weights for each control point, stored in the same
  //O-   order as the control points.  No weights are returned if
  //O-   rational == false
  //O u_knots
  //O-   knot vector in the u direction
  //O v_knots
  //O-   knot vector in the v direction

  void get_hidden_curves( DLIList<Curve*>& curves );
  
  void print_debug_info( const char* line_prefix = 0, bool brief = false );
                         
  int closest_underlying_surface( const CubitVector& position );
  
  static void print_cpt_stats();
  static void reset_cpt_stats();
  
  void read_attributes() { compGeom->read_attributes(); }
  void write_attributes() { compGeom->write_attributes(); }
  
  
  CubitStatus get_graphics( GMem& gmem );

  // handles transformation of composite surface
  // child loops, coedges, etc. must be handled separately
  void notify_transformed();
  
protected:
  
  void update_facet_tool();
  void update_facets_to_ignore();
  CubitStatus closest_trimmed( int underlying_surface,
                               const CubitVector& position, 
                               CubitVector& result );
	
  void update_modified( Surface* modified_surface, 
                        DLIList<CompositeCurve*>& new_curves );
  void update_modified();
  bool     is_dead_coedge( CompositeCoEdge* coedge );
  void remove_dead_coedge( CompositeCoEdge* coedge );

private:
  int HadBridgeRemoved;

  DLIList<Surface*> surfacesToIgnore;	

    // these have no implementation, just private delcarations
    // to prevent the compiler from generating default implementations
  CompositeSurface& operator=(const CompositeSurface&);
  CompositeSurface(const CompositeSurface&);

  CompositeGeom *compGeom;
  CompositeSurface *stitchPartner;
  
  CompositeCoSurf* firstCoSurf;

  CompositeLoop* firstLoop;

  HiddenEntitySet* hiddenSet;
  
  CompSurfFacets* facetTool;
};


inline Surface* CompositeSurface::get_surface( int index )  const
  { return dynamic_cast<Surface*>(compGeom->entity(index)); }

inline CubitSense CompositeSurface::get_sense( int index ) const 
  { return compGeom->sense(index); }

inline int CompositeSurface::num_surfs() const
  { return compGeom->num_entities(); }

inline int CompositeSurface::index_of( Surface* surf ) const
  { return compGeom->index_of( surf ); }

inline void CompositeSurface::update()
  { compGeom->update_cached_data(); }
inline CompositeLoop* CompositeSurface::first_loop( ) const
  { return firstLoop; }

inline CompositeLoop* CompositeSurface::next_loop( CompositeLoop* loop ) const
  { return !loop ? firstLoop : loop->mySurface == this ? loop->loopNext : 0; }

inline CompositeCoSurf* CompositeSurface::next_co_surface( CompositeCoSurf* prev ) const
  { return prev ? prev->surfaceNext : firstCoSurf; }

inline HiddenEntitySet& CompositeSurface::hidden_entities()
{ 
  if( !hiddenSet )
    hiddenSet = new HiddenEntitySet(this);
  return *hiddenSet;
}


#endif

