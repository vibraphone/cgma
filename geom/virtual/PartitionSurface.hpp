//-------------------------------------------------------------------------
// Filename      : PartitionSurface.hpp
//
// Purpose       : 
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 04/21/02
//-------------------------------------------------------------------------

#ifndef PARTITION_SURFACE_HPP
#define PARTITION_SURFACE_HPP

#include <set>
#include "Surface.hpp"
#include "PartitionEntity.hpp"
#include "PartitionCoSurf.hpp"
#include "PartitionLoop.hpp"

class FacetEvalTool;
class CubitFacetData;
class PartitionShell;
class PartitionLump;
class CubitFacetEdge;
class CubitFacet;
class CubitFacetEdgeData;
class CubitPointData;

class PartitionSurface : public Surface, 
                         public PartitionEntity
{

public:

  virtual CubitStatus combine( PartitionSurface* dead_surface );
  
  static PartitionSurface* construct( const CubitSimpleAttrib& attrib, PartitionLump* vol );

  PartitionSurface( PartitionLump* owner );
  
  PartitionSurface* split( DLIList<CubitFacetData*>& facets_to_move );

  virtual ~PartitionSurface();
  
  int num_loops() const;
  PartitionLoop* next_loop( const PartitionLoop* prev = 0 ) const;
  
  CubitStatus add( PartitionLoop* loop );
  CubitStatus remove( PartitionLoop* loop );
  
  CubitStatus add( PartitionCoSurf* cosurf );
  CubitStatus remove( PartitionCoSurf* cosurf );
  
  PartitionCoSurf* next_co_surface( const PartitionCoSurf* prev = 0 ) const;
  PartitionCoSurf* find_first( const PartitionShell* shell ) const;
    // Get the first CoSurface connecting this surface w/ the passed shell
  PartitionCoSurf* find_first( const PartitionLump* lump ) const;
    // Get the first CoSurface connecting this surface w/ a
    // shell owned by the passed volume
  PartitionCoSurf* find_next( const PartitionCoSurf* cosurf ) const;
    // Get the next CoSurface after the passed CoSurface and
    // with the same Shell as the passed CoSurface.
  PartitionCoSurf* find_first( const PartitionShell* shell, CubitSense sense ) const;
    // Get the first CoSurface w/ the passed shell and sense
  
  void get_points( DLIList<PartitionPoint*>& list ) const;
  
  virtual CubitStatus save( CubitSimpleAttrib& attrib );
  
  CubitBox bounding_box() const;

  double measure();
  GeometryType geometry_type();


  virtual CubitStatus move_to_geometry( CubitVector& position );
  virtual void closest_point_trimmed( CubitVector from, CubitVector& closest);
  virtual CubitStatus get_point_normal( CubitVector& point, CubitVector& normal );

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

  virtual CubitStatus closest_point( CubitVector const& pos, 
                             CubitVector* close = 0, CubitVector* norm = 0,
                             CubitVector* curv1 = 0, CubitVector* curv2 = 0);
  virtual CubitStatus principal_curvatures( CubitVector const& loc, 
                                    double& curv1, double& curv2,
                                    CubitVector* closest_location=0 );

  virtual CubitVector position_from_u_v( double u, double v );
  virtual CubitStatus u_v_from_position( CubitVector const& location,
                                 double& u, double& v,
                                 CubitVector* closest_location = NULL );
  
  virtual CubitStatus closest_point_along_vector( CubitVector& from_point, 
                                                  CubitVector& along_vector,
                                                  CubitVector& point_on_surface);
  virtual CubitBoolean is_periodic();
  virtual CubitBoolean is_periodic_in_U( double& period );
  virtual CubitBoolean is_periodic_in_V( double& period );
  virtual CubitBoolean is_singular_in_U( double u_param );
  virtual CubitBoolean is_singular_in_V( double v_param );
  virtual CubitBoolean is_closed_in_U();
  virtual CubitBoolean is_closed_in_V();

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

  virtual CubitStatus uv_derivitives( double u, double v,
                              CubitVector &du, CubitVector &dv );
  virtual CubitBoolean is_parametric();
  virtual CubitBoolean get_param_range_U( double& lo, double& hi );
  virtual CubitBoolean get_param_range_V( double& lo, double& hi );
 
  virtual CubitBoolean is_position_on( CubitVector &test_position );

  virtual CubitPointContainment point_containment( const CubitVector &point );
  virtual CubitPointContainment point_containment( double u, double v );
//  virtual CubitPointContainment point_containment( CubitVector &point, 
//                                           double u, double v );
  CubitPointContainment point_containment( const CubitVector &point, 
                                           PartitionCurve*& boundary_curve );

  CubitSense get_geometry_sense();
  void reverse_sense();

  
  void get_parents_virt( DLIList<TopologyBridge*>& );
  void get_children_virt( DLIList<TopologyBridge*>& );
  int layer() const { return sub_entity_set().get_owner_layer(); }
  GeometryQueryEngine* get_geometry_query_engine() const;
  
  virtual CubitSense get_shell_sense( ShellSM* shell_ptr ) const;

  void append_simple_attribute_virt(const CubitSimpleAttrib&);
  void remove_simple_attribute_virt(const CubitSimpleAttrib&);
  void remove_all_simple_attribute_virt();
  CubitStatus get_simple_attribute(DLIList<CubitSimpleAttrib>&);
  CubitStatus get_simple_attribute( const CubitString& name,
                                    DLIList<CubitSimpleAttrib>& );
  
  void notify_split( FacetEntity* old_tri, FacetEntity* new_tri );
  void notify_destroyed( CubitFacetData* facet );
  void get_facet_data( DLIList<CubitFacetData*>& result_list ) const;
  void set_facet_data( const DLIList<CubitFacetData*>& new_list );
  void replace_facets(DLIList<CubitFacetData*> &dead_facets,
                      DLIList<CubitFacetData*> &new_facets);
  bool has_facets() const { return facetList.size() > 0; }
  void draw_facets( int color ) const;


  CubitStatus init_facet_data();
  

  virtual void transform( const CubitTransformMatrix& xform );

  virtual void print_debug_info( const char* prefix = 0, 
                                 bool print_sub_entity_set = true ) const;
  CubitFacet* closest_facet( const CubitVector& input_position,
                             CubitVector& result_position);
  //- Call the following function for the facets owned by this surface.
                             
  CubitStatus notify_moving_point( CubitPoint* point, 
                                   const CubitVector& new_pos );
  //- Call fix_move_point for facetts of this PartitionSurface, 
  //- and update the facet list for this surface accordingly.

protected:
 
  void interior_facet_points( DLIList<CubitPoint*>& result_list ) const;
  
  void update_facet_tool();

  void reverse_loops();
  
  PartitionSurface( );
  
  virtual PartitionSurface* copy();

  CubitStatus get_save_topology( DLIList<int>& curves );
  
  
  CubitStatus insert_facets( DLIList <PartitionEntity*> points,
                             DLIList <CubitFacetEdgeData *> &boundary_edges,
                             DLIList <CubitFacetData*> &surf_facets);
  
  DLIList<CubitFacetData*> facetList;
  CubitSense geometry_sense;

private:
  
  PartitionSurface( PartitionSurface* split_from );

  PartitionSurface( const PartitionSurface& );

  PartitionLoop* firstLoop;
  PartitionCoSurf* firstCoSurf;

  
};

inline PartitionLoop* 
PartitionSurface::next_loop( const PartitionLoop* prev ) const
{ 
  return !prev                   ? firstLoop           : 
         prev->mySurface == this ? prev->nextInSurface : 
                                   0; 
}
  
inline PartitionCoSurf* 
PartitionSurface::next_co_surface( const PartitionCoSurf* prev ) const
{ 
  return !prev                   ? firstCoSurf       :
         prev->mySurface == this ? prev->surfaceNext :
                                   0;
}

#endif
