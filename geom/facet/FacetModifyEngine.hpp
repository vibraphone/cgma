//-------------------------------------------------------------------------
// Filename      : FacetModifyEngine.hpp
//
// Purpose       : ModifyEngine for faceted geometry
//
// Special Notes : Modeled after GeometryModifyEngine and AcisModifyEngine.
//
// Creator       : John Fowler
//
// Creation Date : 6/02
//
// Owner         : John Fowler
//-------------------------------------------------------------------------

#ifndef FACET_MODIFY_ENGINE_HPP
#define FACET_MODIFY_ENGINE_HPP

#include "GeometryModifyEngine.hpp"
#include <vector>

class Point;
class TopologyBridge;
class TopologyEntity;
class ChollaSurface;
class ChollaCurve;
class ChollaPoint;
class ChollaEngine;
class CurveFacetEvalTool;
class FacetEvalTool;
class CubitPoint;
class CubitFacet;
class CubitFacetEdge;
class CubitQuadFacet;
class Point;
class CoEdgeSM;
class ShellSM;
class FacetBody;
class Body;
class RefFace;
class FacetSurface;
class CubitEvaluatorData;
class SphereEvaluatorData;
class CylinderEvaluatorData;

class FacetModifyEngine : public GeometryModifyEngine
{
  
public:
  
  //HEADER- Constructor and Destructor
private:
  
  FacetModifyEngine();
  
  static FacetModifyEngine* instance_;
  
public:
  
  static inline FacetModifyEngine* instance()
  {
    if( !instance_ )
      instance_ = new FacetModifyEngine;
    return instance_;
  }
  
  static void initialize_settings();
  
  virtual ~FacetModifyEngine();
  //- virtual destructor
  
  virtual Point* make_Point( CubitVector const& point) const ;
  
  virtual Curve* make_Curve(Curve *curve_ptr) const;
  //- creates a curve from an existing curve.  This creates totally
  //- new topology.  This function is useful for constructing geometry
  //- from existing geometry.
  
  virtual Curve* make_Curve( Point const* point1_ptr,
    Point const* point2_ptr,
    Surface* ref_face_ptr,
    const CubitVector *third_point = NULL) const;
  //- Create a curve exactly on the give ref_face.
  //- Make sure the points are on the underlying surface.
  
  virtual Curve* make_Curve( GeometryType curve_type,
    Point const* point1_ptr,
    Point const* point2_ptr,
    DLIList<CubitVector*>& vector_list,
    Surface* ref_face_ptr = NULL) const;
  
  virtual Curve* make_Curve( GeometryType curve_type,
    Point const* point1_ptr,
    Point const* point2_ptr,
    CubitVector const* intermediate_point_ptr,
    CubitSense sense) const;
  
  
  virtual Surface* make_Surface( Surface *old_surface_ptr,
    CubitBoolean extended_from = CUBIT_FALSE) const;
  
  virtual Surface* make_Surface( GeometryType surface_type,
    DLIList<Curve*>& curve_list,
    Surface *old_surface_ptr = NULL,
    bool check_edges = true ) const;
  
  virtual Lump* make_Lump( DLIList<Surface*>& surface_list ) const;
  
  //virtual Body *make_Body(Surface *) const;
  
  virtual BodySM* make_BodySM( Surface * ) const;
  
  virtual BodySM* make_BodySM( DLIList<Lump*>& /*lump_list*/ ) const;
  
  //virtual Chain *make_Chain(Curve *) const;
  
  virtual BodySM* sphere(double radius) const ;
  
  virtual BodySM* brick( double wid, double dep, double hi ) const ;
  
  virtual BodySM* brick( const CubitVector &center, 
                         const CubitVector axes[3],
                         const CubitVector &extension) const ;
  
  virtual BodySM* prism( double height, int sides, double major,
    double minor) const ;
  
  virtual BodySM* pyramid( double height, int sides, double major,
    double minor, double top=0.0) const ;
  
  virtual BodySM* cylinder( double hi, double r1, double r2,
    double r3 ) const ;
  
  virtual BodySM* torus( double r1, double r2 ) const ;
  
  virtual BodySM* planar_sheet ( const CubitVector& p1,
    const CubitVector& p2,
    const CubitVector& p3,
    const CubitVector& p4 ) const;
  
  virtual BodySM* copy_body ( BodySM* bodyPtr) const ;
  
  virtual CubitStatus stitch_surfs(
		      DLIList<BodySM*>& surf_bodies,
		      BodySM*& stitched_body)const;

  virtual CubitStatus subtract(DLIList<BodySM*> &tool_body_list,
    DLIList<BodySM*> &from_bodies,
    DLIList<BodySM*> &new_bodies,
    bool imprint = false,
    bool keep_old = false) const;
  
  virtual CubitStatus imprint(BodySM* BodyPtr1, BodySM* BodyPtr2,
    BodySM*& newBody1, BodySM*& newBody2,
    bool keep_old) const;
  
  virtual CubitStatus imprint(DLIList<BodySM*> &from_body_list,
    DLIList<BodySM*> &new_from_body_list,
    bool keep_old, DLIList<TopologyBridge*> *new_tbs = NULL,
    DLIList<TopologyBridge*> *att_tbs = NULL) const;
  
  virtual CubitStatus imprint( DLIList<BodySM*> &body_list,
    DLIList<Curve*> &ref_edge_list,
    DLIList<BodySM*>& new_body_list,
    bool keep_old_body,
    bool show_messages=CUBIT_TRUE) const;
  
  virtual CubitStatus imprint( DLIList<Surface*> &ref_face_list,
    DLIList<Curve*> &ref_edge_list,
    DLIList<BodySM*>& new_body_list,
    bool keep_old_body ) const;

  virtual CubitStatus imprint( DLIList<Surface*> &surface_list,
    DLIList<DLIList<Curve*>*> &curve_lists_list,
    BodySM*& new_body,
    bool keep_old_body ) const;
  
  virtual CubitStatus imprint( DLIList<BodySM*> &body_list,
    DLIList<CubitVector*> &vector_list,
    DLIList<BodySM*>& new_body_list,
    bool keep_old_body,
    DLIList<TopologyBridge*> *new_tbs = NULL,
    DLIList<TopologyBridge*> *att_tbs = NULL ) const;
  
  virtual CubitStatus imprint_projected_edges( DLIList<Surface*> &ref_face_list,
    DLIList<Curve*> &ref_edge_list,
    DLIList<BodySM*>& new_body_list,
    bool keep_old_body,
    bool keep_free_edges) const;
  
  virtual CubitStatus imprint_projected_edges(DLIList<Surface*> &ref_face_list,
    DLIList<BodySM*> &body_list,
    DLIList<Curve*> &ref_edge_list,
    DLIList<BodySM*>& new_body_list,
    bool keep_old_body,
    bool keep_free_edges) const;
  
  virtual CubitStatus project_edges( DLIList<Surface*> &ref_face_list,
    DLIList<Curve*> &ref_edge_list_in,
    DLIList<Curve*> &ref_edge_list_new,
    bool print_error = true ) const;
  
  virtual CubitStatus intersect(BodySM* tool_body_ptr,
    DLIList<BodySM*> &from_bodies,
    DLIList<BodySM*> &new_bodies,
    bool keep_old = CUBIT_FALSE) const;
  
  virtual CubitStatus chop(DLIList<BodySM*> &bodies, 
                           DLIList<BodySM*> &intersectBodies,
                           DLIList<BodySM*> &outsideBody, 
                           BodySM*& leftoversBody,
                           bool keep_old = CUBIT_FALSE,
                           bool nonreg = CUBIT_FALSE) const;
  
  virtual CubitStatus unite(DLIList<BodySM*> &bodies, 
                            DLIList<BodySM*> &newBodies,
    bool keep_old = CUBIT_FALSE) const;
  
  virtual CubitStatus thicken(DLIList<BodySM*>& bodies, 
    DLIList<BodySM*>& new_bodies,
    double depth,
    bool both = CUBIT_FALSE) const ;

  virtual CubitStatus flip_normals( DLIList<Surface*>& face_list ) const;
  
  virtual CubitStatus  sweep_translational(
    DLIList<GeometryEntity*>& ref_ent_list,
    DLIList<BodySM*>& result_body_list,
    const CubitVector& sweep_vector,
    double draft_angle,
    int draft_type,
    bool switchside,
    bool rigid,
    Surface *stop_surf = NULL,
    BodySM  *stop_body = NULL) const;
  
  virtual CubitStatus  sweep_perpendicular(
    DLIList<GeometryEntity*>& ref_ent_list,
    DLIList<BodySM*>& result_body_list,
    double distance,
    double draft_angle,
    int draft_type,
    bool switchside,
    bool rigid,
    Surface* stop_surf = NULL,
    BodySM  *stop_body = NULL ) const;
  
  virtual CubitStatus  sweep_rotational(
    DLIList<GeometryEntity*>& ref_ent_list,
    DLIList<BodySM*>& result_body_list,
    const CubitVector& point,
    const CubitVector& direction,
    double angle,
    int steps = 0,
    double draft_angle = 0.0,
    int draft_type = 0,
    bool switchside = CUBIT_FALSE,
    bool make_solid = CUBIT_FALSE,
    bool rigid = CUBIT_FALSE , 
    Surface *stop_surf = NULL,
    BodySM  *stop_body = NULL) const;
  
  virtual CubitStatus sweep_along_curve(
    DLIList<GeometryEntity*>& ref_ent_list,
    DLIList<BodySM*>& result_body_list,
    DLIList<Curve*>& ref_edge_list,
    double draft_angle = 0.0,
    int draft_type = 0,
    bool rigid = CUBIT_FALSE,
    Surface *stop_surf = NULL,
    BodySM  *stop_body = NULL ) const;
 
    virtual CubitStatus webcut_with_sweep_surfaces(
                          DLIList<BodySM*> &blank_bodies,
                          DLIList<Surface*> &surfaces,
                          const CubitVector& sweep_vector,
                          bool sweep_perp, 
                          bool through_all,
                          bool outward,
                          bool up_to_next, 
                          Surface *stop_surf, 
                          Curve *curve_to_sweep_along, 
                          DLIList<BodySM*> &results_list,
                          CubitBoolean imprint = false);

    virtual CubitStatus webcut_with_sweep_curves(
                          DLIList<BodySM*> &blank_bodies,
                          DLIList<Curve*> &curves,
                          const CubitVector& sweep_vector,
                          bool through_all, 
                          Surface *stop_surf, 
                          Curve *curve_to_sweep_along, 
                          DLIList<BodySM*> &results_list,
                          CubitBoolean imprint = false);

    virtual CubitStatus webcut_with_sweep_curves_rotated(
                          DLIList<BodySM*> &blank_bodies,
                          DLIList<Curve*> &curves,
                          const CubitVector &point,
                          const CubitVector &sweep_axis,
                          double angle,
                          Surface *stop_surf, 
                          DLIList<BodySM*> &results_list,
                          CubitBoolean imprint = false);

    virtual CubitStatus webcut_with_sweep_surfaces_rotated(
                            DLIList<BodySM*> &blank_bodies,
                            DLIList<Surface*> &surfaces,
                            const CubitVector &point, 
                            const CubitVector &sweep_axis, 
                            double angle, 
                            Surface *stop_surf, 
                            bool up_to_next, 
                            DLIList<BodySM*> &results_list,
                            CubitBoolean imprint = false); 

  //HEADER- Webcut-related functions
  virtual CubitStatus webcut(DLIList<BodySM*>& webcut_body_list,
    const CubitVector &v1,
    const CubitVector &v2,
    const CubitVector &v3,
    DLIList<BodySM*>& results_list,
    bool imprint = false 
    ) const;
  
  virtual CubitStatus webcut(DLIList<BodySM*>& webcut_body_list,
    BodySM const* tool_body,
    DLIList<BodySM*>& results_list,
    bool imprint = false)const;
  
  virtual CubitStatus webcut_across_translate( DLIList<BodySM*>& body_list,
    Surface* plane_surf1,
    Surface* plane_surf2,
    DLIList<BodySM*>& results_list,
    bool imprint = false ) const;
  
  virtual CubitStatus webcut_with_sheet(DLIList<BodySM*> &webcut_body_list,
    
    BodySM *sheet_body,
    DLIList<BodySM*> &new_bodies,
    bool imprint = false );
  
  virtual CubitStatus webcut_with_extended_surf(DLIList<BodySM*> &webcut_body_list,
    Surface *extend_from,
    DLIList<BodySM*> &new_bodies,
    int &num_cut,
    bool imprint = false );
  
  virtual CubitStatus webcut_with_cylinder(DLIList<BodySM*> &webcut_body_list,
    double radius,
    const CubitVector &axis,
    const CubitVector &center,
    DLIList<BodySM*>& results_list,
    bool imprint = false );
  
/*  virtual CubitStatus webcut_with_brick( DLIList<BodySM*>& webcut_body_list, 
    const CubitVector &center,
    const CubitVector axes[3], 
    const CubitVector &extension,
    DLIList<BodySM*> &results_list,
    bool imprint = false );
*/
  
  virtual CubitStatus webcut_with_planar_sheet( DLIList<BodySM*>& webcut_body_list,
    const CubitVector &center,
    const CubitVector axes[2],
    double width, double height,
    DLIList<BodySM*> &results_list,
    bool imprint = false );
  
  virtual CubitStatus webcut_with_curve_loop(DLIList<BodySM*> &webcut_body_list,
    DLIList<Curve*> &ref_edge_list,
    DLIList<BodySM*>& results_list,
    bool imprint = false );
  
  virtual CubitStatus section( DLIList<BodySM*> &section_body_list,
    const CubitVector &point_1,
    const CubitVector &point_2,
    const CubitVector &point_3,
    DLIList<BodySM*>& new_body_list,
    bool keep_normal_side,
    bool keep_old = false,
    bool keep_both_sides = false);
  
  virtual CubitStatus split_body( BodySM *body_ptr,
    DLIList<BodySM*> &new_bodies );
  
  virtual CubitStatus reverse_body( BodySM *body_to_reverse );
  
  virtual CubitStatus split_periodic( BodySM *body_ptr,
    BodySM *&new_body );
  
  virtual CubitStatus regularize_body( BodySM *body_ptr,
    BodySM *&new_body_ptr );
  
  virtual CubitStatus regularize_entity(GeometryEntity *old_entity_ptr,  
                                           BodySM *&new_body_ptr);
  
  virtual CubitStatus offset_curves( DLIList<Curve*>& ref_edge_list, 
                                     DLIList<Curve*>& result_curve_list,
                                     double offset_distance,
                                     const CubitVector& offset_direction, 
                                     int gap_type = 1 );
  
  virtual CubitStatus scale ( BodySM *&body, const CubitVector& factors );

  virtual Curve* trim_curve( Curve* trim_curve, 
                             const CubitVector& trim_vector,
                             const CubitVector& keep_vector,
                             bool keep_old = false );
  
  virtual CubitStatus create_solid_bodies_from_surfs(DLIList<Surface*> &ref_face_list,
                                                     DLIList<BodySM*> &new_bodies,
                                                     bool keep_old = false,
                                                     bool heal = false) const;
  
  virtual Curve* create_arc_three( Point* ref_vertex1, 
                                   Point* ref_vertex2,
                                   Point* ref_vertex3, 
                                   bool full = false );
  
  virtual Curve* create_arc_three( Curve* ref_edge1, 
                                   Curve* ref_edge2,
                                   Curve* ref_edge3, 
                                   bool full = false );
  
  virtual Curve* create_arc_center_edge( 
                                   Point* ref_vertex1, 
                                   Point* ref_vertex2,
                                   Point* ref_vertex3, 
                                   const CubitVector &normal,
                                   double radius = CUBIT_DBL_MAX,
                                   bool full = false );
  
  virtual CubitStatus create_curve_combine( DLIList<Curve*>& curve_list, 
                                    Curve *&new_curve_ptr );
    //-  Uses the solid modeller to create a new RefEdge that is a combination 
    //- of the input chain of edges.  
    //-

  virtual GeometryQueryEngine *get_gqe();
  
  virtual CubitBoolean is_modify_engine(const TopologyBridge *tb_ptr) const;
  
  virtual CubitStatus get_offset_intersections( 
                               Curve* ref_edge1, 
                               Curve* ref_edge2,
                               DLIList<CubitVector*>& intersection_list,
                               double offset,
                               bool ext_first = true );
  
  virtual CubitStatus get_offset_intersections( 
                               Curve* ref_edge_ptr, 
                               Surface* ref_face_ptr,
                               DLIList<CubitVector*> &intersection_list,
                               double offset = 0.0,
                               bool ext_surf = true );
  
  virtual CubitStatus surface_intersection( Surface *surface1_ptr,
    Surface *surface2_ptr,
    DLIList<Curve*> &inter_graph,
    const double tol) const;
  
  virtual CubitStatus get_mid_plane( const CubitVector &point_1,
    const CubitVector &point_2,
    const CubitVector &point_3,
    BodySM *body_to_trim_to,
    BodySM *&midplane_body ) const;

  virtual CubitStatus get_spheric_mid_surface( Surface *surface_ptr1,
					       Surface *surface_ptr2,
					       BodySM *body_to_trim_to,
					       BodySM *&midsurface_body ) const;
  
  virtual CubitStatus get_conic_mid_surface( Surface *surface_ptr1,
					       Surface *surface_ptr2,
					       BodySM *body_to_trim_to,
					       BodySM *&midsurface_body ) const;
  
  virtual CubitStatus get_toric_mid_surface( Surface *surface_ptr1,
					       Surface *surface_ptr2,
					       BodySM *body_to_trim_to,
					       BodySM *&midsurface_body ) const;
  
  virtual CubitStatus tweak_chamfer( DLIList<Curve*> &curve_list, 
                                     double left_offset,
                                     DLIList<BodySM*> &new_bodysm_list,
                                     double right_offset = -1.0,
                                     CubitBoolean keep_old_body = CUBIT_FALSE,
                                     CubitBoolean preview = CUBIT_FALSE ) const;
  /**<  Chamfer curves on solid bodies.  The left and right offsets are with 
    *   respect to the curve direction.  If the given right offset is negative,
    *   the left offset is used.  Users can preview to clarify the meaning of
    *   left and right.
    */

  virtual CubitStatus tweak_chamfer( DLIList<Point*> &point_list, 
                                     double offset1, 
                                     DLIList<BodySM*> &new_bodysm_list,
                                     Curve *edge1 = NULL,
                                     double offset2 = -1.0,
                                     Curve *edge2 = NULL,
                                     double offset3 = -1.0,
                                     Curve *edge3 = NULL,
                                     CubitBoolean keep_old_body = CUBIT_FALSE,
                                     CubitBoolean preview = CUBIT_FALSE ) const;
  /**<  Chamfer vertices on solid or sheet bodies.  On a solid body there can
    *   be up to 3 offsets; on a sheet body up to 2 offsets.  The offsets are
    *   in the direction of the supplied edges.  If multiple vertices are 
    *   supplied, only one offset value is allowed and the edges are not used.
    */

  virtual CubitStatus tweak_fillet( DLIList<Curve*> &curve_list, 
                                    double radius,
                                    DLIList<BodySM*> &new_bodysm_list,
                                    CubitBoolean keep_old_body = CUBIT_FALSE,
                                    CubitBoolean preview = CUBIT_FALSE ) const;
  /**<  Create a round fillet (or blend) at the given curves on solid bodies.
    */

  virtual CubitStatus tweak_fillet( Curve *curve_ptr, 
                                    double start_radius,
                                    double end_radius,
                                    BodySM *&new_body_ptr,
                                    CubitBoolean keep_old_body = CUBIT_FALSE,
                                    CubitBoolean preview = CUBIT_FALSE ) const;
  /**<  Create a round fillet (or blend) at the given curve on a solid body.
    *   The fillet has a variable radius from the start to the end of the curve.
    */

  virtual CubitStatus tweak_fillet( DLIList<Point*> &point_list, 
                                    double radius,
                                    DLIList<BodySM*> &new_bodysm_list,
                                    CubitBoolean keep_old_body = CUBIT_FALSE,
                                    CubitBoolean preview = CUBIT_FALSE ) const;
  /**<  Create a round fillet (or blend) at the given vertices on sheet bodies.
    */

  virtual CubitStatus tweak_move( DLIList<Surface*> &surface_list,
                                  const CubitVector &delta,
                                  DLIList<BodySM*> &new_bodysm_list,
                                  CubitBoolean keep_old_body = CUBIT_FALSE,
                                  CubitBoolean preview = CUBIT_FALSE ) const;
  /**<  Tweak specified faces of a volume or volumes along a vector.
    */

  virtual CubitStatus tweak_move( DLIList<Curve*> &curve_list,
                                  const CubitVector &delta,
                                  DLIList<BodySM*> &new_bodysm_list,
                                  CubitBoolean keep_old_body = CUBIT_FALSE,
                                  CubitBoolean preview = CUBIT_FALSE ) const;
  /**<  Tweak specified curves of a sheet body along a vector.
    */

  virtual CubitStatus tweak_offset( DLIList<Surface*> &surface_list,
                                    double offset_distance,
                                    DLIList<BodySM*> &new_bodysm_list,
                                    CubitBoolean keep_old_body = CUBIT_FALSE,
                                    CubitBoolean preview = CUBIT_FALSE ) const;
  /**<  Tweak specified faces of a volume or volumes by offsetting those faces
    *   by the offset distance.
    */

  virtual CubitStatus tweak_offset( DLIList<Curve*> &curve_list,
                                    double offset_distance,
                                    DLIList<BodySM*> &new_bodysm_list,
                                    CubitBoolean keep_old_body = CUBIT_FALSE,
                                    CubitBoolean preview = CUBIT_FALSE ) const;
  /**<  Tweak specified curves of a sheet body or bodies by offsetting those
    *   curves by the offset distance.
    */

  virtual CubitStatus tweak_remove( DLIList<Surface*> &surface_list,
                                    DLIList<BodySM*> &new_bodysm_list,
                                    CubitBoolean extend_adjoining = CUBIT_TRUE,
                                    CubitBoolean keep_surface = CUBIT_FALSE,
                                    CubitBoolean keep_old_body = CUBIT_FALSE,
                                    CubitBoolean preview = CUBIT_FALSE ) const;
  /**<  Remove surfaces from a body or bodies and then extend the adjoining
    *   surfaces to fill the gap or remove the hole.
    */

  virtual CubitStatus tweak_remove( DLIList<Curve*> &curve_list,
                                    DLIList<BodySM*> &new_bodysm_list, 
                                    CubitBoolean keep_old_body = CUBIT_FALSE,
                                    CubitBoolean preview = CUBIT_FALSE ) const;
  /**<  Remove curves from a sheet body or bodies and then extend the remaining
    *   curves to fill the gap.  If an internal loop of curves is removed the
    *   hole is removed.
    */

  virtual CubitStatus tweak_target( DLIList<Surface*> &surface_list,
                                    DLIList<Surface*> &target_surf_list,
                                    DLIList<BodySM*> &new_bodysm_list,
                                    CubitBoolean reverse_flg = CUBIT_FALSE,
                                    CubitBoolean keep_old_body = CUBIT_FALSE,
                                    CubitBoolean preview = CUBIT_FALSE ) const;
  /**<  Tweak specified faces of a volume or volumes up to target surfaces.
    */

  virtual CubitStatus tweak_target( DLIList<Curve*> &curve_list,
                                    DLIList<Surface*> &target_surf_list, 
                                    DLIList<BodySM*> &new_bodysm_list,
                                    CubitBoolean reverse_flg = CUBIT_FALSE,
                                    CubitBoolean keep_old_body = CUBIT_FALSE,
                                    CubitBoolean preview = CUBIT_FALSE ) const;
  /**<  Tweak specified edges of a surface or set of surfaces (in sheet
    *   bodies) up to a set of target surfaces.  This essentially extends or
    *   trims the attached surfaces of the sheet body.
    */

  virtual CubitStatus tweak_target( DLIList<Curve*> &curve_list,
                                    DLIList<Curve*> &target_curve_list, 
                                    DLIList<BodySM*> &new_bodysm_list,
                                    CubitBoolean reverse_flg = CUBIT_FALSE,
                                    CubitBoolean keep_old_body = CUBIT_FALSE,
                                    CubitBoolean preview = CUBIT_FALSE ) const;
  /**<  Tweak specified edges of a sheet body or bodies up to a list of target
    *   curves that are part of a sheet body.  The target is a surface created by
    *   thickening the owning surfaces of the target curves.
    */
  
  virtual CubitStatus remove_curve_slivers( BodySM *body, double lengthlimit ) const;

  virtual CubitStatus create_net_surface( DLIList<Surface*>& ref_face_list, BodySM *& new_body,
                                          DLIList<DLIList<CubitVector*>*> &vec_lists_u, 
                                          DLIList<DLIList<CubitVector*>*> &vec_lists_v, 
                                          double net_tol = 1e-3,
                                          CubitBoolean heal = CUBIT_TRUE ) const;

  virtual CubitStatus create_net_surface( DLIList<Curve*>& u_curves, DLIList<Curve*>& v_curves,
                                          BodySM *& new_body,
                                          double net_tol = 1e-3, 
                                          CubitBoolean heal = CUBIT_TRUE ) const;

  virtual CubitStatus create_offset_surface( Surface* ref_face_ptr, BodySM*& new_body, double offset_distance ) const;

  virtual CubitStatus create_offset_body( BodySM* body_ptr, BodySM*& new_body, double offset_distance ) const;

  virtual CubitStatus create_skin_surface( DLIList<Curve*>& curves, BodySM*& new_body ) const;

  virtual CubitStatus loft_surfaces( Surface *face1, const double &takeoff1,
                                     Surface *face2, const double &takeoff2,
                                     BodySM*& new_body,
                                     CubitBoolean arc_length_option = CUBIT_FALSE,
                                     CubitBoolean twist_option = CUBIT_FALSE,
                                     CubitBoolean align_direction = CUBIT_TRUE,
                                     CubitBoolean perpendicular = CUBIT_TRUE,
                                     CubitBoolean simplify_option = CUBIT_FALSE) const;

  virtual CubitStatus loft_surfaces_to_body( Surface *face1, const double &takeoff1,
                                             Surface *face2, const double &takeoff2,
                                             BodySM*& new_body,
                                             CubitBoolean arc_length_option,
                                             CubitBoolean twist_option,
                                             CubitBoolean align_direction,
                                             CubitBoolean perpendicular,
                                             CubitBoolean simplify_option) const;
    
  virtual CubitStatus create_surface( DLIList<CubitVector*>& vec_list,
                                      BodySM *&new_body,
                                      Surface *ref_face_ptr, 
			                             CubitBoolean project_points ) const;

  virtual CubitStatus create_weld_surface( CubitVector &root,
                                           Surface *ref_face1, double leg1, Surface *ref_face2, double leg2,
                                           BodySM *&new_body ) const;

  //--------------------------------------------------------------
   //- Methods for building specific facet-based geometry entities
   //--------------------------------------------------------------
  CubitStatus make_facet_point( CubitPoint *thePoint,
                                Point *&new_point_ptr );
  CubitStatus make_facet_point( CubitVector &location,
                                Point *&new_point_ptr );
    //- create a new facet point

  CubitStatus make_facet_curve( Point *start_ptr,
                                Point *end_ptr,
                                Curve *&new_curve_ptr,
                                CurveFacetEvalTool *eval_tool_ptr = NULL);
  CubitStatus make_facet_curve( Point *start_ptr,
                                Point *end_ptr,
                                DLIList<CubitFacetEdge*> &edge_list,
                                DLIList<CubitPoint*> &point_list,
                                Curve *&new_curve_ptr,
                                CurveFacetEvalTool *eval_tool_ptr = NULL);
   //- create a new facet curve 

  CubitStatus make_facet_coedge( Curve *curv_ptr,
                                 CubitSense sense,
                                 CoEdgeSM *&new_coedge_ptr );
   //- create a new facet coedge

  CubitStatus make_facet_loop( DLIList<CoEdgeSM*> &coedge_list,
                               LoopSM *&new_loop_ptr );
   //- create a new facet loop

  CubitStatus make_facet_surface(const CubitEvaluatorData *eval_data,
                                 DLIList<CubitFacet*> &facet_list,
                                 DLIList<CubitPoint*> &point_list,
                                 DLIList<LoopSM*> &my_loops,
                                 int interp_order,
                                 double min_dot,
                                 Surface *&new_surface_ptr,
                                 CubitBoolean use_point_addresses = CUBIT_TRUE,
                                 FacetEvalTool *eval_tool_ptr = NULL);
    //-creates a new FacetSurface given the points and facet list.

  CubitStatus make_facet_surface(DLIList<CubitQuadFacet*> &facet_list,
                                 DLIList<CubitPoint*> &point_list,
                                 DLIList<LoopSM*> &my_loops,
                                 int interp_order,
                                 double min_dot,
                                 Surface *&new_surface_ptr);
    //-creates a new FacetSurface given the points and quad facet list.

  CubitStatus make_facet_shell(DLIList<Surface*> &surface_list,
                               ShellSM *&new_shell_ptr);
    //-creates a new shell, given the list of surfaces.

  CubitStatus make_facet_lump(DLIList<ShellSM*> &shell_list,
                              Lump*& new_lump_ptr);
    //-creates a new lump, given the list of shells.

  CubitStatus make_facet_body(DLIList<Lump*> &lump_list,
                              BodySM *&new_body_ptr);
    //-creates a new body, given the list of lump.

  CubitStatus build_facet_surface( const CubitEvaluatorData **eval_data,
                                   DLIList<CubitFacet *> &facet_list,
                                   DLIList<CubitPoint *> &point_list,
                                   double feature_angle,
                                   int interp_order,
                                   CubitBoolean smooth_non_manifold,
                                   CubitBoolean split_surfaces,
                                   DLIList<Surface *> &surface_list);
  CubitStatus build_facet_surface( DLIList<CubitQuadFacet *> &qfacet_list,
                                   DLIList<CubitFacet *> &tfacet_list,
                                   DLIList<CubitPoint *> &point_list,
                                   double feature_angle,
                                   int interp_order,
                                   CubitBoolean smooth_non_manifold,
                                   CubitBoolean split_surfaces,
                                   DLIList<Surface *> &surface_list);
  CubitStatus build_facet_surface( DLIList<CubitQuadFacet *> &facet_list,
                                   DLIList<CubitPoint *> &point_list,
                                   double feature_angle,
                                   int interp_order,
                                   CubitBoolean smooth_non_manifold,
                                   CubitBoolean split_surfaces,
                                   DLIList<Surface *> &surface_list);
    //- creates one or more new FacetSurfaces with all of its lower
    //- order entities.  Use an optional feature_angle to break
    //- surface

  CubitStatus smooth_facets( RefFace *ref_face_ptr, int niter, CubitBoolean free_laplacian );
    //- attempt to clean up facets by smoothing the points on the surface
  CubitStatus create_shell_offset( BodySM *bodysm_ptr, BodySM *&new_bodysm, double offset );
    // create a shell offset from body
  CubitStatus improve_facets( RefFace *ref_face_ptr );
    // improve the facets by local swaps

  CubitStatus build_cholla_surfaces( DLIList<CubitFacet *> facet_list,
                                     DLIList<CubitPoint *> point_list,
                                     double feature_angle,
                                     int interp_order,
                                     CubitBoolean smooth_non_manifold,
                                     CubitBoolean split_surfaces,
                                     ChollaEngine *&cholla_ptr );

/*  virtual CubitStatus finish_facet_Body(
    GeometryType surface_type,
    const CubitEvaluatorData *eval_data,
    DLIList <CubitFacet *>facet_list,
    DLIList <CubitPoint *>point_list,
    double feature_angle,
    int interp_order,
    CubitBoolean smooth_non_manifold, 
    CubitBoolean split_surfaces,
    BodySM *&body_ptr) const;
*/
  virtual CubitStatus finish_facet_Body( ChollaEngine *&cholla_ptr,
                                         const CubitEvaluatorData **eval_data,
                                         double feature_angle,
                                         int interp_order,
                                         BodySM *&bodysm_ptr) const;

  void set_sphere_eval_data( ChollaEngine *cholla_ptr,
                             double radius,
                             CubitEvaluatorData **&eval_data ) const;

  void set_cylinder_eval_data( ChollaEngine *cholla_ptr,
                               double height,
                               double base_radius_xdir,
                               double base_radius_ydir,
                               double top_radius,
                               CubitEvaluatorData **&eval_data ) const;

  // non-virtual specific functions for building facet-based geometry
  // from Cholla geometry

  CubitStatus build_cholla_geometry(
    const CubitEvaluatorData **eval_data,
    DLIList<ChollaSurface*> &cholla_surface_list,
    DLIList<ChollaCurve*> &cholla_curve_list,
    DLIList<ChollaPoint*> &cholla_point_list,
    CubitBoolean use_feature_angle, 
    double feature_angle, 
    int interp_order,
    DLIList<Surface *> &surface_list);
    // build the CUBIT geometry based on the Cholla entity class lists

  CubitStatus build_cholla_point_geometry(
    DLIList<ChollaPoint*> &cholla_point_list );
    // From the cholla point list, create geometric points for each

  CubitStatus build_cholla_curve_geometry(
    DLIList<ChollaCurve*> &cholla_curve_list );
    // From the cholla curve list, create geometric curves for each

  CubitStatus build_cholla_surface_geometry(
    const CubitEvaluatorData **eval_data,
    DLIList<ChollaSurface*> &cholla_surface_list,
    int interp_order,
    double min_dot,
    DLIList<Surface *> &surface_list);
    // From the facet surface list, create geometric surface,
    // loops and coedges for each surface in the list  

    CubitStatus tolerant_imprint( DLIList<BodySM*> &bodies_in,
                                  DLIList<BodySM*> &new_bodies,
                                   DLIList<TopologyBridge*> *new_tbs = NULL,
                                   DLIList<TopologyBridge*> *att_tbs = NULL ) const;

  static CubitBoolean is_modify_enabled()
    {return modifyEnabled;}

  static void set_modify_enabled(CubitBoolean my_bool)
    {modifyEnabled = my_bool;}
  
      
  
  
protected:
     
private:
  
  static CubitBoolean modifyEnabled;
  
  CubitStatus build_cholla_loop_geometry(
    DLIList<ChollaCurve*> &cholla_curve_list,
    ChollaSurface *chsurf_ptr,
    DLIList<LoopSM*> &loop_list,
    int &ncurves,
    int debug_draw = 0 );
    // From the cholla curve list of a surface, create geometric loops 
 
  void fillinedge( 
    int *edge, 
    int numpointsonanedge, 
    double radius, 
    std::vector<CubitPoint *>& points) const;
    //!  Put points on this edge of a triangle being refined.
    
  void refinetriangle(
    int level, 
    int numpointsonanedge, 
    int *iedge1, 
    int *iedge2, 
    int *iedge3,
    int isign1, 
    int isign2, 
    int isign3, 
    double radius, 
    std::vector<CubitPoint *>& points,
    DLIList<CubitFacet *>& facet_list) const;    
    //! add internal points and make connections for this triangle


} ;

#endif
