#ifndef ACIS_MODIFY_ENGINE_HPP
#define ACIS_MODIFY_ENGINE_HPP

#include "GeometryModifyEngine.hpp"
#include "NoAcisQueryEngine.hpp"

#include <set>
class TopologyBridge;
class GeometryEntity;
class BodySM;
class Lump;
class ShellSM;
class Surface;
class LoopSM;
class Curve;
class Point;
class BodySM;

class CubitBox;
class CubitString;
class GMem;
template <class X> class DLIList;

class AcisModifyEngine : public GeometryModifyEngine
{
public:
  friend class AcisQueryEngine;
private:

  AcisModifyEngine();

  static AcisModifyEngine* instance_;
  
public:

  static inline AcisModifyEngine* instance()
    {
      if( !instance_ )
        instance_ = new AcisModifyEngine;
      return instance_;
    }

  virtual ~AcisModifyEngine();

  CubitString identify();
    //- Return a string identifying this version of ACIS.

  virtual const type_info& entity_type_info() const ;
  virtual BodySM* sphere(double radius) const ;
  virtual BodySM* brick( double wid, double dep, double hi ) const ;
  virtual BodySM* brick( const CubitVector &center, const CubitVector axes[3], const CubitVector &extension) const ;
  virtual BodySM* prism( double height, int sides, double major, double minor) const ;
  virtual BodySM* pyramid( double height, int sides, double major, double minor, double top=0.0) const ;
  virtual BodySM* cylinder( double hi, double r1, double r2, double r3 ) const ;
  virtual BodySM* torus( double r1, double r2 ) const ;
  virtual BodySM* planar_sheet( const CubitVector& p1, const CubitVector& p2, const CubitVector& p3, const CubitVector& p4 ) const;
  virtual CubitStatus subtract( DLIList<BodySM*> &tool_body_ptr, DLIList<BodySM*> &from_bodies, DLIList<BodySM*> &new_bodies, bool keep_old = false) const ;
  virtual CubitStatus imprint( BodySM* BodyPtr1, BodySM* BodyPtr2, BodySM*& newBody1,  BodySM*& newBody2, bool keep_old = false) const ;
  virtual CubitStatus imprint(DLIList<BodySM*> &from_body_list, DLIList<BodySM*> &new_from_body_list, bool keep_old = false) const;
  virtual CubitStatus imprint( DLIList<BodySM*> &body_list, DLIList<Curve*> &ref_edge_list, DLIList<BodySM*>& new_body_list, bool keep_old_body, bool show_messages = true) const;
  virtual CubitStatus imprint( DLIList<Surface*> &surface_list, DLIList<Curve*> &curve_list, DLIList<BodySM*>& new_body_list, bool keep_old_body ) const;
  virtual CubitStatus imprint( DLIList<Surface*> &surface_list, DLIList<DLIList<Curve*>*> &curve_lists_list, BodySM*& new_body, bool keep_old_body ) const;
  virtual CubitStatus imprint( DLIList<BodySM*> &body_list, DLIList<CubitVector*> &vector_list, DLIList<BodySM*>& new_body_list, bool keep_old_body ) const;
  virtual CubitStatus imprint_projected_edges( DLIList<Surface*> &ref_face_list, DLIList<Curve*> &ref_edge_list, DLIList<BodySM*>& new_body_list, bool keep_old_body, bool keep_free_edges ) const;
  virtual CubitStatus imprint_projected_edges(DLIList<Surface*> &ref_face_list, DLIList<BodySM*> &body_list, DLIList<Curve*> &ref_edge_list, DLIList<BodySM*>& new_body_list, bool keep_old_body, bool keep_free_edges ) const;
  virtual CubitStatus tolerant_imprint( DLIList<BodySM*> &bodies_in, DLIList<BodySM*> &new_bodies ) const;
  virtual CubitStatus project_edges( DLIList<Surface*> &ref_face_list, DLIList<Curve*> &ref_edge_list_in, DLIList<Curve*> &ref_edge_list_new, bool ) const;
  virtual CubitBoolean bodies_overlap( BodySM *body_ptr_1, BodySM *body_ptr_2 ) const;
  virtual CubitStatus intersect( BodySM* tool_body_ptr, DLIList<BodySM*> &from_bodies, DLIList<BodySM*> &new_bodies, bool keep_old = false) const ;
  virtual CubitStatus unite( DLIList<BodySM*> &bodies, DLIList<BodySM*> &newBody, bool keep_old = false) const ;
  virtual CubitStatus chop( DLIList<BodySM*> &bodies, DLIList<BodySM*> &intersectBodies, DLIList<BodySM*> &outsideBodies, BodySM*& leftoversBody, bool keep_old = false, bool nonreg = false) const ;
  virtual CubitStatus thicken( DLIList<BodySM*>& bodies, DLIList<BodySM*>& new_bodies, double depth, bool both = false) const;
  virtual CubitStatus flip_normals( DLIList<Surface*>& face_list ) const;
  virtual CubitStatus  sweep_translational(DLIList<GeometryEntity*>& ref_ent_list, DLIList<BodySM*>& result_body_list, const CubitVector& sweep_vector, double draft_angle, int draft_type, bool switchside, bool rigid) const ;
  virtual CubitStatus  sweep_perpendicular(DLIList<GeometryEntity*>& ref_ent_list, DLIList<BodySM*>& result_body_list, double distance, double draft_angle, int draft_type, bool switchside, bool rigid) const ;
  virtual CubitStatus  sweep_rotational(DLIList<GeometryEntity*>& ref_ent_list, DLIList<BodySM*>& result_body_list, const CubitVector& point, const CubitVector& direction, double angle, int steps = 0, double draft_angle = 0.0, int draft_type = 0, bool switchside = false, bool make_solid = false, bool rigid = false ) const ;
  virtual CubitStatus sweep_along_curve(DLIList<GeometryEntity*>& ref_ent_list, DLIList<BodySM*>& result_body_list, DLIList<Curve*>& ref_edge_list, double draft_angle = 0.0, int draft_type = 0, bool rigid = false ) const; 
  virtual CubitStatus webcut(DLIList<BodySM*>& webcut_body_list, const CubitVector &v1, const CubitVector &v2, const CubitVector &v3, DLIList<BodySM*>& results_list, bool imprint = false ) const ;
  virtual CubitStatus webcut( DLIList<BodySM*>& webcut_body_list, BodySM const* tool_body, DLIList<BodySM*>& results_list, bool imprint  = false ) const ;
  virtual CubitStatus webcut_with_cylinder(DLIList<BodySM*> &webcut_body_list, double radius, const CubitVector &axis, const CubitVector &center, DLIList<BodySM*>& results_list, bool imprint = false );
  virtual CubitStatus webcut_with_brick(DLIList<BodySM*>& webcut_body_list, const CubitVector &center, const CubitVector axes[3], const CubitVector &extension, DLIList<BodySM*> &results_list, bool imprint = false );
  virtual CubitStatus webcut_with_planar_sheet(DLIList<BodySM*>& webcut_body_list, const CubitVector &center, const CubitVector axes[2], double width, double height, DLIList<BodySM*> &results_list, bool imprint = false );
  virtual CubitStatus webcut_with_curve_loop(DLIList<BodySM*> &webcut_body_list, DLIList<Curve*> &ref_edge_list, DLIList<BodySM*>& results_list, bool imprint = false);
  virtual CubitStatus webcut_with_sweep_surfaces_rotated(DLIList<BodySM*> &blank_bodies, DLIList<Surface*> &surfaces, const CubitVector& point, const CubitVector& sweep_vector, double angle, Surface *stop_surf, bool up_to_next, DLIList<BodySM*> &results_list, CubitBoolean imprint = false );
  virtual CubitStatus webcut_with_sweep_curves_rotated(DLIList<BodySM*> &blank_bodies, DLIList<Curve*> &curves, const CubitVector& point, const CubitVector& sweep_vector, double angle, Surface *stop_surf, DLIList<BodySM*> &results_list, CubitBoolean imprint = false );
  virtual CubitStatus webcut_with_sweep_surfaces(DLIList<BodySM*> &blank_bodies, DLIList<Surface*> &surfaces, const CubitVector& sweep_vector, bool sweep_perp, bool through_all, bool outward, bool up_to_next,  Surface *stop_surf, Curve *curve_to_sweep_along, DLIList<BodySM*> &results_list, CubitBoolean imprint = false );
  virtual CubitStatus webcut_with_sweep_curves(DLIList<BodySM*> &blank_bodies, DLIList<Curve*> &curves, const CubitVector& sweep_vector, bool through_all, Surface *stop_surf, Curve *curve_to_sweep_along, DLIList<BodySM*> &results_list, CubitBoolean imprint = false );
  virtual CubitStatus section( DLIList<BodySM*> &section_body_list, const CubitVector &point_1, const CubitVector &point_2, const CubitVector &point_3, DLIList<BodySM*>& new_body_list, bool keep_normal_side, bool keep_old = false, bool keep_both_sides = false);
  virtual BodySM* copy_body ( BodySM* bodyPtr) const ;
  virtual Point* make_Point( CubitVector const& point) const ;
  virtual Curve* make_Curve(Curve *curve_ptr) const;
  virtual Curve* make_Curve( Point const* point1_ptr, Point const* point2_ptr, Surface* ref_face_ptr, const CubitVector *third_point = NULL ) const;
  virtual Curve* make_Curve( GeometryType curve_type, Point const* point1_ptr, Point const* point2_ptr, DLIList<CubitVector*>& vector_list, Surface* ref_face_ptr = NULL ) const;
  virtual Curve* make_Curve( GeometryType curve_type, Point const* point1_ptr, Point const* point2_ptr, CubitVector const* intermediate_point_ptr, CubitSense sense) const;
  virtual Surface* make_Surface( Surface *old_surface_ptr, bool extended_from ) const;
  virtual Surface* make_Surface( GeometryType surface_type, DLIList<Curve*>& curve_list, Surface *old_surface_ptr = NULL, bool check_edges = true ) const; 
  virtual Lump* make_Lump( DLIList<Surface*>& surface_list ) const ;
  virtual BodySM* make_BodySM( Surface * ) const;
  virtual BodySM* make_BodySM( DLIList<Lump*>& lump_list ) const ;
  virtual CubitStatus webcut_with_sheet(DLIList<BodySM*>& webcut_body_list, BodySM *sheet_body, DLIList<BodySM*> &new_bodies, bool imprint = false);
  virtual CubitStatus webcut_with_extended_surf(DLIList<BodySM*> &webcut_body, Surface *extend_from, DLIList<BodySM*> &new_bodies, int &num_cut, bool imprint = false);
  virtual CubitStatus split_body( BodySM *body_ptr, DLIList<BodySM*> &new_bodies );
  virtual CubitStatus split_periodic( BodySM *body_ptr, BodySM *&new_body);
  virtual CubitStatus reverse_body( BodySM *body_to_reverse );
  virtual CubitStatus create_body_from_surfs(DLIList<Surface*> &ref_face_list, BodySM *&new_body, bool keep_old = false, bool heal = true) const;
  CubitStatus combine_surfaces( DLIList<Surface*>& face_list, BodySM*& new_body, bool keep_old = false) const;
  CubitStatus combine_sheets( DLIList<BodySM*>& body_list, BodySM*& new_body, bool keep_old = false) const ;
  virtual CubitStatus webcut_across_translate( DLIList<BodySM*>& body_list, Surface* plane_surf1, Surface* plane_surf2, DLIList<BodySM*>& results_list, bool imprint = false) const;
  CubitStatus offset_curves( DLIList<Curve*>& ref_edge_list, DLIList<Curve*>& new_curves, double offset_distance, const CubitVector& offset_direction, int gap_type = 1 );
  Curve* trim_curve( Curve* trim_curve, const CubitVector& trim_vector, const CubitVector& keep_vector, bool keep_old = false );
  Curve* create_arc_three( Point* ref_vertex1, Point* ref_vertex2, Point* ref_vertex3, bool full = false );
  Curve* create_arc_three( Curve* ref_edge1, Curve* ref_edge2, Curve* ref_edge3, bool full = false );
  Curve* create_arc_center_edge( Point* ref_vertex1, Point* ref_vertex2, Point* ref_vertex3, const CubitVector &normal, double radius = CUBIT_DBL_MAX, bool full = false );
  CubitStatus create_curve_combine( DLIList<Curve*>& curve_list, Curve *&new_curve_ptr );
  GeometryQueryEngine *get_gqe();
  virtual CubitBoolean is_modify_engine(const TopologyBridge *tb_ptr) const;
  virtual CubitStatus get_offset_intersections(Curve* ref_edge1, Curve* ref_edge2, DLIList<CubitVector*>& intersection_list, double offset, bool ext_first = true );
  virtual CubitStatus get_offset_intersections(Curve* ref_edge_ptr, Surface* ref_face_ptr, DLIList<CubitVector*> &intersection_list, double offset = 0.0, bool ext_surf = true );
  virtual CubitStatus surface_intersection( Surface *surface1_ptr, Surface *surface2_ptr, DLIList<Curve*> &inter_graph, const double tol) const;
  virtual CubitStatus get_mid_plane( const CubitVector &point_1, const CubitVector &point_2, const CubitVector &point_3, BodySM *body_to_trim_to, BodySM *&midplane_body ) const;
  virtual CubitStatus tweak_chamfer( DLIList<Curve*> &curve_list, double left_offset, DLIList<BodySM*> &new_bodysm_list, double right_offset = -1.0, CubitBoolean keep_old_body = CUBIT_FALSE, CubitBoolean preview = CUBIT_FALSE ) const;
  virtual CubitStatus tweak_chamfer( DLIList<Point*> &point_list, double offset1, DLIList<BodySM*> &new_bodysm_list, Curve *edge1 = NULL, double offset2 = -1.0, Curve *edge2 = NULL, double offset3 = -1.0, Curve *edge3 = NULL, CubitBoolean keep_old_body = CUBIT_FALSE, CubitBoolean preview = CUBIT_FALSE ) const;
  virtual CubitStatus tweak_fillet( DLIList<Curve*> &curve_list, double radius, DLIList<BodySM*> &new_bodysm_list, CubitBoolean keep_old_body = CUBIT_FALSE, CubitBoolean preview = CUBIT_FALSE ) const;
  virtual CubitStatus tweak_fillet( Curve *curve_ptr, double start_radius, double end_radius, BodySM *&new_body_ptr, CubitBoolean keep_old_body = CUBIT_FALSE, CubitBoolean preview = CUBIT_FALSE ) const;
  virtual CubitStatus tweak_fillet( DLIList<Point*> &point_list, double radius, DLIList<BodySM*> &new_bodysm_list, CubitBoolean keep_old_body = CUBIT_FALSE, CubitBoolean preview = CUBIT_FALSE ) const;
  virtual CubitStatus tweak_move(DLIList<Surface*>&, const CubitVector&, DLIList<BodySM*>&, CubitBoolean, CubitBoolean) const;
  virtual CubitStatus tweak_move(DLIList<Curve*>&, const CubitVector&, DLIList<BodySM*>&, CubitBoolean, CubitBoolean) const;

  virtual CubitStatus tweak_offset(DLIList<Surface*>&, double, DLIList<BodySM*>&, CubitBoolean, CubitBoolean) const;
  virtual CubitStatus tweak_offset(DLIList<Curve*>&, double, DLIList<BodySM*>&, CubitBoolean, CubitBoolean) const;
  virtual CubitStatus tweak_remove(DLIList<Surface*>&, DLIList<BodySM*>&, CubitBoolean, CubitBoolean, CubitBoolean, CubitBoolean) const;
  virtual CubitStatus tweak_remove(DLIList<Curve*>&, DLIList<BodySM*>&, CubitBoolean, CubitBoolean) const;
  virtual CubitStatus tweak_target(DLIList<Surface*>&, DLIList<Surface*>&, DLIList<BodySM*>&, CubitBoolean, CubitBoolean, CubitBoolean) const;
  virtual CubitStatus tweak_target(DLIList<Curve*>&, DLIList<Surface*>&, DLIList<BodySM*>&, CubitBoolean, CubitBoolean, CubitBoolean) const;
  virtual CubitStatus tweak_target(DLIList<Curve*>&, DLIList<Curve*>&, DLIList<BodySM*>&, CubitBoolean, CubitBoolean, CubitBoolean) const;

  virtual CubitStatus create_net_surface( DLIList<Surface*>& ref_face_list, BodySM *& new_body, DLIList<DLIList<CubitVector*>*> &vec_lists_u,  DLIList<DLIList<CubitVector*>*> &vec_lists_v, double net_tol = 1e-3, CubitBoolean heal = CUBIT_TRUE ) const;
  virtual CubitStatus create_net_surface( DLIList<Curve*>& u_curves, DLIList<Curve*>& v_curves, BodySM *& new_body, double net_tol = 1e-3, CubitBoolean heal = CUBIT_TRUE ) const;
  virtual CubitStatus create_offset_surface( Surface* ref_face_ptr, BodySM*& new_body, double offset_distance ) const;
  virtual CubitStatus create_offset_body( BodySM* body_ptr, BodySM*& new_body, double offset_distance ) const;
  virtual CubitStatus create_skin_surface( DLIList<Curve*>& curves, BodySM*& new_body ) const;
  virtual CubitStatus loft_surfaces( Surface *face1, const double &takeoff1, Surface *face2, const double &takeoff2, BodySM*& new_body, CubitBoolean arc_length_option = CUBIT_FALSE, CubitBoolean twist_option = CUBIT_FALSE, CubitBoolean align_direction = CUBIT_TRUE, CubitBoolean perpendicular = CUBIT_TRUE, CubitBoolean simplify_option = CUBIT_FALSE) const;
  virtual CubitStatus loft_surfaces_to_body( Surface *face1, const double &takeoff1, Surface *face2, const double &takeoff2, BodySM*& new_body, CubitBoolean arc_length_option, CubitBoolean twist_option, CubitBoolean align_direction, CubitBoolean perpendicular, CubitBoolean simplify_option) const;
  virtual CubitStatus create_surface( DLIList<CubitVector*>& vec_list, BodySM *&new_body, Surface *ref_face_ptr, CubitBoolean project_points ) const;
  virtual CubitStatus create_weld_surface( CubitVector &root, Surface *ref_face1, double leg1, Surface *ref_face2, double leg2, BodySM *&new_body ) const;
  CubitStatus scale( BodySM *&body, const CubitVector& f );

  virtual CubitStatus get_spheric_mid_surface(Surface*, Surface*, BodySM*, BodySM*&) const;
  virtual CubitStatus get_conic_mid_surface(Surface*, Surface*, BodySM*, BodySM*&) const;
  virtual CubitStatus get_toric_mid_surface(Surface*, Surface*, BodySM*, BodySM*&) const;

protected:

private:
  CubitStatus regularize_body( BodySM *body_ptr, BodySM *&new_body_ptr );
  CubitStatus regularize_entity( GeometryEntity *old_refentity_ptr, BodySM *&new_body_ptr );
};

#endif


