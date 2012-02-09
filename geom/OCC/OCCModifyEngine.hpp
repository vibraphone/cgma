//-------------------------------------------------------------------------
// Filename      : OCCModifyEngine.hpp
//
// Purpose       : ModifyEngine for OCC geometry
//
// Special Notes : Modeled after GeometryModifyEngine and AcisModifyEngine.
//
// Creator       : Jane Hu
//
// Creation Date : 6/08
//
//-------------------------------------------------------------------------

#ifndef OCC_MODIFY_ENGINE_HPP
#define OCC_MODIFY_ENGINE_HPP

#include "GeometryModifyEngine.hpp"
#include <vector>
#include <map>
#include "DLIList.hpp"

class Point;
class TopologyBridge;
class CoEdgeSM;
class ShellSM;
class OCCBody;
class OCCSurface;
class OCCCurve;
class OCCPoint;
class TopoDS_Shape;
class TopoDS_Edge;
class TopoDS_Wire;
class TopoDS_Face;
class CubitBox;
class BRepAlgoAPI_BooleanOperation;
class BRepOffsetAPI_ThruSections;
class BRepBuilderAPI_Copy;
class OCCHistory;

class OCCModifyEngine : public GeometryModifyEngine
{
  
public:
  
  //HEADER- Constructor and Destructor
private:
  
  OCCModifyEngine();
  
  static OCCModifyEngine* instance_;

  CubitStatus copy_body_attributes(TopoDS_Shape orig_shape,
                            BRepBuilderAPI_Copy& api_copy)const;
  
public:
  static inline OCCModifyEngine* instance()
  {
    if( !instance_ )
      instance_ = new OCCModifyEngine;
    return instance_;
  }
  
  virtual ~OCCModifyEngine();
  //- virtual destructor
  
  virtual Point* make_Point( CubitVector const& point) const ;
  
  virtual Curve* make_Curve( DLIList<CubitVector*>& point_list,
                             DLIList<CubitVector*>& point_tangents) const;

  virtual Curve* make_Curve(Curve *curve_ptr) const;
  //- creates a curve from an existing curve.  This creates totally
  //- new topology.  This function is useful for constructing geometry
  //- from existing geometry.
  
  virtual Curve* make_Curve( Point const* point1_ptr,
    Point const* point2_ptr,
    Surface* ref_face_ptr = NULL,
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
  
  virtual BodySM* make_extended_sheet( DLIList<Surface*> &surface_list,
                                       CubitBox *clip_box = NULL,
                                       bool preview = false) const;
    //R BodySM*
    //R- Pointer to a newly created BodySM object.
    //I surface_list
    //I- The surface_list from which we want to create an extended sheet.
    //I clip_box
    //I- An optional bounding box to clip the resultant sheet body by.
    //I preview
    //I- If true just draw the sheet instead of creating it
    //- This function creates a sheet body by extending the input surfaces.
    //- The result can be optionally clipped to fit inside of the given
    //- bounding box.
  
  virtual Surface* make_Surface( Surface *old_surface_ptr,
    CubitBoolean extended_from = CUBIT_FALSE) const;
  
  virtual Surface* make_Surface( GeometryType surface_type,
    DLIList<Curve*>& curve_list,
    Surface *old_surface_ptr = NULL,
    bool check_edges = true ) const;
  
  CubitStatus create_rectangle_surface( double width,
                                        double height,
                                        CubitVector plane,
                                        BodySM *&sheet_body) const;

  CubitStatus create_ellipse_surface( Point  *pt1,
                                      Point  *pt3,
                                      CubitVector center_point,
                                      BodySM *&sheet_body) const;

  CubitStatus create_ellipse_surface( double major_radius,
                                      double minor_radius,
                                      CubitVector plane,
                                      BodySM *&sheet_body) const;

  Curve* make_elliptical_Curve( Point const* point1,
                                Point const* point2,
                                CubitVector &center_point,
                                double start_angle,
                                double end_angle,
                                CubitSense sense) const;

  CubitStatus create_circle_surface( Point *pt1,
                                     CubitVector center_point,
                                     Point *pt3,
                                     BodySM *&sheet_body) const;

  CubitStatus create_circle_surface( Point *pt1,
                                     Point *pt3,
                                     CubitVector center_point,
                                     BodySM *&sheet_body) const;

  CubitStatus create_circle_surface( double radius,
                                     CubitVector plane,
                                     BodySM *&sheet_body) const;

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
    bool keep_old,
    DLIList<TopologyBridge*> *new_tbs = NULL,
    DLIList<TopologyBridge*> *att_tbs = NULL) const; 
  
  virtual CubitStatus imprint( DLIList<BodySM*> &body_list,
                               DLIList<Curve*> &ref_edge_list,
                               DLIList<BodySM*>& new_body_list,
                               DLIList<TopologyBridge*> &temporary_bridges,
                               bool keep_old_body,
                               bool show_messages = true) const;
    //- Imprints a list of Bodies with a list of RefEdges.  All entities must
    //- be ACIS entities.  Useful for splitting surfaces.  If edge pierces a
    //- surface a hardpoint will result at the pierce location.
 
  virtual CubitStatus imprint( DLIList<Surface*> &surface_list,
                               DLIList<Curve*> &curve_list,
                               DLIList<TopologyBridge*> &temporary_bridges,
                               DLIList<BodySM*>& new_body_list,
                               bool keep_old_body ) const;
    //- Imprints a list of RefFaces with a list of RefEdges.  This is
    //- useful if the user has a curve which spans several surfaces on
    //- a body and only wants to imprint to selected surfaces.  Algorithm
    //- does not support imprinting to free surfaces.

  virtual CubitStatus imprint( DLIList<Surface*> &surface_list,
                               DLIList<DLIList<Curve*>*> &curve_lists_list,
                               BodySM*& new_body,
                               bool keep_old_body,
                               bool expand = true,
                               DLIList<TopologyBridge*> *new_tbs = NULL,
                               DLIList<TopologyBridge*> *att_tbs = NULL ) const;
    //- Imprints a list of Surfaces with list of Curves, sorted per
    //- Surface (ie., curve_lists_list is same length as surface_list).
  
  virtual CubitStatus imprint( DLIList<BodySM*> &body_list,
                               DLIList<CubitVector*> &vector_list,
                               DLIList<BodySM*>& new_body_list,
                               bool keep_old_body,
                               DLIList<TopologyBridge*> *new_tbs = NULL,
                               DLIList<TopologyBridge*> *att_tbs = NULL,
                               double *tol_in = NULL,
                               bool clean_up_slivers = true) const;
    //- Imprints a list of bodies with a list of vectors.  Useful for
    //- splitting curves and creating hardpoints on surfaces.

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

  virtual CubitStatus hollow( DLIList<BodySM*>& bodies,
                              DLIList<Surface*>& surfs_to_remove,
                              DLIList<BodySM*>& new_bodies,
                              double depth) const ;

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
    BodySM  *stop_body = NULL) const;
  
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
    bool rigid = CUBIT_FALSE,
    Surface *stop_surf = NULL,
    BodySM  *stop_body = NULL ) const;
  
  virtual CubitStatus sweep_along_curve(
    DLIList<GeometryEntity*>& ref_ent_list,
    DLIList<BodySM*>& result_body_list,
    DLIList<Curve*>& ref_edge_list,
    double draft_angle = 0.0,
    int draft_type = 0,
    bool rigid = CUBIT_FALSE,
    Surface *stop_surf = NULL,
    BodySM  *stop_body = NULL ) const;
 

  //HEADER- Webcut-related functions
  virtual CubitStatus webcut(
                      DLIList<BodySM*>& webcut_body_list,
                      const CubitVector &v1,
                      const CubitVector &v2,
                      const CubitVector &v3,
                      DLIList<BodySM*>& neighbor_imprint_list,
                      DLIList<BodySM*>& results_list,
                      ImprintType imprint_type = NO_IMPRINT,
                      bool preview = false) const ;
    //R int
    //R- Number of bodies that were webcut ( >= 0 )
    //I webcut_body_list
    //I- The list of bodies to be webcut
    //I plane
    //I- The plane to be used for webcutting.
    //I merge
    //I- A flag to decide whether the new bodies created by the
    //I- webcutting process should be merged or not.
    //I imprint
    //I- A flag to decide whether the new bodies created by the
    //I- webcutting process should be imprinted or not.
    //- This functions webcuts a list of bodies through a plane.
    //- The newly created bodies are merged and imprinted depeding on
    //- the respective flags.
  
  virtual CubitStatus webcut( DLIList<BodySM*>& webcut_body_list,
                              BodySM const* tool_body,
                              DLIList<BodySM*>& neighbor_imprint_list,
                              DLIList<BodySM*>& results_list,
                              ImprintType imprint_type = NO_IMPRINT,
                              bool preview = false) const ;
    //R int
    //R- Number of bodies that were webcut ( >= 0 )
    //I webcut_body_list
    //I- The list of bodies to be webcut
    //I tool_body
    //I- The body to be used for webcutting.
    //I merge
    //I- A flag to decide whether the new bodies created by the
    //I- webcutting process should be merged or not.
    //I imprint
    //I- A flag to decide whether the new bodies created by the
    //I- webcutting process should be imprinted or not.
    //- This functions webcuts a list of bodies using another body
    //- as the webcutting tool. The newly created bodies are
    //- merged and imprinted depeding on the respective flags.
  
  virtual CubitStatus webcut_across_translate( DLIList<BodySM*>& body_list,
                                               Surface* plane_surf1,
                                               Surface* plane_surf2,
                                               DLIList<BodySM*>& neighbor_imprint_list,
                                               DLIList<BodySM*>& results_list,
                                               ImprintType imprint_type = NO_IMPRINT,
                                               bool preview = false) const;
  // In-process function to webcut a flat plate suitable for singe-single sweeping.
 
  virtual CubitStatus separate_surfaces( DLIList<Surface*> &surf_list,
                                         DLIList<BodySM*> &new_bodies );
      //- Separates surfaces from sheet bodies into separate bodies.  Connected
      //- surfaces will remain connected but be placed in a new body.

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
  
  virtual CubitStatus create_solid_bodies_from_surfs(
                                    DLIList<Surface*> &ref_face_list,
                                    DLIList<BodySM*> &new_bodies,
                                    bool keep_old = false,
                                    bool heal = true,
                                    bool sheet = false) const;
  //- This function assumes that the reffaces sent into
  //- this function are either sheet bodies, or free surfaces.  This
  //- Will have been taken care of in the calling function.  GT?
 
  CubitStatus tweak_bend( DLIList<BodySM*> &bend_bodies,
                          DLIList<BodySM*> &new_bodysm_list,
                          CubitVector& neutral_root,
                          CubitVector& bend_axis,
                          CubitVector& bend_direction,
                          double radius,
                          double angle,
                          DLIList<CubitVector*> bend_regions,
                          double width = -1,
                          CubitBoolean center_bend = CUBIT_FALSE,
                          int num_points = 0,
                          CubitBoolean keep_old_body = CUBIT_FALSE,
                          CubitBoolean preview = CUBIT_FALSE ) const
  /**<  Bend solid bodies based on a bend radius and angle.
    */
  { return CUBIT_FAILURE;}

  virtual Curve* create_arc_three( Point* ref_vertex1, 
                                   Point* ref_vertex2,
                                   Point* ref_vertex3, 
                                   bool full = false )const;
  
  virtual Curve* create_arc_three( Curve* ref_edge1, 
                                   Curve* ref_edge2,
                                   Curve* ref_edge3, 
                                   bool full = false )const;
  
  virtual Curve* create_arc_center_edge( 
                                   Point* ref_vertex1, 
                                   Point* ref_vertex2,
                                   Point* ref_vertex3, 
                                   const CubitVector &normal,
                                   double radius = CUBIT_DBL_MAX,
                                   bool full = false )const;
  
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
  
  CubitStatus get_3_point_plane( const CubitVector & point_1,
                                 const CubitVector & point_2,
                                 const CubitVector & point_3,
                                 TopoDS_Face*& face)const;

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
                                     CubitBoolean preview = CUBIT_FALSE ) ;
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
                                     CubitBoolean preview = CUBIT_FALSE ) ;
  /**<  Chamfer vertices on solid or sheet bodies.  On a solid body there can
    *   be up to 3 offsets; on a sheet body up to 2 offsets.  The offsets are
    *   in the direction of the supplied edges.  If multiple vertices are 
    *   supplied, only one offset value is allowed and the edges are not used.
    */

  virtual CubitStatus tweak_fillet( DLIList<Curve*> &curve_list, 
                                    double radius,
                                    DLIList<BodySM*> &new_bodysm_list,
                                    CubitBoolean keep_old_body = CUBIT_FALSE,
                                    CubitBoolean preview = CUBIT_FALSE ) ;
  /**<  Create a round fillet (or blend) at the given curves on solid bodies.
    */

  virtual CubitStatus tweak_fillet( Curve *curve_ptr, 
                                    double start_radius,
                                    double end_radius,
                                    BodySM *&new_body_ptr,
                                    CubitBoolean keep_old_body = CUBIT_FALSE,
                                    CubitBoolean preview = CUBIT_FALSE ) ;
  /**<  Create a round fillet (or blend) at the given curve on a solid body.
    *   The fillet has a variable radius from the start to the end of the curve.
    */

  virtual CubitStatus tweak_fillet( DLIList<Point*> &point_list, 
                                    double radius,
                                    DLIList<BodySM*> &new_bodysm_list,
                                    CubitBoolean keep_old_body = CUBIT_FALSE,
                                    CubitBoolean preview = CUBIT_FALSE ) ;
  /**<  Create a round fillet (or blend) at the given vertices on sheet bodies.
    */

  virtual CubitStatus tweak_move( DLIList<Surface*> &surface_list,
                                  const CubitVector &delta,
                                  DLIList<BodySM*> &new_bodysm_list,
                                  CubitBoolean keep_old_body = CUBIT_FALSE,
                                  CubitBoolean preview = CUBIT_FALSE  ) const;
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
                                    DLIList<Surface*> *add_surface_list_ptr,
                                    DLIList<double>*, 
                                    DLIList<BodySM*>&,
                                    CubitBoolean keep_old_body = CUBIT_FALSE,
                                    CubitBoolean preview = CUBIT_FALSE ) const;
  /**<  Tweak specified faces of a volume or volumes by offsetting those faces
    *   by the offset distance.
    */

  virtual CubitStatus tweak_offset( DLIList<Curve*> &curve_list,
                                    double offset_distance,
                                    DLIList<Curve*>*, 
                                    DLIList<double>*,
                                    DLIList<BodySM*> &new_bodysm_list,
                                    CubitBoolean keep_old_body = CUBIT_FALSE,
                                    CubitBoolean preview = CUBIT_FALSE ) const;
  /**<  Tweak specified curves of a sheet body or bodies by offsetting those
    *   curves by the offset distance.
    */

  virtual CubitStatus tweak_remove( DLIList<Surface*> &surface_list,
                                    DLIList<BodySM*> &new_bodysm_list,
                                    CubitBoolean extend_adjoining = CUBIT_TRUE,
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
                                    DLIList<Surface*> &target_surfs,
                                    DLIList<BodySM*> &new_bodysm_list,
                                    CubitBoolean extend_flg = CUBIT_TRUE,
                                    CubitPlane *limit_plane = NULL,
                                    CubitBoolean reverse_flg = CUBIT_FALSE,
                                    CubitBoolean keep_old_body = CUBIT_FALSE,
                                    CubitBoolean preview = CUBIT_FALSE ) const;
  /**<  Tweak specified faces of a volume or volumes up to a target surface.
    */

  virtual CubitStatus tweak_target( DLIList<Curve*> &curve_list,
                                    DLIList<Surface*> &target_surf_list,
                                    DLIList<BodySM*> &new_bodysm_list, 
                                    CubitBoolean extend_flg = CUBIT_TRUE,
                                    CubitPlane *limit_plane = NULL,
                                    CubitBoolean reverse_flg = CUBIT_FALSE,
                                    CubitBoolean keep_old_body = CUBIT_FALSE,
                                    CubitBoolean preview = CUBIT_FALSE,
                                    double max_area_increase = 0 ) const;
  /**<  Tweak specified edges of a surface or set of surfaces (in sheet
    *   bodies) up to a target surface.  This essentially extends or
    *   trims the attached surfaces of the sheet body.
    */

  virtual CubitStatus tweak_target( DLIList<Curve*> &curve_list,
                                    DLIList<Curve*> &target_curves, 
                                    DLIList<BodySM*> &new_bodysm_list, 
                                    CubitBoolean extend_flg = CUBIT_TRUE,
                                    CubitPlane *limit_plane = NULL,
                                    CubitBoolean reverse_flg = CUBIT_FALSE,
                                    CubitBoolean keep_old_body = CUBIT_FALSE,
                                    CubitBoolean preview = CUBIT_FALSE,
                                    double max_area_increase = 0 ) const;
  /**<  Tweak specified edges of a sheet body or bodies up to a target curve 
    *   that is part of a sheet body.  The target is a surface created by
    *   thickening the owning surface of the target curve.
    */

  virtual CubitStatus tweak_target( Point *point_ptr,
                                    DLIList<Surface*> &modify_surface_list,
                                    CubitVector &target_loc,
                                    BodySM *&new_bodysm_ptr,
                                    CubitBoolean keep_old_body = CUBIT_FALSE,
                                    CubitBoolean preview = CUBIT_FALSE ) const ;
  /**<  Tweak specified vertex of a sheet body to a given location.  The
    *   given vertex must be part of a planar surface or surfaces attached to
    *   linear curves only.  The user specified which of those surfaces to
    *   actually modify.  The given location will be projected to be on the
    *   given planar surface(s) before being used - this projected location
    *   must be the same on all surfaces.
    */

  virtual CubitStatus  split_curve( Curve* curve_to_split,
                                    const CubitVector& split_location,
                                    DLIList<Curve*>& created_curves ) ;
          //- Splits a curve at the specified location.
          //- the single passed in curve is split into two curves at the split location
          //- the two resulting curves are added to the passed in list


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

  virtual CubitStatus create_skin_surface( DLIList<Curve*>& curves, BodySM*& new_body, DLIList<Curve*>& ) const;

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

  virtual CubitStatus create_surface( DLIList<Point*> &points,
                                      BodySM *&new_body ) const;


  virtual CubitStatus create_weld_surface( CubitVector &root,
                                           Surface *ref_face1, double leg1, Surface *ref_face2, double leg2,
                                           BodySM *&new_body ) const;

  CubitStatus tolerant_imprint( DLIList<BodySM*> &bodies_in,
                                  DLIList<BodySM*> &new_bodies,
                                  DLIList<TopologyBridge*>*,
                                  DLIList<TopologyBridge*>* ) const; 

  virtual CubitStatus tolerant_imprint(DLIList<Surface*> &surfs_in,
                                       DLIList<BodySM*> &new_bodysm_list) const;

  virtual CubitStatus tolerant_imprint_surface_with_curves(
                                         Surface *surface_to_imprint,
                                         DLIList<Curve*> &curves,
                                         DLIList<TopologyBridge*> &temporary_bridges,
                                         BodySM *&new_body,
                                         DLIList<TopologyBridge*> *new_tbs = NULL,
                                         DLIList<TopologyBridge*> *att_tbs = NULL ) const;
  //Imprints a surface with passed-in curves.  Can imprint successfully
  //and expectedly with sloppy/dirty geometry.

  virtual CubitStatus remove_topology(DLIList<Curve*> &curves_to_remove,
                                       DLIList<Surface*> &surfs_to_remove,
                                       double backoff_distance,
                                       double small_edge_size,
                                       DLIList<BodySM*> &new_bodysm_list,
                                       CubitBoolean preview) const;

  virtual CubitStatus curve_surface_intersection( Surface *surface,
                                                  Curve* curve,
                                                  DLIList<Curve*> &new_curves ) const;
  //Intersects input surface with input curve to produce intersection curve(s).
  //If intersection results is nothing or a point, CUBIT_FAILURE is returned.

  virtual void get_possible_invalid_tbs(DLIList<TopologyBridge*> &bridges_in,
                             DLIList<TopologyBridge*> &bridges_out);

protected:

 TopoDS_Face* make_TopoDS_Face( GeometryType surface_type,
  	                        DLIList<DLIList<TopoDS_Edge*>*> topo_edges, 
                                Surface* old_surface_ptr) const;     

 int imprint_toposhapes(TopoDS_Shape*&, TopoDS_Shape*,
                                DLIList<TopoDS_Face*>& on_faces ) const;

 TopoDS_Edge* find_imprinting_edge(TopoDS_Shape& from_shape,
                                   TopoDS_Edge& tool_shape,
                                   DLIList<TopoDS_Face*>& faces)const;

 CubitStatus sort_curves(DLIList<Curve*> curve_list,
                         DLIList<DLIList<TopoDS_Edge*>*>& topo_edges_loops)const;

CubitStatus stitch_surfs(DLIList<BodySM*>& surf_bodies,
                         TopoDS_Shape& stitched_shape) const;

private:
 CubitStatus result_1_imprint(BodySM* from_body,
                              BodySM* tool_body,
                              BodySM*& newBody)const;

 CubitStatus result_3_imprint(BodySM* from_body,
                              BodySM* tool_body,
                              BodySM*& newBody)const;

    //- for periodic surfaces, use webcut first and then unite to get imprints

 virtual bool supports_interoperability() { return true; }
    //- Returns whether intermixing of real and virtual geometry operations
    //- is supported for the current geometry kernel.

 CubitStatus tweak_chamfer_sheet(Point* pnt,
                                 OCCSurface* face,
                                 double d1,
                                 Curve* edge1,
                                 double d2,
                                 Curve* edge2,
                                 DLIList<BodySM*> & new_bodysm_list,
                                 CubitBoolean keep_old_body,
                                 CubitBoolean preview ) ;

 CubitStatus tweak_fillet_chamfer_sheet( DLIList<Point*> & ref_vertex_list,
                               DLIList<OCCSurface*> faces,
                               double radius,
                               CubitBoolean is_fillet,
                               DLIList<BodySM*> & new_bodysm_list,
                               CubitBoolean keep_old_body,
                               CubitBoolean preview ) ;
 
 CubitStatus tweak_chamfer_solid( Point* point_ptr,
                                    OCCBody* body,
                                    double r1,
                                    Curve *c1,
                                    double r2,
                                    Curve *c2,
                                    double r3,
                                    Curve *c3,
                                    DLIList<BodySM *> &new_bodysm_list,
                                    CubitBoolean keep_old_body,
                                    CubitBoolean preview );

 CubitStatus tweak_chamfer_solid( DLIList<Point*> &point_list,
                                    DLIList<OCCBody*> &bodies,
                                    double radius,
                                    DLIList<BodySM*> &new_bodysm_list,
                                    CubitBoolean keep_old_body,
                                    CubitBoolean preview );

 CubitStatus sort_points_by_body_type( DLIList<Point*> &point_list,
                                         DLIList<Point*> &solid_points,
                                         DLIList<Point*> &sheet_points,
                                         DLIList<OCCSurface*> &s_list,
                                         DLIList<OCCBody*> &bodies );

 CubitStatus tweak_fillet( Curve * curve_ptr,
                           double start_radius,
                           double end_radius,
                           BodySM *& new_bodysm_ptr,
                           CubitBoolean keep_old_body,
                           CubitBoolean preview,
                           CubitBoolean if_fillet ) ;

 CubitStatus do_loft(BRepOffsetAPI_ThruSections& loft,
                     Surface * face1,
                     Surface * face2) const;

 CubitStatus do_loft(BRepOffsetAPI_ThruSections& loft,
                     DLIList<DLIList<TopoDS_Edge*>*> loops) const;

 void make_base_for_prim_pyramid(double major,
                                 double minor,
                                 double height,
                                 int sides,
                                 TopoDS_Wire& wire)const;

 int check_intersection(DLIList<TopoDS_Edge*>*& edge_list,
                        TopoDS_Face from_face)const;

 CubitStatus get_shape_list(DLIList<BodySM*>& BodySM_list,
                         DLIList<TopoDS_Shape*>& shape_list,
                         DLIList<CubitBoolean>& is_volume,
                         bool  keep_old,
                         DLIList<CubitBox*>* b_boxes = NULL) const;

 CubitStatus face_edge_imprint( DLIList<Surface*> &ref_face_list,
                                DLIList<Curve*> &edge_list,
                                DLIList<TopoDS_Face*>& face_list,
                                DLIList<TopoDS_Shape*>& shape_list,
                                bool keep_old ) const; 
 void shape_to_bodySM( DLIList<TopoDS_Shape*> shape_list,
                       DLIList<BodySM*>& new_body_list)const;

 void check_operation(TopoDS_Shape& cut_shape,
                      TopoDS_Shape*& from_shape, //output
                      CubitBoolean  is_volume,
                      CubitBoolean& has_changed, //output
                      BRepAlgoAPI_BooleanOperation* op,
                      CubitBoolean keep_old) const;

 CubitStatus get_sweepable_toposhape(OCCSurface*& surface,
                                     const CubitVector* sweep_v_p,
                                     TopoDS_Shape& toposhape)const;

 CubitStatus get_sweepable_toposhape(OCCCurve*& curve,
                                     TopoDS_Shape& toposhape)const;

 CubitStatus do_subtract(DLIList<BodySM*> &from_bodies,
                         DLIList<TopoDS_Shape*> &tool_bodies_copy, 
                         DLIList<CubitBoolean> &is_tool_volume,
                         DLIList<CubitBox*>* tool_boxes,
                         DLIList<BodySM*> &new_bodies,
                         bool keep_old,
                         bool imprint = CUBIT_FALSE) const;

 void get_new_tbs(
         std::map<OCCSurface*, std::pair<CubitVector, int> >& surf_property_map,
         std::map<OCCCurve*, std::pair<CubitVector, int> >& curve_property_map,
         DLIList<OCCPoint*> &points,
         DLIList<OCCSurface*> &new_surfaces,
         DLIList<OCCCurve*> &new_curves,
         DLIList<OCCPoint*> &new_points,
         DLIList<TopologyBridge*> *new_tbs)const; 

 void get_att_tbs(DLIList<OCCSurface*> &new_surfaces,
                  DLIList<OCCCurve*> &new_curves,
                  DLIList<OCCPoint*> &new_points,
                  const CubitString& name,
                  DLIList<TopologyBridge*> *att_tbs)const;

 CubitStatus split_shape_by_location(TopoDS_Shape *&from_shape,
                                     Curve* curve_to_split,
                                     const CubitVector& split_location,
                                     DLIList<Curve*>& created_curves )const;

  //- Removes all all unnessesary faces, curves, vertices and associated
 //- data from a refentity.
 virtual CubitStatus test_regularize_entity( GeometryEntity *old_entity_ptr)
 {return CUBIT_FAILURE; } 

 virtual CubitStatus create_offset_sheet( DLIList<Surface*> &surface_list,
                                    double offset_distance,
                                    DLIList<Surface*> *add_surface_list_ptr,
                                    DLIList<double> *add_offset_list_ptr,
                                    DLIList<BodySM*> &new_body_list,
                                    CubitBoolean preview = CUBIT_FALSE ) const
  /**< Create a sheet body (or bodies) by offsetting the given faces. The
    *  optional additional face list and double list (must be same length)
    *  allow different offset distances for different faces. Adjoining faces
    *  are extended or trimmed to remain joined in the new sheet body.  Radial
    *  faces that cannot be so offset are removed and the resulting wound
    *  healed by the surrounding faces.
    */
  { return CUBIT_FAILURE; }

  virtual CubitBoolean bodies_interfering( BodySM *body1,  BodySM *body2 ) const
  {return CUBIT_FAILURE; }

  virtual CubitStatus stitch( DLIList<BodySM*> &bodies_to_stitch,
                      DLIList<BodySM*> &new_bodies,
                      bool tighten_gaps,
                      double tolerance )const
  {return CUBIT_FAILURE; }

} ;

#endif
