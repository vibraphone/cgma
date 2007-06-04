#ifndef ACIS_HEALER_TOOL_HPP
#define ACIS_HEALER_TOOL_HPP

#include "CubitDefines.h"
#include "GeometryHealerEngine.hpp"

#ifdef UNIX_PORT
#include "undefwin_.h"
#endif
//#include "body.hxx"
#ifdef UNIX_PORT
#include "defwin_.h"
#endif

class RefVertex;
class RefEdge;
class RefEntity;
class CoEdge;
class Loop;
class RefFace;
class Shell;
class RefVolume;
template <class X> class DLIList;
class Body;
class RefFace;
class HealerTool;
class AcisModifyEngine;
class BODY;
class ENTITY;

class AcisHealerTool : public GeometryHealerEngine
{
public:
// ********** BEGIN FRIEND DECLARATIONS        **********

// ********** END FRIEND DECLARATIONS        **********

  ~AcisHealerTool();

  static AcisHealerTool* instance();
  //- Gives access to the singleton object of this class


  CubitStatus init_BODY_for_healing( BODY* BODY_ptr );
  //- Use this function to initialize a body before healing it.
  //- Aggregate attributes (containing such things as tolerances
  //- for healing) are attached to the body in this step.

  CubitStatus end_BODY_for_healing( BODY* BODY_ptr );
  //- Use this function to cleanup ACIS aggregate attributes from
  //- a body after healing.
  
  CubitStatus heal_BODY( BODY* BODY_ptr, int& percent_goodgeom_before,
                         int& percent_goodgeom_after,
                         int& num_splines_simplified,
                         CubitBoolean make_tolerant = CUBIT_FALSE,
                         FILE* logfile_ptr = NULL );
  //- Uses ACIS autohealer to heal the body.  The function 
  //- init_BODY_for_healing must be used on the BODY before calling
  //- this function.  Also, remember to call end_BODY_for_healing
  //- after using the function.  If logfile_ptr is non-NULL, then 
  //- overall healing statistics (before/after) are written to the 
  //- logfile for this body.

  CubitBoolean is_healer_engine( const TopologyBridge* );
  CubitBoolean is_healer_engine( const TopologyEntity* );
  // returns CUBIT_TRUE if the entity belongs to this healer engine

  CubitStatus auto_heal_bodies( DLIList<Body*> &body_list, DLIList<Body*> &new_body_list,
                                DLIList<TopologyEntity*> &bad_geometry,
                                CubitBoolean rebuild = CUBIT_FALSE, 
                                CubitBoolean keep_old = CUBIT_FALSE,
                                CubitBoolean make_tolerant = CUBIT_FALSE,
                                FILE* logfile_ptr = NULL );
  //- Uses the ACIS autohealer to heal the given body list.  The rebuild option
  //- can be used for more rigorous healing, where each surface is pulled off,
  //- healed, and then stitched back into a body.

   CubitStatus heal_bodies( DLIList<Body*> &body_list, DLIList<Body*> &new_body_list,
                            DLIList<TopologyEntity*> &bad_geometry,
                            CubitBoolean rebuild = CUBIT_FALSE,
                            CubitBoolean keep_old = CUBIT_FALSE,
                            CubitBoolean make_tolerant = CUBIT_FALSE,
                            FILE* logfile_ptr = NULL );
   //- Heals the input bodies, using either autoheal or incremental,
   //- depending on the switches that are setup.  Writes to the logfile
   //- if it's been opened.

  CubitStatus analyze_badgeom( DLIList<Body*> &body_list,
                               DLIList<TopologyEntity*> &bad_geometry,
                               FILE* logfile = NULL );
  //- Uses the ACIS healing husk to find bad geometry and give feedback to 
  //- the user.  A lot more could be done here - this is just a quick 
  //- overview for the user (the user will have no idea why the identified 
  //- geometry is bad).  A logfile can give a little more information.

  CubitStatus get_badgeom( DLIList<Body*> &body_list, DLIList<TopologyEntity*> &bad_geometry );
  //- Shows the bad geometry.  The geometry must have been analyzed first.
  //- If body_list is empty, shows only for those bodies that have been
  //- analyzed.

  CubitStatus get_tcurves( DLIList<Body*> &body_list, DLIList<RefEdge*> &t_curves );
  //- Gets tolerant curves.  If the body list is empty, all tolerant 
  //- curves are retrieved.

  CubitStatus heal_incremental( DLIList<Body*> &body_list, 
                                DLIList<Body*> &new_bodies,
                                DLIList<TopologyEntity*> &bad_geometry,
                                double simplify_tol = CUBIT_DBL_MAX,
                                double stitch_min_tol = CUBIT_DBL_MAX,
                                double stitch_max_tol = CUBIT_DBL_MAX,
                                double geombuild_tol = CUBIT_DBL_MAX,
                                double analytic_tol = CUBIT_DBL_MAX,
                                double isospline_tol = CUBIT_DBL_MAX,
                                double reblend_classify_tol = CUBIT_DBL_MAX,
                                double reblend_tol = CUBIT_DBL_MAX,
                                CubitBoolean keep_old = CUBIT_FALSE,
                                CubitBoolean make_tolerant = CUBIT_FALSE,
                                FILE* logfile_ptr = NULL);
  //- Uses the ACIS healing husk to perform one or more healing operations.
  //- Which operations are performed is determined by switches in the
  //- healer.  The user can set default tolerances (separate from this 
  //- function, override the defaults with this function, or just use the 
  //- tolerances calculated by the healer (preferred).
  //
  //- The possible healing steps are:
  //    1) preprocess - removes zero-length edges, sliver faces, duplicate 
  //                    vertices. This is the first step which is always 
  //                    done (user shouldn't suppress this).  Tolerance is SPAresabs.
  //    2) simplify - simplify NURBS into analytic.   
  //                   Default simplify_tol = .0001
  //    3) stitch - stitch geometry together.  Iterative from min to max 
  //                tolerance.
  //                   Default min tol = 10e-5
  //                   Default max tol = 1
  //    4) geombuild - geometry-related healing.  
  //                   Default geombuild_tol = .01
  //                   Default analytic_tangency_tol = .01
  //                   Default isolspline_solver_tol = .01      
  // Individual geombuild steps can be (instead of doing all):
  //          analytic - performs all of the stages of the analytic solver 
  //                     subphase of the geometry building phase. The analytic 
  //                     solver subphase attempts to heal all edges and 
  //                     vertices shared by analytic surfaces. 
  //          isospline - performs all of the stages of the isospline solver
  //                      subphase of the geometry building phase. The 
  //                      isospline solver attempts to heal all edges shared 
  //                      by tangential isoparametric surfaces (e.g., the 
  //                      intersection curve is an isoparametric curve of 
  //                      both splines in the intersection). 
  //          reblend - future option
  //          sharpedge - performs all of the stages of the sharp edge solver 
  //                      subphase of the geometry building phase. The sharp 
  //                      edge solver attempts to heal all edges and vertices 
  //                      that are shared by surfaces that intersect sharply. 
  //                      This includes nontangential surface junctions. 
  //          genericspline - performs all of the stages of the generic spline 
  //                          solver subphase of the geometry building phase. 
  //                          The generic spline solver attempts to heal 
  //                          generic tangential spline junctions, (e.g., the 
  //                          intersection curve is not an isoparametric curve 
  //                          of both splines in the intersection). 
  //          wrapup - handles remaining pcurves, wraps up
  //                   geometry buiilding phase (user shouldn't suppress)
  //    5) postprocess - correction of neg-area faces, duplicate vertices, edge groups
  //                     Last step - always done (user shouldn't suppress)

  void list_incremental();
  //- Lists the current settings for incremental healing

  void list_tolerances( DLIList<Body*> &body_list );
  //- Lists out the tolerances of each body, separately

  double get_default_simplify_tol();
  void set_default_simplify_tol( double tol );
  double get_default_stitch_min_tol();
  void set_default_stitch_min_tol( double tol );
  double get_default_stitch_max_tol();
  void set_default_stitch_max_tol( double tol );
  double get_default_geombuild_tol();
  void set_default_geombuild_tol( double tol );
  double get_default_analytic_tol();
  void set_default_analytic_tol( double tol );
  double get_default_isospline_tol();
  void set_default_isospline_tol( double tol );
  double get_default_reblend_classify_tol();
  void set_default_reblend_classify_tol( double tol );
  double get_default_reblend_tol();
  void set_default_reblend_tol( double tol );
  void reset_default_tolerances();
  void list_default_tolerances();
  //- Functions to set the default tolerances used.  The healer calculates
  //- the default tolerance per body.  These functions allow the user to override
  //- these defaults for all bodies healed.  In incremental healing, the user can 
  //- override these defaults by sending-in tolerances.  For autoheal, these
  //- defaults are used if they are set, otherwise the healer calculates
  //- intelligent defaults.
  
  void clean_attributes( DLIList<Body*>& body_list );
  //- Cleans healer attributes from the bodies.  These can be left if the 
  //- CleanAtt setting was used when doing analysis or healing.

  CubitBoolean get_cleanatt_flg();
  void set_cleanatt_flg( CubitBoolean flg );
  //- Get/set flags which determine whether to clean attributes after
  //- analysis and healing.
  
  int get_show_method(); // 0-none, 1-highlight, 2-draw
  void set_show_method( int method );
  //- Method determines how bad geometry is displayed

  CubitBoolean get_show_summary_flg();
  void set_show_summary_flg( CubitBoolean flg );
  //- Flag to determine whether to list a summary when showing bad geometry

  CubitBoolean get_show_details_flg();
  void set_show_details_flg( CubitBoolean flg );
  //- Flag to determine whether to list details when showing bad geometry

  CubitBoolean get_show_on_heal_flg();
  void set_show_on_heal_flg( CubitBoolean flg );
  //- Flag to determine whether to show bad geometry after healing

  CubitBoolean get_check_vol_on_heal_flg();
  void set_check_vol_on_heal_flg( CubitBoolean flg );
  double get_vol_on_heal_limit();
  void set_vol_on_heal_limit( double limit );
  //- Allow for checking of short curves after healing

  CubitBoolean get_check_surf_on_heal_flg();
  void set_check_surf_on_heal_flg( CubitBoolean flg );
  double get_surf_on_heal_limit();
  void set_surf_on_heal_limit( double limit );
  //- Allow for checking of small surfaces after healing

  CubitBoolean get_check_curve_on_heal_flg();
  void set_check_curve_on_heal_flg( CubitBoolean flg );
  double get_curve_on_heal_limit();
  void set_curve_on_heal_limit( double limit );
  //- Allow for checking of short curves after healing

  CubitBoolean get_show_bad_vertices_flg();
  void set_show_bad_vertices_flg( CubitBoolean flg );
  CubitBoolean get_show_bad_curves_flg();
  void set_show_bad_curves_flg( CubitBoolean flg );
  CubitBoolean get_show_bad_coedges_flg();
  void set_show_bad_coedges_flg( CubitBoolean flg );
  CubitBoolean get_show_bad_loops_flg();
  void set_show_bad_loops_flg( CubitBoolean flg );
  CubitBoolean get_show_bad_surfaces_flg();
  void set_show_bad_surfaces_flg( CubitBoolean flg );
  CubitBoolean get_show_bad_shells_flg();
  void set_show_bad_shells_flg( CubitBoolean flg );
  CubitBoolean get_show_bad_volumes_flg();
  void set_show_bad_volumes_flg( CubitBoolean flg );
  CubitBoolean get_show_bad_bodies_flg();
  void set_show_bad_bodies_flg( CubitBoolean flg );
  //- Get/set for flags for what to show during analysis/show

  void list_onshow_flgs();
  //- Function to list to user what the current onshow flags are

  CubitBoolean get_inc_preprocess_flg();
  void set_inc_preprocess_flg( CubitBoolean flg );
  CubitBoolean get_inc_simplify_flg();
  void set_inc_simplify_flg( CubitBoolean flg );
  CubitBoolean get_inc_stitch_flg();
  void set_inc_stitch_flg( CubitBoolean flg );
  CubitBoolean get_inc_geombuild_flg();
  void set_inc_geombuild_flg( CubitBoolean flg );
  CubitBoolean get_inc_analytic_flg();
  void set_inc_analytic_flg( CubitBoolean flg );
  CubitBoolean get_inc_isospline_flg();
  void set_inc_isospline_flg( CubitBoolean flg );
  CubitBoolean get_inc_reblend_flg();
  void set_inc_reblend_flg( CubitBoolean flg );
  CubitBoolean get_inc_sharpedge_flg();
  void set_inc_sharpedge_flg( CubitBoolean flg );
  CubitBoolean get_inc_genericspline_flg();
  void set_inc_genericspline_flg( CubitBoolean flg );
  CubitBoolean get_inc_wrapup_flg();
  void set_inc_wrapup_flg( CubitBoolean flg );
  CubitBoolean get_inc_postprocess_flg();
  void set_inc_postprocess_flg( CubitBoolean flg );
  //- Functions for controlling incremental healing

  void measure_filter(DLIList<RefEntity*> &ref_entities,
                      double measure_value,
                      double tolerance);
  
  CubitStatus force_simplify_to_plane( DLIList<RefFace*> &ref_face_list, 
                                       DLIList<Body*>& new_body_list, 
                                       CubitBoolean keep = CUBIT_FALSE );
  CubitStatus force_simplify_to_cylinder( DLIList<RefFace*> &ref_face_list, 
                                          DLIList<Body*>& new_body_list, 
                                          CubitBoolean keep = CUBIT_FALSE );
  CubitStatus force_simplify_to_cone( DLIList<RefFace*> &ref_face_list, 
                                      DLIList<Body*>& new_body_list, 
                                      CubitBoolean keep = CUBIT_FALSE );
  CubitStatus force_simplify_to_sphere( DLIList<RefFace*> &ref_face_list, 
                                        DLIList<Body*>& new_body_list, 
                                        CubitBoolean keep = CUBIT_FALSE );
  CubitStatus force_simplify_to_torus( DLIList<RefFace*> &ref_face_list, 
                                       DLIList<Body*>& new_body_list, 
                                       CubitBoolean keep = CUBIT_FALSE );
  //- Forces a spline surface to be an analytical of the type specified.

protected:
   AcisHealerTool();
   //- Class Constructor. (Not callable by user code. Class is constructed
   //- by the {instance()} member function.

private:
   
   CubitStatus analyze_BODY_for_healing( BODY* BODY_ptr );
   //- This function analyzes the BODY, attaching attributes to
   //- the entities in preparation for healing.  Use it to 
   //- determine what's bad in a body.  Be sure to call init_BODY_
   //- for_healing before and end_BODY_for_healing after using
   //- this function.

   BODY* auto_rebuild_BODY( BODY* BODY_ptr, int& percent_goodgeom_before,
                            int& percent_goodgeom_after, 
                            CubitBoolean make_tolerant = CUBIT_FALSE,
                            FILE* logfile_ptr = NULL );
   //- Pulls each surface off the body, heals it, stitches the body back together,
   //- and heals the result.  Function initializes/terminates healing
   //- internally.  If logfile_ptr is non-NULL, then overall healing statistics 
   //- (before/after) are written to the logfile for this body.

   CubitStatus percentage_goodgeom_before( BODY* BODY_ptr, int& percent_goodgeom );
   CubitStatus percentage_goodgeom_after( BODY* BODY_ptr, int& percent_goodgeom );
   //- Retrieves percentage of good geometry from the BODY, before and after
   //- healing.  Body must have been initialized for healing and analyzed, at 
   //- least, before this function is called.

   CubitStatus overall_printout( BODY* BODY_ptr, FILE* file_ptr);
   CubitStatus analysis_printout( BODY* BODY_ptr, FILE* file_ptr);
   CubitStatus calculate_printout( BODY* BODY_ptr, FILE* file_ptr);
   CubitStatus fix_printout( BODY* BODY_ptr, FILE* file_ptr);
   //- Functions to print healing data to a logfile.

   double get_simplify_tol(ENTITY* ent);
   void set_simplify_tol( double tol );
   //- The geometry simplification tolerance is the tolerance at which
   //- spline surfaces get simplified to analytic surfaces. If the tolerance is tight
   //- (as default is), only spline surfaces that are exact analytic surfaces get
   //- simplified. If the tolerance is loosened, then approximate analytic fits to
   //- splines are obtained. In such cases, the gaps between surfaces may increase and
   //- healing in subsequent operations may be more difficult. The need for increasing
   //- tolerances typically arises when analytic surfaces are output as NUBS surfaces
   //- rather than NURBS surfaces. The default tolerance is 0.0001 (length units),
   //- which obtains a very good approximation of analytic surfaces to spline surfaces
   //- in most models.
   
   double get_stitch_min_tol(ENTITY* ent);
   void set_stitch_min_tol( double tol );
   double get_stitch_max_tol(ENTITY* ent);
   void set_stitch_max_tol( double tol );
   //- The minimum and maximum stitching tolerances specify the range
   //- in which stitching between edges is performed. The stitching begins from the
   //- minimum tolerance and increases in steps towards the maximum tolerance.
   
   double get_geombuild_tol(ENTITY* ent);
   void set_geombuild_tol( double tol );
   //- The geometry build tolerance drives the actual geometry building
   //- of the model. This should typically be slightly more (around 3 times) than the
   //- maximum gap size in the model. The maximum gap size calculated during stitching
   //- is used as the geometry building tolerance for automatic healing (autoheal).
   //- However, the user may need to increase the geometry building tolerance if the
   //- healed geometry deviates substantially from the original geometry.
   
   double get_analytic_tol(ENTITY* ent);
   void set_analytic_tol( double tol );
   //- The analytic tolerance is used in the analytic solver subphase
   //- of geometry building, which fixes tangency constraints in the model. This is an
   //- upper bound for deviation of the analytic surfaces to satisfy the constraints.
   //- The default value of 0.01 permits translations of 0.01 to be performed to
   //- surfaces.
   
   double get_isospline_tol(ENTITY* ent);
   void set_isospline_tol( double tol );
   //- The isospline solver attempts to heal all edges shared by
   //- tangential isoparametric surfaces (e.g., the intersection curve is an
   //- isoparametric curve of both splines in the intersection). It calculates
   //- isoparametric junctions of spline geometries intersecting tangentially.

   double get_reblend_classify_tol(ENTITY* ent);
   void set_reblend_classify_tol( double tol );
   double get_reblend_tol(ENTITY* ent);
   void set_reblend_tol( double tol );
   //- Reblend tolerances. Currently there is no documentation on these but
   //- they are included for completeness anyway.

   void reset_switches();
   //- Resets all the tolerances that the user can set to default.

   CubitBoolean is_simplify_initialized(ENTITY* ent);
   CubitBoolean is_stitch_initialized(ENTITY* ent);
   CubitBoolean is_geombuild_initialized(ENTITY* ent);
   CubitBoolean is_initialized(ENTITY* ent);
   //- Query functions to find out if the body has been initialized yet, for each
   //- phase of healing.

   CubitStatus get_bad_vertices( BODY* BODY_ptr, DLIList<RefVertex*>& vertex_list);
   CubitStatus get_bad_edges( BODY* BODY_ptr, DLIList<RefEdge*>& edge_list);
   CubitStatus get_bad_coedges( BODY* BODY_ptr, DLIList<CoEdge*>& coedge_list);
   CubitStatus get_bad_loops( BODY* BODY_ptr, DLIList<Loop*>& loop_list);
   CubitStatus get_bad_faces( BODY* BODY_ptr, DLIList<RefFace*>& face_list);
   CubitStatus get_bad_shells( BODY* BODY_ptr, DLIList<Shell*>& shell_list);
   CubitStatus get_bad_volumes( BODY* BODY_ptr, DLIList<RefVolume*>& volume_list);
   //- Functions to retrieve bad geometry items from the given input BODY.
   //- BODY must have been initialized for healing and analyzed before 
   //- calling these functions.

   CubitStatus get_coedges_not_on_faces(BODY* body, DLIList<CoEdge*> coedge_list);
   // Get coedges that don't lie on surfaces
   CubitStatus get_coedges_no_partner(BODY *body, DLIList<CoEdge*> coedge_list);
   // Get coedges without a partner
   CubitStatus get_vertices_not_on_faces(BODY *body, DLIList<RefVertex*> ref_vertex_list);
   // Get vertices that do not lie on associated faces
   CubitStatus get_vertices_not_on_edges(BODY *body, DLIList<RefVertex*> ref_vertex_list);
   // Get vertices that do not lie on the associated edges
   CubitStatus get_vertices_edges_dont_meet(BODY *body, DLIList<RefVertex*> ref_vertex_list);
   // Get vertices where edges do not meet
   CubitStatus get_discontinuous_curves(BODY *body, DLIList<RefEdge*> ref_edge_list);
   // Get discontinuous curves
   CubitStatus get_degenerate_curves(BODY *body, DLIList<RefEdge*> ref_edge_list);
   // Get degenerate curves
   CubitStatus get_self_intersecting_curves(BODY *body, DLIList<RefEdge*> ref_edge_list);
   // Get self intersecting curves
   CubitStatus get_periodic_curves(BODY *body, DLIList<RefEdge*> ref_edge_list);
   // Get periodic curves
   CubitStatus get_closed_curves(BODY *body, DLIList<RefEdge*> ref_edge_list);
   // Get closed curves
   CubitStatus get_short_edges(BODY *body, DLIList<RefEdge*> ref_edge_list);
   // Get short edges (length less than geometry tolerance)
   CubitStatus get_tangent_edges(BODY *body, DLIList<RefEdge*> ref_edge_list);
   // Get tangent edges (edge joins two tangent faces)
   CubitStatus get_convex_edges(BODY *body, DLIList<RefEdge*> ref_edge_list);
   // Get convex edges 
   CubitStatus get_concave_edges(BODY *body, DLIList<RefEdge*> ref_edge_list);
   // Get concave edges
   CubitStatus get_loops_not_on_faces(BODY *body, DLIList<Loop*> loop_list);
   // Get loops that do not lie on their associated faces
   CubitStatus get_loops_disoriented(BODY *body, DLIList<Loop*> loop_list);
   // Get loops with incorrect orientation
   CubitStatus get_loops_gaps(BODY *body, DLIList<Loop*> loop_list);
   // Get loops that have coedge gaps
   CubitStatus get_loops_open (BODY *body, DLIList<Loop*> loop_list);
   // Get open loops
   CubitStatus get_discontinuous_surfaces(BODY *body, DLIList<RefFace*> ref_face_list);
   // Get discontinuous surfaces
   CubitStatus get_degenerate_surfaces(BODY *body, DLIList<RefFace*> ref_face_list);
   // Get degenerate surfaces
   CubitStatus get_self_intersecting_surfaces(BODY *body, DLIList<RefFace*> ref_face_list);
   // Get self intersecting surfaces
   CubitStatus get_periodic_surfaces(BODY *body, DLIList<RefFace*> ref_face_list);
   // Get periodic surfaces
   CubitStatus get_closed_surfaces (BODY *body, DLIList<RefFace*> ref_face_list);
   // Get closed surfaces

   CubitBoolean is_analyzed( BODY* BODY_ptr );
   //- Determines whether healing analysis attributes exist on the body.

   void get_curves_from_coedges( DLIList<CoEdge*> &coedge_list, DLIList<RefEdge*> &curve_list );
   void get_curves_from_loops( DLIList<Loop*> &loop_list, DLIList<RefEdge*> &curve_list );
   void get_surfaces_from_shells( DLIList<Shell*> &shell_list, DLIList<RefFace*> &surface_list );
   //- Convenience functions to get lower entities from higher ones.

   void print_none( FILE* file_ptr, const char* str );
   void print_vertices( FILE* logfile_ptr, DLIList<RefVertex*> &ref_vertex_list,  
                        const char* str_none, const char* pre_str );
   void print_curves( FILE* logfile_ptr, DLIList<RefEdge*> &curve_list,  
                      const char* str_none, const char* pre_str );
   void print_coedges( FILE* logfile_ptr, DLIList<CoEdge*> &coedge_list,  
                       const char* str_none, const char* pre_str );
   void print_loops( FILE* logfile_ptr, DLIList<Loop*> &loop_list,  
                     const char* str_none, const char* pre_str );
   void print_surfaces( FILE* logfile_ptr, DLIList<RefFace*> &surface_list,  
                        const char* str_none, const char* pre_str );
   void print_shells( FILE* logfile_ptr, DLIList<Shell*> &shell_list,  
                      const char* str_none, const char* pre_str );
   void print_volumes( FILE* logfile_ptr, DLIList<RefVolume*> &volume_list,  
                       const char* str_none, const char* pre_str );
   void print_bodies( FILE* logfile_ptr, DLIList<Body*> &body_list,  
                      const char* str_none, const char* pre_str );
   //- Functions used during analysis to list bad entities for the user.
   //- They can write to the logfile as well as do PRINT_INFO.  The 
   //- str_none is used for the no entities case, whereas the pre_str
   //- is used if entities exist, preceding the entity list.

   CubitStatus get_badgeom( CubitBoolean analyze_flag, DLIList<Body*> &body_list,
                            CubitBoolean after_heal, FILE* logfile_ptr,
                            DLIList<TopologyEntity*> &bad_geometry );
   //- Main workhorse for analyzing and showing bad geometry on bodies
   //- (i.e., doesn't just show, but can analyze too).

   void incremental_presummary( DLIList<Body*>& body_list, 
                                FILE* logfile_ptr = NULL,
                                CubitBoolean keep_old = CUBIT_FALSE );
   //- Gives some feedback to user as to which steps are being done
   //- for incremental healing, as well as tolerances used.

   CubitStatus force_simplify( int type, DLIList<RefFace*> &ref_face_list, 
                               DLIList<Body*>& new_body_list, CubitBoolean keep );
   // Workhorse function for force_simplify methods.  Type is in 
   // order of plane, cylinder, cone, sphere, torus (1-5).

   static AcisHealerTool* instance_;

   CubitBoolean cleanAtt;
   //- Determines whether to clean off attributes after analysis and
   //- healing.

   int showMethod; // 0-none, 1-highlight, 2-draw
   //- Method used for showing bad geometry during analyze, show, healing

   CubitBoolean showSummary;
   //- Flag to determine whether list a summary when showing bad geometry

   CubitBoolean showDetails;
   //- Flag to determine whether to show details when showing bad geometry

   CubitBoolean showOnHeal;
   //- Flag to determine whether to show while healing

   CubitBoolean checkVolOnHeal;
   CubitBoolean checkSurfOnHeal;
   CubitBoolean checkCurveOnHeal;
   double volLimit;
   double surfLimit;
   double curveLimit;

   CubitBoolean showBadVertices;
   CubitBoolean showBadCurves;
   CubitBoolean showBadCoEdges;
   CubitBoolean showBadLoops;
   CubitBoolean showBadSurfaces;
   CubitBoolean showBadShells;
   CubitBoolean showBadVolumes;
   CubitBoolean showBadBodies;
   //- Flags determine what to show during show

   CubitBoolean incPreprocess;  // This is normally always done
   CubitBoolean incSimplify;                        
   CubitBoolean incStitch;
   CubitBoolean incGeombuild;
   CubitBoolean incAnalytic;    // The next 6 are ignored if Geombuild is off
   CubitBoolean incIsospline;
   CubitBoolean incReblend;
   CubitBoolean incSharpedge;
   CubitBoolean incGenericspline;
   CubitBoolean incWrapup;      // Normally always done
   CubitBoolean incPostprocess; // Normally always done
   //- Flags for setting what is done in incremental healing


};

inline
CubitBoolean AcisHealerTool::get_cleanatt_flg()
{ return cleanAtt; }
  
inline 
void AcisHealerTool::set_cleanatt_flg( CubitBoolean flg )
{ cleanAtt = flg; }

inline
int AcisHealerTool::get_show_method()
{ return showMethod; }

inline
void AcisHealerTool::set_show_method( int method )
{ showMethod = method; }

inline
CubitBoolean AcisHealerTool::get_show_summary_flg()
{ return showSummary; }

inline
void AcisHealerTool::set_show_summary_flg( CubitBoolean flg )
{ showSummary = flg; }

inline
CubitBoolean AcisHealerTool::get_show_details_flg()
{ return showDetails; }

inline
void AcisHealerTool::set_show_details_flg( CubitBoolean flg )
{ showDetails = flg; }

inline
CubitBoolean AcisHealerTool::get_show_on_heal_flg()
{ return showOnHeal; }

inline
void AcisHealerTool::set_show_on_heal_flg( CubitBoolean flg )
{ showOnHeal = flg; }

inline
CubitBoolean AcisHealerTool::get_check_vol_on_heal_flg()
{ return checkVolOnHeal; }

inline
void AcisHealerTool::set_check_vol_on_heal_flg( CubitBoolean flg )
{ checkVolOnHeal = flg; }

inline
double AcisHealerTool::get_vol_on_heal_limit()
{ return volLimit; }

inline
void AcisHealerTool::set_vol_on_heal_limit( double limit )
{ volLimit = limit; }

inline
CubitBoolean AcisHealerTool::get_check_surf_on_heal_flg()
{ return checkSurfOnHeal; }

inline
void AcisHealerTool::set_check_surf_on_heal_flg( CubitBoolean flg )
{ checkSurfOnHeal = flg; }

inline
double AcisHealerTool::get_surf_on_heal_limit()
{ return surfLimit; }

inline
void AcisHealerTool::set_surf_on_heal_limit( double limit )
{ surfLimit = limit; }

inline
CubitBoolean AcisHealerTool::get_check_curve_on_heal_flg()
{ return checkCurveOnHeal; }

inline
void AcisHealerTool::set_check_curve_on_heal_flg( CubitBoolean flg )
{ checkCurveOnHeal = flg; }

inline
double AcisHealerTool::get_curve_on_heal_limit()
{ return curveLimit; }

inline
void AcisHealerTool::set_curve_on_heal_limit( double limit )
{ curveLimit = limit; }

inline
CubitBoolean AcisHealerTool::get_show_bad_vertices_flg()
{ return showBadVertices; }

inline
void AcisHealerTool::set_show_bad_vertices_flg( CubitBoolean flg )
{ showBadVertices = flg; }

inline
CubitBoolean AcisHealerTool::get_show_bad_curves_flg()
{ return showBadCurves; }

inline
void AcisHealerTool::set_show_bad_curves_flg( CubitBoolean flg )
{ showBadCurves = flg; }

inline
CubitBoolean AcisHealerTool::get_show_bad_coedges_flg()
{ return showBadCoEdges; }

inline
void AcisHealerTool::set_show_bad_coedges_flg( CubitBoolean flg )
{ showBadCoEdges = flg; }

inline
CubitBoolean AcisHealerTool::get_show_bad_loops_flg()
{ return showBadLoops; }

inline
void AcisHealerTool::set_show_bad_loops_flg( CubitBoolean flg )
{ showBadLoops = flg; }

inline
CubitBoolean AcisHealerTool::get_show_bad_surfaces_flg()
{ return showBadSurfaces; }

inline
void AcisHealerTool::set_show_bad_surfaces_flg( CubitBoolean flg )
{ showBadSurfaces = flg; }

inline
CubitBoolean AcisHealerTool::get_show_bad_shells_flg()
{ return showBadShells; }

inline
void AcisHealerTool::set_show_bad_shells_flg( CubitBoolean flg )
{ showBadShells = flg; }

inline
CubitBoolean AcisHealerTool::get_show_bad_volumes_flg()
{ return showBadVolumes; }

inline
void AcisHealerTool::set_show_bad_volumes_flg( CubitBoolean flg )
{ showBadVolumes = flg; }

inline
CubitBoolean AcisHealerTool::get_show_bad_bodies_flg()
{ return showBadBodies; }

inline
void AcisHealerTool::set_show_bad_bodies_flg( CubitBoolean flg )
{ showBadBodies = flg; }

// Incremental Healing
inline
CubitBoolean AcisHealerTool::get_inc_preprocess_flg()
{ return incPreprocess; }

inline
void AcisHealerTool::set_inc_preprocess_flg( CubitBoolean flg )
{ incPreprocess = flg; }

inline
CubitBoolean AcisHealerTool::get_inc_simplify_flg()
{ return incSimplify; }

inline
void AcisHealerTool::set_inc_simplify_flg( CubitBoolean flg )
{ incSimplify = flg; }

inline
CubitBoolean AcisHealerTool::get_inc_stitch_flg()
{ return incStitch; }

inline
void AcisHealerTool::set_inc_stitch_flg( CubitBoolean flg )
{ incStitch = flg; }

inline
CubitBoolean AcisHealerTool::get_inc_geombuild_flg()
{ return incGeombuild; }

inline
void AcisHealerTool::set_inc_geombuild_flg( CubitBoolean flg )
{ incGeombuild = flg; }

inline
CubitBoolean AcisHealerTool::get_inc_analytic_flg()
{ return incAnalytic; }

inline
void AcisHealerTool::set_inc_analytic_flg( CubitBoolean flg )
{ incAnalytic = flg; }

inline
CubitBoolean AcisHealerTool::get_inc_isospline_flg()
{ return incIsospline; }

inline
void AcisHealerTool::set_inc_isospline_flg( CubitBoolean flg )
{ incIsospline = flg; }

inline
CubitBoolean AcisHealerTool::get_inc_reblend_flg()
{ return incReblend; }

inline
void AcisHealerTool::set_inc_reblend_flg( CubitBoolean flg )
{ incReblend = flg; }

inline
CubitBoolean AcisHealerTool::get_inc_sharpedge_flg()
{ return incSharpedge; }

inline
void AcisHealerTool::set_inc_sharpedge_flg( CubitBoolean flg )
{ incSharpedge = flg; }

inline
CubitBoolean AcisHealerTool::get_inc_genericspline_flg()
{ return incGenericspline; }

inline
void AcisHealerTool::set_inc_genericspline_flg( CubitBoolean flg )
{ incGenericspline = flg; }

inline
CubitBoolean AcisHealerTool::get_inc_wrapup_flg()
{ return incWrapup; }

inline
void AcisHealerTool::set_inc_wrapup_flg( CubitBoolean flg )
{ incWrapup = flg; }

inline
CubitBoolean AcisHealerTool::get_inc_postprocess_flg()
{ return incPostprocess; }

inline
void AcisHealerTool::set_inc_postprocess_flg( CubitBoolean flg )
{ incPostprocess = flg; }

#endif
