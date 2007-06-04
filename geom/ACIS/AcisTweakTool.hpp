#ifndef ACIS_TWEAK_TOOL_HPP
#define ACIS_TWEAK_TOOL_HPP

#include "CubitDefines.h"

class CubitVector;
class Point;
class PointACIS;
class Curve;
class CurveACIS;
class Surface;
class SurfaceACIS;
class BodySM;
class BODY;
class AcisBridge;
class CubitPlane;
class AcisQueryEngine;
class AcisModifyEngine;
template <class X> class DLIList;

class AcisTweakTool
{
public:
// ********** BEGIN FRIEND DECLARATIONS        **********

// ********** END FRIEND DECLARATIONS        **********

  virtual ~AcisTweakTool();

  static AcisTweakTool* instance();
  //- Gives access to the singleton object of this class

  virtual CubitStatus tweak_chamfer( DLIList<Curve*> &curve_list, 
                                     double left_offset,
                                     DLIList<BodySM*> &new_bodysm_list,
                                     double right_offset = -1.0,
                                     CubitBoolean keep_old_body = CUBIT_FALSE,
                                     CubitBoolean preview = CUBIT_FALSE );
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
                                     CubitBoolean preview = CUBIT_FALSE );
  /**<  Chamfer vertices on solid or sheet bodies.  On a solid body there can
    *   be up to 3 offsets; on a sheet body up to 2 offsets.  The offsets are
    *   in the direction of the supplied edges.  If multiple vertices are 
    *   supplied, only one offset value is allowed and the edges are not used.
    */

  virtual CubitStatus tweak_fillet( DLIList<Curve*> &curve_list, 
                                    double radius,
                                    DLIList<BodySM*> &new_bodysm_list,
                                    CubitBoolean keep_old_body = CUBIT_FALSE,
                                    CubitBoolean preview = CUBIT_FALSE );
  /**<  Create a round fillet (or blend) at the given curves on solid bodies.
    */

  virtual CubitStatus tweak_fillet( Curve *curve_ptr, 
                                    double start_radius,
                                    double end_radius,
                                    BodySM *&new_body_ptr,
                                    CubitBoolean keep_old_body = CUBIT_FALSE,
                                    CubitBoolean preview = CUBIT_FALSE );
  /**<  Create a round fillet (or blend) at the given curve on a solid body.
    *   The fillet has a variable radius from the start to the end of the curve.
    */

  virtual CubitStatus tweak_fillet( DLIList<Point*> &point_list, 
                                    double radius,
                                    DLIList<BodySM*> &new_bodysm_list,
                                    CubitBoolean keep_old_body = CUBIT_FALSE,
                                    CubitBoolean preview = CUBIT_FALSE );
  /**<  Create a round fillet (or blend) at the given vertices on sheet bodies.
    */

  virtual CubitStatus tweak_move( DLIList<Surface*> &surface_list,
                                  const CubitVector &delta,
                                  DLIList<BodySM*> &new_bodysm_list,
                                  CubitBoolean keep_old_body = CUBIT_FALSE,
                                  CubitBoolean preview = CUBIT_FALSE );
  /**<  Tweak specified faces of a volume or volumes along a vector.
    */

  virtual CubitStatus tweak_move( DLIList<Curve*> &curve_list,
                                  const CubitVector &delta,
                                  DLIList<BodySM*> &new_bodysm_list,
                                  CubitBoolean keep_old_body = CUBIT_FALSE,
                                  CubitBoolean preview = CUBIT_FALSE );
  /**<  Tweak specified curves of a sheet body along a vector.
    */

  virtual CubitStatus tweak_offset( DLIList<Surface*> &surface_list,
                                    double offset_distance,
                                    DLIList<BodySM*> &new_bodysm_list,
                                    CubitBoolean keep_old_body = CUBIT_FALSE,
                                    CubitBoolean preview = CUBIT_FALSE );
  /**<  Tweak specified faces of a volume or volumes by offsetting those faces
    *   by the offset distance.
    */

  virtual CubitStatus tweak_offset( DLIList<Curve*> &curve_list,
                                    double offset_distance,
                                    DLIList<BodySM*> &new_bodysm_list,
                                    CubitBoolean keep_old_body = CUBIT_FALSE,
                                    CubitBoolean preview = CUBIT_FALSE );
  /**<  Tweak specified curves of a sheet body or bodies by offsetting those
    *   curves by the offset distance.
    */

  virtual CubitStatus tweak_remove( DLIList<Surface*> &surface_list,
                                    DLIList<BodySM*> &new_bodysm_list,
                                    CubitBoolean extend_adjoining = CUBIT_TRUE,
                                    CubitBoolean keep_surface = CUBIT_FALSE,
                                    CubitBoolean keep_old_body = CUBIT_FALSE,
                                    CubitBoolean preview = CUBIT_FALSE );
  /**<  Remove surfaces from a body or bodies and then by default extend
    *   the adjoining surfaces to fill the gap or remove the hole.  This
    *   requires functionality from the ACIS local operations husk.  This
    *   capability can usually successfully remove fillets, chamfers, blind
    *   holes and through holes from the ACIS body.
    */

  virtual CubitStatus tweak_remove( DLIList<Curve*> &curve_list,
                                    DLIList<BodySM*> &new_bodysm_list, 
                                    CubitBoolean keep_old_body = CUBIT_FALSE,
                                    CubitBoolean preview = CUBIT_FALSE );
  /**<  Remove curves from a sheet body or bodies and then extend the remaining
    *   curves to fill the gap.  If an internal loop of curves is removed the
    *   hole is removed.
    */

  virtual CubitStatus tweak_target( DLIList<Surface*> &surface_list,
                                    DLIList<Surface*> &target_surf_list,
                                    DLIList<BodySM*> &new_bodysm_list,
                                    CubitBoolean reverse_flg = CUBIT_FALSE,
                                    CubitBoolean keep_old_body = CUBIT_FALSE,
                                    CubitBoolean preview = CUBIT_FALSE );
  /**<  Tweak specified faces of a volume or volumes up to a set of target
    *   surfaces.  Topology is tweaked to the target surfaces.  Think of this
    *   as extending/trimming the body up past the target surfaces, then
    *   webcutting it off with the extended target surfaces and throwing away
    *   the excess.  The reverse flag should never be needed - if it is, there
    *   may be a bug or a bad normal on a body (i.e., negative volume body).
    */

  virtual CubitStatus tweak_target( DLIList<Curve*> &curve_list,
                                    DLIList<Surface*> &target_surf_list,
                                    DLIList<BodySM*> &new_bodysm_list,
                                    CubitBoolean reverse_flg = CUBIT_FALSE,
                                    CubitBoolean keep_old_body = CUBIT_FALSE,
                                    CubitBoolean preview = CUBIT_FALSE );
  /**<  Tweak specified edges of a surface or set of surfaces (in sheet bodies)
    *   up to a set of target surfaces.  This essentially extends or trims the
    *   attached surfaces of the sheet body.
    */

  virtual CubitStatus tweak_target( DLIList<Curve*> &curve_list,
                                    DLIList<Curve*> &target_curve_list, 
                                    DLIList<BodySM*> &new_bodysm_list,
                                    CubitBoolean reverse_flg = CUBIT_FALSE,
                                    CubitBoolean keep_old_body = CUBIT_FALSE,
                                    CubitBoolean preview = CUBIT_FALSE );
  /**<  Tweak specified edges of a sheet body or bodies up to a set of target
    *   curves that are part of a sheet body.  The target is a set of surfaces
    *   created by thickening the owning surfaces of the target curves.
    */

protected:
   AcisTweakTool();
   //- Class Constructor. (Not callable by user code. Class is constructed
   //- by the {instance()} member function.

private:

  CubitStatus tweak_target_single( DLIList<Surface*> &surface_list,
                                   Surface *target_surf_ptr,
                                   DLIList<BodySM*> &new_bodysm_list,
                                   CubitBoolean reverse_flg = CUBIT_FALSE,
                                   CubitBoolean keep_old_body = CUBIT_FALSE,
                                   CubitBoolean preview = CUBIT_FALSE );
  /**<  Tweak specified faces of a volume or volumes up to a single target
    *   surface.  Topology is tweaked to the target surface.  Think of this
    *   as extending/trimming the body up past the target surface, then
    *   webcutting it off with the extended target surface and throwing away
    *   the excess.  The reverse flag should never be needed - if it is, there
    *   may be a bug or a bad normal on a body (i.e., negative volume body).
    */

  CubitStatus tweak_target_multiple( DLIList<Surface*> &surface_list,
                                     DLIList<Surface*> &target_surf_list,
                                     DLIList<BodySM*> &new_bodysm_list,
                                     CubitBoolean reverse_flg = CUBIT_FALSE,
                                     CubitBoolean keep_old_body = CUBIT_FALSE,
                                     CubitBoolean preview = CUBIT_FALSE );
  /**<  Tweak specified faces of a volume or volumes up to a set of target
    *   surfaces.  Topology is tweaked to the target surfaces.  Think of this
    *   as extending/trimming the body up past the target surfaces, then
    *   webcutting it off with the extended target surfaces and throwing away
    *   the excess.  The reverse flag should never be needed - if it is, there
    *   may be a bug or a bad normal on a body (i.e., negative volume body).
    */

  CubitStatus tweak_target_single( DLIList<Curve*> &curve_list,
                                   Surface *target_surf_ptr,
                                   DLIList<BodySM*> &new_bodysm_list,
                                   CubitBoolean reverse_flg = CUBIT_FALSE,
                                   CubitBoolean keep_old_body = CUBIT_FALSE,
                                   CubitBoolean preview = CUBIT_FALSE );
  /**<  Tweak specified edges of a surface or set of surfaces (in sheet bodies)
    *   up to a single target surface.  This essentially extends or trims the
    *   attached surfaces of the sheet body.
    */

  CubitStatus tweak_target_multiple( DLIList<Curve*> &curve_list,
                                     DLIList<Surface*> &target_surf_list,
                                     DLIList<BodySM*> &new_bodysm_list, 
                                     CubitBoolean reverse_flg = CUBIT_FALSE,
                                     CubitBoolean keep_old_body = CUBIT_FALSE,
                                     CubitBoolean preview = CUBIT_FALSE );
  /**<  Tweak specified edges of a surface or set of surfaces (in sheet bodies)
    *   up to a set of target surfaces.  This essentially extends or trims the
    *   attached surfaces of the sheet body.
    */

  CubitStatus tweak_target_single( DLIList<Curve*> &curve_list,
                                   Curve *target_curve_ptr, 
                                   DLIList<BodySM*> &new_bodysm_list,
                                   CubitBoolean reverse_flg = CUBIT_FALSE,
                                   CubitBoolean keep_old_body = CUBIT_FALSE,
                                   CubitBoolean preview = CUBIT_FALSE );
  /**<  Tweak specified edges of a sheet body or bodies up to a single target
    *   curve that is part of a sheet body.  The target is a surface created
    *   by thickening the owning surface of the target curve.
    */

  CubitStatus tweak_target_multiple( DLIList<Curve*> &curve_list,
                                     DLIList<Curve*> &target_curve_list, 
                                     DLIList<BodySM*> &new_bodysm_list,
                                     CubitBoolean reverse_flg = CUBIT_FALSE,
                                     CubitBoolean keep_old_body = CUBIT_FALSE,
                                     CubitBoolean preview = CUBIT_FALSE );
  /**<  Tweak specified edges of a sheet body or bodies up to a set of target
    *   curves that are part of a sheet body.  The target is a set of surfaces
    *   created by thickening the owning surfaces of the target curves.
    */

  CubitStatus tweak_target_multiple( DLIList<FACE*> &source_FACE_list,
                                     DLIList<FACE*> &target_FACE_list,
                                     BODY *ext_target_BODY_ptr,
                                     DLIList<BodySM*> &debug_BodySM_list,
                                     CubitBoolean reverse_flg = CUBIT_FALSE );
  //- Tweak target multiple workhorse function.  The source FACE list must be
  //- from a single BODY.  The target FACE list can be from multiple BODIES. 
  //- The ext_target_BODY_ptr is passed in - you can create it from the targets
  //- by using the function "create_extended_sheet" (it is passed in instead of
  //- created in this function for efficiency since this function is likely
  //- called inside of a loop for multiple source BODIEs). The debug_BodySM_list
  //- is populated with debug BODIES if the debug flag 167 is turned on.

  CubitStatus get_ACIS_surfaces( DLIList<Surface*> &surface_list,
                                 DLIList<SurfaceACIS*> &acis_list );
   //- Get SurfaceACISs from Surfaces

   CubitStatus get_FACEs( DLIList<Surface*> &surface_list,
                          DLIList<FACE*> &FACE_list );
   //- Get FACEs from Surfaces

   CubitStatus get_EDGEs( DLIList<Curve*> &curve_list,
                          DLIList<EDGE*> &EDGE_list );
   //- Get EDGEs from Curves

   CubitStatus get_VERTICEs( DLIList<Point*> &point_list,
                             DLIList<VERTEX*> &VERTEX_list );
   //- Get VERTICEs from Points 

   CubitStatus tweak_target_single( DLIList<EDGE*> &input_EDGE_list,
                                    FACE *target_FACE_ptr, 
                                    DLIList<BodySM*> &new_body_list,
                                    CubitBoolean reverse_flg = CUBIT_FALSE,
                                    CubitBoolean keep_old_bodies = CUBIT_FALSE,
                                    CubitBoolean preview = CUBIT_FALSE,
                                    DLIList<AcisBridge*> *t_ab_list = NULL );
   //- Tweak EDGES workhorse function.  Note this function will return SUCCESS
   //- if any new bodies were created (i.e., it could fail on some and still
   //- return SUCCESS).  Last argument is needed for preview - if target was
   //- derived from an EDGE send in the AcisBridge of that EDGE - preview 
   //- needs to know about this EDGE if it exists on the input BODY itself.

   CubitStatus tweak_target_multiple( DLIList<EDGE*> &input_EDGE_list,
                                      DLIList<FACE*> &target_FACE_list, 
                                      DLIList<BodySM*> &new_bodysm_list,
                                      CubitBoolean reverse_flg = CUBIT_FALSE,
                                      CubitBoolean keep_old_bodies = CUBIT_FALSE,
                                      CubitBoolean preview = CUBIT_FALSE,
                                      DLIList<AcisBridge*> *t_ab_list = NULL );
   //- Tweak EDGES to multiple FACEs workhorse function.  Note this function
   //- will return SUCCESS if any new bodies were created (i.e., it could fail
   //- on some and still return SUCCESS).  Last argument is needed for preview -
   //- if targets were derived from an EDGE list send in the AcisBridges of
   //- those EDGEs - it is possible (unlikely though) that these EDGEs (if
   //- part of the BODY itself) could exist in the final result and are needed
   //- for the preview.

   CubitStatus get_thickened_BODIES_of_EDGES( const char *command_name,
                                              DLIList<EDGE*> &EDGE_list,
                                              DLIList<EDGE*> &removed_EDGE_list,
                                              BODY *&common_BODY_ptr,
                                              DLIList<BODY*> &thickened_BODY_list,
                                              DLIList<DLIList<EDGE*>*> &output_EDGE_lists,
                                              DLIList<DLIList<FACE*>*> &output_FACE_lists,
                                              DLIList<DLIList<FACE*>*> &conjugate_FACE_lists,
                                              double thickness = 0.2 );
   //- Get thickened BODIES from an input list of EDGEs.  This function can be
   //- called multiple times on the same list of EDGEs - each time, the EDGEs 
   //- from the common body that is thickened are removed from the input EDGE
   //- list.  The thickened BODIES are copied from the sheet the EDGEs are
   //- attached to (multiple BODIES can be returned because the returned
   //- thickened BODIES must be nonmanifold).  For each BODY returned, also
   //- output a list of EDGES (corresponding to the original input EDGEs), a
   //- list of FACEs (corresponding to the original FACEs of the input BODY)
   //- and a list of conjugate FACEs (the "side" FACEs of the thickened BODY -
   //- these are the FACEs we can tweak).

   CubitStatus get_thickened_BODIES_of_VERTICES( const char *command_name,
                                                 DLIList<VERTEX*> &input_VERTEX_list,
                                                 BODY *&common_BODY_ptr,
                                                 DLIList<BODY*> &thickened_BODY_list,
                                                 DLIList<DLIList<FACE*>*> &output_FACE_lists,
                                                 DLIList<DLIList<VERTEX*>*> &output_VERTEX_lists,
                                                 DLIList<DLIList<EDGE*>*> &output_EDGE_lists,
                                                 double thickness = 0.2 );
   //- Get thickened BODIES from an input list of VERTICEs.  This function can
   //- be called multiple times on the same list of VERTICEs - each time, the
   //- VERTICEs from the common body that is thickened are removed from the
   //- input EDGE list.  These thickened BODIES are copied from the sheets the
   //- VERTICEs are attached to (multiple BODIES can be returned because the 
   //- thickened BODIES must be nonmanifold).  For each BODY returned, also
   //- output a list of the FACEs (corresponding to the original FACEs of the
   //- sheet body), the VERTICEs (corresponding to the original VERTICEs), and
   //- a list of EDGEs created by the thickening process (from sweeping the
   //- original VERTEX - note it is possible that no EDGE was created, in which
   //- case a NULL value will exist in the list).

   CubitStatus copy_FACES_from_BODY( DLIList<FACE*> &FACE_list, 
                                     BODY *&copied_BODY_ptr);
   //- Input FACE list must be in a single BODY.  Copy these FACEs off into
   //- a new BODY, while retaining all the Cubit attributes.  The input
   //- BODY (the parent BODY of the FACE_list) is NOT modified.

   CubitStatus remove_FACES_from_BODY( BODY *BODY_ptr, 
                                       DLIList<FACE*> &remove_FACE_list );
   //- Remove the given FACEs from the BODY (modifying the input BODY).

   CubitStatus remove_FACES_from_BODY_except( BODY *BODY_ptr, 
                                              DLIList<FACE*> &keep_FACE_list );
   //- Remove all FACEs from the BODY except those in the keep_FACE_list
   //- (modifying the input BODY).

   CubitStatus thicken_BODY( BODY *BODY_ptr, double thickness );
   //- Thicken a BODY.  Remove CUBIT owner attributes from all VERTICEs, EDGEs,
   //- and FACEs except the original VERTICEs, EDGEs and FACEs.

   VERTEX *find_corresponding_VERTEX( VERTEX *ref_VERTEX_ptr, DLIList<VERTEX*> &VERTEX_list );
   //- Find the corresponding VERTEX in the input VERTEX_list to the ref_VERTEX_ptr
   //- using CUBIT owner attributes.  For example, copy a BODY then use this
   //- function to find a corresponding VERTEX from the parent BODY on the copied
   //- BODY.

   EDGE *find_corresponding_EDGE( EDGE *ref_EDGE_ptr, DLIList<EDGE*> &EDGE_list );
   //- Find the corresponding EDGE in the input EDGE_list to the ref_EDGE_ptr
   //- using CUBIT owner attributes.  For example, copy a BODY then use this
   //- function to find a corresponding EDGE from the parent BODY on the copied
   //- BODY.

   LOOP *find_corresponding_LOOP( LOOP *ref_LOOP_ptr, DLIList<LOOP*> &LOOP_list );
   //- Find the corresponding LOOP in the input LOOP_list to the ref_LOOP_ptr
   //- using CUBIT owner attributes.  For example, copy a BODY then use this
   //- function to find a corresponding LOOP from the parent BODY on the copied
   //- BODY.

   FACE *find_corresponding_FACE( FACE *ref_FACE_ptr, DLIList<FACE*> &FACE_list );
   //- Find the corresponding FACE in the input FACE_list to the ref_FACE_ptr
   //- using CUBIT owner attributes.  For example, copy a BODY then use this
   //- function to find a corresponding FACE from the parent BODY on the copied
   //- BODY.

   FACE *find_corresponding_FACE( AcisBridge *ab_ptr, BODY *BODY_ptr);
   //- Find FACE with given owner in the BODY_ptr

   LOOP *find_corresponding_LOOP( AcisBridge *ab_ptr, BODY *BODY_ptr);
   //- Find LOOP with given owner in the BODY_ptr

   CubitStatus get_owner_list( DLIList<EDGE*> &EDGE_list, 
                               DLIList<AcisBridge*> &owner_list );
   //- Get a list of Cubit owners corresponding to the input EDGE list.
   //- Appends (unique) to the input list.

   CubitStatus get_owner_list( DLIList<FACE*> &FACE_list, 
                               DLIList<AcisBridge*> &owner_list );
   //- Get a list of Cubit owners corresponding to the input FACE list.
   //- Appends (unique) to the input list.

   CubitStatus get_owner_list( DLIList<LOOP*> &LOOP_list, 
                               DLIList<AcisBridge*> &owner_list );
   //- Get a list of Cubit owners corresponding to the input LOOP list.
   //- Appends (unique) to the input list.

   CubitStatus get_corresponding_FACE_list( DLIList<AcisBridge*> &owner_list,
                                            BODY *BODY_ptr,
                                            DLIList<FACE*> &corresponding_FACE_list );
   //- Find the FACEs in BODY_ptr corresponding to the given owners

   CubitStatus get_corresponding_LOOP_list( DLIList<AcisBridge*> &owner_list,
                                            BODY *BODY_ptr,
                                            DLIList<LOOP*> &corresponding_LOOP_list );
   //- Find the LOOPs in BODY_ptr corresponding to the given owners

   CubitStatus prep_for_surface_swap( BODY *thickened_BODY_ptr,
                                      BODY *copied_input_BODY_ptr,
                                      DLIList<AcisBridge*> &owner_FACE_list );
   //- Prepare to swap the non-extended surfaces in the (original)
   //- copied_input_BODY_ptr with the extended surfaces in the
   //- thickened_BODY_ptr.  This function removes ALL BUT the input
   //- FACEs from the thickend_BODY_ptr (so we just have the extended or
   //- trimmed FACEs left), and removes the input FACEs from the 
   //- copied_input_BODY_ptr.  These two BODIEs can be united later for the
   //- end result - the sheet BODY surfaces have been extended/trimmed/etc.!

   CubitStatus unite_BODIES( BODY *copied_input_BODY_ptr,
                             DLIList<BODY*> &thickened_BODY_list,
                             BODY *&output_BODY_ptr );
   //- Unite all the given BODIES.

   CubitStatus get_owner_attribs( BODY *BODY_ptr, AcisBridge *&ab_body_ptr,
     AcisBridge *&ab_lump_ptr, AcisBridge *&ab_shell_ptr);
   //- Get owner attributes on the BODY, LUMP, SHELL

   CubitStatus reset_owner_attribs( BODY *BODY_ptr, AcisBridge *ab_body_ptr,
     AcisBridge *ab_lump_ptr, AcisBridge *ab_shell_ptr);
   //- Reset owner attributes on the BODY, LUMP, SHELL, since some operations
   //- remove them.

   CubitStatus sort_points_by_body_type( DLIList<Point*> &point_list, 
     DLIList<Point*> &solid_points, DLIList<Point*> &sheet_points );
   //- Sort the incoming point list by the type of owning body of the points -
   //- solids or sheet bodies.

   CubitStatus tweak_chamfer_solid( DLIList<Point*> &point_list, 
                                    double radius, 
                                    DLIList<BodySM*> &new_bodysm_list,
                                    CubitBoolean keep_old_body = CUBIT_FALSE,
                                    CubitBoolean preview = CUBIT_FALSE );
  //- Chamfer vertices on solid bodies

   CubitStatus tweak_chamfer_fillet_sheet( DLIList<Point*> &point_list, 
                                           double radius, 
                                           int type,
                                           DLIList<BodySM*> &new_bodysm_list,
                                           CubitBoolean keep_old_body = CUBIT_FALSE,
                                           CubitBoolean preview = CUBIT_FALSE );
   //- Chamfer or fillet vertices on sheet bodies.
   //- Type = 1=chamfer
   //-        2=fillet

   CubitStatus tweak_chamfer_solid( Point* point_ptr, 
                                    double r1,                                      
                                    Curve *c1,
                                    double r2,
                                    Curve *c2,
                                    double r3,
                                    Curve *c3,
                                    BodySM *&new_bodysm_ptr,
                                    CubitBoolean keep_old_body = CUBIT_FALSE,
                                    CubitBoolean preview = CUBIT_FALSE );
   //- Chamfer a vertex on a solid body with variable radii.  Radii correspond
   //- to given curves.

   CubitStatus tweak_chamfer_sheet( Point* point_ptr, 
                                    double r1,                                      
                                    Curve *c1,
                                    double r2,
                                    Curve *c2,
                                    BodySM *&new_bodysm_ptr,
                                    CubitBoolean keep_old_body = CUBIT_FALSE,
                                    CubitBoolean preview = CUBIT_FALSE );
   //- Chamfer a vertex on a sheet body with variable radii.  Radii correspond
   //- to given curves.

   CubitStatus assign_tweak_attribs( BODY *BODY_ptr, const char *att_name,
     DLIList<FACE*> &FACE_list, DLIList<AcisBridge*> &ab_FACE_list, 
     DLIList<EDGE*> &EDGE_list, DLIList<AcisBridge*> &ab_EDGE_list,
     DLIList<VERTEX*> &VERTEX_list, DLIList<AcisBridge*> &ab_VERTEX_list );
   //- Assign named attribs to faces, edges and vertices so we can get back to
   //- them after a fillet or chamfer.  These operations split these entities
   //- to make room for the new surfaces, and the original entities (and
   //- owner atributes) are lost. Thus we copy the owner atts during the split.
   //- This ultimately will allow us to preserve Cubit owners through a chamfer
   //- or fillet operation

   CubitStatus find_corresponding_entity_from_tweak_attrib( BODY *BODY_ptr, 
                                                            const char *att_name,
                                                            const int ent_type, 
                                                            int ent_integer,
                                                            ENTITY *&output_ENTITY_ptr );
   //- Find a corresponding entity to the input ENTITY_ptr using the tweak
   //- attributes assigned by assign_tweak_attribs.
   //- Input is of type ent_type, either FACE_TYPE, EDGE_TYPE, VERTEX_TYPE

   CubitStatus reassign_cubit_owners_from_tweak_attribs( BODY *BODY_ptr,
     const char *att_name,
     DLIList<FACE*> &FACE_list, DLIList<AcisBridge*> &ab_FACE_list, 
     DLIList<EDGE*> &EDGE_list, DLIList<AcisBridge*> &ab_EDGE_list,
     DLIList<VERTEX*> &VERTEX_list, DLIList<AcisBridge*> &ab_VERTEX_list );
   //- Using the tweak attributes assigned by assign_tweak_attribs, reset
   //- the Cubit owners back on the proper entities.

   CubitStatus remove_named_attribs( BODY *BODY_ptr, const char *name );
   //- Remove the named attributes from FACEs, EDGEs and VERTICEs

   CubitStatus blend_edges( DLIList<EDGE*> EDGE_list, double radius );
   //- Perform a blend operation.  This maintains Cubit owners where
   //- possible through the operation.

   CubitStatus chamfer_edges( DLIList<EDGE*> EDGE_list, double r1, double r2 = -1.0 );
   //- Perform a chamfer operation.  This maintains Cubit owners where
   //- possible through the operation.

   CubitStatus chamfer_vertices( DLIList<VERTEX*> VERTEX_list, double radius );
   //- Perform a chamfer operation.  This maintains Cubit owners where
   //- possible through the operation.

   CubitBoolean FACE_surrounded( FACE *ref_FACE_ptr, 
                                 DLIList<FACE*> &FACE_list );
   //- Determine if given (reference) FACE is entirely surrounded by 
   //- the FACEs in FACE_list (reference FACE can be contained in
   //- FACE_list).

   CubitStatus get_outer_EDGEs( FACE *FACE_ptr, DLIList<EDGE*> &EDGE_list );
   //- Get EDGEs from outer LOOP of FACE

   CubitStatus get_EDGES_by_BODY( DLIList<Curve*> &input_curve_list, 
                                  DLIList<BODY*> &BODY_list, 
                                  DLIList<DLIList<EDGE*>*> &BODY_EDGE_lists );
   CubitStatus get_EDGES_by_BODY( DLIList<EDGE*> &input_EDGE_list, 
                                  DLIList<BODY*> &BODY_list, 
                                  DLIList<DLIList<EDGE*>*> &BODY_EDGE_lists );
   //- Get separate lists of EDGEs from the input EDGEs by common BODY.  Be
   //- sure to free the individual lists in BODY_EDGE_lists when you are done
   //- with them, as this function allocates the memory for these lists.

   CubitStatus all_complete_internal_loops( DLIList<EDGE*> &BODY_EDGE_list,
                                            DLIList<LOOP*> &LOOP_list );
   //- Determine if given EDGEs form complete internal LOOPs.  LOOPs must be on
   //- a sheet BODY.  Return those LOOPs.  All given EDGEs must be in the same
   //- BODY.

   CubitStatus remove_LOOPs( DLIList<LOOP*> &LOOP_list );
   //- Remove the given internal LOOPs (from sheet BODIES).

   CubitStatus remove_holes( BODY *sheet_BODY_ptr );
   //- Remove all of the holes (internal LOOPs) in the given sheet BODY.                     

   CubitStatus extrema_pln_BODY( CubitPlane &plane, BODY *BODY_ptr,
                                 double &extrema_dist, int back_side = 0 );
   //- Finds the extrema distance from a plane to an ACIS BODY (perpendicular
   //- distance from plane to farthest extent of BODY), on one side of the
   //- plane (by default on the front of the plane - direction of plane
   //- normal).  If the entire BODY lies on the other side of the plane, the
   //- extrema distance is 0.0.

   CubitStatus tweak_FACEs_to_target( DLIList<FACE*> &tweak_FACE_list, 
     FACE *target_FACE, CubitBoolean reverse_flg = CUBIT_FALSE,
     CubitBoolean skip_self_int_check = CUBIT_FALSE );
   //- Tweak the given FACEs (which must all be part of the same BODY) to the
   //- given target FACE.  The reverse flag should never be needed - if it is,
   //- there may be a bug or a bad normal on a body (i.e., negative volume
   //- body.  However, it is here as a workaround if needed.

   CubitVector surface_normal( FACE *FACE_ptr, CubitVector &location,
                               CubitBoolean calc_loc = CUBIT_TRUE );
   //- Get normal of underlying surface of given FACE.  If calc_loc is 
   //- CUBIT_TRUE, the location vector is calculated as the closest location on
   //- the FACE to the center of the bounding box of the FACE, and the normal
   //- is found there.  If calc_loc is CUBIT_FALSE, the normal is found at the
   //- given location.  Note the normal is NOT adjusted for the FACE sense.

   CubitVector FACE_normal( FACE *FACE_ptr, CubitVector &location,
                            CubitBoolean calc_loc = CUBIT_TRUE );
   //- Get normal of given FACE.  If calc_loc is CUBIT_TRUE, the location
   //- vector is calculated as the closest location on the FACE to the center
   //- of the bounding box of the FACE, and the normal is found there.  If
   //- calc_loc is CUBIT_FALSE, the normal is found at the given location.
   //- Note the normal is adjusted for the FACE sense.

   CubitVector weighted_average_FACE_normal( DLIList<FACE*> &FACE_list );
   //- Get weighted average normal of given FACEs.  The normals are weighted
   //- by a factor related to FACE area utilizing the graphics facets.

   CubitStatus weighted_average_FACE_normal( FACE *FACE_ptr, CubitVector &norm,
                                             double &weight );
   //- Get weighted average normal of give FACE.  The normal is weighted
   //- by a factor related to FACE area utilizing the graphics facets.  Note 
   //- for a surface like a full cylinder, this can result in a 0,0,0 normal
   //- vector.

   CubitStatus create_offset_planar_body( CubitPlane &plane, double offset,
                                          BODY *&planar_BODY_ptr );
   //- Create a planar ACIS BODY from the given CubitPlane and offset.  The
   //- offset is in the direction of the plane's normal.  The plane is 
   //- arbitrarily 10x10x10 in size near the origin.

   CubitStatus create_extended_sheet( DLIList<FACE*> &FACE_list, 
                                      BODY *&ext_BODY_ptr );
   //- Extend out a set of FACEs and place into a new sheet BODY.

   CubitStatus prep_cutting_sheet( BODY *&cutting_sheet, 
                                   BODY *tweaked_BODY_ptr,
                                   BODY *ext_target_BODY_ptr, 
                                   CubitVector &source_norm,
                                   CubitVector &target_norm,
                                   CubitBoolean check_crossing = CUBIT_TRUE );
   //- Prepare the cutting sheet for chopping off the BODY in a
   //- tweak_target_multiple operation.  Inputs include the cutting sheet,
   //- the tweaked BODY (the BODY that will be cut), the extended target BODY
   //- with consistent normals and reference source and target normals.  The
   //- cutting sheet can optionally be trimmed to including only FACEs crossed
   //- by the extended target (the "side" FACEs) that we want to cut.  Note
   //- named "tweak" attributes must exist on the side FACEs prior to calling
   //- this function if "check_crossing" is true.

   CubitStatus chop_off_with_sheet( BODY *BODY_ptr, BODY *sheet_BODY_ptr );
   //- Chop off all the material in a BODY to one side of a sheet body.  The
   //- material is removed on the side of the sheet body according to the
   //- surface normals of the sheet body (material removed in the direction of
   //- the surface normals of the sheet body, which must be consistent).

   CubitStatus get_overlap_area( BODY *BODY_ptr1, BODY *BODY_ptr2, 
                                 double &overlap_area, double accuracy = 1e-5 );
   //- Determine the area of overlap of two BODIES.  Uses an ACIS intersection
   //- boolean so is quite expensive.

   CubitStatus copy_FACEs_into_sheet( DLIList<FACE*> &FACE_list, 
                                      BODY *&sheet_ptr,
                                      CubitBoolean heal = CUBIT_TRUE );
   //- Copies the given FACEs into a sheet BODY.  The FACEs are united together
   //- into the new BODY.  Note a regularized boolean is performed (so the body
   //- is regularized).  By default the sheet BODY is healed, but this can be
   //- turned off.  Note healing will make the normals consistent, unless there
   //- is nonmanifold geometry in the body.

   CubitStatus remove_per_nonmanifold_FACEs( BODY *cutting_tool );
   //- Removes undesired nonmanifold FACEs in the cutting sheet.  See 
   //- diagram in function for one such case.  If any nonmanifold edges exist
   //- in the model that can't be removed, CUBIT_FAILURE is returned.

   CubitStatus tangent_outdir( FACE *FACE_ptr, EDGE *EDGE_ptr, 
                               CubitVector &pos, CubitVector &tangent_outvec );
   //- Finds the direction tangent to and pointing away from the FACE boundary
   //- at the given location on the EDGE.

   CubitStatus remove_aligned_periodic_FACEs( BODY *cutting_tool, 
                                              CubitVector &basis_vec );
   //- Removes cone, sphere or torus FACEs that have normals aligned with the
   //- given basis vector.

   CubitStatus heal_BODY( BODY *BODY_ptr );
   //- Use the ACIS healer to heal the given sheet BODY.

   CubitStatus align_normals( BODY *child, BODY *master );
   //- Align the normals of the child's FACEs to the master FACEs.  It is
   //- assumed the child is a subset of the master (i.e., the master is an
   //- extended BODY).

   FACE * find_overlap_FACE( FACE *FACE_ptr, DLIList<FACE*> &FACE_list );
   //- Find the FACE in FACE_list that overlaps with FACE_ptr.  Returns
   //- NULL if no overlap FACE found.

   CubitStatus make_surf_normals_consistent( BODY *BODY_ptr, 
                                             FACE *seed_FACE = NULL );
   //- Make the shell normals in the given sheet BODY consistent.  The optional
   //- seed FACE must be part of the BODY.  If none given, get_seed_FACE (see
   //- below) is used to select a seed FACE.

   CubitStatus append_neighbors( DLIList<FACE*> &FACE_list );
   //- Helper function for make_shell_normals_consistent

   FACE *get_seed_FACE( DLIList<FACE*> &FACE_list );
   //- Get a seed FACE from the given list of FACEs.  Only a seed FACE without
   //- a named "tweak" attribute on it will be returned.  Also, it may not be
   //- of advantage, but we try to find a FACE that is not of type CONE, TORUS
   //- or SPHERE.  If that is not possible, we just return the first FACE in
   //- the list with no "tweak" attribute on it.
   
   CubitStatus draw_tweak_preview_omt( BODY *BODY_ptr, CubitBoolean flush,
                                       DLIList<AcisBridge*> *prev_ab_list = NULL );
   //- Draw the preview EDGEs for tweak offset, move and target, for solid and
   //- sheetbodies.  The function works by finding EDGEs or FACEs without Cubit
   //- owner attributes on them and drawing those EDGEs as well as the EDGEs
   //- attached to them.  The "prev_ab_list" may be needed for special cases
   //- where the target is on the BODY itself.  I recommend always sending in
   //- the sources and target entities in this list.

   CubitStatus tag_tweak_remove_FACEs_for_preview( BODY *BODY_ptr, 
                                            DLIList<FACE*> &remove_FACE_list );
   //- Add an attribute ("tweak_preview") to the required FACEs for a tweak
   //- remove.  For the preview we draw the EDGEs on the surviving FACEs
   //- adjoining those that are removed.  Input is the copied BODY, along with
   //- the FACEs being removed.

   CubitStatus tag_tweak_remove_FACEs_for_preview( BODY *BODY_ptr, 
                                            DLIList<EDGE*> &remove_EDGE_list );
   //- Add an attribute ("tweak_preview") to the required FACEs for a tweak
   //- remove.  For the preview we draw the EDGEs on the surviving FACEs
   //- adjoining the EDGEs that are removed.  Input is the copied BODY, along 
   //- with the EDGEs being removed.

   CubitStatus draw_tweak_preview_tagged_FACEs( BODY *BODY_ptr, 
                                                CubitBoolean flush );
   //- Draw EDGEs on FACEs with the "tweak_preview" attribute on them.

   static AcisTweakTool* instance_;
   AcisQueryEngine *AQE;  // For convenience
   AcisModifyEngine *AME; // For convenience
};

#endif
