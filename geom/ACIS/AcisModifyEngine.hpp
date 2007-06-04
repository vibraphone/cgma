//-------------------------------------------------------------------------
// Filename      : AcisModifyEngine.hpp
//
// Purpose       : This class creates and modifies Acis-type geometry
//
// Special Notes :
//
// Creator       : Tim Tautges
//
// Creation Date : 2/01
//
// Owner         : Tim Tautges
//-------------------------------------------------------------------------

#ifndef ACIS_MODIFY_ENGINE_HPP
#define ACIS_MODIFY_ENGINE_HPP

#include "GeometryModifyEngine.hpp"
#include "AcisQueryEngine.hpp"

#include <set>

class ENTITY;
class ENTITY_LIST;
class BODY;
class WIRE;
class LUMP;
class SHELL;
class FACE;
class LOOP;
class COEDGE;
class EDGE;
class VERTEX;

class outcome;

class SPAposition;
class SPAtransf;
class SPAbox;

class SURFACE;
class surface;
class PLANE;
class CONE;
class SPHERE;
class TORUS;
class SPLINE;

class TopologyBridge;
class GeometryEntity;
class BodySM;
class Lump;
class ShellSM;
class Surface;
class LoopSM;
class Curve;
class Point;

class AcisBridge;
class BodyACIS;
class LumpACIS;
class ShellACIS;
class SurfaceACIS;
class LoopACIS;
class CoEdgeACIS;
class CurveACIS;
class PointACIS;
class BodySM;

class CubitBox;
class CubitString;
class GMem;
class ProgressTool;
template <class X> class DLIList;

class FeatureCut;
//#ifdef CUBIT_GUI
//class GUITestingDlg; // Testing only (remove when done)
//#endif

class AcisFacetManager;

class AcisModifyEngine : public GeometryModifyEngine
{
public:
// ********** BEGIN FRIEND DECLARATIONS        **********
  friend class AcisHealerTool;
  friend class AcisSurfaceTool;
  friend class AcisTweakTool;
  friend class AcisEdgeTool;
  friend class AcisToolUtil;
  friend class AcisBridge;
  friend class PointACIS;
  friend class CurveACIS;
  friend class CoEdgeACIS;
  friend class LoopACIS;
  friend class SurfaceACIS;
  friend class ShellACIS;
  friend class LumpACIS;
  friend class BodyACIS;
  friend class FeatureCut;
  friend class FeatureCutSurface;
  friend class FeatureExternalCLoop;
  friend class FeatureDecomp;
  friend class FeatureHLoop;
  friend class AcisQueryEngine;
//#ifdef CUBIT_GUI
//  friend class GUITestingDlg; // Temporary - for testing
//#endif

  friend class ATTRIB_CUBIT_OWNER;

// ********** END FRIEND DECLARATIONS        **********

//HEADER- Constructor and Destructor
  private:

  //AcisModifyEngine(AcisQueryEngine *aqe = NULL);
  AcisModifyEngine();

  static AcisModifyEngine* instance_;
  
  std::set<AcisBridge*> deactivatedSet;
  
  bool bridge_deactivated(AcisBridge*) const;
  CubitStatus deactivate_bridge(AcisBridge*) const;
  CubitStatus reactivate_bridge(AcisBridge*) const;
  void cleanout_deactivated_geometry() const;

  CubitStatus cleanup_slivers( BODY *body_to_cleanup );

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

//HEADER- RTTI and safe casting functions.

  virtual const type_info& entity_type_info() const ;
    //R- The geometric modeler type
    //- This function returns the type of the geometric modeler.

//HEADER- Functions for creating geometric primitives

  virtual BodySM* sphere(double radius) const ;
    //R Body*
    //R- A pointer to a newly created Body
    //I radius
    //I- Radius of the sphere
    //- Creates a sphere and assigns it to a Body.
    //- Returns pointer to the newly created body.

  virtual BodySM* brick( double wid, double dep, double hi ) const ;
    //R Body*
    //R- A pointer to a newly created Body
    //I wid
    //I- Width of the brick
    //I dep
    //I- Depth of the brick
    //I hi
    //I- Height of the brick
    //- Creates a cuboid and assigns it to a Body $
    //- Returns pointer to the newly created body.

  virtual BodySM* brick( const CubitVector &center, 
                         const CubitVector axes[3],
                         const CubitVector &extension) const ;
    //R Body*
    //R- A pointer to a newly created Body
    //I center
    //I- Center location of the brick
    //I axes
    //I- XYZ axes of brick
    //I extension
    //I- Size of brick, equivalent to 1/2 width, height, depth
    //- Creates a cuboid and assigns it to a Body.
    //- Returns pointer to the newly created body.

  virtual BodySM* prism( double height, int sides, double major,
                       double minor) const ;
    //- Creates an ACIS prism and assigns it to a Body $
    //- {height, major, minor} input height, major and minor radii. $
    //- {sides} input number of sides. Must be >= 3.
    //- Returns the ID of the new Body or CUBIT_FAILURE


  virtual BodySM* pyramid( double height, int sides, double major,
                         double minor, double top=0.0) const ;
    //- Creates an ACIS pyramid and assigns it to a Body $
    //- {height, major, minor} input height, major and minor radii. $
    //- {sides} input number of sides. Must be >= 3.
    //- {top} radius at top of pyramid.
    //- Returns the ID of the new Body or CUBIT_FAILURE

  virtual BodySM* cylinder( double hi, double r1, double r2,
                          double r3 ) const ;
    //- Creates an ACIS frustum and assigns it to a Body $
    //- {hi} input height $
    //- {r1} input radius in x-direction at base $
    //- {r2} input radius in y-direction at base $
    //- {r3} input radius in x-direction at top
    //- Returns the ID of the new Body or CUBIT_FAILURE

  virtual BodySM* torus( double r1, double r2 ) const ;
    //- Creates an ACIS torus and assigns it to a Body $
    //- {r1} input major_radius $
    //- {r2} input minor_radius
    //- Returns the ID of the new Body or CUBIT_FAILURE

  virtual BodySM* planar_sheet( const CubitVector& p1,
                                const CubitVector& p2,
                                const CubitVector& p3,
                                const CubitVector& p4 ) const;
    //- Creates a solid body consisting of a planar sheet (no volume)
    //- {p1} - 1st corner of the sheet
    //- {p2} - 2nd corner of the sheet
    //- {p3} - 3rd corner of the sheet
    //- {p4} - 4th corner of the sheet

//HEADER- Boolean operator functions.

  virtual CubitStatus subtract( DLIList<BodySM*> &tool_body_ptr,
                                DLIList<BodySM*> &from_bodies,
                                DLIList<BodySM*> &new_bodies,
                                bool imprint = false,
                                bool keep_old = false) const ;
    //- new form; implements new body topology checking and persistent objects

  virtual CubitStatus imprint( BodySM* BodyPtr1, 
                               BodySM* BodyPtr2,
                               BodySM*& newBody1, 
                               BodySM*& newBody2,
                               bool keep_old = false) const ;
  virtual CubitStatus imprint(DLIList<BodySM*> &from_body_list,
                              DLIList<BodySM*> &new_from_body_list,
                               bool keep_old = false,
                               DLIList<TopologyBridge*> *new_tbs = NULL,
                               DLIList<TopologyBridge*> *att_tbs = NULL) const;
    //- imprint a list of bodies, either with itself or with a with_list
    //- performs all imprints on acis bodies before constructing cubit
    //- bodies

  virtual CubitStatus imprint( DLIList<BodySM*> &body_list,
                               DLIList<Curve*> &ref_edge_list,
                               DLIList<BodySM*>& new_body_list,
                               bool keep_old_body,
                               bool show_messages = true) const;
    //- Imprints a list of Bodies with a list of RefEdges.  All entities must
    //- be ACIS entities.  Useful for splitting surfaces.  If edge pierces a
    //- surface a hardpoint will result at the pierce location.

  virtual CubitStatus imprint( DLIList<Surface*> &surface_list,
                               DLIList<Curve*> &curve_list,
                               DLIList<BodySM*>& new_body_list,
                               bool keep_old_body ) const;
    //- Imprints a list of RefFaces with a list of RefEdges.  This is
    //- useful if the user has a curve which spans several surfaces on
    //- a body and only wants to imprint to selected surfaces.  Algorithm
    //- does not support imprinting to free surfaces.

  virtual CubitStatus imprint( DLIList<Surface*> &surface_list,
                               DLIList<DLIList<Curve*>*> &curve_lists_list,
                               BodySM*& new_body,
                               bool keep_old_body ) const;
    //- Imprints a list of Surfaces with list of Curves, sorted per
    //- Surface (ie., curve_lists_list is same length as surface_list).
    //- This version is more efficient than the general-purpose one
    //- above, as we know which curves to imprint with which surfaces.
    //- All input surfaces must be from the same body.

  virtual CubitStatus imprint( DLIList<BodySM*> &body_list,
                               DLIList<CubitVector*> &vector_list,
                               DLIList<BodySM*>& new_body_list,
                               bool keep_old_body,
                               DLIList<TopologyBridge*> *new_tbs = NULL,
                               DLIList<TopologyBridge*> *att_tbs = NULL ) const;
    //- Imprints a list of bodies with a list of vectors.  Useful for
    //- splitting curves and creating hardpoints on surfaces.

  virtual CubitStatus imprint_projected_edges( 
                               DLIList<Surface*> &ref_face_list,
                               DLIList<Curve*> &ref_edge_list,
                               DLIList<BodySM*>& new_body_list,
                               bool keep_old_body,
                               bool keep_free_edges ) const;
    //- Imprints a list of RefFaces with a list of projected RefEdges.

  virtual CubitStatus imprint_projected_edges(
                               DLIList<Surface*> &ref_face_list,
                               DLIList<BodySM*> &body_list,
                               DLIList<Curve*> &ref_edge_list,
                               DLIList<BodySM*>& new_body_list,
                               bool keep_old_body,
                               bool keep_free_edges ) const;
    //- Imprints a list of Bodies with a list of RefEdges which are projected
    //- to a list of RefFaces

  virtual CubitStatus tolerant_imprint( DLIList<BodySM*> &bodies_in,
                                        DLIList<BodySM*> &new_bodies,
                                        DLIList<TopologyBridge*> *new_tbs,
                                        DLIList<TopologyBridge*> *att_tbs ) const;
  //Imprints a list of bodies with each other.  Can imprint successfully
  //and expectedly with sloppy/dirty models.  

  virtual CubitStatus project_edges( DLIList<Surface*> &ref_face_list,
                                     DLIList<Curve*> &ref_edge_list_in,
                                     DLIList<Curve*> &ref_edge_list_new,
                                     bool print_error = true ) const;
    //- Projects list RefEdges to a list of RefFaces

  virtual CubitStatus intersect( BodySM* tool_body_ptr,
                                 DLIList<BodySM*> &from_bodies,
                                 DLIList<BodySM*> &new_bodies,
                                 bool keep_old = false) const ;
    //R CubitStatus
    //R-the result of the Boolean operation: Success or Failure
    //I BodyPtr1, BodyPtr2
    //I-Two Body pointers that will be Booleaned
    //O newBody (or newBody1, newBody2 for imprint)
    //O-the new Body build by boolean operation on two Bodys.
    //O-for imprint, the output is two Bodies.
    //- These functions perform boolean operations on two Bodys and
    //- return the result through the output arguments. If the
    //- boolean operations fails at any stage, a NULL value is assigned
    //- to the output argument(s) and the function returns
    //- CUBIT_FAILURE. If everything goes well, the function returns
    //- CUBIT_SUCCESS. The original Bodys are left untouched during the
    //- boolean operation. The boolean operation is carried out on
    //- copies of the original Bodys.

  virtual CubitStatus unite( DLIList<BodySM*> &bodies, 
                             DLIList<BodySM*> &newBody,
                             bool keep_old = false) const ;
    //R CubitStatus
    //R-the result of the unite operation: Success or Failure
    //I bodies
    //I-DLIList<Body*>: a list of Body pointers that will be united
    //O newBody
    //O-the new Body build by unite operation on the list of  Body pointers.
    //- This function performs a unite of a list of Bodys and returns
    //- the result through the output argument, "newBody". If the unite
    //- operation went through OK, the function returns CUBIT_SUCCESS. If,
    //- for some reason, the unite operation did not go well, the output
    //- argument is assigned a NULL value and the function returns
    //- CUBIT_FAILURE. In either case, the original Bodys are left
    //- untouched.

 virtual CubitStatus chop( DLIList<BodySM*> &bodies, 
                           DLIList<BodySM*> &intersectBodies,
                           DLIList<BodySM*> &outsideBodies,  
                           BodySM*& leftoversBody,
	                         bool keep_old = false, 
                           bool nonreg = false) const ;
    //R CubitStatus
    //R-the result of the chop operation: Success or Failure
    //I bodies
    //I-DLIList<Body*>: a list of Body pointers that will be united
    //O intersectBody, outsideBody, leftoversBody
    //O- new Bodies build by chop operation on the list of  Body pointers.
    //- This function performs a chop of a tool body on a blank body  and
    //- returns the result through the output argument intersectBody,
    //- outsideBody, leftoversBody. If the chop operation went through OK,
    //- the function returns CUBIT_SUCCESS. If, for some reason, the chop
    //- operation did not go well, the output argument is assigned a NULL
    //- value and the function returns CUBIT_FAILURE. In either case, the
    //- original Bodys are left untouched.


 virtual CubitStatus thicken( DLIList<BodySM*>& bodies, 
                              DLIList<BodySM*>& new_bodies,
                              double depth, 
                              bool both = false) const;

    //- Thicken a sheet body into a solid
    //R CubitStatus
    //R-the result of the thicken operation: Success or Failure

 virtual CubitStatus flip_normals( DLIList<Surface*>& face_list ) const;
    //R CubitStatus
    //R-the result of the flip_normals operation: Success or Failure

//HEADER- Functions for geometric sweep operations

  virtual CubitStatus  sweep_translational(
                             DLIList<GeometryEntity*>& ref_ent_list,
                             DLIList<BodySM*>& result_body_list,
                             const CubitVector& sweep_vector,
                             double draft_angle,
                             int draft_type,
                             bool switchside,
                             bool rigid) const ;

  virtual CubitStatus  sweep_perpendicular(
                             DLIList<GeometryEntity*>& ref_ent_list,
                             DLIList<BodySM*>& result_body_list,
                             double distance,
                             double draft_angle,
                             int draft_type,
                             bool switchside,
                             bool rigid) const ;

  virtual CubitStatus  sweep_rotational(
                             DLIList<GeometryEntity*>& ref_ent_list,
                             DLIList<BodySM*>& result_body_list,
                             const CubitVector& point,
                             const CubitVector& direction,
                             double angle,
                             int steps = 0,
                             double draft_angle = 0.0,
                             int draft_type = 0,
                             bool switchside = false,
                             bool make_solid = false,
                             bool rigid = false ) const ;

  virtual CubitStatus sweep_along_curve(
                             DLIList<GeometryEntity*>& ref_ent_list,
                             DLIList<BodySM*>& result_body_list,
                             DLIList<Curve*>& ref_edge_list,
                             double draft_angle = 0.0,
                             int draft_type = 0,
                             bool rigid = false ) const;


//HEADER- Functions for webcut operations

  virtual CubitStatus webcut( 
                      DLIList<BodySM*>& webcut_body_list,
                      const CubitVector &v1,
                      const CubitVector &v2,
                      const CubitVector &v3,
                      DLIList<BodySM*>& results_list,
                      bool imprint = false ) const ;
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
                              DLIList<BodySM*>& results_list,
                              bool imprint  = false ) const ;
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

  virtual CubitStatus webcut_with_cylinder( 
                                    DLIList<BodySM*> &webcut_body_list,
                                    double radius,
                                    const CubitVector &axis,
                                    const CubitVector &center,
                                    DLIList<BodySM*>& results_list,
                                    bool imprint = false );
    //- cores the bodies with a cylinder formed by the given
    //- parameters and the height larger than the model.

  virtual CubitStatus webcut_with_brick( 
                                 DLIList<BodySM*>& webcut_body_list, 
                                 const CubitVector &center,
                                 const CubitVector axes[3], 
                                 const CubitVector &extension,
                                 DLIList<BodySM*> &results_list,
                                 bool imprint = false );
  /**<  Webcuts the bodies in the list with a cutting brick.
    *  The brick is created by the given parameters - center of
    *  brick, xyz axes, and extension.  Extension is 1/2 width,
    *  height and depth. Brick creation is done in the
    *  solid modeling engine to reduce the impact on body ids.
    */

  virtual CubitStatus webcut_with_planar_sheet( 
                                        DLIList<BodySM*>& webcut_body_list,
                                        const CubitVector &center,
                                        const CubitVector axes[2],
                                        double width, double height,
                                        DLIList<BodySM*> &results_list,
                                        bool imprint = false );
  /**<  Webcuts the bodies in the list with a cutting planar sheet.
    *  The sheet is created by the given parameters - center of
    *  sheet, xy axes, and width and height. Sheet creation is done
    *  in the solid modeling engine to reduce the impact on body ids.
    */

  virtual CubitStatus webcut_with_curve_loop( 
                                      DLIList<BodySM*> &webcut_body_list,
                                      DLIList<Curve*> &ref_edge_list,
                                      DLIList<BodySM*>& results_list,
                                      bool imprint = false);
    //- creates a sheet body with the given curve loop
    //- uses the new sheet body to cut the body list
 
  virtual CubitStatus webcut_with_sweep_surfaces_rotated(
                              DLIList<BodySM*> &blank_bodies,
                              DLIList<Surface*> &surfaces,
                              const CubitVector& point,
                              const CubitVector& sweep_vector,
                              double angle, 
                              Surface *stop_surf,
                              bool up_to_next, 
                              DLIList<BodySM*> &results_list,
                              CubitBoolean imprint = false );

  virtual CubitStatus webcut_with_sweep_curves_rotated(
                              DLIList<BodySM*> &blank_bodies,
                              DLIList<Curve*> &curves,
                              const CubitVector& point,
                              const CubitVector& sweep_vector,
                              double angle, 
                              Surface *stop_surf,
                              DLIList<BodySM*> &results_list,
                              CubitBoolean imprint = false );
  //-these 2 functions sweep a surface or curve about an axis, creating a swept
  //volume or surface respectively, which is in turn used to webcut the blank_bodies.
  //stop_surface is a surface where the sweep will terminate.  If stop_surface is
  //specified, there MUST be some intersection between sweep and stop_surface. 

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
                              CubitBoolean imprint = false );
  //-sweeps a surface(s) along a vector, curve or perpendicular to the surface 
  //(last case implied when sweep_vector's length is zero) producing a swept body.
  //If through_all is true, sweep will traverse blank_bodies completely.  If stop_surf
  //is specified, sweep terminates at stop_surface.  Operation fails is there is 
  //NO intersection of stop_surf with swept body. Swept body is
  //then used for webcutting blank_bodies. When multiple surfaces are specified,
  //an attempt is made to stitch them together into a single sheet body.  If this
  //fails, each surface is swept individually and the resultant bodies united.
  //This tool body is then used for the webcutting. 

  virtual CubitStatus webcut_with_sweep_curves(
                              DLIList<BodySM*> &blank_bodies,
                              DLIList<Curve*> &curves,
                              const CubitVector& sweep_vector,
                              bool through_all, 
                              Surface *stop_surf, 
                              Curve *curve_to_sweep_along, 
                              DLIList<BodySM*> &results_list,
                              CubitBoolean imprint = false );
  //-sweeps a curve(s) along a vector, curve. Swept surface(s)are stitched 
  //together and used as cutting tool to webcut blank_bodies. If surface(s) 
  //do not completely cut blank_bodies, webcut fails. If stop_surface is
  //specified, there MUST be some intersection between sweept surface and 
  //stop_surface. 

  CubitStatus trim_up_to_next_surface( BODY *&tool_body, DLIList<BodySM*> &blank_bodies,
                                       VERTEX *surf_VERT );

  virtual CubitStatus section( DLIList<BodySM*> &section_body_list,
                               const CubitVector &point_1,
                               const CubitVector &point_2,
                               const CubitVector &point_3,
                               DLIList<BodySM*>& new_body_list,
                               bool keep_normal_side,
                               bool keep_old = false,
                               bool keep_both_sides = false);
    //- Section will cut a list a bodies and keep a side of the bodies.
    //- The bodies are cut with a planar surface (surface will be extended).

//HEADER- Functions that create GeometryEntities

  virtual BodySM* copy_body ( BodySM* bodyPtr) const ;
    //R Body*
    //R- A pointer to the newly created body
    //I bodyPtr
    //I- A pointer to the body to be copied
    //- This function makes a copy of the input Body and returns a
    //- pointer to the newly created copy. The input Body and the newly
    //- created Body are geometrically identical.

  virtual Point* make_Point( CubitVector const& point) const ;
    //R Point*
    //R- Returned pointer to a Point object.
    //I point_type
    //I- The type of point to be created. At this time, this argument
    //I- is not used and should be set to UNDEFINED_POINT_TYPE.
    //I point
    //I- Input coordinates of the point to be created.
    //- This function creates a Point, given coordinates.  The particular
    //- type of Point object that is created depends on the specific
    //- modeling engine.  For example, if the engine
    //- is AcisModifyEngine, then a PointACIS is created and returned.

  virtual Curve* make_Curve(Curve *curve_ptr) const;
    //- creates a curve from an existing curve.  It also returns
    //- the start and end points.  This creates totally new topology
    //- in ACIS.  This function is useful for constructing geometry from
    //- existing geometry.

  virtual Curve* make_Curve( Point const* point1_ptr,
                             Point const* point2_ptr,
                             Surface* ref_face_ptr,
                             const CubitVector *third_point = NULL ) const;
    //- Create a curve exactly on the give ref_face.
    //- The third point is optional so it can be NULL.  It is
    //- or can be used for periodic curves that result.

  virtual Curve* make_Curve( GeometryType curve_type,
                             Point const* point1_ptr,
                             Point const* point2_ptr,
                             DLIList<CubitVector*>& vector_list,
                             Surface* ref_face_ptr = NULL ) const;
  virtual Curve* make_Curve( GeometryType curve_type,
                             Point const* point1_ptr,
                             Point const* point2_ptr,
                             CubitVector const* intermediate_point_ptr,
                             CubitSense sense) const;
    //R Curve*
    //R- Returned pointer to a Curve object.
    //I curve_type
    //I- The type of curve to be created.
    //I point1_ptr, point2_ptr
    //I- Input end points of the curve to be created.
    //I intermediate_point_ptr
    //I- The coordinates of an intermediate point required for the
    //I- generation of the curve.
    //I vector_list
    //I- Input list of CubitVectors. The new Curve interpolates these
    //I- positions.
    //I ref_face_ptr
    //I- If the optional RefFace pointer is provided, then the input
    //I- locations (CubitVectors) are moved to the surface of the
    //I- RefFace before interpolation is done.
    //I- NOTE: The end Points are *not* moved to the surface. Only the
    //I-       (intermediate) points are.
    //- This function creates a Curve, of type, curve_type, given
    //- the end points of the curve and an intermediate point or a set
    //- of points to interpolate through.
    //- The particular type of Curve object
    //- that is created depends on the specific modeling engine.  For
    //- example, if the engine is AcisModifyEngine, then a CurveACIS
    //- is created and returned.
    //-
    //- NOTE on some hideous programming
    //-   This function not only creates the Curve, but also
    //-   creates the VGI entities Chain and CoVertexes that
    //-   will be connected up with the RefEdge. For the sake
    //-   of expediency we are assuming that this function is
    //-   called when RefEdge needs to be made from RefVertices.
    //-
    //- This function creates non-linear curves (e.g., elliptic and
    //- parabolic curves. curve_type indicates what type of curve is
    //- required.
    //-
    //- Elliptical:
    //- The portion of the ellipse (circle) that is generated goes from
    //- point1 to point2.
    //-
    //- Parabolic:
    //- Construct a parabolic arc from 3 points. intermediate_point is the
    //- peak of the parabola - in this case, the point which is equidistant
    //- from the start and end points of the parabola. The 3 points must form
    //- an isosceles triangle. This definition limits the user to generation
    //- of the tip of parabolic shapes only.

  virtual Surface* make_Surface( Surface *old_surface_ptr,
                                 bool extended_from ) const;
    //R Surface*
    //R- Pointer to a newly created Surface object.
    //I Surface*
    //I- The surface from which we want to create a new one.
    //- This function creates a new surface from an existing one.
    //- The new surface is attached to ACIS geometry.  The acis
    //- goemetry is attached to a full data structure, loop, lump, bodies..

  virtual Surface* make_Surface( GeometryType surface_type,
                                 DLIList<Curve*>& curve_list,
                                 Surface *old_surface_ptr = NULL,
                                 bool check_edges = true ) const; 
    //R Surface*
    //R- Pointer to a newly created Surface object.
    //I surface_type
    //I- The type of surface to be created.
    //I curve_list
    //I- list of curves to be used as the bounds of the Surface
    //- This function creates a Surface, of type, surface_type, given
    //- the list of curves.
    //- The particular type of Surface object created depends on the
    //- specific modeling engine.  For example, if the engine is
    //- AcisModifyEngine, then a SurfaceACIS is created.
    //-
    //- NOTE on some hideous programming
    //-   This function not only creates the Surface, but also
    //-   creates the VGI entities Loop and CoEdges that
    //-   will be connected up with the RefFace. For the sake
    //-   of expediency we are assuming that this function is
    //-   called when RefFace needs to be made from RefEdges.
    //-

  virtual Lump* make_Lump( DLIList<Surface*>& surface_list ) const ;
    //R Lump*
    //R- Pointer to a newly created Lump object.
    //I lump_type
    //I- The type of surface to be created.
    //I surface_list
    //I- list of surfaces to be used as the bounds of the Lump
    //- This function creates a Lump, of type, lump_type, given
    //- the list of surfaces.
    //- The particular type of Lump object created depends on the
    //- specific modeling engine.  For example, if the engine is
    //- AcisModifyEngine, then a LumpACIS is created.
    //-
    //- NOTE on some hideous programming
    //-   This function not only creates the Lump, but also
    //-   creates the VGI entities Shell and CoFaces that
    //-   will be connected up with the RefVolume. For the sake
    //-   of expediency we are assuming that this function is
    //-   called when RefVolume needs to be made from RefFaces.
    //-

    virtual BodySM* make_BodySM( Surface * ) const;
    //- given a Surface, make a BodySM

    virtual BodySM* make_BodySM( DLIList<Lump*>& lump_list ) const ;
    //R BodySM*
    //R- Pointer to a newly created Body object.
    //I lump_list
    //I- list of lumps to be used to create the BodySM
    //- This function creates a BodySM given the list of lumps.
    //- The particular type of BodySM object created depends on the
    //- specific modeling engine.  For example, if the engine is
    //- AcisModifyEngine, then a BodyACIS is created. Non-solid
    //- model based engines can use the default implementation
    //- to return a NULL pointer.
    //-
    //- NOTE on some hideous programming
    //-   This function not only creates the Lump, but also
    //-   creates the VGI entities CoVolumes that
    //-   will be connected up with the Body. For the sake
    //-   of expediency we are assuming that this function is
    //-   called when Body needs to be made from RefVolumes.
    //-

//HEADER- Miscellaneous functions

  virtual CubitStatus webcut_with_sheet(DLIList<BodySM*>& webcut_body_list,
                                        BodySM *sheet_body,
                                        DLIList<BodySM*> &new_bodies,
                                        bool imprint = false);
    //- webcuts a body using a sheet body.
    //- It splits the sheet into two single sided bodies.
    //- it then subtracts this with the webcut body.
    //- The result is splitting the webcut_body into halves.
    //- if the webcut body is a topological torus, this routine
    //- will fail...

  virtual CubitStatus webcut_with_sheet(BODY *webcut_body,
                                        BODY *sheet_body,
                                        BODY *&webcut_body_1,
                                        BODY *&webcut_body_2,
                                        CubitBoolean imprint  = CUBIT_FALSE);
    //- webcuts a body using a sheet body.
    //- It splits the sheet into two single sided bodies.

  virtual CubitStatus webcut_with_extended_surf(DLIList<BodySM*> &webcut_body,
                                                Surface *extend_from,
                                                DLIList<BodySM*> &new_bodies,
                                                int &num_cut,
                                                bool imprint = false);
    //- creates a surface by extending one from the given surface then
    //- webcuts using the this sheet.(see webcut_with_sheet.

  virtual CubitStatus split_body( BodySM *body_ptr,
                                  DLIList<BodySM*> &new_bodies );
    //- splits a body into several bodies (multiple volumes in
    //- a single body go to multiple bodies.

  virtual CubitStatus split_periodic( BodySM *body_ptr,
                                      BodySM *&new_body);
    //- splits the periodic bodies

  virtual CubitStatus reverse_body( BodySM *body_to_reverse );
    //- Reverse body (turn it inside-out).
  
    //- Reflects the body.  This changes acis topology so we need
    //- to create a new body.

  virtual CubitStatus create_solid_bodies_from_surfs(
                                    DLIList<Surface*> &ref_face_list,
                                    DLIList<BodySM*> &new_bodies,
                                    bool keep_old = false,
                                    bool heal = true) const;
  //- This function assumes that the reffaces sent into
  //- this function are either sheet bodies, or free surfaces.  This
  //- Will have been taken care of in the calling function.  GT?
  //- All the surfaces, in ACIS, will be turned into sheet bodies.  They
  //- will be united, then the void enclosed.  The result will be healed
  //- if the option is specified (default is to heal, from experience).
  //- This function will not create new reffaces.

#ifdef BOYD14
  CubitStatus combine_surfaces( DLIList<Surface*>& face_list, 
                                BodySM*& new_body, 
                                bool keep_old = false) const;
  CubitStatus combine_sheets( DLIList<BodySM*>& body_list, 
                              BodySM*& new_body, 
                              bool keep_old = false) const ;
#endif

  virtual CubitStatus webcut_across_translate( DLIList<BodySM*>& body_list,
                                               Surface* plane_surf1,
                                               Surface* plane_surf2,
                                               DLIList<BodySM*>& results_list,
                                               bool imprint = false) const;
  // In-process function to webcut a flat plate suitable for singe-single sweeping.

  CubitStatus offset_curves( DLIList<Curve*>& ref_edge_list, 
                             DLIList<Curve*>& new_curves,
                             double offset_distance,
                             const CubitVector& offset_direction, 
                             int gap_type = 1 );
  //- Creates curves offset from a chain of curves.  The offset direction is
  //- only used if there is one linear curve.  Otherwise, the offset direction
  //- is calculated by ACIS (the cross product of the wires tangent and the
  //- planar normal).  The gap type is 0 - rounded, 1 - extended, 2 - natural.

  Curve* trim_curve( Curve* trim_curve, 
                     const CubitVector& trim_vector,
                     const CubitVector& keep_vector,
                     bool keep_old = false );
    //- Trims or extends a curve, up to the trim_vector.  If trimming, the
    //- keep_vector determines which side of the curve to keep.  If the curve
    //- is not free, the curve is automatically copied before trimming (so
    //- a new curve results).

  Curve* create_arc_three( Point* ref_vertex1, 
                           Point* ref_vertex2,
                           Point* ref_vertex3, 
                           bool full = false );
  Curve* create_arc_three( Curve* ref_edge1,
                           Curve* ref_edge2,
                           Curve* ref_edge3, 
                           bool full = false );
  Curve* create_arc_center_edge( Point* ref_vertex1, 
                                 Point* ref_vertex2,
                                 Point* ref_vertex3, 
                                 const CubitVector &normal,
                                 double radius = CUBIT_DBL_MAX,
                                 bool full = false );
    //- Methods to create arcs.  First uses 3 points on arc, center creates arc
    //- tangent to 3 curves, third creates arc using center and two points on arc.

  CubitStatus create_curve_combine( DLIList<Curve*>& curve_list, 
                                    Curve *&new_curve_ptr );
    //-  Uses the solid modeller to create a new RefEdge that is a combination 
    //- of the input chain of edges.  
    //-

  GeometryQueryEngine *get_gqe();
    //- return the query engine

  AcisQueryEngine *get_acis_query_engine()  {return AcisQueryEngine::instance();}
    //- return the query engine

  virtual CubitBoolean is_modify_engine(const TopologyBridge *tb_ptr) const;
    /**< return CUBIT_TRUE if the tb_ptr belongs to this modify engine
     */

  virtual CubitStatus get_offset_intersections( 
                                      Curve* ref_edge1, 
                                      Curve* ref_edge2,
                                      DLIList<CubitVector*>& intersection_list,
                                      double offset,
                                      bool ext_first = true );
    //- Finds the intersections of a certain distance (offset) between two
    //- curves.  The two curves must lie in a plane.  The first curve is offset
    //- the offset distance in both directions, and the bounded intersections with
    //- the second curve are found.  The first curve and optionally be extended
    //- to infinity for the intersection calculation.  The intent of the function
    //- is so that the user can create a point on a curve a certain distance
    //- from another curve, as in specifying a reference location for a gage
    //- diameter on an arc in an engineering drawing.

  virtual CubitStatus get_offset_intersections( 
                                      Curve* ref_edge_ptr, 
                                      Surface* ref_face_ptr,
                                      DLIList<CubitVector*> &intersection_list,
                                      double offset = 0.0,
                                      bool ext_surf = true );
    //- Finds intersections (points) of the curve and surface.  The surface can
    //- be offset - it is offset to each side and intersections are found.  By
    //- default the surface is extended to infinity (if possible) and the
    //- intersections are found.  The function allocates the CubitVectors in
    //- the returned list, so be sure to free them.

  virtual CubitStatus surface_intersection( Surface *surface1_ptr,
                                            Surface *surface2_ptr,
                                            DLIList<Curve*> &inter_graph,
                                            const double tol) const;
    //- Finds the intersection of two bodies defined by the curves of the
    //- intersection graph.

  virtual CubitStatus get_mid_plane( const CubitVector &point_1,
			                               const CubitVector &point_2,
                                     const CubitVector &point_3,
			                               BodySM *body_to_trim_to,
			                               BodySM *&midplane_body ) const;
  //- Given 3 points that describe the mid_plane, return a body that
  //- contains the trimmed surfaces defining the mid plane.
  //- These surfaces will be part of the midplane_body and will
  //- need to be collected out by the calling code. 

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

  virtual CubitStatus remove_curve_slivers( BodySM *body_sm, double lengthlimit ) const;
  
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

  CubitStatus scale( BodySM *&body, const CubitVector& f );

protected:

private:

//HEADER- Functions that create TopologyEntities


  CubitStatus get_new_Body(DLIList<BodySM*> &old_Body_list,
                           DLIList<BODY*> &old_BODY_list,
                           DLIList<BODY*> &new_BODY_list,
                           DLIList<BodySM*> &new_Body_list,
                           const bool keep_old,
                           const bool delete_old = true) const;

  CubitStatus get_new_Body(DLIList<TopologyBridge*> &old_entity_list,
                           DLIList<BODY*> &new_BODY_list,
                           DLIList<BodySM*> &new_Body_list,
                           const bool keep_old,
                           const bool delete_old = true) const;

  CubitStatus get_new_Body(BodySM *old_Body,
                           BODY *old_BODY,
                           DLIList<BODY*> &new_BODY_list,
                           DLIList<BodySM*> &new_Body_list,
                           const bool keep_old,
                           const bool delete_old = true) const;

  BodySM *get_new_Body(DLIList<BodySM*> &old_Body_list,
                     DLIList<BODY*> &old_BODY_list,
                     BODY *new_BODY,
                     const bool keep_old,
                     const bool delete_old = true) const;

  BodySM *get_new_Body(BodySM *old_Body,
                     BODY *old_BODY,
                     BODY *new_BODY,
                     const bool keep_old,
                     const bool topology_check = false,
                     const bool delete_old = true) const;
    //- function used to create new bodies resulting from booleans; old and
    //- new BODY's are input to this function, along with old Body; new
    //- Body is returned.

  int mark_owners_deactivated_flag(ENTITY *this_ENTITY, bool flag,
                                   bool recurse = CUBIT_TRUE) const;
    //- mark the owners of BODY and all its children with deactivated flag input

  int mark_owners_deactivated_flag(ENTITY_LIST& ENTITIES, bool flag) const;
  
  void print_deleted_reused(DLIList<BODY*> &old_BODY_list,
                            DLIList<BODY*> &new_BODY_list) const;
    //- debugging function for get_new_Body

//HEADER- Functions that create ACIS ENTITYs

  VERTEX* make_VERTEX( CubitVector const& point ) const;
    //R VERTEX*
    //R- Returned VERTEX pointer
    //I VERTEX_type
    //I- The type of VERTEX to be created (not used at this time)
    //I point
    //I- Input (global) coordinate location
    //- This function creates an ACIS VERTEX at the input location.
    //-
    //- All other associated ACIS entities are also created.

  EDGE* make_EDGE( GeometryType EDGE_type,
                   VERTEX* start_VERTEX,
                   VERTEX* end_VERTEX,
                   SPAposition* intermediate_point_ptr,
                   bool sense ) const;
    //R EDGE*
    //R- Returned EDGE pointer
    //I EDGE_type
    //I- The type of geometry entity to be associated with the EDGE
    //I- that is created
    //I start/end_VERTEX
    //I- The start and end VERTEX'es of the EDGE to the created.
    //I intermediate_point_ptr
    //I- This is the intermediate point required to complete the definition
    //I- of the underlying curve. It is not required for creating a
    //I- straight EDGE.
    //I sense
    //I- The "sense" makes sense (no pun intended) only for ellipses in
    //I- the XY plane. For ellipses in the XY plane, the curve can
    //I- defined from point-1 to point-2 in two ways. The direction of
    //I- the curve can be such that the vector going from the center to
    //I- the first point can be rotated in a counter clockwise manner to
    //I- reach the vector going from the center to the second point. Such
    //I- an ellipse is considered to have positive or TRUE sense. When
    //I- the curve is defined in the other direction it is considered to
    //I- have a negative or FALSE sense.
    //I- This argument is required only for some types of quadratic EDGEs
    //I- (such as ellipses and log spirals).
    //- This function creates an ACIS EDGE using start and end VERTEX'es
    //- and an intermediate point.
    //-
    //- Straight, parabolic and elliptical EDGEs are made in this function.
    //-
    //- Elliptical EDGE:
    //- The portion of the ellipse (circle) that is generated goes from
    //- start_VERTEX to end_VERTEX. The sense value determines which part
    //- of the whole ellipse is required. intermediate_point is the
    //- center of the ellipse.  If the radius ratio is not very close
    //- to 1.0 (a circle) then a logarithmic spiral is generated.
    //-
    //- Parabolic EDGE:
    //- Construct a parabolic arc from 3 points. intermediate_point is the
    //- peak of the parabola - in this case, the point which is equidistant
    //- from the start and end points of the parabola. The 3 points must form
    //- an isosceles triangle. This definition limits the user to generation
    //- of the tip of parabolic shapes only.
    //-
    //- All other associated ACIS entities are also created.

  EDGE* make_straight_EDGE( VERTEX* start_VERTEX,
                            VERTEX* end_VERTEX ) const;
    //R EDGE*
    //R- Returned EDGE pointer
    //I start/end_VERTEX
    //I- The start and end VERTEX'es of the EDGE to the created.
    //- This function creates a straight ACIS EDGE using start and
    //- end VERTEX'es.
    //-
    //- All other associated ACIS entities are also created.

  EDGE* make_parabolic_EDGE( VERTEX* from,
                             VERTEX* to,
                             SPAposition* intermediate_point_ptr ) const;
    //R EDGE*
    //R- Returned EDGE pointer
    //I start/end_VERTEX
    //I- The start and end VERTEX'es of the EDGE to the created.
    //I intermediate_point_ptr
    //I- This is the intermediate point required to complete the definition
    //I- of the underlying curve.
    //- This function creates a parabolic ACIS EDGE using start and
    //- end VERTEX'es and the intermediate_point. intermediate_point is the
    //- peak of the parabola - in this case, the point which is equidistant
    //- from the start and end points of the parabola. The 3 points must form
    //- an isosceles triangle. This definition limits the user to generation
    //- of the tip of parabolic shapes only.
    //-
    //- All other associated ACIS entities are also created.

  EDGE* make_elliptical_EDGE( VERTEX* from,
                              VERTEX* to,
                              SPAposition* ctr_ptr,
                              bool sense ) const;
    //R EDGE*
    //R- Returned EDGE pointer
    //I start/end_VERTEX
    //I- The start and end VERTEX'es of the EDGE to the created.
    //I center_point_ptr
    //I- This is the center point of the ellipse.
    //I sense
    //I- The "sense" makes sense (no pun intended) only for ellipses in
    //I- the XY plane. For ellipses in the XY plane, the curve can
    //I- defined from point-1 to point-2 in two ways. The direction of
    //I- the curve can be such that the vector going from the center to
    //I- the first point can be rotated in a counter clockwise manner to
    //I- reach the vector going from the center to the second point. Such
    //I- an ellipse is considered to have positive or TRUE sense. When
    //I- the curve is defined in the other direction it is considered to
    //I- have a negative or FALSE sense.
    //- This function creates an elliptical ACIS EDGE using start and
    //- end VERTEX'es and the center_point. center_point is the center
    //- of the ellipse.  The sense is with respect to the positive Z axis.
    //-
    //- All other associated ACIS entities are also created.

  EDGE* make_circular_EDGE( VERTEX* from,
                            VERTEX* to,
                            SPAposition* third_pt ) const;
    //R EDGE*
    //R- Returned EDGE pointer
    //I start/end_VERTEX
    //I- The start and end VERTEX'es of the EDGE to the created.
    //I third_pt
    //I- This is a third point on the circle.
    //-
    //- All other associated ACIS entities are also created.

  EDGE* make_spiral_EDGE( VERTEX* from,
                          VERTEX* to,
                          SPAposition* ctr_ptr,
                          bool sense) const;

  EDGE* make_surface_EDGE( VERTEX* from,
                           VERTEX* to,
                           FACE *FACE_ptr,
                           const CubitVector &plane_normal,
                           const CubitVector *third_point = NULL) const;

  EDGE* make_spline_EDGE( VERTEX* from,
                          VERTEX* to,
                          SPAposition* pos_array,
                          int number_points)  const;
    //R EDGE*
    //R- Returned EDGE pointer (uses api_mk_ed_cubic)
    //I from, to
    //I- The end VERTEX'es that must be used to create the EDGE.
    //I pos_array
    //I- Input array of ACIS position's used to construct the spline
    //I- (this must include the locations of the end VERTEX'es as the
    //I- first and last entries in the array).
    //I number_points
    //I- The number of entries (points) in the pos_array
    //- These functions create an ACIS EDGE who's geometry is a spline curve.
    //- The spline interpolates the input list of positions.
    //- After the EDGE is created using the input points, the new EDGE's
    //- VERTEX'es are replaced by the input VERTEX'es.

  EDGE* make_spline_EDGE( DLIList<CubitVector*> &vec_list );
    //- Create a spline curve from the input locations. Returns NULL if
    //- unable to create a spline.  Uses api_curve_spline

  FACE* make_type_FACE( CONE *CONE_ptr ) const;
    //R FACE*
    //R- face pointer connected to sheet body.
    //I- Cone geometry where info is to be gleaned.
    //-  This function takes the input cone geometry and creates
    //-  another face using this to create a larger face that
    //-  has a height fitting a factor larger than the model.

  FACE* make_type_FACE( SPHERE *SPHERE_ptr ) const;
    //R FACE*
    //R- face pointer connected to sheet body.
    //I- Sphere geometry where info is to be gleaned.
    //-  This function takes the input sphere geometry and creates
    //-  another face using this to create a full sphere face.
    //-  it is just the shell of the sphere.

  FACE* make_type_FACE( TORUS *TORUS_ptr ) const;
    //R FACE*
    //R- face pointer connected to sheet body.
    //I- Torus geometry where info is to be gleaned.
    //-  This function takes the input torus geometry and creates
    //-  another face using this to create a full torus face.
    //-  it is just the shell of the torus.

  FACE* make_type_FACE( SPLINE *SPLINE_ptr ) const;
    //R FACE*
    //R- face pointer connected to sheet body.
    //I- Spline geometry where info is to be gleaned.
    //-  This function takes the input spline geometry and creates
    //-  another face using this to create a full spline face.
    //-  it is just the shell of the spline.

  FACE* make_type_FACE( PLANE *PLANE_ptr ) const;
    //R FACE*
    //R- face pointer connected to sheet body.
    //I- Plane geometry where info is to be gleaned.
    //- The super bounding box of the model is used to creat the limits of the
    //- new plane.

  FACE* make_FACE( FACE* FACE_ptr, 
                   CubitBoolean extended_from = CUBIT_FALSE
                 ) const;
    //R FACE*
    //R- REturned FACE pointer
    //I- Existing FACE to make a copy from.
    //I- flag to tell whether to extend the face to its full geometric form or not.
    //- This function accepts apointer to a face and makes an api
    //- call to copy this face and construct a sheet body from it.

  FACE* make_FACE( GeometryType FACE_type,
                   DLIList<EDGE*>& EDGE_list,
                   Surface *old_surface_ptr = NULL) const;
    //R FACE*
    //R- Returned FACE pointer
    //I FACE_type
    //I- The type of geometry entity to be associated with the FACE
    //I- that is created
    //I EDGE_list
    //I- Input list of EDGEs that bound the FACE to be made.
    //- This function creates an ACIS FACE using the EDGEs in EDGE_list.
    //- The type of the underlying geometry is determined by the
    //- FACE_type argument.
    //- All other associated ACIS entities are also created.

  BODY* make_planar_quad_BODY ( const CubitVector& Point1,
                                const CubitVector& Point2,
                                const CubitVector& Point3,
                                const CubitVector& Point4 ) const;
    //- Creates an ACIS planar quadrilateral BODY $
    //- {Vertex1...4} the 4 input CubitVertex objects. $
    //- The ACIS FACE is constructed using the vertices in the input order. $
    //- Returns a BODY* == NULL if the 4 input vertices are not planar or $
    //- if the ACIS BODY could not be constructed for any other reason.

  BODY* make_brick_BODY( const CubitVector &center, 
                         const CubitVector axes[3],
                         const CubitVector &extension ) const;
    //- Makes a cuboid at the given location, orientation and size.  The
    //- extension vector is equivalent to 1/2 width, depth, and height.

  EDGE *copy_and_move_EDGE( EDGE *EDGE_ptr, 
                            const CubitVector &dir, 
                            double dist );
  //- Copy an EDGE and move it in the direction given a certain distance.
  //- NULL returned if fails.

  CubitStatus transform_Curve( Curve *curve_ptr, const SPAtransf& tran_vec );
    //-Apply the given transformation to the given Curve.

  CubitStatus offset_EDGES( DLIList<EDGE*> &EDGE_list, 
                            double offset_distance,
                            int gap_type, 
                            DLIList<EDGE*> &offset_EDGE_list );
  //- Offsets the connected chain of EDGES, creating a new list of ACIS EDGES.
  //- The gap type is as follows; 0 = rounded like arcs, 1 = extended like
  //- lines, and 2 = natural like curve extensions. Lines are extend with
  //- lines and circles are extended with circles.

//HEADER- Functions related to webcut operations

  CubitStatus webcut( DLIList<BodySM*>& webcut_body_list,
                      BODY* tool_BODY,
                      DLIList<BodySM*>& results_list,
                      bool imprint = false ) const ;
    //R int
    //R- Number of bodies that were webcut ( >= 0 )
    //I webcut_body_list
    //I- The list of bodies to be webcut
    //I tool_BODY
    //I- The BODY to be used for webcutting.
    //I merge
    //I- A flag to decide whether the new bodies created by the
    //I- webcutting process should be merged or not.
    //I imprint
    //I- A flag to decide whether the new bodies created by the
    //I- webcutting process should be imprinted or not.
    //- This functions webcuts a list of bodies using another body
    //- as the webcutting tool. The newly created bodies are
    //- merged and imprinted depeding on the respective flags.

  CubitStatus FACE_intersection( FACE *FACE1_ptr, FACE* FACE2_ptr,
                                 BODY *&intersect_graph,
                                 bool &interfering) const;
    //- Finds the intersection graph between two ACIS FACE's and returns
    //- the intersection graph.  If the faces do not
    //- intersect, the flag, interfering, is set to CUBIT_FALSE.

  CubitStatus convert_WIRE( BODY *wire_BODY,
                            DLIList <Curve*> &curve_list ) const;
    //- Converts the wire body into a list of curves.  Copies the edges in
    //- the wire into new edges and creates CurveACIS entities from the copies.

  CubitStatus webcut_BODY( BODY* inputAcisBody,
                           BODY* infiniteCuttingTool,
                           BODY*& webcutBody1,
                           BODY*& webcutBody2) const;
    //- Creates 2 ACIS BODYs that are the result of "webcutting" the input ACIS
    //- BODY, inputAcisBody, using the "infinite" solid BODY, infiniteCuttingTool.
    //- The ACIS FACES that overlap (at the cutting tool) are not "merged" -- i.e.,
    //- the resulting BODYs are manifold, possibly with multiple LUMPs some of which
    //- have FACEs that overlap in space. The Cellular Topology Husk is *not* used.
    //- NOTE: After the operation (whether it was successful or not), the input
    //-       BODYs are left intact.

#ifdef BOYD14
  void webcut_failed ( DLIList<BodySM*>& body_list,
                       ENTITY_LIST& BODY_list) const ;
    //- Outputs an error message, clears the echo stream, and calls webcut_cleanup
    //- to delete Bodies and BODYs, as the webcut command has failed
#endif

  void webcut_cleanup ( DLIList<BodySM*>& body_list,
                        ENTITY_LIST& BODY_list) const ;
    //- Deletes the Bodies in body_list and the Acis BODYs in BODY_list


  CubitStatus stitch_FACEs( DLIList<FACE*> &faces,
                            FACE *&result_FACE,
                            BODY *&stitched_BODY );
  //I List of FACEs you want stitched together
  //IO Pointer to FACE...gets set if all FACEs can be stitched into single FACE 
  //IO If stitching of FACEs is successful, stitched_BODY is resultant BODY.  
  //   Means that all faces could be stitched together, forming a single BODY,
  //   but not a single FACE.
  //returns CUBIT_SUCCESS if stitching was successful, other return FAILURE

  void webcut_imprint(BODY* cutting_tool_BODY_ptr,
                      DLIList<BodySM*> &old_body_list,
                      DLIList<BODY*>& new_webcut_BODY_list,
                      DLIList<BODY*>& just_webcut_BODY_list,
                      DLIList<BodySM*>& results_list,
                      DLIList<BodySM*> &imprinted_model_Body_list,
                      DLIList<BODY*> &imprinted_model_BODY_list) const ;
    //I cutting_tool_BODY_ptr
    //I- Pointer to the Cutting Tool BODY. This is required in the case
    //I- where an existing Body - is used as the CT. In that case, its BODY
    //I- is kept out of the imprint operations.
    //IO new_webcut_BODY_list
    //IO- List of BODYs created by webcutting the list of RefBodies
    //IO- specified in the webcut command. New BODYs created as a product
    //IO- of performing imprints on the existing BODYs in the Model are
    //IO- appended to this list.
    //- Performs the required imprint operations to ensure that the topologies
    //- of the boundaries of the new webcut BODYs "match" those of their
    //- neighbouring BODYs (for mesh compatibility reasons).

  BODY* create_infinite_plane_cutting_tool( const CubitVector &vecVertex1,
                                            const CubitVector &vecVertex2,
                                            const CubitVector &vecVertex3,
                                            const SPAbox& super_box,
                                            bool just_face = false ) const;

  BODY* copy_BODY (BODY* body_ptr,
                   bool remove_attribs = true ) const;
    //- Make a copy of the input BODY, body_ptr, and return a pointer
    //- to the new BODY. If the operation fails, a NULL pointer is returned.
    //- If remove_attribs is set to CUBIT_TRUE, then the attributes from the
    //- original BODY will be stripped off the new BODY (these are automatically
    //- copied to the new BODY, by ACIS.

//HEADER- Functions related to geometric sweep operations

  CubitStatus sweep_FACE_along_vector( FACE *& FACE_ptr,
                                       const CubitVector& sweep_vector,
                                       bool& volume_is_negative,
                                       bool primary_side = true,
                                       double draft_angle = 0.0,
                                       int draft_type = 0,
                                       bool rigid = false,
                                       FACE *to_FACE= NULL ) const;
    //- Sweeps the FACE, FACE_ptr, along the sweep_vector, a distance equal
    //- to the length of the vector, using the input draft_angle and draft_type.
    //- The primary side of the FACE is swept unless it is DOUBLE_SIDED
    //- and primary_side is CUBIT_FALSE -- in that case, the face is swept
    //- in the direction of the secondary normal, not the primary normal.
    //- The draft_type is either 0 or 1 (see the ACIS documentation on the
    //- api_sw_face_vec function) and indicates what should be done to "close"
    //- the end FACEs of the sweep operation, if the draft_angle is non-zero.
    //- The draft angle is in radians.
    //- If the sweep failed and if this was due to the creation
    //- of a solid of negative volume, then volume_is_negative is returned as
    //- CUBIT_TRUE.
    //- Returns CUBIT_SUCCESS if sweep successful, CUBIT_FAILURE otherwise.
  CubitStatus sweep_FACE_perpendicular( FACE *& FACE_ptr,
                                        double distance,
                                        bool & volume_is_negative,
                                        bool primary_side = true,
                                        double draft_angle = 0.0,
                                        int draft_type = 0,
                                        bool rigid = false ,
                                        FACE *to_FACE= NULL ) const;
    //- Sweeps the FACE, FACE_ptr, perpendicularly, a distance equal
    //- to the distance given, using the input draft_angle and draft_type.
    //- The primary side of the FACE is swept unless it is DOUBLE_SIDED
    //- and primary_side is CUBIT_FALSE -- in that case, the face is swept
    //- in the direction of the secondary normal, not the primary normal.
    //- The draft_type is either 0 or 1 (see the ACIS documentation on the
    //- api_sw_face_vec function) and indicates what should be done to "close"
    //- the end FACEs of the sweep operation, if the draft_angle is non-zero.
    //- The draft angle is in radians.
    //- If the sweep failed and if this was due to the creation
    //- of a solid of negative volume, then volume_is_negative is returned as
    //- CUBIT_TRUE.
    //- Returns CUBIT_SUCCESS if sweep successful, CUBIT_FAILURE otherwise.

  CubitStatus sweep_EDGE_along_vector( EDGE* EDGE_ptr,
                                       BODY *&new_BODY_ptr,
                                       const CubitVector& sweep_vector,
                                       double draft_angle = 0.0,
                                       int draft_type = 0,
                                       bool rigid = false, 
                                       FACE *to_FACE= NULL ) const;
    //- Sweeps the EDGE, EDGE_ptr, along the sweep_vector, a distance equal
    //- to the length of the vector, using the input draft_angle and draft_type.
    //- The draft_type is either 0 or 1 (see the ACIS documentation on the
    //- api_sw_face_vec function) and indicates what should be done to "close"
    //- the end FACEs of the sweep operation, if the draft_angle is non-zero.
    //- The draft angle is in radians.
    //- Returns CUBIT_SUCCESS if sweep successful, CUBIT_FAILURE otherwise.

  CubitStatus sweep_FACE_about_axis( FACE *& FACE_ptr,
                                     const CubitVector& axis_unit_vector,
                                     const CubitVector& point,
                                     double angle,
                                     bool & volume_is_negative,
                                     bool primary_side = true,
                                     int steps = 0,
                                     double draft_angle = 0.0,
                                     int draft_type = 0,
                                     bool rigid = false,
                                     FACE *to_FACE = NULL ) const;
    //- Revolves the FACE, FACE_ptr, about the unit vector, axis_unit_vector,
    //- passing through the spatial location, point.
    //- The angle of revolution is specified by angle and the input
    //- draft_angle and draft_type specify other parameters of the sweep.
    //- The primary side of the FACE is swept unless it is DOUBLE_SIDED
    //- and primary_side is CUBIT_FALSE -- in that case, the face is swept
    //- in the direction of the secondary normal, not the primary normal.
    //- The draft_type is either 0 or 1 (see the ACIS documentation on the
    //- api_sw_face_axis function) and indicates what should be done to "close"
    //- the end FACEs of the sweep operation, if the draft_angle is non-zero.
    //- The draft angle is in radians.
    //- If the sweep failed and if this was due to the creation
    //- of a solid of negative volume, then volume_is_negative is returned as
    //- CUBIT_TRUE.
    //- Returns CUBIT_SUCCESS if sweep successful, CUBIT_FAILURE otherwise.

  CubitStatus sweep_EDGE_about_axis( EDGE* EDGE_ptr,
                                     BODY *&BODYptr,
                                     const CubitVector& axis_unit_vector,
                                     const CubitVector& point,
                                     double angle,
                                     int steps = 0,
                                     double draft_angle = 0.0,
                                     int draft_type = 0,
                                     bool make_solid = false,
                                     bool rigid = false, 
                                     FACE *to_FACE= NULL ) const;
    //- Revolves the EDGE, EDGE_ptr, about the unit vector, axis_unit_vector,
    //- passing through the spatial location, point.
    //- The angle of revolution is specified by angle and the input
    //- draft_angle and draft_type specify other parameters of the sweep.
    //- The draft_type is either 0 or 1 (see the ACIS documentation on the
    //- api_sw_face_axis function) and indicates what should be done to "close"
    //- the end FACEs of the sweep operation, if the draft_angle is non-zero.
    //- The draft angle is in radians.
    //- Returns CUBIT_SUCCESS if sweep successful, CUBIT_FAILURE otherwise.

  CubitStatus sweep_FACE_along_WIRE( FACE *& FACE_ptr,
                                     BODY* wire_BODY_ptr,
                                     double draft_angle = 0.0,
                                     int draft_type = 0,
                                     bool rigid = false, 
                                     FACE *to_FACE= NULL ) const;
    //- Sweeps the FACE, FACE_ptr, along a WIRE body.
    //- The draft_type is either 0, 1 or 2 (see the ACIS documentation on the
    //- sweep_options) and indicates what should be done to "close"
    //- the end FACEs of the sweep operation, if the draft_angle is non-zero.
    //- The draft angle is in radians.
    //- Returns CUBIT_SUCCESS if sweep successful, CUBIT_FAILURE otherwise.

  CubitStatus sweep_EDGE_along_WIRE( EDGE* EDGE_ptr,
                                     BODY* wire_BODY_ptr,
                                     BODY *&new_BODY_ptr,
                                     double draft_angle = 0.0,
                                     int draft_type = 0,
                                     bool rigid = false, 
                                     FACE *to_FACE= NULL ) const;
    //- Sweeps the EDGE, EDGE_ptr, along a WIRE body.
    //- The draft_type is either 0, 1 or 2 (see the ACIS documentation on the
    //- sweep_options) and indicates what should be done to "close"
    //- the end FACEs of the sweep operation, if the draft_angle is non-zero.
    //- The draft angle is in radians.
    //- Returns CUBIT_SUCCESS if sweep successful, CUBIT_FAILURE otherwise.


//HEADER- Functions for geometric imprint operations

  CubitStatus imprint_BODYs( BODY* body1_ptr, BODY* body2_ptr ) const;
    //- Imprints the 2 input ACIS BODYs. Note, no check is done to see if
    //- the 2 BODYs intersect at all. For efficiency, this check ought
    //- to be done first using the GeometryTool::is_interfering procedure.

  CubitStatus imprint( BODY *BODY_ptr, DLIList<FACE*> &FACE_list,
                       DLIList<EDGE*> &EDGE_list ) const;
    //- Imprints a list of FACEs with a list of EDGEs.  The FACEs must
    //- be from the same BODY that is passed in.  FACEs are modified but
    //- EDGEs are not.

  CubitStatus imprint( BODY *BODY_ptr, DLIList<FACE*> &FACE_list,
                       DLIList<DLIList<EDGE*>*> &EDGE_lists_list ) const;
    //- Same as above function, except more efficient if it is known exactly 
    //- which EDGEs need to be imprinted with which FACEs.  Note as input 
    //- there is a corresponding EDGE_list for each input FACE.  The function
    //- does not attempt to imprint all input EDGEs to all input FACEs (as in
    //- the previous function) - only those in corresponding input list 
    //- positions.

  CubitStatus imprint( BODY *BODY_ptr, FACE *FACE_ptr, EDGE *EDGE_ptr ) const;
    //- Imprint given FACE with given EDGE.  BODY_ptr is BODY that FACE_ptr is on.

  EDGE* project_EDGE( EDGE* EDGE_in_ptr, FACE* FACE_ptr, bool print_error = true ) const;
    //- Projects a EDGE to a FACE returning the newly created EDGE

//HEADER- Functions for geometry and topology queries.

  CubitBoolean BODYs_interfering (BODY* body1_ptr, BODY* body2_ptr ) const;
    //- Checks whether there is any interference (intersection) between the
    //- 2 input ACIS BODYs. The BODYs are not affected by this test. Returns
    //- CUBIT_TRUE or CUBIT_FALSE. The entire intersection boolean is not
    //- performed.

//HEADER- Debug and error-checking functions.

  CubitStatus regularize_body( BodySM *body_ptr,
                               BodySM *&new_body_ptr );

    //- Regularizes the body into a new_body.  The original is destroyed.
    //- This function removes topology that is not essential
    //- to the existance of the body.
  CubitStatus regularize_entity( GeometryEntity *old_refentity_ptr,
                                 BodySM *&new_body_ptr );
   //- clean RefEntity

  CubitStatus get_copied_FACES_of_body( DLIList<SurfaceACIS*>& ref_face_list,
                                        DLIList<FACE*>& FACE_list,
                                        BODY*& copied_BODY_ptr ) const;
  //- Given a list of RefFaces, this function returns a list of Acis FACE's
  //- from a BODY that it copies from the parent BODY of the first RefFace in
  //- the incoming list.  This body is also returned, and still contains
  //- all the CUBIT attributes that the original BODIES had. The referenced
  //- RefFaces are removed from the ref_face_list.  Thus, you can loop
  //- with this function until the ref_face_list is empty.  If the function
  //- fails, the copied BODY will have been removed and the returned list
  //- cleaned out.  See AcisTweakTool for uses of this function.

  CubitStatus get_copied_FACES_of_body( DLIList<SurfaceACIS*>& ref_face_list,
                                        DLIList<FACE*>& FACE_list,
                                        DLIList<SurfaceACIS*>& removed_face_list,
                                        BODY*& copied_BODY_ptr ) const;
  //- Same as above, except also returns RefFaces that were removed.

  CubitStatus get_copied_EDGES_of_body( DLIList<CurveACIS*>& ref_edge_list,
                                        DLIList<EDGE*>& EDGE_list,
                                        BODY*& copied_BODY_ptr ) const;
  //- Given a list of RefEdges, this function returns a list of Acis EDGES's
  //- from a BODY that it copies from the parent BODY of the first RefEdge in
  //- the incoming list.  This body is also returned, and still contains
  //- all the CUBIT attributes that the original BODIES had. The referenced
  //- RefEdges are removed from the ref_edge_list.  Thus, you can loop
  //- with this function until the ref_edge_list is empty.  If the function
  //- fails, the copied BODY will have been removed and the returned list
  //- cleaned out.

  CubitStatus get_copied_EDGES_of_body( DLIList<CurveACIS*>& ref_edge_list,
                                        DLIList<EDGE*>& EDGE_list,
                                        DLIList<CurveACIS*>& removed_edge_list,
                                        BODY*& copied_BODY_ptr ) const;
  //- Same as above, except also returns RefEdges that were removed.

  CubitStatus get_copied_VERTICES_of_body( DLIList<PointACIS*>& ref_vertex_list,
                                           DLIList<VERTEX*>& VERTEX_list,
                                           BODY*& copied_BODY_ptr ) const;
  //- Given a list of RefVertexs, this function returns a list of Acis VERTEX's
  //- from a BODY that it copies from the parent BODY of the first RefVertex in
  //- the incoming list.  This body is also returned, and still contains
  //- all the CUBIT attributes that the original BODIES had. The referenced
  //- RefVertex are removed from the ref_vertex_list.  Thus, you can loop
  //- with this function until the ref_vertex_list is empty.  If the function
  //- fails, the copied BODY will have been removed and the returned list
  //- cleaned out.

  CubitStatus get_copied_VERTICES_of_body( DLIList<PointACIS*>& ref_vertex_list,
                                           DLIList<VERTEX*>& VERTEX_list,
                                           DLIList<PointACIS*>& removed_vertex_list,
                                           BODY*& copied_BODY_ptr ) const;
  //- Same as above, except also returns RefVertices that were removed.
  
  Curve* find_curve_by_end_coord( const CubitVector& coords,
                                  DLIList<Curve*>& curves,
                                  bool& start_flag ) const;

  CubitStatus curve_surface_intersection( Surface *surface, 
                                          Curve* curve,
                                          DLIList<Curve*> &new_curves ) const;
  //Intersects input surface with input curve to produce intersection curve(s).
  //If intersection results is nothing or a point, CUBIT_FAILURE is returned.

  CubitStatus embed_curves_into_surface( Surface *surface, 
                                         DLIList<Curve*> &curves_to_imprint,
                                         BodySM *&new_body ) const;
  //Embeds the given curves into the surface.  Equaivalent to imprinting, but 
  //works well for curves that have tolerance problems, (i.e. might not reach
  //the bounding curves of the surface, or might be almost coincident with a
  //vertex on the surface).  
  //WARNING!!! Only use this function if you KNOW that the curves you want 
  //to embed line on the surface.....otherwise you'll get problems.
  //api_embed_wire_in_faces returns unexpected results if the curve
  //does not lie on the surface.


  CubitStatus imprint_overlapping_curves( DLIList<BodySM*> &body_sms,
                                          DLIList<BodySM*> &new_body_sms,
                                          ProgressTool *progress_tool = NULL,
                               DLIList<TopologyBridge*> *new_tbs = NULL,
                               DLIList<TopologyBridge*> *att_tbs = NULL) const;
  //Given a list of bodies, imprints overlapping curves onto each other.   
  //After a call to this function, all overlapping curves between the 
  //bodies are also mergeable.

  CubitStatus imprint_overlapping_surfaces( DLIList<BodySM*> &body_sms,
                                            DLIList<BodySM*> &new_body_sms,
                                            ProgressTool *progress_tool = NULL) const;
  //Given a list of bodies, if Surface A overlaps with surface B, imprints 
  //appropriate curves of surface A onto surface B then curves of B onto A.
void get_new_ENTITIES(ENTITY *top_ENTITY, DLIList<ENTITY*> &new_ENTITIES) const;
void get_att_ENTITIES(ENTITY *top_ENTITY, DLIList<ENTITY*> &att_ENTITIES, char *att_name) const;

  DLIList<EDGE*> find_new_EDGES(BODY *copied_BODY_ptr);
  // Marks all of the edges without Cubit owners as imprint features

};

#endif
