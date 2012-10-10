/**
 * \file GeometryModifyTool.hpp
 *
 * \brief This class provides the interface for all the geometry
 *        related operations to the outside.
 *
 * This class is implemented using the Singleton pattern.
 *
 * \author Xuechen Liu
 *
 * \date 07/11/96
 *
 */

#ifndef GEOMETRYMODIFYTOOL_HPP 
#define GEOMETRYMODIFYTOOL_HPP

#include <stdio.h>
#include <map>

#include "CubitDefines.h"
#include "GeometryDefines.h"
#include "CubitBox.hpp"
#include "CubitVector.hpp"
#include "DLIList.hpp"
#include "CubitGeomConfigure.h"
#include "FacetShapeDefs.hpp"

class CoEdgeSM;
class LoopSM;
class CubitPlane;
class Surface;
class Lump; 
class Curve;
class Body;
class RefVertex;
class RefEdge;
class RefFace;
class RefVolume;
class RefVertex;
class RefFace;
class RefVolume;
template <class X> class DLIList;
class Loop;
class LoopSM;
class RefEntity;
class CoEdgeSM;
class Surface;
class ShellSM;
class Lump;
class TopologyEntity;
class TopologyBridge;
class ShellSM;
class GeometryModifyEngine;
class TopologyBridge;
class TopologyEntity;
class BodySM;
class Surface;
class Curve;
class TBPoint;
class GeometryEntity;


//! Interface class for modifying geometry. 
class CUBIT_GEOM_EXPORT GeometryModifyTool
{

  friend class MergeTool;
public :

  /*!
    *  \return GeometryModifyTool* - Pointer to the singleton GeometryModifyTool object
    *  \arg SMEPtr
    *   Pointer to a SolidModelingEngine object. The default value
    *   is set to NULL.
    *
    *  Return a pointer to the only instance of the class with the default
    *  solid modeling engine (SME) set to the argument, if it is not NULL.
    *  If the argument is NULL, return the pointer to the existing
    *  instance without modifying the default SME. A valid SMEPtr *must*
    *  be supplied at the time of the first call to this function.
    *  Hence, this instance function should specifically be called at
    *  startup with a valid non-NULL input SMEPtr.
    */
  static GeometryModifyTool* instance( GeometryModifyEngine* GMEPtr = NULL);

  static void delete_instance();

  //! \brief Destructor.
  ~GeometryModifyTool();

  /*! Creates an ACIS sphere and assigns it to a Body $
    *  {radius} input radius
    *  Returns the new Body or NULL
    */
  //! \brief Creates a sphere.
  Body* sphere( double radius,
                int x_shift = 0,
                int y_shift = 0,
                int z_shift = 0,
                double inner_radius = 0,
                bool delete_side = false );
  
  /*!  Creates an ACIS cuboid and assigns it to a Body $
    *  {wid, dep, hi} input width, depth, height
    *  Returns the new Body or NULL
    */
  //! \brief Creates a brick (cube).
  Body* brick( double wid, double dep, double hi );

  /*!  Creates an ACIS cuboid and assigns it to a Body $
    *  {center, axes, extension} input center, xyz axes, extension
    *     (extension vector is 1/2 width, height, depth)
    *  If one of dimensions is zero a planar sheet should be created.
    *  Returns the new Body or NULL
    */
  //! \brief Creates a brick (cube).
  Body* brick( const CubitVector &center,
               const CubitVector axes[3],
               const CubitVector &extension );

  /*!  Creates an ACIS prism and assigns it to a Body $
    *  {height, major, minor} input height, major and minor radii. $
    *  {sides} input number of sides. Must be >= 3.
    *  Returns the new Body or NULL
    */
  //! \brief Creates a brick (cube).
  Body* prism( double height, int sides, double major, double minor );

  /*!  Creates an ACIS pyramid and assigns it to a Body $
    *  {height, major, minor} input height, major and minor radii. $
    *  {sides} input number of sides. Must be >= 3.
    *  {top} radius at top of pyramid.
    *  Returns the new Body or NULL
    */
  //! \brief Creates a pyramid. 
  Body* pyramid( double height, int sides, double major, double minor,
                 double top=0.0 );

  /*!  Creates an ACIS frustum and assigns it to a Body $
    *  {hi} input height $
    *  {r1} input radius in x-direction at base $
    *  {r2} input radius in y-direction at base $
    *  {r3} input radius in x-direction at top
    *  Returns the new Body or NULL
    */
  //! \brief Creates a clyinder. 
  Body* cylinder( double hi, double r1, double r2, double r3 );

  /*!  Creates an ACIS torus and assigns it to a Body $
    *  {r1} input major_radius $
    *  {r2} input minor_radius
    *  Returns the new Body or NULL
    */
  //! \brief Creates a torus. 
  Body* torus( double r1, double r2 );

  /*! Creates a Body consisting of a planar sheet (no volume)
    * p1 - 1st corner of the sheet
    * p2 - 2nd corner of the sheet
    * p3 - 3rd corner of the sheet
    * p4 - 4th corner of the sheet
    */
  //! \brief Creates a planar surface from using four vertices.
  Body* planar_sheet( const CubitVector& p1,
                      const CubitVector& p2,
                      const CubitVector& p3,
                      const CubitVector& p4 );

  /*!  Creates an body consisting of a planar sheet (no volume)
    *  {center, axes, width, height} - input location, orientation and
    *                                  size of sheet
    *  Returns the new Body or NULL
    */
  //! \brief Creates a surface from using four vertices.
  Body* planar_sheet( const CubitVector &center,
                      const CubitVector axes[2],
                      double width, double height );

  /*!  Creates an ACIS body consisting of a planar sheet (no volume)
    *  {plane} - plane it will lie in
    *  {bounding_box} - 3D bounding box it must expand.  Plane will
    *                   have minimal area required to just cut the box.
    *  {extension_type} - 0: no extension, 1: percentage, 2: absolute
    *  {extension} - distance sheet is to extend outside of bounding box
    *  Returns the new Body or NULL
    */
  //! \brief Creates a surface from a plane and bounding box. 
  Body* planar_sheet( const CubitPlane& plane,
                      const CubitBox& bounding_box,
                      int extension_type = 0,
                      double extension = 0.0 );

  /*!  Creates a facet body based on facet information
    *  {volume} - topological and facet information of the volume and 
    *             underlying entities
    *  {points} - a vector of doubles representing points used by
    *             the facets
    *  Returns the new Body or NULL
    */
  //! \brief Creates a volume from facet data. 
  Body* create_body( VolumeFacets& volume, std::map<FacetShapes*, RefEntity*>& entity_map,
                     const FacetPointSet& points, int interp_order);

  //! \brief Scale a body non-uniformly. (x,y,z directions) 
  CubitStatus scale ( Body *&entity, 
                      const CubitVector& point,
                      const CubitVector& factors, 
                      bool check_to_transform = true,
                      bool preview = false);

  /*!  \return  Body*
    *  \return - A pointer to a newly created body.
    *  \arg body_ptr
    *  A pointer to a Body which is to used to make copy.
    *
    *  This function makes a copy of the input Body and returns a
    *  pointer to the newly created copy. The input Body and the newly
    *  created Body are geometrically identical. If the copying
    *  operation fails for some reason, the function returns NULL.
    */
  //! \brief Copy a Body.
  Body* copy_body ( Body* body_ptr );

  /*!  Takes a RefEntity (RefVertex, RefEdge, or RefFace) and creates a
    *  copy of it using the make_RefXxxxx commands.
    */
  //! \brief Copy a RefEntity. 
  RefEntity* copy_refentity( RefEntity *old_entity );

  //! \brief Turn bodies inside-out 
  CubitStatus reverse( DLIList<Body*> &body_list );

  //! \brief Reverse the given faces using the solid modeller (flip normals) 
  CubitStatus reverse( DLIList<RefFace*> &ref_face_list );

  /*!  Aligns the body that contains my_face to match my_face
    *  with the target face..
    */
  //! \brief Align the body, which contains my_face, with the target_face.
  CubitStatus align_body( Body *body_ptr,
                          RefFace *my_face,
                          RefFace *target_face,
                          CubitVector &my_center,
                          CubitVector &axis,
                          CubitVector &target_center,
                          double &angle,
                          bool preview);

  //! \brief Aligns a surface or vertex of the body with a plane defined an axis defined by 2 points.
  CubitStatus align_body( Body *body_ptr,
                          DLIList<RefEntity*>& ref_ent_list,
                          CubitVector first_vector,
                          CubitVector second_vector,
                          CubitVector &my_center_1,
                          CubitVector &my_center_2,
                          CubitVector &axis_of_rot,
                          double &angle,
                          double &angle_2,
                          bool preview );

  
  /*!  This function accepts two bodies in the first list.  The second body
       in the list is subtracted from and the first and the result is put
       into the outsideBodies list.  intersectBodies list is the result of 
       intersecting the original two bodies in the first list. 
  */ 
  //! \brief Intersects and subtracts the two bodies in the list. 
  CubitStatus chop( DLIList<Body*> &bodies,
                    DLIList<Body*> &intersectBodies,
                    DLIList<Body*> &outsideBodies,
                    Body*& leftoversBody,
                    bool keep_old = false,
                    bool nonreg = false );

  CubitStatus hollow( DLIList<Body*>& bodies,
                      DLIList<RefFace*> faces_to_remove,
                      DLIList<Body*>& new_bodies,
                      double depth);

  //! \brief Thickens a sheet body (surface) into a solid body.  Can do both directions. 
  CubitStatus thicken( DLIList<Body*>& bodies,
                       DLIList<Body*>& new_bodies,
                       double depth,
                       bool both = false );
  //! \brief Boolean unite.
  CubitStatus unite( DLIList<Body*> &bodies,
                     DLIList<Body*> &newBodies,
                     bool keep_old = false );

  //! \brief Boolean unite.
  CubitStatus unite( DLIList<BodySM*> &body_sm_list,
                     DLIList<BodySM*> &new_body_sm_list,
                     bool keep_old = false);

  //! \brief Boolean subtract.   
  CubitStatus subtract( DLIList<Body*> &tool_body_list, 
                        DLIList<Body*> &from_bodies,
                        DLIList<Body*> &new_bodies,
                        bool imprint = false,
                        bool keep_old = false );

  //! \brief Boolean subtract.  Subtracts the tool from all the blanks.
  CubitStatus subtract( Body* tool_body, 
                        DLIList<Body*> &from_bodies,
                        DLIList<Body*> &new_bodies,
                        bool imprint = false,
                        bool keep_old = false );

  /*! \brief Checks for consistent normals across body. Uses surf_ref
    as a reference surface.  All other surfaces in volume should have 
    normal point in likewise direction (in or out of volume).
    */
  //! \brief Checks for consistent normals across body.  
  CubitStatus validate_normals( DLIList<Body*>& bodies,
                                RefFace *surf_ref,
                                bool reverse );
  
  //! \brief Boolean intersect.
  CubitStatus intersect( Body *tool_body_ptr, 
                         DLIList<Body*> &from_bodies,
                         DLIList<Body*> &new_bodies,
                         bool keep_old = false,
                         bool preview = false );

  //! \brief Boolean intersect.
  CubitStatus intersect( DLIList<Body*> &from_bodies,
                         DLIList<Body*> &new_bodies,
                         bool keep_old = false,
                         bool preview = false);

  /*!  Cuts the given bodies with the points defining a plane.
    *  This is mainly used for vis purposes to section a body with
    *  a plane and only retain half of it.
    */
  //! \brief Cuts body with plane, throwing away one side.
  CubitStatus section( DLIList<Body*> &section_body_list,
                       const CubitVector &point_1,
                       const CubitVector &point_2,
                       const CubitVector &point_3,
                       DLIList<Body*> &new_body_list,
                       CubitBoolean keep_normal_side,
                       CubitBoolean keep_old = CUBIT_FALSE );

  //! \brief Splits periodic curves and surfaces.
  CubitStatus split_periodic( Body *body_ptr,
                              Body *&new_body_ptr );
  
  //! \brief Stitches sheet bodies at common edges to form larger sheet body or solid. 
  CubitStatus stitch( DLIList<Body*> &bodies_to_stitch,
                      DLIList<Body*> &result_list,
                      bool tighten_gaps,
                      double tolerance = -1.0 );

  //! \brief Beta function.
  CubitStatus discover_topology(RefFace *ref_face_ptr,
                                CubitVector &pos,
                                double &step_size,
                                int num_subdivisions);

#ifdef CGM_KCM  
  //! \brief Convert mesh entities to ACIS geometry.
  CubitStatus mesh2brep(std::vector<double> &xvals,
                        std::vector<double> &yvals,
                        std::vector<double> &zvals,
                        std::vector<unsigned int> &tri_connectivity,
                        DLIList<BodySM*> &new_body_sms);
#endif

  //! \brief Beta function.
  void march_path(CubitVector &start_pos,
                  RefFace *start_face,
                  CubitVector &march_dir,
                  double &step_size);

  //! \brief Beta function.
  void subdivide_pie(CubitVector &dir1, CubitVector &dir2, int num_subdivisions,
                                       DLIList<CubitVector> &all_directions);
  
  /*! Split function for simple surface splitting.  Temporary Curves are
    * created on the RefFace and used to split the RefFace.  The curves are
    * created from the input vec_lists (straight lines, arcs or splines).
    * The input locations are the original input locations from the user,
    * which are displayed for user reference.  If the preview_flg is CUBIT_TRUE,
    * the curves are displayed but no splitting occurs.  If the preview_flg is
    * CUBIT_TRUE and the create_ref_edges_flg is CUBIT_TRUE, then RefEdges are
    * created (but no splitting occurs).  The clear_previous_previews flag
    * controls whether previous previews are cleared from the graphics window
    * (only applicable if preview_flg is CUBIT_TRUE).
    */
  //! \brief Splits a surface with a curve defined by positions.
  CubitStatus split_surface( RefFace *ref_face_ptr,
                             DLIList<CubitVector*> &locations,
                             DLIList<DLIList<CubitVector*>*> &vec_lists,
                             CubitBoolean preview_flg = CUBIT_FALSE,
                             CubitBoolean create_ref_edges_flg = CUBIT_FALSE,
                             CubitBoolean clear_previous_previews = CUBIT_TRUE );


  /*! Same behavior as other split_surface function only it handles multiple surface
      and multiple lists of curves.
    */
  //! \brief Splits MULTIPLE surfaces with a curves defined by positions.
  CubitStatus split_surface( DLIList<RefFace*> &ref_face_list,
                             DLIList<CubitVector*> &locations,
                             DLIList<DLIList<DLIList<CubitVector*>*>*> &list_of_vec_lists,
                             CubitBoolean preview_flg = CUBIT_FALSE,
                             CubitBoolean create_ref_edges_flg = CUBIT_FALSE,
                             CubitBoolean clear_previous_previews = CUBIT_TRUE );

  /*!  Splits surfaces by extending curves at the end of hardlines across the
    *   surface.  Only works on planar surfaces and by extending linear curves.
    *   If specified vertex list is provided, only try to extend curves
    *   attached to those vertices, and error out if unable to perform the
    *   operation using any vertex.  If no vertex list is provided, find
    *   appropriate vertices automatically, and limit the split to gaps less
    *   than or equal to the specified "split surface extend gap threshold".
    *   {ref_face_list} - list of connected surfaces to split
    *   {ref_vertex_list} - (OPTIONAL) - vertices to extend from. Function will
    *                       error out if unable to extend from any vertex. If
    *                       not supplied, function will automatically find
    *                       valid vertices to extend from on the surfaces (AUTO
    *                       mode), and will limit the split to gaps less than
    *                       or equal to the specified extend gap threshold (see
    *                       setting below).
    *   {preview_flg} - routine only displays graphics preview of splits
    *   {create_ref_edges_flg} - valid only if preview_flg=CUBIT_TRUE.  If
    *                            CUBIT_TRUE, create RefEdges *instead* of
    *                            splitting.
    * Note: this function uses SplitSurfaceTool.  Three other settings can be
    *       specified that control this function.  They are (static):
    *         SplitSurfaceTool::set_extend_gap_threshold - for AUTO mode,
    *          only split the surface if distance is less than or equal to
    *          this value.
    *         SplitSurfaceTool::set_extend_tolerance - snaps to vertices along
    *           curves within this tolerance so as to avoid creating short
    *           curves when splitting surfaces.
    *         SplitSurfaceTool:set_extend_normal - if TRUE, if a projection
    *          normal to the curve is less distance than a projection in the
    *          direction of the curve, use the normal projection instead.
    */
  //! \brief  Splits surfaces by extending curves at the end of hardlines across the surface.
  CubitStatus split_surfaces_extend( DLIList<RefFace*> &ref_face_list,
                                     DLIList<RefVertex*> &ref_vertex_list,
                                     CubitBoolean preview_flg = CUBIT_FALSE,
                                     CubitBoolean create_ref_edges_flg = CUBIT_FALSE );

  /*!  Splits a surface or connected set of surfaces in one direction.
    *   Particularly useful for splitting fillets to prepare the model
    *   for sweeping.
    *   {ref_face_list} - list of connected surfaces to split
    *   {num_segs} - number of segments to create (i.e, 2 segments
    *                means to split the surface in the middle).
    *   {fraction} - the fraction along the surfaces to split, not valid if
    *                num_segs > 2.  This value is not used if a distance is
    *                specified instead.
    *   {distance} - if 2 segments, allow the split to be at a
    *                user specified distance across the surface.  Specify
    *                -1.0 to use the fraction instead.
    *   {from_curve_ptr} - (OPTIONAL) if user specified a fraction or distance,
    *                      orient from this curve.  If not specified, a
    *                      default is chosen automatically.
    *   {corner_vertex_list} - (OPTIONAL) user can specify corners of surface
    *                       patch.  Optional, if ommitted routine
    *                       will attempt to use smallest angle criteria
    *                       to pick corners.  User should specify corners
    *                       in ccw order - split direction will be
    *                       perpendicular to edge between first two
    *                       corners, unless curve_dir_ptr is specified
    *                       which overrides this direction.
    *   {curve_dir_ptr} - (OPTIONAL) user specifies split direction - valid only
    *                     for a single input surface (not a chain).
    *                     Can be NULL, in which case the routine will
    *                     split along longest aspect ratio (ie.,
    *                     along a fillet)
    *   {through_vertex_list} - (OPTIONAL) user specifies forced vertices for
    *                           the split to run through (on curves).  Not valid
    *                           with more than 2 segments.
    *   {preview_flg} - routine only displays graphics preview of split
    *   {create_ref_edges_flg} - valid only if preview_flg=CUBIT_TRUE.  If
    *                            CUBIT_TRUE, create RefEdges *instead* of
    *                            splitting.
    * Note: this function uses SplitSurfaceTool.  Five other settings can be
    *       specified that control this function.  They are (static):
    *         SplitSurfaceTool::set_parametric_flg.  If TRUE, find spline
    *           locations in parametric space of surface, otherwise do in 3D
    *           and project  back to surface (handy if parametric space is
    *           ill-formed).  DEFAULT = CUBIT_FALSE (do in 3D).  See
    *           SplitSurfaceTool for more information.
    *         SplitSurfaceTool::set_tolerance  Tolerance for finding split
    *           location (determines how accurately curves are tessellated for
    *           approximating the split location, etc.).  See SplitSurfaceTool
    *           for more information.
    *         SplitSurfaceTool::autoDetectTriangles.  If TRUE, detect triangles
    *           on the end of surface chains automatically (see below).
    *         SplitSurfaceTool::sideAngleThreshold.  The closest corner within
    *           this threshold to 180 is removed (if pointAngleThreshold
    *           criteria is also met).
    *         SplitSurfaceTool::pointAngleThreshold.  The corner with angle
    *           below this becomes the triangle point (if sideAngleThreshold
    *           criteria is also met)
    */
  //! \brief Splits a surface or connected set of surfaces in one direction.
  CubitStatus split_surfaces( DLIList<RefFace*> &ref_face_list,
                              int num_segs,
                              double fraction,
                              double distance,
                              RefEdge *from_curve_ptr,
                              DLIList<RefVertex*> &corner_vertex_list,
                              DLIList<RefVertex*> &through_vertex_list,
                              RefEdge *curve_dir_ptr = NULL,
                              CubitBoolean preview_flg = CUBIT_FALSE,
                              CubitBoolean create_ref_edges_flg = CUBIT_FALSE );

  /*! Splits surfaces a specified distance from curves.
  *   {ref_face_list} - list of surfaces to split
  *   {edge_list} - list of curves to offset from
  *   {num_segs} - number of segments to create (i.e, 2 segments
  *                means to split using two offset curves).
  *   {distance} - distance to offset
  *   {partition_flg} - (OPTIONAL) creates partitions of the offset curves. This typically provides
  *                         surfaces easily map meshed
  *   {blunt_flg} - (OPTIONAL) if true the endings of the input curve graph are right angles
  *                     instead of arcs
  *   {preview_flg} - (OPTIONAL) user specifies forced vertices for
  *                        the split to run through (on curves).  Not valid
  *                        with more than 2 segments.
  *   {create_ref_edges_flg} - (OPTIONAL) Create curves instead of splitting surfaces
  *                                 Requires the preview flag
  */
  //! \brief Splits surfaces a specified distance from curves.
  CubitStatus split_surfaces_offset(
      DLIList<RefFace*> &ref_face_list,
      DLIList<RefEdge*> &edge_list,
      int num_segs,
      double distance,
      CubitBoolean partition_flg = CUBIT_FALSE,
      CubitBoolean blunt_flg = CUBIT_FALSE,
      CubitBoolean preview_flg = CUBIT_FALSE,
      CubitBoolean create_ref_edges_flg = CUBIT_FALSE);

  /*! Automatically create a midsurface from a volume.
  *   {body_list_in} - list of bodies to midsurface
  *   {body_list_out} - resulting midsurfaces
  *   {lower_tol} - lower tolerance for finding surface pairs
  *   {upper_tol} - upper tolerance for finding surface pairs
  *   {preview} - preview midsurfaces
  */
  //! \brief Create a midsurface from a volume.
  CubitStatus auto_mid_surface(
	  DLIList<Body*> &body_list_in,
	  DLIList<Body*> &body_list_out,
      DLIList<Body*> &old_bodies_midsurfaced,
      DLIList<double> &thickness_list,
      double lower_tol,
      double upper_tol,
      CubitBoolean delete_midsurfaced,
      CubitBoolean preview);
	/** Automatically create a midsurface from a volume.
	*   {body_list_in} - list of bodies to midsurface
	*   {body_list_out} - resulting midsurfaces
    *   {thickness_list} - list of matching thickness values for the body out list
	*   {lower_tol} - lower tolerance for finding surface pairs
	*   {upper_tol} - upper tolerance for finding surface pairs
	*   {preview} - preview midsurfaces
	*/

  /*!  Webcuts the bodies in the list with a cutting cylinder.
    *  The cylinder is created by the given parameters.  This
    *  is done in the solid modeling engine to reduce the impact
    *  on body ids.
    */
  //! \brief Webcuts bodies with a cylinder.
  CubitStatus webcut_with_cylinder( DLIList<Body*>& webcut_body_list,
                                    double radius,
                                    const CubitVector &axis,
                                    const CubitVector& center,
                                    DLIList<Body*> &results_list,
                                    DLIList<Body*> &neighboring_bodies,
                                    ImprintType imprint_type = NO_IMPRINT,
                                    CubitBoolean merge = CUBIT_FALSE ,
                                    CubitBoolean preview = CUBIT_FALSE);

#ifdef CAT
  CubitStatus webcut_across_translate( DLIList<Body*>& body_list,
                                       RefFace* plane_surf1,
                                       RefFace* plane_surf2,
                                       DLIList<Body*>& results_list,
                                       ImprintType imprint_type = NO_IMPRINT,
                                       CubitBoolean merge = CUBIT_FALSE,
                                       CubitBoolean preview = CUBIT_FALSE);
  /**<  Webcuts with a flat plate to make a body suitable for single-single
    *   sweeping.  Only experimental at this point.
    */
#endif

  /*!  Webcuts the bodies in the list with a cutting brick.
    *  The brick is created by the given parameters - center of
    *  brick, xyz axes, and extension.  Extension is 1/2 width,
    *  height and depth. If one of the brick dimensions is zero
    *  the resultant planar sheet is used to webcut (webcut_with_
    *  planar_sheet is called).  Brick creation is done in the
    *  solid modeling engine to reduce the impact on body ids.
    */
  //! \brief Webcuts bodies with a brick (cube).
  CubitStatus webcut_with_brick( DLIList<Body*>& webcut_body_list,
                                 const CubitVector &center,
                                 const CubitVector axes[3],
                                 const CubitVector &extension,
                                 DLIList<Body*> &results_list,
                                 DLIList<Body*> &neighbor_list,
                                 ImprintType imprint_type = NO_IMPRINT,
                                 CubitBoolean merge = CUBIT_FALSE,
                                 CubitBoolean preview = CUBIT_FALSE);


  /*!  Webcuts the bodies in the list with a cutting planar sheet.
    *  The sheet is created by the given parameters - center of
    *  sheet, xy axes, and width and height. Sheet creation is done
    *  in the solid modeling engine to reduce the impact on body ids.
    */
  //! \brief Webcuts bodies with a planar sheet.
  CubitStatus webcut_with_planar_sheet( DLIList<Body*>& webcut_body_list,
                                        const CubitVector &center,
                                        const CubitVector axes[2],
                                        double width,
                                        double height,
                                        DLIList<Body*> &results_list,
                                        DLIList<Body*> &neighbor_list,
                                        ImprintType imprint_type = NO_IMPRINT,
                                        CubitBoolean merge = CUBIT_FALSE,
                                        CubitBoolean preview = CUBIT_FALSE);


  /*!  \return  int
    *  \return - Number of bodies that were webcut ( >= 0 )
    *  \arg webcut_body_list
    *  The list of bodies to be webcut
    *  \arg plane
    *  The plane to be used for webcutting.
    *  \arg merge
    *  A flag to decide whether the new bodies created by the
    *  webcutting process should be merged or not.
    *  \arg imprint
    *  A flag to decide whether the new bodies created by the
    *  webcutting process should be imprinted or not.
    *
    *  This function webcuts a list of bodies through a plane.
    *  The newly created bodies are merged and imprinted depeding on
    *  the respective flags.
    */
  //! \brief Webcuts bodies with a plane.
  CubitStatus webcut_with_plane( DLIList<Body*>& webcut_body_list,
                                 const CubitVector &vector1,
                                 const CubitVector &vector2,
                                 const CubitVector &vector3,
                                 DLIList<Body*> &results_list,
                                 DLIList<Body*> &neighbor_list,
                                 ImprintType imprint_type = NO_IMPRINT,
                                 CubitBoolean merge = CUBIT_FALSE,
                                 CubitBoolean preview = CUBIT_FALSE) ;

  /*!  \return  int
    *  \return - Number of bodies that were webcut ( >= 0 )
    *  \arg webcut_body_list
    *  The list of bodies to be webcut
    *  \arg refFace
    *  The refFace to be used to construct the plane for
    *  webcutting.
    *  \arg merge
    *  A flag to decide whether the new bodies created by the
    *  webcutting process should be merged or not.
    *  \arg imprint
    *  A flag to decide whether the new bodies created by the
    *  webcutting process should be imprinted or not.
    *
    *  This function webcuts a list of bodies through a plane
    *  defined by a refFace. The newly created bodies are
    *  merged and imprinted depeding on the respective flags.
    *  It returns the number of bodies that were webcut.
    */
  //! \brief Webcuts bodies with a surface -- must be planar.
  CubitStatus webcut_with_surface( DLIList<Body*>& webcut_body_list,
                                   RefFace* refFace,
                                   DLIList<Body*>& results_list,
                                   DLIList<Body*> &neighbor_list,
                                   ImprintType imprint_type = NO_IMPRINT,
                                   CubitBoolean merge = CUBIT_FALSE,
                                   CubitBoolean preview = CUBIT_FALSE);

   /*!  \return  int
    *  \return - Number of bodies that were webcut ( >= 0 )
    *  \arg webcut_body_list
    *  The list of bodies to be webcut
    *  \arg refedge_list
    *  The refEdge list to be used to construct the sheet body for
    *  webcutting.
    *  \arg merge
    *  A flag to decide whether the new bodies created by the
    *  webcutting process should be merged or not.
    *  \arg imprint
    *  A flag to decide whether the new bodies created by the
    *  webcutting process should be imprinted or not.
    *
    *  This function webcuts a list of bodies through a sheet
    *  body defined by a set of curves. The newly created bodies
    *  are merged and imprinted depending on the respective flags.
    *  It returns the number of bodies that were webcut.
    */
  //! \brief Webcuts bodies with a surface defined by a curve loop.
  CubitStatus webcut_with_curve_loop( DLIList<Body*>& webcut_body_list,
                              DLIList<RefEdge*>& refedge_list,
                              DLIList<Body*>& results_list,
                              DLIList<Body*> &neighbor_list,
                              ImprintType imprint_type = NO_IMPRINT,
                              CubitBoolean merge = CUBIT_FALSE,
                              CubitBoolean preview = CUBIT_FALSE);

  /*!  Webcut the bodies using
    *  the surface as the cutting tool.
    *  This is the real webcut with surface.  The others are
    *  just for planes...  The sheet body is a body with 0 volume.
    */
  //! \brief Webcut bodies using a surface, need to be planar. 
  CubitStatus webcut_with_sheet( DLIList<Body*>& webcut_body_list,
                                 Body *sheet_body,
                                 DLIList<Body*> &new_bodies,
                                 DLIList<Body*> &neighbor_list,
                                 ImprintType imprint_type = NO_IMPRINT,
                                 CubitBoolean merge = CUBIT_FALSE,
                                 CubitBoolean preview = CUBIT_FALSE);

  /*!  \return  int
    *  \return - Number of bodies that were webcut ( >= 0 )
    *  \arg webcut_body_list
    *  The list of bodies to be webcut
    *  \arg plane
    *  The plane to be used for webcutting.
    *  \arg merge
    *  A flag to decide whether the new bodies created by the
    *  webcutting process should be merged or not.
    *  \arg imprint
    *  A flag to decide whether the new bodies created by the
    *  webcutting process should be imprinted or not.
    *
    *  This function webcuts a list of bodies through a plane.
    *  The newly created bodies are merged and imprinted depeding on
    *  the respective flags.
    */
  //! \brief Webcuts bodies with a tool body.
  CubitStatus webcut_with_body( DLIList<Body*>& webcut_body_list,
		                Body* body,
		                DLIList<Body*>& results_list,
                                DLIList<Body*> &neighbor_list,
                                ImprintType imprint_type = NO_IMPRINT,
		                CubitBoolean merge = CUBIT_FALSE,
                                CubitBoolean preview = CUBIT_FALSE);

  /*!  Create a sheet extended from the tool surfaces
    *  and then uses it to webcut the bodies.
    */
  //! \brief Webcuts bodies by with an extended surface(s). 
  CubitStatus webcut_with_extended_sheet( DLIList<Body*> &webcut_body_list,
                                          DLIList<RefFace*> &ref_face_list,
                                          DLIList<Body*> &new_bodies,
                                          DLIList<Body*> &neighbor_list,
                                          int &num_cut,
                                          ImprintType imprint_type = NO_IMPRINT,
                                          CubitBoolean merge = CUBIT_FALSE,
                                          CubitBoolean preview = CUBIT_FALSE);
  
  //! \brief Webcuts bodies with a tool solid created by rotating surfaces about an axis.
  CubitStatus webcut_with_sweep_surfaces_rotated(
                            DLIList<Body*> &webcut_body_list,
                            DLIList<RefFace*> &tool_faces,
                            CubitVector &point,
                            CubitVector &sweep_axis,
                            double angle,
                            RefFace* stop_surf,
                            bool up_to_next,
                            DLIList<Body*> &new_bodies,
                            DLIList<Body*> &neighbor_list,
                            ImprintType imprint_type = NO_IMPRINT,
                            CubitBoolean merge = false,
                            CubitBoolean preview = false);

  //! \brief Webcuts bodies with a tool surface created by rotating curves about an axis.
  CubitStatus webcut_with_sweep_curves_rotated(
                            DLIList<Body*> &webcut_body_list,
                            DLIList<RefEdge*> &tool_curves,
                            CubitVector &point,
                            CubitVector &sweep_axis,
                            double angle,
                            RefFace* stop_surf,
                            DLIList<Body*> &new_bodies,
                            DLIList<Body*> &neighbor_list,
                            ImprintType imprint_type = NO_IMPRINT,
                            CubitBoolean merge = false,
                            CubitBoolean preview = false);

  //! \brief Webcuts bodies with a tool solid created by sweeping surfaces along a vector.
  CubitStatus webcut_with_sweep_surfaces(
                            DLIList<Body*> &webcut_body_list,
                            DLIList<RefFace*> &tool_faces,
                            const CubitVector sweep_vector,
                            bool sweep_perp,
                            bool through_all,
                            bool outward,
                            bool up_to_next,
                            RefFace *stop_surf,
                            RefEdge* edge_to_sweep_along,
                            DLIList<Body*> &new_bodies,
                            DLIList<Body*> &neighbor_list,
                            ImprintType imprint_type = NO_IMPRINT,
                            CubitBoolean merge = false,
                            CubitBoolean preview = false);

  //! \brief Webcuts bodies with a tool surface created by sweeping curves along a vector.
  CubitStatus webcut_with_sweep_curves(
                            DLIList<Body*> &webcut_body_list,
                            DLIList<RefEdge*> &tool_curves,
                            const CubitVector sweep_vector,
                            bool through_all,
                            RefFace *stop_surf,
                            RefEdge* edge_to_sweep_along,
                            DLIList<Body*> &new_bodies,
                            DLIList<Body*> &neighbor_list,
                            ImprintType imprint_type = NO_IMPRINT,
                            CubitBoolean merge = false,
                            CubitBoolean preview = false);
  
  //! \brief Removes small topology; i.e, curves and surfaces.
  CubitStatus remove_topology( DLIList<RefEdge*> &ref_edge_list,
        DLIList<RefFace*> &ref_face_list, double backoff_distance,
        double small_curve_size, DLIList<Body*> &new_body_list,
        CubitBoolean propagate,
        CubitBoolean preview);

  /*!  Creates curves offset from a chain of curves.  The offset direction is
    *  only used if there is one linear curve.  Otherwise, the offset direction
    *  is calculated by ACIS (the cross product of the wires tangent and the
    *  planar normal).  The gap type is 0 - rounded, 1 - extended, 2 - natural.
    */
  //! \brief Creates curves offset from a chain of curves.  The offset direction is
  CubitStatus offset_curves( DLIList<RefEdge*>& ref_edge_list,
                             double offset_distance,
                             const CubitVector& offset_direction,
                             int gap_type = 1 );

  /*!  Trims or extends a curve, up to the trim_vector.  If trimming, the
    *  keep_vector determines which side of the curve to keep.  If the curve
    *  is not free, the curve is automatically copied before trimming (so
    *  a new curve results).
    */
  //! \brief Trims or extends a curve, up to the trim_vector.  If trimming, the
  CubitStatus trim_curve( RefEdge* trim_curve,
                          const CubitVector& trim_vector,
                          const CubitVector& keep_vector );

  //! \brief Imprints a group of bodies with one another.
  CubitStatus imprint( DLIList<Body*> &from_body_list,
                       DLIList<Body*> &new_body_list,
                       CubitBoolean keep_old = CUBIT_FALSE );

  /*!  Imprints a list of Bodies with a list of RefEdges.  Useful for
    *  splitting surfaces.  If edge pierces a surface a hardpoint will
    *  result at the pierce location.  Interface is free of ACIS but
    *  currently only works if all entities are ACIS entities.
    */
  //! \brief Imprints a list of Bodies with a list of RefEdges.  Useful for
  CubitStatus imprint( DLIList<Body*> &body_list,
                       DLIList<RefEdge*> &ref_edge_list,
                       DLIList<Body*>& new_body_list,
                       CubitBoolean keep_old_body = CUBIT_FALSE,
                       CubitBoolean show_messages = CUBIT_TRUE );

  /*!  Imprints a list of RefFaces with a list of RefEdges.  This is
    *  useful if the user has a curve which spans several surfaces on
    *  a body and only wants to imprint to selected surfaces.  Algorithm
    *  does not support imprinting to free surfaces.  This method
    *  is not as reliable as the function to imprint a curve to a
    *  body. Interface is free of ACIS but currently only works if all
    *  entities are ACIS entities.
    */
  //! Imprints a list of RefFaces with a list of RefEdges.  This is
  CubitStatus imprint( DLIList<RefFace*> &ref_face_list,
                       DLIList<RefEdge*> &ref_edge_list,
                       DLIList<Body*>& new_body_list,
                       CubitBoolean keep_old_body = CUBIT_FALSE );

  /**<  Imprints a list of Surfaces with list of Curves, sorted per
    *  Surface (ie., curve_lists_list is same length as surface_list).
    *  This version is more efficient than the general-purpose one
    *  above, as we know which curves to imprint with which surfaces.
    *  Also, the Curves need not be RefEntities.
    *  All input surfaces must be from the same body however.
    */
  //! \brief Imprints a list of Surfaces with list of Curves, sorted per
  CubitStatus imprint( DLIList<Surface*> &surface_list,
                       DLIList<DLIList<Curve*>*> &curve_lists_list,
                       Body*& new_body,
                       CubitBoolean keep_old_body = CUBIT_FALSE,
                       CubitBoolean expand = CUBIT_TRUE);

  /*!  Imprints a list of bodies to a list of locations.  Useful for
    *  splitting curves and creating hardpoints on surfaces.  Interface
    *  is free of ACIS but currently only works if bodies are ACIS bodies.
    */
  //! \brief Imprints a list of bodies to a list of locations.  
  CubitStatus imprint( DLIList<Body*> &body_list,
                       DLIList<CubitVector*> &vector_list,
                       DLIList<Body*>& new_body_list,
                       CubitBoolean keep_old_body = CUBIT_FALSE,
                       CubitBoolean merge = CUBIT_FALSE );

  //HEADER- Sweep-related functions. All of these are implemented
  //HEADER- only for RefEntities whose underlying geometry is
  //HEADER- represented by a solid model such as ACIS. The interface
  //HEADER- itself is free of ACIS.


   //! \brief Imprints a list of RefFaces with a list of projected RefEdges.
   CubitStatus imprint_projected_edges( DLIList<RefFace*> &ref_face_list,
                                        DLIList<RefEdge*> &ref_edge_list,
                                        DLIList<Body*>& new_body_list,
                                        CubitBoolean keep_old_body,
                                        CubitBoolean keep_free_edges );

   //! \brief Imprints a list of Bodies with a list of RefEdges which are projected
   //! to a list of RefFaces
   CubitStatus imprint_projected_edges( DLIList<RefFace*> &ref_face_list,
                                        DLIList<Body*> &body_list,
                                        DLIList<RefEdge*> &ref_edge_list,
                                        DLIList<Body*>& new_body_list,
                                        CubitBoolean keep_old_body,
                                        CubitBoolean keep_free_edges );

  //! \brief Tolerantly imprintings a list of curves on a list of surfaces.
  //! Should be used when you have sloppy/out-of-tolerance geometry.
  CubitStatus tolerant_imprint( DLIList<RefFace*> &ref_faces,
                                DLIList<RefEdge*> &ref_edge_list,
                                DLIList<Body*> &new_bodies,
                                bool merge = false );

  //! \brief Tolerantly imprintings a list of curves onto a single surface.
  //!  Should be used when you have sloppy/out-of-tolerance geometry.
  CubitStatus tolerant_imprint( RefFace *ref_face,
                                DLIList<RefEdge*> &ref_edge_list,
                                Body *&new_body, bool merge = false );

  //! \brief Imprints bodies onto one another.  Should be used when you have sloppy/out-
  //! of-tolerance geometry.
  CubitStatus tolerant_imprint( DLIList<Body*> &bodies, DLIList<Body*> &new_bodies,
                                bool merge = false );

  //! \brief Projects list RefEdges onto a list of RefFaces.
  CubitStatus project_edges( DLIList<RefFace*> &ref_face_list,
                              DLIList<RefEdge*> &ref_edge_list_in,
                             DLIList<RefEdge*> &ref_edge_list_new,
                             CubitBoolean trim_projected = CUBIT_FALSE);

  //! \brief Removes all unnessesary faces, curves, vertices from a body. 
  CubitStatus regularize_body( Body *body_ptr, Body *&new_body );

  //! \brief Removes all unnessesary faces, curves, vertices and associated
  //! data from a refentity.
  CubitStatus regularize_refentity( RefEntity *old_entity_ptr,
                                    Body *&new_body_ptr);
  
  //! \brief Tests if a RefEntity can have any unncessary child entities removed,
  //! or simplified away. 
  CubitStatus test_regularize_refentity( RefEntity *old_entity_ptr);

  CubitStatus split_free_curve( RefEdge *ref_edge,
                                DLIList<CubitVector*> &split_locations );

  //! \brief Separates multi-volume bodies into single-volume bodies.
  CubitStatus split_body( Body *body_ptr,
                          DLIList<Body*> &new_bodies ) const;
    
  /*!  Separates surfaces from sheet bodies into separate bodies.  Connected
    *   surfaces will remain connected but be placed in a new body.
    */
  //! \brief Separates surfaces from sheet bodies into separate bodies.  Connected
  CubitStatus separate_surfaces( DLIList<RefFace*> &ref_face_list,
                                 DLIList<Body*> &new_bodies );


  ///< <HR><H3>GeometryModifyTool options and settings</H3>
  
  //! \brief Sets group imprint flag.
  static void set_group_imprint(CubitBoolean flag) {groupImprint = flag;}
  //! \brief Gets group imprint flag.
  static CubitBoolean get_group_imprint() {return groupImprint;}

  //! \brief Sets all edges imprint flag.
  static void set_all_edges_imprint( CubitBoolean flag );
  //! \brief Gets all edges imprint flag.
  static CubitBoolean get_all_edges_imprint();

  //! \brief Gets boolean regularize flag.
  static void boolean_regularize( CubitBoolean flag );
  //! \brief Sets boolean regularize flag.
  static CubitBoolean boolean_regularize();

  //! \brief Gets unite mixed models flag.
  static CubitBoolean unite_mixed_models();
  //! \brief Sets unite mixed models flag.
  static void unite_mixed_models( CubitBoolean flag );

  //! \brief Gets new ids flag.
  static CubitBoolean get_new_ids() {return newIds;}
  //! \brief Sets new ids flag.
  static void set_new_ids(CubitBoolean flag) {newIds = flag;}

  //! \brief Gets old names flag.
  static CubitBoolean get_old_names() { return oldNames; }
  //! \brief Sets old names flag.
  static void set_old_names(CubitBoolean flag) { oldNames = flag; }

  //! \brief Initializes all the settings in GeometryModifyTool to defaults
  static void initialize_settings();


  #ifdef PROE
  CubitStatus prepare_for_topology_update( BodySM* old_bodysm );

  CubitStatus finish_topology_update( BodySM* new_bodysm,
                                      Body* old_body );
  #endif

  //! \brief Gets the entity being copied, if any. 
  static RefEntity* get_copy_entity() { return copyingEntity; }

  //! \brief Sets the entity being copied.
  static void set_copy_entity( RefEntity *ent) { copyingEntity = ent; }
  
  double get_resultant_angle_score(RefFace *narrow_face,
                                   RefFace *other_face,
                                   RefEdge *common_edge);

  double get_neighbor_type_score(RefEdge *common_edge,
                                 RefFace *other_face,
                                 double small_curve_size,
                                 int &neighbor_is_narrow_or_small);

  double get_dihedral_angle_score(RefFace *f1, RefFace *f2,
                                  RefEdge *common_edge);

  double get_edge_type_score(RefEdge *common_edge, double small_curve_size);
  double get_diff_from_multiple_of_90(double angle);
  void split_surface_with_narrow_region(RefFace *face,
                                        DLIList<CubitVector> &split_pos1_list,
                                        DLIList<CubitVector> &split_pos2_list);

  ///< <HR><H3>Topology/geometry creation functions</H3>

  /*!  \return  RefVertex*
    *  \return - A pointer to a newly created RefVertex
    *  \arg ref_vertex_type
    *  The type of the RefVertex
    *  \arg point
    *  The spatial location to use to create the RefVertex
    *  This function takes a type information and a location to create
    *  a RefVertex. The underlying representation of the RefVertex is
    *  determined by the default GeometryQueryEngine.
    *  \arg color
    *  Allows the color to define during creation
    *  defaults to -1 or CUBIT_DEFAULT_COLOR
    */
  //! \brief Creates a point from an x,y,z location
  RefVertex* make_RefVertex(CubitVector const& point, int color = -1) const;
  RefVertex* make_RefVertex( RefVertex *vertex ) const;

  //! \brief Creates a sheet body from a surface.
  Body *make_Body(Surface *surface) const;

  /*!  \return  RefEdge*
    *  \return - A pointer to a newly created RefEdge
    *  \arg ref_edge_type
    *  The type of the RefEdge
    *  \arg ref_vertex_1
    *  The starting point of the to be created RefEdge.
    *  \arg ref_vertex_2
    *  The end point of the to be created RefEdge.
    *  \arg vector_list
    *  List of CubitVectors.
    *  \arg refface_ptr
    *  Input pointer to a RefFace -- defaults to NULL.
    *
    *  This function creates a RefEdge connected to the input RefVertices
    *  as end points. This RefEdge interpolates between the locations
    *  represented by the input set of CubitVectors. A spline curve
    *  is created. If the input refface_ptr is not NULL, the points are
    *  first moved to the surface before interpolation.
    */
  //! \brief Creates a curve from two or more points.
  RefEdge* make_RefEdge( GeometryType ref_edge_type,
                         RefVertex *ref_vertex_1,
                         RefVertex *ref_vertex_2,
                         DLIList<CubitVector*>& vector_list,
                         RefFace* reffaca_ptr = NULL ) const ;

  /*!  \return  RefEdge*
    *  \return - A pointer to a newly created RefEdge
    *  \arg ref_edge_type
    *  The type of the RefEdge
    *  \arg ref_vertex_1
    *  The starting point of the to be created RefEdge.
    *  \arg ref_vertex_2
    *  \arg refface_ptr
    *
    *  This function creates a RefEdge connected to the input RefVertices
    *  as end points. This is along the ref_face. It is basically a straight
    *  line along the surface.  In periodic surfaces, a third point may
    *  be used to create the give the correct direction.
    */
  //! \brief Creates a curve from two or more points, on a surface.
  RefEdge* make_RefEdge(  RefVertex *ref_vertex_1,
                          RefVertex *ref_vertex_2,
                          RefFace* ref_face_ptr,
                          RefVertex const* ref_vertex_3 = NULL ) const ;
  
  /*!  Give a certain ref edge, create a new one.
    *  This is useful when you are creating surfaces from curves already
    *  in the model.  If you use curves existing in the model (attached
    *  to other surfaces) then acis will mess up.  So we basically
    *  need a copy.  This will also create new vertices...
    */
  //! \brief Create a curve from an existing curve.
  RefEdge* make_RefEdge(RefEdge *ref_edge, bool copy_attribs = true) const;
  
  /*!  \return  RefEdge*
    *  \return - A pointer to a newly created RefEdge
    *  \arg ref_edge_type
    *  The type of the RefEdge
    *  \arg ref_vertex_1
    *  The starting point of the to be created RefEdge.
    *  \arg ref_vertex_2
    *  The end point of the to be created RefEdge.
    *  indermediate_point
    *  The intermediate point to be used to form the RefEdge.
    *  For a parabola, the intermediate point is the tip of the
    *  parabola. For an ellipse, the intermediate point is the
    *  center of the ellipse. If the intermediate point is does
    *  not lie at the center of a circle, an logarithmic spiral
    *  curve is generated. For straight curves this argument is
    *  ignored.
    *  \arg sense
    *
    *  The "sense" makes sense (no pun intended) only for ellipses in
    *  the XY plane. For ellipses in the XY plane, the curve can
    *  defined from point-1 to point-2 in two ways. The direction of
    *  the curve can be such that the vector going from the center to
    *  the first point can be rotated in a counter clockwise manner to
    *  reach the vector going from the center to the second point. Such
    *  an ellipse is considered to have positive or FORWARD sense. When
    *  the curve is defined in the other direction it is considered to
    *  have a negative or REVERSED sense.
    *  This argument is ignored for straight and parabolic curves.
    *  This function takes a type information, two RefVertices, an
    *  intermediate position, and a sense information to create a
    *  RefEdge. The RefVertices must be associated with the same
    *  GeometryQueryEngine. The return value can be a NULL
    *  pointer, if the RefEdge cannot be succesfully made for some reason.
    *
    *  Elliptical:
    *  The portion of the ellipse (circle) that is generated goes from
    *  the first point to the second point.
    *
    *  Parabolic:
    *  Construct a parabolic arc from 3 points. The intermediate point
    *  is the peak of the parabola - in this case, the point which is
    *  equidistant from the start and end points of the parabola. The
    *  3 points must form an isosceles triangle. This definition
    *  limits the user to generation of the tip of parabolic shapes only.
    */
  //! \brief Create a curve, i.e. ellipse, parabola, straight, or arc curves.  
  RefEdge* make_RefEdge(GeometryType ref_edge_type,
                        RefVertex *ref_vertex_1,
                        RefVertex *ref_vertex_2,
                        CubitVector const* intermediate_point = NULL ) const;

  RefEdge* make_elliptical_RefEdge( RefVertex *vert1,
                                    RefVertex *vert2,
                                    CubitVector center_point,
                                    double start_angle,
                                    double end_angle,
                                    CubitSense sense) const;


  //! \brief Create a surface from an existing one. 
  RefFace* make_RefFace(RefFace *from_ref_face,
                        CubitBoolean extended_from = CUBIT_FALSE) const;

  /*!  \return Body*
    *  \return - Pointer to a newly created Body object.
    *  \arg ref_face_list
    *  The RefFace list from which we want to create an extended sheet.
    *  \arg clip_box
    *  An optional bounding box to clip the resultant sheet body by.
    *  \arg preview
    *  If true just draw the extended sheet instead of creating it
    *
    * This function creates a sheet body by extending the input surfaces.
    * The result can be optionally clipped to fit inside of the given
    * bounding box.
    */
  //! \brief This function creates a sheet body by extending the input surfaces.
  Body* make_extended_sheet( DLIList<RefFace*> &ref_face_list,
                             CubitBox *clip_box_ptr = NULL,
                             bool preview = false) const;

  /*!  \return  RefFace*
    *  \return - A pointer to a newly created RefFace
    *  \arg ref_face_type
    *  The type of the RefFace
    *  \arg ref_edge_list
    *  The RefEdges to use to create the RefFace. The RefFace will
    *  have just one loop. The first RefEdge in the list will be
    *  first RefEdge of the loop and so on.
    *  \arg ref_face_ptr
    *  The surface that the curves will "cut" a piece out of.  The
    *  curves must lie on this surface.
    *  \check_edges  (specific to Acis)
    *  Creates duplicate edges if these edges are already in use or
    *  if their vertices are merged. For more info see: AME::make_Surface
    *
    *  This function takes a type information and a list of
    *  RefEdges to create a RefFace. The underlying representation
    *  of the RefFace is determined by the GeometryQueryEngine of
    *  the RefEdges. All the RefEdges in the list must be
    *  associated with the same GeometryQueryEngine. The return
    *  value can be a NULL pointer, if the RefFace cannot be succesfully
    *  made for some reason.
    */
  //! \brief Creates a surface from bounding curves.
  RefFace* make_RefFace(GeometryType ref_face_type,
                        DLIList<RefEdge*>& ref_edge_list,
                        bool is_free_face,
                        RefFace *ref_face_ptr = NULL,
                        bool check_edges = true ) const ;

  /*!  \return  Body*
    *  \return - A pointer to a newly created Body
    *  \arg ref_volume_list
    *  The RefVolumes to use to create the Body
    *
    *  This function takes a list of RefVolumes to create a Body. The
    *  underlying representation of the Body is determined by the
    *  GeometryQueryEngine of the RefVolumes. All the RefVolumes
    *  in the list must be associated with the same
    *  GeometryQueryEngine. The return value can be a NULL
    *  pointer, if the RefFace cannot be succesfully made for some
    *  reason.
    */
  //! \brief Creates a body from a list of volumes. 
  Body* make_Body(DLIList<RefVolume*>& ref_volume_list) const ;

  /*!  Creates a body from a ref_face.  This will always be a
    *  sheet body, with no volume, consisting of a single face.
    */
  //! \brief Creates a body from a surface.
  Body* make_Body(RefFace *from_ref_face,
                  CubitBoolean extended_from = CUBIT_FALSE) const;

  /*!  \return  Body*
    *  \return - A pointer to a newly created Body that has just one RefFace
    *  \arg ref_face_type
    *  The type of the RefFace
    *  \arg ref_edge_list
    *  The RefEdges to use to create the RefFace. The RefFace will
    *  have just one loop. The first RefEdge in the list will be
    *  first RefEdge of the loop and so on.
    *
    *  This function takes a type information and a list of RefEdges
    *  to create a Body that has just one RefFace. The underlying
    *  representation of the Body is determined by the
    *  GeometryQueryEngine of the RefEdges. All the RefEdges in
    *  the list must be associated with the same
    *  GeometryQueryEngine. The return value can be a NULL
    *  pointer, if the Body cannot be succesfully made for some
    *  reason.
    *
    *  Each RefEdge in ref_edge_list MUST be a free edge, i.e., not
    *  attached to a RefFace, for this function to succeed!
    */
  //! \brief Creates a body from a surface created from a list of curves.
  Body* make_Body(GeometryType ref_face_type,
                  DLIList<RefEdge*>& ref_edge_list,
                  RefFace *ref_face_ptr = NULL) const ;

  //! \brief Create bodies by sweeping curves or surfaces along a vector.
  CubitStatus sweep_translational(DLIList<RefEntity*>& ref_ent_list,
                                  const CubitVector& sweep_vector,
                                  double draft_angle,
                                  int draft_type,
                                  CubitBoolean switchside,
                                  CubitBoolean rigid,
                                  CubitBoolean anchor_entity,
                                  CubitBoolean keep_old,
                                  DLIList<Body*>& output_body_list);

  CubitStatus sweep_helical( DLIList<RefEntity*>& ref_ent_list,
                             CubitVector &location,
                             CubitVector &direction,
                             double &thread_distance,
                             double &angle,
                             bool right_handed,
                             CubitBoolean anchor_entity,
                             CubitBoolean keep_old,
                             DLIList<Body*>& output_body_list); 

   CubitStatus sweep_curve_target(CubitPlane ref_plane,
						 DLIList<RefEntity*>& ref_ent_list);

   CubitStatus sweep_curve_target(DLIList<RefEdge*>& curve_list,
                            Body *target_body,
                            DLIList<Body*> &out_bodies,
                            CubitVector direction,
                            CubitPlane stop_plane,
                            bool unite);

   CubitStatus sweep_surface_target(RefFace *face,
                                    Body *target_body,
                                    CubitVector distance,
                                    CubitPlane stop_plane,                                    
                                    double magnitude = 0.0 );                                    

   CubitStatus sweep_surface_target(CubitPlane ref_plane,
						DLIList<RefEntity*>& ref_ent_list);

  CubitStatus sweep_perpendicular(DLIList<RefEntity*>& ref_ent_list,
                                  double distance,
                                  double draft_angle,
                                  int draft_type,
                                  CubitBoolean switchside,
                                  CubitBoolean rigid,
                                  CubitBoolean anchor_entity,
                                  CubitBoolean keep_old,
                                  DLIList<Body*>& output_body_list);

  //! \brief Creates bodies by sweeping surfaces or curves about an axis. 
  CubitStatus sweep_rotational(DLIList<RefEntity*>& ref_ent_list,
                               const CubitVector& point,
                               const CubitVector& direction,
                               double angle,
                               DLIList<Body*>& output_body_list,
                               CubitBoolean anchor_entity,
                               CubitBoolean keep_old,
                               int steps = 0,
                               double draft_angle = 0.0,
                               int draft_type = 0,                               
                               CubitBoolean switchside = CUBIT_FALSE,
                               CubitBoolean make_solid = CUBIT_FALSE,
                               CubitBoolean rigid = CUBIT_FALSE );

  //! \brief Creates bodies by sweeping surfaces or curves along a curve. 
  CubitStatus sweep_along_curve(DLIList<RefEntity*>& ref_ent_list,
                                DLIList<RefEdge*>& ref_edge_list,
                                DLIList<Body*>& output_body_list,
                                CubitBoolean anchor_entity,
                                CubitBoolean keep_old,
                                double draft_angle = 0.0,
                                int draft_type = 0,
                                CubitBoolean rigid = CUBIT_FALSE);

  //! \brief Currently unsupported.
  CubitStatus tweak_bend( DLIList<Body*> &bend_bodies,
	                  DLIList<Body*> &new_body_list,
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
                      CubitBoolean preview = CUBIT_FALSE );

  /*! Chamfer curves on solid and sheet bodies.  The left and right offsets
    * are with respect to the curve direction.  If the given right offset is
    * negative, the left offset is used.  Users can preview to clarify the
    * meaning of left and right.
    */
  //! \brief Chamfer curves on solid and sheet bodies. 
  CubitStatus tweak_chamfer( DLIList<RefEdge*> &ref_edge_list,
                             double left_offset,
                             DLIList<Body*> &new_body_list,
                             double right_offset = -1.0,
                             CubitBoolean keep_old_body = CUBIT_FALSE,
                             CubitBoolean preview = CUBIT_FALSE );

  /*! Chamfer vertices on solid or sheet bodies.  On a solid body there can
    * be up to 3 offsets; on a sheet body up to 2 offsets.  The offsets are
    * in the direction of the supplied edges.  If multiple vertices are
    * supplied, only one offset value is allowed and the edges are not used.
    */
  //! \brief Chamfer vertices on solid or sheet bodies.
  CubitStatus tweak_chamfer( DLIList<RefVertex*> &ref_vertex_list,
                             double offset1,
                             DLIList<Body*> &new_body_list,
                             RefEdge *edge1 = NULL,
                             double offset2 = -1.0,
                             RefEdge *edge2 = NULL,
                             double offset3 = -1.0,
                             RefEdge *edge3 = NULL,
                             CubitBoolean keep_old_body = CUBIT_FALSE,
                             CubitBoolean preview = CUBIT_FALSE );

  //! \brief Creates a round fillet (or blend) at the given curves on solid or
  //!  sheet bodies.
  CubitStatus tweak_fillet( DLIList<RefEdge*> &ref_edge_list,
                            double radius,
                            DLIList<Body*> &new_body_list,
                            CubitBoolean keep_old_body = CUBIT_FALSE,
                            CubitBoolean preview = CUBIT_FALSE );

  /*! Creates a round fillet (or blend) at the given curve on a solid or sheet
    * body.  The fillet has a variable radius from the start to the end of
    * the curve.
    */
  //! \brief Creates a round fillet (or blend) at the given curve. 
  CubitStatus tweak_fillet( RefEdge *ref_edge_ptr,
                            double start_radius,
                            double end_radius,
                            Body *&new_body_ptr,
                            CubitBoolean keep_old_body = CUBIT_FALSE,
                            CubitBoolean preview = CUBIT_FALSE );

  //! \brief Create a round fillet (or blend) at the given vertices on sheet bodies.
  CubitStatus tweak_fillet( DLIList<RefVertex*> &ref_vertex_list,
                            double radius,
                            DLIList<Body*> &new_body_list,
                            CubitBoolean keep_old_body = CUBIT_FALSE,
                            CubitBoolean preview = CUBIT_FALSE );

  //! \brief Tweak specified faces of a volume or volumes along a vector.
  CubitStatus tweak_move( DLIList<RefFace*> &ref_face_list,
                          const CubitVector &delta,
                          DLIList<Body*> &new_body_list,
                          CubitBoolean keep_old_body = CUBIT_FALSE,
                          CubitBoolean preview = CUBIT_FALSE );

  //! \brief Tweak specified curves of a sheet body along a vector.
  CubitStatus tweak_move( DLIList<RefEdge*> &ref_edge_list,
                          const CubitVector &delta,
                          DLIList<Body*> &new_body_list,
                          CubitBoolean keep_old_body = CUBIT_FALSE,
                          CubitBoolean preview = CUBIT_FALSE );

  /*! Tweak specified faces of a volume or volumes by offsetting those faces
    * by the offset distance. Optionally supply additional faces with
    * different offset distances.
    */
  //! \brief Offsets a surface(s) on a volume(s).
  CubitStatus tweak_offset( DLIList<RefFace*> &ref_face_list,
                            double offset_distance,
                            DLIList<RefFace*> *add_ref_face_list_ptr,
                            DLIList<double> *add_offset_list_ptr,
                            DLIList<Body*> &new_body_list,
                            CubitBoolean keep_old_body = CUBIT_FALSE,
                            CubitBoolean preview = CUBIT_FALSE );

  /*! Tweak specified curves of a sheet body or bodies by offsetting those
    * curves by the offset distance.  Optionally supply additional curves with
    * different offset distances.
    */
  //! \brief Offset curves on sheet bodies.
  CubitStatus tweak_offset( DLIList<RefEdge*> &ref_edge_list,
                            double offset_distance,
                            DLIList<RefEdge*> *add_ref_face_list_ptr,
                            DLIList<double> *add_offset_list_ptr,
                            DLIList<Body*> &new_body_list,
                            CubitBoolean keep_old_body = CUBIT_FALSE,
                            CubitBoolean preview = CUBIT_FALSE );

  //! \brief Performs tweak_remove operation on surfaces individually. 
  CubitStatus tweak_remove_individually( DLIList<RefFace*> &ref_face_list,
                                         DLIList<Body*> &new_body_list,
                                         CubitBoolean keep_surface = CUBIT_FALSE,
                                         CubitBoolean keep_old_body = CUBIT_FALSE,
                                         CubitBoolean preview = CUBIT_FALSE );

  //! \brief Performs tweak_remove operation on surfaces collectively. 
  CubitStatus tweak_remove_together( DLIList<RefFace*> &ref_face_list,
                                     DLIList<Body*> &new_body_list,
                                     CubitBoolean extend_adjoining = CUBIT_TRUE,
                                     CubitBoolean keep_surface = CUBIT_FALSE,
                                     CubitBoolean keep_old_body = CUBIT_FALSE,
                                     CubitBoolean preview = CUBIT_FALSE);

  /*! Remove curves from a sheet body or bodies and then extend the remaining
    * curves to fill the gap.  If an internal loop of curves is removed the
    * hole is removed.
  */
  //! \brief Removes a surface from a volume, extending neighboring surfaces.
  CubitStatus tweak_remove( DLIList<RefEdge*> &ref_edge_list,
                            DLIList<Body*> &new_body_list,
                            CubitBoolean keep_old_body = CUBIT_FALSE,
                            CubitBoolean preview = CUBIT_FALSE );

  /*!  Tweak specified faces of a volume or volumes up to target surface.
    *  If extend flag is true, extend out the targets before
    *   tweaking to them (only used for multiple targets; single targets are
    *   always extended).  The optional limit plane is only valid if extend_flg
    *   is TRUE; it will limit the tweak to not go past this plane in the case
    *   where the tweaked body would only partially intersect the extended
    *   targets.The reverse flag should never be needed - if it is there may be
    *   a bug or a bad normal on a body (i.e., negative volume body), and is
    *   only retained for debugging.
    */
  //! \brief Extends (tweaks) surfaces up to a target surface.
  CubitStatus tweak_target( DLIList<RefFace*> &ref_face_list,
                            DLIList<RefFace*> &target_face_list,
                            DLIList<Body*> &new_body_list,
                            CubitBoolean extend_flg = CUBIT_TRUE,
                            CubitPlane *limit_plane = NULL,
                            CubitBoolean reverse_flg = CUBIT_FALSE,
                            CubitBoolean keep_old_body = CUBIT_FALSE,
                            CubitBoolean preview = CUBIT_FALSE );

  /*!  Tweak specified faces of a volume or volumes up to target plane. 
    *   If extend flag is true, extend out the targets before
    *   tweaking to them (only used for multiple targets; single targets are
    *   always extended).  The optional limit plane is only valid if extend_flg
    *   is TRUE; it will limit the tweak to not go past this plane in the case
    *   where the tweaked body would only partially intersect the extended
    *   targets.The reverse flag should never be needed - if it is there may be
    *   a bug or a bad normal on a body (i.e., negative volume body), and is
    *   only retained for debugging.
    */

  //! \brief Extends (tweaks) surfaces up to a target plane.
  CubitStatus tweak_target( DLIList<RefFace*> &ref_face_list,
                            CubitPlane &plane,
                            DLIList<Body*> &new_body_list,
                            CubitBoolean reverse_flg = CUBIT_FALSE,
                            CubitBoolean keep_old_body = CUBIT_FALSE,
                            CubitBoolean preview = CUBIT_FALSE );

  /*!   Tweak specified edges of a surface or set of surfaces (in sheet
    *   bodies) up to target surfaces or plane.  This essentially extends
    *   or trims the attached surfaces of the sheet body.  If extend flag is
    *   true, extend out the targets before tweaking to them (only used for
    *   multiple targets; single targets are always extended).  The optional
    *   limit plane is only valid if extend_flg is TRUE; it will limit the
    *   tweak to not go past this plane in the case where the tweaked body
    *   would only partially intersect the extended targets.The reverse flag
    *   should never be needed - if it is there may be a bug or a bad normal
    *   on a body (i.e., negative volume body), and is only retained for
    *   debugging.
    */
  //! \brief Extends (tweaks) curves up to a target surface.
  CubitStatus tweak_target( DLIList<RefEdge*> &ref_edge_list,
                            DLIList<RefFace*> &target_face_list,
                            DLIList<Body*> &new_body_list,
                            CubitBoolean extend_flg = CUBIT_TRUE,
                            CubitPlane *limit_plane = NULL,
                            CubitBoolean reverse_flg = CUBIT_FALSE,
                            CubitBoolean keep_old_body = CUBIT_FALSE,
                            CubitBoolean preview = CUBIT_FALSE,
                            double max_area_increase = 0 );
  CubitStatus tweak_target( DLIList<RefEdge*> &ref_edge_list,
                            CubitPlane &plane,
                            DLIList<Body*> &new_body_list,
                            CubitBoolean reverse_flg = CUBIT_FALSE,
                            CubitBoolean keep_old_body = CUBIT_FALSE,
                            CubitBoolean preview = CUBIT_FALSE );

  /*!   Tweak specified edges of a surface or set of surfaces (in sheet
    *   bodies) up to target plane.  This essentially extends
    *   or trims the attached surfaces of the sheet body.  If extend flag is
    *   true, extend out the targets before tweaking to them (only used for
    *   multiple targets; single targets are always extended).  The optional
    *   limit plane is only valid if extend_flg is TRUE; it will limit the
    *   tweak to not go past this plane in the case where the tweaked body
    *   would only partially intersect the extended targets.The reverse flag
    *   should never be needed - if it is there may be a bug or a bad normal
    *   on a body (i.e., negative volume body), and is only retained for
    *   debugging.
    */
  //! \brief Extends (tweaks) curves up to a target plane.
  CubitStatus tweak_target( DLIList<RefEdge*> &ref_edge_list,
                            DLIList<RefEdge*> &target_edge_list,
                            DLIList<Body*> &new_body_list,
                            CubitBoolean extend_flg = CUBIT_TRUE,
                            CubitPlane *limit_plane = NULL,
                            CubitBoolean reverse_flg = CUBIT_FALSE,
                            CubitBoolean keep_old_body = CUBIT_FALSE,
                            CubitBoolean preview = CUBIT_FALSE,
                            double max_area_increase = 0 );
  /**<  Tweak specified edges of a sheet body or bodies up to target curves
    *   that are part of a sheet body.  The target is a surface or surfaces
    *   created by thickening the owning surface of the target curve(s).
    *   If extend flag is true, extend out the targets before tweaking to them
    *   (only used for multiple targets; single targets are always extended).
    *   The optional limit plane is only valid if extend_flg is TRUE; it will
    *   limit the tweak to not go past this plane in the case where the tweaked
    *   body would only partially intersect the extended targets. The reverse
    *   flag should never be needed - if it is there may be a bug or a bad normal
    *   on a body (i.e., negative volume body), and is only retained for
    *   debugging.
    */
  //! \brief Tweak specified vertex of a sheet body to a given location.
  CubitStatus tweak_target( RefVertex *ref_vertex_ptr,
                            DLIList<RefFace*> &modify_ref_face_list,
                            CubitVector &target_loc,
                            Body *&new_Body_ptr,
                            CubitBoolean keep_old_body = CUBIT_FALSE,
                            CubitBoolean preview = CUBIT_FALSE );

  //! \brief Converts edges smaller than 'lengthlimit' into tolerant (or fat) vertices.
  //! ACIS only.
  CubitStatus remove_curve_slivers( DLIList<Body*> &bodies, double lengthlimit );

  //! \brief Create a surface using arrays of points in u and v directions. 
  CubitStatus create_net_surface( DLIList<Surface*>& ref_face_list, BodySM *& new_body,
                                  DLIList<DLIList<CubitVector*>*> &vec_lists_u,
                                  DLIList<DLIList<CubitVector*>*> &vec_lists_v,
                                  double net_tol = 1e-3,
                                  CubitBoolean heal = CUBIT_TRUE );

  //! \brief Create a surface using curves in u and v directions. 
  CubitStatus create_net_surface( DLIList<RefEdge*>& u_curves, 
                                  DLIList<RefEdge*>& v_curves,
                                  double net_tol = 1e-3,
                                  CubitBoolean heal = CUBIT_TRUE );

  //! \brief Creates an offset surface some distance from a surface.
  CubitStatus create_offset_surface( RefFace* ref_face_ptr,
                                     double offset_distance );

  /*! Create a sheet body (or bodies) by offsetting the given faces. The
    *  optional additional face list and double list (must be same length)
    *  allow different offset distances for different faces. Adjoining faces
    *  are extended or trimmed to remain joined in the new sheet body.  Radial
    *  faces that cannot be so offset are removed and the resulting wound
    *  healed by the surrounding faces.
    */
  //! \brief Create a body(s) by offsetting the given surfaces. 
  CubitStatus create_offset_sheet( DLIList<RefFace*> &ref_face_list,
                                   double offset_distance,
                                   DLIList<RefFace*> *add_ref_face_list_ptr,
                                   DLIList<double> *add_offset_list_ptr,
                                   DLIList<Body*> &new_body_list,
                                   CubitBoolean preview = CUBIT_FALSE );

  //! \brief Create a body by offsetting all the surface of the body. 
  CubitStatus create_offset_body( Body *body_ptr, Body *&new_body,
                                  double offset_distance );

  //! \brief Create a sheet body skinning (lofting) a series of curves.
  CubitStatus create_skin_surface( DLIList<RefEdge*>& ref_edges, 
                                   Body*& new_body,
                                   DLIList<RefEdge*>& guides);
  
  ////! \brief Create a solid body by lofting between two surfaces.
  //CubitStatus loft_surfaces_to_body( RefFace *face1, const double &takeoff1,
		//		     RefFace *face2, const double &takeoff2,
  //                                   Body*& new_body,
  //                                   CubitBoolean arc_length_option,
  //                                   CubitBoolean twist_option,
  //                                   CubitBoolean align_direction,
  //                                   CubitBoolean perpendicular,
  //                                   CubitBoolean simplify_option);





  ////! \brief Create a solid body by lofting between two surfaces.
  //CubitStatus loft_surfaces_to_body( RefFace *face1, const double &takeoff1,
  //                                   RefFace *face2, const double &takeoff2,
  //                                   DLIList<RefEdge*> &guides,
  //                                   Body*& new_body,
  //                                   CubitBoolean arc_length_option,
  //                                   CubitBoolean twist_option,
  //                                   CubitBoolean align_direction,
  //                                   CubitBoolean perpendicular,
  //                                   CubitBoolean simplify_option);


  CubitStatus loft_surfaces_to_body(DLIList<RefFace*> &surfaces,
                                    DLIList<double> &takeoff_factor_list,
                                    DLIList<RefFace*> &takeoff_vector_surface_list,
                                    DLIList<CubitVector> &surface_takeoff_vector_list,
                                    DLIList<RefEdge*> &takeoff_vector_curve_list,
                                    DLIList<CubitVector> &curve_takeoff_vector_list,
                                    DLIList<RefEdge*> &guides,
                                    DLIList<RefVertex*> &match_vertices_list,
                                    Body*& new_body,
                                    CubitBoolean global_guides,
                                    CubitBoolean closed,
                                    CubitBoolean show_matching_curves,
                                    CubitBoolean preview
                                    );





  CubitStatus create_surface_curve(DLIList<RefEntity*> curve_entity,
                                   DLIList<RefEntity*> target_entity,
				   CubitVector sweep_direction = CubitVector (0,0,0),
				   CubitBoolean distance_flag = CUBIT_FALSE);

  CubitStatus create_rectangle_surface( double width, double height, CubitVector plane );

  CubitStatus create_parallelogram_surface( RefVertex *v1, 
                                            RefVertex *v2,
                                            RefVertex *v3 );

  CubitStatus create_circle_surface( double radius, CubitVector plane ); 

  CubitStatus create_circle_surface( RefVertex *v1, 
                                     RefVertex *v2,
                                     RefVertex *v3 );

  CubitStatus create_circle_surface( RefVertex *v1, 
                                     RefVertex *v2,
                                     CubitVector center_point ); 

  CubitStatus create_ellipse_surface( RefVertex *v1, 
                                      RefVertex *v2,
                                      CubitVector center_point ); 

  CubitStatus create_ellipse_surface( double major_radius, 
                                      double minor_radius,
                                      CubitVector plane );

  CubitStatus idealize_hole_slot_geometry(DLIList<RefEntity*> idealize_entity,
                                     DLIList<RefEntity*> exclude_entity,
                                     double arc_radius,
                                     double slot_arc_radius,
                                     double slot_length,
                                     CubitBoolean preview = CUBIT_FALSE);

  CubitStatus idealize_fillet_geometry(DLIList<RefEntity*> idealize_entity,
                                       DLIList<RefEntity*> exclude_entity,
                                       double fillet_rad,
                                       CubitBoolean internal_flg,
                                       CubitBoolean external_flg,
						               CubitBoolean preview = CUBIT_FALSE);

  CubitStatus create_surface_doubler(DLIList<RefEntity*> doubler_entity,
                                    DLIList<RefEntity*> target_entity,
                                    DLIList<Body*> &body_list_out,
                                    CubitBoolean internal_flg = CUBIT_FALSE,
                                    CubitBoolean extend_flg = CUBIT_TRUE,
                                    CubitPlane *limit_plane = NULL,
                                    CubitVector sweep_direction = CubitVector(0,0,0),
                                    CubitBoolean preview = CUBIT_FALSE);
  
  /*! Fits an analytic spline surface through positions.  If a surface is specified, 
      positions are forced to lie on the surface and restrains the resultant surface
      geometry to match the surface's geometry.  The project parameter will project 
      the nodes or vertices to the specified surface. 
  */
  //! \brief  Fits an analytic spline surface through positions. 
  CubitStatus create_surface( DLIList<CubitVector*>& vec_list,
                              Body *&new_body,
                              RefFace *ref_face_ptr,
			      CubitBoolean project_points );

  CubitStatus create_surface( DLIList<RefVertex*> &vert_list, 
                              Body *&new_body );


  //! \brief Creates a simple triangular weld surface.
  CubitStatus create_weld_surface( CubitVector &root,
                                   RefFace *ref_face1,
                                   double leg1,
                                   RefFace *ref_face2,
                                   double leg2,
                                   Body *&new_body );

  /*!  Creates a body out of the surfaces in the ref_face_list.  The surfaces
    *  can be either free surfaces or sheet bodies or both.
    *  The body should be healed.  The bounding curves and vertices between
    *  the surfaces must be within tolerance.
    */
  //! \brief Creates a body out of the specified surfaces. 
  CubitStatus create_solid_bodies_from_surfs( DLIList<RefFace*> &ref_face_list,
                                      DLIList<Body*> &new_bodies,
                                      CubitBoolean keep_old = CUBIT_FALSE,
                                      CubitBoolean heal = CUBIT_TRUE,
                                      CubitBoolean sheet = CUBIT_FALSE ) const;

  //! \brief Create curves from the intersection of two surfaces.
  CubitStatus surface_intersection( RefFace *ref_face1,
                                    RefFace *ref_face2,
                                    DLIList<RefEdge*> &ref_edge_list );

  RefEdge* create_arc(const CubitVector& position,
                               double radius,
                               double start_angle,
                               double end_angle,
                               CubitVector plane,
                               CubitBoolean preview = CUBIT_FALSE);

  RefEdge* create_arc_radius( RefVertex* ref_vertex1,
                                   RefVertex* ref_vertex2,
                                   const CubitVector &normal,
                                   double radius,
                                   CubitBoolean other_arc = CUBIT_FALSE ,
                                   CubitBoolean full = CUBIT_FALSE,
                                   CubitBoolean preview = CUBIT_FALSE);
  //! \brief Create an arc curve from three points.
  RefEdge* create_arc_three( RefVertex* ref_vertex1,
                             RefVertex* ref_vertex2,
                             RefVertex* ref_vertex3,
                             CubitBoolean full = CUBIT_FALSE,
                             CubitBoolean preview = CUBIT_FALSE );

  //! \brief Create an arc curve tangent to three curves.
  RefEdge* create_arc_three( RefEdge* ref_edge1,
                             RefEdge* ref_edge2,
                             RefEdge* ref_edge3,
                             CubitBoolean full = CUBIT_FALSE,
                             CubitBoolean preview = CUBIT_FALSE );
  
  /*! \brief Create an arc curve from two points and a center point. 
      If full option is specified, a full circle is created.*/
  //! \brief Create an arc curve from two points and a center point. 
  RefEdge* create_arc_center_edge( RefVertex* ref_vertex1,
                                   RefVertex* ref_vertex2,
                                   RefVertex* ref_vertex3,
                                   const CubitVector &normal,
                                   double radius = CUBIT_DBL_MAX,
                                   CubitBoolean full = CUBIT_FALSE,
                                   CubitBoolean preview = CUBIT_FALSE );

  //! \brief  Create a curve that is a combination of the specified curves.
  CubitStatus create_curve_combine( DLIList<RefEdge*>& ref_edge_list,
                                    RefEdge *&new_ref_edge_ptr );


  CubitStatus create_curve_helix( CubitVector &location,
                                  CubitVector &direction,
                                  CubitVector &start_point,
                                  double &thread_distance,
                                  double &angle,
                                  bool right_handed,
                                  RefEdge *&new_ref_edge_ptr);


  //! \brief Sets sepAfterWebcut variable.
  static void set_sep_after_webcut_setting(CubitBoolean val) {sepAfterWebcut = val;}

  //! \brief Gets sepAfterWebcut variable.
  static CubitBoolean get_sep_after_webcut_setting() {return sepAfterWebcut;}

  /*!  Returns CUBIT_TRUE if all the entities have the same geometric query engine and
    *  if that is the same one as the default. */
  //! \brief Check to determine if all entities are of the same geometry engine.
  CubitBoolean same_modify_engine(DLIList<TopologyEntity*> &topo_list) const;

  /*!  Returns CUBIT_TRUE if all the entities have the same geometric query engine and
    *  if that is the same one as the default.  If the check_children parameter is
    *  CUBIT_TRUE, all the child entities will also be checked.
    */
  //! \brief Check to determine if all entities are of the same geometry engine as 
  //! active geometry engine.
  CubitBoolean same_modify_engine(DLIList<RefEntity*> &ref_entity_list,
				  CubitBoolean check_children = CUBIT_FALSE) const;

  
  //! \brief Determines if specified RefFaces and RefEdges are of the same geometry engine.
  //! Returns common engine, if any, and underlying surfaces and curves. 
  GeometryModifyEngine* common_modify_engine( DLIList<RefFace*>& faces,
                                              DLIList<RefEdge*>& edges,
                                              DLIList<Surface*>& surfaces,
                                              DLIList<Curve*>& curves,
                                              CubitBoolean allow_composites = CUBIT_FALSE) const;



  //! \brief Add a geometry modify engine to the list.
  void add_gme(GeometryModifyEngine *gme_ptr);

  //! \brief < remove a geometry modify engine from the list; returns CUBIT_FAILURE
  //!  if it wasn't in the list.
  CubitStatus remove_gme(GeometryModifyEngine *gme_ptr);

  //! \brief Returns a list of GeometryModifyEngines.
  void get_gme_list(DLIList<GeometryModifyEngine*> &gme_list);

  //! \brief Gets the current active GeometryModifyEngine.
  GeometryModifyEngine *get_gme() const;


  //! \brief Gets the GeometryModifyEngine of this entity.
  GeometryModifyEngine *get_engine(TopologyBridge *tb_ptr) const;

  //! \brief Gets the GeometryModifyEngine of this entity.
  GeometryModifyEngine *get_engine(TopologyEntity *te_ptr,
                                   TopologyBridge** bridge = 0) const;

  /*! Finds the intersections of a certain distance (offset) between two
    * curves.  The two curves must lie in a plane.  The first curve is offset
    * the offset distance in both directions, and the bounded intersections with
    * the second curve are found.  The first curve can optionally be extended
    * to infinity for the intersection calculation.  The intent of the function
    * is so that the user can create a point on a curve a certain distance
    * from another curve, as in specifying a reference location for a gage
    * diameter on an arc in an engineering drawing.  The function allocates the
    * CubitVectors in the returned list, so be sure to free them.
    */
  //! \brief Finds the intersections of a certain distance (offset) between two curves.
  CubitStatus get_offset_intersections( RefEdge* ref_edge1, RefEdge* ref_edge2,
                                        DLIList<CubitVector*>& intersection_list,
                                        double offset, CubitBoolean ext_first = CUBIT_TRUE );

  /*! Finds intersections (points) of the curve and surface.  The surface can
    * be offset - it is offset to each side and intersections are found.  By
    * default the surface is extended to infinity (if possible) and the
    * intersections are found.  The function allocates the CubitVectors in
    * the returned list, so be sure to free them.
    */
  //! \brief Finds intersections (points) of the curve and surface.  
  CubitStatus get_offset_intersections( RefEdge* ref_edge_ptr, RefFace* ref_face_ptr,
                                        DLIList<CubitVector*> &intersection_list,
                                        double offset = 0.0,
                                        CubitBoolean ext_surf = CUBIT_TRUE );

  //! \brief From two surface, create a midplane, then trim it with a body.
  CubitStatus get_mid_plane( RefFace *ref_face_1,
                             RefFace *ref_face_2,
                             Body *body_to_trim_to,
                             DLIList<RefFace*> &mid_plane_surfs ) const;

  //! \brief Given 2 surfaces, this returns trimmed surfaces of
  //!  the midsurface (this is an ALPHA feature).
  CubitStatus get_mid_surface( RefFace *ref_face_1,
                             RefFace *ref_face_2,
                             Body *body_to_trim_to,
                             DLIList<RefFace*> &mid_plane_surfs ) const;

  //! \brief Sets the active geometry engine.
  CubitStatus set_default_gme(GeometryModifyEngine* GMEPtr);

  //! Get a subset of bodies from 'remaining_bodies' that share
  //! a common GeometryModifyEngine.  Remove them from
  //! 'remaining_bodies', put the Bodies and corresponding BodySMs in
  //! 'engine_bodies' and 'engine_body_sms', and return the engine.
  //! Returns NULL if all bodies remaining in the list have no
  //! modify engine.
  //! \brief Groups Bodies with the same underlying geometry engine into a list.
  GeometryModifyEngine* group_bodies_by_engine( DLIList<Body*>& remaining_bodies,
                                                DLIList<Body*>& engine_bodies,
                                                DLIList<BodySM*>& engine_body_sms ) const;
protected :

  GeometryModifyTool(GeometryModifyEngine* GMEPtr);
  /**<  Constructor for the (singleton) GeometryModifyTool object
    */

  GeometryModifyEngine* make_RefEdge_common ( RefVertex* start_vertex,
                                              RefVertex* end_vertex,
                                              TBPoint*& start_point,
                                              TBPoint*& end_point,
                                              RefFace* ref_face = 0,
                                              Surface** surface = 0) const;
    //- Common code for misc. make_RefEdge functions.
    //- Input  : start and end vertices and an optional RefFace pointer.
    //- Returns: The modify engine to use to create the curve.
    //- Output : points in the returned modify engine and if a refface
    //-          was given, a surface in the modify engine.
    //- Determines which engine to use to create the curve and returns
    //- TopologyBridges owned by that surface.  New Points are created
    //- if either the RefVertex already has a parent RefEdge or because
    //- a common modify engine could not be found for all the input
    //- entities.

  CubitStatus okay_to_modify( DLIList<Body*>& bodies, const char* op ) const;

  GeometryModifyEngine* common_modify_engine( DLIList<Body*>& bodies,
                                              DLIList<BodySM*>& bodysms ) const;

  GeometryModifyEngine* common_modify_engine( DLIList<RefFace*>& faces,
                                              DLIList<Surface*>& surfaces,
                                              CubitBoolean allow_composites = CUBIT_FALSE ) const;
  GeometryModifyEngine* common_modify_engine( DLIList<RefEdge*>& edges,
                                              DLIList<Curve*>& curves,
                                              CubitBoolean allow_composites = CUBIT_FALSE ) const;
  GeometryModifyEngine* common_modify_engine( DLIList<RefVertex*>& vertices,
                                              DLIList<TBPoint*>& points,
                                              CubitBoolean allow_composites = CUBIT_FALSE ) const;

  GeometryModifyEngine* common_modify_engine( DLIList<TopologyEntity*>& topology_list,
                                              DLIList<TopologyBridge*>& engine_bridges,
                                              CubitBoolean allow_composites
                                              = CUBIT_FALSE ) const;
  /**<  \return GeometryModifyEngine*
    *   A GeometryModifyEngine common at least one
    *   TopologyBridge of each of the passed TopologyEntities, or
    *   NULL if no common geometry engine is found.
    *   \arg topology_list
    *   The input list of TopologyEntities
    *   \arg engine_bridges
    *   Pass back the list of TopolgyBridges associated with each
    *   of the passed TopologyEntities (topology_list) and owned
    *   by the returned geometry engine.
    *   \arg allow_virtual_engine
    *   Return VirtualGeometryEngine::instance() if no common
    *   geometry enginge can be found.
    *
    *   Look for a common geometry engine other than the
    *   VirtualGeometryEngine.  If no common geometry engine other
    *   than VGE can be found and allow_virtual_engine is FALSE,
    *   NULL is returned.  If allow_virtual_engine is TRUE, and no
    *   common geometry engine is found, VGE will be returned, and
    *   engine_bridges will be populated with any virtual geometry
    *   if possible, otherwise with the first topology bridge attached
    *   to each of the passed TopologyEntities.
    */

public:  

  //! \brief Internal use only. 
  CubitStatus finish_sm_op( DLIList<Body*>& input_bodies,
                            DLIList<BodySM*>& new_bodies,
                            DLIList<Body*>& result_bodies,
                            bool print_info = true ) const;
    //- Common code for completion of solid model engine operations
    //- on a set of bodies.
    //I input_bodies
    //I- A list of potentially modified or deleted bodies that were
    //I- were passed as input to the solid modeling operation.
    //I new_bodies
    //I- New bodies created during the solid modeling operation.
    //I result_bodies
    //I- New Bodys created for the BodySMs in new_bodies.
    //I update_input_bodies
    //I- Update DAG for non-deleted input_bodies
    //I delete_old_first
    //I- Clean out all deactivated geometry from input_bodies
    //I- before generating any new entities.  This has effect
    //I- on the resulting topology or IDs, but does affect the
    //I- names on new entities.

protected:

  GeometryModifyEngine* tweak_setup( DLIList<RefFace*>& input_faces,
                                     const char* tweak_function_name,
                                     DLIList<Body*>& old_bodies_out,
                                     DLIList<Surface*>& surfaces_out,
                                     CubitBoolean allow_composites = CUBIT_FALSE);

  GeometryModifyEngine* tweak_setup( DLIList<RefEdge*>& input_edges,
                                     const char* tweak_function_name,
                                     DLIList<Body*>& old_bodies_out,
                                     DLIList<Curve*>& curves_out );

  GeometryModifyEngine* tweak_setup( DLIList<RefVertex*> &input_vertices,
                                     const char* name,
                                     DLIList<Body*> &output_bodies,
                                     DLIList<TBPoint*> &output_points );

  CubitStatus sweep_setup ( const char* sweep_function_name,
                            DLIList<RefEntity*>& entity_list,
                            DLIList<Body*>& body_list,
                            GeometryModifyEngine*& engine,
                            CubitBoolean& changed_new_ids,
                            DLIList<GeometryEntity*>& geom_list,
                            bool keep_old,
                            DLIList<RefEdge*>* edge_list = 0,
                            DLIList<Curve*>* curve_list = 0 );
    //- Common setup code for sweep functions
    //I sweep_function_name
    //I- Name to use in error messages.
    //I entity_list
    //I- Input entity list to be swept
    //O body_list
    //O- Existing bodies that could potentially be modified.
    //O engine
    //O- Engine to call sweep_* on.
    //O change_new_ids
    //O- Save this value and pass into sweep_finish to make sure
    //O- state of 'new ids' is restored.
    //O geom_list
    //O- TopologyBridges to pass to the sweep fuction of the ModifyEngine.
    //I edge_list
    //I- RefEdges describing sweep path
    //O curve_list
    //O- Curves corresponding to edges in edge_list.

  CubitStatus imprint_singly( DLIList<Body*>& body_list,
                              DLIList<Body*>& new_bodies,
                              CubitBoolean keep_old );
    //- Implementation of imprint(..) when group_imprint is false.

  void remove_dead_entity_names( RefEntity* entity ) const;
    //- If the passed entity or any child entities are dead
    //- (have no TopologyBridges), remove entity names.

  Body* update_body( Body* body ) const;
    //- Destroy or update modified body, as appropriate.

private :

  //! \brief Static pointer to the unique instance of this class.
  static GeometryModifyTool* instance_;

  //! \brief Use group imprinting, which is much faster. Default: true.
  static CubitBoolean groupImprint;

  //! \brief An object will be unmeshed if modified or deleted. Default: true.
  static CubitBoolean meshAutodelete;

  //! \brief  An object will be cached if the mesh is autodeleted
  //!  and can be remeshed later.  Default: false.
  static CubitBoolean meshAutodeleteRemesh;

  //! \brief Causes code to reuse old ids, which are more persistent.  Default: false.
  static CubitBoolean newIds;

  //! \brief Attempt entity naming backwards compatible with v8/pre-bool-on-merged.  Default:
  //! false.
  static CubitBoolean oldNames;

  //! \brief Separate the result bodies of webcutting into single-body volumes. 
  //!  Default: true.
  static CubitBoolean sepAfterWebcut;

  /*!  Option for acis to know if we should imprint all of the edges
    *  or only relevant ones.  For now, default is false, but check
    *  the .cc file.  The user can control this by a "set" nonregimprint
    *  command.  Default: true.
    */
  static CubitBoolean allEdgesImprint;

  //! \brief Regularized booleans are performed.  Default: true.
  static CubitBoolean booleanRegularize;

  //! \brief Unite operations can mix sheet bodies with solid bodies.  Default: true.
  static CubitBoolean uniteMixedModels;

  /*! Pointer to the entity being copied.  When copying an entity, used to 
      prevent to copying of unique id and saved-merged-ids attributes. 
  */
  //! \brief This pointer points to the entity that is being copied.  
  static RefEntity *copyingEntity;

  //! \brief The list of geometry modify engines.
  DLIList<GeometryModifyEngine*> gmeList;

  void get_merged_curve_and_surface_ids( DLIList<Body*> &bodies,
                                         DLIList<int> &merged_surface_ids,
                                         DLIList<int> &merged_curve_ids ) const;

  CubitStatus separate_body_after_webcut (DLIList<Body*> &input_list,
                                          DLIList<Body*> &output_list) const;
  /**<  Separates each body in the input list to have one volume per body,
    *   and places the new list of bodies in the output_list. Should
    *   only be called after webcutting. (Checks sepAfterWebcut flag.)
    */

  void fixup_merged_entities( DLIList<int> &merged_surface_ids,
                              DLIList<int> &merged_curve_ids ) const;

  bool contains_intermediate_geom(DLIList<Body*>& list) const;
  bool contains_intermediate_geom(DLIList<TopologyBridge*>& list) const;
  bool contains_composites(DLIList<TopologyBridge*>& bridge_list ) const;
  bool contains_partitions(DLIList<TopologyBridge*>& bridge_list ) const;
  bool contains_partitions( DLIList<Body*>& list ) const;
  bool contains_composites( DLIList<Body*>& list ) const;

  void do_attribute_setup(void);
  void do_attribute_cleanup(void);
  void push_attributes_before_modify(DLIList<BodySM*> &old_sms);
  CubitStatus restore_vg_after_modify(DLIList<BodySM*> &new_sms,
                                      DLIList<Body*> &old_bodies,
                                      GeometryModifyEngine *gme);
  void remove_pushed_attributes(DLIList<BodySM*> &new_sms,
                                      DLIList<Body*> &old_bodies);
  void push_imprint_attributes_before_modify(DLIList<BodySM*> &old_sms);
  void push_named_attributes_to_curves_and_points(DLIList<TopologyBridge*> &tb_list, const char *name_in);
  void remove_imprint_attributes_after_modify(DLIList<BodySM*> &body_sms,
                                              DLIList<BodySM*> &new_sms);

  //! draw a preview of a plane for webcut previews
  static void plane_preview(DLIList<Body*> &body_list,
                            const CubitVector &pt1,
                            const CubitVector &pt2,
                            const CubitVector &pt3);

  void propagate_from_small_edge( RefEdge *edge,
                                  DLIList<RefEdge*> &small_edges,
                                  DLIList<RefFace*> &narrow_faces,
                                  DLIList<RefFace*> &processed_faces,
                                  double small_edge_length);
  void propagate_over_narrow_face( RefFace *narrow_face,
                                   RefEdge *edge,
                                   DLIList<RefFace*> &processed_faces,
                                   DLIList<RefEdge*> &small_edges,
                                   DLIList<RefFace*> &narrow_faces,
                                   double small_edge_length);

  void propagate_merge_tolerance(DLIList<Body*> &body_list);
  void push_tolerance_attribute(DLIList<Body*> &body_list);

  static CubitStatus prepare_for_copy( RefEntity *ref_ents,
                                       TopologyBridge *&top_bridge );

  CubitStatus prepare_bc_for_webcut();
  CubitStatus finish_bc_webcut();

  static CubitStatus finish_copy( TopologyBridge *&new_bridge,
                                  TopologyBridge *old_bridge );

  CubitStatus sweep_finish( const char* const sweep_function_name,
                            DLIList<Body*>& input_body_list,
                            DLIList<BodySM*>& new_body_list,
                            DLIList<Body*>& output_body_list,
                            CubitBoolean restore_newids );
    //- Common cleanup code for sweep functions
    //I sweep_function_name
    //I- Name to use in error messages.
    //I input_body_list
    //I- Bodies that may need to be updated.
    //I new_body_list
    //I- SMbodies returned by GeometryModifyEngine
    //I output_body_list
    //I- ref_bodies returned by GeometryModifyEngine
    //I restore_newids
    //I- Passed back from sweep_setup -- restore setting state.

  CubitStatus finish_webcut( DLIList<Body*>& input_bodies,
                             DLIList<BodySM*>& webcut_results,
                             CubitBoolean merge,
                             CubitStatus webcut_status,
                             DLIList<Body*>& new_bodies,
                             DLIList<int> *merged_surface_ids = NULL,
                             DLIList<int> *merged_curve_ids = NULL,
                             CubitBoolean print_info = CUBIT_TRUE ) const;
    //- Helper function for all webcut functions.
    //- Finish up DAG update, merging, etc. after GME does webcut.
    //I input_bodies
    //I- The list of bodies that were webcut
    //I webcut_results
    //I- The new bodies returned by the GME webcut.
    //I merge
    //I- Merge after webcut
    //I webcut_status
    //I- The return value from the GME webcut function
    //O new_bodies
    //O- The new Bodies constructed from the passed BodySMs.
    //R CubitStatus
    //R- CUBIT_FAILURE on error, webcut_status on success.
    //- Does separate_body_after_webcut for all bodies if
    //- webcut_status == CUBIT_SUCCESS
    //- Merges bodies if webcut_status == CUBIT_SUCCESS AND merge == CUBIT_TRUE

  CubitStatus tweak_remove( DLIList<RefFace*> &ref_face_list,
                            DLIList<Body*> &new_body_list,
                            CubitBoolean extend_adjoining = CUBIT_TRUE,
                            CubitBoolean keep_surface = CUBIT_FALSE,
                            CubitBoolean keep_old_body = CUBIT_FALSE,
                            CubitBoolean preview = CUBIT_FALSE);
  /**<  Remove surfaces from a body or bodies and then extend the adjoining
    *   surfaces to fill the gap or remove the hole.  Only called by
    *   tweak_remove_individually, tweak_remove_together.
    */

  void get_neighboring_bodies( DLIList<Body*> &input_bodies,
                               DLIList<Body*> &neighboring_bodies );


  CubitStatus unite_separately( GeometryModifyEngine *gme_ptr,
                                DLIList<Body*> &body_list,
                                DLIList<Body*> &new_body_list,
                                bool keep_old );
  CubitStatus unite_all( GeometryModifyEngine *gme_ptr,
                         DLIList<Body*> &body_list,
                         DLIList<Body*> &new_body_list,
                         bool keep_old );
  /**<  Private functions called to unite bodies.  First form separates the
    *   body list into sheets & solids and unites each group together 
    *   separately.  Second one unites all bodies together (including sheets
    *   & solids.
    */

  CubitStatus unite_private( GeometryModifyEngine *gme_ptr,
                             DLIList<Body*> &body_list,
                             DLIList<Body*> &new_body_list,
                             bool keep_old );
  /**<  Private lower level function called to unite bodies
    */


  static CubitStatus clean_up_from_copy_failure( TopologyBridge *old_bridge );

  GeometryModifyEngine * pull_common_surfs( DLIList<RefFace*> &ref_face_list,
                                            DLIList<RefFace*> &common_face_list,
                                            DLIList<Surface*> &common_surf_list );
    /**< Pull RefFaces with a common GeometryModifyEngine out of the input
      *  ref_face_list.  Place their surfaces in the output RefFace and Surface
      *  lists, and return the common modify engine.  Note the function returns
      *  a NULL pointer if a RefFace without a modify engine is found in the
      *  input list.
      */

};

inline void GeometryModifyTool::set_all_edges_imprint( CubitBoolean flag )
{allEdgesImprint = flag;}
inline CubitBoolean GeometryModifyTool::get_all_edges_imprint()
{return allEdgesImprint;}

inline CubitBoolean GeometryModifyTool::boolean_regularize()
{return booleanRegularize;}
inline void GeometryModifyTool::boolean_regularize(CubitBoolean flag)
{booleanRegularize = flag;}

inline CubitBoolean GeometryModifyTool::unite_mixed_models()
{return uniteMixedModels;}
inline void GeometryModifyTool::unite_mixed_models(CubitBoolean flag)
{uniteMixedModels = flag;}

inline void GeometryModifyTool::add_gme(GeometryModifyEngine *gme_ptr)
{
  assert(gme_ptr != 0);
  if (!gmeList.move_to(gme_ptr)) gmeList.append(gme_ptr);
}
  /**< add a geometry query engine to the list
   */

inline void GeometryModifyTool::get_gme_list(DLIList<GeometryModifyEngine*> &gme_list)
{
  gme_list += gmeList;
}
  /**< return the list of gme's
   */

inline GeometryModifyEngine *GeometryModifyTool::get_gme() const
{
  GeometryModifyEngine *gme = NULL;
  if (gmeList.size()) {
    const_cast<GeometryModifyTool*>(this)->gmeList.reset();
    gme = gmeList.get();
  }
  return gme;
}

#endif

