/**
 * \file GeometryQueryTool.hpp
 *
 * \brief This class provides the interface for all the query-only geometry
 *        related operations to the outside.
 *
 * This class is implemented using the Singleton pattern.
 *
 * \author Tim Tautges
 *
 * \date 2/2001
 *
 */


#ifndef GEOMETRYQUERYTOOL_HPP
#define GEOMETRYQUERYTOOL_HPP

#define CUBIT_12 3

#include <stdio.h>
#include <typeinfo>
#include <list>
#include <set>
#if !defined(WIN32)
using std::type_info;
#endif

#include "GeometryDefines.h"
#include "DLIList.hpp"
#include "GeometryQueryEngine.hpp"
#include "IntermediateGeomEngine.hpp"
#include "CGMHistory.hpp"

class RefGroup;
class Body;
class BodySM;
class Lump;
class Shell;
class ShellSM;
class Loop;
class LoopSM;
class Chain;
class Surface;
class Curve;
class TopologyBridge;
class CoEdgeSM;

class GeometryEntity;

class RefGroup;

class TopologyEntity;
class TopologyBridge;

class RefEntity ;
class RefVolume ;
class RefFace ;
class RefEdge ;
class CoEdge;
class RefVertex ;

class SenseEntity;
class BasicTopologyEntity;
class RefVolume ;
class RefFace ;
class RefEdge ;
class RefVertex ;
class TBPoint;
class CubitEntity;
class CoEdgeSM;

class GeometryQueryEngine;

class CubitPlane ;
class CubitBox;
class CubitVector;
class CubitTransformMatrix;
class TBOwner;

#ifdef PROE
class RefPart;
class RefAssembly;
#endif

//! Interface class for querying geometry. 
class CUBIT_GEOM_EXPORT GeometryQueryTool
{
public :

  friend class GeometryModifyTool;
   
  bool ige_is_composite(TBOwner *bridge_owner);
  bool ige_is_composite(TopologyBridge *bridge);
  bool ige_is_partition(TBOwner *bridge_owner);


  void ige_remove_modified(DLIList<Surface*> &all_surfs,
                                            DLIList<Curve*> &all_curves,
                                            DLIList<TBPoint*> &all_points);

  static GeometryQueryTool* instance( GeometryQueryEngine* gqePtr = NULL);
  /**<
    *  \return GeometryQueryTool* - Pointer to the singleton GeometryQueryTool object
    *  \arg SMEPtr
    *   Pointer to a GeometryQueryEngine object. The default value
    *   is set to NULL.
    *
    *  Return a pointer to the only instance of the class with the default
    *  geometric query engine set to the argument, if it is not NULL.
    *  If the argument is NULL, return the pointer to the existing
    *  instance without modifying the default SME. A valid GQEPtr *must*
    *  be supplied at the time of the first call to this function.
    *  Hence, this instance function should specifically be called at
    *  startup with a valid non-NULL input GQEPtr.
    */

  ~GeometryQueryTool();
  ///<  Destructor.

  static void delete_instance();
  
  //!
  //! \brief Estimates a good merge tolerance for the volumes passed in.
  //!
  double estimate_merge_tolerance(DLIList<RefVolume*> &vol_list,
                                                    bool accurate_in = false,
                                                    bool report_in = false,
                                                    double lo_val_in = -1.0, 
                                                    double hi_val_in = -1.0,
                                                    int num_calculations_in = 10,
                                                    bool return_calculations_in = false,
                                                    DLIList<double> *merge_tols = NULL,
                                                    DLIList<int> *num_proximities = NULL);

  //!
  //! \brief Find all of the volumes that do not contain any merged entities.
  //!
  void find_floating_volumes(DLIList<RefVolume*> &vol_list,
                             DLIList<RefVolume*> &floating_list);

  //!
  //! \brief Find the nonmanifold curves in the passed-in volumes based on what is merged.
  //!
  void find_nonmanifold_curves(DLIList<RefVolume*> &vol_list, DLIList<RefEdge*> &curve_list);

  //!
  //! \brief Find the nonmanifold vertices in the passed-in volumes based on what is merged.
  //!
  void find_nonmanifold_vertices(DLIList<RefVolume*> &vol_list, DLIList<RefVertex*> &vertex_list);

  CubitStatus register_intermediate_engine( IntermediateGeomEngine* engine );
  void unregister_intermediate_engine( IntermediateGeomEngine* engine );

  void ige_remove_imprint_attributes_after_modify(DLIList<BodySM*> &old_sms,
                                              DLIList<BodySM*> &new_sms);
  void ige_push_imprint_attributes_before_modify
                                ( DLIList<BodySM*> &geom_list );
  void ige_push_named_attributes_to_curves_and_points
                                ( DLIList<TopologyBridge*> &tb_list, const char *name_in );
  void ige_export_geom( DLIList<TopologyBridge*> &geom_list );
  void ige_import_geom( DLIList<TopologyBridge*> &geom_list );
  void ige_remove_attributes( DLIList<TopologyBridge*> &geom_list );
  void ige_attribute_after_imprinting( DLIList<TopologyBridge*> &new_tbs,
                                                    DLIList<TopologyBridge*> &att_tbs,
                                                    DLIList<TopologyBridge*> &tb_list,
                                                        DLIList<Body*> &old_bodies);
  void ige_remove_attributes_from_unmodifed_virtual(DLIList<TopologyBridge*> &bridges);
  
  //Using the source_bridge, finds all bridges that actually have a BridgeManager the owner.
  //This is for obtaining the real TopologyBridge when all you have is a  TopologyBridge 
  //that is Partition Entity.
  void get_tbs_with_bridge_manager_as_owner( TopologyBridge *source_bridge, 
                                             DLIList<TopologyBridge*> &tbs );

  /*! <HR><H1> global-list-functions Global entity list functions </H1>*/

  //
  // Returns ref entities in a generic cubit entity list (overwrites list).
  // returns CUBIT_FAILURE if keyword is not a ref entity name, and
  // optionally prints error message
  //
  //CubitStatus cubit_entity_list( const char* keyword,
  //                               DLIList<CubitEntity*> &entity_list,
  //                               const CubitBoolean print_errors = CUBIT_TRUE);

  //!
  //! \brief Return ref entities in entity_list (overwrites list).
  //! returns CUBIT_FAILURE if keyword is not a ref entity name,
  //! and optionally prints error message
  //!  

  CubitStatus ref_entity_list( char const* keyword,
                               DLIList<RefEntity*> &entity_list,
                               const CubitBoolean print_errors = CUBIT_TRUE);

  //! \brief Returns the bounding box of the model.  Include free entities.
  CubitBox model_bounding_box();

  /*! Returns the bounding box of all bodies in the model. Excludes
  free entities.*/
  //! \brief Returns the bounding box of all bodies in the model. 
  CubitBox bounding_box_of_bodies();


  //! \brief Returns all the bodies in the current session. 
  void bodies (DLIList<Body*> &bodies);

  //! \brief Returns all volumes in the current session. 
  void ref_volumes (DLIList<RefVolume*> &ref_volumes);

  //! \brief Returns all groups in the current session. 
  void ref_groups (DLIList<RefGroup*> &ref_groups);

  //! \brief Returns all surfaces in the current session. 
  void ref_faces (DLIList<RefFace*> &ref_faces);

  //! \brief Returns all curves in the current session. 
  void ref_edges (DLIList<RefEdge*> &ref_edges);

  //! \brief Returns all the vertices in the current session. 
  void ref_vertices (DLIList<RefVertex*> &ref_vertices);

#ifdef PROE
  void ref_parts (DLIList<RefPart*> &ref_parts);
  void ref_assemblies (DLIList<RefAssembly*> &ref_assemblies);
#endif //PROE
  ///< Append global lists to arguments

  //! \brief Number of bodies in current session.
  int num_bodies() const;

  //! \brief Number of volumes in current session.
  int num_ref_volumes() const;

  //! \brief Number of groups in current session.
  int num_ref_groups() const;

  //! \brief Number of surfaces in current session.
  int num_ref_faces() const;

  //! \brief Number of curves in current session.
  int num_ref_edges() const;

  //! \brief Number of vertices in current session.
  int num_ref_vertices() const;

  //! \brief Get RefEntity by type name and id.
  RefEntity *get_ref_entity (const char *type, int id);

  //! \brief Get a RefEntity of the specified type and id.
  RefEntity *get_ref_entity (const type_info& type, int id);

  //! \brief Get entity by id.
  Body *get_body ( int id );

  //! \brief Get entity by id.
  RefVolume *get_ref_volume ( int id );

  //! \brief Get entity by id.
  RefGroup *get_ref_group ( int id );

  //! \brief Get entity by id.
  RefFace *get_ref_face ( int id );

  //! \brief Get entity by id.
  RefEdge *get_ref_edge ( int id );

  //! \brief Get entity by id.
  RefVertex *get_ref_vertex ( int id );
  

  //! \brief Get the first entity in the global list of the specified type
  Body *get_first_body ();
  //! \brief Get the first entity in the global list of the specified type
  RefVolume *get_first_ref_volume ();
  //! \brief Get the first entity in the global list of the specified type
  RefGroup *get_first_ref_group ();
  //! \brief Get the first entity in the global list of the specified type
  RefFace *get_first_ref_face ();
  //! \brief Get the first entity in the global list of the specified type
  RefEdge *get_first_ref_edge ();
  //! \brief Get the first entity in the global list of the specified type
  RefVertex *get_first_ref_vertex ();

  //! \brief Get the next entity in the global list of the specified type
  Body *get_next_body ();
  //! \brief Get the next entity in the global list of the specified type
  RefVolume *get_next_ref_volume ();
  //! \brief Get the next entity in the global list of the specified type
  RefGroup *get_next_ref_group ();
  //! \brief Get the next entity in the global list of the specified type
  RefFace *get_next_ref_face ();
  //! \brief Get the next entity in the global list of the specified type
  RefEdge *get_next_ref_edge ();
  //! \brief Get the next entity in the global list of the specified type
  RefVertex *get_next_ref_vertex ();

  ///! \brief Get the last entity in the global list of the specified type
  Body *get_last_body ();
  ///! \brief Get the last entity in the global list of the specified type
  RefVolume *get_last_ref_volume ();
  ///! \brief Get the last entity in the global list of the specified type
  RefGroup *get_last_ref_group ();
  ///! \brief Get the last entity in the global list of the specified type
  RefFace *get_last_ref_face ();
  ///! \brief Get the last entity in the global list of the specified type
  RefEdge *get_last_ref_edge ();
  ///! \brief Get the last entity in the global list of the specified type
  RefVertex *get_last_ref_vertex ();
  
  
  //! \brief Get all free surfaces, curves, and vertices 
  CubitStatus get_free_ref_entities(DLIList<RefEntity*> &free_entities);

  //! \brief Get the free entities connected to, but not necessarily
  //! vertically related to, the entity.  - If merge_option, then take
  //! into account the fact that the model may be merged (it's slower
  //! that way).
  void get_connected_free_ref_entities(
    RefEntity *entity,
    const int merge_option,
    DLIList<Body*> &body_list,
    DLIList<RefFace*> &ref_face_list,
    DLIList<RefEdge*> &ref_edge_list,
    DLIList<RefVertex*> &ref_vertex_list );

  //! \brief Saves out a temporary geometry file containing specified entities
  //! that are of the same geometry engine. 
  CubitStatus save_temp_geom_files( DLIList<RefEntity*>& ref_entity_list,
                                  const char* filename,
                                  const CubitString &cubit_version,
                                  std::list<CubitString> &files_written,
                                  std::list<CubitString> &types_written);


  //!
  //! * Export entities to a solid model file.
  //!  *  \arg ref_entity_list
  //!  *  A list of RefEntities to be exported or saved to a file.
  //!  *  \arg file_name
  //!  *  The name of the file to write to.
  //!  *  \arg file_type
  //!  *  An optional type of file.
  //!  *  \arg logfile_name
  //!  *  Optional - name of logfile.
  //!  *   \return  CubitStatus - success/failure
  //!  *
  //!  *  Export the current CUBIT geometry (everything in the Model) to a
  //!  *  solid model format. Valid file types are:
  //!  *     "ACIS_SAT"    --  ACIS ASCII (SAT) file format
  //!  *     "ACIS_SAB"    --  ACIS BINARY (SAB) file format
  //!  *     "ACIS_DEBUG"  --  ACIS ASCII debug format
  //!  *     "IGES"        --  IGES file
  //!  *     "STEP"        --  STEP file
  //!  No logfile gets created for SAB/SAT files, but for IGES and
  //!  STEP file a logfile always gets created.  Default filenames
  //!  are assigned if one is not given (iges_export.log, step_export.log).
  //!  
  //!  The function returns CUBIT_FAILURE if anything goes wrong with
  //!  export - improper file type, inaccessible file, mismatch between
  //!  the underlying representation and file type. It returns
  //!  CUBIT_SUCCESS if everything goes well.
  //!  
  //!  NOTE: if the ref_entity_list is empty, GeometryQueryTool gets the list of
  //!  all entities in the current model, including free ref entities
  //! \brief Save a geometry file containing specified entities
  //! that are of the same geometry engine. 
  CubitStatus export_solid_model( DLIList<RefEntity*>& ref_entity_list,
                                  const char* filename,
                                  Model_File_Type filetype,
                                  int &num_ents_exported,
                                  const CubitString &cubit_version,
                                  ModelExportOptions &export_options );

  /*!
    Import all or specified entities in a solid model file.
    *  \arg file_ptr
    *  A pointer to the file to read (can be NULL for IGES and STEP files).
    *  \arg file_type
    *  Type of file.
    *  \arg heal_step - auto-healing of step bodies on import.  This is recommended
    *     because they always need it.
    *  \arg import_bodies (etc...)
    *  Should bodies be import.
    *   \return  CubitStatus - success/failure
    *
    *  Reads in geometry and creates the necessary Reference entities
    *  associated with the input geometry.
    *  Valid file types are:
    *     "ACIS_SAT"    --  ACIS ASCII (SAT) file format
    *     "ACIS_SAB"    --  ACIS BINARY (SAB) file format
    *     "IGES"        --  IGES file
    *     "STEP"        --  STEP file
    *
    *  Function can selectively import solid bodies, free surfaces, free
    *  curves, or free vertices.  For example, the user may not want
    *  to import any free entities.
    *
    *  The function returns CUBIT_FAILURE if anything goes wrong with
    *  import - improper file type, inaccessible file, mismatch between
    *  the underlying representation and file type. It returns
    *  CUBIT_SUCCESS if everything goes well.
    */
  //! \brief Import a geometry file. 
  CubitStatus import_solid_model(const char* file_name,
                                 Model_File_Type file_type,
                                 ModelImportOptions &import_options,
				 DLIList<RefEntity*> *imported_entities = NULL );

  /*!
   * Fire a ray at entities, passing back distances of hits and entities hit
    * \arg origin
    * origin of ray
    * \arg direction
    * direction of ray
    * \arg at_entity_list
    * entities to fire ray at
    * \arg ray_params
    * returned array of parameters (distances) along ray at which entities were hit
    * \arg max_hits
    * maximum number of hits to return, 0 = unlimited (default)
    * \arg ray_radius
    * radius of ray to use for intersecting entities, 0 = use engine default
    * \arg hit_entity_list (pointer)
    * entities hit by ray (list length same as ray_params), default NULL
    * \return - error flag
    *
    *  Fire a ray at specified entities, returning the parameters (distances)
    *  along the ray and optionally the entities hit; return CUBIT_FAILURE if
    *  error.  Returned lists are appended to.
    */

  //! \brief Fire a ray at entities, passing back distances of hits and entities hit
  CubitStatus fire_ray( CubitVector &origin,
                        CubitVector &direction,
                        DLIList<RefEntity*> &at_entity_list,
                        DLIList<double> &ray_params,
                        int max_hits = 0,
                        double ray_radius = 0.0,
                        DLIList<RefEntity*> *hit_entity_list_ptr = 0 );

  /*!
   * Fire a ray at entities, passing back distances of hits and entities hit
    * \arg origin
    * origin of ray
    * \arg direction
    * direction of ray
    * \arg at_entity_list
    * entities to fire ray at
    * \arg ray_params
    * returned array of parameters (distances) along ray at which entities were hit
    * \arg max_hits
    * maximum number of hits to return, 0 = unlimited (default)
    * \arg ray_radius
    * radius of ray to use for intersecting entities, 0 = use engine default
    * \arg hit_entity_list (pointer)
    * entities hit by ray (list length same as ray_params), default NULL
    * \return - error flag
    *
    *  Fire a ray at specified entities, returning the parameters (distances)
    *  along the ray and optionally the entities hit; return CUBIT_FAILURE if
    *  error.  Returned lists are appended to.  NOTE: ALL ENTITIES MUST BE FROM
    *  THE SAME GEOMETRY ENGINE. 
    */
  //! \brief Fire a ray at entities, passing back distances of hits and entities hit
  CubitStatus fire_ray( CubitVector &origin,
                        CubitVector &direction,
                        DLIList<TopologyEntity*> &at_entity_list,
                        DLIList<double> &ray_params,
                        int max_hits = 0,
                        double ray_radius = 0.0,
                        DLIList<TopologyEntity*> *hit_entity_list_ptr = 0 );
  
  //! \brief Debugging function.
  static void geom_debug( DLIList<TopologyEntity*> );

  //! \brief Set facet box flag
  static void set_facet_bbox( CubitBoolean pass_flag )
      {useFacetBBox = pass_flag;}

  //! \brief Get facet box flag
  static CubitBoolean get_facet_bbox()
      {return useFacetBBox;}

  //! \brief Calls engine version of the active geometry engine.
  CubitString get_engine_version_string();

  //! \brief Set the major/minor version of the active geometry engine.
  CubitStatus set_export_allint_version(int version);

  //! \brief Get the major/minor version of the active geometry engine.
  int get_allint_version();

  //! \brief Returns a string with the versions of the active geometry engine.
  CubitStatus list_engine_versions(CubitString &versions);

  //! \brief Gets solid modeler's resolution absolute tolerance
  double get_sme_resabs_tolerance();

  //! \brief Sets solid modeler's resolution absolute tolerance
  double set_sme_resabs_tolerance( double new_resabs );

  //! \brief Set solid modeler integer option.
  CubitStatus set_sme_int_option( const char* opt_name, int val );
  //! \brief Set solid modeler double option.
  CubitStatus set_sme_dbl_option( const char* opt_name, double val );
  //! \brief Set solid modeler string option.
  CubitStatus set_sme_str_option( const char* opt_name, const char* val );


  ///< <HR><H3>Topology/geometry creation functions</H3>
  Body* make_Body(BodySM *bodysm_ptr) const;
  RefFace* make_RefFace(Surface* surface_ptr ) const;
  RefEdge* make_RefEdge(Curve* curve_ptr) const;
  RefVertex* make_RefVertex(TBPoint* point_ptr) const;

  static CubitSense relative_sense( Surface* surface1, Surface* surface2 );

  RefFace* make_free_RefFace(Surface *surface_ptr, bool is_free_surface) const;
  RefEdge* make_free_RefEdge(Curve *curve_ptr ) const;
  RefVertex* make_free_RefVertex(TBPoint *point_ptr) const;
  /**<  These functions can be used to create free ref-entities
    *  from a geometry engine.  Just call populate_topology_
    *  bridges to create the sm_ptr of the desired type, then pass
    *  that to the appropriate function here.
    */


  ///<HR><H3>Topology and geometry deletion</H3>

  /*! \brief  Deletes all Bodies in the input list from the model. Their
    associated Solid Model entities are deleted as well, if they exist,
    and if remove_solid_model_entities is CUBIT_TRUE.
    */
  void delete_Body( DLIList<Body*>& body_list );

  /**<  Deletes the input Body from the model.
    *  Its associated Solid Model entities are deleted as well, if
    *  they exist, and if remove_solid_model_entities is CUBIT_TRUE.
    *  If all went well, the input Body* reference will be set to NULL
    *  as the Body itself has been deleted.
    *  Returns CUBIT_SUCCESS if all went well, otherwise, CUBIT_FAILURE.
    */
  //! \brief Deletes a body.
  CubitStatus delete_Body( Body* body_ptr );

  /*! Behaves exactly as delete_Body, but in addition checks to see if
    *  children of Body are merged.  In some cases, 2 entities can be
    *  forced-merged, where they are not spatially equal.  This regenerates
    *  the graphics on merged entities so they look correct after a partner
    *  has been deleted.
    */
  //! \brief Deletes a body.
  CubitStatus delete_single_Body( Body* body_ptr );

  //! \brief Deletes free RefEnties
  CubitStatus delete_RefEntity( RefEntity* ref_entity_ptr );
  
  //! \brief Deletes the RefFace if it is free
  CubitStatus delete_RefFace( RefFace* ref_face_ptr );

  //! \brief Deletes the RefEdge if it is free
  CubitStatus delete_RefEdge( RefEdge* ref_edge_ptr );

  //! \brief Deletes the RefVertex if it is free
  CubitStatus delete_RefVertex( RefVertex* ref_vertex_ptr );

  void cleanout_deactivated_geometry();
  void cleanout_temporary_geometry ();
  
  //! \brief Deletes all geometry.
  void delete_geometry();

  
  /*! \brief  Creates a list of vectors the length of the number_points that
      interpolate between vector_1 and vector_2.  All of the points
      will lie on the underlying equation of the refface.
  */
  CubitStatus interpolate_along_surface( CubitVector *vector_1,
                                         CubitVector *vector_2,
                                         DLIList<CubitVector*> &vector_list,
                                         RefFace* ref_face_ptr,
                                         int number_points ) const;
  /*!  \return  CubitBoolean
    \return - CUBIT_TRUE/CUBIT_FALSE
    \arg Vec1
    A reference to the first vector.
    \arg Vec2
    A reference to the second vector.
    \arg tolerance_factor
    Factor by which the absolute internal tolerance shall be
    multiplied.
    *  Returns CUBIT_TRUE if the input Vec1 and Vec2 are spatially
    equivalent within a tolerance.  The internal spatial tolerance
    value is multiplied by tolerance_factor before the (spatial)
    test is done. Else, returns CUBIT_FALSE.
    */
  //! \brief Compares two positions for coincidence.
  CubitBoolean about_spatially_equal (const CubitVector& Vec1,
                                      const CubitVector& Vec2,
                                      double tolerance_factor = 1.0);

  //! \brief Compares two vertices for coincidence.
  CubitBoolean about_spatially_equal (RefVertex* refVertex1,
                                      RefVertex* refVertex2,
                                      double tolerance_factor = 1.0);

  /*! \brief Calculates internal surface angles given 2 refedges on the
   * surface.  CoEdge version correctly handles curves in a surface
   * twice. */
  double geometric_angle(RefEdge* ref_edge_1,
                         RefEdge* ref_edge_2,
                         RefFace* ref_face );

  /*! \brief Calculates internal surface angles given 2 refedges on the
   surface.  This version correctly handles curves in a surface twice. */
  double geometric_angle(CoEdge* co_edge_1,
                         CoEdge* co_edge_2 );

  /*! \brief Calculate dihedral angle at curve between two surfaces of the
    volume. */
  double surface_angle( RefFace *ref_face_1, RefFace *ref_face_2,
                        RefEdge *ref_edge = NULL,
                        RefVolume *ref_volume = NULL,
                        double frac = 0.5);

  /*!  Finds the intersections between a curve and a line  If the bounded flag 
    *  is true, it finds only those intersections within the parameter range
    *  of the curve; otherwise it uses the extensions of the curve.  The
    *  closest option is currently valid only if the curve is straight,
    *  in which case it will return the 2 closest intersection locations,
    *  if the straight lines don't actually intersect. The function allocates 
    *  allocates the CubitVectors in the returned list, so be sure to free them.
    */
  //! \brief Finds the intersections of a straight line and a curve.  
  CubitStatus get_intersections( RefEdge* ref_edge1,
                                 CubitVector& point1,
                                 CubitVector& point2,
                                 DLIList<CubitVector*>& intersection_list,
                                 CubitBoolean bounded = CUBIT_FALSE,
                                 CubitBoolean closest = CUBIT_FALSE);
 
  /*!  Finds the intersections of the two curves.  If the bounded flag is
    *  true, it finds only those intersections within the parameter range
    *  of both curves; otherwise it uses the extensions of the curves.  The
    *  closest option is currently valid only if both curves are straight,
    *  in which case it will return the 2 closest intersection locations,
    *  if the straight lines don't actually intersect.  So far, other than
    *  for straight lines, this function only works if both curves are ACIS
    *  curves, unless both curves are linear.  The function allocates the
    *  CubitVectors in the returned list, so be sure to free them.
    */
  //! \brief Finds the intersections of two curves.  
  CubitStatus get_intersections( RefEdge* ref_edge1, RefEdge* ref_edge2,
                                 DLIList<CubitVector*>& intersection_list,
                                 CubitBoolean bounded = CUBIT_FALSE,
                                 CubitBoolean closest = CUBIT_FALSE );

  /*! Finds the intersections of the curve and surface.  The curve is extended
    * if the bounded flag is not false.  The function allocates the CubitVectors
    * in the returned list, so be sure to free them.
  */
  //! \brief Finds the intersections of the curve and surface.  
  CubitStatus get_intersections( RefEdge* ref_edge, RefFace* ref_face,
                                 DLIList<CubitVector*>& intersection_list,
                                 CubitBoolean bounded = CUBIT_FALSE );

  //! Gets the intersection of a curve a plane.  The extended_percent
  //! extends the plane by a percentage value.
  CubitStatus get_intersections( RefEdge* ref_edge, CubitPlane plane,
                                 DLIList<CubitVector*>& intersection_list,
                                 CubitBoolean bounded = CUBIT_FALSE,
                                 double extended_percent = 0.0);


  /*! Gets the extrema position along the first given direction. If there
    * is more than one extrema position, the other directions will be used
    * to determine a unique position.  Directions 2 and 3 can be NULL.
    * Entities supported include bodies, volumes, surfaces, curves and
    * vertices.  
    */
  //! \brief Gets extrema position on an entity.
  CubitStatus entity_extrema( RefEntity *ref_entity_ptr,
                              const CubitVector *dir1,
                              const CubitVector *dir2,
                              const CubitVector *dir3,
                              CubitVector &extrema,
                              RefEntity *&extrema_entity_ptr );

  /*! Gets the extrema position along the first given direction. If there
    * is more than one extrema position, the other directions will be used
    * to determine a unique position.  Directions 2 and 3 can be NULL.
    * Entities supported include bodies, volumes, surfaces, curves and
    * vertices.  The entity the extrema is found on is also returned.
    */
  //! \brief Gets extrema position on a list of entities.
  CubitStatus entity_extrema( DLIList<RefEntity*> &ref_entity_list,
                              const CubitVector *dir1,
                              const CubitVector *dir2,
                              const CubitVector *dir3,
                              CubitVector &extrema,
                              RefEntity *&extrema_entity_ptr );

  /*! Gets the minimum distance between two entities and the closest positions
    * on those entities. Supports vertices, curves, surfaces, volumes and bodies.
  */
  //! \brief Get the minimum distance between two entities.
  CubitStatus entity_entity_distance( GeometryEntity *ge1,
                                      GeometryEntity *ge2,
                                      CubitVector &pos1, CubitVector &pos2,
                                      double &distance );

  /*! Gets the minimum distance between two entities and the closest positions
    on those entities. Supports vertices, curves, surfaces, volumes and bodies.
  */
  //! \brief Get the minimum distance between two entities. (CGM internal use)
  CubitStatus entity_entity_distance( RefEntity *ref_entity_ptr1,
                                      RefEntity *ref_entity_ptr2,
                                      CubitVector &pos1, CubitVector &pos2,
                                      double &distance );


  //! \brief Resets geometry factor back to 500.0
  static void reset_geometry_factor();
  
  //! \brief Sets geometry factor. 
  static void set_geometry_factor( double fac );

  //! \brief Gets geometry factor. 
  static double get_geometry_factor();

  //! \brief Sets bboxMergeTest variable.
  static void set_merge_test_bbox(CubitBoolean tof);

  //! \brief Gets bboxMergeTest variable.
  static CubitBoolean get_merge_test_bbox();

  static void set_merge_test_internal(int tof);
  static int get_merge_test_internal();
  static void set_sliver_curve_cleanup_tolerance( double tol );
  static void set_sliver_surface_cleanup_tolerance( double tol );
  static double get_sliver_curve_cleanup_tolerance(); 
  static double get_sliver_surface_cleanup_tolerance(); 

  //! \brief Initializes all settings of this class
  static void initialize_settings();

  //! \brief Causes attributes hanging on entities to be applied.
  static CubitStatus import_actuate( DLIList<RefEntity*> &entity_list );

  /*! \brief Returns CUBIT_TRUE if all the entities have the same geometric query engine and
    if that is the same one as the default. */
  CubitBoolean same_query_engine(DLIList<TopologyEntity*> &topo_list) const;

  /*!  \return GeometryQueryEngine*
    *   A GeometryQueryEngine common at least one
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
  //! \brief Gets geometry beloning to a common modeling engine. 
  GeometryQueryEngine* common_query_engine(
    DLIList<TopologyEntity*>& topology_list,
    DLIList<TopologyBridge*>& engine_bridges,
    CubitBoolean allow_default_engine
    = CUBIT_FALSE ) const;

  //! \brief Determine if any of the input entities contain the given query engine
  CubitBoolean does_geom_contain_query_engine(DLIList<TopologyEntity*> &topo_list,
                                              GeometryQueryEngine *engine) const;

  //! \brief Determine if any of the input entities contain the given query engine
  CubitBoolean does_geom_contain_query_engine(DLIList<RefEntity*> &ref_entity_list,
                                              GeometryQueryEngine *engine,
                                              CubitBoolean children_too = CUBIT_FALSE) const;

  //! \brief Retrieves the TopologyEntity from the underlying TopologyBridge
  TopologyEntity* entity_from_bridge( TopologyBridge* bridge_ptr ) const;

  //! \brief Adds a geometry query engine to the list
  void add_gqe(GeometryQueryEngine *gqe_ptr);

  /*! \brief Removes a geometry query engine from the list. Returns CUBIT_FAILURE
     if it wasn't in the list */
  CubitStatus remove_gqe(GeometryQueryEngine *gqe_ptr);

  //! \brief Return the list of GeometryQureyEngines.
  void get_gqe_list(DLIList<GeometryQueryEngine*> &gqe_list);

  //! \brief Returns the first gqe on the list.
  GeometryQueryEngine *get_gqe();

  //! \brief Set the default GeometryQueryEngine.
  CubitStatus set_default_gqe(GeometryQueryEngine* gqe);

  bool contains_intermediate_geometry(RefEntity*) const;
  bool contains_intermediate_geometry(DLIList<RefEntity*>& ref_entitylist) const;
  bool is_intermediate_geometry(RefEntity*) const;
  bool is_intermediate_geometry(TopologyBridge*) const;

  /*!- Remove this entity and any dead children where
  a dead child a) has no parent entities and b) has no topology_bridges. */
  CubitStatus destroy_dead_entity( TopologyEntity* topo_ent, bool top = true ) const;

  //! \brief Translate a Body some delta.
  CubitStatus translate( Body* entity, const CubitVector& delta, bool check_to_transform = true,
                         bool preview = false );

  //! \brief Translate a BasicTopologyEntity some delta.
  CubitStatus translate( BasicTopologyEntity* entity, const CubitVector& delta, bool check_to_transform = true,
                         bool preview = false );

  void translate( DLIList<RefEntity*> &entities_to_transform,
          double x, double y, double z, bool check_before_transforming,
          DLIList<RefEntity*> &entities_transformed,
          bool preview = false );


  //! \brief Rotate a Body an angle about an axis. 
  CubitStatus rotate( Body* entity, const CubitVector& axis, double degrees, bool check_to_transform = true,
                      bool preview = false );

  //! \brief Rotate a Body an angle about an axis, defined by a point and a 
  //! direction.
  CubitStatus rotate( Body* entity,
                      const CubitVector& point,
                      const CubitVector& normal,
                      double degrees,
                      bool check_to_transform = true,
                      bool preview = false);

  CubitStatus rotate( DLIList<RefEntity*> &entities_to_transform,  
                      const CubitVector& point,
                      const CubitVector& normal,
                      double degrees,
                      bool check_to_transform, 
                      DLIList<RefEntity*> &entities_transformed,
                      bool preview = false);

  //! \brief Rotate a BacisTopologyEntity an angle about an axis.
  CubitStatus rotate( BasicTopologyEntity* entity, 
                      const CubitVector& axis, 
                      double degrees,
                      bool check_to_transform = true,
                      bool preview = false);
  CubitStatus rotate( BasicTopologyEntity* entity, 
                      const CubitVector& point,
                      const CubitVector& normal,
                      double degrees,
                      bool check_to_transform = true,
                      bool preview = false);

  //! \brief Scale a Body.
  CubitStatus scale( Body* entity,const CubitVector& point, double factor, bool check_to_transform = true, bool preview = false);

  //! \brief Scale a Body different factors in x, y, and z.
  CubitStatus scale( Body* entity,const CubitVector& point, const CubitVector& factors, bool check_to_transform = true, bool preview = false);

  //! \brief Scale a BasicTopologyEntity. 
  CubitStatus scale( BasicTopologyEntity* entity,const CubitVector& point, double factor, bool check_to_transform = true, bool preview = false);

  //! \brief Scale a BasicTopologyEntity different factors in x, y, and z.
  CubitStatus scale( BasicTopologyEntity* entity,const CubitVector& point, const CubitVector& factors,bool check_to_transform = true,
                     bool preview = false);

  void scale( DLIList<RefEntity*> &entities_to_transform, 
              const CubitVector& point,
              double scale_x, double scale_y, double scale_z, 
              bool check_to_transform, 
              DLIList<RefEntity*> &entities_scaled,
              bool preview = false);

  //! \brief Reflect a list of bodies about a plane defined by an axis.
  CubitStatus reflect( DLIList<Body*> bodies,const CubitVector& point, const CubitVector& axis, bool preview = false );

  //! \brief Reflect a BasicTopologyEntity about a plane defined by an axis.
  CubitStatus reflect( BasicTopologyEntity* entity, 
                       const CubitVector& point,
                       const CubitVector& axis,
                       bool check_to_transform = true,
                       bool preview = false);

  void reflect( DLIList<RefEntity*> &entities_to_transform,
                            const CubitVector& point,
                            const CubitVector& axis,
                            bool check_before_transforming,
                            DLIList<RefEntity*> &entities_transformed,
                            bool preview = false);

  //! \brief Need to deprecate.
  CubitStatus restore_transform( Body* body );

  /*! Query to determine if volumes intersect, share common volume.
  Returns CUBIT_TRUE if the two volumes overlap, CUBIT_FALSE if they don't
  overlap.  If the volumes are touching the function should return CUBIT_FALSE.
  volume_1, 
  volume_2
  The two volume pointers that are being tested for overlap.
  The function uses the intersect call to test if the volumes 
  are overlaping.  The full intersect Boolean is needed to see if
  the volumes actually overlap and don't just touch. */
  //! \brief Query to determine if volumes intersect, share common volume.
  CubitBoolean volumes_overlap( RefVolume *volume_1, RefVolume *volume_2);

  /*! Query to determine if bodies intersect, share common volume.
  Returns CUBIT_TRUE if the two bodies overlap, CUBIT_FALSE if they don't
  overlap.  If the bodies are touching the function should return CUBIT_FALSE.
  body_ptr_1, body_ptr_2
  The two body pointers that are being tested for overlap.
  The function uses the intersect call to test if the bodies
  are overlaping.  The full intersect Boolean is needed to see if
  the bodies actually overlap and don't just touch. */
  //! \brief Query to determine if bodies intersect, share common volume.
  CubitBoolean bodies_overlap( Body *body_ptr_1, Body *body_ptr_2 );

  //! Given a list of TB's, construct ref entities for them; if the 2nd list pointer is
  //! non-NULL, pass back the list of ref entities in that list
  CubitStatus construct_refentities(DLIList<TopologyBridge*> &topology_bridges,
                                    DLIList<RefEntity*> *imported_entities = NULL);

  /*! When importing a cub file, embedded in the cub file is how many 
      geometry entities it is supposed to restore.  If geometry that 
      you are improrting is merged with geometry already in the session, 
      you need to keep track of how many geometry entieies get 'merged-away' 
      like this so that import does not fail.  Taking into account 
      the 'merged-away' geometry allows CUBIT to successfully import 
      when you have geometry that will merge-away. */
  /*! \brief Variable needed when importing geometry that will be 
    merged-away with already existing geometry in the session */
  static CubitBoolean trackMergedAwayEnts;
  static CubitBoolean importingSolidModel;

  //Before calling import_solid_model, normally we want to clear out the 
  //map in CAUniqueId, but not when we are importing a cub file.
  static CubitBoolean clearUidMapBeforeImport; 

  //variable defining scope of merge...if true, 
  //merge entities you are importing with any entity (i.e.
  //entities already in the cubit session), otherwise 
  //only merge entities that are importing with other entities 
  //that are importing 
  static CubitBoolean mergeGloballyOnImport;

  static DLIList<int> uidsOfImportingEnts;
  static int entitiesMergedAway; 

  CGMHistory& history();

  CubitStatus get_graphics( Body *body, 
                            GMem *g_mem,
                            std::vector<RefFace*> &face_to_facet_vector,
                            std::vector<RefEntity*> &facet_point_ownership_vector,
                            std::vector<std::pair<RefEntity*, std::pair<int,int> > > &facetedges_on_refedges,
                            unsigned short normal_tolerance, 
                            double distance_tolerance, 
                            double max_edge_length );

  CubitStatus get_graphics( RefFace *ref_face,
                            GMem *gmem,
                            std::vector<RefEntity*> &facet_point_ownership_vector,
                            std::vector<std::pair< RefEntity*, std::pair<int,int> > > &facetedges_on_refedges,
                            unsigned short normal_tolerance = 15, 
                            double distance_tolerance = 0.0, 
                            double max_edge_length = 0.0 );

protected :

  GeometryQueryTool(GeometryQueryEngine* GQEPtr);
  /**<  Constructor for the (singleton) GeometryQueryTool object
    */

private :

  CubitBoolean okay_to_transform( Body* body ) const;
  CubitBoolean okay_to_transform( BasicTopologyEntity* bte ) const;
  CubitStatus notify_observers_of_transform( RefEntity* ref_entity ) const;
  CubitStatus notify_intermediate_of_transform( TopologyEntity* bte,
                                             const CubitTransformMatrix& xform
                                             ) const;

  Shell* make_Shell(ShellSM *shellsm_ptr, bool& shell_modified ) const;
  RefVolume* make_RefVolume( Lump* lump_ptr, bool& vol_modified ) const;

  CubitStatus make_merged_RefFace( Surface* surface_ptr ) const;
    // Helper function for make_RefFace(Surface*)
    // For a merged RefFace, check if the passed surface topology
    // is still compatible.  If so, remerge any unmerged child
    // entities and return CUBIT_SUCCESS.  If not, return CUBIT_FAILURE.

  RefEntity *check_mergeable_refentity(GeometryEntity *bridge) const;
    // check for mergeable ref entity, indicated by a merge attribute on the
    // bridge

  CoEdgeSM *find_merged_coedgesm( Surface* on_this_surface,
                                  CoEdgeSM* merged_with_this ) const;
    // check for mergeable sense entity, which is a child of the entity pointed
    // to by a merge attribute on bridge, which corresponds to re_ptr

  
  void get_merged_away_free_entities( DLIList<RefEntity*> &ref_ents,
                                      DLIList<TopologyBridge*> &free_ents );


  CubitStatus straightline_intersections( RefEdge* ref_edge1,
                                 CubitVector & origin2,
                                 CubitVector & dir2,
                                 DLIList<CubitVector*>& intersection_list,
                                 CubitBoolean bounded = CUBIT_FALSE,
                                 CubitBoolean closest = CUBIT_FALSE);

  //! static pointer to the unique instance of this class.
  static GeometryQueryTool* instance_;

  //! For use in calculating a bounding box, you can do it based on
  //! a set of facets, rather than what default modeling engine  uses.
  static CubitBoolean useFacetBBox;

  //! The list of geometry query engines
  DLIList<GeometryQueryEngine*> gqeList;

  struct IGEComp : public std::binary_function<IntermediateGeomEngine*,
                                              IntermediateGeomEngine*,
                                              bool> {
    bool operator() ( const IntermediateGeomEngine* ptr1,
                      const IntermediateGeomEngine* ptr2 ) const
      { return ptr1->level() < ptr2->level(); }
  };
  typedef std::set<IntermediateGeomEngine*,IGEComp> IGESet;
  IGESet igeSet;

  //! \brief The default geometry query engine
  GeometryQueryEngine *default_gqe;

  //! After imprinting, an attempt at removing sliver curves is made.
  //! Curves less than this tolerance will be removed.
  static double curveSliverCleanUpTolerance;

  //! After imprinting, an attempt at removing sliver surfaces is made.
  //! Removes sliver surfaces whose maximum gap distance among the long 
  //! edges is smaller than the tolerance and who have at most three long edges. 
  static double surfaceSliverCleanUpTolerance;

  //!  This factor is the the multiplier for the ACIS resabs
  //!  when comparingcurves.  ALWAYS when using a multiplier
  //!  use this factor for consistancy.
  static double geometryToleranceFactor;

  static CubitBoolean bboxMergeTest;

  //! Options for refface merging.
  //! 0=off, 1=all, 2=splines only
  static int internalSurfaceMergeTest; 

  int maxPersistentBodyId;
  int maxPersistentRefVolumeId;
  int maxPersistentRefGroupId;
  int maxPersistentRefFaceId;
  int maxPersistentRefEdgeId;
  int maxPersistentRefVertexId;

  static const int CGM_MAJOR_VERSION;
  static const int CGM_MINOR_VERSION;

  CGMHistory mHistory;
};

inline void GeometryQueryTool::reset_geometry_factor()
{
  geometryToleranceFactor = DEFAULT_GEOM_FACTOR;
}

inline void GeometryQueryTool::set_geometry_factor( double factor )
{
  if ( factor < .0099999999999 )
    return;
  else
    geometryToleranceFactor = factor;
}
inline double GeometryQueryTool::get_geometry_factor()
{
  return geometryToleranceFactor;
}

inline void GeometryQueryTool::set_sliver_surface_cleanup_tolerance( double tol )
{
  surfaceSliverCleanUpTolerance = tol;
}

inline void GeometryQueryTool::set_sliver_curve_cleanup_tolerance( double tol )
{
  curveSliverCleanUpTolerance = tol;
}

inline double GeometryQueryTool::get_sliver_curve_cleanup_tolerance()
{
  return curveSliverCleanUpTolerance; 
}

inline double GeometryQueryTool::get_sliver_surface_cleanup_tolerance()
{
  return surfaceSliverCleanUpTolerance; 
}

inline CubitBoolean GeometryQueryTool::get_merge_test_bbox()
{return bboxMergeTest;}
inline int GeometryQueryTool::get_merge_test_internal()
{return internalSurfaceMergeTest;}
inline void GeometryQueryTool::set_merge_test_bbox(CubitBoolean tof)
{bboxMergeTest = tof;}
inline void GeometryQueryTool::set_merge_test_internal(int tof)
{internalSurfaceMergeTest = tof;}

inline void GeometryQueryTool::add_gqe(GeometryQueryEngine *gqe_ptr)
{
  assert(gqe_ptr != 0);
  if (!gqeList.move_to(gqe_ptr)) gqeList.append(gqe_ptr);
}
  /**< add a geometry query engine to the list
   */

inline CubitStatus GeometryQueryTool::remove_gqe(GeometryQueryEngine *gqe_ptr)
{
  assert(gqe_ptr != 0);
  CubitStatus status = CUBIT_FAILURE;
  if (gqeList.move_to(gqe_ptr)) {
    gqeList.remove();
    status = CUBIT_SUCCESS;
  }
  return status;
}

inline void GeometryQueryTool::get_gqe_list(DLIList<GeometryQueryEngine*> &gqe_list)
{
  gqe_list += gqeList;
}
  /**< return the list of gqe's
   */

inline GeometryQueryEngine *GeometryQueryTool::get_gqe()
{
  GeometryQueryEngine *gqe = NULL;
  if (gqeList.size()) {
    gqeList.reset();
    gqe = gqeList.get();
  }
  return gqe;
}

  /**< Return the VirtualGeometryEngine if it exists.  Otherwise return null.
   */

#endif



