//-------------------------------------------------------------------------
// Filename      : AcisQueryEngine.hpp
//
// Purpose       : ACIS version of the geometry query engine (provides geometry
//                 query functions at the solid model
//                 level)
//
// Special Notes :
//
// Creator       : 
//
// Creation Date : 
//
// Owner         : 
//-------------------------------------------------------------------------

#ifndef ACIS_QUERY_ENGINE_HPP
#define ACIS_QUERY_ENGINE_HPP

// ********** BEGIN STANDARD INCLUDES         **********

#include <typeinfo>

#if !defined(NT)
using std::type_info;
#endif

// ********** END STANDARD INCLUDES           **********

// ********** BEGIN ACIS INCLUDES             **********
// ********** END ACIS INCLUDES               **********

// ********** BEGIN CUBIT INCLUDES            **********

#include "GeometryQueryEngine.hpp"
#include "AcisTypes.h"

// ********** END CUBIT INCLUDES              **********

// ********** BEGIN FORWARD DECLARATIONS
class ENTITY;
class ENTITY_LIST;
class BODY;
class SHELL;
class LUMP;
class SURFACE;
class VERTEX;
class WIRE;
class FACE;
class EDGE;
class COEDGE;
class LOOP;
class SPAposition;
class SPAtransf;
class surface;
class outcome;
class PLANE;
class CONE;
class SPHERE;
class TORUS;
class SPLINE;
class SPAbox;

class TopologyBridge;

class GeometryEntity;
class BodySM;
class Lump;
class ShellSM;
class Surface;
class LoopSM;
class CoEdgeSM;
class Curve;
class Point;

class CubitBox;
class CubitString;

class BodyACIS;
class LumpACIS;
class ShellACIS;
class SurfaceACIS;
class LoopACIS;
class CoEdgeACIS;
class CurveACIS;
class PointACIS;

template <class X> class DLIList;
class FACE;
class LOOP;
class EDGE;
class VERTEX;
class FeatureCut;
#ifdef ACIS_HEALER
class AcisHealerTool;
#endif
//#ifdef CUBIT_GUI
//class GUITestingDlg; // Testing only (remove when done)
//#endif

#ifdef ACIS_STEP_TRANSLATOR
class step_header;
#endif

class GMem;
class AcisFacetManager;

// ********** END FORWARD DECLARATIONS        **********

// ********** BEGIN MACRO DEFINITIONS         **********
// ********** END MACRO DEFINITIONS           **********

// ********** BEGIN ENUM DEFINITIONS          **********
// ********** END ENUM DEFINITIONS            **********

class AcisQueryEngine : public GeometryQueryEngine
{
public:
// ********** BEGIN FRIEND DECLARATIONS        **********
  friend class AcisBridge;
  friend class PointACIS;
  friend class CurveACIS;
  friend class CoEdgeACIS;
  friend class LoopACIS;
  friend class SurfaceACIS;
  friend class ShellACIS;
  friend class LumpACIS;
  friend class BodyACIS;
  friend class ATTRIB_CUBIT_OWNER;
  friend class FeatureCut;
  friend class AcisModifyEngine;
  friend class AcisHealerTool;
  friend class AcisSurfaceTool;
  friend class AcisTweakTool;
  friend class AcisEdgeTool;
  friend class AcisToolUtil;
  friend class AcisDrawTool;
//#ifdef CUBIT_GUI
//  friend class GUITestingDlg;
//#endif
// ********** END FRIEND DECLARATIONS        **********

//HEADER- Constructor and Destructor
  
  virtual ~AcisQueryEngine();
  
private:

  AcisQueryEngine();

  static AcisQueryEngine* instance_;

public:
  inline static AcisQueryEngine* instance()
  {
    if( ! instance_ )
      instance_ = new AcisQueryEngine;
    return instance_;
  }

  const char * modeler_type()
     { return "ACIS"; }

  
//HEADER- RTTI and safe casting functions.
  
  virtual const type_info& entity_type_info() const ;
    //R- The geometric modeler type
    //- This function returns the type of the geometric modeler.
  
  
//HEADER- Functions for importing and exporting solid models.
  
  virtual CubitStatus export_solid_model(
    DLIList<TopologyBridge*>& bridge_list,
    const char* file_name,
    const char* file_type,
    const CubitString &cubit_version,
    const char* logfile_name = NULL ) ;
    //R CubitStatus
    //R- CUBIT_SUCCESS/CUBIT_FAILURE
    //I ref_entity_list
    //I- A list of RefEntities to be exported or saved to a file.
    //I file_name
    //I- The name of the file to write to.
    //I file_type
    //I- An optional type of file. 
    //- Export the current CUBIT geometry (everything in the Model) to a 
    //- solid model format. Valid file types are:
    //-    "ACIS_SAT"    --  ACIS ASCII (SAT) file format
    //-    "ACIS_SAB"    --  ACIS BINARY (SAB) file format
    //-    "ACIS_DEBUG"  --  ACIS ASCII debug format
    //-    "IGES"        --  IGES file
    //-    "STEP"        --  STEP file
    //I logfile_name
    //I- Optional - name of logfile.  
    //- No logfile gets created for SAB/SAT files, but for IGES and
    //- STEP file a logfile always gets created.  Default filenames
    //- are assigned if one is not given (iges_export.log, step_export.log).
    //-
    //- The function returns CUBIT_FAILURE if anything goes wrong with 
    //- export - improper file type, inaccessible file, mismatch between
    //- the underlying representation and file type. It returns
    //- CUBIT_SUCCESS if everything goes well.
 
  virtual CubitStatus save_temp_geom_file(DLIList<TopologyBridge*> &ref_entity_list,
                                          const char *filename,
                                          const CubitString &cubit_version,
                                          CubitString &created_file,
                                          CubitString &created_file_type );

  virtual CubitStatus import_temp_geom_file(FILE* file_ptr, 
                                    const char* file_name,
                                    const char* file_type,
                                    DLIList<TopologyBridge*> &bridge_list );

  virtual CubitStatus import_solid_model(
    const char* file_name,
    const char* file_type,
    DLIList<TopologyBridge*>& imported_entities,
    CubitBoolean print_results = CUBIT_TRUE,
    const char* logfile_name = NULL,
    CubitBoolean heal_step = CUBIT_TRUE,
    CubitBoolean import_bodies = CUBIT_TRUE,
    CubitBoolean import_surfaces = CUBIT_TRUE,
    CubitBoolean import_curves = CUBIT_TRUE,
    CubitBoolean import_vertices = CUBIT_TRUE,
    CubitBoolean free_surfaces = CUBIT_TRUE ) ;
    //R CubitStatus
    //R- CUBIT_SUCCESS/CUBIT_FAILURE
    //I file_ptr
    //I- A pointer to the file to read (can be NULL for IGES and STEP files).
    //I file_type
    //I- Type of file. 
    //- Reads in geometry and creates the necessary Reference entities 
    //- associated with the input geometry. 
    //- Valid file types are:
    //-    "ACIS_SAT"    --  ACIS ASCII (SAT) file format
    //-    "ACIS_SAB"    --  ACIS BINARY (SAB) file format
    //-    "IGES"        --  IGES file
    //-    "STEP"        --  STEP file
    //I heal_step - auto-healing of step bodies on import.  This is recommended
    //              because they always need it.
    //I import_bodies (etc...)
    //I- Should bodies be import.
    //- Function can selectively import solid bodies, free surfaces, free
    //- curves, or free vertices.  For example, the user may not want
    //- to import any free entities.
    //-
    //- The function returns CUBIT_FAILURE if anything goes wrong with 
    //- import - improper file type, inaccessible file, mismatch between
    //- the underlying representation and file type. It returns
    //- CUBIT_SUCCESS if everything goes well.
  
  CubitStatus restore_entity_list(ENTITY_LIST &entity_list,
                                  DLIList<TopologyBridge*> &bridge_list,
                                  CubitBoolean print_results = CUBIT_TRUE,
                                  CubitBoolean step_import = CUBIT_FALSE,
                                  CubitBoolean iges_import = CUBIT_FALSE,
                                  CubitBoolean import_bodies = CUBIT_TRUE,
                                  CubitBoolean import_surfaces = CUBIT_TRUE,
                                  CubitBoolean import_curves = CUBIT_TRUE,
                                  CubitBoolean import_vertices = CUBIT_TRUE,
                                  CubitBoolean free_surfaces = CUBIT_TRUE) const;
  
//HEADER- Functions for copying and removing solid model entities
  
  virtual void delete_solid_model_entities(DLIList<BodySM*>& body_list) const ;
    //- Deletes all solid model entities associated with the Bodies in 
    //- the input list. 
  
  virtual CubitStatus delete_solid_model_entities(GeometryEntity* ptr,
                                                  bool delete_lower_order) const;
  virtual CubitStatus delete_solid_model_entities( BodySM* body_ptr ) const;
  virtual CubitStatus delete_solid_model_entities(Surface* surf_ptr ) const;
  virtual CubitStatus delete_solid_model_entities( Curve* curve_ptr ) const;
  virtual CubitStatus delete_solid_model_entities( Point* point_ptr ) const;
    //- Delete specified bridges and associated ACIS topology.
  
//HEADER- Functions related to faceting
  
  virtual CubitStatus get_graphics(Surface* surface_ptr,
                                          int& number_triangles,
                                          int& number_points,
                                          int& number_facets,
                                          GMem* gMem,
                                          unsigned short normal_tolerance,
                                          double distance_tolerance,
                                          double max_edge_length ) const;
  CubitStatus get_graphics( FACE* face_ptr,
                                  int& number_triangles,
                                  int& number_points,
                                  int& number_facets,
                                  GMem* gMem,
                                  unsigned short normal_tolerance,
                                  double distance_tolerance,
                                  double max_edge_length = 0 ) const;

    //R CubitStatus
    //R- CUBIT_SUCCESS/CUBIT_FAILURE
    //I ref_face_ptr
    //I- The RefFace for which hoops facetting information will be
    //I- gathered.
    //I deviationPercentage
    //I- The facetter surface deviation percentage.
    //O number_triangles
    //O- The number of polygons (triangles) needed for facetting.
    //O number_points
    //O- The number of points needed for facetting.
    //O number_facets
    //O- The number of facets needed for facetting.
    //O gMem
    //O- The storage place for facets (and points).
    //- This function gathers and outputs ACIS facet (and point)
    //- information for hoops involved in facetting a RefFace.  If
    //- all goes well, CUBIT_SUCCESS is returned.  Otherwise,
    //- CUBIT_FAILURE is returned.

  virtual CubitStatus get_graphics( Curve* curve_ptr,
                                    int& num_points,
                                    GMem* gMem = NULL,
                                    double tolerance = 0.0 ) const;
    //R CubitStatus
    //R- CUBIT_SUCCESS/CUBIT_FAILURE
    //I ref_edge_ptr
    //I- The RefEdge for which hoops facetting information will be
    //I- gathered.
    //O numSteps
    //O- The number of points in resulting polyline.
    //O gMem
    //O- The storage place for edges, involved in facetting.
    //I tolerance
    //I- The tolerance deviation used when facetting the curve (optional).
    //- This function gathers and outputs ACIS edge information for
    //- hoops involved in facetting a RefEdge.  If all goes well,
    //- CUBIT_SUCCESS is returned.  Otherwise, CUBIT_FAILURE is
    //- returned.
  
  virtual CubitStatus get_isoparametric_points(Surface* ref_face_ptr,
                                               int &nu, int &nv,
                                               GMem *&gMem) const;
  
  virtual CubitStatus get_u_isoparametric_points(Surface* ref_face_ptr,
                                                 double v, int& n,
                                                 GMem *&gMem) const;
  
  virtual CubitStatus get_v_isoparametric_points(Surface* ref_face_ptr,
                                                 double u, int&n,
                                                 GMem *&gMem) const;
  
//HEADER- Miscellaneous functions
  
  virtual CubitStatus get_intersections(
                                 Curve* ref_edge1,
                                 CubitVector& point1,
                                 CubitVector& point2,
                                 DLIList<CubitVector*>& intersection_list,
                                 bool bounded = false,
                                 bool closest = false );
    //- Finds the intersections of the a curve and a straight line defined
    //- by the two points. If the bounded flag is
    //- true, it finds only those intersections within the parameter range
    //- of both curves; otherwise it uses the extensions of the curves.  The
    //- closest option is currently valid only if both curves are straight,
    //- in which case it will return the 2 closest intersection locations,
    //- if the straight lines don't actually intersect. The function allocates
    //- the CubitVectors in the returned list, so be sure to free them.
    //- NOTE: closest option has not been implemented in AGE, it's handled
    //        in GeometryTool (so it errors out if this option is found).

  virtual CubitStatus get_intersections( 
                                 Curve* ref_edge1, 
                                 Curve* ref_edge2,
                                 DLIList<CubitVector*>& intersection_list,
                                 bool bounded = false,
                                 bool closest = false );
    //- Finds the intersections of the two curves.  If the bounded flag is
    //- true, it finds only those intersections within the parameter range
    //- of both curves; otherwise it uses the extensions of the curves.  The
    //- closest option is currently valid only if both curves are straight,
    //- in which case it will return the 2 closest intersection locations,
    //- if the straight lines don't actually intersect. The function allocates 
    //- the CubitVectors in the returned list, so be sure to free them.
    //- NOTE: closest option has not been implemented in AGE, it's handled
    //        in GeometryTool (so it errors out if this option is found).

  virtual CubitStatus get_intersections( 
                                    Curve* ref_edge, 
                                    Surface* ref_face,
                                    DLIList<CubitVector*>& intersection_list,
                                    bool bounded = false );
    //- Finds the intersections of the curve and surface.  The curve is extended
    //- if the bounded flag is not false.  The function allocates the CubitVectors 
    //- in the returned list, so be sure to free them.

  virtual CubitStatus entity_extrema( DLIList<GeometryEntity*> &entity_list, 
                                      const CubitVector *dir1, 
                                      const CubitVector *dir2,
                                      const CubitVector *dir3, 
                                      CubitVector &extrema,
                                      GeometryEntity *&extrema_entity_ptr );
    //- Gets the extrema position along the first given direction. If there 
    //- is more than one extrema position, the other directions will be used 
    //- to determine a unique position.  Directions 2 and 3 can be NULL.
    //- Entities supported include bodies, volumes, surfaces, curves and
    //- vertices.  The entity the extrema is found on is also returned.

  virtual CubitStatus entity_entity_distance( GeometryEntity *ref_entity_ptr1,
                                              GeometryEntity *ref_entity_ptr2,
                                              CubitVector &pos1, 
                                              CubitVector &pos2,
                                              double &distance );
    //- Gets the minimum distance between two entities and the closest positions 
    //- on those entities. Supports vertices, curves, surfaces, volumes and bodies.

  virtual CubitStatus fire_ray( BodySM *body,
                                const CubitVector &ray_point,
                                const CubitVector &unit,
                                DLIList<double>& ray_params,
                                DLIList<GeometryEntity*> *entity_list) const;
    //- Fire a ray at the specified body, returning the entities hit and
    // the parameters along the ray; return non-zero if error

  virtual CubitStatus transform_vec_position( 
    CubitVector const& position_vector,
    BodySM *OSME_ptr,
    CubitVector &transform_vector ) const;
    //R CubitStatus         
    //R- CUBIT_SUCCESS/CUBIT_FAILURE
    //I position_vector on body, body that was transformed
    //O vector showning transform.
    //O-Computes the transform from the transformed body.
  
  virtual int get_major_version();
  virtual int get_minor_version();
  virtual int get_subminor_version();
  virtual int get_allint_version();
  virtual CubitString get_engine_version_string();

    //- pass a string back identifying this query engine
  virtual CubitStatus set_export_version(const int major,
                                         const int minor);
  virtual CubitStatus set_export_allint_version(const int version);
  virtual CubitStatus list_engine_versions(CubitString &versions);
    //- gets, sets, and lists choices for acis save version(s)

  virtual double get_sme_resabs_tolerance() const; 
  virtual double set_sme_resabs_tolerance( double new_resabs );
  // Get/set solid modeler's resolution absolute tolerance

  virtual CubitStatus set_int_option( const char* opt_name, int val );
  virtual CubitStatus set_dbl_option( const char* opt_name, double val );
  virtual CubitStatus set_str_option( const char* opt_name, const char* val );
   //- Set solid modeler options
  
  CubitStatus get_surfs_on_plane( BodySM* body_ptr, 
                                  const CubitVector& pln_orig, 
                                  const CubitVector& pln_norm,
                                  DLIList<Surface*>& ref_face_list ) const;
   
#ifdef ACIS_STEP_TRANSLATOR
  static void step_import_header_data_callback(step_header&);
    //- Used to get data from the step header during read process
#endif

  void delete_ACIS_ENTITY (ENTITY_LIST& to_be_deleted) const;
    //- Deletes the ACIS ENTITYs in the input ENTITY_LIST. The function starts
    //- by removing all lower order ENTITYs related to higher order ENTITYs in the
    //- in the list.  This is done to avoid getting an api error message when deleting
    //- lower order ENTITYs that have been marked for deletion.
  
  CubitStatus delete_ACIS_ENTITY (ENTITY* inputAcisEntity) const;
    //- Deletes the input ACIS ENTITY.  It also removes all references from the
    //- RefEntity datastructure to the input ACIS ENTIY.

  void delete_ACIS_BODY (ENTITY_LIST& BODY_list) const;
    //- Deletes the ACIS BODYs in the input BODY_list. It also removes all 
    //- references (from the RefEntity datastructure) to ACIS ENTITYs belonging 
    //- to the BODY being deleted. Note that pointers to only FACEs, EDGEs 
    //- and VERTEXes are removed from the RefEntity datastructure as these 
    //- are the only non-manifold entities.
    //- NOTE:
    //- Only BODYs are allowed in the input ENTITY_LIST.
  
  CubitStatus delete_ACIS_BODY (BODY* inputAcisBody,
                                CubitBoolean skip_unhook = CUBIT_FALSE) const ;
    //- Carefully deletes the input ACIS BODY. This procedure also removes
    //- all references (from the RefEntity datastructure) to ACIS ENTITYs 
    //- belonging to inputAcisBody, if they exist. Note that 
    //- pointers to only FACEs, EDGEs and VERTEXes are removed from the 
    //- RefEntity datastructure as these are the only non-manifold entities.

  CubitStatus transform( BodySM* body, const SPAtransf& xform );
  CubitStatus transform( BODY *body, const SPAtransf& xform );
  CubitStatus transform( GeometryEntity* ent, const SPAtransf& xform );
  
  CubitStatus translate( BodySM* body, const CubitVector& offset );
  CubitStatus rotate   ( BodySM* body, const CubitVector& axis, double angle );
  CubitStatus scale    ( BodySM* body, double factor );
  CubitStatus scale    ( BodySM* body, const CubitVector& factors );
  CubitStatus reflect  ( BodySM* body, const CubitVector& axis );
  CubitStatus restore_transform( BodySM* body );
  
  CubitStatus translate( GeometryEntity* ent, const CubitVector& offset );
  CubitStatus rotate   ( GeometryEntity* ent, const CubitVector& axis, double angle );
  CubitStatus scale    ( GeometryEntity* ent, double factor );
  CubitStatus scale    ( GeometryEntity* ent, const CubitVector& factors );
  CubitStatus reflect  ( GeometryEntity* ent, const CubitVector& axis );
  
  
  static void listcat( ENTITY_LIST& result, const ENTITY_LIST& other );
    // Append non-deleted entities in "other" to "result".
  
  virtual CubitBoolean volumes_overlap( Lump *lump1,
                                        Lump *lump2) const;

  virtual CubitBoolean bodies_overlap( BodySM *body_ptr_1,
                                       BodySM *body_ptr_2 ) const;
    //R CubitBoolean
    //R- CUBIT_TRUE if the two bodies overlap, CUBIT_FALSE if they don't
    //R- overlap.  If the bodies are touching the function
    //R- should return CUBIT_FALSE.
    //I body_ptr_1, body_ptr_2
    //I- The two body pointers that are being tested for overlap.
    //-  The function uses the intersect call to test if the bodies
    //-  are overlaping.  The full intersect Boolean is needed to see if
    //-  the bodies actually overlap and don't just touch.
  

protected:
  
private:
// HEADER- Private functions to read various file types
  CubitStatus read_acis_file(const char* file_name, const char* file_type, ENTITY_LIST &entity_list);
  CubitStatus read_acis_file(FILE* file_ptr, const char* file_type, ENTITY_LIST &entity_list);
#ifdef ACIS_IGES_TRANSLATOR
  CubitStatus read_iges_file(const char* file_name, const char* logfile_name, 
                             ENTITY_LIST &entity_list, bool heal );
#endif
#ifdef ACIS_CATIA_TRANSLATOR
  CubitStatus read_catia_file(const char* file_name, const char* logfile_name, ENTITY_LIST &entity_list);
#endif
#ifdef ACIS_STEP_TRANSLATOR
  CubitStatus read_step_file(const char* file_name, const char* logfile_name, 
                             ENTITY_LIST &entity_list, bool heal );
  void read_step_part_names(ENTITY_LIST &entity_list);
  void auto_separate_step_body(ENTITY_LIST &pre_separate_list, ENTITY_LIST &post_separate_list);
#endif
#ifdef ACIS_PROE_TRANSLATOR
  CubitStatus read_proe_file(const char* file_name, ENTITY_LIST &entity_list);
#endif

//HEADER- Functions that create TopologyEntities
  BodySM *populate_topology_bridges(BODY *body_ptr) const;
  Lump *populate_topology_bridges(LUMP *lump_ptr) const;
  Surface *populate_topology_bridges(FACE *face_ptr) const;
  Curve *populate_topology_bridges(EDGE *edge_ptr) const;
  Point *populate_topology_bridges(VERTEX *vertex_ptr) const;
    //- go through the entity passed in, make sure it has TopologyBridge
    //- entities (or create them)
  
//HEADER- Topology and geometry query functions
  
#ifdef BOYD14
  void get_EDGEs_of_FACE(FACE* FACEPtr, DLIList<EDGE*>& EDGEList) const ;
    //- Given a FACE, return all its EDGEs
#endif
  
  CubitStatus get_EDGEs_of_Curves( DLIList<Curve*>& ref_edge_list, 
                                   DLIList<EDGE*>& EDGE_list) const ;
    //- Given a list of RefEdges, return all the EDGEs
  
//HEADER- Functions that create ACIS ENTITYs
  
  int save_ENTITY_as_sat_file (ENTITY* entity_ptr, 
                               const char* filename,
                               const char* update_mode) const;
    //- Save the ENTITY, entity_ptr, to the file, filename, in the ACIS
    //- .sat ASCII file format. Does not look for or append the "sat" filetype,
    //- but uses the input string, filename, as is. The update_mode is
    //- the second argument in the "fopen" C/C++ library function.
  
//HEADER- Functions related to creating and deleting BODYs
  
//HEADER- Functions that unhook the double links between VGI
//HEADER- entities and ACIS entities.
  
  void unhook_BODY_from_refentities (BODY* BODY_ptr) const;
    //- Traverses the topology of BODY_ptr and for each ENTITY, calls
    //- remove_ENTITY_from_refentities. In other words, this
    //- call removes all references (from within the RefEntity datastructure)
    //- to all the ACIS ENTITYs in the input ACIS BODY, BODY_ptr.
  
  void remove_ENTITY_from_refentities (ENTITY* acisEntity) const;
    //- Removes all references to the input ACIS ENTITY, acisEntity, from the
    //- entire RefEntity datastructure in Model.
  
  void unhook_ENTITY_from_VGI(ENTITY* ENTITY_ptr)  const;
  void unhook_ENTITY_from_VGI(BODY* BODY_ptr)  const;
  void unhook_ENTITY_from_VGI(LUMP* LUMP_ptr)  const;
  void unhook_ENTITY_from_VGI(SHELL* SHELL_ptr)  const;
  void unhook_ENTITY_from_VGI(FACE* FACE_ptr)  const;
  void unhook_ENTITY_from_VGI(LOOP* LOOP_ptr)  const;
  void unhook_ENTITY_from_VGI(COEDGE* COEDGE_ptr)  const;
  void unhook_ENTITY_from_VGI(EDGE* EDGE_ptr, CubitBoolean remove_lower_entities = CUBIT_TRUE)  const;
  void unhook_ENTITY_from_VGI(VERTEX* VERTEX_ptr)  const;
    //- Unhook the ENTITY's references in VGI. Call the next level
    //- of the series of functions to allow other ENTITYs do their
    //- own unhooking.
  
//HEADER- Functions related to geometric comparison operations
  
  CubitBoolean about_spatially_equal ( const SPAposition &pos1, 
                                       const SPAposition &pos2, 
                                       double tolerance_factor = 1.0 ) const;
    //- Returns CUBIT_TRUE if the input pos1 and pos2 are spatially 
    //- equivalent within a tolerance.  The internal spatial tolerance value is
    //- multiplied by tolerance_factor before the (spatial) test is done.
    //- Else, returns CUBIT_FALSE.
  
  CubitBoolean about_spatially_equal ( VERTEX* VERTEX1, 
                                       VERTEX* VERTEX2, 
                                       double tolerance_factor = 1.0 ) const;
    //- Returns CUBIT_TRUE if the input VERTEX1 and VERTEX2 are spatially 
    //- equivalent within a tolerance.  The internal spatial tolerance value is
    //- multiplied by tolerance_factor before the (spatial) test is done.
    //- Else, returns CUBIT_FALSE
  
//HEADER- Functions for geometry and topology queries.
  
#ifdef BOYD14
  int get_position ( EDGE* EDGE_ptr,
                     double fraction_of_param_range,
                     SPAposition& pos ) const;
    //- Returns the position (pos) at the parameter value represented by
    //- the lowest parameter value + fraction_of_param_range*parameter_range.
    //- Note that this is NOT an arc-length based position -- it is based
    //- on the local parameterization of the underlying curve.
    //- Returns CUBIT_TRUE if all went well; CUBIT_FALSE, otherwise.
#endif
  
  static BODY*   get_BODY  (BodySM * bodysm_ptr );
  static LUMP*   get_LUMP  (Lump   * lump_ptr   );
  static FACE*   get_FACE  (Surface* surface_ptr);
  static EDGE*   get_EDGE  (Curve  * curve_ptr  );
  static VERTEX* get_VERTEX(Point  * point_ptr  );
  
  BodySM* get_body_sm_of_ENTITY(ENTITY* ENTITY_ptr) const ;
  BODY* get_BODY_of_ENTITY(ENTITY* ENTITY_ptr) const ;
    //- get the Body and BODY associated with this ENTITY
  
  BODY* get_BODY_of_ENTITY(BODY* BODY_ptr) const ;
  BODY* get_BODY_of_ENTITY(LUMP* LUMP_ptr) const ;
  BODY* get_BODY_of_ENTITY(WIRE* WIRE_ptr) const ;
  BODY* get_BODY_of_ENTITY(SHELL* SHELL_ptr) const ;
  BODY* get_BODY_of_ENTITY(FACE* FACE_ptr) const ;
  BODY* get_BODY_of_ENTITY(LOOP* LOOP_ptr) const ;
  BODY* get_BODY_of_ENTITY(COEDGE* COEDGE_ptr) const ;
  BODY* get_BODY_of_ENTITY(EDGE* EDGE_ptr) const ;
  BODY* get_BODY_of_ENTITY(VERTEX* VERTEX_ptr) const ;
    //- Return a pointer to the BODY that "owns" the input ENTITY, ENTITY_ptr.
  
  BODY* get_BODY_of_entity(Lump* refvolume_ptr) const ;
    //- Return the BODY that "owns" the LUMP associated with the input 
    //- RefVolume.

  BODY* get_BODY_of_entity(Surface* ref_face_ptr) const ;
    //- Return the BODY that "owns" the FACE associated with the input 
    //- RefFace.

  BODY* get_BODY_of_entity(Curve* ref_edge_ptr) const ;
    //- Return the BODY that "owns" the EDGE associated with the input 
    //- RefEdge.

  BODY* get_BODY_of_entity(Point* ref_vertex_ptr) const ;
    //- Return the BODY that "owns" the VERTEX associated with the input 
    //- RefVertex.

  ENTITY* get_ENTITY_of_entity( TopologyBridge *ref_entity_ptr ) const;
    //- Returns the ENTITY of the given RefEntity
  
  ENTITY* get_ENTITY_of_entity( BODY* BODY_ptr,
                                TopologyBridge* topoentity_ptr) const;
    //- Returns the ENTITY in the input BODY_ptr whose parent attribute 
    //- corresponds to the input TopologyEntity pointer
  
//HEADER- Bounding box functions.
  
  CubitStatus create_super_acis_bounding_box( ENTITY_LIST& entity_list,
                                              SPAbox &super_box ) const ;
    //- Returns an ACIS box object that is the bounding box of the list
    //- of ACIS ENTITIYs, entity_list.
  
  CubitStatus create_super_acis_bounding_box( DLIList<BodySM*>& entity_list,
                                              SPAbox &super_box ) const ;
    //- Returns an ACIS box object that is the bounding box of the list 
    //- of Bodies, body_list.
  
  SPAbox get_acis_entity_bounding_box( ENTITY* entity_ptr) const ;
    //- Returns the bounding box of the input ACIS ENTITY. 
    //- MJP NOTE: 
    //- This member function works only for a limited subset of those ENTITYs $
    //- that support bounding boxes.
  
  double get_max_size_of_acis_box ( const SPAbox& acis_box ) const ;
    //- Returns the length of the largest side of the input ACIS box.
  
  void clear_bounding_box(BODY* BODYPtr) const ;
  void clear_bounding_box(LUMP* LUMPPtr) const ;
  void clear_bounding_box(SHELL* SHELLPtr) const ;
  void clear_bounding_box(FACE* FACEPtr) const ;
  void clear_bounding_box(LOOP* LOOPPtr) const ;
  void clear_bounding_box(EDGE* EDGEPtr) const ;
    //R void
    //I BODYPtr, etc...
    //I- Pointer to a BODY (etc...) whose bounding box and those of the 
    //I- underlying ENTITYs are to be cleared.
    //- This function clears(i.e. sets to NULL) the bounding boxes of
    //- a BODY (etc...) and all the other underlying ENTITYs.

  static CubitBox bounding_box (SPAbox const & box);
  static SPAbox bounding_box(const CubitBox& box);
  CubitBox bounding_box (BodySM* body) const;
    //- Generate a CubitBox from an ACIS box
  
  SPAbox bounding_box (BODY* body) const;
  SPAbox bounding_box (LUMP* lump) const;
  SPAbox bounding_box (FACE* face) const;
  SPAbox bounding_box (EDGE* edge) const;
    //- Return the bounding box for the input ACIS ENTITY.

  CubitStatus bounding_box_from_facets(FACE *face_ptr, 
      SPAposition &bbox_min, SPAposition &bbox_max) const ;
    // - Return the bounding box for the input ACIS ENTITY based on
    // - its graphics facet representation
  
//HEADER- Attribute-related functions
  void get_all_cubit_owners(BODY *BODY_ptr, 
                            DLIList<TopologyBridge*> &tb_list) const;
    //- get all the te's owned by this entity and its children
  
  void remove_cubit_owner_attrib_in_BODY(BODY* BODY_ptr) const ;
    //- Removes the cubit owner attributes from all ENTITYs in the input BODY.
    //- The BODY topology is traversed and cleansed of the cubit owner attribute
    
//HEADER- Debug and error-checking functions.
  
  void ACIS_API_error ( outcome result, char const* message = NULL ) const;
    //- handles the api error messages
  
  void acis_debug_body ( BODY *local_body ) const ;
  void acis_debug_face ( FACE *local_face ) const ;
  void acis_debug_edge ( EDGE *local_edge ) const ;
  void acis_debug_coedge ( COEDGE *local_coedge ) const ;
  void acis_debug_loop ( LOOP *local_loop ) const ;
    //- dumps an acis debug_entity output to standard out for the ENTITY
  
  CubitStatus facet_EDGE(const EDGE* EDGE_ptr, int& num_points,
                         GMem* gMem, double tolerance=0) const;
  
  CubitStatus number_VERTICES( BODY* BODY_ptr, int& num_vertices ) const;
  CubitStatus number_EDGES( BODY* BODY_ptr, int& num_edges ) const;
  CubitStatus number_FACES( BODY* BODY_ptr, int& num_faces ) const;
  CubitStatus number_VOLUMES( BODY* BODY_ptr, int& num_volumes ) const;
  CubitStatus number_ENTITIES( BODY* BODY_ptr, int &num_volumes,
                               int& num_faces, int& num_edges,
                               int& num_vertices ) const;
    //- get number of acis entities in given body
  
  CubitStatus bodysms ( ENTITY *entity, DLIList<BodySM*  > &bodies   ) const;
  CubitStatus lumps   ( ENTITY *entity, DLIList<Lump*    > &lumps    ) const;
  CubitStatus shellsms( ENTITY *entity, DLIList<ShellSM* > &shellsms ) const;
  CubitStatus surfaces( ENTITY *entity, DLIList<Surface* > &surfaces ) const;
  CubitStatus loopsms ( ENTITY *entity, DLIList<LoopSM*  > &loopsms  ) const;
  CubitStatus curves  ( ENTITY *entity, DLIList<Curve*   > &curves   ) const;
  CubitStatus coedgesms(ENTITY *entity, DLIList<CoEdgeSM*> &coedgesms) const;
  CubitStatus points  ( ENTITY *entity, DLIList<Point*   > &points   ) const;
    //- topology traversal of TB's; returns list of TopologyBridge
    //- entities

  CubitBoolean is_Related(ENTITY *entity1, ENTITY *entity2) const;
    //- returns true if two ACIS ENTITYs are the same or related

  CubitStatus get_ENTITY_from_ENTITY(ENTITY *from_entity,
                                     const ENTITY *from_type,
                                     ENTITY_LIST &entities) const;
  CubitStatus get_ENTITY_from_ENTITY(const int ident,
                                     ENTITY *entity,
                                     ENTITY_LIST &entities) const;

  CubitStatus get_child_ENTITYs(ENTITY *from_entity,
                                ENTITY_LIST &entities,
                                CubitBoolean get_all = CUBIT_FALSE) const;

#ifdef BOYD14
  CubitStatus get_parent_ENTITYs(ENTITY *from_entity,
                                ENTITY_LIST &entities,
                                CubitBoolean get_all = CUBIT_FALSE) const;
#endif

    //- topology traversal of TB's; implemented based on native traversal
    //- functions in the modeler; returns lists of acis entities

  CubitStatus get_FACEs( ENTITY* ENTITY_ptr, DLIList<FACE*>& FACE_list ) const;
  CubitStatus get_LOOPs( ENTITY* ENTITY_ptr, DLIList<LOOP*>& LOOP_list ) const;
  CubitStatus get_EDGEs( ENTITY* ENTITY_ptr, DLIList<EDGE*>& EDGE_list ) const;
  CubitStatus get_VERTICEs( ENTITY* ENTITY_ptr, DLIList<VERTEX*>& VERTEX_list ) const;
  //- Traversal functions based on native traversal functions in modeler,
  //- return DLLists of ACIS entities.

  CubitStatus get_intersections( EDGE *EDGE_ptr1, EDGE *EDGE_ptr2, 
                                 DLIList<CubitVector*> &int_list,
                                 bool bounded = false ) const;
  //- Gets all intersections of EDGE and EDGE.  Allocates CubitVectors in
  //- int_list - calling code must free them when done.  By default 
  //- intersections are calculated on infinite curves, but you can limit
  //- them to be bounded.

  CubitStatus get_intersections( EDGE *EDGE_ptr, FACE *FACE_ptr, 
                                 DLIList<CubitVector*> &int_list,
                                 bool bounded = false ) const;
  //- Gets all intersections of EDGE and FACE.  Allocates CubitVectors in
  //- int_list - calling code must free them when done.

  CubitStatus entity_extrema( ENTITY_LIST &ENTITY_list, 
                              const CubitVector *dir1, 
                              const CubitVector *dir2,
                              const CubitVector *dir3, 
                              CubitVector &extrema,
                              ENTITY *&extrema_ENTITY_ptr );
  //- Gets the extrema position along the first given direction. If there 
  //- is more than one extrema position, the other directions will be used 
  //- to determine a unique position.  Directions 2 and 3 can be NULL.
  //- Entities supported include groups, bodies, volumes, surfaces, curves,
  //- and vertices.  The entity the extrema is found on is also returned.

  CubitStatus entity_entity_distance( ENTITY *ENTITY_ptr1,
                                      ENTITY *ENTITY_ptr2,
                                      CubitVector &pos1, 
                                      CubitVector &pos2,
                                      double &distance );
  //- Gets the minimum distance between two entities and the closest positions 
  //- on those entities. Supports vertices, curves, surfaces, volumes and bodies.

  CubitStatus get_EDGE_normal( EDGE *EDGE_ptr, CubitVector &norm );
  //- Gets a normal vector to a planar EDGE.  If the EDGE is linear or
  //- non-planar return CUBIT_FAILURE.

  CubitBoolean is_curve_app_straight( Curve *ref_edge_ptr );
  //- Determines if an edge is approximately straight.  This is useful if
  //- you have splines that should be straight lines.  As of ACIS 6.2, the
  //- IGES reader is still making all edges that are part of faces splines.

  double volume(BODY *this_body) const;
    //- get the volume of the passed in body

#ifdef BOYD14
  CubitSense align_surfaces( SURFACE *right_surf,
                             SURFACE *test_surf,
                             FACE *FACE_ptr );
    //- returns the sense of the two geometric surfaces.
#endif
    
  CubitStatus remove_refinements( ENTITY_LIST& list );
    //- Clear any refinements set on objects in .SAT file after
    //- import.
  
  static const int MAX_NUM_CURVE_POINTS;
  static const int AQE_SUBMINOR_VERSION;
  static SPAposition *curvePoints;
  static double   *curveParams;
  AcisFacetManager* facetManager;
  
  bool stepInitialized;

  int exportVersion;
};

// ********** BEGIN INLINE FUNCTIONS          **********
// ********** END INLINE FUNCTIONS            **********

// ********** BEGIN FRIEND FUNCTIONS          **********
// ********** END FRIEND FUNCTIONS            **********

// ********** BEGIN EXTERN FUNCTIONS          **********
// ********** END EXTERN FUNCTIONS            **********

// ********** BEGIN HELPER CLASS DECLARATIONS **********
// ********** END HELPER CLASS DECLARATIONS   **********

#endif
