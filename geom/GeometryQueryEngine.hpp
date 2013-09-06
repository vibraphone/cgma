//-------------------------------------------------------------------------
// Filename      : GeoemtricQueryEngine.hpp
//
// Purpose       : Define the interface for all solid model query
//                 engines.
//
// Special Notes : This is an abstract base class.
//
// Creator       : Tim Tautges
//
// Creation Date : 2/01
//
// Owner         : Tim Tautges
//-------------------------------------------------------------------------

#ifndef GEOMETRIC_QUERY_ENGINE_HPP
#define GEOMETRIC_QUERY_ENGINE_HPP

#include "CubitObserver.hpp"
#include "CubitDefines.h"
#include "GeometryDefines.h"
#include "CubitCompat.h"
#include "CubitGeomConfigure.h"

class CubitString;
class CubitVector;
template <class X> class DLIList;
class Body;
class GMem;

class TopologyBridge;
class TopologyEntity;
class GeometryEntity;
class BodySM;
class Lump; 
class Surface;
class Curve;
class TBPoint;

enum Model_File_Type {
  MFT_NOT_DEFINED = 0,
  ACIS_TYPE,
  ACIS_SAT_TYPE,
  ACIS_SAB_TYPE,
  ACIS_DEBUG_TYPE,
  IGES_TYPE,
  CATIA_TYPE,
  STEP_TYPE,
  PROE_TYPE,
  GRANITE_TYPE,
  GRANITE_G_TYPE,
  GRANITE_SAT_TYPE,
  GRANITE_PROE_PART_TYPE,
  GRANITE_PROE_ASM_TYPE,
  GRANITE_NEUTRAL_TYPE,
  NCGM_TYPE,
  CATIA_NCGM_TYPE,
  CATPART_TYPE,
  CATPRODUCT_TYPE,
  FACET_TYPE, 
  SOLIDWORKS_TYPE,
  OCC_TYPE
};


extern const char*
get_MFT_string(Model_File_Type);

inline const char*
get_MFT_string(Model_File_Type type)
{
  switch (type)
  {
    case MFT_NOT_DEFINED :
      return "Model File Type NOT DEFINED";
    case ACIS_TYPE :
      return "ACIS";
    case ACIS_SAT_TYPE :
      return "ACIS_SAT";
    case ACIS_SAB_TYPE :
      return "ACIS_SAB";
    case ACIS_DEBUG_TYPE :
      return "ACIS_DEBUG";
    case IGES_TYPE :
      return "IGES";
    case CATIA_TYPE :
      return "CATIA";
    case STEP_TYPE :
      return "STEP";
    case PROE_TYPE :
      return "PROE";
    case GRANITE_TYPE :
      return "GRANITE";
    case GRANITE_G_TYPE :
      return "GRANITE_G";
    case GRANITE_SAT_TYPE :
      return "GRANITE_SAT";
    case GRANITE_PROE_PART_TYPE :
      return "GRANITE_PROE_PART";
    case GRANITE_PROE_ASM_TYPE :
      return "GRANITE_PROE_ASM";
    case GRANITE_NEUTRAL_TYPE :
      return "GRANITE_NEUTRAL";
    case NCGM_TYPE :
      return "NCGM";
    case CATIA_NCGM_TYPE :
      return "CATIA_NCGM";
    case CATPART_TYPE :
      return "CATPART";
    case CATPRODUCT_TYPE :
      return "CATPRODUCT";
    case FACET_TYPE :
      return "FACET";
    case SOLIDWORKS_TYPE :
      return "SOLIDWORKS";
    default :
      PRINT_ERROR("Model_File_Type index %i is not handled properly in get_MFT_string()\n", type);
  }
  
  return "";
}  //  get_MFT_string()


static const char* const ACIS_SAT = ".sat";
static const char* const ACIS_SAB = ".sab";
static const char* const GRANITE_G = ".g";
static const char* const GRANITE_PROE_PART = ".prt";
static const char* const GRANITE_PROE_ASM = ".asm";
static const char* const GRANITE_PROE_NEUTRAL = ".neu";
static const char* const STEP_IMPORT_DEFAULT_LOG  = "step_import.log";
static const char* const IGES_IMPORT_DEFAULT_LOG = "iges_import.log";
static const char* const CATIA_IMPORT_DEFAULT_LOG = "catia_import.log";
static const char* const ACIS_SAT_IMPORT_DEFAULT_LOG = "sat_import.log";
static const char* const GRANITE_EXPORT_DEFAULT_LOG = "granite_export.log";
static const char* const STEP_EXPORT_DEFAULT_LOG = "step_export.log";
static const char* const IGES_EXPORT_DEFAULT_LOG = "iges_export.log";
static const char* const STEP_IGES_TRANSLATOR = "step_iges_translator";
static const char* const TRANSLATOR_DIR =  "/translator/";

class CUBIT_GEOM_EXPORT GeometryQueryEngine
{

public:

    //!- virtual destructor
  virtual ~GeometryQueryEngine() {}

     //!R CubitStatus
     //!R- CUBIT_SUCCESS/CUBIT_FAILURE
     //!I ref_entity_list
     //!I- A list of RefEntities to be exported or saved to a file.
     //!I file_name
     //!I- The name of the file to write to.
     //!I file_type
     //!I- An optional type of file.
     //!- Export the current CUBIT geometry (everything in the Model) to a
     //!- solid model format. Valid file types are:
     //!-    "ACIS_SAT"    --  ACIS ASCII (SAT) file format
     //!-    "ACIS_SAB"    --  ACIS BINARY (SAB) file format
     //!-    "ACIS_DEBUG"  --  ACIS ASCII debug format
     //!-    "IGES"        --  IGES file
     //!-    "STEP"        --  STEP file
     //!I logfile_name
     //!I- Optional - name of logfile.
     //!- No logfile gets created for SAB/SAT files, but for IGES and
     //!- STEP file a logfile always gets created.  Default filenames
     //!- are assigned if one is not given (iges_export.log, step_export.log).
     //!-
     //!- The function returns CUBIT_FAILURE if anything goes wrong with
     //!- export - improper file type, inaccessible file, mismatch between
     //!- the underlying representation and file type. It returns
     //!- CUBIT_SUCCESS if everything goes well.
      virtual CubitStatus export_solid_model(
                                   DLIList<TopologyBridge*>& bridge_list,
                                   const char* file_name,
                                   Model_File_Type file_type,
                                   const CubitString &cubit_version,
                                   ModelExportOptions &export_options ) = 0;

      virtual CubitStatus export_solid_model(
                                   DLIList<TopologyBridge*>& bridge_list,
                                   char*& p_buffer,
                                   int& n_buffer_size,
                                   bool b_export_buffer) = 0;
    
     //! Saves out a temporary geometry file.  Entities in list must all be 
     //! of same modeling engine.
     virtual CubitStatus save_temp_geom_file(
                                 DLIList<TopologyBridge*> &ref_entity_list,
                                 const char *filename,
                                 const CubitString &cubit_version,
                                 CubitString &created_file,
                                 CubitString &created_file_type ) = 0;

     //! Imports entities into CGM from file that was embedded in a cub file.  
     //! This file is could an ACIS, granite, catia, or another supported
     //! solid modeling engine file.
     virtual CubitStatus import_temp_geom_file(
                                 FILE* file_ptr,
                                 const char* file_name,
                                 Model_File_Type file_type,
                                 DLIList<TopologyBridge*> &bridge_list  ) = 0;

     //!R CubitStatus
     //!R- CUBIT_SUCCESS/CUBIT_FAILURE
     //!I file_ptr
     //!I- A pointer to the file to read (can be NULL for IGES and STEP files).
     //!I file_type
     //!I- Type of file.
     //!- Reads in geometry and creates the necessary Reference entities
     //!- associated with the input geometry.
     //!- Valid file types are:
     //!-    "ACIS_SAT"    --  ACIS ASCII (SAT) file format
     //!-    "ACIS_SAB"    --  ACIS BINARY (SAB) file format
     //!-    "IGES"        --  IGES file
     //!-    "STEP"        --  STEP file
     //!I heal_step - auto-healing of step bodies on import.  This is recommended
     //!              because they always need it.
     //!I import_bodies (etc...)
     //!I- Should bodies be import.
     //!- Function can selectively import solid bodies, free surfaces, free
     //!- curves, or free vertices.  For example, the user may not want
     //!- to import any free entities.
     //!-
     //!- The function returns CUBIT_FAILURE if anything goes wrong with
     //!- import - improper file type, inaccessible file, mismatch between
     //!- the underlying representation and file type. It returns
     //!- CUBIT_SUCCESS if everything goes well.
     virtual CubitStatus import_solid_model(
                                const char* file_name,
                                Model_File_Type file_type,
                                DLIList<TopologyBridge*>& imported_entities,
                                ModelImportOptions &import_options ) = 0;

     virtual CubitStatus import_solid_model(DLIList<TopologyBridge*> &imported_entities,
                                            const char* pBuffer,
                                            const int n_buffer_size) = 0;

     //O imported_entities
     //O- List of top-level entities read from file
     //I print_results
     //I- If false, fail silently (don't write error messages to stdout or stderr)

      virtual CubitStatus get_underlying_curves(Curve * curve_ptr,
                                 DLIList<TopologyBridge*>& curve_list)  ;
      // -currently used only by VirtualGeometryEngine to get the underlying
      // - curves from virtual curves.
      virtual CubitStatus get_underlying_surfaces(Surface * surf_ptr,
                                 DLIList<TopologyBridge*>& surf_list)  ;
      virtual CubitStatus get_underlying_bridges(TopologyBridge* bridge_ptr,
                                 DLIList<TopologyBridge*>& bridge_list);

      virtual CubitStatus get_intersections(
                                    Curve* ref_edge1, CubitVector& point1,
                                    CubitVector& point2,
                                    DLIList<CubitVector*>& intersection_list,
                                    bool bounded = false,
                                    bool closest = false ) = 0;

      virtual CubitStatus get_intersections(
                                    Curve* ref_edge1, Curve* ref_edge2,
                                    DLIList<CubitVector*>& intersection_list,
                                    bool bounded = false,
                                    bool closest = false ) = 0;
      //- Finds the intersections of the two curves.  If the bounded flag is
      //- true, it finds only those intersections within the parameter range
      //- of both curves; otherwise it uses the extensions of the curves.  The
      //- closest option is currently valid only if both curves are straight,
      //- in which case it will return the 2 closest intersection locations,
      //- if the straight lines don't actually intersect. The function allocates
      //- the CubitVectors in the returned list, so be sure to free them.

      virtual CubitStatus get_intersections(
                                   Curve* curve,
                                   Surface* surface,
                                   DLIList<CubitVector*>& intersection_list,
                                   bool bounded = false ) = 0;
      //- Finds the intersections of the curve and surface.  The curve is extended
      //- if the bounded flag is not false.  The function allocates the CubitVectors
      //- in the returned list, so be sure to free them.

      virtual CubitStatus entity_extrema(
                                DLIList<GeometryEntity*> &ref_entity_list,
                                const CubitVector *dir1,
                                const CubitVector *dir2,
                                const CubitVector *dir3,
                                CubitVector &extrema,
                                GeometryEntity *&extrema_entity_ptr ) = 0;
      //- Gets the extrema position along the first given direction. If there
      //- is more than one extrema position, the other directions will be used
      //- to determine a unique position.  Directions 2 and 3 can be NULL.
      //- Entities supported include bodies, volumes, surfaces, curves and
      //- vertices.  The entity the extrema is found on is also returned.

      virtual CubitStatus entity_entity_distance(
                                         GeometryEntity *ref_entity_ptr1,
                                         GeometryEntity *ref_entity_ptr2,
                                         CubitVector &pos1,
                                         CubitVector &pos2,
                                         double &distance ) = 0;
      //- Gets the minimum distance between two entities and the closest positions
      //- on those entities. Supports vertices, curves, surfaces, volumes and bodies.

      virtual void delete_solid_model_entities(DLIList<BodySM*>& body_list) const = 0;
      virtual CubitStatus delete_solid_model_entities( BodySM* body_ptr ) const = 0;
      virtual CubitStatus delete_solid_model_entities(Surface* surf_ptr ) const = 0;
      virtual CubitStatus delete_solid_model_entities( Curve* curve_ptr ) const = 0;
      virtual CubitStatus delete_solid_model_entities( TBPoint* point_ptr ) const = 0;
      //- Deletes the solid model entities associcated with the input
      //- free-floating RefEntity. If the input RefEntity is not free-floating,
      //- then its underlying ACIS ENTITYs are not deleted and CUBIT_FAILURE
      //- is returned.
      //- Returns CUBIT_SUCCESS if all went well, otherwise, CUBIT_FAILURE.

  virtual CubitStatus fire_ray( CubitVector &origin,
                                CubitVector &direction,
                                DLIList<TopologyBridge*> &at_entity_list,
                                DLIList<double> &ray_params,
                                int max_hits = 0,
                                double ray_radius = 0.0,
                                DLIList<TopologyBridge*> *hit_entity_list=0 ) const = 0;
    //- Fire a ray at specified entities, returning the parameters (distances)
    //- along the ray and optionally the entities hit.  Returned lists are
    //- appended to.  Input entities can be any of bodies, volumes, faces,
    //- edges or vertices.  Optionally you can specify the maximum number of
    //- hits to return (default = 0 = unlimited), and the ray radius to use for
    //- intersecting the entities (default = 0.0 = use modeller default).

  virtual CubitStatus get_isoparametric_points(Surface* ref_face_ptr,
                                               int &nu, int &nv,
                                               GMem *&gMem) const = 0;

  virtual CubitStatus get_u_isoparametric_points(Surface* ref_face_ptr,
                                                 double v, int& n,
                                                 GMem *&gMem) const = 0;

  virtual CubitStatus get_v_isoparametric_points(Surface* ref_face_ptr,
                                                 double u, int&n,
                                                 GMem *&gMem) const = 0;

  virtual CubitStatus transform_vec_position(
    CubitVector const& position_vector,
    BodySM *OSME_ptr,
    CubitVector &transform_vector ) const = 0;
    //R CubitStatus
    //R- CUBIT_SUCCESS/CUBIT_FAILURE
    //I position_vector on body, body that was transformed
    //O vector showning transform.
    //O-Computes the transform from the transformed body.

  virtual int curve_is_on_ignored_surface(Curve* /* curve */,
                    Surface* /* surf */) { return 0; }

  virtual const char* modeler_type() = 0;
  virtual int get_major_version() = 0;
  virtual int get_minor_version() = 0;
  virtual int get_subminor_version() = 0;
  virtual int get_allint_version();
  virtual CubitString get_engine_version_string() = 0;

  virtual bool is_intermediate_engine() {return FALSE;}

    //- pass a string back identifying this query engine
  virtual CubitStatus set_export_version(const int major,
                                         const int minor);
  virtual CubitStatus set_export_allint_version(const int version);
  virtual CubitStatus list_engine_versions(CubitString &versions);
    //- gets, sets, and lists choices for solid modeling engine version(s)

  virtual double get_sme_resabs_tolerance() const = 0; // Gets solid modeler's resolution absolute tolerance
  virtual double set_sme_resabs_tolerance( double new_resabs ) = 0;

  virtual CubitStatus set_int_option( const char* opt_name, int val ) = 0;
  virtual CubitStatus set_dbl_option( const char* opt_name, double val ) = 0;
  virtual CubitStatus set_str_option( const char* opt_name, const char* val ) = 0;
   //- Set solid modeler options

  virtual CubitStatus get_graphics( BodySM *bodysm, GMem* g_mem,
                                    std::vector<Surface*> &surface_to_facets_vector,
                                    std::vector<TopologyBridge*> &vertex_edge_to_point_vector,
                                    std::vector<std::pair<TopologyBridge*, std::pair<int,int> > > &facet_edges_on_curves,
                                    unsigned short normal_tolerance, 
                                    double distance_tolerance, double max_edge_length ) const;
  
  virtual CubitStatus get_graphics( Surface* surface_ptr,
                                    GMem* gMem,
                                    unsigned short normal_tolerance=15,
                                    double distance_tolerance=0,
                                    double longest_edge = 0 ) const = 0;

  virtual CubitStatus get_graphics( Surface* surface_ptr,
                            GMem *gmem,
                            std::vector<TopologyBridge*> &vertex_edge_to_point_vector,
                            std::vector<std::pair<TopologyBridge*, std::pair<int,int> > > &facet_edges_on_curves,
                            unsigned short normal_tolerance = 15, 
                            double distance_tolerance = 0.0, 
                            double max_edge_length = 0.0 ) const;

    //R CubitStatus
    //R- CUBIT_SUCCESS/CUBIT_FAILURE
    //I ref_face_ptr
    //I- The RefFAce for which hoops facetting information will be
    //I- gathered.
    //O gMem
    //O- The sturage place for facets (and points).
    //= This function gathers and outputs ACIS facet (and point)
    //- information for hoops involved in facetting an RefFace.  If
    //- all goes well, CUBIT_Success is retuned.

  virtual CubitStatus get_graphics( Curve* curve_ptr,
                                    GMem* gMem,
                                    double angle_tolerance=0,
                                    double distance_tolerance=0,
                                    double max_edge_length=0 ) const = 0;
    //R CubitStatus
    //R- CUBIT_SUCCESS/CUBIT_FAILURE
    //I ref_edge_ptr
    //I- The RefEdge for which hoops facetting information will be
    //I- gathered/
    //O numSteps
    //O- The number of edges gathered.
    //O gMem
    //O- The sturage place for edges involved in facetting.
    //I tolerance
    //I- The tolerance deviation used when facetting the curve.
    //- This function tessellates the edge which is typically needed
    //- for graphically displaying the edge. If all goes well,
    //- CUBIT_SUCCESS is returned. Otherwise, CUBIT_FAILURE is returned.


  virtual CubitStatus translate( BodySM* body, const CubitVector& offset ) = 0;
  virtual CubitStatus rotate   ( BodySM* body, const CubitVector& axis, double angle ) = 0;
  virtual CubitStatus scale    ( BodySM* body, double factor ) = 0;
  virtual CubitStatus scale    ( BodySM* body, const CubitVector& factors ) = 0;
  virtual CubitStatus reflect  ( BodySM* body, const CubitVector& axis ) = 0;
  virtual CubitStatus restore_transform( BodySM* body ) = 0;

  virtual CubitStatus translate( GeometryEntity* ent, const CubitVector& offset ) = 0;
  virtual CubitStatus rotate   ( GeometryEntity* ent, const CubitVector& axis, double degrees ) = 0;
  virtual CubitStatus scale    ( GeometryEntity* ent, double factor ) = 0;
  virtual CubitStatus scale    ( GeometryEntity* ent, const CubitVector& factors ) = 0;
  virtual CubitStatus reflect  ( GeometryEntity* ent, const CubitVector& axis ) = 0;

  virtual CubitBoolean bodies_overlap (BodySM *body_ptr_1,
                                       BodySM *body_ptr_2 ) const = 0;

  virtual CubitBoolean volumes_overlap (Lump *lump1, Lump *lump2 ) const = 0; 

  virtual TopologyBridge* get_visible_entity_at_point(TopologyBridge* /* hidden_tb */, CubitVector* /* point */){return NULL;};

  virtual CubitStatus get_visible_entities( TopologyBridge *hidden_tb, DLIList<TopologyBridge*> &real_tbs );

   protected:

   private:

} ;

#endif
