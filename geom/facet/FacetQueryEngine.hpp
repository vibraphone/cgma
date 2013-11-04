//-------------------------------------------------------------------------
// Filename      : FacetQueryEngine.hpp
//
// Purpose       : Facet geometry engine.
//
//                 This class is implemented as a Singleton pattern. Only
//                 one instance is created and it is accessed through the 
//                 {instance()} static member function.
//
// Special Notes :
//
// Creator       : David R. White
//
// Creation Date : 6/29/00
//
//-------------------------------------------------------------------------

#ifndef FACET_GEOMETRY_ENGINE_HPP
#define FACET_GEOMETRY_ENGINE_HPP

// ********** BEGIN STANDARD INCLUDES         **********

#include <typeinfo>
#if !defined(WIN32)
using std::type_info;
#endif

// ********** END STANDARD INCLUDES           **********

// ********** BEGIN CUBIT INCLUDES            **********

#include "CubitFileIOWrapper.hpp"
#include "GeometryQueryEngine.hpp"

// ********** END CUBIT INCLUDES              **********

// ********** BEGIN FORWARD DECLARATIONS
class TopologyEntity;
class TopologyBridge;
class RefEntity;
class Body;
class Shell;
class ShellSM;
class Loop;
class Chain;
class CoEdgeSM;
class LoopSM;
class RefVolume;
class RefFace;
class RefEdge;
class RefVertex;
class TBPoint;
class Curve;
class Surface;
class Lump;
class BodySM;

class GeometryEntity;
class CubitBox;
class CubitString;

class FacetLump;
class FacetShell;
class FacetLoop;
class FacetSurface;
class FacetBody;
class FacetCoEdge;
class FacetCurve;
class FacetPoint;

class CubitFacet;
class CubitQuadFacet;
class CubitFacetEdge;
class CubitPoint;
class FacetEntity;
class CurveFacetEvalTool;
class FacetEvalTool;

// ********** END FORWARD DECLARATIONS        **********

// ********** BEGIN MACRO DEFINITIONS         **********
// ********** END MACRO DEFINITIONS           **********

// ********** BEGIN ENUM DEFINITIONS          **********
typedef enum {
  CUBIT_FACET_FILE,
  AVS_FILE,
  CHOLLA_FILE,
  FROM_FACET_LIST, 
  STL_FILE
} FacetFileFormat;

// ********** END ENUM DEFINITIONS            **********

class FacetQueryEngine : public GeometryQueryEngine
{
public:
// ********** BEGIN FRIEND DECLARATIONS        **********
  friend class FacetSurface;
  
// ********** END FRIEND DECLARATIONS        **********

//HEADER- Constructor and Destructor
  
  static FacetQueryEngine* instance();
    //- Singleton pattern
    //- Controlled access and creation of the sole instance of this class.

  virtual ~FacetQueryEngine();

  static void delete_instance();
  
  const char* modeler_type()
     { return "facet"; }

  int get_major_version();

  int get_minor_version();

  int get_subminor_version();

  CubitString get_engine_version_string();
  
//HEADER- RTTI and safe casting functions.
  
  virtual const type_info& entity_type_info() const ;
    //R- The geometric modeler type
    //- This function returns the type of the geometric modeler.
  
  virtual CubitBoolean is_solid_modeler_type() const 
    {return CUBIT_FALSE;}
    //R CubitBoolean
    //R- This  is not a solid modeling engine.
//HEADER- Functions for importing and exporting solid models.

  CubitBoolean can_delete_bodies(DLIList<Body*>body_list);
  
  virtual TBPoint* make_Point( GeometryType point_type,
                             CubitVector const& point) const ;
  virtual Curve* make_Curve(Curve *) const;
  virtual Curve* make_Curve( TBPoint const* ,
                             TBPoint const* ,
                             RefFace* ,
                             CubitVector * ) const;
  virtual Curve* make_Curve( GeometryType ,
                             TBPoint const* ,
                             TBPoint const* ,
                             DLIList<CubitVector*>& ,
                             RefFace*  ) const;
  virtual Curve* make_Curve( GeometryType ,
                             TBPoint const* ,
                             TBPoint const* ,
                             CubitVector const* ,
                             CubitSense ) const;
  virtual Surface* make_Surface( Surface *,
                                 DLIList<Loop*> &,
                                 CubitBoolean  ) const;


  virtual Surface* make_Surface( GeometryType , 
                                 DLIList<Curve*>& ,
                                 DLIList<Loop*> &,
                                 Surface *) const ;
  virtual Lump* make_Lump( GeometryType , 
                           DLIList<Surface*>&  ) const ;
  virtual BodySM* make_BodySM( Surface * ) const;
    virtual BodySM* make_BodySM( DLIList<Lump*>&  ) const ;

  Body* copy_body( Body *body_ptr );

  virtual CubitStatus get_graphics( BodySM *bodysm, GMem* g_mem,
                         std::vector<Surface*> &surface_to_facets_vector,
                         std::vector<TopologyBridge*> &vertex_edge_to_point_vector,
                         std::vector<std::pair<TopologyBridge*, std::pair<int,int> > > &facet_edges_on_curves,
                         unsigned short normal_tolerance, 
                         double distance_tolerance, double max_edge_length ) const;

  virtual CubitStatus get_graphics( Surface* surface_ptr,
                                          GMem* gMem,
                                          unsigned short normal_tolerance,
                                          double distance_tolerance,
                                          double max_edge_length ) const;

  virtual CubitStatus get_graphics( Curve* curve_ptr,
                                    GMem* gMem,
                                    double angle_tolerance=0,
                                    double distance_tolerance=0,
                                    double max_edge_length=0 ) const;
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
    //I- The tolerance deviation used when facetting the curve (optional
    //I- and currently IGNORED by this engine).
    //- This function gathers and outputs ACIS edge information for
    //- hoops involved in facetting a RefEdge.  If all goes well,
    //- CUBIT_SUCCESS is returned.  Otherwise, CUBIT_FAILURE is
    //- returned.
 
  virtual CubitStatus get_isoparametric_points(Surface* ,
                                               int &, int &,
                                               GMem*&) const;
  
  virtual CubitStatus get_u_isoparametric_points(Surface* ref_face_ptr,
                                                 double v, int &n,
                                                 GMem *&gMem) const;
  
  virtual CubitStatus get_v_isoparametric_points(Surface* ref_face_ptr,
                                                 double u, int &n,
                                                 GMem *&gMem) const;
  
  virtual CubitStatus transform_vec_position( 
    CubitVector const& ,
    BodySM *,
    CubitVector & ) const;
  
  virtual CubitStatus get_intersections( Curve*, CubitVector& point1,
                                         CubitVector&,
                                         DLIList<CubitVector*>& ,
                                         CubitBoolean,
                                         CubitBoolean );

  virtual CubitStatus get_intersections( Curve* , Curve* ,
                                         DLIList<CubitVector*>& ,
                                         CubitBoolean,
                                         CubitBoolean );
  virtual CubitStatus get_intersections( Curve* ref_edge1, Curve* ref_edge2,
                                         DLIList<CubitVector*>& intersection_list,
                                         double offset, 
                                         CubitBoolean ext_first = CUBIT_TRUE );
  virtual CubitStatus get_intersections( Curve* ref_edge, Surface* ref_face,
                                        DLIList<CubitVector*>& intersection_list,
                                        CubitBoolean bounded = CUBIT_FALSE );

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
                                              CubitVector &pos1, CubitVector &pos2,
                                              double &distance );
  //- Gets the minimum distance between two entities and the closest positions 
  //- on those entities. Supports vertices, curves, surfaces, volumes and bodies.

  virtual CubitStatus export_solid_model( DLIList<TopologyBridge*>& bridge_list,
                                          const char* file_name,
                                          Model_File_Type file_type,
                                          const CubitString &cubit_version,
                                          ModelExportOptions &export_options );

  virtual CubitStatus export_solid_model( DLIList<TopologyBridge*>& bridge_list,
                                          char*& p_buffer,
                                          int& n_buffer_size,
                                          bool b_export_buffer)
{return CUBIT_FAILURE;}

  virtual CubitStatus save_temp_geom_file( DLIList<TopologyBridge*>& ref_entity_list,
                                          const char *file_name,
                                          const CubitString &cubit_version,
                                          CubitString &created_file,
                                          CubitString &created_file_type ); 

 virtual CubitStatus import_temp_geom_file(FILE* file_ptr,
                                      const char* file_name,
                                      Model_File_Type file_type,
                                      DLIList<TopologyBridge*> &bridge_list );

 virtual CubitStatus import_solid_model(const char* file_name,
                                         Model_File_Type file_type,
                                         DLIList<TopologyBridge*>& imported_entities,
                                         ModelImportOptions &import_options );

 virtual CubitStatus import_solid_model(DLIList<TopologyBridge*> &imported_entities,
                                        const char* pBuffer,
                                        const int n_buffer_size)
{return CUBIT_FAILURE;}

private:
  CubitStatus import_solid_model(FILE *file_ptr,
                                 DLIList<TopologyBridge*> &imported_entities );

public:
  virtual void delete_solid_model_entities(DLIList<BodySM*>& body_list) const;
    //- Deletes all solid model entities associated with the Bodies in 
    //- the input list. 
      
  virtual CubitStatus delete_solid_model_entities( BodySM* body_ptr ) const;
  virtual CubitStatus delete_solid_model_entities(Surface* surf_ptr ) const;
  virtual CubitStatus delete_solid_model_entities( Curve* curve_ptr ) const;
  virtual CubitStatus delete_solid_model_entities( TBPoint* point_ptr ) const;

  virtual CubitStatus fire_ray( CubitVector &origin,
                                CubitVector &direction,
                                DLIList<TopologyBridge*> &at_entity_list,
                                DLIList<double> &ray_params,
                                int max_hits,
                                double ray_radius,
                                DLIList<TopologyBridge*> *hit_entity_list=0 ) const;
    //- Fire a ray at specified entities, returning the parameters (distances)
    //- along the ray and optionally the entities hit.  Returned lists are
    //- appended to.  Input entities can be any of bodies, volumes, faces,
    //- edges or vertices.  Optionally you can specify the maximum number of
    //- hits to return (default = 0 = unlimited), and the ray radius to use for
    //- intersecting the entities (default = 0.0 = use modeller default).
    //- NOTE: returned entities hit might be "hidden" beneath virtual entities.
    //-       To resolve to visible entities, use "get_visible_ents_for_hits"
    //-       in GeometryQueryTool.

  virtual double get_sme_resabs_tolerance() const; // Gets solid modeler's resolution absolute tolerance
  virtual double set_sme_resabs_tolerance( double new_resabs );

  virtual CubitStatus set_int_option( const char* opt_name, int val );
  virtual CubitStatus set_dbl_option( const char* opt_name, double val );
  virtual CubitStatus set_str_option( const char* opt_name, const char* val );
    //- Set solid modeler options

  static CubitStatus make_facets( int *conn, int nfacets,
                                  DLIList<CubitQuadFacet *> &facet_list );
  static CubitStatus make_facets( int *conn, int nfacets,
                                  DLIList<CubitFacet *> &facet_list );
    //- create facets from a list of points and connectivity

  CubitStatus ensure_is_ascii_stl_file(FILE * fp, CubitBoolean &is_ascii);
  //- returns true in is_ascii if fp points to an ascii stl file

  CubitStatus read_facets_stl_tolerance(   
                                              DLIList<CubitFacet *> &tfacet_list,
                                              DLIList<CubitPoint *> &point_list,
                                              const char * file_name,
                                              int &npoints, 
                                              int &ntri,
                                              long& seek_address,
                                              double tolerance);
  //- read facets from an stl file and combine vertices within tolerance distance

    CubitStatus read_facets_stl(   
                                              DLIList<CubitFacet *> &tfacet_list,
                                              DLIList<CubitPoint *> &point_list,
                                              const char * file_name,
                                              int &npoints, 
                                              int &ntri,
                                              long& seek_address);
  //- read facets from an stl file
  
  CubitStatus import_facets( const char *file_name, 
                             CubitBoolean use_feature_angle, 
                             double feature_angle,
                             double tolerance,
                             int interp_order,
                             CubitBoolean smooth_non_manifold,
                             CubitBoolean split_surfaces,
                             CubitBoolean stitch,
                             CubitBoolean improve,
                             DLIList <CubitQuadFacet *>&quad_facet_list,
                             DLIList <CubitFacet *> &tri_facet_list,
                             DLIList<Surface *> &surface_list,
                             FacetFileFormat file_format = CUBIT_FACET_FILE );
    //- import facets from a file and create a geometry model

  static CubitStatus read_facets( const char * file_name,
                                  int *&conn,
                                  int &npoints, 
                                  int &nquad, int &ntri, 
                                  FacetFileFormat file_format = CUBIT_FACET_FILE );
  static CubitStatus read_cholla_file( const char *file_name, 
                                       double &feature_angle,
                                       DLIList<CubitPoint *> &point_list, 
                                       DLIList<CubitFacet *> &facet_list);
    //- read points and facets from a file

  CubitBoolean is_close(CubitVector &this_point,
                        DLIList<CubitFacet *>&facet_list,
                        CubitFacet *&lastFacet,
                        double tol);
    //- determine if one of the facets in the list is within a
    //- certain distance of the point.

  static CubitStatus export_facets(DLIList<CubitFacet*> &facet_list,
                                  char *filename);
    //-  export a list of facets to a facet file for debugging purposes

  CubitStatus gather_all_facet_entities( DLIList<FacetBody*> &facet_bodies,
                                         DLIList<FacetLump*> &facet_lumps,
                                         DLIList<FacetShell*> &facet_shells,
                                         DLIList<FacetSurface*> &facet_surfaces,
                                         DLIList<FacetLoop*> &facet_loops,
                                         DLIList<FacetCoEdge*> &facet_coedges,
                                         DLIList<FacetCurve*> &facet_curves,
                                         DLIList<FacetPoint*> &facet_points );

  CubitStatus save_facets( FILE *fp, DLIList<FacetSurface*> facet_surfaces,
                                     DLIList<FacetCurve*>   facet_curves, 
                                     DLIList<FacetPoint*>   facet_points ); 

  CubitStatus save_eval_tools( FILE *fp, DLIList<FacetSurface*> facet_surfaces,
                                         DLIList<FacetCurve*> facet_curves );

  CubitStatus dump_facets( FILE *fp,
                           DLIList<CubitFacet *> &facet_list,
                           DLIList<CubitFacetEdge *> &edge_list, 
                           DLIList<CubitPoint *> &point_list );
  CubitStatus gather_facets( DLIList<FacetSurface *> facet_surfaces,
                             DLIList<FacetCurve *> facet_curves,
                             DLIList<FacetPoint *> facet_points,
                             DLIList<CubitFacet *> &facet_list,
                             DLIList<CubitFacetEdge *> &edge_list, 
                             DLIList<CubitPoint *> &point_list );
    //- functions for saving the facet geometry representation to a cubit file
  
  CubitStatus restore_eval_tools( FILE *fp,
                                  unsigned int endian, 
                                  int num_facets,
                                  int num_edges,
                                  int num_points,
                                  CubitFacet **facets,     
                                  CubitFacetEdge **edges,     
                                  CubitPoint **points,
                                  int &num_cfet,
                                  int &num_fet,
                                  CurveFacetEvalTool **&cfeval_tools,
                                  FacetEvalTool **&feval_tools );
    //- restore facets from CUB file

CubitStatus create_super_facet_bounding_box(
                                DLIList<BodySM*>& body_list,
                                CubitBox& super_box );
CubitStatus create_facet_bounding_box(
                                BodySM* bodySM,
                                CubitBox& bbox );

  CubitStatus restore_transform( BodySM* body );

  CubitStatus translate( BodySM* body, const CubitVector& offset );
  CubitStatus rotate   ( BodySM* body, const CubitVector& axis, double angle );
  CubitStatus scale    ( BodySM* body, double factor );
  CubitStatus scale    ( BodySM* body, const CubitVector& factors );
  CubitStatus reflect  ( BodySM* body, const CubitVector& axis );

  CubitStatus translate( GeometryEntity* ent, const CubitVector& offset );
  CubitStatus rotate   ( GeometryEntity* ent, const CubitVector& axis, double degrees );
  CubitStatus scale    ( GeometryEntity* ent, double factor );
  CubitStatus scale    ( GeometryEntity* ent, const CubitVector& factors );
  CubitStatus reflect  ( GeometryEntity* ent, const CubitVector& axis );

  CubitStatus get_connected_patch( DLIList<FacetSurface*>& remaining_surfs,
                                   DLIList<FacetSurface*>& output_patch );
  virtual CubitBoolean bodies_overlap (BodySM *body_ptr_1,
                                       BodySM *body_ptr_2 ) const;
  virtual CubitBoolean volumes_overlap (Lump *lump_ptr_1,
                                        Lump *lump_ptr_2 ) const;
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
  
  FacetQueryEngine();
  
private:

  CubitStatus write_topology( FILE *file_ptr, 
                              DLIList<FacetBody*> &facet_bodies,
                              DLIList<FacetLump*> &facet_lumps,
                              DLIList<FacetShell*> &facet_shells,
                              DLIList<FacetSurface*> &facet_surfaces,
                              DLIList<FacetLoop*> &facet_loops,
                              DLIList<FacetCoEdge*> &facet_coedges,
                              DLIList<FacetCurve*> &facet_curves,
                              DLIList<FacetPoint*> &facet_points );

  CubitStatus restore_topology( FILE *file_ptr, 
                                unsigned int endian, 
                                int num_points,
                                CubitPoint **points_array,
                                int num_cfet,
                                CurveFacetEvalTool **cfev_array,
                                int num_fet,
                                FacetEvalTool **fev_array,
                                DLIList<TopologyBridge*> &imported_entities );

  CubitStatus restore_facets( FILE *file_ptr,
                              unsigned int endian, 
                              int &num_facets,
                              int &num_edges,
                              int &num_points,
                              CubitPoint**&points_array,
                              int &num_cfet,
                              int &num_fet,
                              CurveFacetEvalTool **&cfet_array,
                              FacetEvalTool **&fet_array);

  CubitStatus read_facets( FILE *fp, 
                           unsigned int endian, 
                           int &num_facets, int &num_edges, int &num_points, 
                           CubitFacet **&facets, CubitFacetEdge **&edges,     
                           CubitPoint **&points );

  static FacetQueryEngine* instance_;
    //- static pointer to unique instance of this class

  static CubitStatus init_hash_points( int num_points );
  static CubitStatus add_hash_point( CubitPoint *point_ptr );
  static CubitPoint *get_hash_point( int id );
  static void delete_hash_points( );
  static int get_hash_key( int id );
  static CubitStatus get_all_hash_points(DLIList<CubitPoint *> &point_list);
  static int hashPointSize;
  static DLIList<CubitPoint *> *hashPointArray;
    //- hash table functions used for reading the facet file

  static const int FQE_MAJOR_VERSION;
  static const int FQE_MINOR_VERSION;
  static const int FQE_SUBMINOR_VERSION;
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
