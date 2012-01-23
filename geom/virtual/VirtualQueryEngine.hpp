//-------------------------------------------------------------------------
// Filename      : VirtualQueryEngine.hpp
//
// Purpose       : Interface for creation/destruction of CompositeEntities.
//
// Special Notes : This is not a complete modeling engine.  This class is
//                 a container of functions for interaction with 
//                 CompositeEntities
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 7/31/97
//-------------------------------------------------------------------------

#ifndef VIRTUAL_GEOMETRY_ENGINE_HPP
#define VIRTUAL_GEOMETRY_ENGINE_HPP

// ********** BEGIN STANDARD INCLUDES         **********
// ********** END STANDARD INCLUDES           **********

// ********** BEGIN MOTIF INCLUDES            **********
// ********** END MOTIF INCLUDES              **********

// ********** BEGIN OPEN INVENTOR INCLUDES    **********
// ********** END OPEN INVENTOR INCLUDES      **********

// ********** BEGIN ACIS INCLUDES             **********
// ********** END ACIS INCLUDES               **********

// ********** BEGIN CUBIT INCLUDES            **********

#include "GeometryQueryEngine.hpp"
#include "CubitDefines.h"
#include "GeometryDefines.h"
#include "CubitString.hpp"
#include "CubitVector.hpp"

// ********** END CUBIT INCLUDES              **********

// ********** BEGIN FORWARD DECLARATIONS      **********
class BasicTopologyEntity;
class Body;
class RefVolume;
class RefEdge;
class RefFace;
class RefVertex;
class Point;
class Curve;
class GMem;
template <class X> class DLIList;
class RefEntity;
class Body;
class CubitVector;
class TopologyEntity;
class CoEdge;

class CompositeCurve;
class CompositeSurface;
class PartitionCurve;
class PartitionSurface;

// ********** END FORWARD DECLARATIONS        **********

// ********** BEGIN MACRO DEFINITIONS         **********
// ********** END MACRO DEFINITIONS           **********

// ********** BEGIN ENUM DEFINITIONS          **********

class VirtualQueryEngine : public GeometryQueryEngine
{
  
public:
  
  virtual int curve_is_on_ignored_surface(Curve *curve,
                    Surface *surf);
  void remove_virtual_geometry( RefEntity* entity_ptr,
                                CubitBoolean all_children );
  void remove_virtual_geometry( Body* body, bool all_children );
  void remove_virtual_geometry( RefVolume* vol, bool all_children );
  void remove_virtual_geometry( RefFace* face );
    //- Remove all possible virtual geometry from the passed
    //- entity and its children.
    
  
  static const double param_epsilon_fraction;
    //- This value is used to determine epsilon for
    //- compensating for rounding error in parameter
    //- conversion between virtaul geometry and real
    //- geometry.
  
#ifdef BOYD15
  static CubitBoolean test_and_round( double min, double max,
                                      double& value, double epsilon );
    //R CubitBoolean
    //R- CUBIT_TRUE/CUBIT_FALSE
    //I min/max
    //I- The range for which to test value
    //I value
    //I- The value to test
    //O value
    //O- The possible adjustment made to value
    //I epsilon
    //I- The allowable distance outside the range
    //- This method is used by various VirtualEntities to
    //- compensate for rounding errors in conversions.  The
    //- method returns CUBIT_TRUE if the passed value is within
    //- the inclusive range specified by min and max.  If the
    //- passed value is outside of the range, but within epsilon
    //- of min or max, the value will be set to min or max, and
    //- CUBIT_TRUE will be returned.  Otherwise CUBIT_FALSE will
    //- be returned.
#endif
  
  inline static VirtualQueryEngine* instance();
    //R CompositeModelingeEngine
    //R- A pointer to the single instance of this class
    //- Access to VirtualQueryEngine must be through this function.
    //- This forces single instance of the class. (constructor is protected)
  
  static void delete_instance();

  void register_attributes();
   //-  VGE registers attributes required for virtual geometry
   // -- only call this once per session

  ~VirtualQueryEngine( );
    //- Destructor
  
  virtual const type_info& entity_type_info( ) const
    { return typeid(VirtualQueryEngine); }
    //R- The type of this entity
    //- This function returns the type of this class.

const char * modeler_type()
   { return "virtual"; }

  
  virtual int get_major_version();
  virtual int get_minor_version();
  virtual int get_subminor_version();
  virtual CubitString get_engine_version_string();

  virtual bool is_intermediate_engine() {return TRUE;}
  
  virtual CubitStatus get_graphics( Surface* surface_ptr,
                                           GMem* gMem,
                                           unsigned short normal_tolerance,
                                           double distance_tolerance,
                                           double max_edge_length ) const;
    //R CubitStatus
    //R- CUBIT_SUCCESS/CUBIT_FAILURE
    //I ref_face_ptr
    //I- The RefFAce for which hoops facetting information will be 
    //I- gathered.
    //O- The number of polygons (triangles) needed for facetting.
    //O number_points
    //O- The number of points needed for facetting
    //O number_Facets
    //O- The number of facets needed for facetting.
    //O gMem
    //O- The sturage place for facets (and points).
    //= This function gathersw and outputs ACIS facet (and point)
    //- information for hoops involved in facetting an RefFace.  If
    //- all goes well, CUBIT_Success is retuned.
  
  virtual CubitStatus get_graphics( Curve* curve_ptr,
                                    GMem* gMem = NULL,
                                    double tolerance = 0.0 ) const;
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
    //I- The tolerance deviation used when facetting the curve (optional
    //I- and currently IGNORED by this engine).
    //- This function gathers and outputs ACIS edge information for
    //- hoops involved in facetting a RefEdge.  If all goes well,
    //- CUBIT_SUCCESS is returned. Otherwise, CUBIT_FAILURE is retuned.
    

  virtual CubitStatus get_isoparametric_points(Surface* ref_face_ptr,
                                               int &nu, int &nv,
                                               GMem *&gMem) const;
  
  virtual CubitStatus get_u_isoparametric_points(Surface* ref_face_ptr,
                                                 double v, int& n,
                                                 GMem *&gMem) const;
  
  virtual CubitStatus get_v_isoparametric_points(Surface* ref_face_ptr,
                                                 double u, int&n,
                                                 GMem *&gMem) const;
  
  virtual CubitStatus transform_vec_position( 
    CubitVector const& position_vector,
    BodySM *OSME_ptr,
    CubitVector &transformed_vector ) const;


  CubitStatus translate( BodySM* body, const CubitVector& offset );
  CubitStatus rotate   ( BodySM* body, const CubitVector& axis, double angle );
  CubitStatus scale    ( BodySM* body, double factor );
  CubitStatus scale    ( BodySM* body, const CubitVector& factors );
  CubitStatus reflect  ( BodySM* body, const CubitVector& axis );
  CubitStatus restore_transform( BodySM* body );

  CubitStatus translate( GeometryEntity* ent, const CubitVector& offset );
  CubitStatus rotate   ( GeometryEntity* ent, const CubitVector& axis, double degrees );
  CubitStatus scale    ( GeometryEntity* ent, double factor );
  CubitStatus scale    ( GeometryEntity* ent, const CubitVector& factors );
  CubitStatus reflect  ( GeometryEntity* ent, const CubitVector& axis );

  static void get_VEs( RefVolume* volume_ptr,
                DLIList<TopologyBridge*>& ve_list,
                CubitBoolean visible = CUBIT_FALSE,
                const CubitBoolean children_too = CUBIT_TRUE);
  static void get_VEs( RefFace* face_ptr, 
                DLIList<TopologyBridge*>& ve_list,
                CubitBoolean visible = CUBIT_FALSE,
                const CubitBoolean children_too = CUBIT_TRUE );
  static void get_VEs( RefEdge* edge_ptr,
                DLIList<TopologyBridge*>& ve_list,
                CubitBoolean visible = CUBIT_FALSE,
                const CubitBoolean children_too = CUBIT_TRUE );
  static void get_VEs( RefVertex* vertex_ptr,
                DLIList<TopologyBridge*>& ve_list,
                CubitBoolean visible = CUBIT_FALSE,
                const CubitBoolean children_too = CUBIT_TRUE );
  static void get_VEs( TopologyEntity* te_ptr,
                DLIList<TopologyBridge*>& ve_list,
                CubitBoolean visible = CUBIT_FALSE,
                const CubitBoolean children_too = CUBIT_TRUE );
    //R void
    //I te_ptr, volume_ptr, face_ptr, edge_ptr
    //I- Entity to query.
    //O ve_list
    //O- result set.
    //I visible
    //I- Return only entities visible in the model.
    //- Get any virtual geometry of the passed entity or its children.
  
  virtual TopologyBridge* get_visible_entity_at_point(TopologyBridge* hidden_tb, CubitVector* point);
   //R TopologyBridge*
   //I hidden_tb, point
   //I- Returns the lowest level visible entity that contains 'point.' (used for fire_ray())

  CubitBoolean virtuals_created() const;
    //- returns CUBIT_TRUE if there are virtual entities in the model
  
#ifdef BOYD15
  static CubitBoolean is_real( TopologyEntity* entity );
#endif
  static CubitBoolean has_virtual( TopologyEntity* entity);
  
  static CubitBoolean is_virtual(TopologyEntity *entity,
                                 const CubitBoolean children_too = CUBIT_FALSE);
  static CubitBoolean is_virtual(DLIList<Body*> &entity_list,
                                 const CubitBoolean children_too = CUBIT_FALSE);
    //- returns CUBIT_TRUE if the entities passed in are virtual or
    //- contain any virtual entities
  
  static CubitBoolean is_virtual(DLIList<RefEntity*> &entity_list,
                          const CubitBoolean children_too = CUBIT_FALSE);
    //- returns CUBIT_TRUE if the entities passed in are virtual or
    //- contain any virtual entities

  CubitBoolean is_partition( RefEntity *ref_entity );

  CubitStatus get_sister_partitions( RefEntity *ref_entity,
                                     DLIList<RefEntity*> &sisters);
  
  CubitStatus sort_edges( DLIList<RefEdge*>& edge_list ) const;
    //R CubitStatus
    //R- CUBIT_SUCCESS/CUBIT_FAILURE
    //I edge_list
    //I- Input: The list of RefEdges to sort.
    //I- Output: The list of sorted RefEdges.
    //O composite_sense_wrt_first
    //O- The sense of the CompositeCurve relative to the first RefEdge's
    //O- Curve geometry object.  NOT w.r.t. the RefEdge.
    //O- in the SORTED list of RefEdges RETURNED.
    //- This function sorts a list of RefEdges using topological
    //- information.  The RefEdges are sorted such that they form a
    //- consecutive chain of edges.  CUBIT_SUCCESS is returned if the
    //- RefEdges are C0 continous, and CUBIT_FAILURE if they are not.
    //- The C0 continuity is checked topologically, which assumes the
    //- model is manifold.
		
  static CubitStatus associate_curves( DLIList<Curve*>& target_set,  
                                       RefEdge* edge_ptr,
                                       DLIList<Curve*>& reset_set );
  //- Find the subset of target_set that correspond to the curves contained
  //- in the passed edge_ptr OR any hosts of the passed edge_ptr.

  virtual CubitStatus get_underlying_curves(Curve *curve_ptr, 
                             DLIList<TopologyBridge*>& curve_list);
  virtual CubitStatus get_underlying_surfaces(Surface * surf_ptr,
                                 DLIList<TopologyBridge*>& surf_list)  ;
  virtual CubitStatus get_underlying_bridges(TopologyBridge* bridge_ptr,
                                        DLIList<TopologyBridge*>& bridge_list);

  virtual CubitBoolean bodies_overlap (BodySM *body_ptr_1, BodySM *body_ptr_2 ) const;
  virtual CubitBoolean volumes_overlap (Lump *lump1, Lump *lump2 ) const;

protected:
  
  VirtualQueryEngine(); 
    //- Constructor.
  
  static VirtualQueryEngine* instance_;

//**** Display related methods
  
  CubitStatus get_composite_curve_facetting( CompositeCurve* ccurve_ptr,
                                             GMem* gMem ) const;
  CubitStatus get_partition_curve_facetting( PartitionCurve* pcurve_ptr,
                                             GMem* gMem ) const;
  
  CubitStatus get_composite_surface_facetting( CompositeSurface* surf_ptr,
                                               GMem* gMem,
                                               unsigned short normal_tol,
                                               double absolute_tol,
                                               double longest_edge ) const;

  CubitStatus get_partition_surface_facetting( PartitionSurface* surf_ptr,
                                               GMem* gMem,
                                               unsigned short normal_tol,
                                               double absolute_tol,
                                               double longest_edge ) const;
  
#ifdef BOYD15
  CubitVector gmem_to_vector( GMem* gmem_ptr, int point_index ) const;
  void vector_to_gmem( CubitVector vect, GMem* gmem_ptr, 
                       int point_index ) const;
    //R CubitVector
    //R- The corresponding CubitVector
    //I vect
    //I- The corresponding CubitVector
    //I gmem_ptr
    //I- A pointer to a GMem object
    //I point_index
    //I- The index of a point in the GPoints list of the GMem object.
    //- These methods convert a GMem object's GPoint to a 
    //- CubitVector, and vise versa.
  
//	CubitStatus get_virtual_curve_facetting( ParasiteCurve* vcurve_ptr,
//						 int& num_steps,
//						 GMem* gMem ) const;
    //- Type-specific methods for generating display information.	
#endif

  virtual CubitStatus export_solid_model( DLIList<TopologyBridge*>& bridge_list,
                                          const char* file_name,
                                          const char* file_type,
                                          const CubitString &cubit_version,
                                          const char* logfile_name = NULL );

  virtual CubitStatus export_solid_model( DLIList<TopologyBridge*>& bridge_list,
					  char*& p_buffer,
					  int& n_buffer_size,
					  bool b_export_buffer);
  
  virtual CubitStatus save_temp_geom_file(DLIList<TopologyBridge*>& bridge_list,
                                          const char *file_name,
                                          const CubitString &cubit_version,
                                          CubitString &created_file,
                                          CubitString &created_file_type);

  virtual CubitStatus import_temp_geom_file(FILE* file_ptr, 
                                 const char* file_name,
                                 const char* file_type,
                                 DLIList<TopologyBridge*> &bridge_list);

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
                             CubitBoolean free_surfaces = CUBIT_TRUE );

  virtual CubitStatus import_solid_model(DLIList<TopologyBridge*> &imported_entities,
					 const char* pBuffer,
					 const int n_buffer_size);

  virtual void delete_solid_model_entities(DLIList<BodySM*>&) const;
  virtual CubitStatus delete_solid_model_entities( BodySM* body_ptr ) const;
  virtual CubitStatus delete_solid_model_entities(Surface* surf_ptr ) const;
  virtual CubitStatus delete_solid_model_entities( Curve* curve_ptr ) const;
  virtual CubitStatus delete_solid_model_entities( Point* point_ptr ) const;

  virtual CubitStatus fire_ray( const CubitVector &origin,
                                const CubitVector &direction,
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

  virtual double get_sme_resabs_tolerance() const;
  virtual double set_sme_resabs_tolerance( double new_resabs );

  virtual CubitStatus set_int_option( const char* opt_name, int val );
  virtual CubitStatus set_dbl_option( const char* opt_name, double val );
  virtual CubitStatus set_str_option( const char* opt_name, const char* val );

private:
  
  void default_error_message( const char callers_name[] ) const;
    //- Error message for inhereted functions which are not implemented.
  
  static const int VGE_MAJOR_VERSION;
  static const int VGE_MINOR_VERSION;
  static const int VGE_SUBMINOR_VERSION;

public: 
  virtual CubitStatus get_intersections(Curve* curve, CubitVector &point1,
                                         CubitVector &point2,
                                         DLIList<CubitVector*>& intersection_list,
                                         CubitBoolean bounded,
                                         CubitBoolean closest );

  virtual CubitStatus get_intersections( Curve* ref_edge1, Curve* ref_edge2,
                                         DLIList<CubitVector*>& intersection_list,
                                         bool bounded = CUBIT_FALSE,
                                         bool closest = CUBIT_FALSE );
  //- Finds the intersections of the two curves.  If the bounded flag is
  //- true, it finds only those intersections within the parameter range
  //- of both curves; otherwise it uses the extensions of the curves.  The
  //- closest option is currently valid only if both curves are straight,
  //- in which case it will return the 2 closest intersection locations,
  //- if the straight lines don't actually intersect. The function allocates 
  //- the CubitVectors in the returned list, so be sure to free them.
  //- NOT DEFINED YET IN THIS GEOMETRY ENGINE.

  virtual CubitStatus get_intersections( Curve* ref_edge, Surface* ref_face,
                                         DLIList<CubitVector*>& intersection_list,
                                         bool bounded = CUBIT_FALSE );
  //- Finds the intersections of the curve and surface.  The curve is extended
  //- if the bounded flag is not false.  The function allocates the CubitVectors 
  //- in the returned list, so be sure to free them.
  //- NOT DEFINED YET IN THIS GEOMETRY ENGINE.

  virtual CubitStatus entity_extrema( DLIList<GeometryEntity*> &ref_entity_list, 
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
  //- NOT DEFINED YET IN THIS GEOMETRY ENGINE.

  virtual CubitStatus entity_entity_distance( GeometryEntity *ref_entity_ptr1,
                                              GeometryEntity *ref_entity_ptr2,
                                              CubitVector &pos1, CubitVector &pos2,
                                              double &distance );
  //- Gets the minimum distance between two entities and the closest positions 
  //- on those entities. Supports vertices, curves, surfaces, volumes and bodies.
  //- NOT DEFINED YET IN THIS GEOMETRY ENGINE.

};


// ********** BEGIN INLINE FUNCTIONS          **********
inline VirtualQueryEngine* VirtualQueryEngine::instance()
{
  if( instance_ == NULL )
  {
    instance_ = new VirtualQueryEngine;
    assert( instance_ != NULL );
  }
  
  return instance_;
}

// ********** END INLINE FUNCTIONS            **********

// ********** BEGIN FRIEND FUNCTIONS          **********
// ********** END FRIEND FUNCTIONS            **********

// ********** BEGIN EXTERN FUNCTIONS          **********
// ********** END EXTERN FUNCTIONS            **********

// ********** BEGIN HELPER CLASS DECLARATIONS **********
// ********** END HELPER CLASS DECLARATIONS   **********

#endif
