//-------------------------------------------------------------------------
// Filename      : SWQueryEngine.hpp
//
// Purpose       : SolidWorks solid model query engine
//
// Special Notes : 
//
// Creator       : Byron Hanks
//
// Creation Date : 4/01
//
// Owner         : Byron Hanks
//-------------------------------------------------------------------------

#ifndef SOLIDWORKS_QUERY_ENGINE_HPP
#define SOLIDWORKS_QUERY_ENGINE_HPP

#include "GeometryQueryEngine.hpp"

//#include "amapp.h"
#import "C:\Program Files\SolidWorks\sldworks.tlb" raw_interfaces_only, raw_native_types, no_namespace, named_guids     //Import the SolidWorks type library

#include "SWImport.hpp"

//class BodySM;
class Lump;
class ShellSM;
class LoopSM;
class CoEdgeSM;
//class Surface;
//class Curve;
//class Point;

//class SWPart;

class CubitVector;


class SWQueryEngine : public GeometryQueryEngine
{

private:

  SWQueryEngine();

  static SWQueryEngine* instance_;

public:
  inline static SWQueryEngine* instance()
  {
    if( ! instance_ )
    {
      instance_ = new SWQueryEngine;
      instance_->initialize();
    }
    return instance_;
  }

  virtual ~SWQueryEngine();

  CubitStatus initialize();

  ISldWorks *GetSldWorks();
//  LPATTRIBUTEDEF GetAttributeDef();

      virtual CubitStatus export_solid_model( DLIList<TopologyBridge*>& bridge_list,
                                              const char* file_name,
                                              const char* file_type,
                                              const CubitString &cubit_version,
                                              const char* logfile_name = NULL );


     virtual CubitStatus save_temp_geom_file(DLIList<TopologyBridge*> &ref_entity_list,
                                             const char *filename,
                                             const CubitString &cubit_version,
                                             CubitString &created_file,
                                             CubitString &created_file_type );

     virtual CubitStatus import_temp_geom_file(FILE* file_ptr, 
                                            const char* file_name,
                                            const char* file_type,
                                            DLIList<TopologyBridge*> &bridge_list);

     virtual CubitStatus import_solid_model(const char* file_name,
                                            const char* file_type,
                                            DLIList<TopologyBridge*>& imported_entities,
                                            bool print_results = CUBIT_TRUE,
                                            const char* logfile_name = NULL,
                                            bool heal_step = CUBIT_TRUE,
                                            bool import_bodies = CUBIT_TRUE,
                                            bool import_surfaces = CUBIT_TRUE,
                                            bool import_curves = CUBIT_TRUE,
                                            bool import_vertices = CUBIT_TRUE,
                                            bool free_surfaces = CUBIT_TRUE
                                            );

      virtual CubitStatus get_intersections(
                                    Curve* ref_edge1, CubitVector& point1,
                                    CubitVector& point2,
                                    DLIList<CubitVector*>& intersection_list,
                                    bool bounded = false,
                                    bool closest = false );

      virtual CubitStatus get_intersections( Curve* ref_edge1, Curve* ref_edge2,
                                             DLIList<CubitVector*>& intersection_list,
                                             bool bounded = false,
                                             bool closest = false );

      virtual CubitStatus get_intersections( Curve* curve, Surface* surface,
                                             DLIList<CubitVector*>& intersection_list,
                                             bool bounded = false );
      //- Finds the intersections of the curve and surface.  The curve is extended
      //- if the bounded flag is not false.  The function allocates the CubitVectors 
      //- in the returned list, so be sure to free them.

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

      virtual CubitStatus entity_entity_distance( GeometryEntity *ref_entity_ptr1,
                                                  GeometryEntity *ref_entity_ptr2,
                                                  CubitVector &pos1, 
                                                  CubitVector &pos2,
                                                  double &distance );
      //- Gets the minimum distance between two entities and the closest positions 
      //- on those entities. Supports vertices, curves, surfaces, volumes and bodies.

      virtual void delete_solid_model_entities(DLIList<BodySM*>& body_list) const;
  virtual CubitStatus delete_solid_model_entities( BodySM* body_ptr ) const;
  virtual CubitStatus delete_solid_model_entities(Surface* surf_ptr ) const;
  virtual CubitStatus delete_solid_model_entities( Curve* curve_ptr ) const;
  virtual CubitStatus delete_solid_model_entities( Point* point_ptr ) const;
      //- Deletes the solid model entities associcated with the input 
      //- free-floating RefEntity. If the input RefEntity is not free-floating,
      //- then its underlying ACIS ENTITYs are not deleted and CUBIT_FAILURE
      //- is returned.
      //- Returns CUBIT_SUCCESS if all went well, otherwise, CUBIT_FAILURE.

      virtual CubitStatus fire_ray(BodySM *body,
                                   const CubitVector &ray_start_point,
                                   const CubitVector &unit_direction,
                                   DLIList<double> &ray_params,
                                   DLIList<GeometryEntity*> *entity_list = NULL) const;

//  virtual CubitStatus is_point_in_body(Body *body,
//                                       const CubitVector &point_coords,
//                                       CubitPointContainment &is_point_in) const;
    //- Find whether a point is found in a given body, is_point_in equals
    //- -1 for undetermined, 1 for inside, 2 for boundary, and 0 for outside
  
//  virtual CubitBoolean get_mass_props(Body* body, CubitVector &cofg) const;
    //- Get certain mass properties, cofg is center of gravity
  
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
    CubitVector &transform_vector ) const;

  virtual const char* modeler_type();
  virtual int get_major_version();
  virtual int get_minor_version();
  virtual int get_subminor_version();
  virtual int get_allint_version();
  virtual CubitString get_engine_version_string();

  virtual bool is_intermediate_engine() {return FALSE;}

    //- pass a string back identifying this query engine
  virtual CubitStatus set_export_version(const int major,
                                         const int minor);
  virtual CubitStatus set_export_allint_version(const int version);
  virtual CubitStatus list_engine_versions(CubitString &versions);
    //- gets, sets, and lists choices for solid modeling engine version(s)

  virtual double get_sme_resabs_tolerance() const; // Gets solid modeler's resolution absolute tolerance
  virtual double set_sme_resabs_tolerance( double new_resabs );

  virtual CubitStatus set_int_option( const char* opt_name, int val );
  virtual CubitStatus set_dbl_option( const char* opt_name, double val );
  virtual CubitStatus set_str_option( const char* opt_name, const char* val );

//  virtual CubitString identify();

  virtual CubitStatus get_graphics( Surface* surface_ptr,
                                    int& number_of_triangles,
                                    int& number_of_points,
                                    int& number_of_facets,
                                    GMem* gMem,
                                    unsigned short normal_tolerance=15,
                                    double distance_tolerance=0,
                                    double longest_edge = 0 ) const;

  virtual CubitStatus get_graphics( Curve* curve_ptr,
                                    int& num_points,
                                    GMem* gMem = NULL,
                                    double tolerance = 0.0 ) const;

  virtual void reset();

  virtual CubitStatus translate( BodySM* body, const CubitVector& offset );
  virtual CubitStatus rotate   ( BodySM* body, const CubitVector& axis, double angle );
  virtual CubitStatus scale    ( BodySM* body, double factor );
  virtual CubitStatus scale    ( BodySM* body, const CubitVector& factors );
  virtual CubitStatus reflect  ( BodySM* body, const CubitVector& axis );
  virtual CubitStatus restore_transform( BodySM* body );

  virtual CubitStatus translate( GeometryEntity* ent, const CubitVector& offset );
  virtual CubitStatus rotate   ( GeometryEntity* ent, const CubitVector& axis, double degrees );
  virtual CubitStatus scale    ( GeometryEntity* ent, double factor );
  virtual CubitStatus scale    ( GeometryEntity* ent, const CubitVector& factors );
  virtual CubitStatus reflect  ( GeometryEntity* ent, const CubitVector& axis );

  virtual CubitBoolean bodies_overlap (BodySM *body_ptr_1,
                                       BodySM *body_ptr_2 ) const;

public:
  CubitStatus bodysms(IUnknown *entity, DLIList<BodySM*> &bodies) const;
  CubitStatus lumps(SWPart *pSWPart, IUnknown *entity, DLIList<Lump*> &lumps) const;
  CubitStatus shellsms(SWPart *pSWPart, IUnknown *entity, DLIList<ShellSM*> &shellsms) const;
  CubitStatus surfaces(SWPart *pSWPart, IUnknown *entity, DLIList<Surface*> &surfaces) const;
  CubitStatus loopsms(SWPart *pSWPart, IUnknown *entity, DLIList<LoopSM*> &loopsms) const;
  CubitStatus curves(SWPart *pSWPart, IUnknown *entity, DLIList<Curve*> &curves) const;
  CubitStatus coedgesms(SWPart *pSWPart, IUnknown *entity, DLIList<CoEdgeSM*> &coedgesms) const;
  CubitStatus points(SWPart *pSWPart, IUnknown *entity, DLIList<Point*> &points) const;

private:
  //- topology traversal of TB's; implemented based on native traversal
  //- functions in the modeler; returns lists of SW entities - interface
  //- need to be released by the calling function.
  CubitStatus get_ENTITY_from_ENTITY(IID queryIID,
                                     IUnknown *entity,
                                     DLIList<IUnknown *> &entities) const;

private:
    ISldWorks *m_pSWApp;
    IAttributeDef *m_pAttDef; // attribute definition

    SWImport *m_pImport;
} ;

#endif

