#ifndef ACIS_QUERY_ENGINE_HPP
#define ACIS_QUERY_ENGINE_HPP
#include <typeinfo>

#if !defined(NT) && !defined(CANT_USE_STD)
using std::type_info;
#endif
#include "GeometryQueryEngine.hpp"
class TopologyBridge;

class GeometryEntity;
class AssemblySM;
class PartSM;
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

template <class X> class DLIList;
class GMem;

class AcisQueryEngine : public GeometryQueryEngine
{
public:
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
  virtual const type_info& entity_type_info() const ;
  virtual CubitStatus export_solid_model(DLIList<TopologyBridge*>& bridge_list, const char* file_name, const char* file_type, const CubitString &cubit_version, const char* logfile_name = NULL ) ;
  virtual CubitStatus save_temp_geom_file(DLIList<TopologyBridge*> &ref_entity_list, const char *filename, const CubitString &cubit_version, CubitString &created_file, CubitString &created_file_type );
  virtual CubitStatus import_temp_geom_file(FILE* file_ptr, const char* file_name, const char* file_type, DLIList<TopologyBridge*> &bridge_list );
//  virtual CubitStatus import_solid_model(FILE* file_ptr, const char* file_name, const char* file_type, DLIList<TopologyBridge*>& imported_entities, CubitBoolean print_results = CUBIT_TRUE, const char* logfile_name = NULL, CubitBoolean heal_step = CUBIT_TRUE, CubitBoolean import_bodies = CUBIT_TRUE, CubitBoolean import_surfaces = CUBIT_TRUE, CubitBoolean import_curves = CUBIT_TRUE, CubitBoolean import_vertices = CUBIT_TRUE, CubitBoolean free_surfaces = CUBIT_TRUE ) ;
  virtual CubitStatus import_solid_model(SolidModelerType model_type) const;
  virtual void delete_solid_model_entities(DLIList<BodySM*>& body_list) const ;
  virtual CubitStatus delete_solid_model_entities(AssemblySM* ref_assembly_ptr) const;
  virtual CubitStatus delete_solid_model_entities(PartSM* ref_part_ptr) const;
  virtual CubitStatus delete_solid_model_entities(GeometryEntity* ptr, bool delete_lower_order) const;
  virtual CubitStatus delete_solid_model_entities( BodySM* body_ptr ) const;
  virtual CubitStatus delete_solid_model_entities(Surface* surf_ptr ) const;
  virtual CubitStatus delete_solid_model_entities( Curve* curve_ptr ) const;
  virtual CubitStatus delete_solid_model_entities( Point* point_ptr ) const;
  virtual CubitStatus get_graphics(Surface*, int&, int&, int&, GMem*, short unsigned int, double, double) const;
  
  virtual CubitStatus get_graphics(Surface*, int&, int&, int&, GMem*, unsigned short, double) const;
  virtual CubitStatus get_graphics( Curve* curve_ptr, int& num_points, GMem* gMem = NULL, double tolerance = 0.0 ) const;
  virtual CubitStatus get_isoparametric_points(Surface* ref_face_ptr, int &nu, int &nv, GMem *&gMem) const;
  virtual CubitStatus get_u_isoparametric_points(Surface* ref_face_ptr, double v, int& n, GMem *&gMem) const;
  virtual CubitStatus get_v_isoparametric_points(Surface* ref_face_ptr, double u, int&n, GMem *&gMem) const;
  virtual CubitStatus get_intersections(Curve* ref_edge1, CubitVector& point1, CubitVector& point2, DLIList<CubitVector*>& intersection_list, bool bounded = false, bool closest = false );
  virtual CubitStatus get_intersections(Curve* ref_edge1, Curve* ref_edge2, DLIList<CubitVector*>& intersection_list, bool bounded = false, bool closest = false );
  virtual CubitStatus get_intersections(Curve* ref_edge, Surface* ref_face, DLIList<CubitVector*>& intersection_list, bool bounded = false );
  virtual CubitStatus entity_extrema( DLIList<GeometryEntity*> &entity_list, const CubitVector *dir1, const CubitVector *dir2, const CubitVector *dir3, CubitVector &extrema, GeometryEntity *&extrema_entity_ptr );
  virtual CubitStatus entity_entity_distance( GeometryEntity *ref_entity_ptr1, GeometryEntity *ref_entity_ptr2, CubitVector &pos1, CubitVector &pos2, double &distance );
  virtual CubitStatus fire_ray( BodySM *body, const CubitVector &ray_point, const CubitVector &unit, DLIList<double>& ray_params, DLIList<GeometryEntity*> *entity_list) const;
  virtual CubitStatus transform_vec_position(CubitVector const& position_vector, BodySM *OSME_ptr, CubitVector &transform_vector ) const;
  virtual int get_major_version();
  virtual int get_minor_version();
  virtual int get_subminor_version();
  virtual int get_allint_version();
  virtual CubitString get_engine_version_string();
  virtual CubitStatus set_export_version(const int major, const int minor);
  virtual CubitStatus set_export_allint_version(const int version);
  virtual CubitStatus list_engine_versions(CubitString &versions);
  virtual double get_sme_resabs_tolerance() const; 
  virtual double set_sme_resabs_tolerance( double new_resabs );
  virtual CubitStatus set_steptools_path( const char* path );
  virtual CubitStatus set_igestools_path( const char* path );
  virtual CubitStatus set_int_option( const char* opt_name, int val );
  virtual CubitStatus set_dbl_option( const char* opt_name, double val );
  virtual CubitStatus set_str_option( const char* opt_name, const char* val );
  CubitStatus import_datum_curves(FILE* inputFile,  const char* filetype ) const;
  CubitStatus get_surfs_on_plane( BodySM* body_ptr, const CubitVector& pln_orig, const CubitVector& pln_norm, DLIList<Surface*>& ref_face_list ) const;
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

  CubitStatus import_solid_model(const char*, const char*, DLIList<TopologyBridge*>&, bool, const char*, bool, bool, bool, bool, bool, bool);
  virtual CubitBoolean bodies_overlap(BodySM*, BodySM*) const;
  virtual CubitBoolean volumes_overlap (Lump *lump1, Lump *lump2 ) const;
};

#endif
