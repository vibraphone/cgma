#ifndef ACIS_SURFACE_TOOL_HPP
#define ACIS_SURFACE_TOOL_HPP

//class RefEntity;

#include "CubitDefines.h"
#include "DLIList.hpp"

class CubitVector;
class EDGE;
template <class X> class DLIList;

class Curve;
class Surface;
class BodySM;
class BodyACIS;
class SurfaceACIS;
class CurveACIS;
class TopologyBridge;

class TopologyEntity;
class Body;
class RefFace;
class RefEdge;

class AcisSurfaceTool
{
public:
// ********** BEGIN FRIEND DECLARATIONS        **********

// ********** END FRIEND DECLARATIONS        **********

  ~AcisSurfaceTool();

  static AcisSurfaceTool* instance();
  //- Gives access to the singleton object of this class

  CubitStatus create_net_surface( BodySM *& new_body,
                                  DLIList<DLIList<CubitVector*>*> &vec_lists_u, 
                                  DLIList<DLIList<CubitVector*>*> &vec_lists_v, 
                                  double net_tol = 1e-3,
                                  CubitBoolean heal = CUBIT_TRUE,
                                  CubitBoolean verbose = CUBIT_TRUE);
  CubitStatus create_net_surface( Body *& new_body,
                                  DLIList<DLIList<CubitVector*>*> &vec_lists_u, 
                                  DLIList<DLIList<CubitVector*>*> &vec_lists_v, 
                                  double net_tol = 1e-3,
                                  CubitBoolean heal = CUBIT_TRUE,
                                  CubitBoolean verbose = CUBIT_TRUE);
  //- Creates a net surface using an existing MAPPED surfaces.  Surface
  //- must be meshed or failure occurs.  Useful for approximating a
  //- surface through a set of meshed composite surfaces, or just a group
  //- of mapped surfaces.  If using a set of mapped surfaces, the collection
  //- must form a logical rectangle.  The created RefFace will be in a sheet body.  
  //- Curves must intersect each other within the tolerance given.  From experience, 
  //- bad coedges result, so we heal the resultant body.  Hopefully this will be 
  //- fixed sometime...(SRS, 6-22-99)

  CubitStatus create_net_surface( DLIList<Curve*>& u_curves, 
                                  DLIList<Curve*>& v_curves,
                                  BodySM *& new_body,
                                  double net_tol = 1e-3, 
                                  CubitBoolean heal = CUBIT_TRUE,
                                  CubitBoolean verbose = CUBIT_TRUE);
  CubitStatus create_net_surface( DLIList<RefEdge*>& u_curves,
                                  DLIList<RefEdge*>& v_curves,
                                  Body *& new_body,
                                  double net_tol = 1e-3, 
                                  CubitBoolean heal = CUBIT_TRUE,
                                  CubitBoolean verbose = CUBIT_TRUE);
  //- Creates a net surface through all the u and v curves given.  Curves 
  //- must intersect each other within the tolerance given.  From experience, 
  //- bad coedges result, so we heal the resultant body.  Hopefully this will 
  //- be fixed sometime... (SRS, 6-22-99)

  CubitStatus create_offset_surface( Surface* ref_face_ptr, 
                                     BodySM*& new_body, 
                                     double offset_distance );
  CubitStatus create_offset_surface( RefFace* ref_face_ptr, Body*& new_body, 
                                     double offset_distance );
  //- Offsets one surface, without regard to any others.  Returned surface
  //- is in a sheet body.

  CubitStatus create_offset_body( BodySM* body_ptr, BodySM*& new_body, 
                                  double offset_distance );
  CubitStatus create_offset_body( Body* body_ptr, Body*& new_body, 
                                  double offset_distance );
  //- Creates a new body by offsetting from another body

  CubitStatus create_skin_surface( DLIList<Curve*>& curves, BodySM*& new_body );
  CubitStatus create_skin_surface( DLIList<RefEdge*>& curves, Body*& new_body );
  //- Skinning puts a surface through a set of curves that are somewhat
  //- parallel to each other.

  CubitStatus loft_surfaces( Surface *face1, const double &takeoff1,
                             Surface *face2, const double &takeoff2,
                             BodySM*& new_body,
                             CubitBoolean arc_length_option = CUBIT_FALSE,
                             CubitBoolean twist_option = CUBIT_FALSE,
                             CubitBoolean align_direction = CUBIT_TRUE,
                             CubitBoolean perpendicular = CUBIT_TRUE,
                             CubitBoolean simplify_option = CUBIT_FALSE);
  CubitStatus loft_surfaces( RefFace *face1, const double &takeoff1,
                             RefFace *face2, const double &takeoff2,
                             Body*& new_body,
                             CubitBoolean arc_length_option = CUBIT_FALSE,
                             CubitBoolean twist_option = CUBIT_FALSE,
                             CubitBoolean align_direction = CUBIT_TRUE,
                             CubitBoolean perpendicular = CUBIT_TRUE,
                             CubitBoolean simplify_option = CUBIT_FALSE);

	// loft face1 to face2
  CubitStatus loft_surfaces_to_body( 
                   Surface *face1, const double &takeoff1,
									 Surface *face2, const double &takeoff2,
									 BodySM*& new_body,
									 CubitBoolean arc_length_option,
									 CubitBoolean twist_option,
									 CubitBoolean align_direction,
									 CubitBoolean perpendicular,
									 CubitBoolean simplify_option);
   CubitStatus loft_surfaces_to_body( RefFace *face1, const double &takeoff1,
									 RefFace *face2, const double &takeoff2,
									 Body*& new_body,
									 CubitBoolean arc_length_option,
									 CubitBoolean twist_option,
									 CubitBoolean align_direction,
									 CubitBoolean perpendicular,
									 CubitBoolean simplify_option);
   
	// loft surface to surface, creating new body
  CubitStatus create_surface( DLIList<CubitVector*>& vec_list,
                              BodySM *&new_body,
                              Surface *surface_ptr, 
			                        CubitBoolean project_points );
  CubitStatus create_surface( DLIList<CubitVector*>& vec_list,
                              Body *&new_body,
                              RefFace *ref_face_ptr, 
			                        CubitBoolean project_points );
  //- Creates a surface out of a list of points.  The points must form
  //- a polygon (there must be at least 3 points), but do not need to be 
  //- absolutely planar.  The resultant surface will have one loop with
  //- straight curves connecting the points.  If an existing RefFace is 
  //- passed-in the created surface will use it's geometry.  If project_
  //- points is true the points will be projected to the surface then the
  //- new surface will be created, otherwise the existing points must lie
  //- on the surface.

  CubitStatus create_weld_surface( CubitVector &root, 
                                   Surface *surface_ptr1, double leg1, 
                                   Surface *surface_ptr2, double leg2, 
                                   BodySM *&new_body );
  CubitStatus create_weld_surface( CubitVector &root, RefFace *ref_face1, 
    double leg1, RefFace *ref_face2, double leg2, Body *&new_body );
  //- Creates a triangular surface utilizing a root point, two directions
  //- based on surface normal directions and two lengths.  

protected:
   AcisSurfaceTool();
   //- Class Constructor. (Not callable by user code. Class is constructed
   //- by the {instance()} member function.

private:

  static BodyACIS* get_acis_body( Body* );
  static SurfaceACIS* get_acis_surface( RefFace* );
  static CurveACIS* get_acis_curve( RefEdge* );
  static TopologyBridge* get_acis_bridge( TopologyEntity* );

  CubitStatus create_surface( DLIList<CubitVector*>& vec_list,
                              BodySM *&new_body,
                              Surface *surface_ptr = NULL );
  //- Creates a surface out of a list of points.  The points must form
  //- a polygon (there must be at least 3 points), but do not need to be 
  //- absolutely planar.  The resultant surface will have one loop with
  //- straight curves connecting the points.  If an existing RefFace is 
  //- passed-in the created surface will use it's geometry (the calling 
  //- code should make sure the points like on the surface).

  CubitStatus create_EDGE_on_Surface( const CubitVector &start_in, 
                                      const CubitVector &end_in,
                                      SurfaceACIS* surface_ptr, 
                                      EDGE *&EDGE_ptr );
  //- Creates an EDGE that lies on the given RefFace.

  void extend_weld_EDGE( EDGE *&EDGE_ptr, double desired_length );
  //- Extends or contracts the EDGE to be the desired length.  Note that the 
  //- EDGE_ptr may be modified.  Assumes edge parameter space is well-defined -
  //- if not, the desired length may not be achieved.

  static AcisSurfaceTool* instance_;
};

#endif
