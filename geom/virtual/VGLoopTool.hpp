//-------------------------------------------------------------------------
// Filename      : VGLoopTool.hpp
//
// Purpose       : Loop Modification Code
//
// Special Notes : Provides common functionality for composite and 
//                 partition geometry.
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 08/08/02
//-------------------------------------------------------------------------

#ifndef VG_LOOP_TOOL
#define VG_LOOP_TOOL

 
#include "CubitDefines.h"
#include "DLIList.hpp"
#include <vector>

class CubitVector;

template <class SURFACE, class LOOP, class COEDGE,class CURVE, class POINT>
class VGLoopTool 
{
public:
/*  
  static CubitStatus remove_curve( COEDGE* coedge1, 
                                   COEDGE* coedge2,
                                   DLIList<COEDGE*>& new_loop_coedges );
    //- Remove a curve from loop topology
    //- coedge1/2:  The two coedges of a curve
    //- If the passed coedges are in different loops, the loops
    //- will be combined into a single loop and new_loop_coedges will
    //- be empty.  
    //- If the passed coedges are in the same loop, unless the coedges
    //- are adjacent in that loop, the loop will be split into two
    //- loops and the coedges for the new loop passed back in 
    //- new_loop_coedges.

  static CubitStatus insert_curve( COEDGE* new_coedge1, 
                                   COEDGE* new_coedge2,
                                   SURFACE* surface,
                                   LOOP*& new_loop1,
                                   LOOP*& new_loop2 );
    //- Insert a curve into the loop topology
    //- new_coedge1/2:  The two new coedges of the curve to be inserted
    //- loop1/2_coedge: Any coedges of the existing loops in which the
    //-                 new coedges are to be inserted.  (May be NULL).
    //- new_loop_coedges: The ordered list of coedges from which the
    //-                   caller should create a new loop.  The loops
    //-                   are specified with a coedge to avoid the need
    //-                   for this class to also be a template of the
    //-                   loop type.
    //- Note:  In the comments below, when it is stated that a new loop
    //-        is created, this means that the coedges for the new loop
    //-        are passed back in new_loop_coedges.  This function does
    //-        not do the actual loop creation.
    //-
    //- If both loop1_coedge and loop2_coedge are null, both coedges
    //- are passed back in new_loop_coedges without any topology
    //- modification.  The caller should create either one or two
    //- new loops (a curve imprinted on the interior of the surface
    //- or a hole in the surface and a new surface in the hole)
    //- depending on the curve having the same start and end vertices.
    //-
    //- If one of loop1_coedge and loop2_coedge is NULL, or both
    //- are part of the same loop, both new coedges are inserted in 
    //- the same loop.  This will result in either a sipe or the 
    //- loop being split into two loops.  If the loop is split, the 
    //- coedges for the new loop are passed back in new_loop_coedges.
    //-
    //- If both loop1_coedge and loop2_coedge are non-null and part
    //- of different loops, the loops are joined using the new 
    //- coedges as a bridge.
*/
  static double loop_angle_metric( COEDGE* first_coedge );
/*  
  static bool is_position_within_loop( const CubitVector& position,
                                       COEDGE* first_coedge,
                                       CubitVector* closest_on_loop = 0,
                                       CubitVector* normal = 0 );
                                       
  static double loop_area( COEDGE* first_coedge );

  static CubitStatus closest_point_trimmed( 
                                      DLIList<COEDGE*>& surface_coedges,
                                      const CubitVector& from_pt,
                                      CubitVector& result_pt,
                                      CubitPointContainment* cont = 0 );
*/
  static void get_loop_polyline( COEDGE* first_coedge,
                                 std::vector<CubitVector>& result_set );

private:
/*  
  static CubitStatus previous_coedge( COEDGE* coedge, SURFACE* surface,
                                      COEDGE*& result );
    //- Find the coedge in the passed surface after which the passed
    //- coedge will be inserted.  The returned coedge will share its
    //- start point (end point of the curve if the coedge sense is 
    //- reversed) with the end point (start point of the curve if the
    //- coedge sense is reversed) of the passed coedge.  If the coedge
    //- has no points in common with any loops in the passed surface,
    //- the result is NULL and CUBIT_SUCCESS is returned.  CUBIT_FAILURE
    //- is returned if something is wrong with the topology.
  
  static COEDGE* closest_loop_coedge( const CubitVector& from_pt,
                                      COEDGE* first_coedge,
                                      COEDGE*& other_coedge,
                                      CubitVector* closest_pt = 0 );

  static CubitStatus closest_surface_coedge( SURFACE* surface,
                                         const CubitVector& from_pt,
                                         COEDGE*& coedge1,
                                         COEDGE*& coedge2,
                                         CubitVector* pt_on_curve = 0 );

  static bool inside_of_curve( const CubitVector& curve_tangent,
                               const CubitVector& curve_position,
                               const CubitVector& surface_position,
                               const CubitVector& surface_normal );                                     

  static bool odd_coedge_count( COEDGE* coedge );
*/  
  static CubitStatus get_loop_polygon( COEDGE* first_coedge,
                                       DLIList<CubitVector*>& polygon_points );
                                       
};

#include "VGLoopTool.cpp"

#endif

