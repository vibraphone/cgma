#ifndef REF_EDGE_CCAPI_HPP
#define REF_EDGE_CCAPI_HPP

#include "EntityType.h"
#include "CubitDefines.h"
#include "CubitVectorStruct.h"

#ifdef __cplusplus
extern "C" {
#endif

  void *RefEdge_get_address(void *this_ref_edge,
                            enum EntityType inputEntityType);
    /* R void **/
    /* R- Returned void pointer */
    /* I inputenum EntityType */
    /* I- The input type to get the address of. */
    /* - This function returns a void pointer that points to the */
    /* - "appropriate" portion of this object.  The appropriate */
    /* - portion is determined by the input enum EntityType variable. */
    /* - Returns NULL if the input type and the type of "this" */
    /* - are not related by inheritance. */
    /* - Note: The RTTI capabilities encoded in these functions */
    /* -       are designed to work with any form of multiple  */
    /* -       inheritance, as well.  Multiple inheritance is what */
    /* -       necessitates having this function defined in every  */
    /* -       class in the hierarchy. */
    /* - Note: This function can also be used to merely check if */
    /* -       an object of one type is related to another type */
    /* -       through inheritance.  If a non-NULL pointer is */
    /* -       returned, then this is true. */
      
  enum CubitStatus RefEdge_get_point_direction(void *this_ref_edge,
                                                 /* CubitVector & */ struct CubitVectorStruct *origin, 
                                                 /* CubitVector & */ struct CubitVectorStruct *direction );
    /* - Only valid for straight lines */
    /* - Finds the underlying line's origin and direction unit vector */
    /* - Returns CUBIT_FAILURE if curve is not a line */

  enum CubitStatus RefEdge_get_center_radius(void *this_ref_edge,
                                               /* CubitVector & */ struct CubitVectorStruct *center,
                                               /* double& */ double *radius );
    /* - Only valid for arcs */
    /* - Finds the underlying arc's center point and radius */
    /* - Returns CUBIT_FAILURE if curve is not an arc */
  
  enum EntityType RefEdge_entity_type(void *this_ref_edge); 

  enum EntityType RefEdge_get_child_ref_entity_type(void *this_ref_edge);
    /* R enum EntityType */
    /* R- A type value. */
    /* - This function returns the type of the child RefEntity of   */
    /* - RefEdge, which is the type value of the RefVertex class. */
      
  enum EntityType RefEdge_get_parent_ref_entity_type(void *this_ref_edge);
    /* R enum EntityType */
    /* R- A type value. */
    /* - This function returns the type of the parent RefEntity of   */
    /* - RefEdge, which is the type value of the RefFace class. */
      
  enum EntityType RefEdge_get_topology_bridge_type(void *this_ref_edge); 
    /* R enum EntityType */
    /* R- The type of GeometryEntity associated with this object */
    /* - This function returns the type of GeometryEntity associated with  */
    /* - this BasicTopologyEntity.  */
      
    /*  ********* TOPOLOGY ******************** */
  
    enum EntityType RefEdge_get_grouping_entity_type(void *this_ref_edge);
    /* R enum EntityType */
    /* R- The type of the corresponding GroupingEntity type. */
    /* - This function returns the type of the corresponding */
    /* - GroupingEntity . */

  enum EntityType RefEdge_get_sense_entity_type(void *this_ref_edge);  
    /* R enum EntityType */
    /* R- The type of SenseEntity associated with this object */
    /* - This function returns the type of SenseEntity associated with  */
    /* - this BasicTopologyEntity.  */

    /* RefVertex * */ void *RefEdge_startRefVertex(void *this_ref_edge);
    /* RefVertex * */ void *RefEdge_endRefVertex(void *this_ref_edge);
    
  void RefEdge_switch_vertices(void *this_ref_edge);  
    /* R RefVertex **/
    /* R- Returned RefVertex pointer */
    /* - These functions get the start and end RefVertex'es of this */
    /* - RefEdge. */
    /* - */
    /* - MJP Notes: */
    /* - The start and end RefVertices are cached in the RefEdge (see */
    /* - private member data, start/end_RefVertex). These need to be */
    /* - updated if the end RefVertices could have changed during a */
    /* - geometric operation.  See the routine, RefEdge::update to */
    /* - see how this is done. */

    /* Chain * */ void *RefEdge_get_first_chain_ptr(void *this_ref_edge);
    /* R Chain **/
    /* R- Pointer to the first Chain */
    /* - This function returns a pointer to the first Chain that this */
    /* - RefEdge points to. */

  enum CubitStatus RefEdge_get_co_edges(void *this_ref_edge,
                                          /* DLCoEdgeList & */ void ***co_edges_found_list, int *co_edges_found_list_size,
                                          /* RefFace * */ void *input_ref_face_ptr);
    /* R enum CubitStatus  */
    /* R-CubitSucces/CubitFailure */
    /* O co_edges_found_list and optional input_ref_face_ptr */
    /* O- Gets the coedges that */
    /* - correspond to "this" RefEdge, they are populated into the */
    /* - co_edges_found_list. If an input_ref_face_ptr is sent it */
    /* - gets the co_edges that are just associated with that ref_face. */
    /* - Remember that usually there will be only one but for hard lines */
    /* - and sipes there may be two co-edges/ref-edge/ref-face. */
    /* - Note: This function will just blindly append the co-edges it */
    /* - finds into the co_edges_found_list.  A merge-unique might be */
    /* - more expensive than we want. */

  double RefEdge_angle_between(void *this_ref_edge,
                                 /* RefEdge * */ void *other_edge_ptr, /* RefFace * */ void *face_ptr );
    /* - Returns the "inside" angle between this edge and */
    /*  other_edge_ptr.  Figures out which edge comes first.  Only */
    /*  requirement is that both edges must be on face_ptr and on */
    /*  the same loop.  Angle is in radians. */
  
  enum CubitSense RefEdge_sense(void *this_ref_edge,
                             /* RefFace * */ void *face );
    /* - Sense of this edge wrt the given face. */
    /* - Return CUBIT_FORWARD or CUBIT_REVERSE. */
    /* - Returns CUBIT_UNKNOWN if there is more than one coedge, and */
    /* - senses don't agree, or if there are no coedges. */
    
  int RefEdge_dimension(void *this_ref_edge); 
    /* - returns dimension of the actual entity. */

    /* RefVertex * */ void *RefEdge_common_ref_vertex_1(void *this_ref_edge,
                                                          /* RefEdge * */ void *otherRefEdgePtr );
    /* - Find common RefVertex */

    /* RefVertex * */ void *RefEdge_common_ref_vertex_2(void *this_ref_edge,
                                                          /* RefEdge * */ void *next_ref_edge,
                                                          /* RefFace * */ void *ref_face_ptr );
    /* - Finds the common vertex between the next_ref_edge and */
    /* - this edge, assuming that this edge comes before next_ref_edge */
    /* - in a the loop around the surface. This is needed because */
    /* - two edges may share two vertices... */

  enum CubitBoolean RefEdge_common_vertices(void *this_ref_edge,
                                              /* RefEdge * */ void *otherRefEdgePtr, 
                                              /* DLRefVertexList & */ void ***common_verts, int *common_verts_size);
    /* -Populates the common_verts with the vertices that are common */
    /* -between the two curves. */
    /* -Returns CUBIT_TRUE if there are common vertices/ false otherwise. */

    /* RefEdge * */ void *RefEdge_get_other_curve(void *this_ref_edge,
                                                    /* RefVertex * */ void *common_vertex,
                                                    /* RefFace * */ void *ref_face_ptr);
    /* - Finds the other curve on the ref_face_ptr that shares */
    /* - the common vertex */
  
  enum CubitStatus RefEdge_get_two_co_edges(void *this_ref_edge,
                                              /* RefEdge * */ void *next_ref_edge,
                                              /* RefFace * */ void *ref_face_ptr,
                                              /* CoEdge *& */ void **co_edge_this,
                                              /* CoEdge *& */ void **co_edge_next );
    /* R enum CubitStatus  */
    /* R- CUBIT_SUCCESS/CUBIT_FAILURE */
    /* O next_ref_edge, ref_face_ptr, co_edge_this, co_edge_next */
    /* O-Returns the co_edge that corrisponds to 'this' ref_edge */
    /* - and the one that corrisponds to next_ref_edge, with */
    /* - respect to the ref_face_ptr. */
    /* - Special Note : next_ref_edge must follow 'this' ref_edge in a Loop */
    /* - on the ref_face_ptr, this is assumed so the function */
    /* - will assert if this is not done...   */

    /* RefFace * */ void *RefEdge_other_face(void *this_ref_edge,
/* RefFace * */ void *not_face, /* RefVolume * */ void *ref_volume);
    /* - return the (an) other face sharing this edge, which also borders */
    /* - ref_volume if non-NULL */

    /* RefVertex * */ void *RefEdge_other_vertex(void *this_ref_edge,
                                                   /* RefVertex * */ void *refVertexPtr );
    /* - Returns the vertex at the other end of the edge. */

    /* RefVertex * */ void *RefEdge_closest_vertex(void *this_ref_edge,
                                                   const /* CubitVector & */ struct CubitVectorStruct *point3);
  
    /*  ********* GEOMETRY ******************** */
  
    /* Curve * */ void *RefEdge_get_curve_ptr(void *this_ref_edge) ;
    /* R Curve **/
    /* R- A pointer to the Curve to which this RefEdge points.  */
    /* - This function returns a pointer to the Curve */
    /* - to which the current RefEdge points. */
      
  struct CubitVectorStruct RefEdge_start_coordinates(void *this_ref_edge);
  struct CubitVectorStruct RefEdge_end_coordinates(void *this_ref_edge);
    /* R CubitVector */
    /* R- Returned location. */
    /* - These functions return the start and end global coordinate */
    /* - locations of the RefEdge. */
    /* - */
    /* - NOTE: */
    /* - These coordinates remain consistent throughout the life of the */
    /* - RefEdge, regardless of the fact that the actual start and */
    /* - end RefVerex'es may change (e.g., after a merge operation). */
      
  void RefEdge_move_to_curve (void *this_ref_edge,
                                /* CubitVector & */ struct CubitVectorStruct *vector );
    /* - Moves the given location (CubitVector or CubitNode) to the closest  */
    /* - point on the Curve */

  enum CubitStatus RefEdge_get_interior_extrema(void *this_ref_edge,
                                                  /* DLCubitVectorList& */ struct CubitVectorStruct **interior_points,
                                                enum CubitSense *return_sense);

  enum CubitStatus RefEdge_closest_point_trimmed(void *this_ref_edge,
                                                   /* CubitVector const& */ struct CubitVectorStruct *location, 
                                                   /* CubitVector & */ struct CubitVectorStruct *closest_location);
    /* R void     */
    /* I location */
    /* I- The point to which the closest point on the Curve is desired. */
    /* O closest_location */
    /* O- The point on the Curve, closest to the input location which */
    /* O- will be on the Curve.  This is input as a reference  */
    /* O- so that the function can modify its contents. */

  
  enum CubitStatus RefEdge_closest_point(void *this_ref_edge,
                                           /* CubitVector const& */ struct CubitVectorStruct *location, 
                                           /* CubitVector & */ struct CubitVectorStruct *closest_location,
                                           /* CubitVector * */ struct CubitVectorStruct *tangent_ptr,
                                           /* CubitVector * */ struct CubitVectorStruct *curvature_ptr);
    /* R void     */
    /* I location */
    /* I- The point to which the closest point on the Curve is desired. */
    /* O closest_location */
    /* O- The point on the Curve, closest to the input location which */
    /* O- might not be on the Curve.  This is input as a reference  */
    /* O- so that the function can modify its contents. */
    /* O tangent_ptr */
    /* O- The tangent to the Curve (output as a unit vector) at the  */
    /* O- closest_location. */
    /* O curvature_ptr */
    /* O- The curvature of the Curve at the closest_location. */
    /* - This function computes the point on the Curve closest to the input  */
    /* - location. */
    /* - */
    /* - If the tangent and/or curvature is required, then the calling code */
    /* - is responsible for allocating space for the CubitVector(s) and */
    /* - sending in the relevant non-NULL pointers.  If either of these */
    /* - pointers is NULL, the related quantity is not computed. */
    /* - */
    /* - Notes: */
    /* - The tangent direction is always in the positive direction of the  */
    /* - *owning RefEdge*, regardless of the positive direction of the */
    /* - underlying solid model entities, if any. */
      
  void RefEdge_tangent_1(void *this_ref_edge,
                         const /* CubitVector & */ struct CubitVectorStruct *point, /* CubitVector & */ struct CubitVectorStruct *tangent_vec );
    /* - Return the tangent for the point on this RefEdge that is closest */
    /* - to the input "point".  The tangent direction is always in the */
    /* - positive direction of the RefEdge.  The positive direction of */
    /* - the RefEdge is an invariant through its lifecycle. */


  enum CubitStatus RefEdge_tangent_2(void *this_ref_edge,
                                     const /* CubitVector & */ struct CubitVectorStruct *point,
                                       /* CubitVector & */ struct CubitVectorStruct *tangent_vec,
                                       /* RefFace * */ void *ref_face_ptr );
    /* R enum CubitStatus  */
    /* R- CUBIT_SUCCESS/CUBIT_FAILURE   */
    /* O point, tangent_vec, ref_face_ptr */
    /* O- Get the correct tangent with respect to */
    /* - the ref_face_ptr. */
    /* - This tangent function is not the safest method for getting */
    /* - the tangent on a surface.  In face if there could be */
    /* - another direction, this function could fail...(assert). */
    /* - This function assumes that there is only 1 co-edge per */
    /* - this ref_face for this REfEdge. */
  
  enum CubitStatus RefEdge_tangent_3(void *this_ref_edge,
                                     const /* CubitVector & */ struct CubitVectorStruct *point, /* CubitVector & */ struct CubitVectorStruct *tangent_vec,
                                       /* RefEdge * */ void *next_ref_edge,
                                       /* RefFace * */ void *ref_face_ptr );
    /* - Retruns the tangent for the point on this RefEdge that is closest */
    /* - to the input "point".  Also, because the next_ref_edge and */
    /* - the ref_face_ptr are given, the tangent is oriented correctly. */
    /* - NOTE: It is assumed that 'this' RefEdge and next_ref_edge are */
    /* - in order in a "Loop" sense on the given ref_face_ptr.  The next */
    /* - ref_edge obviously must follow 'this' edge in this loop. */
  
  double RefEdge_measure(void *this_ref_edge);
    /* - A generic geometric extent function. */
    /* - Returns volume for RefVolumes, area for RefFaces, length for RefEdge, */
    /* - and 1.0 for RefVertices */
    /* - A RefGroup calculates the maximum dimension of its contained */
    /* - entities and returns the sum of the measures() of all entities */
    /* - of that dimension. */
    /* - Default return value is 0.0 for all other RefEntities. */

    /* CubitString */ const char *measure_label();
    /* - Returns the type of measure: (volume, area, length, or N/A) */

  double RefEdge_get_arc_length_1(void *this_ref_edge);
  double RefEdge_get_arc_length_2(void *this_ref_edge,
                                  const /* CubitVector & */ struct CubitVectorStruct *point1,
                                  const /* CubitVector & */ struct CubitVectorStruct *point2 );
  double RefEdge_get_arc_length_3(void *this_ref_edge,
                                  const /* CubitVector & */ struct CubitVectorStruct *point1, const int whichEnd );
    /* - Various arc length calculations, some are redundant */

  double RefEdge_get_chord_length(void *this_ref_edge);
    /* - Calculates and returns the straight-line distance */
    /* - between the startRefVertex and the endRefVertex */
  
    /* CubitVector */ struct CubitVectorStruct RefEdge_center_point(void *this_ref_edge);
    /* - Returns location at the "actual" center of this RefEdge (along the  */
    /* - arc of the RefEdge)  */

  enum CubitStatus RefEdge_mid_point_1(void *this_ref_edge,
                                       const /* CubitVector & */ struct CubitVectorStruct *point1, const /* CubitVector & */ struct CubitVectorStruct *point2,
                                         /* CubitVector & */ struct CubitVectorStruct *midPoint );
    /* - Calculate midpoint between the 2 points on this RefEdge */

  enum CubitStatus RefEdge_mid_point_2(void *this_ref_edge,
                                         /* CubitVector & */ struct CubitVectorStruct *mid_point);
    /* - Calculate midpoint on this RefEdge */

  enum CubitStatus RefEdge_position_from_fraction(void *this_ref_edge,
                                                  const double fraction_along_curve,
                                                    /* CubitVector & */ struct CubitVectorStruct *ouput_position );
    /* R enum CubitStatus  */
    /* I fraction in parameter space along refedge. (1/3,2/3,4/5...) */
    /* O- position where percent in parameter space lies. */
    /* -This function takes the given fraction, finds the parameter value, */
    /* -and calculates this position. This is based off the vgi curve. */

  double RefEdge_start_param(void *this_ref_edge);
    /* R double */
    /* - this function returns the starting parameter of the underlying curve */
  
  double RefEdge_end_param(void *this_ref_edge);
    /* R double */
    /* - this function returns the ending parameter of the underlying curve */
  
  enum CubitBoolean RefEdge_get_param_range(void *this_ref_edge,
                                              /* double& */ double *start_param,
                                              /* double& */ double *end_param );
    /* R enum CubitBoolean */
    /* R- CUBIT_TRUE/FALSE */
    /* O start_param, end_param */
    /* O- The "lower" and "upper" parameter values of the RefEdge. */
    /* - This function returns the parameter bounds for the RefEdge. */
    /* - The start_param represents the parameter value of the  */
    /* - start location of the RefEdge and the end_param that of the */
    /* - end location. */
    /* - */
    /* - CUBIT_TRUE is returned if the RefEdge is defined parametrically */
    /* - and CUBIT_FALSE if it is not.  In the latter case, the  */
    /* - output values of start_ and end_param are undetermined. */
    /* - */
    /* - MJP Note: */
    /* - The numercial value of the start_param could be higher */
    /* - than that of the end_param. */

  double RefEdge_u_from_position (void *this_ref_edge,
                                  const /* CubitVector & */ struct CubitVectorStruct *input_position);
    /* R double */
    /* R- The returned "u" parameter value in local parametric space */
    /* I input_position */
    /* I- The input position for which "u" is to be computed. */
    /* - This function returns the coordinate of a point in the local */
    /* - parametric (u) space that corresponds to the input position in */
    /* - global (world) space.  The input point is first moved to the */
    /* - closest point on the Curve and the parameter value of that */
    /* - point is determined. */

  enum CubitStatus RefEdge_position_from_u (void *this_ref_edge,
                                            double u_value,
                                              /* CubitVector & */ struct CubitVectorStruct *output_position);
    /* R enum CubitStatus  */
    /* R- CUBIT_SUCCESS/FAILURE */
    /* I u_value */
    /* I- The input u parameter value */
    /* O output_position */
    /* O- The output position */
    /* - This function returns the coordinates of a point in the global */
    /* - (world) space that corresponds to the input parametric position  */
    /* - in the local space. */
    /* - */
    /* - If the input parameter value is not defined for the Curve, then  */
    /* - the input CubitVector is not modified and CUBIT_FAILURE is */
    /* - returned. Otherwise, position is appropriately modified and */
    /* - CUBIT_SUCCESS is returned. */
    /* - */
    /* - If the curve is periodic, the input u_value is first "normalized" */
    /* - to the fundamental period of the Curve before its position */
    /* - in global space is determined. */

  double RefEdge_u_from_arc_length (void *this_ref_edge,
                                    double root_param,
                                    double arc_length );
    /* R double */
    /* R- Returned parameter value */
    /* I root_param */
    /* I- The parameter value of the "root point" */
    /* I arc_length */
    /* I- A distance along the Curve */
    /* - This function returns the parameter value of the point that is */
    /* - "arc_length" away from the root point in the */
    /* - positive sense direction of the owning RefEdge. */
    /* -  */
    /* - A negative value for distance would force the search to go in the  */
    /* - negative (sense) direction of the RefEdge. */
    /* - */
    /* - NOTE: */
    /* - The important assumption that is made in this routine is that */
    /* - the end points of the RefEdge that owns this CurveACIS are the same */
    /* - as the end points of the first solid model entity in the list of  */
    /* - solid model entities associated with this Curve. */
      
  enum CubitStatus RefEdge_point_from_arc_length (void *this_ref_edge,
                                                  const /* CubitVector & */ struct CubitVectorStruct *root_point, 
                                                  double const arc_length,
                                                    /* CubitVector & */ struct CubitVectorStruct *new_point );
    /* - Return a point arc_length distance from root_point */
    /* - on this RefEdge.  */
    /* - */
    /* - If arc_length is negative, the new point */
    /* - is in the negative sense direction (along the RefEdge) from */
    /* - the root point. */
    /* -  */
    /* - If the curve is not periodic and the point arc_length away */
    /* - from root_point in the appropriate direction goes beyond */
    /* - the end point of the RefEdge, that end point is returned */
    /* - as new_point. */
    /* - */
    /* - If the curve is periodic and the point arc_length away */
    /* - from root_point in the appropriate direction goes beyond */
    /* - the end point of the RefEdge, wrap around is done. */
    /* - */
    /* - NOTE: I have had problems with this function if the root point */
    /* - is the start vertex, for some reason by a factor of 1e-7 the */
    /* - parameter value is different between the start_param and the */
    /* - parameter found from the start_vertex location.  So I will switch */
    /* - EdgeMeshTool to go off the start parameter instead of this function. */
  
  double RefEdge_length_from_u(void *this_ref_edge,
                               double parameter1,
                               double parameter2 );

  enum CubitBoolean RefEdge_is_periodic_1(void *this_ref_edge);
    /* R enum CubitBoolean */
    /* R- CUBIT_TRUE/CUBIT_FALSE */
    /* - This function determines whether the underlying geometry of the */
    /* - Curve is periodic or not.  Returns CUBIT_TRUE if it is and  */
    /* - CUBIT_FALSE if it is not. */
  
  enum CubitBoolean RefEdge_is_periodic_2(void *this_ref_edge,
                                            /* double& */ double *period);
    /* R enum CubitBoolean */
    /* R- CUBIT_TRUE/CUBIT_FALSE */
    /* O period */
    /* O- Returned period value */
    /* - This function determines whether the underlying geometry of the */
    /* - Curve is periodic or not.  Returns CUBIT_TRUE if it is and  */
    /* - CUBIT_FALSE if it is not. */
    /* - */
    /* - If it is periodic, then it returns the period in the input */
    /* - reference variable, "period". This value is set to 0.0 if */
    /* - the Curve is not periodic. */

  int RefEdge_get_mark(void *this_ref_edge);
  void RefEdge_set_mark(void *this_ref_edge,
                        int set_value );
    /* - get/set generic 2-bit mark. */

  enum CubitStatus RefEdge_relative_sense(void *this_ref_edge,
                                            /* RefEdge * */ void *ref_edge_ptr_2,
                                          double tolerance_factor,
                                          enum CubitSense *sense,
                                            /* CubitBoolean & */ enum CubitBoolean *spatially_equal,
                                            /* CubitBoolean & */ enum CubitBoolean *tangent_warning,
                                          enum CubitBoolean force_merge);
    /* - calculates the relative sense of two refedges and tries to */
    /* - see if there are two points that are close on each edge (1/3 and 2/3 */
    /* - param value). */

  enum CubitBoolean RefEdge_about_spatially_equal(void *this_ref_edge,
                                                    /* RefEdge * */ void *ref_edge_ptr_2,
                                                  double tolerance_factor,
                                                  enum CubitSense *sensePtr,
                                                  enum CubitBoolean notify_refEntity,
                                                  enum CubitBoolean force_merge);

    /* R enum CubitBoolean */
    /* R-CUBIT_TRUE/CUBIT_FALSE */
    /* I RefEdge, double, CubitSense*, enum CubitBoolean */
    /* I- Second RefEdge to compare, Tolerance factor to for GEOMETRY_RESABS, */
    /* I- and flag for notifying compared RefEntities. */
    /* O CubitSense*, and returned enum CubitBoolean. */
    /* O- if the two refEdges are spatially equal within the GEOMETRY_RESABS* */
    /* - the tolerance_factor, then CUBIT_TRUE will be returned.  Otherwise */
    /* - CUBIT_FALSE is returned.  If there is a non-null sesnePtr sent in, */
    /* - the sense pointer value will be assinged the relative sense of the */
    /* - RefEdges, ie, CUBIT_FORWARD/CUBIT_REVERSED. */
    /* - To do the comparison, the end points, 1/3 and 2/3 world */
    /* - points are spatially compared. */

  int RefEdge_validate(void *this_ref_edge);
    /* - Check that entity is valid. Returns number of problems detected. */

/* =========  Add Code by SRS of Cat,  3/3/99 2:28:31 PM  ========= */
  enum CubitBoolean RefEdge_is_tolerant(void *this_ref_edge);
    /* - Returns CUBIT_TRUE if refedge is a tolerant edge.  Tolerant */
    /* - edges can be created in the ACIS healer if it is unable to */
    /* - heal the edge. */
/* =========  Code End by SRS of Cat,  3/3/99 2:28:31 PM  ========= */

    /* RefVertex * */ void *RefEdge_get_startRefVertex(void *this_ref_edge);
    /* RefVertex * */ void *RefEdge_get_endRefVertex(void *this_ref_edge);
    /* - retrieve the start/end vertices from the chain */

 
#ifdef __cplusplus
}
#endif

#endif
