#ifndef REFFACE_CCAPI
#define REFFACE_CCAPI

#include "CubitDefines.h"
#include "CubitVectorStruct.h"
#include "LocalStart.h"
#include "EntityType.h"

/* class CubitVectorStruct; */

#ifdef __cplusplus
extern "C" {
#endif

enum CubitStatus RefFace_get_co_faces(void *RefFace_ptr,
			    /* DLCoFaceList & */ void ***co_faces_found_list,
                            int *co_faces_found_list_size,
			    /* RefVolume * */ void *input_ref_volume_ptr);
/*R enum CubitStatus */
/*R- CubitSuccess/CubitFailure */
/*O co_faces_found_list */
/*O-Populates the co_faces_found_list with the CoFaces that are */
/*-associated with this RefFace.  Because the input RefVolume is sent */
/*-in it fills the list just with the CoFaces that are associated with */
/*-this RefVolume.  Note that usually there will be only one CoFace with */
/*-the RefVolume but with HardSurfaces there may be more (2). */
  
enum CubitStatus RefFace_ordered_loops(void *RefFace_ptr,
			     /* DLLoopList & */ void ***loop_list,
                             int *loop_list_size);
/*- Gets the loops in order from outside to inside. */
/*- This function is used for all the loop extracting from the refface, */
/*- i.e., nodes, ref-edges, ref-vertices.  This function orders the */
/*- loops based on the angle metric calculation.  I believe this metric */
/*- is actually currently calculated in AcisGeometryEngine but I see */
/*- no reason that it should be done there.  I believe it could be */
/*- done in the Loop class. */

int RefFace_co_edge_loops(void *RefFace_ptr,
		      /* DLCoEdgeLoopList & */ void ***co_edge_loops,
		     int *co_edge_loops_size);
/*- Returns a list of lists.  Each of the included lists contains a list */
/*- of CoEdges and represents an ordered list of CoEdges associated */
/*- with each of the Loops of this RefFace. */
  
int RefFace_ref_edge_loops(void *RefFace_ptr,
		      /* DLRefEdgeLoopList & */ void ***ref_edge_loops,
		      int *ref_edge_loops_size);
/*- Returns a list of lists.  Each of the included lists contains a list */
/*- of RefEdges and represents an ordered list of RefEdges associated */
/*- with each of the Loops of this RefFace. */
/*- NOTE: All of the ref_edge_lists in ref_edge_loops will  */
/*-       need to be deleted by the *calling* function.   */
  
void RefFace_ref_vertex_loops(void *RefFace_ptr,
			 /* DLRefVertLoopList & */ void ***ref_vert_loop_list,
                         int *ref_vert_loop_list_size);
/*- Returns a list of lists.  Each of the included lists contains a list */
/*- of RefVertex'es and represents an ordered list of RefVertex'es  */
/*- associated with each of the Loops of this RefFace.   */
/*- NOTE: All of the ref_vertex_lists in the ref_vert_loop_list will  */
/*-       need to be deleted by the *calling* function. */
  
int RefFace_number_of_Loops(void *RefFace_ptr);
/*- Returns the number of Loops associated with this RefFace */
  
  
/* RefVolume* */ void * RefFace_ref_volume(void *RefFace_ptr); 
/*- Return the first RefVolume pointer to the volume which owns */ 
/*- this RefFace */
/*- Note: There may be more than one RefVolume that owns this RefFace. */
/*-       This method just gets the first in the list. */

void RefFace_hard_points(void *RefFace_ptr,
		    /* DLRefVertexList & */ void ***new_hard_point_list,
                    int *new_hard_point_list_size);
/*- Populate the input DLRefVertexList with the list of hard points  */
/*- that are defined for this RefFace */
  
int RefFace_adjoins(void *RefFace_ptr,
	       /* RefFace * */ void *input_face_ptr);
/*- Returns CUBIT_TRUE if this RefFace adjoins (is connected via a RefEdge) */
/*- the input RefFace */
  
/* RefEdge* */ void *RefFace_common_ref_edge(void *RefFace_ptr,
					/* RefFace * */ void *input_face_ptr);
/*- Returns a common RefEdge* if this RefFace shares one with the */
/*- input RefFace */

enum EntityType RefFace_get_grouping_entity_type(void *RefFace_ptr);
/*R EntityType */
/*R- The type of the corresponding GroupingEntity type. */
/*- This function returns the type of the corresponding */
/*- GroupingEntity . */
  
enum EntityType RefFace_get_sense_entity_type(void *RefFace_ptr);
/*R EntityType */
/*R- The type of SenseEntity associated with this object */
/*- This function returns the type of SenseEntity associated with  */
/*- this BasicTopologyEntity.  */
  
enum EntityType RefFace_get_child_ref_entity_type(void *RefFace_ptr);
/*R EntityType */
/*R- A type value. */
/*- This function returns the type of the child RefEntity of   */
/*- RefFace, which is the type value of the RefEdge class. */
  
enum EntityType RefFace_get_parent_ref_entity_type(void *RefFace_ptr);
/*R EntityType */
/*R- A type value. */
/*- This function returns the type of the parent RefEntity of   */
/*- RefFace, which is the type value of the RefVolume class. */
  
enum EntityType RefFace_get_topology_bridge_type(void *RefFace_ptr);
/*R EntityType */
/*R- The type of GeometryEntity associated with this object */
/*- This function returns the type of GeometryEntity associated with  */
/*- this BasicTopologyEntity.  */

/* CoFace* */ void *RefFace_get_matching_CoFace(void *RefFace_ptr,
					   /* RefVolume * */ void *ref_volume_ptr);
/*R CoFace* */
/*R- Returned CoFace pointer */
/*I ref_volume_ptr */
/*I- The RefVolume to which matching is done. */
/*- This function returns the CoFace that is associated with both "this" */
/*- RefFace as well as with the input RefVolume.  The function is useful */
/*- when a merge operation has resulted in a RefFace that is shared by */
/*- more than 1 RefVolume (most often, two). In this case, the RefFace */
/*- would be associated with more than one CoFace, each belonging to a  */
/*- different RefVolume. */
/*- If there is no match (i.e., this RefFace is not associated with the */
/*- input RefVolume) then a NULL pointer is returned. */
  
/* CubitVector */ struct CubitVectorStruct RefFace_center_point(void *RefFace_ptr);
/*- Return the approximate (spatial) center of this RefFace */
  
enum CubitSense RefFace_sense_1(void *RefFace_ptr,
		    /* RefVolume * */ void *volume);
/*-Determines the sense of "this" with respect to the passed-in volume. */
		
enum CubitSense RefFace_sense_2(void *RefFace_ptr,
		    /* RefFace * */ void *face_ptr);
/*-Determine the relative sense of the passed face with respect to */
/*-this face using the senses of common RefEdges.  i.e. if  */
/*-CUBIT_REVERSED is returned, than the sense of the passed face */
/*-should be the opposite of that of this face with respect to */
/*-any volume.  CUBIT_UNKNOWN is returned if there are no */
/*-common RefEdges between RefFaces or there is more than one */
/*-common RefEdge. */

enum CubitSense RefFace_get_geometry_sense(void *RefFace_ptr);
/*- Gets the sense of the reface with respect to the underlying */
/*- geometry engines representation of the surface. */

enum CubitBoolean RefFace_about_spatially_equal (void *RefFace_ptr,
                                         /* RefFace * */ void *ref_face_ptr_2,
                                       double tolerance_factor,
                                                 enum CubitBoolean notify_refEntity,
                                                 enum CubitBoolean test_bbox,
                                       int test_internal,
                                                 enum CubitBoolean force_merge);
/*R CubitBoolean */
/*R-CUBIT_TRUE/CUBIT_FALSE */
/*I RefFace*, double, CubitBoolean */
/*I- Second RefFace to compare, Tolerance factor to for GEOMETRY_RESABS, */
/*I- and flag for notifying compared RefEntities. */
/*O CubitBoolean */
/*O- If the two RefFaces are spatially equal within the GEOMETRY_RESABS* */
/*- the tolerance_factor, then CUBIT_TRUE will be returned.  Otherwise */
/*- CUBIT_FALSE is returned. */
/*- The comparison is done by first checking the bounding boxes of the */
/*- RefFaces.  If this test is passed then the ref_edges of each face */
/*- are looped through and compared.  A bounding box check for each */
/*- edge is also done first before a comparison, for speed. */
  
enum CubitSense RefFace_compare_alignment(void *RefFace_ptr,
				/* RefFace * */ void *second_ref_face_ptr);
/*R CubitSense */
/*R- Sense of this reface with respect to the second one passed in. */
/*I RefFace * */
/*I- pointer to second ref face with which the alignment is compared. */
/*- This function will compare the sense of the two ref-faces, or */
/*- rather their normals. */
/*- NOTE: It is ASSUMED that BOTH reffaces are SPATIALLY EQUAL. */
/*- If this is not followed this could explode. */
  

/* CubitVector */ struct CubitVectorStruct RefFace_normal_at(void *RefFace_ptr,
                                                               /* const CubitVector& */ const struct CubitVectorStruct location, 
                                                               /* RefVolume * */ void *volume);

/*- Calculate normal for input location (optional input RefVolume to  */
/*- allow for feature consolidation).  */
/*- Note that the input location is modified to the coordinates */
/*- of the closest point on the surface. */
/*- */
/*- MJP NOTE: */
/*- In the previous implementation, the result of this function call  */
/*- would not only be the returned unit vector which is the normal */
/*- at the location, but the function would also fill in the myPosition */
/*- and myParametricPosition data members of RefFace.  These data */
/*- members have been removed in this implementation of RefFace. */
/*- However, the Surface::normal_at function that gets called returns */
/*- an additional parameter which is the location on the underlying */
/*- surface that is closest to the input location. */

void RefFace_reverse_normal(void *RefFace_ptr);
/*- switch the sense of this face with respect to the underlying */
/*- geometry, so that all normals point in the opposite */
/*- direction. The orientation of quads on the surface are switched */
/*- to agree with this normal. */

/*========  Change Code by DZ of Cat,  on 10/29/98 8:46:59 AM  ======== */
enum CubitBoolean RefFace_set_outward_normal(void *RefFace_ptr,
				   /* RefVolume * */ void *volume);
/*- Set the normal of this face to point outward wrt to given volume. */
/*- Assumes there is only one coface of the volume for this RefFace. */
/*- Uses the above "reverse_normal" function. */
/*- return true only if reverse_normal function called. */
/*========  Change End by DZ of Cat,  on 10/29/98 8:46:59 AM  ======== */

void RefFace_move_to_surface (void *RefFace_ptr,
                           /* CubitVector& */ struct CubitVectorStruct *location, 
                         LocalStart start);
/*- Moves the given node back onto its surface */
  
void RefFace_find_closest_point_trimmed(void *RefFace_ptr,
				   /* CubitVector */ struct CubitVectorStruct from_point, 
				   /* CubitVector& */ struct CubitVectorStruct *point_on_surface);
/*R void */
/*I CubitVector */
/*I- point from which to find closest point on trimmed surface */
/*O CubitVector */
/*O- point on trimmed surface closest to passed-in point  */
/*- This function finds the closest point on a TRIMMED surface to the */
/*- passed-in point.   */
  
enum CubitStatus RefFace_get_principal_curvatures(void *RefFace_ptr,
					/* const CubitVector& */ const struct CubitVectorStruct point,
                                        /* double& */ double *curvature1,
					   /* double& */ double *curvature2,
					/* RefVolume * */ void *ref_volume_ptr);
/*R enum CubitStatus  */
/*R- CUBIT_SUCCESS/FAILURE */
/*I point */
/*I- Input location. The coordinates of this input point are  */
/*I- modified to those of the point closest to this one, on the */
/*I- surface of this RefFace. */
/*O curvature1/curvature2 */
/*O- Output principal curvature values. */
/*I ref_volume_ptr */
/*I- Input RefVolume pointer */
/*- This function first computes the point on the surface closest to the */
/*- input point and sets the values of "point" to this closest */
/*- location.  Then, the principal curvatures of the surface at this  */
/*- new point are computed and returned. If the input RefVolume pointer */
/*- is not NULL, it is used when computing the curvatures. */
  
struct CubitVectorStruct RefFace_position_from_u_v(void *RefFace_ptr,
                                       double u, double v);

/*- Return a CubitVector (representing a position vector corresponding  */
/*- to the input point in {u,v} space */
  
enum CubitStatus RefFace_u_v_from_position(void *RefFace_ptr,
				 /* CubitVector const& */ const struct CubitVectorStruct location,
                                 /* double& */ double *u, 
				 /* double& */ double *v,
				 /* CubitVector * */ struct CubitVectorStruct *closest_location);
/*R enum CubitStatus  */
/*R- CUBIT_SUCCESS/FAILURE */
/*I location */
/*I- The input point in global space */
/*O u, v */
/*O- The returned u, v coordinate values (in local parametric space) */
/*O- of the closest_point */
/*O closest_location */
/*O- The point on the Surface closest to the input location */
/*I refvolume_ptr */
/*I- The reference RefVolume with respect to which, the normal */
/*I- is to be computed.  If the pointer is NULL, then the  */
/*I- first underlying solid model entity is used to compute the */
/*I- normal. */
/*- This function returns the {u, v} coordinates of the point  */
/*- on the Surface closest to the input point (specified in global */
/*- space). The closest_location is also returned. */
  
enum CubitBoolean RefFace_is_parametric(void *RefFace_ptr);
/*R CubitBoolean */
/*R- CUBIT_TRUE/CUBIT_FALSE */
/*- This function determines whether the underlying geometry of the */
/*- Surface is parametrically defined or not.  Returns CUBIT_TRUE if  */
/*- it is and CUBIT_FALSE if it is not. */
  
enum CubitBoolean RefFace_get_param_range_U(void *RefFace_ptr,
                                  /* double& */ double *lower_bound,
				  /* double& */ double *upper_bound);
/*R CubitBoolean */
/*R- CUBIT_TRUE/CUBIT_FALSE */
/*O lower_bound */
/*O- The lower bound of the parametric range in the U direction. */
/*O- This is set to 0.0 if the surface is not parametric. */
/*O upper_bound */
/*O- The upper bound of the parametric range in the U direction. */
/*O- This is set to 0.0 if the surface is not parametric. */
/*- Returns the lower and upper parametric bounds of the  */
/*- surface in U, if it is parametric.  Otherwise, it returns */
/*- CUBIT_FALSE and zeroes for the upper and lower parametric */
/*- bounds. */
  
enum CubitBoolean RefFace_get_param_range_V(void *RefFace_ptr,
                                  /* double& */ double *lower_bound,
				  /* double& */ double *upper_bound);
/*R CubitBoolean */
/*R- CUBIT_TRUE/CUBIT_FALSE */
/*O lower_bound */
/*O- The lower bound of the parametric range in the V direction. */
/*O- This is set to 0.0 if the surface is not parametric. */
/*O upper_bound */
/*O- The upper bound of the parametric range in the V direction. */
/*O- This is set to 0.0 if the surface is not parametric. */
/*- Returns the lower and upper parametric bounds of the  */
/*- surface in V, if it is parametric.  Otherwise, it returns */
/*- CUBIT_FALSE and zeroes for the upper and lower parametric */
/*- bounds. */
  
enum CubitBoolean RefFace_is_periodic(void *RefFace_ptr);
/*R CubitBoolean */
/*R- CUBIT_TRUE/CUBIT_FALSE */
/*- This function determines whether the underlying geometry of the */
/*- Surface is periodic or not.  Returns CUBIT_TRUE if it is and  */
/*- CUBIT_FALSE if it is not. */
  
enum CubitBoolean RefFace_is_periodic_in_U(void *RefFace_ptr,
                                 /* double& */ double *period);
/*R CubitBoolean */
/*R- CUBIT_TRUE/CUBIT_FALSE */
/*O period */
/*O- The value of the period in the U direction. */
/*- Determines whether the surface object is  */
/*- periodic in the U direction or not.  If it is, it */
/*- returns CUBIT_TRUE and the value of the period. Otherwise, */
/*- it returns CUBIT_FALSE and a value of 0.0 or the period. */
  
enum CubitBoolean RefFace_is_periodic_in_V(void *RefFace_ptr,
                                 /* double& */ double *period);
/*R CubitBoolean */
/*R- CUBIT_TRUE/CUBIT_FALSE */
/*O period */
/*O- The value of the period in the V direction. */
/*- Determines whether the surface object is  */
/*- periodic in the V direction or not.  If it is, it */
/*- returns CUBIT_TRUE and the value of the period. Otherwise, */
/*- it returns CUBIT_FALSE and a value of 0.0 or the period. */

enum CubitBoolean RefFace_is_singular_in_U(void *RefFace_ptr,
                                 double u_param);

enum CubitBoolean RefFace_is_singular_in_V(void *RefFace_ptr,
                                 double v_param);
/*R CubitBoolean */
/*R- CUBIT_TRUE/CUBIT_FALSE */
/*I double u/v parameter value. */
/*- Determines if the surface is singular in a given direction */
/*- at a given parameter value. */

enum CubitBoolean RefFace_is_closed_in_U(void *RefFace_ptr);

enum CubitBoolean RefFace_is_closed_in_V(void *RefFace_ptr);
/*R CubitBoolean */
/*R- CUBIT_TRUE/CUBIT_FALSE */
/*- Determines if the surface is closed, smoothly or not in the */
/*- given parameter direction. */
/*- A periodic surface is always closed but a closed surface is */
/*- is not always periodic. */

enum CubitStatus RefFace_uv_derivitives(void *RefFace_ptr,
                              double u_param,
                              double v_param,
                              /* CubitVector & */ struct CubitVectorStruct *du,
                              /* CubitVector & */ struct CubitVectorStruct *dv);
/*R enum CubitStatus */
/*R- CUBIT_SUCCESS/CUBIT_FAILURE */
/*O- du, dv */
/*- Determines the u and v derivitives from the given parameter */
/*- values. */
   
int RefFace_dimension(void *RefFace_ptr);
/*- Returns the geometric dimension of RefFace entities.  */

double RefFace_area(void *RefFace_ptr);
/*- get the area of the underlying surface */

double RefFace_measure(void *RefFace_ptr);

/* CubitString */ const char *RefFace_measure_label(void *RefFace_ptr);

enum CubitBoolean RefFace_is_planar(void *RefFace_ptr);
/*R CubitBoolean */
/*R CUBIT_TRUE/CUBIT_FALSE */
/*- This function returns CUBIT_TRUE if the underlying geometry */
/*- of the face is planar. CUBIT_FALSE otherwise. */

enum CubitStatus RefFace_get_point_normal(void *RefFace_ptr,
				/* CubitVector& */ struct CubitVectorStruct *origin, 
				/* CubitVector& */ struct CubitVectorStruct *normal);
/*- Only valid for planar surfaces */
/*- Finds the underlying plane's origin and normal (unit) vector */
/*- Returns CUBIT_FAILURE if surface is not a plane */
  
int RefFace_validate(void *RefFace_ptr);
/*- Check that entity is valid. Returns number of problems detected. */

double RefFace_get_crack_length(void *RefFace_ptr);
/*- return the length of the periodic crack, or 0.0 if non-periodic. */

void RefFace_add_hard_point_1(void *RefFace_ptr,
                         double x, double y, double z);
/*- Add a hard point to this RefFace.  The input "point" is first */
/*- moved to the underlying surface, and then a hard point is added */
/*- at this, possibly different, location on the surface. */
/*-*************************************************************** */
/*- MJP Note: */
/*- Currently, the new RefVertex that is created is NOT PART OF THE */
/*- main DAG datastructure. The new RefVertex, however, has its */
/*- own little mini-DAG which consists of a single CDODAGNode. */
/*- Discuss this with the team before making it a part of the main DAG.  */
/*- It is, however, deleted appropriately when the RefFace is deleted. */
/*-*************************************************************** */
      
void RefFace_add_hard_point_2(void *RefFace_ptr,
                           /* RefVertex * */ void *ref_vertex_ptr);
/*- Add a hard point to this RefFace.  */
/*-*************************************************************** */
/*- MJP Note: */
/*- Currently, the new RefVertex that is created is NOT PART OF THE */
/*- main DAG datastructure. The new RefVertex, however, has its */
/*- own little mini-DAG which consists of a single CDODAGNode. */
/*- Discuss this with the team before making it a part of the main DAG.  */
/*- It is, however, deleted appropriately when the RefFace is deleted. */
/*-*************************************************************** */
    
/***********Graphics Related Functions********* */
void RefFace_draw_my_edges (void *RefFace_ptr,
                       int color);
  
void RefFace_draw_normal(void *RefFace_ptr,
		    int color,
		    double length,
		    /* RefVolume * */ void *wrt_vol);

/* Surface* */ void *RefFace_get_surface_ptr(void *RefFace_ptr);
/*R Surface* */
/*R- A pointer to the Surface to which the current  */
/*R- face points.  */
/*- This function returns a pointer to the Surface */
/*- to which the current face points. */
  
enum EntityType RefFace_entity_type(void *RefFace_ptr);
  
void* RefFace_get_address(void *RefFace_ptr,
                     enum EntityType inputEntityType);
/*R void* */
/*R- Returned void pointer */
/*I inputEntityType */
/*I- The input type to get the address of. */
/*- This function returns a void pointer that points to the */
/*- "appropriate" portion of this object.  The appropriate */
/*- portion is determined by the input EntityType variable. */
/*- Returns NULL if the input type and the type of "this" */
/*- are not related by inheritance. */
/*- Note: The RTTI capabilities encoded in these functions */
/*-       are designed to work with any form of multiple  */
/*-       inheritance, as well.  Multiple inheritance is what */
/*-       necessitates having this function defined in every  */
/*-       class in the hierarchy. */
/*- Note: This function can also be used to merely check if */
/*-       an object of one type is related to another type */
/*-       through inheritance.  If a non-NULL pointer is */
/*-       returned, then this is true. */
  
  void *RefFace_TopologyEntity(void *RefFace_ptr);
    /* cast to a TopologyEntity, then return a void * */

  void *RefFace_RefEntity(void *RefFace_ptr);
    /* cast to a RefEntity, then return a void * */

  enum CubitStatus RefFace_setup_use_facets( void *RefFace_ptr,
                                               /* DLCubitFacetList &*/ void ***facet_list, int *facet_list_size,
                                             int interp_order);
    /* set up the underlying geometry to use a faceted base representation. */

  void RefFace_stop_use_facets(void *RefFace_ptr);
    /* stop using the underlying geometry to use a faceted base representation. */

  enum CubitBoolean RefFace_get_use_facets(void *RefFace_ptr);
    /* find out if the surface is using the faceted rep. */

#ifdef __cplusplus
}
#endif

#endif
