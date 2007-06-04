/**
 * \file RefFace_ccapi.cpp
 *
 * \brief C api to classes in RefFace class
 *
 * \author Tim Tautges
 *
 * \date 7/2000
 *
 */

#include "RefFace_ccapi.h"
#include "RefFace.hpp"
#include "CubitDefines.h"
#include "DLLoopList.hpp"
#include "DLCoEdgeLoopList.hpp"
#include "DLRefEdgeLoopList.hpp"
#include "DLRefVertLoopList.hpp"
#include "DLCoFaceList.hpp"
#include "DLRefVertexList.hpp"
#include "DLRefEdgeList.hpp"
#include "DLRefVolumeList.hpp"
#include "DLCubitFacetList.hpp"
#include "copy_defines.h"

    /* topology */

  CubitStatus RefFace_get_co_faces(void *this_ref_face,  /* DLCoFaceList& */ void ***co_faces_found_list, int *co_faces_found_list_size,
                                     /* RefVolume * */ void *input_ref_volume_ptr )
{
  RefFace *temp_ref_face = (RefFace *) this_ref_face;
  RefVolume *temp_input_ref_volume_ptr = (RefVolume *) input_ref_volume_ptr;

  DLCoFaceList temp_list;

  CubitStatus status = temp_ref_face->get_co_faces(temp_list, temp_input_ref_volume_ptr);
  COPY_LIST_TO_ARRAY(temp_list, *co_faces_found_list, *co_faces_found_list_size);

  return status;  
}
    //R CubitStatus
    //R- CubitSuccess/CubitFailure
    //O co_faces_found_list
    //O-Populates the co_faces_found_list with the CoFaces that are
    //-associated with this RefFace.  Because the input RefVolume is sent
    //-in it fills the list just with the CoFaces that are associated with
    //-this RefVolume.  Note that usually there will be only one CoFace with
    //-the RefVolume but with HardSurfaces there may be more (2).

  CubitStatus RefFace_ordered_loops(void *this_ref_face,  /* DLLoopList& */ void ***loop_list, int *loop_list_size )
{
  RefFace *temp_ref_face = (RefFace *) this_ref_face;

  DLLoopList temp_list;

  CubitStatus status = temp_ref_face->ordered_loops(temp_list);
  COPY_LIST_TO_ARRAY(temp_list, *loop_list, *loop_list_size);

  return status;  
}
    //- Gets the loops in order from outside to inside.
    //- This function is used for all the loop extracting from the refface,
    //- i.e., nodes, ref-edges, ref-vertices.  This function orders the
    //- loops based on the angle metric calculation.  I believe this metric
    //- is actually currently calculated in AcisGeometryEngine but I see
    //- no reason that it should be done there.  I believe it could be
    //- done in the Loop class.

  int RefFace_co_edge_loops (void *this_ref_face,  /* DLCoEdgeLoopList& */ void ***co_edge_loops, int *co_edge_loops_size )
{
  RefFace *temp_ref_face = (RefFace *) this_ref_face;

  DLCoEdgeLoopList temp_list;

  int status = temp_ref_face->co_edge_loops(temp_list);
  COPY_LIST_TO_ARRAY(temp_list, *co_edge_loops, *co_edge_loops_size);

  return status;  
}
    //- Returns a list of lists.  Each of the included lists contains a list
    //- of CoEdges and represents an ordered list of CoEdges associated
    //- with each of the Loops of this RefFace.
  
  int RefFace_ref_edge_loops (void *this_ref_face,  /* DLRefEdgeLoopList& */ void ***ref_edge_loops, int *ref_edge_loops_size )
{
  RefFace *temp_ref_face = (RefFace *) this_ref_face;

  DLRefEdgeLoopList temp_list;

  int status = temp_ref_face->ref_edge_loops(temp_list);
  COPY_LIST_TO_ARRAY(temp_list, *ref_edge_loops, *ref_edge_loops_size);

  return status;  
}
    //- Returns a list of lists.  Each of the included lists contains a list
    //- of RefEdges and represents an ordered list of RefEdges associated
    //- with each of the Loops of this RefFace.
    //- NOTE: All of the ref_edge_lists in ref_edge_loops will 
    //-       need to be deleted by the *calling* function.  
  
  void RefFace_ref_vertex_loops(void *this_ref_face,  /* DLRefVertLoopList& */ void ***ref_vert_loop_list, int *ref_vert_loop_list_size)
{
  RefFace *temp_ref_face = (RefFace *) this_ref_face;

  DLRefVertLoopList temp_list;

  temp_ref_face->ref_vertex_loops(temp_list);
  COPY_LIST_TO_ARRAY(temp_list, *ref_vert_loop_list, *ref_vert_loop_list_size);
}
    //- Returns a list of lists.  Each of the included lists contains a list
    //- of RefVertex'es and represents an ordered list of RefVertex'es 
    //- associated with each of the Loops of this RefFace.  
    //- NOTE: All of the ref_vertex_lists in the ref_vert_loop_list will 
    //-       need to be deleted by the *calling* function.
  
  int RefFace_number_of_Loops (void *this_ref_face)
{
  RefFace *temp_ref_face = (RefFace *) this_ref_face;

  return temp_ref_face->number_of_Loops();
}
    //- Returns the number of Loops associated with this RefFace
  
  
    /* RefVolume* */ void *RefFace_ref_volume(void *this_ref_face)
{
  RefFace *temp_ref_face = (RefFace *) this_ref_face;

  return temp_ref_face->ref_volume();
}
    //- Return the first RefVolume pointer to the volume which owns
    //- this RefFace
    //- Note: There may be more than one RefVolume that owns this RefFace.
    //-       This method just gets the first in the list.

  void RefFace_hard_points(void *this_ref_face,  /* DLRefVertexList& */ void ***new_hard_point_list, int *new_hard_point_list_size )
{
  RefFace *temp_ref_face = (RefFace *) this_ref_face;

  DLRefVertexList temp_list;
  
  temp_ref_face->hard_points(temp_list);
  COPY_LIST_TO_ARRAY(temp_list, *new_hard_point_list, *new_hard_point_list_size);
}
    //- Populate the input DLRefVertexList with the list of hard points 
    //- that are defined for this RefFace
  
  int RefFace_adjoins (void *this_ref_face,  /* RefFace* */ void *input_face_ptr )
{
  RefFace *temp_ref_face = (RefFace *) this_ref_face;

  RefFace *temp_input_face_ptr = (RefFace *) input_face_ptr;
  
  return temp_ref_face->adjoins(temp_input_face_ptr);
}
    //- Returns CUBIT_TRUE if this RefFace adjoins (is connected via a RefEdge)
    //- the input RefFace
  
/* RefEdge* */ void *RefFace_common_ref_edge (void *this_ref_face,  /* RefFace* */ void *input_face_ptr )
{
  RefFace *temp_ref_face = (RefFace *) this_ref_face;

  RefFace *temp_input_face_ptr = (RefFace *) input_face_ptr;
  
  return temp_ref_face->common_ref_edge(temp_input_face_ptr);
}
    //- Returns a common RefEdge* if this RefFace shares one with the
    //- input RefFace

  EntityType RefFace_get_grouping_entity_type(void *) 
{ return Loop_TYPE;}
    //R EntityType
    //R- The type of the corresponding GroupingEntity type.
    //- This function returns the type of the corresponding
    //- GroupingEntity .
  
  EntityType RefFace_get_sense_entity_type(void *) 
{ return CoFace_TYPE;}
    //R EntityType
    //R- The type of SenseEntity associated with this object
    //- This function returns the type of SenseEntity associated with 
    //- this BasicTopologyEntity. 
  
  EntityType RefFace_get_child_ref_entity_type(void *) 
{ return RefEdge_TYPE;}
    //R EntityType
    //R- A type value.
    //- This function returns the type of the child RefEntity of  
    //- RefFace, which is the type value of the RefEdge class.
  
  EntityType RefFace_get_parent_ref_entity_type(void *) 
{ return RefVolume_TYPE;}
    //R EntityType
    //R- A type value.
    //- This function returns the type of the parent RefEntity of  
    //- RefFace, which is the type value of the RefVolume class.
  
  EntityType RefFace_get_topology_bridge_type(void *) 
{ return Surface_TYPE;}
    //R EntityType
    //R- The type of GeometryEntity associated with this object
    //- This function returns the type of GeometryEntity associated with 
    //- this BasicTopologyEntity. 

/* CoFace* */ void *RefFace_get_matching_CoFace(void *this_ref_face, /* RefVolume* */ void *ref_volume_ptr)
{
  RefFace *temp_ref_face = (RefFace *) this_ref_face;

  RefVolume *temp_ref_volume_ptr = (RefVolume *) ref_volume_ptr;
  
  return temp_ref_face->get_matching_CoFace(temp_ref_volume_ptr);
}
    //R CoFace*
    //R- Returned CoFace pointer
    //I ref_volume_ptr
    //I- The RefVolume to which matching is done.
    //- This function returns the CoFace that is associated with both "this"
    //- RefFace as well as with the input RefVolume.  The function is useful
    //- when a merge operation has resulted in a RefFace that is shared by
    //- more than 1 RefVolume (most often, two). In this case, the RefFace
    //- would be associated with more than one CoFace, each belonging to a 
    //- different RefVolume.
    //- If there is no match (i.e., this RefFace is not associated with the
    //- input RefVolume) then a NULL pointer is returned.
  
    /*  geometry */

  CubitVectorStruct RefFace_center_point(void *this_ref_face)
{
  RefFace *temp_ref_face = (RefFace *) this_ref_face;

  return temp_ref_face->center_point();
}
    //- Return the approximate (spatial) center of this RefFace
  
  CubitSense RefFace_sense_1(void *this_ref_face, /* RefVolume* */ void *volume)
{
  RefFace *temp_ref_face = (RefFace *) this_ref_face;

  RefVolume *temp_volume = (RefVolume *) volume;
  
  return temp_ref_face->sense(temp_volume);
}
    //-Determines the sense of "this" with respect to the passed-in volume.
		
  CubitSense RefFace_sense_2(void *this_ref_face, /* RefFace* */ void *face_ptr)
{
  RefFace *temp_ref_face = (RefFace *) this_ref_face;

  RefFace *temp_face_ptr = (RefFace *) face_ptr;
  
  return temp_ref_face->sense(temp_face_ptr);
}
    //-Determine the relative sense of the passed face with respect to
    //-this face using the senses of common RefEdges.  i.e. if 
    //-CUBIT_REVERSED is returned, than the sense of the passed face
    //-should be the opposite of that of this face with respect to
    //-any volume.  CUBIT_UNKNOWN is returned if there are no
    //-common RefEdges between RefFaces or there is more than one
    //-common RefEdge.

  CubitSense RefFace_get_geometry_sense(void *this_ref_face)
{
  RefFace *temp_ref_face = (RefFace *) this_ref_face;

  return temp_ref_face->get_geometry_sense();
}
    //- Gets the sense of the reface with respect to the underlying
    //- geometry engines representation of the surface.

  CubitBoolean RefFace_about_spatially_equal (void *this_ref_face,  /* RefFace* */ void *ref_face_ptr_2,
                                              double tolerance_factor,
                                              CubitBoolean notify_refEntity,
                                              CubitBoolean test_bbox,
                                              int test_internal,
                                              CubitBoolean force_merge)
{
  RefFace *temp_ref_face = (RefFace *) this_ref_face;

  RefFace *temp_ref_face_ptr_2 = (RefFace *) ref_face_ptr_2;
  
  return temp_ref_face->about_spatially_equal(temp_ref_face_ptr_2, tolerance_factor, notify_refEntity,
                                              test_bbox, test_internal, force_merge);
}
  
    //R CubitBoolean
    //R-CUBIT_TRUE/CUBIT_FALSE
    //I RefFace*, double, CubitBoolean
    //I- Second RefFace to compare, Tolerance factor to for GEOMETRY_RESABS,
    //I- and flag for notifying compared RefEntities.
    //O CubitBoolean
    //O- If the two RefFaces are spatially equal within the GEOMETRY_RESABS*
    //- the tolerance_factor, then CUBIT_TRUE will be returned.  Otherwise
    //- CUBIT_FALSE is returned.
    //- The comparison is done by first checking the bounding boxes of the
    //- RefFaces.  If this test is passed then the ref_edges of each face
    //- are looped through and compared.  A bounding box check for each
    //- edge is also done first before a comparison, for speed.
  
  CubitSense RefFace_compare_alignment(void *this_ref_face,  /* RefFace* */ void *second_ref_face_ptr )
{
  RefFace *temp_ref_face = (RefFace *) this_ref_face;

  RefFace *temp_second_ref_face_ptr = (RefFace *) second_ref_face_ptr;
  
  return temp_ref_face->compare_alignment(temp_second_ref_face_ptr);
}
    //R CubitSense
    //R- Sense of this reface with respect to the second one passed in.
    //I RefFace *
    //I- pointer to second ref face with which the alignment is compared.
    //- This function will compare the sense of the two ref-faces, or
    //- rather their normals.
    //- NOTE: It is ASSUMED that BOTH reffaces are SPATIALLY EQUAL.
    //- If this is not followed this could explode.
  

  CubitVectorStruct RefFace_normal_at(void *this_ref_face, const CubitVectorStruct location, /* RefVolume* */ void *volume)
{
  RefFace *temp_ref_face = (RefFace *) this_ref_face;

  RefVolume *temp_volume = (RefVolume *) volume;
  
  return temp_ref_face->normal_at(location, temp_volume);
}
    //- Calculate normal for input location (optional input RefVolume to 
    //- allow for feature consolidation). 
    //- Note that the input location is modified to the coordinates
    //- of the closest point on the surface.
    //-
    //- MJP NOTE:
    //- In the previous implementation, the result of this function call 
    //- would not only be the returned unit vector which is the normal
    //- at the location, but the function would also fill in the myPosition
    //- and myParametricPosition data members of RefFace.  These data
    //- members have been removed in this implementation of RefFace.
    //- However, the Surface::normal_at function that gets called returns
    //- an additional parameter which is the location on the underlying
    //- surface that is closest to the input location.
  void RefFace_reverse_normal(void *this_ref_face)
{
  RefFace *temp_ref_face = (RefFace *) this_ref_face;

  temp_ref_face->reverse_normal();
}
    //- switch the sense of this face with respect to the underlying
    //- geometry, so that all normals point in the opposite
    //- direction. The orientation of quads on the surface are switched
    //- to agree with this normal.

//========  Change Code by DZ of Cat,  on 10/29/98 8:46:59 AM  ========
  CubitBoolean RefFace_set_outward_normal(void *this_ref_face, /* RefVolume * */ void *volume )
{
  RefFace *temp_ref_face = (RefFace *) this_ref_face;

  RefVolume *temp_volume = (RefVolume *) volume;
  
  return temp_ref_face->set_outward_normal(temp_volume);
}
    //- Set the normal of this face to point outward wrt to given volume.
    //- Assumes there is only one coface of the volume for this RefFace.
    //- Uses the above "reverse_normal" function.
    //- return true only if reverse_normal function called.
//========  Change End by DZ of Cat,  on 10/29/98 8:46:59 AM  ========

  void RefFace_move_to_surface (void *this_ref_face,  /* CubitVector& */ CubitVectorStruct *location, LocalStart start)
{
  RefFace *temp_ref_face = (RefFace *) this_ref_face;

  CubitVector temp_location = *location;
  
  temp_ref_face->move_to_surface(temp_location, start);

  *location = temp_location;
}
    //- Moves the given node back onto its surface
  
  void RefFace_find_closest_point_trimmed(void *this_ref_face, CubitVectorStruct from_point, 
                                            /* CubitVector& */ CubitVectorStruct *point_on_surface)
{
  RefFace *temp_ref_face = (RefFace *) this_ref_face;

  CubitVector temp_point_on_surface = *point_on_surface;
  
  temp_ref_face->find_closest_point_trimmed(from_point, temp_point_on_surface);

  *point_on_surface = temp_point_on_surface;
}
    //R void
    //I CubitVector
    //I- point from which to find closest point on trimmed surface
    //O CubitVector
    //O- point on trimmed surface closest to passed-in point 
    //- This function finds the closest point on a TRIMMED surface to the
    //- passed-in point.  
  
  CubitStatus RefFace_get_principal_curvatures(void *this_ref_face,  const CubitVectorStruct point,
                                                 /* double& */ double *curvature1,
                                                 /* double& */ double *curvature2,
                                                 /* RefVolume* */ void *ref_volume_ptr)
{
  RefFace *temp_ref_face = (RefFace *) this_ref_face;

  CubitVector temp_point = point;

  RefVolume *temp_ref_volume_ptr = (RefVolume *) ref_volume_ptr;
  
  return temp_ref_face->get_principal_curvatures(temp_point, *curvature1, *curvature2, temp_ref_volume_ptr);
}
    //R CubitStatus 
    //R- CUBIT_SUCCESS/FAILURE
    //I point
    //I- Input location. The coordinates of this input point are 
    //I- modified to those of the point closest to this one, on the
    //I- surface of this RefFace.
    //O curvature1/curvature2
    //O- Output principal curvature values.
    //I ref_volume_ptr
    //I- Input RefVolume pointer
    //- This function first computes the point on the surface closest to the
    //- input point and sets the values of "point" to this closest
    //- location.  Then, the principal curvatures of the surface at this 
    //- new point are computed and returned. If the input RefVolume pointer
    //- is not NULL, it is used when computing the curvatures.
  
  CubitVectorStruct RefFace_position_from_u_v (void *this_ref_face, double u, double v)
{
  RefFace *temp_ref_face = (RefFace *) this_ref_face;

  return temp_ref_face->position_from_u_v(u, v);
}
    //- Return a CubitVector (representing a position vector corresponding 
    //- to the input point in {u,v} space
  
  CubitStatus RefFace_u_v_from_position (void *this_ref_face, const CubitVectorStruct location,
                                           /* double& */ double *u, 
                                           /* double& */ double *v,
                                           /* CubitVector* */ CubitVectorStruct *closest_location)
{
  RefFace *temp_ref_face = (RefFace *) this_ref_face;

  CubitVector temp_location(location);

  CubitVector *temp_closest_location =
    (closest_location ? new CubitVector(*closest_location) : NULL);
    
  CubitStatus status = temp_ref_face->u_v_from_position(temp_location, *u, *v, temp_closest_location);

  if (temp_closest_location != NULL) {
    *closest_location = *temp_closest_location;
    delete temp_closest_location;
  }

  return status;
}
    //R CubitStatus 
    //R- CUBIT_SUCCESS/FAILURE
    //I location
    //I- The input point in global space
    //O u, v
    //O- The returned u, v coordinate values (in local parametric space)
    //O- of the closest_point
    //O closest_location
    //O- The point on the Surface closest to the input location
    //I refvolume_ptr
    //I- The reference RefVolume with respect to which, the normal
    //I- is to be computed.  If the pointer is NULL, then the 
    //I- first underlying solid model entity is used to compute the
    //I- normal.
    //- This function returns the {u, v} coordinates of the point 
    //- on the Surface closest to the input point (specified in global
    //- space). The closest_location is also returned.
  
  CubitBoolean RefFace_is_parametric(void *this_ref_face)
{
  RefFace *temp_ref_face = (RefFace *) this_ref_face;

  return temp_ref_face->is_parametric();
}
    //R CubitBoolean
    //R- CUBIT_TRUE/CUBIT_FALSE
    //- This function determines whether the underlying geometry of the
    //- Surface is parametrically defined or not.  Returns CUBIT_TRUE if 
    //- it is and CUBIT_FALSE if it is not.
  
  CubitBoolean RefFace_get_param_range_U(void *this_ref_face,  /* double& */ double *lower_bound,
                                           /* double& */ double *upper_bound )
{
  RefFace *temp_ref_face = (RefFace *) this_ref_face;

  return temp_ref_face->get_param_range_U(*lower_bound, *upper_bound);
}
    //R CubitBoolean
    //R- CUBIT_TRUE/CUBIT_FALSE
    //O lower_bound
    //O- The lower bound of the parametric range in the U direction.
    //O- This is set to 0.0 if the surface is not parametric.
    //O upper_bound
    //O- The upper bound of the parametric range in the U direction.
    //O- This is set to 0.0 if the surface is not parametric.
    //- Returns the lower and upper parametric bounds of the 
    //- surface in U, if it is parametric.  Otherwise, it returns
    //- CUBIT_FALSE and zeroes for the upper and lower parametric
    //- bounds.
  
  CubitBoolean RefFace_get_param_range_V(void *this_ref_face,  /* double& */ double *lower_bound,
                                           /* double& */ double *upper_bound )
{
  RefFace *temp_ref_face = (RefFace *) this_ref_face;

  return temp_ref_face->get_param_range_V(*lower_bound, *upper_bound);
}
    //R CubitBoolean
    //R- CUBIT_TRUE/CUBIT_FALSE
    //O lower_bound
    //O- The lower bound of the parametric range in the V direction.
    //O- This is set to 0.0 if the surface is not parametric.
    //O upper_bound
    //O- The upper bound of the parametric range in the V direction.
    //O- This is set to 0.0 if the surface is not parametric.
    //- Returns the lower and upper parametric bounds of the 
    //- surface in V, if it is parametric.  Otherwise, it returns
    //- CUBIT_FALSE and zeroes for the upper and lower parametric
    //- bounds.
  
  CubitBoolean RefFace_is_periodic(void *this_ref_face)
{
  RefFace *temp_ref_face = (RefFace *) this_ref_face;

  return temp_ref_face->is_periodic();
}
    //R CubitBoolean
    //R- CUBIT_TRUE/CUBIT_FALSE
    //- This function determines whether the underlying geometry of the
    //- Surface is periodic or not.  Returns CUBIT_TRUE if it is and 
    //- CUBIT_FALSE if it is not.
  
  CubitBoolean RefFace_is_periodic_in_U(void *this_ref_face,  /* double& */ double *period )
{
  RefFace *temp_ref_face = (RefFace *) this_ref_face;

  return temp_ref_face->is_periodic_in_U(*period);
}
    //R CubitBoolean
    //R- CUBIT_TRUE/CUBIT_FALSE
    //O period
    //O- The value of the period in the U direction.
    //- Determines whether the surface object is 
    //- periodic in the U direction or not.  If it is, it
    //- returns CUBIT_TRUE and the value of the period. Otherwise,
    //- it returns CUBIT_FALSE and a value of 0.0 or the period.
  
  CubitBoolean RefFace_is_periodic_in_V(void *this_ref_face,  /* double& */ double *period )
{
  RefFace *temp_ref_face = (RefFace *) this_ref_face;

  return temp_ref_face->is_periodic_in_V(*period);
}
    //R CubitBoolean
    //R- CUBIT_TRUE/CUBIT_FALSE
    //O period
    //O- The value of the period in the V direction.
    //- Determines whether the surface object is 
    //- periodic in the V direction or not.  If it is, it
    //- returns CUBIT_TRUE and the value of the period. Otherwise,
    //- it returns CUBIT_FALSE and a value of 0.0 or the period.

  CubitBoolean RefFace_is_singular_in_U(void *this_ref_face,  double u_param )
{
  RefFace *temp_ref_face = (RefFace *) this_ref_face;

  return temp_ref_face->is_singular_in_U(u_param);
}
  CubitBoolean RefFace_is_singular_in_V(void *this_ref_face,  double v_param )
{
  RefFace *temp_ref_face = (RefFace *) this_ref_face;

  return temp_ref_face->is_singular_in_V(v_param);
}
    //R CubitBoolean
    //R- CUBIT_TRUE/CUBIT_FALSE
    //I double u/v parameter value.
    //- Determines if the surface is singular in a given direction
    //- at a given parameter value.

  CubitBoolean RefFace_is_closed_in_U(void *this_ref_face)
{
  RefFace *temp_ref_face = (RefFace *) this_ref_face;

  return temp_ref_face->is_closed_in_U();
}
  CubitBoolean RefFace_is_closed_in_V(void *this_ref_face)
{
  RefFace *temp_ref_face = (RefFace *) this_ref_face;

  return temp_ref_face->is_closed_in_V();
}
    //R CubitBoolean
    //R- CUBIT_TRUE/CUBIT_FALSE
    //- Determines if the surface is closed, smoothly or not in the
    //- given parameter direction.
    //- A periodic surface is always closed but a closed surface is
    //- is not always periodic.

  CubitStatus RefFace_uv_derivitives(void *this_ref_face,  double u_param,
                                     double v_param,
                                       /* CubitVector & */ CubitVectorStruct *du,
                                       /* CubitVector & */ CubitVectorStruct *dv )
{
  RefFace *temp_ref_face = (RefFace *) this_ref_face;

  CubitVector temp_du, temp_dv;
  
  CubitStatus status = temp_ref_face->uv_derivitives(u_param, v_param, temp_du, temp_dv);

  *du = temp_du;
  *dv = temp_dv;

  return status;
}
    //R CubitStatus
    //R- CUBIT_SUCCESS/CUBIT_FAILURE
    //O- du, dv
    //- Determines the u and v derivitives from the given parameter
    //- values.

  int RefFace_dimension(void *this_ref_face)
{
  RefFace *temp_ref_face = (RefFace *) this_ref_face;

  return temp_ref_face->dimension();
}
    //- Returns the geometric dimension of RefFace entities. 

  double RefFace_area(void *this_ref_face)
{
  RefFace *temp_ref_face = (RefFace *) this_ref_face;

  return temp_ref_face->area();
}
    //- get the area of the underlying surface

  double RefFace_measure(void *this_ref_face)
{
  RefFace *temp_ref_face = (RefFace *) this_ref_face;

  return temp_ref_face->measure();
}
    /* CubitString */ const char *RefFace_measure_label(void *this_ref_face)
{
  RefFace *temp_ref_face = (RefFace *) this_ref_face;

  return temp_ref_face->measure_label().c_str();
}
  
  CubitBoolean RefFace_is_planar(void *this_ref_face)
{
  RefFace *temp_ref_face = (RefFace *) this_ref_face;

  return temp_ref_face->is_planar();
}
    //R CubitBoolean
    //R CUBIT_TRUE/CUBIT_FALSE
    //- This function returns CUBIT_TRUE if the underlying geometry
    //- of the face is planar. CUBIT_FALSE otherwise.

  CubitStatus RefFace_get_point_normal(void *this_ref_face,  /* CubitVector& */ CubitVectorStruct *origin, /* CubitVector& */ CubitVectorStruct *normal )
{
  RefFace *temp_ref_face = (RefFace *) this_ref_face;

  CubitVector temp_origin, temp_normal;
  
  CubitStatus status = temp_ref_face->get_point_normal(temp_origin, temp_normal);

  *origin = temp_origin;

  *normal = temp_normal;

  return status;
}
    //- Only valid for planar surfaces
    //- Finds the underlying plane's origin and normal (unit) vector
    //- Returns CUBIT_FAILURE if surface is not a plane
  
  int RefFace_validate(void *this_ref_face)
{
  RefFace *temp_ref_face = (RefFace *) this_ref_face;

  return temp_ref_face->validate();
}
    //- Check that entity is valid. Returns number of problems detected.

  double RefFace_get_crack_length(void *this_ref_face)
{
  RefFace *temp_ref_face = (RefFace *) this_ref_face;

  return temp_ref_face->get_crack_length();
}
    //- return the length of the periodic crack, or 0.0 if non-periodic.

    /* geometry modification */

  void RefFace_add_hard_point(void *this_ref_face,  double x, double y, double z )
{
  RefFace *temp_ref_face = (RefFace *) this_ref_face;

  temp_ref_face->add_hard_point(x, y, z);
}
    //- Add a hard point to this RefFace.  The input "point" is first
    //- moved to the underlying surface, and then a hard point is added
    //- at this, possibly different, location on the surface.
    //-***************************************************************
    //- MJP Note:
    //- Currently, the new RefVertex that is created is NOT PART OF THE
    //- main DAG datastructure. The new RefVertex, however, has its
    //- own little mini-DAG which consists of a single CDODAGNode.
    //- Discuss this with the team before making it a part of the main DAG. 
    //- It is, however, deleted appropriately when the RefFace is deleted.
    //-***************************************************************
      
  void RefFace_add_hard_point(void *this_ref_face, /* RefVertex* */ void *ref_vertex_ptr )
{
  RefFace *temp_ref_face = (RefFace *) this_ref_face;

  RefVertex *temp_ref_vertex_ptr = (RefVertex*) ref_vertex_ptr;
  
  temp_ref_face->add_hard_point(temp_ref_vertex_ptr);
}
    //- Add a hard point to this RefFace. 
    //-***************************************************************
    //- MJP Note:
    //- Currently, the new RefVertex that is created is NOT PART OF THE
    //- main DAG datastructure. The new RefVertex, however, has its
    //- own little mini-DAG which consists of a single CDODAGNode.
    //- Discuss this with the team before making it a part of the main DAG. 
    //- It is, however, deleted appropriately when the RefFace is deleted.
    //-***************************************************************

    /* other functions */

///**********Graphics Related Functions*********//
  void RefFace_draw_my_edges (void *this_ref_face,  int color)
{
  RefFace *temp_ref_face = (RefFace *) this_ref_face;

  temp_ref_face->draw_my_edges(color);
}
  
  void RefFace_draw_normal   (void *this_ref_face, int color,
                              double length,
                              /* RefVolume * */ void *wrt_vol)
{
  RefFace *temp_ref_face = (RefFace *) this_ref_face;

  RefVolume *temp_wrt_vol = (RefVolume *) wrt_vol;
  
  temp_ref_face->draw_normal(color, length, temp_wrt_vol);
}

    /* Surface* */ void *RefFace_get_surface_ptr(void *this_ref_face)
{
  RefFace *temp_ref_face = (RefFace *) this_ref_face;

  return temp_ref_face->get_surface_ptr();
}
    //R Surface*
    //R- A pointer to the Surface to which the current 
    //R- face points. 
    //- This function returns a pointer to the Surface
    //- to which the current face points.
  
  EntityType RefFace_entity_type(void *) 
{ return RefFace_TYPE;}

  
  void* RefFace_get_address(void *this_ref_face, EntityType inputEntityType)
{
  RefFace *temp_ref_face = (RefFace *) this_ref_face;

  return temp_ref_face->get_address(inputEntityType);
}
    //R void*
    //R- Returned void pointer
    //I inputEntityType
    //I- The input type to get the address of.
    //- This function returns a void pointer that points to the
    //- "appropriate" portion of this object.  The appropriate
    //- portion is determined by the input EntityType variable.
    //- Returns NULL if the input type and the type of "this"
    //- are not related by inheritance.
    //- Note: The RTTI capabilities encoded in these functions
    //-       are designed to work with any form of multiple 
    //-       inheritance, as well.  Multiple inheritance is what
    //-       necessitates having this function defined in every 
    //-       class in the hierarchy.
    //- Note: This function can also be used to merely check if
    //-       an object of one type is related to another type
    //-       through inheritance.  If a non-NULL pointer is
    //-       returned, then this is true.

enum CubitStatus RefFace_setup_use_facets( void *this_ref_face, 
                                             /* DLCubitFacetList &*/ void ***facet_list, 
                                           int *facet_list_size,
                                           int interp_order) 
{
  RefFace *temp_ref_face = (RefFace *) this_ref_face;

  DLCubitFacetList temp_list;
  COPY_ARRAY_TO_LIST(*facet_list, *facet_list_size, temp_list);

  return temp_ref_face->setup_use_facets(temp_list, interp_order);
}
  
    /* set up the underlying geometry to use a faceted base representation. */

void RefFace_stop_use_facets(void *RefFace_ptr) 
{
  RefFace *temp_ref_face = (RefFace *) RefFace_ptr;
  temp_ref_face->stop_use_facets();
}

enum CubitBoolean RefFace_get_use_facets(void *RefFace_ptr) 
{
  RefFace *temp_ref_face = (RefFace *) RefFace_ptr;
  return temp_ref_face->get_use_facets();
}
  
