/*-------------------------------------------------------------------------
 * Filename      : Body.hpp
 *
 * Purpose       : This class represents the topological entity, Body,
 *                 which is the highest level abstraction for a complete
 *                 geometric model.
 *
 * Special Notes : Body is a GroupingEntity in the TopologyEntity hierarchy.
 *                 It is also a RefEntity and get a lot of its functionality
 *                 and data from the RefEntity base class for meshing
 *                 purposes.
 *
 *                 The user is provided a handle to Bodies.
 *
 * Creator       : 
 *
 * Creation Date : 
 *
 * Owner         : 
 *-------------------------------------------------------------------------
 */

#ifndef BODY_CCAPI_HPP
#define BODY_CCAPI_HPP

#include "EntityType.h"
#include "CubitVectorStruct.h"
#include "CubitBoxStruct.h"
#include "CubitDefines.h"

#ifdef __cplusplus
extern "C" {
#endif

  /* topology */

enum EntityType Body_get_sense_entity_type(void *this_body);
  /*R EntityType
   *R- The type of the specific sense entity.
   *- This function returns the type of the specific sense entity.
   */

enum EntityType Body_get_basic_topology_entity_type(void *this_body);
  /*R EntityType
   *R- The type of the corresponding BasicTopologyEntity type.
   *- This function returns the type of the corresponding 
   *- BasicTopologyEntity .
   */

enum EntityType Body_get_child_basic_topology_entity_type(void *this_body);
  /*R EntityType
   *R- A type value.
   *- This function returns the type of the child RefEntity of  
   *- Body, which is the type value of the RefVolume class.
   */
      
enum EntityType Body_get_child_ref_entity_type(void *this_body);
  /*R EntityType
   *R- A type value.
   *- This function returns the type of the child RefEntity of  
   *- Body, which is the type value of the RefVolume class.
   */

enum EntityType Body_get_parent_ref_entity_type(void *this_body);
  /*R EntityType
   *R- A type value.
   *- This function returns the type of the child RefEntity of  
   *- Body, which is the type value of the RefVolume class.
   */

    /* geometry */

struct CubitBoxStruct Body_bounding_box(void *this_body);
  /*- Returns the bounding box of this body
   */

struct CubitVectorStruct Body_center_point(void *this_body);
  /*- Return a CubitVector set to the centroid of this body
   */

double Body_measure(void *this_body);

  /*- Heading: geometry modification functions */

/* Body* */ void * Body_copy(void *this_body);
  /*- Copy this Body to another body
   */
  
enum CubitStatus Body_move(void *this_body, double dx, double dy, double dz,
                 enum CubitBoolean update_now);
  /*- Move this Body by dx, dy and dz
   */
  
enum CubitStatus Body_rotate(void *this_body, double x, double y, double z, double angle,
                   enum CubitBoolean update_now);
  /*- Rotates this Body about a vector defined by x, y and z
   */
  
enum CubitStatus Body_scale(void *this_body, double scaling_factor,
                  enum CubitBoolean update_now);
  /*- Scales this Body by the factor, scaling_factor
   */

enum CubitStatus Body_restore(void *this_body, enum CubitBoolean update_now);
  /*- Restores this Body by replacing the transformation matrix 
   *- associated with it with a unit matrix
   */

enum CubitStatus Body_reverse(void *this_body, enum CubitBoolean update_now);
  /*- reverses the orientation of the faces on this body
   */
  
enum CubitStatus Body_transform_position(void *this_body, /* CubitVector const& */ struct CubitVectorStruct position_vector,
                                      /* CubitVector & */ struct CubitVectorStruct *transformed_vector);
  /*R enum CubitStatus
   *R- CUBIT_SUCCESS/FAILURE
   *I position_vector
   *I- The input vector that needs to be transformed.
   *O transformed_vector
   *O- The transformed CubitVector -- sent in by reference.
   *- This function will transform the input position_vector using
   *- the transformation matrix associated with the underlying geometry
   *- of this Body.  If there is none, then CUBIT_FAILURE is returned
   *- and transformed_vector is left untouched.
   *- Note: The input vector *must* have 3 components.
   */

int Body_validate(void *this_body);
  /*- do a measure and api entity check.
   */
  
/* CubitString */ const char *Body_measure_label(void *this_body);
  /*- Functions related to "measuring" the Body
   */

enum EntityType Body_entity_type(void *this_body);

void* Body_get_address(void *this_body, enum EntityType inputEntityType);
  /*R void*
   *R- Returned void pointer
   *I inputEntityType
   *I- The input type to get the address of.
   *- This function returns a void pointer that points to the
   *- "appropriate" portion of this object.  The appropriate
   *- portion is determined by the input EntityType variable.
   *- Returns NULL if the input type and the type of "this"
   *- are not related by inheritance.
   *- Note: The RTTI capabilities encoded in these functions
   *-       are designed to work with any form of multiple 
   *-       inheritance, as well.  Multiple inheritance is what
   *-       necessitates having this function defined in every 
   *-       class in the hierarchy.
   *- Note: This function can also be used to merely check if
   *-       an object of one type is related to another type
   *-       through inheritance.  If a non-NULL pointer is
   *-       returned, then this is true.
   */
      
void Body_copiedFromId_1(void *this_body,  int Id );
int Body_copiedFromId_2(void *this_body);
  /*- Sets/Gets the id of the body that it was copied from */

#ifdef __cplusplus
}
#endif

#endif
