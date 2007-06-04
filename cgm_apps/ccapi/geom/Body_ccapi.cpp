#include "Body_ccapi.h"

#include "Body.hpp"
#include "CubitBox.hpp"

  /* topology */

EntityType Body_get_sense_entity_type(void *) 
{ return CoVolume_TYPE; }
  /*R EntityType
   *R- The type of the specific sense entity.
   *- This function returns the type of the specific sense entity.
   */

EntityType Body_get_basic_topology_entity_type(void *) 
{ return InvalidEntity_TYPE; }
  /*R EntityType
   *R- The type of the corresponding BasicTopologyEntity type.
   *- This function returns the type of the corresponding 
   *- BasicTopologyEntity .
   */

EntityType Body_get_child_basic_topology_entity_type(void *) 
{ return RefVolume_TYPE; }
  /*R EntityType
   *R- A type value.
   *- This function returns the type of the child RefEntity of  
   *- Body, which is the type value of the RefVolume class.
   */
      
EntityType Body_get_child_ref_entity_type(void *) 
{ return RefVolume_TYPE; }
  /*R EntityType
   *R- A type value.
   *- This function returns the type of the child RefEntity of  
   *- Body, which is the type value of the RefVolume class.
   */

EntityType Body_get_parent_ref_entity_type(void *) 
{ return InvalidEntity_TYPE; }
  /*R EntityType
   *R- A type value.
   *- This function returns the type of the child RefEntity of  
   *- Body, which is the type value of the RefVolume class.
   */

    /* geometry */

CubitBoxStruct Body_bounding_box(void *this_body) 
{
  Body *temp_body = (Body *) this_body;
  return temp_body->bounding_box();
}
  /*- Returns the bounding box of this body
   */

CubitVectorStruct Body_center_point(void *this_body)
{
  Body *temp_body = (Body *) this_body;
  return temp_body->center_point();
}
  
  /*- Return a CubitVector set to the centroid of this body
   */

double Body_measure(void *this_body)
{
  Body *temp_body = (Body *) this_body;
  return temp_body->measure();
}

  /*- Heading: geometry modification functions */

/* Body* */ void * Body_copy(void *this_body)
{
  Body *temp_body = (Body *) this_body;
  return temp_body->copy();
}
  /*- Copy this Body to another body
   */
  
CubitStatus Body_move(void *this_body, double dx, double dy, double dz,
                 CubitBoolean update_now)
{
  Body *temp_body = (Body *) this_body;
  return temp_body->move(dx, dy, dz, update_now);
}
  /*- Move this Body by dx, dy and dz
   */
  
CubitStatus Body_rotate(void *this_body, double x, double y, double z, double angle,
                   CubitBoolean update_now)
{
  Body *temp_body = (Body *) this_body;
  return temp_body->rotate(x, y, z, angle, update_now);
}
  /*- Rotates this Body about a vector defined by x, y and z
   */
  
CubitStatus Body_scale(void *this_body, double scaling_factor,
                  CubitBoolean update_now)
{
  Body *temp_body = (Body *) this_body;
  return temp_body->scale(scaling_factor, update_now);
}
  /*- Scales this Body by the factor, scaling_factor
   */

CubitStatus Body_restore(void *this_body, CubitBoolean update_now)
{
  Body *temp_body = (Body *) this_body;
  return temp_body->restore(update_now);
}
  /*- Restores this Body by replacing the transformation matrix 
   *- associated with it with a unit matrix
   */

CubitStatus Body_reverse(void *this_body, CubitBoolean update_now)
{
  Body *temp_body = (Body *) this_body;
  return temp_body->reverse(update_now);
}
  /*- reverses the orientation of the faces on this body
   */
  
CubitStatus Body_transform_position(void *this_body, /* CubitVector const& */ CubitVectorStruct position_vector,
                                      /* CubitVector & */ CubitVectorStruct *transformed_vector)
{
  Body *temp_body = (Body *) this_body;

  CubitVector temp_position_vector(position_vector);
  CubitVector temp_transformed_vector;
  
  CubitStatus status = temp_body->transform_position(temp_position_vector, temp_transformed_vector);

  *transformed_vector = temp_transformed_vector;

  return status;
}
  /*R CubitStatus
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

int Body_validate(void *this_body)
{
  Body *temp_body = (Body *) this_body;
  return temp_body->validate();
}
  /*- do a measure and api entity check.
   */
  
/* CubitString */ const char *Body_measure_label(void *this_body)
{
  Body *temp_body = (Body *) this_body;
  return temp_body->measure_label().c_str();
}
  /*- Functions related to "measuring" the Body
   */

EntityType Body_entity_type(void *)
{ return Body_TYPE ; }

void* Body_get_address(void *this_body, EntityType inputEntityType)
{
  Body *temp_body = (Body *) this_body;
  return temp_body->get_address(inputEntityType);
}
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
      
void Body_copiedFromId_1(void *this_body,  int Id )
{
  Body *temp_body = (Body *) this_body;
  temp_body->copiedFromId(Id);
}

int Body_copiedFromId_2 (void *this_body)
{
  Body *temp_body = (Body *) this_body;
  return temp_body->copiedFromId();
}
  /*- Sets/Gets the id of the body that it was copied from */
