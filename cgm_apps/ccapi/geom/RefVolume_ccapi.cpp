#include "RefVolume_ccapi.h"
#include "RefVolume.hpp"
#include "RefFace.hpp"
#include "Lump.hpp"

    /* topology */

enum EntityType RefVolume_get_grouping_entity_type(void *) 
{ return Shell_TYPE; }
  /* R EntityType */
  /* R- The type of the corresponding GroupingEntity type. */
  /* - This function returns the type of the corresponding */
  /* - GroupingEntity . */

enum EntityType RefVolume_get_sense_entity_type(void *)
{ return CoVolume_TYPE; }
  /* R EntityType */
  /* R- The type of SenseEntity associated with this object */
  /* - This function returns the type of SenseEntity associated with  */
  /* - this BasicTopologyEntity.  */

/* Lump* */ void * RefVolume_get_lump_ptr(void *this_ref_volume) 
{
  RefVolume *temp_ref_volume = (RefVolume*) this_ref_volume;
  
  return temp_ref_volume->get_lump_ptr();
}

  /* R Lump* */
  /* R- A pointer to the Lump to which the current  */
  /* R- volume  points.  */
  /* - This function returns a pointer to the Lump */
  /* - to which the current volume points. */

enum EntityType RefVolume_get_child_ref_entity_type(void *) 
{ return RefFace_TYPE; }
  /* R EntityType */
  /* R- A type value. */
  /* - This function returns the type of the child RefEntity of   */
  /* - RefVolume, which is the type value of the RefFace class. */
      
enum EntityType RefVolume_get_parent_ref_entity_type(void *) 
{ return InvalidEntity_TYPE; }
  /* R EntityType */
  /* R- A type value. */
  /* - This function returns the type of the parent RefEntity of   */
  /* - RefVolume -- and there isn't one. */

enum EntityType RefVolume_get_topology_bridge_type(void *)
{ return Lump_TYPE ;}
  /* R EntityType */
  /* R- The type of GeometryEntity associated with this object */
  /* - This function returns the type of GeometryEntity associated with  */
  /* - this BasicTopologyEntity.  */
      
int RefVolume_num_boundary_components(void *this_ref_volume)
{
  RefVolume *temp_ref_volume = (RefVolume*) this_ref_volume;
  
  return temp_ref_volume->num_boundary_components();
}

  /* - Return the number of connected components bounding this volume. */
  /* - Note: this counts the number of ref_edge-connected ref_face  */
  /* - components, which may not be the same as the number of acis shells. */

  /* geometry */

CubitVectorStruct RefVolume_center_point(void *this_ref_volume) 
{
  RefVolume *temp_ref_volume = (RefVolume*) this_ref_volume;
  
  return temp_ref_volume->center_point();
}
  /* - Returns centroid of the RefVolume */
  
/* **********Graphics Related Functions********** */
int RefVolume_dimension(void *this_ref_volume)  
{
  RefVolume *temp_ref_volume = (RefVolume*) this_ref_volume;
  
  return temp_ref_volume->dimension();
}
  /* - returns dimension of the actual entity.  */

    /* other functions */

enum EntityType RefVolume_entity_type(void *) 
{ return RefVolume_TYPE ; }

void * RefVolume_get_address(void *this_ref_volume, enum EntityType inputEntityType) 
{
  RefVolume *temp_ref_volume = (RefVolume*) this_ref_volume;
  
  return temp_ref_volume->get_address(inputEntityType);
}
  /* R void * */
  /* R- Returned void pointer */
  /* I inputEntityType */
  /* I- The input type to get the address of. */
  /* - This function returns a void pointer that points to the */
  /* - "appropriate" portion of this object.  The appropriate */
  /* - portion is determined by the input EntityType variable. */
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
      
/* CubitString */ const char *RefVolume_measure_label(void *this_ref_volume) 
{
  RefVolume *temp_ref_volume = (RefVolume*) this_ref_volume;
  
  return temp_ref_volume->measure_label().c_str();
}

void RefVolume_draw_my_faces(void *this_ref_volume, int color) 
{
  RefVolume *temp_ref_volume = (RefVolume*) this_ref_volume;
  
  temp_ref_volume->draw_my_faces(color);
}

