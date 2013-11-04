//- File:           CAEntityColor.cpp
//- Owner:          Joel Kopp
//- Description:    Cubit Attribute for entity colors.
//- Checked By:
//- Version:

#include "CAEntityColor.hpp"
#include "RefEntity.hpp"
#include "RefVolume.hpp"
#include "CastTo.hpp"

CubitAttrib* CAEntityColor_creator(RefEntity* entity, const CubitSimpleAttrib &p_csa)
{
  return new CAEntityColor(entity, p_csa);
}

CAEntityColor::CAEntityColor(RefEntity* new_attrib_owner,
                       const CubitSimpleAttrib &csa_ptr)
        : CubitAttrib(new_attrib_owner)
{
  entityColor = 0;

  if(!csa_ptr.isEmpty())
  {
    const std::vector<int>& i_list = csa_ptr.int_data_list();
    assert(i_list.size() > 0);
    entityColor = i_list[0];
  }
}

CAEntityColor::~CAEntityColor()
{
}


CubitStatus CAEntityColor::actuate()
{
   if ( hasActuated)
      return CUBIT_SUCCESS;
   
   if ( !attribOwnerEntity )
      return CUBIT_FAILURE;
   
   deleteAttrib = CUBIT_FALSE;
   
   //- If actuating after import, change the color.  Else (auto actuated) color is already changed...
   if( attribOwnerEntity->color() == CUBIT_DEFAULT_COLOR || 
	   dynamic_cast<RefVolume*>(attribOwnerEntity) )
   {
      attribOwnerEntity->color(entityColor);
   }
   
   hasActuated = CUBIT_TRUE;
   return CUBIT_SUCCESS;
}

CubitStatus CAEntityColor::update()
{
      delete_attrib(CUBIT_TRUE);
      return CUBIT_SUCCESS;

/*
   if ( hasUpdated ) 
      return CUBIT_SUCCESS;
   
   if ( !attribOwnerEntity)
      return CUBIT_FAILURE;
   
   entityColor = attribOwnerEntity->color();
   
   // set the updated flag
   hasUpdated = CUBIT_TRUE;

   if( entityColor == -1 )
   {
      delete_attrib(CUBIT_TRUE);
      return CUBIT_SUCCESS;
   }
   
   return CUBIT_SUCCESS;
   */
}

void CAEntityColor::merge_owner(CubitAttrib *deletable_attrib)
{
    // take the id with the lowest value
  CAEntityColor *other_caecolor = CAST_TO(deletable_attrib, CAEntityColor);
  
  if (other_caecolor)
    entityColor = other_caecolor->color();
}

CubitSimpleAttrib CAEntityColor::cubit_simple_attrib()
{
  std::vector<CubitString> cs_list;
  std::vector<int> i_list;

  i_list.push_back ( entityColor );
    
  cs_list.push_back(att_internal_name());

  return CubitSimpleAttrib(&cs_list, NULL, &i_list);
}

void CAEntityColor::print()
{
    // print info on this attribute
  
  PRINT_INFO("CAEntityColor: owner = %s %d:   color ",
             attribOwnerEntity->class_name(), attribOwnerEntity->id());
  if ( entityColor > -1 )
     PRINT_INFO("                             %d\n",
     entityColor);
  else
     PRINT_INFO("\n");
}
