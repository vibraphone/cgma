//- File:           CAEntitySense.cpp
//- Owner:          Corey Ernst
//- Description:    Cubit Attribute for entity colors.
//- Checked By:
//- Version:

#include "CAEntitySense.hpp"
#include "RefEntity.hpp"
#include "RefEdge.hpp"
#include "RefFace.hpp"
#include "BridgeManager.hpp"
#include "BasicTopologyEntity.hpp"
#include "CastTo.hpp"


CubitAttrib* CAEntitySense_creator(RefEntity* entity, CubitSimpleAttrib *p_csa)
{
  CAEntitySense *new_attrib = NULL;
  if (NULL == p_csa)
  {
    new_attrib = new CAEntitySense(entity);
  }
  else
  {
    new_attrib = new CAEntitySense(entity, p_csa);
  }

  return new_attrib;
}

CAEntitySense::CAEntitySense(RefEntity* new_attrib_owner,
                       CubitSimpleAttrib *csa_ptr)
        : CubitAttrib(new_attrib_owner)
{
   assert ( csa_ptr != 0 );
   
   std::vector<int, std::allocator<int> > i_list = csa_ptr->int_data_list();
   
   assert( i_list.size() > 0);
   int i =  i_list[0];
   if( i == -1 ) 
     entitySense = CUBIT_UNKNOWN;
   else if( i == 0 )
     entitySense = CUBIT_FORWARD;
   else if( i == 1 )
     entitySense = CUBIT_REVERSED;
}

CAEntitySense::CAEntitySense(RefEntity* new_attrib_owner)
        : CubitAttrib(new_attrib_owner)
{
  entitySense = CUBIT_FORWARD;
}

CAEntitySense::~CAEntitySense()
{
}


CubitStatus CAEntitySense::actuate()
{
   if ( hasActuated)
      return CUBIT_SUCCESS;
   
   if ( !attribOwnerEntity )
      return CUBIT_FAILURE;
   
   deleteAttrib = CUBIT_FALSE;  

   int dimension = attribOwnerEntity->dimension();
   if( 1 == dimension || 2 == dimension )
   {
     BasicTopologyEntity* bte_ptr = CAST_TO(attribOwnerEntity, BasicTopologyEntity);
     CubitSense tmp_sense = CUBIT_UNKNOWN;
     if( bte_ptr )     
     {
       tmp_sense = bte_ptr->bridge_manager()->topology_bridge()->bridge_sense();
       if( tmp_sense != entitySense )
       {
         if( 1==dimension )
         {
           RefEdge *ref_edge = dynamic_cast<RefEdge*>( bte_ptr );
           ref_edge->reverse_tangent();
         }
         else if( 2 == dimension )
         {
           RefFace *ref_face = dynamic_cast<RefFace*>( bte_ptr );
           ref_face->reverse_normal();
         }
       }
     }
   }
   
   hasActuated = CUBIT_TRUE;
   return CUBIT_SUCCESS;
}

CubitStatus CAEntitySense::update()
{
  int dimension = attribOwnerEntity->dimension();
  if( 1 == dimension || 2 == dimension )
  {
    BasicTopologyEntity* bte_ptr = CAST_TO(attribOwnerEntity, BasicTopologyEntity);
    entitySense = CUBIT_UNKNOWN;
    if( bte_ptr )     
      entitySense = bte_ptr->bridge_manager()->topology_bridge()->bridge_sense();
  }  

  if( entitySense == CUBIT_FORWARD || entitySense == CUBIT_UNKNOWN )  
    delete_attrib(CUBIT_TRUE);

  return CUBIT_SUCCESS;
}

CubitSimpleAttrib CAEntitySense::cubit_simple_attrib()
{
  std::vector<CubitString> cs_list;
  std::vector<int, std::allocator<int> >  i_list;

  i_list.push_back( entitySense );
    
  cs_list.push_back(att_internal_name());

  CubitSimpleAttrib csattrib_ptr(&cs_list, NULL, &i_list);
  
  return csattrib_ptr;
}

void CAEntitySense::print()
{
    // print info on this attribute
  
  PRINT_INFO("CAEntitySense: owner = %s %d:   color ",
             attribOwnerEntity->class_name(), attribOwnerEntity->id());
  if( entitySense == CUBIT_UNKNOWN )
     PRINT_INFO("                             CUBIT_UNKNOWN\n");
  else if( entitySense == CUBIT_FORWARD )
    PRINT_INFO("                             CUBIT_FORWARD\n");
  else if( entitySense == CUBIT_REVERSED )
    PRINT_INFO("                             CUBIT_REVERSED\n");  
}
