//- File:           CAEntityTol.cpp
//- Owner:          W. Roshan Quadros
//- Description:    Cubit Attribute for entity ids.
//- Checked By:
//- Version:


#include "CAEntityTol.hpp"
#include "CAMergePartner.hpp"
#include "TDUniqueId.hpp"
#include "TopologyBridge.hpp"
#include "RefEntity.hpp"
#include "RefVertex.hpp"
#include "RefEdge.hpp"
#include "RefFace.hpp"
#include "CastTo.hpp"
#include "MergeTool.hpp"
#include "RefEntityFactory.hpp"
#include "GeometryQueryTool.hpp"
#include "GSaveOpen.hpp"
#include "CADeferredAttrib.hpp"
#include "BasicTopologyEntity.hpp"
#include "GeometryEntity.hpp"

CubitAttrib* CAEntityTol_creator(RefEntity* entity, const CubitSimpleAttrib &p_csa)
{
    return new CAEntityTol(entity, p_csa);
}

CAEntityTol::CAEntityTol(RefEntity* new_attrib_owner,
                       const CubitSimpleAttrib &csa_ptr)
        : CubitAttrib(new_attrib_owner)
{
  entityTol = 0.0;
  if(!csa_ptr.isEmpty())
  {
   PRINT_DEBUG_95( "Creating ENTITY_TOL attribute from CSA for %s %d\n",
      (attribOwnerEntity ? attribOwnerEntity->class_name() : "(none)"),
      (attribOwnerEntity ? attribOwnerEntity->id() : 0));
   
   const std::vector<double>& d_list = csa_ptr.double_data_list();

   assert(d_list.size() == 1);
   entityTol =  d_list[0];
  }
}

CAEntityTol::CAEntityTol(RefEntity* new_attrib_owner)
        : CubitAttrib(new_attrib_owner)
{
  entityTol = 0.0;
    
  PRINT_DEBUG_95( "Creating ENTITY_TOL attribute for %s %d\n",
              (attribOwnerEntity ? attribOwnerEntity->class_name() : "(none)"),
              (attribOwnerEntity ? attribOwnerEntity->id() : 0));
}

CAEntityTol::~CAEntityTol()
{
}

CubitStatus CAEntityTol::reset()
{

  return CUBIT_SUCCESS;
}

CubitStatus CAEntityTol::actuate()
{
   if ( hasActuated)
      return CUBIT_SUCCESS;
   
   if ( !attribOwnerEntity )
      return CUBIT_FAILURE;
   
   deleteAttrib = CUBIT_FALSE;
   
   //- If actuating after import, change the tol.  Else (auto actuated) tol is already changed...
   if( attribOwnerEntity->local_tolerance() == 0.0 )
   {
      attribOwnerEntity->local_tolerance(entityTol);
   }
   
   hasActuated = CUBIT_TRUE;
   return CUBIT_SUCCESS;
}

CubitStatus CAEntityTol::update()
{
  //delete_attrib(CUBIT_TRUE);
  //return CUBIT_SUCCESS;

  if (hasUpdated) return CUBIT_SUCCESS;
  
  PRINT_DEBUG_95( "Updating ENTITY_TOL attribute for %s %d\n",
              attribOwnerEntity->class_name(), attribOwnerEntity->id());

    // set the updated flag
  hasUpdated = CUBIT_TRUE;

    // first, remove this attrib in its old form from the geometry entity
  if (hasWritten == CUBIT_TRUE) 
  {
    attribOwnerEntity->remove_attrib_geometry_entity(this);
    hasWritten = CUBIT_FALSE;
  }
  
  double local_tol = attribOwnerEntity->local_tolerance();
  if( local_tol == 0.0)
  {
    delete_attrib(CUBIT_TRUE);
    return CUBIT_SUCCESS;
  }
  else
  {
         
    // reset the delete flag if it was set before
    delete_attrib(CUBIT_FALSE);
    
      // now, write to geometry entity
    entityTol = local_tol;
    attribOwnerEntity->write_specific_cubit_attrib(this);
  }
  return CUBIT_SUCCESS;
}

void CAEntityTol::merge_owner(CubitAttrib *deletable_attrib)
{
    // take the id with the lowest value
  CAEntityTol *other_ca_tol = CAST_TO(deletable_attrib, CAEntityTol);
  
  if (other_ca_tol && ( entityTol == 0.0 || entityTol < other_ca_tol->tolerance() ))
    entityTol = other_ca_tol->tolerance();
}

CubitSimpleAttrib CAEntityTol::cubit_simple_attrib()
{
  std::vector<CubitString> cs_list;
  std::vector<double> d_list;
  std::vector<int, std::allocator<int> > i_list;

  d_list.push_back ( entityTol );

  cs_list.push_back(att_internal_name());

  CubitSimpleAttrib csattrib_ptr(&cs_list, &d_list, &i_list);
 
  return csattrib_ptr;
}

void CAEntityTol::print()
{
    // print info on this attribute
  
  PRINT_INFO("CAEntityTol: owner = %s %d:  tolerance =%f\n",
             attribOwnerEntity->class_name(), attribOwnerEntity->id(),
             entityTol);
}


CubitSimpleAttrib CAEntityTol::split_owner()
{
    // if this entity is to be split, pass back a simple attribute with
    // duplicate tol data to be put on new entity
  PRINT_DEBUG_95("CAEntityName::split_owner()\n");
  update();
  return cubit_simple_attrib();
}

