//- Class:          CABodies
//- Description:    Cubit Attribute for a list of bodies the entity is part of.
//- Author: Hong-Jun Kim
//- Checked By: Tim Tautges
//- Version:

#include "CABodies.hpp"
#include "BasicTopologyEntity.hpp"
#include "Body.hpp"
#include "RefEntityName.hpp"
#include "CastTo.hpp"
#include "TDParallel.hpp"

CubitAttrib* CABodies_creator(RefEntity* entity, CubitSimpleAttrib *p_csa)
{
  CABodies *new_attrib = NULL;
  if (NULL == p_csa)
  {
    new_attrib = new CABodies(entity);
  }
  else
  {
    new_attrib = new CABodies(entity, p_csa);
  }

  return new_attrib;
}


CABodies::CABodies(RefEntity* new_attrib_owner,
                           CubitSimpleAttrib *csa_ptr)
        : CubitAttrib(new_attrib_owner)
{
  assert ( csa_ptr != NULL );
  if (DEBUG_FLAG(138))
  {
    PRINT_DEBUG_138( "Creating BODIES attribute from CSA for %s %d\n",
		     (attribOwnerEntity ? attribOwnerEntity->class_name() : "(none)"),
		     (attribOwnerEntity ? attribOwnerEntity->id() : 0));
  }

  //initialize();
  
  DLIList<int*> *i_list = csa_ptr->int_data_list();

  // first, the ints
  i_list->reset();
  
  int num_bodies = *(i_list->get_and_step());
  
  // bodyID
  int i;
  for (i = num_bodies; i > 0; i--)
    bodyUniqueId.append(*(i_list->get_and_step()));
}

CABodies::CABodies(RefEntity* new_attrib_owner)
  : CubitAttrib(new_attrib_owner)
{
  if (DEBUG_FLAG(138))
  {
  PRINT_DEBUG_138( "Creating BODIES attribute for %s %d\n",
              (attribOwnerEntity ? attribOwnerEntity->class_name() : "(none)"),
              (attribOwnerEntity ? attribOwnerEntity->id() : 0));
  }
}

CABodies::~CABodies()
{
}

CubitStatus CABodies::actuate()
{
  CubitStatus status = CUBIT_SUCCESS;

  if (hasActuated == CUBIT_TRUE) return CUBIT_SUCCESS;

  if (DEBUG_FLAG(138))
  {
    PRINT_DEBUG_138( "Actuating BODIES attribute for %s %d\n",
		     attribOwnerEntity->class_name(), attribOwnerEntity->id());
  }

  // create a TDParallel for the entity, if it doesn't already exist
  TDParallel *par = (TDParallel *) attrib_owner()->get_TD(&TDParallel::is_parallel);

  if (par != NULL) {
    // check to make sure it's the same body list
    par->body_unique_id_list()->reset();
    bodyUniqueId.reset();
    int size = par->body_unique_id_list()->size();

    for (int i = 0; i < size; i++) {
      if (par->body_unique_id_list()->get_and_step() != bodyUniqueId.get_and_step()) {
	PRINT_ERROR("Different body found for %s %d.\n",
		    attrib_owner()->class_name(), attrib_owner()->id());
	return CUBIT_FAILURE;
      }
    }
  }
  else {
    // else make a new one
    par = new TDParallel(attrib_owner(), &bodyUniqueId);
  }

  //status = par->set_local_non_local_list();
  delete_attrib(CUBIT_TRUE);
  hasActuated = CUBIT_TRUE;

  return status;
}

CubitStatus CABodies::update()
{
  if (hasUpdated) return CUBIT_SUCCESS;
  
  if (DEBUG_FLAG(138))
  {
    PRINT_DEBUG_138( "Updating BODIES attribute for %s %d\n",
              attribOwnerEntity->class_name(), attribOwnerEntity->id());
  }

  // set the updated flag
  hasUpdated = CUBIT_TRUE;
  
  // if the owner has a body list, save it, otherwise delete this one
  TDParallel *td_par = (TDParallel *) attrib_owner()->get_TD(&TDParallel::is_parallel);
  
  if (td_par == NULL) {
    delete_attrib(CUBIT_TRUE);
  }
  else {
    int size = td_par->body_unique_id_list()->size();
    td_par->body_unique_id_list()->reset();
    bodyUniqueId.clean_out();

    for (int i = 0; i < size; i++) {
      bodyUniqueId.append(td_par->body_unique_id_list()->get_and_step());
    }
	
    if (delete_attrib() == CUBIT_TRUE) delete_attrib(CUBIT_FALSE);
  }

  return CUBIT_SUCCESS;
}

CubitSimpleAttrib* CABodies::cubit_simple_attrib()
{
  DLIList<CubitString*> cs_list;
  DLIList<int> i_list;

  // attribute internal name
  cs_list.append(new CubitString(att_internal_name()));

  // bodyID
  bodyUniqueId.reset();
  i_list.append(bodyUniqueId.size());
  int i;
  for (i = bodyUniqueId.size(); i > 0; i--) {
    i_list.append(bodyUniqueId.get_and_step());
  }
  
  CubitSimpleAttrib* csattrib_ptr = new CubitSimpleAttrib(&cs_list, NULL, &i_list);
  
  for (i = cs_list.size(); i--;) delete cs_list.get_and_step();

  return csattrib_ptr;
}

CubitStatus CABodies::reset()
{
  bodyUniqueId.clean_out();
  return CUBIT_SUCCESS;
}







