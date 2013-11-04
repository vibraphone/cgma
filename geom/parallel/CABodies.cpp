//- Class:          CABodies
//- Description:    Cubit Attribute for a list of bodies the entity is part of.
//- Author: Hong-Jun Kim
//- Version:

#include "CABodies.hpp"
#include "BasicTopologyEntity.hpp"
#include "Body.hpp"
#include "RefEntityName.hpp"
#include "CastTo.hpp"
#include "TDParallel.hpp"

CubitAttrib* CABodies_creator(RefEntity* entity, const CubitSimpleAttrib& p_csa)
{
  CABodies *new_attrib = NULL;
  CubitSimpleAttrib csa = p_csa;
  new_attrib = new CABodies(entity, &csa);
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

  std::vector<int> i_list = csa_ptr->int_data_list();

  // first, the ints
  if (i_list.size() > 0)
  {
    m_interface = i_list[0]; // is interface

    if (i_list.size() > 1)
    {
      m_uniqueID = i_list[1]; // unique ID

      if(i_list.size() > 2)
      {
        // shared bodies
        int num_list = i_list[2];
        for (int i = 0; i < num_list; i++) 
          m_sharedBodies.append(i_list[3+i]);

        // shared procs
        int new_start = 3 + num_list;
        if(i_list.size() > new_start + 1)
        {
          num_list = i_list[new_start];
          for (int i = 1; i <= num_list; i++) 
           m_sharedProcs.append(i_list[new_start + i]);

          // ghost procs
          new_start += num_list+1;
          if(i_list.size() > new_start + 1)
          {
           num_list = i_list[new_start];
           for (int i = 1; i <= num_list; i++) 
             m_ghostProcs.append(i_list[new_start + i]);
          } 
        }
      }
    }
  }
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
    if (par->is_interface() != m_interface) {
      PRINT_ERROR("TDParallel interface check is failed for %s %d.\n",
		  attrib_owner()->class_name(), attrib_owner()->id());
      return CUBIT_FAILURE;
    }

    // check to make sure it's the same body list
    par->get_shared_body_list()->reset();
    m_sharedBodies.reset();
    int size = par->get_shared_body_list()->size();

    for (int i = 0; i < size; i++) {
      if (par->get_shared_body_list()->get_and_step() != m_sharedBodies.get_and_step()) {
	PRINT_ERROR("Different body found for %s %d.\n",
		    attrib_owner()->class_name(), attrib_owner()->id());
	return CUBIT_FAILURE;
      }
    }

    par->get_shared_proc_list()->reset();
    m_sharedProcs.reset();
    size = par->get_shared_proc_list()->size();

    for (int i = 0; i < size; i++) {
      if (par->get_shared_proc_list()->get_and_step() != m_sharedProcs.get_and_step()) {
	PRINT_ERROR("Different processor found for %s %d.\n",
		    attrib_owner()->class_name(), attrib_owner()->id());
	return CUBIT_FAILURE;
      }
    }

    par->get_ghost_proc_list()->reset();
    m_ghostProcs.reset();
    size = par->get_ghost_proc_list()->size();

    for (int i = 0; i < size; i++) {
      if (par->get_ghost_proc_list()->get_and_step() != m_ghostProcs.get_and_step()) {
	PRINT_ERROR("Different ghost processor found for %s %d.\n",
		    attrib_owner()->class_name(), attrib_owner()->id());
	return CUBIT_FAILURE;
      }
    }

    if ((int)par->get_unique_id() != m_uniqueID) {
      PRINT_ERROR("Different unique ID found for %s %d.\n",
		  attrib_owner()->class_name(), attrib_owner()->id());
      return CUBIT_FAILURE;
    }
  }
  else {
    // else make a new one
    par = new TDParallel(attrib_owner(), &m_sharedBodies, &m_sharedProcs,
			 &m_ghostProcs, m_uniqueID, m_interface);
  }

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
    m_interface = td_par->is_interface();
    m_uniqueID = td_par->get_unique_id();
    int size = td_par->get_shared_body_list()->size();
    td_par->get_shared_body_list()->reset();
    m_sharedBodies.clean_out();

    for (int i = 0; i < size; i++) {
      m_sharedBodies.append(td_par->get_shared_body_list()->get_and_step());
    }

    size = td_par->get_shared_proc_list()->size();
    td_par->get_shared_proc_list()->reset();
    m_sharedProcs.clean_out();

    for (int i = 0; i < size; i++) {
      m_sharedProcs.append(td_par->get_shared_proc_list()->get_and_step());
    }

    size = td_par->get_ghost_proc_list()->size();
    td_par->get_ghost_proc_list()->reset();
    m_ghostProcs.clean_out();

    for (int i = 0; i < size; i++) {
      m_ghostProcs.append(td_par->get_ghost_proc_list()->get_and_step());
    }
	
    if (delete_attrib() == CUBIT_TRUE) delete_attrib(CUBIT_FALSE);
  }

  return CUBIT_SUCCESS;
}

CubitSimpleAttrib CABodies::cubit_simple_attrib()
{
  std::vector<CubitString> cs_list;
  std::vector<int> i_list;

  // attribute internal name
  cs_list.push_back(att_internal_name());

  // is interface
  i_list.push_back(m_interface);

  // unique ID
  i_list.push_back(m_uniqueID);

  // shared bodies
  i_list.push_back(m_sharedBodies.size());
  int i;
  for (i = m_sharedBodies.size(); i > 0; i--) {
    i_list.push_back(m_sharedBodies.get_and_step());
  }

  // shared procs
  i_list.push_back(m_sharedProcs.size());
  for (i = m_sharedProcs.size(); i > 0; i--) {
    i_list.push_back(m_sharedProcs.get_and_step());
  }

  // ghost procs
  i_list.push_back(m_ghostProcs.size());
  for (i = m_ghostProcs.size(); i > 0; i--) {
    i_list.push_back(m_ghostProcs.get_and_step());
  }
  
  CubitSimpleAttrib csattrib =  CubitSimpleAttrib(&cs_list, NULL, &i_list);
  
  return csattrib;
}

CubitStatus CABodies::reset()
{
  m_sharedBodies.clean_out();
  return CUBIT_SUCCESS;
}







