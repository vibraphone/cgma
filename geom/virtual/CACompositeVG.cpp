// CACompositeVG class

#include "CACompositeVG.hpp"
#include "CADeferredAttrib.hpp"
#include "CAUniqueId.hpp"
#include "TDUniqueId.hpp"
#include "CastTo.hpp"
#include "RefFace.hpp"
#include "RefEdge.hpp"
#include "BasicTopologyEntity.hpp"
#include "PartitionEntity.hpp"
#include "TopologyBridge.hpp"

#include "CompositeTool.hpp"
#include "CompositeSurface.hpp"
#include "CompositeCurve.hpp"

CubitAttrib* CACompositeVG_creator(RefEntity* entity, const CubitSimpleAttrib &p_csa)
{
  return new CACompositeVG(entity, p_csa);
}

CACompositeVG::CACompositeVG(RefEntity *owner, const CubitSimpleAttrib &simple_attrib)
        : CubitAttrib(owner)
{

  compositeId = -1;

  if(!simple_attrib.isEmpty())
  {

      // generate a simple attribute containing the data in this CA
    const std::vector<CubitString>& cs_list = simple_attrib.string_data_list();
    const std::vector<int>& i_list = simple_attrib.int_data_list();

      // (no string)

      // now the integers
      // compositeId, numSubEntities
    compositeId = i_list[0];
    int numSubEntities = i_list[1];
    int i;
    for (i = 0; i < numSubEntities; i++)
    {
      subEntityIds.append(i_list[i+2]);
    }

      // If we are constructing from a CubitSimpleAttrib,
      // then this attrib is already written
    has_written(CUBIT_TRUE);
  }
}

CubitStatus CACompositeVG::update()
{
  hasUpdated = CUBIT_TRUE;
  delete_attrib(CUBIT_TRUE);
  return CUBIT_SUCCESS;
}

CubitStatus CACompositeVG::reset()
{

    // reset the info on this CACVG; don't reset info on underlying
    // entities, since the lists on this CACVG get *assigned* to that
    // entity and not *appended*
  compositeId = -1;
  subEntityIds.clean_out();
  return CUBIT_SUCCESS;
}
  
CubitSimpleAttrib CACompositeVG::cubit_simple_attrib()
{
    // generate a simple attribute containing the data in this CA
  std::vector<CubitString> cs_list;
  std::vector<int> i_list;

    // first the string
  cs_list.push_back(att_internal_name());

    // now the integers
    // compositeId, numSubEntities
  i_list.push_back(compositeId);
  i_list.push_back(subEntityIds.size());
  int i;
  subEntityIds.reset();
  for (i = subEntityIds.size(); i > 0; i--)
    i_list.push_back(subEntityIds.get_and_step());
    

  return CubitSimpleAttrib(&cs_list, NULL, &i_list);
}

CubitStatus CACompositeVG::actuate()
{
    // actuate this CA
  if (has_actuated() == CUBIT_TRUE) return CUBIT_SUCCESS;

    // actuate all the unique id attribute first, since CA_COMPOSITE_VG
    // depends on CA_UNIQUEID
  CAUniqueId::actuate_all();
  
    // now, gather up all the sub entities; if they're not all around,
    // exit without setting the actuate flag
  DLIList<ToolDataUser*> tdu_list;
  int i;
  subEntityIds.reset();
  CubitBoolean all_actuated = CUBIT_TRUE;

  DLIList<CACompositeVG*> comp_attribs;

  for (i = subEntityIds.size(); i > 0; i--) {
    ToolDataUser *tdu = TDUniqueId::find_td_unique_id(subEntityIds.get_and_step());
    if (tdu != NULL) tdu_list.append(tdu);
    else {

      // put this CA on a list for future actuation
      //has_actuated( CUBIT_FALSE );
      bool result = CADeferredAttrib::add_unactuated_ca(this);

      if (true == result)
        PRINT_DEBUG_90("Can't actuate composite vg attrib on %s %d; adding to unactuated "
                       "list.\n",
                      attrib_owner()->class_name(), attrib_owner()->id());
      return CUBIT_FAILURE;
    }

    //The portion of code below examines the entities that this entity will be 
    //composited with.  If the CACompositeVG attribute on these entities-to-be-composited 
    //do not have their 'has_actuated' variable set to TRUE, then this Composite attribute
    //is not ready to be actuated.  All entities-to-be-composited must have their 
    //'has_actuated' flag set to TRUE to be ready to be composited
    
    RefEntity *ref_entity = CAST_TO( tdu, RefEntity );
    if( ref_entity == NULL )
      return CUBIT_FAILURE;
    else
    {
      CACompositeVG *comp_vg_attrib = (CACompositeVG *) ref_entity->get_cubit_attrib(CA_COMPOSITE_VG);
      if( comp_vg_attrib ) 
      {
        comp_attribs.append( comp_vg_attrib );
        if( ref_entity == attrib_owner() ) //skip owner of this attribute..it will be actuated in a minute
          continue;
        if( !comp_vg_attrib->has_actuated() )
          all_actuated = CUBIT_FALSE;
        continue;
      }
      else
      {
        all_actuated = CUBIT_FALSE;
        continue;  
      }
    }
  }

  has_actuated( CUBIT_TRUE );

  //if all sub entities CACompositeVG attribs haven't been here, 
  //set this one and get out 
  if( all_actuated == CUBIT_FALSE )
    return CUBIT_SUCCESS;

    // ok, we've got a list of tdu's; cast to list of ref edges or faces (currently
    // the only entities we know how to composite)
  DLIList<RefEdge*> ref_edges;
  CAST_LIST(tdu_list, ref_edges, RefEdge);
  DLIList<RefFace*> ref_faces;
  CAST_LIST(tdu_list, ref_faces, RefFace);

    // do a little error checking: the entities should all be the same type, and
    // should be either ref edges or ref faces
  if ((ref_edges.size() > 0 && ref_faces.size() > 0) ||
      (ref_edges.size() > 0 && ref_edges.size() != tdu_list.size()) ||
      (ref_faces.size() > 0 && ref_faces.size() != tdu_list.size())) {
    PRINT_ERROR("Entities for composite containing %s %d not all the same.\n",
                attrib_owner()->class_name(), attrib_owner()->id());
    return CUBIT_FAILURE;
  }
  else if (ref_edges.size() == 0 && ref_faces.size() == 0) {
    PRINT_ERROR("Can't find any faces or edges for composite containing %s %d.\n",
                attrib_owner()->class_name(), attrib_owner()->id());
    return CUBIT_FAILURE;
  }

    // first, check to see that we can composite these entities; if
    // not, return without setting actuated flag; this CA will get 
    // re-actuated if a parent gets composited or partitioned.
  DLIList<BasicTopologyEntity*> bte_list;
  if (ref_edges.size() > 0) {
    CAST_LIST_TO_PARENT(ref_edges, bte_list);
  }
  else {
    CAST_LIST_TO_PARENT(ref_faces, bte_list);
  }
  
  if (!CompositeTool::instance()->okayToComposite(bte_list, NULL, NULL, CUBIT_FALSE))
  {
    for( int k = comp_attribs.size(); k--; )
      comp_attribs.get_and_step()->has_actuated( CUBIT_FALSE );
    has_actuated( CUBIT_FALSE );
    return CUBIT_FAILURE;
  }
    
// ok, we can build the composite
  RefEntity *new_entity = NULL;
  if (ref_edges.size() > 0) {
    PRINT_DEBUG_90("Creating composite edge from %d edges\n",
                   ref_edges.size());

    //need to get the RefEdge that owns this attribute
    RefEdge *edge_to_keep= NULL;
    edge_to_keep = CAST_TO( attrib_owner(), RefEdge );

    if( edge_to_keep == NULL )
      return CUBIT_FAILURE;

    new_entity = CompositeTool::instance()->composite(ref_edges, NULL, edge_to_keep);

    new TDUniqueId(new_entity, compositeId);
    CADeferredAttrib::owner_created(new_entity, compositeId);
  }
  else if (ref_faces.size() > 0) {
    PRINT_DEBUG_90("Creating composite face from %d faces\n",
                   ref_faces.size());

    //need to get the RefFace that owns this attribute
    RefFace *face_to_keep = NULL;
    face_to_keep = CAST_TO( attrib_owner(), RefFace );
    
    if( face_to_keep == NULL )
      return CUBIT_FAILURE;

    new_entity = CompositeTool::instance()->composite(ref_faces, face_to_keep);

    new TDUniqueId(new_entity, compositeId);
    CADeferredAttrib::owner_created(new_entity, compositeId);
  }
  else {
      // shouldn't get here
    assert(CUBIT_FALSE);
  }
/*
    // otherwise, we're done
    // set the actuated flag for all the CACVG's for this composite
  for (i = tdu_list.size(); i > 0; i--) {
    RefEntity *entity;
    if (ref_faces.size() > 0) entity = ref_faces.get_and_step();
    else entity = ref_edges.get_and_step();
    CACompositeVG *cacvg = (CACompositeVG *) entity->get_cubit_attrib(CA_COMPOSITE_VG);
    assert(cacvg != 0);
    cacvg->has_actuated(CUBIT_TRUE);
  }
*/
    // ok, we've composited; now check for CACVG's on any children, and
    // call actuate functions on those if they exist
  check_child_cacvgs(new_entity);

  return CUBIT_SUCCESS;
}

void CACompositeVG::check_child_cacvgs(RefEntity *new_entity) 
{
    // check the new entity's children for CACVG's, and actuate those
    // if they exist
  DLIList<RefEntity*> children;
  new_entity->get_child_ref_entities(children);

  DLIList<CubitAttrib*> ca_list;
  
    // get the CA list first, then actuate; if we looped over entities,
    // actuating as we went, some entities might eventually be gone
    // when we got to them
  int i;
  for (i = children.size(); i > 0; i--) {
    RefEntity *child = children.get_and_step();
    child->find_cubit_attrib_type(CA_COMPOSITE_VG, ca_list);
  }

  //get unactuated deferred attribs
  DLIList<CubitAttrib*> unact_deferred_attribs;
  unact_deferred_attribs = CADeferredAttrib::get_unactuated_deferred_attribs();


  for (i = ca_list.size(); i > 0; i--) 
  {
    CubitAttrib *ca_ptr = ca_list.get();
    if (ca_ptr->has_actuated() == CUBIT_TRUE)
      ca_list.remove();
    else
      ca_list.step();
  }

    // for the same reason, don't delete until we've gone through the list
  for (i = ca_list.size(); i > 0; i--) {
    CubitAttrib *ca_ptr = ca_list.get();

    //don't try and actuate unactuated deferred attribs.....wait till later
    if( unact_deferred_attribs.move_to( ca_ptr ) )
    {
      ca_list.remove();
      continue;
    }

    if (ca_ptr->delete_attrib() == CUBIT_FALSE && 
        ca_ptr->has_actuated() == CUBIT_FALSE) 
    {
      ca_ptr->actuate();
    }

      // remove undeletable attribs from the list so we don't delete them later
    if (ca_ptr->delete_attrib() == CUBIT_FALSE)
      ca_list.remove();
    else
      ca_list.step();
  }
  
  for (i = ca_list.size(); i > 0; i--)
    delete ca_list.get_and_step();
}


