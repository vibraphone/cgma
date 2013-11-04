// CAVirtualVG class

#include "CAVirtualVG.hpp"
#include "CubitSimpleAttrib.hpp"
#include "TopologyEntity.hpp"
#include "RefEntity.hpp"
#include "BasicTopologyEntity.hpp"
#include "RefEdge.hpp"
#include "RefVertex.hpp"
#include "TDUniqueId.hpp"
#include "CADeferredAttrib.hpp"

CubitAttrib* CAVirtualVG_creator(RefEntity* entity, const CubitSimpleAttrib &p_csa)
{
  return new CAVirtualVG(entity, p_csa);
}


CAVirtualVG::CAVirtualVG(RefEntity *owner, const CubitSimpleAttrib &simple_attrib)
        : CubitAttrib(owner)
{
  numVC = 0;
  numVV = 0;

  if(!simple_attrib.isEmpty())
  {
     // generate a simple attribute containing the data in this CA
    const std::vector<double> &d_list = simple_attrib.double_data_list();
    const std::vector<int> &i_list = simple_attrib.int_data_list();

      // (no string)

      // now the integers
      // numVV, numVC
    int offset = 0;
    numVV = i_list[offset++];
    numVC = i_list[offset++];

      // vgUIDs
    int i;
    for (i = numVV+(3*numVC); i > 0; i--)
      vgUIDs.append(i_list[offset++]);

      // numVCPoints
    int sum = 0;
    for (i = numVC; i > 0; i--) {
      int temp_int = i_list[offset++];
      numVCPoints.append(temp_int);
      sum += temp_int;
    }

    offset = 0;
      // pointUIDList
    for (i = numVV+sum; i > 0; i--) {
      double x = d_list[offset++];
      double y = d_list[offset++];
      double z = d_list[offset++];
      posVector.append(new CubitVector( x, y, z ) );
    }

      // If we are constructing from a CubitSimpleAttrib,
      // then this attrib is already written
    has_written(CUBIT_TRUE);
  }
}

CubitStatus CAVirtualVG::update()
{

    // the owner entity is virtual; put a CAVirtualVG on the underlying entity to represent
    // this entity

  if (hasUpdated) return CUBIT_SUCCESS;
/*  
  assert(attrib_owner() != 0);

  TopologyEntity *topo_entity = CAST_TO(attrib_owner(), TopologyEntity);
  assert(topo_entity != 0);
  
  DLIList<VirtualEntity*> ve_list;
  VirtualGeometryEngine::instance()->get_VEs(topo_entity, ve_list, CUBIT_FALSE, CUBIT_FALSE);
  
  for( int i = ve_list.size(); i--; )
  {
    ParasiteEntity *vge = dynamic_cast<ParasiteEntity*>(ve_list.get_and_step());

  if (vge == NULL) {
      // this entity isn't a virtual entity - if this entity doesn't have any virtual
      // entities registered, set delete flag, then exit
    if (numVV == 0 && numVC == 0)
      delete_attrib(CUBIT_TRUE);
    else {
      PRINT_DEBUG_90("Keeping CA_VIRTUAL_VG for %s %d\n",
                     attrib_owner()->class_name(), attrib_owner()->id());
      hasUpdated = CUBIT_TRUE;
    }
    
    continue;
  }

    // ok, we have a virtual entity; first get the underlying entity, and a CAVVG for that
    // entity
  BasicTopologyEntity* bte_ptr = vge->bte_bound_to();
  
  if (!bte_ptr) {
    PRINT_ERROR("Couldn't find underlying BTE for virtual entity %s %d.\n",
                attrib_owner()->class_name(), attrib_owner()->id());
    return CUBIT_FAILURE;
  }
  
  CAVirtualVG *other_CAVVG = (CAVirtualVG *) bte_ptr->get_cubit_attrib(CA_VIRTUAL_VG);

    // if that other CAPVG's written flag is set, it's an old one from a
    // previous write and needs to be reset
  if (other_CAVVG->has_written() == CUBIT_TRUE) {
    other_CAVVG->reset();
    other_CAVVG->has_written(CUBIT_FALSE);
  }

    // now put virtual geometry-specific data on the attribute
  ParasitePoint *vvertex = CAST_TO(vge, ParasitePoint);
  ParasiteCurve *vcurve = CAST_TO(vge, ParasiteCurve);
  if (vvertex != NULL) {
      // put uids and position on ca
    other_CAVVG->add_vvertex(vvertex);
    other_CAVVG->delete_attrib(CUBIT_FALSE);
  }

  else if (vcurve != NULL) {
    other_CAVVG->add_vcurve(vcurve);
    other_CAVVG->delete_attrib(CUBIT_FALSE);
  }

  else {
    PRINT_ERROR("Shouldn't get here in CAVirtualVG::update.\n");
    return CUBIT_FAILURE;
  }
  }

  hasUpdated = CUBIT_TRUE;
  if (numVV == 0 && numVV == 0) delete_attrib(CUBIT_TRUE);
  
  return CUBIT_SUCCESS;
*/ 

  delete_attrib(CUBIT_TRUE);
  return CUBIT_SUCCESS;
}


CubitStatus CAVirtualVG::reset()
{
    //- reset info; called from CAU and also from update!
  numVV = 0;
  numVC = 0;
    //- the number of virtual points and curves contained in this attr

  vgUIDs.clean_out();
    //- unique ids of virtual points and curves contained in this attr

  posVector.clean_out();
    //- position vectors for virtual curves

  numVCPoints.clean_out();
    //- for each virtual curve, the number of virtual points in that curve

  return CUBIT_SUCCESS;
}

CubitSimpleAttrib CAVirtualVG::cubit_simple_attrib()
{
    // generate a simple attribute containing the data in this CA
  std::vector<CubitString> cs_list;
  std::vector<double> d_list;
  std::vector<int> i_list;

    // first the string
  cs_list.push_back(att_internal_name());

    // now the integers
    // numVV, numVC
  i_list.push_back(numVV);
  i_list.push_back(numVC);

    // vgUIDs
  int i;
  vgUIDs.reset();
  for (i = vgUIDs.size(); i > 0; i--)
    i_list.push_back(vgUIDs.get_and_step());

    // numVCPoints
  numVCPoints.reset();
  for (i = numVCPoints.size(); i > 0; i--)
    i_list.push_back(numVCPoints.get_and_step());
  
    // now the doubles
  posVector.reset();
  for (i = posVector.size(); i > 0; i--) {
    d_list.push_back(posVector.get()->x());
    d_list.push_back(posVector.get()->y());
    d_list.push_back(posVector.get_and_step()->z());
  }
  
  return CubitSimpleAttrib(&cs_list, &d_list, &i_list);
}

CubitStatus CAVirtualVG::actuate()
{
    // actuate this CA
/*
    // first, actuate the points, with this CA's owner as their owners
  int i, j;
  vgUIDs.reset();
  posVector.reset();

  DLIList<int> leftover_uids, leftover_vcpoints;
  DLIList<CubitVector*> leftover_posvector;
  int leftover_vvs = 0, leftover_vcs = 0;

  DLIList<ParasitePoint*> vpoint_list;
  for (i = 0; i < numVV; i++) {
    CubitVector pos_vec(*posVector.get_and_step());
      // make the virtual vertex
    BasicTopologyEntity *bte_ptr = CAST_TO(attrib_owner(), BasicTopologyEntity);
    RefVertex *new_vertex = 
      VirtualGeometryEngine::instance()->create_VirtualVertex(pos_vec, bte_ptr);
      // add a unique id to it, for referencing later
    TDUniqueId *uid_ptr = new TDUniqueId(new_vertex, vgUIDs.get());
    CADeferredAttrib::owner_created( new_vertex, vgUIDs.get_and_step() );

    if (new_vertex == NULL) {
      PRINT_WARNING("Problems making new vertex with uid = %d.\n", vgUIDs.prev());
      leftover_uids.append(vgUIDs.prev());
      leftover_posvector.append(posVector.prev());
      leftover_vvs++;
    }
  }

    // now actuate CAVVG's on lower order entities
  DLIList<RefEntity*> children;
  attrib_owner()->get_child_ref_entities(children);
  for (i = children.size(); i > 0; i--) {
    children.get_and_step()->actuate_cubit_attrib(CA_VIRTUAL_VG);
  }

    // now handle creating virtual curves
  numVCPoints.reset();
  DLIList<RefEdge*> vcurve_list;
  for (i = 0; i < numVC; i++) {
      // start by grabbing all the stuff off the attribute's lists

    int start_uid = vgUIDs.get_and_step(),
      end_uid = vgUIDs.get_and_step();

      // now get the intermediate points; these should all be virtual, and should
      // be in this entity's virtual point list
    DLIList<CubitVector*> vec_list;
      // get the number of virtual points in the list
    int num_points = numVCPoints.get_and_step();
    for (j = 0; j < num_points; j++)
      vec_list.append(posVector.get_and_step());

    int virtual_curve_uid = vgUIDs.get_and_step();
    
      // the first two are start and end points, and may not be virtual
    ToolDataUser *tdu = TDUniqueId::find_td_unique_id(start_uid);
    RefVertex *start_vertex = CAST_TO(tdu, RefVertex);
    tdu = TDUniqueId::find_td_unique_id(end_uid);
    RefVertex *end_vertex = CAST_TO(tdu, RefVertex);
    
    if (!start_vertex || !end_vertex) {
      PRINT_DEBUG_90("Couldn't restore virtual curve with uid = %d.\n",
                     virtual_curve_uid);
        // cache leftover data for restoring later
      leftover_uids.append(start_uid);
      leftover_uids.append(end_uid);
      leftover_uids.append(virtual_curve_uid);
      leftover_vcpoints.append(num_points);
      leftover_posvector += vec_list;
      leftover_vcs++;
      continue;
    }
    
      // make the virtual curve 
    RefEdge *virtual_edge = 
      VirtualGeometryEngine::instance()->create_VirtualEdge(start_vertex,
                                                            end_vertex,
                                                            vec_list);
    if (!virtual_edge) {
      PRINT_DEBUG_90("Couldn't restore virtual curve with uid = %d.\n",
                     virtual_curve_uid);
        // cache leftover data for restoring later
      leftover_uids.append(start_uid);
      leftover_uids.append(end_uid);
      leftover_uids.append(virtual_curve_uid);
      leftover_vcpoints.append(num_points);
      leftover_posvector += vec_list;
      leftover_vcs++;
      continue;
    }

    ParasiteEntity* curve = 
      dynamic_cast<ParasiteEntity*>(virtual_edge->get_geometry_entity_ptr());
    assert(curve != NULL);
    curve->bind_to( dynamic_cast<BasicTopologyEntity*>(attrib_owner())
      ->get_geometry_entity_ptr());

      // save the curve's unique id
    TDUniqueId *uid_ptr = new TDUniqueId(virtual_edge, virtual_curve_uid);
    
    CADeferredAttrib::owner_created( virtual_edge, virtual_curve_uid );
    //virtual_edge->actuate_cubit_attrib(CA_VIRTUAL_VG);
      //virtual_edge->actuate_cubit_attrib(CA_PARTITION_VG);
  }

  if (0 == leftover_vvs && 0 == leftover_vcs)
    hasActuated = CUBIT_TRUE;
  else {
      // have some leftover data - reset data in this attribute
    numVV = leftover_vvs;
    numVC = leftover_vcs;
    vgUIDs = leftover_uids;
    numVCPoints = leftover_vcpoints;
    posVector = leftover_posvector;
    hasActuated = CUBIT_FALSE;

      // now add this attribute to the unactuated list in CADA
    CADeferredAttrib::add_unactuated_ca(this);
  }

   // we're done
  return (CUBIT_FALSE == hasActuated ?
          CUBIT_FAILURE : CUBIT_SUCCESS);
*/ return CUBIT_FAILURE;
}

/*
void CAVirtualVG::add_vvertex(ParasitePoint *vpoint)
{

  TopologyEntity *vpoint_topo_entity = vpoint->topology_entity();
  RefVertex *owning_vertex = CAST_TO(vpoint_topo_entity, RefVertex);
  assert(owning_vertex != 0);
  int vvertex_uid = TDUniqueId::get_unique_id(owning_vertex);
  
    //- adds data for this vpoint to this CA
  numVV++;
  vgUIDs.insert_first(vvertex_uid);
  posVector.insert_first(new CubitVector(vpoint->coordinates()));
}
*/
/*
void CAVirtualVG::add_vcurve(ParasiteCurve *vcurve)
{
    // need to get list of vpoint uids defining the virtual curve
    // the owner should be a RefEdge
  TopologyEntity *vcurve_topo_entity = vcurve->topology_entity();
  RefEdge *virtual_edge = CAST_TO(vcurve_topo_entity, RefEdge);
  assert(virtual_edge != 0);

    // first get start and end points, and their uids
  DLIList<int> uid_list;
  RefVertex *temp_vertex = virtual_edge->start_vertex();
  int start_uid = TDUniqueId::get_unique_id(temp_vertex);
  vgUIDs.append(start_uid);
  temp_vertex = virtual_edge->end_vertex();
  int end_uid = TDUniqueId::get_unique_id(temp_vertex);
  vgUIDs.append(end_uid);

    // now get the other points and their uids
  DLIList<ParasitePoint*> vpoint_list;
  vcurve->getSegments(vpoint_list);
  vpoint_list.reset();

  for (int i = vpoint_list.size(); i > 0; i--) {
    ParasitePoint *temp_point = vpoint_list.get_and_step();
    posVector.append(new CubitVector(temp_point->coordinates()));
  }

  int owner_uid = TDUniqueId::get_unique_id(virtual_edge);
  
    //- adds data for this vcurve to this CA
  vgUIDs.append(owner_uid);
  numVCPoints.append(vpoint_list.size());
  numVC++;
}

*/

