// CAPartitionVG class

#include <vector>
#include <utility>
#include <algorithm>

#include "CAPartitionVG.hpp"
#include "TDUniqueId.hpp"
#include "CastTo.hpp"
#include "RefFace.hpp"
#include "RefEdge.hpp"
#include "RefVertex.hpp"
#include "PartitionEntity.hpp"
#include "PartitionTool.hpp"
#include "BasicTopologyEntity.hpp"
#include "PartitionCurve.hpp"
#include "PartitionSurface.hpp"
#include "CAVirtualVG.hpp"
#include "MergeTool.hpp"

struct my_sort : public std::binary_function< std::pair<double, int>, std::pair<double, int>, bool> {
  bool operator() (std::pair<double, int> const& x, std::pair<double, int> const & y ) { return x.first < y.first; }
};

CubitAttrib* CAPartitionVG_creator(RefEntity* entity, CubitSimpleAttrib *p_csa)
{
  CAPartitionVG *new_attrib = NULL;
  if (NULL == p_csa)
  {
    new_attrib = new CAPartitionVG(entity);
  }
  else
  {
    new_attrib = new CAPartitionVG(entity, p_csa);
  }

  return new_attrib;
}

CAPartitionVG::CAPartitionVG(RefEntity *owner) 
        : CubitAttrib(owner)
{
  numPC = 0;
  numPS = 0;
}

CAPartitionVG::CAPartitionVG(RefEntity *owner, CubitSimpleAttrib *simple_attrib) 
        : CubitAttrib(owner)
{
    // generate a simple attribute containing the data in this CA
  DLIList<CubitString*> *cs_list = simple_attrib->string_data_list();
  DLIList<double*> *d_list = simple_attrib->double_data_list();
  DLIList<int*> *i_list = simple_attrib->int_data_list();

  cs_list->reset();
  d_list->reset();
  i_list->reset();
  
    // (no string)

    // now the integers
    // numVP, numVC
  numPC = *(i_list->get_and_step());
  numPS = *(i_list->get_and_step());

    // numBdyCurves
  int temp, i, sum = 0;
  for (i = numPS; i > 0; i--) {
    temp = *(i_list->get_and_step());
    numBdyCurves.append(temp);
    sum += temp;
  }

    // vgUIDs: 3 for each PC, numPS+sum for PS
  for (i = 3*numPC+sum+numPS; i > 0; i--)
    vgUIDs.append(*(i_list->get_and_step()));


   // If the CubitSimpleAttrib already exists,
   // then this attribute is already written
  has_written(CUBIT_TRUE);
}

CubitStatus CAPartitionVG::update()
{
/*
    // this attribute behaves in a peculiar way: it detects whether the owner
    // is itself a partition entity, and if so, adds data about this partition
    // entity to the underlying entity

  if (hasUpdated) return CUBIT_SUCCESS;
  
  assert(attrib_owner() != 0);

  TopologyEntity *topo_entity = CAST_TO(attrib_owner(), TopologyEntity);
  assert(topo_entity != 0);
  DLIList<TopologyBridge*> bridge_list;
  topo_entity->bridge_manager()->get_bridge_list( bridge_list );
  for( int i = bridge_list.size(); i--; )
  {
    TopologyBridge *topo_bridge = bridge_list.get_and_step();
    PartitionEntity *partition_entity = CAST_TO(topo_bridge, PartitionEntity);

  if (partition_entity == NULL) {
      // this entity isn't a partition entity - if this entity doesn't have any virtual
      // entities registered, set delete flag, then exit
    if (numPC == 0 && numPS == 0)
      delete_attrib(CUBIT_TRUE);
    else {
      PRINT_INFO("Keeping CA_PARTITION_VG for %s %d\n",
                 attrib_owner()->class_name(), attrib_owner()->id());
      hasUpdated = CUBIT_TRUE;
    }
    
    continue;
  }
  
    // ok, we have a partition entity; first get the underlying entity, and a CAPVG 
    // for that entity
  BasicTopologyEntity* bte_ptr = partition_entity->get_underlying_BTE_ptr();
  if (!bte_ptr) {
    PRINT_ERROR("Couldn't find bound_to\n");
    return CUBIT_FAILURE;
  }
  
  CAPartitionVG *other_CAPVG = (CAPartitionVG *) bte_ptr->get_cubit_attrib(CA_PARTITION_VG);

    // if that other CAPVG's written flag is set, it's an old one from a
    // previous write and needs to be reset
  if (other_CAPVG->has_written() == CUBIT_TRUE) {
    other_CAPVG->reset();
    other_CAPVG->has_written(CUBIT_FALSE);
  }

    // now put virtual geometry-specific data on the attribute
  PartitionCurve *partition_curve = CAST_TO(partition_entity, PartitionCurve);
  PartitionSurface *partition_surface = CAST_TO(partition_entity, PartitionSurface);
  
  if (partition_curve != NULL) {
    other_CAPVG->add_pcurve(partition_curve);
    other_CAPVG->delete_attrib(CUBIT_FALSE);
  }

  else if (partition_surface != NULL) {
    other_CAPVG->add_psurface(partition_surface);
    other_CAPVG->delete_attrib(CUBIT_FALSE);
  }

  else {
    PRINT_ERROR("Shouldn't get here in CAPartitionVG::update.\n");
    return CUBIT_FAILURE;
  }
  }

  hasUpdated = CUBIT_TRUE;
  if (numPC == 0 && numPS == 0) delete_attrib(CUBIT_TRUE);
  
  return CUBIT_SUCCESS;
*/ 
  delete_attrib(CUBIT_TRUE);
  return CUBIT_SUCCESS;
}


CubitStatus CAPartitionVG::reset()
{
  numPC = 0;
  numPS = 0;
  
  vgUIDs.clean_out();
  numBdyCurves.clean_out();
  
  return CUBIT_SUCCESS;
}

CubitSimpleAttrib *CAPartitionVG::cubit_simple_attrib()
{
    // generate a simple attribute containing the data in this CA
  DLIList<CubitString*> cs_list;
  DLIList<double> d_list;
  DLIList<int> i_list;

    // first the string
  cs_list.append(new CubitString(att_internal_name()));

    // now the integers
    // numVP, numVC
  i_list.append(numPC);
  i_list.append(numPS);

    // numBdyCurves
  int i;
  for (i = numBdyCurves.size(); i > 0; i--)
    i_list.append(numBdyCurves.get_and_step());

    // vgUIDs
  vgUIDs.reset();
  for (i = vgUIDs.size(); i > 0; i--)
    i_list.append(vgUIDs.get_and_step());

  CubitSimpleAttrib* csattrib_ptr = new CubitSimpleAttrib(&cs_list,
                                                          &d_list,
                                                          &i_list);

  for( i=cs_list.size(); i--;) delete cs_list.get_and_step();

  return csattrib_ptr;
}

CubitStatus CAPartitionVG::actuate()
{
    // actuate this CA

    // actuate partition VG attributes on next-lower order entities
  RefEntity *owner = CAST_TO(attrib_owner(), RefEntity);
  if (owner->dimension() > 1) {
    DLIList<RefEntity*> lower_entities;
    owner->get_child_ref_entities(lower_entities);

      // don't check return values here - there may be other CA's hanging
      // around unactuated, but we may still be able to actuate later
    CubitAttribUser::actuate_cubit_attrib(lower_entities, CA_PARTITION_VG);
    CubitAttribUser::actuate_cubit_attrib(lower_entities, CA_VIRTUAL_VG);
  }

  //actuate it

    // if this is an edge, now partition it
  RefEdge *owner_edge = CAST_TO(owner, RefEdge);
  RefFace *owner_face = CAST_TO(owner, RefFace);

  //have to get some of the data from CAVirtualVG in order to partition
  //get CA_VIRTUAL_VG attrib associated with this entity
  DLIList<CubitAttrib*> vg_attribs;
  attrib_owner()->find_cubit_attrib_type( CA_VIRTUAL_VG, vg_attribs );
  CAVirtualVG *ca_vg_ptr = NULL;
  if( vg_attribs.size() != 0 )
    ca_vg_ptr = CAST_TO( vg_attribs.get(), CAVirtualVG );

  DLIList<RefFace*> new_faces;
  DLIList<RefEdge*> new_edges;

  int i,j,k;
  if (owner_edge) //if an edge has been partitioned 
  {
    DLIList<RefVertex*> new_vertices;
    DLIList<RefEdge*> new_edges;
    DLIList<CubitVector*> split_points;
      
    if( !ca_vg_ptr ) //if NO virtual geometry has been used to partiton this curve
    {
      //get vertices of curve
      DLIList<RefEntity*> lower_entities;
      owner->get_child_ref_entities(lower_entities);

      //get the points that are not owned by a vertex  
      vgUIDs.reset();
      DLIList<RefVertex*> vertices;
      for(i=vgUIDs.size()/3; i--;)
      {
        ToolDataUser *tdu = TDUniqueId::find_td_unique_id(vgUIDs.get_and_step());
        RefVertex *s_vert = CAST_TO( tdu, RefVertex ); 
        tdu = TDUniqueId::find_td_unique_id(vgUIDs.get_and_step());
        RefVertex *e_vert = CAST_TO( tdu, RefVertex ); 
      
        //get vertices owned by this curve 
        if( !lower_entities.move_to( s_vert ) )
          vertices.append_unique( s_vert );
        if( !lower_entities.move_to( e_vert ) )
          vertices.append_unique( e_vert );

        vgUIDs.step();
      }
      //convert vertices to vectors to split curve
      for(i=vertices.size(); i--;)
      {
        RefVertex *cur_vertex = vertices.get_and_step();
        split_points.append( new CubitVector( cur_vertex->coordinates() ) );
      }

      //partition the curve with these split points 
      split_points.reset(); 
      PartitionTool::instance()->partition( owner_edge, split_points,
                                            new_vertices, new_edges );
      //may need to merge some vertices
      vertices += new_vertices;
      MergeTool::instance()->merge_refvertices( vertices );

    }
    else //virtual geometry HAS been used to partiton this curve
    {
      split_points = ca_vg_ptr->posVector;
      DLIList<int>vertex_unique_ids = ca_vg_ptr->vgUIDs;
      std::vector< std::pair<double, int> > list_of_pairs;

      //before partitioning edge, reorder split points, from lowest u to 
      //highest u; we use a pair so that the corresponding uids are reordered as well 
      for( i=split_points.size(); i--;)
      {
        CubitVector *split_point = split_points.get_and_step();
        double u_param = owner_edge->u_from_position( *split_point );
        int uuid = vertex_unique_ids.get_and_step();
        std::pair<double, int> my_pair;
        my_pair.first = u_param;
        my_pair.second = uuid;
        list_of_pairs.push_back( my_pair );
      }
      
      std::sort(list_of_pairs.begin(), list_of_pairs.end(), my_sort() ); 

      //partition the curve with these split points 
      split_points.reset(); 
      PartitionTool::instance()->partition( owner_edge, split_points,
                                            new_vertices, new_edges );
         
      assert( list_of_pairs.size() == (unsigned int)new_vertices.size() );
      //associate the vertex uuids with the new vertices 
      new_vertices.reset();
      vertex_unique_ids.reset();
      std::vector< std::pair<double, int> >::iterator iter = list_of_pairs.begin();
      for( i=new_vertices.size(); i--; )
      {
        new TDUniqueId( new_vertices.get_and_step(), (*iter).second );
        iter++;
      }
    }

    //associate the curve uuids with the new curve 
    new_edges.reset();
    RefEdge* ref_edge;
    vgUIDs.reset();
    for( i=new_edges.size(); i--; )
    {
      ref_edge = new_edges.get_and_step();
      //get uuids of start and end vertices on curve
      RefVertex *start_vertex = ref_edge->start_vertex();
      RefVertex *end_vertex = ref_edge->end_vertex();

      int s_vert_uuid = TDUniqueId::get_unique_id( start_vertex );
      int e_vert_uuid = TDUniqueId::get_unique_id( end_vertex );
      
      for( j=vgUIDs.size(); j--;)
      { 
        int s_uuid = vgUIDs.get_and_step();
        int e_uuid = vgUIDs.get_and_step();
        if( (s_vert_uuid == s_uuid && e_vert_uuid == e_uuid ) ||
            (e_vert_uuid == s_uuid && s_vert_uuid == e_uuid ) )
        {
         new TDUniqueId( ref_edge, vgUIDs.get_and_step() );
         break;
        }
        else
          vgUIDs.step();
      }
    }
  }
  else if (owner_face)  //partition a surface
  {
    DLIList<RefFace*> input_faces;
    input_faces.append( owner_face );
    DLIList<CubitVector*> segments;

    //for each virtual curve used to partition surface
    ca_vg_ptr->numVCPoints.reset();
    ca_vg_ptr->posVector.reset();
    ca_vg_ptr->vgUIDs.step( ca_vg_ptr->numVV );
    for( i=ca_vg_ptr->numVC; i--;)
    {
      //get coordinates of start/end vertices of virtual curve
      ToolDataUser *tdu = TDUniqueId::find_td_unique_id(ca_vg_ptr->vgUIDs.get_and_step());
      RefVertex *s_vert = CAST_TO( tdu, RefVertex ); 
      tdu = TDUniqueId::find_td_unique_id(ca_vg_ptr->vgUIDs.get_and_step());
      RefVertex *e_vert = CAST_TO( tdu, RefVertex ); 
      DLIList<RefVertex*> vertices;

      CubitVector *vec1, *vec2;
      if( s_vert)
      {
        vec1 = new CubitVector(s_vert->coordinates());
        vertices.append( s_vert );
      }

      else
        vec1 = new CubitVector( *(ca_vg_ptr->posVector.get() ));
      
      segments.append( vec1 );

      if( e_vert)
      {
        vec2 = new CubitVector(e_vert->coordinates());
        vertices.append( e_vert );
      }
      else
        vec2 = new CubitVector( *(ca_vg_ptr->posVector.get() ));

      //append any intermediate segments
      for( j=ca_vg_ptr->numVCPoints.get_and_step(); j--; ) 
        segments.append( ca_vg_ptr->posVector.get_and_step() );

      segments.append( vec2 );
      //partition the surf
      new_edges.clean_out();
      new_faces.clean_out();
      PartitionTool::instance()->insert_edge( input_faces, segments,
                                              new_faces, new_edges);

      //may need to merge some vertices
      DLIList<RefVertex*> verts_to_merge;
      for( j=new_edges.size(); j--;)
      {
        verts_to_merge.append( new_edges.get()->start_vertex() );
        verts_to_merge.append( new_edges.get()->end_vertex() );
      }
      verts_to_merge += vertices;
      verts_to_merge.uniquify_unordered();

      MergeTool::instance()->merge_refvertices( verts_to_merge  );

      //give new edge uuids
      new_edges.reset();
      for( j=new_edges.size(); j--;)
        new TDUniqueId( new_edges.get_and_step(), ca_vg_ptr->vgUIDs.get_and_step() );
      
      delete vec1;
      delete vec2;
      segments.clean_out();
    
    
    //associate the face uuids with the new faces
    new_faces.reset();
    RefFace* ref_face;
    vgUIDs.reset();
    for( j=new_faces.size(); j--; )
    {
      ref_face = new_faces.get_and_step();

      //get uuids of each edge in surface   
      DLIList<RefEdge*> edges;
      ref_face->ref_edges( edges );
      DLIList<int> edge_uuids;

      edges.reset();
      for( k=edges.size(); k--;)
        edge_uuids.append( TDUniqueId::get_unique_id( edges.get_and_step() ));

      //look at all the groups of boundary edges of each surface
      vgUIDs.reset();
      for( k=numBdyCurves.size(); k--;)
      {
        int kk;
        DLIList<int> bdy_curve_uuids;
        int num_bdy_curves = numBdyCurves.get_and_step();
        if( edge_uuids.size() != num_bdy_curves )
          continue;
        for( kk=num_bdy_curves; kk--;)
          bdy_curve_uuids.append( vgUIDs.get_and_step());
              
        CubitBoolean match = CUBIT_TRUE;
        kk=num_bdy_curves;
        while( kk && match )
        {
          kk--;
          if( !bdy_curve_uuids.move_to( edge_uuids.get_and_step() ) )
          {
            match = CUBIT_FALSE; 
            break;
          }
        }
        if( match )
        {
          new TDUniqueId( ref_face, vgUIDs.get() );
          break;
        }
        vgUIDs.step();
      }
    }

    input_faces.clean_out();
    input_faces += new_faces;

    }
  }
 
  hasActuated = CUBIT_TRUE;
 
   // otherwise, we're done
  return CUBIT_SUCCESS;
}



