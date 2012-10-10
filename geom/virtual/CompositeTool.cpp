#include <assert.h>
#include <vector>

#include "CompositeTool.hpp"

#include "DLIList.hpp"
#include "CastTo.hpp"
#include "CubitUtil.hpp"

#include "RefEntityFactory.hpp"
#include "RefEdge.hpp"
#include "RefFace.hpp"
#include "RefVertex.hpp"
#include "RefVolume.hpp"
#include "Loop.hpp"
#include "Shell.hpp"
#include "Chain.hpp"
#include "CoEdge.hpp"
#include "CoFace.hpp"
#include "CoVertex.hpp"
#include "Body.hpp"

#include "CompositeCurve.hpp"
#include "CompositeSurface.hpp"
#include "CompositeLump.hpp"
#include "CompositeBody.hpp"

//#include "VirtualQueryEngine.hpp"
#include "PartitionTool.hpp"
#include "BasicTopologyEntity.hpp"
#include "ModelQueryEngine.hpp"
#include "RefEntityName.hpp"

#include "CompositeEngine.hpp"
#include "Point.hpp"
#include "GeometryQueryTool.hpp"
#include "DAG.hpp"
//#include "CAMergePartner.hpp"
//#include "TDUniqueId.hpp"
#include "CubitAttrib.hpp"
#include "CADefines.hpp"

#include "PartitionEngine.hpp"
#include "PartitionSurface.hpp"
#include "SegmentedCurve.hpp"
#include "PartPTCurve.hpp"
#include "PartitionPoint.hpp"
#include "CompositeCombineEvent.hpp"
#include "GMem.hpp"
#include "GfxDebug.hpp"

CompositeTool* CompositeTool::instance_ = NULL;


//-------------------------------------------------------------------------
// Purpose       : Constructor 
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 10/17/03
//-------------------------------------------------------------------------
CompositeTool::CompositeTool()
{
}


//-------------------------------------------------------------------------
// Purpose       : Destructor
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 
//-------------------------------------------------------------------------
CompositeTool::~CompositeTool()
{
  assert( instance_ == this );
  instance_ = NULL;
}

//-------------------------------------------------------------------------
// Purpose       : Create a composite curve
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 
//-------------------------------------------------------------------------
RefEdge* CompositeTool::composite( DLIList<RefEdge*>& edge_list,
                                   RefVertex* keep,
                                   RefEdge* survivor  )
{
  int i;
  DLIList<RefVertex*> vertex_list;
  DLIList<RefEdge*> modified_edges;
  
  for( i = edge_list.size(); i--; )
  {
    RefEdge* edge = edge_list.get_and_step();
    RefVertex* start_vtx = edge->start_vertex();
    RefVertex* end_vtx = edge->end_vertex();
    start_vtx->marked(0);
    end_vtx->marked(0);
    vertex_list.append( start_vtx );
    vertex_list.append( end_vtx );
  }
  for( i = vertex_list.size(); i--; )
  {
    RefVertex* vtx = vertex_list.get_and_step();
    vtx->marked( vtx->marked() + 1 );
  }
  // Check for edges that will have both vertices removed.
  // If any exist then we want to set the flag to have the
  // DAG updated after each vertex is removed.  Otherwise
  // things get left in a bad state.
  bool update_dag = false;
  for(i=edge_list.size(); i>0; i--)
  {
    RefEdge* edge = edge_list.get_and_step();
    if(edge->start_vertex()->marked() == 2 &&
       edge->end_vertex()->marked() == 2)
    {
      update_dag = true;
      i=0;
    }
  }
  for( i = vertex_list.size(); i--; )
    if( vertex_list.step_and_get()->marked() == 2 )
      vertex_list.get()->marked(0);
    else
      vertex_list.change_to( 0 );
  vertex_list.remove_all_with_value(0);
  
  if ( keep && vertex_list.move_to(keep) )
    vertex_list.remove();
  
  if( vertex_list.size() == 0 )
  {
    PRINT_ERROR("Cannot create any composites from passed curves.\n");
    return 0;
  }
  vertex_list.reset();
  RefEdge* result = 0;
  DLIList<RefEdge*> vtx_edges;
  for( i = vertex_list.size(); i--; )
  {
    RefVertex* vtx = vertex_list.get_and_step();
    vtx_edges.clean_out();
    vtx->ref_edges(vtx_edges);
    //check to make sure that the two edges resides on the same surface(s)
    if ( vtx_edges.size() == 2 )
    {
      RefEdge* first_edge = vtx_edges.get();
      RefEdge* second_edge = vtx_edges.step_and_get();
      DLIList <RefFace*> common_faces;
      int size = first_edge->common_ref_faces(second_edge, common_faces);
      if (first_edge->num_ref_faces() != second_edge->num_ref_faces() ||
          size != first_edge->num_ref_faces())
      {
         PRINT_INFO("Curve %d and curve %d don't reside on the same surfaces.\n", first_edge->id(), second_edge->id());
         continue;
      }
    
      result = remove_vertex( vtx, false, update_dag, survivor );
      if (result)
      {
        // First remove any edges from the list that got destroyed in the
        // last remove_vertex operation.
        if(result == first_edge && modified_edges.move_to(second_edge))
        {
          modified_edges.remove();
        }
        else if(result == second_edge && modified_edges.move_to(first_edge))
        {
          modified_edges.remove();
        }
        modified_edges.append_unique(result);
      }
    }  
  }
  
  for ( i = modified_edges.size(); i--; )
    if ( ! modified_edges.step_and_get()->get_curve_ptr() )
      modified_edges.change_to(0);
  modified_edges.remove_all_with_value(0);
  
  result = modified_edges.size() ? modified_edges.get() : 0;

  DLIList<Surface*> update_surfaces, curve_surfaces;
  DLIList<TopologyBridge*> curve_bridges;
  for ( i = modified_edges.size(); i--; )
  {
    RefEdge* edge = modified_edges.get_and_step();
    edge->set_id( RefEntityFactory::instance()->next_ref_edge_id() );
    edge->bridge_manager()->get_bridge_list(curve_bridges);
    curve_bridges.reset();
    for (int j = curve_bridges.size(); j--; )
    {
      curve_surfaces.clean_out();
      curve_bridges.get_and_step()->surfaces(curve_surfaces);
      update_surfaces.merge_unique(curve_surfaces);
    }
  }

  update_surfaces.reset();
  for ( i = update_surfaces.size(); i--; )
  {
    GeometryQueryTool::instance()->make_RefFace(update_surfaces.get_and_step());
  }
  assert( !(survivor && survivor->deactivated()) );
  DAG::instance()->cleanout_deactivated_DAG_nodes();

  return result;
}


//-------------------------------------------------------------------------
// Purpose       : Composite over a vertex
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 03/12/03
//-------------------------------------------------------------------------
RefEdge* CompositeTool::remove_vertex( RefVertex* vertex,
                                       bool remove_partitions, /* = false */
                                       bool update_dag ,       /* = true  */
                                       RefEdge* keep_edge      /* = NULL */ )
{
  DLIList<RefEdge*> vertex_edges;
  int i;
  
    // Get two ref-edges
  vertex->ref_edges( vertex_edges );
  if ( vertex_edges.size() != 2 )
  {
    PRINT_ERROR("Cannot composite over valence-%d vertex.\n", vertex_edges.size());
    return 0;
  }
  
  DLIList<TopologyBridge*> point_bridges;
  DLIList<CoEdgeSM*> coedgesms;
    
  RefEdge *edge = 0, *refedge1 = 0, *refedge2 = 0;
  point_bridges.clean_out();
  vertex->bridge_manager()->get_bridge_list( point_bridges );
    
  std::vector<CubitString> names_to_add;
  Curve* end_result = 0;
  point_bridges.reset();
  while( point_bridges.size() )
  {
    Curve* result_curve = 0;
    TBPoint* point = dynamic_cast<TBPoint*>(point_bridges.pop());

      // get the RefEdges that will be joined when the vertex is removed
    DLIList<Curve*> curves;
    point->curves(curves);
    curves.reset();
    
    if (curves.size() == 0)  // merged free vertex -- ignore it
      continue;

      // If only one curve, need to stitch point with some other
      // point in same vertex.
    if (curves.size() == 1) {
    
      /* Stitching code not complete.  Fail for now. */
      return 0;
      
      /*
      Curve* first_curve = curves.pop();
    
      TBPoint* other = 0;
      for ( i = point_bridges.size(); i--; ) {
        TBPoint* temp = dynamic_cast<TBPoint*>(point_bridges.step_and_get());
        curves.clean_out();
        temp->curves(curves);
        if ( curves.size() == 1 && curves.get()->owner() != first_curve->owner() ) 
        {
          point_bridges.extract();
          other = temp;
          break;
        }
      }
    
      if (!other)
        return 0;
      
      point = CompositeEngine::instance().stitch_points(point, other);
      if (!point)
        return 0;
        
      curves.append(first_curve);
      
      */
    } // end if (curves.size() == 1)
    

   // assert(curves.size() == 2);
    if(curves.size() != 2)
    {
      break;
    }


    TopologyEntity* topo = curves.get()->topology_entity();
    refedge1 = CAST_TO(topo, RefEdge);
    assert(refedge1 != NULL);
    curves.get()->set_saved_id(refedge1->id());
    DLIList<CubitString*> tmp_names;
    std::vector<CubitString> edge1_names;
    refedge1->entity_names( tmp_names );
    curves.get()->set_saved_names( tmp_names );  
    for( int i=tmp_names.size(); i--; )
      edge1_names.push_back( *(tmp_names.get_and_step()) );
    tmp_names.clean_out();
    refedge1->remove_entity_names();

    topo = curves.next()->topology_entity();
    refedge2 = CAST_TO(topo, RefEdge);
    assert(refedge2 != NULL);   
    curves.next()->set_saved_id(refedge2->id());
    std::vector<CubitString> edge2_names;
    refedge2->entity_names( tmp_names );
    curves.next()->set_saved_names( tmp_names );    
    for( int i=tmp_names.size(); i--; )
      edge2_names.push_back( *(tmp_names.get_and_step()) );   
    refedge2->remove_entity_names();

    Curve* keep = 0;
    if (keep_edge) {
      if (keep_edge->bridge_manager() == curves.get()->owner())
        keep = curves.get();
      else if (keep_edge->bridge_manager() == curves.next()->owner())
        keep = curves.next();
    }

      // remove the point to make the composite curve
    PartitionPoint* part = dynamic_cast<PartitionPoint*>(point);
    if( remove_partitions && part )
      result_curve = PartitionEngine::instance().remove_point( part );
    if ( !result_curve )
      result_curve = CompositeEngine::instance().remove_point( point, remove_partitions, keep );

    if( !result_curve )
    {
      //reapply names
      for( int i=edge1_names.size(); i--; )
        refedge1->entity_name( edge1_names[i] );
      for( int i=edge2_names.size(); i--; )
        refedge2->entity_name( edge2_names[i] );

      PRINT_ERROR("Failed to remove vertex %d\n", vertex->id() );
      break;
    }
    else if( !edge )
    {
      edge = dynamic_cast<RefEdge*>(result_curve->topology_entity());     
    }
    else
    {
      if( result_curve->owner() )
        result_curve->owner()->remove_bridge( result_curve );
      edge->bridge_manager()->add_bridge( result_curve );
        
      Curve* curve = edge->get_curve_ptr();
      bool geom_sense = curve->relative_sense( result_curve ) == CUBIT_FORWARD;
      bool bridge_sense = curve->bridge_sense() == result_curve->bridge_sense();
      if( geom_sense != bridge_sense )
        result_curve->reverse_bridge_sense();

      coedgesms.clean_out();
      result_curve->coedgesms( coedgesms );
      while( coedgesms.size() )
      { 
        CoEdgeSM* coe = coedgesms.pop();
        if( coe->owner() )
          coe->owner()->remove_bridge( coe );
      }
    }

    names_to_add.insert(names_to_add.end(), edge1_names.begin(), edge1_names.end() );
    names_to_add.insert(names_to_add.end(), edge2_names.begin(), edge2_names.end() );
    
    if ( result_curve )
      end_result = result_curve;
  }
  
  RefEdge* result = NULL;
  if(end_result)
    result = dynamic_cast<RefEdge*>(end_result->topology_entity());
    
  RefEdge* dead = NULL;
  if(result)
  {
    // notify observers that one edge is being composited into another
    // TODO - make a simple function for this notification since it is  times????
    dead = result != refedge1 ? refedge1 : result != refedge2 ? refedge2 : 0;
    update_combined_edges( result, dead );

    //append all names to this ref entity
    int k;
    for(k=0; k<names_to_add.size(); k++ )      
      result->entity_name( names_to_add[k] );
    
  }

  if ( result && update_dag )
  {
    CubitBoolean flag = CGMApp::instance()->attrib_manager()->auto_actuate_flag(CA_MERGE_PARTNER);
    CGMApp::instance()->attrib_manager()->set_auto_actuate_flag(CA_MERGE_PARTNER, CUBIT_TRUE);
    
    DLIList<RefFace*> parent_faces;
    result->ref_faces( parent_faces );
    
    for( i = parent_faces.size(); i--; )
    {
      Surface* surf = parent_faces.get_and_step()->get_surface_ptr();
      GeometryQueryTool::instance()->make_RefFace( surf );
    }
    if ( !parent_faces.size() && end_result) // free curves
    {
      RefEdge* edge = GeometryQueryTool::instance()->make_RefEdge(end_result);
      RefEdge* dead = edge == refedge1 ? refedge2 : edge == refedge2 ? refedge1 : 0;
      assert( dead && !dead->get_curve_ptr() );
      dead->remove_from_DAG();
    }
    assert( !(keep_edge && keep_edge->deactivated()) );
    DAG::instance()->cleanout_deactivated_DAG_nodes();
    CGMApp::instance()->attrib_manager()->set_auto_actuate_flag(CA_MERGE_PARTNER, flag);
  }

  
  return result;
}


/*
CubitStatus CompositeTool::make_mergeable( GeometryEntity* bridge1,
                                           GeometryEntity* bridge2 )
{
  BasicTopologyEntity *bte1, *bte2;
  bte1 = dynamic_cast<BasicTopologyEntity*>(bridge1->topology_entity());
  bte2 = dynamic_cast<BasicTopologyEntity*>(bridge2->topology_entity());
  if( (!bte1) == (!bte2) ) // 1 must be null and one non-null
    return CUBIT_FAILURE;
  
  if( bte2 )
  {
    bte1 = bte2;
    GeometryEntity* tmp = bridge1;
    bridge1 = bridge2;
    bridge2 = tmp;
  }
  
  DLIList<CubitSimpleAttrib*> csa_list;
  bridge1->get_simple_attribute( csa_list );
  CubitSimpleAttrib* csa = 0;
  for( int i = csa_list.size(); !csa && i--; )
    if( CubitAttrib::attrib_type(csa_list.step_and_get()) == CA_MERGE_PARTNER )
      csa = csa_list.get();
  
  if( !csa )
  {
    int merge_id = TDUniqueId::get_unique_id(bte1);
    csa = CAMergePartner::cubit_simple_attrib(merge_id);
    bridge1->append_simple_attribute( csa );
  }
 
  assert(csa);
  bridge2->append_simple_attribute( csa );
  delete csa;
  return CUBIT_SUCCESS;
}
*/
  

//-------------------------------------------------------------------------
// Purpose       : Create a composite surface
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 
//-------------------------------------------------------------------------
RefFace* CompositeTool::composite( DLIList<RefFace*>& face_list,
                                   RefFace* keep_face /* = 0 */ )
{
  int i;
  DLIList<RefFace*> modified_faces;

    // Get all edges in the faces in face_list
  DLIList<ModelEntity*> query_input(face_list.size()), query_output;
  CAST_LIST_TO_PARENT(face_list, query_input);
  ModelQueryEngine::instance()
    ->query_model( query_input, DagType::ref_edge_type(), query_output );
  DLIList<RefEdge*> edge_list(query_output.size());
  CAST_LIST(query_output, edge_list, RefEdge);
  
    // Remove all edges that don't occur in exactly two faces in list
  edge_list.last();
  DLIList<CoEdge*> coedges;
  for (i = edge_list.size(); i--; )
  {
    RefEdge* edge = edge_list.step_and_get();
    coedges.clean_out();
    edge->co_edges(coedges);
    if (coedges.size() == 2)
    {
      RefFace* face1 = coedges.get()->get_ref_face();
      RefFace* face2 = coedges.next()->get_ref_face();
      // Check to make sure both faces are in the faces we are compositing
      // and also make sure the faces are not the same otherwise we
      // will composite out hardlines in one of the faces.
      if (face1 != face2 &&
        face_list.is_in_list(face1) && face_list.is_in_list(face2))
        continue;
    }
    
    edge_list.change_to(0);
  }
  edge_list.remove_all_with_value(0);
  
  if( edge_list.size() == 0 )
  {
    PRINT_ERROR("Cannot create composites from the specified surfaces.\n");
    return 0;
  }

  edge_list.reset();
  for( i = edge_list.size(); i--; )
  {
    RefEdge* edge = edge_list.get_and_step();
    RefFace* face = remove_edge( edge, false, false, keep_face );
    if (face)
      modified_faces.append_unique(face);
  }
  
  for ( i = modified_faces.size(); i--; )
    if( ! modified_faces.step_and_get()->get_surface_ptr() )
      modified_faces.change_to(0);
  modified_faces.remove_all_with_value(0);  
  
  RefFace* result = modified_faces.size() ? modified_faces.get() : 0;

#ifdef ALPHA_TREADSWEEP
  // Make sure names from all of the faces get added to the
  // list of names in the result.
  if(result)
  {
    for(i=face_list.size(); i--;)
    {
      RefFace *cur_face = face_list.get_and_step();
      if(cur_face != result)
      {
        result->merge_entity_names(cur_face);
      }
    }
  }
#endif

  DLIList<RefVolume*> vol_list, surf_vols;
  for ( i = modified_faces.size(); i--; )
  {
    RefFace* face = modified_faces.get_and_step();
    face->set_id( RefEntityFactory::instance()->next_ref_face_id() );
    surf_vols.clean_out();
    face->ref_volumes(surf_vols);
    vol_list.merge_unique(surf_vols);
  }
  
  vol_list.reset();
  for ( i = vol_list.size(); i--; )
  {
    RefVolume* vol = vol_list.get_and_step();
    Lump* lump = vol->get_lump_ptr();
    GeometryQueryTool::instance()->make_Body(lump->bodysm());
  }
  assert( !(keep_face && keep_face->deactivated()) );
  DAG::instance()->cleanout_deactivated_DAG_nodes();
    
  return result;
}

//-------------------------------------------------------------------------
// Purpose       : Composite over a curve
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 03/12/03
//-------------------------------------------------------------------------
RefFace* CompositeTool::remove_edge( RefEdge* edge,
                                     bool remove_partitions, /* = false */
                                     bool update_dag,        /* = true  */
                                     RefFace* keep_face      /* = NULL  */ )
{
  int i;

  DLIList<CoEdge*> edge_coedges;
  edge->co_edges(edge_coedges);
  if (edge_coedges.size() != 2 && 
     (edge->get_curve_ptr()->geometry_type() != POINT_CURVE_TYPE ||
      edge_coedges.size() != 1))
  {
    PRINT_ERROR("Cannot composite over %d-valent curve %d (%s)\n",
      edge_coedges.size(), edge->id(), edge->entity_name().c_str());
    return NULL;
  }

  DLIList<TopologyBridge*> curves(edge->bridge_manager()->number_of_bridges());
  edge->bridge_manager()->get_bridge_list( curves );

  DLIList<Surface*> surfaces;
  assert(curves.size());
  curves.get()->surfaces( surfaces );
  if ( surfaces.size() > 2  || surfaces.size() < 1 )
  {
    PRINT_ERROR("Cannot composite over %d-valence curve %d\n",
      surfaces.size(), edge->id() );
    return 0;
  }
    
  RefFace* face = 0;

  // These will be the two RefFaces associated with the two 
  // surfaces being composited.  We need to store them before we
  // get to the call to remove_curve because after that we won't
  // have access to both of them and we need to pass both of them
  // down to update_combined_faces.
  RefFace *face1 = NULL;
  RefFace *face2 = NULL;
    
  std::vector<CubitString> names_to_add;
  curves.reset();
  Surface* surface = 0;
  for( i = curves.size(); i--; )
  {
    Curve* curve = dynamic_cast<Curve*>(curves.get_and_step());

      // get the two Surfaces that will be involved in the composite operation
      // when the curve is removed
    DLIList<Surface*> surfs;
    curve->surfaces(surfs);
    assert(surfs.size() > 0  && surfs.size() <= 2);  // should be one or two surfs
    face1 = dynamic_cast<RefFace*>(surfs.get()->topology_entity());
    face2 = dynamic_cast<RefFace*>(surfs.next()->topology_entity());

    //Save the names and ids
    surfs.get()->set_saved_id(face1->id());
    DLIList<CubitString*> tmp_names;
    face1->entity_names( tmp_names );        
    surfs.get()->set_saved_names( tmp_names );
    std::vector<CubitString> face1_names;
    for( int i=tmp_names.size(); i--; )
      face1_names.push_back( *(tmp_names.get_and_step()) );
    tmp_names.clean_out();       
    face1->remove_entity_names();
    
    surfs.next()->set_saved_id(face2->id());
    face2->entity_names( tmp_names );        
    surfs.next()->set_saved_names( tmp_names );
    std::vector<CubitString> face2_names;
    for( int i=tmp_names.size(); i--; )
      face2_names.push_back( *(tmp_names.get_and_step()) );
    tmp_names.clean_out();       
    face2->remove_entity_names();

      // if only one surface, the edge should be internal to a composite surface
    surfs.reset();

      // Composite the surfaces
    Surface* result = 0;
    PartitionCurve* partcurve = 0;
    if ( !(partcurve = dynamic_cast<SegmentedCurve*>(curve)) )
      partcurve = dynamic_cast<PartPTCurve*>(curve);
    
    Surface* keep = 0;
    if (keep_face) {
      if (surfs.get()->owner() == keep_face->bridge_manager())
        keep = surfs.get();
      else if (surfs.next()->owner() == keep_face->bridge_manager())
        keep = surfs.next();
    }
    if( remove_partitions && partcurve )
      result = PartitionEngine::instance().remove_curve( partcurve );
    if( !result )
      result = CompositeEngine::instance().remove_curve( curve, remove_partitions, keep );

    if( !result )
    {
      //reapply names
      for( int i=face1_names.size(); i--; )
        face1->entity_name( face1_names[i] );
      for( int i=face2_names.size(); i--; )
        face2->entity_name( face2_names[i] );

      PRINT_ERROR("Failed to remove curve %d\n",edge->id());
      break;
    }
    else if( !face )  // get surviving refface
    {
      face = dynamic_cast<RefFace*>(result->topology_entity());
    }
    else  // merge resulting composites into surviving refface
    {
      RefFace* old_face = dynamic_cast<RefFace*>(result->topology_entity());
      if( old_face != face )
      {
        if( old_face )
          old_face->bridge_manager()->remove_bridge( result );
        face->bridge_manager()->add_bridge( result );
        if( GeometryQueryTool::relative_sense( face->get_surface_ptr(), result )
            == CUBIT_REVERSED && face->get_surface_ptr()->bridge_sense() ==
            result->bridge_sense() )
          result->reverse_bridge_sense();
            
      }
    }
    
    names_to_add.insert(names_to_add.end(), face1_names.begin(), face1_names.end() );
    names_to_add.insert(names_to_add.end(), face2_names.begin(), face2_names.end() );
      
    if( result )
      surface = result;
  }
  
  if( !surface )
  {
    PRINT_ERROR("Failed to remove curve %d\n", edge->id() );
    return 0;
  }
  
  RefFace* result = dynamic_cast<RefFace*>(surface->topology_entity());
  assert( result == face );

  //append all names to this ref entity
  int k;
  for(k=0; k<names_to_add.size(); k++ )      
    result->entity_name( names_to_add[k] );

    // notify observers that one face is being composited into another

  // We need to pass the result face and also the one that wasn't chosen
  // as the result face.
  RefFace *delete_face = face1;
  if(result == delete_face)
    delete_face = face2;
  update_combined_faces( result, edge, delete_face );
  
  if ( update_dag )
  {
    DLIList<RefVolume*> modified_volumes;
    result->ref_volumes(modified_volumes);
    for( i = modified_volumes.size(); i--; )
    {
      RefVolume* volume = modified_volumes.get_and_step();
      Lump* lump = volume->get_lump_ptr();
      GeometryQueryTool::instance()->make_Body(lump->bodysm());
    }

    assert( !(keep_face && keep_face->deactivated()) );
    DAG::instance()->cleanout_deactivated_DAG_nodes();
  }
  
  return result;
}

//-------------------------------------------------------------------------
// Purpose       : Create a composite volume
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 09/25/01
//-------------------------------------------------------------------------
RefVolume* CompositeTool::composite( DLIList<RefVolume*>& volumes,
                                     Body** body_ptr /* = 0 */ )
{
  int i;
  DLIList<RefVolume*> modified_vols;
  DLIList<Body*> input_bodies;
  for (i = volumes.size(); i--; )
    input_bodies.append_unique( volumes.get_and_step()->get_body_ptr() );

    // Get all faces in the volumes
  DLIList<ModelEntity*> query_input(volumes.size()), query_output;
  CAST_LIST_TO_PARENT(volumes, query_input);
  ModelQueryEngine::instance()
    ->query_model( query_input, DagType::ref_face_type(), query_output );
  DLIList<RefFace*> faces(query_output.size());
  CAST_LIST(query_output, faces, RefFace);
  
    // Remove all faces that don't occur in exactly two volumes in list
  faces.last();
  DLIList<CoFace*> cofaces;
  for (i = faces.size(); i--; )
  {
    RefFace* face = faces.step_and_get();
    cofaces.clean_out();
    face->co_faces(cofaces);
    if (cofaces.size() == 2)
    {
      RefVolume* vol1 = cofaces.get()->get_ref_volume();
      RefVolume* vol2 = cofaces.next()->get_ref_volume();
      if (volumes.is_in_list(vol1) && volumes.is_in_list(vol2))
        continue;
    }
    
    faces.change_to(0);
  }
  faces.remove_all_with_value(0);
  
  if( faces.size() == 0 )
  {
    PRINT_ERROR("Cannot create composites from the specified volumes.\n");
    return 0;
  }
  
  faces.reset();
  for( i = faces.size(); i--; )
  {
    RefFace* face = faces.get_and_step();
    RefVolume* vol = remove_face( face, false, false );
    if (vol)
      modified_vols.append_unique(vol);
  }
  
  for ( i = modified_vols.size(); i--; )
    if( ! modified_vols.step_and_get()->get_lump_ptr() )
      modified_vols.change_to(0);
  modified_vols.remove_all_with_value(0);  
  
  RefVolume* result = modified_vols.size() ? modified_vols.get() : 0;
  for ( i = modified_vols.size(); i--; )
  {
    RefVolume* vol = modified_vols.get_and_step();
    vol->set_id( RefEntityFactory::instance()->next_ref_volume_id() );
  }
  
  input_bodies.reset();
  for ( i = input_bodies.size(); i--; )
  {
    Body* body = input_bodies.get_and_step();
    BodySM* sm = body->get_body_sm_ptr();
    if (sm)
      GeometryQueryTool::instance()->make_Body(sm);
    else
      GeometryQueryTool::instance()->destroy_dead_entity( body );
  }
  DAG::instance()->cleanout_deactivated_DAG_nodes();
   
  if (body_ptr)
    *body_ptr = result->get_body_ptr();
  return result;
}

//-------------------------------------------------------------------------
// Purpose       : Composite over a surface
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 03/12/03
//-------------------------------------------------------------------------
RefVolume* CompositeTool::remove_face( RefFace* face,
                                  bool remove_partitions, /* = false */
                                  bool update_dag         /* = true  */)
{
  int i;

  DLIList<RefVolume*> vol_list(2);
  face->ref_volumes( vol_list );

  DLIList<CoFace*> face_cofaces;
  face->co_faces(face_cofaces);
  if (face_cofaces.size() != 2)
  {
    PRINT_ERROR("Cannot composite over %d-valent surface %d (%s)\n",
      face_cofaces.size(), face->id(), face->entity_name().c_str());
    return NULL;
  }

  DLIList<TopologyBridge*> surfs(face->bridge_manager()->number_of_bridges());
  face->bridge_manager()->get_bridge_list( surfs );
  
    // Seam topology
  Surface *surf1 = 0, *surf2 = 0;
  if (surfs.size()  == 2)
  {
    surfs.reset();
    surf1 = dynamic_cast<Surface*>(surfs.get());
    surf2 = dynamic_cast<Surface*>(surfs.next());
    
    DLIList<RefEdge*> edge_list;
    DLIList<TopologyBridge*> bridges, ptcurves;
    DLIList<Surface*> surfaces;
    face->ref_edges( edge_list );
    while( edge_list.size() )
    {
      RefEdge* edge = edge_list.pop();
      edge->bridge_manager()->get_bridge_list( bridges );
      TopologyBridge *curve1 = 0, *curve2 = 0;
      while (bridges.size())
      {
        TopologyBridge* curve = bridges.pop();
        surfaces.clean_out();
        curve->surfaces( surfaces );
        if (surfaces.is_in_list(surf1))
          curve1 = curve;
        if (surfaces.is_in_list(surf2))
          curve2 = curve;
      }
      
      assert(curve1 && curve2);
      if (curve1 == curve2)
        continue;
      
      for (int i = 0; i < 2; i++)
      {
        RefVertex* vtx = i ? edge->end_vertex() : edge->start_vertex();
        vtx->bridge_manager()->get_bridge_list( bridges );
        TopologyBridge *point1 = 0, *point2 = 0;
        while (bridges.size())
        {
          TopologyBridge* point = bridges.pop();
          ptcurves.clean_out();
          point->get_parents(ptcurves);
          if (ptcurves.is_in_list( curve1 ))
            point1 = point;
          if (ptcurves.is_in_list( curve2 ))
            point2 = point;
        }
        
        if (point1 != point2)
          CompositeEngine::instance().stitch_points( 
            dynamic_cast<TBPoint*>(point1), 
            dynamic_cast<TBPoint*>(point2) );
      }
      
      CompositeEngine::instance().stitch_curves( 
            dynamic_cast<Curve*>(curve1), 
            dynamic_cast<Curve*>(curve2) );
    }
  }
  else if (surfs.size() == 1)
  {
     surf1 = dynamic_cast<Surface*>(surfs.get());
     surf2 = 0;
  }
  else
  {
    return 0;
  }

  DLIList<Lump*> lumps, tmp;
  surf1->lumps( lumps );
  if (surf2)
  {
    surf2->lumps(tmp);
    lumps.merge_unique(tmp);
  }
  if ( lumps.size() > 2  || lumps.size() < 1 )
  {
    PRINT_ERROR("Cannot composite over %d-valence surface %d\n",
      lumps.size(), face->id() );
    return 0;
  }
    
    // get the two Lumps that will be involved in the composite operation
    // when the surface is removed
  Lump* lump = 0;
  RefVolume *tmp_vol = dynamic_cast<RefVolume*>(lumps.get()->topology_entity());
  lumps.get()->set_saved_id( tmp_vol->id() );
  DLIList<CubitString*> names;
  tmp_vol->entity_names( names );
  lumps.get()->set_saved_names( names );
  names.clean_out();  

  tmp_vol = dynamic_cast<RefVolume*>(lumps.get()->topology_entity());
  lumps.next()->set_saved_id( tmp_vol->id() );
  tmp_vol->entity_names( names );
  lumps.next()->set_saved_names( names );

    // if only one lump, the surface should be internal to a composite lump
  lumps.reset();

    // Composite the lumps
  PartitionSurface* partsurf = dynamic_cast<PartitionSurface*>(surf1);
  if( surf1 == surf2 && remove_partitions && partsurf )
    lump = PartitionEngine::instance().remove_surface( partsurf );
  if( !lump )
    lump = CompositeEngine::instance().remove_surface( surf1, surf2, remove_partitions );

  if( !lump )
  {
    PRINT_ERROR("Failed to remove surface %d\n",face->id());
    return 0;
  }
  
  RefVolume* result = dynamic_cast<RefVolume*>(lump->topology_entity());
  
    // notify observers that one face is being composited into another
  RefVolume* dead_vol =0;
  if (vol_list.size() > 1)
  {
    vol_list.move_to( result );
    dead_vol = vol_list.next();
  }
  update_combined_vols( result, dead_vol );
  
  if ( update_dag )
  {
    DLIList<Body*> modified_bodies;
    result->bodies(modified_bodies);
    for( i = modified_bodies.size(); i--; )
    {
      Body* body = modified_bodies.get_and_step();
      BodySM* sm = body->get_body_sm_ptr();
      if (sm)
        GeometryQueryTool::instance()->make_Body(sm);
      else
        GeometryQueryTool::instance()->destroy_dead_entity( body );
    }

    DAG::instance()->cleanout_deactivated_DAG_nodes();
  }
  
  return result;
}


//-------------------------------------------------------------------------
// Purpose       : Create a composite body
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 10/02/01
//-------------------------------------------------------------------------
Body* CompositeTool::composite( DLIList<Body*>& )
{
  return 0;
}


//-------------------------------------------------------------------------
// Purpose       : Check for the topological validity of a possible composite
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 
//-------------------------------------------------------------------------
CubitBoolean CompositeTool::okayToComposite( 
                        DLIList<BasicTopologyEntity*>& bte_list,
                        DLIList<BasicTopologyEntity*>* boundary,
                        DLIList<BasicTopologyEntity*>* internal,
                        bool print_errors,
                        bool force_same_parents ) const
{
  //check that all entities are of the same type
  bte_list.reset();
  DagType type = bte_list.get_and_step()->dag_type();
  if( (type != DagType::ref_face_type()) && 
      (type != DagType::ref_edge_type()) && 
      (type != DagType::ref_volume_type()) )
  {
    if (print_errors) 
      PRINT_ERROR("Invalid entities passed to CompositeTool::okayToComposite().\n");
    return CUBIT_FALSE;
  }
  int i;
  for( i = bte_list.size(); i > 1; i-- )
  {
    if( bte_list.get_and_step()->dag_type() != type )
    {
      if (print_errors) PRINT_ERROR("Cannot combine entities of different dimensionality.\n");
      return CUBIT_FALSE;
    }
  }
  
  //check that all entities are C0-continous
  DLIList<RefEntity*> unchecked_entities, current_entities,
                  next_entities, child_entities, temp_list ;
  CAST_LIST_TO_PARENT( bte_list, unchecked_entities);
  current_entities.append( unchecked_entities.pop() );

  while( current_entities.size() > 0 )
  {
    next_entities.clean_out();
    child_entities.clean_out();

    
    // Get all adjacent entities
    for( i = current_entities.size(); i > 0; i-- )
    {
      temp_list.clean_out();
      current_entities.get_and_step()->get_child_ref_entities( temp_list );
      child_entities.merge_unique( temp_list );
    }

    for( i = child_entities.size(); i > 0; i-- )
    {
      temp_list.clean_out();
      child_entities.get_and_step()->get_parent_ref_entities( temp_list );
      next_entities.merge_unique( temp_list );
    }
    
    temp_list.clean_out();
    for( i = next_entities.size(); i > 0; i-- )
    {
      RefEntity* entity_ptr = next_entities.get_and_step();

      if( unchecked_entities.move_to( entity_ptr ) )
      {
        temp_list.append( entity_ptr );
        unchecked_entities.extract();
      }
    }
    current_entities = temp_list;
  }
  if( unchecked_entities.size() == 1 )
  {
    if (print_errors) 
    {
      PRINT_ERROR("%s is not connected to any others.\n",
        unchecked_entities.get()->entity_name().c_str() );
      if( unchecked_entities.get()->dag_type() == DagType::ref_volume_type() )
        PRINT_ERROR("\tTry merging first.\n");
    }
    return CUBIT_FALSE;
  }
  else if( unchecked_entities.size() > 1 )
  {
    if (print_errors) {
      PRINT_ERROR("Entities \n");
      for( i = unchecked_entities.size(); i > 0; i-- )
        PRINT_INFO("%d, ",unchecked_entities.get_and_step()->id());
      PRINT_INFO("do not appear to be connected to the others.\n");
    }
    
    return CUBIT_FALSE;
  }
  
  //Check that all common children have at most two
  //parents (and those parents are in the list of entities to
  //composite).
  unchecked_entities.clean_out();
  child_entities.clean_out();
  CAST_LIST_TO_PARENT( bte_list, unchecked_entities);
  for( i = unchecked_entities.size(); i > 0; i-- )
  {
    temp_list.clean_out();
    unchecked_entities.get_and_step()->get_child_ref_entities( temp_list );
    child_entities.merge_unique( temp_list );
  }
  
  DLIList<SenseEntity*> sense_entity_list;
  for( i = child_entities.size(); i > 0; i-- )
  {
    int parent_count = 0;
    current_entities.clean_out();
    RefEntity* child_ptr = child_entities.step_and_get();
    sense_entity_list.clean_out();
    CAST_TO(child_ptr,BasicTopologyEntity)->get_sense_entity_list( sense_entity_list );
    for( int k = sense_entity_list.size(); k > 0; k-- )
      current_entities.append( sense_entity_list.get_and_step()->
        get_grouping_entity_ptr()->get_basic_topology_entity_ptr() );
    for( int j = current_entities.size(); j> 0; j-- )
      if( unchecked_entities.move_to( current_entities.get_and_step() ) )
        parent_count++;
    switch( parent_count )
    {
      case 0:
//        assert( parent_count != 0 ); break;
        return CUBIT_FALSE;
      case 1:
        if( boundary != NULL ) boundary->
          append( CAST_TO( child_entities.get(), BasicTopologyEntity) );
        break;
      case 2:
        if( current_entities.size() != 2 )
        {
          if (print_errors) PRINT_ERROR("Child entity %d cannot be removed.\n",
                                        child_entities.get()->id() );
          return CUBIT_FALSE;
        }
        if( internal != NULL ) internal->
          append( CAST_TO( child_entities.get(), BasicTopologyEntity) );
        break;
      default:
        if (print_errors) PRINT_ERROR("Entities are not simply connected.\n");
        return CUBIT_FALSE;
    }
  }
  
    //Check that all entities to merge are common to all parent grouping entities.
  if( force_same_parents )
  {
    assert( type.functional_type() == DagType::BasicTopologyEntity_TYPE );
    const DagType type_2 = DagType( type.dimension(), DagType::GroupingEntity_TYPE );
    DLIList<ModelEntity*> first_entity_parents, current_entity_parents;
    bte_list.reset();
    ModelQueryEngine::instance()->query_model( 
      *(bte_list.get_and_step()), type_2, first_entity_parents ); 

    for( i = bte_list.size(); i > 1; i-- )
    {
      current_entity_parents.clean_out();
      ModelQueryEngine::instance()->query_model(
        *(bte_list.get_and_step()), type_2, current_entity_parents );

      if( current_entity_parents != first_entity_parents )  
      {
        if (print_errors) PRINT_ERROR("All entities to be combined must have common "
                                      "topological parents.\n");
        return CUBIT_FALSE;
      }
    }
  }
    
  return CUBIT_TRUE;
}

//-------------------------------------------------------------------------
// Purpose       : Classify child topology
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 
//-------------------------------------------------------------------------
CubitStatus CompositeTool::classify_children( 
                                 DLIList<BasicTopologyEntity*>& bte_list,
                                 DLIList<BasicTopologyEntity*>& boundary_children,
                                 DLIList<BasicTopologyEntity*>& interior_children,
                                 DLIList<BasicTopologyEntity*>& unknown_children )
{
  DLIList<RefEntity*> ref_entity_list, all_children;
  RefEntity *child_ptr, *parent_ptr;
  BasicTopologyEntity* bte_ptr;
  int i, j;
  for( i = bte_list.size(); i > 0; i-- )
  {
    ref_entity_list.clean_out();
    bte_list.get_and_step()->get_child_ref_entities( ref_entity_list );
    all_children.merge_unique( ref_entity_list );
  }
  for( i = all_children.size(); i > 0; i-- )
  {
    ref_entity_list.clean_out();
    child_ptr = all_children.get_and_step();
    child_ptr->get_parent_ref_entities( ref_entity_list );
    
    int parent_count = 0;
    for( j = ref_entity_list.size(); j > 0; j-- )
    {
      parent_ptr = ref_entity_list.get_and_step();
      bte_ptr = CAST_TO( parent_ptr, BasicTopologyEntity );
      if( bte_list.move_to( bte_ptr ) ) parent_count++;
    }
    
    bte_ptr = CAST_TO( child_ptr, BasicTopologyEntity );
    if( parent_count == 1 )
      boundary_children.append( bte_ptr );
    else if( parent_count == ref_entity_list.size() )
      interior_children.append( bte_ptr );
    else
      unknown_children.append( bte_ptr );
  }
  return CUBIT_SUCCESS;
}
  
  
  

//-------------------------------------------------------------------------
// Purpose       : remove a composite edge, restoring the orignial edges.
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 
//-------------------------------------------------------------------------
CubitStatus CompositeTool::uncomposite( RefEdge* composite_edge, 
                                        DLIList<RefEdge*>* restored_edges)
{
  int i;
  DLIList<TBPoint*> hidden_points;
  DLIList<RefFace*> face_list;
  DLIList<Curve*> curve_list;
  
  composite_edge->ref_faces( face_list );
  
  if( composite_edge->bridge_manager()->number_of_bridges() != 1 )
  {
    PRINT_ERROR("Curve %d is a merged curve.  Unmerge before removing composite.\n",
      composite_edge->id() );
    return CUBIT_FAILURE;
  }
  
  CompositeCurve* curve_ptr = dynamic_cast<CompositeCurve*>(composite_edge->get_curve_ptr());
  if( !curve_ptr )
    return CUBIT_FAILURE;
  
  for( i = 0; i < curve_ptr->num_curves(); i++ )  
    curve_list.append( curve_ptr->get_curve(i) );

  //remove all names off this ref entity
  composite_edge->remove_entity_names();
  
  curve_ptr->get_hidden_points( hidden_points );
  CubitStatus result = CUBIT_SUCCESS;
  hidden_points.reset();
  for( i = hidden_points.size(); i--; )
  {
    //remove attributes off the hidden points
    CompositeEngine::strip_attributes( hidden_points.get() );

    if( ! CompositeEngine::instance().
      restore_point( hidden_points.get_and_step() ) )
      result = CUBIT_FAILURE;
  }
  
  DLIList<TopologyBridge*> bridge_list;
  composite_edge->bridge_manager()->get_bridge_list(bridge_list);
  int smallest_id = 0;
  for ( i = bridge_list.size(); i--; )
  {
    GeometryEntity* geom_ptr = dynamic_cast<GeometryEntity*>(bridge_list.get_and_step());
    int saved_id = geom_ptr->get_saved_id();
    if ( saved_id && (!smallest_id || saved_id < smallest_id) &&
         !RefEntityFactory::instance()->get_ref_edge(saved_id) )
      smallest_id = saved_id;
  }
  if ( smallest_id && composite_edge->id() != smallest_id )
    composite_edge->set_id( smallest_id );

  
  for( i = face_list.size(); i--; ) 
  {
    Surface* surf = face_list.get_and_step()->get_surface_ptr();
    GeometryQueryTool::instance()->make_RefFace( surf );
  }
  for( i = face_list.size(); i--; )
    GeometryQueryTool::instance()->destroy_dead_entity( face_list.get_and_step() );
  
  for( i = curve_list.size(); i--; )
  {
    //remove attributes off the underlying curves
    Curve *tmp_curve = curve_list.get_and_step();
    CompositeEngine::strip_attributes( tmp_curve );

    BridgeManager* bm = dynamic_cast<BridgeManager*>(tmp_curve->owner());
    RefEdge* edge = bm ? dynamic_cast<RefEdge*>(bm->topology_entity()) : 0;
    if( restored_edges )
      restored_edges->append( edge );
    DLIList<CubitString*> underlying_names;
    tmp_curve->get_saved_names( underlying_names );
    int k;
    for( k=underlying_names.size(); k--; )    
      edge->entity_name( *(underlying_names.get_and_step()) );    
  }
    
  return result;
}

//-------------------------------------------------------------------------
// Purpose       : Restore a hidden point, remerging points and
//                 parent curves if possible.
//
// Special Notes : Helper function for uncomposite(RefFace*,...)
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 08/06/03
//-------------------------------------------------------------------------
CubitStatus CompositeTool::restore_merged_point( TBPoint* hidden_pt,
                                      DLIList<RefFace*>& modified,
                                      bool force )
{
  int i;
  
  HiddenEntitySet* hs = dynamic_cast<HiddenEntitySet*>(hidden_pt->owner());
  if (!hs)
    return CUBIT_FAILURE;
  CompositeCurve* owner = dynamic_cast<CompositeCurve*>(hs->owner());
  if (!owner)
    return CUBIT_FAILURE;
  
  RefEdge* old_edge = dynamic_cast<RefEdge*>(owner->topology_entity());
  int num_curves;
  DLIList<TopologyBridge*> curve_bridges;
  if (old_edge)
  {
    num_curves = old_edge->bridge_manager()->number_of_bridges();
    old_edge->bridge_manager()->get_bridge_list(curve_bridges);
  }
  else
  {
    num_curves = 1;
    curve_bridges.append(owner);
  }
     
  DLIList<CompositeCurve*> curves(num_curves);
  CAST_LIST(curve_bridges, curves, CompositeCurve);
  DLIList<TBPoint*> points(curves.size());
  
  assert(curves.is_in_list(owner));
  CubitVector position = hidden_pt->coordinates();
  
  DLIList<TBPoint*> curve_pts;
  curves.reset();
  GeometryQueryTool* gqt = GeometryQueryTool::instance();
  double factor = gqt->get_geometry_factor();
  for ( i = curves.size(); i--; )
  {
    curve_pts.clean_out();
    curves.get_and_step()->get_hidden_points(curve_pts);
    for (int j = curve_pts.size(); j--; )
    {
      TBPoint* curve_pt = curve_pts.get_and_step();
      if (gqt->about_spatially_equal(position, curve_pt->coordinates(), factor))
      {
        points.append(curve_pt);
        break;
      }
    }
  }
  
  if (!force && points.size() != curve_bridges.size())
  {
    RefEdge* edge = dynamic_cast<RefEdge*>(owner->topology_entity());
    PRINT_ERROR("Cannot restore vertex hidden by curve %d without unmerging "
                "curve and possibly adjacent surfaces.  Use the 'force' "
                "option to force unmerging of geometry.\n", edge ? edge->id() : 0);
    return CUBIT_FAILURE;
  }
  
  points.last();
  for ( i = points.size(); i--; )
  {
    TBPoint* pt = points.step_and_get();
    if (CompositePoint* cp = dynamic_cast<CompositePoint*>(pt))
      pt = cp->get_point();
    
    if (!CompositeEngine::instance().restore_point(pt))
      points.change_to(0);
    else if (CompositePoint* cp = dynamic_cast<CompositePoint*>(pt->owner()))
      points.change_to(cp);
    else
      points.change_to(pt);
  }
  points.remove_all_with_value(0);
  
  if (!points.size())
    return CUBIT_FAILURE;
  
  points.reset();
  RefVertex* new_vtx = gqt->make_RefVertex(points.get());
  for ( i = points.size(); i--; )
  {
    TBPoint* pt = points.get_and_step();
    if (pt->owner())
      assert(pt->topology_entity() == new_vtx);
    else
      new_vtx->bridge_manager()->add_bridge(pt);
  }
          
  CubitVector center_1, center_2;
  DLIList<TopologyBridge*> pt_curves(2);
  points.reset();
  points.get()->get_parents(pt_curves);
  
    // Ignore composite point-curves.  They are transient and will be
    // deleted later.  Don't try to create RefEdges for them now.
  if (pt_curves.size() != 2)
  {
    assert(pt_curves.size() == 1);
    CompositeCurve* ccurve = dynamic_cast<CompositeCurve*>(pt_curves.get());
    assert(ccurve->num_curves() == 0 && ccurve->geometry_type() == POINT_CURVE_TYPE );
    return CUBIT_SUCCESS;
  }
  
  dynamic_cast<Curve*>(pt_curves.get())->position_from_fraction(0.5, center_1);
  dynamic_cast<Curve*>(pt_curves.next())->position_from_fraction(0.5, center_2);
  
  DLIList<Curve*> new_edge_1_curves(points.size()), 
                  new_edge_2_curves(points.size());
  
  points.reset();
  for ( i = points.size(); i--; )
  {
    pt_curves.clean_out();
    points.get_and_step()->get_parents(pt_curves);
    assert(pt_curves.size() == 2);
    Curve* curve_1 = dynamic_cast<Curve*>(pt_curves.get());
    Curve* curve_2 = dynamic_cast<Curve*>(pt_curves.next());
    CubitVector ptc1, ptc2;
    curve_1->position_from_fraction( 0.5, ptc1 );
    curve_2->position_from_fraction( 0.5, ptc2 );
    bool close1 = (center_1 - ptc1).length_squared() < (center_2 - ptc1).length_squared();
    bool close2 = (center_1 - ptc2).length_squared() < (center_2 - ptc2).length_squared();
    if (close1 && !close2)
    {
      new_edge_1_curves.append(curve_1);
      new_edge_2_curves.append(curve_2);
    }
    else if (!close1 && close2)
    {
      new_edge_2_curves.append(curve_1);
      new_edge_1_curves.append(curve_2);
    }
    else
    {
      assert(close1 != close2);
    }
  }

  for ( i = 0; i < 2; i++ )
  {
    int j;
    DLIList<Curve*>& new_edge_curves = i ? new_edge_2_curves : new_edge_1_curves;
    DLIList<Curve*>& otr_edge_curves = i ? new_edge_1_curves : new_edge_2_curves;
    RefEdge* new_edge = 0;
    new_edge_curves.reset();
    for (j = new_edge_curves.size(); j--; )
    {
      Curve* curve = new_edge_curves.get_and_step();
      if (curve->owner())
      {
        RefEdge* edge = dynamic_cast<RefEdge*>(curve->topology_entity());
        assert(!!edge);
        if (!otr_edge_curves.is_in_list(edge->get_curve_ptr()))
          new_edge = edge;
      }
    }
    if (!new_edge)
    {
      Curve* curve = new_edge_curves.get();
      if (curve->owner())
        curve->owner()->remove_bridge(curve);
      new_edge = gqt->make_RefEdge(curve);
    }
  
    Curve* first = new_edge->get_curve_ptr();
    for (j = new_edge_curves.size(); j--; )
    {
      Curve* curve = new_edge_curves.get_and_step();
      if (curve->owner())
      {
        if (curve->topology_entity() == new_edge)
          continue;
        else
          curve->owner()->remove_bridge(curve);
      }
      
      bool reversed = first->relative_sense(curve) == CUBIT_REVERSED;
      bool saved = first->bridge_sense() == curve->bridge_sense();
      if (reversed == saved)
        curve->reverse_bridge_sense();
      new_edge->bridge_manager()->add_bridge(curve);
    }
  }
  
  if (old_edge)
  {
    DLIList<RefFace*> face_list;
    old_edge->ref_faces(face_list);
    modified += face_list;
  }
  return CUBIT_SUCCESS;
}
  
  
  
  

  


//-------------------------------------------------------------------------
// Purpose       : remove a composite face, restoring the orignial faces.
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 
//-------------------------------------------------------------------------
CubitStatus CompositeTool::uncomposite( RefFace* composite_face,
                                        DLIList<RefFace*>* restored_faces,
                                        bool force_unmerge )
{
  int i;
  DLIList<Curve*> hidden_curves;
  DLIList<RefVolume*> vol_list;
  DLIList<Surface*> surface_list;
  
  composite_face->ref_volumes( vol_list );
  
  if( composite_face->bridge_manager()->number_of_bridges() != 1 )
  {
    PRINT_ERROR("Surface %d is a merged curve.  Unmerge before removing composite.\n",
      composite_face->id() );
    return CUBIT_FAILURE;
  }
  
  CompositeSurface* surf_ptr = dynamic_cast<CompositeSurface*>(composite_face->get_surface_ptr());
  if( !surf_ptr )
    return CUBIT_FAILURE;

  //remove all names off this ref entity
  composite_face->remove_entity_names();
  
  for( i = 0; i < surf_ptr->num_surfs(); i++ )
    surface_list.append( surf_ptr->get_surface(i) );
  
  surf_ptr->get_hidden_curves( hidden_curves );
  CubitStatus result = CUBIT_SUCCESS;
  hidden_curves.reset();
  DLIList<TopologyBridge*> points(2);
  DLIList<RefFace*> modified_faces;
  for( i = hidden_curves.size(); i--; )
  {
    //remove attributes off the hidden curves 
    CompositeEngine::strip_attributes( hidden_curves.get() );

    Curve* curve = hidden_curves.get_and_step();

    if(DEBUG_FLAG(87))
    {
      GMem gmem;
      curve->get_geometry_query_engine()->get_graphics(curve,&gmem);
      GfxDebug::draw_polyline(gmem.point_list(),gmem.pointListCount,CUBIT_RED);
      GfxDebug::flush();
    }

    curve->get_children(points, true);
    while (points.size())
    {
      TBPoint* pt = dynamic_cast<TBPoint*>(points.pop());

      //remove attributes off the underlying points 
      CompositeEngine::strip_attributes( pt ); 

      HiddenEntitySet* hs = dynamic_cast<HiddenEntitySet*>(pt->owner());
      if (hs && dynamic_cast<CompositeCurve*>(hs->owner()))
      {
        if (!restore_merged_point(pt, modified_faces, force_unmerge))
        {
          PRINT_ERROR("Failed to restore hidden vertex.\n");
          result = CUBIT_FAILURE;
          break;
        }
      }
    }
    
    if (!result)
      break;
    
    if( ! CompositeEngine::instance().restore_curve( curve ) )
      result = CUBIT_FAILURE;
  }
  
  if (modified_faces.size())
  {
    DLIList<ModelEntity*> query_input(modified_faces.size());
    CAST_LIST_TO_PARENT(modified_faces, query_input);
    DLIList<ModelEntity*> query_output;
    ModelQueryEngine::instance()->
      query_model(query_input,DagType::ref_volume_type(),query_output);
    while(query_output.size())
      vol_list.append_unique(dynamic_cast<RefVolume*>(query_output.pop()));
  }
                                 
  
  DLIList<TopologyBridge*> bridge_list;
  composite_face->bridge_manager()->get_bridge_list(bridge_list);
  int smallest_id = 0;
  for ( i = bridge_list.size(); i--; )
  {
    GeometryEntity* geom_ptr = dynamic_cast<GeometryEntity*>(bridge_list.get_and_step());
    int saved_id = geom_ptr->get_saved_id();
    if ( saved_id && (!smallest_id || saved_id < smallest_id) &&
         !RefEntityFactory::instance()->get_ref_face(saved_id) )
      smallest_id = saved_id;
  }
  if ( smallest_id && composite_face->id() != smallest_id )
    composite_face->set_id( smallest_id );
  
  for( i = vol_list.size(); i--; )
  {
    Lump* lump= dynamic_cast<Lump*>(
      vol_list.get_and_step()->bridge_manager()->topology_bridge() );
    GeometryQueryTool::instance()->make_Body( lump->bodysm() );
  }
  for( i = vol_list.size(); i--; )
    GeometryQueryTool::instance()->destroy_dead_entity( vol_list.get_and_step() );
  
  for( i = surface_list.size(); i--; )
  {
    Surface *tmp_surf = surface_list.get_and_step();
    //remove attributes off the surfaces 
    CompositeEngine::strip_attributes( tmp_surf ); 

    BridgeManager* bm = dynamic_cast<BridgeManager*>(tmp_surf->owner());
    RefFace* face = bm ? dynamic_cast<RefFace*>(bm->topology_entity()) : 0;
    if( restored_faces )
      restored_faces->append( face );
    DLIList<CubitString*> underlying_names;
    tmp_surf->get_saved_names( underlying_names );
    int k;
    for( k=underlying_names.size(); k--; )
      face->entity_name( *(underlying_names.get_and_step()) );

  }
  
  return result;
}


//-------------------------------------------------------------------------
// Purpose       : Remove a composite volume
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 09/25/01
//-------------------------------------------------------------------------
CubitStatus CompositeTool::uncomposite( RefVolume* composite_vol,
                                        DLIList<RefVolume*>* restored_vols )
{
  int i;
  DLIList<Surface*> hidden_surfs;
  DLIList<Body*> body_list;
  DLIList<Lump*> lump_list;
  
  composite_vol->bodies( body_list );
  
  CompositeLump* lump_ptr = dynamic_cast<CompositeLump*>(composite_vol->get_lump_ptr());
  if( !lump_ptr )
    return CUBIT_FAILURE;
  
  for( i = 0; i < lump_ptr->num_lumps(); i++ )
    lump_list.append( lump_ptr->get_lump(i) );
  
  lump_ptr->get_hidden_surfaces( hidden_surfs );
  CubitStatus result = CUBIT_SUCCESS;
  hidden_surfs.reverse();
  while( hidden_surfs.size() )
  {
    Surface* surf = hidden_surfs.pop();
    Surface* stitch = 0;
    if( ! CompositeEngine::instance().restore_surface( surf, stitch ) )
      result = CUBIT_FAILURE;
    if (stitch)
      hidden_surfs.remove( stitch );
  }
  
  int saved_id = composite_vol->get_lump_ptr()->get_saved_id();
  if (saved_id && saved_id < composite_vol->id() &&
      !RefEntityFactory::instance()->get_ref_volume(saved_id))
    composite_vol->set_id( saved_id );
  
  for( i = body_list.size(); i--; )
    GeometryQueryTool::instance()->make_Body( body_list.get_and_step()->get_body_sm_ptr() );
  for( i = body_list.size(); i--; )
    GeometryQueryTool::instance()->destroy_dead_entity( body_list.get_and_step() );
  
  if( restored_vols )
  {
    for( i = lump_list.size(); i--; )
    {
      TBOwner* owner = lump_list.get_and_step()->owner();
      if (BridgeManager* bm = dynamic_cast<BridgeManager*>(owner))
        restored_vols->append( dynamic_cast<RefVolume*>(bm->topology_entity()) );
      else if (CompositeLump* lump = dynamic_cast<CompositeLump*>(owner))
        restored_vols->append_unique( dynamic_cast<RefVolume*>(lump->topology_entity() ) );
    }
  }
  
  return result;
}

//-------------------------------------------------------------------------
// Purpose       : Remove composite body
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 10/02/01
//-------------------------------------------------------------------------
CubitStatus CompositeTool::uncomposite( Body* , 
                                        DLIList<Body*>* )
{
  return CUBIT_FAILURE;
}




//-------------------------------------------------------------------------
// Purpose       : Create multiple composite curves
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 12/03/99
//-------------------------------------------------------------------------
CubitStatus CompositeTool::composite( DLIList<RefEdge*>& edges, 
                                      DLIList<RefEdge*>& new_edge_list,
                                      DLIList<RefVertex*>* vertices_to_keep, 
                                      double max_curve_angle /*radians*/ )
{
  new_edge_list.clean_out();
  double min_cos_angle = cos( max_curve_angle );
  int i;
  
  //check for duplicates
  DLIList<RefEdge*> edges_to_composite;
  for( i = edges.size(); i > 0; i-- )
  {
    RefEdge* edge_ptr = edges.get_and_step();
    if( ! edges_to_composite.append_unique( edge_ptr ) )
    {
      PRINT_WARNING("Duplicate curve %d in list of curves to composite.\n",
        edge_ptr->id());
    }
  }
  
  
  if( edges_to_composite.size() < 2 ) return CUBIT_SUCCESS;
  
  fast_edge_sort( edges_to_composite, true );
  
  //look for a place to start (a point that we cannot 
  //composite across.)

  RefVertex *vtx_ptr;
  RefEdge* first_edge;
  DLIList<RefEdge*> vertex_edges;

  edges_to_composite.last();
  RefEdge* prev_edge = edges_to_composite.get();
  edges_to_composite.reset();
  
  for( i = edges_to_composite.size(); i >= 0; i-- )
  {
    first_edge = edges_to_composite.get_and_step();
    vtx_ptr = first_edge->common_ref_vertex( prev_edge );
    if( ! vtx_ptr ) break;
    
    vertex_edges.clean_out();
    vtx_ptr->ref_edges( vertex_edges );
    vertex_edges.remove( first_edge );
    vertex_edges.remove( prev_edge );
    if( vertex_edges.size() > 0 ) break;
    
    prev_edge = first_edge;
  }
  
  if( i < 0 ) //if all the edges get composited into one
  {
    edges_to_composite.reset();
    first_edge = edges_to_composite.get();
  }
    
  //composite chains of consecutive edges
  
  edges_to_composite.move_to( first_edge );
  DLIList<RefEdge*> edge_list, sublist;
  int edge_count = edges_to_composite.size();
  while( edge_count > 0 )
  {
    edge_list.clean_out();
    RefEdge* edge_ptr = edges_to_composite.get_and_step();
    edge_list.append( edge_ptr );
    edge_count--;
    
    //Find a chain of connected edges that can form a 
    //topologically valid composite.
    
    while( edge_count > 0 )
    {
      RefEdge* next_edge = edges_to_composite.get();
      vtx_ptr = edge_ptr->common_ref_vertex( next_edge );
      if( ! vtx_ptr ) break;
      
      vertex_edges.clean_out();
      vtx_ptr->ref_edges( vertex_edges );
      vertex_edges.remove( edge_ptr );
      vertex_edges.remove( next_edge );
      if( vertex_edges.size() > 0 ) break;
      
      edge_ptr = next_edge;
      edges_to_composite.step();
      edge_list.append( edge_ptr );
      edge_count--;
    }
      
    if( edge_list.size() == 1 ) continue;
    assert( edge_list.size() );
    
    //Now split that chain of edges further based on
    //user-specified constraints (vertices and/or angle.)
    //This is done seperatly from the above loop to handle
    //cases where the edge_list forms a closed loop, and the
    //user specified one position on that loop to keep.
    
    //Find a start point in the chain.
    edge_list.last();
    edge_ptr = edge_list.get();
    edge_list.reset();
    RefVertex* keep_vertex = 0;
    if( edge_ptr->common_ref_vertex( edge_list.get() ) )
    //if we have a closed loop
    {
      edge_ptr = edge_list.get_and_step();
      for( i = edge_list.size(); i > 0; i-- )
      {
        RefEdge* next_edge = edge_list.get();
        vtx_ptr = edge_ptr->common_ref_vertex( next_edge );
        assert( vtx_ptr != 0 );
        
        if( vertices_to_keep && vertices_to_keep->is_in_list( vtx_ptr ) )
        {
          keep_vertex = vtx_ptr;
          break;
        }
        CubitVector cur_tan = tangent( edge_ptr, vtx_ptr ) * -1.0;
        CubitVector next_tan = tangent( next_edge, vtx_ptr );
        double cos_angle = (cur_tan % next_tan) / 
                       (cur_tan.length() * next_tan.length());
        if( cos_angle < min_cos_angle )
        {
          keep_vertex = vtx_ptr;
          break;
        }
        
        edge_ptr = next_edge;
        edge_list.step();
      }
    }
    
    //If we didn't find any split locations, that's okay.
    //We'll let the compositing code choose which vertex to keep.
    
    int subcount = edge_list.size();
    while( subcount > 0 )
    {
      sublist.clean_out();
      edge_ptr = edge_list.get_and_step();
      subcount--;
      sublist.append( edge_ptr );
      
      while( subcount > 0 )
      {
        RefEdge* next_edge = edge_list.get();
        vtx_ptr = edge_ptr->common_ref_vertex( next_edge );
        assert( vtx_ptr != 0 );
        
        if( vertices_to_keep && vertices_to_keep->is_in_list( vtx_ptr ) )
        {
          keep_vertex = vtx_ptr;
          break;
        }
        CubitVector cur_tan = tangent( edge_ptr, vtx_ptr ) * -1.0;
        CubitVector next_tan = tangent( next_edge, vtx_ptr );
        double cos_angle = (cur_tan % next_tan) / 
                        (cur_tan.length() * next_tan.length());
        if( cos_angle < min_cos_angle )
        {
          keep_vertex = vtx_ptr;
          break;
        }
        
        edge_ptr = next_edge;
        edge_list.step();
        subcount--;
        sublist.append( edge_ptr );
      }
      
    
      if( sublist.size() > 1 )
      {
        edge_ptr = composite( sublist, keep_vertex );
        if( ! edge_ptr ) return CUBIT_FAILURE;
        new_edge_list.append( edge_ptr );
      }
    }
  }
  
  return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Purpose       : Get the tangent along an edge away from a vertex.
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 12/03/99
//-------------------------------------------------------------------------
CubitVector CompositeTool::tangent( RefEdge* edge_ptr, RefVertex* vtx_ptr ) const
{
  CubitVector result;
  edge_ptr->tangent( vtx_ptr->coordinates(), result );
  if( vtx_ptr == edge_ptr->end_vertex() )
    result *= -1;
  return result;
} 

//-------------------------------------------------------------------------
// Purpose       : Sort edges topologically.  
//
// Special Notes : Not safe with secondary links
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 12/03/99
//-------------------------------------------------------------------------
CubitStatus CompositeTool::fast_edge_sort( DLIList<RefEdge*>& edge_list,
                                           bool valence2_vertices ) const
{
  int i;
  if( edge_list.size() < 2 ) return CUBIT_SUCCESS;
  
  //make sure all RefEdges are initially marked as zero.
  DLIList<RefVertex*> vtx_list(2), all_vertices(edge_list.size() * 2);
  DLIList<RefEdge*> vertex_edges;
  RefEdge* edge_ptr;
  RefVertex* vtx_ptr;
  for( i = edge_list.size(); i > 0; i-- )
  {
    edge_ptr = edge_list.get_and_step();
    vtx_list.clean_out();
    edge_ptr->ref_vertices( vtx_list );
    for( int j = vtx_list.size(); j > 0; j-- )
    {
      vtx_ptr = vtx_list.get_and_step();
      vtx_ptr->marked( 0 );
    }
  }
  for( i = edge_list.size(); i > 0; i-- )
  {
    edge_ptr = edge_list.get_and_step();
    vtx_list.clean_out();
    edge_ptr->ref_vertices( vtx_list );
    for( int j = vtx_list.size(); j > 0; j-- )
    {
      vtx_ptr = vtx_list.get_and_step();
      if( vtx_ptr->marked() ) continue;
    
      all_vertices.append( vtx_ptr );
      vtx_ptr->marked( 1 );

      vertex_edges.clean_out();
      vtx_ptr->ref_edges( vertex_edges );
      for( int k = vertex_edges.size(); k > 0; k-- )
        vertex_edges.get_and_step()->marked( 0 );
    }
  }
  
  //now mark all of the edges in the edge list with a 1
  for( i = edge_list.size(); i > 0; i-- )
    edge_list.get_and_step()->marked( 1 );
    
  DLIList<RefEdge*> sorted_edge_list( edge_list.size() );
  
  int passes = 0;
  while( sorted_edge_list.size() < edge_list.size() ) 
  {
    passes++;
    RefVertex* start_vtx = 0;
    RefEdge* start_edge = 0;
    
    //choose a start_vtx and start_edge
    for( i = all_vertices.size(); i > 0; i-- )
    {
      vtx_ptr = all_vertices.get_and_step();
      //if( ! vtx_ptr->marked() ) continue;
      
      int edge_count = 0;
      vertex_edges.clean_out();
      vtx_ptr->ref_edges( vertex_edges );
      for( int j = vertex_edges.size(); j > 0; j-- )
      {
        edge_ptr = vertex_edges.get_and_step();
        if( edge_ptr->marked() )
        {
          start_edge = edge_ptr;
          edge_count++;
          if( (edge_count > 1) || (valence2_vertices && vertex_edges.size() > 2) ) 
            break;
          // If the vertex has more than two edges, consider it a 
          // valid starting point for the search because it will not
          // be possible to create a composite across that vertex
          // anyway.
        }
      }
      if( edge_count == 1 )
      {
        start_vtx = vtx_ptr;
        break;
      }
    }
    
    //If we have a circular loop, choose any vertex as the start vertex
    if( ! start_vtx )
    {
      start_edge = 0;
      all_vertices.reset();
      for( i = all_vertices.size(); ! start_edge && (i > 0); i-- )
      {
        start_vtx = all_vertices.get_and_step();
        vertex_edges.clean_out();
        start_vtx->ref_edges( vertex_edges );
        for( int j = vertex_edges.size(); j > 0; j-- )
        {
          edge_ptr = vertex_edges.get_and_step();
          if( edge_ptr->marked() )
          {
            start_edge = edge_ptr;
            break;
          }
        }
        if( start_vtx->marked() ) break;
      }
    }

    assert( start_edge != 0 );
    assert( start_vtx != 0 );


    start_vtx->marked( 0 );
    edge_ptr = start_edge;
    vtx_ptr = start_vtx;
    
    do
    {
      sorted_edge_list.append( edge_ptr );
      edge_ptr->marked( 0 );
      
      RefVertex* next_vtx = edge_ptr->other_vertex( vtx_ptr );
      next_vtx->marked( 0 );
      
      vertex_edges.clean_out();
      next_vtx->ref_edges( vertex_edges );
      RefEdge* next_edge = 0;
      
      if( !valence2_vertices || vertex_edges.size() < 3 )
      {
        for( i = vertex_edges.size(); i > 0; i-- )
        {
          RefEdge* edge2_ptr = vertex_edges.get_and_step();
          if( edge2_ptr->marked() )
          {
            next_edge = edge2_ptr;
            break;
          }
        }
      }
      
      edge_ptr = next_edge;
      vtx_ptr = next_vtx;
    } while( edge_ptr );
  }
  
  assert( sorted_edge_list.size() == edge_list.size() );
  
  edge_list.clean_out();
  edge_list = sorted_edge_list;
  return (passes == 1) ? CUBIT_SUCCESS : CUBIT_FAILURE;
}


CubitStatus CompositeTool::composite( DLIList<RefFace*>& faces_to_composite,
                                      DLIList<RefFace*>& new_composites,
                                      double max_angle )
{
  int i;
  DLIList<CoEdge*> coedges(2);
  double cos_max_angle = cos( max_angle );
  CubitStatus result = CUBIT_SUCCESS;
  
  DLIList<BasicTopologyEntity*> input, boundary, internal, other;
  CAST_LIST_TO_PARENT( faces_to_composite, input );
  classify_children( input, boundary, internal, other );
  DLIList<RefEdge*> internal_edges(internal.size());
  CAST_LIST( internal, internal_edges, RefEdge );
  if ( !internal_edges.size() )
    return CUBIT_FAILURE;
    
    // Group faces into sets to be composited
  
  std::vector< DLIList<RefFace*> > face_sets;
  for ( i = faces_to_composite.size(); i--; )
    faces_to_composite.get_and_step()->marked(0);
  
  while ( internal_edges.size() )
  {
    RefEdge* edge = internal_edges.pop();
    coedges.clean_out();
    edge->co_edges(coedges);
    if ( coedges.size() != 2 )
      continue;
    
    CoEdge* coedge1 = coedges.get_and_step();
    CoEdge* coedge2 = coedges.get_and_step();
    RefFace* face1 = coedge1->get_ref_face();
    RefFace* face2 = coedge2->get_ref_face();
    
    if ( max_angle < 4.0 )
    {
      CubitVector point, norm1, norm2;
      bool reversed = coedge1->get_sense() == coedge2->get_sense();

      edge->position_from_fraction( 0.33, point );
      norm1 = face1->normal_at( point );
      norm2 = face2->normal_at( point );
      if ( reversed ) norm2 *= -1.0;
      double cos_angle = (norm1 % norm2)/(norm1.length()*norm2.length());
      if ( cos_angle < cos_max_angle )
        continue;

      edge->position_from_fraction( 0.67, point );
      norm1 = face1->normal_at( point );
      norm2 = face2->normal_at( point );
      if ( reversed ) norm2 *= -1.0;
      cos_angle = (norm1 % norm2)/(norm1.length()*norm2.length());
      if ( cos_angle < cos_max_angle )
        continue;
    }
    
    if ( face1->marked() && face2->marked() )
    {
      if ( face1->marked() != face2->marked() )
      {
        DLIList<RefFace*>& merge_set = face_sets[face2->marked()-1];
        face_sets[face1->marked()-1] += merge_set;
        while( merge_set.size() )
          merge_set.pop()->marked( face1->marked() );
      }
    }
    else if ( face2->marked() )
    {
      face1->marked(face2->marked());
      face_sets[face2->marked()-1].append(face1);
    }
    else
    {
      if ( !face1->marked() )
      {
        face_sets.push_back(DLIList<RefFace*>());
        face1->marked(face_sets.size());
        face_sets[face1->marked()-1].append(face1);
      }
      face2->marked(face1->marked());
      face_sets[face1->marked()-1].append(face2);
    }
  }
  
    // clear all marks
  std::vector< DLIList<RefFace*> >::iterator s_itor;
  for ( s_itor = face_sets.begin(); s_itor != face_sets.end(); ++s_itor )
  {
    for ( i = s_itor->size(); i--; )
      s_itor->get_and_step()->marked(0);
  }

    // composite faces
  for ( s_itor = face_sets.begin(); s_itor != face_sets.end(); ++s_itor )
  {
    if ( s_itor->size() < 2 )
      continue;
    
    RefFace* comp = composite(*s_itor);
    if ( comp )
      new_composites.append(comp);
    else
      result = CUBIT_FAILURE;
  }

  return result;
}

//-------------------------------------------------------------------------
// Purpose       : Composite surfaces and curves.
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 04/13/00
//-------------------------------------------------------------------------
RefFace* CompositeTool::composite( DLIList<RefFace*>& faces_to_composite,
                                   DLIList<RefEdge*>& result_edges,
                                   DLIList<RefVertex*>* vertices_to_keep,
                                   double max_vertex_angle )
{
  result_edges.clean_out();
  RefFace* result = composite( faces_to_composite );
  if( ! result ) return 0;
  
  DLIList<RefEdge*> face_edges;
  result->ref_edges( face_edges );
  composite( face_edges, result_edges, vertices_to_keep, max_vertex_angle );
  
  return result;
}

//-------------------------------------------------------------------------
// Purpose       : Composite surfaces and curves.
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 04/13/00
//-------------------------------------------------------------------------
CubitStatus CompositeTool::composite( DLIList<RefFace*>& faces_to_composite,
                                      DLIList<RefFace*>& result_faces,
                                      DLIList<RefEdge*>& result_edges,
                                      double max_curve_angle,
                                      DLIList<RefVertex*>* vertices_to_keep,
                                      double max_vtx_angle  )
{
  result_faces.clean_out();
  result_edges.clean_out();
  
  CubitStatus status = 
    composite( faces_to_composite, result_faces, max_curve_angle );
  if( ! status ) return status;
  
  DLIList<ModelEntity*> query_input, query_output;
  CAST_LIST_TO_PARENT( result_faces, query_input );
  ModelQueryEngine::instance()->
    query_model( query_input, DagType::ref_edge_type(), query_output );
  
  DLIList<RefEdge*> edge_list;
  CAST_LIST( query_output, edge_list, RefEdge );
  
  return composite( edge_list, result_edges, vertices_to_keep, max_vtx_angle );
}

  

  
CubitStatus CompositeTool::stitch_points( TBPoint* pt1, TBPoint* pt2 )
{
  if ( pt1->owner() && pt2->owner() && pt1->owner() != pt2->owner() )
    return CUBIT_FAILURE;
  
  CompositePoint* result = CompositeEngine::instance().stitch_points( pt1, pt2 );
  return result ? CUBIT_SUCCESS : CUBIT_FAILURE;
}

CubitStatus CompositeTool::stitch_curves( Curve* cv1, Curve* cv2 )
{
  if ( cv1->owner() && cv2->owner() && cv1->owner() != cv2->owner() )
    return CUBIT_FAILURE;
  
  DLIList<TopologyBridge*> points;
  TBPoint *start1, *end1, *start2, *end2;

  points.clean_out();
  cv1->get_children( points );
  points.reset();
  start1 = dynamic_cast<TBPoint*>(points.get_and_step());
  end1 = dynamic_cast<TBPoint*>(points.get_and_step());

  points.clean_out();
  cv2->get_children( points );
  points.reset();
  start2 = dynamic_cast<TBPoint*>(points.get_and_step());
  end2 = dynamic_cast<TBPoint*>(points.get_and_step());
  
  assert( start1 && start2 && end1 && end2 );
  
  if ( (start1 == end1) != (start2 == end2) )
    return CUBIT_FAILURE;
  
  if ( start1 == end2 || start1->owner() == end2->owner() )
    std::swap(start2, end2);
  
  if ( (start1->owner() != start2->owner()) || (end1->owner() != end2->owner()) )
    return CUBIT_FAILURE;
  
  if ( start1 != start2 && !stitch_points( start1, start2 ) )
    return CUBIT_FAILURE;
  
  if ( start1 != end1 && end1 != end2 && !stitch_points( end1, end2 ) )
    return CUBIT_FAILURE;
  
  CompositeCurve* result = CompositeEngine::instance().stitch_curves( cv1, cv2 );
  return result ? CUBIT_SUCCESS : CUBIT_FAILURE;
}

        
//-------------------------------------------------------------------------
// Purpose       : Functions overloaded by CompositeToolMesh
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 05/10/04
//-------------------------------------------------------------------------
CubitStatus CompositeTool::update_combined_vols( RefVolume*, RefVolume* ) 
  {return CUBIT_SUCCESS;}
CubitStatus CompositeTool::update_combined_faces( RefFace*, RefEdge*, RefFace* ) 
  {return CUBIT_SUCCESS;}
CubitStatus CompositeTool::update_combined_edges( RefEdge*, RefEdge* ) 
  {return CUBIT_SUCCESS;}
