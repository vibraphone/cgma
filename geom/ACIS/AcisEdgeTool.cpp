//- Class:       AcisEdgeTool
//- Description: Functions to create surfaces
// 
//- Owner: Eric Nielsen      
//- Checked by:
//- Version: $Id:

// ********** BEGIN ACIS INCLUDES             **********

#if CUBIT_ACIS_VERSION < 1100
#include "kernel/kernapi/api/kernapi.hxx"
#include "kernel/kerndata/top/body.hxx"
#include "kernel/kerndata/top/vertex.hxx"
#include "kernel/kerndata/geom/point.hxx"
#include "intersct/sg_husk/split/esplit.hxx"
#else
#include "kernapi.hxx" 
#include "vertex.hxx" 
#include "point.hxx"
#include "coverapi.hxx"
#include "esplit.hxx"
#endif

// ********** END ACIS INCLUDES               **********

// ********** BEGIN CUBIT INCLUDES            **********
#include "AcisEdgeTool.hpp"
#include "RefEdge.hpp"
#include "RefVertex.hpp"

#include "AcisQueryEngine.hpp"
#include "AcisModifyEngine.hpp"
#include "GeometryQueryTool.hpp"
#include "GeometryModifyTool.hpp"
#include "ModelQueryEngine.hpp"
#include "Body.hpp"
#include "BodyACIS.hpp"
#include "CurveACIS.hpp"

// ********** END CUBIT INCLUDES              **********

AcisEdgeTool* AcisEdgeTool::instance_ = 0;

// Method: instance
// provides access to the unique model for this execution.
// sets up this instance on first access
AcisEdgeTool* AcisEdgeTool::instance()
{
  if (instance_ == 0) {
    instance_ = new AcisEdgeTool();
  }
  return instance_;
}

AcisEdgeTool::AcisEdgeTool() 
{
}

AcisEdgeTool::~AcisEdgeTool()
{
}

CubitStatus 
AcisEdgeTool::create_curve_combine( DLIList<Curve*>& curve_list, 
                                    Curve *&new_curve_ptr )
{
   AcisQueryEngine *const AQE = AcisQueryEngine::instance();
   GeometryQueryTool *const GQT = GeometryQueryTool::instance();
   
   new_curve_ptr = create_curve_combine( curve_list );
   if (!new_curve_ptr)
     return CUBIT_FAILURE;
   return CUBIT_SUCCESS;
}   

Curve* AcisEdgeTool::create_curve_combine( DLIList<Curve*>& curve_list )
{
   int i;
   
   DLIList<CurveACIS*> acis_curves(curve_list.size());
   CAST_LIST( curve_list, acis_curves, CurveACIS );
   if (curve_list.size() != acis_curves.size())
   {
     PRINT_ERROR("In AcisEdgeTool::create_curve_combine\n"
                 "       Not all input curves are ACIS Curves.\n");
     return 0;
   }

   AcisQueryEngine *const AQE = AcisQueryEngine::instance();
   
   CurveACIS *ref_edge_ptr1, *ref_edge_ptr2;
   EDGE *EDGE_ptr1, *EDGE_ptr2;
   EDGE *combined_EDGE_ptr = NULL;
   EDGE *new_EDGE_ptr = NULL;

   outcome result;

   acis_curves.reset();
   ref_edge_ptr1 = acis_curves.get_and_step();
   ref_edge_ptr2 = acis_curves.get_and_step();
   
   EDGE_ptr1 = ref_edge_ptr1->get_EDGE_ptr();;
   if( !EDGE_ptr1 )
   {
      PRINT_ERROR( "Unable to get ACIS edge from CUBIT curve -- aborting\n" );
      return 0;
   }
   EDGE_ptr2 = ref_edge_ptr2->get_EDGE_ptr();
   if( !EDGE_ptr2 )
   {
      PRINT_ERROR( "Unable to get ACIS edge from CUBIT curve -- aborting\n" );
      return 0;
   }
   
   result = api_combine_edges( EDGE_ptr1, EDGE_ptr2, combined_EDGE_ptr );

   if( !result.ok() )
   {
      AQE->ACIS_API_error(result);
      PRINT_ERROR( "Error combining Curves\n"  );
      return 0;
   }

   // Loop on remaining edges
   for( i=2; i<acis_curves.size(); i++ )
   {
      new_EDGE_ptr = combined_EDGE_ptr;

      ref_edge_ptr2 = acis_curves.get_and_step();
      EDGE_ptr2 = ref_edge_ptr2->get_EDGE_ptr();;
      if( !EDGE_ptr2 )
      {
         PRINT_ERROR( "Unable to get ACIS edge from CUBIT curve -- aborting\n" );
         api_delent( new_EDGE_ptr );
         return 0;
      }

      result = api_combine_edges( new_EDGE_ptr, EDGE_ptr2, combined_EDGE_ptr );
      if( !result.ok() )
      {
         AQE->ACIS_API_error(result);
         PRINT_ERROR( "Error combining Curves\n"  );
         return 0;
      }

      api_delent( new_EDGE_ptr );
   }

   return AQE->populate_topology_bridges( combined_EDGE_ptr );
}


// KGM
#if 0
CubitStatus 
AcisEdgeTool::split_edges_at_vertices( DLIList<RefEdge*>& ref_edge_list, 
                                       DLIList<RefVertex*>& ref_vertex_list,
                                       DLIList<Body*> &new_body_list,
                                       bool keep_old_body)
{
   DLIList<CubitVector*> positions;
   DLIList<Curve*> curve_list(ref_edge_list.size());
   DLIList<BodySM*> new_sms;
   DLIList<Body*> old_bodies;
   CubitStatus result;
   
   if (!common_setup("split_edges_at_vertices", ref_edge_list, curve_list, old_bodies))
     return CUBIT_FAILURE;
  
   ref_vertex_list.reset();
   for (int i = ref_vertex_list.size(); i--; )
     positions.append( new CubitVector( ref_vertex_list.get_and_step()->coordinates() ) );
   
   result = split_edges_at_vertices( curve_list, positions, new_sms, keep_old_body );  
   while (positions.size())
     delete positions.pop();
   
   if (!GeometryModifyTool::instance()->finish_sm_op(old_bodies, new_sms, new_body_list))
     result = CUBIT_FAILURE;
  
   return result;
}


CubitStatus 
AcisEdgeTool::split_edges_at_vertices( DLIList<Curve*>& ref_edge_list, 
                                       DLIList<CubitVector*>& position_list,
                                       DLIList<BodySM*> &new_body_list,
                                       bool keep_old_body)
{
    AcisQueryEngine  *const AQE = AcisQueryEngine::instance();
    AcisModifyEngine *const AME = AcisModifyEngine::instance();
    outcome result;
    
    BodySM *body_ptr;
    
    BODY *BODY_ptr;
    BODY *copied_BODY_ptr;
    
    int delete_attribs =
        (GeometryModifyTool::instance()->get_new_ids() || keep_old_body);
    
    // we may be able to make this work with mutiple vertices
    // but for now allow just 1
    if (position_list.size() != 1){
        PRINT_ERROR( "Currently only 1 vertex is supported aborting\n");
        return CUBIT_FAILURE;
    }
    CubitVector pos = *position_list.get();
    

    // Copy the incoming ref_edge_list since we will be pulling
    // edges out of it.
    DLIList<CurveACIS*> copied_ref_edge_list(ref_edge_list.size());
    CAST_LIST( ref_edge_list, copied_ref_edge_list, CurveACIS );
    if (ref_edge_list.size() != copied_ref_edge_list.size())
    {
        PRINT_ERROR( "Non-ACIS Curve(s) in AcisEdgeTool::split_edges_at_vertices\n");
        return CUBIT_FAILURE;
    }
    
    copied_ref_edge_list.reset();
    while( copied_ref_edge_list.size() )
    {
        DLIList<EDGE*> EDGE_list;
        if( AME->get_copied_EDGES_of_body( copied_ref_edge_list, EDGE_list, 
            copied_BODY_ptr ) == CUBIT_FAILURE ){
            PRINT_ERROR("Edge Has No Body Associated With It\n");
            break;
        }
        
        // Get original Body and BODY
        body_ptr = AQE->get_body_sm_of_ENTITY( copied_BODY_ptr );
        BODY_ptr = dynamic_cast<BodyACIS*>(body_ptr)->get_BODY_ptr();
        
        // Now cleanout the owner attributes from the copied BODY, if required
        if( delete_attribs )
            AQE->remove_cubit_owner_attrib_in_BODY(copied_BODY_ptr);
        
        EDGE *EDGE_ptr = EDGE_list.get();

        // Make a copy of the vertex coming in to use in the split
        // edge, it will become part of the edge?
        APOINT* P1 = new APOINT (pos.x(),pos.y(),pos.z());
        VERTEX* VERTEX_copy = new VERTEX (P1);

        // Now, split the edges on this body
        // void function no result returned, nice..
        sg_split_edge_at_vertex( EDGE_ptr, VERTEX_copy );

        //Check to see if the vertex has an owner to determine if the operation 
        // was a success?
        if (!VERTEX_copy->owner()){
            PRINT_ERROR( "Unable to split curve \n");
            api_delent(VERTEX_copy);
            // APOINT removed too?
            continue;
        }
        // If we've made it this far, the copied_BODY has been
        // modified and we can update it in CUBIT
        BodySM* new_body_ptr = AME->get_new_Body( body_ptr, BODY_ptr, 
                                      copied_BODY_ptr, keep_old_body );
        
        if (new_body_ptr)
          new_body_list.append( new_body_ptr );
    }
    
    return CUBIT_SUCCESS;
}


CubitStatus AcisEdgeTool::common_setup( const char* name, 
                                   DLIList<RefEdge*>& ref_edge_list,
                                   DLIList<Curve*>& curve_list,
                                   DLIList<Body*>& old_bodies )
{
  AcisQueryEngine *const AQE = AcisQueryEngine::instance();
  ModelQueryEngine *const MQE = ModelQueryEngine::instance();
  GeometryModifyTool *const GMT = GeometryModifyTool::instance();
  
    // Get Body(s) of RefEdges
  DLIList<ModelEntity*> query_input(ref_edge_list.size()), query_output;
  CAST_LIST_TO_PARENT( ref_edge_list, query_input );
  MQE->query_model( query_input, DagType::body_type(), query_output );
  CAST_LIST( query_output, old_bodies, Body );
  
    // Check for virtual geometry
  if ( GMT->contains_intermediate_geometry(old_bodies))
  {
      PRINT_ERROR("%s edges on bodies containing virtual geometry\n"
         "       is not allowed.\n"
         "       Delete virtual geometry on these bodies before operation.\n",
         name);
      return CUBIT_FAILURE;
  }
  
    // Get ACIS Curves
  ref_edge_list.reset();
  for (int i = ref_edge_list.size(); i--; )
  {
    RefEdge* edge_ptr = ref_edge_list.get_and_step();
    TopologyBridge* geom_ptr = edge_ptr->bridge_manager()->topology_bridge(AQE);
    Curve* curve_ptr = dynamic_cast<Curve*>(geom_ptr);
    if (curve_ptr)
      curve_list.append(curve_ptr);
    else
      PRINT_ERROR("Curve %d is not an ACIS Curve.\n", edge_ptr->id());
  }
  
  if (curve_list.size() != ref_edge_list.size())
  {
    PRINT_ERROR("Non-ACIS Curves in AcisEdgeTool::%s\n", name);
    return CUBIT_FAILURE;
  }
  
  return CUBIT_SUCCESS;
}
#endif
