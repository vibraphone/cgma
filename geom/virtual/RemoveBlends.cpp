#include "RemoveBlends.hpp"
#include "CubitUtil.hpp"
//#include "Model.hpp"
#include "GeometryQueryTool.hpp"

//#include "SVDrawTool.hpp"
#include "GMem.hpp"
#include "GfxDebug.hpp"

#include "VirtualQueryEngine.hpp"
#include "PartitionEngine.hpp"


#include "CompositeTool.hpp"
#include "PartitionTool.hpp"
#include "PartSurfFacetTool.hpp"

#include "RefFace.hpp"
#include "RefEdge.hpp"
#include "RefVertex.hpp"

#include "CompositeLump.hpp"
#include "PartitionCurve.hpp"
#include "SubSurface.hpp"
#include "SplitSurfaceTool.hpp"

CubitStatus RemoveBlends::remove_blend(RefFace* ref_face,
                                          DLIList<CubitVector*>& locations,
                                          DLIList<RefFace*>& composite_faces)
{
   CubitStatus result = CUBIT_SUCCESS;
   int i;

   locations.reset();
   RefFace* the_face = ref_face;
   if (the_face && VirtualQueryEngine::is_virtual(the_face))
   {
      PRINT_WARNING ("Can not collapse a virtual surface.\n");
      return CUBIT_FAILURE;
   }

   // Find the curve to split along.
   SplitSurfaceTool                sst;
   DLIList<DLIList<CubitVector*>*> vec_lists;
   DLIList<Curve*>                 curve_list;

   DLIList<CubitVector*> *vec_list = new DLIList<CubitVector*>;

   vec_list->append( new CubitVector( *locations.get_and_step() ) );
   vec_list->append( new CubitVector( *locations.get_and_step() ) );

   vec_lists.append( vec_list );

   sst.calculate_split_curves(the_face, locations, vec_lists, curve_list);

   //int num_segs = 2; 
   //double fraction =.5;
   //double distance = -1.;
   //RefEdge* from_curve_ptr = 0;
   //DLIList<RefVertex*> corner_vertex_list;
   //DLIList<RefVertex*> through_vertex_list;
   //RefEdge *curve_dir_ptr = 0;
   //CubitBoolean preview_flg = false;
   //CubitBoolean create_ref_edges_flg = false;
   //CubitBoolean just_curves_flg = false;
   //DLIList<DLIList<Curve*>*> curve_lists_list; 

//   sst.calculate_split_curves(ref_face_list, num_segs, fraction, distance, 
//                          from_curve_ptr, corner_vertex_list, through_vertex_list,
//                          curve_dir_ptr,
//                          preview_flg, create_ref_edges_flg, just_curves_flg,
//                          curve_lists_list ); 

   //partition the_face
   DLIList<RefEdge*> new_edges;
   DLIList<RefFace*> result_faces;
   result = PartitionTool::instance()->
          partition_face_by_curves( the_face, curve_list, result_faces, CUBIT_TRUE, &new_edges );
   if( result == CUBIT_FAILURE )
   {
      PRINT_ERROR("Failed to partition surface %d into 2 distinct pieces\n", the_face->id());
      return result;
   }

   // We MUST clean up an ambiguous partition.
   if (new_edges.size() > 1)
   {
      PRINT_ERROR("Attempted to partition surface %d more than 2 parts,\n"
                    "Ambiguous surface collapse. Collapse failed.\n",
                    the_face->id());
      DLIList<RefFace*> orig_faces;
      PartitionTool::instance()->unpartitionAll(result_faces, orig_faces);
      
      DLIList<RefEdge*> edges, edge_results;
      edges.clean_out();
      // We started out with one face, we should the same face back
      assert( the_face == orig_faces.next(0) );
      the_face->ref_edges( edges );

      // Do a quick sort to find the PartitionCurves
      for ( i = edges.size(); i--; ) 
      {
        if ( !dynamic_cast<PartitionCurve*>(edges.step_and_get()->get_curve_ptr()))
          edges.change_to(0);
      }
      edges.remove_all_with_value(0);

      // Now remove the partitioned edges
      // TODO: what if we are collapsing with an adjacent partitioned surface.
      // What happens to the partitioned edges?
      if ( edges.size() )
        PartitionTool::instance()->unpartitionAll(edges, edge_results);

      return CUBIT_FAILURE;
   }

   if (result_faces.size() < 2)
   {
      PRINT_WARNING("Surface %d has been improperly partitioned,\n"
                    "Ambiguous surface collapse. Collapse failed.\n",
                    the_face->id());
      return result;
   }
 
   //if(new_face == the_face)
   //   new_face = GeometryQueryTool::instance()->get_last_ref_face ();   

   new_edges.clean_out();
   //Re-composite surfaces
   DLIList<RefFace*> composite_result_faces;
   DLIList<RefFace*> remaining_faces;
   DLIList<RefFace*> surf_list;
   DLIList<int> id_list;
   RefFace*  face_to_composite = NULL;
   RefFace*  temp_face = NULL;
   RefFace*  second_face = NULL;
   RefFace*  composite_face = NULL;
   int  num_common_edge1, num_common_edge2;
   int  size;
   int  composite_count = 0;
   while (0 < (composite_faces.size() + remaining_faces.size()))
   {
      if(composite_count == 2)
        break;

      size = composite_faces.size();
      if( size > 0)
      {
        composite_faces.reset();
        face_to_composite = composite_faces.get();
        composite_faces.remove(face_to_composite);
      }

      else
      {
        remaining_faces.reset();
        face_to_composite = remaining_faces.get();
        remaining_faces.remove(face_to_composite);
      }
      
      if (NULL != second_face)
      {
         num_common_edge1 = second_face->common_ref_edges(face_to_composite,
                                                       new_edges); 
         new_edges.clean_out();
         if (num_common_edge1 == 0)
           continue;
         
         if (NULL == composite_face)
           num_common_edge2 = 0;
         else
           num_common_edge2 = composite_face->
                   common_ref_edges(face_to_composite, new_edges);
         new_edges.clean_out();
      }

      else 
      {
         if (NULL == result_faces[0])
           num_common_edge1 = 0;
         else 
         {
            num_common_edge1 = result_faces[0]->common_ref_edges(face_to_composite, 
                                                       new_edges);  
            new_edges.clean_out();
         }

         if(NULL == result_faces[1])
            num_common_edge2 = 0;
         else
         {
            num_common_edge2 = result_faces[1]->common_ref_edges(face_to_composite,
                                                          new_edges);
            new_edges.clean_out();
         }
      }

      if (num_common_edge1 == 0 && num_common_edge2 == 0)
         continue;

      else if(num_common_edge1 == 0 || num_common_edge2 == 0 || size == 0)
      {  
         surf_list.clean_out();
         if (NULL != second_face)
            surf_list.append(second_face);

         else if (num_common_edge2 == 0 || size == 0)        
         {
            temp_face = result_faces[0];
            second_face = result_faces[1];
            surf_list.append(temp_face);
         }

         else
         {
            temp_face = result_faces[1];
            second_face = result_faces[0];
            surf_list.append(temp_face);
         }
         surf_list.append(face_to_composite);
     
         // composite the surfaces together
         CompositeTool::instance()->composite(surf_list, composite_result_faces);
         composite_count++;

         if( composite_result_faces.size() > 0 )
         {
            DLIList<RefEdge*> surface_edges, new_edge_list;
            remaining_faces.clean_out();
            for(  i = 0; i < composite_result_faces.size(); i++)
            {
               // now composite together the edges of the new surface
               surface_edges.clean_out();
               composite_result_faces.get_and_step()->ref_edges( surface_edges );
               new_edge_list.clean_out();
               CompositeTool::instance()->composite( surface_edges, new_edge_list);

               // keep the lists coherent
               composite_face = composite_result_faces.pop();
               composite_faces.remove(composite_face);
               id_list.append( composite_face->id() );
            }
         }
                                                           
         else
            PRINT_ERROR("Composite surface %d and surface %d failed.\n",
                temp_face->id(), face_to_composite->id());
      }

      else 
         remaining_faces.append(face_to_composite);      
   }
   if (id_list.size() > 0)
      CubitUtil::list_entity_ids("Created composite surfaces ", id_list );
   return result;
}

CubitStatus RemoveBlends::remove_blends(DLIList<RefFace*>& ref_face_list,
                                        int num_segs, double fraction, double distance,
                                        RefEdge* from_curve_ptr,
                                        DLIList<RefVertex*>& corner_vertex_list,
                                        DLIList<RefVertex*>& through_vertex_list,
                                        RefEdge *curve_dir_ptr,
                                        CubitBoolean preview_flg,
                                        DLIList<CubitVector*>& locations,
                                        DLIList<RefFace*>& composite_faces)
{
   CubitStatus result = CUBIT_SUCCESS;
   int i;

   if (ref_face_list.size() <= 1)
   {
      PRINT_WARNING ("Must specify a chain of surfaces.\n");
      return CUBIT_FAILURE;
   }
   locations.reset();

   // check and make sure that we aren't splitting virtual surfaces (yet)
   for (i = 0; i < ref_face_list.size(); i++)
   {
     RefFace* the_face = ref_face_list.get_and_step();
     if (the_face && VirtualQueryEngine::is_virtual(the_face))
     {
        PRINT_WARNING ("Can not collapse a virtual surface.\n");
        return CUBIT_FAILURE;
     }
   }

   // Find the curve to split along.
   SplitSurfaceTool                sst;
   DLIList<DLIList<CubitVector*>*> vec_lists;
   DLIList<Curve*>*                curve_list;
   DLIList<DLIList<Curve*>*> curve_lists_list; 

   DLIList<CubitVector*> *vec_list = new DLIList<CubitVector*>;

   vec_list->append( new CubitVector( *locations.get_and_step() ) );
   vec_list->append( new CubitVector( *locations.get_and_step() ) );

   vec_lists.append( vec_list );

   CubitBoolean create_ref_edges_flg = false;
   CubitBoolean just_curves_flg = false;

   sst.calculate_split_curves(ref_face_list, num_segs, fraction, distance, 
                          from_curve_ptr, corner_vertex_list, through_vertex_list,
                          curve_dir_ptr,
                          preview_flg, create_ref_edges_flg, just_curves_flg,
                          curve_lists_list ); 

   //partition the_face
   DLIList<RefEdge*> new_edges;
   DLIList<RefFace*> result_faces;
   RefFace* the_face;
   DLIList<int> id_list;

   for (i = 0; i < ref_face_list.size(); i++)
   {
     // get the face and the curves that split this face
     the_face = ref_face_list.get_and_step();
     curve_list = curve_lists_list.get_and_step();

     result = PartitionTool::instance()->
            partition_face_by_curves( the_face, *curve_list, result_faces, CUBIT_FALSE, &new_edges );
     if( result == CUBIT_FAILURE )
        return result;

     if (new_edges.size() > 1)
     {
        PRINT_WARNING("Surface %d has been partitioned into more than 2 parts,\n"
                      "Ambiguous surface collapse.\n",
                      the_face->id());
        return result;
     }
 
     new_edges.clean_out();
     //Re-composite surfaces
     DLIList<RefFace*> composite_result_faces;
     DLIList<RefFace*> remaining_faces;
     DLIList<RefFace*> surf_list;
     RefFace*  face_to_composite = NULL;
     RefFace*  temp_face = NULL;
     RefFace*  second_face = NULL;
     RefFace*  composite_face = NULL;
     int  num_common_edge1, num_common_edge2;
     int  size;
     int  composite_count = 0;
     while (0 < (composite_faces.size() + remaining_faces.size()))
     {
        if(composite_count == 2)
          break;

        size = composite_faces.size();
        if( size > 0)
        {
          composite_faces.reset();
          face_to_composite = composite_faces.get();
          composite_faces.remove(face_to_composite);
        }
        else
        {
          remaining_faces.reset();
          face_to_composite = remaining_faces.get();
          remaining_faces.remove(face_to_composite);
        }
        
        if (NULL != second_face)
        {
           num_common_edge1 = second_face->common_ref_edges(face_to_composite,
                                                         new_edges); 
           new_edges.clean_out();
           if (num_common_edge1 == 0)
             continue;
           
           if (NULL == composite_face)
             num_common_edge2 = 0;
           else
             num_common_edge2 = composite_face->
                     common_ref_edges(face_to_composite, new_edges);
           new_edges.clean_out();
        }
        else 
        {
           if (NULL == result_faces[0])
             num_common_edge1 = 0;
           else 
           {
              num_common_edge1 = result_faces[0]->common_ref_edges(face_to_composite, 
                                                         new_edges);  
              new_edges.clean_out();
           }

           if(NULL == result_faces[1])
              num_common_edge2 = 0;
           else
           {
              num_common_edge2 = result_faces[1]->common_ref_edges(face_to_composite,
                                                            new_edges);
              new_edges.clean_out();
           }
        }

        if (num_common_edge1 == 0 && num_common_edge2 == 0)
           continue;

        else if(num_common_edge1 == 0 || num_common_edge2 == 0 || size == 0)
        {  
           surf_list.clean_out();
           if (NULL != second_face)
              surf_list.append(second_face);

           else if (num_common_edge2 == 0 || size == 0)        
           {
              temp_face = result_faces[0];
              second_face = result_faces[1];
              surf_list.append(temp_face);
           }

           else
           {
              temp_face = result_faces[1];
              second_face = result_faces[0];
              surf_list.append(temp_face);
           }
           surf_list.append(face_to_composite);
       
           // composite the surfaces together
           CompositeTool::instance()->composite(surf_list, composite_result_faces);
           composite_count++;

           if( composite_result_faces.size() > 0 )
           {
              DLIList<RefEdge*> surface_edges, new_edge_list;
              remaining_faces.clean_out();
              for(  i = 0; i < composite_result_faces.size(); i++)
              {
                 // now composite together the edges of the new surface
                 surface_edges.clean_out();
                 composite_result_faces.get_and_step()->ref_edges( surface_edges );
                 new_edge_list.clean_out();
                 CompositeTool::instance()->composite( surface_edges, new_edge_list);

                 // keep the lists coherent
                 composite_face = composite_result_faces.pop();
                 composite_faces.remove(composite_face);
                 id_list.append( composite_face->id() );
              }
           }
           PRINT_ERROR("Composite surface %d and surface %d failed.\n",
                temp_face->id(), face_to_composite->id());
       }
     }
   }
   if (id_list.size() > 0)
      CubitUtil::list_entity_ids("Created composite surfaces ", id_list );
   return result;
}

