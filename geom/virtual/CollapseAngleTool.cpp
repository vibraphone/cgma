//-------------------------------------------------------------------------
//- Filename:       CollapseAngleTool
//- Purpose:  To collapse small angle for preparing for mesh
//-      "collapse angle less than <angle> mesh_size <length> [preview]\n",
//-      "collapse angle at vertex <id> curve <id> [arc_length <length>] curve
//-       <id> [arc_length <length> | perpendicular | same_size | tangent] 
//-      [preview]\n"
//-
//- Creator:       Jiangtao Hu
//- Creation date: 02/09/2005
//-------------------------------------------------------------------------

// ********** BEGIN STANDARD INCLUDES         **********
// ********** END STANDARD INCLUDES           **********
                                                                                
// ********** BEGIN MOTIF INCLUDES            **********
// ********** END MOTIF INCLUDES              **********
                                                                                
// ********** BEGIN ACIS INCLUDES             **********
// ********** END ACIS INCLUDES               **********
                                                                                
// ********** BEGIN CUBIT INCLUDES            **********

#include "GMem.hpp"                                                           
#include "GeometryQueryTool.hpp"
#include "GeomMeasureTool.hpp"
#include "RefVolume.hpp"
#include "RefEdge.hpp"
#include "RefFace.hpp"
#include "RefVertex.hpp"
#include "CollapseAngleTool.hpp"
#include "Point.hpp"
#include "CubitBox.hpp"
#include "GfxDebug.hpp"
#include "PartitionTool.hpp"
#include "CubitUtil.hpp"
#include "CompositeTool.hpp"
// ********** END CUBIT INCLUDES              **********

CollapseAngleTool *CollapseAngleTool::instance_ = NULL;
// ********** BEGIN PUBLIC FUNCTIONS          **********
CollapseAngleTool *CollapseAngleTool::instance()
{
  if (instance_ == NULL) instance_ = new CollapseAngleTool();
                                                                                
  return instance_;
}

void CollapseAngleTool::collapse_one_angle(RefVertex *vex_ptr,
                                           RefEdge   *edge_to_remove,
                                           RefEdge   *the_other_edge,
                                           double     length1,
                                           double     length2,
                                           CubitBoolean get_position,
                                           CubitVector  &position,
                                           CubitBoolean preview,
                                           CubitBoolean if_comp_vertex,
                                           double       angle)
{
    if (CUBIT_TRUE == preview)
       draw_preview(vex_ptr, edge_to_remove, the_other_edge, length1,
                    length2, get_position, position);
    else
    {
       collapse_angle( vex_ptr, edge_to_remove, the_other_edge,
                       length1, length2, get_position, position);
       if (if_comp_vertex)
       {
          DLIList<RefEdge*> edge_list;
          vex_ptr->ref_edges(edge_list);
          if (edge_list.size() == 2)
          {
             DLIList<RefEdge*> composite_edges;
             angle *= (CUBIT_PI / 180);
             int id = vex_ptr->id();
             CompositeTool::instance()->composite(edge_list, composite_edges,
                                                  NULL, angle);   
             if (composite_edges.size() > 0)
             {
                DLIList<int> id_list;
                id_list.append(id);
                CubitUtil::list_entity_ids("Composite out of vertex ", id_list);
             }
          }
       }
    }
}

CubitStatus CollapseAngleTool::auto_collapse(double length,
                                             double angle,
                                             CubitBoolean preview,
                                             CubitBoolean if_comp_vertex,
                                             double       max_angle)
{
    //clear out the partly_drawn_curve list
    partly_drawn_curve.clean_out();
    //finding all angles less than 'angle'
    DLIList<RefEntity*> entity_list;
    DLIList<RefVolume*> volume_list;
#ifdef BOYD17 
    DLIList<RefVolume*> final_volume_list;
#endif
    RefVolume * ref_volume;
    DLIList <RefEdge *> large_edge_angles_list;
    DLIList <RefEdge *> small_edge_angles_list;
    DLIList <double> large_angles_list;
    DLIList <double> small_angles_list;
    int total_interior;
    int total_fuzzy;
    CubitStatus  result = CUBIT_SUCCESS;

    // get an instance of GeometryQueryTool
    GeometryQueryTool *gqt = GeometryQueryTool::instance();
                                                                                
    volume_list.clean_out();
    result = gqt->ref_entity_list("volume", entity_list);
    if (result == CUBIT_FAILURE || entity_list.size() < 1)
    {
      PRINT_ERROR("This command is prepared for volume mesh only.\n"
                  "The angles must reside on volumes.\n");
      return result;
    }
                                                                                
    CAST_LIST(entity_list, volume_list, RefVolume);
                                                                                
    int i;
    small_edge_angles_list.clean_out();
    for (i=0; i < volume_list.size(); i++)
    {
       ref_volume = volume_list.get_and_step();
       assert (ref_volume != NULL);
       GeomMeasureTool::find_interior_curve_angles(ref_volume,
                                                   360,
                                                   angle,
                                                   large_edge_angles_list,
                                                   small_edge_angles_list,
                                                   large_angles_list,
                                                   small_angles_list,
                                                   total_interior,
                                                   total_fuzzy);
    }
                                                                                
    //collapse small angles, use perpendicular to the second curve.
    int size = small_edge_angles_list.size();
    int init_size = size+1;
    int index = 0;
    RefVertex *vex_ptr = NULL;
    RefEdge   *edge_to_remove = NULL;
    RefEdge   *the_other_edge = NULL;
    double    length1 = -1.0, length2 = -1.0;
    CubitVector  position;

    CubitBoolean if_step_over = CUBIT_FALSE;
    while (size != 0)
    {
       if (init_size <= size)
          break;
       init_size = size;
       small_edge_angles_list.reset();
       assert (size%2 == 0);
       double angle1, angle2;
#ifdef BOYD17 
       DLIList<RefVertex *> vertex_list;
#endif
       if (CUBIT_TRUE == preview || CUBIT_TRUE == if_step_over)
       {
          if_step_over = CUBIT_FALSE;
          for (i=0; i<index; i++)
            small_edge_angles_list.step();
       }
                                                                                
       //find the vertex for each edge pair.
       edge_to_remove = small_edge_angles_list.get_and_step();
       the_other_edge = small_edge_angles_list.get_and_step();
       index += 2;
       assert (edge_to_remove != NULL && the_other_edge != NULL);
       vex_ptr = edge_to_remove->common_ref_vertex(the_other_edge);
       if (vex_ptr == NULL)
       {
          size -= 2;
          if_step_over = CUBIT_TRUE;
          continue;
       }
                                                            
       //check to make sure there's one common refface between the edges.
       if(edge_to_remove->num_of_common_ref_face(the_other_edge) != 1)
       {
          PRINT_INFO("Curves %d and %d are not candidates for collapsing because\n",
                     edge_to_remove->id(), the_other_edge->id());
          PRINT_INFO("the two curves need to have exactly one common face.\n");
          size -= 2;
          if_step_over = CUBIT_TRUE;
          continue;
       }

       //find the true to_be_removed edge
       angle1 = the_surface_angle(edge_to_remove, vex_ptr, length);
       angle2 = the_surface_angle(the_other_edge, vex_ptr, length);
       if (angle1 != -1.0 && angle2 != -1.0 &&
          (fabs(angle1 - CUBIT_PI) > fabs(angle2 - CUBIT_PI)))
          exchange_edges(edge_to_remove, the_other_edge);
       else if (angle1 == -1.0 || angle2 == -1.0 ||
               (fabs(angle1 - CUBIT_PI) == fabs(angle2 - CUBIT_PI)))
       {
          //curved curve has higher priority than straight curve.
          if ((edge_to_remove->geometry_type() == STRAIGHT_CURVE_TYPE)
              &&(the_other_edge->geometry_type() != STRAIGHT_CURVE_TYPE))
             exchange_edges(edge_to_remove, the_other_edge);
                                                                                
          //shorter curve has higher priority than longer curve.
          else if (((edge_to_remove->geometry_type() != STRAIGHT_CURVE_TYPE)
             &&(the_other_edge->geometry_type() != STRAIGHT_CURVE_TYPE)) ||
              ((edge_to_remove->geometry_type() == STRAIGHT_CURVE_TYPE)
              &&(the_other_edge->geometry_type() == STRAIGHT_CURVE_TYPE)))
          {
             if (edge_to_remove->get_arc_length() >
                the_other_edge->get_arc_length())
                exchange_edges(edge_to_remove, the_other_edge);
          }
       }//else if
                                                                                
       // determine the length1
       if (length >= (edge_to_remove->get_arc_length() - GEOMETRY_RESABS))
          length1 = -1.0;
                                                                                
       else
          length1 = length;
       //determine the position on curve2
       RefVertex * tmp_vertex = edge_to_remove->other_vertex(vex_ptr);
       CubitVector tmp_vec = tmp_vertex->coordinates();
       if (length1 != -1.0)
       {
         result = position_from_length(edge_to_remove, vex_ptr, length1, tmp_vec);
       }
                                                                                
       the_other_edge->closest_point(tmp_vec, position);
                                                                                
       //check if position is within the_other_edge
       CubitBox position_box(position);
       CubitBox edge_box = the_other_edge->bounding_box();
                                                                                
       // check if it's a preview command.
       if (CUBIT_TRUE == preview)
       {
          if (!position_box.overlap(0, edge_box))
             draw_preview(vex_ptr, edge_to_remove, the_other_edge,
                          length1, -1.0, CUBIT_FALSE, position);
                                                                                
          else
             draw_preview(vex_ptr, edge_to_remove, the_other_edge,
                          length1, -1.0, CUBIT_TRUE, position);
          size -= 2;
          continue;
       }
                                                                                
       //collapse angle
       else if (!position_box.overlap(0, edge_box))
          result = collapse_angle( vex_ptr, edge_to_remove, the_other_edge,
                                   length1, -1.0, CUBIT_FALSE, position);
       else
          result = collapse_angle( vex_ptr, edge_to_remove, the_other_edge,
                                   length1, length2, CUBIT_TRUE, position);
       if (result != CUBIT_SUCCESS)
       {
          result = CUBIT_SUCCESS;
          continue ;
       }
                                                                                
       if (if_comp_vertex)
       {
          DLIList<RefEdge*> edge_list;
          vex_ptr->ref_edges(edge_list);
          int id = vex_ptr->id();
          if (edge_list.size() == 2)
          {
             DLIList<RefEdge*> composite_edges;
             max_angle *= (CUBIT_PI / 180);
             CompositeTool::instance()->composite(edge_list, composite_edges,
                                                  NULL, max_angle);
             if (composite_edges.size() > 0)
             {
                DLIList<int> id_list;
                id_list.append(id);
                CubitUtil::list_entity_ids("Composite out of vertex ", id_list);
             }
          }
       }

       small_edge_angles_list.clean_out();
       for (i=0; i < volume_list.size(); i++)
       {
          ref_volume = volume_list.get_and_step();
          assert (ref_volume != NULL);
          GeomMeasureTool::find_interior_curve_angles(ref_volume,
                                                      360,
                                                      angle,
                                                      large_edge_angles_list,
                                                      small_edge_angles_list,
                                                      large_angles_list,
                                                      small_angles_list,
                                                      total_interior,
                                                      total_fuzzy);
       }
       size = small_edge_angles_list.size();
   }//while
   return result;
}

// ********** END PUBLIC FUNCTIONS            **********
                                                                                
// ********** BEGIN PROTECTED FUNCTIONS       **********
// ********** END PROTECTED FUNCTIONS         **********
                                                                                
// ********** BEGIN PRIVATE FUNCTIONS         **********
CollapseAngleTool::CollapseAngleTool()
{
}

CollapseAngleTool::~CollapseAngleTool()
{
}

void  CollapseAngleTool::exchange_edges(RefEdge *&edge,
                                        RefEdge *&the_other_edge)
{
   RefEdge* temp_edge;
   temp_edge = edge;
   edge = the_other_edge;
   the_other_edge = temp_edge;
}

double CollapseAngleTool::the_surface_angle(RefEdge *edge,
                                           RefVertex *root_vertex,
                                           double  length)
{
   DLIList <RefFace *> ref_faces;
   double arc_length = length/2;
   double fraction;
                                                                                
   if (length > edge->get_arc_length())
     arc_length = edge->get_arc_length()/2;
                                                                                
   fraction = edge->fraction_from_arc_length(root_vertex, arc_length);
                                                                                
   edge->ref_faces(ref_faces);
   if (ref_faces.size() != 2)
      return -1.0;
                                                                                
   double angle =  GeometryQueryTool::instance()->
                 surface_angle(ref_faces.get_and_step(),ref_faces.get(),
                 edge, NULL, fraction);
   return angle;
}

CubitStatus CollapseAngleTool::position_from_length(RefEdge *edge,
                                                  RefVertex *root_vertex,
                                                  double  arc_length,
                                                  CubitVector& v_new)
{
   CubitStatus result;
   double sense = 1.0;
   if (root_vertex != edge->start_vertex())
      sense = -1.0;
   CubitVector v_root = root_vertex->get_point_ptr()->coordinates();
   result = edge->point_from_arc_length(v_root, arc_length*sense, v_new);
   return result;
}
                                                                                 
//Only change new_vertex when the partition succeeds.
CubitStatus CollapseAngleTool::partition_curve(RefEdge *&edge,
                                              RefVertex *root_vertex,
                                              double  arc_length,
                                              RefVertex *&new_vertex)
{
   CubitStatus  result = CUBIT_SUCCESS;
   CubitVector v_new;
   result = position_from_length(edge, root_vertex, arc_length, v_new);
   if (result == CUBIT_SUCCESS)
      result = partition_curve(edge, root_vertex, v_new, new_vertex);
   return result;
}
                                                                                
                                                                                
//Only change new_vertex when the partition succeeds.
CubitStatus CollapseAngleTool::partition_curve(RefEdge *&edge,
                                              RefVertex *root_vertex,
                                              CubitVector  &position,
                                              RefVertex *&new_vertex)
{
   CubitStatus  result = CUBIT_SUCCESS;
   DLIList<CubitVector*> vect_list;
   vect_list.append(  new CubitVector( position ) );
   DLIList<RefVertex*> vertex_list;
   RefVertex* the_vertex = NULL;
   edge->ref_vertices(vertex_list);
   if (vertex_list.size() == 1)
     the_vertex = vertex_list.get();

   DLIList<RefEdge*> tmp_list;
   result = PartitionTool::instance()->
          partition( edge, vect_list, tmp_list );
   delete vect_list.pop();
   if (!result)
   {
      PRINT_ERROR("Split of curve %d failed.\n", edge->id());
      return result;
   }
//Only change new_vertex when the partition succeeds.
   RefEdge *temp_edge = tmp_list.get_and_step();
   if (NULL != the_vertex)
   {
     new_vertex = temp_edge->other_vertex(the_vertex);
     double len1 = temp_edge->get_arc_length();
     double len2 = tmp_list.get()->get_arc_length();
     if (len1 > len2)
       edge =  tmp_list.get(); 
     else   
       edge = temp_edge;
   }
   else if (temp_edge != edge)
   {
     new_vertex = edge->common_ref_vertex(temp_edge);
     if (root_vertex && edge->start_vertex() != root_vertex && 
         edge->end_vertex() != root_vertex)
       edge = temp_edge;
   }
   else
   {
     new_vertex = edge->common_ref_vertex(tmp_list.get());
     if (root_vertex && edge->start_vertex() != root_vertex &&
         edge->end_vertex() != root_vertex)
       edge = tmp_list.get();
   } 

   return result;
}
                                                                                
CubitStatus CollapseAngleTool::partition_surface(RefFace *common_face,
                                                RefEdge *&edge_to_remove,
                                                RefVertex *root,
                                                RefVertex *vertex1,
                                                RefVertex *vertex2,
                                                RefFace  *&result_face)
{
   CubitStatus  result = CUBIT_SUCCESS;
   DLIList<CubitVector*> positions;
   positions.append(new CubitVector(vertex1->get_point_ptr()->coordinates()));
   positions.append(new CubitVector(vertex2->get_point_ptr()->coordinates()));
                                                                                
   DLIList<RefEdge*> new_edges;
   result_face = PartitionTool::instance()->
      insert_edge( common_face, positions, CUBIT_FALSE, new_edges );
                                                                                
   RefEdge *new_edge = new_edges.get();
   if(new_edge)
   {
     RefVertex * start = new_edge->start_vertex();
     RefVertex * end = new_edge->end_vertex();
     if (start != vertex1 && end != vertex1)
       edge_to_remove = root->common_ref_edge(start == vertex2 ? end : start);
   }

   delete positions.get_and_step();
   delete positions.get();
   if( result_face == NULL )
       result = CUBIT_FAILURE;
                                                                                
   return result;
}

CubitStatus CollapseAngleTool::draw_preview(RefVertex *vex_ptr,
                                           RefEdge   *edge_to_remove,
                                           RefEdge   *the_other_edge,
                                           double     length1,
                                           double     length2,
                                           CubitBoolean get_position,
                                           CubitVector &position)
{
   CubitStatus result = CUBIT_SUCCESS;
 
   //flush out previous drawing
   GfxDebug::display_all();
                                                                                
   //wireframe display of future composite surface per design review.
   DLIList<RefFace*> ref_faces;
   edge_to_remove-> ref_faces(ref_faces);
   RefFace *common_face = edge_to_remove->common_ref_face(the_other_edge);
   assert(ref_faces.size() == 2);
                                                                                
   RefFace *highlight_surf = ref_faces.get_and_step();
   if (highlight_surf == common_face)
       highlight_surf = ref_faces.get();
                                                                                
   DLIList<RefEdge*> ref_edges;
   highlight_surf->ref_edges(ref_edges);
   ref_edges.remove(edge_to_remove);
                                                                                
   int i;
   for (i = 0; i < partly_drawn_curve.size(); i++)
      ref_edges.remove(partly_drawn_curve.get_and_step());

   for(i = 0; i < ref_edges.size(); i++)
      GfxDebug::draw_ref_entity(ref_edges.get_and_step(), CUBIT_BLUE);
                                                                                
   //check if there will be a cutting edge
   CubitBoolean if_partition = if_partition_surf(vex_ptr, edge_to_remove,
                                               common_face);
   if (length1 == -1.0 && !get_position && length2 == -1.0 && !if_partition)
   {
      ref_edges.clean_out();
      common_face->ref_edges(ref_edges);
      ref_edges.remove(edge_to_remove);
                                                                                
      for (i = 0; i < partly_drawn_curve.size(); i++)
         ref_edges.remove(partly_drawn_curve.get_and_step());

      for(i = 0; i < ref_edges.size(); i++)
        GfxDebug::draw_ref_entity(ref_edges.get_and_step(), CUBIT_BLUE);
      return result;
   }
                                                                                
   CubitVector location, location2, location3;
   if (length1 == -1.0)
       location = edge_to_remove->other_vertex(vex_ptr)->coordinates();
                                                                                
   if (length1 != -1.0)
       result = position_from_length(edge_to_remove, vex_ptr,
                                       length1, location);
                                                                                
   if ( !get_position && length2 == -1.0 &&
        (length1 != -1.0 ||
         (length1 == -1.0 && if_partition)))
      location2=the_other_edge->other_vertex(vex_ptr)->coordinates();
                                                                                
   else if(get_position)
      location2 = position;
                                                                                
   else if (length2 != -1.0)
      result = position_from_length(the_other_edge, vex_ptr, length2, location2);
                                                                                
   GfxDebug::draw_line(location, location2, CUBIT_BLUE);
                                                                                
   //Draw part of curves
   double length;
   CubitVector temp_vex;
   GPoint pts[5];
   location3 = vex_ptr->coordinates();

   partly_drawn_curve.append(the_other_edge);
   partly_drawn_curve.append(edge_to_remove);
   if (the_other_edge->geometry_type() == STRAIGHT_CURVE_TYPE)
      GfxDebug::draw_line(location3, location2, CUBIT_BLUE);
                                                                                
   else
   {
      length = the_other_edge->get_arc_length(location2, location3);
      pts[0].x = location3.x();
      pts[0].y = location3.y();
      pts[0].z = location3.z();
      pts[4].x = location2.x();
      pts[4].y = location2.y();
      pts[4].z = location2.z();
      position_from_length(the_other_edge, vex_ptr, length/4, temp_vex);
      pts[1].x = temp_vex.x();
      pts[1].y = temp_vex.y();
      pts[1].z = temp_vex.z();
      position_from_length(the_other_edge, vex_ptr, length/2, temp_vex);
      pts[2].x = temp_vex.x();
      pts[2].y = temp_vex.y();
      pts[2].z = temp_vex.z();
      position_from_length(the_other_edge, vex_ptr, length*0.75, temp_vex);
      pts[3].x = temp_vex.x();
      pts[3].y = temp_vex.y();
      pts[3].z = temp_vex.z();
      GfxDebug::draw_polyline(pts, 5, CUBIT_BLUE);
   }
                                                                                
   RefVertex *vertex = edge_to_remove->other_vertex(vex_ptr);
   location3 = vertex->coordinates();
   if (edge_to_remove->geometry_type() == STRAIGHT_CURVE_TYPE)
      GfxDebug::draw_line(location3, location, CUBIT_BLUE);
                                                                                
   else
   {
      length = edge_to_remove->get_arc_length(location, location3);
      pts[0].x = location3.x();
      pts[0].y = location3.y();
      pts[0].z = location3.z();
      pts[4].x = location.x();
      pts[4].y = location.y();
      pts[4].z = location.z();
                                                                                
      position_from_length(edge_to_remove, vertex, length/4, temp_vex);
      pts[1].x = temp_vex.x();
      pts[1].y = temp_vex.y();
      pts[1].z = temp_vex.z();
      position_from_length(edge_to_remove, vertex, length/2, temp_vex);
      pts[2].x = temp_vex.x();
      pts[2].y = temp_vex.y();
      pts[2].z = temp_vex.z();
      position_from_length(edge_to_remove, vertex, length*0.75, temp_vex);
      pts[3].x = temp_vex.x();
      pts[3].y = temp_vex.y();
      pts[3].z = temp_vex.z();
      GfxDebug::draw_polyline(pts, 5, CUBIT_BLUE);
   }
                                                                                
   GfxDebug::flush();
                                                                                
   return result;
}
                                                                                
CubitBoolean CollapseAngleTool::if_partition_surf(RefVertex *vex_ptr,
                                                  RefEdge  *edge_to_remove,
                                                  RefFace  *common_face)
{
   CubitBoolean partition_surface_ =CUBIT_FALSE;
   RefVertex *temp_vtx = NULL;
   RefEdge * temp_edge = NULL;
   if (edge_to_remove->start_vertex() == vex_ptr)
       temp_vtx = edge_to_remove->end_vertex();
   else
       temp_vtx = edge_to_remove->start_vertex();
                                                                                
   temp_edge = edge_to_remove->get_other_curve(temp_vtx, common_face);
                                                                                
   if (temp_edge->num_of_common_ref_face(edge_to_remove) > 1)
       partition_surface_ = CUBIT_TRUE;
                                                                                
   return partition_surface_;
}

CubitStatus CollapseAngleTool::collapse_angle(RefVertex *vex_ptr,
                                             RefEdge   *edge_to_remove,
                                             RefEdge   *the_other_edge,
                                             double     length1,
                                             double     length2,
                                             CubitBoolean get_position,
                                             CubitVector &position)
{
  CubitStatus result = CUBIT_SUCCESS;
  RefVertex *new_vertex1 = edge_to_remove->other_vertex(vex_ptr);
  RefVertex *new_vertex2 = the_other_edge->other_vertex(vex_ptr);
  DLIList<RefFace*> ref_faces;
  edge_to_remove-> ref_faces(ref_faces);
  RefFace *common_face = NULL;
  DLIList<RefFace*> common_face_list;
  int nfaces = edge_to_remove->num_of_common_ref_face(the_other_edge);

  if (nfaces == 1)
    common_face = edge_to_remove->common_ref_face(the_other_edge);
  else
  {
    assert (nfaces == 2);
    edge_to_remove->common_ref_faces(the_other_edge, common_face_list);
  }

  if (length1 != -1.0)
  {
    result = partition_curve(edge_to_remove, vex_ptr, length1, new_vertex1);
    if (result != CUBIT_SUCCESS)
    {
       PRINT_ERROR("Operation failed at vertex %d, curve %d and curve %d.\n",
                vex_ptr->id(), edge_to_remove->id(), the_other_edge->id());
       return result;
    }
  }
                                                                                
  if (get_position)
  {
    result = partition_curve(the_other_edge, NULL, position, new_vertex2);
    if (result != CUBIT_SUCCESS)
    {
       PRINT_ERROR("Operation failed at vertex %d, curve %d and curve %d.\n",
                vex_ptr->id(), edge_to_remove->id(), the_other_edge->id());
       return result;
    }
  }
                                                                                
  else if (!get_position && length2 != -1.0)
  {
    result = partition_curve(the_other_edge, vex_ptr, length2, new_vertex2);
    if (result != CUBIT_SUCCESS)
    {
       PRINT_ERROR("Operation failed at vertex %d, curve %d and curve %d.\n",
                vex_ptr->id(), edge_to_remove->id(), the_other_edge->id());
       return result;
    }
  }
                                                                                
  if (nfaces > 1) //determine the common face.
  {
     CubitBox edge_box(new_vertex1->coordinates());
     edge_box |= new_vertex2->bounding_box();
     CubitBox face_box1, face_box2;
     RefFace *face1, *face2;
     face1 = common_face_list.get_and_step();
     face2 = common_face_list.get();
     face_box1 = face1->bounding_box();
     face_box2 = face2->bounding_box();
     if(edge_box.overlap(0, face_box1))
     {
       // Face 1 overlaps with the partition points.
       if(edge_box.overlap(0, face_box2))
       {
         // Both face 1 and face 2 overlap with the parition points.
         // We must do a further check to try to choose the best one.
         // Take the mid point between the two vertices and calculate
         // the distance from both faces.  Choose the one with the 
         // smaller distance.  This approach was chosen because there
         // was a case where one of the potential common faces was
         // actually a composite that had a 90 degree turn in it
         // so that the end points of the edge that would be used
         // to partition it fell on the composite surface but the
         // rest of the interior of the edge didn't lie on the surface.
         // However, the edge completely lied on the other potential
         // surface and so it was the better candidate.
         CubitVector mid_pos = (new_vertex1->coordinates() +
                new_vertex2->coordinates())/2;
         Surface *surf1 = face1->get_surface_ptr();
         Surface *surf2 = face2->get_surface_ptr();
         if(surf1 && surf2)
         {
           CubitVector pos1, pos2;
           surf1->closest_point(mid_pos, &pos1);
           surf2->closest_point(mid_pos, &pos2);
           if((mid_pos - pos1).length_squared() <
             (mid_pos - pos2).length_squared())
           {
             common_face = face1;
           }
           else
           {
             common_face = face2;
           }
         }
         else
         {
            PRINT_ERROR("Operation failed because one of the potential common faces doesn't have a surface.\n");
            return CUBIT_FAILURE;
         }
       }
       else
       {
         // Only face 1 overlaps so use it as the common face.
         common_face = face1;
       }
     }
     else if(edge_box.overlap(0, face_box2))
     {
       // Only face 2 overlaps so use it as the common face.
       common_face = face2;
     }
     else
     {
       // Neither face overlaps so throw an error.
       PRINT_ERROR("Operation failed to find a valid common face for partitioning.\n");
       return CUBIT_FAILURE;
     }
  }

  CubitBoolean partition_surface_ = CUBIT_FALSE;
  // check if edge_to_remove is the only one that shares the two surfaces.
  if (length1 == -1.0)
     partition_surface_ = if_partition_surf(vex_ptr, edge_to_remove,
                                            common_face);
                                                                                
  RefFace * face_to_remove = common_face;
  if (length1 != -1.0 || get_position || length2 != -1.0 || partition_surface_)
  {
     RefFace *result_face;
     result = partition_surface(common_face, edge_to_remove, vex_ptr, 
                                new_vertex1, new_vertex2, result_face);
     if (result != CUBIT_SUCCESS)
     {
       PRINT_ERROR("Split of surface %d failed.\n", common_face->id());
       PRINT_ERROR("Operation failed at vertex %d, curve %d and curve %d.\n",
                vex_ptr->id(), edge_to_remove->id(), the_other_edge->id());
       return result;
     }
                                                                                
     DLIList<RefEdge*> ref_edges;
     result_face->ref_edges(ref_edges);
     int i;
     for (i = 0; i <ref_edges.size(); i++)
        if (ref_edges.get_and_step() == edge_to_remove)
        {
           face_to_remove = result_face;
           break;
        }
  }
                                                                                
  if (face_to_remove != common_face)
  {
    ref_faces.remove(common_face);
    ref_faces.append(face_to_remove);
  }
                                                                                
  DLIList<RefFace*> result_faces;
  CompositeTool::instance()->composite(ref_faces, result_faces);
                                                                                
  if( result_faces.size() > 0 )
  {
     DLIList<int> id_list;
     int e;
     for(  e = 0; e < result_faces.size(); e++)
     {
        id_list.append( result_faces.get_and_step()->id() );
        CubitUtil::list_entity_ids("Created composite surfaces ", id_list );
     }
  }
                                                                                
  else
  {
     PRINT_ERROR("Operation failed at vertex %d, curve %d and curve %d.\n",
                vex_ptr->id(), edge_to_remove->id(), the_other_edge->id());
     result = CUBIT_FAILURE;
  }

  //composite out new_vertex1
  DLIList<RefEdge*> edge_list;
  RefEdge* composite_edge = NULL;
  int id = new_vertex1->id();
  new_vertex1->ref_edges(edge_list);
  if (edge_list.size() == 2)
    composite_edge = CompositeTool::instance()->composite(edge_list);
  if (composite_edge)
  {
     DLIList<int> id_list;
     id_list.append(id);
     CubitUtil::list_entity_ids("Composite out of vertex ", id_list);
  }

  return result;
}             
