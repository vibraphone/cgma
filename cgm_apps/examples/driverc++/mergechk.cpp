/**
 * \file mergechk.cpp
 *
 * \brief mergechk, another simple C++ driver for CGM
 *
 * This program acts as a simple driver for CGM.  It reads in a geometry,
 * performs imprints between all the bodies, merges them, and writes information
 * on the results.  It also performs pairwise intersections between the
 * bodies to check for overlaps.  Results are written to stardard output.
 *
 */

#include "GeometryModifyTool.hpp"
#include "GeometryQueryTool.hpp"
#include "AcisQueryEngine.hpp"
#include "MergeTool.hpp"
#include "CubitUtil.hpp"
#include "CubitMessage.hpp"
#include "CubitDefines.h"
#include "RefEntity.hpp"
#include "Body.hpp"
#include "RefVolume.hpp"
#include "RefFace.hpp"
#include "RefEdge.hpp"
#include "RefVertex.hpp"
#include "CubitObserver.hpp"
#include "CastTo.hpp"
#include "AcisQueryEngine.hpp"
#include "AcisModifyEngine.hpp"
#include "AppUtil.hpp"
#include "RefEntityFactory.hpp"

// forward declare some functions used and defined later
CubitStatus read_geometry(int, char **);
CubitStatus evaluate_overlaps();
CubitStatus imprint_bodies();
CubitStatus print_unmerged_surfaces();

// macro for printing a separator line
#define PRINT_SEPARATOR   PRINT_INFO("=======================================\n");


/// main program - initialize, then send to proper function
int main (int argc, char **argv)
{

  CubitObserver::init_static_observers();
    // Initialize the GeometryTool
  
  CGMApp::instance()->startup( argc, argv );
  GeometryQueryTool *gti = GeometryQueryTool::instance();
  AcisQueryEngine::instance();
  AcisModifyEngine::instance();

    // If there aren't any file arguments, print usage and exit
  if (argc == 1) {
    PRINT_INFO("Usage: mergechk <geom_file> [<geom_file> ...]\n");
    exit(0);
  }
  
    // Read in the geometry from files specified on the command line
  CubitStatus status = read_geometry(argc, argv);
  if (status == CUBIT_FAILURE) exit(1);
  else if (gti->num_bodies() == 0) {
    PRINT_WARNING("No bodies read; exiting.\n");
    int ret_val = ( CubitMessage::instance()->error_count() );

    exit(ret_val);
  }

    // Check for overlaps
  status = evaluate_overlaps();
  if (status == CUBIT_FAILURE) exit(1);

    // Imprint bodies together, reporting on results
  status = imprint_bodies();
  if (status == CUBIT_FAILURE) exit(1);
  
    // Merge bodies
  status = MergeTool::instance()->merge_all_bodies();
  if (status == CUBIT_FAILURE) exit(1);
  
    // Print number and ids of non-shared surfaces
  status = print_unmerged_surfaces();
  if (status == CUBIT_FAILURE) exit(1);

  int ret_val = ( CubitMessage::instance()->error_count() );
  if ( ret_val > 0 )
  {
    PRINT_ERROR("Errors found during Mergechk session.\n");
  }
  return ret_val;
  
}

/// attribs module: list, modify attributes in a give model or models
/// 
/// Arguments: file name(s) of geometry files in which to look
///
CubitStatus read_geometry(int num_files, char **argv) 
{
  CubitStatus status = CUBIT_SUCCESS;
  GeometryQueryTool *gti = GeometryQueryTool::instance();
  assert(gti);
  int i;
  
    // For each file, open and read the geometry
  FILE *file_ptr;

  PRINT_SEPARATOR;

  for (i = 1; i < num_files; i++) {
    status = gti->import_solid_model(argv[i], "ACIS_SAT");
    if (status != CUBIT_SUCCESS) {
      PRINT_ERROR("Problems reading geometry file %s.\n", argv[i]);
    }
  }
  PRINT_SEPARATOR;

  return CUBIT_SUCCESS;
}

CubitStatus evaluate_overlaps() 
{
    // evaluate overlaps by intersecting bodies pairwise
  int i;
  GeometryQueryTool *gqti = GeometryQueryTool::instance();
  assert(gqti);
  GeometryModifyTool *gmti = GeometryModifyTool::instance();
  assert(gmti);
  
    // make a copy of the body list for use in this function
  DLIList<Body*> all_bodies, all_new_bodies;
  gqti->bodies(all_bodies);
  
    // step backward on this list, extracting the last body and using it as a tool
    // for remaining bodies
  
  for (i = all_bodies.size(); i > 1; i--) {
    all_bodies.last();
    Body *tool = all_bodies.remove();
    
      // intersect the tool with remaining bodies; make sure and keep old bodies,
      // since we're not using copies; save new bodies for evaluation later
    DLIList<Body*> new_bodies;
    DLIList<Body*> temp_bodies = all_bodies;
    
    CubitStatus status = gmti->intersect(tool, temp_bodies, new_bodies, CUBIT_TRUE);
    if (CUBIT_FAILURE == status) return status;
    all_new_bodies += new_bodies;
  }
  
    // count number of geometric entities in new bodies; if there are no overlaps,
    // this number should be zero
  DLIList<RefEntity*> child_entities, temp_children;

    // first get all child entities of the new bodies
  for (i = all_new_bodies.size(); i > 0; i--) {
    all_new_bodies.get_and_step()->get_all_child_ref_entities(temp_children);
    child_entities += temp_children;
    temp_children.clean_out();
  }
  
    // then filter the list, keeping only unique entities
  temp_children.clean_out();
  temp_children.merge_unique(child_entities);
  
    // now report
  PRINT_SEPARATOR;
  if (temp_children.size() == 0)
    PRINT_INFO("No body overlaps.\n");

  else {
      // check for various types of entities; check by casting the list
      // to that type
    PRINT_INFO("Overlaps found: \n");
    DLIList<RefEntity*> temp_entities;
    
      // check vertices
      // temp_entities is cleaned out in the CAST_LIST macro, so no need
      // to do that here
    CAST_LIST(temp_children, temp_entities, RefVertex);
    if (temp_entities.size() > 0) 
      PRINT_INFO("   Vertices: %d\n", temp_entities.size());
    
      // check edges
    CAST_LIST(temp_children, temp_entities, RefEdge);
    if (temp_entities.size() > 0) 
      PRINT_INFO("   Edges: %d\n", temp_entities.size());
    
      // check faces
    CAST_LIST(temp_children, temp_entities, RefFace);
    if (temp_entities.size() > 0) 
      PRINT_INFO("   Faces: %d\n", temp_entities.size());
    
      // check volumes
    CAST_LIST(temp_children, temp_entities, RefVolume);
    if (temp_entities.size() > 0) 
      PRINT_INFO("   Volumes: %d\n", temp_entities.size());

      // now delete all the bodies produced by the intersections
    DLIList<Body*> new_bodies;
    CAST_LIST(temp_children, new_bodies, Body);
    for (i = new_bodies.size(); i > 0; i--) 
      gqti->delete_Body(new_bodies.get_and_step());
  }
  PRINT_SEPARATOR;
  
    // we're done
  return CUBIT_SUCCESS;
}

CubitStatus imprint_bodies() 
{
    // imprint all the bodies together, and report number of new 
    // entities formed
    //
  GeometryQueryTool *gti = GeometryQueryTool::instance();
  GeometryModifyTool *gmti = GeometryModifyTool::instance();

    // first, count old entities
  int num_vertices = gti->num_ref_vertices();
  int num_edges = gti->num_ref_edges();
  int num_faces = gti->num_ref_faces();
  int num_volumes = gti->num_ref_volumes();
  int num_bodies = gti->num_bodies();
  
    // imprint the bodies together, discarding old bodies
  DLIList<Body*> old_bodies, new_bodies;
  gti->bodies(old_bodies);
  gmti->imprint(old_bodies, new_bodies);
  
    // now count new numbers of entities, subtracting old numbers
  num_vertices = gti->num_ref_vertices() - num_vertices;
  num_edges = gti->num_ref_edges() - num_edges;
  num_faces = gti->num_ref_faces() - num_faces;
  num_volumes = gti->num_ref_volumes() - num_volumes;
  num_bodies = gti->num_bodies() - num_bodies;

    // report results
  PRINT_SEPARATOR;
  if (!num_vertices && !num_edges && !num_faces && !num_volumes && !num_bodies)
    PRINT_INFO("Imprinting resulted in no new entities.\n");
  else {
    PRINT_INFO("Imprinting resulted in the following numbers of new entities:\n");
    if (num_vertices) PRINT_INFO("   %d vertices.\n", num_vertices);
    if (num_edges)    PRINT_INFO("   %d edges.\n", num_edges);
    if (num_faces)    PRINT_INFO("   %d faces.\n", num_faces);
    if (num_volumes)  PRINT_INFO("   %d volumes.\n", num_volumes);
    if (num_bodies)   PRINT_INFO("   %d bodies.\n", num_bodies);
  }
  PRINT_SEPARATOR;
  
    // ok, we're done
  return CUBIT_SUCCESS;
}

CubitStatus print_unmerged_surfaces()
{
    // Print number and ids of non-shared surfaces
    // Non-shared surfaces are one with < 2 volumes connected to them,
    // or less than two parent entities
  RefFace *face = RefEntityFactory::instance()->get_first_ref_face();

  DLIList<RefEntity*> parents;
  DLIList<RefFace*> unmerged_faces;
  int i;

  for (i = 0; i < GeometryQueryTool::instance()->num_ref_faces(); i++) {
    
      // get the parent volumes
    parents.clean_out();
    face->get_parent_ref_entities(parents);
    if (parents.size() < 2) unmerged_faces.append(face);
    
    face = GeometryQueryTool::instance()->get_next_ref_face();
  }
  
    // we should be at the beginning of the face list again
  assert(face == RefEntityFactory::instance()->get_first_ref_face());

    // now print information on unmerged surfaces
  
    // first cast faces to a cubit entity list; use cast_list_to_parent,
    // since it's much more efficient
  DLIList<CubitEntity*> temp_entities;
  CAST_LIST_TO_PARENT(unmerged_faces, temp_entities);
  PRINT_SEPARATOR;
  PRINT_INFO("There were %d unmerged surfaces; their ids are:\n",
             unmerged_faces.size());
  CubitUtil::list_entity_ids("\0", temp_entities);
  PRINT_SEPARATOR;
  
    // now we're done
  return CUBIT_SUCCESS;
}

  
