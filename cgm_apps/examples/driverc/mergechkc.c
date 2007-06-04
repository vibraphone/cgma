/**
 * \file mergechk.c
 *
 * \brief mergechk, a C implementation of the like-named C++ driver for CGM
 *
 * This program acts as a simple driver for CGM.  It reads in a geometry,
 * performs imprints between all the bodies, merges them, and writes information
 * on the results.  It also performs pairwise intersections between the
 * bodies to check for overlaps.  Results are written to stardard output.
 *
 */

#include "GeometryTool_ccapi.h"
#include "MergeTool_ccapi.h"
#include "CubitUtil_ccapi.h"
#include "CubitDefines.h"
#include "RefEntity_ccapi.h"
#include "Body_ccapi.h"
#include "RefVolume_ccapi.h"
#include "RefFace_ccapi.h"
#include "RefEdge_ccapi.h"
#include "RefVertex_ccapi.h"
#include "CubitObserver_ccapi.h"

/*  forward declare some functions used and defined later */
enum CubitStatus read_geometry(int, char **);
enum CubitStatus evaluate_overlaps();
enum CubitStatus imprint_bodies();
enum CubitStatus print_unmerged_surfaces();

/*  macro for printing a separator line */
#define PRINT_SEPARATOR  printf("=======================================\n");

/* / main program - initialize, then send to proper function */
int c_main (int argc, char **argv)
{
  void *dum = NULL;
  enum CubitStatus status;
  
  CubitObserver_init_static_observers(dum);

    /*  Initialize the GeometryTool */
  GeometryTool_instance(dum);

    /*  If there aren't any file arguments, print usage and exit */
  if (argc == 1) {
    printf("Usage: mergechk <geom_file> [<geom_file> ...]\n");
    exit(0);
  }
  
    /*  Read in the geometry from files specified on the command line */
  status = read_geometry(argc, argv);
  if (status == CUBIT_FAILURE) exit(1);
  else if (GeometryTool_num_bodies() == 0) {
    printf("No bodies read; exiting.\n");
    exit(0);
  }

    /*  Check for overlaps */
  status = evaluate_overlaps();
  if (status == CUBIT_FAILURE) exit(1);

    /*  Imprint bodies together, reporting on results */
  status = imprint_bodies();
  if (status == CUBIT_FAILURE) exit(1);
  
    /*  Merge bodies */
  status = MergeTool_merge_all_bodies();
  if (status == CUBIT_FAILURE) exit(1);
  
    /*  Print number and ids of non-shared surfaces */
  status = print_unmerged_surfaces();
  if (status == CUBIT_FAILURE) exit(1);

}

/* / attribs module: list, modify attributes in a give model or models */
/* /  */
/* / Arguments: file name(s) of geometry files in which to look */
/* / */
enum CubitStatus read_geometry(int num_files, char **argv) 
{
  enum CubitStatus status = CUBIT_SUCCESS;
  int i;
  
    /*  For each file, open and read the geometry */
  FILE *file_ptr;

  PRINT_SEPARATOR;

  for (i = 1; i < num_files; i++) {
    file_ptr = fopen(argv[i], "r");
    if (file_ptr == NULL) printf("Could not open file %s\n", argv[i]);
    else {
      status = GeometryTool_import_solid_model_1(file_ptr, argv[i], "ACIS_SAT",
                                                 NULL, CUBIT_TRUE, CUBIT_TRUE,
                                                 CUBIT_TRUE, CUBIT_TRUE, CUBIT_TRUE,
                                                 CUBIT_TRUE);
      if (status != CUBIT_SUCCESS) {
        printf("Problems reading geometry file %s.\n", argv[i]);
      }
      fclose(file_ptr);
    }
  }
  PRINT_SEPARATOR;

  return CUBIT_SUCCESS;
}

enum CubitStatus evaluate_overlaps() 
{
    /*  evaluate overlaps by intersecting bodies pairwise */
  int i, j;
  
    /*  make a copy of the body list for use in this function */
  void **all_bodies = NULL;
  void **all_new_bodies = NULL;
  void **temp_all_bodies = NULL;
  int all_bodies_size, all_new_bodies_size;
  void *tool;
  void **new_bodies = NULL;
  int new_bodies_size;
  enum CubitStatus status;
  int int_bodies_size;
  void **child_entities = NULL;
  void **temp_children = NULL;
  int child_entities_size, temp_children_size;
  void *this_entity;
  void **temp_entities = NULL;
  int num_entities;

  GeometryTool_bodies(&all_bodies, &all_bodies_size);
  
    /*  step backward on this list, extracting the last body and using it as a tool */
    /*  for remaining bodies */

  all_new_bodies_size = 0;
  all_new_bodies = malloc(all_bodies_size * all_bodies_size * sizeof(void *));
  int_bodies_size = all_bodies_size-1;
  
  for (i = 0; i < all_bodies_size-1; i++) {
    tool = all_bodies[i];
    
      /*  intersect the tool with remaining bodies; make sure and keep old bodies, */
      /*  since we're not using copies; save new bodies for evaluation later */

    temp_all_bodies = &all_bodies[i+1];
    status = GeometryTool_intersect(tool, 
                                    &temp_all_bodies, &int_bodies_size,
                                    &new_bodies, &new_bodies_size,
                                    CUBIT_TRUE);

    for (j = 0; j < new_bodies_size; j++)
      all_new_bodies[all_new_bodies_size++] = new_bodies[j];

    int_bodies_size--;
  }
  
    /*  count number of geometric entities in new bodies; if there are no overlaps, */
    /*  this number should be zero */
    /*  first get all child entities of the new bodies */
  child_entities = malloc(1000 * all_new_bodies_size * sizeof(void *));
  child_entities_size = 0;
  
  for (i = 0; i < all_new_bodies_size; i++) {
    this_entity = Body_get_address(all_new_bodies[i], RefEntity_TYPE);
    RefEntity_get_all_child_ref_entities(this_entity, &temp_children, &temp_children_size);

    for (j = 0; j < temp_children_size; j++)
      child_entities[child_entities_size++] = temp_children[j];
  }
  
    /*  then filter the list, keeping only unique entities */
  free(temp_children);
  temp_children = malloc(1000 * all_new_bodies_size * sizeof(void *));
  temp_children_size = 0;
  
  DLList_merge_unique(&temp_children, &temp_children_size, &child_entities,
                      &child_entities_size, CUBIT_FALSE);
  
    /*  now report */
  PRINT_SEPARATOR;
  if (temp_children_size == 0)
    printf("No body overlaps.\n");

  else {
      /*  check for various types of entities; check by casting the list */
      /*  to that type */
    printf("Overlaps found: \n");
    
      /*  check vertices */
      /*  temp_entities is cleaned out in the CAST_LIST macro, so no need */
      /*  to do that here */
    num_entities = 0;
    for (i = 0; i < temp_children_size; i++) {
      this_entity = RefEntity_get_address(temp_children[i], RefVertex_TYPE);
      if (this_entity != NULL) num_entities++;
    }

    if (num_entities > 0) 
      printf("   Vertices: %d\n", num_entities);
    
    num_entities = 0;
    for (i = 0; i < temp_children_size; i++) {
      this_entity = RefEntity_get_address(temp_children[i], RefEdge_TYPE);
      if (this_entity != NULL) num_entities++;
    }

    if (num_entities > 0) 
      printf("   Edges: %d\n", num_entities);
    
    num_entities = 0;
    for (i = 0; i < temp_children_size; i++) {
      this_entity = RefEntity_get_address(temp_children[i], RefFace_TYPE);
      if (this_entity != NULL) num_entities++;
    }

    if (num_entities > 0) 
      printf("   Faces: %d\n", num_entities);
    
    num_entities = 0;
    for (i = 0; i < temp_children_size; i++) {
      this_entity = RefEntity_get_address(temp_children[i], RefVolume_TYPE);
      if (this_entity != NULL) num_entities++;
    }

    if (num_entities > 0) 
      printf("   Volumes: %d\n", num_entities);
    
      /*  now delete all the bodies produced by the intersections */
    num_entities = 0;
    for (i = 0; i < all_new_bodies_size; i++) {
      GeometryTool_delete_Body_2(all_new_bodies[i], CUBIT_TRUE);
    }
  }
  
  PRINT_SEPARATOR;
  
    /*  we're done */
  return CUBIT_SUCCESS;
}

enum CubitStatus imprint_bodies() 
{
    /*  imprint all the bodies together, and report number of new  */
    /*  entities formed */
    /*  */

  void **old_bodies = NULL;
  void **new_bodies = NULL;
  int old_bodies_size, new_bodies_size;
  
    /*  first, count old entities */
  int num_vertices = GeometryTool_num_ref_vertices();
  int num_edges = GeometryTool_num_ref_edges();
  int num_faces = GeometryTool_num_ref_faces();
  int num_volumes = GeometryTool_num_ref_volumes();
  int num_bodies = GeometryTool_num_bodies();
  
    /*  imprint the bodies together, discarding old bodies */
  GeometryTool_bodies(&old_bodies, &old_bodies_size);
  GeometryTool_imprint_1(&old_bodies, &old_bodies_size, 
                         &new_bodies, &new_bodies_size, CUBIT_FALSE);
  
    /*  now count new numbers of entities, subtracting old numbers */
  num_vertices = GeometryTool_num_ref_vertices() - num_vertices;
  num_edges = GeometryTool_num_ref_edges() - num_edges;
  num_faces = GeometryTool_num_ref_faces() - num_faces;
  num_volumes = GeometryTool_num_ref_volumes() - num_volumes;
  num_bodies = GeometryTool_num_bodies() - num_bodies;

    /*  report results */
  PRINT_SEPARATOR;
  if (!num_vertices && !num_edges && !num_faces && !num_volumes && !num_bodies)
    printf("Imprinting resulted in no new entities.\n");
  else {
    printf("Imprinting resulted in the following numbers of new entities:\n");
    if (num_vertices) printf("   %d vertices.\n", num_vertices);
    if (num_edges)    printf("   %d edges.\n", num_edges);
    if (num_faces)    printf("   %d faces.\n", num_faces);
    if (num_volumes)  printf("   %d volumes.\n", num_volumes);
    if (num_bodies)   printf("   %d bodies.\n", num_bodies);
  }
  PRINT_SEPARATOR;
  
    /*  ok, we're done */
  return CUBIT_SUCCESS;
}

enum CubitStatus print_unmerged_surfaces()
{
    /*  Print number and ids of non-shared surfaces */
    /*  Non-shared surfaces are one with < 2 volumes connected to them, */
    /*  or less than two parent entities */
  void *face;
  void **parents = NULL;
  int parents_size;
  void **unmerged_faces = NULL;
  int unmerged_faces_size;
  int unmerged_faces_dim;
  int i;
  void **temp_entities = NULL;
  void *face_entity;

  unmerged_faces_dim = GeometryTool_num_ref_faces();
  unmerged_faces = malloc(unmerged_faces_dim * sizeof(void *));
  unmerged_faces_size = 0;
  
  face = GeometryTool_get_first_ref_face();
  face_entity = RefFace_get_address(face, RefEntity_TYPE);

  for (i = 0; i < GeometryTool_num_ref_faces(); i++) {
    
      /*  get the parent volumes */
    RefEntity_get_parent_ref_entities(face_entity, &parents, &parents_size);
    if (parents_size < 2) unmerged_faces[unmerged_faces_size++] = face;
    
    face = GeometryTool_get_next_ref_face();
  face_entity = RefFace_get_address(face, RefEntity_TYPE);
  }
  
    /*  we should be at the beginning of the face list again */
  if (face != GeometryTool_get_first_ref_face())
    printf("Warning: first face not reached after loop over all faces!!!\n");

    /*  now print information on unmerged surfaces */
  
    /*  first cast faces to a cubit entity list */
  temp_entities = malloc(unmerged_faces_size*sizeof(void *));
  
  for (i = 0; i < unmerged_faces_size; i++)
    temp_entities[i] = RefFace_get_address(unmerged_faces[i], CubitEntity_TYPE);
  
  PRINT_SEPARATOR;
  printf("There were %d unmerged surfaces; their ids are:\n",
         unmerged_faces_size);
  CubitUtil_list_entity_ids_1("\0", &temp_entities, &unmerged_faces_size, 80, "\n", 0, 1,
                              8, ",", "(none)");
  PRINT_SEPARATOR;
  
    /*  now we're done */
  return CUBIT_SUCCESS;
}

  
