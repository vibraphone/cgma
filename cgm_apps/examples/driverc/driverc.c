/* test program */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "GeometryTool_ccapi.h"
#include "CubitApp_ccapi.h"
#include "TopologyEntity_ccapi.h"
#include "CubitVectorStruct.h"
#include "Body_ccapi.h"
#include "RefFace_ccapi.h"

enum CubitStatus read_geom_file(int argc, char **argv);
enum CubitStatus evaluate_geometry();

int c_main(int argc, char **argv) 
{
  enum CubitStatus result;

  void *gti;
  void *cai;
  
  cai = CubitApp_instance();
  gti = GeometryTool_instance(NULL);

  if (argc == 1) {
    printf("Usage: %s <geom_file1.sat> [<geom_file2.sat> ...]\n", argv[0]);
    return 1;
  }
  
  printf("Reading geometry file(s)\n");
  result = read_geom_file(argc, argv);
  if (result == CUBIT_FAILURE) return 1;
  
    /* now do some fun geometry stuff */
  printf("Evaluating geometry\n");
  result = evaluate_geometry();
  if (result == CUBIT_FAILURE) return 1;

  return result;
}

enum CubitStatus read_geom_file(int argc, char **argv)
{
    /* For each file, open and read the geometry */
  FILE *file_ptr;
  enum CubitStatus status = CUBIT_FAILURE;

  int i;
  for (i = 1; i < argc; i++) {
    file_ptr = fopen(argv[i], "r");
    if (file_ptr == NULL) printf("Could not open file %s\n", argv[i]);
    else {
      status = GeometryTool_import_solid_model_1(file_ptr, argv[i], "ACIS_SAT",
                                                 NULL, CUBIT_FALSE,
                                                 CUBIT_TRUE, CUBIT_TRUE, CUBIT_TRUE,
                                                 CUBIT_TRUE, CUBIT_TRUE);
      if (status != CUBIT_SUCCESS) {
        printf("Problems reading geometry file %s.\n", argv[i]);
      }
      fclose(file_ptr);
    }
  }
  
  return status;
}

enum CubitStatus evaluate_geometry() 
{

    /* Evaluate the model; for each surface, on a 5x5 grid in parameter space, print:
     * - normal
     * - tangents in u and v directions
     * - curvatures in u and v
     */

  enum CubitStatus result;
  int body_num;
  void *body_ptr;
  void **ref_faces = NULL;
  int num_ref_faces;
  void *current_surf;
  double lower_u, upper_u, lower_v, upper_v;
  enum CubitBoolean status_u, status_v;
  const int NUM_STEPS = 5;
  const double PI = acos(-1.0);
  const double deg_per_radian = 360.0 / (2.0 * PI);
  int iu, iv;
  double range_u, range_v;
  int face_num;
  
    /* get all the surfaces in the model */
  GeometryTool_ref_faces(&ref_faces, &num_ref_faces);

  for (face_num = 0; face_num < num_ref_faces; face_num++) {

      /* get the face */
    current_surf = ref_faces[face_num];
    
      /* get the parameter space of the surface */
    status_u = RefFace_get_param_range_U(current_surf, &lower_u, &upper_u);
    status_v = RefFace_get_param_range_V(current_surf, &lower_v, &upper_v);

    if (status_u == CUBIT_FALSE || status_v == CUBIT_FALSE) {
      printf("Surface %d is not parametric in %s; skipping\n",
             face_num, (status_v ? "U" : "V"));
      continue;
    }
    else {
      printf("Surface %d is parametric in u and v; Bounds: %f < U < %f; %f < V < %f\n",
             face_num, lower_u, upper_u, lower_v, upper_v);
      printf("Surface data: u, v, normal(x,y,z), tangentu(x,y,z), tangentv(x,y,z), du, dv = \n");
    }
  
      /* compute the ranges of steps */
    range_u = (upper_u - lower_u) / (double) NUM_STEPS;
    range_v = (upper_v - lower_v) / (double) NUM_STEPS;

      /* now compute normal, tangents and curvatures at various points on the
         surface */
    for (iu = 0; iu <= NUM_STEPS; iu++) {
      for (iv = 0; iv <= NUM_STEPS; iv++) {
        double posu, posv;
        struct CubitVectorStruct new_pos, surf_normal;
        double curvu, curvv;
        struct CubitVectorStruct du, dv;
      
          /* compute the parameter positions */
        posu = lower_u + iu*range_u;
        posv = lower_v + iv*range_v;
      
          /* get the coordinates of the point on the surface */
        new_pos = RefFace_position_from_u_v(current_surf, posu, posv);
      
          /* get the normal; don't need a wrt_volume */
        surf_normal = RefFace_normal_at(current_surf, new_pos, NULL);
      
          /* get the tangent vectors, one for each parametric direction */
        result = RefFace_uv_derivitives(current_surf, posu, posv, &du, &dv);
      
          /* get the principle curvatures for the surface */
        result = RefFace_get_principal_curvatures(current_surf, new_pos, &curvu, &curvv, NULL);
      
          /* print the results */
        printf("%f %f %f %f %f %f %f %f %f %f %f %f %f \n", posu, posv,
               surf_normal.xVal, surf_normal.yVal, surf_normal.zVal,
               du.xVal, du.yVal, du.zVal,
               dv.xVal, dv.yVal, dv.zVal,
               curvu, curvv);
      }
    }
  
  }

    /* ok, we're done printing fun geometry stuff */
  return CUBIT_SUCCESS;
}
