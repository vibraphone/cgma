/**
 * \file facets.cpp
 *
 * \brief facets, a simple facet-based C++ driver for CGM
 *
 */

#ifdef NT
//#include "stdafx.h"
#endif

#include <math.h>
#include "AppUtil.hpp"
#include "CGMApp.hpp"
#include "GeometryQueryTool.hpp"
#include "FacetModifyEngine.hpp"
#include "CubitObserver.hpp"
#include "CubitPointData.hpp"
#include "CubitFacetData.hpp"
#include "DLIList.hpp"
#include "Surface.hpp"
#include "ShellSM.hpp"
#include "Lump.hpp"
#include "Body.hpp"
#include "RefFace.hpp"
#include "RefVertex.hpp"

extern "C" void gl_cleanup()
{}

/// main program - initialize, then send to proper function
int main (int argc, char **argv)
{

  if (argc == 1) {
    PRINT_INFO("Usage: facets <angle_val>\n"
               "  where <angle_val> is the facet dihedral angle below which a\n"
               "  model edge will be constructed.\n");
    return 0;
  }
  
    // first argument should be an angle
  double angle = atof(argv[1]);
    
    // Initialize the GeometryTool
  GeometryQueryTool *gqt = GeometryQueryTool::instance();
  FacetModifyEngine *fme = FacetModifyEngine::instance();

  DLIList<CubitFacet*> f_list;
  DLIList <CubitPoint*> p_list;

    // define some really simple facets, corresponding to the
    // facets for a brick
  typedef CubitPointData CPD; typedef CubitFacetData CFD;
  CPD *p1 = new CPD(0, 0, 0); CPD *p2 = new CPD(1, 0, 0);
  CPD *p3 = new CPD(1, 1, 0); CPD *p4 = new CPD(0, 1, 0);
  CPD *p5 = new CPD(0, 0, 1); CPD *p6 = new CPD(1, 0, 1);
  CPD *p7 = new CPD(1, 1, 1); CPD *p8 = new CPD(0, 1, 1);
  p_list.append(p1); p_list.append(p2);
  p_list.append(p3); p_list.append(p4);
  p_list.append(p5); p_list.append(p6);
  p_list.append(p7); p_list.append(p8);
  CFD *f1 = new CFD(p1,p3,p2); CFD *f2 = new CFD(p1,p4,p3); // xy-
  CFD *f3 = new CFD(p5,p6,p8); CFD *f4 = new CFD(p6,p7,p8); // xy+
  CFD *f5 = new CFD(p1,p8,p4); CFD *f6 = new CFD(p1,p5,p8); // yz-
  CFD *f7 = new CFD(p3,p7,p6); CFD *f8 = new CFD(p3,p6,p2); // yz+
  CFD *f9 = new CFD(p1,p6,p5); CFD *f10 = new CFD(p1,p2,p6); // xz-
  CFD *f11 = new CFD(p3,p4,p8); CFD *f12 = new CFD(p3,p8,p7); // xz+
  f_list.append(f1); f_list.append(f2); f_list.append(f3); 
  f_list.append(f4); f_list.append(f5); f_list.append(f6); 
  f_list.append(f7); f_list.append(f8); f_list.append(f9); 
  f_list.append(f10); f_list.append(f11); f_list.append(f12); 
  DLIList<LoopSM*> my_loops;

  DLIList<Surface*> surf_list;
  CubitStatus result;
  result = fme->build_facet_surface(NULL, f_list, p_list, angle, 4, 
                                    false, false, surf_list);
  
  if ( surf_list.size() == 0 || result != CUBIT_SUCCESS )
  {
    PRINT_ERROR("Problems building mesh based surfaces.\n");
    return result;
  }
  else
    PRINT_INFO("Constructed %d surfaces.\n", surf_list.size());

    //Now build the shell.  If we had it set up right this would be
    //in a loop.  We need to store list of DLBlockSurfaceLists on each
    //blockvolumemesh to store the shell information.  But that will
    //be saved for later.
  ShellSM *shell_ptr;
  result = fme->make_facet_shell(surf_list, shell_ptr);
                                 
  if ( shell_ptr == NULL || result != CUBIT_SUCCESS )
  {
    PRINT_ERROR("Problems building mesh based shell entity.\n");
    return result;
  }
  DLIList<ShellSM*> shell_list;
  shell_list.append(shell_ptr);
  Lump *lump_ptr;
  result = fme->make_facet_lump(shell_list, lump_ptr);
                                
  if ( lump_ptr == NULL || result != CUBIT_SUCCESS )
  {
    PRINT_ERROR("Problems building mesh based lump entity.\n");
    return result;
  }
  DLIList<Lump*> lump_list;
  lump_list.append(lump_ptr);

  BodySM *bodysm_ptr;
  Body *body_ptr;
  result = fme->make_facet_body(lump_list, bodysm_ptr);
                                
  body_ptr = GeometryQueryTool::instance()->make_Body(bodysm_ptr);
  
  if ( body_ptr == NULL || result != CUBIT_SUCCESS )
  {
    PRINT_ERROR("Problems building mesh based body entity.\n");
    return result;
  }

  if (!body_ptr) {
    exit(1);
  }

  PRINT_INFO("Body successfully created.\n");

  PRINT_INFO("Number of vertices = %d\n", gqt->num_ref_vertices());
  PRINT_INFO("Number of edges = %d\n", gqt->num_ref_edges());
  PRINT_INFO("Number of faces = %d\n", gqt->num_ref_faces());
  PRINT_INFO("Number of volumes = %d\n", gqt->num_ref_volumes());
  PRINT_INFO("Number of bodies = %d\n", gqt->num_bodies());

    // print vertex positions
  DLIList<RefVertex*> verts;
  gqt->ref_vertices(verts);
  int i;
  for (i = verts.size(); i > 0; i--) {
    CubitVector coords = verts.get_and_step()->coordinates();
    PRINT_INFO("Vertex %d: %4.2f, %4.2f, %4.2f.\n", 
               8-i, coords.x(), coords.y(), coords.z());
  }
  
  RefFace *face = gqt->get_first_ref_face();

    // find closest point on each coordinate plane face center
  CubitVector test_position, result_position, normal;
  
  test_position.set(.5, .5, 0);
  face->find_closest_point_trimmed(test_position, result_position);
  normal = face->normal_at(result_position);
  PRINT_INFO("Point (%4.2f, %4.2f, %4.2f):\n  "
             "closest=(%4.2f, %4.2f, %4.2f), normal=(%4.2f, %4.2f, %4.2f).\n",
             test_position.x(), test_position.y(), test_position.z(), 
             result_position.x(), result_position.y(), result_position.z(), 
             normal.x(), normal.y(), normal.z());
  
  test_position.set(.5, .5, 1);
  face->find_closest_point_trimmed(test_position, result_position);
  normal = face->normal_at(result_position);
  PRINT_INFO("Point (%4.2f, %4.2f, %4.2f):\n  "
             "closest=(%4.2f, %4.2f, %4.2f), normal=(%4.2f, %4.2f, %4.2f).\n",
             test_position.x(), test_position.y(), test_position.z(), 
             result_position.x(), result_position.y(), result_position.z(), 
             normal.x(), normal.y(), normal.z());
  
  test_position.set(.5, 0, .5);
  face->find_closest_point_trimmed(test_position, result_position);
  normal = face->normal_at(result_position);
  PRINT_INFO("Point (%4.2f, %4.2f, %4.2f):\n  "
             "closest=(%4.2f, %4.2f, %4.2f), normal=(%4.2f, %4.2f, %4.2f).\n",
             test_position.x(), test_position.y(), test_position.z(), 
             result_position.x(), result_position.y(), result_position.z(), 
             normal.x(), normal.y(), normal.z());
  
  test_position.set(.5, 1, .5);
  face->find_closest_point_trimmed(test_position, result_position);
  normal = face->normal_at(result_position);
  PRINT_INFO("Point (%4.2f, %4.2f, %4.2f):\n  "
             "closest=(%4.2f, %4.2f, %4.2f), normal=(%4.2f, %4.2f, %4.2f).\n",
             test_position.x(), test_position.y(), test_position.z(), 
             result_position.x(), result_position.y(), result_position.z(), 
             normal.x(), normal.y(), normal.z());

  test_position.set(0, .5, .5);
  face->find_closest_point_trimmed(test_position, result_position);
  normal = face->normal_at(result_position);
  PRINT_INFO("Point (%4.2f, %4.2f, %4.2f):\n  "
             "closest=(%4.2f, %4.2f, %4.2f), normal=(%4.2f, %4.2f, %4.2f).\n",
             test_position.x(), test_position.y(), test_position.z(), 
             result_position.x(), result_position.y(), result_position.z(), 
             normal.x(), normal.y(), normal.z());
  
  test_position.set(1, .5, .5);
  face->find_closest_point_trimmed(test_position, result_position);
  normal = face->normal_at(result_position);
  PRINT_INFO("Point (%4.2f, %4.2f, %4.2f):\n  "
             "closest=(%4.2f, %4.2f, %4.2f), normal=(%4.2f, %4.2f, %4.2f).\n",
             test_position.x(), test_position.y(), test_position.z(), 
             result_position.x(), result_position.y(), result_position.z(), 
             normal.x(), normal.y(), normal.z());

  double a = 1.0/3.0;
  double b = 1.0/3.0;
  double c = 1.0/3.0;
  CubitVector temp_point(a,b,c);
  CubitVector eval_point;

  f1->evaluate( temp_point, &eval_point );
  PRINT_INFO("Evaluation of facet 1 at (%4.2f, %4.2f, %4.2f) is "
             "(%4.2f, %4.2f, %4.2f).\n",
              a, b, c, eval_point.x(), eval_point.y(), eval_point.z() );

  a = 0.0; b = 1.0; c = 0.0;
  temp_point.set(a,b,c);
  f2->evaluate( temp_point, &eval_point );
  PRINT_INFO("Evaluation of facet 2 at (%4.2f, %4.2f, %4.2f) is (%4.2f, %4.2f, %4.2f).\n",
              a, b, c, eval_point.x(), eval_point.y(), eval_point.z() );
  a = 0.5; b = 0.5; c = 0.0;
  temp_point.set(a,b,c);
  f3->evaluate( temp_point, &eval_point );
  PRINT_INFO("Evaluation of facet 3 at (%4.2f, %4.2f, %4.2f) is (%4.2f, %4.2f, %4.2f).\n",
              a, b, c, eval_point.x(), eval_point.y(), eval_point.z() );
  a = 0.0; b = 0.5; c = 0.5;
  temp_point.set(a,b,c);
  f4->evaluate( temp_point, &eval_point );
  PRINT_INFO("Evaluation of facet 4 at (%4.2f, %4.2f, %4.2f) is (%4.2f, %4.2f, %4.2f).\n",
              a, b, c, eval_point.x(), eval_point.y(), eval_point.z() );

  int ret_val = ( CubitMessage::instance()->error_count() );

  return ret_val;
}
