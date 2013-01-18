//- Class:       Cholla
//- Description: C Interface for the Cholla module
//-              Facet-based geometry definition and evaluation
//- Owner:       Steven J. Owen
//- Checked by:
//- Version:

#include "time.h"
#include "Cholla.h"
#include "ChollaEngine.hpp"
#include "CubitPoint.hpp"
#include "CubitPointData.hpp"
#include "CubitFacetEdge.hpp"
#include "CubitFacetEdgeData.hpp"
#include "CubitFacet.hpp"
#include "CubitFacetData.hpp"
#include "CubitQuadFacet.hpp"
#include "CubitQuadFacetData.hpp"
#include "CastTo.hpp"
#include "TDFacetBoundaryPoint.hpp"
#include "FacetEvalTool.hpp"
#include "FacetDataUtil.hpp"

#include <iostream>
#include <fstream>
using std::ifstream;

#define END_OF_BLOCKS "END_OF_BLOCKS"

// Internal Static function prototypes

static void time_stamp( FILE *fp );
static CubitStatus build_facets(int numTri, int numQuad, int numEdge, 
                                int numVert, int* triEdge, 
                                int *quadEdge, 
                                int* edgeVert, double* vert,
                                CubitPoint **point_array, 
                                CubitFacetEdge **edge_array, 
                                CubitFacet **facet_array,
                                CubitQuadFacet **qfacet_array,
                                DLIList <FacetEntity *> &facet_list);
static CubitStatus extractEdgeTangsFromFile( int num_file_edges,
                                             int *edge_nodes_file,
                                             double *edge_node_tangs_file,
                                             int numTri,
                                             int *edgeVert,
                                             double *edgeVertTang );
static CubitStatus extractTriNormalsFromFile( int num_file_tris,
                                              int *tri_nodes_file,
                                              double *tri_node_norms_file,
                                              int numTri,
                                              int *edgeVert,
                                              int *triEdge,
                                              double *triVertNorms );
static CubitStatus extractQuadNormalsFromFile( int num_file_quads,
                                               int *quad_nodes_file,
                                               double *quad_node_norms_file,
                                               int numQuad,
                                               int *edgeVert,
                                               int *quadEdge,
                                               double *quadVertNorms );
static void constructTriVerts( int triEdge[3],
                               int *edgeVert,
                               int this_tri[3] );
static void constructQuadVerts( int quadEdge[3],
                                int *edgeVert,
                                int this_quad[4] );
static void checkMemoryAllocations( int num_nodes_per_elem,
                                    int additional_num_elems,
                                    int *num_elems,
                                    int **nodes,
                                    double **node_normals );

static int classify_local_convexity_at_edge(CubitFacetEdge *edge_ptr);

//===========================================================================
//  Function: constructBezier
//  Purpose:  Return control points for a quartic Bezier for each face, using
//            a given feature angle tolerance to distinguish C0 and C1 edges.
//            Separate arrays of tris and quad faces should be supplied
//  Date:     10/28/2002
//  Author:   sjowen
//===========================================================================
void constructBezier( double angle, int numTri, int numQuad, int numEdge,
                      int numVert, int* triEdge, int* quadEdge,
                      int* edgeVert, double* vert,
                      double* edgeCtrlPts, double* triCtrlPts,
                      double* quadCtrlPts )
{

  // create arrays of facet entities from the arrays

  CubitPoint **point_array = new CubitPoint * [numVert];
  CubitFacetEdge **edge_array = new CubitFacetEdge * [numEdge];
  CubitFacet **facet_array = NULL;
  if (numTri > 0)
    facet_array = new CubitFacet * [numTri];
  CubitQuadFacet **qfacet_array = NULL;
  if (numQuad > 0)
    qfacet_array = new CubitQuadFacet * [numQuad];
  DLIList<FacetEntity *> facet_list;
  CubitStatus stat = CUBIT_SUCCESS;
  stat = build_facets(numTri, numQuad, numEdge, numVert, triEdge, quadEdge, 
                      edgeVert, vert, point_array, edge_array, facet_array,
                      qfacet_array, facet_list);

  // create lists of points and edges

  int ii;
  DLIList<CubitPoint *> point_list;
  for (ii=0; ii<numVert; ii++)  
    point_list.append(point_array[ii]);  

  DLIList<CubitFacetEdge *> edge_list;
  for (ii=0; ii<numEdge; ii++)
    edge_list.append(edge_array[ii]);

  // prepare the ChollaEngine

  DLIList<FacetEntity *> fpoint_list, fedge_list;
  CAST_LIST_TO_PARENT( point_list, fpoint_list );
  CAST_LIST_TO_PARENT( edge_list,  fedge_list );
  ChollaEngine *ceng = new ChollaEngine( facet_list, fedge_list, fpoint_list );

  CubitBoolean use_feature_angle;
  if (angle < 0.00001 || angle > 179.9999)
    use_feature_angle = CUBIT_FALSE;
  else
    use_feature_angle = CUBIT_TRUE;
  int interp_order = 4;  // assume we are creating beziers
  CubitBoolean smooth_non_manifold = CUBIT_TRUE;
  CubitBoolean split_surfaces = CUBIT_FALSE;

  // create the geometry
    
  //CubitStatus rv = 
  ceng->create_geometry(use_feature_angle,angle,
                        interp_order,smooth_non_manifold,
                        split_surfaces);

  // get the control points on edges

  int jj;
  int index;
  CubitVector *control_points;
  CubitFacetEdge *edge_ptr;
  for (ii=0; ii<numEdge; ii++)
  {
    edge_ptr = edge_array[ii];
    control_points = edge_ptr->control_points();
    assert(control_points != NULL);
    for (jj=0; jj<NUM_EDGE_CPTS; jj++)
    {
      index = 3 * (ii * NUM_EDGE_CPTS + jj); 
      edgeCtrlPts[index] = control_points[jj].x();
      edgeCtrlPts[index+1] = control_points[jj].y();
      edgeCtrlPts[index+2] = control_points[jj].z();
    }
  }

  // get the control points on triangles

  CubitFacet *facet_ptr;
  for (ii=0; ii<numTri; ii++)
  {
    facet_ptr = facet_array[ii];
    control_points = facet_ptr->control_points();
    assert(control_points != NULL);
    for (jj=0; jj<NUM_TRI_CPTS; jj++)
    {
      index = 3 * (ii * NUM_TRI_CPTS + jj); 
      triCtrlPts[index] = control_points[jj].x();
      triCtrlPts[index+1] = control_points[jj].y();
      triCtrlPts[index+2] = control_points[jj].z();
    }
  }

    // get the control points on quads

  //CubitVector qctrl_points[NUM_QUAD_CPTS];
  CubitQuadFacet *qfacet_ptr;
  for (ii=0; ii<numQuad; ii++)
  {
    qfacet_ptr = qfacet_array[ii];
    index = 3 * (ii * NUM_QUAD_CPTS);
    qfacet_ptr->get_control_points( &(quadCtrlPts[index]) );
  }

  // delete the ChollaEngine

  ceng->delete_eval_tools();
  ceng->delete_me();
  delete ceng;

  // delete temp arrays (facets are deleted with the eval tool)

  for (ii=0; ii<numQuad; ii++)
  {
    qfacet_array[ii]->remove_tri_facets();
    delete qfacet_array[ii];
  }
  delete [] point_array;
  delete [] edge_array;
  if (facet_array)
    delete [] facet_array; 
  if (qfacet_array)
    delete [] qfacet_array; 

}

//===========================================================================
//  Function: constructTriNormals
//  Purpose:  Return normals and tangents for a set of triangles using
//            a given feature angle tolerance to distinguish C0 and C1 edges.
//  Date:     10/10/2003
//  Author:   sjowen
//===========================================================================
void constructTriNormals( double angle, int numTri, int numEdge,
                          int numVert, int* triEdge,
                          int* edgeVert, double* vert,
                          double* triVertNorms, double* edgeVertTang,
                          int* edgeTypeFlags)
{

  // create arrays of facet entities from the arrays

  int* quadEdge = NULL;
  CubitPoint **point_array = new CubitPoint * [numVert];
  CubitFacetEdge **edge_array = new CubitFacetEdge * [numEdge];
  CubitFacet **facet_array = NULL;
  if (numTri > 0)
    facet_array = new CubitFacet * [numTri];
  int numQuad = 0;
  CubitQuadFacet **qfacet_array = NULL;
  DLIList<FacetEntity *> facet_list;
  CubitStatus stat = CUBIT_SUCCESS;
  stat = build_facets(numTri, numQuad, numEdge, numVert, triEdge, quadEdge, 
                      edgeVert, vert, point_array, edge_array, facet_array,
                      qfacet_array, facet_list);

  // create lists of points and edges

  int ii;
  DLIList<CubitPoint *> point_list;
  for (ii=0; ii<numVert; ii++)  
    point_list.append(point_array[ii]);  

  DLIList<CubitFacetEdge *> edge_list;
  for (ii=0; ii<numEdge; ii++){
    edge_list.append(edge_array[ii]);
  }
  

  // prepare the ChollaEngine

  DLIList<FacetEntity *> fpoint_list, fedge_list;
  CAST_LIST_TO_PARENT( point_list, fpoint_list );
  CAST_LIST_TO_PARENT( edge_list,  fedge_list );
  ChollaEngine *ceng = new ChollaEngine( facet_list, fedge_list, fpoint_list );

  CubitBoolean use_feature_angle;
  if (angle < 0.00001 || angle > 179.9999)
    use_feature_angle = CUBIT_FALSE;
  else
    use_feature_angle = CUBIT_TRUE;
  int interp_order = 4;  // assume we are creating beziers
  CubitBoolean smooth_non_manifold = CUBIT_TRUE;
  CubitBoolean split_surfaces = CUBIT_FALSE;

  // create the geometry
    
  //CubitStatus rv = 
  ceng->create_geometry(use_feature_angle,angle,
                        interp_order,smooth_non_manifold,
                        split_surfaces);

  // get the control points on edges and infer the tangents

  int index = 0;
  CubitVector *control_points;
  CubitFacetEdge *edge_ptr;
  CubitVector v0, v1, t0, t1;
  for (ii=0; ii<numEdge; ii++)
  {
    edge_ptr = edge_array[ii];
    v0 = edge_ptr->point(0)->coordinates();
    v1 = edge_ptr->point(1)->coordinates();
    control_points = edge_ptr->control_points();
    assert(control_points != NULL);
    index = 6 * ii;
    t0 = control_points[0] - v0;
    t1 = v1 - control_points[2];
    if(edge_ptr->is_flipped()){
      t0 *= (-1.0);
      t1 *= (-1.0);
    }
    t0.normalize();
    t1.normalize();
    edgeVertTang[index  ] = t0.x();
    edgeVertTang[index+1] = t0.y();
    edgeVertTang[index+2] = t0.z();
    edgeVertTang[index+3] = t1.x();
    edgeVertTang[index+4] = t1.y();
    edgeVertTang[index+5] = t1.z();
    edgeTypeFlags[ii]=classify_local_convexity_at_edge(edge_ptr);
  }

  // get the normal on triangles

  CubitVector n0, n1, n2;
  CubitFacet *facet_ptr;
  for (ii=0; ii<numTri; ii++)
  {
    facet_ptr = facet_array[ii];
    index = 9 * ii;
    n0 = facet_ptr->point(0)->normal(facet_ptr);
    n1 = facet_ptr->point(1)->normal(facet_ptr);
    n2 = facet_ptr->point(2)->normal(facet_ptr);
    triVertNorms[index    ] = n0.x();
    triVertNorms[index + 1] = n0.y();
    triVertNorms[index + 2] = n0.z();
    triVertNorms[index + 3] = n1.x();
    triVertNorms[index + 4] = n1.y();
    triVertNorms[index + 5] = n1.z();
    triVertNorms[index + 6] = n2.x();
    triVertNorms[index + 7] = n2.y();
    triVertNorms[index + 8] = n2.z();
  }

  // delete the ChollaEngine

  ceng->delete_eval_tools();
  ceng->delete_me();
  delete ceng;

  // delete temp arrays (facets are deleted with the eval tool)

  for (ii=0; ii<numQuad; ii++)
  {
    qfacet_array[ii]->remove_tri_facets();
    delete qfacet_array[ii];
  }
  delete [] point_array;
  delete [] edge_array;
  if (facet_array)
    delete [] facet_array; 
  if (qfacet_array)
    delete [] qfacet_array; 
}

//===========================================================================
//  Function: constructQuadNormals
//  Purpose:  Return normals and tangents for a set of quads using
//            a given feature angle tolerance to distinguish C0 and C1 edges.
//  Date:     10/10/2003
//  Author:   sjowen
//===========================================================================
void constructQuadNormals( double angle, int numQuad, int numEdge,
                           int numVert, int* quadEdge,
                           int* edgeVert, double* vert,
                           double* quadVertNorms, double* edgeVertTang,
                           int* edgeTypeFlags )
{

  // create arrays of facet entities from the arrays

  int *triEdge = NULL;
  CubitPoint **point_array = new CubitPoint * [numVert];
  CubitFacetEdge **edge_array = new CubitFacetEdge * [numEdge];
  CubitFacet **facet_array = NULL;
  CubitQuadFacet **qfacet_array = NULL; 
  int numTri = 0;
  if (numQuad > 0)
    qfacet_array = new CubitQuadFacet * [numQuad];
  DLIList<FacetEntity *> facet_list;
  CubitStatus stat = CUBIT_SUCCESS;
  stat = build_facets(numTri, numQuad, numEdge, numVert, triEdge, quadEdge, 
                      edgeVert, vert, point_array, edge_array, facet_array,
                      qfacet_array, facet_list);
  // create lists of points and edges

  int ii;
  
  DLIList<CubitPoint *> point_list;
  for (ii=0; ii<numVert; ii++)  
    point_list.append(point_array[ii]);  

  DLIList<CubitFacetEdge *> edge_list;
  for (ii=0; ii<numEdge; ii++)
    edge_list.append(edge_array[ii]);

  // prepare the ChollaEngine

  DLIList<FacetEntity *> fpoint_list, fedge_list;
  CAST_LIST_TO_PARENT( point_list, fpoint_list );
  CAST_LIST_TO_PARENT( edge_list,  fedge_list );
  ChollaEngine *ceng = new ChollaEngine( facet_list, fedge_list, fpoint_list );

  CubitBoolean use_feature_angle;
  if (angle < 0.00001 || angle > 179.9999)
    use_feature_angle = CUBIT_FALSE;
  else
    use_feature_angle = CUBIT_TRUE;
  int interp_order = 4;  // assume we are creating beziers
  CubitBoolean smooth_non_manifold = CUBIT_TRUE;
  CubitBoolean split_surfaces = CUBIT_FALSE;

  // create the geometry
 
  ceng->create_geometry(use_feature_angle,angle,
                        interp_order,smooth_non_manifold,
                        split_surfaces);

  // get the control points on edges and infer the tangents

  int mydebug = 0;
  FILE *fp = NULL;
  if (mydebug)
  {
    fp = fopen("edges.txt", "w");
  }
  int index = 0;
  CubitVector *control_points;
  CubitFacetEdge *edge_ptr;
  CubitVector v0, v1, t0, t1;
  for (ii=0; ii<numEdge; ii++)
  {
    edge_ptr = edge_array[ii];
    v0 = edge_ptr->point(0)->coordinates();
    v1 = edge_ptr->point(1)->coordinates();
    control_points = edge_ptr->control_points();
    assert(control_points != NULL);
    if (mydebug)
    {
      int kk;
      for(kk=0; kk<3; kk++);
      fprintf(fp, "(%d) %.6lf %.6lf %.6lf\n", edge_ptr->id(), v0.x(), v0.y(), v0.z()); 
      fprintf(fp, "    %.6lf %.6lf %.6lf\n", control_points[0].x(), control_points[0].y(), control_points[0].z());
      fprintf(fp, "    %.6lf %.6lf %.6lf\n", control_points[1].x(), control_points[1].y(), control_points[1].z());
      fprintf(fp, "    %.6lf %.6lf %.6lf\n", control_points[2].x(), control_points[2].y(), control_points[2].z());
      fprintf(fp, "    %.6lf %.6lf %.6lf\n", v1.x(), v1.y(), v1.z()); 
    }
    index = 6 * ii; 
    t0 = control_points[0] - v0;
    t1 = v1 - control_points[2];
    if(edge_ptr->is_flipped()){
      t0 *= (-1.0);
      t1 *= (-1.0);
    }
    t0.normalize();
    t1.normalize();
    edgeVertTang[index  ] = t0.x();
    edgeVertTang[index+1] = t0.y();
    edgeVertTang[index+2] = t0.z();
    edgeVertTang[index+3] = t1.x();
    edgeVertTang[index+4] = t1.y();
    edgeVertTang[index+5] = t1.z();
    edgeTypeFlags[ii]=classify_local_convexity_at_edge(edge_ptr);
  }
  if (mydebug)
  {
    fclose(fp);
  }

  // get the normal on quad

  CubitVector n0, n1, n2, n3;
  CubitQuadFacet *qfacet_ptr;
  for (ii=0; ii<numQuad; ii++)
  {
    qfacet_ptr = qfacet_array[ii];
    index = 12 * ii;
    n0 = qfacet_ptr->point(0)->normal(qfacet_ptr);
    n1 = qfacet_ptr->point(1)->normal(qfacet_ptr);
    n2 = qfacet_ptr->point(2)->normal(qfacet_ptr);
    n3 = qfacet_ptr->point(3)->normal(qfacet_ptr);
    quadVertNorms[index     ] = n0.x();
    quadVertNorms[index + 1 ] = n0.y();
    quadVertNorms[index + 2 ] = n0.z();
    quadVertNorms[index + 3 ] = n1.x();
    quadVertNorms[index + 4 ] = n1.y();
    quadVertNorms[index + 5 ] = n1.z();
    quadVertNorms[index + 6 ] = n2.x();
    quadVertNorms[index + 7 ] = n2.y();
    quadVertNorms[index + 8 ] = n2.z();
    quadVertNorms[index + 9 ] = n3.x();
    quadVertNorms[index + 10] = n3.y();
    quadVertNorms[index + 11] = n3.z();
  }
  
  // delete the ChollaEngine

  ceng->delete_eval_tools();
  ceng->delete_me();
  delete ceng;

  // delete temp arrays

  for (ii=0; ii<numQuad; ii++)
  {
    qfacet_array[ii]->remove_tri_facets();
    delete qfacet_array[ii];
  }
  delete [] point_array;
  delete [] edge_array;
  if (facet_array)
    delete [] facet_array; 
  if (qfacet_array)
    delete [] qfacet_array; 
}

//===========================================================================
//  Function: build_facets
//  Purpose:  static function that creates facet entities
//            from initial the arrays
//  Date:     10/28/2002
//  Author:   sjowen
//===========================================================================
static CubitStatus build_facets(int numTri, int numQuad, int numEdge, 
                                int numVert, int* triEdge, int *quadEdge, 
                                int* edgeVert, double* vert,
                                CubitPoint **point_array, 
                                CubitFacetEdge **edge_array, 
                                CubitFacet **facet_array,
                                CubitQuadFacet **qfacet_array,
                                DLIList <FacetEntity *> &facet_list)
{
  // create the points

  double x, y, z;
  int ii;
  
  for (ii=0; ii<numVert; ii++)
  {
    x = vert[ii*3];
    y = vert[ii*3 + 1];
    z = vert[ii*3 + 2];
    point_array[ii] = (CubitPoint *) new CubitPointData( x, y, z );
    point_array[ii]->set_id(ii);
  }

  // create the edges

  int ip, jp;
  for (ii=0; ii<numEdge; ii++)
  {
    ip = edgeVert[ii*2];
    jp = edgeVert[ii*2 + 1];
    assert(ip < numVert && jp < numVert && ip >= 0 && jp >= 0);
    edge_array[ii] = (CubitFacetEdge *) new CubitFacetEdgeData( point_array[ip],
                                                                point_array[jp] );
    edge_array[ii]->set_id(ii);
  }

  // create the tri facets

  int begin_face;
  int jj, iedge;
  CubitFacetEdge *edges[4];
  CubitFacet *facet_ptr = NULL;

  for (ii=0; ii<numTri; ii++)
  {
    begin_face = 3 * ii;
    for(jj=0; jj<3; jj++)
    {
      iedge = triEdge[begin_face+jj];
      edges[jj] = edge_array[iedge];
    }
    facet_ptr = (CubitFacet *)
      new CubitFacetData(edges[0], edges[1], edges[2]);
    facet_list.append(facet_ptr);
    facet_array[ii] = facet_ptr;
  }

  // create the quad facets

  CubitQuadFacet *qfacet_ptr = NULL;
  for (ii=0; ii<numQuad; ii++)
  {
    begin_face = 4 * ii;
    for(jj=0; jj<4; jj++)
    {
      iedge = quadEdge[begin_face+jj];
      edges[jj] = edge_array[iedge];
    }
    qfacet_ptr = (CubitQuadFacet *)
      new CubitQuadFacetData(edges[0], edges[1], edges[2], edges[3]);
    facet_list.append(qfacet_ptr->get_tri_facet(0));
    facet_list.append(qfacet_ptr->get_tri_facet(1));
    qfacet_array[ii] = qfacet_ptr;
  }

  int mydebug = 0;
  if (mydebug)
  {
    DLIList<CubitFacet*> myfacet_list;
    CAST_LIST(facet_list, myfacet_list, CubitFacet);
    FacetDataUtil::write_facets("C:\\CUBIT\\cholla_apps\\examples\\sierra\\test.facets",myfacet_list);
  }

  return CUBIT_SUCCESS;

}

//===========================================================================
//  Function: evalBezierFace
//  Purpose:  Evaluate a set of quartic Bezier patches at a specified 
//            parametric location given its control points. Evaluates 
//            location, normal and/or derivative. Expects list of either 
//            quads or tris (not both)
//
//  Date:     10/28/2002
//  Author:   sjowen
//===========================================================================
void evalBezierFace( int numFace, int numEdge, int numVert, int numVertPerFace, 
                     int* faceEdge, int* edgeVert, 
                     double* vert,  double* faceCtrlPts,
                     double* edgeCtrlPts, 
                     int numLocs, double* paramLocation, 
                     double* location, double* normal, double* deriv )
{
    // create arrays of facet entities from the arrays

  int numTri = 0;
  int numQuad = 0;
  CubitPoint **point_array = new CubitPoint * [numVert];
  CubitFacetEdge **edge_array = new CubitFacetEdge * [numEdge];
  CubitFacet **facet_array = NULL;
  int *tri_edge = (numVertPerFace == 3) ? faceEdge : NULL;
  int *quad_edge = (numVertPerFace == 4) ? faceEdge : NULL;
  if (numVertPerFace == 3)
  {
    facet_array = new CubitFacet * [numFace];
    numTri = numFace;
  }
  CubitQuadFacet **qfacet_array = NULL;
  if (numVertPerFace == 4)
  {
    qfacet_array = new CubitQuadFacet * [numFace];
    numQuad = numFace;
  }
  DLIList<FacetEntity *> facet_list;
  CubitStatus stat = CUBIT_SUCCESS;
  stat = build_facets(numTri, numQuad, numEdge, numVert, tri_edge, quad_edge, 
                      edgeVert, vert, point_array, edge_array, facet_array,
                      qfacet_array, facet_list);

  // set the control points on edges

  int index, ii, jj;
  CubitFacetEdge *edge_ptr;
  for (ii=0; ii<numEdge; ii++)
  {
    edge_ptr = edge_array[ii];
    index = 3 * (ii * NUM_EDGE_CPTS);
    edge_ptr->set_control_points( &(edgeCtrlPts[index]) );
  }

  // set the control points on triangles

  CubitFacet *facet_ptr;
  for (ii=0; ii<numTri; ii++)
  {
    facet_ptr = facet_array[ii];
    index = 3 * (ii * NUM_TRI_CPTS); 
    facet_ptr->set_control_points( &(faceCtrlPts[index]) );
  }

    // set the control points on quads

  CubitQuadFacet *qfacet_ptr;
  for (ii=0; ii<numQuad; ii++)
  {
    qfacet_ptr = qfacet_array[ii];
    index = 3 * (ii * NUM_QUAD_CPTS);
    qfacet_ptr->set_control_points( &(faceCtrlPts[index]) );
  }

  // Do the evaluations

  for (ii=0; ii<numTri; ii++)
  {
    for (jj=0; jj<numLocs; jj++)
    {
      facet_ptr = facet_array[ii];
      CubitVector areacoord(paramLocation[jj*3],
                            paramLocation[jj*3+1],
                            paramLocation[jj*3+2]);
      CubitVector eval_point;
      CubitVector *eval_point_ptr = &eval_point;
      CubitVector eval_normal;
      CubitVector *eval_normal_ptr = NULL;
      if (normal != NULL)
      {
        eval_normal_ptr = &eval_normal;
      }
      facet_ptr->evaluate( areacoord, eval_point_ptr, eval_normal_ptr );
      index = 3 * ii * jj;
      location[index] = eval_point.x();
      location[index+1] = eval_point.y();
      location[index+2] = eval_point.z();
      if (normal != NULL)
      {
        normal[index] = eval_normal.x();
        normal[index+1] = eval_normal.y();
        normal[index+2] = eval_normal.z();
      }
    }
  }

  for (ii=0; ii<numQuad; ii++)
  {
    for (jj=0; jj<numLocs; jj++)
    {
      qfacet_ptr = qfacet_array[ii];
      CubitVector eval_point;
      CubitVector *eval_point_ptr = &eval_point;
      CubitVector eval_normal;
      CubitVector *eval_normal_ptr = NULL;
      if (normal != NULL)
        eval_normal_ptr = &eval_normal;
      CubitVector eval_du, eval_dv;
      CubitVector *eval_du_ptr = NULL;
      CubitVector *eval_dv_ptr = NULL;
      if (deriv != NULL)
      {
        eval_du_ptr = &eval_du;
        eval_dv_ptr = &eval_dv;
      }

      qfacet_ptr->evaluate( paramLocation[2*jj], paramLocation[2*jj+1],
                            eval_point_ptr, eval_normal_ptr,
                            eval_du_ptr, eval_dv_ptr);
      index = 3 * ii * jj;
      location[index] = eval_point.x();
      location[index+1] = eval_point.y();
      location[index+2] = eval_point.z();
      if (normal != NULL)
      {
        normal[index] = eval_normal.x();
        normal[index+1] = eval_normal.y();
        normal[index+2] = eval_normal.z();
      }
      if (deriv != NULL)
      {
        index = 6 * ii * jj;
        deriv[index] = eval_du.x();
        deriv[index+1] = eval_du.y();
        deriv[index+2] = eval_du.z();
        deriv[index+3] = eval_dv.x();
        deriv[index+4] = eval_dv.y();
        deriv[index+5] = eval_dv.z();
      }
    }
  }
  
  // delete the temp arrays

  for (ii=0; ii<numQuad; ii++)
    delete qfacet_array[ii];
  for (ii=0; ii<numTri; ii++)
    delete facet_array[ii];
  for (ii=0; ii<numEdge; ii++)
    delete edge_array[ii];
  for (ii=0; ii<numVert; ii++)
    delete point_array[ii];
  delete [] point_array;
  delete [] edge_array;
  if (facet_array)
    delete [] facet_array; 
  if (qfacet_array)
    delete [] qfacet_array; 
}

//===========================================================================
//  Function: evalBezierFaceFromNorms
//  Purpose:  This is the same as the evalBezierFace except it differes by 
//            the input arguments.  This function takes a list of normals 
//            and tangents at the face vertices and computes bezier control 
//            points internally.  Normals and tangents should have been 
//            computed previously in constructTriNormals or constructQuadNormals.  
//            This function is not as computationally efficient as evalBezierFace 
//            since it requires Bezier control points to be computed
//            as part of the call - however, it is more memory efficient, 
//            requiring fewer variables to be stored with the calling application. 
//            The following argument list describes only those arguments that 
//            differ from evalBezierFace above.
//
//  Date:     10/15/2003
//  Author:   sjowen
//===========================================================================
void evalBezierFaceFromNorms( int numFace, int numEdge, int numVert, 
                              int numVertPerFace, int* faceEdge, int* edgeVert, 
                              double* vert, double* vertNorms, double* edgeVertTang,
                              int numLocs, double* paramLocation, 
                              double* location, double* normal, double* deriv )
{
    // create arrays of facet entities from the arrays

  int numTri = 0;
  int numQuad = 0;
  CubitPoint **point_array = new CubitPoint * [numVert];
  CubitFacetEdge **edge_array = new CubitFacetEdge * [numEdge];
  CubitFacet **facet_array = NULL;
  int *triEdge = (numVertPerFace == 3) ? faceEdge : NULL;
  int *quadEdge = (numVertPerFace == 4) ? faceEdge : NULL;
  if (numVertPerFace == 3)
  {
    facet_array = new CubitFacet * [numFace];
    numTri = numFace;
  }
  CubitQuadFacet **qfacet_array = NULL;
  if (numVertPerFace == 4)
  {
    qfacet_array = new CubitQuadFacet * [numFace];
    numQuad = numFace;
  }
  DLIList<FacetEntity *> facet_list;
  CubitStatus stat = CUBIT_SUCCESS;
  stat = build_facets(numTri, numQuad, numEdge, numVert, triEdge, quadEdge, 
                      edgeVert, vert, point_array, edge_array, facet_array,
                      qfacet_array, facet_list);

  // set the normals on triangle vertices

  int index, ii, jj;
  CubitPoint *pt;
  CubitFacet *facet_ptr;
  CubitVector pt_normal;
  for (ii=0; ii<numTri; ii++)
  {
    facet_ptr = facet_array[ii];
    for(jj=0; jj<3; jj++)
    {
      index = ii * 9 + (3 * jj);
      pt_normal.x( vertNorms[ index ] );
      pt_normal.y( vertNorms[ index + 1 ] );
      pt_normal.z( vertNorms[ index + 2 ] );
      pt = facet_ptr->point(jj);
      TDFacetBoundaryPoint::add_facet_boundary_point( pt, facet_ptr, pt_normal );
    }
  }

  // set the normals on quad vertices

  CubitQuadFacet *qfacet_ptr;
  for (ii=0; ii<numQuad; ii++)
  {
    qfacet_ptr = qfacet_array[ii];
    for(jj=0; jj<4; jj++)
    {
      index = ii * 12 + (3 * jj);
      pt_normal.x( vertNorms[ index ] );
      pt_normal.y( vertNorms[ index + 1 ] );
      pt_normal.z( vertNorms[ index + 2 ] );
      pt = qfacet_ptr->point(jj);
      TDFacetBoundaryPoint::add_facet_boundary_point( pt, qfacet_ptr, pt_normal );
    }
  }

  // set the control points on edges

  int mydebug = 0;
  FILE *fp = NULL;
  if (mydebug)
    fp = fopen("edges.txt", "w");
  CubitFacetEdge *edge_ptr;
  CubitVector N0, N1, T0, T1, P0, P1;
  CubitVector Pi[3];
  CubitPoint *pt0, *pt1;
  for (ii=0; ii<numEdge; ii++)
  {
    index = 6 * ii;
    edge_ptr = edge_array[ii];
    pt0 = edge_ptr->point(0);
    pt1 = edge_ptr->point(1);
    P0 = pt0->coordinates();
    P1 = pt1->coordinates();
    N0 = pt0->normal( edge_ptr );
    N1 = pt1->normal( edge_ptr );
    T0.x( edgeVertTang[ index ] );
    T0.y( edgeVertTang[ index + 1] );
    T0.z( edgeVertTang[ index + 2] );
    T1.x( edgeVertTang[ index + 3] );
    T1.y( edgeVertTang[ index + 4] );
    T1.z( edgeVertTang[ index + 5] );
    FacetEvalTool::init_edge_control_points(P0, P1, N0, N1, T0, T1, Pi);
    if (mydebug)
    {
      int kk;
      for(kk=0; kk<3; kk++);
      fprintf(fp, "(%d) %.6lf %.6lf %.6lf\n", edge_ptr->id(), P0.x(), P0.y(), P0.z()); 
      fprintf(fp, "    %.6lf %.6lf %.6lf\n", Pi[0].x(), Pi[0].y(), Pi[0].z());
      fprintf(fp, "    %.6lf %.6lf %.6lf\n", Pi[1].x(), Pi[1].y(), Pi[1].z());
      fprintf(fp, "    %.6lf %.6lf %.6lf\n", Pi[2].x(), Pi[2].y(), Pi[2].z());
      fprintf(fp, "    %.6lf %.6lf %.6lf\n", P1.x(), P1.y(), P1.z()); 
    }
    edge_ptr->control_points( Pi, 4 );
  }
  if (mydebug)
    fclose(fp);

  // set the control points on triangles

  for (ii=0; ii<numTri; ii++)
  {
    facet_ptr = facet_array[ii];
    facet_ptr->init_patch();
  }

    // set the control points on quads

  for (ii=0; ii<numQuad; ii++)
  {
    qfacet_ptr = qfacet_array[ii];
    qfacet_ptr->init_patch();
  }

  // Do the evaluations

  index = 0;
  int nindex = 0;
  for (ii=0; ii<numTri; ii++)
  {
    facet_ptr = facet_array[ii];
    for (jj=0; jj<numLocs; jj++)
    {  
      CubitVector areacoord(paramLocation[jj*3],
                            paramLocation[jj*3+1],
                            paramLocation[jj*3+2]);
      CubitVector eval_point;
      CubitVector *eval_point_ptr = &eval_point;
      CubitVector eval_normal;
      CubitVector *eval_normal_ptr = NULL;
      if (normal != NULL)
      {
        eval_normal_ptr = &eval_normal;
      }
      facet_ptr->evaluate( areacoord, eval_point_ptr, eval_normal_ptr );
      location[index++] = eval_point.x();
      location[index++] = eval_point.y();
      location[index++] = eval_point.z();
      if (normal != NULL)
      {
        normal[nindex++] = eval_normal.x();
        normal[nindex++] = eval_normal.y();
        normal[nindex++] = eval_normal.z();
      }
    }
  }

  int dindex = 0;
  nindex = index = 0;
  for (ii=0; ii<numQuad; ii++)
  {
    for (jj=0; jj<numLocs; jj++)
    {
      qfacet_ptr = qfacet_array[ii];
      CubitVector eval_point;
      CubitVector *eval_point_ptr = &eval_point;
      CubitVector eval_normal;
      CubitVector *eval_normal_ptr = NULL;
      if (normal != NULL)
        eval_normal_ptr = &eval_normal;
      CubitVector eval_du, eval_dv;
      CubitVector *eval_du_ptr = NULL;
      CubitVector *eval_dv_ptr = NULL;
      if (deriv != NULL)
      {
        eval_du_ptr = &eval_du;
        eval_dv_ptr = &eval_dv;
      }

      qfacet_ptr->evaluate( paramLocation[2*jj], paramLocation[2*jj+1],
                            eval_point_ptr, eval_normal_ptr,
                            eval_du_ptr, eval_dv_ptr);
      location[index++] = eval_point.x();
      location[index++] = eval_point.y();
      location[index++] = eval_point.z();
      if (normal != NULL)
      {
        normal[nindex++] = eval_normal.x();
        normal[nindex++] = eval_normal.y();
        normal[nindex++] = eval_normal.z();
      }
      if (deriv != NULL)
      {
        index = 6 * ii * jj;
        deriv[dindex++] = eval_du.x();
        deriv[dindex++] = eval_du.y();
        deriv[dindex++] = eval_du.z();
        deriv[dindex++] = eval_dv.x();
        deriv[dindex++] = eval_dv.y();
        deriv[dindex++] = eval_dv.z();
      }
    }
  }
  
  // delete the temp arrays

  for (ii=0; ii<numQuad; ii++)
    delete qfacet_array[ii];
  for (ii=0; ii<numTri; ii++)
    delete facet_array[ii];
  for (ii=0; ii<numEdge; ii++)
    delete edge_array[ii];
  for (ii=0; ii<numVert; ii++)
    delete point_array[ii];
  delete [] point_array;
  delete [] edge_array;
  if (facet_array)
    delete [] facet_array; 
  if (qfacet_array)
    delete [] qfacet_array; 
}

//===========================================================================
//  Function: projToBezierFace
//  Purpose:  Project a set of x-y-z locations to Bezier patches.  Finds the 
//            closest point on one of the patches.  Evaluates location, 
//            normal and/or derivative.  Expects list of either quads or 
//            tris (not both)
//
//  Date:     12/08/2003
//  Author:   sjowen
//===========================================================================
void projToBezierFace( int numFace, int numEdge, int numVert, int numVertPerFace, 
                     int* faceEdge, int* edgeVert, 
                     double* vert,  double* faceCtrlPts,
                     double* edgeCtrlPts, 
                     int numLocs, double* xyz, 
                     int specify_tol, double converge_tol,
                     double* xyzOnFace, double* normal, double* deriv )
{
    // create arrays of facet entities from the arrays

  int numTri = 0;
  int numQuad = 0;
  CubitPoint **point_array = new CubitPoint * [numVert];
  CubitFacetEdge **edge_array = new CubitFacetEdge * [numEdge];
  CubitFacet **facet_array = NULL;
  int *tri_edge = (numVertPerFace == 3) ? faceEdge : NULL;
  int *quad_edge = (numVertPerFace == 4) ? faceEdge : NULL;
  if (numVertPerFace == 3)
  {
    facet_array = new CubitFacet * [numFace];
    numTri = numFace;
  }
  CubitQuadFacet **qfacet_array = NULL;
  if (numVertPerFace == 4)
  {
    qfacet_array = new CubitQuadFacet * [numFace];
    numQuad = numFace;
  }
  DLIList<FacetEntity *> facetentity_list;
  CubitStatus stat = CUBIT_SUCCESS;
  stat = build_facets(numTri, numQuad, numEdge, numVert, quad_edge, tri_edge, 
                      edgeVert, vert, point_array, edge_array, facet_array,
                      qfacet_array, facetentity_list);
  DLIList<CubitFacet *> facet_list;
  CAST_LIST(facetentity_list, facet_list, CubitFacet);

  // set the control points on edges

  double edgelen;
  double compare_tol = 0.0;
  int index, ii;
  CubitFacetEdge *edge_ptr;
  for (ii=0; ii<numEdge; ii++)
  {
    edge_ptr = edge_array[ii];
    index = 3 * (ii * NUM_EDGE_CPTS);
    edge_ptr->set_control_points( &(edgeCtrlPts[index]) );
    if ( !specify_tol )
    {
        edgelen = edge_ptr->length();
        compare_tol += edgelen; 
    }
  }
  if ( specify_tol )
  {
      compare_tol = converge_tol;
  }
  else
  {
      compare_tol /= numEdge;
      compare_tol *= 1.0e-3;
  }

  // set the control points on triangles

  CubitFacet *facet_ptr;
  for (ii=0; ii<numTri; ii++)
  {
    facet_ptr = facet_array[ii];
    index = 3 * (ii * NUM_TRI_CPTS); 
    facet_ptr->set_control_points( &(faceCtrlPts[index]) );
  }

  // set the control points on quads

  CubitQuadFacet *qfacet_ptr;
  for (ii=0; ii<numQuad; ii++)
  {
    qfacet_ptr = qfacet_array[ii];
    index = 3 * (ii * NUM_QUAD_CPTS);
    qfacet_ptr->set_control_points( &(faceCtrlPts[index]) );
  }

  // Do the projections

  CubitFacet *last_facet = NULL;
  int interp_order = 4;
  CubitBoolean trim = CUBIT_TRUE;
  CubitBoolean outside;
  for (ii=0; ii<numLocs; ii++)
  {
    CubitVector cur_location(xyz[ii*3],
                             xyz[ii*3+1],
                             xyz[ii*3+2]);
    CubitVector eval_point;
    CubitVector *eval_point_ptr = &eval_point;
    CubitVector eval_normal;
    CubitVector *eval_normal_ptr = NULL;
    if (normal != NULL)
    {
      eval_normal_ptr = &eval_normal;
    }
    FacetEvalTool::project_to_facets(facet_list, last_facet,
                                     interp_order, compare_tol,
                                     cur_location, trim,
                                     &outside, eval_point_ptr,
                                     eval_normal_ptr);
    index = 3 * ii;
    xyzOnFace[index] = eval_point.x();
    xyzOnFace[index+1] = eval_point.y();
    xyzOnFace[index+2] = eval_point.z();
    if (normal != NULL)
    {
      normal[index] = eval_normal.x();
      normal[index+1] = eval_normal.y();
      normal[index+2] = eval_normal.z();
    }
  }

  // no derivatives yet!!

  
  // delete the temp arrays

  for (ii=0; ii<numQuad; ii++)
    delete qfacet_array[ii];
  for (ii=0; ii<numTri; ii++)
    delete facet_array[ii];
  for (ii=0; ii<numEdge; ii++)
    delete edge_array[ii];
  for (ii=0; ii<numVert; ii++)
    delete point_array[ii];
  delete [] point_array;
  delete [] edge_array;
  if (facet_array)
    delete [] facet_array; 
  if (qfacet_array)
    delete [] qfacet_array; 
}

//===========================================================================
//  Function: projToBezierFaceFromNorms
//  Purpose:  This is the same as the projToBezierFace except it differes by 
//  the input arguments.  This function takes a list of normals and tangents 
//  at the face vertices and computes bezier control points internally.  
//  Normals and tangents should have been computed previously in 
//  constructTriNormals or constructQuadNormals.  This function is not as 
//  computationally efficient as evalBezierFace since it requires Bezier 
//  control points to be computed as part of the call - however, it is more 
//  memory efficient, requiring fewer variables to be stored with the calling 
//  application. The following argument list describes only those arguments 
//  that differ from projToBezierFace above.
//
//  Date:     12/08/2003
//  Author:   sjowen
//===========================================================================
void projToBezierFaceFromNorms( int numFace, int numEdge, int numVert, 
                              int numVertPerFace, int* faceEdge, int* edgeVert, 
                              double* vert, double* vertNorms, double* edgeVertTang,
                              int numLocs, double* xyz, 
                              int specify_tol, double converge_tol,
                              double* xyzOnFace, double* normal, double* deriv )
{
      // create arrays of facet entities from the arrays

  int numTri = 0;
  int numQuad = 0;
  CubitPoint **point_array = new CubitPoint * [numVert];
  CubitFacetEdge **edge_array = new CubitFacetEdge * [numEdge];
  CubitFacet **facet_array = NULL;
  int *triEdge = (numVertPerFace == 3) ? faceEdge : NULL;
  int *quadEdge = (numVertPerFace == 4) ? faceEdge : NULL;
  if (numVertPerFace == 3)
  {
    facet_array = new CubitFacet * [numFace];
    numTri = numFace;
  }
  CubitQuadFacet **qfacet_array = NULL;
  if (numVertPerFace == 4)
  {
    qfacet_array = new CubitQuadFacet * [numFace];
    numQuad = numFace;
  }
  DLIList<FacetEntity *> facetentity_list;
  CubitStatus stat = CUBIT_SUCCESS;
  stat = build_facets(numTri, numQuad, numEdge, numVert, triEdge, quadEdge, 
                      edgeVert, vert, point_array, edge_array, facet_array,
                      qfacet_array, facetentity_list);
  DLIList<CubitFacet *> facet_list;
  CAST_LIST(facetentity_list, facet_list, CubitFacet);

  // set the normals on triangle vertices

  int index, ii, jj;
  CubitPoint *pt;
  CubitFacet *facet_ptr;
  CubitVector pt_normal;
  for (ii=0; ii<numTri; ii++)
  {
    facet_ptr = facet_array[ii];
    for(jj=0; jj<3; jj++)
    {
      index = ii * 9 + (3 * jj);
      pt_normal.x( vertNorms[ index ] );
      pt_normal.y( vertNorms[ index + 1 ] );
      pt_normal.z( vertNorms[ index + 2 ] );
      pt = facet_ptr->point(jj);
      TDFacetBoundaryPoint::add_facet_boundary_point( pt, facet_ptr, pt_normal );
    }
  }

  // set the normals on quad vertices

  CubitQuadFacet *qfacet_ptr;
  for (ii=0; ii<numQuad; ii++)
  {
    qfacet_ptr = qfacet_array[ii];
    for(jj=0; jj<4; jj++)
    {
      index = ii * 12 + (3 * jj);
      pt_normal.x( vertNorms[ index ] );
      pt_normal.y( vertNorms[ index + 1 ] );
      pt_normal.z( vertNorms[ index + 2 ] );
      pt = qfacet_ptr->point(jj);
      TDFacetBoundaryPoint::add_facet_boundary_point( pt, qfacet_ptr, pt_normal );
    }
  }

  // set the control points on edges

  double edgelen;
  double compare_tol = 0.0;
  CubitFacetEdge *edge_ptr;
  CubitVector N0, N1, T0, T1, P0, P1;
  CubitVector Pi[3];
  CubitPoint *pt0, *pt1;
  for (ii=0; ii<numEdge; ii++)
  {
    index = 6 * ii;
    edge_ptr = edge_array[ii];
    pt0 = edge_ptr->point(0);
    pt1 = edge_ptr->point(1);
    P0 = pt0->coordinates();
    P1 = pt1->coordinates();
    N0 = pt0->normal( edge_ptr );
    N1 = pt1->normal( edge_ptr );
    T0.x( edgeVertTang[ index ] );
    T0.y( edgeVertTang[ index + 1] );
    T0.z( edgeVertTang[ index + 2] );
    T1.x( edgeVertTang[ index + 3] );
    T1.y( edgeVertTang[ index + 4] );
    T1.z( edgeVertTang[ index + 5] );
    FacetEvalTool::init_edge_control_points(P0, P1, N0, N1, T0, T1, Pi);
    edge_ptr->control_points( Pi, 4 );
    if ( !specify_tol )
    {
        edgelen = edge_ptr->length();
        compare_tol += edgelen; 
    }
  }
  if ( specify_tol )
  {
      compare_tol = converge_tol;
  }
  else
  {
      compare_tol /= numEdge;
      compare_tol *= 1.0e-3;
  }

  // set the control points on triangles

  for (ii=0; ii<numTri; ii++)
  {
    facet_ptr = facet_array[ii];
    facet_ptr->init_patch();
  }

    // set the control points on quads

  for (ii=0; ii<numQuad; ii++)
  {
    qfacet_ptr = qfacet_array[ii];
    qfacet_ptr->init_patch();
  }

    // Do the projections

  CubitFacet *last_facet = NULL;
  int interp_order = 4;
  CubitBoolean trim = CUBIT_TRUE;
  CubitBoolean outside;
  
  for (ii=0; ii<numLocs; ii++)
  {
    CubitVector cur_location(xyz[ii*3],
                             xyz[ii*3+1],
                             xyz[ii*3+2]);
    CubitVector eval_point;
    CubitVector *eval_point_ptr = &eval_point;
    CubitVector eval_normal;
    CubitVector *eval_normal_ptr = NULL;
    if (normal != NULL)
    {
      eval_normal_ptr = &eval_normal;
    }
    FacetEvalTool::project_to_facets(facet_list, last_facet,
                                     interp_order, compare_tol,
                                     cur_location, trim,
                                     &outside, eval_point_ptr,
                                     eval_normal_ptr);
      //double dist = cur_location.distance_between(eval_point);
    index = 3 * ii;
    xyzOnFace[index] = eval_point.x();
    xyzOnFace[index+1] = eval_point.y();
    xyzOnFace[index+2] = eval_point.z();
    if (normal != NULL)
    {
      normal[index] = eval_normal.x();
      normal[index+1] = eval_normal.y();
      normal[index+2] = eval_normal.z();
    }
  }

  // no derivatives yet!!

  // delete the temp arrays

  for (ii=0; ii<numQuad; ii++)
    delete qfacet_array[ii];
  for (ii=0; ii<numTri; ii++)
    delete facet_array[ii];
  for (ii=0; ii<numEdge; ii++)
    delete edge_array[ii];
  for (ii=0; ii<numVert; ii++)
    delete point_array[ii];
  delete [] point_array;
  delete [] edge_array;
  if (facet_array)
    delete [] facet_array; 
  if (qfacet_array)
    delete [] qfacet_array; 
}

//---------------------------------------------------------------------------
// Functions for debugging
//---------------------------------------------------------------------------

//===========================================================================
//  Function: dumpMesh
//  Purpose:  dump the face mesh to a file that can be read by CUBIT
//            Use the same definition of parameters as used with 
//            resolveFaceVectors
//  Date:     10/28/2002
//  Author:   sjowen
//===========================================================================
void dumpMesh(const char *fileName, int includeResults, double angle,
              int numTri, int numQuad, int numEdge,
              int numVert, int* triEdge, int* quadEdge,
              int* edgeVert, double* vert,
              double* edgeCtrlPts, double* triCtrlPts,
              double* quadCtrlPts)
{
  FILE *fp = fopen(fileName, "w");
  assert(fp != NULL); // couldn't open file for writing
  
  // write the header info

  fprintf(fp, "Cholla version 1.0\n");
  time_stamp( fp );
  fprintf(fp, "%d %d %d %d %d\n", numTri, numQuad, numEdge, numVert, includeResults);
  fprintf(fp, "%f\n", angle);
  
  // faceEdgeBegin

  int ii;

  // faceEdge

  fprintf(fp, "triEdge\n");
  for(ii=0; ii<numTri; ii++)
  {
    fprintf(fp, "%d %d %d\n", triEdge[ii*3], triEdge[ii*3+1], triEdge[ii*3+2]);
  }

  fprintf(fp, "quadEdge\n");
  for(ii=0; ii<numQuad; ii++)
  {
    fprintf(fp, "%d %d %d %d\n", quadEdge[ii*4], quadEdge[ii*4+1],
                                 quadEdge[ii*4+2], quadEdge[ii*4+3]);
  }

  // edgeVert

  fprintf(fp, "edgeVert\n");
  for (ii=0; ii<numEdge; ii++)
  {
    fprintf(fp, "%d %d\n", edgeVert[ii*2], edgeVert[ii*2+1]);
  }

  // vert

  fprintf(fp, "vert\n");
  for (ii=0; ii<numVert; ii++)
  {
    fprintf(fp, "%f %f %f\n", vert[ii*3], vert[ii*3+1], vert[ii*3+2]);
  }

  if (includeResults != 0)
  {
     // edge control points

    int jj;
    fprintf(fp, "edgeCtrlPts\n");
    for(ii=0; ii<numEdge; ii++)
    {
      for (jj=0; jj<NUM_EDGE_CPTS; jj++)
      {
        fprintf(fp, "%f %f %f\n", edgeCtrlPts[ii*jj*3], 
                                     edgeCtrlPts[ii*jj*3+1], 
                                     edgeCtrlPts[ii*jj*3+2]);
      }
    }

    // triangle control points

    fprintf(fp, "triCtrlPts\n");
    for(ii=0; ii<numTri; ii++)
    {
      for (jj=0; jj<NUM_TRI_CPTS; jj++)
      {
        fprintf(fp, "%f %f %f\n", triCtrlPts[ii*jj*3], 
                                     triCtrlPts[ii*jj*3+1], 
                                     triCtrlPts[ii*jj*3+2]);
      }
    }

    // quad control points

    fprintf(fp, "quadCtrlPts\n");
    for(ii=0; ii<numQuad; ii++)
    {
      for (jj=0; jj<NUM_QUAD_CPTS; jj++)
      {
        fprintf(fp, "%f %f %f\n", quadCtrlPts[ii*jj*3], 
                                     quadCtrlPts[ii*jj*3+1], 
                                     quadCtrlPts[ii*jj*3+2]);
      }
    }
  }
  fclose(fp);
}

//===========================================================================
//  Function: dumpResults
//  Purpose:  dump only the results to a file
//  Date:     10/28/2002
//  Author:   sjowen
//===========================================================================
void dumpResults(const char *fileName, int numEdge, int numTri, int numQuad, 
                 double* edgeCtrlPts, double* triCtrlPts, 
                 double* quadCtrlPts )
{
  FILE *fp = fopen( fileName, "w" );
  assert(fp != NULL);

    // write the header info

  fprintf(fp, "Cholla version 1.0 Results\n");
  time_stamp( fp );
  fprintf(fp, "%d %d %d\n", numEdge, numTri, numQuad );

  int ii,jj;
  fprintf(fp, "edgeCtrlPts\n");
  for(ii=0; ii<numEdge; ii++)
  {
    for (jj=0; jj<NUM_EDGE_CPTS; jj++)
    {
      fprintf(fp, "%f %f %f\n", edgeCtrlPts[ii*jj*3], 
                                   edgeCtrlPts[ii*jj*3+1], 
                                   edgeCtrlPts[ii*jj*3+2]);
    }
  }

  // triangle control points

  fprintf(fp, "triCtrlPts\n");
  for(ii=0; ii<numTri; ii++)
  {
    for (jj=0; jj<NUM_TRI_CPTS; jj++)
    {
      fprintf(fp, "%f %f %f\n", triCtrlPts[ii*jj*3], 
                                   triCtrlPts[ii*jj*3+1], 
                                   triCtrlPts[ii*jj*3+2]);
    }
  }

  // quad control points

  fprintf(fp, "quadCtrlPts\n");
  for(ii=0; ii<numQuad; ii++)
  {
    for (jj=0; jj<NUM_QUAD_CPTS; jj++)
    {
      fprintf(fp, "%f %f %f\n", quadCtrlPts[ii*jj*3], 
                                   quadCtrlPts[ii*jj*3+1], 
                                   quadCtrlPts[ii*jj*3+2]);
    }
  }

  fclose(fp);
}

//===========================================================================
//  Function: importMesh
//  Purpose:  import the face mesh from the specified file
//            allocate the arrays and pass them back
//  Date:     11/28/2002
//  Author:   sjowen
//===========================================================================
void importMesh(const char *fileName, int *includeResults,
                double *angle, int *numTri, int *numQuad, int *numEdge,
                int *numVert, int** triEdge, int** quadEdge,
                int** edgeVert, double** vert,
                double** edgeCtrlPts, double** triCtrlPts,
                double** quadCtrlPts)
{
  FILE *fp = fopen(fileName, "r");
  assert(fp != NULL);  // couldn't open file for reading
  
  // read the header info

  char version[128];
  char* s = fgets(version, 128, fp);
  assert(strcmp(version, "Cholla version 1.0\n") == 0);  // don't recognize file
  char filetime[128];
  s = fgets(filetime, 128, fp); 
  
  char line[256];
  int num_tri;
  int num_quad;
  int num_edge;
  int num_vert;
  s = fgets(line,256,fp);
  sscanf(line, "%d %d %d %d %d\n", &num_tri, &num_quad, &num_edge, &num_vert, includeResults); 
  *numTri = num_tri;
  *numQuad = num_quad;
  *numEdge = num_edge;
  *numVert = num_vert;
  int nread = fscanf(fp, "%lf\n", angle);

  // triEdge

  *triEdge = new int [num_tri*3];
  int *tri_edge = *triEdge;
  char array_name[128];
  nread = fscanf(fp, "%s\n", array_name);
  assert(strcmp(array_name, "triEdge") == 0);  // check start of array
  int ii;
  for(ii=0; ii<num_tri; ii++)
  {
    nread = fscanf(fp, "%d %d %d\n", &tri_edge[ii*3], &tri_edge[ii*3+1], &tri_edge[ii*3+2]);
  }

  // quadEdge

  *quadEdge = new int [num_quad*4];
  int *quad_edge = *quadEdge;
  nread = fscanf(fp, "%s\n", array_name);
  assert(strcmp(array_name, "quadEdge") == 0);  // check start of array
  for(ii=0; ii<num_quad; ii++)
  {
    nread = fscanf(fp, "%d %d %d %d\n", &quad_edge[ii*4], &quad_edge[ii*4+1], 
                                &quad_edge[ii*4+2], &quad_edge[ii*4+3]);
  }

  // edgeVert

  *edgeVert = new int [2*num_edge];
  int *edge_vert = *edgeVert;
  nread = fscanf(fp, "%s\n", array_name);
  assert(strcmp(array_name, "edgeVert") == 0);  // check start of array
  for (ii=0; ii<num_edge; ii++)
  {
    nread = fscanf(fp, "%d %d\n", &edge_vert[ii*2], &edge_vert[ii*2+1]);
  }

  // vert

  *vert = new double [3*num_vert];
  double *my_vert = *vert;
  nread = fscanf(fp, "%s\n", array_name);
  assert(strcmp(array_name, "vert") == 0);  // check start of array
  for (ii=0; ii<num_vert; ii++)
  {
    nread = fscanf(fp, "%lf %lf %lf\n", &my_vert[ii*3], &my_vert[ii*3+1], &my_vert[ii*3+2]);
  }

  *edgeCtrlPts = new double [3*num_edge*NUM_EDGE_CPTS];
  if (numTri != 0)
    *triCtrlPts = new double [3*num_tri*NUM_TRI_CPTS];
  if (numQuad != 0)
    *quadCtrlPts = new double [3*num_quad*NUM_QUAD_CPTS];

  double *edge_ctrl_pts = *edgeCtrlPts;
  double *tri_ctrl_pts = *triCtrlPts;
  double *quad_ctrl_pts = *quadCtrlPts;
  if (*includeResults != 0)
  {
     // edge control points
    
    nread = fscanf(fp, "%s\n", array_name);
    assert(strcmp(array_name, "edgeCtrlPts") == 0);  // check start of array
    for(ii=0; ii<num_edge*NUM_EDGE_CPTS; ii++)
    {
      nread = fscanf(fp, "%lf %lf %lf\n", &edge_ctrl_pts[ii*3], 
                                  &edge_ctrl_pts[ii*3+1], 
                                  &edge_ctrl_pts[ii*3+2]);
    }

     // triangle control points
    
    nread = fscanf(fp, "%s\n", array_name);
    assert(strcmp(array_name, "triCtrlPts") == 0);  // check start of array
    for(ii=0; ii<num_tri*NUM_TRI_CPTS; ii++)
    {
      nread = fscanf(fp, "%lf %lf %lf\n", &tri_ctrl_pts[ii*3], 
                                  &tri_ctrl_pts[ii*3+1], 
                                  &tri_ctrl_pts[ii*3+2]);
    }

    // quad control points
    
    nread = fscanf(fp, "%s\n", array_name);
    assert(strcmp(array_name, "quadCtrlPts") == 0);  // check start of array
    for(ii=0; ii<num_quad*NUM_QUAD_CPTS; ii++)
    {
      nread = fscanf(fp, "%lf %lf %lf\n", &quad_ctrl_pts[ii*3], 
                                  &quad_ctrl_pts[ii*3+1], 
                                  &quad_ctrl_pts[ii*3+2]);
    }

  }
  fclose(fp);
}

//===========================================================================
//  Function: importResults
//  Purpose:  import only the results to from a file
//  Note:     allocates the normal and tangent arrays
//  Date:     10/28/2002
//  Author:   sjowen
//===========================================================================
void importResults(const char *fileName, int *numEdge, int *numTri, 
                   int *numQuad, double** edgeCtrlPts, double** triCtrlPts, 
                   double** quadCtrlPts )
{
  FILE *fp = fopen(fileName, "r");
  assert(fp != NULL);  // couldn't open file for reading
  
  // read the header info

  char version[128];
  char* s = fgets(version, 128, fp);
  assert(strcmp(version, "Cholla version 1.0 Results\n") == 0);  // don't recognize file
  char filetime[128];
  s = fgets(filetime, 128, fp);
  int nread = fscanf(fp, "%d %d %d\n", numEdge, numTri, numQuad);
  int num_edge = *numEdge;
  int num_tri = *numTri;
  int num_quad = *numQuad;
  
  *edgeCtrlPts = new double [3*num_edge*NUM_EDGE_CPTS];
  if (numTri != 0)
    *triCtrlPts = new double [3*num_tri*NUM_TRI_CPTS];
  if (numQuad != 0)
    *quadCtrlPts = new double [3*num_quad*NUM_QUAD_CPTS];

  double *edge_ctrl_pts = *edgeCtrlPts;
  double *tri_ctrl_pts = *triCtrlPts;
  double *quad_ctrl_pts = *quadCtrlPts;

   // edge control points
  
  int ii;
  char array_name[128];
  nread = fscanf(fp, "%s\n", array_name);
  assert(strcmp(array_name, "edgeCtrlPts") == 0);  // check start of array
  for(ii=0; ii<num_edge*NUM_EDGE_CPTS; ii++)
  {
    nread = fscanf(fp, "%lf %lf %lf\n", &edge_ctrl_pts[ii*3], 
                                &edge_ctrl_pts[ii*3+1], 
                                &edge_ctrl_pts[ii*3+2]);
  }

   // triangle control points
  
  nread = fscanf(fp, "%s\n", array_name);
  assert(strcmp(array_name, "triCtrlPts") == 0);  // check start of array
  for(ii=0; ii<num_tri*NUM_TRI_CPTS; ii++)
  {
    nread = fscanf(fp, "%lf %lf %lf\n", &tri_ctrl_pts[ii*3], 
                                &tri_ctrl_pts[ii*3+1], 
                                &tri_ctrl_pts[ii*3+2]);
  }

  // quad control points
  
  nread = fscanf(fp, "%s\n", array_name);
  assert(strcmp(array_name, "quadCtrlPts") == 0);  // check start of array
  for(ii=0; ii<num_tri*NUM_QUAD_CPTS; ii++)
  {
    nread = fscanf(fp, "%lf %lf %lf\n", &quad_ctrl_pts[ii*3], 
                                &quad_ctrl_pts[ii*3+1], 
                                &quad_ctrl_pts[ii*3+2]);
  }

  fclose(fp);

}

//===========================================================================
//  Function: time_stamp
//  Purpose:  write the current time to a file
//  Date:     11/28/2002
//  Author:   sjowen
//===========================================================================
static void time_stamp( FILE *fp )
{
  struct tm *newtime;
  bool am = true;
  time_t long_time;

  time( &long_time );                /* Get time as long integer. */
  newtime = localtime( &long_time ); /* Convert to local time. */

  if( newtime->tm_hour > 12 )        /* Set up extension. */
    am = false;
  if( newtime->tm_hour > 12 )        /* Convert from 24-hour */
          newtime->tm_hour -= 12;    /*   to 12-hour clock.  */
  if( newtime->tm_hour == 0 )        /*Set hour to 12 if midnight. */
          newtime->tm_hour = 12;

  fprintf( fp, "%.19s %s\n", asctime( newtime ), am ? "AM" : "PM" );
}

//===========================================================================
//  Function: evalBezierEdge
//  Purpose:  Evaluate a quartic Bezier curve at a specified parametric
//  location given its control points. Evaluates location, and/or tangent.

//  Date:     01/20/2004
//  Author:   jfowler
//===========================================================================
void evalBezierEdge( int numEdge, int numVert, int* edgeVert,
                     double* vert, double* edgeCtrlPts, 
                     int numLocs, double* paramLocation, 
                     double* location, double* tangent )
{
CubitPoint **point_array = new CubitPoint * [numVert];
CubitFacetEdge **edge_array = new CubitFacetEdge * [numEdge];
CubitFacet **facet_array = NULL;
CubitQuadFacet **qfacet_array = NULL;
DLIList<FacetEntity *> facet_list;
int ii, index, numTri, numQuad;
int *triEdge, *quadEdge;
  triEdge = quadEdge = NULL;
  numTri = numQuad = 0;
  CubitStatus stat = CUBIT_SUCCESS;
  stat = build_facets(numTri, numQuad, numEdge, numVert, triEdge, quadEdge, 
                      edgeVert, vert, point_array, edge_array, facet_array,
                      qfacet_array, facet_list);

  CubitFacetEdge *edge_ptr;
  //  Set the control points on the edges.
  for (ii=0; ii<numEdge; ii++)
  {
    edge_ptr = edge_array[ii];
    index = 3 * (ii * NUM_EDGE_CPTS);
    edge_ptr->set_control_points( &(edgeCtrlPts[index]) );
  } 

  // do the evaluations
  int jj, tindex;
  index = tindex = 0;
    
    for ( ii = 0; ii < numEdge; ii++ ) {
    for ( jj = 0; jj < numLocs; jj++ ) {
      
      CubitVector eval_point;
      CubitVector *eval_point_ptr = &eval_point;
      CubitVector eval_tangent;
      CubitVector *eval_tangent_ptr = NULL;
      if ( tangent != NULL ) 
        eval_tangent_ptr = &eval_tangent;
      edge_ptr = edge_array[ii];
      edge_ptr->evaluate_single(paramLocation[jj],eval_point_ptr);      
      location[index++] = eval_point.x();
      location[index++] = eval_point.y();
      location[index++] = eval_point.z();
      if ( tangent != NULL ) {  
        edge_ptr->evaluate_single_tangent(paramLocation[jj],eval_tangent_ptr);      
        eval_tangent.normalize(); 
        tangent[tindex++] = eval_tangent.x();
        tangent[tindex++] = eval_tangent.y();
        tangent[tindex++] = eval_tangent.z();             
      }      
    }  
  }
  //  delete the temp arrays

  for (ii=0; ii<numEdge; ii++)
    delete edge_array[ii];
  for (ii=0; ii<numVert; ii++)
    delete point_array[ii];
  delete [] point_array;
  delete [] edge_array;

}

//===========================================================================
//  Function: evalBezierEdgeFromTans
//  Purpose:  Evaluate a quartic Bezier curve at a specified parametric
//  location given its tangents at the end-points. Evaluates location, 
//  and/or tangent.
//  Date:     01/20/2004
//  Author:   jfowler
//===========================================================================
void evalBezierEdgeFromTans( int numEdge, int numVert, int* edgeVert, 
                              double* vert, double* edgeVertTang,
                              int numLocs, double* paramLocation, 
                              double* location, double* tangent )
{
CubitPoint **point_array = new CubitPoint * [numVert];
CubitFacetEdge **edge_array = new CubitFacetEdge * [numEdge];
CubitFacet **facet_array = NULL;
CubitQuadFacet **qfacet_array = NULL;
DLIList<FacetEntity *> facet_list;
int ii, index, numTri, numQuad;
int *triEdge, *quadEdge;
  triEdge = quadEdge = NULL;
  numTri = numQuad = 0;
  CubitStatus stat = CUBIT_SUCCESS;
  stat = build_facets(numTri, numQuad, numEdge, numVert, triEdge, quadEdge, 
                      edgeVert, vert, point_array, edge_array, facet_array,
                      qfacet_array, facet_list);

  // set the control points on edges

  int mydebug = 0;
  FILE *fp = NULL;
  if (mydebug)
    fp = fopen("edges.txt", "w");
  CubitFacetEdge *edge_ptr;
  CubitVector T0, T1, P0, P1;
  CubitVector Pi[3];
  CubitPoint *pt0, *pt1;
  for (ii=0; ii<numEdge; ii++)
  {
    index = 6 * ii;
    edge_ptr = edge_array[ii];
    pt0 = edge_ptr->point(0);
    pt1 = edge_ptr->point(1);
    P0 = pt0->coordinates();
    P1 = pt1->coordinates();
    T0.x( edgeVertTang[ index ] );
    T0.y( edgeVertTang[ index + 1] );
    T0.z( edgeVertTang[ index + 2] );
    T1.x( edgeVertTang[ index + 3] );
    T1.y( edgeVertTang[ index + 4] );
    T1.z( edgeVertTang[ index + 5] );
    FacetEvalTool::init_edge_control_points_single(P0, P1, T0, T1, Pi);
    if (mydebug)
    {
      int kk;
      for(kk=0; kk<3; kk++);
      fprintf(fp, "(%d) %.6lf %.6lf %.6lf\n", edge_ptr->id(), P0.x(), P0.y(), P0.z()); 
      fprintf(fp, "    %.6lf %.6lf %.6lf\n", Pi[0].x(), Pi[0].y(), Pi[0].z());
      fprintf(fp, "    %.6lf %.6lf %.6lf\n", Pi[1].x(), Pi[1].y(), Pi[1].z());
      fprintf(fp, "    %.6lf %.6lf %.6lf\n", Pi[2].x(), Pi[2].y(), Pi[2].z());
      fprintf(fp, "    %.6lf %.6lf %.6lf\n", P1.x(), P1.y(), P1.z()); 
    }
    edge_ptr->control_points( Pi, 4 );
  } 
  if (mydebug)
    fclose(fp);

  // do the evaluations
  int jj, tindex;
  index = tindex = 0;
    
  for ( ii = 0; ii < numEdge; ii++ ) {
    for ( jj = 0; jj < numLocs; jj++ ) {
      
      CubitVector eval_point;
      CubitVector *eval_point_ptr = &eval_point;
      CubitVector eval_tangent;
      CubitVector *eval_tangent_ptr = NULL;
      if ( tangent != NULL ) 
        eval_tangent_ptr = &eval_tangent;
      edge_ptr = edge_array[ii];
      edge_ptr->evaluate_single(paramLocation[jj],eval_point_ptr);      
      location[index++] = eval_point.x();
      location[index++] = eval_point.y();
      location[index++] = eval_point.z();
      if ( tangent != NULL ) {  
        edge_ptr->evaluate_single_tangent(paramLocation[jj],eval_tangent_ptr);      
        eval_tangent.normalize(); 
        tangent[tindex++] = eval_tangent.x();
        tangent[tindex++] = eval_tangent.y();
        tangent[tindex++] = eval_tangent.z();             
      }      
    }  
  }
  //  delete the temp arrays

  for (ii=0; ii<numEdge; ii++)
    delete edge_array[ii];
  for (ii=0; ii<numVert; ii++)
    delete point_array[ii];
  delete [] point_array;
  delete [] edge_array;
}

//===========================================================================
//  Function: projBezierEdgeFromTans
//  Purpose:  Project a point to an edge given tangents at the endpoints.
//            The routine only finds one point.  Rarely there might be
//            two or more closest points.
//  Date:     01/20/2004             
//  Author:   jfowler
//===========================================================================
void projToBezierEdgeFromTans( int numEdge, int numVert, int* edgeVert, 
                              double* vert, double* edgeVertTang,
                              int numLocs, double* xyz, 
                              double* xyzOnEdge, double* tangent )
{
CubitPoint **point_array = new CubitPoint * [numVert];
CubitFacetEdge **edge_array = new CubitFacetEdge * [numEdge];
CubitFacet **facet_array = NULL;
CubitQuadFacet **qfacet_array = NULL;
DLIList<FacetEntity *> facet_list;
int ii, index, numTri, numQuad;
int *triEdge, *quadEdge;
  triEdge = quadEdge = NULL;
  numTri = numQuad = 0;
  CubitStatus stat = CUBIT_SUCCESS;
  stat = build_facets(numTri, numQuad, numEdge, numVert, triEdge, quadEdge, 
                      edgeVert, vert, point_array, edge_array, facet_array,
                      qfacet_array, facet_list);

  // set the control points on edges

  CubitFacetEdge *edge_ptr;
  CubitVector T0, T1, P0, P1;
  CubitVector Pi[3];
  CubitPoint *pt0, *pt1;
  for (ii=0; ii<numEdge; ii++)
  {
    index = 6 * ii;
    edge_ptr = edge_array[ii];
    pt0 = edge_ptr->point(0);
    pt1 = edge_ptr->point(1);
    P0 = pt0->coordinates();
    P1 = pt1->coordinates();
    T0.x( edgeVertTang[ index ] );
    T0.y( edgeVertTang[ index + 1] );
    T0.z( edgeVertTang[ index + 2] );
    T1.x( edgeVertTang[ index + 3] );
    T1.y( edgeVertTang[ index + 4] );
    T1.z( edgeVertTang[ index + 5] );
    FacetEvalTool::init_edge_control_points_single(P0, P1, T0, T1, Pi);

    edge_ptr->control_points( Pi, 4 );
  } 

  // do the projections
  int jj, kk, xyzindex, tindex;
  index = tindex = 0;
  double dsquared, dsquared_test, dderiv, dderiv2, t, t_min = 0.0;
  CubitVector xyz_pt;//, second_d;
  t = 0.0;
  dsquared_test =  CUBIT_DBL_MAX;  
  for ( ii = 0; ii < numEdge; ii++ ) {
    edge_ptr = edge_array[ii];
    for ( jj = 0; jj < numLocs; jj++ ) {
      xyzindex = 3*jj;
      xyz_pt.x(xyz[xyzindex]);
      xyz_pt.y(xyz[xyzindex+1]);
      xyz_pt.z(xyz[xyzindex+2]);    
      t = -1.0;
      dsquared_test =  CUBIT_DBL_MAX;  
      //  Get square of distance from xyz to curve at intervals of 0.2
      //  in the parameter t.  Save the closest distance found as a 
      //  starting point for Newton-Raphson iterations.  It was found 
      //  that an interval of 0.5 was too coarse and that the N-R
      //  sometimes converged on a false value.
      CubitVector eval_point;
      CubitVector *eval_point_ptr = &eval_point;
      CubitVector eval_tangent;
      CubitVector *eval_tangent_ptr = &eval_tangent;
      CubitVector second_d;
      CubitVector *second_d_ptr = &second_d;
      for ( kk = 0; kk < 11; kk++ ) {
        edge_ptr->evaluate_single(t,eval_point_ptr);      
        dsquared = ( (xyz_pt.x()-eval_point.x())*(xyz_pt.x()-eval_point.x()) +
                     (xyz_pt.y()-eval_point.y())*(xyz_pt.y()-eval_point.y()) +
                     (xyz_pt.z()-eval_point.z())*(xyz_pt.z()-eval_point.z()) );             
        if ( fabs(dsquared) < fabs(dsquared_test) ) {
          t_min = t;
          dsquared_test = dsquared;
        }      
        t += 0.2;
      }
      double dderiva, dderivb;
      if ( t_min == -1.00 ) {
      //  Check whether the slope changed signs between -1.0 and -0.8.
      //  If so, the min must have occurred in this interval -- set
      //  t_min to -0.9.
        t = -1.0;
        edge_ptr->evaluate_single(t,eval_point_ptr);      
        edge_ptr->evaluate_single_tangent(t,eval_tangent_ptr);      
        dderiva = 2.0*( (eval_point.x()-xyz_pt.x())*eval_tangent.x() +
                       (eval_point.y()-xyz_pt.y())*eval_tangent.y() +
                       (eval_point.z()-xyz_pt.z())*eval_tangent.z() );
        t = -0.8;
        edge_ptr->evaluate_single(t,eval_point_ptr);      
        edge_ptr->evaluate_single_tangent(t,eval_tangent_ptr);      
        dderivb = 2.0*( (eval_point.x()-xyz_pt.x())*eval_tangent.x() +
                       (eval_point.y()-xyz_pt.y())*eval_tangent.y() +
                       (eval_point.z()-xyz_pt.z())*eval_tangent.z() );
        if ( dderiva*dderivb < 0.0 ) t_min = -0.9;
   
      } else if ( t_min == 1.00 ) {
      //  Check the other end of the range, similarly.
        t = 1.0;
        edge_ptr->evaluate_single(t,eval_point_ptr);      
        edge_ptr->evaluate_single_tangent(t,eval_tangent_ptr);      
        dderiva = 2.0*( (eval_point.x()-xyz_pt.x())*eval_tangent.x() +
                        (eval_point.y()-xyz_pt.y())*eval_tangent.y() +
                        (eval_point.z()-xyz_pt.z())*eval_tangent.z() );
        t = 0.8;
        edge_ptr->evaluate_single(t,eval_point_ptr);      
        edge_ptr->evaluate_single_tangent(t,eval_tangent_ptr);      
        dderivb = 2.0*( (eval_point.x()-xyz_pt.x())*eval_tangent.x() +
                        (eval_point.y()-xyz_pt.y())*eval_tangent.y() +
                        (eval_point.z()-xyz_pt.z())*eval_tangent.z() );
        if ( dderiva*dderivb < 0.0 ) t_min = 0.9;  
      }      
      t = t_min;
      if ( (t > -1.0) && (t < 1.0) ) {
        int mm;
        mm = 0;
        while ( mm < 10 ) {  //  to avoid possible infinite loop
          mm++;
          edge_ptr->evaluate_single(t,eval_point_ptr);      
          edge_ptr->evaluate_single_tangent(t,eval_tangent_ptr);      
          dderiv = 2.0*( (eval_point.x()-xyz_pt.x())*eval_tangent.x() +
                         (eval_point.y()-xyz_pt.y())*eval_tangent.y() +
                         (eval_point.z()-xyz_pt.z())*eval_tangent.z() );

          edge_ptr->evaluate_2nd_derivative(t, second_d_ptr);            
  
          dderiv2 = (eval_point.x()-xyz_pt.x())*second_d.x() + 
                     eval_tangent.x()*eval_tangent.x() +
                    (eval_point.y()-xyz_pt.y())*second_d.y() + 
                     eval_tangent.y()*eval_tangent.y() +
                    (eval_point.z()-xyz_pt.z())*second_d.z() + 
                     eval_tangent.z()*eval_tangent.z();

           t = t - dderiv/dderiv2; // Newton-Raphson

           if ( t < -1.0 ) {
             t = -1.0;
             edge_ptr->evaluate_single(t,eval_point_ptr);      
             edge_ptr->evaluate_single_tangent(t,eval_tangent_ptr);      
             break;
           }
           if ( t > 1.0 ) {
             t = 1.0;
             edge_ptr->evaluate_single(t,eval_point_ptr);      
             edge_ptr->evaluate_single_tangent(t,eval_tangent_ptr);      
             break;
           }
           if ( fabs(dderiv) < 1.e-11 ) break;
        }
      } else {  //  At an endpoint of the paramaterization.
          edge_ptr->evaluate_single(t,eval_point_ptr);       
          edge_ptr->evaluate_single_tangent(t,eval_tangent_ptr);      
      }      
      xyzOnEdge[index++] = eval_point.x();
      xyzOnEdge[index++] = eval_point.y();
      xyzOnEdge[index++] = eval_point.z();
      if ( tangent != NULL ) {  
        eval_tangent.normalize(); 
        tangent[tindex++] = eval_tangent.x();
        tangent[tindex++] = eval_tangent.y();
        tangent[tindex++] = eval_tangent.z();             
      }      
    }  
  }  

  //  delete the temp arrays

  for (ii=0; ii<numEdge; ii++)
    delete edge_array[ii];
  for (ii=0; ii<numVert; ii++)
    delete point_array[ii];
  delete [] point_array;
  delete [] edge_array;

}
                       

//===========================================================================
//  Function: projBezierEdge
//  Purpose:  Project a point to an edge.
//            The routine only finds one point.  Rarely there might be
//            two or more closest points.
//  Date:     01/20/2004             
//  Author:   jfowler
//===========================================================================
void projToBezierEdge( int numEdge, int numVert, int* edgeVert, 
                              double* vert, double* edgeCtrlPts,
                              int numLocs, double* xyz, 
                              double* xyzOnEdge, double* tangent, double* t_value )
{
CubitPoint **point_array = new CubitPoint * [numVert];
CubitFacetEdge **edge_array = new CubitFacetEdge * [numEdge];
CubitFacet **facet_array = NULL;
CubitQuadFacet **qfacet_array = NULL;
DLIList<FacetEntity *> facet_list;
int ii, index, numTri, numQuad;
int *triEdge, *quadEdge;
  triEdge = quadEdge = NULL;
  numTri = numQuad = 0;
  CubitStatus stat = CUBIT_SUCCESS;
  stat = build_facets(numTri, numQuad, numEdge, numVert, triEdge, quadEdge, 
                      edgeVert, vert, point_array, edge_array, facet_array,
                      qfacet_array, facet_list);

  // set the control points on edges

  CubitFacetEdge *edge_ptr;
    //CubitVector T0, T1, P0, P1;
    //CubitVector Pi[3];

  for (ii=0; ii<numEdge; ii++)
  {
    edge_ptr = edge_array[ii];
    index = 3 * (ii * NUM_EDGE_CPTS);
    edge_ptr->set_control_points( &(edgeCtrlPts[index]) );
  } 

  // do the projections
  int jj, kk, xyzindex, tindex;
  index = tindex = 0;
  double dsquared, dsquared_test, dderiv, dderiv2, t, t_min = 0.0;
  CubitVector xyz_pt;//, second2142_d;
  t = 0.0;
  dsquared_test =  CUBIT_DBL_MAX;  
  for ( ii = 0; ii < numEdge; ii++ ) {
    edge_ptr = edge_array[ii];
    for ( jj = 0; jj < numLocs; jj++ ) {
      xyzindex = 3*jj;
      xyz_pt.x(xyz[xyzindex]);
      xyz_pt.y(xyz[xyzindex+1]);
      xyz_pt.z(xyz[xyzindex+2]);    
      t = -1.0;
      dsquared_test =  CUBIT_DBL_MAX;  
      //  Get square of distance from xyz to curve at intervals of 0.2
      //  in the parameter t.  Save the closest distance found as a 
      //  starting point for Newton-Raphson iterations.  It was found 
      //  that an interval of 0.5 was too coarse and that the N-R
      //  sometimes converged on a false value.
      CubitVector eval_point;
      CubitVector *eval_point_ptr = &eval_point;
      CubitVector eval_tangent;
      CubitVector *eval_tangent_ptr = &eval_tangent;
      CubitVector second_d;
      CubitVector *second_d_ptr = &second_d;
      for ( kk = 0; kk < 11; kk++ ) {
        edge_ptr->evaluate_single(t,eval_point_ptr);      
        dsquared = ( (xyz_pt.x()-eval_point.x())*(xyz_pt.x()-eval_point.x()) +
                     (xyz_pt.y()-eval_point.y())*(xyz_pt.y()-eval_point.y()) +
                     (xyz_pt.z()-eval_point.z())*(xyz_pt.z()-eval_point.z()) );             
        if ( fabs(dsquared) < fabs(dsquared_test) ) {
          t_min = t;
          dsquared_test = dsquared;
        }      
        t += 0.2;
      }
      double dderiva, dderivb;
      if ( t_min == -1.00 ) {
      //  Check whether the slope changed signs between -1.0 and -0.8.
      //  If so, the min must have occurred in this interval -- set
      //  t_min to -0.9.
        t = -1.0;
        edge_ptr->evaluate_single(t,eval_point_ptr);      
        edge_ptr->evaluate_single_tangent(t,eval_tangent_ptr);      
        dderiva = 2.0*( (eval_point.x()-xyz_pt.x())*eval_tangent.x() +
                       (eval_point.y()-xyz_pt.y())*eval_tangent.y() +
                       (eval_point.z()-xyz_pt.z())*eval_tangent.z() );
        t = -0.8;
        edge_ptr->evaluate_single(t,eval_point_ptr);      
        edge_ptr->evaluate_single_tangent(t,eval_tangent_ptr);      
        dderivb = 2.0*( (eval_point.x()-xyz_pt.x())*eval_tangent.x() +
                       (eval_point.y()-xyz_pt.y())*eval_tangent.y() +
                       (eval_point.z()-xyz_pt.z())*eval_tangent.z() );
        if ( dderiva*dderivb < 0.0 ) t_min = -0.9;
   
      } else if ( t_min == 1.00 ) {
      //  Check the other end of the range, similarly.
        t = 1.0;
        edge_ptr->evaluate_single(t,eval_point_ptr);      
        edge_ptr->evaluate_single_tangent(t,eval_tangent_ptr);      
        dderiva = 2.0*( (eval_point.x()-xyz_pt.x())*eval_tangent.x() +
                        (eval_point.y()-xyz_pt.y())*eval_tangent.y() +
                        (eval_point.z()-xyz_pt.z())*eval_tangent.z() );
        t = 0.8;
        edge_ptr->evaluate_single(t,eval_point_ptr);      
        edge_ptr->evaluate_single_tangent(t,eval_tangent_ptr);      
        dderivb = 2.0*( (eval_point.x()-xyz_pt.x())*eval_tangent.x() +
                        (eval_point.y()-xyz_pt.y())*eval_tangent.y() +
                        (eval_point.z()-xyz_pt.z())*eval_tangent.z() );
        if ( dderiva*dderivb < 0.0 ) t_min = 0.9;  
      }      
      t = t_min;
      if ( (t > -1.0) && (t < 1.0) ) {
        int mm;
        mm = 0;
        while ( mm < 10 ) {  //  to avoid possible infinite loop
          mm++;
          edge_ptr->evaluate_single(t,eval_point_ptr);      
          edge_ptr->evaluate_single_tangent(t,eval_tangent_ptr);      
          dderiv = 2.0*( (eval_point.x()-xyz_pt.x())*eval_tangent.x() +
                         (eval_point.y()-xyz_pt.y())*eval_tangent.y() +
                         (eval_point.z()-xyz_pt.z())*eval_tangent.z() );

          edge_ptr->evaluate_2nd_derivative(t, second_d_ptr);            
  
          dderiv2 = (eval_point.x()-xyz_pt.x())*second_d.x() + 
                     eval_tangent.x()*eval_tangent.x() +
                    (eval_point.y()-xyz_pt.y())*second_d.y() + 
                     eval_tangent.y()*eval_tangent.y() +
                    (eval_point.z()-xyz_pt.z())*second_d.z() + 
                     eval_tangent.z()*eval_tangent.z();

           t = t - dderiv/dderiv2; // Newton-Raphson

           if ( t < -1.0 ) {
             t = -1.0;
             edge_ptr->evaluate_single(t,eval_point_ptr);      
             edge_ptr->evaluate_single_tangent(t,eval_tangent_ptr);      
             break;
           }
           if ( t > 1.0 ) {
             t = 1.0;
             edge_ptr->evaluate_single(t,eval_point_ptr);      
             edge_ptr->evaluate_single_tangent(t,eval_tangent_ptr);      
             break;
           }
           if ( fabs(dderiv) < 1.e-11 ) break;
        }
      } else {  //  At an endpoint of the paramaterization.
          edge_ptr->evaluate_single(t,eval_point_ptr);       
          edge_ptr->evaluate_single_tangent(t,eval_tangent_ptr);      
      }      
      xyzOnEdge[index++] = eval_point.x();
      xyzOnEdge[index++] = eval_point.y();
      xyzOnEdge[index++] = eval_point.z();
      if ( tangent != NULL ) {  
        eval_tangent.normalize(); 
        tangent[tindex++] = eval_tangent.x();
        tangent[tindex++] = eval_tangent.y();
        tangent[tindex++] = eval_tangent.z();             
      }      
    }  
  } 
   
  *t_value = t;
  
  //  delete the temp arrays

  for (ii=0; ii<numEdge; ii++)
    delete edge_array[ii];
  for (ii=0; ii<numVert; ii++)
    delete point_array[ii];
  delete [] point_array;
  delete [] edge_array;

}
                       
//===========================================================================
//  Function: constructQuadNormalsFromFile
//  Purpose:  Similar to constructQuadNormals, except that instead of
//            approximating the surface normals and tangents from the input
//            facets, the normals and derivatives are extracted from a file
//            that was exported from Cubit along with the exodus file.  The
//            purpose for this is that Cubit will have the exact surface
//            normals and tangents, instead of the approximations from the
//            geometry.
//  Date:     08/31/2004
//  Author:   mlstate
//===========================================================================
CubitStatus constructQuadNormalsFromFile( const char *filename,
                                         int numQuad,
                                         int numEdge,
                                         int* quadEdge,
                                         int* edgeVert,
                                         double* quadVertNorms,
                                         double* edgeVertTang )
{
    CubitStatus status = CUBIT_SUCCESS;
    int num_file_quads = 0,
        num_file_edges = 0,
        *edge_nodes_file = NULL,
        *quad_nodes_file = NULL;
    double *edge_node_tangs_file = NULL,
           *quad_node_norms_file = NULL;

    status = readMBGNormalFile( filename,
                                NULL,
                                &num_file_quads,
                                &num_file_edges,
                                &edge_nodes_file,
                                &edge_node_tangs_file,
                                NULL,
                                NULL,
                                &quad_nodes_file,
                                &quad_node_norms_file );

    if ( status == CUBIT_SUCCESS && numEdge > 0 )
    {
        status = extractEdgeTangsFromFile( num_file_edges,
                                           edge_nodes_file,
                                           edge_node_tangs_file,
                                           numEdge,
                                           edgeVert,
                                           edgeVertTang );
    }
    if ( status == CUBIT_SUCCESS && numQuad > 0 )
    {
        status = extractQuadNormalsFromFile( num_file_quads,
                                             quad_nodes_file,
                                             quad_node_norms_file,
                                             numQuad,
                                             edgeVert,
                                             quadEdge,
                                             quadVertNorms );
    }

    delete [] edge_nodes_file;
    delete [] quad_nodes_file;
    delete [] edge_node_tangs_file;
    delete [] quad_node_norms_file;
    return status;
}
                       
//===========================================================================
//  Function: constructTriNormalsFromFile
//  Purpose:  Similar to constructTriNormals, except that instead of
//            approximating the surface normals and tangents from the input
//            facets, the normals and derivatives are extracted from a file
//            that was exported from Cubit along with the exodus file.  The
//            purpose for this is that Cubit will have the exact surface
//            normals and tangents, instead of the approximations from the
//            geometry.
//  Date:     08/31/2004
//  Author:   mlstate
//===========================================================================
CubitStatus constructTriNormalsFromFile( const char *filename,
                                         int numTri,
                                         int numEdge,
                                         int* triEdge,
                                         int* edgeVert,
                                         double* triVertNorms,
                                         double* edgeVertTang )
{
    CubitStatus status = CUBIT_SUCCESS;
    int num_file_tris = 0,
        num_file_edges = 0,
        *edge_nodes_file = NULL,
        *tri_nodes_file = NULL;
    double *edge_node_tangs_file = NULL,
           *tri_node_norms_file = NULL;

    status = readMBGNormalFile( filename,
                                &num_file_tris,
                                NULL,
                                &num_file_edges,
                                &edge_nodes_file,
                                &edge_node_tangs_file,
                                &tri_nodes_file,
                                &tri_node_norms_file,
                                NULL,
                                NULL );

    if ( status == CUBIT_SUCCESS && numEdge > 0 )
    {
        status = extractEdgeTangsFromFile( num_file_edges,
                                           edge_nodes_file,
                                           edge_node_tangs_file,
                                           numEdge,
                                           edgeVert,
                                           edgeVertTang );
    }
    if ( status == CUBIT_SUCCESS && numTri > 0 )
    {
        status = extractTriNormalsFromFile( num_file_tris,
                                            tri_nodes_file,
                                            tri_node_norms_file,
                                            numTri,
                                            edgeVert,
                                            triEdge,
                                            triVertNorms );
    }

    delete [] edge_nodes_file;
    delete [] tri_nodes_file;
    delete [] edge_node_tangs_file;
    delete [] tri_node_norms_file;
    return status;
}

//===========================================================================
//  Function: extractEdgeTangsFromFile
//  Purpose:  extract from the file the edge tangents.
//  Date:     08/31/2004
//  Author:   mlstate
//===========================================================================
static CubitStatus extractEdgeTangsFromFile
(
    int num_file_edges,
    int *edge_nodes_file,
    double *edge_node_tangs_file,
    int numEdge,
    int *edgeVert,
    double *edgeVertTang
)
{
#define MIN_EDGE_NODE_ID( n1, n2 ) n1 < n2 ? n1 : n2
    CubitStatus status = CUBIT_SUCCESS;
    DLIList<int *> *edge_hash = new DLIList<int*> [num_file_edges];
    int * edge_data = new int [num_file_edges * 3 ];

    int i;
    
    for ( i = 0; i < num_file_edges; i++ )
    {
        int *data = &(edge_data[i*3]);

        data[0] = edge_nodes_file[i*2];
        data[1] = edge_nodes_file[i*2 + 1];
        data[2] = i;

        int key = MIN_EDGE_NODE_ID( data[0], data[1] );
        key %= num_file_edges;
        edge_hash[key].append( data );
    }

    for ( i = 0; i < numEdge; i++ )
    {
        CubitBoolean found = CUBIT_FALSE;
        int *this_edge = &(edgeVert[i*2]);
        int key = MIN_EDGE_NODE_ID( this_edge[0], this_edge[1] );
        key %= num_file_edges;
        for ( int j = 0; j < edge_hash[key].size(); j++ )
        {
            int *data = edge_hash[key].get_and_step();
            if ( data[0] == this_edge[0] &&
                 data[1] == this_edge[1] )
            {
                edgeVertTang[ i*2*3     ] = edge_node_tangs_file[ data[2] * 2 * 3 ];
                edgeVertTang[ i*2*3 + 1 ] = edge_node_tangs_file[ data[2] * 2 * 3 + 1 ];
                edgeVertTang[ i*2*3 + 2 ] = edge_node_tangs_file[ data[2] * 2 * 3 + 2 ];
                edgeVertTang[ i*2*3 + 3 ] = edge_node_tangs_file[ data[2] * 2 * 3 + 3 ];
                edgeVertTang[ i*2*3 + 4 ] = edge_node_tangs_file[ data[2] * 2 * 3 + 4 ];
                edgeVertTang[ i*2*3 + 5 ] = edge_node_tangs_file[ data[2] * 2 * 3 + 5 ];
                found = CUBIT_TRUE;
                break;
            }
            else if ( data[0] == this_edge[1] &&
                      data[1] == this_edge[0] )
            {
                edgeVertTang[ i*2*3     ] = edge_node_tangs_file[ data[2] * 2 * 3 + 3 ];
                edgeVertTang[ i*2*3 + 1 ] = edge_node_tangs_file[ data[2] * 2 * 3 + 4 ];
                edgeVertTang[ i*2*3 + 2 ] = edge_node_tangs_file[ data[2] * 2 * 3 + 5 ];
                edgeVertTang[ i*2*3 + 3 ] = edge_node_tangs_file[ data[2] * 2 * 3  ];
                edgeVertTang[ i*2*3 + 4 ] = edge_node_tangs_file[ data[2] * 2 * 3 + 1 ];
                edgeVertTang[ i*2*3 + 5 ] = edge_node_tangs_file[ data[2] * 2 * 3 + 2 ];
                found = CUBIT_TRUE;
                break;
            }
        }
        if ( found == CUBIT_FALSE )
        {
            status = CUBIT_FAILURE;
            break;
        }
    }
    delete [] edge_data;
    delete [] edge_hash;

    return status;
 }

//===========================================================================
//  Function: extractTriNormalsFromFile
//  Purpose:  extract from the file the normals for a block of tri elements.
//  Date:     08/31/2004
//  Author:   mlstate
//===========================================================================
static CubitStatus extractTriNormalsFromFile
(
    int num_file_tris,
    int *tri_nodes_file,
    double *tri_node_norms_file,
    int numTri,
    int *edgeVert,
    int *triEdge,
    double *triVertNorms
)
{
#define MIN_TRI_NODE_ID( n1, n2, n3 ) n1 < n2 ?  n1 < n3 ? n1 : n3 : n2 < n3 ? n2 : n3

    CubitStatus status = CUBIT_SUCCESS;
    DLIList< int* > *tri_hash = new DLIList< int*>[num_file_tris];
    int * tri_data = new int [num_file_tris * 4 ];

    int i;
    
    for ( i = 0; i < num_file_tris; i++ )
    {
        int *data = &(tri_data[i*4]);

        data[0] = tri_nodes_file[i*3];
        data[1] = tri_nodes_file[i*3 + 1];
        data[2] = tri_nodes_file[i*3 + 2];
        data[3] = i;

        int key = MIN_TRI_NODE_ID( data[0], data[1], data[2] );
        key %= num_file_tris;
        tri_hash[key].append( data );
    }

    for ( i = 0; i < numTri; i++ )
    {
        CubitBoolean found = CUBIT_FALSE;
        int this_tri[3];
        
        constructTriVerts( &triEdge[i*3], edgeVert, this_tri );

        int key = MIN_TRI_NODE_ID( this_tri[0], this_tri[1], this_tri[2] );
        key %= num_file_tris;
        for ( int j = 0; j < tri_hash[key].size(); j++ )
        {
            int *data = tri_hash[key].get_and_step();

            if ( this_tri[0] == data[0] &&
                 this_tri[1] == data[1] &&
                 this_tri[2] == data[2] )
            {
                triVertNorms[ i*3*3     ] = tri_node_norms_file[ data[3]*3*3     ];
                triVertNorms[ i*3*3 + 1 ] = tri_node_norms_file[ data[3]*3*3 + 1 ];
                triVertNorms[ i*3*3 + 2 ] = tri_node_norms_file[ data[3]*3*3 + 2 ];
                triVertNorms[ i*3*3 + 3 ] = tri_node_norms_file[ data[3]*3*3 + 3 ];
                triVertNorms[ i*3*3 + 4 ] = tri_node_norms_file[ data[3]*3*3 + 4 ];
                triVertNorms[ i*3*3 + 5 ] = tri_node_norms_file[ data[3]*3*3 + 5 ];
                triVertNorms[ i*3*3 + 6 ] = tri_node_norms_file[ data[3]*3*3 + 6 ];
                triVertNorms[ i*3*3 + 7 ] = tri_node_norms_file[ data[3]*3*3 + 7 ];
                triVertNorms[ i*3*3 + 8 ] = tri_node_norms_file[ data[3]*3*3 + 8 ];
                found = CUBIT_TRUE;
                break;
            }
            else if ( this_tri[0] == data[1] &&
                      this_tri[1] == data[2] &&
                      this_tri[2] == data[0] )
            {
                triVertNorms[ i*3*3     ] = tri_node_norms_file[ data[3]*3*3 + 3 ];
                triVertNorms[ i*3*3 + 1 ] = tri_node_norms_file[ data[3]*3*3 + 4 ];
                triVertNorms[ i*3*3 + 2 ] = tri_node_norms_file[ data[3]*3*3 + 5 ];
                triVertNorms[ i*3*3 + 3 ] = tri_node_norms_file[ data[3]*3*3 + 6 ];
                triVertNorms[ i*3*3 + 4 ] = tri_node_norms_file[ data[3]*3*3 + 7 ];
                triVertNorms[ i*3*3 + 5 ] = tri_node_norms_file[ data[3]*3*3 + 8 ];
                triVertNorms[ i*3*3 + 6 ] = tri_node_norms_file[ data[3]*3*3     ];
                triVertNorms[ i*3*3 + 7 ] = tri_node_norms_file[ data[3]*3*3 + 1 ];
                triVertNorms[ i*3*3 + 8 ] = tri_node_norms_file[ data[3]*3*3 + 2 ];
                found = CUBIT_TRUE;
                break;
            }
            else if ( this_tri[0] == data[2] &&
                      this_tri[1] == data[0] &&
                      this_tri[2] == data[1] )
            {
                triVertNorms[ i*3*3     ] = tri_node_norms_file[ data[3]*3*3 + 6 ];
                triVertNorms[ i*3*3 + 1 ] = tri_node_norms_file[ data[3]*3*3 + 7 ];
                triVertNorms[ i*3*3 + 2 ] = tri_node_norms_file[ data[3]*3*3 + 8 ];
                triVertNorms[ i*3*3 + 3 ] = tri_node_norms_file[ data[3]*3*3     ];
                triVertNorms[ i*3*3 + 4 ] = tri_node_norms_file[ data[3]*3*3 + 1 ];
                triVertNorms[ i*3*3 + 5 ] = tri_node_norms_file[ data[3]*3*3 + 2 ];
                triVertNorms[ i*3*3 + 6 ] = tri_node_norms_file[ data[3]*3*3 + 3 ];
                triVertNorms[ i*3*3 + 7 ] = tri_node_norms_file[ data[3]*3*3 + 4 ];
                triVertNorms[ i*3*3 + 8 ] = tri_node_norms_file[ data[3]*3*3 + 5 ];
                found = CUBIT_TRUE;
                break;
            }
        }
        if ( found == CUBIT_FALSE )
        {
            status = CUBIT_FAILURE;
            break;
        }
    }

    delete [] tri_data;
    delete [] tri_hash;
    return status;
}

//===========================================================================
//  Function: extractQuadNormalsFromFile
//  Purpose:  extract from the file the normals for a block of quad elements.
//  Date:     08/31/2004
//  Author:   mlstate
//===========================================================================
static CubitStatus extractQuadNormalsFromFile
(
    int num_file_quads,
    int *quad_nodes_file,
    double *quad_node_norms_file,
    int numQuad,
    int *edgeVert,
    int *quadEdge,
    double *quadVertNorms
)
{
#define MIN_QUAD_NODE_ID( n1, n2, n3, n4 ) MIN_EDGE_NODE_ID( (MIN_EDGE_NODE_ID( n1, n2 )), (MIN_EDGE_NODE_ID( n3, n4 )) )

    CubitStatus status = CUBIT_SUCCESS;
    DLIList< int* > *quad_hash = new DLIList< int*>[num_file_quads];
    int * quad_data = new int [num_file_quads * 5 ];

    int i;
    for ( i = 0; i < num_file_quads; i++ )
    {
        int *data = &(quad_data[i*5]);

        data[0] = quad_nodes_file[i*4];
        data[1] = quad_nodes_file[i*4 + 1];
        data[2] = quad_nodes_file[i*4 + 2];
        data[3] = quad_nodes_file[i*4 + 3];
        data[4] = i;

        int key = MIN_QUAD_NODE_ID( data[0], data[1], data[2], data[3] );
        key %= num_file_quads;
        quad_hash[key].append( data );
    }

    for ( i = 0; i < numQuad; i++ )
    {
        CubitBoolean found = CUBIT_FALSE;
        int this_quad[4];
        
        constructQuadVerts( &quadEdge[i*4], edgeVert, this_quad );

        int key = MIN_QUAD_NODE_ID( this_quad[0], this_quad[1], this_quad[2], this_quad[3] );
        key %= num_file_quads;
        for ( int j = 0; j < quad_hash[key].size(); j++ )
        {
            int *data = quad_hash[key].get_and_step();

            if ( this_quad[0] == data[0] &&
                 this_quad[1] == data[1] &&
                 this_quad[2] == data[2] &&
                 this_quad[3] == data[3] )
            {
                quadVertNorms[ i*4*3     ] = quad_node_norms_file[ data[4]*4*3     ];
                quadVertNorms[ i*4*3 + 1 ] = quad_node_norms_file[ data[4]*4*3 + 1 ];
                quadVertNorms[ i*4*3 + 2 ] = quad_node_norms_file[ data[4]*4*3 + 2 ];
                quadVertNorms[ i*4*3 + 3 ] = quad_node_norms_file[ data[4]*4*3 + 3 ];
                quadVertNorms[ i*4*3 + 4 ] = quad_node_norms_file[ data[4]*4*3 + 4 ];
                quadVertNorms[ i*4*3 + 5 ] = quad_node_norms_file[ data[4]*4*3 + 5 ];
                quadVertNorms[ i*4*3 + 6 ] = quad_node_norms_file[ data[4]*4*3 + 6 ];
                quadVertNorms[ i*4*3 + 7 ] = quad_node_norms_file[ data[4]*4*3 + 7 ];
                quadVertNorms[ i*4*3 + 8 ] = quad_node_norms_file[ data[4]*4*3 + 8 ];
                quadVertNorms[ i*4*3 + 9 ] = quad_node_norms_file[ data[4]*4*3 + 9 ];
                quadVertNorms[ i*4*3 + 10] = quad_node_norms_file[ data[4]*4*3 + 10];
                quadVertNorms[ i*4*3 + 11] = quad_node_norms_file[ data[4]*4*3 + 11];
                found = CUBIT_TRUE;
                break;
            }
            else if ( this_quad[0] == data[1] &&
                      this_quad[1] == data[2] &&
                      this_quad[2] == data[3] &&
                      this_quad[3] == data[0] )
            {
                quadVertNorms[ i*4*3     ] = quad_node_norms_file[ data[4]*4*3 + 3 ];
                quadVertNorms[ i*4*3 + 1 ] = quad_node_norms_file[ data[4]*4*3 + 4 ];
                quadVertNorms[ i*4*3 + 2 ] = quad_node_norms_file[ data[4]*4*3 + 5 ];
                quadVertNorms[ i*4*3 + 3 ] = quad_node_norms_file[ data[4]*4*3 + 6 ];
                quadVertNorms[ i*4*3 + 4 ] = quad_node_norms_file[ data[4]*4*3 + 7 ];
                quadVertNorms[ i*4*3 + 5 ] = quad_node_norms_file[ data[4]*4*3 + 8 ];
                quadVertNorms[ i*4*3 + 6 ] = quad_node_norms_file[ data[4]*4*3 + 9 ];
                quadVertNorms[ i*4*3 + 7 ] = quad_node_norms_file[ data[4]*4*3 + 10 ];
                quadVertNorms[ i*4*3 + 8 ] = quad_node_norms_file[ data[4]*4*3 + 11 ];
                quadVertNorms[ i*4*3 + 9 ] = quad_node_norms_file[ data[4]*4*3 + 0 ];
                quadVertNorms[ i*4*3 + 10] = quad_node_norms_file[ data[4]*4*3 + 1 ];
                quadVertNorms[ i*4*3 + 11] = quad_node_norms_file[ data[4]*4*3 + 2 ];
                found = CUBIT_TRUE;
                break;
            }
            else if ( this_quad[0] == data[2] &&
                      this_quad[1] == data[3] &&
                      this_quad[2] == data[0] &&
                      this_quad[3] == data[1] )
            {
                quadVertNorms[ i*4*3     ] = quad_node_norms_file[ data[4]*4*3 + 6 ];
                quadVertNorms[ i*4*3 + 1 ] = quad_node_norms_file[ data[4]*4*3 + 7 ];
                quadVertNorms[ i*4*3 + 2 ] = quad_node_norms_file[ data[4]*4*3 + 8 ];
                quadVertNorms[ i*4*3 + 3 ] = quad_node_norms_file[ data[4]*4*3 + 9 ];
                quadVertNorms[ i*4*3 + 4 ] = quad_node_norms_file[ data[4]*4*3 + 10 ];
                quadVertNorms[ i*4*3 + 5 ] = quad_node_norms_file[ data[4]*4*3 + 11 ];
                quadVertNorms[ i*4*3 + 6 ] = quad_node_norms_file[ data[4]*4*3 + 0 ];
                quadVertNorms[ i*4*3 + 7 ] = quad_node_norms_file[ data[4]*4*3 + 1 ];
                quadVertNorms[ i*4*3 + 8 ] = quad_node_norms_file[ data[4]*4*3 + 2 ];
                quadVertNorms[ i*4*3 + 9 ] = quad_node_norms_file[ data[4]*4*3 + 3 ];
                quadVertNorms[ i*4*3 + 10] = quad_node_norms_file[ data[4]*4*3 + 4 ];
                quadVertNorms[ i*4*3 + 11] = quad_node_norms_file[ data[4]*4*3 + 5 ];
                found = CUBIT_TRUE;
                break;
            }
            else if (  this_quad[0] == data[3] &&
                       this_quad[1] == data[0] &&
                       this_quad[2] == data[1] &&
                       this_quad[3] == data[2] )
            {
                quadVertNorms[ i*4*3     ] = quad_node_norms_file[ data[4]*4*3 + 9 ];
                quadVertNorms[ i*4*3 + 1 ] = quad_node_norms_file[ data[4]*4*3 + 10 ];
                quadVertNorms[ i*4*3 + 2 ] = quad_node_norms_file[ data[4]*4*3 + 11 ];
                quadVertNorms[ i*4*3 + 3 ] = quad_node_norms_file[ data[4]*4*3 + 0 ];
                quadVertNorms[ i*4*3 + 4 ] = quad_node_norms_file[ data[4]*4*3 + 1 ];
                quadVertNorms[ i*4*3 + 5 ] = quad_node_norms_file[ data[4]*4*3 + 2 ];
                quadVertNorms[ i*4*3 + 6 ] = quad_node_norms_file[ data[4]*4*3 + 3 ];
                quadVertNorms[ i*4*3 + 7 ] = quad_node_norms_file[ data[4]*4*3 + 4 ];
                quadVertNorms[ i*4*3 + 8 ] = quad_node_norms_file[ data[4]*4*3 + 5 ];
                quadVertNorms[ i*4*3 + 9 ] = quad_node_norms_file[ data[4]*4*3 + 6 ];
                quadVertNorms[ i*4*3 + 10] = quad_node_norms_file[ data[4]*4*3 + 7 ];
                quadVertNorms[ i*4*3 + 11] = quad_node_norms_file[ data[4]*4*3 + 8 ];
                found = CUBIT_TRUE;
                break;
            }
        }
        if ( found == CUBIT_FALSE )
        {
            status = CUBIT_FAILURE;
            break;
        }
    }

    delete [] quad_data;
    delete [] quad_hash;
    return status;
}

//===========================================================================
//  Function: constructTriVerts
//  Purpose:  Given an array of edge indices for a tri element and an array
//            of edge to vertex info, construct an array containing the vertices
//            on the input tri.
//  Date:     08/31/2004
//  Author:   mlstate
//===========================================================================
static void constructTriVerts
(
    int triEdge[3],  // <I> Array of edges on this tri.
    int *edgeVert,   // <I> Array of edge to vertex info.
    int triVert[3]   // <O> The vertices on the input tri.
)
{
    int *verts_edge1 = &(edgeVert[triEdge[0]*2]);
    int *verts_edge2 = &(edgeVert[triEdge[1]*2]);
    int *verts_edge3 = &(edgeVert[triEdge[2]*2]);

    if ( verts_edge1[0] != verts_edge2[0] &&
         verts_edge1[0] != verts_edge2[1] )
    {
        triVert[0] = verts_edge1[0];
        triVert[1] = verts_edge1[1];
    }
    else
    {
        triVert[0] = verts_edge1[1];
        triVert[1] = verts_edge1[0];
    }
    if ( verts_edge3[0] != verts_edge2[0] &&
         verts_edge3[0] != verts_edge2[1] )
    {
        triVert[2] = verts_edge3[1];
    }
    else
    {
        triVert[2] = verts_edge3[0];
    }
}

//===========================================================================
//  Function: constructQuadVerts
//  Purpose:  Given an array of edge indices for a quad element and an array
//            of edge to vertex info, construct an array containing the vertices
//            on the input quad.
//  Date:     08/31/2004
//  Author:   mlstate
//===========================================================================
static void constructQuadVerts
(
    int QuadEdge[4],  // <I> Array of edges on this quad.
    int *edgeVert,    // <I> Array of edge to vertex info.
    int QuadVert[4]   // <O> The vertices on the input quad.
)
{
    int *verts_edge1 = &(edgeVert[QuadEdge[0]*2]);
    int *verts_edge2 = &(edgeVert[QuadEdge[1]*2]);
    int *verts_edge3 = &(edgeVert[QuadEdge[2]*2]);
    int *verts_edge4 = &(edgeVert[QuadEdge[3]*2]);

    if ( verts_edge1[0] != verts_edge2[0] &&
         verts_edge1[0] != verts_edge2[1] )
    {
        QuadVert[0] = verts_edge1[0];
        QuadVert[1] = verts_edge1[1];
    }
    else
    {
        QuadVert[0] = verts_edge1[1];
        QuadVert[1] = verts_edge1[0];
    }

    if ( verts_edge3[0] != verts_edge4[0] &&
         verts_edge3[0] != verts_edge4[1] )
    {
        QuadVert[2] = verts_edge3[0];
        QuadVert[3] = verts_edge3[1];
    }
    else
    {
        QuadVert[2] = verts_edge3[1];
        QuadVert[3] = verts_edge3[0];
    }
}

//===========================================================================
//  Function: readMBGNormalFile
//  Purpose:  Read the raw data out of the Normal/Tangents file.
//  Date:     08/31/2004
//  Author:   mlstate
//===========================================================================
CubitStatus readMBGNormalFile
(
    const char *filename,
    int *num_tris,                // <O>
    int *num_quads,               // <O>
    int *num_edges,               // <O>
    int **edge_nodes,             // <OF> len = 2 * num_edges
    double **edge_node_tangents,  // <OF> len = 3 * 2 * num_edges
    int **tri_nodes,              // <OF> len = 3 * num_tris
    double **tri_node_normals,    // <OF> len = 3 * 3 * num_tris
    int **quad_nodes,             // <OF> len = 4 * num_quads
    double **quad_node_normals    // <OF> len = 3 * 4 * num_quads
)
{
#define MAX_FILE_LINE 512

    int num_t = 0;
    int num_q = 0;
    int curr_quad = 0;
    int curr_tri = 0;
    int num_e = 0;
    int n[4];

    if ( num_tris           ) *num_tris = 0;
    if ( num_quads          ) *num_quads = 0;
    if ( num_edges          ) *num_edges = 0;
    if ( edge_nodes         ) *edge_nodes = NULL;
    if ( edge_node_tangents ) *edge_node_tangents = NULL;
    if ( tri_nodes          ) *tri_nodes = NULL;
    if ( quad_nodes         ) *quad_nodes = NULL;
    if ( tri_node_normals   ) *tri_node_normals = NULL;
    if ( quad_node_normals  ) *quad_node_normals = NULL;

    if ( strlen( filename ) == 0 )
    {
        PRINT_ERROR(" No filename specified\n");
        return CUBIT_FAILURE;
    }

    ifstream geom_file( filename );

    char fileline[MAX_FILE_LINE];

    while ( geom_file.getline( fileline, MAX_FILE_LINE ) )
    {
        int block_id,
            num_corners,
            num_elems;

        if ( !strcmp( fileline, END_OF_BLOCKS ) ) break;
        if ( 3 != sscanf( fileline, "BLK %d %d %d",
                          &block_id, &num_corners, &num_elems ) )
        {
            PRINT_ERROR( "MBG Normal file has incorrect file format\n" );
            geom_file.close();
            return CUBIT_FAILURE;
        }

        if ( num_corners == 3 )
        {
            checkMemoryAllocations( num_corners, num_elems, &num_t, tri_nodes,
                                    tri_node_normals );
        }
        else if ( num_corners == 4 )
        {
            checkMemoryAllocations( num_corners, num_elems, &num_q, quad_nodes,
                                    quad_node_normals );
        }
        else
        {
            PRINT_ERROR( "MBG Normal file can only support quads and tris.\n" );
            geom_file.close();
            return CUBIT_FAILURE;
        }

        for ( int ielem = 0; ielem < num_elems; ielem++ )
        {
            int *this_nodes = NULL;
            double *this_norms = NULL;

            if ( !geom_file.getline( fileline, MAX_FILE_LINE ) )
            {
                PRINT_ERROR( "MBG Normal file has incorrect file format\n" );
                geom_file.close();
                return CUBIT_FAILURE;
            }

            if ( num_corners == 4 )
            {
                if ( quad_nodes        ) this_nodes = &((*quad_nodes       )[4     * curr_quad]);
                if ( quad_node_normals ) this_norms = &((*quad_node_normals)[4 * 3 * curr_quad]);
                curr_quad++;
            }
            else if ( num_corners == 3 )
            {
                if ( tri_nodes        ) this_nodes = &((*tri_nodes       )[3     * curr_tri]);
                if ( tri_node_normals ) this_norms = &((*tri_node_normals)[3 * 3 * curr_tri]);
                curr_tri++;
            }
            if ( ( num_corners == 4 && 4 != sscanf( fileline, "%d %d %d %d",
                                                    &n[0], &n[1], &n[2], &n[3]  ) ) ||
                 ( num_corners == 3 && 3 != sscanf( fileline, "%d %d %d",
                                                    &n[0], &n[1], &n[2] ) ) )
            {
                PRINT_ERROR( "MBG Normal file has incorrect file format\n" );
                geom_file.close();
                return CUBIT_FAILURE;
            }

            if ( this_nodes )
                memcpy( this_nodes, n, num_corners * sizeof( int ) );

            // Read in the surface normals at the vertices of this element.

            for ( int icorner = 0; icorner < num_corners; icorner++ )
            {
                if ( !geom_file.getline( fileline, MAX_FILE_LINE ) )
                {
                    PRINT_ERROR( "MBG Normal file has incorrect file format\n" );
                    geom_file.close();
                    return CUBIT_FAILURE;
                }
                if ( this_norms )
                {
                    double norm[3];

                    if ( 3 != sscanf( fileline, "%lf %lf %lf", &norm[0], &norm[1], &norm[2] ) )
                    {
                        PRINT_ERROR( "MBG Normal file has incorrect file format\n" );
                        geom_file.close();
                        return CUBIT_FAILURE;
                    }
                    memcpy( &(this_norms[icorner * 3]), norm, 3 * sizeof( double ) );
                }
            }
        }
    }

    if ( !geom_file.getline( fileline, MAX_FILE_LINE ) )
    {
        PRINT_ERROR( "MBG Normal file has incorrect file format\n" );
        geom_file.close();
        return CUBIT_FAILURE;
    }

    if ( 1 != sscanf( fileline, "%d", &num_e ) )
    {
        PRINT_ERROR( "MBG Normal file has incorrect file format\n" );
        geom_file.close();
        return CUBIT_FAILURE;
    }

    if ( edge_nodes || edge_node_tangents )
    {
        if ( edge_nodes         ) *edge_nodes         = new int [ 2 * num_e ];
        if ( edge_node_tangents ) *edge_node_tangents = new double [3 * 2 * num_e];

        for ( int iedge = 0; iedge < num_e; iedge++ )
        {
            if ( !geom_file.getline( fileline, MAX_FILE_LINE ) )
            {
                PRINT_ERROR( "MBG Normal file has incorrect file format\n" );
                geom_file.close();
                return CUBIT_FAILURE;
            }
            if ( 2 != sscanf( fileline, "%d %d", &n[0], &n[1] ) )
            {
                PRINT_ERROR( "MBG Normal file has incorrect file format\n" );
                geom_file.close();
                return CUBIT_FAILURE;
            }

            for ( int iend = 0; iend < 2; iend++ )
            {
                double d[3];

                if ( edge_nodes )
                {
                    (*edge_nodes)[iedge * 2 + iend ] = n[iend];
                }

                if ( !geom_file.getline( fileline, MAX_FILE_LINE ) )
                {
                    PRINT_ERROR( "MBG Normal file has incorrect file format\n" );
                    geom_file.close();
                    return CUBIT_FAILURE;
                }
                if ( 3 != sscanf( fileline, "%lf %lf %lf", &d[0], &d[1], &d[2] ) )
                {
                    PRINT_ERROR( "MBG Normal file has incorrect file format\n" );
                    geom_file.close();
                    return CUBIT_FAILURE;
                }

                int index = ((iedge*2) + iend) * 3;
                if ( edge_node_tangents )
                {
                    (*edge_node_tangents)[ index     ] = d[0];
                    (*edge_node_tangents)[ index + 1 ] = d[1];
                    (*edge_node_tangents)[ index + 2 ] = d[2];
                }
            }
        }
    }

    if ( num_quads ) *num_quads = num_q;
    if ( num_tris  ) *num_tris  = num_t;
    if ( num_edges ) *num_edges = num_e;

    geom_file.close();
    return CUBIT_SUCCESS;
}

//===========================================================================
//  Function: checkMemoryAllocations
//  Purpose:  While reading data out of the Normal/Tangents file, reallocate
//            Arrays to ensure that we have enough memory.
//  Date:     08/31/2004
//  Author:   mlstate
//===========================================================================
static void checkMemoryAllocations
(
    int num_nodes_per_elem,
    int additional_num_elems,
    int *num_elems,
    int **nodes,
    double **node_normals
)
{
    int old_num = *num_elems;

    *num_elems += additional_num_elems;
    if ( node_normals && *node_normals )
    {
        double *new_mem = new double [ *num_elems * 3 * num_nodes_per_elem ];
        memcpy( new_mem, *node_normals,  old_num * 3 * num_nodes_per_elem * sizeof( double ) );
        delete *node_normals;
        *node_normals = new_mem;
    }
    else if ( node_normals )
    {
        *node_normals = new double [ *num_elems * num_nodes_per_elem * 3 ];
    }

    if ( nodes && *nodes )
    {
        int *new_mem = new int [ *num_elems * num_nodes_per_elem ];
        memcpy( new_mem, *nodes,  old_num * num_nodes_per_elem * sizeof( int ) );
        delete *nodes;
        *nodes = new_mem;
    }
    else if ( nodes )
    {
        *nodes = new int [ *num_elems * num_nodes_per_elem ];
    }
}

// This code classify the type of edge with respect to convexity,
// C0 or C1 continuity, and topology.  The returned int uses the following
// number code:
// -1 = non-convex, manifold edge
//  0 = non-feature, manifold  edge (ie, C1 edge)
//  1 = convex, manifold edge
//  2 = non-manifold edge
//  Numbers higher than two are error codes
//  3 = edge only attached to a single facet (boundary of 2-d sheet).
static int classify_local_convexity_at_edge(CubitFacetEdge *edge_ptr)
{
    //easy check... if not a feature return 0
  if(!edge_ptr->is_feature())
    return 0;
  
    //else we know that we are at feature... determine the type
  CubitFacet *adj_facet_1 = NULL, *adj_facet_2 = NULL;
  CubitPoint *point_1 = NULL, *point_2 = NULL, *point_3 = NULL;
  CubitVector v2, v3;
  CubitVector v0 = edge_ptr->point(0)->coordinates();
  CubitVector v1 = edge_ptr->point(1)->coordinates();
  DLIList<CubitFacet*> adj_facets;
  edge_ptr->facets(adj_facets);
  
    //if non-manifold edge, mark as two
  if(adj_facets.size() > 2){
    return 2;
  }
    //if boundary on 2-d (shouldn't happen) mark as 3
    //anything above 2 is essentially an error code...
  else if(adj_facets.size() < 2){
    return 3;
  }
    //otherwise, determine whether it is 1 or -1
    //this else should be a separate function:  classify_edge_convexity()
  adj_facet_1 = adj_facets.get_and_step();
  int edge_index = adj_facet_1->edge_index(edge_ptr);
    //order such that *_1 is the facet that uses
    // this edge as the forward edge...
  if(adj_facet_1->edge_use(edge_index) == 1)
    adj_facet_2 = adj_facets.get();
  else{
    adj_facet_1 = adj_facets.get();
    adj_facets.reset();
    adj_facet_2 = adj_facets.get();
  }
    // order the points so that the Jacobian of the matrix
    // will tell use whether we are convex or not at this edge.
  adj_facet_1->points(point_1, point_2, point_3);
  if(point_1 == edge_ptr->point(0)){
    v2 = point_3->coordinates();
  }
  else if(point_1 == edge_ptr->point(1)){
    v2 = point_2->coordinates();
  }
  else{
    v2 = point_1->coordinates();
  }
  adj_facet_2->points(point_1, point_2, point_3);
  if(point_1 == edge_ptr->point(0)){
    v3 = point_2->coordinates();
  }
  else if(point_1 == edge_ptr->point(1)){
    v3 = point_3->coordinates();
  }
  else{
    v3 = point_1->coordinates();
  }

  v1-=v0;
  v2-=v0;
  v3-=v0;
    //if Jacobian is negative, we are convex here.
  if( ( v1 % (v2* v3) ) < 0.0){
    return 1;
  }
  else{
    return -1;
  }

}
/*
int test_edges(int num_edge, int* edge_vert, double* my_vert, CubitFacetEdge** edge_array){
    int ii;
    int return_val = 0.0;
  CubitFacetEdgeData *cfed_ptr = NULL;
  int my_v1, my_v2;
  CubitVector vert1, vert2;
  CubitPoint *poi1, *poi2;
  CubitVector coord1, coord2;
  for(ii=0; ii<num_edge; ++ii){
    cfed_ptr = CAST_TO(edge_array[ii], CubitFacetEdgeData);
    poi1 = cfed_ptr->start_node();
    poi2 = cfed_ptr->end_node();
    coord1.set(poi1->x(),poi1->y(),poi1->z());
    coord2.set(poi2->x(),poi2->y(),poi2->z());
    my_v1 = edge_vert[ii*2];
    my_v2 = edge_vert[ii*2+1];
    vert1.set(my_vert[my_v1*3],my_vert[my_v1*3+1],my_vert[my_v1*3+2]);
    vert2.set(my_vert[my_v2*3],my_vert[my_v2*3+1],my_vert[my_v2*3+2]);
    double my_dist = vert1.distance_between(coord1);
    if(my_dist > .01){
        PRINT_INFO("Edge (%i)(1), distance (%f)\n",ii,my_dist);
        return_val++;
    }
    my_dist = vert2.distance_between(coord2);
    if(my_dist > .01){
        PRINT_INFO("Edge (%i)(2), distance (%f)\n",ii,my_dist);
        return_val++;
    }
  }
  return return_val;
} 

  
*/
