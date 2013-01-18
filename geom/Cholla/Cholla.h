/** \file Cholla.h
     
//- Class:       Cholla
//- Description: C-style Interface for the Cholla module
//-              Facet-based geometry definition and evaluation
//- Owner:       Steven J. Owen
//- Checked by:
//- Version:
*/

#ifndef CHOLLA____H
#define CHOLLA____H

#define NUM_EDGE_CPTS 3
#define NUM_TRI_CPTS  6
#define NUM_QUAD_CPTS 15

#include "CubitDefines.h"
/**
 * Return control points for a quartic Bezier for each face, using
 * a given feature angle tolerance to distinguish C0 and C1 edges.
 * Separate arrays of tris and quad faces should be supplied
 *
 * @param angle         - (IN) feature angle tolerance (degrees) 
 * @param numTri        - (IN) number of quad faces
 * @param numQuad       - (IN) number of triangle faces
 * @param numEdge       - (IN) number of edges
 * @param numVert       - (IN) number of vertices
 * @param triEdge       - (IN) array of face to edge numbering,
 *                             length 3 * numTri
 * @param quadEdge      - (IN) array of face to edge numbering,
 *                             length 4 * numQuad
 * @param edgeVert      - (IN) array of edge to vertex numbering,
 *                             length 2*numEdge
 * @param vert          - (IN) array of vertex coordinates,
 *                             length 3*numVert
 * @param edgeCtrlPts   - (OUT) array of control points computed at the edges
 *                              length 3*(3*numEdge)
 * @param triCtrlPts    - (OUT) array of control points computed at the triangles
 *                              length 3*(6*numTri)
 * @param quadCtrlPts   - (OUT) array of control points computed at the quads
 *                              length 3*(15*numQuad)
 */
void constructBezier( double angle, int numTri, int numQuad, int numEdge,
                      int numVert, int* triEdge, int* quadEdge,
                      int* edgeVert, double* vert,
                      double* edgeCtrlPts, double* triCtrlPts,
                      double* quadCtrlPts );

/**
 * Return normal and tangents for on a list of triangles
 * a given feature angle tolerance to distinguish C0 and C1 edges.
 * Separate arrays of tris and quad faces should be supplied
 *
 * @param angle         - (IN) feature angle tolerance (degrees) 
 * @param numTri        - (IN) number of tri faces
 * @param numEdge       - (IN) number of edges
 * @param numVert       - (IN) number of vertices
 * @param triEdge       - (IN) array of face to edge numbering,
 *                             length 3 * numTri
 * @param edgeVert      - (IN) array of edge to vertex numbering,
 *                             length 2*numEdge
 * @param vert          - (IN) array of vertex coordinates,
 *                             length 3*numVert
 * @param triVertNorms  - (OUT) array of normals at the triangle vertices
 *                              length 3*(3*numTri)
 * @param edgeVertTang  - (OUT) array of tangents at the edge vertices
 *                              length 3*(2*numEdge)
 * @param edgeTypeFlags - (OUT) array of flags labeling the type of each edge,
 *                              length numEdge \n
 *                              2 = non-manifold edge \n
 *                              1 = concave feature edge \n
 *                              0 = non-feature edge \n
 *                              -1= convex feature edge \n
 */
void constructTriNormals( double angle, int numTri, int numEdge,
                          int numVert, int* triEdge,
                          int* edgeVert, double* vert,
                          double* triVertNorms, double* edgeVertTang,
                          int* edgeTypeFlags);

/**
 * Return normal and tangents for on a list of quad faces
 * a given feature angle tolerance to distinguish C0 and C1 edges.
 * Separate arrays of tris and quad faces should be supplied
 *
 * @param angle         - (IN) feature angle tolerance (degrees) 
 * @param numQuad       - (IN) number of quad faces
 * @param numEdge       - (IN) number of edges
 * @param numVert       - (IN) number of vertices
 * @param quadEdge      - (IN) array of face to edge numbering,
 *                             length 4 * numQuad
 * @param edgeVert      - (IN) array of edge to vertex numbering,
 *                             length 2*numEdge
 * @param vert          - (IN) array of vertex coordinates,
 *                             length 3*numVert
 * @param quadVertNorms - (OUT) array of normals at the quad vertices
 *                              length 3*(4*numQuad)
 * @param edgeVertTang  - (OUT) array of tangents at the edge vertices
 *                              length 3*(2*numEdge)
 * @param edgeTypeFlags - (OUT) array of flags labeling the type of each edge,
 *                              length numEdge \n
 *                              2 = non-manifold edge \n
 *                              1 = concave feature edge \n
 *                              0 = non-feature edge \n
 *                              -1= convex feature edge \n
 */
void constructQuadNormals( double angle, int numQuad, int numEdge,
                           int numVert, int* quadEdge,
                           int* edgeVert, double* vert,
                           double* quadVertNorms, double* edgeVertTang,
                           int* edgeTypeFlags);

/**
 * Return normal and tangents for on a list of triangles
 * a given  a file containing the geometric normals and tangents.
 * Separate arrays of tris and quad faces should be supplied
 *
 * @param filename      - (IN) full path and name of the text file containing
 *                             the geometric surface normals and tangents.
 *                             This file is output from Cubit when exporting
 *                             an exodus file.
 * @param numTri        - (IN) number of tri faces
 * @param triVert       - (IN) array of face to vertex numbering,
 *                             length 3 * numTri
 * @param triVertNorms  - (OUT) array of normals at the triangle vertices
 *                              length 3*(3*numTri)
 * @param edgeVertTang  - (OUT) array of tangents at the edge vertices
 *                              length 3*(2*numEdge)
 */
CubitStatus constructTriNormalsFromFile( const char *filename,
                                         int numTri,
                                         int numEdge,
                                         int* triEdge,
                                         int* edgeVert,
                                         double* triVertNorms,
                                         double* edgeVertTang );

/**
 * Return normal and tangents for on a list of quad faces
 * a given a file containing the geometric normals and tangents.
 * Separate arrays of tris and quad faces should be supplied
 *
 * @param filename      - (IN) full path and name of the text file containing
 *                             the geometric surface normals and tangents.
 *                             This file is output from Cubit when exporting
 *                             an exodus file.
 * @param numQuad       - (IN) number of quad faces
 * @param numEdge       - (IN) number of edges
 * @param quadEdge      - (IN) array of face to edge numbering,
 *                             length 4 * numQuad
 * @param edgeVert      - (IN) array of edge to vertex numbering,
 *                             length 2*numEdge
 * @param quadVertNorms - (OUT) array of normals at the quad vertices
 *                              length 3*(4*numQuad)
 * @param edgeVertTang  - (OUT) array of tangents at the edge vertices
 *                              length 3*(2*numEdge)
 */
CubitStatus constructQuadNormalsFromFile( const char *filename,
                                          int numQuad,
                                          int numEdge,
                                          int* quadEdge,
                                          int* edgeVert,
                                          double* quadVertNorms,
                                          double* edgeVertTang );

/**
 * Return normal and tangents for on a list of quad faces
 * a given a file containing the geometric normals and tangents.
 * This returns all of the data in the file rather than just for the
 * specified input elements and edges like constructQuadNormalsFromFile
 * and constructTriNormalsFromFile
 *
 * @param filename      - (IN) full path and name of the text file containing
 *                             the geometric surface normals and tangents.
 *                             This file is output from Cubit when exporting
 *                             an exodus file.
 * @param numTri        - (OUT) number of quad faces
 * @param numQuad       - (OUT) number of quad faces
 * @param numEdge       - (OUT) number of edges
 * @param edgeVert      - (OUT) array of vertices on each edge.
 *                              length 2*numEdge (Must be freed by caller)
 * @param edgeVertTang  - (OUT) array of tangents at the edge vertices
 *                              length 3*(2*numEdge) (Must be freed by caller)
 * @param TriVert       - (OUT) array of vertices on each tri.
 *                              length 3*numTri
 * @param TriVertNorms  - (OUT) array of normals at the tri vertices
 *                              length 3*(4*numTri) (Must be freed by caller)
 * @param QuadVert      - (OUT) array of vertices on each quad.
 *                              length 3*numQuad
 * @param QuadVertNorms - (OUT) array of normals at the quad vertices
 *                              length 3*(4*numQuad) (Must be freed by caller)
 */
CubitStatus readMBGNormalFile( const char *filename,
                               int *numTri,
                               int *numQuad,
                               int *numEdge,
                               int **edgeVert,
                               double **edgeVertTang,
                               int **triVert,
                               double **triVertNorms,
                               int **quadVert,
                               double **quadVertNorms );

/**
 * Evaluate a set of quartic Bezier patches at a specified parametric location
 * given its control points. Evaluates location, normal and/or derivative.
 * Expects list of either quads or tris (not both)
 *
 * @param numFace       - (IN) number of faces to evaluate
 * @param numEdge       - (IN) number of edges on the faces
 * @param numVert       - (IN) number of vertices
 * @param numVertPerFace- (IN) has value of of 3 or 4, for tris or quads
 * @param faceEdge      - (IN) array of face to edge indecies (into edgeVert)
 *                             length numVertPerFace * numFace
 * @param edgeVert      - (IN) array of edge to vertex indecies (into vert)
 *                             length 2 * numEdge
 * @param vert          - (IN) array of vertex locations
 *                             length 3 * numVert
 * @param faceCtrlPts   - (IN) array of control points for the faces (same
 *                             order as faceEdge) 
 *                             if (numFacePerVert == 3) ncp = 6
 *                             if (numFacePerVert == 4) ncp = 15
 *                             length = 3 * (ncp * numFace)
 * @param edgeCtrlPts   - (IN) array of control points for the edges (same
 *                             order as edgeVert array)
 *                             length 3 * (3 * numEdge)
 * @param numLocs       - (IN) number of locations to evaluate per face
 * @param paramLocation - (IN) parametric location to evaluate on each face
 *                             if (numFacePerVert == 3) (u,v,w)
 *                               length = 3 * numLocs  
 *                             if (numFacePerVert == 4) (u,v)
 *                               length = 2 * numLocs
 * @param location      - (OUT) evaluated point on surface at paramLocation
 *                              length = numFace * 3 * numLocs
 * @param normal        - (IN/OUT) evaluated normal at parmLocation. If NULL
 *                              then will not be evaluated. Return vector 
 *                              is normalized
 *                              length = numFace * 3 * numLocs
 * @param deriv         - (IN/OUT) evaluated derivative at parmLocation with 
 *                              respect to uv or uvw system. If NULL
 *                              then will not be evaluated.
 *                              if (numFacePerVert == 3) 
 *                                length = numFace * numLocs * 3*3 
 *                                             (du/dx,du/dy,du/dz, 
 *                                              dv/dx,dv/dy,dv/dz,
 *                                              dw/dx,dw/dz,dw/dz)
 *                              if (numFacePerVert == 4) 
 *                                length = numFace * numLocs * 3*2 
 *                                             (du/dx,du/dy,du/dz, 
 *                                              dv/dx,dv/dy,dv/dz)                          
 */
void evalBezierFace( int numFace, int numEdge, int numVert, int numVertPerFace, 
                     int* faceEdge, int* edgeVert, 
                     double* vert,  double* faceCtrlPts,
                     double* edgeCtrlPts, 
                     int numLocs, double* paramLocation, 
                     double* location, double* normal, double* deriv );

/**
 * This is the same as the evalBezierFace except it differes by the input
 * arguments.  This function takes a list of normals and tangents at the face
 * vertices and computes bezier control points internally.  Normals and tangents
 * should have been computed previously in constructTriNormals
 * or constructQuad Normals.  This function is not as computationally efficient
 * as evalBezierFace since it requires Bezier control points to be computed
 * as part of the call - however, it is more memory efficient, requiring fewer
 * variables to be stored with the calling application. The following argument 
 * list describes only those arguments that differ from evalBezierFace above.
 *
 * @param vertNorms -     (IN) array of normals at the face vertices
 *                              length 3*(numVertPerFace*numFace)
 * @param edgeVertTang  - (OUT) array of tangents at the edge vertices
 *                              length 3*(2*numEdge)
 */
void evalBezierFaceFromNorms( int numFace, int numEdge, int numVert, 
                              int numVertPerFace, int* faceEdge, int* edgeVert, 
                              double* vert, double* vertNorms, double* edgeVertTang,
                              int numLocs, double* paramLocation, 
                              double* location, double* normal, double* deriv );

/**
 * Project a set of x-y-z locations to Bezier patches.  Finds the closest point on
 * one of the patches.  Evaluates location, normal and/or derivative.
 * Expects list of either quads or tris (not both)
 *
 * @param numFace       - (IN) number of faces to evaluate
 * @param numEdge       - (IN) number of edges on the faces
 * @param numVert       - (IN) number of vertices
 * @param numVertPerFace- (IN) has value of of 3 or 4, for tris or quads
 * @param faceEdge      - (IN) array of face to edge indecies (into edgeVert)
 *                             length numVertPerFace * numFace
 * @param edgeVert      - (IN) array of edge to vertex indecies (into vert)
 *                             length 2 * numEdge
 * @param vert          - (IN) array of vertex locations
 *                             length 3 * numVert
 * @param faceCtrlPts   - (IN) array of control points for the faces (same
 *                             order as faceEdge) 
 *                             if (numFacePerVert == 3) ncp = 6
 *                             if (numFacePerVert == 4) ncp = 15
 *                             length = 3 * (ncp * numFace)
 * @param edgeCtrlPts   - (IN) array of control points for the edges (same
 *                             order as edgeVert array)
 *                             length 3 * (3 * numEdge)
 * @param numLocs       - (IN) number of locations to evaluate per face
 * @param xyz           - (IN) xyz locations to evaluate on each face                             
 *                               length = 3 * numLocs                             
 * @param specify_tol   - (IN) 1 = a converge_tol is specified,
 *                             0 = a converge_tol is not specified.  A default
 *                                 tol will be computed and used.
 * @param converge_tol  - (IN) The convergance tolerance for this projection.
 *                             Will be ignored if specify_tol == 0
 * @param xyzOnFace     - (OUT) locations on patches closest to input points xyz.
 *                              length = 3 * numLocs
 * @param normal        - (IN/OUT) evaluated normals at xyz. If NULL
 *                              then will not be evaluated. Return vector 
 *                              is normalized
 *                              length = 3 * numLocs
 * @param deriv         - (IN/OUT) evaluated derivatives at xyz with 
 *                              respect to uv or uvw system. If NULL
 *                              then will not be evaluated.
 *                              if (numFacePerVert == 3) 
 *                                length = numLocs * 3*3 
 *                                             (du/dx,du/dy,du/dz, 
 *                                              dv/dx,dv/dy,dv/dz,
 *                                              dw/dx,dw/dz,dw/dz)
 *                              if (numFacePerVert == 4) 
 *                                length = numLocs * 3*2 
 *                                             (du/dx,du/dy,du/dz, 
 *                                              dv/dx,dv/dy,dv/dz)                          
 */
void projToBezierFace( int numFace, int numEdge, int numVert, int numVertPerFace, 
                     int* faceEdge, int* edgeVert, 
                     double* vert,  double* faceCtrlPts,
                     double* edgeCtrlPts, 
                     int numLocs, double* xyz, 
                     int specify_tol, double converge_tol,
                     double* xyzOnFace, double* normal, double* deriv );

/**
 * This is the same as the projToBezierFace except it differes by the input
 * arguments.  This function takes a list of normals and tangents at the face
 * vertices and computes bezier control points internally.  Normals and tangents
 * should have been computed previously in constructTriNormals
 * or constructQuadNormals.  This function is not as computationally efficient
 * as evalBezierFace since it requires Bezier control points to be computed
 * as part of the call - however, it is more memory efficient, requiring fewer
 * variables to be stored with the calling application. The following argument 
 * list describes only those arguments that differ from projToBezierFace above.
 *
 * @param vertNorms -     (IN) array of normals at the face vertices
 *                              length 3*(numVertPerFace*numFace)
 * @param edgeVertTang  - (OUT) array of tangents at the edge vertices
 *                              length 3*(2*numEdge)
 */
void projToBezierFaceFromNorms( int numFace, int numEdge, int numVert, int numVertPerFace, 
                     int* faceEdge, int* edgeVert, 
                     double* vert,  double* vertNorms,
                     double* edgeVertTang, 
                     int numLocs, double* xyz, 
                     int specify_tol, double converge_tol,
                     double* xyzOnFace, double* normal, double* deriv );

/**
 * Evaluate a quartic Bezier curve at a specified parametric location
 * given its tangents at the end-points. Evaluates location, 
 * and/or tangent.
 *
 * @param numEdge       - (IN) number of edges on the faces
 * @param numVert       - (IN) number of vertices
 * @param edgeVert      - (IN) array of edge to vertex indecies (into vert)
 *                             length 2 * numEdge
 * @param vert          - (IN) array of vertex locations
 *
 * @param edgeVertTang  - (IN) array of tangents at the edge vertices
 *                              length 3*(2*numEdge)
 * @param numLocs       - (IN) number of locations to evaluate per edge
 * @param paramLocation - (IN) parametric location to evaluate on each edge
 *                             length = numLocs
 * @param location      - (OUT) evaluated point on edge at paramLocation
 *                              length = numEdge * 3 * numLocs
 * @param tangent       - (IN/OUT) evaluated tangent at parmLocation. If NULL
 *                              then will not be evaluated. Return vector 
 *                              is normalized
 *                              length = numEdge * 3 * numLocs                
 */
void evalBezierEdgeFromTans( int numEdge, int numVert, int* edgeVert, 
                              double* vert, double* edgeVertTang,
                              int numLocs, double* paramLocation, 
                              double* location, double* tangent );
/**
 * Evaluate a quartic Bezier curve at a specified parametric location
 * given its control points. Evaluates location, and/or tangent.
 *
 * @param numEdge       - (IN) number of edges on the faces
 * @param numVert       - (IN) number of vertices
 * @param edgeVert      - (IN) array of edge to vertex indecies (into vert)
 *                             length 2 * numEdge
 * @param vert          - (IN) array of vertex locations
 *                             length 3 * numVert
 * @param edgeCtrlPts   - (IN) array of control pt locations for each edge in
 *                              the bezier.  The beziers are quartic, so there
 *                              are 3 control pts per curve in addition to the
 *                              2 end pts for each curve for a total of 5 pts
 *                              per edge.
 *                              length = 3 * ( 3 * numEdge )
 * @param numLocs       - (IN) number of locations to evaluate per edge
 * @param paramLocation - (IN) parametric location to evaluate on each edge
 *                             length = numLocs
 * @param location      - (OUT) evaluated point on edge at paramLocation
 *                              length = numEdge * 3 * numLocs
 * @param tangent       - (IN/OUT) evaluated tangent at parmLocation. If NULL
 *                              then will not be evaluated. Return vector 
 *                              is normalized
 *                              length = numEdge * 3 * numLocs                
 */
void evalBezierEdge( int numEdge, int numVert, int* edgeVert,
                     double* vert, double* edgeCtrlPts, 
                     int numLocs, double* paramLocation, 
                     double* location, double* tangent );

/**
 * Project to a quartic Bezier curve at a specified parametric location
 * given its tangents at the end-points. Evaluates location, 
 * and/or tangent.  Finds only one point; in rare cases, more than one
 * closest point might exist.
 *
 * @param numEdge       - (IN) number of edges on the faces
 * @param numVert       - (IN) number of vertices
 * @param edgeVert      - (IN) array of edge to vertex indecies (into vert)
 *                             length 2 * numEdge
 * @param vert          - (IN) array of vertex locations
 *                              length = 3 * NumVerts
 * @param edgeVertTang  - (IN) array of tangents at the edge vertices
 *                              length 3*(2*numEdge)
 * @param numLocs       - (IN) number of locations to evaluate per edge
 * @param xyz           - (IN) xyz locations to project to each edge        
 *                               length = 3 * numLocs                             
 * @param xyzOnEdge     - (OUT) locations on edges closest to input points xyz.
 *                              length = numEdge * 3 * numLocs
 * @param tangent       - (IN/OUT) evaluated tangent at parmLocation. If NULL
 *                              then will not be evaluated. Return vector 
 *                              is normalized
 *                              length = numEdge * 3 * numLocs                
 */
void projToBezierEdgeFromTans( int numEdge, int numVert, int* edgeVert, 
                              double* vert, double* edgeVertTang,
                              int numLocs, double* xyz, 
                              double* xyzOnEdge, double* tangent );

/**
 * Project a 3d point to a quartic Bezier curve at a specified parametric
 * location given its control points.
 *
 * @param numEdge       - (IN) number of edges on the faces
 * @param numVert       - (IN) number of vertices
 * @param edgeVert      - (IN) array of edge to vertex indecies (into vert)
 *                             length 2 * numEdge
 * @param vert          - (IN) array of vertex locations
 *                              length = 3 * numVert
 * @param edgeCtrlPts   - (IN) array of control pt locations for each edge in
 *                              the bezier.  The beziers are quartic, so there
 *                              are 3 control pts per curve in addition to the
 *                              2 end pts for each curve for a total of 5 pts
 *                              per edge.
 *                              length = 3 * ( 3 * numEdge )
 * @param numLocs       - (IN) number of locations to evaluate per edge
 * @param xyz           - (IN) xyz locations to evaluate on each edge        
 *                               length = 3 * numLocs                             
 * @param xyzOnEdge     - (OUT) locations on edges closest to input points xyz.
 *                              length = numEdge * 3 * numLocs
 * @param tangent       - (IN/OUT) evaluated tangent at parmLocation. If NULL
 *                              then will not be evaluated. Return vector 
 *                              is normalized
 *                              length = numEdge * 3 * numLocs 
 * @param t_value       - (OUT) value of the parameter at the projected point               
 *                              length = numEdge * numLocs
 */
void projToBezierEdge( int numEdge, int numVert, int* edgeVert, 
                       double* vert, double* edgeCtrlPts,
                       int numLocs, double* xyz, 
                       double* xyzOnEdge, double* tangent, double* t_value );

/**
 * Functions for debugging
 **/

/* 
 * 
 *
 * 
 * 
 */
/* dump the face mesh to a file that can be read by CUBIT
 * Use the same definition of parameters as used with resolveFaceVectors
 *
 * @param fileName      - name of file to write to (will be overwritten)
 * @param includeResults- if (0) then don't write face/edge vectors to file
 */
void dumpMesh(const char *fileName, int includeResults, double angle,
              int numTri, int numQuad, int numEdge,
              int numVert, int* triEdge, int* quadEdge,
              int* edgeVert, double* vert,
              double* edgeCtrlPts, double* triCtrlPts,
              double* quadCtrlPts);

/* dump only the results to a file
 */
void dumpResults(const char *fileName, int numEdge, int numTri, int numQuad, 
                 double* edgeCtrlPts, double* triCtrlPts, 
                 double* quadCtrlPts );

/* import the face mesh from the specified file
 * This is used by Cubit for debugging purposes
 *
 * @param fileName      - name of file to read from.  
 *                        (file typically written by dumpMesh)
 */
void importMesh(const char *fileName, int *resultsIncluded,
                double *angle, int *numTri, int *numQuad, int *numEdge,
                int *numVert, int** triEdge, int** quadEdge,
                int** edgeVert, double** vert,
                double** edgeCtrlPts, double** triCtrlPts,
                double** quadCtrlPts);

/* read only the results from a file
 */
void importResults(const char *fileName, int *numEdge, int *numTri, 
                   int *numQuad, double** edgeCtrlPts, double** triCtrlPts, 
                   double** quadCtrlPts );


#endif
