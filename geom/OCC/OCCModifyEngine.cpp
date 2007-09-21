//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

// Filename      : OCCModifyEngine.cpp
//
// Purpose       : ModifyEngine for faceted geometry
//
// Special Notes : Modeled after GeometryModifyEngine and AcisModifyEngine.
//
// Creator       : John Fowler
//
// Creation Date : 6/02
//
// Owner         : John Fowler
//-------------------------------------------------------------------------

#include "OCCModifyEngine.hpp"
#include "OCCQueryEngine.hpp"

#include "CubitMessage.hpp"
#include "CubitDefines.h"

#include "CubitUtil.hpp"
#include "CubitPoint.hpp"
#include "CubitPointData.hpp"
#include "CubitFacet.hpp"
#include "CubitQuadFacet.hpp"
#include "CubitFacetData.hpp"
#include "CubitFacetEdge.hpp"
#include "CubitFacetEdgeData.hpp"
#include "GeometryQueryTool.hpp"
#include "GeometryModifyTool.hpp"
#include "ChollaSurface.hpp"
#include "ChollaCurve.hpp"
#include "ChollaPoint.hpp"
#include "OCCPoint.hpp"
#include "OCCCurve.hpp"
#include "CurveFacetEvalTool.hpp"
#include "FacetEvalTool.hpp"
#include "OCCCoEdge.hpp"
#include "OCCLoop.hpp"
#include "OCCSurface.hpp"
#include "OCCShell.hpp"
#include "OCCLump.hpp"
#include "OCCBody.hpp"
#include "ChollaEngine.hpp"
#include "TDGeomFacet.hpp"
#include "CubitFileIOWrapper.hpp"
#include "TDFacetBoundaryPoint.hpp"
#include "Cholla.h"
#include "Body.hpp"
#include "GfxDebug.hpp"
#include "RefFace.hpp"
#include "FacetDataUtil.hpp"
#include "FBDataUtil.hpp"
#include "FBIntersect.hpp"
#include "IntegerHash.hpp"
#include "FacetboolInterface.hpp"
#include "CpuTimer.hpp"
#include "ProgressTool.hpp"
#include "AppUtil.hpp"
#include "SphereEvaluator.hpp"
#include "CylinderEvaluator.hpp"

OCCModifyEngine* OCCModifyEngine::instance_ = 0;

//===============================================================================
// Function   : OCCModifyEngine
// Member Type: PUBLIC
// Description: constructor
// Author     : John Fowler
// Date       : 10/02
//===============================================================================
OCCModifyEngine::OCCModifyEngine()
{
//  assert( !instance_ );

    // add this modify engine to geometrymodifytool
  GeometryModifyTool::instance()->add_gme(this);
}


//===============================================================================
// Function   : ~OCCModifyEngine
// Member Type: PUBLIC
// Description: destructor
// Author     : John Fowler
// Date       : 10/02
//===============================================================================
OCCModifyEngine::~OCCModifyEngine() 
{
        instance_ = 0;
}

//===============================================================================
// Function   : make_Point
// Member Type: PUBLIC
// Description: make a geometric entity point
// Author     : John Fowler
// Date       : 10/02
//===============================================================================
Point* OCCModifyEngine::make_Point( CubitVector const& /*point*/) const
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return (Point*) NULL;
}

//===============================================================================
// Function   : make_Curve
// Member Type: PUBLIC
// Description: make a curve
// Author     : John Fowler
// Date       : 10/02
//===============================================================================
Curve* OCCModifyEngine::make_Curve(Curve * /*curve_ptr*/) const
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return (Curve*) NULL;
}

//===============================================================================
// Function   : make_Curve
// Member Type: PUBLIC
// Description: make a curve
// Author     : John Fowler
// Date       : 10/02
//===============================================================================
Curve* OCCModifyEngine::make_Curve( Point const* /*point1_ptr*/,
                             Point const* /*point2_ptr*/,
                             Surface* /*ref_face_ptr*/,
                             const CubitVector * /*third_point*/) const
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return (Curve*) NULL;
}

//===============================================================================
// Function   : make_Curve
// Member Type: PUBLIC
// Description: make a curve
// Author     : John Fowler
// Date       : 10/02
//===============================================================================
Curve* OCCModifyEngine::make_Curve( GeometryType /*curve_type*/,
                             Point const* /*point1_ptr*/,
                             Point const* /*point2_ptr*/,
                             DLIList<CubitVector*>& /*vector_list*/,
                             Surface* /*ref_face_ptr*/) const
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return (Curve*) NULL;
}

//===============================================================================
// Function   : make_Curve
// Member Type: PUBLIC
// Description: make a curve
// Author     : John Fowler
// Date       : 10/02
//===============================================================================
Curve* OCCModifyEngine::make_Curve( GeometryType /*curve_type*/,
                             Point const* /*point1_ptr*/,
                             Point const* /*point2_ptr*/,
                             CubitVector const* /*intermediate_point_ptr*/,
                             CubitSense /*sense*/) const
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return (Curve*) NULL;
}

//===============================================================================
// Function   : make_Surface
// Member Type: PUBLIC
// Description: make a surface
// Author     : John Fowler
// Date       : 10/02
//===============================================================================
Surface* OCCModifyEngine::make_Surface( Surface * /*old_surface_ptr*/,
                                 CubitBoolean /*extended_from*/) const
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return (Surface*) NULL;
}

//===============================================================================
// Function   : make_Surface
// Member Type: PUBLIC
// Description: make a surface
// Author     : John Fowler
// Date       : 10/02
//===============================================================================
Surface* OCCModifyEngine::make_Surface( GeometryType /*surface_type*/,
                                 DLIList<Curve*>& /*curve_list*/,
                                 Surface * /*old_surface_ptr*/,
                                 bool /*check_edges*/) const
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return (Surface*) NULL;
}

//===============================================================================
// Function   : make_Lump
// Member Type: PUBLIC
// Description: make a lump
// Author     : John Fowler
// Date       : 10/02
//===============================================================================
Lump* OCCModifyEngine::make_Lump( DLIList<Surface*>& /*surface_list*/ ) const
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return (Lump*) NULL;
}

//===============================================================================
// Function   : make_BodySM
// Member Type: PUBLIC
// Description: make a BodySM
// Author     : John Fowler
// Date       : 10/02
//===============================================================================
BodySM* OCCModifyEngine::make_BodySM( Surface * ) const
    {return NULL ;}

//===============================================================================
// Function   : make_BodySM
// Member Type: PUBLIC
// Description: make a BodySM
// Author     : John Fowler
// Date       : 10/02
//===============================================================================
BodySM* OCCModifyEngine::make_BodySM( DLIList<Lump*>& /*lump_list*/ ) const
    {return NULL ;}

//===============================================================================
// Function   : fillinedge
// Member Type: PRIVATE
// Description: put points on this edge of triangle to be refined
// Author     : John Fowler
// Date       : 11/02
//===============================================================================
void OCCModifyEngine::fillinedge( 
  int *edge, 
  int numpointsonanedge, 
  double radius, 
  STD(vector)<CubitPoint *>& points) const
{
  int numintervals, i;
  double xbegin, ybegin, zbegin, xend, yend, zend, xincr, yincr, zincr;
  double dist;
  CubitPoint *new_point;
  
  numintervals = numpointsonanedge - 1;
  if ( numintervals < 1 ) return;
  xbegin = points[edge[0]]->x();
  ybegin = points[edge[0]]->y();  
  zbegin = points[edge[0]]->z();
  xend = points[edge[numintervals]]->x();
  yend = points[edge[numintervals]]->y();  
  zend = points[edge[numintervals]]->z();  
  xincr = (xend-xbegin)/(double)(numintervals);  
  yincr = (yend-ybegin)/(double)(numintervals);  
  zincr = (zend-zbegin)/(double)(numintervals);  
  for ( i = 0; i < numintervals - 1; i++ ) {
    xbegin += xincr; 
    ybegin += yincr;
    zbegin += zincr;

    dist = sqrt(xbegin*xbegin + ybegin*ybegin + zbegin*zbegin);
    new_point = (CubitPoint *) new CubitPointData( xbegin*radius/dist, ybegin*radius/dist, 
    							zbegin*radius/dist ); 
    edge[i+1] = points.size(); // get the point number
    points.push_back(new_point);   
  }

}

//===============================================================================
// Function   : refinetriangle
// Member Type: PRIVATE
// Description: add internal points and make connections for this triangle
// Author     : John Fowler
// Date       : 11/02
//===============================================================================
void OCCModifyEngine::refinetriangle(
  int level, 
  int numpointsonanedge, 
  int *iedge1, 
  int *iedge2, 
  int *iedge3,
  int isign1, 
  int isign2, 
  int isign3, 
  double radius,
  STD(vector)<CubitPoint *>& points,
  DLIList<CubitFacet *>& facet_list) const
{
int i, numintervals, icount, j, jcount;
int ntris;
int iedge1cnt, iedge2cnt, iedge3cnt;

int increment, trissofar, i1s, i2s, i1e;  
int nverts1, nverts;
double dlev;
int *vertnumarray, index, iskip;
CubitFacet *facet_ptr;
CubitPoint *new_point;
		double x1inc, y1inc, z1inc, x2inc, y2inc, z2inc;
		double xstart, ystart, zstart, xend, yend, zend, dist;
			double xinternalinc, yinternalinc, zinternalinc;
  
  numintervals = numpointsonanedge - 1;
  iedge1cnt = ( isign1 == 1 ) ? 0 : numintervals;
  iedge2cnt = ( isign2 == 1 ) ? 0 : numintervals;
  iedge3cnt = ( isign3 == 1 ) ? 0 : numintervals;

	index = points.size();
	dlev = 1. + (double)level;
	nverts1 = (int)(0.5*(dlev+1.)*dlev);

	ntris = (level+1)*(level+1);
	dlev += 1.;	
	nverts = (int)(0.5*(dlev+1.)*dlev);
	vertnumarray = new int[nverts];

//  Put the point numbers for the interior points into vertnumarray.
	if ( numintervals > 2 ) {  // numintervals must be at least 3 to make interior points
		jcount = 1;
		icount = 2;
		for ( i = 2; i < numintervals; i++ ) {
			icount += 2;
			for ( j = 0; j < jcount; j++ ) {
				vertnumarray[icount] = index++;
				icount += 1;
			}
			jcount += 1;		
		}	
	}	
	i = 3;

	iskip = 2;
	vertnumarray[0] = iedge1[iedge1cnt]; iedge1cnt += isign1;
	vertnumarray[1] = iedge1[iedge1cnt]; iedge1cnt += isign1;
	iedge2cnt += isign2;
	vertnumarray[2] = iedge2[iedge2cnt]; iedge2cnt += isign2;

	while ( i < nverts1 ) {
		vertnumarray[i] = iedge1[iedge1cnt]; iedge1cnt += isign1;
		vertnumarray[i+iskip] = iedge2[iedge2cnt]; iedge2cnt += isign2;
		i += iskip+1;
		iskip += 1;	
	}
	for ( i = nverts1; i < nverts; i++ ) {
		vertnumarray[i] = iedge3[iedge3cnt]; iedge3cnt += isign3;
	}	
	

//!  Make the internal points, and put them on the sphere.

	if ( numintervals > 2 ) {
		int i1first, i1last, i2first, i2last;
		if ( isign1 == 1 ) {
			i1first = 0; i1last = numintervals;
		} else {
			i1last = 0; i1first = numintervals;
		}
		if ( isign2 == 1 ) {
			i2first = 0; i2last = numintervals;
		} else {
			i2last = 0; i2first = numintervals;
		}
		x1inc = (points[iedge1[i1last]]->x() - points[iedge1[i1first]]->x())/numintervals;	
		y1inc = (points[iedge1[i1last]]->y() - points[iedge1[i1first]]->y())/numintervals;	
		z1inc = (points[iedge1[i1last]]->z() - points[iedge1[i1first]]->z())/numintervals;	
		x2inc = (points[iedge2[i2last]]->x() - points[iedge2[i2first]]->x())/numintervals;	
		y2inc = (points[iedge2[i2last]]->y() - points[iedge2[i2first]]->y())/numintervals;	
		z2inc = (points[iedge2[i2last]]->z() - points[iedge2[i2first]]->z())/numintervals;	

		icount = 2;
		jcount = 1;
		for ( i = 2; i < numintervals; i++ ) {
		xstart = points[iedge1[i1first]]->x() + (double)i*x1inc;
		ystart = points[iedge1[i1first]]->y() + (double)i*y1inc;
		zstart = points[iedge1[i1first]]->z() + (double)i*z1inc;
		xend = points[iedge2[i2first]]->x() + (double)i*x2inc;
		yend = points[iedge2[i2first]]->y() + (double)i*y2inc;
		zend = points[iedge2[i2first]]->z() + (double)i*z2inc;
			xinternalinc = (xend-xstart)/(icount);    
			yinternalinc = (yend-ystart)/(icount);    
			zinternalinc = (zend-zstart)/(icount);    
				for ( j = 0; j < jcount; j++ ) {
					xstart += xinternalinc;
					ystart += yinternalinc;
					zstart += zinternalinc;

    					dist = sqrt(xstart*xstart + ystart*ystart + zstart*zstart);
  					new_point = (CubitPoint *) new CubitPointData( xstart*radius/dist, 
  											ystart*radius/dist,
   											zstart*radius/dist );
  					points.push_back(new_point);					
				}
			icount += 1;
			jcount += 1;
		}	
	
	}

//!  Make the connections.
	increment = trissofar = 0;
	i1s = 0;
	i2s = 1;
	while ( trissofar < ntris ) {
		i1e = i1s + increment;
		while ( i1s < i1e) {
  			facet_ptr = new CubitFacetData( points[vertnumarray[i1s]],
  							points[vertnumarray[i2s]], 
							points[vertnumarray[i2s+1]] );
  			facet_list.append( facet_ptr );
  			facet_ptr = new CubitFacetData( points[vertnumarray[i1s]],
  							points[vertnumarray[i2s+1]], 
							points[vertnumarray[i1s+1]] );
  			facet_list.append( facet_ptr );
			i1s++;
			i2s++;
			trissofar += 2;
		}
  			facet_ptr = new CubitFacetData( points[vertnumarray[i1s]],
  							points[vertnumarray[i2s]], 
							points[vertnumarray[i2s+1]] );
  			facet_list.append( facet_ptr );
		increment++;
		trissofar++;
		i1s++;
		i2s += 2;
	}
	delete [] vertnumarray;
}

//===============================================================================
// Function   : sphere
// Member Type: PUBLIC
// Description: build a sphere with facets
// Author     : John Fowler
// Date       : 10/02
//===============================================================================
BodySM* OCCModifyEngine::sphere(double radius) const
{
  CubitStatus rv = CUBIT_SUCCESS;
  DLIList <CubitFacet *>facet_list;
  DLIList <CubitPoint *>point_list;
  CubitPoint *new_point;
  double feature_angle;
  int i, interp_order;
  CubitBoolean smooth_non_manifold, split_surfaces;
  BodySM *body_ptr = NULL;
  STD(vector)<CubitPoint *> points;
//! We use a vector for the points because when refining the triangles, we will
//! need fast random access to them.  After all of the points have been made,
//! they will go onto the DLIList point_list.  Then the points vector will be
//! deleted.
//! xp and yp are the fundamental units for an icosahedron.
  const double xp = 0.525731112119133606*radius, zp = 0.850650808352039932*radius;
//!  Makes an icosahedron and then refines the triangles as they are made.  Number
//!  of levels is given by the variable "level". 
  int level, numpointsonanedge;
  int *edge[30];
//! 20 faces in an icosahedron.  const int basefacets[][] holds the connections.
//! 30 edges in an icosahedron.  const int edgeendverts[][] holds the edge endpoints.
/*
  const int basefacets[20][3] = { {0,1,2},{0,2,3},{3,2,4},{2,5,4},
  				  {2,1,5},{5,1,6},{5,6,7},{4,5,7},
				  {4,7,8},{8,7,9},{9,7,6},{9,6,10},
				  {9,10,11},{11,10,0},{0,10,1},{10,6,1},
				  {3,11,0},{3,8,11},{3,4,8},{9,11,8} };
*/
  const int edgeendverts[30][2] = { {0,1},{0,2},{0,3},{0,11},{0,10},{1,2},
  			            {1,5},{1,6},{1,10},{2,3},{2,4},{2,5},
			            {3,4},{3,8},{3,11},{4,5},{4,7},{4,8},
			            {5,6},{5,7},{6,7},{6,9},{6,10},{7,8},
			            {7,9},{8,9},{8,11},{9,10},{9,11},{10,11} };
//!  triedges[][] holds the three edges for each of the 20 triangles.  A minus sign means
//!  that the edge will be traversed backward.  To accommodate a zero, one has been added
//!  to the magnitude of the edge.  So, for example, the first triangle has edges 0, 5, and 1,
//!  with the last being traversed from end to beginning.
  const int triedges[20][3] = { {1,2,6},{2,3,10},{-10,13,11},{12,11,-16,},
  				{-6,12,7,},{-7,19,8,},{19,20,21},{16,17,20},
				{17,18,24,},{-24,26,25},{-25,-22,-21},{-22,28,23},
				{28,29,30},{-30,-4,-5},{5,1,-9},{-23,-9,-8},
				{15,-3,-4},{14,15,27},{13,14,18},{29,-26,-27} };
//!  12 points in an icosahedron.  svert[][] holds these. 
  const double svert[12][3] = { {-xp,0.,zp}, {xp,0.,zp}, {0.,zp,xp}, {-zp, xp, 0.},
  				{0.,zp,-xp}, {zp,xp,0.}, {zp,-xp,0.}, {xp,0.,-zp},
				{-xp,0.,-zp}, {0.,-zp,-xp}, {0.,-zp,xp}, {-zp,-xp,0.} };
    
  level = 7;  // gives a sphere with 642 vertices and 1280 triangles.
  	      // Eventually this should be user-selectable.
  numpointsonanedge = 2 + level;  

  for ( i = 0; i < 30; i++ ) { // make the edges
    edge[i] = new int[numpointsonanedge];
    edge[i][0] = edgeendverts[i][0];
    edge[i][numpointsonanedge-1] = edgeendverts[i][1];
  }
  
  for ( i = 0; i < 12; i++ ) { // make the icosahedron vertices
    new_point = (CubitPoint *) new CubitPointData( svert[i][0],svert[i][1],svert[i][2] );
    points.push_back(new_point);
  }
   
  for ( i = 0; i < 30; i++ ) { // put points on the edges
    fillinedge(edge[i], numpointsonanedge, radius, points); 
  }

 int sign1, sign2, sign3, edg1, edg2, edg3; 
  for ( i = 0; i < 20; i++ ) { // refine the 20 triangles
    edg1 = ( triedges[i][0] > 0 ) ? triedges[i][0] - 1 : -triedges[i][0] - 1;
    edg2 = ( triedges[i][1] > 0 ) ? triedges[i][1] - 1 : -triedges[i][1] - 1;
    edg3 = ( triedges[i][2] > 0 ) ? triedges[i][2] - 1 : -triedges[i][2] - 1;
//! sign1, etc., says in which direction to traverse an edge
    sign1 = ( triedges[i][0] > 0 ) ? 1 : -1;
    sign2 = ( triedges[i][1] > 0 ) ? 1 : -1;
    sign3 = ( triedges[i][2] > 0 ) ? 1 : -1;
    
    refinetriangle(level,numpointsonanedge,edge[edg1],edge[edg2],edge[edg3],
    		sign1,sign2,sign3,radius,points,facet_list);
  }

//! Put the points in point_list and then delete the points vector.	 
  for ( unsigned int z = 0; z < points.size(); z++ ) {
    point_list.append(points[z]);
  }

  points.clear();

  for ( i = 0; i < 30; i++ ) delete[] edge[i];  
  feature_angle = -1.0;
  interp_order = 0;
  smooth_non_manifold = CUBIT_TRUE;
  split_surfaces = CUBIT_FALSE;

  ChollaEngine *cholla_ptr = NULL;
  OCCModifyEngine *fme = const_cast<OCCModifyEngine *> (this);
  rv = fme->build_cholla_surfaces( facet_list,
                                   point_list,
                                   feature_angle,
                                   interp_order,
                                   smooth_non_manifold,
                                   split_surfaces,
                                   cholla_ptr );
  if ( rv == CUBIT_SUCCESS )
  {
      CubitEvaluatorData **sphere_data;

      set_sphere_eval_data( cholla_ptr, radius, sphere_data );

      finish_facet_Body( cholla_ptr,
                         (const CubitEvaluatorData **)sphere_data,
                         feature_angle,
                         interp_order,
                         body_ptr);

      if ( cholla_ptr )
      {
         cholla_ptr->delete_me();
         delete cholla_ptr;
      }
      if ( sphere_data[0] ) delete sphere_data[0];
      delete sphere_data;
  }
  return body_ptr;
}

void OCCModifyEngine::set_sphere_eval_data
(
    ChollaEngine *cholla_ptr,
    double radius,
    CubitEvaluatorData **&eval_data ) const
{
    eval_data = new CubitEvaluatorData*;

    SphereEvaluatorData *data = new SphereEvaluatorData;
    data->radius = radius;
    data->center.set( 0.0, 0.0, 0.0 );
    eval_data[0] = data;
}

//===============================================================================
// Function   : brick
// Member Type: PUBLIC
// Description: build a brick with facets
// Author     : John Fowler
// Date       : 10/02
//===============================================================================
BodySM* OCCModifyEngine::brick( double wid, double dep, double hi ) const
{
  CubitStatus rv = CUBIT_SUCCESS;
  double xmin, xmax, ymin, ymax, zmin, zmax;
  DLIList <CubitFacet *>facet_list;
  DLIList <CubitPoint *>point_list;
  CubitPoint *new_point;
  CubitFacet *facet_ptr;
  int i, numpoints, numtris;
  double feature_angle;
  int interp_order;
  CubitBoolean smooth_non_manifold, split_surfaces;
  BodySM *body_ptr = NULL;
  STD(vector)<CubitPoint *> points;
  
  numpoints = 14;
  numtris = 24;
  
  xmin = -0.5*wid;
  xmax = 0.5*wid;
  ymin = -0.5*dep;
  ymax = 0.5*dep;
  zmin = -0.5*hi;
  zmax = 0.5*hi;
  
  new_point = (CubitPoint *) new CubitPointData( xmin,ymin,zmin ); 
  points.push_back(new_point);
  new_point = (CubitPoint *) new CubitPointData( xmax,ymin,zmin ); 
  points.push_back(new_point);
  new_point = (CubitPoint *) new CubitPointData( xmax,ymin,zmax ); 
  points.push_back(new_point);  
  new_point = (CubitPoint *) new CubitPointData( xmin,ymin,zmax ); 
  points.push_back(new_point); 
  new_point = (CubitPoint *) new CubitPointData( xmin,ymax,zmin ); 
  points.push_back(new_point);
  new_point = (CubitPoint *) new CubitPointData( xmax,ymax,zmin ); 
  points.push_back(new_point);
  new_point = (CubitPoint *) new CubitPointData( xmax,ymax,zmax ); 
  points.push_back(new_point);  
  new_point = (CubitPoint *) new CubitPointData( xmin,ymax,zmax ); 
  points.push_back(new_point);   
  new_point = (CubitPoint *) new CubitPointData( 0.5*(xmin+xmax),0.5*(ymin+ymax),zmin ); 
  points.push_back(new_point);  
  new_point = (CubitPoint *) new CubitPointData( xmax,0.5*(ymin+ymax),0.5*(zmin+zmax) ); 
  points.push_back(new_point);    
  new_point = (CubitPoint *) new CubitPointData( 0.5*(xmin+xmax),0.5*(ymin+ymax),zmax ); 
  points.push_back(new_point);  
  new_point = (CubitPoint *) new CubitPointData( xmin,0.5*(ymin+ymax),0.5*(zmin+zmax) ); 
  points.push_back(new_point);   
  new_point = (CubitPoint *) new CubitPointData( 0.5*(xmin+xmax),ymin,0.5*(zmin+zmax) ); 
  points.push_back(new_point);   
  new_point = (CubitPoint *) new CubitPointData( 0.5*(xmin+xmax),ymax,0.5*(zmin+zmax) ); 
  points.push_back(new_point);   
  
  for ( i = 0; i < numpoints; i++ ) {
    point_list.append(points[i]);
  }    
  
  // bottom face      
  facet_ptr = new CubitFacetData( points[0],points[1], points[12] );
  facet_list.append( facet_ptr );
  facet_ptr = new CubitFacetData( points[1],points[2], points[12] );
  facet_list.append( facet_ptr );
  facet_ptr = new CubitFacetData( points[2],points[3], points[12] );
  facet_list.append( facet_ptr );
  facet_ptr = new CubitFacetData( points[3],points[0], points[12] );
  // back face
  facet_list.append( facet_ptr );
  facet_ptr = new CubitFacetData( points[1],points[0], points[8] );
  facet_list.append( facet_ptr );
  facet_ptr = new CubitFacetData( points[5],points[1], points[8] );
  facet_list.append( facet_ptr );  
  facet_ptr = new CubitFacetData( points[4],points[5], points[8] );
  facet_list.append( facet_ptr );
  facet_ptr = new CubitFacetData( points[0],points[4], points[8] );
  facet_list.append( facet_ptr );
  // left face
  facet_ptr = new CubitFacetData( points[0],points[3], points[11] );
  facet_list.append( facet_ptr );
  facet_ptr = new CubitFacetData( points[3],points[7], points[11] );
  facet_list.append( facet_ptr );
  facet_ptr = new CubitFacetData( points[7],points[4], points[11] );
  facet_list.append( facet_ptr );
  facet_ptr = new CubitFacetData( points[4],points[0], points[11] );
  facet_list.append( facet_ptr ); 
  // top face
  facet_ptr = new CubitFacetData( points[7],points[6], points[13] );
  facet_list.append( facet_ptr );
  facet_ptr = new CubitFacetData( points[6],points[5], points[13] );
  facet_list.append( facet_ptr );
  facet_ptr = new CubitFacetData( points[5],points[4], points[13] );
  facet_list.append( facet_ptr );
  facet_ptr = new CubitFacetData( points[4],points[7], points[13] );
  facet_list.append( facet_ptr ); 
  // right face
  facet_ptr = new CubitFacetData( points[1],points[5], points[9] );
  facet_list.append( facet_ptr );
  facet_ptr = new CubitFacetData( points[5],points[6], points[9] );
  facet_list.append( facet_ptr );
  facet_ptr = new CubitFacetData( points[6],points[2], points[9] );
  facet_list.append( facet_ptr );
  facet_ptr = new CubitFacetData( points[2],points[1], points[9] );
  facet_list.append( facet_ptr ); 
  // front face  
  facet_ptr = new CubitFacetData( points[3],points[2], points[10] );
  facet_list.append( facet_ptr );
  facet_ptr = new CubitFacetData( points[2],points[6], points[10] );
  facet_list.append( facet_ptr );
  facet_ptr = new CubitFacetData( points[6],points[7], points[10] );
  facet_list.append( facet_ptr );
  facet_ptr = new CubitFacetData( points[7],points[3], points[10] );
  facet_list.append( facet_ptr );   
  
  points.clear(); //  clear out the points vector since we are through with it.
  
  feature_angle = 100.0;
  interp_order = 0;
  smooth_non_manifold = CUBIT_TRUE;
  split_surfaces = CUBIT_FALSE;
  
  ChollaEngine *cholla_ptr = NULL;
  OCCModifyEngine *fme = const_cast<OCCModifyEngine *> (this);
  rv = fme->build_cholla_surfaces( facet_list,
                                   point_list,
                                   feature_angle,
                                   interp_order,
                                   smooth_non_manifold,
                                   split_surfaces,
                                   cholla_ptr );
  if ( rv == CUBIT_SUCCESS )
  {
      finish_facet_Body( cholla_ptr,
                         NULL,
                         feature_angle,
                         interp_order,
                         body_ptr);
      if ( cholla_ptr )
      {
         cholla_ptr->delete_me();
         delete cholla_ptr;
      }

  }
  return body_ptr;
}

//===============================================================================
// Function   : finish_facet_Body
// Member Type: PRIVATE
// Description: general function for creating a facet-based Body given a closed
//              set of facets
// Author     : John Fowler
// Date       : 10/02
//===============================================================================
CubitStatus OCCModifyEngine::finish_facet_Body( ChollaEngine *&cholla_ptr,
                                                  const CubitEvaluatorData **eval_data,
                                                  double feature_angle,
                                                  int interp_order,
                                                  BodySM *&bodysm_ptr) const
{
  DLIList<Surface *> surface_list;
  ShellSM *shell_ptr;
  DLIList<ShellSM*> shell_list;
  Lump *lump_ptr;
  DLIList<Lump*> lump_list;
  CubitStatus rv;
  DLIList<ChollaSurface *> cholla_surface_list;
  DLIList<ChollaCurve *> cholla_curve_list;
  DLIList<ChollaPoint *> cholla_point_list;

  cholla_ptr->get_curves( cholla_curve_list );
  cholla_ptr->get_surfaces( cholla_surface_list );
  cholla_ptr->get_points( cholla_point_list );

  CubitBoolean use_feature_angle;
  if (feature_angle < 0.0)
    use_feature_angle = CUBIT_FALSE;
  else
    use_feature_angle = CUBIT_TRUE;
  
  GeometryQueryTool *gti = GeometryQueryTool::instance();  
  
  OCCModifyEngine *fme = const_cast<OCCModifyEngine *> (this);

  rv = fme->build_cholla_geometry( eval_data,
                                   cholla_surface_list,
                                   cholla_curve_list,
                                   cholla_point_list,
                                   use_feature_angle, 
                                   feature_angle,
                                   interp_order, surface_list );
  
  if ( surface_list.size() == 0 || rv != CUBIT_SUCCESS )
  {
    PRINT_ERROR("Problems building facet based surfaces.\n");
    rv = CUBIT_FAILURE;
    goto end_brick;
  }

  // make a body out of it
  rv = fme->make_facet_shell(surface_list,shell_ptr);
  if ( shell_ptr == NULL || rv != CUBIT_SUCCESS )
  {
    PRINT_ERROR("Problems building facet based shell entity.\n");
    rv = CUBIT_FAILURE;
    goto end_brick;
  } 
      //Set the sense for the surfaces (will be cofaces) on this shell.
      //Assumption: The sense is always forward when creating geom from facets.
      // (This may not be correct -especially with multiple shells in a body)
  int ii;
  OCCShell* facet_shell;
  facet_shell = CAST_TO( shell_ptr, OCCShell );
  for( ii = surface_list.size(); ii > 0; ii-- )
  {
    Surface* surf = surface_list.get_and_step();
    OCCSurface* facet_surf = CAST_TO( surf, OCCSurface );
    facet_surf->set_shell_sense( facet_shell, CUBIT_FORWARD );
  }
   			     
  shell_list.append(shell_ptr);
  rv = fme->make_facet_lump(shell_list,lump_ptr);
  if ( lump_ptr == NULL || rv != CUBIT_SUCCESS )
  {
    PRINT_ERROR("Problems building facet based lump entity.\n");
    rv = CUBIT_FAILURE;
    goto end_brick;
  }
  lump_list.append(lump_ptr);
  rv = fme->make_facet_body(lump_list,bodysm_ptr);
  
  if ( rv != CUBIT_SUCCESS )
  {
    PRINT_ERROR("Problems building facet based body entity.\n");
    rv = CUBIT_FAILURE;
    goto end_brick;
  }
   
  PRINT_INFO("Body successfully created.\n");
  PRINT_INFO("  Number of vertices = %d\n", gti->num_ref_vertices());
  PRINT_INFO("  Number of edges = %d\n", gti->num_ref_edges());
  PRINT_INFO("  Number of faces = %d\n", gti->num_ref_faces());
  PRINT_INFO("  Number of volumes = %d\n", gti->num_ref_volumes());
  PRINT_INFO("  Number of bodies = %d\n", gti->num_bodies());
   
  return rv;
  
end_brick:
  bodysm_ptr = (BodySM *)NULL;
  return rv;
  
}

//===============================================================================
// Function   : brick
// Member Type: PUBLIC
// Description: create a brick with facets given center axes and extension
// Author     : John Fowler
// Date       : 10/02
//===============================================================================
BodySM* OCCModifyEngine::brick( const CubitVector &/*center*/, 
                                  const CubitVector* /*axes[3]*/,
                                  const CubitVector &/*extension*/) const
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return (BodySM*) NULL;
}

//===============================================================================
// Function   : prism
// Member Type: PUBLIC
// Description: create a prism with facets
// Author     : John Fowler
// Date       : 10/02
//===============================================================================
BodySM* OCCModifyEngine::prism( double /*height*/, int /*sides*/, double /*major*/,
                               double /*minor*/) const
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return (BodySM*) NULL;
}

//===============================================================================
// Function   : pyramid
// Member Type: PUBLIC
// Description: create a pyramid with facets
// Author     : John Fowler
// Date       : 10/02
//===============================================================================
BodySM* OCCModifyEngine::pyramid( double /*height*/, int /*sides*/, double /*major*/,
                                 double /*minor*/, double /*top*/) const
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return (BodySM*) NULL;
}

//===============================================================================
// Function   : cylinder
// Member Type: PUBLIC
// Description: create a cylinder with facets
// Author     : John Fowler
// Date       : 10/02
//===============================================================================
BodySM* OCCModifyEngine::cylinder( double hi, double r1, double r2, double r3 ) const
{
  CubitStatus rv = CUBIT_SUCCESS;
  DLIList <CubitFacet *>facet_list;
  DLIList <CubitPoint *>point_list;
  CubitPoint *new_point;
  CubitFacet *facet_ptr;
  int i, j, numpoints;
  double feature_angle;
  int interp_order;
  CubitBoolean smooth_non_manifold, split_surfaces;
  BodySM *body_ptr = NULL;
  STD(vector)<CubitPoint *> points;
  int reslevel;  // relative level of resolution
  int nl;  // number of axial divisions
  int nr;  // number of radil divisions
  double cfac, rinc, linc;
  double x, y, z;
  int istart, iend, V3, pend;
  double zoffset, lpos, rpos, xrad, yrad;
  double rad_ratio = 1.0;

  reslevel = 4;
  nl = reslevel;
  nr = 8*reslevel;
    
  rinc = 360.0/(double)nr;
  linc = hi/(double)nl;
  cfac = CUBIT_PI/180.;
  
  if ( r3 > 0.0 ) {
    //  Cylinder:
    numpoints = (nl + 1)*nr + 2;
    istart = 0; iend = nl+1;
    V3 = (nl+1)*nr;
    pend = nl;
    rad_ratio = r2/r1;
  } else {
    //  Cone:
   	numpoints = nl*nr + 2;
    istart = 0; iend = nl;
    V3 = nl*nr;
    pend = nl-1;
  }
  
  //  Make the points.
  
  zoffset = 0.0;
  lpos = -0.5*hi; 
  xrad = r1;
  yrad = r2;
  for ( i = istart; i < iend; i++ ) {
    rpos = 10.0;
    xrad = zoffset*r3/hi + (hi - zoffset)*r1/hi;
    yrad = zoffset*r3*rad_ratio/hi + (hi - zoffset)*r2/hi;
    for ( j = 0; j < nr; j++ ) {
      x = xrad*cos(cfac*rpos);
      y = yrad*sin(cfac*rpos);
      z = lpos;
      new_point = (CubitPoint *) new CubitPointData( x,y,z );
      points.push_back(new_point);
      rpos += rinc;   
    }
    lpos += linc;
    zoffset += linc;
  } 
  //  Add the two apoint on the axis at the ends.
  new_point = (CubitPoint *) new CubitPointData( 0.,0.,-0.5*hi );
  points.push_back(new_point);
  new_point = (CubitPoint *) new CubitPointData( 0.,0.,0.5*hi );
  points.push_back(new_point);
  
  for ( i = 0; i < numpoints; i++ ) {
    point_list.append(points[i]);
  }  
  
  //  Make the triangles.
  int vertnum;
  vertnum = 0;
  for ( i = 0; i < pend; i++ ) {
    for ( j = 0; j < nr-1; j++ ) {
      facet_ptr = new CubitFacetData( points[vertnum+j],points[vertnum+j+1], points[vertnum+j+nr] );
      facet_list.append( facet_ptr );     
      facet_ptr = new CubitFacetData( points[vertnum+j+1],points[vertnum+j+1+nr], points[vertnum+j+nr] );
      facet_list.append( facet_ptr ); 
    }
    facet_ptr = new CubitFacetData( points[vertnum],points[vertnum+nr], points[vertnum+2*nr-1] );
    facet_list.append( facet_ptr );     
    facet_ptr = new CubitFacetData( points[vertnum+nr-1],points[vertnum], points[vertnum+2*nr-1] );
    facet_list.append( facet_ptr );     
    vertnum += nr;
  }
  
  //  Endcap(s)
  for ( i = 0; i < nr-1; i++ ) { // top cap
    facet_ptr = new CubitFacetData( points[vertnum+i],points[vertnum+i+1], points[V3+1] );
    facet_list.append( facet_ptr );  
  }   
  facet_ptr = new CubitFacetData( points[nr-1+vertnum],points[vertnum], points[V3+1] );
  facet_list.append( facet_ptr );    
  
  for ( i = 0; i < nr-1; i++ ) { // bottom cap
    facet_ptr = new CubitFacetData( points[i+1],points[i], points[V3] );
    facet_list.append( facet_ptr );  
  }   
  facet_ptr = new CubitFacetData( points[0],points[nr-1], points[V3] );
  facet_list.append( facet_ptr );    
  
  points.clear(); //  clear out the points vector since we are through with it.
  
  feature_angle = 135.0;
  interp_order = 0;
  smooth_non_manifold = CUBIT_TRUE;
  split_surfaces = CUBIT_FALSE;
  
  ChollaEngine *cholla_ptr = NULL;
  OCCModifyEngine *fme = const_cast<OCCModifyEngine *> (this);
  rv = fme->build_cholla_surfaces( facet_list,
                                   point_list,
                                   feature_angle,
                                   interp_order,
                                   smooth_non_manifold,
                                   split_surfaces,
                                   cholla_ptr );

  if ( rv == CUBIT_SUCCESS )
  {
      CubitEvaluatorData **cyl_data = NULL;

      set_cylinder_eval_data( cholla_ptr, hi, r1, r2, r3, cyl_data );
      finish_facet_Body( cholla_ptr,
                         (const CubitEvaluatorData**)cyl_data,
                         feature_angle,
                         interp_order,
                         body_ptr);
      if ( cholla_ptr )
      {
         cholla_ptr->delete_me();
         delete cholla_ptr;
      }
      if ( cyl_data )
      {
          if ( cyl_data[0] ) delete cyl_data[0];
          if ( cyl_data[1] ) delete cyl_data[1];
          if ( cyl_data[2] ) delete cyl_data[2];
          delete []cyl_data;
      }
  }
  return body_ptr;
}

void OCCModifyEngine::set_cylinder_eval_data
(
    ChollaEngine *cholla_ptr,
    double height,
    double base_radius_xdir,
    double base_radius_ydir,
    double top_radius,
    CubitEvaluatorData **&eval_data ) const
{
    DLIList<ChollaSurface *> cholla_surface_list;
    cholla_ptr->get_surfaces( cholla_surface_list );
    cholla_surface_list.reset();

    eval_data = NULL;

    if ( cholla_surface_list.size() != 3 )
    {
        // This cylinder/cone is shaped quite unusually such that the cholla engine
        // found more than 3 faces on it.  This is likely because it is quite smashed
        // in the minor radius direction.  As such, there is not a cylindrical
        // surface so exit.

        return;
    }

    eval_data = new CubitEvaluatorData* [3];
    eval_data[0] =
    eval_data[1] =
    eval_data[2] = NULL;

    int isurf;
    for ( isurf = 0; isurf < cholla_surface_list.size(); isurf++ )
    {
        eval_data[isurf] = NULL;
        DLIList<ChollaCurve*> cholla_curve_list;
        ChollaSurface *csurf = cholla_surface_list.get_and_step();
        csurf->get_curves( cholla_curve_list );
        if ( cholla_curve_list.size() == 2 )
        {
            // this is the cylindrical face around the cylinder.
            CylinderEvaluatorData *data = new CylinderEvaluatorData;
            data->height = height;
            data->height = height;
            data->base_radius_x = base_radius_xdir;
            data->base_radius_y = base_radius_ydir;
            data->top_radius = top_radius;

            eval_data[isurf] = data;
        }
    }
}

//===============================================================================
// Function   : torus
// Member Type: PUBLIC
// Description: create a torus with facets
// Author     : John Fowler
// Date       : 10/02
//===============================================================================
BodySM* OCCModifyEngine::torus( double r1, double r2 ) const
{
  CubitStatus rv = CUBIT_SUCCESS;
  DLIList <CubitFacet *>facet_list;
  DLIList <CubitPoint *>point_list;
  CubitPoint *new_point;
  CubitFacet *facet_ptr;
  int numpoints;
  double feature_angle;
int interp_order;
CubitBoolean smooth_non_manifold, split_surfaces;
  BodySM *body_ptr = NULL;
  STD(vector)<CubitPoint *> points;
  int reslevel;  // relative level of resolution
  double theta, thetainc, phi, phiinc, x, z, xp, yp, zp, rmajor, rminor;
  int numtheta, numphi, i, j;
 
  reslevel = 4;
  numtheta = 8*reslevel;
  numphi = 8*reslevel;  
  numpoints = numtheta*numphi;
  rmajor = r1;
  rminor = r2;  
  thetainc = 2.*CUBIT_PI/(double)numtheta;
  phiinc = 2.*CUBIT_PI/(double)numphi;
  phi = 0.;

//  Make the points in the y=0 plane
  for ( j = 0; j < numphi; j++ ) {
    theta = 0.;
    for ( i = 0; i < numtheta; i++ ) {
      x = rmajor + rminor*cos(theta);
      z = rminor*sin(theta);
//  Rotate around the z axis
      xp = x*cos(phi);
      zp = z;
      yp = x*sin(phi);
      new_point = (CubitPoint *) new CubitPointData( xp,yp,zp );
      points.push_back(new_point);
      theta += thetainc;
    }
    phi += phiinc;	
  }

  for ( i = 0; i < numpoints; i++ ) {
    point_list.append(points[i]);
  } 
//  Make the triangles  
  int m, k, m2, numtris;
  m = numtheta;
  numtris = 0;
  for ( j = 0; j < numphi; j++ ) {
    if ( j == numphi-1 ) m2 = 0;
    else m2 = m;
    for ( i = 0; i < numtheta; i++ ) {
      k = (i+1)%numtheta;
      facet_ptr = new CubitFacetData( points[i+m-numtheta], points[m2+i], points[m2+k] );
      facet_list.append( facet_ptr );     
      facet_ptr = new CubitFacetData( points[i+m-numtheta], points[m2+k], points[m-numtheta+k] );
      facet_list.append( facet_ptr );
      numtris += 2; 	
    }
    m += numtheta;
  }

  points.clear(); //  clear out the points vector since we are through with it.

  feature_angle = 135.0;
  interp_order = 0;
  smooth_non_manifold = CUBIT_TRUE;
  split_surfaces = CUBIT_FALSE;
  
  ChollaEngine *cholla_ptr = NULL;
  OCCModifyEngine *fme = const_cast<OCCModifyEngine *> (this);
  rv = fme->build_cholla_surfaces( facet_list,
                                   point_list,
                                   feature_angle,
                                   interp_order,
                                   smooth_non_manifold,
                                   split_surfaces,
                                   cholla_ptr );


  if ( rv == CUBIT_SUCCESS )
  {
      finish_facet_Body( cholla_ptr,
                         NULL,
                         feature_angle,
                         interp_order,
                         body_ptr);
      if ( cholla_ptr )
      {
         cholla_ptr->delete_me();
         delete cholla_ptr;
      }
  }
  return body_ptr;
}

//===============================================================================
// Function   : planar_sheet
// Member Type: PUBLIC
// Description: create a planar_sheet with facets
// Author     : John Fowler
// Date       : 10/02
//===============================================================================
BodySM* OCCModifyEngine::planar_sheet ( const CubitVector& /*p1*/,
                                       const CubitVector& /*p2*/,
                                       const CubitVector& /*p3*/,
                                       const CubitVector& /*p4*/ ) const
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return (BodySM*) NULL;
}

//===============================================================================
// Function   : copy_body
// Member Type: PUBLIC
// Description: copy a facet-based body
// Author     : John Fowler
// Date       : 10/02
//===============================================================================
BodySM* OCCModifyEngine::copy_body ( BodySM* /*bodyPtr*/) const
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return (BodySM*) NULL;
}

//===============================================================================
// Function   : subtract
// Member Type: PUBLIC
// Description: subtract boolean operation on facet-based bodies
// Author     : John Fowler
// Date       : 10/02
//===============================================================================
CubitStatus     OCCModifyEngine::subtract(DLIList<BodySM*> &tool_body_list,
                                            DLIList<BodySM*> &from_bodies,
                                            DLIList<BodySM*> &new_bodies,
                                            bool keep_old) const
{

  CubitStatus status = CUBIT_FAILURE;
  int i;
  BodySM *tool_body, *from_body;
  FacetboolInterface *fbint;
  CubitFacetboolOp op;

  bool *to_be_deleted = new bool[from_bodies.size()];

  op = CUBIT_FB_SUBTRACTION;

  from_bodies.reset();
  tool_body_list.reset();

  for ( i = 0; i < from_bodies.size(); i++ ) to_be_deleted[i] = false;

  for ( i = tool_body_list.size(); i > 0; i-- ) { 
    tool_body = tool_body_list.get_and_step(); 
    fbint = new FacetboolInterface;
    status = fbint->dofacetboolean_subtract(tool_body,from_bodies,new_bodies,
                                            keep_old,to_be_deleted,op);
    delete fbint;
    OCCQueryEngine::instance()->delete_solid_model_entities(tool_body);    
  }

  for ( i = 0; i < from_bodies.size(); i++ ) {
    from_body = from_bodies.get_and_step();
    if ( to_be_deleted[i] == true ) 
      OCCQueryEngine::instance()->delete_solid_model_entities(from_body);  
  }

  delete [] to_be_deleted;
    
  return status; 
     
}

//===============================================================================
// Function   : imprint
// Member Type: PUBLIC
// Description: imprint boolean operation on facet-based bodies
// Author     : John Fowler
// Date       : 10/02
//===============================================================================
CubitStatus     OCCModifyEngine::imprint(BodySM* /*BodyPtr1*/, BodySM* /*BodyPtr2*/,
                                           BodySM*& /*newBody1*/, BodySM*& /*newBody2*/,
                                           bool  /*keep_old*/) const
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return CUBIT_FAILURE;
}

//===============================================================================
// Function   : imprint
// Member Type: PUBLIC
// Description: imprint boolean operation on facet-based bodies
// Author     : John Fowler
// Date       : 10/02
//===============================================================================
CubitStatus     OCCModifyEngine::imprint(DLIList<BodySM*> &from_body_list ,
                                           DLIList<BodySM*> &new_from_body_list,
                                           bool keep_old) const
{
  CubitStatus success = CUBIT_SUCCESS;

    // total number of imprints to be done
  const int num_bodies = from_body_list.size();
//  int total_imprints = (num_bodies *(num_bodies-1))/2;
  int i, j;
  CubitBox bbox1, bbox2; 
  STD(vector)<BodySM*> bodies_vector;
  for ( i = 0; i < from_body_list.size(); i++ ) 
    bodies_vector.push_back(from_body_list.get_and_step());
    
    // pass in keep_old to the delete_owner_attrib flag; if we're not keeping
    // old bodies, we want to keep the owner attrib, so we can pick up entities
    // that didn't change
//  bool delete_attribs =
      (GeometryModifyTool::instance()->get_new_ids() || keep_old);

  from_body_list.reset();
  for (i = 0; success && i < num_bodies - 1; i++) {
    for (j = i + 1; j < num_bodies; j++) {
      BodySM *Body_Ptr1 = bodies_vector[i];
      success = OCCQueryEngine::instance()->create_facet_bounding_box(
                                             Body_Ptr1,
                                             bbox1);
      
      if (CubitMessage::instance()->Interrupt())
      {
        success = CUBIT_FAILURE;
        break;
      }    
      BodySM *Body_Ptr2 = bodies_vector[j];
      success = OCCQueryEngine::instance()->create_facet_bounding_box(
                                               Body_Ptr2,
                                               bbox2);
      if (bbox1.overlap(GEOMETRY_RESABS, bbox2)) {
        FacetboolInterface *FBInt = new FacetboolInterface;
        BodySM *out_Body_Ptr1, *out_Body_Ptr2;
        FBInt->dofacetboolean_2bodies_imprint(Body_Ptr1,Body_Ptr2,
                                              out_Body_Ptr1,out_Body_Ptr2,
                                              keep_old);
        delete FBInt;
        instance()->get_gqe()->delete_solid_model_entities(bodies_vector[i]);
        instance()->get_gqe()->delete_solid_model_entities(bodies_vector[j]);
        bodies_vector[i] = out_Body_Ptr1;
        bodies_vector[j] = out_Body_Ptr2;                                         
      }

    }
  }

  from_body_list.reset();
  for ( i = 0; i < from_body_list.size(); i++ ) {  
      new_from_body_list.append(bodies_vector[i]);
  }
   
  return success;
}

//===============================================================================
// Function   : imprint
// Member Type: PUBLIC
// Description: imprint boolean operation on facet-based bodies
// Author     : John Fowler
// Date       : 10/02
//===============================================================================
CubitStatus     OCCModifyEngine::imprint( DLIList<BodySM*> &body_list,
                                           DLIList<Curve*> &ref_edge_list,
                                           DLIList<BodySM*>& new_body_list,
                                           bool keep_old,
                                           bool show_messages) const
{
  CubitStatus success = CUBIT_SUCCESS;
//  bool delete_attribs =
      (GeometryModifyTool::instance()->get_new_ids() || keep_old);
  int i;
  CubitBox edge_list_bbox, bbox2; 
  const int num_bodies = body_list.size();

  FacetboolInterface *FBInt = new FacetboolInterface;
      
  FBInt->make_FB_edge_list(ref_edge_list);
  FBInt->get_edge_list_bbox(edge_list_bbox);

  for (i = 0; success && i < num_bodies; i++) {
    BodySM *Body_Ptr = body_list.next(i);
    success = OCCQueryEngine::instance()->create_facet_bounding_box(
                                             Body_Ptr,
                                             bbox2);
    if (edge_list_bbox.overlap(GEOMETRY_RESABS, bbox2)) {
    BodySM *new_body;
      FBInt->FB_imprint_with_curves(Body_Ptr,new_body,keep_old);
      new_body_list.append(new_body);      
    }
  }

  delete FBInt;

  body_list.reset();
  for ( i = 0; i < body_list.size(); i++ ) {  
      instance()->get_gqe()->delete_solid_model_entities(body_list.get_and_step());
  }

  return success;
}

//===============================================================================
// Function   : imprint
// Member Type: PUBLIC
// Description: imprint boolean operation on facet-based bodies
// Author     : John Fowler
// Date       : 10/02
//===============================================================================
CubitStatus     OCCModifyEngine::imprint( DLIList<Surface*> &/*ref_face_list*/,
                                           DLIList<Curve*> &/*ref_edge_list*/,
                                           DLIList<BodySM*>& /*new_body_list*/,
                                           bool /*keep_old_body*/ ) const
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return CUBIT_FAILURE;
}

//===============================================================================
// Function   : imprint
// Member Type: PUBLIC
// Description: imprint boolean operation on facet-based bodies
// Author     : John Fowler
// Date       : 10/02
//===============================================================================
CubitStatus     OCCModifyEngine::imprint( DLIList<Surface*> &/*surface_list*/,
                                           DLIList<DLIList<Curve*>*> &/*curve_lists_list*/,
                                           BodySM*& /*new_body*/,
                                           bool /*keep_old_body*/ ) const
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return CUBIT_FAILURE;
}

//===============================================================================
// Function   : imprint
// Member Type: PUBLIC
// Description: imprint boolean operation on facet-based bodies
// Author     : John Fowler
// Date       : 10/02
//===============================================================================
CubitStatus     OCCModifyEngine::imprint( DLIList<BodySM*> &/*body_list*/,
                                           DLIList<CubitVector*> &/*vector_list*/,
                                           DLIList<BodySM*>& /*new_body_list*/,
                                           bool keep_old /*keep_old_body*/ ) const
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return CUBIT_FAILURE;
}

//===============================================================================
// Function   : imprint_projected_edges
// Member Type: PUBLIC
// Description: 
// Author     : John Fowler
// Date       : 10/02
//===============================================================================
CubitStatus     OCCModifyEngine::imprint_projected_edges( DLIList<Surface*> &/*ref_face_list*/,
                                                           DLIList<Curve*> &/*ref_edge_list*/,
                                                           DLIList<BodySM*>& /*new_body_list*/,
                                                           bool /*keep_old_body*/,
                                                           bool /*keep_free_edges*/) const
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return CUBIT_FAILURE;
}

//===============================================================================
// Function   : imprint_projected_edges
// Member Type: PUBLIC
// Description: 
// Author     : John Fowler
// Date       : 10/02
//===============================================================================
CubitStatus     OCCModifyEngine::imprint_projected_edges(DLIList<Surface*> &/*ref_face_list*/,
                                                           DLIList<BodySM*> &/*body_list*/,
                                                           DLIList<Curve*> &/*ref_edge_list*/,
                                                           DLIList<BodySM*>& /*new_body_list*/,
                                                           bool /*keep_old_body*/,
                                                           bool /*keep_free_edges*/) const
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return CUBIT_FAILURE;
}

//===============================================================================
// Function   : project_edges
// Member Type: PUBLIC
// Description: 
// Author     : John Fowler
// Date       : 10/02
//===============================================================================
CubitStatus     OCCModifyEngine::project_edges( DLIList<Surface*> &/*ref_face_list*/,
                                                 DLIList<Curve*> &/*ref_edge_list_in*/,
                                                 DLIList<Curve*> &/*ref_edge_list_new*/,
                                                 bool /*print_error*/ ) const
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return CUBIT_FAILURE;
}


//===============================================================================
// Function   : intersect
// Member Type: PUBLIC
// Description: intersect boolean operation between facet-based bodies
// Author     : John Fowler
// Date       : 10/02
//===============================================================================
CubitStatus     OCCModifyEngine::intersect(BodySM*  tool_body_ptr,
                                             DLIList<BodySM*>  &from_bodies,
                                             DLIList<BodySM*>  &new_bodies,
                                             bool  keep_old) const
{

  CubitStatus status = CUBIT_FAILURE;
  int i;
  BodySM *from_body, *newBody;
  FacetboolInterface *fbint;
  CubitFacetboolOp op;
  DLIList<BodySM*> bodies;

  op = CUBIT_FB_INTERSECTION;
  fbint = new FacetboolInterface;

  from_bodies.reset();
  for ( i = from_bodies.size(); i > 0; i-- ) {
    bodies.clean_out();
    bodies.append(tool_body_ptr);    
    from_body = from_bodies.get_and_step();
    bodies.append(from_body);
    fbint = new FacetboolInterface;
    status = fbint->dofacetboolean(bodies,newBody,keep_old,op);
    delete fbint;     
    if ( status == CUBIT_SUCCESS && newBody) new_bodies.append(newBody);
  }
    
  return status;
  
/*
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return CUBIT_FAILURE;
*/  
}

//===============================================================================
// Function   : chop
// Member Type: PUBLIC
// Description: chop boolean operation between facet-based bodies
// Author     : John Fowler
// Date       : 10/02
//===============================================================================
CubitStatus      OCCModifyEngine::chop(DLIList<BodySM*>& bodies, 
                                         DLIList<BodySM*> &intersectBodies, 
                                         DLIList<BodySM*> &outsideBodies,
                                         BodySM*& leftoversBody,
                                         bool keep_old ,
                                         bool nonreg) const
{
//Note:  there is never any leftover body.
FacetboolInterface *fbint;
CubitFacetboolOp op;
bool keep_old_this = keep_old;
CubitStatus status;
BodySM *body_1, *body_2;
bool intersection_found;
BodySM* intersectBody = NULL;
BodySM* outsideBody = NULL;

  if ( bodies.size() > 2 ) {
    PRINT_ERROR("Chop not yet supported for more than two bodies.\n");
    return CUBIT_FAILURE;
  }
  body_1 = bodies.get_and_step();
  body_2 = bodies.get();
  op = CUBIT_FB_INTERSECTION;
  fbint = new FacetboolInterface;  
  status = fbint->dofacetboolean_2bodies(body_2,body_1,intersectBody,
                                         keep_old_this,intersection_found,op);
  intersectBodies.append( intersectBody );
 
  delete fbint;                                     
  op = CUBIT_FB_SUBTRACTION;
  fbint = new FacetboolInterface;  
  status = fbint->dofacetboolean_2bodies(body_1,body_2,outsideBody,
                                       keep_old_this,intersection_found,op);
  outsideBodies.append( outsideBody );

  delete fbint;                                       

  if ( keep_old == false ) {
      OCCQueryEngine::instance()->delete_solid_model_entities(body_1);
      OCCQueryEngine::instance()->delete_solid_model_entities(body_2);        
  }

  return status;
  
}

//===============================================================================
// Function   : unite
// Member Type: PUBLIC
// Description: unite boolean operation between facet-based bodies
// Author     : John Fowler
// Date       : 10/02
//===============================================================================
CubitStatus     OCCModifyEngine::unite(DLIList<BodySM*> &bodies, 
                                         DLIList<BodySM*> &newBodies,
                                         bool keep_old) const
{

CubitStatus status;
FacetboolInterface *fbint;
CubitFacetboolOp op;

  op = CUBIT_FB_UNION;
  fbint = new FacetboolInterface;

  BodySM *newBody = NULL;
  status = fbint->dofacetboolean(bodies,newBody,keep_old,op);

  newBodies.append( newBody );
  delete fbint;
  return status;
/*
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return CUBIT_FAILURE;
*/
  
}


//===============================================================================
// Function   : thicken
// Member Type: PUBLIC
// Description: 
// Author     : John Fowler
// Date       : 10/02
//===============================================================================
CubitStatus OCCModifyEngine::thicken(DLIList<BodySM*>& /*bodies*/, 
                                       DLIList<BodySM*>& /*new_bodies*/,
                                       double /*depth*/,
                                       bool /*both*/) const
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return CUBIT_FAILURE;
}


//===============================================================================
// Function   : flip_normals
// Member Type: PUBLIC
// Description: 
// Author     : John Fowler
// Date       : 10/02
//===============================================================================
CubitStatus OCCModifyEngine :: flip_normals( DLIList<Surface*>& face_list ) const
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return CUBIT_FAILURE;
}


//===============================================================================
// Function   : sweep_translational
// Member Type: PUBLIC
// Description: 
// Author     : John Fowler
// Date       : 10/02
//===============================================================================
CubitStatus OCCModifyEngine:: sweep_translational(
  DLIList<GeometryEntity*>& /*ref_ent_list*/,
  DLIList<BodySM*>& /*result_body_list*/,
  const CubitVector& /*sweep_vector*/,
  double /*draft_angle*/,
  int /*draft_type*/,
  bool /*switchside*/,
  bool /*rigid*/) const
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return CUBIT_FAILURE;
}

//===============================================================================
// Function   : sweep_perpendicular
// Member Type: PUBLIC
// Description: 
// Author     : John Fowler
// Date       : 10/02
//===============================================================================
CubitStatus OCCModifyEngine:: sweep_perpendicular(
  DLIList<GeometryEntity*>& /*ref_ent_list*/,
  DLIList<BodySM*>& /*result_body_list*/,
  double /*distance*/,
  double /*draft_angle*/,
  int /*draft_type*/,
  bool /*switchside*/,
  bool /*rigid*/) const
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return CUBIT_FAILURE;
}

//===============================================================================
// Function   : sweep_rotational
// Member Type: PUBLIC
// Description: 
// Author     : John Fowler
// Date       : 10/02
//===============================================================================
CubitStatus OCCModifyEngine:: sweep_rotational(
  DLIList<GeometryEntity*>& /*ref_ent_list*/,
  DLIList<BodySM*>& /*result_body_list*/,
  const CubitVector& /*point*/,
  const CubitVector& /*direction*/,
  double /*angle*/,
  int /*steps*/,
  double /*draft_angle*/,
  int /*draft_type*/,
  bool /*switchside*/,
  bool /*make_solid*/,
  bool /*rigid*/) const
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return CUBIT_FAILURE;
}

//===============================================================================
// Function   : sweep_along_curve
// Member Type: PUBLIC
// Description: 
// Author     : John Fowler
// Date       : 10/02
//===============================================================================
CubitStatus OCCModifyEngine::sweep_along_curve(
  DLIList<GeometryEntity*>& /*ref_ent_list*/,
  DLIList<BodySM*>& /*result_body_list*/,
  DLIList<Curve*>& /*ref_edge_list*/,
  double /*draft_angle*/,
  int /*draft_type*/,
  bool /*rigid*/) const
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return CUBIT_FAILURE;
}

//HEADER- Webcut-related functions

//===============================================================================
// Function   : webcut
// Member Type: PUBLIC
// Description: 
// Author     : John Fowler
// Date       : 10/02
//===============================================================================
CubitStatus OCCModifyEngine::webcut(DLIList<BodySM*>& webcut_body_list,
                              const CubitVector &v1,
                              const CubitVector &v2,
                              const CubitVector &v3,
                              DLIList<BodySM*>& results_list,
                              bool imprint ) const
{

CubitBox super_bbox;
CubitStatus status;
int i;

  CubitBoolean delete_bodies = (GeometryModifyTool::instance()->get_new_ids() ?
                                CUBIT_FALSE : CUBIT_TRUE);

//  Find the bounding box of all of the bodies.  This will be used to make
//  the cutting plane.


  status = OCCQueryEngine::instance()->create_super_facet_bounding_box(webcut_body_list,
                                             super_bbox);
  
//  Find the size of the cutting plane (x,y) in terms of super_bbox.
  DLIList<CubitVector> intersection_points;
  FBDataUtil::intersect_plane_with_boundingbox(super_bbox,
                                  v1,v2,v3,intersection_points); 
  int numpts = intersection_points.size();;
  if ( numpts < 3 ) {
      PRINT_INFO("INFO: Cutting Tool overlaps the original volume,\n"
                 "      Or cutting plane does not pass through volume.\n"
                 "         The original volume is unaffected.\n" );

    return CUBIT_FAILURE;
  }                                            
double xsize, ysize, xcen, ycen, zcen;
double xmin, ymin, zmin, xmax, ymax, zmax, xx, yy, zz;
  xmin = ymin = zmin = CUBIT_DBL_MAX;
  xmax = ymax = zmax = -xmin + 1;
  xcen = ycen = zcen = 0.0;
  for ( i = 0; i < numpts; i++ ) {
    xx = intersection_points[i].x();
    yy = intersection_points[i].y();
    zz = intersection_points[i].z();
    
    xcen += xx;
    ycen += yy;
    zcen += zz;

    xmin = (xmin < xx ) ? xmin : xx; 
    ymin = (ymin < yy ) ? ymin : yy; 
    zmin = (zmin < zz ) ? zmin : zz; 
    xmax = (xmax > xx ) ? xmax : xx; 
    ymax = (ymax > yy ) ? ymax : yy; 
    zmax = (zmax > zz ) ? zmax : zz;        
  }
  xcen /= numpts; ycen /= numpts; zcen /= numpts;  
  
  //  Could do this better by rotating the intersection points into the
  //  x-y plane and then getting xsize and ysize.  Factor of 1.3 is just
  //  to make sure that the plane extends beyond the bodies.
  xsize = ysize = 1.3*sqrt( (xmax-xmin)*(xmax-xmin) + 
                            (ymax-ymin)*(ymax-ymin) +
                            (zmax-zmin)*(zmax-zmin) );
  STD(vector)<double> cutter_verts;
  STD(vector)<int> cutter_connections;
  int numx, numy;
  numx = 20;
  numy = 20;
  //  Make the cutter surface.
  status = FBDataUtil::FBmake_xy_plane(cutter_verts, cutter_connections, 
                  xsize, ysize, numx, numy);
  CubitVector va, vb;
  va = v1 - v2;
  vb = v3 - v2;               
  CubitVector rotate_to;
 
  rotate_to = vb*va;
  rotate_to.normalize();
  CubitVector center_pt(xcen,ycen,zcen);
  status = FBDataUtil::rotate_FB_object(cutter_verts,rotate_to,center_pt);
  FacetboolInterface *fbint;
  bool cutter_is_plane = true;
  //  Now make the facetbool objects for the bodies. 
  webcut_body_list.reset();                
BodySM *body_sm;

  for ( i = webcut_body_list.size(); i > 0; i-- ) {
    CubitBoolean intersects;
    body_sm = webcut_body_list.get_and_step();
    fbint = new FacetboolInterface;
    status = fbint->webcut_FB(body_sm,cutter_verts,cutter_connections,
                              cutter_is_plane,delete_bodies,
                              intersects,results_list); 
    delete fbint;    
    if ( status == CUBIT_FAILURE )
    {
        PRINT_ERROR(" Unable to perform webcut.\n" );
        return CUBIT_FAILURE;
    }

    if ( status == CUBIT_SUCCESS && intersects )  
    {
      instance()->get_gqe()->delete_solid_model_entities(body_sm);
    }
    else
    {
      PRINT_INFO("INFO: Cutting Tool overlaps the original volume,\n"
                 "      Or cutting plane does not pass through volume.\n"
                 "         The original volume is unaffected.\n" );
    }
  }
                                 
  return status;
    
}

//===============================================================================
// Function   : webcut
// Member Type: PUBLIC
// Description: 
// Author     : John Fowler
// Date       : 10/02
//===============================================================================
CubitStatus    OCCModifyEngine::webcut(DLIList<BodySM*>& /*webcut_body_list*/,
                                 BodySM const* /*tool_body*/,
                                 DLIList<BodySM*>& /*results_list*/,
                                 bool /*imprint*/ ) const
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return CUBIT_FAILURE;
}

//===============================================================================
// Function   : webcut_across_translate
// Member Type: PUBLIC
// Description: 
// Author     : John Fowler
// Date       : 10/02
//===============================================================================
CubitStatus    OCCModifyEngine::webcut_across_translate( DLIList<BodySM*>& /*body_list*/,
                                                          Surface* /*plane_surf1*/,
                                                          Surface* /*plane_surf2*/,
                                                          DLIList<BodySM*>& /*results_list*/,
                                                          bool /*imprint*/ ) const
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return CUBIT_FAILURE;
}

//===============================================================================
// Function   : webcut_with_sheet
// Member Type: PUBLIC
// Description: 
// Author     : John Fowler
// Date       : 10/02
//===============================================================================
CubitStatus OCCModifyEngine::webcut_with_sheet(DLIList<BodySM*> & /*webcut_body_list*/,
                                                 BodySM * /*sheet_body*/,
                                                 DLIList<BodySM*> & /*new_bodies*/,
                                                 bool /*imprint*/ )
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return CUBIT_FAILURE;
}

//===============================================================================
// Function   : webcut_with_extended_surf
// Member Type: PUBLIC
// Description: 
// Author     : John Fowler
// Date       : 10/02
//===============================================================================
CubitStatus OCCModifyEngine::webcut_with_extended_surf(DLIList<BodySM*> & /*webcut_body_list*/,
                                                         Surface * /*extend_from*/,
                                                         DLIList<BodySM*> & /*new_bodies*/,
                                                         int & /*num_cut*/,
                                                         bool /*imprint*/ )
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return CUBIT_FAILURE;
}

//===============================================================================
// Function   : webcut_with_cylinder
// Member Type: PUBLIC
// Description: 
// Author     : John Fowler
// Date       : 10/02
//===============================================================================
CubitStatus OCCModifyEngine::webcut_with_cylinder(DLIList<BodySM*> &webcut_body_list,
                                            double radius,
                                            const CubitVector &axis,
                                            const CubitVector &center,
                                            DLIList<BodySM*>& results_list,
                                            bool imprint )
{

CubitBox super_bbox;
CubitStatus status;
int i;
CubitVector bodies_center, my_center, diagonal, my_axis;

  CubitBoolean delete_bodies = (GeometryModifyTool::instance()->get_new_ids() ?
                                CUBIT_FALSE : CUBIT_TRUE);

  status = OCCQueryEngine::instance()->create_super_facet_bounding_box(
    webcut_body_list,super_bbox);
  STD(vector)<double> cutter_verts;
  STD(vector)<int> cutter_connections;
  int nr, nz;
  double length;
  
  diagonal = super_bbox.diagonal();
//  length = 2.3*diagonal.length() + 2.0*center.length();
  bodies_center = super_bbox.center();
  length = 3.*sqrt((bodies_center.x() - center.x())*(bodies_center.x() - center.x()) +
                   (bodies_center.y() - center.y())*(bodies_center.y() - center.y()) +
                   (bodies_center.z() - center.z())*(bodies_center.z() - center.z()) );
  length += 3.*diagonal.length();
    //length = sqrt(length*length + radius*radius);
    //  bodies_center += center;

  nr = 30;
  nz = 5;

  //  Make the cutter surface.
  status = FBDataUtil::FBmake_cylinder(cutter_verts, cutter_connections, 
                  radius, length, nr, nz);
  my_center = center;
  my_axis = axis;
  status = FBDataUtil::rotate_FB_object(cutter_verts,my_axis,my_center);

  FacetboolInterface *fbint;
  bool cutter_is_plane = false;
  //  Now make the facetbool objects for the bodies.                  
  webcut_body_list.reset();                
  BodySM* body_sm;
  for ( i = webcut_body_list.size(); i > 0; i-- ) {
    CubitBoolean intersects;
    body_sm = webcut_body_list.get_and_step();
    fbint = new FacetboolInterface;
     status = fbint->webcut_FB(body_sm,cutter_verts,cutter_connections,
                               cutter_is_plane,delete_bodies,
                               intersects,results_list);   
    delete fbint;    
    if ( status == CUBIT_FAILURE )
    {
        PRINT_ERROR(" Unable to perform webcut.\n" );
        return CUBIT_FAILURE;
    }

    if ( status == CUBIT_SUCCESS && intersects )  
    {
      instance()->get_gqe()->delete_solid_model_entities(body_sm);
    }
    else
    {
      PRINT_INFO("INFO: Cutting Tool overlaps the original volume,\n"
                 "      Or cutting plane does not pass through volume.\n"
                 "         The original volume is unaffected.\n" );
    }
  }

  return status;
  
}

//===============================================================================
// Function   : webcut_with_brick
// Member Type: PUBLIC
// Description: 
// Author     : John Fowler
// Date       : 10/02
//===============================================================================
CubitStatus OCCModifyEngine::webcut_with_brick( 
                                      DLIList<BodySM*>& /*webcut_body_list*/, 
                                      const CubitVector &/*center*/,
                                      const CubitVector* /*axes[3]*/, 
                                      const CubitVector &/*extension*/,
                                      DLIList<BodySM*> &/*results_list*/,
                                      bool /*imprint*/ )
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return CUBIT_FAILURE;
}

//===============================================================================
// Function   : webcut_with_planar_sheet
// Member Type: PUBLIC
// Description: 
// Author     : John Fowler
// Date       : 10/02
//===============================================================================
CubitStatus OCCModifyEngine::webcut_with_planar_sheet( 
                                          DLIList<BodySM*>& /*webcut_body_list*/,
                                          const CubitVector &/*center*/,
                                          const CubitVector* /*axes[2]*/,
                                          double /*width*/, 
                                          double /*height*/,
                                          DLIList<BodySM*> &/*results_list*/,
                                          bool /*imprint*/ )
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return CUBIT_FAILURE;
}

//===============================================================================
// Function   : webcut_with_curve_loop
// Member Type: PUBLIC
// Description: 
// Author     : John Fowler
// Date       : 10/02
//===============================================================================
CubitStatus OCCModifyEngine::webcut_with_curve_loop(
                                              DLIList<BodySM*> &/*webcut_body_list*/,
                                              DLIList<Curve*> &/*ref_edge_list*/,
                                              DLIList<BodySM*>& /*results_list*/,
                                              bool /*imprint*/)
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return CUBIT_FAILURE;
}

//===============================================================================
// Function   : section
// Member Type: PUBLIC
// Description: 
// Author     : John Fowler
// Date       : 10/02
//===============================================================================
CubitStatus OCCModifyEngine::section( DLIList<BodySM*> &/*section_body_list*/,
                                        const CubitVector &/*point_1*/,
                                        const CubitVector &/*point_2*/,
                                        const CubitVector &/*point_3*/,
                                        DLIList<BodySM*>& /*new_body_list*/,
                                        bool /*keep_normal_side*/,
                                        bool /*keep_old*/,
                                        bool /*keep_both_sides*/)
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return CUBIT_FAILURE;
}

//===============================================================================
// Function   : split_body
// Member Type: PUBLIC
// Description: Splits multiple lumps in one body into separate bodies
// Author     : Corey Ernst 
// Date       : 08/04
//===============================================================================
CubitStatus OCCModifyEngine::split_body( BodySM *body_ptr,
                                           DLIList<BodySM*> &new_bodies )
{
  //get the all lumps of input body
  DLIList<Lump*> lumps;
  body_ptr->lumps( lumps );

  if( lumps.size() == 1 )
  {
    new_bodies.append( body_ptr );
    return CUBIT_SUCCESS;
  }

  //for each lump except one first one, create a new body
  DLIList<Lump*> single_lump;
  lumps.reset();
  lumps.step();
  int i;
  OCCBody *tmp_facet_body = static_cast<OCCBody*>( body_ptr );
  for( i=lumps.size()-1; i--; )
  {
    BodySM *bodysm_ptr;
    single_lump.clean_out();
    tmp_facet_body->remove_lump( static_cast<OCCLump*>(lumps.get())); 
    single_lump.append( lumps.get_and_step() );
    make_facet_body(single_lump, bodysm_ptr);
    if( bodysm_ptr )
      new_bodies.append(bodysm_ptr);
  }

  new_bodies.append( body_ptr );

  return CUBIT_SUCCESS;
}


//===============================================================================
// Function   : reverse_body
// Member Type: PUBLIC
// Description: Turn body inside-out
// Author     : Jason Kraftcheck
// Date       : 05/25/04
//===============================================================================
CubitStatus OCCModifyEngine::reverse_body( BodySM* body_ptr )
{
  OCCBody* body = dynamic_cast<OCCBody*>(body_ptr);
  if (0 == body)
  {
    PRINT_ERROR("Non-facet body in FME::reverse.\n");
    return CUBIT_FAILURE;
  }
  
    // Flip CoFace senses
  DLIList<OCCShell*> shells;
  body->get_shells( shells );
  while (shells.size())
    shells.pop()->reverse();
  
  return CUBIT_SUCCESS;
}
    


//===============================================================================
// Function   : split_periodic
// Member Type: PUBLIC
// Description: 
// Author     : John Fowler
// Date       : 10/02
//===============================================================================
CubitStatus OCCModifyEngine::split_periodic( BodySM * /*body_ptr*/,
                                               BodySM *& /*new_body*/ )
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return CUBIT_FAILURE;
}

//===============================================================================
// Function   : regularize_body
// Member Type: PUBLIC
// Description: 
// Author     : John Fowler
// Date       : 10/02
//===============================================================================
CubitStatus    OCCModifyEngine::regularize_body( BodySM * /*body_ptr*/,
                                                   BodySM *& /*new_body_ptr*/ )
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return CUBIT_FAILURE;
}

//===============================================================================
// Function   : regularize_refentity
// Member Type: PUBLIC
// Description: 
// Author     : John Fowler
// Date       : 10/02
//===============================================================================
CubitStatus  OCCModifyEngine::regularize_entity( GeometryEntity * /*old_entity_ptr*/,  
                                                      BodySM *& /*new_body_ptr*/ )
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return CUBIT_FAILURE;
}

//===============================================================================
// Function   : offset_curves
// Member Type: PUBLIC
// Description: 
// Author     : John Fowler
// Date       : 10/02
//===============================================================================
CubitStatus OCCModifyEngine::offset_curves( DLIList<Curve*>& /*ref_edge_list*/, 
                                              DLIList<Curve*>&,
                                              double /*offset_distance*/,
                                              const CubitVector& /*offset_direction*/, 
                                              int /*gap_type*/ )
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return CUBIT_FAILURE;
}

//===============================================================================
// Function   : trim_curve
// Member Type: PUBLIC
// Description: 
// Author     : John Fowler
// Date       : 10/02
//===============================================================================
Curve* OCCModifyEngine::trim_curve( Curve* /*trim_curve*/, 
                                      const CubitVector& /*trim_vector*/,
                                      const CubitVector& /*keep_vector*/,
                                      bool )
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return 0;
}

//===============================================================================
// Function   : create_body_from_surfs
// Member Type: PUBLIC
// Description: 
// Author     : Steve Owen
// Date       : 9/11/03
//===============================================================================
CubitStatus OCCModifyEngine::create_body_from_surfs(DLIList<Surface*> & ref_face_list,
                                                      BodySM *& new_body,
                                                      bool keep_old,
                                                      bool heal) const
{
  // get the facets from the faces

  int ii;
  Surface *surf_ptr;
  OCCSurface *fsurf_ptr;
  DLIList<CubitFacet *> facet_list;
  DLIList<CubitPoint *> point_list;
#ifdef BOYD17 
  DLIList<CubitFacetEdge *>edge_list;
#endif

  if (heal)
  {
    PRINT_WARNING("\"heal\" option for facet-based geometry is not supported.\n"); 
  }

  ref_face_list.reset();
  for(ii=0; ii<ref_face_list.size(); ii++)
  {
    surf_ptr = ref_face_list.get_and_step();
    fsurf_ptr = dynamic_cast<OCCSurface *>(surf_ptr);
    assert(fsurf_ptr != NULL);
    if (fsurf_ptr == NULL)
      return CUBIT_FAILURE;
    point_list.clean_out();
    fsurf_ptr->get_my_facets( facet_list, point_list );
  }

  // copy the facets

  DLIList<CubitFacet *> new_facet_list;
  DLIList<CubitPoint *> new_point_list;
  DLIList<CubitFacetEdge *> new_edge_list;
  FacetDataUtil::copy_facets( facet_list, new_facet_list, new_point_list, new_edge_list );

  // generate new geometry

  const char *file_name = NULL; 
  CubitBoolean use_feature_angle = CUBIT_FALSE; 
  double feature_angle = 0.0;
  int interp_order = 0;
  double tol=0.0;
  CubitBoolean smooth_non_manifold = CUBIT_FALSE;
  CubitBoolean split_surfaces = CUBIT_FALSE;
  CubitBoolean stitch = CUBIT_TRUE;
  CubitBoolean improve = CUBIT_TRUE;
  DLIList <CubitQuadFacet *> quad_facet_list;
  DLIList<Surface *> surface_list;
  FacetFileFormat file_format = FROM_FACET_LIST;
  CubitStatus rv = 
    OCCQueryEngine::instance()->import_facets( file_name, 
                                               use_feature_angle,
                                               feature_angle,
                                               tol,
                                               interp_order,
                                               smooth_non_manifold,
                                               split_surfaces,
                                               stitch,
                                               improve,
                                               quad_facet_list,
                                               new_facet_list,
                                               surface_list,
                                               file_format );

  if (rv == CUBIT_SUCCESS)
  {
    surf_ptr = surface_list.get();
    new_body = surf_ptr->bodysm();
  
    // delete the old model
  
    if (!keep_old)
    {
      DLIList<BodySM*> mybody_list;
      DLIList<BodySM*> body_list;
      for(ii=0; ii<ref_face_list.size(); ii++)
      {
        mybody_list.clean_out();
        surf_ptr = ref_face_list.get_and_step();
        surf_ptr->bodysms(mybody_list);
        if (mybody_list.size() == 1)
        {
          body_list.append_unique( mybody_list.get() );
        }
      }
      if (body_list.size() > 0)
      {
        OCCQueryEngine::instance()->delete_solid_model_entities(body_list);
        //GeometryQueryTool::instance()->delete_Body( body_list );
      }
    }
  }
  else
  {
    new_body = NULL;
  }

  return CUBIT_SUCCESS;
}

//===============================================================================
// Function   : create_arc_three
// Member Type: PUBLIC
// Description: 
// Author     : John Fowler
// Date       : 10/02
//===============================================================================
Curve* OCCModifyEngine::create_arc_three( Point* /*ref_vertex1*/, 
                                            Point* /*ref_vertex2*/,
                                            Point* /*ref_vertex3*/, 
                                            bool /*full*/ )
{ return NULL;

}

//===============================================================================
// Function   : create_arc_three
// Member Type: PUBLIC
// Description: 
// Author     : John Fowler
// Date       : 10/02
//===============================================================================
Curve* OCCModifyEngine::create_arc_three( Curve* /*ref_edge1*/, 
                                            Curve* /*ref_edge2*/,
                                            Curve* /*ref_edge3*/, 
                                            bool /*full*/  )
{ return NULL;
}

//===============================================================================
// Function   : create_arc_center_edge
// Member Type: PUBLIC
// Description: 
// Author     : John Fowler
// Date       : 10/02
//===============================================================================
Curve* OCCModifyEngine::create_arc_center_edge( Point* /*ref_vertex1*/, 
                                                  Point* /*ref_vertex2*/,
                                                  Point* /*ref_vertex3*/,
                                                  const CubitVector& /*normal*/, 
                                                  double /*radius*/,
                                                  bool /*full*/ ) 
{ 
  return NULL; 
}

CubitStatus 
OCCModifyEngine::create_curve_combine( DLIList<Curve*>& curve_list, 
                                    Curve *&new_curve_ptr )
{
  PRINT_ERROR("Curve combine is not implemented for facet based models\n");
  return CUBIT_FAILURE;
}

//===============================================================================
// Function   : get_gqe
// Member Type: PUBLIC
// Description: get the facet geometry query engince instance pointer
// Author     : John Fowler
// Date       : 10/02
//===============================================================================
GeometryQueryEngine *OCCModifyEngine::get_gqe()
{
  return OCCQueryEngine::instance();
}

//===============================================================================
// Function   : is_modify_engine
// Member Type: PUBLIC
// Description: return CUBIT_TRUE if the tb_ptr belongs to this modify engine
// Author     : John Fowler
// Date       : 10/02
//===============================================================================
CubitBoolean OCCModifyEngine::is_modify_engine(const TopologyBridge *tb_ptr) const 
{
  return tb_ptr->get_geometry_query_engine() == OCCQueryEngine::instance();
}

//===============================================================================
// Function   : get_offset_intersections
// Member Type: PUBLIC
// Description: 
// Author     : John Fowler
// Date       : 10/02
//===============================================================================
CubitStatus OCCModifyEngine::get_offset_intersections( Curve* /*ref_edge1*/, 
                                                         Curve* /*ref_edge2*/,
                                                         DLIList<CubitVector*>& /*intersection_list*/,
                                                         double /*offset*/,
                                                         CubitBoolean /*ext_first*/ )
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return CUBIT_FAILURE;
}

//===============================================================================
// Function   : get_offset_intersections
// Member Type: PUBLIC
// Description: 
// Author     : John Fowler
// Date       : 10/02
//===============================================================================
CubitStatus OCCModifyEngine::get_offset_intersections( Curve* /*ref_edge_ptr*/, 
                                                         Surface* /*ref_face_ptr*/,
                                                         DLIList<CubitVector*> & /*intersection_list*/,
                                                         double /*offset*/,
                                                         CubitBoolean /*ext_surf*/ )
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return CUBIT_FAILURE;
}

//===============================================================================
// Function   : surface_intersection
// Member Type: PUBLIC
// Description: 
// Author     : John Fowler
// Date       : 10/02
//===============================================================================
CubitStatus OCCModifyEngine::surface_intersection( Surface * /*surface1_ptr*/,
                                                     Surface * /*surface2_ptr*/,
                                                     DLIList<Curve*> &/*inter_graph*/,
                                                     const double /*tol*/) const
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return CUBIT_FAILURE;
}

//===============================================================================
// Function   : get_mid_plane
// Member Type: PUBLIC
// Description: 
// Author     : John Fowler
// Date       : 10/02
//===============================================================================
CubitStatus OCCModifyEngine::get_mid_plane( const CubitVector & /*point_1*/,
                                              const CubitVector & /*point_2*/,
                                              const CubitVector & /*point_3*/,
                                              BodySM * /*body_to_trim_to*/,
                                              BodySM *& /*midplane_body*/ ) const
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return CUBIT_FAILURE;
}

//=============================================================================
// Function   : tweak_chamfer
// Member Type: PUBLIC
// Description: Chamfer curves on solid bodies.  The left and right offsets are
//              with respect to the curve direction.  If the given right offset
//              is negative, the left offset is used.  Users can preview to
//              clarify the meaning of left and right.
// Author     : 
// Date       : 
//=============================================================================
CubitStatus OCCModifyEngine::tweak_chamfer( DLIList<Curve*> & /*curve_list*/, 
                                              double /*left_offset*/,
                                              DLIList<BodySM*> & /*new_bodysm_list*/,
                                              double /*right_offset*/,
                                              CubitBoolean /*keep_old_body*/,
                                              CubitBoolean /*preview*/ ) const
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return CUBIT_FAILURE;
}

//=============================================================================
// Function   : tweak_chamfer
// Member Type: PUBLIC
// Description: Chamfer vertices on solid or sheet bodies.  On a solid body 
//              there can be up to 3 offsets; on a sheet body up to 2 offsets.
//              The offsets are in the direction of the supplied edges.  If 
//              multiple vertices are supplied, only one offset value is 
//              allowed and the edges are not used.
// Author     : 
// Date       : 
//=============================================================================
CubitStatus
OCCModifyEngine::tweak_chamfer( DLIList<Point*> & /*point_list*/, 
                                  double /*offset1*/,
                                  DLIList<BodySM*> & /*new_bodysm_list*/,
                                  Curve * /*edge1*/,
                                  double /*offset2*/,
                                  Curve * /*edge2*/,
                                  double /*offset3*/,
                                  Curve * /*edge3*/,
                                  CubitBoolean /*keep_old_body*/,
                                  CubitBoolean /*preview*/ ) const
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return CUBIT_FAILURE;
}

//=============================================================================
// Function   : tweak_fillet
// Member Type: PUBLIC
// Description: Create a round fillet (or blend) at the given curves on solid 
//              bodies.
// Author     : 
// Date       : 
//=============================================================================
CubitStatus OCCModifyEngine::tweak_fillet( DLIList<Curve*> & /*curve_list*/, 
                                             double /*radius*/,
                                             DLIList<BodySM*> & /*new_bodysm_list*/,
                                             CubitBoolean /*keep_old_body*/,
                                             CubitBoolean /*preview*/ ) const
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return CUBIT_FAILURE;
}

//=============================================================================
// Function   : tweak_fillet
// Member Type: PUBLIC
// Description: Create a round fillet (or blend) at the given curves on a solid 
//              body.  The fillet has a variable radius from the start to the
//              end of the curve.
// Author     : 
// Date       : 
//=============================================================================
CubitStatus OCCModifyEngine::tweak_fillet( Curve * /*curve_ptr*/, 
                                             double /*start_radius*/,
                                             double /*end_radius*/,
                                             BodySM *& /*new_bodysm_ptr*/,
                                             CubitBoolean /*keep_old_body*/,
                                             CubitBoolean /*preview*/ ) const
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return CUBIT_FAILURE;
}

//=============================================================================
// Function   : tweak_fillet
// Member Type: PUBLIC
// Description: Create a round fillet (or blend) at the given vertices on sheet
//              bodies.
// Author     : 
// Date       : 
//=============================================================================
CubitStatus
OCCModifyEngine::tweak_fillet( DLIList<Point*> & /*ref_vertex_list*/, 
                                 double /*radius*/,
                                 DLIList<BodySM*> & /*new_bodysm_list*/,
                                 CubitBoolean /*keep_old_body*/,
                                 CubitBoolean /*preview*/ ) const
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return CUBIT_FAILURE;
}

//=============================================================================
// Function   : tweak_move
// Member Type: PUBLIC
// Description: Tweak specified faces of a volume or volumes along a vector.
// Author     : 
// Date       : 
//=============================================================================
CubitStatus OCCModifyEngine::tweak_move( DLIList<Surface*> & /*surface_list*/, 
                                           const CubitVector & /*delta*/,
                                           DLIList<BodySM*> & /*new_bodysm_list*/, 
                                           CubitBoolean /*keep_old_body*/ ) const
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return CUBIT_FAILURE;
}

//=============================================================================
// Function   : tweak_move
// Member Type: PUBLIC
// Description: Tweak specified curves of a sheet body along a vector.
// Author     : 
// Date       : 
//=============================================================================
CubitStatus OCCModifyEngine::tweak_move( DLIList<Curve*> & /*curve_list*/,
                                           const CubitVector & /*delta*/,
                                           DLIList<BodySM*> & /*new_bodysm_list*/, 
                                           CubitBoolean /*keep_old_body*/ ) const
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return CUBIT_FAILURE;
}

//=============================================================================
// Function   : tweak_offset
// Member Type: PUBLIC
// Description: Tweak specified faces of a volume or volumes by offsetting
//              those faces by the offset distance.
// Author     : 
// Date       : 
//=============================================================================
CubitStatus OCCModifyEngine::tweak_offset( DLIList<Surface*> & /*surface_list*/, 
                                             double /*offset_distance*/,
                                             DLIList<BodySM*> & /*new_bodysm_list*/,
                                             CubitBoolean /*keep_old_body*/ ) const
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return CUBIT_FAILURE;
}

//=============================================================================
// Function   : tweak_offset
// Member Type: PUBLIC
// Description: Tweak specified curves of a sheet body or bodies by offsetting
//              those curves by the offset distance.
// Author     : 
// Date       : 
//=============================================================================
CubitStatus OCCModifyEngine::tweak_offset( DLIList<Curve*> & /*curve_list*/,  
                                             double /*offset_distance*/,
                                             DLIList<BodySM*> & /*new_bodysm_list*/,
                                             CubitBoolean /*keep_old_body*/ ) const
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return CUBIT_FAILURE;
}

//=============================================================================
// Function   : tweak_remove
// Member Type: PUBLIC
// Description: Function to remove surfaces from a body and then extend the 
//              remaining surfaces to fill the gap or hole.
// Author     : 
// Date       : 
//=============================================================================
CubitStatus OCCModifyEngine::tweak_remove( DLIList<Surface*> & /*surface_list*/,
                                             DLIList<BodySM*> & /*new_bodysm_list*/,
                                             CubitBoolean /*extend_adjoining*/,
                                             CubitBoolean /*keep_surface*/,
                                             CubitBoolean /*keep_old_body*/ ) const
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return CUBIT_FAILURE;
}

//=============================================================================
// Function   : tweak_remove
// Member Type: PUBLIC
// Description: Function to remove curves from a sheet body and then extend the 
//              remaining curves or fill the gap or hole.
// Author     : 
// Date       : 
//=============================================================================
CubitStatus OCCModifyEngine::tweak_remove( DLIList<Curve*> & /*curve_list*/,
                                             DLIList<BodySM*> & /*new_bodysm_list*/, 
                                             CubitBoolean /*keep_old_body*/ ) const
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return CUBIT_FAILURE;
}

//=============================================================================
// Function   : tweak_target
// Member Type: PUBLIC
// Description: Tweak specified faces of a volume or volumes up to a target 
//              surface.
// Author     : 
// Date       : 
//=============================================================================
CubitStatus OCCModifyEngine::tweak_target( DLIList<Surface*> & /*surface_list*/,
                                             Surface * /*target_surf_ptr*/,
                                             DLIList<BodySM*> & /*new_bodysm_list*/,
                                             CubitBoolean /*reverse_flg*/,
                                             CubitBoolean /*keep_old_body*/ ) const
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return CUBIT_FAILURE;
}

//=============================================================================
// Function   : tweak_target
// Member Type: PUBLIC
// Description: Tweak specified edges of a surface or set of surfaces (in sheet
//              bodies) up to a target surface.
// Author     : 
// Date       : 
//=============================================================================
CubitStatus OCCModifyEngine::tweak_target( DLIList<Curve*> & /*curve_list*/,
                                             Surface * /*target_surf_ptr*/, 
                                             DLIList<BodySM*> & /*new_bodysm_list*/, 
                                             CubitBoolean /*keep_old_body*/ ) const
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return CUBIT_FAILURE;
}

//=============================================================================
// Function   : tweak_target
// Member Type: PUBLIC
// Description: Tweak specified edges of a sheet body or bodies up to a target
//              curve that is part of a sheet body.  The target is a surface 
//              created by thickening the owning surface of the target curve.
// Author     : 
// Date       : 
//=============================================================================
CubitStatus OCCModifyEngine::tweak_target( DLIList<Curve*> & /*curve_list*/,
                                             Curve * /*target_curve_ptr*/, 
                                             DLIList<BodySM*> & /*new_bodysm_list*/, 
                                             CubitBoolean /*keep_old_body*/ ) const
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return CUBIT_FAILURE;
}

//================================================================================
// Description: Creates a net surface.
// Author     : Tyronne Lim
// Date       : 08/18/03
//================================================================================
CubitStatus OCCModifyEngine::create_net_surface( DLIList<Surface*>& /*ref_face_list*/, 
                                                   BodySM *& /*new_body*/,
                                                   DLIList<DLIList<CubitVector*>*> & /*vec_lists_u*/, 
                                                   DLIList<DLIList<CubitVector*>*> & /*vec_lists_v*/, 
                                                   double /*net_tol*/, 
                                                   CubitBoolean /*heal*/ ) const
{
   PRINT_ERROR("Function not implemented in this engine.\n");
   return CUBIT_FAILURE;
}

//================================================================================
// Description: Creates a net surface.
// Author     : Tyronne Lim
// Date       : 08/18/03
//================================================================================
CubitStatus OCCModifyEngine::create_net_surface( DLIList<Curve*>& /*u_curves*/, 
                                                   DLIList<Curve*>& /*v_curves*/,
                                                   BodySM *& /*new_body*/, 
                                                   double /*net_tol*/, 
                                                   CubitBoolean /*heal*/ ) const
{
   PRINT_ERROR("Function not implemented in this engine.\n");
   return CUBIT_FAILURE;
}

//================================================================================
// Description: Creates an offset surface.
// Author     : Tyronne Lim
// Date       : 08/18/03
//================================================================================
CubitStatus OCCModifyEngine::create_offset_surface( Surface* /*ref_face_ptr*/, 
                                                      BodySM*& /*new_body*/, 
                                                      double /*offset_distance*/ ) const
{
   PRINT_ERROR("Function not implemented in this engine.\n");
   return CUBIT_FAILURE;
}

//================================================================================
// Description: Creates an offset body.
// Author     : Tyronne Lim
// Date       : 08/18/03
//================================================================================
CubitStatus OCCModifyEngine::create_offset_body( BodySM* body_ptr, 
                                                   BodySM*& new_bodysm, 
                                                   double offset_distance ) const
{
  return OCCModifyEngine::instance()->
    create_shell_offset( body_ptr, new_bodysm, offset_distance );
}

//================================================================================
// Description: Creates a skin surface.
// Author     : Tyronne Lim
// Date       : 08/18/03
//================================================================================
CubitStatus OCCModifyEngine::create_skin_surface( DLIList<Curve*>& /*curves*/, 
                                                    BodySM*& /*new_body*/ ) const
{
   PRINT_ERROR("Function not implemented in this engine.\n");
   return CUBIT_FAILURE;
}

//================================================================================
// Description: Creates a body from lofting surfaces.
// Author     : Tyronne Lim
// Date       : 08/18/03
//================================================================================
CubitStatus OCCModifyEngine::loft_surfaces( Surface * /*face1*/, 
                                              const double & /*takeoff1*/,
                                              Surface * /*face2*/, 
                                              const double & /*takeoff2*/,
                                              BodySM*& /*new_body*/,
                                              CubitBoolean /*arc_length_option*/, 
                                              CubitBoolean /*twist_option*/,
                                              CubitBoolean /*align_direction*/, 
                                              CubitBoolean /*perpendicular*/,
                                              CubitBoolean /*simplify_option*/ ) const
{
   PRINT_ERROR("Function not implemented in this engine.\n");
   return CUBIT_FAILURE;
}

//================================================================================
// Description: Creates a body by lofting surfaces between bodies.
// Author     : Tyronne Lim
// Date       : 08/18/03
//================================================================================
CubitStatus OCCModifyEngine::loft_surfaces_to_body( Surface * /*face1*/, 
                                                      const double & /*takeoff1*/,
                                                      Surface * /*face2*/, 
                                                      const double & /*takeoff2*/,
                                                      BodySM*& /*new_body*/,
                                                      CubitBoolean /*arc_length_option*/, 
                                                      CubitBoolean /*twist_option*/,
                                                      CubitBoolean /*align_direction*/, 
                                                      CubitBoolean /*perpendicular*/,
                                                      CubitBoolean /*simplify_option*/ ) const
{
   PRINT_ERROR("Function not implemented in this engine.\n");
   return CUBIT_FAILURE;
}
 
//================================================================================
// Description: Creates a surface.
// Author     : Tyronne Lim
// Date       : 08/18/03
//================================================================================
CubitStatus OCCModifyEngine::create_surface( DLIList<CubitVector*>& /*vec_list*/, 
                                               BodySM *& /*new_body*/, 
                                               Surface * /*ref_face_ptr*/,
                                               CubitBoolean /*project_points*/ ) const
{
   PRINT_ERROR("Function not implemented in this engine.\n");
   return CUBIT_FAILURE;
}

//================================================================================
// Description: Creates a weld surface.
// Author     : Tyronne Lim
// Date       : 08/18/03
//================================================================================
CubitStatus OCCModifyEngine::create_weld_surface( CubitVector & /*root*/,
                                                    Surface * /*ref_face1*/, 
                                                    double /*leg1*/, 
                                                    Surface * /*ref_face2*/, 
                                                    double /*leg2*/,
                                                    BodySM *& /*new_body*/ ) const
{
   PRINT_ERROR("Function not implemented in this engine.\n");
   return CUBIT_FAILURE;
}

//================================================================================
//    Facet-based geometry entities
//    Methods for building specific facet-based geometry entities
//================================================================================

//================================================================================
// Function   : make_facet_point
// Member Type: PUBLIC
// Description: create a new facet point given a pointer to the 
//              associated CubitPoint
// Author     : sjowen
// Date       : 12/28/00
//================================================================================
CubitStatus OCCModifyEngine::make_facet_point( CubitPoint *thePoint,
                                                 Point *&new_point_ptr )
{
  //We don't know in this function what to attach it to, so do that later.
  DLIList<Curve*> curves;

  OCCPoint *new_facet_point = new OCCPoint( thePoint, curves );

  new_point_ptr = (Point*) new_facet_point;

  return CUBIT_SUCCESS;
}

//================================================================================
// Function   : make_facet_point
// Member Type: PUBLIC
// Description: create a new facet point given its location
// Author     : sjowen
// Date       : 12/6/00
//================================================================================
CubitStatus OCCModifyEngine::make_facet_point( CubitVector &location,
                                                 Point *&new_point_ptr )
{
  //We don't know in this function what to attach it to, so do that later.
  DLIList<Curve*> curves;

  OCCPoint *new_facet_point = new OCCPoint( location, curves );

  new_point_ptr = (Point*) new_facet_point;

  return CUBIT_SUCCESS;
}

//================================================================================
// Function   : make_facet_curve
// Member Type: PUBLIC
// Description: create a new facet curve
//              (Assumes the curve will have an associated parent surface.
//               Evaluations will be done on the surface -- uses the same
//               order evaluations as the surface)
// Author     : sjowen
// Date       : 12/6/00
//================================================================================
CubitStatus OCCModifyEngine::make_facet_curve( Point *start_ptr,
                                                 Point *end_ptr, 
                                                 Curve *&new_curve_ptr,
                                                 CurveFacetEvalTool *eval_tool_ptr)
{

  //We don't know in this function what to attach it to, so do that later.
  DLIList<CoEdgeSM*> coedgesms; 

  OCCCurve *new_facet_curv = new OCCCurve( eval_tool_ptr, start_ptr, 
                                               end_ptr, coedgesms );
  new_curve_ptr = (Curve*) new_facet_curv;

  OCCPoint *facet_start_point = CAST_TO( start_ptr, OCCPoint );
  OCCPoint *facet_end_point = CAST_TO( end_ptr, OCCPoint );
  facet_start_point->add_curve( new_curve_ptr );
  facet_end_point->add_curve( new_curve_ptr );

  return CUBIT_SUCCESS;
}

//================================================================================
// Function   : make_facet_curve
// Member Type: PUBLIC
// Description: create a new facet curve given its edge facets
//              (this can be used for generating a facet curve representation
//               that does not have a parent surface associated -- currently
//               only linear evaluations are possible)
// Author     : sjowen
// Date       : 2/17/01
//================================================================================
CubitStatus OCCModifyEngine::make_facet_curve(
  Point *start_ptr,              // endpoints on the curve
  Point *end_ptr,
  DLIList<CubitFacetEdge*> &edge_list, // the ordered facet edges on this curve
  DLIList<CubitPoint*> &point_list,    // the ordered points on this curve
  Curve *&new_curve_ptr,          // return the new curve pointer
  CurveFacetEvalTool *curve_facet_tool )
{
  if ( curve_facet_tool == NULL )
  {
    curve_facet_tool = new CurveFacetEvalTool;
    CubitStatus status = curve_facet_tool->initialize( edge_list, point_list );
    if ( status != CUBIT_SUCCESS )
    {
        return status;
    }
  }
  if (!curve_facet_tool)
    return CUBIT_FAILURE;

  //We don't know in this function what to attach it to, so do that later.
  DLIList<CoEdgeSM*> coedgesms; 

  OCCCurve *new_facet_curv = new OCCCurve( curve_facet_tool, start_ptr, 
                                               end_ptr, coedgesms );
  new_curve_ptr = (Curve*) new_facet_curv;

  OCCPoint *facet_start_point = CAST_TO( start_ptr, OCCPoint );
  OCCPoint *facet_end_point = CAST_TO( end_ptr, OCCPoint );
  facet_start_point->add_curve( new_curve_ptr );
  facet_end_point->add_curve( new_curve_ptr );

  return CUBIT_SUCCESS;
}

//================================================================================
// Function   : make_facet_coedge
// Member Type: PUBLIC
// Description: create a new facet coedge
// Author     : sjowen
// Date       : 12/6/00
//================================================================================
CubitStatus OCCModifyEngine::make_facet_coedge( Curve *curv_ptr,
                                                    CubitSense sense,
                                                    CoEdgeSM *&new_coedge_ptr )
{
    //We don't know in this function what to attach it to, so do that later.
  LoopSM *loop_ptr = NULL;

    //Now create a new loop.
 
  OCCCoEdge *new_facet_coedge = new OCCCoEdge(curv_ptr, loop_ptr, sense);

  new_coedge_ptr = (CoEdgeSM*) new_facet_coedge;

  OCCCurve *facet_curve = CAST_TO( curv_ptr, OCCCurve );
  facet_curve->add_coedge( new_coedge_ptr );

  return CUBIT_SUCCESS;
}

//================================================================================
// Function   : make_facet_loop
// Member Type: PUBLIC
// Description: create a new facet loop
// Author     : sjowen
// Date       : 12/6/00
//================================================================================
CubitStatus OCCModifyEngine::make_facet_loop( DLIList<CoEdgeSM*> &coedge_list,
                                                  LoopSM *&new_loop_ptr )
{
    //We don't know in this function what to attach it to, so do that later.
  Surface *surf_ptr = NULL;

    //Now create a new loop.
 
  OCCLoop *new_facet_loop = new OCCLoop(surf_ptr, coedge_list );

  new_loop_ptr = (LoopSM*) new_facet_loop;
  int ii;
  for ( ii = coedge_list.size(); ii > 0; ii-- )
  {
    CoEdgeSM* coedge_ptr = coedge_list.get_and_step();
    OCCCoEdge *facet_coedge = CAST_TO(coedge_ptr, OCCCoEdge);
    facet_coedge->add_loop(new_loop_ptr);
  }
  return CUBIT_SUCCESS;
}

//================================================================================
// Function   : make_facet_surface
// Member Type: PUBLIC
// Description: create a new facet surface from Quad facets
// Author     : sjowen
// Date       : 12/6/00
//================================================================================
CubitStatus OCCModifyEngine::make_facet_surface(
  DLIList<CubitQuadFacet*> &quad_facet_list,
  DLIList<CubitPoint*> &point_list,
  DLIList<LoopSM*> &my_loops,
  int interp_order,
  double min_dot,
  Surface *&new_surface_ptr)
{
  CubitQuadFacet *qfacet;
  DLIList<CubitFacet*> facet_list;
  int ii;

  quad_facet_list.reset();
  for (ii=0; ii<quad_facet_list.size(); ii++)
  {
    qfacet = quad_facet_list.get_and_step();
    facet_list.append( qfacet->get_tri_facet( 0 ) );
    facet_list.append( qfacet->get_tri_facet( 1 ) );
  }
  return make_facet_surface( NULL, facet_list, point_list, my_loops, 
                             interp_order, min_dot,  new_surface_ptr );
}

//================================================================================
// Function   : make_facet_surface
// Member Type: PUBLIC
// Description: create a new facet surface
// Author     : sjowen
// Date       : 12/6/00
//================================================================================
CubitStatus OCCModifyEngine::make_facet_surface( const CubitEvaluatorData *eval_data,
                                                   DLIList<CubitFacet*> &facet_list,
                                                   DLIList<CubitPoint*> &point_list,
                                                   DLIList<LoopSM*> &my_loops,
                                                   int interp_order,
                                                   double min_dot,
                                                   Surface *&new_surface_ptr,
                                                   CubitBoolean use_point_addresses,
                                                   FacetEvalTool *facet_eval_tool)
{
    //Create a new surface given the facets and point list.
  if (facet_eval_tool == NULL)
  {
    facet_eval_tool = new FacetEvalTool(facet_list, point_list, 
                                         interp_order, min_dot );
  }
    
    //We don't know in this function what to attach it to, so do that later.
  DLIList<ShellSM*> shellsms;
  OCCSurface *new_facet_surf;

    //Now create a new surface.
  if ( eval_data && eval_data->ask_type() == SPHERE_SURFACE_TYPE )
  {
      new_facet_surf = new OCCSurface( (const SphereEvaluatorData*)eval_data,
                                         facet_eval_tool,
                                         shellsms,
                                         my_loops );
  }
  else if ( eval_data && eval_data->ask_type() == CONE_SURFACE_TYPE )
  {
      new_facet_surf = new OCCSurface( (const CylinderEvaluatorData*)eval_data,
                                         facet_eval_tool,
                                         shellsms,
                                         my_loops );
  }
  else
  {
      new_facet_surf = new OCCSurface( facet_eval_tool,
                                         shellsms, my_loops );
  }

  new_surface_ptr = (Surface*) new_facet_surf;
  int ii;
  for ( ii = my_loops.size(); ii > 0; ii-- )
  {
    LoopSM* loop_ptr = my_loops.get_and_step();
    OCCLoop *facet_loop = CAST_TO(loop_ptr, OCCLoop);
    facet_loop->add_surface(new_surface_ptr);
#ifdef BOYD17 
    DLIList<CoEdgeSM*> coedges;
#endif
  }

   // generate the curve facet evaluaion tools for each of the curves on this
   // surface - they are based on the surface evaluation tool

  DLIList<OCCCoEdge*> coedges;
  new_facet_surf->get_coedges( coedges );
  for (ii=0; ii<coedges.size(); ii++)
  {
    OCCCoEdge *coedgesm_ptr = coedges.get_and_step();
    DLIList<OCCCurve*> curve_list;
    coedgesm_ptr->get_curves( curve_list );
    OCCCurve *facet_curve =curve_list.get();
    CurveFacetEvalTool *eval_tool = facet_curve->get_eval_tool(); 
    if (eval_tool == NULL) 
    {
      Point *start_point = facet_curve->start_point();
      Point *end_point = facet_curve->end_point();
      OCCCoEdge *facet_coedge_ptr = CAST_TO( coedgesm_ptr, OCCCoEdge );
      CubitSense sense = facet_coedge_ptr->get_sense();

      CubitStatus status = CUBIT_SUCCESS;
      eval_tool = new CurveFacetEvalTool;
      if (!use_point_addresses)
      {
        CubitVector temp_vec1 = start_point->coordinates();
        CubitVector temp_vec2 = end_point->coordinates();
        status = eval_tool->initialize( facet_eval_tool,
                                        temp_vec1,
                                        temp_vec2,
                                        sense );
        if ( status != CUBIT_SUCCESS )
        {
            return status;
        }
      }
      else
      {
        OCCPoint *start_facet_pt_ptr = CAST_TO( start_point, OCCPoint );
        OCCPoint *end_facet_pt_ptr = CAST_TO( end_point, OCCPoint );
        CubitPoint *start_pt = start_facet_pt_ptr->get_cubit_point();
        CubitPoint *end_pt = end_facet_pt_ptr->get_cubit_point();
        status = eval_tool->initialize( facet_eval_tool,
                                        start_pt,
                                        end_pt,
                                        sense );
        if ( status != CUBIT_SUCCESS )
        {
            return status;
        }
      }
      if( !eval_tool->has_good_curve_data() )
        return CUBIT_FAILURE;
      facet_curve->set_eval_tool( eval_tool );
    }
  }
  return CUBIT_SUCCESS;
}

//================================================================================
// Function   : make_facet_shell
// Member Type: PUBLIC
// Description: create a new facet shell
// Author     : sjowen
// Date       : 12/6/00
//================================================================================
CubitStatus OCCModifyEngine::make_facet_shell(DLIList<Surface*> &surface_list,
                                                  ShellSM *&new_shell_ptr)
{
  
  Lump* my_lump = NULL;
  OCCShell *new_facet_shell = new OCCShell(my_lump, surface_list);
  new_shell_ptr = (ShellSM*) new_facet_shell;
  int ii;
  for ( ii = surface_list.size(); ii > 0; ii-- )
  {
    Surface* surface_ptr = surface_list.get_and_step();
    OCCSurface *facet_surface = CAST_TO(surface_ptr, OCCSurface);
    facet_surface->add_shell(new_shell_ptr);
  }
  return CUBIT_SUCCESS;
}

//================================================================================
// Function   : make_facet_lump
// Member Type: PUBLIC
// Description: create a new facet lump
// Author     : sjowen
// Date       : 12/6/00
//================================================================================
CubitStatus OCCModifyEngine::make_facet_lump(DLIList<ShellSM*> &shell_list,
                                                 Lump *&new_lump_ptr)
{
  
  OCCLump *new_facet_lump = new OCCLump(shell_list);
  new_lump_ptr = (Lump*) new_facet_lump;
  int ii;
  for ( ii = shell_list.size(); ii > 0; ii-- )
  {
    ShellSM* shell_ptr = shell_list.get_and_step();
    OCCShell *facet_shell = CAST_TO(shell_ptr, OCCShell);
    facet_shell->add_lump(new_lump_ptr);
  }
  return CUBIT_SUCCESS;
}

//================================================================================
// Function   : make_facet_body
// Member Type: PUBLIC
// Description: create a new facet body
// Author     : sjowen
// Date       : 12/6/00
//================================================================================
CubitStatus OCCModifyEngine::make_facet_body(DLIList<Lump*> &lump_list,
                                                 BodySM *&new_body_ptr)
{
  
  OCCBody *new_facet_body = new OCCBody(lump_list);
  new_body_ptr = (BodySM*) new_facet_body;
    //Now attach the lower entites.

  int ii;
  for ( ii = lump_list.size(); ii > 0; ii-- )
  {
    Lump* lump_ptr = lump_list.get_and_step();
    OCCLump *facet_lump = CAST_TO(lump_ptr, OCCLump);
    facet_lump->add_body(new_body_ptr);
  }
  return CUBIT_SUCCESS;
}

//================================================================================
// Function   : build_facet_surface
// Member Type: PUBLIC
// Description: same as build_facet_surface below, but takes a list of quad facets
// Author     : sjowen
// Date       : 4/30/01
//================================================================================
CubitStatus OCCModifyEngine::build_facet_surface(
  DLIList<CubitQuadFacet *> &quad_facet_list,
  DLIList<CubitPoint *>&point_list,  // fills this in if not provided
  double feature_angle,
  int interp_order,
  CubitBoolean smooth_non_manifold,
  CubitBoolean split_surfaces,
  DLIList<Surface *> &surface_list) 
{
  CubitQuadFacet *qfacet;
  DLIList<CubitFacet*> facet_list;
  int ii;

  quad_facet_list.reset();
  for (ii=0; ii<quad_facet_list.size(); ii++)
  {
    qfacet = quad_facet_list.get_and_step();
    facet_list.append( qfacet->get_tri_facet( 0 ) );
    facet_list.append( qfacet->get_tri_facet( 1 ) );
  }
  return build_facet_surface( NULL, facet_list, point_list, feature_angle, 
                             interp_order, smooth_non_manifold,
                             split_surfaces, surface_list );
}

//================================================================================
// Function   : build_facet_surface
// Member Type: PUBLIC
// Description: same as build_facet_surface below, both quad and triangle facets
// Author     : sjowen
// Date       : 4/30/01
//================================================================================
CubitStatus OCCModifyEngine::build_facet_surface(
  DLIList<CubitQuadFacet *> &quad_facet_list,
  DLIList<CubitFacet *> &tri_facet_list,
  DLIList<CubitPoint *>&point_list,  // fills this in if not provided
  double feature_angle,
  int interp_order,
  CubitBoolean smooth_non_manifold,
  CubitBoolean split_surfaces,
  DLIList<Surface *> &surface_list) 
{
  CubitStatus rv = CUBIT_SUCCESS;
  CubitQuadFacet *qfacet;
  DLIList<CubitFacet*> facet_list;
  int ii;
  quad_facet_list.reset();
  for (ii=0; ii<quad_facet_list.size(); ii++)
  {
    qfacet = quad_facet_list.get_and_step();
    facet_list.append( qfacet->get_tri_facet( 0 ) );
    facet_list.append( qfacet->get_tri_facet( 1 ) );
  }
  facet_list += tri_facet_list;

  rv = build_facet_surface( NULL, facet_list, point_list, feature_angle, 
                            interp_order, smooth_non_manifold,
                            split_surfaces, surface_list );
  return rv; 
}

//================================================================================
// Function   : build_cholla_surfaces
// Member Type: PUBLIC
// Description: create one or more surfaces from the list of facets.  Also creates
//              the lower order entities.  Generate curves (split surfaces), based
//              on feature angle.  If feature_angle < 0 then ignore feature_angle
// Note       : This function is the main interface for Cubit to the Cholla library
// Author     : sjowen
// Date       : 4/26/01
//================================================================================
CubitStatus OCCModifyEngine::build_cholla_surfaces
(
  DLIList<CubitFacet *> facet_list,
  DLIList<CubitPoint *> point_list,  // fills this in if not provided
  double feature_angle,
  int interp_order,
  CubitBoolean smooth_non_manifold,
  CubitBoolean split_surfaces,
  ChollaEngine *&cholla_ptr
)
{
  CubitStatus rv = CUBIT_SUCCESS;

  // generate the edge and point lists from the facets

  DLIList<CubitFacetEdge *> edge_list;
  DLIList<FacetEntity*> pe_list;
  DLIList<FacetEntity *> pf_list;
  FacetDataUtil::get_edges( facet_list, edge_list );
  if (facet_list.size() > 0 && point_list.size() == 0)
    FacetDataUtil::get_points( facet_list, point_list );

  CAST_LIST( point_list, pf_list, FacetEntity );
  CAST_LIST( edge_list, pe_list, FacetEntity );

  // create the surfaces

  if (rv == CUBIT_SUCCESS)
  {
    DLIList<FacetEntity*> facet_entity_list;
    CAST_LIST( facet_list, facet_entity_list, FacetEntity );
    cholla_ptr = new ChollaEngine( facet_entity_list, 
                                   pe_list, pf_list );
    CubitBoolean use_feature_angle;
    if (feature_angle < 0.0)
      use_feature_angle = CUBIT_FALSE;
    else
      use_feature_angle = CUBIT_TRUE;
    rv = cholla_ptr->create_geometry( use_feature_angle, feature_angle, 
                                      interp_order, smooth_non_manifold, 
                                      split_surfaces );
  }
  return rv;
}
//================================================================================
// Function   : build_facet_surface
// Member Type: PUBLIC
// Description: create one or more surfaces from the list of facets.  Also creates
//              the lower order entities.  Generate curves (split surfaces), based
//              on feature angle.  If feature_angle < 0 then ignore feature_angle
// Note       : This function is the main interface for Cubit to the Cholla library
// Author     : sjowen
// Date       : 4/26/01
//================================================================================
CubitStatus OCCModifyEngine::build_facet_surface(
  const CubitEvaluatorData **eval_data,
  DLIList<CubitFacet *> &facet_list,
  DLIList<CubitPoint *> &point_list,  // fills this in if not provided
  double feature_angle,
  int interp_order,
  CubitBoolean smooth_non_manifold,
  CubitBoolean split_surfaces,
  DLIList<Surface *> &surface_list) 
{
  CubitStatus rv = CUBIT_SUCCESS;

  // generate the edge and point lists from the facets

  DLIList<CubitFacetEdge *> edge_list;
  DLIList<FacetEntity*> pe_list;
  DLIList<FacetEntity *> pf_list;
  FacetDataUtil::get_edges( facet_list, edge_list );
  if (point_list.size() == 0)
    FacetDataUtil::get_points( facet_list, point_list );

  CAST_LIST( point_list, pf_list, FacetEntity );
  CAST_LIST( edge_list, pe_list, FacetEntity );

  // create the surfaces

  if (rv == CUBIT_SUCCESS)
  {
    DLIList<FacetEntity*> facet_entity_list;
    CAST_LIST( facet_list, facet_entity_list, FacetEntity );
    ChollaEngine *cholla_ptr = new ChollaEngine( facet_entity_list, 
                                                 pe_list, pf_list );
    cholla_ptr->set_flip_flag(CUBIT_TRUE);
    CubitBoolean use_feature_angle;
    if (feature_angle < 0.0)
      use_feature_angle = CUBIT_FALSE;
    else
      use_feature_angle = CUBIT_TRUE;
    rv = cholla_ptr->create_geometry( use_feature_angle, feature_angle, 
                                      interp_order, smooth_non_manifold, 
                                      split_surfaces );
    
    // make CUBIT geometry out of the Cholla entities

    if (rv == CUBIT_SUCCESS)
    {
      DLIList<ChollaSurface *> cholla_surface_list;
      cholla_ptr->get_surfaces( cholla_surface_list );
      DLIList<ChollaCurve *> cholla_curve_list;
      cholla_ptr->get_curves( cholla_curve_list );
      DLIList<ChollaPoint *> cholla_point_list;
      cholla_ptr->get_points( cholla_point_list );
      rv = build_cholla_geometry(eval_data, cholla_surface_list, cholla_curve_list,
                                 cholla_point_list, use_feature_angle, 
                                 feature_angle, interp_order, surface_list);
    }
    cholla_ptr->delete_me();
    delete cholla_ptr;
  }
  return rv;
}

//===============================================================================
// Function   : build_cholla_geometry
// Member Type: PRIVATE
// Description: build the CUBIT geometry based on the Facet entity class lists
// Author     : sjowen
// Date       : 10/02
//===============================================================================
CubitStatus OCCModifyEngine::build_cholla_geometry(
  const CubitEvaluatorData **eval_data,
  DLIList<ChollaSurface*> &cholla_surface_list,
  DLIList<ChollaCurve*> &cholla_curve_list,
  DLIList<ChollaPoint*> &cholla_point_list,
  CubitBoolean use_feature_angle, 
  double feature_angle, 
  int interp_order,
  DLIList<Surface *> &surface_list)
{
  CubitStatus stat = CUBIT_SUCCESS;

  //- convert feature angle to a dot product

  double min_dot = 0.0;
  if (use_feature_angle)
  {
    if (feature_angle > 180.0) feature_angle = 180.0;
    if (feature_angle < 0.0) feature_angle = 0.0;
    double rad_angle = (180.0 - feature_angle) * CUBIT_PI / 180.0;
    min_dot = cos( rad_angle );
  }

  // build the CUBIT geometry based on Cholla entities

  stat = build_cholla_point_geometry( cholla_point_list );
  
  if (stat == CUBIT_SUCCESS)
    stat = build_cholla_curve_geometry( cholla_curve_list );

  if (stat == CUBIT_SUCCESS)
    stat = build_cholla_surface_geometry( eval_data,
                                          cholla_surface_list,
                                          interp_order, min_dot,
                                          surface_list);

  return stat;
}

//===============================================================================
// Function   : build_cholla_point_geometry
// Member Type: PRIVATE
// Description: From the cholla point list, create geometric points for each
// Author     : sjowen
// Date       : 10/02
//===============================================================================
CubitStatus OCCModifyEngine::build_cholla_point_geometry(
  DLIList<ChollaPoint*> &cholla_point_list )
{
  CubitStatus stat = CUBIT_SUCCESS;
  int kk;
  //int mydebug = 0;
  for ( kk = cholla_point_list.size(); kk > 0; kk-- )
  {
    ChollaPoint *chpnt_ptr = cholla_point_list.get_and_step();
    Point *point_ptr = (Point *)chpnt_ptr->get_geometric_point();
    if (point_ptr == NULL)
    {
      FacetEntity *facet_ptr = chpnt_ptr->get_facets();
      CubitPoint *cp_ptr = CAST_TO( facet_ptr, CubitPoint );
      Point *point = NULL;
      stat = make_facet_point( cp_ptr, point );
      if ( point == NULL || stat != CUBIT_SUCCESS )
      {
        PRINT_ERROR("Problems building mesh based points.\n");
        return stat;
      }
      point_ptr = (Point *)point;
      chpnt_ptr->assign_geometric_point((void *)point_ptr);
    }
  }
  return stat;
}

//===============================================================================
// Function   : build_cholla_curve_geometry
// Member Type: PRIVATE
// Description: From the cholla curve list, create geometric curves for each
// Author     : sjowen
// Date       : 10/02
//===============================================================================
CubitStatus OCCModifyEngine::build_cholla_curve_geometry(
  DLIList<ChollaCurve*> &cholla_curve_list )
{
  CubitStatus stat = CUBIT_SUCCESS;
  int kk;
  for ( kk = cholla_curve_list.size(); kk > 0; kk-- )
  {
    ChollaCurve *chcurv_ptr = cholla_curve_list.get_and_step();
    Curve *curv_ptr = (Curve *)chcurv_ptr->get_geometric_curve();
    if (curv_ptr == NULL)
    {
      DLIList<ChollaPoint*> fpoint_list = chcurv_ptr->get_points( );
      ChollaPoint *fpm0_ptr = fpoint_list.get_and_step();
      ChollaPoint *fpm1_ptr = fpoint_list.get_and_step();
      CubitPoint *start_point, *end_point;
      chcurv_ptr->get_ends( start_point, end_point );
      if (fpm0_ptr->get_facets() != start_point)
      {
        ChollaPoint *temp_ptr;
        temp_ptr = fpm0_ptr;
        fpm0_ptr = fpm1_ptr;
        fpm1_ptr = temp_ptr;
      }
      Point *start_ptr = (Point *) fpm0_ptr->get_geometric_point();
      Point *end_ptr = (Point *) fpm1_ptr->get_geometric_point();
      
      // if this is a curve without a parent surface then handle it 
      // differently.  (Curves with parents use the surface to evaluate to
      // With only a curve, it must evaluate to the curve)
      
      DLIList<ChollaSurface*> fsm_list = chcurv_ptr->get_surfaces();
      if (fsm_list.size() == 0)
      {
        DLIList<FacetEntity*> facet_list;
        DLIList<CubitPoint*> point_list;  // needs to be filled in
        facet_list = chcurv_ptr->get_facet_list();
        DLIList<CubitFacetEdge*> edge_list;
        CAST_LIST( facet_list, edge_list, CubitFacetEdge );
        if (stat != CUBIT_SUCCESS)
          return stat;

        // the CurveFacetEvalTool was computed in Cholla.  Retreive it
        // and pass it to make_facet_curve so it can store it with the 
        // new curve
        CurveFacetEvalTool *eval_tool_ptr = chcurv_ptr->get_eval_tool();
        if (eval_tool_ptr == NULL)
          return CUBIT_FAILURE;
        stat = make_facet_curve( start_ptr, end_ptr,
                                 edge_list, point_list,
                                 curv_ptr, eval_tool_ptr );
        if (stat == CUBIT_SUCCESS)
        {
          Curve *curvsm_ptr = CAST_TO( curv_ptr, Curve );
          GeometryQueryTool::instance()->make_free_RefEdge( curvsm_ptr );
        }
      }
      else
      {
        CurveFacetEvalTool *eval_tool_ptr = chcurv_ptr->get_eval_tool();
        if (eval_tool_ptr == NULL)
          return CUBIT_FAILURE;
        stat = make_facet_curve( start_ptr, end_ptr,
                                 curv_ptr, eval_tool_ptr );
        if ( curv_ptr == NULL || stat != CUBIT_SUCCESS )
        {
          PRINT_ERROR("Problems building mesh based curves.\n");
          return stat;
        }
      }
      chcurv_ptr->assign_geometric_curve((void *)curv_ptr);
    }
  }
  return stat;
}

//===============================================================================
// Function   : build_cholla_loop_geometry
// Member Type: PRIVATE
// Description: From the cholla curve list of a surface, create geometric loops 
// Author     : sjowen
// Date       : 10/02
//===============================================================================
CubitStatus OCCModifyEngine::build_cholla_loop_geometry(
  DLIList<ChollaCurve*> &cholla_curve_list,   // curves on a surface
  ChollaSurface *chsurf_ptr,                    // the surface
  DLIList<LoopSM*> &loop_list,                  // append to this loop list
  int &ncurves,                                 // number of curves used
  int debug_draw )
{

  // find the first unused curve on the list.  Use that as the starting curve
  // in the loop

  ncurves = 0;
  CubitStatus stat = CUBIT_SUCCESS;
  int ii, jj;
  int found = 0;
  ChollaCurve *chcurv_ptr = NULL;
  for(ii=0; ii<cholla_curve_list.size() && !found; ii++)
  {
    chcurv_ptr = cholla_curve_list.get_and_step();
    if (chcurv_ptr->get_flag() == 0)
      found = 1;
  }
  if (!found)
    return CUBIT_FAILURE;  // all the curves have already been used in this surface

  // get the orientation of this curve w.r.t. this surface
  
  CubitSense orientation;
  stat = ChollaEngine::determine_curve_orientation( chsurf_ptr, 
                                   chcurv_ptr, orientation );
  if (orientation == CUBIT_UNKNOWN  || stat == CUBIT_FAILURE)
    return CUBIT_FAILURE;

  CubitPoint *start_ptr, *end_ptr;
  chcurv_ptr->get_ends( start_ptr, end_ptr );

  // create a new coedge

  CoEdgeSM *coedge_ptr = NULL;
  Curve *curv_ptr = (Curve *)chcurv_ptr->get_geometric_curve();
  stat = make_facet_coedge( curv_ptr, orientation, coedge_ptr );
  if ( coedge_ptr == NULL || stat != CUBIT_SUCCESS )
  {
    PRINT_ERROR("Problems building mesh based coedges.\n");
    return stat;
  }

  if (debug_draw)
  {
    CubitVector s = start_ptr->coordinates();
    CubitVector e = end_ptr->coordinates();
    CubitSense sense = coedge_ptr->sense();
    if (sense == CUBIT_FORWARD)
      GfxDebug::draw_line(s, e, CUBIT_GREEN);
    else
      GfxDebug::draw_line(s, e, CUBIT_RED);

    GfxDebug::draw_point(s, CUBIT_GREEN);
    GfxDebug::draw_label(0, s.x(), s.y(), s.z(), CUBIT_GREEN);
    
    GfxDebug::draw_point(e, CUBIT_RED);
    GfxDebug::draw_label(1, e.x(), e.y(), e.z(), CUBIT_RED);

    GfxDebug::flush();
    debug_draw=debug_draw;
  }

  // start a list of coedges for this loop

  DLIList<CoEdgeSM*> coedge_list;
  coedge_list.append( coedge_ptr );
  chcurv_ptr->set_flag(1);
  ChollaCurve *last_chcurv_ptr = chcurv_ptr;

  // loop through and determine the rest of the coedges based on the
  // orientation of the first curve in the loop
   
  chcurv_ptr->get_ends( start_ptr, end_ptr );
  FacetEntity *n0_ptr, *n1_ptr;
  if (orientation == CUBIT_FORWARD)
  {
    n0_ptr = start_ptr;
    n1_ptr = end_ptr;
  }
  else
  {
    n1_ptr = start_ptr;
    n0_ptr = end_ptr;
  }

  while (n0_ptr != n1_ptr)
  {

    // find the next curve in the loop

    int ntimes = 0;
    found = 0;
    while(!found)
    {
      for ( jj = 0; jj < cholla_curve_list.size() && !found; jj++)
      {
        chcurv_ptr = cholla_curve_list.get_and_step();
        if (last_chcurv_ptr != chcurv_ptr && (chcurv_ptr->get_flag() == 0 || ntimes > 0))
        {
          chcurv_ptr->get_ends( start_ptr, end_ptr );
          if (start_ptr == n1_ptr)
          {
            found = 1;
            orientation = CUBIT_FORWARD;
            n1_ptr = end_ptr;
          }
          else if (end_ptr == n1_ptr)
          {
            found = 1;
            orientation = CUBIT_REVERSED;
            n1_ptr = start_ptr;
          }
        }
      }
      if (!found)
      {

        // in this special case - we have gone through the whole list without
        // finding the next curve on the loop.  In the first pass we have assumed
        // that each curve is only used once.  This may not be the case.  Go
        // back and see of we can find a curve that has already been used that
        // will work.  If we fail the second time through, then we fail!

        ntimes++;
        if (ntimes > 1)
          // return CUBIT_FAILURE;  // didn't find the next curve in the loop
          // DEBUG ONLY -- THIS NEEDS TO BE FIXED
          goto create_loop;
      }
    }

    // create a new coedge

    coedge_ptr = NULL;
    curv_ptr = (Curve *)chcurv_ptr->get_geometric_curve();
    stat = make_facet_coedge( curv_ptr, orientation,coedge_ptr );
    if ( coedge_ptr == NULL || stat != CUBIT_SUCCESS )
    {
      PRINT_ERROR("Problems building mesh based coedges.\n");
      return stat;
    }

    if (debug_draw)
    {
      CubitVector s = start_ptr->coordinates();
      CubitVector e = end_ptr->coordinates();
      CubitSense sense = coedge_ptr->sense();
      if (sense == CUBIT_FORWARD)
        GfxDebug::draw_line(s, e, CUBIT_GREEN);
      else
        GfxDebug::draw_line(s, e, CUBIT_RED);

      GfxDebug::draw_point(s, CUBIT_GREEN);
      GfxDebug::draw_label(0, s.x(), s.y(), s.z(), CUBIT_GREEN);
    
      GfxDebug::draw_point(e, CUBIT_RED);
      GfxDebug::draw_label(1, e.x(), e.y(), e.z(), CUBIT_RED);

      GfxDebug::flush();
      debug_draw=debug_draw;
    }
    coedge_list.append( coedge_ptr );
    if (coedge_list.size() > 2*cholla_curve_list.size())
      return CUBIT_FAILURE;  // more curves in the loop than their are curves ????
    chcurv_ptr->set_flag(1);
    last_chcurv_ptr = chcurv_ptr;
  } 

  // create the loop

create_loop:
  LoopSM *loop_ptr = NULL;
  stat = make_facet_loop( coedge_list, loop_ptr );
  if ( coedge_ptr == NULL || stat != CUBIT_SUCCESS )
  {
    PRINT_ERROR("Problems building mesh based loops.\n");
    return stat;
  }
  loop_list.append( loop_ptr );
  ncurves = coedge_list.size();
      
  return stat;
}

//===============================================================================
// Function   : build_cholla_surface_geometry (PRIVATE)
// Member Type: PRIVATE
// Description: From the facet surface list, create geometric surface,
//              loops and coedges for each surface in the list
// Author     : sjowen
// Date       : 10/02
//===============================================================================
CubitStatus OCCModifyEngine::build_cholla_surface_geometry(
  const CubitEvaluatorData **eval_data,
  DLIList<ChollaSurface*> &cholla_surface_list,
  int interp_order,
  double min_dot,
  DLIList<Surface *> &surface_list)
{
  CubitStatus stat = CUBIT_SUCCESS;
  int ii, isurface;

  // make sure the facet flags have been reset

  ChollaCurve *chcurv_ptr;
  cholla_surface_list.reset();
  for ( isurface = 0; isurface < cholla_surface_list.size(); isurface++ )
  {
    ChollaSurface *chsurf_ptr = cholla_surface_list.get_and_step();
    chsurf_ptr->reset_facet_flags();
  }

  // now loop through surfaces and create them

  int mydebug = 0;
  cholla_surface_list.reset();
  for ( isurface = 0; isurface < cholla_surface_list.size(); isurface++ )
  {
    const CubitEvaluatorData *surf_eval_data = eval_data == NULL ? NULL : eval_data[isurface];
    int num_curves = 0;
    ChollaSurface *chsurf_ptr = cholla_surface_list.get_and_step();

    DLIList<ChollaCurve*> chcurve_list;
    chsurf_ptr->get_curves( chcurve_list );

    // init flags for curves on this surface

    for (ii=0; ii< chcurve_list.size(); ii++)
    {
      chcurv_ptr = chcurve_list.get_and_step();
      chcurv_ptr->set_flag(0);
      if (mydebug)
      {
        chcurv_ptr->debug_draw();
      }
    }

    // generate loops - keep going until we have used all the curves in the list

    DLIList<LoopSM*> loop_list;
    int debug_draw = 0;
    while (num_curves < chcurve_list.size())
    {
      int ncurves = 0;
      stat = build_cholla_loop_geometry( chcurve_list, chsurf_ptr, 
                                         loop_list, ncurves, debug_draw );
      if (stat != CUBIT_SUCCESS)
      {
        PRINT_ERROR("Can't build geometric loop from facets\n");
        return stat;
      }
      num_curves += ncurves;
    }

    // create the surface

    Surface *surf_ptr = (Surface *)chsurf_ptr->get_geometric_surface();
    if (surf_ptr == NULL)
    {
      DLIList<FacetEntity*> facet_entity_list;
      DLIList<CubitPoint*> point_list;
      chsurf_ptr->get_points(point_list);
      chsurf_ptr->get_facets(facet_entity_list);
      DLIList<CubitFacet*> facet_list;
      CAST_LIST( facet_entity_list, facet_list, CubitFacet );
      FacetEvalTool *eval_tool_ptr = chsurf_ptr->get_eval_tool();
      if (eval_tool_ptr == NULL)
        return CUBIT_FAILURE;
      stat = make_facet_surface(surf_eval_data, facet_list,point_list,loop_list,
                                interp_order,min_dot,surf_ptr, CUBIT_TRUE, eval_tool_ptr);
      if ( surf_ptr == NULL || stat != CUBIT_SUCCESS )
      {
        PRINT_ERROR("Problems building mesh based surfaces.\n");
        return stat;
      }
      chsurf_ptr->assign_geometric_surface((void *)surf_ptr);
      surface_list.append(surf_ptr);
    }
  }
  return stat;
}

//================================================================================
// Description: improve the facets on a surface
// Author     : Steve Owen
// Date       : 8/28/2003
//================================================================================
CubitStatus OCCModifyEngine::improve_facets( RefFace *ref_face_ptr )
{
  //CubitStatus rv = CUBIT_SUCCESS;

  // This must be a facet-based surface - return now if it isn't

  Surface *surf_ptr = ref_face_ptr->get_surface_ptr();
  OCCSurface *fsurf_ptr = dynamic_cast<OCCSurface *>(surf_ptr);
  if (fsurf_ptr == NULL)
  {
    PRINT_ERROR("Couldn't improve facets on surface %d. It is not a facet-based geometry.\n",
                 ref_face_ptr->id());
    return CUBIT_FAILURE;
  }

  CubitPoint *point_ptr;
  DLIList<CubitPoint *> point_list;
  DLIList<CubitFacet *> facet_list;
  DLIList<CubitFacetEdge *> edge_list;

  // get the facets from the facet-based surface

  fsurf_ptr->get_my_facets(facet_list, point_list);
  FacetDataUtil::get_edges(facet_list, edge_list);

  // compute normals if necessary  

  int ii;
  for (ii=0; ii<point_list.size(); ii++)
  {
    point_ptr = point_list.get_and_step();
    point_ptr->normal();
  }

  // main clean loop

  CubitFacet *facet0, *facet1;
  CubitFacetData *cfacet0, *cfacet1;
  CubitPoint *pt0, *pt1, *pt2, *pt3;
  CubitFacetEdge *edge;
  CubitFacetEdgeData *cedge;
  int nflip = 0;
  CubitVector n0, n1, n2, n3;
  CubitVector c0, c1, c2, c3;
  double q0, q1, q2, q3;
  CubitVector nt0, nt1, nt2, nt3;
  DLIList<CubitFacet *>adj_facets;
  int nedges = edge_list.size();
  for(ii=0; ii<nedges; ii++)
  {
    edge = edge_list.get();

    // can't be a feature (edge is on a curve)
    if(edge->is_feature())
    {
      edge_list.step();
      continue;
    }
    
    pt0 = edge->point(0);
    pt1 = edge->point(1);
    adj_facets.clean_out();
    edge->facets( adj_facets );

    // must be 2 adjacent facets.
    if (adj_facets.size() != 2)
    {
      edge_list.step();
      continue;
    }
    facet0 = adj_facets.get_and_step();
    facet1 = adj_facets.get();

    // get the adjacent vertices, make sure we are oriented correctly 
    // - swap pt0 and pt1 if necessary
    pt2 = facet0->next_node(pt1);
    if(pt2 == pt0)
    {
      pt2 = pt0;
      pt0 = pt1;
      pt1 = pt2;
      pt2 = facet0->next_node(pt1);
    }
    pt3 = facet1->next_node(pt0);

    // get the normals from the vertices
    n0 = pt0->normal( facet0 );
    n1 = pt1->normal( facet0 );
    n2 = pt2->normal( facet0 );
    n3 = pt3->normal( facet1 );
    c0 = pt0->coordinates();
    c1 = pt1->coordinates();
    c2 = pt2->coordinates();
    c3 = pt3->coordinates();

    // compute triangle normals by averaging vertices
    nt0 = n0 + n1 + n2;
    nt0.normalize();
    nt1 = n0 + n1 + n3;
    nt1.normalize();
    nt2 = n0 + n2 + n3;
    nt2.normalize();
    nt3 = n1 + n2 + n3;
    nt3.normalize();

    // compute quality for each of the four tris... 2 existing tris
    // and the two potential tris if we flip
    q0 = FacetDataUtil::quality(c0, c1, c2, nt0);
    q1 = FacetDataUtil::quality(c1, c0, c3, nt1);
    q2 = FacetDataUtil::quality(c0, c3, c2, nt2);
    q3 = FacetDataUtil::quality(c1, c2, c3, nt3);
    q0 = CUBIT_MIN( q0, q1 );
    q2 = CUBIT_MIN( q2, q3 );
    if (q2 > q0)
    {
      // flip

      facet_list.remove( facet0 );  // these take too long - find another way to remove them
      facet_list.remove( facet1 );
      edge_list.change_to(NULL);

      cfacet0 = dynamic_cast<CubitFacetData *>( facet0 );
      cfacet1 = dynamic_cast<CubitFacetData *>( facet1 );
      cedge = dynamic_cast<CubitFacetEdgeData *>( edge );
      
      delete cfacet0;
      delete cfacet1;
      delete cedge;

      // get edges

      CubitFacetEdge *e[4];
      e[0] = pt2->get_edge( pt0 );
      e[1] = pt0->get_edge( pt3 );
      e[2] = pt3->get_edge( pt1 );
      e[3] = pt1->get_edge( pt2 );
      cedge = new CubitFacetEdgeData( pt2, pt3 );
      edge = dynamic_cast<CubitFacetEdge *>( cedge ); 

      cfacet0 = new CubitFacetData( edge, e[0], e[1] );
      cfacet1 = new CubitFacetData( edge, e[2], e[3] );
      facet0 = dynamic_cast<CubitFacet *>( cfacet0 );
      facet1 = dynamic_cast<CubitFacet *>( cfacet1 );
       
      edge_list.append( edge );
      facet_list.append( facet0 );
      facet_list.append( facet1 );
      
      nflip++;
    }
    edge_list.step();
  }
  edge_list.remove_all_with_value(NULL);
  assert(edge_list.size() == nedges);
  if (nflip > 0)
  {
    fsurf_ptr->get_eval_tool()->replace_facets(facet_list);
  }
  return CUBIT_SUCCESS;
}


//================================================================================
// Description: smooth the facets on a surface
// Author     : Steve Owen
// Date       : 8/28/2003
// Note: free_laplacian = TRUE will not project the node back to a tangent plane
//       free_laplacian = FALSE will project to local tangent plane
//================================================================================
CubitStatus OCCModifyEngine::smooth_facets( RefFace *ref_face_ptr, int niter,
                                              CubitBoolean free_laplacian )
{
  CubitStatus rv = CUBIT_SUCCESS;

  // This must be a facet-based surface - return now if it isn't

  Surface *surf_ptr = ref_face_ptr->get_surface_ptr();
  OCCSurface *fsurf_ptr = dynamic_cast<OCCSurface *>(surf_ptr);
  if (fsurf_ptr == NULL)
  {
    PRINT_ERROR("Couldn't smooth facets on surface %d. It is not a facet-based geometry.\n",
                 ref_face_ptr->id());
    return CUBIT_FAILURE;
  }

  CubitPoint *point_ptr;
  DLIList<CubitPoint *>point_list;

  // get the points from the facet-basec surface

  fsurf_ptr->get_my_points(point_list);

  // compute normals if necessary  

  int ii, jj;
  for (ii=0; ii<point_list.size(); ii++)
  {
    point_ptr = point_list.get_and_step();
    point_ptr->normal();
    point_ptr->marked(0);
  }

  // get the points on the facet curves associated with this surface

  DLIList<CubitPoint *>cpoint_list;
  OCCCurve *fcurv_ptr;
  Curve *curv_ptr;

  int idx;
  DLIList<Curve *> curves;
  fsurf_ptr->curves(curves);
  OCCCurve **fcurve_array = new OCCCurve * [curves.size()];
  for (ii=0; ii<curves.size(); ii++)
  {
    idx = ii+1;
    curv_ptr = curves.get_and_step();
    fcurv_ptr = dynamic_cast<OCCCurve *>(curv_ptr);
    if (!fcurv_ptr)
    {
      PRINT_ERROR("Couldn't smooth facets on surface %d. At least one of it curves is not a facet-based geometry.\n",
                   ref_face_ptr->id());
      return CUBIT_FAILURE;
    }

    // mark the points on the curve with the index of the curve in the array
    // this will allow us to project the point back to the owning curve

    cpoint_list.clean_out();
    fcurv_ptr->get_points(cpoint_list);
    for (jj=0; jj<cpoint_list.size(); jj++)
    {
      point_ptr = cpoint_list.get_and_step();
      point_ptr->marked(  idx  );
    }
    fcurve_array[idx-1] = fcurv_ptr;

    // mark the vertices with a negative value so we don't touch them

    DLIList<Point *> verts;
    fcurv_ptr->points(verts);
    for(jj=0; jj<verts.size(); jj++)
    {
      Point *vert = verts.get_and_step();
      OCCPoint *fvert = dynamic_cast<OCCPoint *> (vert);
      fvert->get_cubit_point()->marked(-1);
    }
  }
  
  // laplacian smoothing

  DLIList<CubitPoint *> adj_points;
  int iter, nadj;
  CubitPoint *adj_point_ptr;
  CubitVector new_location; 

  for (iter = 0; iter < niter; iter++)
  {
    for (ii=0; ii<point_list.size(); ii++)
    {
      point_ptr = point_list.get_and_step();

      // don't smooth points at vertices

      if (point_ptr->marked() >= 0)
      {
        adj_points.clean_out();
        point_ptr->adjacent_points( adj_points );
        new_location.set(0.0, 0.0, 0.0);
        nadj=0;
        for (jj = 0; jj<adj_points.size(); jj++)
        {
          adj_point_ptr = adj_points.get_and_step();
          if (point_ptr->marked()) // is on a curve
          {
            // only use points on the same curve to influence the new location
            if (point_ptr->marked() == adj_point_ptr->marked()) 
            {
              new_location += adj_point_ptr->coordinates();
              nadj++;
            }
          }

          // interior nodes use all adjacent points
          else
          {
            new_location += adj_point_ptr->coordinates();
            nadj++;
          }
        }
        new_location /= nadj;

        // project it to a tangent plane defined at the point.
        // except if we are on a curve or we're using free-laplacian option
        if (!free_laplacian && point_ptr->marked() == 0)
          new_location = point_ptr->project_to_tangent_plane( new_location );

        // for points on a curve project to the curve definition
        if(point_ptr->marked() != 0)
        {
          const CubitVector temp = new_location;
          fcurve_array[point_ptr->marked()-1]->closest_point( temp, 
                                                            new_location );
        }
      
        if (point_ptr->check_inverted_facets(new_location) == CUBIT_SUCCESS)
        {
          point_ptr->set(new_location);
        }
      }
    }

    // update the normals

    for (ii=0; ii<point_list.size(); ii++)
    {
      point_ptr = point_list.get_and_step();
      point_ptr->compute_avg_normal();
    }
  }

  return rv;
}

//================================================================================
// Description: create a shell that is the normal offset of the given shell
// Author     : Steve Owen
// Date       : 8/28/2003
//================================================================================
CubitStatus OCCModifyEngine::create_shell_offset( BodySM *bodysm_ptr, 
                                                    BodySM *&new_bodysm, double offset )
{
  // check that this is a facet body

  OCCBody *fbody_ptr = dynamic_cast<OCCBody *>(bodysm_ptr);
  if (fbody_ptr == NULL)
  {
    PRINT_ERROR("Body is not a facet-based geometry.  This command is intended for use with facetted models\n");
    return CUBIT_FAILURE;
  }

  // The facet and point lists

  DLIList<CubitPoint *> point_list;
  DLIList<CubitFacet *> facet_list;
  
  // get the surfaces.  Check that they are also all facet bodies

  int ii;
  Surface *surf_ptr;
  OCCSurface *fsurf_ptr;
  DLIList<Surface *> surface_list;
  fbody_ptr->surfaces( surface_list );
  for (ii=0; ii<surface_list.size(); ii++)
  {
    surf_ptr = surface_list.get_and_step();
    fsurf_ptr = dynamic_cast<OCCSurface *>(surf_ptr);
    if (fsurf_ptr == NULL)
    {
      PRINT_ERROR("At least one of the surfaces in Body is not facet-based geometry.\n"  
                  "This command is intended for use with facetted models\n");
      return CUBIT_FAILURE;
    }
    fsurf_ptr->get_my_facets(facet_list, point_list);
  }

  // create an array of points so we can reference them quickly when generating new facets

  CubitPoint **points = new CubitPoint * [point_list.size()];

  // make a unique set of points

  int id = 0;
  CubitPoint *point_ptr;
  for (ii=0; ii < point_list.size(); ii++)
  {
    point_list.get_and_step()->marked(1);
  }
  for (ii=0; ii<point_list.size(); ii++)
  {
    point_ptr = point_list.get_and_step();
    if (point_ptr->marked() == 1)
    {
      point_ptr->marked(0);
      points[id] = point_ptr;
      point_ptr->set_id( id );
      id ++;
    }
  }
  int npoints = id;

  // get the curves and mark any point on a curve as non-smoothable.

  int jj;
  Curve *curv_ptr;
  OCCCurve *fcurv_ptr;
  DLIList<Curve *> curve_list;
  fbody_ptr->curves( curve_list );
  for(ii=0; ii<curve_list.size(); ii++)
  {
    curv_ptr = curve_list.get_and_step();
    fcurv_ptr = dynamic_cast<OCCCurve *>(curv_ptr);
    if (fcurv_ptr == NULL)
    {
      PRINT_ERROR("At least one of the curves in Body is not a facet-based geometry.\n"  
                  "This command is intended for use with facetted models\n"); 
      return CUBIT_FAILURE;
    }
    point_list.clean_out();
    fcurv_ptr->get_points( point_list );
    for(jj=0; jj<point_list.size(); jj++)
    {
      point_ptr = point_list.get_and_step();
      point_ptr->marked(1);
    }
  }

  // smooth the normals to reduce chance of overlap on concave regions

  CubitVector new_normal, normal;
  const int niter = 3;
  int iter, nadj;
  CubitPoint *adj_point_ptr;
  DLIList<CubitPoint *> adj_points;

  for (iter = 0; iter < niter; iter++)
  {
    for (ii=0; ii<npoints; ii++)
    {
      point_ptr = points[ii];

      // don't smooth points on curves

      if (point_ptr->marked() != 1)
      {
        adj_points.clean_out();
        point_ptr->adjacent_points( adj_points );
        new_normal.set(0.0, 0.0, 0.0);
        nadj = adj_points.size();
        for (jj = 0; jj<nadj; jj++)
        {
          adj_point_ptr = adj_points.get_and_step();
          new_normal += adj_point_ptr->normal();
        }
        new_normal.normalize();
        point_ptr->normal(new_normal);
      }
    }
  }

  // generate a new set of points offset from the original

  CubitPointData *new_point;
  CubitVector new_location;
  CubitPoint **new_points = new CubitPoint * [npoints];
  for (ii=0; ii<npoints; ii++)
  {
    point_ptr = points[ii];
    new_location = point_ptr->coordinates() + point_ptr->normal() * offset;
    new_point = new CubitPointData( new_location );
    new_point->set_id( ii );
    new_points[ii] = (CubitPoint *) new_point;
  }

  // generate triangles for the new shell

  DLIList<CubitFacet *> fix_list;
  DLIList<CubitFacet *> orig_list;
  double dot;
  double mindot = 1.0;
  CubitFacet *facet_ptr, *new_facet_ptr;
  DLIList<CubitFacet *> new_facet_list;
  CubitPoint *p0, *p1, *p2;
  CubitPoint *np0, *np1, *np2;
  for(ii=0; ii<facet_list.size(); ii++)
  {
    facet_ptr = facet_list.get_and_step();
    facet_ptr->points( p0, p1, p2 );
    np0 = new_points[p0->id()];
    np1 = new_points[p1->id()];
    np2 = new_points[p2->id()];
    new_facet_ptr = (CubitFacet *) new CubitFacetData( np0, np1, np2 );
    new_facet_list.append( new_facet_ptr );

    // compare normals to see if we inverted anything

    normal = facet_ptr->normal();
    new_normal = new_facet_ptr->normal();
    dot = normal % new_normal;
    if (dot < mindot) mindot = dot;
    if (dot < 0.707)
    {
      new_facet_ptr->marked(1);
      fix_list.append( new_facet_ptr );
      orig_list.append( facet_ptr );
    }
  }

  // go back and try to uninvert tris if needed

  int kk;
  DLIList<CubitFacet *> problem_list;
  if (fix_list.size() > 0)
  {
    CubitVector new_location; 
    for (iter = 0; iter < niter; iter++)
    {
      for (ii=0; ii<fix_list.size(); ii++)
      {
        facet_ptr = fix_list.get_and_step();
        point_list.clean_out();
        facet_ptr->points( point_list );
        for (jj=0; jj<point_list.size(); jj++)
        {
          point_ptr = point_list.get_and_step();

          // don't smooth marked points

          if (point_ptr->marked() != 1)
          {
            adj_points.clean_out();
            point_ptr->adjacent_points( adj_points );
            new_location.set(0.0, 0.0, 0.0);
            nadj = adj_points.size();
            for (kk = 0; kk<nadj; kk++)
            {
              adj_point_ptr = adj_points.get_and_step();
              new_location += adj_point_ptr->coordinates();
            }
            new_location /= nadj;
            point_ptr->set(new_location);
          }
        }
      }
    }

    // see if it worked

    fix_list.reset();
    orig_list.reset();
    for(ii=0; ii<fix_list.size(); ii++)
    {
      new_facet_ptr = fix_list.get_and_step();
      facet_ptr = orig_list.get_and_step();
      new_facet_ptr->update_plane();
      new_normal = new_facet_ptr->normal();
      normal = facet_ptr->normal();
      dot = normal % new_normal;
      if (dot < 0.707)
      {
        problem_list.append( new_facet_ptr );
      }
    }
  }

  // set the normals back the way they were and clean up

  for (ii=0; ii<npoints; ii++)
  {
    point_ptr = points[ii];
    point_ptr->compute_avg_normal();
  }
  delete [] points;
  delete [] new_points;

  // build the geometry

  const char *file_name = NULL;
  CubitBoolean use_feature_angle = CUBIT_TRUE;
  double feature_angle = 0.0;
  int interp_order = 0;
  CubitBoolean smooth_non_manifold = CUBIT_TRUE;
  CubitBoolean split_surfaces = CUBIT_FALSE;
  CubitBoolean stitch = CUBIT_FALSE;
  CubitBoolean improve = CUBIT_FALSE;
  DLIList<CubitQuadFacet *> quad_facet_list;
  DLIList<Surface *> new_surface_list;
  FacetFileFormat file_format = FROM_FACET_LIST;

  CubitStatus rv =
  OCCQueryEngine::instance()->import_facets( file_name, use_feature_angle, feature_angle, 0.0, 
                                   interp_order, smooth_non_manifold,
                                   split_surfaces, stitch, improve, quad_facet_list,
                                   new_facet_list, new_surface_list, file_format );
  
  if(problem_list.size() > 0)
  {
    PRINT_WARNING("Offset process on facet body, may have resulted in %d overalapping facets.\n"
                  "Check results carefully. Problem facet(s) are displayed in red.)\n", 
                  problem_list.size());
    for (ii=0; ii<problem_list.size(); ii++)
    {
      facet_ptr = problem_list.get_and_step();
      GfxDebug::draw_facet(facet_ptr, CUBIT_RED);
    }
    GfxDebug::flush();
  }

  if( new_surface_list.size() )
  {
    new_bodysm = new_surface_list.get()->bodysm();
    if( !new_bodysm ) 
    {
      PRINT_ERROR("Problem offsetting Body\n");
      return CUBIT_FAILURE;
    }
  }
  else
  {
    PRINT_ERROR("Problem offsetting Body\n");
    return CUBIT_FAILURE;
  }
  return rv;
}

CubitStatus OCCModifyEngine::webcut_with_sweep_surfaces(
                                 DLIList<BodySM*> &blank_bodies,
                                 DLIList<Surface*> &surfaces,
                                 const CubitVector& sweep_vector,
                                 bool sweep_perp, 
                                 bool through_all,
                                 bool outward,
                                 bool up_to_next, 
                                 Surface *stop_surf, 
                                 Curve *curve_to_sweep_along, 
                                 DLIList<BodySM*> &results_list,
                                 CubitBoolean imprint)
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return CUBIT_FAILURE;
}

CubitStatus OCCModifyEngine::webcut_with_sweep_curves(
                                 DLIList<BodySM*> &blank_bodies,
                                 DLIList<Curve*> &curves,
                                 const CubitVector& sweep_vector,
                                 bool through_all, 
                                 Surface *stop_surf, 
                                 Curve *curve_to_sweep_along, 
                                 DLIList<BodySM*> &results_list,
                                 CubitBoolean imprint)
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return CUBIT_FAILURE;
}

CubitStatus OCCModifyEngine::webcut_with_sweep_surfaces_rotated(
                                 DLIList<BodySM*> &blank_bodies,
                                 DLIList<Surface*> &surfaces,
                                 const CubitVector &point, 
                                 const CubitVector &sweep_axis, 
                                 double angle, 
                                 Surface *stop_surf, 
                                 bool up_to_next, 
                                 DLIList<BodySM*> &results_list,
                                 CubitBoolean imprint)
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return CUBIT_FAILURE;
}

CubitStatus OCCModifyEngine::webcut_with_sweep_curves_rotated(
                                 DLIList<BodySM*> &blank_bodies,
                                 DLIList<Curve*> &curves,
                                 const CubitVector &point, 
                                 const CubitVector &sweep_axis, 
                                 double angle, 
                                 Surface *stop_surf, 
                                 DLIList<BodySM*> &results_list,
                                 CubitBoolean imprint)
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return CUBIT_FAILURE;
}

CubitStatus OCCModifyEngine::scale( BodySM *&body, const CubitVector& factors )
{
  return OCCQueryEngine::instance()->scale( body, factors );
}

CubitStatus OCCModifyEngine::tolerant_imprint( DLIList<BodySM*> &bodies_in,
                                                 DLIList<BodySM*> &new_bodies )  const
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return CUBIT_FAILURE;
}

// EOF
