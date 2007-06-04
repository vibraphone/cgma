//-----------------------------------------------------------------------------
//
//  File: CubitPoint.hpp
//
//  Purpose:  Point Class used for mesh based geometry and other tools.  It is
//            primarily used in conjuction with CubitFacet and CubitFacetEdge 
//            classes.
//
//  Notes:    Note that this class does not contain any private data.
//            All data should be defined within the child classes inherited
//            from this class.  The current Cubit data class that inherits
//            from CubitPoint is CubitPointData.  This is done so that 
//            other applications using CubitPoints can use their own
//            point data, but take advantage of the CGM/Cubit functionality.
//            Please do not add private data to this class; instead add the
//            data to the children and access through virtual functions.
// 
//            Do not create a CubitPoint directly.  For example, don't do: 
//              CubitPoint *cp = new CubitPoint(...);
//            You should instead create the appropriate child class, and 
//            cast it to a CubitPoint for use.  For example:
//              CubitPoint *cp = (CubitPoint *) new CubitPointData(...);
//                        
//-----------------------------------------------------------------------------

#ifndef CUBITPOINT_HPP
#define CUBITPOINT_HPP

// Include for CubitBoolean
#include "CubitDefines.h"
#include "CubitVector.hpp"
#include "DLIList.hpp"
#include "MemoryManager.hpp"
#include "ToolDataUser.hpp"
#include "CubitMatrix.hpp"
#include "FacetEntity.hpp"
class CubitFacet;
class CubitQuadFacet;
class CubitFacetEdge;
class CubitPoint;
class CubitTransformMatrix;
class Surface;
class CubitBox;

class CubitPoint : public FacetEntity
{
private:

  static double boxTol;

protected:

  int markedFlag;
    //-generic marking flag.

  CubitVector *surfNormal;
    //- The normalized surface normal (if required)

  double dCoef;
    //- D coeficient in tangent plane equation Ax+By+Cz+D=0
    //- where A,B,C is surfNormal

  double uVal, vVal, sizeVal;
    //- parametric location of point

  CubitVector *surfU, *surfV;
    //- Orthogonal surface tangent vectors at the point

  double *coefVector;
    //- coefficients for quadric surface approximation

  IttyBit isFeature;
    //- set if this point is a feature

public:

  CubitPoint();
  virtual ~CubitPoint();

  virtual int id() = 0; 
  virtual void set_id( int ) {};

  virtual double x() = 0;
  virtual double y() = 0;
  virtual double z() = 0;
  virtual void set( const CubitVector &pos ) = 0;

  virtual void marked(int marked) {markedFlag = marked;};
  virtual int marked() {return markedFlag;};
      
  virtual CubitVector coordinates() const = 0;
  virtual void coordinates(double point_array[3]) = 0;
  
  virtual void add_facet( CubitFacet *facet) = 0;
  virtual void remove_facet( CubitFacet *facet) = 0;
  virtual int num_adj_facets() = 0;

  virtual void facets( DLIList<CubitFacet*> &facet_list ) = 0;
  virtual void edges( DLIList<CubitFacetEdge*> &edge_list ) = 0;
  virtual void points( DLIList<CubitPoint*> &point_list ) = 0;

  virtual void tris( DLIList<CubitFacet*> &facet_list ) { facets(facet_list); }

  virtual void normal( CubitVector &surf_norm );
  virtual CubitVector normal();
  virtual CubitVector* normal_ptr();

  virtual void d_coef( const double d_coefficient ) {dCoef = d_coefficient;};
  virtual double d_coef() {return dCoef;};

  virtual double u() { return uVal; };
  virtual double v() { return vVal; };
  virtual double size() {return sizeVal; };
  virtual void set_uv( double u, double v ) { uVal = u; vVal = v; };
  virtual void set_uvs(double u, double v, double s) {uVal = u; vVal = v; sizeVal = s; }
    //- get and set the u-v coordinates (careful with internal boundaries - see TDFacetBoundarPoint)
  virtual CubitStatus get_uv( CubitFacet *facet, double &u, double &v );
  virtual CubitStatus get_uvs( CubitFacet *facet, double &u, double &v, double &s);
    //- return the u-v coords with respect to the surface that the facet is on
  
  virtual void du( CubitVector &duvec );
  virtual CubitVector du();

  virtual void dv( CubitVector &dvvec );
  virtual CubitVector dv();

  virtual double *coef_vector( );
  virtual void coef_vector( const CubitMatrix& coef );

  virtual CubitStatus merge_points( CubitPoint *cp, CubitBoolean keep_point = CUBIT_FALSE );

  void shared_facets( CubitPoint* other_point, 
                      CubitFacet*& facet1, 
                      CubitFacet*& facet2 );

  void shared_facets( CubitPoint* other_point, 
                      DLIList<CubitFacet*>& result_list );

  CubitFacetEdge *shared_edge( CubitPoint* other_point );

  
  void adjacent_points( CubitPoint **adj_points,
                        int &num_adj_points );
  void adjacent_points( DLIList<CubitPoint*>& result_list );
    //- return array of points sharing adjacent facets

  CubitBox bounding_box();
    //-return a box around this point

  void facets_on_surf( int surf_id, DLIList<CubitFacet *> &facet_list,
                       CubitBoolean &on_internal_boundary );
    //- return facets adjacent this point that are on the given surface
    //- surf_id is the FacetEvalTool ToolID

  CubitVector normal( CubitFacet *facet_ptr );
  CubitVector normal( CubitQuadFacet *qfacet_ptr );
  CubitVector normal( CubitFacetEdge *edge_ptr );
    //- return normal with respect to the facet and edge

  CubitVector tangent( CubitFacetEdge *edge_ptr, double mindot );
    //- return tangent with respect to the edge
  void next_feature_edges( CubitFacetEdge *this_edge_ptr,
                           DLIList <CubitFacetEdge *> feature_edge_list );
    //- given a facet boundary edge and this point, get a list
    //- of the next fetaure edges at this point

  CubitVector project_to_tangent_plane( CubitVector &pt );
    //- Project a point to the tangent plane defined at the CubitPoint

  void transform_to_local( CubitVector &glob_vec,
                           CubitVector &loc_vec );
  void transform_to_global( CubitVector &loc_vec,
                            CubitVector &glob_vec );
    //- transform a point between local and global system

  void define_tangent_vectors();
    //- set up the tangent vectors at the point for quadratic interpolation

  void get_parents( DLIList<FacetEntity*> &facet_list );

  void debug_draw( int color = -1, int flush = 1, int draw_uv=0 );
    // debug drawing

  void compute_avg_normal();
    //- compute the avg normal at this point based on adjacent facets

  CubitFacetEdge *get_edge( CubitPoint *other_point );
    //- return the edge between two points if one exists

  void transform(CubitTransformMatrix &tfmat);
  void rotate_normal(CubitTransformMatrix &rotmat);
    //- apply transformations to the point
    
  CubitStatus check_inverted_facets( const CubitVector& new_position );
    //- Check if moving this CubitPoint to the passed position
    //- will invert any adjacent facets.  Returns true if no facets
    //- will be inverted, false if one or more facets will be 
    //- inerted when the point is moved.

  void set_as_feature() { isFeature = 1; }
  CubitBoolean is_feature( ){return (isFeature ? CUBIT_TRUE : CUBIT_FALSE); }
    // set and get the isFeature bit

  static void set_box_tol( double tol ) {boxTol = tol;}
};

inline void CubitPoint::coef_vector(const CubitMatrix &input_matrix) 
{
  if (!coefVector) coefVector = new double[5];
  for (int i=0; i<5; i++) {
    coefVector[i] = input_matrix.get(i,0);
  }
}

inline void CubitPoint::normal( CubitVector &surf_norm )
{
  if(!surfNormal) surfNormal = new CubitVector (surf_norm);
  else *surfNormal = surf_norm;
}

inline CubitVector CubitPoint::normal()
{
  if (surfNormal==NULL) compute_avg_normal();
  return *surfNormal;
}

inline CubitVector* CubitPoint::normal_ptr()
{
  return surfNormal;
}


inline void CubitPoint::du( CubitVector &duvec )
{
  if(!surfU) surfU = new CubitVector (duvec);
  else *surfU = duvec;
}

inline void CubitPoint::dv( CubitVector &dvvec )
{
  if(!surfV) surfV = new CubitVector (dvvec);
  else *surfV = dvvec;
}

inline CubitVector CubitPoint::du()
{
  assert(surfU != NULL);
  return *surfU;
}

inline CubitVector CubitPoint::dv()
{
  assert(surfV != NULL);
  return *surfV;
}

inline double *CubitPoint::coef_vector( )
{
  assert (coefVector != NULL);
  return coefVector;
}

class CubitPointComparator {
public:	
  bool operator () (CubitPoint * a, CubitPoint * b) const
  {
    return ( a->x() < b->x() ) || (a->x()==b->x() && a->y()<b->y()) || (a->x() == b->x() && a->y()==b->y() && a->z() < b->z() );
  }
};


#endif







