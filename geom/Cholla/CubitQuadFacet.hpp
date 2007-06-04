#ifndef CUBITQUADFACET_HPP
#define CUBITQUADFACET_HPP

#include "CubitDefines.h"
#include "FacetEntity.hpp"
class CubitFacet;
class CubitVector;
class CubitPoint;

class CubitQuadFacet : public FacetEntity
{
private:
  void map_to_tri_system(double u, double v, int &tri_index,
                         double &A, double &B, double &C);
    //- map the u,v quad system into one of its triangle's area coordinate system
  CubitStatus eval_derivatives( double u, double v,
                                CubitVector &point,
                                CubitVector &du,
                                CubitVector &dv );
    //- evaluate derivatives on the quad

public:
  CubitQuadFacet( );
  virtual ~CubitQuadFacet();

  virtual int tri_to_quad_index( int tri_index, int pt_index ) = 0;
    // return the point index on the quad facet gicen one of its triangles
    // and the point index on the triangle (indecies start at 0)

  virtual CubitFacet *get_tri_facet( int index ) = 0;
    // return one of the underlying triangles (index should be 0 or 1)

  virtual CubitFacet *get_tri_facet_at_point( CubitPoint *point_ptr ) = 0;
    // return the underlying triangle of which the point_ptr is a vertex

  virtual void remove_tri_facets(  ) = 0;
    // used prior to delete if tri facets have already been deleted by
    // another mechanism.

  virtual void points (CubitPoint *points[4]) = 0;
    // return the four vertices of the quad facet

  virtual CubitPoint *point( int index ) = 0;
    // return the specified point
  virtual CubitFacetEdge *edge( int index ) = 0;
    // return the specified edge

  void get_control_points( double *ctrl_pts );
  void set_control_points( double *ctrl_pts );
    // get and set the control points

  CubitStatus evaluate( double u, double v, 
                        CubitVector *eval_point,
                        CubitVector *eval_normal = NULL,
                        CubitVector *eval_du = NULL,
                        CubitVector *eval_dv = NULL );
    // evaluate u, v location and normal on facet. Pass in NULL if 
    // point or normal are not to be evaluated

  virtual void debug_draw(int color=-1, int flush_it = 1, int draw_uv=0);
    // draw the underlying tri facets

  virtual void get_parents(DLIList<FacetEntity *> &){}; 
    // dummy for this class 
  virtual void points(DLIList<CubitPoint*> &point_list ); 
  virtual void facets(DLIList<CubitFacet*> &facet_list );
  virtual void edges(DLIList<CubitFacetEdge*> &edge_list );

  CubitStatus init_patch(  );
    // computes the interior control points for the quad facet and
    // stores them with the quad facet.  Assumes edge control points
    // on the four boundary edges have already been computed

};

#endif


