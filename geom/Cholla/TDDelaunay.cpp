//-------------------------------------------------------------------------
// File:        TDDelaunay.cpp
// Description: Support for FacetorTool.  Maintains circumcircle
//              info at the triangles.  Do it up template style.
// Author:      chynes
// Date:        6/3/2002
//-------------------------------------------------------------------------
#ifdef INLINE_TEMPLATES
#define MY_INLINE inline
#else
#define MY_INLINE
#endif

#include "CubitMessage.hpp"
#include "TDDelaunay.hpp"
#ifndef SQR
#define SQR(x) ((x) * (x))
#endif

template< class TRIA, class NODE > MY_INLINE
TDDelaunay< TRIA, NODE >::TDDelaunay()
{
  mRadius = -1.0;
  visitFlag = INT_MIN;
  sortIndex = -1;
}

template< class TRIA, class NODE > MY_INLINE
TDDelaunay< TRIA, NODE >::~TDDelaunay()
{
}

//-------------------------------------------------------------------------
// Function:    reset
// Description: reset the circumradius stuff
// Note:        
// Author:      sjowen
// Date:        9/19/2003
//-------------------------------------------------------------------------
template< class TRIA, class NODE > MY_INLINE
void TDDelaunay< TRIA, NODE >::reset( )
{
  mRadius = -1.0;
  visitFlag = INT_MIN;
  sortIndex = -1;
}

//-------------------------------------------------------------------------
// Function:    circumcenter
// Description: return the circumcenter of a triangle in 3D.
//              computes and stores with the tooldata if it has not
//              yet been initialized
// Note:        This function was copied from NodeTri.cpp
// Author:      chynes
// Date:        6/3/2002
//-------------------------------------------------------------------------
template< class TRIA, class NODE > MY_INLINE
CubitVector& TDDelaunay< TRIA, NODE >::circumcenter( TRIA *tri_ptr )
{
  if (mRadius < 0.0) {

    NODE *corner[3];
    tri_ptr->tri_nodes( corner[0], corner[1], corner[2] );

      // Use coordinates relative to point `a' of the triangle. 
    CubitVector vec_ba(corner[0]->coordinates(), 
                       corner[1]->coordinates());
    CubitVector vec_ca(corner[0]->coordinates(), 
                       corner[2]->coordinates());

      // Squares of lengths of the edges incident to `a'. 
    double ba_length = vec_ba.length_squared();
    double ca_length = vec_ca.length_squared();
  
      // Cross product of these edges. 
      // (Take your chances with floating-point roundoff.)
    CubitVector cross_bc = vec_ba * vec_ca;

      // Calculate the denominator of the formulae. 
    double denominator = 0.5 / (cross_bc % cross_bc);
    assert(denominator != 0.0);

      // Calculate offset (from `a') of circumcenter. 
    mCenter  = (ba_length * vec_ca - ca_length * vec_ba) * cross_bc;
    mCenter *= denominator;

      // radius is length from point `a' to center
    mRadius = mCenter.length();

      // Add point `a' to get global coordinate of center
    mCenter += corner[0]->coordinates();
  }
  return mCenter;
}

//-------------------------------------------------------------------------
// Function:    circumcenter2d
// Description: same as circumcenter, but assumes 2d triangle (x,y,0)
// Author:      chynes
// Date:        6/3/2002
//-------------------------------------------------------------------------
template< class TRIA, class NODE > MY_INLINE
CubitVector& TDDelaunay< TRIA, NODE >::circumcenter2d( TRIA *tri_ptr )
{
  if (mRadius < 0.0) {

    NODE *corner[3];
    tri_ptr->tri_nodes( corner[0], corner[1], corner[2] );

    double A11 = corner[1]->coordinates().x() - 
                 corner[0]->coordinates().x();
    double A12 = corner[1]->coordinates().y() - 
                 corner[0]->coordinates().y();
    double A21 = corner[2]->coordinates().x() - 
                 corner[0]->coordinates().x();
    double A22 = corner[2]->coordinates().y() - 
                 corner[0]->coordinates().y();
    double det = A11*A22 - A21*A12;
    if(det == 0.0)
    {
      PRINT_ERROR("Co-linear points can not be used to determine circle.\n");
      mCenter.x(0.0);//assert(0);
      mCenter.y(0.0);
      mCenter.z(0.0);
      return mCenter;
    }
    
    double B1 = (A11*A11) + (A12*A12);
    double B2 = (A21*A21) + (A22*A22);
    det += det;
    double xc = (B1*A22 - B2*A12)/det;
    double yc = (B2*A11 - B1*A21)/det;
    mRadius = xc * xc + yc * yc;
    mCenter.x( corner[0]->coordinates().x() + xc );
    mCenter.y( corner[0]->coordinates().y() + yc );
  }
  return mCenter;
}

//-------------------------------------------------------------------------
// Function:    radius
// Description: return radius of circumcircle for 3D triangle.  Computes
//              and stores value if necessary 
// Author:      chynes
// Date:        6/3/2002
//-------------------------------------------------------------------------
template< class TRIA, class NODE > MY_INLINE
double TDDelaunay< TRIA, NODE >::radius( TRIA *tri_ptr )
{
  if (mRadius < 0.0) {
    this->circumcenter( tri_ptr );
  }
  return mRadius;
}

//-------------------------------------------------------------------------
// Function:    radius2d
// Description: same as radius, but assumes 2d triangle (x,y,0)
// Author:      chynes
// Date:        6/3/2002
//-------------------------------------------------------------------------
template< class TRIA, class NODE > MY_INLINE
double TDDelaunay< TRIA, NODE >::radius2d( TRIA *tri_ptr )
{
  if (mRadius < 0.0) {
    this->circumcenter2d( tri_ptr );
  }
  return mRadius;
}

//-------------------------------------------------------------------------
// Function:    circumsphere
// Description: define radius (squared) and circumsphere center for a tet
// Author:      sjowen
// Date:        8/4/2003
//-------------------------------------------------------------------------
template< class TRIA, class NODE > MY_INLINE 
int TDDelaunay< TRIA, NODE >:: 
circumsphere( TRIA *tet_ptr, CubitVector &center, double &radsq )
{
  if (mRadius != -1)
  {
    radsq = mRadius;
    center = mCenter;
    return CUBIT_SUCCESS;
  }
  double reltol = DBL_EPSILON * 100.0;

  NODE *na, *nb, *nc, *nd;
  tet_ptr->tet_nodes(nc, nb, na, nd);
  CubitVector a = na->coordinates();
  CubitVector b = nb->coordinates();
  CubitVector c = nc->coordinates();
  CubitVector d = nd->coordinates();

  CubitVector da = a - d;
  CubitVector db = b - d;
  CubitVector dc = c - d;

  double rhsa = 0.5*(SQR(da.x()) + SQR(da.y()) + SQR(da.z()));
  double rhsb = 0.5*(SQR(db.x()) + SQR(db.y()) + SQR(db.z()));
  double rhsc = 0.5*(SQR(dc.x()) + SQR(dc.y()) + SQR(dc.z()));

  double cpa = db.y()*dc.z() - dc.y()*db.z();
  double cpb = dc.y()*da.z() - da.y()*dc.z();
  double cpc = da.y()*db.z() - db.y()*da.z();

  double det = da.x()*cpa + db.x()*cpb + dc.x()*cpc;

  double xmax = CUBIT_MAX(fabs(a.x()),fabs(b.x()));
         xmax = CUBIT_MAX(xmax,fabs(c.x())); 
         xmax = CUBIT_MAX(xmax,fabs(d.x()));
  double ymax = CUBIT_MAX(fabs(a.y()),fabs(b.y()));
         ymax = CUBIT_MAX(ymax,fabs(c.y())); 
         ymax = CUBIT_MAX(ymax,fabs(d.y()));
  double zmax = CUBIT_MAX(fabs(a.z()),fabs(b.z()));
         zmax = CUBIT_MAX(zmax,fabs(c.z())); 
         zmax = CUBIT_MAX(zmax,fabs(d.z()));
  double tolabs = reltol*xmax*ymax*zmax;
  if (fabs(det) <= tolabs) {
    radsq = mRadius = -1.0; 
    return CUBIT_FAILURE;
  }
  center.x( (rhsa*cpa + rhsb*cpb + rhsc*cpc)/det );
  cpa = db.x()*rhsc - dc.x()*rhsb;
  cpb = dc.x()*rhsa - da.x()*rhsc;
  cpc = da.x()*rhsb - db.x()*rhsa;
  center.y( (da.z()*cpa + db.z()*cpb + dc.z()*cpc)/det );
  center.z( -(da.y()*cpa + db.y()*cpb + dc.y()*cpc)/det );
  radsq = SQR(center.x()) + SQR(center.y()) + SQR(center.z());
  center += d;

  mRadius = radsq;
  mCenter = center;

  return CUBIT_SUCCESS;
}






