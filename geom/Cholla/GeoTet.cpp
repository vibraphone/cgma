//-------------------------------------------------------------------------
// Filename      : GeoTet.cpp
//
// Purpose       : Tet class used for geometry operations.  See also GeoNode.
//
// Description   : light-weight entity used for computational geometry 
//                 tools.  Don't load this up with extra stuff!
//
// Creator       : Steve Owen
//
// Creation Date : 8/13/2003
//
// Owner         : Steve Owen
//-------------------------------------------------------------------------

#include "GeoTet.hpp"
#include "GeoNode.hpp"
#include "GfxDebug.hpp"

#ifndef SQR
#define SQR(x) ((x) * (x))
#endif

//--------------------------------------------------------------------------
// Function: GeoTet
// Description: constructor
// Author: sjowen
// Date: 8/13/2003
//---------------------------------------------------------------------------
GeoTet::GeoTet( GeoNode *nodes[4] ) : isMarked(0), isInside(0)
{
  mNodes[0] = nodes[0];
  mNodes[1] = nodes[1];
  mNodes[2] = nodes[2];
  mNodes[3] = nodes[3];
}

//--------------------------------------------------------------------------
// Function: GeoTet
// Description: destructor
// Author: sjowen
// Date: 8/13/2003
//---------------------------------------------------------------------------
GeoTet::~GeoTet()
{
}

//--------------------------------------------------------------------------
// Function: tet_nodes
// Description: return the nodes on this tet
// Author: sjowen
// Date: 8/13/2003
//---------------------------------------------------------------------------
void GeoTet::tet_nodes( GeoNode *&na, GeoNode *&nb, GeoNode *&nc, GeoNode *&nd )
{
  na = mNodes[0];
  nb = mNodes[1];
  nc = mNodes[2];
  nd = mNodes[3];
}

//--------------------------------------------------------------------------
// Function: tet_face_nodes
// Description: return the on a face of this tet.
// Author: mlstate
// Date: 9/08/2004
//---------------------------------------------------------------------------
void GeoTet::tet_face_nodes( int face_index, GeoNode *&na, GeoNode *&nb, GeoNode *&nc )
{
    int node_face_idx[4][3] = { { 1, 2, 3 },
                                { 3, 2, 0 },
                                { 0, 1, 3 },
                                { 2, 1, 0 } };
    na = mNodes[node_face_idx[face_index][0]];
    nb = mNodes[node_face_idx[face_index][1]];
    nc = mNodes[node_face_idx[face_index][2]];

}

//--------------------------------------------------------------------------
// Function: get_connected_tet
// Description: return the adjacent tet at the face of face_indx.
// Author: mlstate
// Date: 9/08/2004
//---------------------------------------------------------------------------
GeoTet *GeoTet::get_connected_tet( int face_indx )
{
    GeoNode *n1, *n2, *n3;
    tet_face_nodes( face_indx, n1, n2, n3 );
    return get_connected_tet( n1, n2, n3 );
}

//--------------------------------------------------------------------------
// Function: get_connected_tet
// Description: return the adjacent tet at the face defined by the three nodes
// Author: sjowen
// Date: 8/13/2003
//---------------------------------------------------------------------------
GeoTet *GeoTet::get_connected_tet( GeoNode *na, GeoNode *nb, GeoNode *nc )
{
  GeoTet *adj_tet = NULL;
  int ii;

  // get the tets adjacent 
  DLIList<GeoTet *> *tet_list_ptr = na->tet_list();
  GeoTet *tet_ptr;
  for (ii=0; ii<tet_list_ptr->size() && !adj_tet; ii++)
  {
    tet_ptr = tet_list_ptr->get_and_step();
    if (tet_ptr != this)
    {
      if (tet_ptr->node_index( nb ) >= 0 &&
          tet_ptr->node_index( nc ) >= 0)
      {
        adj_tet = tet_ptr;
      }
    }
  }

  return adj_tet;
}

//--------------------------------------------------------------------------
// Function: node_index
// Description: return the index of the node in the tet
// Author: sjowen
// Date: 8/13/2003
//---------------------------------------------------------------------------
int GeoTet::node_index( GeoNode *node_ptr )
{
  int ii;
  for (ii=0; ii<4; ii++)
  {
    if (mNodes[ii] == node_ptr)
      return ii;
  }
  return -1;
}

//--------------------------------------------------------------------------
// Function: opposite_edge
// Description: return nodes on the opposite edge of the tet.  a_node and 
//              b_node must be on this tet.  c_node and d_node are in
//              no particular order
// Author: sjowen
// Date: 01/29/2004
//---------------------------------------------------------------------------
void GeoTet::opposite_edge( GeoNode *a_node, GeoNode *b_node, 
                            GeoNode *&c_node, GeoNode *&d_node )
{
  int ii;
  c_node = d_node = NULL;
  for (ii=0; ii<4; ii++)
  {
    if (mNodes[ii] != a_node && mNodes[ii] != b_node)
      if (!c_node)
        c_node = mNodes[ii];
      else if(!d_node)
        d_node = mNodes[ii];
      else
        assert(0);  // a_node or b_node are not on this tet
  }
}

//--------------------------------------------------------------------------
// Function: dgtet
// Description: global debug function to draw a GeoTet 
// Author: sjowen
// Date: 8/13/2003
//---------------------------------------------------------------------------
void dgtet( GeoTet *tet )
{
  GfxDebug::draw_geotet( tet, CUBIT_YELLOW );
  GfxDebug::flush();
}

//--------------------------------------------------------------------------
// Function: dgnode
// Description: global debug function to draw a GeoNode 
// Author: sjowen
// Date: 8/13/2003
//---------------------------------------------------------------------------
void dgnode( GeoNode *node )
{
  GfxDebug::draw_geonode( node, CUBIT_RED );
  GfxDebug::flush();
}

CubitStatus GeoTet::circumsphere( CubitVector &center, double *radius )
{
  double reltol = DBL_EPSILON * 100.0;

  CubitVector a = mNodes[0]->coordinates();
  CubitVector b = mNodes[1]->coordinates();
  CubitVector c = mNodes[2]->coordinates();
  CubitVector d = mNodes[3]->coordinates();

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
    return CUBIT_FAILURE;
  }
  center.x( (rhsa*cpa + rhsb*cpb + rhsc*cpc)/det );
  cpa = db.x()*rhsc - dc.x()*rhsb;
  cpb = dc.x()*rhsa - da.x()*rhsc;
  cpc = da.x()*rhsb - db.x()*rhsa;
  center.y( (da.z()*cpa + db.z()*cpb + dc.z()*cpc)/det );
  center.z( -(da.y()*cpa + db.y()*cpb + dc.y()*cpc)/det );
  center += d;

  if ( radius )
  {
      double radsq = SQR(center.x()) + SQR(center.y()) + SQR(center.z());
      *radius = sqrt( radsq );
  }

  return CUBIT_SUCCESS;
}

// EOF

