//-----------------------------------------------------------------------------
//
//  File: CubitFacetEdgeData.hpp
//
//  Purpose:  Child class of CubitFacetEdge.  It is the Cubit-specific version of
//            of the CubitFacetEdge.
//
//  Notes:    Note that this class contains data which is accessed
//            virtually from CubitFacetEdge.  In most cases, you can create
//            a CubitFacetEdgeData and treat it as if it is a CubitFacetEdge.
//            For example:
//               CubitFacetEdge *cfe = (CubitFacetEdge *) new CubitFacetEdgeData(...);
//            There should be no reason to reference the CubitFacetEdgeData 
//            directly.  This is done to allow different data representations
//            of a facet edge in addition to that used by Cubit.
//                        
//-----------------------------------------------------------------------------

#ifndef CUBITFACETEDGEDATA_HPP
#define CUBITFACETEDGEDATA_HPP

// Include for CubitBoolean
#include "CubitDefines.h"
#include "CubitPoint.hpp"
#include "DLIList.hpp"
#include "CubitFacetEdge.hpp"

class CubitVector;
class CubitFacet;
class CubitPoint;

class CubitFacetEdgeData : public CubitFacetEdge
{
private:
  CubitPoint *pointArray[2];
  DLIList<CubitFacet *> adjFacetList;
  int entityId;
    //- for debug tracking...
public:
  CubitFacetEdgeData(CubitPoint *p1, CubitPoint *p2);
  CubitFacetEdgeData(CubitPoint *p1, CubitPoint *p2,
                 CubitFacet *facet1, CubitFacet *facet2,
                 int edge_index1, int edge_index2);
    //- constructors
  ~CubitFacetEdgeData();
    //- destructor
  int id(){return entityId;}
  void set_id( int ent_id ) {entityId=ent_id;};
  
  CubitPoint *point( int index) { return pointArray[index]; };
  CubitPoint *start_node() { return point(0); }
  CubitPoint *end_node() { return point(1); }
    //- get one of its points
  void set_point( CubitPoint *pt, int index ) { pointArray[index] = pt; }
    //- set one of the points
  CubitFacet *adj_facet( int index ); 
    //- get the list of adjacent facets
  void facets(DLIList<CubitFacet*> &facet_list );

  void edges(DLIList<CubitFacetEdge*> &edge_list );
  void points(DLIList<CubitPoint*> &point_list );
  int num_adj_facets()
    { return adjFacetList.size(); }
  int number_tris() { return num_adj_facets(); }
  void marked (int my_flag ) { set_flag(my_flag); }
  int marked() { return get_flag(); }
  
  void add_facet(CubitFacet *facet_ptr){ adjFacetList.append(facet_ptr); }
  CubitStatus remove_facet( CubitFacet *facet_ptr );
  
  CubitStatus merge_edges( CubitFacetEdgeData* other_edge );
  
  void flip();

};


#endif


