//-----------------------------------------------------------------------------
//
//  File: CubitPointData.hpp
//
//  Purpose:  Child class of CubitPoint.  It is the Cubit-specific version of
//            of the CubitPoint.
//
//  Notes:    Note that this class contains data which is accessed
//            virtually from CubitPoint.  In most cases, you can create
//            a CubitPointData and treat it as if it is a CubitPoint.
//            For example:
//               CubitPoint *cp = (CubitPoint *) new CubitPointData(...);
//            There should be no reason to reference the CubitFacetPoint 
//            directly.  This is done to allow different data representations
//            of a point in addition to that used by Cubit.
//                        
//-----------------------------------------------------------------------------


#ifndef CUBITPOINTDATA_HPP
#define CUBITPOINTDATA_HPP

// Include for CubitBoolean
#include "CubitDefines.h"
#include "CubitVector.hpp"
#include "DLIList.hpp"
#include "MemoryManager.hpp"
#include "ToolDataUser.hpp"
#include "CubitMatrix.hpp"
#include "CubitPoint.hpp"
class CubitFacet;

class CubitPointData : public CubitPoint
{
private:

  CubitVector coords;
  DLIList<CubitFacet*> *attachedFacets;

  static MemoryManager memoryManager;
    //- memory management object

  int entityId;

public:
  
  CubitPointData(double x_val, double y_val, double z_val );
  CubitPointData(double x_val, double y_val, double z_val,int *);
  CubitPointData( const CubitVector &new_point );
  ~CubitPointData();

  SetDynamicMemoryAllocation(memoryManager)
      //- class specific new and delete operators

  int id(){ return entityId;}
  void set_id( int new_id ) { entityId = new_id; }

  double x(){return coords.x();}
  double y(){return coords.y();}
  double z(){return coords.z();}
  void set( const CubitVector &pos ) { coords = pos; }

  void marked(int marked){ markedFlag = marked;}
  int marked(){return markedFlag;}
    //- generic marker for efficient sorting.
      
  CubitVector coordinates() const { return coords; }
  void coordinates(double point_array[3]);

  void add_facet( CubitFacet *facet);
  void remove_facet( CubitFacet *facet );
  int num_adj_facets();

  void facets( DLIList<CubitFacet*> &facet_list )
    { if (attachedFacets) facet_list += *attachedFacets; }
  void tris( DLIList<CubitFacet*> &facet_list ) 
		{ facets(facet_list); }
  void edges( DLIList<CubitFacetEdge*> &edge_list );
  void points( DLIList<CubitPoint*> &point_list )
    { point_list.append( this ); }

  void compute_avg_normal();
  CubitStatus merge_points( CubitPoint *dead_point, CubitBoolean keep_point = CUBIT_FALSE );
  CubitStatus collapse_edge( CubitPointData *dead_point );
};

#endif


