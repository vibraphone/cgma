//-----------------------------------------------------------------------------
//
//  File: FaceterPointData.hpp
//
//  Purpose:  Child class of CubitPoint.  It is the Cubit-specific version of
//            of the CubitPoint.
//
//  Notes:    Note that this class contains data which is accessed
//            virtually from CubitPoint.  In most cases, you can create
//            a FaceterPointData and treat it as if it is a CubitPoint.
//            For example:
//               CubitPoint *cp = (CubitPoint *) new FaceterPointData(...);
//            There should be no reason to reference the CubitFacetPoint 
//            directly.  This is done to allow different data representations
//            of a point in addition to that used by Cubit.
//                        
//-----------------------------------------------------------------------------


#ifndef FACETERPOINTDATA_HPP
#define FACETERPOINTDATA_HPP

// Include for CubitBoolean
#include "CubitDefines.h"
#include "CubitVector.hpp"
#include "DLIList.hpp"
#include "MemoryManager.hpp"
#include "ToolDataUser.hpp"
#include "CubitMatrix.hpp"
#include "CubitPoint.hpp"
class CubitFacet;
class RefEntity;

class FaceterPointData : public CubitPoint
{
private:

  CubitVector coords;
  DLIList<CubitFacet*> *attachedFacets;
  RefEntity *myOwner;
  FaceterPointData *mNext, *mPrev;
  double myInteriorAngle;

  static MemoryManager memoryManager;
    //- memory management object

  int entityId;

public:
  
  FaceterPointData(double x_val, double y_val, double z_val );
  FaceterPointData( const CubitVector &new_point );
  ~FaceterPointData();

  SetDynamicMemoryAllocation(memoryManager)
      //- class specific new and delete operators

  int id(){ return entityId;}

    //The following functions are used for sorting FaceterPointData
    //based on their interior angle.
  static int sort_by_angle ( FaceterPointData *&pt_1,
                             FaceterPointData *&pt_2);
  void set_interior_angle(double angle)
    {myInteriorAngle = angle;}
  double get_interior_angle()
    {return myInteriorAngle;}
  void set_prev(FaceterPointData *prev)
    {mPrev = prev;}
  FaceterPointData* get_prev()
    {return mPrev;}
  void set_next(FaceterPointData *next)
    {mNext = next;}
  FaceterPointData* get_next()
    {return mNext;}
  
    

  double x(){return coords.x();}
  double y(){return coords.y();}
  double z(){return coords.z();}
  void set( const CubitVector &pos ) { coords = pos; }

  void marked(int marked){ markedFlag = marked;}
  int marked(){return markedFlag;}
    //- generic marker for efficient sorting.
      
  CubitVector coordinates() const { return coords; }
  void coordinates(double point_array[3]);
  void owner(RefEntity* my_owner)
    {myOwner = my_owner;}
  RefEntity* owner()
    {return myOwner;}
  void add_facet( CubitFacet *facet);
  void remove_facet( CubitFacet *facet );
  int num_adj_facets();

  void facets( DLIList<CubitFacet*> &facet_list)
    { if (attachedFacets) facet_list += *attachedFacets; }
  void edges( DLIList<CubitFacetEdge*> &edge_list);
  void points( DLIList<CubitPoint*> &point_list )
    { point_list.append( this ); }

  void compute_avg_normal();
};

#endif
