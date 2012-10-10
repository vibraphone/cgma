//-----------------------------------------------------------------------------
//
//  File: ImprintPointData.hpp
//
//  Purpose:  Child class of CubitPoint.  It is the Cubit-specific version of
//            of the CubitPoint.
//
//  Notes:    Note that this class contains data which is accessed
//            virtually from CubitPoint.  In most cases, you can create
//            a ImprintPointData and treat it as if it is a CubitPoint.
//            For example:
//               CubitPoint *cp = (CubitPoint *) new ImprintPointData(...);
//            There should be no reason to reference the CubitFacetPoint 
//            directly.  This is done to allow different data representations
//            of a point in addition to that used by Cubit.
//                        
//-----------------------------------------------------------------------------


#ifndef IMPRINTPOINTDATA_HPP
#define IMPRINTPOINTDATA_HPP

// Include for CubitBoolean
#include "CubitDefines.h"
#include "CubitVector.hpp"
#include "DLIList.hpp"
#include "CubitPoint.hpp"
#include "CastTo.hpp"
#include "RefEntity.hpp"
#include "RefVertex.hpp"
#include "MemoryManager.hpp"
#include "ImprintMatchData.hpp"
#include "ImprintLineSegment.hpp"

enum PointType { UNSET_POINT = -1,
                 VERTEX_ON_BOUNDARY_1,
                 VERTEX_ON_BOUNDARY_2,
                 CREATE_NEW_VERTEX,
                 VERTEX_ON_BOTH_BOUNDARIES,
                 ON_BOUNDARY_1,
                 ON_BOUNDARY_2,
                 ON_BOTH_BOUNDARIES,
                 ON_SURFACE_1,
                 ON_SURFACE_2};



class ImprintPointData : public CubitPoint
{
private:

  CubitVector coords;
  DLIList<CubitFacet*> *attachedFacets;
  RefEntity *myOwner;
  PointType myPointType;
  CubitBoolean ownerRefVertex;
  ImprintPointData *matchingPoint;
  CubitBoolean startPartition;
  CubitBoolean endPartition;
  DLIList<ImprintMatchData*> *pointMatches;
  CubitBoolean isMatched;
  ImprintLineSegment *prevSeg;
  ImprintLineSegment *nextSeg;
  
  int listLoopPos;
  int loopPos;
  int loopSize;
  
  int entityId;
  void initialize_imprint_data();
    //- Initialize just the data for this class.

  static MemoryManager memoryManager;
    //- memory management object

public:
  
  ImprintPointData(double x_val, double y_val, double z_val );
  ImprintPointData( const CubitVector &new_point );
  ImprintPointData( ImprintPointData *copy, CubitVector &new_coords );
  
  ~ImprintPointData();

  SetDynamicMemoryAllocation(memoryManager)
      //- class specific new and delete operators

  int id(){ return entityId;}

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

  void facets( DLIList<CubitFacet*> &facet_list)
    { if (attachedFacets) facet_list += *attachedFacets; }
  void edges( DLIList<CubitFacetEdge*> &edge_list);
  void points( DLIList<CubitPoint*> &point_list)
    { point_list.append( this ); }

  void compute_avg_normal();

    //-Functions just for imprint point data.
  void owner(RefEntity* my_owner)
    {
      myOwner = my_owner;
    }
  RefEntity* owner()
    {return myOwner;}
  CubitBoolean is_owner_vertex()
    {if (myPointType == CREATE_NEW_VERTEX ||
         myPointType == VERTEX_ON_BOUNDARY_1 ||
         myPointType == VERTEX_ON_BOUNDARY_2 ||
         myPointType == VERTEX_ON_BOTH_BOUNDARIES )
      return CUBIT_TRUE;
    else
      return CUBIT_FALSE;
    }
  void set_point_type(PointType type)
    {myPointType = type;}
  PointType get_point_type()
    {return myPointType;}
  CubitStatus set_matching_point(ImprintPointData *other);
  ImprintPointData* get_matching_point()
    {return matchingPoint;}
  void set_start_partition()
    {startPartition = CUBIT_TRUE;}
  void set_end_partition()
    {endPartition = CUBIT_TRUE;}
  CubitBoolean get_start_partition()
    {return startPartition;}
  CubitBoolean get_end_partition()
    {return endPartition;}
  void set_loop_pos( int ii, int jj, int loop_size )
  {
    listLoopPos = ii;
    loopPos = jj;
    loopSize = loop_size;
  }
  int get_list_loop_pos()
    {return listLoopPos;}
  int get_loop_pos()
    {return loopPos;}
  int get_loop_size()
    {return loopSize;}
  void add_match_data(ImprintMatchData *match_data);
  void get_match_list(DLIList <ImprintMatchData*> &match_data_list);
  void set_matched(ImprintMatchData *match_data);
  void set_unmatched();
  
  CubitBoolean is_matched()
    {return isMatched;}
  ImprintLineSegment* get_prev_seg()
    {return prevSeg;}
  void set_prev_seg(ImprintLineSegment *prev_seg)
    {prevSeg = prev_seg;}
  
  ImprintLineSegment* get_next_seg()
    {return nextSeg;}
  void set_next_seg(ImprintLineSegment *next_seg)
    {nextSeg = next_seg;}
  
      
};

#endif
