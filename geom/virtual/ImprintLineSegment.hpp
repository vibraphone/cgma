//---------------------------------------------------------------------------
//  Class: ImprintLineSegment
//  Description: Utility class for boundary imprinting (ImprintBoundaryTool).
//               Stores a line segment and its bounding box.  Also holds on
//               to the owner of the edge.
//  Owner: David R. White
//  Creation Date: 4/15/2002
//---------------------------------------------------------------------------
#ifndef IMPRINTLINESEGMENT_HPP
#define IMPRINTLINESEGMENT_HPP

#include "MemoryManager.hpp"
template <class X> class DLIList;
#include "CubitBox.hpp"
class CubitVector;
class RefEntity;
class ImprintPointData;
class ImprintMatchData;

class ImprintLineSegment 
{
private:
  ImprintPointData *startPoint;
  ImprintPointData *endPoint;
  ImprintLineSegment *nextSeg, *prevSeg;
  CubitBox boundingBox;
  RefEntity *myOwner;
  DLIList <ImprintMatchData*> *segMatches;
  
  CubitBoolean inactiveFlag;
  static MemoryManager memoryManager;
    //- memory management object
public:
  ImprintLineSegment(ImprintPointData *start_point,
                     ImprintPointData *end_point,
                     RefEntity *owner_v);
    
  ~ImprintLineSegment();
  
  SetDynamicMemoryAllocation(memoryManager)
      //- class specific new and delete operators
  CubitBox& bounding_box()
    {return boundingBox;}
    //- Bounding box for the line segment.

  ImprintPointData *get_start()
    {return startPoint;}
  ImprintPointData *get_end()
    {return endPoint;}
  RefEntity* owner()
    {return myOwner;}
    //- Entity that this segment lies on.

  CubitBoolean get_inactive()
    {return inactiveFlag;}
  void set_inactive(CubitBoolean m)
    {inactiveFlag=m;}

  ImprintLineSegment* get_next()
    {return nextSeg;}
  ImprintLineSegment* get_prev()
    {return prevSeg;}

  void set_next(ImprintLineSegment *v)
    {nextSeg = v;}
  void set_prev(ImprintLineSegment *v)
    {prevSeg = v;}

  void add_match_data(ImprintMatchData* m_d);
  
  void get_match_data(DLIList<ImprintMatchData*> &match_data);
  
};
#endif
