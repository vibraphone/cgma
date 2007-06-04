//----------------------------------------------------------------------------
// Class: ImprintLineSegment
// Description: A class to represent a linesegment for boundary imprinting.
//
// Creator: David R. White
//----------------------------------------------------------------------------
#include "ImprintLineSegment.hpp"
#include "ImprintPointData.hpp"
#include "ImprintMatchData.hpp"
#include "CubitVector.hpp"
#include "CubitBox.hpp"
#include "RefEntity.hpp"
#include "DLIList.hpp"

ImprintLineSegment::ImprintLineSegment(ImprintPointData *start_point,
                                       ImprintPointData *end_point,
                                       RefEntity *owner_v)
{
  CubitVector min, max;
  CubitVector p1 = start_point->coordinates();
  CubitVector p2 = end_point->coordinates();
  if ( p1.x() < p2.x() ) {
     min.x(p1.x());
     max.x(p2.x());
  } 
  else {
    min.x(p2.x());
    max.x(p1.x());
  }
  if ( p1.y() < p2.y() ) {
    min.y(p1.y());
    max.y(p2.y());
  } 
  else {
    min.y(p2.y());
    max.y(p1.y());
  }
  if ( p1.z() < p2.z() ) {
    min.z(p1.z());
    max.z(p2.z());
  } 
  else {
    min.z(p2.z());
    max.z(p1.z());
  }
  boundingBox = CubitBox(min, max);
  startPoint = start_point;
  endPoint = end_point;
  nextSeg = NULL;
  prevSeg = NULL;
  myOwner = owner_v;
  inactiveFlag = CUBIT_FALSE;
  segMatches = NULL;
  startPoint->set_next_seg(this);
  endPoint->set_prev_seg(this);
}
ImprintLineSegment::~ImprintLineSegment()
{
  if ( segMatches != NULL )
    delete segMatches;
}
void ImprintLineSegment::add_match_data(ImprintMatchData* m_d)
{
  if ( segMatches == NULL )
    segMatches = new DLIList<ImprintMatchData*>;
  segMatches->append(m_d);
}
void ImprintLineSegment::get_match_data(DLIList<ImprintMatchData*> &match_data)
{
  if ( segMatches == NULL )
    return;
  else
    match_data += (*segMatches);
}
