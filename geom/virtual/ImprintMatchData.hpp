//------------------------------------------------------------------------------
// Class: ImprintMatchData
// Description:  Holds data for the intersecting process of boundary imprinting.
// Owner: David R. White
// Date: 5/31/02
//------------------------------------------------------------------------------
#ifndef IMPRINTMATCHDATA_HPP
#define IMPRINTMATCHDATA_HPP

#include "CubitDefines.h"
#include "CubitVector.hpp"
//#include "ImprintPointData.hpp"
class ImprintPointData;
class ImprintLineSegment;

class ImprintMatchData
{
private:
  ImprintPointData *closestPoint;
  ImprintLineSegment *closestSeg;
  double closestDist;
  CubitVector *pointOnSeg;
  
public :
  ImprintMatchData ()
    {
      closestPoint = NULL;
      closestSeg = NULL;
      closestDist = CUBIT_DBL_MAX;
      pointOnSeg = NULL;
    }
  ~ImprintMatchData ()
    {
      if ( pointOnSeg )
        delete pointOnSeg;
    }
  void set_closest_point(ImprintPointData *point)
    {closestPoint = point;}
  ImprintPointData* get_closest_point()
    {return closestPoint;}

  void set_closest_seg(ImprintLineSegment *seg)
    {closestSeg = seg;}
  ImprintLineSegment* get_closest_seg()
    {return closestSeg;}

  void set_closest_dist(double dist)
    {closestDist = dist;}
  double get_closest_dist()
    {return closestDist;}

  void set_point_on(CubitVector &point_on)
    {
      if ( pointOnSeg == NULL )
        pointOnSeg = new CubitVector(point_on);
      else
        pointOnSeg->set(point_on);
    }
  CubitVector* get_point_on()
    {
      return pointOnSeg;
    }
  
};

#endif
