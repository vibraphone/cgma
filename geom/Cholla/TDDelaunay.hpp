//-------------------------------------------------------------------------
// Class:       TDDelaunay
// Description: Support for TriDelaunayTool.  Maintains circumcircle
//              info at the triangles.  Do it up template style
// Author:      chynes
// Date:        6/3/2002
//-------------------------------------------------------------------------

#ifndef TD_DELAUNAY_HPP
#define TD_DELAUNAY_HPP

#include "ToolData.hpp"
#include "CubitVector.hpp"
#include "MemoryManager.hpp"
#include "CastTo.hpp"

template< class TRIA, class NODE > 
class TDDelaunay : public virtual ToolData
{
private:

//  static MemoryManager memoryManager;
//   //- memory management object

  CubitVector mCenter;
  double mRadius;
  int visitFlag;
  int sortIndex;

public:

  TDDelaunay<TRIA, NODE>();
  virtual ~TDDelaunay<TRIA, NODE>();
   //-constructor and destructor

  static int is_delaunay(const ToolData* td)
     {return (CAST_TO(td, const TDDelaunay) != NULL);}

  void reset();
    // reset members to default

  CubitVector &circumcenter2d( TRIA *tri_ptr );
  double radius2d( TRIA *tri_ptr );
   //- compute radius and circumcircle info for a 2D triangle 

  CubitVector &circumcenter( TRIA *tri_ptr );
  double radius( TRIA *tri_ptr );
   //- compute radius and circumcircle info for a 2D triangle 

  int circumsphere( TRIA *tet_ptr, CubitVector &center, double &rad );
   //- compute the radius and circumsphere of a tet

  int visit_flag( ) {return visitFlag;};
  void visit_flag( int visit ) {visitFlag = visit;};
   //- get and set the visites flag

  int tri_sort_list( ) {return sortIndex;};
  void tri_sort_list( int index ) {sortIndex = index;};
   //- get and set the index of the sorting array used for 
   //- prioritizing which tris will be processed first

//  SetDynamicMemoryAllocation(memoryManager)
//   //- class specific new and delete operators
//
//  static void set_memory_allocation_increment(int increment = 0)
//                {memoryManager.set_memory_allocation_increment(increment);}
//   //- set block memory size increment
//
//  static void destroy_memory()
//                {memoryManager.destroy_memory();}
//   //- destroy all memory allocted to this object*/
};

#include "TDDelaunay.cpp"

#endif // TD_DELAUNAY_HPP

