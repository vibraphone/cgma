//-------------------------------------------------------------------------
// Class:       TDChordal
// Description: Support for Chordal Axis.
// Author:      jitken
// Date:        1/20/2002
//-------------------------------------------------------------------------

#ifndef TD_CHORDAL_HPP
#define TD_CHORDAL_HPP

#include "ToolData.hpp"
#include "CubitVector.hpp"
#include "MemoryManager.hpp"
#include "CastTo.hpp"
#include "DLIList.hpp"

enum TriType {UNDEFINED, JUNCTION, SLEEVE, TERMINATED, DISCARDED};

class CubitPoint;

class TDChordal : public virtual ToolData
{
private:
  //bool computedMidpoint[3];
  bool boundaryEdge[3];
  //CubitPoint *midPoint[3];
  TriType triGenre;
  int numOfBoundaryEdges;
  bool visited;
 
public:

  TDChordal();
  virtual ~TDChordal();
  //-constructor and destructor

  static int is_chordal(const ToolData* td)
  { return ((dynamic_cast<const TDChordal*> (td)) != NULL); }

  CubitStatus flag_boundary_edge(int index);

  void set_tritype(TriType new_type){ triGenre = new_type; }
  
  TriType get_tritype(){ return triGenre; }

  CubitStatus determine_tritype();

  void get_non_boundary_edges(DLIList <int> &edge_index);
  
  void mark_visited(){ visited = TRUE; }

  void unmark_visited(){ visited = FALSE; }

  bool get_visited(){ return visited; }

 //  CubitPoint *get_midpoint(int index);
  
};

inline CubitStatus TDChordal::flag_boundary_edge(int index){ 
  assert(index >=0 && index < 3);
  
  //should only be flag once.
  if(boundaryEdge[index])
    return CUBIT_FAILURE;
  
  boundaryEdge[index] = true;
  numOfBoundaryEdges++;
  
  return CUBIT_SUCCESS;
}

#endif 

