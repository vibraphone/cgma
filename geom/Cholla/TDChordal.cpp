//-------------------------------------------------------------------------
// Class:       TDChordal
// Description: Support for Chordal Axis.
// Author:      jitken
// Date:        1/20/2002
//-------------------------------------------------------------------------

#include "TDChordal.hpp"
#include "DLIList.hpp"

TDChordal::TDChordal(){
  int ii;
  for (ii = 0; ii < 3; ii++){
    boundaryEdge[ii] = false;
  }
  triGenre = UNDEFINED;
  numOfBoundaryEdges = 0;
  visited = FALSE;

}

TDChordal::~TDChordal(){

}


CubitStatus TDChordal::determine_tritype(){
  
  switch(numOfBoundaryEdges){
  case 0: 
    triGenre = JUNCTION;
    break;
  case 1:
    triGenre = SLEEVE;
    break;
  case 2:
    triGenre = TERMINATED;
    break;
  case 3:
    triGenre = DISCARDED;
    break;
  default:
    return CUBIT_FAILURE;
    break;
  }
  
  return CUBIT_SUCCESS;

}

void TDChordal::get_non_boundary_edges(DLIList <int> &edge_index){
  int ii = 0;
  for(ii = 0; ii<3; ii++){
    if(!boundaryEdge[ii]){
      edge_index.append(ii);
    }
  }

}

// bool TDChordal::computed_midpoint(int index){
//   assert(index >= 0 && index < 3);
//   return computedMidpoint[index];
// }

// CubitPoint *TDChordal::get_midpoint(int index){
//   assert(index >= 0 && index < 3);
//   return midPoint[index];
// }
