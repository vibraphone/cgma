#ifndef GSTREENODE_HPP
#define GSTREENODE_HPP

#include "DLIList.hpp"
#include "CubitPoint.hpp"

class GridSearchTreeNode {
private:
  DLIList <CubitPoint *> plist;

public:
  long icoord;
  long jcoord;
  long kcoord;

  GridSearchTreeNode (long i, long j, long k) {
    icoord=i;
    jcoord=j;
    kcoord=k;
  }
  ~GridSearchTreeNode() {}

  DLIList <CubitPoint *> get_list(){
    return plist;
  }

  void add(CubitPoint * data) {
    plist.append(data);
  }


class GSTNodeComparator {
public:	
  bool operator () (GridSearchTreeNode * a, GridSearchTreeNode * b) const
  {
    return ( a->icoord < b->icoord ) || (a->icoord==b->icoord && a->jcoord<b->jcoord) || (a->icoord == b->icoord && a->jcoord==b->jcoord && a->kcoord < b->kcoord );
  }
};	

};


#endif


