
#ifndef GSTREE_HPP
#define GSTREE_HPP

#include "GridSearchTreeNode.hpp"

#include <typeinfo>
#if !defined(NT) && !defined(CANT_USE_STD)
using std::type_info;
#endif

#include <map>

typedef STD(map)< GridSearchTreeNode* , int, GridSearchTreeNode::GSTNodeComparator > gmap;

class GridSearchTree {

private:
  double epsilon;
  gmap nodemap;
  gmap::iterator pos;
public:

  GridSearchTree(double tolerance) {
    epsilon = tolerance;
  }

  ~GridSearchTree() {}
  

  CubitPoint * fix (CubitPoint * data);

};

#endif


