
#ifndef GSTREE_HPP
#define GSTREE_HPP

#include "GridSearchTreeNode.hpp"

#include <typeinfo>
#if !defined(WIN32)
using std::type_info;
#endif

#include <map>

typedef std::map< GridSearchTreeNode* , int, GridSearchTreeNode::GSTNodeComparator > gmap;

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


