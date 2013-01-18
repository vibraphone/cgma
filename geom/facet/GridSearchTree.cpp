#include <math.h>
#include "GridSearchTree.hpp"

double dist (CubitPoint * a, CubitPoint * b){
  return sqrt((a->x() - b->x())*(a->x() - b->x()) + (a->y() - b->y())*(a->y() - b->y()) + (a->z() - b->z())*(a->z() - b->z()) );
}

CubitPoint * GridSearchTree::fix(CubitPoint * data) {
  
  long i, j, k;
  double x, y, z;

  x=data->x();
  y=data->y();
  z=data->z();
  
  // get the coordinates of the box containing the point
  i=(long)(x/(2*epsilon));
  j=(long)(y/(2*epsilon));
  k=(long)(z/(2*epsilon));
  

  if (x<0) i--;
  if (y<0) j--;
  if (z<0) k--;


  int ofi;
  int ofj;
  int ofk;
  int ii;

  // calculate the i, j, k offset for the seven neighboring boxes to be searched
  if (fabs(x-i*2*epsilon) < epsilon) ofi = -1; else ofi = 1;
  if (fabs(y-j*2*epsilon) < epsilon) ofj = -1; else ofj = 1;
  if (fabs(z-k*2*epsilon) < epsilon) ofk = -1; else ofk = 1;

  // mindist holds the distance to the current closest point
  double mindist=2*epsilon;
  // curdist holds the distance to the current point
  double curdist;
  // closest is the current closest point
  CubitPoint * closest = NULL;
  // curpoint is the current point
  CubitPoint * curpoint;

  // construct grid cell ( box ) 
  GridSearchTreeNode * curnode = new GridSearchTreeNode(i+ofi,j+ofj,k+ofk);
  // attempt to find it in the tree
  pos = nodemap.find(curnode);
  if (pos!=nodemap.end()) 
  {
    // if cell is in the tree search its list for close points
    DLIList<CubitPoint*> curlist = (*pos).first->get_list();
    for (ii= curlist.size(); ii>0; ii--) 
    {
      curpoint = curlist.get_and_step();
      curdist = dist(data, curpoint);
      if (curdist<mindist) 
      {
        closest = curpoint;
        mindist=curdist;
      }
    }
  }
  delete curnode;
  
  // construct grid cell ( box ) 
  curnode = new GridSearchTreeNode(i+ofi,j+ofj,k);
  // attempt to find it in the tree
  pos = nodemap.find(curnode);
  if (pos!=nodemap.end()) 
  {
    // if cell is in the tree search its list for close points
    DLIList<CubitPoint*> curlist = (*pos).first->get_list();
    for (ii= curlist.size(); ii>0; ii--) 
    {
      curpoint = curlist.get_and_step();
      curdist = dist(data, curpoint);
      if (curdist<mindist) 
      {
        closest = curpoint;
        mindist=curdist;
      }
    }
  }
  delete curnode;
  
  curnode = new GridSearchTreeNode(i+ofi,j,k+ofk);
  pos = nodemap.find(curnode);
  if (pos!=nodemap.end()) 
  {
    DLIList<CubitPoint*> curlist = (*pos).first->get_list();
    for (ii= curlist.size(); ii>0; ii--) 
    {
      curpoint = curlist.get_and_step();
      curdist = dist(data, curpoint);
      if (curdist<mindist) 
      {
        closest = curpoint;
        mindist=curdist;
      }
    }
  }
  delete curnode;
  
  curnode = new GridSearchTreeNode(i+ofi,j,k);
  pos = nodemap.find(curnode);
  if (pos!=nodemap.end()) 
  {
    DLIList<CubitPoint*> curlist = (*pos).first->get_list();
    for (ii= curlist.size(); ii>0; ii--) 
    {
      curpoint = curlist.get_and_step();
      curdist = dist(data, curpoint);
      if (curdist<mindist) 
      {
        closest = curpoint;
        mindist=curdist;
      }
    }
  }
  delete curnode;
  
  curnode = new GridSearchTreeNode(i,j+ofj,k+ofk);
  pos = nodemap.find(curnode);
  if (pos!=nodemap.end()) 
  {
    DLIList<CubitPoint*> curlist = (*pos).first->get_list();
    for (ii= curlist.size(); ii>0; ii--) 
    {
      curpoint = curlist.get_and_step();
      curdist = dist(data, curpoint);
      if (curdist<mindist) 
      {
        closest = curpoint;
        mindist=curdist;
      }
    }
  }
  delete curnode;
  
  curnode = new GridSearchTreeNode(i,j+ofj,k);
  pos = nodemap.find(curnode);
  if (pos!=nodemap.end()) 
  {
    DLIList<CubitPoint*> curlist = (*pos).first->get_list();
    for (ii= curlist.size(); ii>0; ii--) 
    {
      curpoint = curlist.get_and_step();
      curdist = dist(data, curpoint);
      if (curdist<mindist) 
      {
        closest = curpoint;
        mindist=curdist;
      }
    }
  }
  delete curnode;
  
  curnode = new GridSearchTreeNode(i,j,k+ofk);
  pos = nodemap.find(curnode);
  if (pos!=nodemap.end()) 
  {
    DLIList<CubitPoint*> curlist = (*pos).first->get_list();
    for (ii= curlist.size(); ii>0; ii--) 
    {
      curpoint = curlist.get_and_step();
      curdist = dist(data, curpoint);
      if (curdist<mindist) 
      {
        closest = curpoint;
        mindist=curdist;
      }
    }
  }
  delete curnode;
  

  curnode = new GridSearchTreeNode(i,j,k);
  pos = nodemap.find(curnode);
  if (pos!=nodemap.end()) 
  { 
    DLIList<CubitPoint*> curlist = (*pos).first->get_list();
    for (ii= curlist.size(); ii>0; ii--) 
    {
      curpoint = curlist.get_and_step();
      curdist = dist(data, curpoint);
      if (curdist<mindist) 
      {
        closest = curpoint;
        mindist=curdist;
      }
    }
  }

  // if closest point is within epsilon distance return the closest point
  if (mindist<=epsilon) 
  {
    delete curnode;
    return closest;
  }
  else 
  {          
    // add current point and cell to tree
    if (pos==nodemap.end()) 
    {
      curnode->add(data);
      nodemap.insert(gmap::value_type(curnode, 1));
    }
    else 
    {
      (*pos).first->add(data);
      delete curnode;
    }
    
    return data;
  }
}

