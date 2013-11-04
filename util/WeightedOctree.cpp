//-------------------------------------------------------------------------
// Filename      : WOctree.cpp
//
// Purpose       : Class for O(ln n) search for nodes within a specified
//                 tolerance of a specified position.
//
// Creator       : Corey McBride 
//
// Creation Date : 02/18/10
//-------------------------------------------------------------------------

#include "WeightedOctree.hpp"
#include "CubitDefines.h"
#include "GeometryDefines.h"
#include <vector>
#include <algorithm>

const int DEFAULT_MAX_WEIGHT_PER_BOX = 100;
// The default value to use for the maximum weight of each box
// , if none is specified in the constructr.

//-------------------------------------------------------------------------
// Purpose       : Constructor
//
// Creator       : Corey McBride 
//
// Creation Date : 02/18/10
//-------------------------------------------------------------------------
template<class X, class C, class E>
WeightedOctree<X,C,E>::WeightedOctree(int max_weight_per_box,double tolerance)
{

  this->Tolerance = tolerance < GEOMETRY_RESABS ? GEOMETRY_RESABS : tolerance;
  this->ToleranceSquared=this->Tolerance*this->Tolerance;

  this->MaxBoxWeight= max_weight_per_box < 0 ? DEFAULT_MAX_WEIGHT_PER_BOX : max_weight_per_box;
 
  this->Root = new Cell();

}

//-------------------------------------------------------------------------
// Purpose       : Destructor
//
// Creator       : Corey McBride 
//
// Creation Date : 02/18/10
//-------------------------------------------------------------------------
template<class X, class C, class E>
WeightedOctree<X,C,E>::~WeightedOctree()
{
  // Release all memory
  delete this->Root;


  this->Root = 0;
}
template<class X, class C, class E>
bool WeightedOctree<X,C,E>::addNode(const X& node,
                                    std::vector<Cell*>& oldCells,
                                    std::vector<Cell*>& newCells,
                                    std::vector<Cell*>& modifiedCells)
{
  newCells.clear();
  oldCells.clear();
  modifiedCells.clear();

  Cell* currentCell=this->Root;

  CubitVector coords = E::coordinates(node);
  int weight=E::weight(node);
  while(!currentCell->leaf())
  {
    int x = coords.x() < currentCell->center().x() ? Cell::X_MIN : Cell::X_MAX;
    int y = coords.y() < currentCell->center().y() ? Cell::Y_MIN : Cell::Y_MAX;
    int z = coords.z() < currentCell->center().z() ? Cell::Z_MIN : Cell::Z_MAX;

    currentCell=currentCell->child(x|y|z);
  }


  currentCell->addNode(node);
  currentCell->addWeight(weight);

  this->NodeToCellMap[node]=currentCell;

  Cell* currentParent=currentCell->parent();
  while(currentParent)
  {
    currentParent->addWeight(weight);
    currentParent=currentParent->parent();
  }


  typename std::map<X,Cell*>::iterator nodeIter;
  if(currentCell->weight()>this->MaxBoxWeight && currentCell->nodeCount()>=2)
  {
    //Try to rebalance first.
    if(currentCell->parent())
    {
      if(this->rebalanceBranch(currentCell->parent(),oldCells,newCells))
      {
        for(unsigned int i=0;i<newCells.size();i++)
        {
          Cell* cell=newCells[i];
          std::vector<X> cellNodes;
          cell->appendNodes(cellNodes);

          for(unsigned int j=0;j<cellNodes.size();j++)
          {
            X& localNode=cellNodes[j];

            nodeIter=this->NodeToCellMap.find(localNode);
            if(nodeIter!=this->NodeToCellMap.end())
            {
              E::setLastCell(localNode,nodeIter->second);
            }
            this->NodeToCellMap[localNode]=cell;
          }
        }
        return true;
      }
    }

    //didn't rebalance so
    //Try to split
    if(this->splitCell(currentCell,newCells))
    {
      //The cell was split. We need to update the node map.
      for(unsigned int i=0;i<newCells.size();i++)
      {
        Cell* cell=newCells[i];
        std::vector<X> cellNodes;
        cell->appendNodes(cellNodes);

        for(unsigned int j=0;j<cellNodes.size();j++)
        {
          X& localNode=cellNodes[j];
          nodeIter=this->NodeToCellMap.find(localNode);
          if(nodeIter!=this->NodeToCellMap.end())
          {
            E::setLastCell(localNode,nodeIter->second);
          }
          this->NodeToCellMap[localNode]=cell;
        }
      }
      oldCells.push_back(currentCell);
    }
    else
    {
      //the cells wasn't split
      if(currentCell->nodeCount()==1)
      {
        newCells.push_back(currentCell);
      }
      else
      {

        modifiedCells.push_back(currentCell);
      }
    }
  }
  else
  {
    if(currentCell->nodeCount()==1)
    {
      newCells.push_back(currentCell);
    }
    else
    {
      modifiedCells.push_back(currentCell);
    }
  }
  return true;
}
template<class X, class C, class E>
bool WeightedOctree<X,C,E>::removeNode(X& node,
                                       std::vector<Cell*>& oldCells,
                                       std::vector<Cell*>& newCells,
                                       std::vector<Cell*>& modifiedCells)
{
  typename std::map<X,Cell*>::iterator iter;
  iter=this->NodeToCellMap.find(node);
  if(iter==this->NodeToCellMap.end())
  {
    return false;
  }
  Cell* currentCell=iter->second;


  E::setLastCell(node,iter->second);

  //remove node from cell.
  if(currentCell->removeNode(node))
  {
    this->NodeToCellMap.erase(node);
    int weight=E::weight(node);
    currentCell->addWeight(-weight);

    //update weight of parents.
    Cell* currentParent=currentCell->parent();
    while(currentParent)
    {
      currentParent->addWeight(-weight);
      currentParent=currentParent->parent();
    }

    //Determine if the parent can be collapsed.
    currentParent=currentCell->parent();
    Cell* lastParent=NULL;
    while(currentParent && currentParent->weight()<this->MaxBoxWeight)
    {
      //Traverse up the tree until the total weight is less that the max weight
      lastParent=currentParent;
      currentParent=currentParent->parent();
    }


    if(lastParent)
    {
      //Parent can be collapsed
      this->mergeCellChildren(lastParent,oldCells);
      newCells.push_back(lastParent);

      std::vector<X> cellNodes;
      lastParent->appendNodes(cellNodes);

      for(unsigned int j=0;j<cellNodes.size();j++)
      {
        X& localNode=cellNodes[j];
        iter=this->NodeToCellMap.find(localNode);
        if(iter!=this->NodeToCellMap.end())
        {
          E::setLastCell(localNode,iter->second);
        }

        this->NodeToCellMap[localNode]=lastParent;
      }
    }
    else
    {
      //Parent cannot be collapsed.
      if(currentCell->nodeCount()==0)
      {
        oldCells.push_back(currentCell);

      }
      else
      {
        modifiedCells.push_back(currentCell);
      }
    }


    return true;
  }
  return false;
}
template <class X, class C, class E>
bool WeightedOctree<X,C,E>::calculateCellCenter(Cell* cell)
{
  if(cell->nodeCount()==0)
  {
    return false;
  }
  std::vector<X> nodes;
  cell->appendNodes(nodes);

  //Calculate Center of cell.
  double wx=0.0;
  double wy=0.0;
  double wz=0.0;
  for(unsigned int i=0;i<nodes.size();i++)
  {
    const X&  node = nodes[i];
    CubitVector coords=E::coordinates(node);
    int weight=E::weight(node);

    wx+=coords.x()*weight;
    wy+=coords.y()*weight;
    wz+=coords.z()*weight;
  }
  wx/=cell->weight();
  wy/=cell->weight();
  wz/=cell->weight();


  cell->setCenter(wx,wy,wz);
  return true;
}
template <class X, class C, class E>
bool WeightedOctree<X,C,E>::splitCell(Cell* cell,std::vector<Cell*>& newLeafs)
{
  if(!cell->leaf())
  {
    return false;
  }

  if(cell->nodeCount()>=2 && cell->weight()>this->MaxBoxWeight)
  {
    for(int q=0;q<8;q++)
    {
      Cell* newCell= new Cell();
      newCell->setParent(cell);
      cell->setChild(q,newCell);
    }
    this->calculateCellCenter(cell);

    std::vector<X> nodes;
    cell->appendNodes(nodes);
    //Assign nodes to children based on center.
    int lastQuadrant=0;
    for(unsigned int i=0;i<nodes.size();i++)
    {
      const X&  node= nodes[i];
      CubitVector coords = E::coordinates(node);
      CubitVector center=cell->center();
      int quadrant=0;
      double distance=center.distance_between_squared(coords);
      if(distance<this->ToleranceSquared)
      {
        lastQuadrant++;
        if(lastQuadrant>=8)
        {
          lastQuadrant=0;
        }
        quadrant=lastQuadrant;
      }
      else
      {
        int x = coords.x() < center.x() ? Cell::X_MIN : Cell::X_MAX;
        int y = coords.y() < center.y() ? Cell::Y_MIN : Cell::Y_MAX;
        int z = coords.z() < center.z() ? Cell::Z_MIN : Cell::Z_MAX;
        quadrant=x|y|z;
      }

      cell->child(quadrant)->addNode(node);
      int weight=E::weight(node);
      cell->child(quadrant)->addWeight(weight);
    }
    int weight=cell->weight();
    cell->removeAllNodes();
    cell->setWeight(weight);

    for(unsigned int c=0;c<8;c++)
    {
      if(cell->child(c)->nodeCount()>0)
      {
        newLeafs.push_back(cell->child(c));
      }
    }


    return true;

  }

  return false;
}  


template<class X, class C, class E>
void WeightedOctree<X,C,E>::mergeCellChildren(Cell* cell,std::vector<Cell*>& oldLeafs)
{
  if(!cell->leaf())
  {
    std::vector<X> nodes;
    for(unsigned int i=0;i<8;i++)
    {
      this->mergeCellChildren(cell->child(i),oldLeafs);
      cell->child(i)->appendNodes(nodes);
      delete cell->child(i);
      cell->setChild(i,NULL);
    }
    for(unsigned int j=0;j<nodes.size();j++)
    {
      int nodeWeight=E::weight(nodes[j]);
      cell->addNode(nodes[j]);
      cell->addWeight(nodeWeight);
    }
  }
  else if(cell->nodeCount()!=0)
  {
    oldLeafs.push_back(cell);
  }
}



template <class X, class C, class E>
bool WeightedOctree<X,C,E>::rebalanceBranch(Cell* cell,
                                            std::vector<Cell*>& oldLeafs,
                                            std::vector<Cell*>& newLeafs)
{
  if(cell->leaf())
    return false;

  for(unsigned int i=0;i<8;i++)
  {
    if(!cell->child(i)->leaf())
    {
      //Found a child that is not a leaf;
      //Currently we only do 1 level rebalancing.
      return false;
    }
  }
  //All children are leafs;

  std::vector<X> nodes;
  cell->appendAllNodes(nodes);

  cell->setWeight(0);

  for(unsigned int i=0;i<nodes.size();i++)
  {
    const X& node = nodes[i];
    cell->addNode(node);

    int weight=E::weight(node);
    cell->addWeight(weight);
  }

  //First find all cells with nodes;
  if(cell->nodeCount()>=2 && cell->weight()>this->MaxBoxWeight)
  {
    // cell can be split. 
    //So we don't delete the children so they can be reused.
    for(unsigned int i=0;i<8;i++)
    {
      cell->child(i)->removeAllNodes();
    }
    this->recursiveSplitCell(cell,newLeafs);
  }
  else
  {
    //cell won't be split so we delete all of the children.
    for(unsigned int i=0;i<8;i++)
    {
      if(cell->child(i)->nodeCount()>0)
      {
        oldLeafs.push_back(cell->child(i));
      }
      delete cell->child(i);
      cell->setChild(i,NULL);
    }

    newLeafs.push_back(cell);
  }

  return true;
}


template <class X, class C, class E>
bool WeightedOctree<X,C,E>::recursiveSplitCell(Cell* cell,std::vector<Cell*>& newLeafs)
{
  if(!cell->leaf())
  {
    return false;
  }
  if(cell->nodeCount()>=2 && cell->weight()>this->MaxBoxWeight)
  {
    if(!cell->child(0))
    {
      for(int q=0;q<8;q++)
      {
        Cell* newCell= new Cell();
        newCell->setParent(cell);
        cell->setChild(q,newCell);
      }
    }
    this->calculateCellCenter(cell);

    std::vector<X> nodes;
    cell->appendNodes(nodes);
    //Assign nodes to children based on center.
    int lastQuadrant=0;
    for(unsigned int i=0;i<nodes.size();i++)
    {
      const X&  node= nodes[i];
      CubitVector coords = E::coordinates(node);
      CubitVector center=cell->center();
      int quadrant=0;
      double distance=center.distance_between_squared(coords);
      if(distance<this->ToleranceSquared)
      {
        lastQuadrant++;
        if(lastQuadrant>=8)
        {
          lastQuadrant=0;
        }
        quadrant=lastQuadrant;
      }
      else
      {
        int x = coords.x() < center.x() ? Cell::X_MIN : Cell::X_MAX;
        int y = coords.y() < center.y() ? Cell::Y_MIN : Cell::Y_MAX;
        int z = coords.z() < center.z() ? Cell::Z_MIN : Cell::Z_MAX;
        quadrant=x|y|z;
      }

      cell->child(quadrant)->addNode(node);
      int weight=E::weight(node);
      cell->child(quadrant)->addWeight(weight);
   }

    int weight=cell->weight();
    cell->removeAllNodes();
    cell->setWeight(weight);

    for(unsigned int c=0;c<8;c++)
    {
      if(cell->child(c)->nodeCount()>0)
      {
        if(!this->recursiveSplitCell(cell->child(c),newLeafs))
        {
          newLeafs.push_back(cell->child(c));
        }
      }
    }

    return true;

  }

  return false;


}

template <class X, class C, class E>
void WeightedOctree<X,C,E>::clear()
{
  delete this->Root;
  this->Root = new Cell;
}

template <class X, class C, class E>
bool WeightedOctree<X,C,E>::initilizeTree(std::vector<X>& nodes,std::vector<Cell*>& newCells)
{
  if(this->Root->leaf() && this->Root->nodeCount()==0)
  {
    for(unsigned int i=0;i<nodes.size();i++)
    {
      const X&  node= nodes[i];
      this->Root->addNode(node);
      int weight=E::weight(node);
      this->Root->addWeight(weight);
    }

    if(!this->recursiveSplitCell(this->Root,newCells))
    {
      //this cell wasn't split;
      newCells.push_back(this->Root);
    } 
    return true;
  }

  return false;
}
template <class X, class C, class E>
WeightedOctree<X,C,E>::Cell::Cell()
{
  this->TotalWeight=0;
  this->IsLeaf=true;


  this->Parent=NULL;
  for(int i=0;i<8;i++)
  {
    this->Children[i]=NULL;
  }
}

template <class X, class C, class E>
void WeightedOctree<X,C,E>::Cell::setParent(Cell* parent)
{
  this->Parent=parent;
}


template <class X, class C, class E>
void WeightedOctree<X,C,E>::Cell::setChild(int quadrant,Cell* c)
{
  assert( quadrant >= 0 && quadrant < 8);
  this->Children[quadrant]=c;
}

template <class X, class C, class E>
void WeightedOctree<X,C,E>::Cell::appendAllNodes( std::vector<X>& list )
{
  if(this->Nodes.size())
  {
    appendNodes( list );
  }
 
  for( int i = 0; i < 8; i++ )
  {
    if(this->Children[i])
      this->Children[i]->appendAllNodes(list);
  }

}

template <class X, class C, class E>
bool WeightedOctree<X,C,E>::Cell::leaf()
{
 return !this->Children[0];
}
template <class X, class C, class E>
void WeightedOctree<X,C,E>::Cell::appendNodes( std::vector<X>& list )
{
  list.insert(list.end(),this->Nodes.begin(),this->Nodes.end());
}

template <class X, class C, class E>
    bool WeightedOctree<X,C,E>::Cell::removeNode(const X&  node)
{
  typename std::vector<X>::iterator iter;
  iter=std::remove(this->Nodes.begin(),this->Nodes.end(),node);
  if(iter!=this->Nodes.end())
  {
    //The nodes was found;
    this->Nodes.erase(iter,this->Nodes.end());
    return true;
  }

  return false;
}
template <class X, class C, class E>
void WeightedOctree<X,C,E>::Cell::addNode(const X&  node )
{
  
  this->Nodes.push_back(node);
}
template <class X, class C, class E>
int WeightedOctree<X,C,E>::Cell::weight()
  {
  return this->TotalWeight;
  }

template <class X, class C, class E>
void WeightedOctree<X,C,E>::Cell::setWeight(int w)
{
   this->TotalWeight=w;
}
template <class X, class C, class E>
void WeightedOctree<X,C,E>::Cell::addWeight(int w)
{
   this->TotalWeight+=w;
}

template <class X, class C, class E>
    void WeightedOctree<X,C,E>::Cell::setCenter(CubitVector& center)
{
  this->Center=center;
}
template <class X, class C, class E>
    void WeightedOctree<X,C,E>::Cell::removeAllNodes()
{
  this->Nodes.clear();
  this->TotalWeight=0;
  
}

template <class X, class C, class E>
    void WeightedOctree<X,C,E>::Cell::setCenter(double& x,double& y,double& z)
{
  this->Center.set(x,y,z);
}

