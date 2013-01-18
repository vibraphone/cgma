//-------------------------------------------------------------------------
// Filename      : WeightedOctree.hpp
//
// Purpose       : Class for O(ln n) search for nodes within a specified
//                 tolerance of a specified position.
//
//
// Creator       : Corey McBride 
//
// Creation Date : 02/18/10
//-------------------------------------------------------------------------

#ifndef WEIGHTEDOCT_TREE_HPP
#define WEIGHTEDOCT_TREE_HPP

#include "CubitVector.hpp"
#include <vector>
#include <map>


template <class X> class DefaultWeightedNodeQuery
{
  public:
    static inline CubitVector coordinates(const X& x)
    {
      return x.coordinates();
    }
    static inline int weight(const X& x)
    {
      return x.weight();
    }

    template<class C>
    static inline void setLastCell( X& x,C* c)
    {
       x.setLastCell(c);
    }

};
struct DefaultOctreeCell
{
};


template < class X, class C = DefaultOctreeCell, class E = DefaultWeightedNodeQuery<X> >
class WeightedOctree
{
  public:
   class Cell: public C
    {

    public:
      Cell();
       
    
    bool removeNode(const X& node);
    void addNode( const X&  node );
      // asserts if not a leaf node.
      
    bool leaf(); 
      // Is this box a leaf in the tree?

    void appendNodes( std::vector<X>& list );
    const std::vector<X>& nodes() const { return Nodes; }

      // If this node is a leaf node, append all contained CubitNodes


    int nodeCount() { return Nodes.size(); }

    void setWeight(int w);
    void addWeight(int w);
    int weight();

    void setCenter(double& x,double& y,double& z);
    void setCenter(CubitVector& center);
    CubitVector& center() {return Center;}

    void removeAllNodes();



      enum {
        X_MIN = 0,
        Y_MIN = 0,
        Z_MIN = 0,
        X_MAX = 1,
        Y_MAX = 2,
        Z_MAX = 4
      };
      // Constants for selecting one of 8 child boxes.
      // Can be either bitwise-or'ed (|) or added (+).
      // MAX values take precedence over MIN values 
      // if both are specified.  (MIN is just the
      // absense of MAX).
      //
      Cell* child( int quadrant )
      {
        assert( quadrant >= 0 && quadrant < 8);
        return this->Children[quadrant];
      }
      // quadrant should be a bitwise OR (or sum) of 
      // the constants above.  If both a MIN and MAX 
      // value is specified for a given principal 
      // direction, the M  : children_(0),

      void setChild(int quadrant,Cell* c);

      void setParent(Cell*);
      Cell* parent() {return this->Parent;}

      void appendAllNodes( std::vector<X>& list );
      // If this node is a leaf node, append all contained CubitNodes.
      // Otherwise descend oct-tree getting CubitNodes contained in
      // all child leaf nodes.


    protected:
      Cell* Parent;
      Cell* Children[8]; 

      std::vector<X> Nodes;
      int TotalWeight;
      CubitVector Center;
      bool IsLeaf;

    };


   
  public:

    WeightedOctree(int max_weight_per_box,double tolerance=1e-6);
      //- Constructor
      //- tolerance - tolerance for comparison of node positions
      //- max_weight_per_box - subdivide box if total weight exceeds this 
      //- min_box_dimension - don't subdivide boxes in tree if
      //-                     all three dimensions are less than this

    ~WeightedOctree();
      // - Destructor
    
    bool initilizeTree(std::vector<X>& nodes,std::vector<Cell*>& newCells);
    void clear();

    bool addNode(const X& node,
      std::vector<Cell*>& oldCells,
      std::vector<Cell*>& newCells,
      std::vector<Cell*>& modifiedCells);
    bool removeNode(X& node,
      std::vector<Cell*>& oldCells,
      std::vector<Cell*>& newCells,
      std::vector<Cell*>& modifiedCells);



  protected:
    bool calculateCellCenter(Cell* cell); 

    void mergeCellChildren(Cell* cell,std::vector<Cell*>& oldLeafs);

    bool splitCell(Cell* cell,std::vector<Cell*>& newLeafs);
      // Change from leaf cell to non-leaf.  Splits nodes amongst 8 new
      // cell boxes using the value of centroid of the cell
      

     bool rebalanceBranch(Cell* cell,
         std::vector<Cell*>& oldLeafs,
         std::vector<Cell*>& newLeafs);

     bool recursiveSplitCell(Cell* cell,
         std::vector<Cell*>& newLeafs);



    int MaxBoxWeight;
      // values passed to constructor
  
    Cell* Root;
      // root of oct-tree

    int TotalWeight;
    
    std::map<X,Cell*> NodeToCellMap;

    double Tolerance;
    double ToleranceSquared;
  

 

};

#include "WeightedOctree.cpp"

#endif

