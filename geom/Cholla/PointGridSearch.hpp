//- Class: PointGridSearch
//- Description:  The PointGridSearch class is used to maintain a "bucket sort" of
//-               all points used for rapid high performance nearest
//-               neighbor searches.  The object contains a 3-dimensional
//-               array of point lists for each grid cell (i.e. box) containing
//-               points that lie within the grid cell.  appropriate calls that
//-               return point lists that lie in the
//-               neighborhood of an input mesh entity are provided.
//-
//-               PointGridSearch is a simple neighborhood searcher that stores
//-               CubitPoints in a data base broken up into small cubes about
//-               the meshing region of interest.  PointGridSearch provides
//-               functionality to return all points,
//-               that lie within the proximity of a user specified
//-               neighborhood.  The neighborhood is a cartesian bounding box
//-               that can be defined by the user in one of many ways.  These
//-               include:
//-
//-               1) a point (defined by a vector or a points coordinates) in
//-                  which case the neighborhood is the grid cell that
//-                  contains the point,
//-
//-               2) a bounding box that encloses two points (defined by
//-                  vectors or nodal coordinates),
//-
//-               3) a bounding box that encloses a list of points
//-
//-               4) a bounding box that encloses a CubitFacet
//-
//-
//-
//- Owner: David White
//- Checked by: 
//- Modified From: GridSearch

#ifndef POINT_GRID_SEARCH_HPP
#define POINT_GRID_SEARCH_HPP

#include "CubitPoint.hpp"
#include "CubitVector.hpp"

template <class X> class DLIList;
class MemoryManager;

class PointGridSearch
{
private:

  DLIList<CubitPoint*>** neighborhoodList;
    //
    //- array of CubitPoint lists for each grid cell
  double maxEdgeLength;
    //- max EdgeLength as found in constructor.

  void add_point_to_cell(CubitPoint* point, int i);
  void remove_point_from_cell(CubitPoint* point, int i);
    //
    //- functions that add or remove a point from neighborhoodList[i]
  
protected:
   
  CubitVector gridRangeMinimum;
  CubitVector gridRangeMaximum;
  CubitVector gridCellWidth;
  CubitVector gridCellWidthInverse;
  int         numberGridCellsX;
  int         numberGridCellsY;
  int         numberGridCellsZ;
  int         numberGridCells;
    //
    //- search grid parameters

  CubitVector boundingRangeMinimum;
  CubitVector boundingRangeMaximum;
  int         boundingCellMinimumX;
  int         boundingCellMinimumY;
  int         boundingCellMinimumZ;
  int         boundingCellMaximumX;
  int         boundingCellMaximumY;
  int         boundingCellMaximumZ;
    //
    //- bounding box range and cell parameters

  int grid_cell_index(CubitPoint* point);
    //
    //- returns the cell index of the point


  virtual void cell_from_range();
    //
    //- evaluates the integer cell range given the double entity range
  
  void bounding_range(DLIList<CubitPoint*>& point_list);
  void bounding_range(DLIList<CubitFacet*>& facet_list);
  void bounding_range(CubitPoint* point);
  void bounding_range(const CubitVector& vec);
    //
    //- evaluates the bounding box entity range from the MRefEntity,
    //- point_list,
    //- point, or vector

  public:

    PointGridSearch( DLIList<CubitPoint*> *point_list,
                    DLIList<CubitFacet*> *facet_list,
                    float grid_scale = 2.0);
    PointGridSearch(DLIList <CubitPoint*> &point_list,
                    double grid_cell_size,
                    double grid_scale);
  

   virtual ~PointGridSearch();
    //
    //- grid search constructors and destructor

    CubitVector grid_range_minimum();
    CubitVector grid_range_maximum();
    CubitVector grid_cell_width();
    int         number_grid_cells_x();
    int         number_grid_cells_y();
    int         number_grid_cells_z();
    int         number_grid_cells();
    //
    //- returns grid parameters

    double      fraction_empty_cells();
    double      average_points_per_occupied_cell();
    int         total_points_in_grid_cells();
    //
    //- returns grid search data base information

    CubitVector bounding_range_minimum();
    CubitVector bounding_range_maximum();
    int         bounding_cell_minimum_x();
    int         bounding_cell_minimum_y();
    int         bounding_cell_minimum_z();
    int         bounding_cell_maximum_x();
    int         bounding_cell_maximum_y();
    int         bounding_cell_maximum_z();
    //
    //- returns bounding box range and cell parameters

    void add_point(CubitPoint* point);
    void remove_point(CubitPoint* point);
    //
    //- adds or removes points from the neighbor lists

    int in_grid(const CubitVector& position);
    //
    //- Determines whether the passed-in position is within the
    //- current grid. Returns CUBIT_TRUE or CUBIT_FALSE.

    void change_cell(CubitPoint* point);
    //
    //- changes a points cell (removes from old list and appends to new one)
    //- if it has changed

    void set_neighborhood_bounds(CubitVector& vec);
    void set_neighborhood_bounds(const CubitVector& center, double size);
    void set_neighborhood_bounds(CubitVector& center, double x,
				 double y, double z);
    void set_neighborhood_bounds(CubitVector& vec_1, CubitVector& vec_2);
    void set_neighborhood_bounds(CubitPoint* point, CubitVector& vec);
    void set_neighborhood_bounds(CubitPoint* point);
    void set_neighborhood_bounds(CubitPoint* point_1, CubitPoint* point_2);
    void set_neighborhood_bounds(DLIList<CubitPoint*>& point_list);
    void set_neighborhood_bounds(CubitFacet* facet);
    void set_neighborhood_bounds(DLIList<CubitFacet*>& facet_list);
    void set_neighborhood_bounds_max_edge(const CubitVector &center,
                                          double factor = 1.0);
    //
    //- evaluates the bounding box (entity range and cell range) from the
    //- input entity(ies)

    void get_neighborhood_points(DLIList<CubitPoint*> &point_list);
    void get_neighborhood_facets(DLIList<CubitFacet*> &facet_list);
    void get_neighborhood_points_sorted(DLIList<CubitPoint*> &point_list,
                                       const CubitVector& center,
                                       double cut_off);
    //
    //- returns lists of entities that occupy the cell range set by a previous
    //- bounding box call


};

inline CubitVector PointGridSearch::grid_range_minimum()  {return gridRangeMinimum;}
inline CubitVector PointGridSearch::grid_range_maximum()  {return gridRangeMaximum;}
inline CubitVector PointGridSearch::grid_cell_width()     {return gridCellWidth;}
inline int         PointGridSearch::number_grid_cells_x() {return numberGridCellsX;}
inline int         PointGridSearch::number_grid_cells_y() {return numberGridCellsY;}
inline int         PointGridSearch::number_grid_cells_z() {return numberGridCellsZ;}
inline int         PointGridSearch::number_grid_cells()   {return numberGridCells;}

inline CubitVector PointGridSearch::bounding_range_minimum()
                               {return boundingRangeMinimum;}
inline CubitVector PointGridSearch::bounding_range_maximum()
                               {return boundingRangeMaximum;}
inline int         PointGridSearch::bounding_cell_minimum_x()
                               {return boundingCellMinimumX;} 
inline int         PointGridSearch::bounding_cell_minimum_y()
                               {return boundingCellMinimumY;} 
inline int         PointGridSearch::bounding_cell_minimum_z()
                               {return boundingCellMinimumZ;} 
inline int         PointGridSearch::bounding_cell_maximum_x()
                               {return boundingCellMaximumX;} 
inline int         PointGridSearch::bounding_cell_maximum_y()
                               {return boundingCellMaximumY;} 
inline int         PointGridSearch::bounding_cell_maximum_z()
                               {return boundingCellMaximumZ;} 

inline int         PointGridSearch::grid_cell_index(CubitPoint* point)
{
  CubitVector range_vec = point->coordinates();
  range_vec -= gridRangeMinimum;

  int i = (int) (range_vec.x() * gridCellWidthInverse.x());
  int j = (int) (range_vec.y() * gridCellWidthInverse.y());
  int k = (int) (range_vec.z() * gridCellWidthInverse.z());

  return numberGridCellsX * (numberGridCellsY * k + j) + i;
}

inline void PointGridSearch::bounding_range(const CubitVector& vec)
{
  if (vec.x() < boundingRangeMinimum.x())
    boundingRangeMinimum.x(vec.x());
  else if (vec.x() > boundingRangeMaximum.x())
    boundingRangeMaximum.x(vec.x());

  if (vec.y() < boundingRangeMinimum.y())
    boundingRangeMinimum.y(vec.y());
  else if (vec.y() > boundingRangeMaximum.y())
    boundingRangeMaximum.y(vec.y());

  if (vec.z() < boundingRangeMinimum.z())
    boundingRangeMinimum.z(vec.z());
  else if (vec.z() > boundingRangeMaximum.z())
    boundingRangeMaximum.z(vec.z());
}

inline void PointGridSearch::bounding_range(CubitPoint* point)
{
  double xpos = point->x();
  double ypos = point->y();
  double zpos = point->z();
  
  if (xpos < boundingRangeMinimum.x())
    boundingRangeMinimum.x(xpos);
  else if (xpos > boundingRangeMaximum.x())
    boundingRangeMaximum.x(xpos);

  if (ypos < boundingRangeMinimum.y())
    boundingRangeMinimum.y(ypos);
  else if (ypos > boundingRangeMaximum.y())
    boundingRangeMaximum.y(ypos);

  if (zpos < boundingRangeMinimum.z())
    boundingRangeMinimum.z(zpos);
  else if (zpos > boundingRangeMaximum.z())
    boundingRangeMaximum.z(zpos);
}

#endif // GRID_SEARCH_BASE_HPP

    
