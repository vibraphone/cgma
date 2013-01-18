//- Class: PointGridSearch
//- Description:  The PointGridSearch class is used to maintain a "bucket sort" of
//-               all mesh points used for rapid high performance nearest
//-               neighbor searches.  The object contains a 3-dimensional
//-               array of point lists for each grid cell (i.e. box) containing
//-               points that lie within the grid cell.  appropriate calls that
//-               return point, edge, face, and hex lists that lie in the
//-               neighborhood of an input mesh entity are provided.
//- Owner: Jim Hipp
//- Checked by: 

#include <math.h>
#include <assert.h>

#include "PointGridSearch.hpp"
#include "MemoryManager.hpp"
#include "DLIList.hpp"
#include "CastTo.hpp"
#include "CubitVector.hpp"
#include "CubitPoint.hpp"
#include "CubitFacet.hpp"
#include "CpuTimer.hpp"
#include "CubitMessage.hpp"
#include "CubitBox.hpp"
#include "TDCellIndex.hpp"
#include "SDLList.hpp"
#include "SDLIndexedDoubleList.hpp"
#include "GeometryDefines.h"

const double GRID_EXTENSION   =  0.1;

//
// Constructor / Destructor ...................................................
//
PointGridSearch::PointGridSearch( DLIList<CubitPoint*> *point_list,
                                  DLIList<CubitFacet*> *facet_list,
                                  float grid_scale)
{
    // find the geometric bounding range mesh
  assert( point_list && point_list->size() );
  bounding_range( *point_list );

    //find an estimate of the cell size based on the average square distance
    //of the edges on the surface of the geometry

  double grid_cell_size = 1.;
  
    // edge lengths
  double sum = 0.0;
  int i;
  maxEdgeLength = CUBIT_DBL_MIN;
  assert( facet_list || facet_list->size());
  for ( i = 0; i < facet_list->size(); i++)
  {
    CubitPoint *p1, *p2;
    facet_list->next(i)->get_edge_1(p1, p2);
    CubitVector temp = p2->coordinates() - p1->coordinates();
    double tmp_lsqrd =  temp.length_squared();
    sum += tmp_lsqrd;
  }
  grid_cell_size = sqrt( sum / facet_list->size() );

  maxEdgeLength = grid_cell_size;
    // evaluate grid parameters
  CubitVector cell(1.0, 1.0, 1.0);
  double grid_cell_width = grid_scale * grid_cell_size;
  cell *= ( GRID_EXTENSION + 5.0 ) * grid_cell_width;
  gridRangeMinimum = boundingRangeMinimum - cell;
  gridRangeMaximum = boundingRangeMaximum + cell;
  
  CubitVector range_width = gridRangeMaximum - gridRangeMinimum;
  
  numberGridCellsX = (int) ceil(range_width.x() / grid_cell_width + 0.5);
  numberGridCellsY = (int) ceil(range_width.y() / grid_cell_width + 0.5);
  numberGridCellsZ = (int) ceil(range_width.z() / grid_cell_width + 0.5);
  numberGridCells  = numberGridCellsX * numberGridCellsY * numberGridCellsZ;
  
  gridCellWidth.x(range_width.x() / ((double) numberGridCellsX));
  gridCellWidth.y(range_width.y() / ((double) numberGridCellsY));
  gridCellWidth.z(range_width.z() / ((double) numberGridCellsZ));
  
  gridCellWidthInverse.x(1.0 / gridCellWidth.x());
  gridCellWidthInverse.y(1.0 / gridCellWidth.y());
  gridCellWidthInverse.z(1.0 / gridCellWidth.z());
  
    // allocate neighborhood list array
  neighborhoodList = new DLIList<CubitPoint*>* [numberGridCells];
  assert(neighborhoodList != NULL);
  for ( i = 0; i < numberGridCells; i++)
  {
    neighborhoodList[i] = NULL;
  }
}
PointGridSearch::PointGridSearch(DLIList <CubitPoint*> &point_list,
                                 double grid_cell_size,
                                 double grid_scale)
{
    // find the bounding range of the reference geometry
  bounding_range(point_list);
  if ( grid_cell_size < CUBIT_RESABS )
    grid_cell_size = 
      (boundingRangeMaximum - boundingRangeMinimum).length() / 5.;    
  
  CubitVector cell(1.0, 1.0, 1.0);
  double grid_cell_width = grid_scale * grid_cell_size;
  cell *= grid_cell_width;
  gridRangeMinimum = boundingRangeMinimum - cell;
  gridRangeMaximum = boundingRangeMaximum + cell;

  CubitVector range_width = gridRangeMaximum - gridRangeMinimum;

  numberGridCellsX = (int) ceil(range_width.x() / grid_cell_width + 0.5);
  numberGridCellsY = (int) ceil(range_width.y() / grid_cell_width + 0.5);
  numberGridCellsZ = (int) ceil(range_width.z() / grid_cell_width + 0.5);
  numberGridCells  = numberGridCellsX * numberGridCellsY * numberGridCellsZ;

  gridCellWidth.x(range_width.x() / ((double) numberGridCellsX));
  gridCellWidth.y(range_width.y() / ((double) numberGridCellsY));
  gridCellWidth.z(range_width.z() / ((double) numberGridCellsZ));

  gridCellWidthInverse.x(1.0 / gridCellWidth.x());
  gridCellWidthInverse.y(1.0 / gridCellWidth.y());
  gridCellWidthInverse.z(1.0 / gridCellWidth.z());

    // allocate neighborhood list array

  neighborhoodList = new DLIList<CubitPoint*>* [numberGridCells];
  assert(neighborhoodList != NULL);
  int i;
  for ( i = 0; i < numberGridCells; i++)
  {
    neighborhoodList[i] = NULL;
  }
}
PointGridSearch::~PointGridSearch()
{
    // delete each DLIList<CubitPoint*> and the array of pointers containing them

    // first, delete the tdcellindex for each point
  int i;
  for ( i = 0; i < numberGridCells; i++) {
    if (!neighborhoodList[i]) continue;
    for (int j = neighborhoodList[i]->size(); j > 0; j--)
      neighborhoodList[i]->get_and_step()->delete_TD(&TDCellIndex::is_cell_index);
  }

    for (i = 0; i < numberGridCells; i++)
      delete neighborhoodList[i];

    delete [] neighborhoodList;
}

//
// Private Member Functions ...................................................
//

void PointGridSearch::cell_from_range()
{
    // Evaluate boundingCellMin and Max from the boundingRangeMin and Max
    // parameters

    CubitVector range_vec = boundingRangeMinimum - gridRangeMinimum;
    boundingCellMinimumX  = (int) (range_vec.x() * gridCellWidthInverse.x());
    boundingCellMinimumY  = (int) (range_vec.y() * gridCellWidthInverse.y());
    boundingCellMinimumZ  = (int) (range_vec.z() * gridCellWidthInverse.z());

    if (boundingCellMinimumX < 0) boundingCellMinimumX = 0;
    if (boundingCellMinimumY < 0) boundingCellMinimumY = 0;
    if (boundingCellMinimumZ < 0) boundingCellMinimumZ = 0;

    range_vec = boundingRangeMaximum - gridRangeMinimum;
    boundingCellMaximumX  = (int) (range_vec.x() * gridCellWidthInverse.x());
    boundingCellMaximumY  = (int) (range_vec.y() * gridCellWidthInverse.y());
    boundingCellMaximumZ  = (int) (range_vec.z() * gridCellWidthInverse.z());

    if (boundingCellMaximumX >= numberGridCellsX)
        boundingCellMaximumX = numberGridCellsX - 1;
    if (boundingCellMaximumY >= numberGridCellsY)
        boundingCellMaximumY = numberGridCellsY - 1;
    if (boundingCellMaximumZ >= numberGridCellsZ)
        boundingCellMaximumZ = numberGridCellsZ - 1;
}

void PointGridSearch::bounding_range(DLIList<CubitPoint*>& point_list)
{

  if ( !point_list.size() )
    return;
  
    // initialize min and max range values to the first point coordinates
  boundingRangeMinimum = point_list.get()->coordinates();
  boundingRangeMaximum = point_list.get_and_step()->coordinates();
  
    // find the min and max coordinates that completely encloses the
    // point list
  
  for (int i = 1; i < point_list.size(); i++)
    bounding_range(point_list.get_and_step());

}
void PointGridSearch::bounding_range(DLIList<CubitFacet*>& facet_list)
{
  int i, j;
  DLIList<CubitPoint*> point_list;
  for ( i = 0; i < facet_list.size(); i++)
  {
    CubitFacet* facet = facet_list.get_and_step();
    
    DLIList<CubitPoint*> temp_point_list;
    facet->points(temp_point_list);
    for (j = 0; j < temp_point_list.size(); j++)
    {
      CubitPoint* point = temp_point_list.get_and_step();
      
      if (!point->marked())
      {
        point->marked(CUBIT_TRUE);
        point_list.append(point);
      }
    }
  }
    // unmark the found points
  
  for (i = 0; i < point_list.size(); i++)
  {
    point_list.get_and_step()->marked(CUBIT_FALSE);
  }
  if ( !point_list.size() )
    return;
  
    // initialize min and max range values to the first point coordinates
  boundingRangeMinimum = point_list.get()->coordinates();
  boundingRangeMaximum = point_list.get_and_step()->coordinates();
  
    // find the min and max coordinates that completely encloses the
    // point list
  
  for (i = 1; i < point_list.size(); i++)
    bounding_range(point_list.get_and_step());

}

//
// Public Member Functions ....................................................
//

void PointGridSearch::change_cell(CubitPoint* point)
{
  int old_cell, current_cell;
  // if cell index has changed then swap lists
  ToolData* td_temp = point->get_TD( &TDCellIndex::is_cell_index);
  TDCellIndex* td_index = CAST_TO(td_temp, TDCellIndex);
  if (td_temp == NULL) {
    td_index = new TDCellIndex;
    assert( td_index != NULL );
    point->add_TD(td_index);
  }
  old_cell = td_index->cell_index();
  
  current_cell = grid_cell_index(point);
  if (current_cell >= 0 && current_cell < numberGridCells)
  {
    if (current_cell != old_cell)
    {
      remove_point_from_cell(point, old_cell);
      add_point_to_cell(point, current_cell);
      td_index->cell_index(current_cell);
    }
  }
}

int PointGridSearch::in_grid(const CubitVector& position)
{
  CubitVector range_vec;
  int i, j, k, index;

  range_vec = position - gridRangeMinimum;
  i = (int) (range_vec.x() * gridCellWidthInverse.x());
  j = (int) (range_vec.y() * gridCellWidthInverse.y());
  k = (int) (range_vec.z() * gridCellWidthInverse.z());
  index = numberGridCellsX * (numberGridCellsY * k + j) + i;
  if(index >= 0 && index < numberGridCells)
  {
    return CUBIT_TRUE;
  }
  else
  {
    return CUBIT_FALSE;
  }
}

void PointGridSearch::add_point(CubitPoint* point)
{
    // append point to the appropriate list

    int i = grid_cell_index(point);
    if (i < 0 || i > numberGridCells){
      return;
    }

    ToolData* td_temp = point->get_TD(&TDCellIndex::is_cell_index);
    TDCellIndex* td_index;
    
    if(td_temp == NULL){
      td_index = new TDCellIndex;
      assert(td_index != NULL);
      point->add_TD(td_index);
      td_index->cell_index(i); 
      add_point_to_cell(point, i);
    }
    else{  // the tool data already exists, avoid multiple listings of point
      td_index = CAST_TO(td_temp, TDCellIndex);
      assert( td_index != NULL);
      int iold = td_index->cell_index();
      
      if(iold != i){
	remove_point_from_cell(point,iold);
	td_index->cell_index(i);
	add_point_to_cell(point, i);
      }
      else{  // iold==i, so check whether the point is already on the list
	if (neighborhoodList[i] && neighborhoodList[i]->move_to(point)) return;
	else add_point_to_cell(point, i);
      }
    }
}

void PointGridSearch::remove_point(CubitPoint* point)
{
    // remove point from the appropriate list
    int index = -1;
    ToolData* td_temp = point->get_TD(&TDCellIndex::is_cell_index);
    if (!td_temp)
      return;
    TDCellIndex* td_index = CAST_TO(td_temp, TDCellIndex);
    if( td_index ){
      index = td_index->cell_index();
      point->delete_TD( &TDCellIndex::is_cell_index );
    }
    remove_point_from_cell(point, index);
}

void PointGridSearch::set_neighborhood_bounds(CubitVector& vec)
{
    // set the range and cell min and max values for the vector

    boundingRangeMinimum = vec;
    boundingRangeMaximum = vec;

    cell_from_range();
}

void PointGridSearch::set_neighborhood_bounds(const CubitVector& center,
                                              double size)
{
    // set the range and cell min and max values as a cube of breadth =
    // 2 * size centered at center

    CubitVector search_extent_box(size, size, size);
    boundingRangeMinimum = center - search_extent_box;
    boundingRangeMaximum = center + search_extent_box;

    cell_from_range();
}
void PointGridSearch::set_neighborhood_bounds_max_edge(const CubitVector& center,
                                                       double factor)
{
    // set the range and cell min and max values as a cube of breadth =
    // 2 * size centered at center

    CubitVector search_extent_box(maxEdgeLength*factor, 
                                  maxEdgeLength*factor, 
                                  maxEdgeLength*factor);
    boundingRangeMinimum = center - search_extent_box;
    boundingRangeMaximum = center + search_extent_box;

    cell_from_range();
}

void PointGridSearch::set_neighborhood_bounds(CubitVector& center, double x,
					 double y, double z)
{
    // set the range and cell min and max values as a box of size 2x, 2y, and
    // 2z, centered at center

    CubitVector search_extent_box(x, y, z);
    boundingRangeMinimum = center - search_extent_box;
    boundingRangeMaximum = center + search_extent_box;

    cell_from_range();
}

void PointGridSearch::set_neighborhood_bounds(CubitVector& vec_1,
					 CubitVector& vec_2)
{
    // initialize min and max range values to the first vectors coordinates
   
    boundingRangeMinimum = vec_1;
    boundingRangeMaximum = vec_1;

    // find the range and cell min and max values for the two vectors

    bounding_range(vec_2);
    cell_from_range();
}

void PointGridSearch::set_neighborhood_bounds(CubitPoint* point, CubitVector& vec)
{
    // initialize min and max range values to the point coordinates

    boundingRangeMinimum = point->coordinates();
    boundingRangeMaximum = point->coordinates();

    // find the range and cell min and max values for the point and vector

    bounding_range(vec);
    cell_from_range();
}

void PointGridSearch::set_neighborhood_bounds(CubitPoint* point)
{
    // initialize min and max range values to the point coordinates

    boundingRangeMinimum = point->coordinates();
    boundingRangeMaximum = point->coordinates();

    // find the cell min and max values

    cell_from_range();
}

void PointGridSearch::set_neighborhood_bounds(CubitPoint* point_1, CubitPoint* point_2)
{
    // initialize min and max range values to the first point coordinates

    boundingRangeMinimum = point_1->coordinates();
    boundingRangeMaximum = point_1->coordinates();

    // find the range and cell min and max values for the two points

    bounding_range(point_2);
    cell_from_range();
}

void PointGridSearch::set_neighborhood_bounds(DLIList<CubitPoint*>& point_list)
{
    // find the range and cell min and max values for the point list

    bounding_range(point_list);
    cell_from_range();
}

void PointGridSearch::set_neighborhood_bounds(CubitFacet* facet)
{
    // find the points belonging to the face

    DLIList<CubitPoint*> point_list;
    facet->points(point_list);

    // find the bounding range and cell min and max values for the point list

    bounding_range(point_list);
    cell_from_range();
}

void PointGridSearch::set_neighborhood_bounds(DLIList<CubitFacet*>& facet_list)
{
    DLIList<CubitPoint*> point_list;

    // find all points belonging to the faces in face_list
    int i;
    for ( i = 0; i < facet_list.size(); i++)
    {
      CubitFacet* facet = facet_list.get_and_step();

      DLIList<CubitPoint*> temp_point_list;
      facet->points(temp_point_list);
      for (int j = 0; j < temp_point_list.size(); j++)
      {
	CubitPoint* point = temp_point_list.get_and_step();

	if (!point->marked())
	{
	  point->marked(CUBIT_TRUE);
	  point_list.append(point);
	}
      }
    }

    // unmark the found points

    for (i = 0; i < point_list.size(); i++)
    {
      point_list.get_and_step()->marked(CUBIT_FALSE);
    }

    // find the bounding range and cell min and max values for the point list

    bounding_range(point_list);
    cell_from_range();
}

void PointGridSearch::get_neighborhood_points( DLIList<CubitPoint*> &point_list )
{
  // retrieve points over the current bounding box range

  point_list.clean_out();

  for (int k = boundingCellMinimumZ; k <= boundingCellMaximumZ; k++)
  {
    int kn = numberGridCellsY * k;
    for (int j = boundingCellMinimumY; j <= boundingCellMaximumY; j++)
    {
      int jn = numberGridCellsX * (kn + j);
      for (int i = boundingCellMinimumX; i <= boundingCellMaximumX; i++)
      {
        int in = jn + i;
        assert( in >= 0 && in < numberGridCells );
        if (neighborhoodList[in])
        {
           point_list += *(neighborhoodList[in]);
        }
      }
    }
  }
}

void PointGridSearch::get_neighborhood_points_sorted(
  DLIList<CubitPoint*> &point_list,
  const CubitVector& center, double cut_off)
{
  point_list.clean_out();
  DLIList<CubitPoint*> temp_point_list;
  int i;

  for (int k = boundingCellMinimumZ; k <= boundingCellMaximumZ; k++)
  {
    int kn = numberGridCellsY * k;
    for (int j = boundingCellMinimumY; j <= boundingCellMaximumY; j++)
    {
      int jn = numberGridCellsX * (kn + j);
      for ( i = boundingCellMinimumX; i <= boundingCellMaximumX; i++)
      {
	int in = jn + i;
        if (neighborhoodList[in])
        {
          temp_point_list += *(neighborhoodList[in]);
        }
      }
    }
  }

  // evaluate point distance to center ... remove those larger than cut_off
  SDLByDouble sorted_index_list;
  IndexedDouble *ID;  
  CubitVector vec;
  cut_off *= cut_off;
  temp_point_list.reset();
  for ( i = 0; i < temp_point_list.size(); i++)
  {
    vec = center - temp_point_list.get_and_step()->coordinates();
    double distance = vec.length_squared();
    if (distance < cut_off)
    {
      ID = new IndexedDouble( i, distance );
      sorted_index_list.append( ID );
    }
  }

  sorted_index_list.sort();
  temp_point_list.reset();
  for ( i = 0; i < sorted_index_list.size(); i++ )
  {
    ID = sorted_index_list.get_and_step();
    point_list.append( temp_point_list.next( ID->index() ) );
    delete ID;
  }
}

void PointGridSearch::get_neighborhood_facets( DLIList<CubitFacet*> &facet_list )
{
  // retrieve points over the current bounding box range

  facet_list.clean_out();
  DLIList<CubitPoint*> point_list;
  get_neighborhood_points( point_list );

  // retrieve all faces attached to the points in point_list

  for (int i = 0; i < point_list.size(); i++)
  {
    CubitPoint* point = point_list.get_and_step();

    DLIList<CubitFacet*> temp_facet_list;
    point->facets(temp_facet_list);

    for (int j = 0; j < temp_facet_list.size(); j++)
    {
      CubitFacet* facet = temp_facet_list.get_and_step();

      if (!facet->marked())
      {
        facet->marked(CUBIT_TRUE);
        facet_list.append(facet);
      }
    }
  }

  // unmark the found faces and return face_list

  for (int m = 0; m < facet_list.size(); m++)
  {
    facet_list.get_and_step()->marked(CUBIT_FALSE);
  }

}

double PointGridSearch::fraction_empty_cells()
{
  int empty_cell = 0;

  for (int i = 0; i < numberGridCells; i++)
  {
    if (!neighborhoodList[i]) empty_cell++;
  }

  return ((double) (empty_cell) / (double) (numberGridCells));
}

double PointGridSearch::average_points_per_occupied_cell()
{
  int total_points    = 0;
  int occupied_cells = 0;

  for (int i = 0; i < numberGridCells; i++)
  {
    if (neighborhoodList[i])
    {
      total_points += neighborhoodList[i]->size();
      occupied_cells++;
    }
  }

  return ((double) (total_points) / (double) (occupied_cells));
}

int PointGridSearch::total_points_in_grid_cells()
{
  int total_points = 0;

  for (int i = 0; i < numberGridCells; i++)
  {
    if (neighborhoodList[i])
    {
      total_points += neighborhoodList[i]->size();
    }
  }

  return total_points;
}


void PointGridSearch::add_point_to_cell(CubitPoint* point, int i)
{
  if (neighborhoodList[i])
  {
    neighborhoodList[i]->append(point);
  }
  else
  {
    neighborhoodList[i] = new DLIList<CubitPoint*>;
    neighborhoodList[i]->append(point);
  }
}
void PointGridSearch::remove_point_from_cell(CubitPoint* point, int i)
{
  if( i < 0 ) return;
  if (neighborhoodList[i])
  { 
    if (neighborhoodList[i]->move_to(point))
    {
      neighborhoodList[i]->extract();
    }
    if (!neighborhoodList[i]->size())
    {
      delete neighborhoodList[i];
      neighborhoodList[i] = NULL;
    }
  }
}

