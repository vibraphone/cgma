//-------------------------------------------------------------------------
// Filename      : CompSurfFacets.cpp
//
// Purpose       : Encapsulate facet date used to speed up composite surfaces
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 03/20/03
//-------------------------------------------------------------------------

#include "CompSurfFacets.hpp"
#include "Surface.hpp"
#include "GeometryQueryEngine.hpp"
#include "GMem.hpp"
#include "OctTree.hpp"
#include "GfxDebug.hpp"

#include "CubitFacetData.hpp"
#include "CubitFacetEdgeData.hpp"
#include "CubitPointData.hpp"
#include "FacetEvalTool.hpp"
#include "FacetDataUtil.hpp"
#include "DLIList.hpp"


CompSurfFacets::CompSurfFacets()
{
  pointsConsolidated = false;
  facetEvalTool = NULL;
}

CompSurfFacets::~CompSurfFacets()
{
  //delete the facetEvalTool...deletes the facets also
  if (facetEvalTool)
    delete facetEvalTool;
  facetEvalTool = NULL;
  
  //delete the facets not in the tool
  FacetDataUtil::delete_facets( facetsToIgnore );
}

CubitStatus CompSurfFacets::setup( const SurfPtrList& surface_data )
{
  GMem gmem;
  DLIList<CubitFacet*> facet_list;
  DLIList<CubitPoint*> point_list;

  ignoreFlags.resize(surface_data.size());

  for ( unsigned int index = 0; index < surface_data.size(); index++ )
  {
    Surface* surface = surface_data[index];
    ignoreFlags[index] = 0;
    
    CubitStatus result = surface->get_geometry_query_engine()
      ->get_graphics( surface, &gmem, 2);
    if ( !result )
      return CUBIT_FAILURE;
 

    //make the points
    DLIList<CubitPoint*> tmp_points;
    GPoint* g_itor = gmem.point_list();
    GPoint* g_end = g_itor + gmem.pointListCount;
    for ( ; g_itor != g_end; ++g_itor )
    {
      CubitPoint *new_point_ptr =
           new CubitPointData(g_itor->x, g_itor->y, g_itor->z); 
      tmp_points.append( new_point_ptr );
    }

    //make the facets
    int* f_itor = gmem.facet_list();
    int* f_end = f_itor + gmem.fListCount;
    int num_facets = 0;
    while ( f_itor != f_end )
    {
      int count = *f_itor++;
      assert( count == 3 );
      
      CubitPoint *facet_points[3];

      for (int i=0; i<count; i++ )
      {
        int index2 = *f_itor++;
        facet_points[i] = tmp_points[index2];
      }
      
      CubitFacet *facet = (CubitFacet*) new CubitFacetData( 
        facet_points[0], facet_points[1], facet_points[2] );

      facet_list.append( facet );
      allFacets.append( facet );

      facetToSurfaceMap.insert( std::map<CubitFacet*, int>::value_type( facet, index ) ); 
      num_facets++;
    }
    point_list += tmp_points;
    numFacetsPerSurface.append( num_facets );
  }

  //create a FacetEvalTool with the points and facets
  facetEvalTool = new FacetEvalTool( facet_list, point_list, -1, 0.707106 );
  
  pointsConsolidated = false;
  return CUBIT_SUCCESS;
}

void CompSurfFacets::set_ignore_flag(DLIList<int> &indicies, int flag)
{
  if( indicies.size() == 0 )
    return;

  DLIList<CubitFacet*> facets_to_use;
  facetsToIgnore.clean_out();

  int where_to_start = 0;
  int i;
  for( i=0; i<ignoreFlags.size(); i++ ) 
  {
    ignoreFlags[i] = 0;
    
    //should we ignore this one?
    int j;
    if( flag == 1 )
    {
      for( j=indicies.size(); j--; )
      {
        if( indicies.get_and_step() == i )
          ignoreFlags[i] = 1;
      }
    }

    int num_facets = numFacetsPerSurface[i];

    //put the ignored ones in a list to prevent memory leak 
    for( j=where_to_start; j<where_to_start+num_facets; j++ )
    {
      if( ignoreFlags[i] == 0 )
        facets_to_use.append( allFacets[j] );
      else
        facetsToIgnore.append( allFacets[j] );
    }

    where_to_start += numFacetsPerSurface[i];
  }

  //give these facets to the tool
  facetEvalTool->replace_facets( facets_to_use ); 

  //now reset the bounding box on the 
  facetEvalTool->reset_bounding_box();

}

int CompSurfFacets::get_ignore_flag(int index)
{
   return ignoreFlags[index];
}

int CompSurfFacets::closest_index( const CubitVector& pos,
                                   CubitVector* closest_pos ) const
{
  CubitFacet *closest_facet = facetEvalTool->closest_facet( pos );
 
  std::map<CubitFacet*, int>::const_iterator my_iter;
  my_iter = facetToSurfaceMap.find( closest_facet );
  
  int closest_index = -1;

  if( my_iter != facetToSurfaceMap.end() )
    closest_index = my_iter->second;

  if( closest_pos )
  {
    CubitPoint *pt1, *pt2;
    closest_facet->closest_point_trimmed( pos, *closest_pos, pt1, pt2 );
  }

  return closest_index;
}

int CompSurfFacets::closest_index( const CubitVector& pos,
                                   DLIList<int>& index_list,
                                   CubitVector* closest_pos ) 
{
  //facets must be consolidated before calling this function
  //if not, you won't know if a position is on multiple
  //surfaces (on seam between 2 surfaces)
  if (!pointsConsolidated)
  {
    consolidate_points(GEOMETRY_RESABS);
    pointsConsolidated = true;
  }

  CubitFacet *closest_facet = facetEvalTool->closest_facet( pos );
  
  std::map<CubitFacet*, int>::const_iterator iter;
  iter = facetToSurfaceMap.find( closest_facet );
  
  int closest_index = -1;

  if( iter != facetToSurfaceMap.end() )
    closest_index = iter->second;

  index_list.append( closest_index );

  CubitVector tmp_closest_pos;
  CubitPoint *pt1, *pt2;
  closest_facet->closest_point_trimmed( pos, tmp_closest_pos, pt1, pt2 );

  if( closest_pos )
    *closest_pos = tmp_closest_pos;
 
  double tol_sq = GEOMETRY_RESABS*GEOMETRY_RESABS;
  double dist_sq = pos.distance_between_squared( tmp_closest_pos );  
  
  //if the closest distance is less than resabs, there could be more facets 
  //less then resabs away too....find them
  if( dist_sq < tol_sq )
  {
    //examine neighboring facets
    CubitPoint *point0, *point1, *point2;
    closest_facet->points( point0, point1, point2 ); 
   
    DLIList<CubitFacet*> neighbor_facets;
    point0->tris( neighbor_facets );
    point1->tris( neighbor_facets );
    point2->tris( neighbor_facets );

    neighbor_facets.uniquify_unordered();
    
    int k;
    for( k=neighbor_facets.size(); k--; )
    {
      CubitFacet *tmp_facet = neighbor_facets.get_and_step();
      tmp_facet->closest_point_trimmed( pos, tmp_closest_pos, pt1, pt2 );
     
      dist_sq = pos.distance_between_squared( tmp_closest_pos );  
      if( dist_sq < tol_sq )
      {
        iter = facetToSurfaceMap.find( tmp_facet );
        if( iter != facetToSurfaceMap.end() )
          index_list.append( iter->second );
      }
    }
  }

  index_list.uniquify_unordered();
  
  return index_list.size();
}

void CompSurfFacets::debug_draw_facets() const
{
  int i;
  for( i=0; i<allFacets.size(); i++ )
    allFacets[i]->debug_draw( 4 );

  GfxDebug::flush();
}

void CompSurfFacets::consolidate_points( double tolerance )
{
  DLIList<CubitFacet*> facet_list;
  facetEvalTool->get_facets( facet_list );

  //group the points by surface
  DLIList< DLIList<CubitFacet*>* > list_of_lists;
  int i;
  numFacetsPerSurface.reset();
  int index_into_facet_array = 0;
  for( i=0; i<numFacetsPerSurface.size(); i++ )
  {
    int num_facets = numFacetsPerSurface.get_and_step();

    if( ignoreFlags[i] == 0 )
    {
      DLIList<CubitFacet*>* tmp_facet_list = new DLIList<CubitFacet*>;
      int k;
      for(k=0; k<num_facets; k++ )
        tmp_facet_list->append( allFacets[index_into_facet_array+k] ); 
      list_of_lists.append( tmp_facet_list );
    }
    index_into_facet_array += num_facets;
  }
  
  int num_points_merged, num_edges_merged;
  DLIList<CubitPoint*> unmerged_points;
  FacetDataUtil::merge_coincident_vertices( list_of_lists, tolerance,
                                            num_points_merged, num_edges_merged,
                                            unmerged_points );

  //replace facets so points list get updates
  facetEvalTool->replace_facets( facet_list );

  //clean up
  for( i=list_of_lists.size(); i--; )
    delete list_of_lists.get_and_step();

  pointsConsolidated = true;
}

class CubitVectOctTreeEval {
public:
static inline const CubitVector& coordinates(CubitVector* p)
  { return *p; };
};

void CompSurfFacets::graphics( double tolerance, GMem& gmem )
{
  if (!pointsConsolidated)
  {
    consolidate_points(tolerance);
    pointsConsolidated = true;
  }
  
  DLIList<CubitPoint*> point_list;
  DLIList<CubitFacet*> facet_list;
  facetEvalTool->get_points( point_list );
  facetEvalTool->get_facets( facet_list );
  point_list.reset();
  facet_list.reset();

  //mark all the points to -1
  int i;
  for( i=point_list.size(); i--; )
    point_list.get_and_step()->marked(-1);

  gmem.allocate_tri( facet_list.size() );
 
  int* f_itor = gmem.facet_list();
  DLIList<CubitPoint*> tmp_points;  
  int index = 0;

  for( i=facet_list.size(); i--; )
  {
    CubitFacet *tmp_facet = facet_list.get_and_step();
    CubitPoint *p0, *p1, *p2;
    tmp_facet->points(p0, p1, p2 );

    *(f_itor++) = 3;

    if( p0->marked() == -1 )
    {
      tmp_points.append( p0 ); 
      p0->marked( index ); 
      index++;
    }
    *(f_itor++) = p0->marked();

    if( p1->marked() == -1 )
    {
      tmp_points.append( p1 ); 
      p1->marked( index ); 
      index++;
    }
    *(f_itor++) = p1->marked();

    if( p2->marked() == -1 )
    {
      tmp_points.append( p2 ); 
      p2->marked( index ); 
      index++;
    }
    *(f_itor++) = p2->marked();
  }


  GPoint* g_itor = gmem.point_list();
  for( i=tmp_points.size(); i--; )
  {
    CubitPoint *tmp_point = tmp_points.get_and_step();
    g_itor->x = (float)tmp_point->x();
    g_itor->y = (float)tmp_point->y();
    g_itor->z = (float)tmp_point->z();
    g_itor++;
  }
  
  gmem.fListCount = (facet_list.size() * 4);
  gmem.pointListCount = tmp_points.size();
}

