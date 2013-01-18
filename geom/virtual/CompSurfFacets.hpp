//-------------------------------------------------------------------------
// Filename      : CompSurfFacets.hpp
//
// Purpose       : Encapsulate facet date used to speed up composite surfaces
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 03/20/03
//-------------------------------------------------------------------------

#ifndef COMP_SURF_FACETS_HPP
#define COMP_SURF_FACETS_HPP

#include "CubitVector.hpp"
#include <vector>
#include <map>
#include "DLIList.hpp"

class Surface;
class GMem;
class CubitFacet;
class FacetEvalTool;


template <class X> class DLIList;

class CompSurfFacets
{
  typedef std::vector<CubitVector> PointList;
  typedef std::vector<int> IntegerList;
  typedef std::vector<Surface*> SurfPtrList;

  public:

    CompSurfFacets();
    ~CompSurfFacets();
    
    CubitStatus setup( const SurfPtrList& surface_data );
    
    int closest_index( const CubitVector& from_point,
                       CubitVector* point_on_facet = 0 ) const;
                       
    int closest_index( const CubitVector& from_point,
                       DLIList<int>& index_list,
                       CubitVector* point_on_facet = 0 );
                       
    void debug_draw_facets() const;
    
    void graphics( double tolerance, GMem& gmem );
    void set_ignore_flag(DLIList<int> &indicies, int flag);
    int get_ignore_flag(int index);
  
  protected:
    
    void closest_tri_point( IntegerList::const_iterator facet, 
                            const CubitVector& from_pt,
                            CubitVector& result_pt,
                            bool& interior ) const;
  
    void consolidate_points( double tolerance );
    void consolidate_few_points(double tolerance);
    void consolidate_many_points(double tolerance);
  
  private:
  
    bool pointsConsolidated;
    IntegerList ignoreFlags;
    DLIList<CubitFacet*> facetsToIgnore;
    
    std::map<CubitFacet*, int> facetToSurfaceMap;
    DLIList<CubitFacet*> allFacets;
    FacetEvalTool *facetEvalTool;
    DLIList<int> numFacetsPerSurface;

};

#endif

    
