/** \class TDFacetboolData
//- Class: TDFacetboolData
//- Owner: John Fowler
//- Description: Data for maintaining persistent IDs.
//-              surfaceIndex and edgeIndex[] are indices into a list of entities.
//- Checked By: 
//- Version:
*/

#ifndef TDFACETBOOLDATA_HPP
#define TDFACETBOOLDATA_HPP

#include "CubitDefines.h"
#include "ToolData.hpp"
class FacetEntity;
class CubitFacet;

class TDFacetboolData : public ToolData
{
private:
    //!  An index into an array of FacetSurfaces denoting which FacetSurface
    //!  this triangle came from.
  int surfaceIndex;
    //!  An index into an array of CubitFacetEdges denoting which
    //!  CubitFacetEdges (if any) that the three triangle edges came from. 
  int edgeIndex[3]; 
    //!  true if this triangle came from body 1; false if from body 2. 
  bool body_1_flag;
    //!Is the associated facet reversed or not.  The determines whether
    //! we need to adjust the edge array before returning it.
  bool isReversed;
    //!To try to avoid a problem with the order of the edge indices,
    //! we force the calling code to provide whether the associated facet
    //! is reversed or not.  The public version is below.
  int *get_edge_indices() { return edgeIndex; }
  
public:

  TDFacetboolData();
    //- constructor

  ~TDFacetboolData();

  static int is_facetbool_facet(const ToolData* td);

  void set(int sv, int e0v, int e1v, int e2v, bool parent, bool is_reversed); 

  static TDFacetboolData* get(CubitFacet *facet_ptr);    

  static CubitStatus add_facetbool_facet(FacetEntity *facet_entity);
  
    //!Get the edge indices array.  Give a flag telling whether the facet
    //! is backwards.  If it has been marked as reversed (or forward) since
    //! the last time this function was called, the array will be reversed.
    //! This is a little dangerous since the calling code can have a pointer
    //! to this array, so please be careful.
  int *get_edge_indices(bool rev_flag);

  int get_surf_index() { return surfaceIndex; }
  
  bool parent_is_body_1() { return body_1_flag; }
    //! Returns the isReversed.  This allows to determine whether the facets
    //! orientation has changed since the egde array was last modified.
  bool is_reversed()
    {return isReversed;}

    //!  Manually set the is_reversed flag.  This shouldn't be needed.
    //! By calling get_edge_indices(bool rev_flag), the isReversed flag
    //! will be automatically modified if necessary.
  void is_reversed(bool rev);
  
};
    

#endif // TDFACETBOOLDATA_HPP


