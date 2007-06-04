//- Class:       TDFacetboolData
//- Description: Tool data for holding surface and curve index information.
//- Owner:       John Fowler
//- Checked by:
//- Version:
#include "TDFacetboolData.hpp"
#include "CastTo.hpp"
#include "TDFacetboolData.hpp"
#include "CubitFacet.hpp"

TDFacetboolData::TDFacetboolData()
{
  isReversed=false;
}

TDFacetboolData::~TDFacetboolData()
{

}

int TDFacetboolData::is_facetbool_facet(const ToolData* td)
{
  return (CAST_TO(const_cast<ToolData*>(td), TDFacetboolData) != NULL);
}

CubitStatus TDFacetboolData::add_facetbool_facet( FacetEntity *facet_ptr )
{
  TDFacetboolData* td = (TDFacetboolData*) 
                           facet_ptr->get_TD( &TDFacetboolData::is_facetbool_facet );
  if ( td == NULL )
  {
    td = new TDFacetboolData;
    facet_ptr->add_TD( td );
  }

  return CUBIT_SUCCESS;
}

void TDFacetboolData::set(int sv, int e0v, int e1v, int e2v, bool parent, bool is_reversed)
    { surfaceIndex = sv; 
      edgeIndex[0] = e0v; 
      edgeIndex[1] = e1v; 
      edgeIndex[2] = e2v;
      body_1_flag = parent;
      isReversed=is_reversed;
    } 

TDFacetboolData* TDFacetboolData::get(CubitFacet *facet_ptr)
{
  TDFacetboolData *td = (TDFacetboolData*) 
                           facet_ptr->get_TD(&TDFacetboolData::is_facetbool_facet);
  if ( td != NULL )
  {
    return td;
  }
  return (TDFacetboolData*) NULL;
}

//manually set the isReversed flag
void TDFacetboolData::is_reversed(bool rev_flag)
{
  if(rev_flag != isReversed){
    isReversed = rev_flag;
    int temp_int;
    temp_int = edgeIndex[0];
    edgeIndex[0]=edgeIndex[2];
    edgeIndex[2]=temp_int;
  }
}
//get indices.  We force the caller to pass in the rev_flag of the
//associated facet so that we can be sure to return the edges in the
//correct order.
int *TDFacetboolData::get_edge_indices(bool rev_flag)
{
  is_reversed(rev_flag);
  return get_edge_indices();
}
