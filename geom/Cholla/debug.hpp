#ifndef FACETDEBUG_HPP
#define FACETDEBUG_HPP

#include "CubitColorConstants.hpp"

class FacetEntity;
class CubitFacet;
class CubitFacetEdge;
class CubitPoint;
class CubitVector;
class CubitBox;
template <class X> class DLIList;

void dcolor(int icol);

void ddraw( FacetEntity *facet_ptr );

void dfdraw( CubitFacet *facet_ptr );

void dedraw( CubitFacetEdge *facet_ptr );

void dpdraw( CubitPoint *facet_ptr );

void dview();

void dzoom(CubitBox &box);

void dldraw( DLIList<FacetEntity *>&facet_list );

void dfldraw( DLIList<CubitFacet *>&facet_list);

void deldraw( DLIList<CubitFacetEdge *>&edge_list);

void dpldraw( DLIList<CubitPoint *>&point_list);

int dflcheck( DLIList<CubitFacet *>&facet_list);

int dcheck( DLIList<FacetEntity *>&facet_list);

int dfcheck( CubitFacet *facet_ptr );

void dray( const CubitVector &start, const CubitVector &vec, double length=0.1 );

void dpoint( const CubitVector &pt );

int get_color();

#endif


