//- Class:       FacetEntity
//- Description: FacetEntity class - the base class in the Facet Entity Tree.
//- Owner:       Steve Owen
//- Checked by:
//- Version: $Id: 

#ifndef FACETENTITY_HPP
#define FACETENTITY_HPP
 
#include "ToolDataUser.hpp"
#include "DLIList.hpp"
class CubitQuadFacet;
class CubitFacet;
class CubitFacetEdge;
class CubitPoint;

class FacetEntity: public ToolDataUser
{   
private:
   
public:

  FacetEntity();
  virtual ~FacetEntity();
  virtual void get_parents( DLIList<FacetEntity*> &facet_list ) = 0;
  virtual void edges( DLIList<CubitFacetEdge*> &edge_list) = 0;
  virtual void facets( DLIList<CubitFacet*> &facet_list) = 0;
  virtual void points( DLIList<CubitPoint*> &point_list) = 0;
  virtual void debug_draw(int =-1, int = 1, int = 0) = 0;
};

   

#endif // FACETENTITY_HPP


