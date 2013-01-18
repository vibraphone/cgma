#ifndef FACET_SHAPE_DEFS_HPP
#define FACET_SHAPE_DEFS_HPP

#include <vector>

typedef std::vector<double> FacetPointSet;

struct FacetShapes
{};

struct VertexFacets : public FacetShapes
{
  VertexFacets() {}
  int point;
};

struct CurveFacets : public FacetShapes
{
  CurveFacets() {}
  VertexFacets* vertexTopology[2];
  std::vector<int> points; //inclusive of end points
};

struct SurfaceFacets : public FacetShapes
{
  SurfaceFacets() {}
  std::vector<std::vector<std::pair<CurveFacets*, CubitSense> > > curveTopology;
  std::vector<int> facetConnectivity; 
};

struct VolumeFacets : public FacetShapes
{
  VolumeFacets() {}
  std::vector<std::pair<SurfaceFacets*, CubitSense> > surfaceTopology;
};

#endif //FACET_SHAPE_DEFS_HPP
