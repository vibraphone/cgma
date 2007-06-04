#include "CubitFacetData.hpp"
#include "CubitPointData.hpp"
#include "TDGeomFacet.hpp"
#include "TDFacetBoundaryPoint.hpp"
#include "TDFacetBoundaryEdge.hpp"


const int NODE_ALLOC_SIZE           = 1024;
const int TRI_ALLOC_SIZE            = 4*NODE_ALLOC_SIZE;


MemoryManager CubitFacetData::memoryManager("CubitFacetData", sizeof(CubitFacetData),
                                      TRI_ALLOC_SIZE,
                                      STATIC_MEMORY_MANAGER);

MemoryManager CubitPointData::memoryManager("CubitPointData", sizeof(CubitPointData),
                                        NODE_ALLOC_SIZE,
                                      STATIC_MEMORY_MANAGER);

MemoryManager TDGeomFacet::memoryManager("TDGeomFacet", sizeof(TDGeomFacet),
                                        NODE_ALLOC_SIZE,
                                      STATIC_MEMORY_MANAGER);

MemoryManager TDFacetBoundaryEdge::memoryManager("TDFacetBoundaryEdge", sizeof(TDFacetBoundaryEdge),
                                        NODE_ALLOC_SIZE,
                                      STATIC_MEMORY_MANAGER);

MemoryManager TDFacetBoundaryPoint::memoryManager("TDFacetBoundaryPoint", sizeof(TDFacetBoundaryPoint),
                                        NODE_ALLOC_SIZE,
                                      STATIC_MEMORY_MANAGER);


