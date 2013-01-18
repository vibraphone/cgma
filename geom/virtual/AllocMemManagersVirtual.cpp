#include "FaceterPointData.hpp"
#include "FaceterFacetData.hpp"
#include "ImprintPointData.hpp"
#include "ImprintLineSegment.hpp"

const int NODE_ALLOC_SIZE           = 1024;
const int TRI_ALLOC_SIZE            = 4*NODE_ALLOC_SIZE;

MemoryManager FaceterFacetData::memoryManager("FaceterFacetData", sizeof(FaceterFacetData),
                                      TRI_ALLOC_SIZE,
                                      STATIC_MEMORY_MANAGER);
MemoryManager FaceterPointData::memoryManager("FaceterPointData", sizeof(FaceterPointData),
                                        NODE_ALLOC_SIZE,
                                      STATIC_MEMORY_MANAGER);
MemoryManager ImprintPointData::memoryManager("ImprintPointData", sizeof(ImprintPointData),
                                              NODE_ALLOC_SIZE,
                                              STATIC_MEMORY_MANAGER);
MemoryManager ImprintLineSegment::memoryManager("ImprintLineSegment", sizeof(ImprintLineSegment),
                                                NODE_ALLOC_SIZE,
                                                STATIC_MEMORY_MANAGER);

