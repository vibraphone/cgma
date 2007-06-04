//-------------------------------------------------------------------------
// Filename      : PartitionLumpImprint.hpp
//
// Purpose       : Imprint a lump with a polyline-loop
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 02/05/03
//-------------------------------------------------------------------------

#ifndef PARTITION_LUMP_IMPRINT
#define PARTITION_LUMP_IMPRINT

#include "DLIList.hpp"
#include "RTree.hpp"
#include "CubitVector.hpp"
#include <map>

class PartitionCurve;
class PartitionCoEdge;
class PartitionLoop;
class PartitionSurface;
class PartitionLump;
class PartitionEntity;

class CubitPoint;
class CubitFacet;
class CubitFacetData;
class CubitFacetEdge;
class CubitFacetEdgeData;

class PartitionLumpImprint
{
  public:
  
    PartitionLumpImprint( PartitionLump* lump );
    
    PartitionSurface* imprint( DLIList<CubitFacet*>& facets,
                               DLIList<PartitionEntity*>& new_entities );
    
    PartitionSurface* imprint( DLIList<CubitFacetData*>& facets,
                               DLIList<CubitVector*>& vtx_points,
                               DLIList<PartitionEntity*>& new_entities );
  
  private:
  
    PartitionSurface* imprint( DLIList<CubitVector*>& vtx_points,
                               DLIList<PartitionEntity*>& new_entities );
  
    PartitionLoop* imprint( DLIList<CubitPoint*>& loop,
                            DLIList<CubitPoint*>& vtx_points );
  
    bool add( PartitionEntity* entity );
    
    PartitionEntity* point_owner( CubitPoint* point ) ;
    
    void set_point_owner( CubitPoint* point, PartitionEntity* owner );
    
    void get_owned_points( PartitionEntity*, DLIList<CubitPoint*>& );
    
    inline CubitPoint* point( int point_id ) const
      { return loopPoints.next(point_id); }
    
    int num_points() const
      { return loopPoints.size(); }
    
    void init( DLIList<CubitFacet*>* facets );
    void begin_loop( DLIList<CubitPoint*>& loop );
    
    void clean_up_loop();
    CubitStatus abort_imprint();
    
    CubitStatus do_imprint();
    
    CubitStatus make_vertices( DLIList<CubitPoint*>& vtx_points );
    
    CubitStatus make_volume_curves();
    
    CubitStatus get_curves( DLIList<PartitionCoEdge*>& result_list );
    
    CubitStatus partitionCurve( CubitPoint* );
    CubitStatus partitionSurface( int first_point_id, int last_point_id,
                                  PartitionSurface* surf = 0 );
    CubitStatus makePointCurve( CubitPoint* );
    CubitStatus makeFreePoint( CubitPoint* );
    CubitStatus makeFreeCurve( int start_vert_point, int end_vert_point );

    PartitionEntity* find_closest( const CubitVector& pos,
                                   DLIList<PartitionEntity*>& list,
                                   bool use_tolerance = true );

    CubitStatus seam_curve( DLIList<CubitFacetEdgeData*>& edges,
                            PartitionCurve* curve,
                            DLIList<CubitFacetData*>& facets );

    DLIList<CubitFacetData*> facetList;
    DLIList<CubitFacetEdge*> boundaryEdges;
    DLIList<CubitPoint*> loopPoints;
    std::map<CubitPoint*,PartitionEntity*> pointAssoc;
    
    PartitionLump* lump;
    PartitionSurface* newSurface;
    RTree<PartitionEntity*> rTree;
    DLIList<PartitionEntity*> entityList;
    
    DLIList<PartitionEntity*> newEntities;
};

#endif

