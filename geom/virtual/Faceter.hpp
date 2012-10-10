//-------------------------------------------------------------------------
// Filename      : Faceter.hpp
//
// Purpose       : Facets a surface.
//
// Special Notes : 
//
// Creator       : David White
//
// Creation Date : 03/01/02
//-------------------------------------------------------------------------

#ifndef FACETER_HPP
#define FACETER_HPP

#include "CubitDefines.h"
#include "GeometryDefines.h" //need GEOMETRY_RESABS
#include "DLIList.hpp"
#include "CubitVector.hpp"

class CubitFacet;
class CubitPoint;
class CubitVector;
class FaceterPointData;
class RefFace;
class RefEdge;
class PointGridSearch;
class GMem;


template <class X> class DLIList;

class Faceter
{

  public:

  Faceter( RefFace *face );
  ~Faceter();

  CubitStatus facet_surface( DLIList <CubitFacet*> &results,
                             DLIList <CubitPoint*> &point_list);
    //- Facets the surface and returns the list of facets that
    //- approximate the surface.

  protected:

  private:

  RefFace *thisRefFacePtr;
  PointGridSearch *gridSearchPtr;
  const double gridCellScale;
  CubitBoolean avoidedOverlap;
  DLIList <CubitPoint*> *globalPointList;
//  GeomIntersectionTool* iToolPtr;

  CubitStatus facet_loop(DLIList <CubitPoint*> *loop_ptr,
                         DLIList <CubitFacet*> &results);
    //- Facets' the loop of CubitPoints.

  CubitStatus facet_factory(FaceterPointData* curr_faceter_data,
                            CubitFacet *&new_facet, 
                            DLIList <FaceterPointData*> &order_list);
    //- given the current data, create a new trianglular facet.
    //- Also update the data associated with the curr_data (next and prev data).
  CubitBoolean avoid_facet(FaceterPointData* curr_faceter_data,
                           DLIList<FaceterPointData*> &order_list);
    //- Determine if the facet should be avoided.  Also adds curr_faceter_data
    //- back to the order_list at the begining of the list (since we
    //- pop through the list, this is really the end).  Also modify
    //- the angle by angle + 2pi so we don't hit this point for a while and
    //- hopefully serendiputiously we avoid having to create this facet...
    //- Set the class variable avoidedOverlap to be true if the
    //- facet is avoided.  Returns FALSE if nothing is done and the
    //- facet is okay to create.
                           
  
  CubitStatus get_boundary_points( DLIList <DLIList<CubitPoint*>*> &boundary_point_loops ) const;
    //- Get the loops of boundary facet points.
  
  CubitStatus get_curve_facets( RefEdge* curve, DLIList<CubitPoint*>& segments ) const;
    //- Get the facet points in order from start to end for the curve.

  void max_min_edge_ratio( DLIList <DLIList <CubitPoint*>*> &boundary_point_loops,
	                       double &ratio,
						   double &cell_size);
    //- Gets the ratio for the gridcell size for the gridsearch tool.

  CubitStatus interior_angle(DLIList <CubitPoint*> *loop_ptr,
                             double &my_angle);
    //- Gets the interior angle of the loop_ptr measured at the
    //- current point in the list.  The angle is taken between the prev
    //- point and next point in the list.


};


#endif
