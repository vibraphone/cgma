//- Class: ChollaPoint
//- Owner: Steven J. Owen
//- Description: Maintains a list of mesh.  This is used to store
//-        the exterior skin while creating a geometry from a mesh.
//- Created: 4/27/01
//- Checked By:
//- Version:

#ifndef ChollaPoint_HPP
#define ChollaPoint_HPP

#include "DLIList.hpp"
#include "ChollaEntity.hpp"

class Curve;
class FacetEntity;
class ChollaCurve;
class ChollaSurface;

class ChollaPoint : public ChollaEntity
{
private:

  DLIList<ChollaCurve*> curveList;
  FacetEntity *myCubitPoint;
  void *myPoint;
  int id;

public:

  ChollaPoint( );
    //- default constructor

  ~ChollaPoint();
    // destructor

  void add_facet(FacetEntity *exterior_node)
    {myCubitPoint = exterior_node; }
    //- define the node associated with this point

  FacetEntity *get_facets()
    {return myCubitPoint;}
    //- get the point

  void add_curve( ChollaCurve *fcm_ptr )
    {curveList.append_unique( fcm_ptr );}
    //- associate a curve with this point

  DLIList<ChollaCurve*> &get_curves()
    {return curveList;}
    //- get the list of curves attached to this point

  DLIList<ChollaCurve*> *get_curve_list_ptr()
    {return &curveList;}
    //- get the pointer to the list of curves

  void assign_geometric_point(void *point)
    {myPoint = point;}
    //- set the geometric point associated with the ChollaPoint

  void* get_geometric_point()
    {return myPoint;}
    //- return the geometric point associated with the ChollaPoint

  int get_id(){return id;}

  void get_surfaces(DLIList<ChollaSurface *> &surf_list);
    //- get list of associated cholla surfaces
};

#endif


