//- Class:       ChollaPoint
//- Description: Temporary class for constructing the facet-based geometry
//-
//- Owner:       Steven J. Owen
//- Checked by:
//- Version:

#include "ChollaPoint.hpp"
#include "ChollaCurve.hpp"
#include "ChollaSurface.hpp"

//===============================================================================
//Function:  ChollaPoint (PUBLIC) (constructor)
//===============================================================================
ChollaPoint::ChollaPoint( )
{
  static int count = 1; 
  id = count++; 
  myCubitPoint = NULL; 
  myPoint = NULL; 
}
//===============================================================================
//Function:  ~ChollaPoint (PUBLIC) (destructor)
//===============================================================================
ChollaPoint::~ChollaPoint()
{
}

//===============================================================================
//Function:  get_surfaces (PUBLIC) 
//===============================================================================
void ChollaPoint::get_surfaces(DLIList<ChollaSurface *> &surf_list)
{
  int ii, jj;
  ChollaCurve *curv_ptr;
  DLIList<ChollaSurface*> *surf_list_ptr;
  ChollaSurface *surf_ptr;

  curveList.reset();
  for(ii=0; ii<curveList.size(); ii++)
  {
    curv_ptr = curveList.get_and_step();
    surf_list_ptr = curv_ptr->get_surface_list_ptr();
    for(jj=0; jj<surf_list_ptr->size(); jj++)
    {
      surf_ptr = surf_list_ptr->get_and_step();
      surf_list.append_unique( surf_ptr );
    }
  }
}



//EOF

