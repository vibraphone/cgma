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
  myMergePartner = NULL;
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

//===============================================================================
//Function:  is_in_surface (PUBLIC) 
//===============================================================================
CubitBoolean ChollaPoint::is_in_surface( ChollaSurface *cholla_surf )
{
  int ii;
  ChollaCurve *curv_ptr;
  DLIList<ChollaSurface*> *surf_list_ptr;

  curveList.reset();
  for(ii=0; ii<curveList.size(); ii++) 
  {
    curv_ptr = curveList.get_and_step();
    surf_list_ptr = curv_ptr->get_surface_list_ptr(); 

    if( surf_list_ptr->is_in_list( cholla_surf ) )
      return CUBIT_TRUE;
  }

  return CUBIT_FALSE;
}

//===============================================================================
//Function:  is_in_volume (PUBLIC) 
//===============================================================================
CubitBoolean ChollaPoint::is_in_volume( ChollaVolume *cholla_vol )
{
  int ii;
  ChollaCurve *curv_ptr;  

  curveList.reset();
  for(ii=0; ii<curveList.size(); ii++) 
  {
    curv_ptr = curveList.get_and_step();

    if( curv_ptr->is_in_volume( cholla_vol ) )
      return CUBIT_TRUE;
  }

  return CUBIT_FALSE;
}

//===============================================================================
//Function:  is_in_curve (PUBLIC) 
//===============================================================================
CubitBoolean ChollaPoint::is_in_curve( ChollaCurve *chcurve )
{
  for (int ii=0; ii<curveList.size(); ii++)
  {
    ChollaCurve *curv = curveList.get_and_step();
    if (curv == chcurve)
      return CUBIT_TRUE;
  }
  return CUBIT_FALSE;
}

//=============================================================================
//Function:  verify_curves (PUBLIC)
//Description:  verify that all curves at this point have this point as an adjacency
//Notes: 
//=============================================================================
CubitStatus ChollaPoint::verify_curves()
{
  for(int ii=0; ii<curveList.size(); ii++)
  {
    ChollaCurve *crv = curveList.get_and_step();
    if (!crv->has_point(this))
      return CUBIT_FAILURE;
  }
  return CUBIT_SUCCESS;
}

//EOF

