//- Class: LoopParamTool
//-------------------------------------------------------------------------
// Filename      : LoopParamTool.hpp
//
// Purpose       : Surface Parameterization for mesh by triangulation flattening specific for Planar surfaces
//
// Creator       : Shiraj Khan
//
// Creation Date : 1/21/2003
//
// Owner         : Shiraj Khan
//-------------------------------------------------------------------------

#ifndef LOOP_PARAM_TOOL_HPP
#define LOOP_PARAM_TOOL_HPP

#include "ParamTool.hpp"
#include "CubitVector.hpp"
#include "CubitPoint.hpp"
#include "TDVector.hpp"

class LoopParamTool : public ParamTool
{
public:

    //- constructor
  LoopParamTool();

    //- deconstructor
  ~LoopParamTool();
	
  virtual CubitStatus set_up_space();
  virtual CubitStatus transform_to_uv(const CubitVector &, CubitVector &);
  virtual CubitStatus transform_to_xyz(CubitVector &, const CubitVector &);
    //These functions are required because this is a child of ParamTool.
    //They do nothing.

  CubitStatus transform_loopspoints_to_uv(DLIList<DLIList<CubitPoint *>*> &points);
    ///
    ///This function is the function that changes the CubitPoints locations to the XY plane,
    ///and should be called only after calling new_space_LoopParam.
    ///Also, the original locations of the points are stored on a TDVector on each of the
    ///points.  So to reverse the transformation, the calling code needs to get this
    ///ToolData off the CubitPoint, and reset its location based on the vector stored there.
    ///The calling code is reponsible for removing and deleting the memory of the TDVectors.
    ///
		                                        
  CubitStatus new_space_LoopParam( DLIList<DLIList<CubitPoint *>*> &loop_cubit_points,
                                   CubitVector* normal = NULL);
    ///
    ///This function must be called first for this tool.  This sets up the plane and
    ///the transformation to the XY plane.
    ///
  
		                              
private:
	
    //- reference surface
  CubitVector Du;
  CubitVector Dv;
  CubitVector uvCenter;
  double a, b, c, d;
    //- transformation info - used with set_up_space

  CubitStatus check_selfintersecting_coincident_edges(DLIList<DLIList<CubitPoint *>*> &loop_cubit_points);
    ///
    /// Checks to see if the tranformation of the points caused self intersections.
    /// Returns CUBIT_FAILURE if there are any found.
    ///

  CubitStatus transform_to_bestfit_plane(DLIList<DLIList<CubitPoint *>*> &loop_cubit_points);
    ///
    /// Transforms the points to a least-squares-fit plane.
    /// Creates and stores the original locations on TDVectors that are added to the CubitPoints.
    ///

  CubitStatus transform_to_uv_local(CubitVector &xyz_location, CubitVector &uv_location);
    ///
    /// Transforms the vector to the xy plane.
    ///

  CubitStatus transform_to_xyz_local(CubitVector &xyz_location, CubitVector &uv_location);
    ///
    /// Transforms from the XY plane to the bestfit plane.
    ///
	
  bool double_equal(double val, double equal_to);
    ///
    /// Returns true if val is within roundoff of equal_to, otherwise false.
    ///
};

#endif // LOOP_PARAM_TOOL_HPP

