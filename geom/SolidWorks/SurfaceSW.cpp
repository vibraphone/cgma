//-------------------------------------------------------------------------
// Filename      : SurfaceSW.cpp
//
// Purpose       : 
//
// Special Notes :
//
// Creator       : Joel Kopp, Malcolm J. Panthaki
//
// Creation Date : 8/10/00
//
// Owner         : Joel Kopp, Malcolm J. Panthaki
//-------------------------------------------------------------------------


// Precompiled header
#include "stdafx.h"

#include <math.h>

//#include <amapp.h>
#import "C:\Program Files\SolidWorks\sldworks.tlb" raw_interfaces_only, raw_native_types, no_namespace, named_guids     //Import the SolidWorks type library
#include <swconst.h>

#include "CubitVector.hpp"
#include "SurfaceSW.hpp"

#include "CubitUtil.hpp"
#include "SWQueryEngine.hpp"

#include "GMem.hpp"
#include "SWPart.hpp"


extern HRESULT ExtractVARIANTArrayData(VARIANT &v, long lArraySize, double *dArray);



//-------------------------------------------------------------------------
// Purpose       : The constructor
//
// Special Notes :
//
// Creator       : Joel Kopp, Xuechen Liu
//
// Creation Date : 8/10/00
//-------------------------------------------------------------------------
SurfaceSW::SurfaceSW(SWPart *pPart)
    : sense_(CUBIT_FORWARD)
{
	m_pSWFace = NULL;
	m_pSWPart = pPart;
}

//-------------------------------------------------------------------------
// Purpose       : The default destructor. 
//
// Special Notes :
//
// Creator       : Joel Kopp, Raikanta Sahu
//
// Creation Date : 8/10/00
//-------------------------------------------------------------------------
SurfaceSW::~SurfaceSW() 
{
	set_FACE_ptr(NULL);
}

IFace2 *SurfaceSW::get_FACE_ptr() const
{
	return m_pSWFace;
}

void SurfaceSW::set_FACE_ptr(IFace2 *face)
{
	if (face == m_pSWFace)
		return;

	if (m_pSWFace)
	{
		m_pSWPart->remove_cubit_owner(m_pSWFace);
        m_pSWFace->Release();
	}

	m_pSWFace = face;

	if (m_pSWFace)
	{
		m_pSWPart->set_cubit_owner(m_pSWFace, this);
        m_pSWFace->AddRef();
	}

}

//-------------------------------------------------------------------------
// Purpose       : Get geometry modeling engine: SWGeometryEngine
//
// Special Notes :
//
// Creator       : Joel Kopp, Stephen J. Verzi
//
// Creation Date : 8/3/00
//-------------------------------------------------------------------------
GeometryQueryEngine* SurfaceSW::get_geometry_query_engine() const
{
  return SWQueryEngine::instance();   
}                 

//-------------------------------------------------------------------------
// Purpose       : The purpose of this function is to append a
//                 attribute to the GE. The name is attached to the 
//                 underlying solid model entity this one points to.
//
//
// Special Notes : 
//
// Creator       : Joel Kopp, Malcolm J. Panthaki
//
// Creation Date : 8/17/00
//-------------------------------------------------------------------------
void SurfaceSW::append_simple_attribute_virt(CubitSimpleAttrib* csattrib_ptr)
{
    assert(m_pSWFace);
	m_pSWPart->append_simple_attribute_virt(m_pSWFace, csattrib_ptr);
}

//-------------------------------------------------------------------------
// Purpose       : The purpose of this function is to remove a simple 
//                 attribute attached to this geometry entity. The name is 
//                 removed from the underlying BODY this points to.
//
// Special Notes : 
//
// Creator       : Joel Kopp, David R. White
//
// Creation Date : 8/17/00
//-------------------------------------------------------------------------
void SurfaceSW::remove_simple_attribute_virt(CubitSimpleAttrib* csattrib_ptr)
{
    assert(m_pSWFace);
	m_pSWPart->remove_simple_attribute_virt(m_pSWFace, csattrib_ptr);
}

//-------------------------------------------------------------------------
// Purpose       : The purpose of this function is to remove all simple 
//                 attributes attached to this geometry entity.  Also
//                 removes lingering GTC attributes.
//
//
// Special Notes : 
//
// Creator       : Joel Kopp, Greg Nielson
//
// Creation Date : 8/17/00
//-------------------------------------------------------------------------
void SurfaceSW::remove_all_simple_attribute_virt()
{
    assert(m_pSWFace);
	m_pSWPart->remove_all_simple_attribute_virt(m_pSWFace);
}

//-------------------------------------------------------------------------
// Purpose       : The purpose of this function is to get the  
//                 attributes attached to this geometry entity. The name is 
//                 attached to the underlying BODY this points to.
//
// Special Notes : 
//
// Creator       : Joel Kopp, Malcolm J. Panthaki
//
// Creation Date : 8/17/00
//-------------------------------------------------------------------------
CubitStatus SurfaceSW::get_simple_attribute(DLIList<CubitSimpleAttrib*>&
                                           cubit_simple_attrib_list)
{
    assert(m_pSWFace);
	return m_pSWPart->get_simple_attribute(m_pSWFace, cubit_simple_attrib_list);
}
CubitStatus SurfaceSW::get_simple_attribute(const CubitString& name,
                                   DLIList<CubitSimpleAttrib*>& list)
  { return m_pSWPart->get_simple_attribute(m_pSWFace, name, list); }

//-------------------------------------------------------------------------
// Purpose       : Get the bounding box of the object.
//
// Special Notes :
//
// Creator       : Joel Kopp, Raikanta Sahu
//
// Creation Date : 8/17/00
//-------------------------------------------------------------------------
CubitBox SurfaceSW::bounding_box() const 
{
    HRESULT hr = NOERROR;

	// Get the bounding box
    VARIANT vCoords;
    IFace2 *pSWFace = get_FACE_ptr();
    hr = pSWFace->GetBox(&vCoords);

	double SW_box[6];
    hr = ExtractVARIANTArrayData(vCoords, 6, SW_box);
    hr = VariantClear(&vCoords);

    CubitVector corner1(SW_box[0], SW_box[1], SW_box[2]);
    CubitVector corner2(SW_box[3], SW_box[4], SW_box[5]);

    CubitVector asmCorner1;
    CubitVector asmCorner2;
    m_pSWPart->transformPartToAssembly(corner1, asmCorner1);
    m_pSWPart->transformPartToAssembly(corner2, asmCorner2);

	// Convert to a CubitBox and return it
    return CubitBox(asmCorner1, asmCorner2);
}

CubitStatus SurfaceSW::get_point_normal( CubitVector &point, CubitVector &normal )
{
  HRESULT hr = NOERROR;
  VARIANT vParams;
  if( geometry_type() != PLANE_SURFACE_TYPE )
  {
    assert(false);
    return CUBIT_FAILURE;
  }

  IFace2 *face = get_FACE_ptr();

  ISurface *surface = NULL;
  hr = face->IGetSurface ( &surface );

  hr = surface->get_PlaneParams( &vParams );
    surface->Release();

  assert(SUCCEEDED(hr));

  double surfParams[6];
  hr = ExtractVARIANTArrayData(vParams, 6, surfParams);
  hr = VariantClear(&vParams);

  CubitVector partPoint(surfParams[3], surfParams[4], surfParams[5]);
  CubitVector partNormal(surfParams[0], surfParams[1], surfParams[2]);
  CubitVector partOrigin(0.0, 0.0, 0.0);
  CubitVector asmNormal;
  CubitVector asmOrigin;

  m_pSWPart->transformPartToAssembly(partPoint, point);
  m_pSWPart->transformPartToAssembly(partNormal, asmNormal);
  m_pSWPart->transformPartToAssembly(partOrigin, asmOrigin);

  normal = asmNormal - asmOrigin;

  return CUBIT_SUCCESS;
}   

//-------------------------------------------------------------------------
// Purpose       : Computes the closest_point on the surface to the input 
//                 location.  Optionally, it also computes and returns
//                 the normal to the surface and the principal curvatures
//                 at closest_location.
//
// Special Notes : The querying is done on the *first* FACE.
//
//                 If the normal and/or the principal curvatures are 
//                 needed by the calling code, it must allocate space
//                 for these CubitVectors and pass the relevant non_NULL
//                 pointers in.  These are optional.
//
// Creator       : Joel Kopp, Malcolm J. Panthaki
//
// Creation Date : 8/17/00
//-------------------------------------------------------------------------
CubitStatus SurfaceSW::closest_point( CubitVector const& location, 
                                        CubitVector* closest_location,
                                        CubitVector* unit_normal_ptr,
                                        CubitVector* curvature1_ptr,
                                        CubitVector* curvature2_ptr)
{
  VARIANT vSurfData;
  HRESULT hr = NOERROR;
  IFace2 *face = NULL;
  ISurface *surf = NULL;
  CubitVector partOrigin;
  CubitVector asmOrigin;

  // Now compute the point on the FACE closest to the input location,
  // and any other information that is required (normal and/or curvatures)
  // Use the first FACE in the FACE list.
  face = get_FACE_ptr();
  assert(face);

  hr = face->IGetSurface(&surf);
  assert(surf);

  // get closest point on surface - x,y,z and u,v
  CubitVector tempLocation = location;
  CubitVector partLocation;
  m_pSWPart->transformAssemblyToPart(tempLocation, partLocation);
  hr = surf->GetClosestPointOn ( partLocation.x(), partLocation.y(), 
                partLocation.z(), &vSurfData );

  double closest_point[5];
  hr = ExtractVARIANTArrayData(vSurfData, 5, closest_point);
  hr = VariantClear(&vSurfData);

  // Fill in the closest_location object if it was passed in
  if (closest_location)
  {
    CubitVector partClosest(closest_point[0], closest_point[1], closest_point[2]);
    m_pSWPart->transformPartToAssembly(partClosest, *closest_location);
  }

  // compute the transformed part origin in the assembly coord sys
  // - this is needed later for transforming any vector output
  partOrigin.set(0.0, 0.0, 0.0);
  m_pSWPart->transformPartToAssembly(partOrigin, asmOrigin);

    // Need to compute the normal
  if (unit_normal_ptr != NULL)
  {
    double u, v;
    u = closest_point[3];
    v = closest_point[4];

    hr = surf->Evaluate(u, v, 0, 0, &vSurfData); // no u or v derivatives

    double dEvaluate[6];
    hr = ExtractVARIANTArrayData(vSurfData, 6, dEvaluate);
    hr = VariantClear(&vSurfData);

    // Reverse the normal, if the sense of the FACE is reversed wrt the 
    // underlying SW surface object.
    CubitVector partUnitNormal;
    if ( get_relative_surface_sense() == CUBIT_REVERSED )
    {
      partUnitNormal.set(-dEvaluate[3], -dEvaluate[4], -dEvaluate[5]);
    }
    else
    {
      partUnitNormal.set(dEvaluate[3], dEvaluate[4], dEvaluate[5]);
    }

    CubitVector asmUnitNormal;
    m_pSWPart->transformPartToAssembly(partUnitNormal, asmUnitNormal);
    *unit_normal_ptr = asmUnitNormal - asmOrigin;
  }

  // need curvature
  if (curvature1_ptr || curvature2_ptr)
  {
    hr = surf->EvaluateAtPoint ( closest_point[0], closest_point[1],
                                      closest_point[2], &vSurfData );

    double surfEval[11];
    hr = ExtractVARIANTArrayData(vSurfData, 11, surfEval);
    hr = VariantClear(&vSurfData);
	
	  // Fill in the first principal curvature, if necessary
	  if (curvature1_ptr != NULL)
	  {
      CubitVector partCurvature1;
      // First set curvature1 using the direction of this curvature
      partCurvature1.set( surfEval[3], surfEval[4], surfEval[5] );
      if( !(fabs(surfEval[3]) < CUBIT_RESABS &&
            fabs(surfEval[4]) < CUBIT_RESABS && 
            fabs(surfEval[5]) < CUBIT_RESABS ) )
      {
        // Now multiply it with the magnitude of the curvature
        partCurvature1 = (surfEval[9]) * partCurvature1;
        // Reverse the normal, if the sense of the FACE is reversed wrt the 
        // underlying SW surface object.
      }
      if ( get_relative_surface_sense() == CUBIT_REVERSED )
      {
	      partCurvature1 = partCurvature1 * -1.0;
      }

      CubitVector asmCurvature1;
      m_pSWPart->transformPartToAssembly(partCurvature1, asmCurvature1);
      *curvature1_ptr = asmCurvature1 - asmOrigin;
    }

    // Fill in the second principal curvature, if necessary
    if (curvature2_ptr != NULL)
    {
      CubitVector partCurvature2;
	      // First set curvature2 using the direction of this curvature
      partCurvature2.set( surfEval[6], surfEval[7], surfEval[8] );

      // Now multiply it with the magnitude of the curvature
      if( !(fabs(surfEval[6]) < CUBIT_RESABS &&
            fabs(surfEval[7]) < CUBIT_RESABS && 
            fabs(surfEval[8]) < CUBIT_RESABS ) )
      {
        // Now multiply it with the magnitude of the curvature
        partCurvature2 = (surfEval[10]) * partCurvature2;
        // Reverse the normal, if the sense of the FACE is reversed wrt the 
        // underlying SW surface object.
      }
      if ( get_relative_surface_sense() == CUBIT_REVERSED )
      {
        partCurvature2 = partCurvature2 * -1.0;
      }

      CubitVector asmCurvature2;
      m_pSWPart->transformPartToAssembly(partCurvature2, asmCurvature2);
      *curvature2_ptr = asmCurvature2 - asmOrigin;
    }
  }

  surf->Release();
  return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Purpose       : Computes the closest_point on the trimmed surface to the 
//                 input location. 
//
// Special Notes : 
// Creator       : Joel Kopp, Brett W. Clark
//
// Creation Date : 8/18/00
//-------------------------------------------------------------------------
void SurfaceSW::closest_point_trimmed( CubitVector from_point, 
                                         CubitVector& point_on_surface)
{
  HRESULT hr = NOERROR;

  // Get the FACE to be used for computing the closest location.
  IFace2 *face = get_FACE_ptr();

  CubitVector partFromPoint;
  m_pSWPart->transformAssemblyToPart(from_point, partFromPoint);

  VARIANT vPoint;
  hr = face->GetClosestPointOn ( partFromPoint.x(), partFromPoint.y(), 
                                 partFromPoint.z(), &vPoint );
  if ( !SUCCEEDED(hr) )
  {
    PRINT_ERROR("Problems with SW closest point trimmed.\n");
  }

  double closest_point[5];
  hr = ExtractVARIANTArrayData(vPoint, 5, closest_point);
  hr = VariantClear(&vPoint);

  CubitVector partClosest(closest_point[0], closest_point[1], closest_point[2]);
  m_pSWPart->transformPartToAssembly(partClosest, point_on_surface);
}

//-------------------------------------------------------------------------
// Purpose       : This functions computes the point on the surface that is 
//                 closest to the input location and then calculates the 
//                 magnitudes of the principal curvatures at this (possibly, 
//                 new) point on the surface. Specifying the RefVolume for 
//                 reference is optional.
//
// Special Notes :
//
// Creator       : Joel Kopp, Malcolm J. Panthaki
//
// Creation Date : 8/18/00
//-------------------------------------------------------------------------
CubitStatus SurfaceSW::principal_curvatures(
	CubitVector const& location, 
	double& curvature_1,
	double& curvature_2,
	CubitVector* closest_location )
{
	// Get the principal curvature vectors
	CubitVector curvature1_vector;
	CubitVector curvature2_vector;
	CubitStatus result = this->closest_point( location,
											closest_location,
											(CubitVector *)NULL,
											&curvature1_vector,
											&curvature2_vector );

	if ( result == CUBIT_FAILURE )
	{
		PRINT_ERROR("In SurfaceSW::principal_curvatures\n"
					"       Could not compute the principal curvature vectors.\n");
		return CUBIT_FAILURE;
	}

	// Extract the magnitudes of the principal curvatures
	curvature_1 = curvature1_vector.length();
	curvature_2 = curvature2_vector.length();

	return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Purpose       : Given values of the two parameters, get the position.
//
// Special Notes :
//
// Creator       : Joel Kopp, Raikanta Sahu
//
// Creation Date : 8/18/00
//-------------------------------------------------------------------------
CubitVector SurfaceSW::position_from_u_v (double u, double v)
{
  HRESULT hr = NOERROR;

  // Get the first FACE
  IFace2 *face = get_FACE_ptr();
  ISurface *surface = NULL;
  hr = face->IGetSurface ( &surface );

  // Use the parametric position to get the position of the 
  // parametric point
  VARIANT vEval;
  hr = surface->Evaluate ( u, v, 0, 0, &vEval );
  surface->Release();
  assert(SUCCEEDED(hr));

  double surfEval[6];
  hr = ExtractVARIANTArrayData(vEval, 6, surfEval);
  hr = VariantClear(&vEval);

  CubitVector partPosition(surfEval[0], surfEval[1], surfEval[2]);
  CubitVector asmPosition;
  m_pSWPart->transformPartToAssembly(partPosition, asmPosition);
  return asmPosition;
}

//-------------------------------------------------------------------------
// Purpose       : This function returns the {u, v} coordinates of the point 
//                 on the Surface closest to the input point (specified in 
//                 global space). The closest_location is also returned.
//
// Special Notes :
//
// Creator       : Joel Kopp, Malcolm J. Panthaki
//
// Creation Date : 8/18/00
//-------------------------------------------------------------------------
CubitStatus SurfaceSW::u_v_from_position (CubitVector const& location,
                                            double& u, 
                                            double& v,
                                            CubitVector* closest_location )
{
  HRESULT hr = NOERROR;

  // A RefVolume was provided for reference and the calling code needs 
  // at least one of the following to be generated: normal or the
  // principal curvatures (i.e., at least one of those pointers is
  // non-NULL)
  IFace2 *face = get_FACE_ptr();

  // If no valid FACE was found, crash and burn :)
  if (face == NULL)
  {
    PRINT_ERROR("In SurfaceSW::u_v_from_position\n"
                "  THIS IS A BUG - PLEASE REPORT IT!\n");
    return CUBIT_FAILURE;
  }

  // Get the parameter values of the closest point, based on the FACE 
  // just extracted, above
  ISurface *surface = NULL;
  hr = face->IGetSurface ( &surface );

  CubitVector tempLocation(location);
  CubitVector partLocation;
  m_pSWPart->transformAssemblyToPart(tempLocation, partLocation);

  VARIANT vPoint;
  hr = surface->GetClosestPointOn(partLocation.x(), partLocation.y(), partLocation.z(),
                                  &vPoint);
  surface->Release();

  double closest_point[5];
  hr = ExtractVARIANTArrayData(vPoint, 5, closest_point);
  hr = VariantClear(&vPoint);

  u = closest_point[3];
  v = closest_point[4];

  if ( closest_location )
  {
    CubitVector partClosest(closest_point[0], closest_point[1], closest_point[2]);
    m_pSWPart->transformPartToAssembly(partClosest, *closest_location);
  }

  return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Purpose       : Determines whether the SW surface object associated
//                 with one of the FACEs of this SurfaceSW object is 
//                 periodic or not.
//
//                 The querying is done on the *first* FACE.
//
// Creator       : Joel Kopp, Malcolm J. Panthaki
//
// Creation Date : 8/18/00
//-------------------------------------------------------------------------
CubitBoolean SurfaceSW::is_periodic()
{
    CubitBoolean bPeriodicU;
    CubitBoolean bPeriodicV;
    double dPeriodU;
    double dPeriodV;

    periodicInfo(bPeriodicU, dPeriodU, bPeriodicV, dPeriodV);

    if ((CUBIT_TRUE == bPeriodicU) || (CUBIT_TRUE == bPeriodicV))
        return CUBIT_TRUE;
    else
        return CUBIT_FALSE;
}

//
// utility routine to determine periodic properties of a face
//
void SurfaceSW::periodicInfo(CubitBoolean &bPeriodicU, double &dPeriodU,
                             CubitBoolean &bPeriodicV, double &dPeriodV)
{
  HRESULT hr = NOERROR;

  dPeriodU = 0.0;
  dPeriodV = 0.0;

  // Get the first FACE in the list of FACEs associated with this
  // SurfaceSW object
  IFace2 *face = get_FACE_ptr();

    // get the number of loops on this face
    // if the number is less than 2, the face is not periodic
// NOTE - I am basing this loop count check on the code in SurfMapTool::crack_periodic_surface 
// which crashes if there
  long nLoops = 0;
  hr = face->GetLoopCount(&nLoops);
  if (2 > nLoops)
  {
    bPeriodicU = CUBIT_FALSE;
    bPeriodicV = CUBIT_FALSE;
    return;
  }

  // get the uv bounds of the face
  // NOTE - cubit only cares whether the face is periodic.  If the surface is periodic
  // but the face is not (eg. partial cylinder) we need to return false
  VARIANT vBounds;
  hr = face->GetUVBounds(&vBounds);

  double dFaceUVBounds[4];
  hr = ExtractVARIANTArrayData(vBounds, 4, dFaceUVBounds);
  hr = VariantClear(&vBounds);

  // find the face range in u and v
  double uRange = dFaceUVBounds[1] - dFaceUVBounds[0];
  double vRange = dFaceUVBounds[3] - dFaceUVBounds[2];

  // Get corresponding surface for current face
  ISurface *surface = NULL;
  hr = face->IGetSurface ( &surface );

  // Now ask SW whether its underlying surface is periodic
  VARIANT vParams;
  hr = surface->Parameterization ( &vParams );
  surface->Release();
  assert(SUCCEEDED(hr));

  double surfParams[11];
  hr = ExtractVARIANTArrayData(vParams, 11, surfParams);
  hr = VariantClear(&vParams);

  union PackedInts        //unpack packed integers
  {
    double value;
    int intData[2];
  } packedIntData;

  double dDiff;

  // check each parameter direction to see if it is periodic
  packedIntData.value = surfParams[4]; // u direction
  int isPeriodic_start = packedIntData.intData[0];
  int isPeriodic_end = packedIntData.intData[1];

  if (isPeriodic_start == 13701 || isPeriodic_start == 13736 ||
      isPeriodic_end == 13701 || isPeriodic_end == 13736)
  {
    // the underlying surface is periodic, but check to see if the face
    // spans the period
    double period_U = surfParams[1] - surfParams[0];

//TODO - what is the correct tolerance to use here?
    dDiff = uRange - period_U;
    if (fabs(dDiff) > CUBIT_RESABS)
    {
      bPeriodicU = CUBIT_FALSE;
    }
    else
    {
      bPeriodicU = CUBIT_TRUE;
      dPeriodU = period_U;
    }
	}
  else
  {
    bPeriodicU = CUBIT_FALSE;
  }

  packedIntData.value = surfParams[5]; // v direction
  isPeriodic_start = packedIntData.intData[0];
  isPeriodic_end = packedIntData.intData[1];

  if (isPeriodic_start == 13701 || isPeriodic_start == 13736 ||
      isPeriodic_end == 13701 || isPeriodic_end == 13736)
  {
    double period_V = surfParams[3] - surfParams[2];
    
    // the underlying surface is periodic, but check to see if the face
    // spans the period

//TODO - what is the correct tolerance to use here?
    dDiff = vRange - period_V;
    if (fabs(dDiff) > CUBIT_RESABS)
    {
      bPeriodicV = CUBIT_FALSE;
    }
    else
    {
      bPeriodicV = CUBIT_TRUE;
      dPeriodV = period_V;
    }
  }
  else
  {
    bPeriodicV = CUBIT_FALSE;
  }

  return;
}

//-------------------------------------------------------------------------
// Purpose       : Determines whether the SW surface object associated
//                 with one of the FACEs of this SurfaceSW object is 
//                 periodic in the U direction or not.  If it is, it
//                 returns CUBIT_TRUE and the value of the period. Otherwise,
//                 it returns CUBIT_FALSE and a value of 0.0 or the period.
//
//                 The querying is done on the *first* FACE.
//
// Creator       : Joel Kopp, Malcolm J. Panthaki
//
// Creation Date : 8/18/00
//-------------------------------------------------------------------------
CubitBoolean SurfaceSW::is_periodic_in_U( double& period_U ) 
{
    CubitBoolean bPeriodicU;
    CubitBoolean bPeriodicV;
    double dPeriodU;
    double dPeriodV;

    periodicInfo(bPeriodicU, dPeriodU, bPeriodicV, dPeriodV);

    if (bPeriodicU)
    {
        period_U = dPeriodU;
        return CUBIT_TRUE;
    }
    else
    {
        period_U = 0.0;
        return CUBIT_FALSE;
    }
}

//-------------------------------------------------------------------------
// Purpose       : Determines whether the SW surface object associated
//                 with one of the FACEs of this SurfaceSW object is 
//                 periodic in the V direction or not. If it is, it
//                 returns CUBIT_TRUE and the value of the period. Otherwise,
//                 it returns CUBIT_FALSE and a value of 0.0 or the period.
//
//                 The querying is done on the *first* FACE.
//
// Creator       : Joel Kopp, Malcolm J. Panthaki
//
// Creation Date : 8/21/00
//-------------------------------------------------------------------------
CubitBoolean SurfaceSW::is_periodic_in_V( double& period_V ) 
{
    CubitBoolean bPeriodicU;
    CubitBoolean bPeriodicV;
    double dPeriodU;
    double dPeriodV;

    periodicInfo(bPeriodicU, dPeriodU, bPeriodicV, dPeriodV);

    if (bPeriodicV)
    {
        period_V = dPeriodV;
        return CUBIT_TRUE;
    }
    else
    {
        period_V = 0.0;
        return CUBIT_FALSE;
    }
}

//-------------------------------------------------------------------------
// Purpose       : Determines whether the SW surface object associated
//                 with one of the FACEs of this SurfaceSW object is 
//                 singular in the U direction or not at a given u parameter.
//
// Note          : The assumption is made that the u_param is in the
//                  parameter bounds of the surface.
//
// Creator       : Joel Kopp, David White
//
// Creation Date : 8/21/00
//-------------------------------------------------------------------------
CubitBoolean SurfaceSW::is_singular_in_U( double u_param )
{
    assert(false);
	PRINT_ERROR("Function Currently Unavailable.\n");
	return CUBIT_FALSE;
}  

//-------------------------------------------------------------------------
// Purpose       : Determines whether the SW surface object associated
//                 with one of the FACEs of this SurfaceSW object is 
//                 singular in the V direction or not at a given v parameter.
//
// Creator       : Joel Kopp, David White
//
// Creation Date : 8/21/00
//-------------------------------------------------------------------------
CubitBoolean SurfaceSW::is_singular_in_V( double v_param )
{
    assert(false);
	PRINT_ERROR("Function Currently Unavailable.\n");
	return CUBIT_FALSE;
}

//-------------------------------------------------------------------------
// Purpose       : Determines whether the SW surface object associated
//                 with one of the FACEs of this SurfaceSW object is 
//                 closed in the U direction or not.
//
// Creator       : Joel Kopp, David White
//
// Creation Date : 8/21/00
//-------------------------------------------------------------------------
CubitBoolean SurfaceSW::is_closed_in_U()
{
    assert(false);
	PRINT_ERROR("Function Currently Unavailable.\n");
	return CUBIT_FALSE;
}

//-------------------------------------------------------------------------
// Purpose       : Determines whether the SW surface object associated
//                 with one of the FACEs of this SurfaceSW object is 
//                 closed in the V direction or not.
//
// Creator       : Joel Kopp, David White
//
// Creation Date : 8/21/00
//-------------------------------------------------------------------------
CubitBoolean SurfaceSW::is_closed_in_V()
{
    assert(false);
	PRINT_ERROR("Function Currently Unavailable.\n");
	return CUBIT_FALSE;
}

//-------------------------------------------------------------------------
// Purpose       : Calculates the derivitives at a given parameter location.
//
// Creator       : Joel Kopp, David White
//
// Creation Date : 8/21/00
//-------------------------------------------------------------------------
CubitStatus SurfaceSW::uv_derivitives( double u_param,
                                         double v_param,
                                         CubitVector &du,
                                         CubitVector &dv )
{
  //Get the first face and work of it...
  IFace2 *face = get_FACE_ptr();

  ISurface *surface = NULL;
  HRESULT hr = face->IGetSurface ( &surface );
  assert(SUCCEEDED(hr));

  VARIANT vEval;
  hr = surface->Evaluate ( u_param, v_param, 1, 1, &vEval ); // 1's represent # of derivatives
  surface->Release();
  assert(SUCCEEDED(hr));

  double surfEval[15];
  hr = ExtractVARIANTArrayData(vEval, 15, surfEval);
  hr = VariantClear(&vEval);

  CubitVector partDU(surfEval[3], surfEval[4], surfEval[5]);
  CubitVector partDV(surfEval[6], surfEval[7], surfEval[8]);
  CubitVector asmDU;
  CubitVector asmDV;
  CubitVector asmOrigin;

  m_pSWPart->transformPartToAssembly(partDU, asmDU);
  m_pSWPart->transformPartToAssembly(partDV, asmDV);
  m_pSWPart->transformPartToAssembly(CubitVector(0.0, 0.0, 0.0), asmOrigin);
  du = asmDU - asmOrigin;
  dv = asmDV - asmOrigin;

  return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Purpose       : Determines whether the SW surface object associated
//                 with one of the FACEs of this SurfaceSW object is 
//                 parametrically defined or not.
//
//                 The querying is done on the *first* FACE.
//
// Creator       : Joel Kopp, Malcolm J. Panthaki
//
// Creation Date : 8/21/00
//-------------------------------------------------------------------------
CubitBoolean SurfaceSW::is_parametric() 
{
    // all SolidWorks surfaces can be evaluated with parameters
    return CUBIT_TRUE;
}

//-------------------------------------------------------------------------
// Purpose       : Returns the lower and upper parametric bounds of the 
//                 surface in U, if it is parametric.  Otherwise, it returns
//                 CUBIT_FALSE and zeroes for the upper and lower parametric
//                 bounds.
//
//                 The querying is done on the *first* FACE.
//
// Creator       : Joel Kopp, Malcolm J. Panthaki
//
// Creation Date : 8/21/00
//-------------------------------------------------------------------------
CubitBoolean SurfaceSW::get_param_range_U( double& lower_bound,
                                             double& upper_bound )
{
  HRESULT hr = NOERROR;

  // Get the first FACE in the list of FACEs associated with this
  // SurfaceSW object
  IFace2 *face = get_FACE_ptr();

  VARIANT vBounds;
  hr = face->GetUVBounds(&vBounds);
// NOTE - calling surface->Parameterization instead of face->GetUVBounds returns the
// parameter range of the underlying surface.  This can be very large (eg. cylinder with
// infinite extent) and cause problems in calculations.  Byron - 7/2/01

  double dBounds[4];
  hr = ExtractVARIANTArrayData(vBounds, 4, dBounds);
  hr = VariantClear(&vBounds);

  lower_bound = dBounds[0];
  upper_bound = dBounds[1];

  return CUBIT_TRUE;
}

//-------------------------------------------------------------------------
// Purpose       : Returns the lower and upper parametric bounds of the 
//                 surface in V, if it is parametric.  Otherwise, it returns
//                 CUBIT_FALSE and zeroes for the upper and lower parametric
//                 bounds.
//
//                 The querying is done on the *first* FACE.
//
// Creator       : Joel Kopp, Malcolm J. Panthaki
//
// Creation Date : 12/17/96
//-------------------------------------------------------------------------
CubitBoolean SurfaceSW::get_param_range_V( double& lower_bound,
                                             double& upper_bound )
{
  HRESULT hr = NOERROR;

  // Get the first FACE in the list of FACEs associated with this
  // SurfaceSW object
  IFace2 *face = get_FACE_ptr();

  VARIANT vBounds;
  hr = face->GetUVBounds(&vBounds);
// NOTE - calling surface->Parameterization instead of face->GetUVBounds returns the
// parameter range of the underlying surface.  This can be very large (eg. cylinder with
// infinite extent) and cause problems in calculations.  Byron - 7/2/01

  double dBounds[4];
  hr = ExtractVARIANTArrayData(vBounds, 4, dBounds);
  hr = VariantClear(&vBounds);

  lower_bound = dBounds[2];
  upper_bound = dBounds[3];

  return CUBIT_TRUE;
}

//-------------------------------------------------------------------------
// Purpose       : Returns a surface type ID -- the values of these are
//                 determined by SW.
//
//                 The querying is done on the *first* FACE.
//
// Creator       : Joel Kopp, Malcolm J. Panthaki
//
// Creation Date : 8/21/00
//-------------------------------------------------------------------------
GeometryType SurfaceSW::geometry_type()
{
  HRESULT hr = NOERROR;

  // Get the first FACE associated with this SurfaceSW object
  IFace2 *face = get_FACE_ptr();

  ISurface *surface = NULL;
  hr = face->IGetSurface ( &surface );

  // Ask SW for the type of the underlying surface of that FACE
  long surface_type = 0;
  hr = surface->Identity ( &surface_type );
  surface->Release();

  GeometryType local_type;

  switch (surface_type)
  {
  case CONE_TYPE:
    local_type = CONE_SURFACE_TYPE;
    break;
  case PLANE_TYPE:
    local_type = PLANE_SURFACE_TYPE;
    break;
  case SPHERE_TYPE:
    local_type = SPHERE_SURFACE_TYPE;
    break;
  //case CYLINDER_TYPE:
    //local_type = CYLINDER_SURFACE_TYPE;  // no such type in CUBIT GeometryDefines.h
    //break;
  case BSURF_TYPE:
    local_type = SPLINE_SURFACE_TYPE;
    break;
  case TORUS_TYPE:
    local_type = TORUS_SURFACE_TYPE;
    break;
  case BLEND_TYPE:
			local_type = BEST_FIT_SURFACE_TYPE;
			break;
		/*case OFFSET_TYPE:
			local_type = OFFSET_SURFACE_TYPE;    // no such type in CUBIT GeometryDefines.h
			break;*/
		/*case EXTRU_TYPE:
			local_type = EXTRU_SURFACE_TYPE;    // no such type in CUBIT GeometryDefines.h
			break;*/
		/*case SREV_TYPE:
			local_type = SREV_SURFACE_TYPE;    // no such type in CUBIT GeometryDefines.h
			break;*/
  default:
            //assert(false);
			local_type = UNDEFINED_SURFACE_TYPE;
			break;
	}

	return local_type;
}

//-------------------------------------------------------------------------
// Purpose       :This function will compare the alignment of the first and
//                second SW FACES.  This is done by computing the normals
//                for the faces at a point determined by the bounding box
//                of the first FACE. 
//                Once the normals are found, they are put into the
//                correct sense orientation. They are then dotted to
//                determine the sense. This function is currently used
//                for determining if the sense_ should change when the
//                first_face is removed. 
//
// Special Notes: The center point is used because of two assumptions:
//                1) the two faces are spacially equivalent, so the center 
//                of one will be near the center of the other.
//                2) the point_normal function takes the input point and
//                finds the closes point on the actual surface to do the 
//                normal calculation.
//
// Creator       : Joel Kopp, David White
//
// Creation Date : 8/21/00
//-------------------------------------------------------------------------
CubitStatus SurfaceSW::compare_alignment( IFace2 *first_face,
                                          IFace2 *second_face,
                                          CubitSense &sense )
{
    assert(false);
    return CUBIT_FAILURE;
}

//-------------------------------------------------------------------------
// Purpose       : Returns the area of the Surface
//
// Special Notes :
//
// Creator       : Joel Kopp, Raikanta Sahu
//
// Creation Date : 8/22/00
//-------------------------------------------------------------------------
double SurfaceSW::measure() 
{
    HRESULT hr = NOERROR;

	// Get a FACE from the list of FACES
	IFace2 *face = get_FACE_ptr();

	double area = 0.0;
	hr = face->GetArea ( &area );
		
	return area;
}

//-------------------------------------------------------------------------
// Purpose       : This function tests the passed in position to see if
//                 it's on the underlying surface.
//
// Special Notes :
//
// Creator       : Joel Kopp, David R. White
//
// Creation Date : 8/22/00
//-------------------------------------------------------------------------
CubitBoolean SurfaceSW::is_position_on( CubitVector &test_position )
{
    // transform test position from assembly to part coordinate system
    CubitVector partTestPoint;
    m_pSWPart->transformAssemblyToPart(test_position, partTestPoint);

    // get the projection of the point on the face
	IFace2 *face = get_FACE_ptr();

    VARIANT vClose;
	HRESULT hr = face->GetClosestPointOn ( test_position.x(), test_position.y(),
                                            test_position.z(), &vClose );
    assert(SUCCEEDED(hr));

	double closePoint[5];
    hr = ::ExtractVARIANTArrayData(vClose, 5, closePoint);
    hr = VariantClear(&vClose);

	CubitVector partClosePoint (closePoint[0], closePoint[1], closePoint[2]);

	if ( partTestPoint.within_tolerance(partClosePoint, GEOMETRY_RESABS) )
	{
		return CUBIT_TRUE;
	}
	else
	{
		return CUBIT_FALSE;
	}
}

CubitSense SurfaceSW::get_geometry_sense()
{
	CubitSense sense = get_relative_surface_sense();

	return sense;
}

void SurfaceSW::bodysms(DLIList<BodySM*> &bodies) 
{
	if (m_pSWFace)
		SWQueryEngine::instance()->bodysms(m_pSWFace, bodies);
}

void SurfaceSW::lumps(DLIList<Lump*> &lumps)
{
	if (m_pSWFace)
		SWQueryEngine::instance()->lumps(m_pSWPart, m_pSWFace, lumps);
}

void SurfaceSW::shellsms(DLIList<ShellSM*> &shellsms)
{
	if (m_pSWFace)
		SWQueryEngine::instance()->shellsms(m_pSWPart, m_pSWFace, shellsms);
}

void SurfaceSW::surfaces(DLIList<Surface*> &surfaces)
{
	if (m_pSWFace)
		SWQueryEngine::instance()->surfaces(m_pSWPart, m_pSWFace, surfaces);
}

void SurfaceSW::loopsms(DLIList<LoopSM*> &loopsms)
{
	if (m_pSWFace)
		SWQueryEngine::instance()->loopsms(m_pSWPart, m_pSWFace, loopsms);
}

void SurfaceSW::curves(DLIList<Curve*> &curves)
{
	if (m_pSWFace)
		SWQueryEngine::instance()->curves(m_pSWPart, m_pSWFace, curves);
}

void SurfaceSW::coedgesms(DLIList<CoEdgeSM*> &coedgesms)
{
	if (m_pSWFace)
		SWQueryEngine::instance()->coedgesms(m_pSWPart, m_pSWFace, coedgesms);
}

void SurfaceSW::points(DLIList<Point*> &points)
{
	if (m_pSWFace)
        SWQueryEngine::instance()->points(m_pSWPart, m_pSWFace, points);
}

//-------------------------------------------------------------------------
// Purpose       : Return the sense of the first FACE (in the list of FACEs
//                 in this Surface) wrt its underlying SW surface.
//
// Special Notes : 
//
// Creator       : Joel Kopp, Malcolm J. Panthaki
//
// Creation Date : 8/22/00
//-------------------------------------------------------------------------
CubitSense SurfaceSW::get_FACE_sense()
{
    HRESULT hr = NOERROR;

	// Get the first FACE
	IFace2 *face = this->get_FACE_ptr();

	// Return the sense value
	VARIANT_BOOL faceSense = FALSE;
	hr = face->FaceInSurfaceSense ( &faceSense );
	
	if (faceSense == FALSE)
	{
		return CUBIT_FORWARD;
	}
	else if (faceSense == TRUE)
	{
		return CUBIT_REVERSED;
	}
	else
	{
		PRINT_ERROR("In SurfaceSW::get_FACE_sense\n"
					"       Could not get the sense of the FACE "
					"wrt its surface.\n");
		return CUBIT_FORWARD;
	}
}

CubitSense SurfaceSW::get_relative_surface_sense()
{
	CubitSense FACEsense = this->get_FACE_sense();
	CubitSense sense = CUBIT_FORWARD;
	
	if ((FACEsense == CUBIT_REVERSED && sense_ == CUBIT_FORWARD) ||
	  (FACEsense == CUBIT_FORWARD && sense_ == CUBIT_REVERSED))
	{
		sense = CUBIT_REVERSED;
	}

	return sense;
}

void SurfaceSW::reverse_sense()
{
	sense_ = CubitUtil::opposite_sense( sense_ );
}

CubitPointContainment
SurfaceSW::point_containment( const CubitVector &point )
{
    assert(false);
    return CUBIT_PNT_UNKNOWN;
}

CubitPointContainment
SurfaceSW::point_containment( double u, double v )
{
    assert(false);
    return CUBIT_PNT_UNKNOWN;
}

CubitPointContainment
SurfaceSW::point_containment( CubitVector &point, double u, double v )
{
    assert(false);
    return CUBIT_PNT_UNKNOWN;
}

CubitStatus SurfaceSW::facet_face(
	int& number_triangles, int& number_points, int& facet_list_size,
	GMem* g_mem, unsigned short normal_tolerance, float distance_tolerance)
{
    HRESULT hr = NOERROR;

    IFace2 *pSWFace = get_FACE_ptr();
    assert(pSWFace);

	// Because this may be unnecessarily called twice,
	// say there is one triangle.
	if (!g_mem)
	{
		number_triangles = 1;
		number_points = 3;
		facet_list_size = 4;
		return CUBIT_SUCCESS;
	}

	// First, tell the mesh manager about g_mem
	//facetManager->set_gmem(g_mem);

	// Update the default refinement to reflect tolerances
	long numTol = 2;
	long toleranceType[2];
	double toleranceVal[2];
    double oldToleranceVal[2];
	toleranceType[0] = swSurfChordTessellationTol;
	toleranceType[1] = swSurfAngularTessellationTol;
	toleranceVal[0] = (double)normal_tolerance; // TODO - are these the correct tolerances in the correct order?
	toleranceVal[1] = (double)distance_tolerance;
	VARIANT_BOOL res = FALSE;

	IModeler *modeler;
	hr = SWQueryEngine::instance()->GetSldWorks()->IGetModeler(&modeler);

    hr = modeler->GetToleranceValue(toleranceType[0], &oldToleranceVal[0]);
    hr = modeler->GetToleranceValue(toleranceType[1], &oldToleranceVal[1]);

    // TODO - use tolerances in faceting faces

//    double dTemp;
//	hr = modeler->SetToleranceValue(toleranceType[0], toleranceVal[0], &dTemp);
//	hr = modeler->SetToleranceValue(toleranceType[1], toleranceVal[1], &dTemp);

	// Now facet the face (this will add the facets to g_mem).
	long tricount = 0;
	hr = pSWFace->GetTessTriangleCount ( &tricount );
	long pointcount = tricount*3;
	long coordcount = tricount*9;

    VARIANT vTessTris;
	hr = pSWFace->GetTessTriangles ( TRUE, &vTessTris );

//	hr = modeler->SetToleranceValue(toleranceType[0], oldToleranceVal[0], &dTemp);
//	hr = modeler->SetToleranceValue(toleranceType[1], oldToleranceVal[1], &dTemp);

    SAFEARRAY* psa = V_ARRAY(&vTessTris);
    float *dTessPoints;
    hr = SafeArrayAccessData(psa, (void**)&dTessPoints);

    long lUBound;
    long lLBound;
    hr = SafeArrayGetUBound(psa, 1, &lUBound);
    hr = SafeArrayGetLBound(psa, 1, &lLBound);

    long iArraySize = lUBound - lLBound + 1;
    assert(iArraySize == coordcount);

    // Make sure there is enough space to store points
	g_mem->allocate_tri(tricount);

    CubitVector tessPoint;
    CubitVector transformed;
	for (int i=0; i<pointcount; i++) // number of points 
    {
        tessPoint.set(dTessPoints[(i*3)], dTessPoints[1+(i*3)],dTessPoints[2+(i*3)]);

        m_pSWPart->transformPartToAssembly(tessPoint, transformed);

		g_mem->point_list()[i].x = transformed.x();
		g_mem->point_list()[i].y = transformed.y();
		g_mem->point_list()[i].z = transformed.z();
	}
    // TODO - make this private and keep track of the number of points in the GMem class
	g_mem->pointListCount = pointcount;
    g_mem->fListCount = tricount * 4;

    hr = SafeArrayUnaccessData(psa);
    hr = VariantClear(&vTessTris);

	//int facetList[tricount];
	for(int j=0; j<tricount; j++)
	{
		// since all triangles have 3 points.
		g_mem->facet_list()[4*j] = 3;
        g_mem->facet_list()[4*j+1] = 3*j;
        g_mem->facet_list()[4*j+2] = 3*j+1;
        g_mem->facet_list()[4*j+3] = 3*j+2;
	}

	// Fill in the return variables and return.
	// These variables really aren't needed, because they
	// are returned as part of g_mem.
	number_points    = g_mem->pointListCount;
	facet_list_size  = g_mem->fListCount;
	number_triangles = (int)tricount;

	return CUBIT_SUCCESS;
}

void SurfaceSW::get_parents_virt(DLIList<TopologyBridge*> &parents )
{
  m_pSWPart->shells(parents);
}

void SurfaceSW::get_children_virt(DLIList<TopologyBridge*> &children )
{
    ILoop2        *pSWLoop = NULL;
    IEnumLoops2   *pLoopEnum = NULL;
    long nFetched = 0;

    HRESULT hr = m_pSWFace->EnumLoops(&pLoopEnum);
    if (SUCCEEDED(hr) && pLoopEnum)
    {
        while (S_OK == pLoopEnum->Next(1, &pSWLoop, &nFetched))
        {
            TopologyBridge *owner = m_pSWPart->cubit_owner(pSWLoop);
            if (owner)
              children.append(owner);
        }
        pLoopEnum->Release();
    }
}

CubitSense SurfaceSW::get_shell_sense( ShellSM* shell_ptr ) const
{
  assert(false);
  return CUBIT_UNKNOWN;
}

CubitStatus SurfaceSW::closest_point_uv_guess( CubitVector const& location,
                                              double& u_guess, double& v_guess,
                                              CubitVector* closest_location,
                                              CubitVector* unit_normal )
{
  assert(false);
  return CUBIT_FAILURE;
}
