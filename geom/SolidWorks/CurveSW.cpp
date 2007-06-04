//-------------------------------------------------------------------------
// Filename      : CurveSW.cpp
//
// Purpose       : 
//
// Special Notes :
//
// Creator       : Joel Kopp, Malcolm J. Panthaki
//
// Creation Date : 07/19/00
//
// Owner         : Joel Kopp
//-------------------------------------------------------------------------

// Precompiled header
#include "stdafx.h"

#include <vector>

//#include <amapp.h>
#import "C:\Program Files\SolidWorks\sldworks.tlb" raw_interfaces_only, raw_native_types, no_namespace, named_guids     //Import the SolidWorks type library

#include <swconst.h>

#include <assert.h>

#include "CubitVector.hpp"
#include "SWQueryEngine.hpp"
#include "CurveSW.hpp"
//#include "CurveFacetEvalTool.hpp"
#include "DLIList.hpp"
#include "GMem.hpp"
#include "SWPart.hpp"
#include "gtfordm.h"


extern HRESULT ExtractVARIANTArrayData(VARIANT &v, long lArraySize, double *dArray);


//-------------------------------------------------------------------------
// Purpose       : The constructor
//
// Special Notes :
//
// Creator       : Joel Kopp, Xuechen Liu
//
// Creation Date : 07/19/00
//-------------------------------------------------------------------------
CurveSW::CurveSW(SWPart *pPart)
    : sense_(CUBIT_FORWARD)
{
	m_pSWEdge = NULL;
	m_pSWPart = pPart;
}

//-------------------------------------------------------------------------
// Purpose       : The destructor. 
//
// Special Notes :
//
// Creator       : Joel Kopp, Raikanta Sahu
//
// Creation Date : 07/19/00
//-------------------------------------------------------------------------
CurveSW::~CurveSW() 
{
    if (m_extrema.size())
    {
        m_extrema.reset();
        for (int iPt=0; iPt<m_extrema.size(); iPt++)
        {
            delete m_extrema.get_and_step();
        }
    }

	set_EDGE_ptr(NULL);
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
GeometryQueryEngine* CurveSW::get_geometry_query_engine() const
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
// Creation Date : 07/19/00
//-------------------------------------------------------------------------
void CurveSW::append_simple_attribute_virt(CubitSimpleAttrib* csattrib_ptr)
{
    m_pSWPart->append_simple_attribute_virt(m_pSWEdge, csattrib_ptr);
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
// Creation Date : 07/19/00
//-------------------------------------------------------------------------
void CurveSW::remove_simple_attribute_virt(CubitSimpleAttrib* csattrib_ptr)
{
    m_pSWPart->remove_simple_attribute_virt(m_pSWEdge, csattrib_ptr);
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
// Creation Date : 07/19/00
//-------------------------------------------------------------------------
void CurveSW::remove_all_simple_attribute_virt()
{
    m_pSWPart->remove_all_simple_attribute_virt(m_pSWEdge);
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
// Creation Date : 07/19/00
//-------------------------------------------------------------------------
CubitStatus CurveSW::get_simple_attribute(DLIList<CubitSimpleAttrib*>&
                                            cubit_simple_attrib_list)
{
    return m_pSWPart->get_simple_attribute(m_pSWEdge, cubit_simple_attrib_list);
}
CubitStatus CurveSW::get_simple_attribute(const CubitString& name,
                                   DLIList<CubitSimpleAttrib*>& list)
  { return m_pSWPart->get_simple_attribute(m_pSWEdge, name, list); }

//-------------------------------------------------------------------------
// Purpose       : Get the bounding box of the object.
//
// Special Notes :
//
// Creator       : Joel Kopp, Raikanta Sahu
//
// Creation Date : 07/19/00
//-------------------------------------------------------------------------
CubitBox CurveSW::bounding_box() const 
{
    HRESULT hr = NOERROR;
    IEdge *edge = get_EDGE_ptr();

    int i;
	double dTemp;
    double oldToleranceVal;

    double xmin;
    double xmax;
    double ymin;
    double ymax;
    double zmin;
    double zmax;

    // get the length of the edge
    VARIANT vParams;
    hr = edge->GetCurveParams2(&vParams);

	double params[11];
    hr = ExtractVARIANTArrayData(vParams, 11, params);
    hr = VariantClear(&vParams);

	double tLow  = params[6];
	double tHigh = params[7];

    ICurve *curve = NULL;
    hr = edge->IGetCurve(&curve);

    double dLength;
    hr = curve->GetLength(tLow, tHigh, &dLength);

    // fit with a tolerance relative to the length of the edge
    double dCurveTol = dLength / .01;

	IModeler *modeler;
	hr = SWQueryEngine::instance()->GetSldWorks()->IGetModeler(&modeler);

    hr = modeler->SetToleranceValue(swBSCurveOutputTol, dCurveTol, &oldToleranceVal);

    VARIANT vBCurve;
    hr = curve->GetBCurveParams(FALSE, &vBCurve);

    hr = modeler->SetToleranceValue(swBSCurveOutputTol, oldToleranceVal, &dTemp);
    modeler->Release();

    double *dData;
    SAFEARRAY* psa = V_ARRAY(&vBCurve);
    hr = SafeArrayAccessData(psa, (void**)&dData);
    assert(SUCCEEDED(hr));

	union PackedInts        //unpack packed integers
	{
		double value;
		int intData[2];
	} packedIntData;

    // the first two doubles are really packed integers
    long index;
    hr = SafeArrayGetLBound(psa, 1, &index);

	packedIntData.value = dData[index];
    index++;
	int iDimension = packedIntData.intData[0];
	int iOrder = packedIntData.intData[1];

    packedIntData.value = dData[index];
    index++;
    int iNumCtrlPts = packedIntData.intData[0];

    // skip the knots to get the control point data
    index += (iNumCtrlPts + iOrder);

    // find the max and min x, y, z values
    xmin = xmax = dData[index];
    ymin = ymax = dData[index+1];
    zmin = zmax = dData[index+2];
    index += iDimension;

    for (i=1; i<iNumCtrlPts; i++)
    {
        if (dData[index] < xmin) xmin = dData[index];
        if (dData[index] > xmax) xmax = dData[index];

        if (dData[index] < ymin) ymin = dData[index+1];
        if (dData[index] > ymax) ymax = dData[index+1];

        if (dData[index] < zmin) zmin = dData[index+2];
        if (dData[index] > zmax) zmax = dData[index+2];

        index += iDimension;
    }

    hr = SafeArrayUnaccessData(psa);
    hr = VariantClear(&vBCurve);

    CubitVector corner1(xmin, ymin, zmin);
    CubitVector corner2(xmax, ymax, zmax);

    CubitVector asmCorner1;
    CubitVector asmCorner2;
    m_pSWPart->transformPartToAssembly(corner1, asmCorner1);
    m_pSWPart->transformPartToAssembly(corner2, asmCorner2);

	// Convert to a CubitBox and return it
    return CubitBox(asmCorner1, asmCorner2);
}

//-------------------------------------------------------------------------
// Purpose       : merges "this" with input GeometryEntity
//
// Special Notes :
//
// Creator       : Joel Kopp, Jihong Ma
//
// Creation Date : 07/19/00
//-------------------------------------------------------------------------
CubitStatus CurveSW::merge(GeometryEntity* /*deadGEPtr*/)
{ 
    assert(false);
    return CUBIT_FAILURE; 
}

TopologyEntity* CurveSW::unmerge( DLIList<RefVolume*> /*volumes*/ )
{
    assert(false);
	PRINT_ERROR("Unmerging of SW-based entities currently disabled.\n");
	return NULL;
}

//-------------------------------------------------------------------------
// Purpose       : Return the length of the curve.
//
// Special Notes :
//
// Creator       : Joel Kopp, Raikanta Sahu
//
// Creation Date : 07/19/00
//-------------------------------------------------------------------------
double CurveSW::measure()
{
    HRESULT hr = NOERROR;
    double startParam;
    double endParam;

	IEdge *edge = get_EDGE_ptr();
	
	if ( edge == NULL )
	{
		return 0.0;
	}

	// First get the start and end parameter values for this SW EDGE
    VARIANT vParams;
    hr = edge->GetCurveParams2(&vParams);

    double pdParams[11];
    hr = ExtractVARIANTArrayData(vParams, 11, pdParams);
    hr = VariantClear(&vParams);

    startParam = pdParams[6];
    endParam = pdParams[7];

	// Switch startParam and endParam if startParam is greater than endParam

	if ( startParam > endParam )
	{
		double temp;
		temp = endParam;
		endParam = startParam;
		startParam = temp;
	}

	else if ( startParam == endParam )
	{
		PRINT_ERROR("In CurveSW::measure\n"
					"       Couldn't resolve parameterization for "
					"the input edge\n");
		return 0.0;
	}

	double curveLength = 0.0;
	ICurve *curve = NULL;
	hr = edge->IGetCurve ( &curve );

    if (curve)
    {
    	hr = curve->GetLength ( startParam, endParam, &curveLength );
        curve->Release();
    }
    else
        assert(false);

	return curveLength;
}
//-------------------------------------------------------------------------
// Purpose       : Return the arc length along the Curve starting from
//                 the point represented by the parameter1 going to the 
//                 point represented by parameter2.
//
// Special Notes : The sign of the returned length value is always positive.
//                 Parameter1 and parameter2 are with respect to the EDGE.
//
// Creator       : Joel Kopp, Malcolm J. Panthaki
//
// Creation Date : 07/19/00
//-------------------------------------------------------------------------
double CurveSW::length_from_u( double parameter1, double parameter2 )
{
    HRESULT hr = NOERROR;

    IEdge *edge = this->get_EDGE_ptr();
     // Get the SW curve associated with the first EDGE
	ICurve *curve = NULL;

    hr = edge->IGetCurve(&curve);
	
	if ( curve == NULL )
	{
		return 0.0;
	}

	if ( parameter1 > parameter2 )
	{
		double temp = 0.0;
		temp = parameter2;
		parameter2 = parameter1;
		parameter1 = temp;
	}

	else if ( parameter1 == parameter2 )
	{
		PRINT_ERROR("In CurveSW::length_from_u\n"
					"       Identical start and end parameters\n" );
		return 0.0;
	}

	double curveLength = 0.0;
	hr = curve->GetLength ( parameter1, parameter2, &curveLength );
    curve->Release();

	return curveLength;
}

//-------------------------------------------------------------------------
// Purpose       : Returns CUBIT_TRUE and the associated period value, if 
//                 the SW curve associated with the first EDGE is periodic.
//                 Otherwise returns CUBIT_FALSE and a value of 0.0 for
//                 the period.
//
// Special Notes : Can use BCurve parameters since all curves generated
//                 in SolidWorks are BCurves.
//
// Creator       : Joel Kopp, Malcolm J. Panthaki
//
// Creation Date : 7/28/00
//-------------------------------------------------------------------------
CubitBoolean CurveSW::is_periodic(double& period)
{
    HRESULT hr = NOERROR;
    double uStart;
    double uEnd;
    VARIANT_BOOL bClosed = FALSE;
    VARIANT_BOOL bPeriodic = FALSE;
    VARIANT_BOOL bRetval = FALSE;

	// Get the SW curve associated with the first EDGE
	IEdge *edge = get_EDGE_ptr();

	ICurve *curve = NULL;
	hr = edge->IGetCurve ( &curve );
	
	if (curve == NULL)
		return CUBIT_FALSE;

    hr = curve->GetEndParams(&uStart, &uEnd, &bClosed, &bPeriodic, &bRetval);
    curve->Release();

    if (bClosed && bPeriodic)
    {
        period = fabs(uEnd - uStart);
        return CUBIT_TRUE;
    }
    else
    {
        period = 0.0;
        return CUBIT_FALSE;
    }
}

//------------------------------------------------------------------
// Purpose: Returns CUBIT_TRUE and the associated parametric values, 
//          if the SW curve associated with the first EDGE is 
//          parametric.
//          Otherwise returns CUBIT_FALSE and the values of 
//          the lower and upper parametric bounds are undetermined.
//
//
// Creator       : Joel Kopp, Malcolm J. Panthaki
//
// Creation Date : 8/29/00
//-------------------------------------------------------------------
CubitBoolean CurveSW::get_param_range(double& lower_bound,
                                      double& upper_bound)
{
    HRESULT hr = NOERROR;

	// Get the EDGE
	IEdge *edge = get_EDGE_ptr();

	// Get the start and end parameters of this EDGE
    VARIANT vParams;
	hr = edge->GetCurveParams2( &vParams );

	double curveParams[11];
    hr = ExtractVARIANTArrayData(vParams, 11, curveParams);
    hr = VariantClear(&vParams);

	double start_param = curveParams[6];
	double end_param = curveParams[7];

    // get sense from curve
    union struct_PackedData
    {
        double doubleval;
        int intvals[2];
    } packedData;

    packedData.doubleval = curveParams[10];
    int iSense = packedData.intvals[1];

    // get sense from edge and adjust parameters
	if ( CUBIT_FORWARD == get_EDGE_sense())
	{
        assert(iSense); // make sure sense from curve agrees with edge
		lower_bound = start_param;
		upper_bound = end_param;
	}

	else
	{
	  // The sense is reversed so switcherooo!!!
        assert(0 == iSense); // make sure sense from curve agrees with edge
		lower_bound = end_param;
		upper_bound = start_param;
	}
	// All SW curves are parametrically defined so...
	return CUBIT_TRUE;
}

	
	// Finds the extrema along this RefEdge.  It is the responsibility of the
	// calling code to delete the CubitVectors added to interior_points!
CubitStatus CurveSW::get_interior_extrema(DLIList<CubitVector*>& interior_points,
                                            CubitSense& return_sense)
{
    return_sense = sense_;

    // if data hasn't been cached, compute it
    if (!m_extrema.size())
    {
        // set the tolerance to 1 percent of the curve length
        double dLength = measure();
        double dTol = dLength * 0.01;

        // get a tessellation for the curve
        std::vector<double> coords;
        HRESULT hr = tessellation(dTol, coords);
        if (!SUCCEEDED(hr))
            return CUBIT_FAILURE;

        assert(0 == coords.size()%3);
        long nPoints = coords.size() / 3;

        // if the curve is a straight line there are no interior extrema
        if (3 > nPoints)
            return CUBIT_SUCCESS;

        CubitVector prev_pt(coords[0], coords[1], coords[2]);
        CubitVector curr_pt(coords[3], coords[4], coords[5]);
        CubitVector prev_vct = curr_pt - prev_pt;
        CubitVector next_vct;
        CubitVector next_pt;
        CubitVector transformed;

        for( int iPt=2; iPt<nPoints; iPt++ )
        {
            // Get a vector between the next two points
            next_pt.set(coords[3*iPt], coords[3*iPt+1], coords[3*iPt+2]);
            next_vct = next_pt - curr_pt;

            // In CurveACIS::get_interior_extrema, the extrema seem to
            // be evaluated with respect to the principle axes, so do
            // the same here.  The extrema are points at which the
            // derivitive in the specified direction (principle axis)
            // is zero.  So look for a sign change in the slope across
            // a point wrt each principle direction.
            if( (prev_vct.x() * next_vct.x() < 0.) ||  // x extrema
                (prev_vct.y() * next_vct.y() < 0.) ||  // y extrema
                (prev_vct.z() * next_vct.z() < 0.)  )  // z extrema
            {
                // transform the point from the part to the assembly coordinate system
                //
                // TODO - should the points be transformed before checking
                //        for extrema so the extrema are found in the assembly
                //        coordinate system rather than the part coordinate system?
                //
                m_pSWPart->transformPartToAssembly(curr_pt, transformed);

                m_extrema.append( new CubitVector( transformed ) );
            }

            // Advance to next point.
            prev_vct = next_vct;
            curr_pt = next_pt;
        }
    }
    // return a copy of cached data
    if (m_extrema.size())
    {
        CubitVector *pVect;
        m_extrema.reset();
        for (int iPt=0; iPt<m_extrema.size(); iPt++)
        {
            pVect = m_extrema.get_and_step();
            interior_points.append(new CubitVector(pVect->x(), pVect->y(), pVect->z()));
        }
        return CUBIT_SUCCESS;
    }

    return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Purpose       : This function computes the point on the curve closest 
//                 to the input location.  Optionally, it can also compute
//                 the tangent and curvature on the Curve at the point on
//                 on the Curve closest to the input location.
//
// Special Notes : The tangent direction is always in the positive direction of the 
//                 owning RefEdge, regardless of the positive direction of the
//                 underlying solid model entities.
//
//                 If the calling code needs the tangent and/or the curvature,
//                 it is responsible for allocating the memory for these
//                 CubitVector(s) and sending in the relevant non-NULL
//                 pointers to this routine.
//
// Creator       : Joel Kopp, Malcolm J. Panthaki
//
// Creation Date : 7/31/00
//-------------------------------------------------------------------------
CubitStatus CurveSW::closest_point( CubitVector const& location, 
                                      CubitVector& closest_location,
                                      CubitVector* tangent_ptr,
                                      CubitVector* curvature_ptr,
                                      double* param)
{
    HRESULT hr = NOERROR;

    // transform the input point from the assembly to the part coordinate system
    CubitVector tempLocation(location);
    CubitVector partLocation;
    m_pSWPart->transformAssemblyToPart(tempLocation, partLocation);

#ifdef BOYD17 
	double new_point_x = 0.0;
	double new_point_y = 0.0;
	double new_point_z = 0.0;
#endif

	double tangent_x = 0.0;
	double tangent_y = 0.0;
	double tangent_z = 0.0;

    double temp_param;

	IEdge *edge = get_EDGE_ptr();

	// Only the closest point is required
    VARIANT vPoint;
	hr = edge->GetClosestPointOn ( partLocation.x(), partLocation.y(), partLocation.z(), &vPoint );

	double closePoint[5];
    hr = ExtractVARIANTArrayData(vPoint, 5, closePoint);
    hr = VariantClear(&vPoint);

    //  - transform back to the assembly coord system
    CubitVector partClosest(closePoint[0], closePoint[1], closePoint[2]);
    m_pSWPart->transformPartToAssembly(partClosest, closest_location);

    temp_param = closePoint[3];

    if (param)
	{
		if (get_EDGE_sense() == CUBIT_REVERSED)
			temp_param = -temp_param;
		adjust_periodic_parameter(temp_param);
        *param = temp_param;
	}

	// return the tangent
	if (tangent_ptr)
	{
        VARIANT vTangent;
		hr = edge->Evaluate ( temp_param, &vTangent );

		double tangent[7];
        hr = ExtractVARIANTArrayData(vTangent, 7, tangent);
        hr = VariantClear(&vTangent);

		tangent_x = tangent[3];
		tangent_y = tangent[4];
		tangent_z = tangent[5];

		if ((sense_ == CUBIT_FORWARD && get_EDGE_sense() == CUBIT_REVERSED) ||
			(sense_ == CUBIT_REVERSED && get_EDGE_sense() == CUBIT_FORWARD))
		{
			tangent_x = -tangent_x;
			tangent_y = -tangent_y;
			tangent_z = -tangent_z;
		}

        CubitVector partTangent(tangent_x, tangent_y, tangent_z);
        CubitVector partOrigin(0.0, 0.0, 0.0);
        CubitVector asmTangent;
        CubitVector asmOrigin;

        m_pSWPart->transformPartToAssembly(partTangent, asmTangent);
        m_pSWPart->transformPartToAssembly(partOrigin, asmOrigin);

        *tangent_ptr = asmTangent - asmOrigin;
	}

	// NOTE: SolidWorks has no way to compute the curvature of a curve.
    // Use the OLE for D&M interfaces, supported by SolidWorks
	if (curvature_ptr)
	{
        LPDMCURVE pDMCurve = NULL;
        hr = edge->QueryInterface(IID_IDMCurve, (LPVOID*)&pDMCurve);
        if (SUCCEEDED(hr) && pDMCurve)
        {
            long nParams = 1;
            double dFirstDeriv[3];
            double dSecondDeriv[3];
            hr = pDMCurve->GetDerivatives(nParams, &temp_param, dFirstDeriv, dSecondDeriv, NULL);
            pDMCurve->Release();

            CubitVector first(dFirstDeriv);
            CubitVector second(dSecondDeriv);

            // curvature = |f'(u) X f''(u)| / (|f'(u)|^3)
            CubitVector tempVector = first * second; // cross product
            double dFirstLength = first.length();
            double dCurvature = tempVector.length() / (dFirstLength * dFirstLength * dFirstLength);

            second.normalize(); // unit normal vector points to center of curvature
            second*=dCurvature;

            CubitVector partOrigin(0.0, 0.0, 0.0);
            CubitVector asmOrigin;
            CubitVector asmVector;
            m_pSWPart->transformPartToAssembly(second, asmVector);
            m_pSWPart->transformPartToAssembly(partOrigin, asmOrigin);

            *curvature_ptr = asmVector - asmOrigin;
        }
        else
        {
            assert(false);
        }

//		// get the sense wrt the SW curve, consists of 2 sense values
//		if ((sense_ == CUBIT_FORWARD && get_EDGE_sense() == CUBIT_REVERSED) ||
//			(sense_ == CUBIT_REVERSED && get_EDGE_sense() == CUBIT_FORWARD))
//			curvature = -(curvature);
//
//		curvature_ptr->set( curvature.x(), curvature.y(), curvature.z() );
	}


	return CUBIT_SUCCESS;
}

//------------------------------------------------------------------
// Purpose       : return the tangent to this Curve
//                 
//
// Special Notes :
//
// Creator       : Joel Kopp, Tim Tautges
//
// Creation Date : 7/31/00
//------------------------------------------------------------------
void CurveSW::get_tangent( CubitVector const& location,
                             CubitVector& tangent)
{
	CubitVector closest_location;
	this->closest_point(location, closest_location, &tangent);
}

//------------------------------------------------------------------
// Purpose       : Return the curvature of this Curve - currently disabled
//
// Special Notes :
//
// Creator       : Joel Kopp, Tim Tautges
//
// Creation Date : 7/31/00
//------------------------------------------------------------------
void CurveSW::get_curvature( CubitVector const& location, 
                               CubitVector& curvature)
{
    assert(false);
	PRINT_ERROR("get_curvature currently disabled.\n");   

	/*CubitVector closest_location, tangent;
	this->closest_point(location, closest_location, &tangent, &curvature);*/
}

//------------------------------------------------------------------
// Purpose       : Return a pointer to the SW EDGE
//                 in this Curve.
//
// Special Notes :
//
// Creator       : Joel Kopp, Malcolm J. Panthaki
//
// Creation Date : 7/31/00
//------------------------------------------------------------------
IEdge *CurveSW::get_EDGE_ptr() const
{
	return m_pSWEdge;
}

void CurveSW::set_EDGE_ptr(IEdge *edge)
{
	if (edge == m_pSWEdge)
		return;

	if (m_pSWEdge)
	{
		m_pSWPart->remove_cubit_owner(m_pSWEdge);
        m_pSWEdge->Release();
	}

	m_pSWEdge = edge;

	if (m_pSWEdge)
	{
		m_pSWPart->set_cubit_owner(m_pSWEdge, this);
        m_pSWEdge->AddRef();
	}
}

//------------------------------------------------------------------
// Purpose: This function returns the coordinate of a point in the local
//          parametric (u) space that corresponds to the input position 
//          in global (world) space.  The input point is first moved to 
//          the closest point on the Curve and the parameter value of 
//          that point is determined. 
//
// Creator       : Joel Kopp, Malcolm J. Panthaki
//
// Creation Date : 8/1/00
//-------------------------------------------------------------------
CubitStatus CurveSW::position_from_u (double u_value,
                                        CubitVector& output_position)
{
    HRESULT hr = NOERROR;
#ifdef BOYD17 
	ICurve *curve = NULL;
#endif
	CubitStatus status = CUBIT_SUCCESS;

	// Get the parameter range of this Curve
	double lower_bound = 0.0;
	double upper_bound = 0.0;
	CubitBoolean paramBool = get_param_range(lower_bound, upper_bound);

	// Make sure the requested u_value is either within the range or, if
	// the Curve is periodic, then reduce the input value down to the
	// fundamental range

	if (u_value > upper_bound || u_value < lower_bound)
	{
		adjust_periodic_parameter(u_value);
		if (u_value > upper_bound || u_value < lower_bound)
		{
			PRINT_ERROR("In CurveSW::position_from_u\n"
					  "       Input parameter value %lf is not within the"
					  " parameter range of this Curve (%lf to %lf)\n",
					  u_value, lower_bound, upper_bound);
			//      assert(param_interval >> u_value == TRUE);
			// find the closest endpoint and return failure
			if (lower_bound < upper_bound)
			{
				if (u_value < lower_bound)
					u_value = lower_bound;
				else
					u_value = upper_bound;
			}
			else
			{
				if (u_value < upper_bound)
					u_value = upper_bound;
				else
					u_value = lower_bound;
			}
			status = CUBIT_FAILURE;
		}
	}

    // Get the first SW curve
	IEdge *edge = this->get_EDGE_ptr();

	// Now that we have a "valid" parameter value, get its global location
	// on the Curve
    VARIANT vPosition;
    hr = edge->Evaluate ( u_value, &vPosition );

	double position[7];
    hr = ExtractVARIANTArrayData(vPosition, 7, position);
    hr = VariantClear(&vPosition);

    // transform output to assembly coordinate system
    CubitVector partPosition(position[0], position[1], position[2]);
    m_pSWPart->transformPartToAssembly(partPosition, output_position);
	
	return status;
}

//-------------------------------------------------------------------------
// Purpose       : This function returns the coordinate of a point in the local
//                 parametric (u) space that corresponds to the input position 
//                 in global (world) space.  The input point is first moved to 
//                 the closest point on the Curve and the parameter value of 
//                 that point is determined. 
//
// Special Notes : 
//
// Creator       : Joel Kopp, Malcolm J. Panthaki
//
// Creation Date : 8/2/00
//-------------------------------------------------------------------------
double CurveSW::u_from_position (const CubitVector& input_position)
{
	// Get the closest point on the Curve to the input position
	CubitVector closest_point;
	double u_val;
	this->closest_point(input_position, closest_point,
					  NULL, NULL, &u_val);
	// closest_point already makes adjustments for sense and periodicity

	return u_val;
}

//------------------------------------------------------------------
// Purpose: This function returns the parameter value of the point 
//          that is "arc_length" away from the root point, in the
//          positive sense direction of the owning RefEdge.
//
// Special Notes : 
//   If arc_length is negative, the new point (whose parameter value
//   is being computed) is in the negative sense direction (along
//   the RefEdge) from the root point (whose parameter value is
//   root_param).
//
//   If the curve is not periodic and the new point, "arc_length"
//   away from the root point in the appropriate direction, goes
//   beyond the end point of the first EDGE, that end point is used
//   to generate the returned parameter value.
//
// If the curve is periodic and the new point, "arc_length" away
// from the root point in the appropriate direction, goes beyond
// the end point of the first EDGE, wrap around is done.  After
// wrap around, the point is treated as with other curves
//
// NOTE:
// The important assumption that is made in this routine is that
// the end points of the RefEdge that owns this CurveSW are the
// same as the end points of the first SW EDGE in the list of EDGEs
// associated with this CurveSW.
//
// Assume that the parameter root_param is with respect to the
// RefEdge as well as arc_length.  Before calling the SW "curve",
// we need to get them with respect to the curve.   
//
// Creator       : Joel Kopp, Malcolm J. Panthaki
//
// Creation Date : 8/2/00
//------------------------------------------------------------------
double CurveSW::u_from_arc_length ( double root_param,
                                      double arc_length )
{
    HRESULT hr = NOERROR;

	// Sanity check
	if (arc_length == 0.0)
	{
		return root_param;
	}
	// get the EDGE and curve
	IEdge *edge = this->get_EDGE_ptr();
	ICurve *curve = NULL;
    
    hr = edge->IGetCurve(&curve);

	// Check for NULL geometry
	if ( curve == NULL )
		return start_param();

    curve->Release();
    curve = NULL;

	// Get the lower and upper parameter values for the first EDGE.
	// Note that the lower and upper bounds are with respect to the start and
	// end locations of the EDGE, not the RefEdge.
	double EDGE_start_param = start_param();
	double EDGE_end_param = end_param();

	// If the sense of the RefEdge with respect to this Curve is REVERSED,
	// then change the sign of the input arc_length and then work on the
	// first EDGE associated with this Curve instead of constantly worrying 
	// about the relative sense value.
	if ( sense_ == CUBIT_REVERSED )
	{
		arc_length = -(arc_length);
	}

	// If the Curve is periodic, then let SW do the work.  Keep in mind 
	// that SW EDGEs also have a sense with respect to their underlying 
	// SW curves.
	double period = 0.0;
	if ( this->is_periodic(period) == CUBIT_TRUE ) 
	{
		double curveLength = measure();
		double param = root_param + (arc_length / curveLength * (EDGE_end_param - EDGE_start_param));
		return param;
	}

	// Now, let's deal with non-periodic Curves...

	// Find out how much "headroom" there is between the root_point and the 
	// "relevant" end of the EDGE.
	//
	// If arc_length is negative, then we have to find the 
	// distance between root point and the start of the EDGE. If arc_length
	// is positive, then we have to find the distance between root_point and
	// the end of the EDGE.
	double headroom = CUBIT_DBL_MAX;
	if(arc_length < 0.0)
		headroom = this->length_from_u( root_param, EDGE_start_param );
	else
		headroom = this->length_from_u( root_param, EDGE_end_param );

	// If the arc_length specifies a point past the "relevant" end of the 
	// EDGE, choose that end -- i.e., don't go beyond the end of the
	// EDGE.
	if (fabs(arc_length) >  headroom) 
	{
	  // This EDGE is not "periodic", so stop at the relevant endpoint
	  // and return its parameter value.
	  // NOTE: We have already taken into account the sense of the EDGE
	  //       with respect to its underlying curve...whew!!! :-)
		if (arc_length > 0.0)
		{
			return EDGE_end_param;
		}
		else
		{
			return EDGE_start_param;
		}
	}

	// The EDGE is non-periodic but arc_length doesn't put us past the end of 
	// the EDGE.
	// NOTE: The following function does not work correctly if the
	//       curve is bounded and the required arc_length puts you over
	//       the end of the curve.  The result, in that case, is *NOT 
	//       DEFINED* and there is no error reported!!  We can use it
	//       as we've already checked to make sure we don't go off the
	//       edge of the EDGE :-)
	else
	{
	  // Assert that root_param lies on the bounded EDGE
		assert( EDGE_end_param > EDGE_start_param ?
				root_param >= (EDGE_start_param - GEOMETRY_RESABS*10) &&
				root_param <= (EDGE_end_param + GEOMETRY_RESABS*10) :
				root_param <= (EDGE_start_param + GEOMETRY_RESABS*10) &&
				root_param >= (EDGE_end_param - GEOMETRY_RESABS*10) );

		double curveLength = measure();
		double param = root_param + (arc_length / curveLength) * (EDGE_end_param - EDGE_start_param);

        return param;
	}
}

//-------------------------------------------------------------------------
// Purpose       : This function tests the passed in position to see if
//                 is on the underlying curve. 
//
// Special Notes :
//
// Creator       : Joel Kopp, David R. White
//
// Creation Date : 8/2/00
//-------------------------------------------------------------------------
CubitBoolean CurveSW::is_position_on( const CubitVector &test_position )
{
    assert(false);
    return CUBIT_FALSE;
}

//-------------------------------------------------------------------------
// Purpose       : This function returns the type of underlying curve. 
//
// Special Notes : It checks to see if *any* of the SW curves associated
//                 with the EDGEs in the list of EDGEs of this Curve is of
//                 a particular type and returns the appropriate value
//                 of the enum, CurveType.
//
// Creator       : Joel Kopp, Malcolm J. Panthaki
//
// Creation Date : 8/2/00
//-------------------------------------------------------------------------
GeometryType CurveSW::geometry_type()
{
  HRESULT hr = NOERROR;

  IEdge *edge = this->get_EDGE_ptr();
  ICurve *curve = NULL;

  hr = edge->IGetCurve(&curve);
    
  // Get its type 
  if ( curve == NULL )
    return POINT_CURVE_TYPE;

  GeometryType local_type;
  VARIANT_BOOL is_circle = FALSE;
  VARIANT_BOOL is_line = FALSE;
  VARIANT_BOOL is_ellipse = FALSE;

  hr = curve->IsCircle ( &is_circle );
  hr = curve->IsLine ( &is_line );
  hr = curve->IsEllipse ( &is_ellipse );
  curve->Release();
  curve = NULL;

  if (is_circle)
    local_type = ARC_CURVE_TYPE;
  else if (is_line)
    local_type = STRAIGHT_CURVE_TYPE;
  else if (is_ellipse)
    local_type = ELLIPSE_CURVE_TYPE;
  else
    local_type = UNDEFINED_CURVE_TYPE;

  return local_type;
}

CubitStatus CurveSW::get_point_direction( CubitVector& point, 
                                            CubitVector& direction )
{
  HRESULT hr = NOERROR;

  if( geometry_type() != STRAIGHT_CURVE_TYPE )
  {
    assert(false);
    return CUBIT_FAILURE;
  }

  IEdge *edge = get_EDGE_ptr();
  ICurve *curve = NULL;
  hr = edge->IGetCurve ( &curve );

  VARIANT vParams;
  double lineParams[6];
  hr = curve->get_LineParams( &vParams );
  curve->Release();
  curve = NULL;

  hr = ExtractVARIANTArrayData(vParams, 6, lineParams);
  hr = VariantClear(&vParams);

  // transform from part to assembly coord system
  CubitVector partPoint(lineParams[0], lineParams[1], lineParams[2]);
  CubitVector partOrigin(0.0, 0.0, 0.0);
  CubitVector partDirection(lineParams[3], lineParams[4], lineParams[5]);

  CubitVector asmOrigin;
  CubitVector asmDirection;
  m_pSWPart->transformPartToAssembly(partPoint, point);
  m_pSWPart->transformPartToAssembly(partOrigin, asmOrigin);
  m_pSWPart->transformPartToAssembly(partDirection, asmDirection);

  direction = asmDirection - asmOrigin;

  return CUBIT_SUCCESS;
}

CubitStatus CurveSW::get_center_radius( CubitVector& center, 
                                         double& radius )
{
    assert(false);
    return CUBIT_FAILURE;
/*
	HRESULT res;
	
	if( geometry_type() != ELLIPSE_CURVE_TYPE ||
		geometry_type() != CIRCLE_CURVE_TYPE )
		return CUBIT_FAILURE;

	LPEDGE edge = get_EDGE_ptr();
	LPCURVE curve = NULL;
	edge->IGetCurve ( &curve );

	if(geometry_type() == CIRCLE_CURVE_TYPE)
	{
		double* circleParams;
		circleParams = (double*)malloc(7 * sizeof(double));
		res = curve->get_ICircleParams( circleParams );

		if(res)
		{
			center.set( circleParams[0], circleParams[1], circleParams[2] );
			radius = circleParams[6];

			return CUBIT_SUCCESS;
		}
		else
			return CUBIT_FAILURE;
	}

	else if(geometry_type() == ELLIPSE_CURVE_TYPE)
	{
		double* ellipseParams;
		ellipseParams = (double*)malloc(11 * sizeof(double));
		res = curve->IGetEllipseParams ( ellipseParams );

		if(res)
		{
			center.set( ellipseParams[0], ellipseParams[1], ellipseParams[2] );
			radius = ellipseParams[3]; // major axis radius

			return CUBIT_SUCCESS;
		}
		else
			return CUBIT_FAILURE;
	}

	else
		return CUBIT_FAILURE;
*/
}

//-------------------------------------------------------------------------
// Purpose       : This function returns the start parameter.
//
// Special Notes : The start param is with respect to the ref_edge.
//
// Creator       : Joel Kopp, David White
//
// Creation Date : 8/2/00
//-------------------------------------------------------------------------
double CurveSW::start_param()
{
	double start = 0.0, end = 0.0;

	get_param_range( start, end );
	return start;
}
//-------------------------------------------------------------------------
// Purpose       : This function returns the end parameter.
//
// Special Notes : The end param is with respect to the ref_edge.
//
// Creator       : Joel Kopp, David White
//
// Creation Date : 8/2/00
//-------------------------------------------------------------------------
double CurveSW::end_param()
{
	double start = 0.0, end = 0.0;

	get_param_range( start, end );
	return end;
}

void CurveSW::bodysms(DLIList<BodySM*> &bodies) 
{
	if (m_pSWEdge)
		SWQueryEngine::instance()->bodysms(m_pSWEdge, bodies);
}

void CurveSW::lumps(DLIList<Lump*> &lumps)
{
	if (m_pSWEdge)
		SWQueryEngine::instance()->lumps(m_pSWPart, m_pSWEdge, lumps);
}

void CurveSW::shellsms(DLIList<ShellSM*> &shellsms)
{
	if (m_pSWEdge)
		SWQueryEngine::instance()->shellsms(m_pSWPart, m_pSWEdge, shellsms);
}

void CurveSW::surfaces(DLIList<Surface*> &surfaces)
{
	if (m_pSWEdge)
		SWQueryEngine::instance()->surfaces(m_pSWPart, m_pSWEdge, surfaces);
}

void CurveSW::loopsms(DLIList<LoopSM*> &loopsms)
{
	if (m_pSWEdge)
		SWQueryEngine::instance()->loopsms(m_pSWPart, m_pSWEdge, loopsms);
}

void CurveSW::curves(DLIList<Curve*> &curves)
{
	if (m_pSWEdge)
		SWQueryEngine::instance()->curves(m_pSWPart, m_pSWEdge, curves);
}

void CurveSW::coedgesms(DLIList<CoEdgeSM*> &coedgesms)
{
	if (m_pSWEdge)
		SWQueryEngine::instance()->coedgesms(m_pSWPart, m_pSWEdge, coedgesms);
}

void CurveSW::points(DLIList<Point*> &points)
{
	if (m_pSWEdge)
        SWQueryEngine::instance()->points(m_pSWPart, m_pSWEdge, points);
}

//-------------------------------------------------------------------------
// Purpose       : Check for G1 discontinuity
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 05/16/00
//-------------------------------------------------------------------------
CubitBoolean CurveSW::G1_discontinuous( 
      double /*param*/, CubitVector* /*mtan*/, CubitVector* /*ptan*/ )
{
    // TODO - this function gets called, but does this need to be implemented?
    // It seems like there should be a better way to do this, like breaking up
    // the curve into G1 continuous pieces for processing so the check can be
    // made once.
  return CUBIT_FALSE;
}

//-------------------------------------------------------------------------
// Purpose       : define a tool for evaluation of curve based on
//                 underlying facet representation
//
// Special Notes :
//
// Creator       : Joel Kopp, Steve J. Owen
//
// Creation Date : 8/2/00
//-------------------------------------------------------------------------
//CubitStatus CurveSW::setup_use_facets(
//	FacetEvalTool *surf_eval_tool,
//	CubitVector &start, CubitVector &end,
//	CubitSense orientation )
//{
//    assert(false);
//	if ( facetEvalTool )
//	{
//		  //get rid of the old facets.
//		delete facetEvalTool;
//	}
//	//create a facet eval tool.

//	facetEvalTool = new CurveFacetEvalTool( surf_eval_tool,
//										  start, end, orientation );

//	useFacets = CUBIT_TRUE;
//	return CUBIT_SUCCESS;
//}

// ********** END PUBLIC FUNCTIONS         **********

// ********** BEGIN PROTECTED FUNCTIONS    **********
// ********** END PROTECTED FUNCTIONS      **********

// ********** BEGIN PRIVATE FUNCTIONS      **********

//-------------------------------------------------------------------------
// Purpose       : Return the sense of the first EDGE (in the list of EDGEs
//                 in this Curve) wrt its underlying SW curve.
//
// Special Notes : In SW, the curves always share the sense of the edge.
//
// Creator       : Joel Kopp, Malcolm J. Panthaki
//
// Creation Date : 8/29/00
//-------------------------------------------------------------------------
CubitSense CurveSW::get_EDGE_sense()
{
  HRESULT hr = NOERROR;

  IEdge *edge = this->get_EDGE_ptr();

  VARIANT_BOOL isReversed;
  hr = edge->IsParamReversed(&isReversed);

  if (isReversed)
    return CUBIT_REVERSED;
  else
    return CUBIT_FORWARD;
}

//----------------------------------------------------------------
// Adjusts the input parameter so that it falls within the
// parameter range of this Curve, if possible.  Necessary for
// periodic curves.
//----------------------------------------------------------------
void CurveSW::adjust_periodic_parameter(double& param)
{
    // Adjustment only legal if this is a periodic curve.
	CubitBoolean result;

	double period;
	double lower_bound = 0.0;
	double upper_bound = 0.0;
	if ( this->is_periodic(period) && (fabs(period) > CUBIT_RESABS))
	{
		result = get_param_range(lower_bound, upper_bound);
		assert( (upper_bound - lower_bound) > CUBIT_RESABS * 100 );

	  // Make sure period is positive
	if (period < 0)
		period = -period;

	  // Move the parameter above the low param
	while (param < (lower_bound - CUBIT_RESABS))
		param += period;
	  // Move the parameter below the high param
	while (param > (upper_bound + CUBIT_RESABS))
		param -= period;
	}
}

CubitPointContainment 
CurveSW::point_containment( const CubitVector &point )
{
    assert(false);
    return CUBIT_PNT_UNKNOWN;
}

CubitStatus
CurveSW::facet_edge(double tolerance, int& num_points, GMem* g_mem)
{
    // initialize output
    num_points = 0;

    // get a tessellation for the curve
    std::vector<double> coords;
    HRESULT hr = tessellation(tolerance, coords);
    if (!SUCCEEDED(hr))
        return CUBIT_FAILURE;

    assert(0 == coords.size()%3);
    int nPts = coords.size() / 3;
    int nPoly = nPts - 1;

    // Make sure there is enough space to store points
    g_mem->allocate_polylines(nPoly);

    CubitVector tessPoint;
    CubitVector transformed;
    for (int i=0; i<nPts; i++) // number of points 
    {
        tessPoint.set(coords[(i*3)], coords[1+(i*3)],coords[2+(i*3)]);

        m_pSWPart->transformPartToAssembly(tessPoint, transformed);

        g_mem->point_list()[i].x = (float)(transformed.x());
        g_mem->point_list()[i].y = (float)(transformed.y());
        g_mem->point_list()[i].z = (float)(transformed.z());
	}
    // TODO - make this private and keep track of the number of points in the GMem class
    g_mem->pointListCount = nPts;
    num_points = nPts;

    return CUBIT_SUCCESS;
}

HRESULT CurveSW::tessellation(const double &dChordTol, std::vector<double> &coordinates)
{
  HRESULT hr = NOERROR;

  // get the SolidWorks Edge
  IEdge *pSWEdge = get_EDGE_ptr();
  assert(pSWEdge);

    // bet the curve for the edge
  ICurve *curve = NULL;
  hr = pSWEdge->IGetCurve ( &curve );

    // Get the length to determine length tolerance.
    //   get the start and end parameters for the edge and compute the
    //   length between them
  double params[11];
  VARIANT vParams;
  hr = pSWEdge->GetCurveParams2(&vParams);
  hr = ExtractVARIANTArrayData(vParams, 11, params);
  hr = VariantClear(&vParams);

  double tLow  = params[6]; // start param
  double tHigh = params[7]; // end param

  double length = 0.0;
  hr = curve->GetLength ( tLow, tHigh, &length );

  // The length tolerance is the shortest length the tesselation chord
  // may have.  In order to keep within MAX_NUM_CURVE_POINTS, must make
  // the shortest chord length the total curve length divided by the 
  // max number of tesselation points.
  const int MAX_NUM_CURVE_POINTS = 750;
  double lengthTol = length / ((double)MAX_NUM_CURVE_POINTS);

  VARIANT vStartPoint;
  VARIANT vEndPoint;
  VariantInit(&vStartPoint);
  VariantInit(&vEndPoint);

  SAFEARRAY *psaStartPt;
  SAFEARRAY *psaEndPt;

  // To continue, must place the start and end points in a SafeArray
  // to pass to SW functions.
  IVertex *startVertex = NULL;
  hr = pSWEdge->IGetStartVertex ( &startVertex );
  IVertex *endVertex = NULL;
  hr = pSWEdge->IGetEndVertex ( &endVertex );

  if (startVertex)
  {
    assert(endVertex);
    hr = startVertex->GetPoint ( &vStartPoint );
    hr = endVertex->GetPoint ( &vEndPoint );
  }
  else
  {
    // fill safearrays with start and end points from curve param data
    psaStartPt = SafeArrayCreateVector(VT_R8, 0, 3);
    psaEndPt = SafeArrayCreateVector(VT_R8, 0, 3);

    long isa;
    for (isa=0; isa<3; isa++)
      SafeArrayPutElement(psaStartPt, &isa, (void *)&params[isa]);

    for (isa=3; isa<6; isa++)
      SafeArrayPutElement(psaEndPt, &isa, (void *)&params[isa]);

    vStartPoint.vt = VT_ARRAY;
    vStartPoint.parray = psaStartPt;      
    vEndPoint.vt = VT_ARRAY;
    vEndPoint.parray = psaEndPt;
  }

  // Now get the tesselation points
  long lTessSize;
  VARIANT vTessPts;
  hr = curve->GetTessPts ( (double)dChordTol, lengthTol, vStartPoint, vEndPoint, &vTessPts );
  curve->Release();
  // NOTE - hr is not returning an HRESULT - Byron

  hr = VariantClear(&vStartPoint);
  hr = VariantClear(&vEndPoint);


  SAFEARRAY* psa = V_ARRAY(&vTessPts);
  double *dTessPoints;
  hr = SafeArrayAccessData(psa, (void**)&dTessPoints);

  long lUBound;
  long lLBound;
  hr = SafeArrayGetUBound(psa, 1, &lUBound);
  hr = SafeArrayGetLBound(psa, 1, &lLBound);

  lTessSize = lUBound - lLBound + 1;
//    int nPoly = (tesssize / 3) - 1;
//    int nPts = nPoly + 1;
  for (int i=0; i<lTessSize; i++)
    coordinates.push_back(dTessPoints[i]);

  hr = SafeArrayUnaccessData(psa);
  hr = VariantClear(&vTessPts);

  return S_OK;
}

void CurveSW::get_parents_virt(DLIList<TopologyBridge*> &parents )
{
    ICoEdge       *pSWCoEdge = NULL;
    IEnumCoEdges  *pCoEdgeEnum = NULL;
    long nFetched = 0;

    HRESULT hr = m_pSWEdge->EnumCoEdges(&pCoEdgeEnum);
    if (SUCCEEDED(hr) && pCoEdgeEnum)
    {
        while (S_OK == pCoEdgeEnum->Next(1, &pSWCoEdge, &nFetched))
        {
            TopologyBridge *owner = m_pSWPart->cubit_owner(pSWCoEdge);
            if (owner)
              parents.append(owner);
        }
        pCoEdgeEnum->Release();
    }
}

void CurveSW::get_children_virt(DLIList<TopologyBridge*> &children )
{
  IVertex *swvertex;
  int num_vertex = 0;

  // start vertex, if any
  HRESULT hr = m_pSWEdge->IGetStartVertex(&swvertex);
  if (SUCCEEDED(hr) && swvertex)
  {
    num_vertex++;

    TopologyBridge *owner = m_pSWPart->cubit_owner(swvertex);
    if (owner)
      children.append(owner);
  }

  hr = m_pSWEdge->IGetEndVertex(&swvertex);
  if (SUCCEEDED(hr) && swvertex)
  {
    num_vertex++;

    TopologyBridge *owner = m_pSWPart->cubit_owner(swvertex);
    if (owner)
      children.append(owner);
  }

  if (0 == num_vertex)
  {
    TopologyBridge *pseudo_vert = m_pSWPart->findPseudoVertex(m_pSWEdge);
    if (pseudo_vert)
      children.append(pseudo_vert);
  }
}
