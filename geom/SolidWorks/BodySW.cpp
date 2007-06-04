//-------------------------------------------------------------------------
// Filename      : BodySW.cpp
//
// Purpose       : 
//
// Special Notes : COM implementation.
//
// Creator       : Joel Kopp, Xuechen Liu
//
// Creation Date : 8/3/00
//
// Owner         : Joel Kopp, Malcolm J. Panthaki
//-------------------------------------------------------------------------

// Precompiled header
#include "stdafx.h"

//#include <amapp.h>
#import "C:\Program Files\SolidWorks\sldworks.tlb" raw_interfaces_only, raw_native_types, no_namespace, named_guids     //Import the SolidWorks type library

#include <assert.h>


#include "BodySW.hpp"

#include "SWQueryEngine.hpp"
#include "SWPart.hpp"

extern HRESULT ExtractVARIANTArrayData(VARIANT &v, long lArraySize, double *dArray);

//-------------------------------------------------------------------------
// Purpose       : 
//
// Special Notes :
//
// Creator       : Joel Kopp, Xuechen Liu
//
// Creation Date : 8/3/00
//-------------------------------------------------------------------------
BodySW::BodySW(SWPart* pPart)
{
	m_pSWBody = NULL;
    m_pSWPart = pPart;
}

//-------------------------------------------------------------------------
// Purpose       : The destructor.
//
// Special Notes :
//
// Creator       : Joel Kopp, Raikanta Sahu
//
// Creation Date : 8/3/00
//-------------------------------------------------------------------------
BodySW::~BodySW()
{
    set_BODY_ptr(NULL);
}

//-------------------------------------------------------------------------
// Purpose       : 
//
// Special Notes : 
//
// Creator       : 
//
// Creation Date : 
//-------------------------------------------------------------------------

IBody2 *BodySW::get_BODY_ptr() const
{
	assert (m_pSWBody);
	return m_pSWBody;
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
GeometryQueryEngine* BodySW::get_geometry_query_engine() const
{
    return (GeometryQueryEngine*)SWQueryEngine::instance();   
}                 

//-------------------------------------------------------------------------
// Purpose       : The purpose of this function is to append a
//                 attribute to the OSME. The name is attached to the 
//                 underlying solid model entity this one points to.
//
//
// Special Notes : 
//
// Creator       : Joel Kopp, Malcolm J. Panthaki
//
// Creation Date : 8/3/00
//-------------------------------------------------------------------------
void BodySW::append_simple_attribute_virt(CubitSimpleAttrib* csattrib_ptr)
{
    m_pSWPart->append_simple_attribute_virt(this, csattrib_ptr);
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
// Creation Date : 8/3/00
//-------------------------------------------------------------------------
void BodySW::remove_simple_attribute_virt(CubitSimpleAttrib* csattrib_ptr)
{
    m_pSWPart->remove_simple_attribute_virt(this, csattrib_ptr);
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
// Creation Date : 8/3/00
//-------------------------------------------------------------------------
void BodySW::remove_all_simple_attribute_virt()
{
    m_pSWPart->remove_all_simple_attribute_virt(this);
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
// Creation Date : 8/3/00
//-------------------------------------------------------------------------
CubitStatus BodySW::get_simple_attribute(DLIList<CubitSimpleAttrib*>&
                                           cubit_simple_attrib_list)
{
	return m_pSWPart->get_simple_attribute(this, cubit_simple_attrib_list);
}
CubitStatus BodySW::get_simple_attribute(const CubitString& name,
                                   DLIList<CubitSimpleAttrib*>& list)
  { return m_pSWPart->get_simple_attribute(this, name, list); }

//-------------------------------------------------------------------------
// Purpose       : Makes a copy of this BodySW object
//
// Special Notes : Also copies the SW BODY that is contained in this
//                 BodySW object.
//
// Creator       : Joel Kopp, Malcolm Panthaki
//
// Creation Date : 8/3/00
//-------------------------------------------------------------------------
//OtherSolidModelEntity* BodySW::copy()
//{
//  LPBODY new_BODY_ptr = 
//    SWGeometryEngine::instance()->copy_BODY(get_BODY_ptr());
  
//  if (new_BODY_ptr != NULL)
//  {
//      // Make a new BodySW object using the new BODY
//    BodySW* new_BodySW_ptr = new BodySW(new_BODY_ptr);
    
      // Return the new beast if all's well :-)
//    if (new_BodySW_ptr != NULL)
//    {
//      return new_BodySW_ptr;
//    }
    
//    else
//    {
//      PRINT_ERROR("In BodySW::copy\n"
//                  "       Problems creating a new BodySW object\n");
//      return (OtherSolidModelEntity *)NULL;
//    }
//  }
//  else
//  {
//    PRINT_ERROR("In BodySW::copy\n"
//                "       Problems copying the SW BODY\n");
//    return (OtherSolidModelEntity *)NULL;
//  }
//}

//-------------------------------------------------------------------------
// Purpose       : Move the SW BODY by dx, dy and dz
//
// Special Notes : 
//
// Creator       : Joel Kopp, Malcolm Panthaki
//
// Creation Date : 8/7/00
//-------------------------------------------------------------------------
CubitStatus BodySW::move(double x_offset, double y_offset, double z_offset)
{
    assert(false);
    return CUBIT_FAILURE;
}

//-------------------------------------------------------------------------
// Purpose       : Rotate the SW BODY by theta_x, theta_y, and theta_z
//
// Special Notes : 
//
// Creator       : Joel Kopp, Malcolm Panthaki
//
// Creation Date : 8/7/00
//-------------------------------------------------------------------------
CubitStatus BodySW::rotate( double theta_x, double theta_y, double theta_z,
						    double theta)
{
    assert(false);
	return CUBIT_FAILURE;
}

//-------------------------------------------------------------------------
// Purpose       : Scale the SW BODY by the factor, scaling_factor
//
// Special Notes : 
//
// Creator       : Joel Kopp, Malcolm Panthaki
//
// Creation Date : 8/7/00
//-------------------------------------------------------------------------
CubitStatus BodySW::scale(double scaling_factor)
{
    assert(false);
    return CUBIT_FAILURE;
}

CubitStatus BodySW::scale(double scaling_factor_x,
                            double scaling_factor_y,
                            double scaling_factor_z)
{
    assert(false);
    return CUBIT_FAILURE;
}


//-------------------------------------------------------------------------
// Purpose       : Restore the SW BODY by replacing the current 
//                 transformation matrix with a unit matrix
//
// Special Notes : 
//
// Creator       : Joel Kopp, Malcolm Panthaki
//
// Creation Date : 8/8/00
//-------------------------------------------------------------------------
CubitStatus BodySW::restore() // must pass component 
												 // pointer
{
    assert(false);
    return CUBIT_FAILURE;
}

//-------------------------------------------------------------------------
// Purpose       : Reverse the SW BODY (so far, the only objects reversed
//				   are the curves)
//
// Special Notes : 
//
// Creator       : Joel Kopp, Tim Tautges
//
// Creation Date : 8/8/00
//-------------------------------------------------------------------------
CubitStatus BodySW::reverse()
{
    assert(false);
    return CUBIT_FAILURE;
}

//-------------------------------------------------------------------------
// Purpose       : Reverse the SW BODY (so far, the only objects reversed
//				   are the curves)
//
// Special Notes : 
//
// Creator       : Joel Kopp, Tim Tautges
//
// Creation Date : 8/8/00
//-------------------------------------------------------------------------
CubitStatus BodySW::reverse(IBody2 *body)
{
    assert(false);
    return CUBIT_FAILURE;
}

void BodySW::set_BODY_ptr(IBody2 *body)
{
	if (body == m_pSWBody)
		return;

	if (m_pSWBody)
	{
		m_pSWPart->remove_cubit_owner(m_pSWBody);
		m_pSWBody->Release();
	}

	m_pSWBody = body;

	if (m_pSWBody)
	{
		m_pSWPart->set_cubit_owner(m_pSWBody, this);
		m_pSWBody->AddRef();
	}
}

void BodySW::bodysms(DLIList<BodySM*> &bodies) 
{
	if (m_pSWBody)
		SWQueryEngine::instance()->bodysms(m_pSWBody, bodies);
}

void BodySW::lumps(DLIList<Lump*> &lumps)
{
	if (m_pSWBody)
		SWQueryEngine::instance()->lumps(m_pSWPart, m_pSWBody, lumps);
}

void BodySW::shellsms(DLIList<ShellSM*> &shellsms)
{
	if (m_pSWBody)
		SWQueryEngine::instance()->shellsms(m_pSWPart, m_pSWBody, shellsms);
}

void BodySW::surfaces(DLIList<Surface*> &surfaces)
{
	if (m_pSWBody)
		SWQueryEngine::instance()->surfaces(m_pSWPart, m_pSWBody, surfaces);
}

void BodySW::loopsms(DLIList<LoopSM*> &loopsms)
{
	if (m_pSWBody)
		SWQueryEngine::instance()->loopsms(m_pSWPart, m_pSWBody, loopsms);
}

void BodySW::curves(DLIList<Curve*> &curves)
{
	if (m_pSWBody)
		SWQueryEngine::instance()->curves(m_pSWPart, m_pSWBody, curves);
}

void BodySW::coedgesms(DLIList<CoEdgeSM*> &coedgesms)
{
	if (m_pSWBody)
		SWQueryEngine::instance()->coedgesms(m_pSWPart, m_pSWBody, coedgesms);
}

void BodySW::points(DLIList<Point*> &points)
{
	if (m_pSWBody)
        SWQueryEngine::instance()->points(m_pSWPart, m_pSWBody, points);
}

// ********** END PUBLIC FUNCTIONS         **********

// ********** BEGIN PROTECTED FUNCTIONS    **********
// ********** END PROTECTED FUNCTIONS      **********

// ********** BEGIN PRIVATE FUNCTIONS      **********

//CubitStatus BodySW::transform_BODY(double* transform, LPCOMPONENT comp)
//{
//    assert(false);
//    return CUBIT_FAILURE;
//}


// ********** END PRIVATE FUNCTIONS        **********

// ********** BEGIN HELPER CLASSES         **********
// ********** END HELPER CLASSES           **********

// ********** BEGIN EXTERN FUNCTIONS       **********
// ********** END EXTERN FUNCTIONS         **********

// ********** BEGIN STATIC FUNCTIONS       **********
// ********** END STATIC FUNCTIONS         **********

CubitStatus BodySW::get_transforms( CubitTransformMatrix &tfm )
{
    assert(false);
    return CUBIT_FAILURE;
}

//void BodySW::get_mass_props(CubitVector &cofg)
//{
//    LPBODY pSWBody = get_BODY_ptr();
//    assert(pSWBody);

//    VARIANT vProps;
//    double density = 1.0;
//    HRESULT hr = pSWBody->GetMassProperties(density, &vProps);

//    double dProps[12];
//    hr = ExtractVARIANTArrayData(vProps, 12, dProps);
//    hr = VariantClear(&vProps);

//    CubitVector partCG(dProps[0], dProps[1], dProps[2]);

//    m_pSWPart->transformPartToAssembly(partCG, cofg);
//}
CubitStatus BodySW::mass_properties( CubitVector& centroid, double& volume )
{
  //see above for partial implementation
  assert(false);
  return CUBIT_FAILURE;
}

CubitPointContainment BodySW::point_containment( const CubitVector& pos )
{
  assert(false);
  return CUBIT_PNT_UNKNOWN;
}

void BodySW::get_parents_virt(DLIList<TopologyBridge*> &parents )
{
  // no parents to add
}

void BodySW::get_children_virt(DLIList<TopologyBridge*> &children )
{
  this->m_pSWPart->lumps(children);
}
