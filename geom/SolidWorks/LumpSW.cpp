//-------------------------------------------------------------------------
// Filename      : LumpSW.cpp
//
// Purpose       : 
//
// Special Notes :
//
// Creator       : Joel Kopp, Xuechen Liu
//
// Creation Date : 8/22/00
//
// Owner         : Joel Kopp, Malcolm J. Panthaki
//-------------------------------------------------------------------------

// Precompiled header
#include "stdafx.h"

//#include <amapp.h>
#import "C:\Program Files\SolidWorks\sldworks.tlb" raw_interfaces_only, raw_native_types, no_namespace, named_guids     //Import the SolidWorks type library
//#include <swconst.h>

#include <assert.h>

#include "SWQueryEngine.hpp"
#include "LumpSW.hpp"
#include "SWPart.hpp"
#include "BodySM.hpp"


extern HRESULT ExtractVARIANTArrayData(VARIANT &v, long lArraySize, double *dArray);


//-------------------------------------------------------------------------
// Purpose       : The constructor with a pointer to a body. 
//
// Special Notes :
//
// Creator       : Joel Kopp, Xuechen Liu
//
// Creation Date : 8/22/00
//-------------------------------------------------------------------------
LumpSW::LumpSW(SWPart *pPart)
{
    // Calculate a bounding box if there isn't one already
    //SWGeometryEngine::instance()->bounding_box(LUMPPtr);
	m_pSWBody = NULL;
	m_pSWPart = pPart;
}

//-------------------------------------------------------------------------
// Purpose       : The destructor
//
// Special Notes :
//
// Creator       : Joel Kopp, Raikanta Sahu
//
// Creation Date : 8/22/00
//-------------------------------------------------------------------------
LumpSW::~LumpSW()
{
	set_BODY_ptr(NULL);
}

//-------------------------------------------------------------------------
// Purpose       : Return a pointer to the body associated with
//                 the object.
//
// Special Notes :
//
// Creator       : Joel Kopp, Raikanta Sahu
//
// Creation Date : 8/22/00
//-------------------------------------------------------------------------
IBody2 *LumpSW::get_BODY_ptr() const
{
	return m_pSWBody;
}

void LumpSW::set_BODY_ptr(IBody2 *lump)
{
	if (lump == m_pSWBody)
		return;

	if (m_pSWBody)
	{
		m_pSWPart->remove_cubit_owner(m_pSWBody);
        m_pSWBody->Release();
	}

	m_pSWBody = lump;

	if (m_pSWBody)
	{
		m_pSWPart->set_cubit_owner(m_pSWBody, this);
        m_pSWBody->AddRef();
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
GeometryQueryEngine* LumpSW::get_geometry_query_engine() const
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
// Creation Date : 8/22/00
//-------------------------------------------------------------------------
void LumpSW::append_simple_attribute_virt(CubitSimpleAttrib* csattrib_ptr)
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
// Creation Date : 8/22/00
//-------------------------------------------------------------------------
void LumpSW::remove_simple_attribute_virt(CubitSimpleAttrib* csattrib_ptr)
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
// Creation Date : 8/23/00
//-------------------------------------------------------------------------
void LumpSW::remove_all_simple_attribute_virt()
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
// Creation Date : 8/23/00
//-------------------------------------------------------------------------
CubitStatus LumpSW::get_simple_attribute(DLIList<CubitSimpleAttrib*>&
                                           cubit_simple_attrib_list)
{
    return m_pSWPart->get_simple_attribute(this, cubit_simple_attrib_list);
}
CubitStatus LumpSW::get_simple_attribute(const CubitString& name,
                                   DLIList<CubitSimpleAttrib*>& list)
  { return m_pSWPart->get_simple_attribute(this, name, list); }

//-------------------------------------------------------------------------
// Purpose       : Get the bounding box of the object.
//
// Special Notes :
//
// Creator       : Joel Kopp, jihong Ma
//
// Creation Date : 8/23/00
//-------------------------------------------------------------------------
CubitBox LumpSW::bounding_box() const 
{
    HRESULT hr = NOERROR;
    IBody2 *pSWBody = get_BODY_ptr();

    VARIANT vCoords;
    hr = pSWBody->GetBodyBox(&vCoords);

	double dSWbox[6];
    hr = ExtractVARIANTArrayData(vCoords, 6, dSWbox);
    hr = VariantClear(&vCoords);

    // transform the box to the assembly coordinate space
    CubitVector corner1(dSWbox[0], dSWbox[1], dSWbox[2]);
    CubitVector corner2(dSWbox[3], dSWbox[4], dSWbox[5]);

    CubitVector asmCorner1;
    CubitVector asmCorner2;
    m_pSWPart->transformPartToAssembly(corner1, asmCorner1);
    m_pSWPart->transformPartToAssembly(corner2, asmCorner2);

	// Convert to a CubitBox and return it
    return CubitBox(asmCorner1, asmCorner2);
}

//-------------------------------------------------------------------------
// Purpose       : Returns the volume of the Lump
//
// Special Notes :
//
// Creator       : Joel Kopp, Raikanta Sahu
//
// Creation Date : 8/23/00
//-------------------------------------------------------------------------
double LumpSW::measure()
{
	IBody2 *lump = get_BODY_ptr();

	// Get the mass properties.  The 1.0 is the density, which is
	// unimportant since all we need is the volume
    VARIANT vMassProps;
	HRESULT hr = lump->GetMassProperties ( 1.0, &vMassProps );
    assert(SUCCEEDED(hr));

	double massProps[12];
    hr = ExtractVARIANTArrayData(vMassProps, 12, massProps);
    hr = VariantClear(&vMassProps);

	double volume = massProps[3];

	return volume;
}

void LumpSW::bodysms(DLIList<BodySM*> &bodies) 
{
	if (m_pSWBody)
		SWQueryEngine::instance()->bodysms(m_pSWBody, bodies);
}

void LumpSW::lumps(DLIList<Lump*> &lumps)
{
	if (m_pSWBody)
		SWQueryEngine::instance()->lumps(m_pSWPart, m_pSWBody, lumps);
}

void LumpSW::shellsms(DLIList<ShellSM*> &shellsms)
{
	if (m_pSWBody)
		SWQueryEngine::instance()->shellsms(m_pSWPart, m_pSWBody, shellsms);
}

void LumpSW::surfaces(DLIList<Surface*> &surfaces)
{
	if (m_pSWBody)
		SWQueryEngine::instance()->surfaces(m_pSWPart, m_pSWBody, surfaces);
}

void LumpSW::loopsms(DLIList<LoopSM*> &loopsms)
{
	if (m_pSWBody)
		SWQueryEngine::instance()->loopsms(m_pSWPart, m_pSWBody, loopsms);
}

void LumpSW::curves(DLIList<Curve*> &curves)
{
	if (m_pSWBody)
		SWQueryEngine::instance()->curves(m_pSWPart, m_pSWBody, curves);
}

void LumpSW::coedgesms(DLIList<CoEdgeSM*> &coedgesms)
{
	if (m_pSWBody)
		SWQueryEngine::instance()->coedgesms(m_pSWPart, m_pSWBody, coedgesms);
}

void LumpSW::points(DLIList<Point*> &points)
{
	if (m_pSWBody)
        SWQueryEngine::instance()->points(m_pSWPart, m_pSWBody, points);
}

void LumpSW::get_parents_virt(DLIList<TopologyBridge*> &parents )
{
  BodySM* bodysm = m_pSWPart->getBody();

  parents.append(bodysm);
}

void LumpSW::get_children_virt(DLIList<TopologyBridge*> &children )
{
  this->m_pSWPart->shells(children);
}

CubitStatus LumpSW::mass_properties( CubitVector &centroid, double &volume )
{
  assert(false);
  return CUBIT_FAILURE;
}
