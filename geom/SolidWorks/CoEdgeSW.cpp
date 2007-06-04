//-------------------------------------------------------------------------
// Filename      : CoEdgeSW.cpp
//
// Purpose       : 
//
// Special Notes :
//
// Creator       : Joel Kopp, Xuechen Liu
//
// Creation Date : 8/9/00
//
// Owner         : Joel Kopp, Malcolm J. Panthaki
//-------------------------------------------------------------------------


// Precompiled header
#include "stdafx.h"

//#include <amapp.h>
#import "C:\Program Files\SolidWorks\sldworks.tlb" raw_interfaces_only, raw_native_types, no_namespace, named_guids     //Import the SolidWorks type library

#include "CoEdgeSW.hpp"
#include "SWQueryEngine.hpp"
#include "SWPart.hpp"


//-------------------------------------------------------------------------
// Purpose       : A constructor with a pointer to a SW COEDGE.
//
// Special Notes :
//
// Creator       : Joel Kopp, Xuechen Liu
//
// Creation Date : 8/9/00
//-------------------------------------------------------------------------
CoEdgeSW::CoEdgeSW(SWPart *pPart)
{
	m_pSWCoEdge = NULL;
	m_pSWPart = pPart;
}

//-------------------------------------------------------------------------
// Purpose       : The destructor
//
// Special Notes :
//
// Creator       : Joel Kopp, Raikanta Sahu
//
// Creation Date : 8/9/00
//-------------------------------------------------------------------------
CoEdgeSW::~CoEdgeSW()
{
    set_COEDGE_ptr(NULL);
}

ICoEdge *CoEdgeSW::get_COEDGE_ptr() const
{
	return m_pSWCoEdge;
}

void CoEdgeSW::set_COEDGE_ptr(ICoEdge *coedge)
{
	if (coedge == m_pSWCoEdge)
		return;

	if (m_pSWCoEdge)
	{
		m_pSWPart->remove_cubit_owner(m_pSWCoEdge);
        m_pSWCoEdge->Release();
	}

	m_pSWCoEdge = coedge;

	if (m_pSWCoEdge)
	{
		m_pSWPart->set_cubit_owner(m_pSWCoEdge, this);
        m_pSWCoEdge->AddRef();
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
GeometryQueryEngine* CoEdgeSW::get_geometry_query_engine() const
{
	return SWQueryEngine::instance();   
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
// Creation Date : 8/10/00
//-------------------------------------------------------------------------
void CoEdgeSW::append_simple_attribute_virt(CubitSimpleAttrib* csattrib_ptr)
{
	m_pSWPart->append_simple_attribute_virt(m_pSWCoEdge, csattrib_ptr);
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
// Creation Date : 8/10/00
//-------------------------------------------------------------------------
void CoEdgeSW::remove_simple_attribute_virt(CubitSimpleAttrib* csattrib_ptr)
{
	m_pSWPart->remove_simple_attribute_virt(m_pSWCoEdge, csattrib_ptr);
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
// Creation Date : 8/10/00
//-------------------------------------------------------------------------
void CoEdgeSW::remove_all_simple_attribute_virt()
{
    m_pSWPart->remove_all_simple_attribute_virt(m_pSWCoEdge);
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
// Creation Date : 8/10/00
//-------------------------------------------------------------------------
CubitStatus CoEdgeSW::get_simple_attribute(DLIList<CubitSimpleAttrib*>&
                                             cubit_simple_attrib_list)
{
    return m_pSWPart->get_simple_attribute(m_pSWCoEdge, cubit_simple_attrib_list);
}
CubitStatus CoEdgeSW::get_simple_attribute(const CubitString& name,
                                   DLIList<CubitSimpleAttrib*>& list)
{
  return m_pSWPart->get_simple_attribute(m_pSWCoEdge, name, list);
}

CubitSense CoEdgeSW::sense()
{
    HRESULT hr = NOERROR;
	ICoEdge *coedge = get_COEDGE_ptr();

	VARIANT_BOOL sense = FALSE;
	hr = coedge->GetSense( &sense );
	if (sense == TRUE) // if coedge and edge have same sense
		return CUBIT_FORWARD;
	else if (sense == FALSE)
		return CUBIT_REVERSED;
	else
    {
        assert(false);
		return CUBIT_UNKNOWN;
    }
}

void CoEdgeSW::bodysms(DLIList<BodySM*> &bodies) 
{
	if (m_pSWCoEdge)
		SWQueryEngine::instance()->bodysms(m_pSWCoEdge, bodies);
}

void CoEdgeSW::lumps(DLIList<Lump*> &lumps)
{
	if (m_pSWCoEdge)
		SWQueryEngine::instance()->lumps(m_pSWPart, m_pSWCoEdge, lumps);
}

void CoEdgeSW::shellsms(DLIList<ShellSM*> &shellsms)
{
	if (m_pSWCoEdge)
		SWQueryEngine::instance()->shellsms(m_pSWPart, m_pSWCoEdge, shellsms);
}

void CoEdgeSW::surfaces(DLIList<Surface*> &surfaces)
{
	if (m_pSWCoEdge)
		SWQueryEngine::instance()->surfaces(m_pSWPart, m_pSWCoEdge, surfaces);
}

void CoEdgeSW::loopsms(DLIList<LoopSM*> &loopsms)
{
	if (m_pSWCoEdge)
		SWQueryEngine::instance()->loopsms(m_pSWPart, m_pSWCoEdge, loopsms);
}

void CoEdgeSW::curves(DLIList<Curve*> &curves)
{
	if (m_pSWCoEdge)
		SWQueryEngine::instance()->curves(m_pSWPart, m_pSWCoEdge, curves);
}

void CoEdgeSW::coedgesms(DLIList<CoEdgeSM*> &coedgesms)
{
	if (m_pSWCoEdge)
		SWQueryEngine::instance()->coedgesms(m_pSWPart, m_pSWCoEdge, coedgesms);
}

void CoEdgeSW::points(DLIList<Point*> &points)
{
	if (m_pSWCoEdge)
        SWQueryEngine::instance()->points(m_pSWPart, m_pSWCoEdge, points);
}

void CoEdgeSW::get_parents_virt(DLIList<TopologyBridge*> &parents )
{
  ILoop2 *swloop;
  HRESULT hr = m_pSWCoEdge->IGetLoop2(&swloop);
  if (SUCCEEDED(hr) && swloop)
  {
    TopologyBridge *owner = m_pSWPart->cubit_owner(swloop);
    if (owner)
      parents.append(owner);
  }
}

void CoEdgeSW::get_children_virt(DLIList<TopologyBridge*> &children )
{
  IEdge *swedge;
  HRESULT hr = m_pSWCoEdge->IGetEdge(&swedge);
  if (SUCCEEDED(hr) && swedge)
  {
    TopologyBridge *owner = m_pSWPart->cubit_owner(swedge);
    if (owner)
      children.append(owner);
  }
}
