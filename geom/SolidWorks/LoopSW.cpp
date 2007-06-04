//-------------------------------------------------------------------------
// Filename      : LoopSW.cpp
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

#include "SWQueryEngine.hpp"
#include "LoopSW.hpp"
#include "SWPart.hpp"


//-------------------------------------------------------------------------
// Purpose       : A constructor with a pointer to a SW LOOP.
//
// Special Notes :
//
// Creator       : Joel Kopp, Xuechen Liu
//
// Creation Date : 8/9/00
//-------------------------------------------------------------------------
LoopSW::LoopSW(SWPart *pPart)
{
	m_pSWLoop = NULL;
	m_pSWPart = pPart;
}

//-------------------------------------------------------------------------
// Purpose       : The destructor.
//
// Special Notes :
//
// Creator       : Joel Kopp, Raikanta Sahu
//
// Creation Date : 8/9/00
//-------------------------------------------------------------------------
LoopSW::~LoopSW()
{
	set_LOOP_ptr(NULL);
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
GeometryQueryEngine* LoopSW::get_geometry_query_engine() const
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
// Creation Date : 8/9/00
//-------------------------------------------------------------------------
void LoopSW::append_simple_attribute_virt(CubitSimpleAttrib* csattrib_ptr)
{
    m_pSWPart->append_simple_attribute_virt(m_pSWLoop, csattrib_ptr);
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
// Creation Date : 8/9/00
//-------------------------------------------------------------------------
void LoopSW::remove_simple_attribute_virt(CubitSimpleAttrib* csattrib_ptr)
{
    m_pSWPart->remove_simple_attribute_virt(m_pSWLoop, csattrib_ptr);
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
// Creation Date : 8/9/00
//-------------------------------------------------------------------------
void LoopSW::remove_all_simple_attribute_virt()
{
	m_pSWPart->remove_all_simple_attribute_virt(m_pSWLoop);
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
// Creation Date : 8/9/00
//-------------------------------------------------------------------------
CubitStatus LoopSW::get_simple_attribute(DLIList<CubitSimpleAttrib*>&
                                           cubit_simple_attrib_list)
{
	return m_pSWPart->get_simple_attribute(m_pSWLoop, cubit_simple_attrib_list);
}
CubitStatus LoopSW::get_simple_attribute(const CubitString& name,
                                   DLIList<CubitSimpleAttrib*>& list)
  { return m_pSWPart->get_simple_attribute(m_pSWLoop, name, list); }

ILoop2 *LoopSW::get_LOOP_ptr() const
{
	return m_pSWLoop;
}

void LoopSW::set_LOOP_ptr(ILoop2 *loop)
{
	if (loop == m_pSWLoop)
		return;

	if (m_pSWLoop)
	{
		m_pSWPart->remove_cubit_owner(m_pSWLoop);
        m_pSWLoop->Release();
	}

	m_pSWLoop = loop;

	if (m_pSWLoop)
	{
		m_pSWPart->set_cubit_owner(m_pSWLoop, this);
        m_pSWLoop->AddRef();
	}
}

//-------------------------------------------------------------------------
// Purpose       : merges "this" with input dead_OSMEPtr
//
// Special Notes :
//
// Creator       : Joel Kopp, Jihong Ma
//
// Creation Date : 8/9/00
//-------------------------------------------------------------------------
//CubitStatus LoopSW::merge( OtherSolidModelEntity* dead_OSMEPtr)
//{
//    assert(false);
//	return CUBIT_FAILURE;
//}

//TopologyEntity* LoopSW::unmerge( DLIList<RefVolume*> /*volumes*/ )
//{
//    assert(false);
//  PRINT_ERROR("Unmerging of SW-based entities is currently disabled.\n");
//  return NULL;
//}

void LoopSW::bodysms(DLIList<BodySM*> &bodies) 
{
	if (m_pSWLoop)
		SWQueryEngine::instance()->bodysms(m_pSWLoop, bodies);
}

void LoopSW::lumps(DLIList<Lump*> &lumps)
{
	if (m_pSWLoop)
		SWQueryEngine::instance()->lumps(m_pSWPart, m_pSWLoop, lumps);
}

void LoopSW::shellsms(DLIList<ShellSM*> &shellsms)
{
	if (m_pSWLoop)
		SWQueryEngine::instance()->shellsms(m_pSWPart, m_pSWLoop, shellsms);
}

void LoopSW::surfaces(DLIList<Surface*> &surfaces)
{
	if (m_pSWLoop)
		SWQueryEngine::instance()->surfaces(m_pSWPart, m_pSWLoop, surfaces);
}

void LoopSW::loopsms(DLIList<LoopSM*> &loopsms)
{
	if (m_pSWLoop)
		SWQueryEngine::instance()->loopsms(m_pSWPart, m_pSWLoop, loopsms);
}

void LoopSW::curves(DLIList<Curve*> &curves)
{
	if (m_pSWLoop)
		SWQueryEngine::instance()->curves(m_pSWPart, m_pSWLoop, curves);
}

void LoopSW::coedgesms(DLIList<CoEdgeSM*> &coedgesms)
{
	if (m_pSWLoop)
		SWQueryEngine::instance()->coedgesms(m_pSWPart, m_pSWLoop, coedgesms);
}

void LoopSW::points(DLIList<Point*> &points)
{
	if (m_pSWLoop)
        SWQueryEngine::instance()->points(m_pSWPart, m_pSWLoop, points);
}

void LoopSW::get_parents_virt(DLIList<TopologyBridge*> &parents )
{
  IFace2 *swface;
  HRESULT hr = m_pSWLoop->IGetFace(&swface);
  if (SUCCEEDED(hr) && swface)
  {
    TopologyBridge *owner = m_pSWPart->cubit_owner(swface);
    if (owner)
      parents.append(owner);
  }
}

void LoopSW::get_children_virt(DLIList<TopologyBridge*> &children )
{
    ICoEdge       *pSWCoEdge = NULL;
    IEnumCoEdges  *pCoEdgeEnum = NULL;
    long nFetched = 0;

    HRESULT hr = m_pSWLoop->EnumCoEdges(&pCoEdgeEnum);
    if (SUCCEEDED(hr) && pCoEdgeEnum)
    {
        while (S_OK == pCoEdgeEnum->Next(1, &pSWCoEdge, &nFetched))
        {
            TopologyBridge *owner = m_pSWPart->cubit_owner(pSWCoEdge);
            if (owner)
              children.append(owner);
        }
        pCoEdgeEnum->Release();
    }
}

