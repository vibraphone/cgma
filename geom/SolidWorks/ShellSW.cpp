//-------------------------------------------------------------------------
// Filename      : ShellSW.cpp
//
// Purpose       : 
//
// Special Notes :
//
// Creator       : Joel Kopp, Xuechen Liu
//
// Creation Date : 8/23/00
//
// Owner         : Joel Kopp, Malcolm J. Panthaki
//-------------------------------------------------------------------------

// Precompiled header
#include "stdafx.h"

//#include <amapp.h>
#import "C:\Program Files\SolidWorks\sldworks.tlb" raw_interfaces_only, raw_native_types, no_namespace, named_guids     //Import the SolidWorks type library

#include "SWQueryEngine.hpp"
#include "ShellSW.hpp"
#include "SWPart.hpp"

#include "Surface.hpp"


//-------------------------------------------------------------------------
// Purpose       : A constructor with a pointer to a SW SHELL (from faces).
//                 There are no shells in SW.
//
// Special Notes :
//
// Creator       : Joel Kopp, Xuechen Liu
//
// Creation Date : 8/23/00
//-------------------------------------------------------------------------
ShellSW::ShellSW(SWPart *pPart)
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
// Creation Date : 8/23/00
//-------------------------------------------------------------------------
ShellSW::~ShellSW()
{
	set_BODY_ptr(NULL);
}

void ShellSW::set_BODY_ptr(IBody2 *shell)
{
	if (shell == m_pSWBody)
		return;

	if (m_pSWBody)
	{
		m_pSWPart->remove_cubit_owner(m_pSWBody);
        m_pSWBody->Release();
	}

	m_pSWBody = shell;

	if (m_pSWBody)
	{
		m_pSWPart->set_cubit_owner(m_pSWBody, this);
        m_pSWBody->AddRef();
	}
}

IBody2 *ShellSW::get_BODY_ptr() const
{
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
GeometryQueryEngine* ShellSW::get_geometry_query_engine() const
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
// Creation Date : 8/23/00
//-------------------------------------------------------------------------
void ShellSW::append_simple_attribute_virt(CubitSimpleAttrib* csattrib_ptr)
{
    assert(false);
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
// Creation Date : 8/23/00
//-------------------------------------------------------------------------
void ShellSW::remove_simple_attribute_virt(CubitSimpleAttrib* csattrib_ptr)
{
    assert(false);
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
void ShellSW::remove_all_simple_attribute_virt()
{
    assert(false);
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
CubitStatus ShellSW::get_simple_attribute(DLIList<CubitSimpleAttrib*>&
                                           cubit_simple_attrib_list)
{
    assert(false);
    return CUBIT_FAILURE;
}
CubitStatus ShellSW::get_simple_attribute(const CubitString& ,
                                   DLIList<CubitSimpleAttrib*>& )
  { assert(0); return CUBIT_FAILURE; }

void ShellSW::bodysms(DLIList<BodySM*> &bodies) 
{
	if (m_pSWBody)
		SWQueryEngine::instance()->bodysms(m_pSWBody, bodies);
}

void ShellSW::lumps(DLIList<Lump*> &lumps)
{
	if (m_pSWBody)
		SWQueryEngine::instance()->lumps(m_pSWPart, m_pSWBody, lumps);
}

void ShellSW::shellsms(DLIList<ShellSM*> &shellsms)
{
	if (m_pSWBody)
		SWQueryEngine::instance()->shellsms(m_pSWPart, m_pSWBody, shellsms);
}

void ShellSW::surfaces(DLIList<Surface*> &surfaces)
{
	if (m_pSWBody)
		SWQueryEngine::instance()->surfaces(m_pSWPart, m_pSWBody, surfaces);
}

void ShellSW::loopsms(DLIList<LoopSM*> &loopsms)
{
	if (m_pSWBody)
		SWQueryEngine::instance()->loopsms(m_pSWPart, m_pSWBody, loopsms);
}

void ShellSW::curves(DLIList<Curve*> &curves)
{
	if (m_pSWBody)
		SWQueryEngine::instance()->curves(m_pSWPart, m_pSWBody, curves);
}

void ShellSW::coedgesms(DLIList<CoEdgeSM*> &coedgesms)
{
	if (m_pSWBody)
		SWQueryEngine::instance()->coedgesms(m_pSWPart, m_pSWBody, coedgesms);
}

void ShellSW::points(DLIList<Point*> &points)
{
	if (m_pSWBody)
        SWQueryEngine::instance()->points(m_pSWPart, m_pSWBody, points);
}

void ShellSW::get_parents_virt(DLIList<TopologyBridge*> &parents )
{
  this->m_pSWPart->lumps(parents);
}

void ShellSW::get_children_virt(DLIList<TopologyBridge*> &children )
{
  DLIList<Surface*> surfs;
  surfaces(surfs);

  surfs.reset();
  int i;
  for (i=0; i<surfs.size(); i++)
    children.append(surfs.get_and_step());
}

