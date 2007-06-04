//-------------------------------------------------------------------------
// Filename      : SWBridge.cpp
//
// Purpose       : Many functions are identical for each SW-specific
//                 TopologyBridge.  This class implements those functions.
//
// Creator       : Joel Kopp, Darryl Melander
//
// Creation Date : 8/30/00
//
// Owner         : Joel Kopp, Darryl Melander
//-------------------------------------------------------------------------

#include "stdafx.h"

//#include <amapp.h>
#import "C:\Program Files\SolidWorks\sldworks.tlb" raw_interfaces_only, raw_native_types, no_namespace, named_guids     //Import the SolidWorks type library

#include "SWBridge.hpp"
#include "SWQueryEngine.hpp"
#include "SWPart.hpp"
#include "CubitSimpleAttrib.hpp"

SWBridge::SWBridge()
{
}

SWBridge::~SWBridge()
{
    remove_all_simple_attribute_virt();
}

void SWBridge::append_simple_attribute_virt(CubitSimpleAttrib* attrib_ptr)
{
	// First, remove any like-typed attributes
	remove_simple_attribute_virt(attrib_ptr);


    // store the attribute locally
    CubitSimpleAttrib *pSimpleAttrib = new CubitSimpleAttrib(*attrib_ptr);

    m_simpleAttribList.append(pSimpleAttrib);

    return;
}

void SWBridge::remove_simple_attribute_virt(CubitSimpleAttrib* attrib_ptr)
{
    CubitSimpleAttrib *pAttrib;
    CubitString typeString;

    if (0 == m_simpleAttribList.size())
        return;

    typeString = attrib_ptr->string_data_list()->get()->c_str();


    int nRemoved = 0;
    m_simpleAttribList.reset();
    while (pAttrib = m_simpleAttribList.get_and_step())
    {
        if (0 == strcmp(typeString.c_str(), pAttrib->string_data_list()->get()->c_str()))
        {
            delete pAttrib;
            m_simpleAttribList.remove();
            nRemoved++;
        }
    }
    // there should never be more than one of the same type in the list
    // since we check when adding
    assert(1 >= nRemoved);
    return;
}

void SWBridge::remove_all_simple_attribute_virt()
{
    while( m_simpleAttribList.size() )
      delete m_simpleAttribList.pop();
}

CubitStatus SWBridge::get_simple_attribute(DLIList<CubitSimpleAttrib*>&
                                           cubit_simple_attrib_list)
{
  m_simpleAttribList.reset();
  for ( int i = m_simpleAttribList.size(); i--; )
    cubit_simple_attrib_list.append( 
      new CubitSimpleAttrib(m_simpleAttribList.get_and_step() ));

	return CUBIT_SUCCESS;
}

CubitStatus SWBridge::get_simple_attribute(const CubitString& name,
                       DLIList<CubitSimpleAttrib*>& cubit_simple_attrib_list)
{
  m_simpleAttribList.reset();
  for ( int i = m_simpleAttribList.size(); i--; )
  {
    CubitSimpleAttrib* csa_ptr = m_simpleAttribList.get_and_step();
    if ( csa_ptr->character_type() == name )
      cubit_simple_attrib_list.append( new CubitSimpleAttrib(csa_ptr) );
  }

	return CUBIT_SUCCESS;
}

/*void SWBridge::print_attribs(LPENTITY entity) 
{
	//- print all the attributes info for this entity

	DLCubitSimpleAttribList csa_list;

	SWBridge::get_simple_attribute(entity, csa_list);

	int i;
	PRINT_INFO("Entity attributes:\n");
	for (i = csa_list.size(); i > 0; i--)
		csa_list.get_and_step()->print();
}*/
/*
SWQueryEngine *SWBridge::get_SW_query_engine() const
{
  SWQueryEngine *sqe = SWQueryEngine::instance();

    // if this function gets called, it needs to return non-NULL
  assert(sqe);

  return sqe;
}
*/
