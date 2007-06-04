#include "stdafx.h"

#include <algorithm>

#include "SWQueryEngine.hpp"
#include "SWPart.hpp"
#include "CastTo.hpp"

#include "BodySW.hpp"
#include "LumpSW.hpp"
#include "ShellSW.hpp"

#include "gtfordm.h"

extern HRESULT ExtractVARIANTArrayData(VARIANT &v, long lArraySize, double *dArray);

static void rebindKeys(LPDMREF pRef, ATTRIBMAP &oldAttribMap, KEYMAP &oldKeyMap,
                       ATTRIBMAP &newAttribMap, KEYMAP &newKeyMap, ATTRIBMAP &newDanglingMap);

SWPart::SWPart()
{
    m_pSWBody = 0;

    m_body = 0;
	m_lump = 0;
	m_shell = 0;
}

SWPart::~SWPart()
{
    CubitSimpleAttrib* pAttrib = NULL;

    // get the attributes associated with the SW entity passed in
    m_AttribMapIter = m_AttribMap.begin();

    // remove all attributes from the map and delete the simple attrib objects
    int nRemoved = 0;
    while (m_AttribMapIter != m_AttribMap.end())
    {
        pAttrib = (*m_AttribMapIter).second;
        m_AttribMapIter = m_AttribMap.erase(m_AttribMapIter); // returns iterator to next or end
        delete pAttrib;
        nRemoved++;
    }
}

BodySM *SWPart::getBody()
{
  BodySM *p_body = dynamic_cast<BodySM*>(m_body);
  assert(!!p_body);
  return p_body;
}

//void SWPart::setBody(BodySM *pBody)
//{
//    m_body = pBody;
//}

//LumpSM *SWPart::getLump()
//{
//    return m_lump;
//}

//void SWPart::setLump(LumpSM *pLump)
//{
//    m_lump = pLump;
//}

//ShellSM* SWPart::getShell()
//{
//    return m_shell;
//}

//void SWPart::setShell(ShellSM *pShell)
//{
//    m_shell = pShell;
//}

TopologyBridge* SWPart::cubit_owner(IUnknown *entity)
{
    HRESULT hr = NOERROR;

    IBody2 *pSWBody;
    IFace2 *pSWFace;
    ILoop2 *pSWLoop;
    ICoEdge *pSWCoEdge;
    IEdge *pSWEdge;
    IVertex *pSWVertex;

	TopologyBridge *pBridge;

    hr = entity->QueryInterface(IID_IFace2, (LPVOID*)&pSWFace);
    if (SUCCEEDED(hr))
    {
		pBridge = cubit_owner(pSWFace);
		pSWFace->Release();
		return pBridge;
    }
    else
    {
        hr = entity->QueryInterface(IID_ILoop2, (LPVOID*)&pSWLoop);
        if (SUCCEEDED(hr))
        {
			pBridge = cubit_owner(pSWLoop);
			pSWLoop->Release();
			return pBridge;
        }
        else
        {
            hr = entity->QueryInterface(IID_ICoEdge, (LPVOID*)&pSWCoEdge);
            if (SUCCEEDED(hr))
            {
				pBridge = cubit_owner(pSWCoEdge);
				pSWCoEdge->Release();
				return pBridge;
            }
            else
            {
                hr = entity->QueryInterface(IID_IEdge, (LPVOID*)&pSWEdge);
                if (SUCCEEDED(hr))
                {
					pBridge = cubit_owner(pSWEdge);
					pSWEdge->Release();
					return pBridge;
                }
                else
                {
                    hr = entity->QueryInterface(IID_IVertex, (LPVOID*)&pSWVertex);
                    if (SUCCEEDED(hr))
                    {
						pBridge = cubit_owner(pSWVertex);
						pSWVertex->Release();
						return pBridge;
                    }
                    else
                    {
                        hr = entity->QueryInterface(IID_IBody2, (LPVOID*)&pSWBody);
                        if (SUCCEEDED(hr))
                        {
							assert(false);
                            // multiple CGM owners, so return none
                            return NULL;
                        }
                        else
                        {
                            assert(false); // a case not yet handled
                        }
                    }
                }// IEdge *
            }// ICoEdge *
        }// ILoop2 *
    }// IFace2 *

    return NULL;
}

TopologyBridge* SWPart::cubit_owner(IFace2 *entity)
{
    HRESULT hr = NOERROR;

    // quick check for exact match
	m_FaceMapIter = m_FaceMap.find(entity);

	// if not found - check each face using SolidWorks compare function
	// TODO - this is only necessary since going from a loop to a face gives
	//        a different interface pointer.  This expensive check is the only
	//        way to find the match.  Report this to SolidWorks and see if they
	//        will fix it.
	if (m_FaceMapIter == m_FaceMap.end())
	{
        // double check
        VARIANT_BOOL bSame = 0;
		for (m_FaceMapIter = m_FaceMap.begin(); m_FaceMapIter != m_FaceMap.end(); ++m_FaceMapIter)
		{
			IFace2 *pMapFace = (*m_FaceMapIter).first;
			hr = pMapFace->IsSame(entity, &bSame);
			if (bSame)
				break;
		}
	}

	if (m_FaceMapIter == m_FaceMap.end())
        return NULL;
    else
    {
        return (*m_FaceMapIter).second;
    }
}

TopologyBridge* SWPart::cubit_owner(ILoop2 *entity)
{
	m_LoopMapIter = m_LoopMap.find(entity);
	if (m_LoopMapIter == m_LoopMap.end())
		return NULL;
	else
		return (*m_LoopMapIter).second;
}

TopologyBridge* SWPart::cubit_owner(ICoEdge *entity)
{
	m_CoEdgeMapIter = m_CoEdgeMap.find(entity);
	if (m_CoEdgeMapIter == m_CoEdgeMap.end())
		return NULL;
	else
		return (*m_CoEdgeMapIter).second;
}

TopologyBridge* SWPart::cubit_owner(IEdge *entity)
{
	m_EdgeMapIter = m_EdgeMap.find(entity);
	if (m_EdgeMapIter == m_EdgeMap.end())
		return NULL;
	else
		return (*m_EdgeMapIter).second;
}

TopologyBridge* SWPart::cubit_owner(IVertex *entity)
{
	m_VertexMapIter = m_VertexMap.find(entity);
	if (m_VertexMapIter == m_VertexMap.end())
		return NULL;
	else
		return (*m_VertexMapIter).second;
}



void SWPart::lumps(DLIList<TopologyBridge*> &tb_list)
{
    tb_list.append(CAST_TO(m_lump, TopologyBridge));
    return;
}

void SWPart::shells(DLIList<TopologyBridge*> &tb_list)
{
    tb_list.append(CAST_TO(m_shell, TopologyBridge));
    return;
}
//void SWPart::cubit_owner( enum EntityType ownerType,
//                          DLIList<TopologyBridge*> &tb_list )
//{
//    HRESULT hr = NOERROR;
//    LPBODY pSWBody = NULL;

//    switch (ownerType)
//    {
//    case Lump_TYPE:

//        tb_list.append(CAST_TO(m_lump, TopologyBridge));
//        return;

//    case Shell_TYPE:
//        tb_list.append(CAST_TO(m_shell, TopologyBridge));
//        return;

//    default:
//        assert(false); // not handled
//    }

//    return;
//}

void SWPart::cubit_owner(DLIList<IUnknown *> &entity_list,
                 DLIList<TopologyBridge*> &tb_list)
{
	for (int i = entity_list.size(); i > 0; i--)
	{
		IUnknown *entity = entity_list.get_and_step();
		TopologyBridge *tb = cubit_owner(entity);
		if (tb != NULL)
			tb_list.append(CAST_TO(tb, TopologyBridge));
	}
}

void SWPart::set_cubit_owner(IUnknown *entity, TopologyBridge *cubit_entity)
{
    HRESULT hr = NOERROR;

    // Some checks
	assert(cubit_entity);
	assert(entity);

    IBody2 *pSWBody = NULL;
    hr = entity->QueryInterface(IID_IBody2, (LPVOID*)&pSWBody);
    if (SUCCEEDED(hr) && pSWBody)
    {

        BodySW *bodySW = CAST_TO(cubit_entity, BodySW);
        if (bodySW)
        {
            assert(NULL == m_pSWBody);
            m_pSWBody = pSWBody;
            m_body = cubit_entity;

            return;
        }

        LumpSW *lump = CAST_TO(cubit_entity, LumpSW);
        if (lump)
        {
            assert(m_pSWBody);
            m_lump = cubit_entity;

            return;
        }

        ShellSW *shell = CAST_TO(cubit_entity, ShellSW);
        if (shell)
        {
            assert(m_pSWBody);
            m_shell = cubit_entity;

            return;
        }

        pSWBody->Release(); // not holding ref count from part, only from BodySW
    }

    IFace2 *pSWFace;
    ILoop2 *pSWLoop;
    ICoEdge *pSWCoEdge;
    IEdge *pSWEdge;
    IVertex *pSWVertex;

    hr = entity->QueryInterface(IID_IFace2, (LPVOID*)&pSWFace);
    if (SUCCEEDED(hr))
    {
        std::pair<FACEMAP::iterator, bool> res;
        res = m_FaceMap.insert(FACEMAP::value_type(pSWFace, cubit_entity));
        assert(res.second);
        pSWFace->Release(); // map not holding ref count
    }
    else
    {
        hr = entity->QueryInterface(IID_ILoop2, (LPVOID*)&pSWLoop);
        if (SUCCEEDED(hr))
        {
            std::pair<LOOPMAP::iterator, bool> res;
            res = m_LoopMap.insert(LOOPMAP::value_type(pSWLoop, cubit_entity));
            assert(res.second);
            pSWLoop->Release(); // map not holding ref count
        }
        else
        {
            hr = entity->QueryInterface(IID_ICoEdge, (LPVOID*)&pSWCoEdge);
            if (SUCCEEDED(hr))
            {
                std::pair<COEDGEMAP::iterator, bool> res;
                res = m_CoEdgeMap.insert(COEDGEMAP::value_type(pSWCoEdge, cubit_entity));
                assert(res.second);
                pSWCoEdge->Release(); // map not holding ref count
            }
            else
            {
                hr = entity->QueryInterface(IID_IEdge, (LPVOID*)&pSWEdge);
                if (SUCCEEDED(hr))
                {
                    std::pair<EDGEMAP::iterator, bool> res;
                    res = m_EdgeMap.insert(EDGEMAP::value_type(pSWEdge, cubit_entity));
                    assert(res.second);
                    pSWEdge->Release(); // map not holding ref count
                }
                else
                {
                    hr = entity->QueryInterface(IID_IVertex, (LPVOID*)&pSWVertex);
                    if (SUCCEEDED(hr))
                    {
                        std::pair<VERTEXMAP::iterator, bool> res;
                        res = m_VertexMap.insert(VERTEXMAP::value_type(pSWVertex, cubit_entity));
                        assert(res.second);
                        pSWVertex->Release(); // map not holding ref count
                    }
                    else
                    {
                            assert(false);
                    } // IVertex *
                }// IEdge *
            }// ICoEdge *
        }// ILoop2 *
    }// IFace2 *

    return;
}

void
SWPart::remove_cubit_owner(IUnknown *entity)
{
    HRESULT hr = NOERROR;
    IBody2 *pSWBody = NULL;
    IFace2 *pSWFace = NULL;
    ILoop2 *pSWLoop = NULL;
    ICoEdge *pSWCoEdge = NULL;
    IEdge *pSWEdge = NULL;
    IVertex *pSWVertex = NULL;

    hr = entity->QueryInterface(IID_IBody2, (LPVOID*)&pSWBody);
    if (SUCCEEDED(hr))
    {
        // the body has multiple cubit owners - body, lump, and shell
        // so set the SW body and all owners to NULL

        if(m_pSWBody) // only do this once
        {
            m_pSWBody = NULL;
            m_body = NULL;
            m_lump = NULL;
            m_shell = NULL;
        }

        pSWBody->Release();
    }
    else
    {
        TopologyBridge* pBridge = NULL;
        pBridge = cubit_owner(entity);
        if (pBridge)
        {
            hr = entity->QueryInterface(IID_IFace2, (LPVOID*)&pSWFace);
            if (SUCCEEDED(hr))
            {
                FACEMAP::size_type nRemoved;
                nRemoved = m_FaceMap.erase(pSWFace);
                assert(1 == nRemoved);
                pSWFace->Release();
            }
            else
            {
                hr = entity->QueryInterface(IID_ILoop2, (LPVOID*)&pSWLoop);
                if (SUCCEEDED(hr))
                {
                    LOOPMAP::size_type nRemoved;
                    nRemoved = m_LoopMap.erase(pSWLoop);
                    assert(1 == nRemoved);
                    pSWLoop->Release();
                }
                else
                {
                    hr = entity->QueryInterface(IID_ICoEdge, (LPVOID*)&pSWCoEdge);
                    if (SUCCEEDED(hr))
                    {
                        COEDGEMAP::size_type nRemoved;
                        nRemoved = m_CoEdgeMap.erase(pSWCoEdge);
                        assert(1 == nRemoved);
                        pSWCoEdge->Release();
                    }
                    else
                    {
                        hr = entity->QueryInterface(IID_IEdge, (LPVOID*)&pSWEdge);
                        if (SUCCEEDED(hr))
                        {
                            EDGEMAP::size_type nRemoved;
                            nRemoved = m_EdgeMap.erase(pSWEdge);
                            assert(1 == nRemoved);
                            pSWEdge->Release();
                        }
                        else
                        {
                            hr = entity->QueryInterface(IID_IVertex, (LPVOID*)&pSWVertex);
                            if (SUCCEEDED(hr))
                            {
                                VERTEXMAP::size_type nRemoved;
                                nRemoved = m_VertexMap.erase(pSWVertex);
                                assert(1 == nRemoved);
                                pSWVertex->Release();
                            }
                            else
                            {
                                assert(false); // type not supported
                            } // vertex
                        } // edge
                    } // coedge
                } // loop
            } // face

    //        pBridge->ENTITY_ptr(NULL);
        }
    }
}

void
SWPart::addPseudoVertex(TopologyBridge* pPseudoVertex, IEdge *pClosedSWEdge)
{
    std::pair<EDGEMAP::iterator, bool> res;
    res = m_PseudoVertexMap.insert(EDGEMAP::value_type(pClosedSWEdge, pPseudoVertex));
    assert(res.second);
    pClosedSWEdge->Release(); // map not holding ref count
}

TopologyBridge *SWPart::findPseudoVertex(IEdge *pClosedSWEdge)
{
	m_EdgeMapIter = m_PseudoVertexMap.find(pClosedSWEdge);
	if (m_EdgeMapIter == m_EdgeMap.end())
		return NULL;
	else
		return (*m_EdgeMapIter).second;
}

void SWPart::removePseudoVertex(IEdge *pClosedSWEdge)
{
    EDGEMAP::size_type nRemoved;
    nRemoved = this->m_PseudoVertexMap.erase(pClosedSWEdge);
    assert(1 == nRemoved);
}

HRESULT SWPart::setTransforms(IComponent2 *pSWComponent)
{
    HRESULT hr = NOERROR;
    VARIANT vTransform;

    hr = pSWComponent->GetXform(&vTransform);

    double pdTransform[16];
    hr = ExtractVARIANTArrayData(vTransform, 16, pdTransform);
    hr = VariantClear(&vTransform);

    m_PartToAsm.set_to_identity();

    // copy the transform from SolidWorks - first 9 doubles are a 3X3 matrix
    int ii,jj;
    for (ii=0; ii<3; ii++)
    {
        for (jj=0; jj<3; jj++)
        {
            m_PartToAsm.set( ii, jj, pdTransform[3*jj+ii] );
        }
    }

    // next 3 doubles are the transformation
    m_PartToAsm.translate(pdTransform[9], pdTransform[10], pdTransform[11]);

    // next double is the scale factor
    assert(1.0 == pdTransform[12]);

    // invert the transformation to get the assembly to part transformation
    CubitTransformMatrix matCopy = m_PartToAsm;// copy
    m_AsmToPart = matCopy.inverse(); // invert

    return hr;
}

CubitTransformMatrix const&
SWPart::getPartToAssemblyTransform()
{
    return m_PartToAsm;
}

CubitTransformMatrix const&
SWPart::getAssemblyToPartTransform()
{
    return m_AsmToPart;
}

void
SWPart::transformPartToAssembly(CubitVector &partVect, CubitVector &assemblyVect)
{
    assemblyVect = m_PartToAsm * partVect;
}

void
SWPart::transformAssemblyToPart(CubitVector &assemblyVect, CubitVector &partVect)
{
    partVect = m_AsmToPart * assemblyVect;
}

void
SWPart::append_simple_attribute_virt(IUnknown *pSWEntity, CubitSimpleAttrib* attrib_ptr)
{
	// First, remove any like-typed attributes
	remove_simple_attribute_virt(pSWEntity, attrib_ptr);


    // copy the attribute passed in
    CubitSimpleAttrib *pSimpleAttrib = new CubitSimpleAttrib();
    *pSimpleAttrib = *attrib_ptr;

    m_AttribMap.insert(ATTRIBMAP::value_type(pSWEntity, pSimpleAttrib));

    LPDMREFKEY pRefKey = NULL;
    HRESULT hr = pSWEntity->QueryInterface(IID_IDMReferenceKey, (LPVOID*)&pRefKey);
    if (pRefKey)
    {
        unsigned long nByte;
        hr = pRefKey->GetKeySize(&nByte);
        BYTE* pData = new BYTE[nByte];

        hr = pRefKey->GetKey(nByte, pData);
        SWRefKey refKey;
        refKey.putData(nByte, pData);
        delete [] pData;

        m_KeyMap.insert(KEYMAP::value_type(pSWEntity, refKey));
        pRefKey->Release();
    }

    return;
}

void
SWPart::remove_simple_attribute_virt(IUnknown *pSWEntity, CubitSimpleAttrib* attrib_ptr)
{
    CubitSimpleAttrib *pAttrib;
    CubitString typeString;

    if (0 == m_AttribMap.size()) // no attributes in map
        return;

    // get the attributes associated with the SW entity passed in
    m_AttribMapIter = m_AttribMap.lower_bound(pSWEntity);

    if (m_AttribMapIter == m_AttribMap.end()) // no attributes on this entity
        return;

    // get the attribute type string
    typeString = attrib_ptr->string_data_list()->get()->c_str();

    // remove any attributes of the same type - should be at most one, since adding an
    // attribute removes any other like typed attributes before adding.
    int nRemoved = 0;
    while (m_AttribMapIter != m_AttribMap.upper_bound(pSWEntity))
    {
        pAttrib = (*m_AttribMapIter).second;
        if (0 == strcmp(typeString.c_str(), pAttrib->string_data_list()->get()->c_str()))
        {
            m_AttribMapIter = m_AttribMap.erase(m_AttribMapIter); // returns iterator to next or end
            delete pAttrib;
            nRemoved++;
        }
        else
        {
            // move to next in map
            ++m_AttribMapIter;
        }
    }
    // there should never be more than one of the same type in the list
    // since we check when adding
    assert(1 >= nRemoved);
    return;
}

void
SWPart::remove_all_simple_attribute_virt(IUnknown *pSWEntity)
{
    CubitSimpleAttrib *pAttrib = NULL;

    if (0 == m_AttribMap.size()) // no attributes in map
        return;

    // get the attributes associated with the SW entity passed in
    m_AttribMapIter = m_AttribMap.lower_bound(pSWEntity);

    if (m_AttribMapIter == m_AttribMap.end()) // no attributes on this entity
        return;

    // remove all attributes associated with the SW entity
    int nRemoved = 0;
    while (m_AttribMapIter != m_AttribMap.upper_bound(pSWEntity))
    {
        pAttrib = (*m_AttribMapIter).second;
        m_AttribMapIter = m_AttribMap.erase(m_AttribMapIter); // returns iterator to next or end
        delete pAttrib;
        nRemoved++;
    }

    return;
}

CubitStatus
SWPart::get_simple_attribute(IUnknown *pSWEntity,
                             DLIList<CubitSimpleAttrib*> &cubit_simple_attrib_list)
{
    DLIList<CubitSimpleAttrib*> simpleAttribList;

    m_AttribMapIter = m_AttribMap.lower_bound(pSWEntity);

    if (m_AttribMapIter == m_AttribMap.end()) // no attributes for this entity
        return CUBIT_SUCCESS;

    while (m_AttribMapIter != m_AttribMap.upper_bound(pSWEntity))
    {
        CubitSimpleAttrib *csa = (*m_AttribMapIter).second;
        CubitSimpleAttrib *copy = new CubitSimpleAttrib(csa);
        simpleAttribList.append(copy);
        m_AttribMapIter++;
    }

    // return the list
    cubit_simple_attrib_list = simpleAttribList;

	return CUBIT_SUCCESS;
}

CubitStatus
SWPart::get_simple_attribute(IUnknown *pSWEntity, const CubitString& name,
                             DLIList<CubitSimpleAttrib*> &cubit_simple_attrib_list)
{
    DLIList<CubitSimpleAttrib*> simpleAttribList;

    m_AttribMapIter = m_AttribMap.lower_bound(pSWEntity);

    if (m_AttribMapIter == m_AttribMap.end()) // no attributes for this entity
        return CUBIT_SUCCESS;

    while (m_AttribMapIter != m_AttribMap.upper_bound(pSWEntity))
    {
        CubitSimpleAttrib *csa = (*m_AttribMapIter).second;
        if( csa->character_type() == name )
        {
          CubitSimpleAttrib *copy = new CubitSimpleAttrib(csa);
          simpleAttribList.append(copy);
        }
        m_AttribMapIter++;
    }

    // return the list
    cubit_simple_attrib_list = simpleAttribList;

	return CUBIT_SUCCESS;
}

void
SWPart::append_simple_attribute_virt(TopologyBridge* pBridge, CubitSimpleAttrib* attrib_ptr)
{
    std::vector<CubitSimpleAttrib*> *pvect = 0;
    BodySW *pBodySW = CAST_TO(pBridge, BodySW);
    if (pBodySW)
    {
        pvect = &m_BodyAttribs;
    }
    else
    {
        LumpSW *pLumpSW = CAST_TO(pBridge, LumpSW);
        if (pLumpSW)
        {
            pvect = &m_LumpAttribs;
        }
    }
    if (0 == pvect)
    {
        assert(false); // invalid case
        return;
    }

	// First, remove any like-typed attributes
	remove_simple_attribute_virt(pBridge, attrib_ptr);

    // copy the attribute passed in
    CubitSimpleAttrib *pSimpleAttrib = new CubitSimpleAttrib();
    *pSimpleAttrib = *attrib_ptr;

    pvect->push_back(pSimpleAttrib);

    return;
}

void
SWPart::remove_simple_attribute_virt(TopologyBridge* pBridge, CubitSimpleAttrib* attrib_ptr)
{
    std::vector<CubitSimpleAttrib*> *pvect = 0;
    BodySW *pBodySW = CAST_TO(pBridge, BodySW);
    if (pBodySW)
    {
        pvect = &m_BodyAttribs;
    }
    else
    {
        LumpSW *pLumpSW = CAST_TO(pBridge, LumpSW);
        if (pLumpSW)
        {
            pvect = &m_LumpAttribs;
        }
    }
    if (0 == pvect)
    {
        assert(false); // invalid case
        return;
    }

    CubitSimpleAttrib *pAttrib;
    CubitString typeString;

    if (pvect->empty()) // no attributes in map
        return;

    // get the attributes associated with the SW entity passed in
    std::vector<CubitSimpleAttrib*>::iterator it;
    it = pvect->begin();

    // get the attribute type string
    typeString = attrib_ptr->string_data_list()->get()->c_str();

    // remove any attributes of the same type - should be at most one, since adding an
    // attribute removes any other like typed attributes before adding.
    int nRemoved = 0;
    while (it != pvect->end())
    {
        pAttrib = *it;
        if (0 == strcmp(typeString.c_str(), pAttrib->string_data_list()->get()->c_str()))
        {
            it = pvect->erase(it); // returns iterator to next or end
            delete pAttrib;
            nRemoved++;
        }
        else
        {
            // move to next in map
            ++it;
        }
    }
    // there should never be more than one of the same type in the list
    // since we check when adding
    assert(1 >= nRemoved);
    return;
}

void
SWPart::remove_all_simple_attribute_virt(TopologyBridge* pBridge)
{
    std::vector<CubitSimpleAttrib*> *pvect = 0;
    BodySW *pBodySW = CAST_TO(pBridge, BodySW);
    if (pBodySW)
    {
        pvect = &m_BodyAttribs;
    }
    else
    {
        LumpSW *pLumpSW = CAST_TO(pBridge, LumpSW);
        if (pLumpSW)
        {
            pvect = &m_LumpAttribs;
        }
    }
    if (0 == pvect)
    {
        assert(false); // invalid case
        return;
    }

    if (pvect->empty()) // no attributes in map
        return;

    CubitSimpleAttrib *pAttrib = NULL;

    // get the attributes in the list
    std::vector<CubitSimpleAttrib*>::iterator it;
    it = pvect->begin();

    // remove all attributes associated with the SW entity
    int nRemoved = 0;
    while (it != pvect->end())
    {
        pAttrib = *it;
        it = pvect->erase(it); // returns iterator to next or end
        delete pAttrib;
        nRemoved++;
    }

    return;
}

CubitStatus
SWPart::get_simple_attribute(TopologyBridge* pBridge,
                             DLIList<CubitSimpleAttrib*> &cubit_simple_attrib_list)
{
    std::vector<CubitSimpleAttrib*> *pvect = 0;
    BodySW *pBodySW = CAST_TO(pBridge, BodySW);
    if (pBodySW)
    {
        pvect = &m_BodyAttribs;
    }
    else
    {
        LumpSW *pLumpSW = CAST_TO(pBridge, LumpSW);
        if (pLumpSW)
        {
            pvect = &m_LumpAttribs;
        }
    }
    if (0 == pvect)
    {
        assert(false); // invalid case
        return CUBIT_FAILURE;
    }

    if (pvect->empty()) // no attributes in map
        return CUBIT_SUCCESS;

    DLIList<CubitSimpleAttrib*> simpleAttribList;

    std::vector<CubitSimpleAttrib*>::iterator it;
    it = pvect->begin();

    while (it != pvect->end())
    {
        CubitSimpleAttrib *csa = (*it);
        CubitSimpleAttrib *copy = new CubitSimpleAttrib(csa);
        simpleAttribList.append(copy);
        it++;
    }

    // return the list
    cubit_simple_attrib_list = simpleAttribList;

	return CUBIT_SUCCESS;
}

CubitStatus
SWPart::get_simple_attribute(TopologyBridge* pBridge, const CubitString& name,
                             DLIList<CubitSimpleAttrib*> &cubit_simple_attrib_list)
{
    std::vector<CubitSimpleAttrib*> *pvect = 0;
    BodySW *pBodySW = CAST_TO(pBridge, BodySW);
    if (pBodySW)
    {
        pvect = &m_BodyAttribs;
    }
    else
    {
        LumpSW *pLumpSW = CAST_TO(pBridge, LumpSW);
        if (pLumpSW)
        {
            pvect = &m_LumpAttribs;
        }
    }
    if (0 == pvect)
    {
        assert(false); // invalid case
        return CUBIT_FAILURE;
    }

    if (pvect->empty()) // no attributes in map
        return CUBIT_SUCCESS;

    DLIList<CubitSimpleAttrib*> simpleAttribList;

    std::vector<CubitSimpleAttrib*>::iterator it;
    it = pvect->begin();

    while (it != pvect->end())
    {
        CubitSimpleAttrib *csa = (*it);
        if( csa->character_type() == name )
        {
          simpleAttribList.append( new CubitSimpleAttrib(csa) );
        }
        it++;
    }

    // return the list
    cubit_simple_attrib_list = simpleAttribList;

	return CUBIT_SUCCESS;
}

void
SWPart::updateAttribs()
{
    // Body and lump attributes don't need to be updated.
    // Since there is only one body and one lump these attributes
    // are always associated to the existing body and lump, respectively

    HRESULT hr = NOERROR;
    LPDMREF pRef = NULL;

    // traverse the old maps and fill new maps
    // then we can clear the old maps and replace them with the new maps
    ATTRIBMAP newAttribMap;
    ATTRIBMAP newDanglingMap;
    KEYMAP newKeyMap;


// update other attributes


    // if there are no attributes just return
    if (m_AttribMap.empty() && m_DanglingAttribMap.empty())
        return;

    hr = m_pImport->document()->QueryInterface(IID_IDMReference, (LPVOID*)&pRef);
    if (!pRef)
    {
        assert(false); // this interface should be supported

        // make all attributes dangling

        // combine the old AttribMap and old DanglingAttribMap to get the new Dangling map
//        newDanglingAttribMap = m_AttribMap;
//        m_AttribMap.clear();
    }

    rebindKeys(pRef, m_AttribMap, m_KeyMap, newAttribMap, newKeyMap, newDanglingMap);
    rebindKeys(pRef, m_DanglingAttribMap, m_KeyMap, newAttribMap, newKeyMap, newDanglingMap);


    // TODO - report attributes that did not bind

    // store the new maps
    m_DanglingAttribMap.clear();
    m_DanglingAttribMap = newDanglingMap;

    m_AttribMap.clear();
    m_AttribMap = newAttribMap;

    m_KeyMap.clear();
    m_KeyMap = newKeyMap;


    pRef->Release();

	return;
}

void
SWPart::setImport(SWImport *pImport)
{
    m_pImport = pImport;
}

static void rebindKeys(LPDMREF pRef, ATTRIBMAP &oldAttribMap, KEYMAP &oldKeyMap,
                       ATTRIBMAP &newAttribMap, KEYMAP &newKeyMap, ATTRIBMAP &newDanglingMap)
{
    HRESULT hr = NOERROR;
    IUnknown *pOldUnk = NULL;
    IUnknown *pNewUnk = NULL;
    LPDMREFKEY pRefKey = NULL;
    IFace2 *pSWFace = NULL;
    IEdge *pSWEdge = NULL;
    IVertex *pSWVertex = NULL;
    IID objIID;

    SWRefKey key;
    CubitSimpleAttrib* pCSA;


    ATTRIBMAP::iterator AttribIter;
    KEYMAP::iterator KeyIter;

    AttribIter = oldAttribMap.begin(); // first attribute
    while (AttribIter != oldAttribMap.end())
    {
        // get the entity pointer (no longer a valid pointer!) and the attribute pointer
        pOldUnk = (*AttribIter).first;
        pCSA = (*AttribIter).second;

        // TODO - check for multiple map entries for a given key
        // so we only bind each key once.

        // look up the key corresponding to the entity pointer in the key map
        KeyIter = oldKeyMap.find(pOldUnk);
        assert(oldKeyMap.end() != KeyIter); // key should always exist in the map
        key = (*KeyIter).second;

        // use the key to get the new entity pointer
        hr = pRef->BindKeyToInterface(IID_IUnknown, key.size(),
                                      key.keyData(), (LPVOID*)&pNewUnk);

        if (pNewUnk)
        {
            hr = pNewUnk->QueryInterface(IID_IDMReferenceKey, (LPVOID*)&pRefKey);
            pNewUnk->Release();
        }


        // if the key doesn't bind, this attribute is dangling
        if (!pRefKey)
        {
            // need to make this attribute dangling
            newDanglingMap.insert(ATTRIBMAP::value_type(pOldUnk, pCSA));

            // hang on to the key, since it could be resolved in the future
            newKeyMap.insert(KEYMAP::value_type(pOldUnk, key));

        }
        else
        {
            // insert the new entity pointer and attribute in the new attrib map
            // insert the new entity pointer and key in the new key map

            hr = pRefKey->GetObjectType(&objIID);
            if(IID_IDMFace == objIID)
            {
                hr = pRefKey->QueryInterface(IID_IFace2, (LPVOID*)&pSWFace);
                newAttribMap.insert(ATTRIBMAP::value_type(pSWFace, pCSA));
                newKeyMap.insert(KEYMAP::value_type(pSWFace, key));
                pSWFace->Release(); // reference is held by CGM entity
            }
            else if(IID_IDMEdge == objIID)
            {
                hr = pRefKey->QueryInterface(IID_IEdge, (LPVOID*)&pSWEdge);
                newAttribMap.insert(ATTRIBMAP::value_type(pSWEdge, pCSA));
                newKeyMap.insert(KEYMAP::value_type(pSWEdge, key));
                pSWEdge->Release(); // reference is held by CGM entity
            }
            else if(IID_IDMVertex == objIID)
            {
                hr = pRefKey->QueryInterface(IID_IVertex, (LPVOID*)&pSWVertex);
                newAttribMap.insert(ATTRIBMAP::value_type(pSWVertex, pCSA));
                newKeyMap.insert(KEYMAP::value_type(pSWVertex, key));
                pSWVertex->Release(); // reference is held by CGM entity
            }
            else
            {
                assert(false);
            }
            pRefKey->Release();
        }

        AttribIter++; // next attribute
    }

}

void
SWPart::setName(IComponent2 *pSWComp)
{
    USES_CONVERSION;

    HRESULT hr = NOERROR;

    BSTR bsName = ::SysAllocString(A2OLE(""));
    hr = pSWComp->get_Name2(&bsName);
    char *pszName = OLE2A(bsName);

    m_csPartName = pszName;
    ::SysFreeString(bsName);
}

CubitString
SWPart::getName()
{
    return m_csPartName;
}
