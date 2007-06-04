#include "stdafx.h"

#include "SWImport.hpp"
#include "SWPart.hpp"

#include "GeometryQueryTool.hpp"

#include "BodySM.hpp"

#include "BodySW.hpp"
#include "LumpSW.hpp"
#include "ShellSW.hpp"
#include "SurfaceSW.hpp"
#include "LoopSW.hpp"
#include "CoEdgeSW.hpp"
#include "CurveSW.hpp"
#include "PointSW.hpp"

#include "CastTo.hpp"


static HRESULT GetSWLeafComponents(IComponent2 *pRootComponent, DLIList<IComponent2 *> &compList);

extern HRESULT ExtractVARIANTArrayData(VARIANT &v, long lArraySize, LPDISPATCH *pDispArray);


SWImport::SWImport()
{
    m_SWDocType = swDocNONE;
    m_pSWModelDoc = NULL;
}

SWImport::~SWImport()
{
    SWPart* pSWPart;
    IComponent2 *pSWComp = NULL;
    int i;

    m_SWPartList.reset();
    for (i=0; i<m_SWPartList.size(); i++)
    {
        pSWPart = m_SWPartList.get_and_step();
        delete pSWPart;
    }

    m_SWComponentList.reset();
    for (i=0; i<m_SWComponentList.size(); i++)
    {
        pSWComp = m_SWComponentList.get_and_step();
        pSWComp->Release();
    }

    if (m_pSWModelDoc)
        m_pSWModelDoc->Release();
    m_pSWModelDoc = NULL;
}

CubitStatus SWImport::setDocument(IModelDoc2 *pSWDoc)
{
    HRESULT hr = NOERROR;
    assert(pSWDoc);

    // check the document type - only handle part and assembly documents
    long lDocType;
    hr = pSWDoc->GetType(&lDocType);
    if ((swDocPART == lDocType) || (swDocASSEMBLY == lDocType))
    {
        pSWDoc->AddRef();
        m_pSWModelDoc = pSWDoc;

        m_SWDocType = (enum swDocumentTypes_e) lDocType;

        return CUBIT_SUCCESS;
    }
    else
    {
        return CUBIT_FAILURE;
    }
}

IModelDoc2 *SWImport::document()
{
    return m_pSWModelDoc;
}

void SWImport::addPart(SWPart* &pSWPart)
{
    pSWPart = NULL;
    assert(m_pSWModelDoc);            // make sure document already set
    assert(swDocPART == m_SWDocType); // this method is only for part documents
    assert(0 == m_SWPartList.size()); // this function should only be called once

    // create a new SWPart
    SWPart *pNewSWPart = new SWPart();
    if (NULL == pNewSWPart)
        return;

    // add the part to the list
    m_SWPartList.append(pNewSWPart);

    pNewSWPart->setImport(this);

    // return the part
    pSWPart = pNewSWPart;
    return;
}

void SWImport::addPart(IComponent2 *pSWComp, SWPart* &pSWPart)
{
    pSWPart = NULL;
    assert(pSWComp);
    assert(m_pSWModelDoc);                // make sure document already set
    assert(swDocASSEMBLY == m_SWDocType); // this method is only for assembly documents

    // create a new SWPart
    SWPart *pNewSWPart = new SWPart();
    if (NULL == pNewSWPart)
        return;

    pNewSWPart->setImport(this);

    updateComponentAndPartLists(pSWComp, pNewSWPart);

    // return the part
    pSWPart = pNewSWPart;
    return;
}

void SWImport::updateComponentAndPartLists(IComponent2 *pSWComp, SWPart* &pNewSWPart)
{
    // add the part to the list
    m_SWPartList.append(pNewSWPart);

    // add the component to the list
    pSWComp->AddRef();
    m_SWComponentList.append(pSWComp);

    // part and component lists should be the same size
    assert(m_SWPartList.size() == m_SWComponentList.size());

    return;
}

CubitStatus 
SWImport::import_solid_model(IModelDoc2 *pSWModelDoc,
                             DLIList<TopologyBridge*>& imported_entities)
{
//    USES_CONVERSION;

#ifdef BOYD17 
    IModelDoc2 *pSWDoc = NULL;
#endif
    IBody2 *pSWBody = NULL;
    HRESULT hr;
    long lDocType = swDocNONE;
    int i;
    SWPart *pSWPart = NULL;
    int bUpdating = FALSE;

    // see if we are importing a new doc or updating one already imported
    if (m_pSWModelDoc)
    {
        // only allowing import of a single document
        // TODO - do we need to check doc name?
        if (pSWModelDoc != m_pSWModelDoc)
            return CUBIT_FAILURE;

        bUpdating = TRUE;

        // delete the existing geometry
        GeometryQueryTool::instance()->delete_geometry();
    }
    else
    {
        setDocument(pSWModelDoc);
    }

    // check for and handle assembly and part documents
    hr = pSWModelDoc->GetType(&lDocType);
    assert(SUCCEEDED(hr));

    switch (lDocType)
    {

    case swDocPART:
        {
            IPartDoc *pSWPartDoc = NULL;
            hr = pSWModelDoc->QueryInterface(IID_IPartDoc, (LPVOID*)&pSWPartDoc);
            pSWModelDoc->Release();
            pSWModelDoc = NULL;

            if (FAILED(hr))
            {
                return CUBIT_FAILURE;
            }

            // get the body pointer
            hr = pSWPartDoc->IBodyObject2(&pSWBody);
            pSWPartDoc->Release();
            if (NULL == pSWBody)
            {
                return CUBIT_FAILURE;
            }

            // use the existing part if updating, otherwise create a new part
            if (bUpdating)
            {
                assert(1 == m_SWPartList.size()); // only one part in a part doc
                m_SWPartList.reset();
                pSWPart = m_SWPartList.get();
            }
            else
            {
                addPart(pSWPart);
            }

            if (pSWPart)
            {
                // Build a VGI Body (i.e, the entire VGI structure) from this 
                // SolidWorks part
                BodySM *this_bodysm = populate_topology_bridges(pSWPart, pSWBody);
                assert(this_bodysm);
                imported_entities.append(this_bodysm);

                // update the attributes in the part
                pSWPart->updateAttribs();
            }

            pSWBody->Release();
        }
        break;

    case swDocASSEMBLY:
        {
            IConfiguration *pConfig = NULL;
            IComponent2 *pRootComponent = NULL;
            IComponent2 *pLeafComponent = NULL;

            DLIList<IComponent2 *> compList;

            hr = pSWModelDoc->IGetActiveConfiguration(&pConfig);
            if (SUCCEEDED(hr) && pConfig)
            {
                hr = pConfig->IGetRootComponent2(&pRootComponent);
                if (SUCCEEDED(hr) && pRootComponent)
                {
                    // save old part list
                    DLIList<SWPart*> oldPartList;
                    oldPartList = m_SWPartList;

                    m_SWPartList.clean_out();
                    m_SWComponentList.clean_out();

                    // recursively traverse the assembly and get the leaf components (parts)
                    hr = GetSWLeafComponents(pRootComponent, compList);

                    int nComp = 0;
                    compList.reset();
                	for (i=0; i<compList.size(); i++)
                    {
                        pLeafComponent = compList.get_and_step();
                        nComp++;
                        hr = pLeafComponent->IGetBody(&pSWBody);
                        if (pSWBody)
                        {
                            // find an existing part that matches this component
                            pSWPart = NULL;
                            if (bUpdating)
                                findComponentPart(pLeafComponent, oldPartList, pSWPart);

                             // update the lists or add a new part
                            if (pSWPart)
                            {
                                oldPartList.remove(pSWPart);
                                updateComponentAndPartLists(pLeafComponent, pSWPart);
                            }
                            else
                            {
                                addPart(pLeafComponent, pSWPart);
                            }

                            if (pSWPart)
                            {

                                // set the part transformation
                                hr = pSWPart->setTransforms(pLeafComponent);
                                pSWPart->setName(pLeafComponent);

                                // Build a VGI Body (i.e, the entire VGI structure) from this 
                                // SolidWorks part
                                BodySM *this_bodysm = populate_topology_bridges(pSWPart, pSWBody);
                                assert(this_bodysm);
                                imported_entities.append(this_bodysm);

                                // update the attributes in the part
                                pSWPart->updateAttribs();
                            }
                            pSWBody->Release();
                        }
                        pLeafComponent->Release();

                    }
                    // keep pointers to parts that may have gone away in case they come back
                    // in a future import
                    m_SWPartList += oldPartList;
                }
            }
        }
        break;

    default:
        assert(false);
        return CUBIT_FAILURE;
    }

    return CUBIT_SUCCESS;
}

BodySM * SWImport::populate_topology_bridges(SWPart *pSWPart, IBody2 *pSWBody)
{
    USES_CONVERSION;
    assert ( pSWBody != NULL);

    BodySW* pBodySW = NULL;
    LumpSW* pLumpSW = NULL;
    ShellSW* pShellSW = NULL;
    IFace2 *pSWFace = NULL;
    HRESULT hr = NOERROR;

// TODO - update body/bodies when reading in geometry

    // First check to make sure we haven't already created a BodySW
    // from this SolidWorks Body.
//  LPENTITY pEntity = NULL;
//  hr = pSWBody->QueryInterface(IID_IEntity, (LPVOID*)&pEntity);
//  if (FAILED(hr))
//      return NULL;

//  SWBridge* pBridge = ATTRIB_CUBIT_OWNER::cubit_owner(pEntity);
//  pEntity->Release();
//  pBodySW = CAST_TO(pBridge, BodySW);

    if (pBodySW == NULL)
    {
        pBodySW = new BodySW(pSWPart);
        pLumpSW = new LumpSW(pSWPart);
        pShellSW = new ShellSW(pSWPart);
    }
    else
        assert(FALSE);

    assert(NULL != pBodySW);
    assert(NULL != pLumpSW);
    assert(NULL != pShellSW);

    pBodySW->set_BODY_ptr(pSWBody);
    pLumpSW->set_BODY_ptr(pSWBody);
    pShellSW->set_BODY_ptr(pSWBody);

    IEnumFaces2 *pFaceEnum = NULL;
    hr = pSWBody->EnumFaces( &pFaceEnum);

    pSWFace = NULL;
    long nFetched = 0;
    long nFace = 0;
    while (S_OK == pFaceEnum->Next(1, &pSWFace, &nFetched))
    {
        assert(pSWFace);
        nFace++;

        populate_topology_bridges(pSWPart, pSWFace);
        pSWFace->Release();
        pSWFace = NULL;
    }
    pFaceEnum->Release();

    return pBodySW;
}

Surface * SWImport::populate_topology_bridges(SWPart *pSWPart, IFace2 *pSWFace)
{
    USES_CONVERSION;
    DIAGNOSTIC("processing SW Face 0x%x\n", pSWFace);

    HRESULT hr = NOERROR;
    assert ( pSWFace != NULL);

    SurfaceSW* pSurfaceSW = NULL;

    // First check to make sure we haven't already created a SurfaceSW
    // from this SW FACE.
    TopologyBridge* pBridge = pSWPart->cubit_owner(pSWFace);
    pSurfaceSW = CAST_TO(pBridge, SurfaceSW);
    if (pSurfaceSW == NULL)
    {
//        LPENTITY pEnt;
//        LPDMREFKEY pKey;
//        unsigned long lKeySize;
//        hr = pSWFace->QueryInterface(IID_IEntity, (LPVOID*)&pEnt);
//        if (pEnt)
//        {
//            char attname[32];
//            sprintf(attname, "CubitFace%d", S_nAttrib++);
//            BSTR bsAttName = ::SysAllocString(A2OLE(attname));
//            LPMODELDOC pModelDoc = NULL;
//            LPATTRIBUTE pAttr = NULL;
//            hr = m_pSWApp->get_IActiveDoc(&pModelDoc);
//            hr = m_pAttDef->ICreateInstance2(pModelDoc, pEnt, bsAttName, 0, &pAttr);
//
//            ::SysFreeString(bsAttName);
//            pModelDoc->Release();
//            pEnt->Release();
//        }
//        hr = pSWFace->QueryInterface(IID_IDMReferenceKey, (LPVOID*)&pKey);
//        if (pKey)
//        {
//            pKey->GetKeySize(&lKeySize);
//            pKey->Release();
//        }

        assert(pBridge == NULL);
        pSurfaceSW = new SurfaceSW(pSWPart) ;
    }

    assert(pSurfaceSW != NULL);
    pSurfaceSW->set_FACE_ptr(pSWFace);


    // now loops
    IEnumLoops2 *pLoopEnum = NULL;
    hr = pSWFace->EnumLoops(&pLoopEnum);

    ILoop2 *pSWLoop = NULL;
    ICoEdge *pSWCoedge = NULL;
    IEdge *pSWEdge = NULL;
    long nFetched = 0;
    int nLoop = 0;
    while (S_OK == pLoopEnum->Next(1, &pSWLoop, &nFetched))
    {
        nLoop++;

// This looks like a SW bug -- immediately after getting the loop from the face
// if I get the face from the loop I get a different object.  The only way to tell
// they are the same is using pFace->IsSame.  This is only happening in assemblies
// where a part is inserted more than once.  Byron 10/4/01
//        LPFACE pTestFace = NULL;
//        hr = pSWLoop->IGetFace(&pTestFace);
//        assert(pSWFace == pTestFace);
//        pTestFace->Release();

//        LPENTITY pEnt;
//        LPDMREFKEY pKey;
//        unsigned long lKeySize;
//        hr = pSWLoop->QueryInterface(IID_IEntity, (LPVOID*)&pEnt);
//        if (pEnt)
//        {
//            char attname[32];
//            sprintf(attname, "CubitLoop%d", S_nAttrib++);
//            BSTR bsAttName = ::SysAllocString(A2OLE(attname));
//            LPMODELDOC pModelDoc = NULL;
//            LPATTRIBUTE pAttr = NULL;
//            hr = m_pSWApp->get_IActiveDoc(&pModelDoc);
//            hr = m_pAttDef->ICreateInstance2(pModelDoc, pEnt, bsAttName, 0, &pAttr);
//
//            ::SysFreeString(bsAttName);
//            pModelDoc->Release();
//            pEnt->Release();
//        }
//        hr = pSWLoop->QueryInterface(IID_IDMReferenceKey, (LPVOID*)&pKey);
//        if (pKey)
//        {
//            hr = pKey->GetKeySize(&lKeySize);
//            pKey->Release();
//        }

        // First check to make sure we haven't already created a LoopSW
        // from this SW LOOP.
        pBridge = pSWPart->cubit_owner(pSWLoop);
        LoopSW *pLoopSW = CAST_TO(pBridge, LoopSW);
        if (pLoopSW == NULL)
        {
            assert(pBridge == NULL);
            pLoopSW = new LoopSW(pSWPart) ;
        }

        pLoopSW->set_LOOP_ptr(pSWLoop);

        // now coedges
        IEnumCoEdges *pCoEdgeEnum = NULL;
        hr = pSWLoop->EnumCoEdges(&pCoEdgeEnum);

        nFetched = 0;
        int nCoedge = 0;
        while (S_OK == pCoEdgeEnum->Next(1, &pSWCoedge, &nFetched))
        {
            nCoedge++;
            assert(pSWCoedge);

//        hr = pSWCoedge->QueryInterface(IID_IEntity, (LPVOID*)&pEnt);
//        if (pEnt)
//        {
//            char attname[32];
//            sprintf(attname, "CubitCoEdge%d", S_nAttrib++);
//            BSTR bsAttName = ::SysAllocString(A2OLE(attname));
//            LPMODELDOC pModelDoc = NULL;
//            LPATTRIBUTE pAttr = NULL;
//            hr = m_pSWApp->get_IActiveDoc(&pModelDoc);
//            hr = m_pAttDef->ICreateInstance2(pModelDoc, pEnt, bsAttName, 0, &pAttr);
//
//            ::SysFreeString(bsAttName);
//            pModelDoc->Release();
//            pEnt->Release();
//        }
//        hr = pSWCoedge->QueryInterface(IID_IDMReferenceKey, (LPVOID*)&pKey);
//        if (pKey)
//        {
//            hr = pKey->GetKeySize(&lKeySize);
//            pKey->Release();
//        }

            // First check to make sure we haven't already created a CoEdgeSW
            // from this SW CoEdge.
            pBridge = pSWPart->cubit_owner(pSWCoedge);
            CoEdgeSW *pCoEdgeSW = CAST_TO(pBridge, CoEdgeSW);
            if (pCoEdgeSW == NULL)
            {
                assert(pBridge == NULL);
                pCoEdgeSW = new CoEdgeSW(pSWPart) ;
            }

            pCoEdgeSW->set_COEDGE_ptr(pSWCoedge);

            pSWCoedge->IGetEdge(&pSWEdge);
            pSWCoedge->Release();

            populate_topology_bridges(pSWPart, pSWEdge);
            pSWEdge->Release();
        }
        if (0 >= nCoedge)
        {
            // loop with no edges found - delete this loop
            // this appears to be a problem in an old SW file in the samples directory???
            delete pLoopSW;
            assert(false);
        }
        pCoEdgeEnum->Release();
        pSWLoop->Release();
    }
    pLoopEnum->Release();

    return pSurfaceSW;
}

Curve * SWImport::populate_topology_bridges(SWPart *pSWPart, IEdge *pSWEdge)
{
    USES_CONVERSION;
    DIAGNOSTIC("processing SW Edge 0x%x\n", pSWEdge);
    assert (pSWEdge != NULL);

    HRESULT hr = NOERROR;
    CurveSW* pCurveSW = NULL;

    // First check to make sure we haven't already created a CurveSW
    // from this SW EDGE.
    TopologyBridge* pBridge = pSWPart->cubit_owner(pSWEdge);
    pCurveSW = CAST_TO(pBridge, CurveSW);
    if (pCurveSW == NULL)
    {
//        LPENTITY pEnt;
//        LPDMREFKEY pKey;
//        unsigned long lKeySize;
//        hr = pSWEdge->QueryInterface(IID_IEntity, (LPVOID*)&pEnt);
//        if (pEnt)
//        {
//            char attname[32];
//            sprintf(attname, "CubitEdge%d", S_nAttrib++);
//            BSTR bsAttName = ::SysAllocString(A2OLE(attname));
//            LPMODELDOC pModelDoc = NULL;
//            LPATTRIBUTE pAttr = NULL;
//            hr = m_pSWApp->get_IActiveDoc(&pModelDoc);
//            hr = m_pAttDef->ICreateInstance2(pModelDoc, pEnt, bsAttName, 0, &pAttr);
//
//            ::SysFreeString(bsAttName);
//            pModelDoc->Release();
//            pEnt->Release();
//        }
//        hr = pSWEdge->QueryInterface(IID_IDMReferenceKey, (LPVOID*)&pKey);
//        if (pKey)
//        {
//            pKey->GetKeySize(&lKeySize);
//            pKey->Release();
//        }

        assert(pBridge == NULL);
        pCurveSW = new CurveSW(pSWPart);

        assert(pCurveSW != NULL);
        pCurveSW->set_EDGE_ptr(pSWEdge);

        //VARIANT_BOOL bReversed;
        //hr = pSWEdge->IsParamReversed(&bReversed);

        // get the start and end vertices
        IVertex *pStartVertex = NULL;
        hr = pSWEdge->IGetStartVertex(&pStartVertex);
        IVertex *pEndVertex = NULL;
        hr = pSWEdge->IGetEndVertex(&pEndVertex);

        // check for closed
        if (pStartVertex)
        {
            populate_topology_bridges(pSWPart, pStartVertex);
            pStartVertex->Release();

            assert(pEndVertex); // should not be NULL
            populate_topology_bridges(pSWPart, pEndVertex);
            pEndVertex->Release();
        }
        else // periodic edge
        {
            assert(NULL == pEndVertex);

            // this edge has no vertices
            // so create a PointSW for CGM
            PointSW *pPseudoPoint = new PointSW(pSWPart);
            assert(pPseudoPoint);

            pPseudoPoint->set_closed_edge_ptr(pSWEdge);
        }
    }


    return pCurveSW;
}

Point * SWImport::populate_topology_bridges(SWPart *pSWPart, IVertex *pSWVertex)
{
    USES_CONVERSION;
    DIAGNOSTIC("processing SW Vertex 0x%x\n", pSWVertex);
    assert (pSWVertex != NULL);

    PointSW* pPointSW = NULL;

    // First check to make sure we haven't already created a PointSW
    // from this SW VERTEX.
    TopologyBridge* pBridge = pSWPart->cubit_owner(pSWVertex);
    pPointSW = CAST_TO(pBridge, PointSW);
    if (pPointSW == NULL)
    {
//        LPENTITY pEnt;
//        LPDMREFKEY pKey;
//        unsigned long lKeySize;
//        HRESULT hr = pSWVertex->QueryInterface(IID_IEntity, (LPVOID*)&pEnt);
//        if (pEnt)
//        {
//            char attname[32];
//            sprintf(attname, "CubitVertex%d", S_nAttrib++);
//            BSTR bsAttName = ::SysAllocString(A2OLE(attname));
//            LPMODELDOC pModelDoc = NULL;
//            LPATTRIBUTE pAttr = NULL;
//            hr = m_pSWApp->get_IActiveDoc(&pModelDoc);
//            hr = m_pAttDef->ICreateInstance2(pModelDoc, pEnt, bsAttName, 0, &pAttr);
//
//            ::SysFreeString(bsAttName);
//            pModelDoc->Release();
//            pEnt->Release();
//        }
//        hr = pSWVertex->QueryInterface(IID_IDMReferenceKey, (LPVOID*)&pKey);
//        if (pKey)
//        {
//            pKey->GetKeySize(&lKeySize);
//            pKey->Release();
//        }

        assert(pBridge == NULL);
        pPointSW = new PointSW(pSWPart) ;
        pPointSW->set_VERTEX_ptr(pSWVertex);
    }

    assert(pPointSW != NULL);

    return pPointSW;
}

static HRESULT
GetSWLeafComponents(IComponent2 *pRootComponent, DLIList<IComponent2 *> &compList)
{
    HRESULT hr = NOERROR;
    VARIANT vChildren;
    IComponent2 *pComponent = NULL;

    int lNumChildren;
    hr = pRootComponent->IGetChildrenCount(&lNumChildren);

    if (lNumChildren)
    {
        LPDISPATCH *pDispChildren = new LPDISPATCH[lNumChildren];

        hr = pRootComponent->GetChildren(&vChildren);

        hr = ExtractVARIANTArrayData(vChildren, lNumChildren, pDispChildren);

        for (int i=0; i<lNumChildren; i++)
        {
            hr = pDispChildren[i]->QueryInterface(IID_IComponent2, (LPVOID*)&pComponent);
            if (SUCCEEDED(hr) && pComponent)
            {
                hr = GetSWLeafComponents(pComponent, compList);
                pComponent->Release();
            }
        }
        delete [] pDispChildren;
        hr = VariantClear(&vChildren);
    }
    else
    {
        // no children means this is a leaf component (part document)
        // add to the leaf component list
        compList.append(pRootComponent);
        pRootComponent->AddRef(); // ref count for list entry
    }

    return hr;
}

void SWImport::findComponentPart(IComponent2 *pLeafComponent, DLIList<SWPart*> &oldPartList, SWPart * &pSWPart)
{
    USES_CONVERSION;

    HRESULT hr = NOERROR;
    SWPart* pTest = 0;
    int bFound = 0;

    oldPartList.reset();

    // get the component name
    char* pszName;
    BSTR bsName = ::SysAllocString(A2OLE(""));
    hr = pLeafComponent->get_Name2(&bsName);
    pszName = OLE2A(bsName);

    CubitString csName;
    for (int i=0; i<oldPartList.size() && !bFound; i++)
    {
        pTest = oldPartList.get();
        csName = pTest->getName();
        if (csName == pszName )
        {
            bFound = TRUE;
            pSWPart = pTest;
        }
    }
    ::SysFreeString(bsName);
}
