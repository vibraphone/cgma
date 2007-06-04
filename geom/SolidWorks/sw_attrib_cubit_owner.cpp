//-------------------------------------------------------------------------
// Filename      : attrib_cubit_owner.cpp
//
// Purpose       : COM implementation.  This attribute represents a pointer 
//                 from an SW ENTITY to a Cubit SWBridge
//
// Creator       : Joel Kopp, Greg Neilson
//
// Owner         : Tim Tautges
//-------------------------------------------------------------------------

// ********** BEGIN SOLIDWORKS INCLUDES       **********
// Precompiled header
#include "stdafx.h"

//#include <amapp.h>
#import "C:\Program Files\SolidWorks\sldworks.tlb" raw_interfaces_only, raw_native_types, no_namespace, named_guids     //Import the SolidWorks type library
#include <swconst.h>
//#import "c:\program files\solidworks\sldworks.tlb"
// ********** END SOLIDWORKS INCLUDES         **********

// ********** BEGIN STANDARD INCLUDES         **********
#include <stdio.h>
#include <memory.h>
// ********** END STANDARD INCLUDES           **********

// ********** BEGIN CUBIT INCLUDES            **********
//#include "attrib_cubit_owner.hpp"
#include "SWBridge.hpp"
#include "TopologyEntity.hpp"
#include "TopologyBridge.hpp"
#include "CubitMessage.hpp"
#include "CubitAttribUser.hpp"
#include "CastTo.hpp"
#include "BodySW.hpp"
#include "LumpSW.hpp"
#include "ShellSW.hpp"
#include "SurfaceSW.hpp"
#include "LoopSW.hpp"
#include "CoEdgeSW.hpp"
#include "CurveSW.hpp"
#include "PointSW.hpp"

//#include "attrib_snl_simple.hpp"
#include "GeometryQueryTool.hpp"
#include "SWQueryEngine.hpp"
#include "DLIList.hpp"
// ********** END CUBIT INCLUDES              **********

// ********** BEGIN MACRO DEFINITIONS         **********
/*#define THIS() SW_ATTRIB_CUBIT_OWNER
#define THIS_LIB NONE
#define PARENT() ATTRIB_SNL
#define PARENT_LIB NONE

#define SW_ATTRIB_CUBIT_OWNER_NAME "cubit_owner"

ATTRIB_DEF( "cubit_owner_attribute")

SAVE_DEF
// Don't save

RESTORE_DEF
// Don't restore

COPY_DEF
   set_cubit_owner ( from->cubit_owner() );
//    PRINT_INFO("COPY_DEF called on attrib_cubit_owner.\n");

SCAN_DEF
// (no specific pointer data)

FIX_POINTER_DEF
// (no specific pointer data)

TERMINATE_DEF
// Don't do anything special*/
// ********** END MACRO DEFINITIONS           **********
int SW_ATTRIB_CUBIT_OWNER::ownerAttCount = 0;
IAttributeDef *SW_ATTRIB_CUBIT_OWNER::attDef = NULL;
CubitBoolean SW_ATTRIB_CUBIT_OWNER::initialized = CUBIT_FALSE;


  static DLIList<IFace2*> S_SWFaceList;
  static DLIList<ILoop2*> S_SWLoopList;
  static DLIList<ICoEdge*> S_SWCoEdgeList;
  static DLIList<IEdge*> S_SWEdgeList;
  static DLIList<IVertex*> S_SWVertexList;

  static SWBridge* S_body;
  static SWBridge* S_lump;
  static SWBridge* S_shell;
  static DLIList<SWBridge*> S_surfaceList;
  static DLIList<SWBridge*> S_loopList;
  static DLIList<SWBridge*> S_coedgeList;
  static DLIList<SWBridge*> S_curveList;
  static DLIList<SWBridge*> S_pointList;
  



// TODO - make a simple constructor and make a separate method to
// attach an attribute to the SolidWorks entity

// make a cubit_owner attribute
SW_ATTRIB_CUBIT_OWNER::SW_ATTRIB_CUBIT_OWNER(IUnknown *owner,
                                       SWBridge *cubit_owner)
    : /*ATTRIB_SNL(owner),*/ cubitOwnerData(NULL)
{
    assert(false);
    USES_CONVERSION;

	HRESULT hr;
	set_cubit_owner(cubit_owner);

    // used to name attribute (uses attribute number, since attrib names
	// must be unique.)
    char pszAttribName[128] = "CUBIT Entity ";
    char pszUnique[33];


    IEntity *pEntity = NULL;
    hr = owner->QueryInterface(IID_IEntity, (LPVOID*)&pEntity);

    // if QueryInterface for IEntity fails, make sure the owner is the body
    //  pEntity == NULL will attach an attribute to the document
    if (FAILED(hr))
    {
        IBody2 *pBody = NULL;
        hr = owner->QueryInterface(IID_IBody2, (LPVOID*)&pBody);
        if (FAILED(hr))
            return;
        assert(pBody);
        pBody->Release();

        //BodySW *pBodySW = CAST_TO(cubit_owner, BodySW);
        BodySW *pBodySW = dynamic_cast<BodySW*> (cubit_owner);
        LumpSW *pLumpSW = CAST_TO(cubit_owner, LumpSW);
        if (pBodySW)
            strcpy(pszUnique, "Body");
        else if (pLumpSW)
            strcpy(pszUnique, "Lump");
        else
            return; // TODO - return an error
    }
    else
    {
        ShellSW *pShellSW = CAST_TO(cubit_owner, ShellSW);
        if (pShellSW)
            strcpy(pszUnique, "Shell");
        else
        {
            ownerAttCount++;
            itoa(ownerAttCount, pszUnique, 10); // base 10
        }
    }

    // append unique part to attribute base name
    strcat(pszAttribName, pszUnique);
    bsAttribName = ::SysAllocString(A2OLE(pszAttribName));

	IModelDoc2 *modDoc = NULL;

	hr = SWQueryEngine::instance()->GetSldWorks()->get_IActiveDoc2(&modDoc);
    assert(modDoc);

    hr = attDef->ICreateInstance3(modDoc, NULL, pEntity, bsAttribName,
                                  0, 0, &attInstance );
    if (pEntity)
        pEntity->Release();
    pEntity = NULL;

    if (modDoc)
        modDoc->Release();


	IParameter *paramInstance = NULL;
    BSTR bsParamName = ::SysAllocString(A2OLE("Cubit_Owner"));
	hr = attInstance->IGetParameter ( bsParamName, &paramInstance );
    ::SysFreeString(bsParamName);
    attInstance->Release();

	double* cubitEntityNum = (double*)this; //cast to double
	VARIANT_BOOL paramInstanceBool = FALSE;
	hr = paramInstance->SetDoubleValue ( *cubitEntityNum, &paramInstanceBool );
    paramInstance->Release();

	if(paramInstanceBool == FALSE)
	{
		//PRINT_ERROR("Problems setting the double value for the parameter.\n");
		//return;
	}
}

SW_ATTRIB_CUBIT_OWNER::~SW_ATTRIB_CUBIT_OWNER()
{
    assert(false);
    if (bsAttribName)
    {
        ::SysFreeString(bsAttribName);
        bsAttribName = NULL;
    }
    if (attInstance)
        attInstance->Release();
    attInstance = NULL;
    if (attDef)
        attDef->Release();
    attDef = NULL;
}


// static function to initialize attribute definition in SolidWorks
HRESULT
SW_ATTRIB_CUBIT_OWNER::initialize()
{
    assert(false);
    USES_CONVERSION;
    HRESULT hr = NOERROR;

	assert(!initialized);
    if (initialized)
        return S_OK;

    BSTR bsAttDefName = ::SysAllocString(A2OLE("Cubit_Owner_Attribute"));

    // create an attribute definition in SolidWorks
    hr = SWQueryEngine::instance()->GetSldWorks()->IDefineAttribute ( bsAttDefName, &attDef );
    ::SysFreeString(bsAttDefName);
    if (FAILED(hr))
        return hr;

	initialized = CUBIT_TRUE;

    // add a parameter to the attribute definition
    BSTR bsParamName = ::SysAllocString(A2OLE("Cubit_Owner"));
	VARIANT_BOOL attDefParam = FALSE;
	hr = attDef->AddParameter ( bsParamName, swParamTypeDouble, 0.0, 0, &attDefParam );
    ::SysFreeString(bsParamName);
	if(FAILED(attDefParam == FALSE))
	{
		PRINT_ERROR("Problems adding the Cubit_Owner parameter.\n");
		return E_FAIL;
	}
	
	VARIANT_BOOL did_register = FALSE;
	hr = attDef->Register ( &did_register );
	if(FAILED(hr) || (did_register == FALSE))
	{
		PRINT_ERROR("Problems registering the Attribute Definition.\n");
		return E_FAIL;
	}

    return hr;
}

// static
SWBridge* SW_ATTRIB_CUBIT_OWNER::cubit_owner(IUnknown *entity)
{
    HRESULT hr = NOERROR;

    IBody2 *pSWBody;
    IFace2 *pSWFace;
    ILoop2 *pSWLoop;
    ICoEdge *pSWCoEdge;
    IEdge *pSWEdge;
    IVertex *pSWVertex;

    // TODO - DLIList uses an int for an index - why not a long?
    int indexSW;

    // TODO - make overloads for this method which take specific interfaces
    // rather than IUnknown??

    hr = entity->QueryInterface(IID_IFace2, (LPVOID*)&pSWFace);
    if (SUCCEEDED(hr))
    {
        indexSW = S_SWFaceList.where_is_item(pSWFace);
        if (indexSW < 0)
            return NULL;
        else
        {
            S_surfaceList.reset();
            S_surfaceList.step(indexSW);
            return S_surfaceList.get();
        }
    }
    else
    {
        hr = entity->QueryInterface(IID_ILoop2, (LPVOID*)&pSWLoop);
        if (SUCCEEDED(hr))
        {
            indexSW = S_SWLoopList.where_is_item(pSWLoop);
            if (indexSW < 0)
                return NULL;
            else
            {
                S_loopList.reset();
                S_loopList.step(indexSW);
                return S_loopList.get();
            }
        }
        else
        {
            hr = entity->QueryInterface(IID_ICoEdge, (LPVOID*)&pSWCoEdge);
            if (SUCCEEDED(hr))
            {
                indexSW = S_SWCoEdgeList.where_is_item(pSWCoEdge);
                if (indexSW < 0)
                    return NULL;
                else
                {
                    S_coedgeList.reset();
                    S_coedgeList.step(indexSW);
                    return S_coedgeList.get();
                }
            }
            else
            {
                hr = entity->QueryInterface(IID_IEdge, (LPVOID*)&pSWEdge);
                if (SUCCEEDED(hr))
                {
                    indexSW = S_SWEdgeList.where_is_item(pSWEdge);
                    if (indexSW < 0)
                        return NULL;
                    else
                    {
                        S_curveList.reset();
                        S_curveList.step(indexSW);
                        return S_curveList.get();
                    }
                }
                else
                {
                    hr = entity->QueryInterface(IID_IVertex, (LPVOID*)&pSWVertex);
                    if (SUCCEEDED(hr))
                    {
                        indexSW = S_SWVertexList.where_is_item(pSWVertex);
                        if (indexSW < 0)
                            return NULL;
                        else
                        {
                            S_pointList.reset();
                            S_pointList.step(indexSW);
                            return S_pointList.get();
                        }
                    }
                    else
                    {
                        hr = entity->QueryInterface(IID_IBody2, (LPVOID*)&pSWBody);
                        if (SUCCEEDED(hr))
                        {
                            // multiple CGM owners, so return none
                            return NULL;
                        }
                        else
                        {
                            assert(false);
                        }
                    }
                }// IEdge *
            }// ICoEdge *
        }// ILoop2 *
    }// IFace2 *

    return NULL;

//    USES_CONVERSION;

//    HRESULT hr = NOERROR;

//    LPENTITY pEntity = NULL;
//    LPATTRIBUTE attribute = NULL;
//    hr = entity->QueryInterface(IID_IEntity, (LPVOID*)&pEntity);
//    if (SUCCEEDED(hr) && (NULL != pEntity))
//    {
        // Get the attribute if there is one
//        hr = pEntity->IFindAttribute ( attDef, 1, &attribute );
//        // (uses static definition assigned in initialize function; only
//        // one instance per entity, so instance # is 1)
//    }
//    else
    // Since the SolidWorks Body object is not an "entity" the
    // attribute is placed directly on the document.  Look on the document
    // for a body attribute.
//    {
////        assert(false); //TODO - look on document for body attribute or lump attribute?
////        hr = entity->QueryInterface(IID_IBody, (LPVOID*)&pBody);
////        if (SUCCEEDED(hr))
////        {
////            LPMODELDOC modDoc = NULL;
////            GeometryQueryEngine *gqe = 
////              GeometryQueryTool::instance()->get_gqe(SWQueryEngine_TYPE);
////            SWQueryEngine *sqe = CAST_TO(gqe, SWQueryEngine);
////            assert(sqe);

////            hr = sqe->GetSldWorks()->get_IActiveDoc(&modDoc);
////            LPPARTDOC pPartDoc = NULL;
////            hr = modDoc->QueryInterface(IID_IPartDoc, (LPVOID*)&pPartDoc);
////
////            LPFEATURE pFeat = NULL;
////            char pszAttribName[] = "CUBIT Entity Body";
////            BSTR bsAttribName = ::SysAllocString(A2OLE(pszAttribName));
////            hr = pPartDoc->IFeatureByName(bsAttribName, &pFeat);
////        }

//    }


	// If there is an attribute and an owner, return it, else return NULL.
//	if (attribute == NULL)
//		return NULL;
//	else
//	{
//        BSTR bsParamName = ::SysAllocString(A2OLE("Cubit_Owner"));
//		LPPARAMETER param = NULL;
//		attribute->IGetParameter ( bsParamName, &param );
//        ::SysFreeString(bsParamName);

//		double doubleval = 0.0;
//		param->GetDoubleValue ( &doubleval );
//		double* doubleval_ptr = &doubleval;
//		SW_ATTRIB_CUBIT_OWNER *ACOattribute = 
//			(SW_ATTRIB_CUBIT_OWNER *)doubleval_ptr;
//		return ACOattribute->cubit_owner();
//	}
}

//static
void SW_ATTRIB_CUBIT_OWNER::cubit_owner(DLIList<IUnknown *> &entity_list,
                                     DLIList<SWBridge*> &tb_list)
{
    assert(false);
/*
	for (int i = entity_list.size(); i > 0; i--)
	{
		LPENTITY entity = entity_list.get_and_step();
		SWBridge *tb = cubit_owner(entity);
		if (tb != NULL)
			tb_list.append(tb);
	}
*/
}

// static
void SW_ATTRIB_CUBIT_OWNER::cubit_owner(DLIList<IUnknown *> &entity_list,
                                     DLIList<TopologyBridge*> &tb_list)
{
	for (int i = entity_list.size(); i > 0; i--)
	{
		IUnknown *entity = entity_list.get_and_step();
		SWBridge *tb = cubit_owner(entity);
		if (tb != NULL)
			tb_list.append(CAST_TO(tb, TopologyBridge));
	}
}

//static
void SW_ATTRIB_CUBIT_OWNER::cubit_owner( enum EntityType ownerType,
                           DLIList<TopologyBridge*> &tb_list )
{
    // TODO - this is currently implemented for part document so it
    // returns the single body, lump, or shell object
    switch (ownerType)
    {
    case Lump_TYPE:

        tb_list.append(CAST_TO(S_lump, TopologyBridge));
        return;

    case Shell_TYPE:
        tb_list.append(CAST_TO(S_shell, TopologyBridge));
        return;

    default:
        assert(false); // not handled
    }

    return;
}


// static
TopologyEntity* SW_ATTRIB_CUBIT_OWNER::get_topology_entity(IUnknown *entity)
{
    USES_CONVERSION;
    HRESULT hr = NOERROR;

    IEntity *pEntity;
    hr = entity->QueryInterface(IID_IEntity, (LPVOID*)&pEntity);
    if (FAILED(hr) || (NULL == pEntity))
        return NULL;

	// Get the attribute if there is one
	IAttribute *attribute = NULL;
	hr = pEntity->IFindAttribute ( attDef, 1, &attribute );
    hr = pEntity->Release();

	// If there is an attribute and an owner, return it, else return NULL.
	if (attribute == NULL)
		return NULL;
	else
	{
        BSTR bsParamName = ::SysAllocString(A2OLE("Cubit_Owner"));
		IParameter *param = NULL;
		hr = attribute->IGetParameter ( bsParamName, &param );
        ::SysFreeString(bsParamName);

		double doubleval = 0.0;
		param->GetDoubleValue ( &doubleval );
		double* doubleval_ptr = &doubleval;
		SW_ATTRIB_CUBIT_OWNER *ACOattribute = 
			(SW_ATTRIB_CUBIT_OWNER *)doubleval_ptr;

		if (attribute == NULL || ACOattribute->cubit_owner() == NULL)
			return NULL;
		else
		{
			TopologyBridge* tb = CAST_TO(ACOattribute->cubit_owner(), TopologyBridge);
			return tb->topology_entity();
		}
	}
}

// set the member data.
void SW_ATTRIB_CUBIT_OWNER::set_cubit_owner (SWBridge* new_cubit_owner)
{
    assert(false);
	//if (new_cubit_owner != NULL)
		//backup();
	cubitOwnerData = new_cubit_owner;
}

// static
void SW_ATTRIB_CUBIT_OWNER::set_cubit_owner(IUnknown *entity,
                                         SWBridge *cubit_entity)
{
    USES_CONVERSION;
    HRESULT hr = NOERROR;

    // Some checks
	if (cubit_entity == NULL)
	{
		PRINT_ERROR("Trying to set the CUBIT_OWNER attrribute of "
					"a SW entity to be NULL in set_cubit_owner\n");
		return;
	}

	if (entity == NULL)
	{
		PRINT_ERROR("Trying to set the CUBIT_OWNER attribute of a "
					" NULL SW entity in set_cubit_owner\n");
		return;
	}

    BodySW *body = CAST_TO(cubit_entity, BodySW);
    if (body)
    {
        S_body = cubit_entity;
        return;
    }

    LumpSW *lump = CAST_TO(cubit_entity, LumpSW);
    if (lump)
    {
        S_lump = cubit_entity;
        return;
    }

    ShellSW *shell = CAST_TO(cubit_entity, ShellSW);
    if (shell)
    {
        S_shell = cubit_entity;
        return;
    }

    // see if we've already got this entity
    // TODO - we may not need to do append_unique here if we already
    // check the list before calling set_cubit_owner

    int swStat;
    int stat;

    IFace2 *pSWFace;
    ILoop2 *pSWLoop;
    ICoEdge *pSWCoEdge;
    IEdge *pSWEdge;
    IVertex *pSWVertex;

    hr = entity->QueryInterface(IID_IFace2, (LPVOID*)&pSWFace);
    if (SUCCEEDED(hr))
    {
        swStat = S_SWFaceList.append_unique(pSWFace);
        assert(CUBIT_TRUE == swStat);
        stat = S_surfaceList.append_unique(cubit_entity);
        assert(CUBIT_TRUE == stat);
    }
    else
    {
        hr = entity->QueryInterface(IID_ILoop2, (LPVOID*)&pSWLoop);
        if (SUCCEEDED(hr))
        {
            swStat = S_SWLoopList.append_unique(pSWLoop);
            assert(CUBIT_TRUE == swStat);
            stat = S_loopList.append_unique(cubit_entity);
            assert(CUBIT_TRUE == stat);
        }
        else
        {
            hr = entity->QueryInterface(IID_ICoEdge, (LPVOID*)&pSWCoEdge);
            if (SUCCEEDED(hr))
            {
                swStat = S_SWCoEdgeList.append_unique(pSWCoEdge);
                assert(CUBIT_TRUE == swStat);
                stat = S_coedgeList.append_unique(cubit_entity);
                assert(CUBIT_TRUE == stat);
            }
            else
            {
                hr = entity->QueryInterface(IID_IEdge, (LPVOID*)&pSWEdge);
                if (SUCCEEDED(hr))
                {
                    swStat = S_SWEdgeList.append_unique(pSWEdge);
                    assert(CUBIT_TRUE == swStat);
                    stat = S_curveList.append_unique(cubit_entity);
                    assert(CUBIT_TRUE == stat);
                }
                else
                {
                    hr = entity->QueryInterface(IID_IVertex, (LPVOID*)&pSWVertex);
                    if (SUCCEEDED(hr))
                    {
                        swStat = S_SWVertexList.append_unique(pSWVertex);
                        assert(CUBIT_TRUE == swStat);
                        stat = S_pointList.append_unique(cubit_entity);
                        assert(CUBIT_TRUE == stat);
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

/*
	// Search for the attribute
	LPATTRIBUTE temp_attribute = NULL;
    LPENTITY pEntity = NULL;
    hr = entity->QueryInterface(IID_IEntity, (LPVOID*)&pEntity);
    if (SUCCEEDED(hr) && (NULL != pEntity))
    {
	    hr = pEntity->IFindAttribute ( attDef, 1, &temp_attribute );
        pEntity->Release();
    }

	// If found, reset it's value
	if (temp_attribute != NULL)
	{
        BSTR bsParamName = ::SysAllocString(A2OLE("Cubit_Owner"));
		LPPARAMETER temp_param = NULL;
		temp_attribute->IGetParameter ( bsParamName, &temp_param );
        ::SysFreeString(bsParamName);

		double doubleval = 0.0;
		temp_param->GetDoubleValue ( &doubleval );
		double* doubleval_ptr = &doubleval;
		SW_ATTRIB_CUBIT_OWNER *temp_attrib = 
			(SW_ATTRIB_CUBIT_OWNER *)doubleval_ptr;
		temp_attrib->set_cubit_owner(cubit_entity);
	}

	// Otherwise, create a new attribute for SW entity
	else
	{
              // create SW attribute
              // assign double parameter in attribute to cubit_entity

		new SW_ATTRIB_CUBIT_OWNER(entity, cubit_entity);
	}
*/
}

// static
void
SW_ATTRIB_CUBIT_OWNER::remove_cubit_owner(IUnknown *entity,
                                       CubitBoolean children_too)
{
    HRESULT hr = NOERROR;
    IBody2 *pSWBody = NULL;
    IFace2 *pSWFace = NULL;
    ILoop2 *pSWLoop = NULL;
    ICoEdge *pSWCoEdge = NULL;
    IEdge *pSWEdge = NULL;
    IVertex *pSWVertex = NULL;

    SWBridge* pBridge = NULL;
    pBridge = SW_ATTRIB_CUBIT_OWNER::cubit_owner(entity);
    if (pBridge)
    {
        hr = entity->QueryInterface(IID_IFace2, (LPVOID*)&pSWFace);
        if (SUCCEEDED(hr))
        {
            S_surfaceList.remove(pBridge);
            S_SWFaceList.remove((IFace2*)entity);
            pSWFace->Release();
        }

        hr = entity->QueryInterface(IID_ILoop2, (LPVOID*)&pSWLoop);
        if (SUCCEEDED(hr))
        {
            S_loopList.remove(pBridge);
            S_SWLoopList.remove((ILoop2 *)entity);
            pSWLoop->Release();
        }

        hr = entity->QueryInterface(IID_ICoEdge, (LPVOID*)&pSWCoEdge);
        if (SUCCEEDED(hr))
        {
            S_coedgeList.remove(pBridge);
            S_SWCoEdgeList.remove((ICoEdge *)entity);
            pSWCoEdge->Release();
        }

        hr = entity->QueryInterface(IID_IEdge, (LPVOID*)&pSWEdge);
        if (SUCCEEDED(hr))
        {
            S_curveList.remove(pBridge);
            S_SWEdgeList.remove((IEdge *)entity);
            pSWEdge->Release();
        }

        hr = entity->QueryInterface(IID_IVertex, (LPVOID*)&pSWVertex);
        if (SUCCEEDED(hr))
        {
            S_pointList.remove(pBridge);
            S_SWVertexList.remove((IVertex *)entity);
            pSWVertex->Release();
        }

//        pBridge->ENTITY_ptr(NULL);
    }
    else
    {
        hr = entity->QueryInterface(IID_IBody2, (LPVOID*)&pSWBody);
        if (SUCCEEDED(hr))
        {
            BodySW* pBodySW = CAST_TO(pBridge, BodySW);
            if (pBodySW)
            {
                S_body = NULL;
            }

            LumpSW* pLumpSW = CAST_TO(pBridge, LumpSW);
            if (pLumpSW)
            {
                S_lump = NULL;
            }

            ShellSW* pShellSW = CAST_TO(pBridge, ShellSW);
            if (pShellSW)
            {
                S_shell = NULL;
            }

            pSWBody->Release();
        }
    }
}

// Virtual function called when an owner entity is being split,
// such as from a subtract operation.
// Also called when an ENTITY is copied in part, such as when
// a FACE is copied without also copying a LUMP to attach it to.
/*void SW_ATTRIB_CUBIT_OWNER::split_owner(LPDISPATCH entityDisp)
{
	//   PRINT_INFO("split_owner called for attrib_cubit_owner on ");

	// must get simple attributes for new entity from CubitAttribUser
	if (!cubitOwnerData)
		return;

	TopologyBridge* tb = CAST_TO(cubitOwnerData, TopologyBridge);
	CubitEntity* cep = CAST_TO(tb->topology_entity(), CubitEntity);
	//   if (cep)
	//     PRINT_INFO("%s %d\n",
	//                cep->class_name(),
	//                cep->id());
	//   else
	//     PRINT_INFO("unknown ENTITY.\n");
	CubitAttribUser *cau = CAST_TO(tb->topology_entity(), CubitAttribUser);

	if (cau != NULL && cau->num_cubit_attrib() != 0)
	{
		DLCubitSimpleAttribList csa_list;
		cau->split_owner(csa_list);

		  // now put those csa's on the new entity, if there are any
		if (csa_list.size() != 0)
		{
			for (int i = csa_list.size(); i > 0; i--)
			{
				CubitSimpleAttrib *csa_ptr = csa_list.get_and_step();
				assert(csa_ptr != NULL);
				new ATTRIB_SNL_SIMPLE ( entityDisp,
										csa_ptr->string_data_list(),
										csa_ptr->double_data_list(),
										csa_ptr->int_data_list());
			}
		}
	}

	// in the case of a split, lose this attribute
	//unhook();
	//lose();
}*/

// Virtual function called when two entities are to be merged,
// as may happen with 'unite'
void SW_ATTRIB_CUBIT_OWNER::merge_owner(IUnknown *entity, 
									 VARIANT_BOOL delete_this)
{
    assert(false);
/*
    USES_CONVERSION;
    HRESULT hr = NOERROR;

	//   PRINT_INFO("merge_owner called for attrib_cubit_owner.\n");
	// Handle merging via the survivor
	if (delete_this)
		return;

	// Get CAU
	if (!cubitOwnerData)
		return;
	TopologyBridge* tb = CAST_TO(cubitOwnerData, TopologyBridge);
	CubitAttribUser *cau =
		(tb ? CAST_TO(tb->topology_entity(), CubitAttribUser) :
		NULL);
	if (cau == NULL || cau->num_cubit_attrib() == 0)
		return;
	
	// get the owner attribute on other_entity
	LPATTRIBUTE other_owner_att = NULL;
    LPENTITY pEntity = NULL;
    hr = entity->QueryInterface(IID_IEntity, (LPVOID*)&pEntity);
    if (SUCCEEDED(hr) && (NULL != pEntity))
    {
        hr = pEntity->IFindAttribute ( attDef, 1, &other_owner_att );
        pEntity->Release();
    }
	RefEntity *ref_entity = NULL;

	if (other_owner_att == NULL)
		return;

    BSTR bsParamName = ::SysAllocString(A2OLE("Cubit_Owner"));
	LPPARAMETER param = NULL;
	other_owner_att->IGetParameter ( bsParamName, &param );
    ::SysFreeString(bsParamName);

	double doubleval = 0.0;
	param->GetDoubleValue ( &doubleval );
	double* doubleval_ptr = &doubleval;
	SW_ATTRIB_CUBIT_OWNER *ACO_other_owner_att = 
		(SW_ATTRIB_CUBIT_OWNER *)doubleval_ptr;

	if (ACO_other_owner_att && ACO_other_owner_att->cubit_owner()) 
	{
		tb = CAST_TO(ACO_other_owner_att->cubit_owner(), TopologyBridge);
		ref_entity = CAST_TO(tb->topology_entity(), RefEntity);
		assert(ref_entity);
	}

	cau->merge_owner(ref_entity);
*/
}

// Virtual function called when the owner is being translated
void SW_ATTRIB_CUBIT_OWNER::trans_owner(double transform[16])
{
    assert(false);
/*
	//   PRINT_INFO("trans_owner called for attrib_cubit_owner.\n");
	if (!cubitOwnerData)
		return;
	TopologyBridge* tb = CAST_TO(cubitOwnerData, TopologyBridge);
	if (tb == NULL)
		return;
	CubitAttribUser *cau = CAST_TO(tb->topology_entity(), CubitAttribUser);
	if (cau == NULL || cau->num_cubit_attrib() == 0)
		return;

	// dummy variables representing transform
	CubitVector translate_vec(transform[9],
							  transform[10],
							  transform[11]);

	CubitVector matrow1(transform[0],
						transform[1],
						transform[2]);

	CubitVector matrow2(transform[3],
						transform[4],
						transform[5]);

	CubitVector matrow3(transform[6],
						transform[7],
						transform[8]);

	double scale_factor = transform[12];

	cau->transf_owner(matrow1, matrow2, matrow3, translate_vec, scale_factor);
*/
}

