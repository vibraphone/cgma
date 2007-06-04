//-------------------------------------------------------------------------
// Filename      : SWQueryEngine.cpp
//
// Purpose       : Performs all query-type operations on SolidWorks geometry
//
// Special Notes :
//
// Creator       : Byron Hanks
//
// Creation Date : 4/01
//
// Owner         : Byron Hanks
//-------------------------------------------------------------------------

#include "stdafx.h"

#include "swconst.h"
//#include "amapp.h"
#import "C:\Program Files\SolidWorks\sldworks.tlb" raw_interfaces_only, raw_native_types, no_namespace, named_guids     //Import the SolidWorks type library


#include "gtfordm.h"

#include "SWQueryEngine.hpp"

#include "GeometryQueryTool.hpp"

#include "BodySW.hpp"
#include "LumpSW.hpp"
#include "ShellSW.hpp"
#include "SurfaceSW.hpp"
#include "LoopSW.hpp"
#include "CoEdgeSW.hpp"
#include "CurveSW.hpp"
#include "PointSW.hpp"

#include "SWPart.hpp"
#include "RefEdge.hpp"
#include "RefFace.hpp"
#include "Body.hpp"

#include "GMem.hpp"

#include "CGMEngineDynamicLoader.hpp"

CGM_ENGINE_EXPORT_CREATE_GQE(SolidWorks)
{
  return SWQueryEngine::instance();
}

SWQueryEngine* SWQueryEngine::instance_ = 0;

SWQueryEngine::SWQueryEngine()
{
    assert(NULL == instance_);

    // add this engine to geometryquerytool
    GeometryQueryTool::instance()->add_gqe(this);

    m_pSWApp = NULL;
    m_pAttDef = NULL;

    m_pImport = NULL;
}

SWQueryEngine::~SWQueryEngine()
{
    if (m_pImport)
        delete m_pImport;
    m_pImport = NULL;

    if (m_pAttDef)
        m_pAttDef->Release();
    m_pAttDef = NULL;

    if (m_pSWApp)
        m_pSWApp->Release();
    m_pSWApp = NULL;

    instance_ = 0;
}

CubitStatus
SWQueryEngine::initialize()
{
    USES_CONVERSION;

    assert(NULL == m_pImport);
    m_pImport = new SWImport();
    if (NULL == m_pImport)
        return CUBIT_FAILURE;

    // get a pointer to the SolidWorks application object
    CLSID appCLSID;
    char pszAppName[] = "SldWorks.Application";
    HRESULT hr = CLSIDFromProgID(A2OLE(pszAppName), &appCLSID);
    hr = ::CoCreateInstance(appCLSID, NULL, CLSCTX_LOCAL_SERVER,
                            IID_ISldWorks, (LPVOID*)&m_pSWApp);
    if (FAILED(hr))
        return CUBIT_FAILURE;

    hr = m_pSWApp->put_Visible(TRUE);

    // initialize attributes

//    BSTR bsAttDefName = ::SysAllocString(A2OLE("CubitAttribute"));

    // create an attribute definition in SolidWorks
//    LPATTRIBUTEDEF attDef = NULL;
//    hr = m_pSWApp->IDefineAttribute ( bsAttDefName, &attDef );
//    ::SysFreeString(bsAttDefName);
//    if (SUCCEEDED(hr) && attDef)
//    {
        // register the attribute definition
//        VARIANT_BOOL did_register = FALSE;
//        hr = attDef->Register ( &did_register );
//        if (did_register)
//            m_pAttDef = attDef;
//        else
//            attDef->Release();
//    }

    return CUBIT_SUCCESS;
}

ISldWorks *SWQueryEngine::GetSldWorks()
{
    if (NULL == m_pSWApp)
        this->initialize();
	assert(m_pSWApp);
	return m_pSWApp;
}
/*
LPATTRIBUTEDEF SWQueryEngine::GetAttributeDef()
{
	assert(m_pAttDef);
	return m_pAttDef;
}
*/
CubitStatus 
SWQueryEngine::export_solid_model( DLIList<TopologyBridge*>& bridge_list,
                                   const char* file_name,
                                   const char* file_type,
                                   const CubitString &cubit_version,
                                   const char* logfile_name ) 
{
    assert(false);
   return CUBIT_FAILURE;
}

CubitStatus SWQueryEngine::save_temp_geom_file(DLIList<TopologyBridge*> &ref_entity_list,
                                               const char *filename,
                                               const CubitString &cubit_version,
                                               CubitString &created_file,
                                               CubitString &created_file_type ) 
{
   return CUBIT_FAILURE;
}




CubitStatus 
SWQueryEngine::import_temp_geom_file(FILE* file_ptr, 
                                     const char* file_name,
                                     const char* file_type,
                                     DLIList<TopologyBridge*> &bridge_list)
{
  return CUBIT_FAILURE;
}

CubitStatus 
SWQueryEngine::import_solid_model(const char* file_name,
                                  const char* file_type,
                                  DLIList<TopologyBridge*>& imported_entities,
                                  CubitBoolean print_results,
                                  const char* logfile_name,
                                  CubitBoolean heal_step,
                                  CubitBoolean import_bodies,
                                  CubitBoolean import_surfaces,
                                  CubitBoolean import_curves,
                                  CubitBoolean import_vertices,
                                  CubitBoolean free_surfaces) 
{
    USES_CONVERSION;

    IModelDoc2 *pSWDoc = NULL;
    IBody2 *pSWBody = NULL;
    HRESULT hr;
    long lDocType = swDocNONE;

// in order to do this and initialize only when needed, import_solid_model must be non const
// Why is it a const member function?????????
//    if (NULL == m_pSWApp)
//        this->initialize();
    
    // make sure SolidWorks application has been started
    if (NULL == m_pSWApp)
        return CUBIT_FAILURE;

    // make sure import object has been created
    if (NULL == m_pImport)
        return CUBIT_FAILURE;

    // if a file name was passed in, open the file in SolidWorks
    if (file_name && strlen(file_name))
    {
        BSTR bsFileName = A2BSTR(file_name);

        long lOpenErr;
        hr = m_pSWApp->IOpenDoc5(bsFileName, swDocPART, swOpenDocOptions_Silent,
                                 NULL, &lOpenErr, &pSWDoc);
        ::SysFreeString(bsFileName);
        if (FAILED(hr) || (NULL == pSWDoc))
        {
            return CUBIT_FAILURE;
        }
    }
    // if no file name, then get the active SolidWorks document
    else
    {
        hr = m_pSWApp->get_IActiveDoc2(&pSWDoc);
        if (FAILED(hr) || (NULL == pSWDoc))
        {
            return CUBIT_FAILURE;
        }
    }

    return m_pImport->import_solid_model(pSWDoc, imported_entities);
}

CubitStatus SWQueryEngine::get_intersections(Curve* ref_edge1, CubitVector& point1,
                              CubitVector& point2,
                              DLIList<CubitVector*>& intersection_list,
                              bool bounded,
                              bool closest )
{
    assert(false);
    return CUBIT_FAILURE;
}

CubitStatus SWQueryEngine::get_intersections( Curve* ref_edge1, Curve* ref_edge2,
                                         DLIList<CubitVector*>& intersection_list,
                                         CubitBoolean bounded,
                                         CubitBoolean closest)
{
    assert(false);
    return CUBIT_FAILURE;
}

CubitStatus SWQueryEngine::get_intersections( Curve* ref_edge, Surface* ref_face,
                                              DLIList<CubitVector*>& intersection_list,
                                              CubitBoolean bounded )
{
    assert(false);
    return CUBIT_FAILURE;
}

//================================================================================
// Description: Find extrema position on an entity list
// Author     : 
// Date       : 
//================================================================================
CubitStatus 
SWQueryEngine::entity_extrema( DLIList<GeometryEntity*> &, 
                               const CubitVector *, const CubitVector *,
                               const CubitVector *, CubitVector &,
                               GeometryEntity *& )
{
  PRINT_ERROR("Option not supported for SolidWorks based geometry.\n");
  return CUBIT_FAILURE;
}
 
//================================================================================
// Description: Find distance between two entities and closest positions.
// Author     : 
// Date       : 
//================================================================================
CubitStatus
SWQueryEngine::entity_entity_distance( GeometryEntity *,
                                       GeometryEntity *,
                                       CubitVector &, CubitVector &,
                                       double & )
{
  PRINT_ERROR("Option not supported for SolidWorks based geometry.\n");
  return CUBIT_FAILURE;
}

void SWQueryEngine::delete_solid_model_entities(DLIList<BodySM*>& body_list) const
{
    assert(false);
}

CubitStatus SWQueryEngine::delete_solid_model_entities( BodySM* body_ptr ) const
{
    assert(false);
    return CUBIT_SUCCESS;
}

CubitStatus SWQueryEngine::delete_solid_model_entities(Surface* surf_ptr ) const
{
    assert(false);
    return CUBIT_SUCCESS;
}

CubitStatus SWQueryEngine::delete_solid_model_entities( Curve* curve_ptr ) const
{
    assert(false);
    return CUBIT_SUCCESS;
}

CubitStatus SWQueryEngine::delete_solid_model_entities( Point* point_ptr ) const
{
    assert(false);
    return CUBIT_SUCCESS;
}


CubitStatus SWQueryEngine::fire_ray(BodySM *body,
                                   const CubitVector &ray_point,
                                   const CubitVector &unit_direction,
                                   DLIList<double> &ray_params,
                                   DLIList<GeometryEntity*> *entity_list) const
{
    assert(false);
    return CUBIT_SUCCESS;
}

CubitStatus SWQueryEngine::get_isoparametric_points(Surface* ref_face_ptr,
                                               int &nu, int &nv,
                                               GMem *&gMem) const
{
    assert(false);
    return CUBIT_SUCCESS;
}

CubitStatus SWQueryEngine::get_u_isoparametric_points(Surface* ref_face_ptr,
                                                 double v, int& n,
                                                 GMem *&gMem) const
{
    assert(false);
    return CUBIT_SUCCESS;
}

CubitStatus SWQueryEngine::get_v_isoparametric_points(Surface* ref_face_ptr,
                                                 double u, int&n,
                                                 GMem *&gMem) const
{
    assert(false);
    return CUBIT_SUCCESS;
}
  
CubitStatus SWQueryEngine::transform_vec_position( 
    CubitVector const& position_vector,
    BodySM *OSME_ptr,
    CubitVector &transform_vector ) const
{
    assert(false);
    return CUBIT_SUCCESS;
}

double SWQueryEngine::get_sme_resabs_tolerance() const // Gets solid modeler's resolution absolute tolerance
{
    assert(false);
    return GEOMETRY_RESABS;
}

double SWQueryEngine::set_sme_resabs_tolerance( double new_resabs )
{
    assert(false);
    return GEOMETRY_RESABS;
}

CubitStatus SWQueryEngine::set_int_option( const char* opt_name, int val )
{
    assert(false);
    return CUBIT_SUCCESS;
}

CubitStatus SWQueryEngine::set_dbl_option( const char* opt_name, double val )
{
    assert(false);
    return CUBIT_SUCCESS;
}

CubitStatus SWQueryEngine::set_str_option( const char* opt_name, const char* val )
{
    assert(false);
    return CUBIT_SUCCESS;
}

//CubitString SWQueryEngine::identify()
//{
//   CubitString version_string = "SolidWorks Version ";
//   return version_string;
//}

CubitStatus SWQueryEngine::get_graphics( Surface* surface_ptr,
                                         int& number_of_triangles,
                                         int& number_of_points,
                                         int& number_of_facets,
                                         GMem* gMem,
                                         unsigned short normal_tolerance,
                                         double distance_tolerance,
                                         double longest_edge) const
{
  CubitStatus rv;

  //CurveSW const* pCurveSW = CAST_TO(pCurve, const CurveSW);
  SurfaceSW *pSurfSW = dynamic_cast<SurfaceSW *>(surface_ptr);
  assert(pSurfSW);

  rv = pSurfSW->facet_face(number_of_triangles, number_of_points,
                           number_of_facets, gMem, normal_tolerance,
                           distance_tolerance);

  if (rv == CUBIT_FAILURE)
  PRINT_ERROR("Unable to facet surface\n");
  return rv;
}

CubitStatus SWQueryEngine::get_graphics( Curve * curve_ptr,
                                         int& num_points,
                                         GMem* gMem,
                                         double tolerance ) const
{
  double resabs = CUBIT_RESABS;
  double resfit = 0.0001; // TODO - make tolerance relative

  CubitStatus rv;

  CurveSW *pCurveSW = dynamic_cast<CurveSW *>(curve_ptr);
  assert(pCurveSW);

  if( tolerance != 0.0 )
  {
    if( resabs >= resfit )
      tolerance = resabs * 2;
    else
      tolerance = resfit;
  }

  rv = pCurveSW->facet_edge(tolerance, num_points, gMem);

  if (rv == CUBIT_FAILURE)
		PRINT_ERROR("Unable to facet curve\n");
	return rv;
}

// This method performs topology traversal given a starting SolidWorks entity and the type
// to find in the traversal.
// NOTE - release the interface pointers after calling this routine
CubitStatus SWQueryEngine::get_ENTITY_from_ENTITY(IID queryIID,
                                                     IUnknown *entity,
                                                     DLIList<IUnknown *> &entities) const
{
    HRESULT       hr = NOERROR;
    IBody2        *pSWBody = NULL;
    IFace2        *pSWFace = NULL;
    IFace2        *pSWFace2 = NULL;
    IEnumFaces2   *pFaceEnum = NULL;
    ILoop2        *pSWLoop = NULL;
    IEnumLoops2   *pLoopEnum = NULL;
    ICoEdge       *pSWCoEdge = NULL;
    IEnumCoEdges  *pCoEdgeEnum = NULL;
    IEdge         *pSWEdge = NULL;
    IEnumEdges    *pEdgeEnum = NULL;
    IVertex       *pSWVertex = NULL;
    IVertex       *pSWStartVtx = NULL;
    IVertex       *pSWEndVtx = NULL;
    long nFetched = 0;
    CubitStatus status = CUBIT_FAILURE;

    hr = entity->QueryInterface(IID_IBody2, (LPVOID*)&pSWBody);
    if (SUCCEEDED(hr))
    {
        if (IID_IFace2 == queryIID)
        {
            hr = pSWBody->EnumFaces( &pFaceEnum);
            if (SUCCEEDED(hr) && pFaceEnum)
            {
                while (S_OK == pFaceEnum->Next(1, &pSWFace, &nFetched))
                {
                    entities.append(pSWFace);
                }
                pFaceEnum->Release();
                status = CUBIT_SUCCESS; // body must have at least one face
            }
        }
        else
        {
            assert(false); // unhandled case
        }
        pSWBody->Release();
    }
    else
    {
        hr = entity->QueryInterface(IID_IFace2, (LPVOID*)&pSWFace);
        if (SUCCEEDED(hr))
        {
            if (IID_ILoop2 == queryIID)
            {
                hr = pSWFace->EnumLoops(&pLoopEnum);
                if (SUCCEEDED(hr) && pLoopEnum)
                {
                    while (S_OK == pLoopEnum->Next(1, &pSWLoop, &nFetched))
                    {
                        entities.append(pSWLoop);
                    }
                    pLoopEnum->Release();
                }
                status = CUBIT_SUCCESS; // face can have zero or more loops
            }
            else if (IID_IEdge == queryIID)
            {
                hr = pSWFace->EnumEdges(&pEdgeEnum);
                if (SUCCEEDED(hr) && pEdgeEnum)
                {
                    while (S_OK == pEdgeEnum->Next(1, &pSWEdge, &nFetched))
                    {
                        entities.append(pSWEdge);
                    }
                    pEdgeEnum->Release();
                }
                status = CUBIT_SUCCESS; // face can have zero or more edges
            }
            else
            {
                assert(false); // unhandled case
            }
            pSWFace->Release();
        }
        else
        {
            hr = entity->QueryInterface(IID_ILoop2, (LPVOID*)&pSWLoop);
            if (SUCCEEDED(hr))
            {
                if (IID_ICoEdge == queryIID)
                {
                    hr = pSWLoop->EnumCoEdges(&pCoEdgeEnum);
                    if (SUCCEEDED(hr) && pCoEdgeEnum)
                    {
                        while (S_OK == pCoEdgeEnum->Next(1, &pSWCoEdge, &nFetched))
                        {
                            entities.append(pSWCoEdge);
                        }
                        pCoEdgeEnum->Release();
                        status = CUBIT_SUCCESS; // loop must have 1 or more coedges
                    }
                }
                else
                {
                    assert(false);
                }
                pSWLoop->Release();
            }
            else
            {
                hr = entity->QueryInterface(IID_ICoEdge, (LPVOID*)&pSWCoEdge);
                if (SUCCEEDED(hr))
                {
                    if (IID_IEdge == queryIID)
                    {
                        hr = pSWCoEdge->IGetEdge(&pSWEdge);
                        entities.append(pSWEdge);
                        status = CUBIT_SUCCESS;
                    }
                    else if (IID_IFace2 == queryIID)
                    {
                        hr = pSWCoEdge->IGetLoop2(&pSWLoop);
                        hr = pSWLoop->IGetFace(&pSWFace);
                        pSWLoop->Release();
                        entities.append(pSWFace);
                        status = CUBIT_SUCCESS;
                    }
                    else if (IID_ILoop2 == queryIID)
                    {
                        hr = pSWCoEdge->IGetLoop2(&pSWLoop);
                        entities.append(pSWLoop);
                        status = CUBIT_SUCCESS;
                    }
                    else
                    {
                        assert(false);
                    }
                    pSWCoEdge->Release();
                }
                else
                {
                    hr = entity->QueryInterface(IID_IEdge, (LPVOID*)&pSWEdge);
                    if (SUCCEEDED(hr))
                    {
                        if (IID_IVertex == queryIID)
                        {
                            hr = pSWEdge->IGetStartVertex(&pSWStartVtx);
                            hr = pSWEdge->IGetEndVertex(&pSWEndVtx);

                            // check for periodic edge (no vertex)
                            if (pSWStartVtx)
                            {
                                entities.append(pSWStartVtx);

                                assert(pSWEndVtx);
                                if (pSWEndVtx != pSWStartVtx)
                                {
                                    entities.append(pSWEndVtx);
                                }
                                else
                                {
                                    pSWEndVtx->Release();
                                }
                            }
                            status = CUBIT_SUCCESS;
                        }
                        else if (IID_IFace2 == queryIID)
                        {
                            hr = pSWEdge->IGetTwoAdjacentFaces2(&pSWFace, &pSWFace2);
                            entities.append(pSWFace);
                            entities.append(pSWFace2);
                            status = CUBIT_SUCCESS;
                        }
                        else
                        {
                            assert(false);
                        }
                        pSWEdge->Release();
                    }
                    else
                    {
                        hr = entity->QueryInterface(IID_IVertex, (LPVOID*)&pSWVertex);
                        if (SUCCEEDED(hr))
                        {
                            if (IID_IEdge == queryIID)
                            {
                                hr = pSWVertex->EnumEdges(&pEdgeEnum);
                                if (SUCCEEDED(hr) && pEdgeEnum)
                                {
                                    while (S_OK == pEdgeEnum->Next(1, &pSWEdge, &nFetched))
                                    {
                                        entities.append(pSWEdge);
                                    }
                                    pEdgeEnum->Release();
                                    status = CUBIT_SUCCESS; // vertex must have at least one edge
                                }
                            }
                            else
                            {
                                assert(false);
                            }
                            pSWVertex->Release();
                        }
                        else
                        {
                            assert(false);
                        } // IVertex *
                    }// IEdge *
                }// ICoEdge *
            }// ILoop2 *
        }// IFace2 *
    }// LPBODY

    return status;
}


CubitStatus SWQueryEngine::bodysms(IUnknown *entity, DLIList<BodySM*> &bodies) const
{
    assert(false);
    return CUBIT_FAILURE;
}

CubitStatus SWQueryEngine::lumps(SWPart *pSWPart, IUnknown *entity, DLIList<Lump*> &lumps) const
{
    // there is exactly one lump per part in SolidWorks, so no matter what entity asks for it
    // just return the one stored in the SWpart
    DLIList<TopologyBridge*> tb_list;
    pSWPart->lumps(tb_list);
    assert(1 == tb_list.size());
    CAST_LIST(tb_list, lumps, Lump);
    return CUBIT_SUCCESS;
}

CubitStatus SWQueryEngine::shellsms(SWPart *pSWPart, IUnknown *entity, DLIList<ShellSM*> &shellsms) const
{
    // solidworks doesn't have the concept of a shell in their API.  For now, just return the
    // shell from the SWPart - this only works for parts with a single shell
	DLIList<TopologyBridge*> tb_list;
    pSWPart->shells(tb_list);
	CAST_LIST(tb_list, shellsms, ShellSM);
    return CUBIT_SUCCESS;
}

CubitStatus SWQueryEngine::surfaces(SWPart *pSWPart, IUnknown *entity, DLIList<Surface*> &surfaces) const
{
    DLIList<IUnknown *> entities;
    CubitStatus status = get_ENTITY_from_ENTITY(IID_IFace2, entity, entities);
    if (status == CUBIT_FAILURE) 
        return status;

    DLIList<TopologyBridge*> tb_list;
    pSWPart->cubit_owner(entities, tb_list);
    CAST_LIST(tb_list, surfaces, Surface);

    IUnknown *pUnk;
	for (int i=entities.size(); i>0; i--)
	{
		pUnk = entities.get_and_step();
        pUnk->Release();
	}

    return CUBIT_SUCCESS;
}

CubitStatus SWQueryEngine::loopsms(SWPart *pSWPart, IUnknown *entity, DLIList<LoopSM*> &loopsms) const
{
	DLIList<IUnknown *> entities;
	CubitStatus status = get_ENTITY_from_ENTITY(IID_ILoop2, entity, entities);
	if (status == CUBIT_FAILURE) 
		return status;

	DLIList<TopologyBridge*> tb_list;
    pSWPart->cubit_owner(entities, tb_list);
	CAST_LIST(tb_list, loopsms, LoopSM);

    IUnknown *pUnk;
	for (int i=entities.size(); i>0; i--)
	{
		pUnk = entities.get_and_step();
        pUnk->Release();
	}

	return CUBIT_SUCCESS;
}

CubitStatus SWQueryEngine::curves(SWPart *pSWPart, IUnknown *entity, DLIList<Curve*> &curves) const
{
	DLIList<IUnknown *> entities;
	CubitStatus status = get_ENTITY_from_ENTITY(IID_IEdge, entity, entities);
	if (status == CUBIT_FAILURE) 
		return status;

	DLIList<TopologyBridge*> tb_list;
    pSWPart->cubit_owner(entities, tb_list);
	CAST_LIST(tb_list, curves, Curve);

    IUnknown *pUnk;
	for (int i=entities.size(); i>0; i--)
	{
		pUnk = entities.get_and_step();
        pUnk->Release();
	}

	return CUBIT_SUCCESS;
}

CubitStatus SWQueryEngine::coedgesms(SWPart *pSWPart, IUnknown *entity, DLIList<CoEdgeSM*> &coedgesms) const
{
	DLIList<IUnknown *> entities;
	CubitStatus status = get_ENTITY_from_ENTITY(IID_ICoEdge, entity, entities);
	if (status == CUBIT_FAILURE) 
		return status;

	DLIList<TopologyBridge*> tb_list;
    pSWPart->cubit_owner(entities, tb_list);
	CAST_LIST(tb_list, coedgesms, CoEdgeSM);

    IUnknown *pUnk;
	for (int i=entities.size(); i>0; i--)
	{
		pUnk = entities.get_and_step();
        pUnk->Release();
	}

	return CUBIT_SUCCESS;
}

CubitStatus SWQueryEngine::points(SWPart *pSWPart, IUnknown *entity, DLIList<Point*> &points) const
{
	DLIList<IUnknown *> entities;
	CubitStatus status = get_ENTITY_from_ENTITY(IID_IVertex, entity, entities);
	if (status == CUBIT_FAILURE)
		return status;

    // if none found, this could be a closed edge with a pseudo vertex
    if (0 == entities.size())
    {
        HRESULT hr = NOERROR;
        IEdge *pSWEdge = NULL;
        hr = entity->QueryInterface(IID_IEdge, (LPVOID*)&pSWEdge);
        if (SUCCEEDED(hr))
        {
            TopologyBridge *pPoint = pSWPart->findPseudoVertex(pSWEdge);

            Point *point;
            point = CAST_TO(pPoint, Point);

            points.append(point);
            return CUBIT_SUCCESS;
        }
    }


	DLIList<TopologyBridge*> tb_list;
    pSWPart->cubit_owner(entities, tb_list);
	CAST_LIST(tb_list, points, Point);

    IUnknown *pUnk;
	for (int i=entities.size(); i>0; i--)
	{
		pUnk = entities.get_and_step();
        pUnk->Release();
	}

	return CUBIT_SUCCESS;
}

CubitBoolean SWQueryEngine::bodies_overlap (BodySM *body_ptr_1, BodySM *body_ptr_2 ) const
{
  assert(false);
  return CUBIT_FALSE;
}

CubitStatus SWQueryEngine::translate( BodySM* body, const CubitVector& offset )
{
  assert(false);
  return CUBIT_FAILURE;
}

CubitStatus SWQueryEngine::rotate   ( BodySM* body, const CubitVector& axis, double angle )
{
  assert(false);
  return CUBIT_FAILURE;
}

CubitStatus SWQueryEngine::scale    ( BodySM* body, double factor )
{
  assert(false);
  return CUBIT_FAILURE;
}

CubitStatus SWQueryEngine::scale    ( BodySM* body, const CubitVector& factors )
{
  assert(false);
  return CUBIT_FAILURE;
}

CubitStatus SWQueryEngine::reflect  ( BodySM* body, const CubitVector& axis )
{
  assert(false);
  return CUBIT_FAILURE;
}

CubitStatus SWQueryEngine::restore_transform( BodySM* body )
{
  assert(false);
  return CUBIT_FAILURE;
}

CubitStatus SWQueryEngine::translate( GeometryEntity* ent, const CubitVector& offset )
{
  assert(false);
  return CUBIT_FAILURE;
}

CubitStatus SWQueryEngine::rotate   ( GeometryEntity* ent, const CubitVector& axis, double degrees )
{
  assert(false);
  return CUBIT_FAILURE;
}

CubitStatus SWQueryEngine::scale    ( GeometryEntity* ent, double factor )
{
  assert(false);
  return CUBIT_FAILURE;
}

CubitStatus SWQueryEngine::scale    ( GeometryEntity* ent, const CubitVector& factors )
{
  assert(false);
  return CUBIT_FAILURE;
}

CubitStatus SWQueryEngine::reflect  ( GeometryEntity* ent, const CubitVector& axis )
{
  assert(false);
  return CUBIT_FAILURE;
}

const char* SWQueryEngine::modeler_type()
{
  return "SolidWorks";
}

int SWQueryEngine::get_major_version()
{
  assert(false);
  return 0;
}

int SWQueryEngine::get_minor_version()
{
  assert(false);
  return 0;
}

int SWQueryEngine::get_subminor_version()
{
  assert(false);
  return 0;
}

int SWQueryEngine::get_allint_version()
{
  assert(false);
  return 0;
}

CubitString SWQueryEngine::get_engine_version_string()
{
  CubitString number_string = CubitString(get_major_version());
  number_string += ".";
  number_string += CubitString(get_minor_version());
  number_string += ".";
  number_string += CubitString(get_subminor_version());
  CubitString version_string = "ACIS Version ";
  version_string += number_string;
  return version_string;
}

CubitStatus SWQueryEngine::set_export_version(const int major, const int minor)
{
  PRINT_ERROR("Unable to change export version for SolidWorks\n");
  return CUBIT_FAILURE;
}

CubitStatus SWQueryEngine::set_export_allint_version(const int version)
{
  PRINT_ERROR("Unable to change export version for SolidWorks\n");
  return CUBIT_FAILURE;
}

CubitStatus SWQueryEngine::list_engine_versions(CubitString &versions)
{
  assert(false);
  return CUBIT_FAILURE;
}

void SWQueryEngine::reset()
{
  assert(false);
}
