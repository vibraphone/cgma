//-------------------------------------------------------------------------
// Filename      : PointSW.cpp
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

#include "stdafx.h"

//#include <amapp.h>
#import "C:\Program Files\SolidWorks\sldworks.tlb" raw_interfaces_only, raw_native_types, no_namespace, named_guids     //Import the SolidWorks type library

#include <assert.h>

#include "PointSW.hpp"
#include "SWQueryEngine.hpp"
#include "SWPart.hpp"


extern HRESULT ExtractVARIANTArrayData(VARIANT &v, long lArraySize, double *dArray);


//-------------------------------------------------------------------------
// Purpose       : The constructor with a pointer to the associated VERTEX . 
//
// Special Notes :
//
// Creator       : Joel Kopp, Xuechen Liu
//
// Creation Date : 8/22/00
//-------------------------------------------------------------------------
PointSW::PointSW(SWPart *pPart)
{
    m_pClosedEdge = NULL;
	m_pSWVertex = NULL;
	m_pSWPart = pPart;
}

//-------------------------------------------------------------------------
// Purpose       : The destructor. 
//
// Special Notes :
//
// Creator       : Joel Kopp, Raikanta Sahu
//
// Creation Date : 8/22/00
//-------------------------------------------------------------------------
PointSW::~PointSW() 
{
    if (m_pSWVertex)
        set_VERTEX_ptr(NULL);

    if (m_pClosedEdge)
        set_closed_edge_ptr(NULL);
}

//-------------------------------------------------------------------------
// Purpose       : This function adds a VERTEX pointer to the list of the
//                 VERTEX pointers associated with this object.  
//
// Special Notes :
//
// Creator       : Joel Kopp, Xuechen Liu
//
// Creation Date : 8/22/00
//-------------------------------------------------------------------------
void PointSW::set_VERTEX_ptr(IVertex *vertex)
{
    // this PointSW should only refer to one SolidWorks entity -- either a SW vertex
    // or a closed SW edge
    assert(NULL == m_pClosedEdge);

	if (vertex == m_pSWVertex)
		return;

	if (m_pSWVertex)
	{
		m_pSWPart->remove_cubit_owner(m_pSWVertex);
        m_pSWVertex->Release();
	}

	m_pSWVertex = vertex;

	if (m_pSWVertex)
	{
		m_pSWPart->set_cubit_owner(m_pSWVertex, this);
        m_pSWVertex->AddRef();
	}
}


IVertex *PointSW::get_VERTEX_ptr() const
{
	return m_pSWVertex;
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
GeometryQueryEngine* PointSW::get_geometry_query_engine() const
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
void PointSW::append_simple_attribute_virt(CubitSimpleAttrib* csattrib_ptr)
{
    if (m_pSWVertex)
        m_pSWPart->append_simple_attribute_virt(m_pSWVertex, csattrib_ptr);
    else
        assert(false); // attrib on pseudo vertex not supported
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
void PointSW::remove_simple_attribute_virt(CubitSimpleAttrib* csattrib_ptr)
{
    if (m_pSWVertex)
        m_pSWPart->remove_simple_attribute_virt(m_pSWVertex, csattrib_ptr);
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
// Creation Date : 8/22/00
//-------------------------------------------------------------------------
void PointSW::remove_all_simple_attribute_virt()
{
    if (m_pSWVertex)
        m_pSWPart->remove_all_simple_attribute_virt(m_pSWVertex);
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
// Creation Date : 8/22/00
//-------------------------------------------------------------------------
CubitStatus PointSW::get_simple_attribute(DLIList<CubitSimpleAttrib*>&
                                            cubit_simple_attrib_list)
{
    if (m_pSWVertex)
        return m_pSWPart->get_simple_attribute(m_pSWVertex, cubit_simple_attrib_list);
    else
        return CUBIT_FAILURE; // attrib on pseudo vertex not supported
}
CubitStatus PointSW::get_simple_attribute(const CubitString& name,
                                   DLIList<CubitSimpleAttrib*>& list)
{
  if (m_pSWVertex)
    return CUBIT_FAILURE;
    
 return m_pSWPart->get_simple_attribute(m_pSWVertex, name, list); 
}

//-------------------------------------------------------------------------
// Purpose       : Returns the coordinates of this Point. 
//
// Special Notes :
//
// Creator       : Joel Kopp, Malcolm J. Panthaki
//
// Creation Date : 8/22/00
//-------------------------------------------------------------------------
CubitVector PointSW::coordinates() const
{
    HRESULT hr = NOERROR;
	double point[3];

	// Get the first SW VERTEX associated with this PointSW object.
	IVertex *vertex = get_VERTEX_ptr();

	// Get the coordinates and return them
	if (vertex != NULL)
	{
		// Get the coordinates of the VERTEX
        VARIANT vPoint;
		hr = vertex->GetPoint ( &vPoint );

        double dTemp[3];
        hr = ExtractVARIANTArrayData(vPoint, 3, dTemp);
        hr = VariantClear(&vPoint);

        memcpy(point, dTemp, 3*sizeof(double));
	}
	else
	{
        // for a closed edge return the start point of the edge as the vertex coords
        if (m_pClosedEdge)
        {
            VARIANT vParams;
            hr = m_pClosedEdge->GetCurveParams2(&vParams);

            double dCurveParams[11];
            hr = ExtractVARIANTArrayData(vParams, 11, dCurveParams);
            hr = VariantClear(&vParams);

            memcpy(point, dCurveParams, 3*sizeof(double));
        }
        else
        {
            assert(false);
        }
	}

    CubitVector partCoords(point[0], point[1], point[2]);
    CubitVector coords;
    m_pSWPart->transformPartToAssembly(partCoords, coords);

	return coords;
}

//-------------------------------------------------------------------------
// Purpose       : Get the bounding box of the object.
//
// Special Notes :
//
// Creator       : Joel Kopp, Raikanta Sahu
//
// Creation Date : 8/22/00
//-------------------------------------------------------------------------
CubitBox PointSW::bounding_box() const 
{
	CubitVector temp_vector = this->coordinates();
	CubitBox temp_box(temp_vector);
	return temp_box;
}

TopologyEntity* PointSW::unmerge( DLIList<RefVolume*> )
{
    assert(false);
	PRINT_ERROR("Unmerging of SW-based entities is currently disabled.\n");
	return NULL;
}

CubitStatus PointSW::move(CubitVector &delta )
{
    assert(false);
	PRINT_ERROR("Translation of SW-based entities is currently disabled.\n");
	return CUBIT_FAILURE;
}

void PointSW::bodysms(DLIList<BodySM*> &bodies) 
{
	if (m_pSWVertex)
		SWQueryEngine::instance()->bodysms(m_pSWVertex, bodies);
}

void PointSW::lumps(DLIList<Lump*> &lumps)
{
	if (m_pSWVertex)
		SWQueryEngine::instance()->lumps(m_pSWPart, m_pSWVertex, lumps);
}

void PointSW::shellsms(DLIList<ShellSM*> &shellsms)
{
	if (m_pSWVertex)
		SWQueryEngine::instance()->shellsms(m_pSWPart, m_pSWVertex, shellsms);
}

void PointSW::surfaces(DLIList<Surface*> &surfaces)
{
	if (m_pSWVertex)
		SWQueryEngine::instance()->surfaces(m_pSWPart, m_pSWVertex, surfaces);
}

void PointSW::loopsms(DLIList<LoopSM*> &loopsms)
{
	if (m_pSWVertex)
		SWQueryEngine::instance()->loopsms(m_pSWPart, m_pSWVertex, loopsms);
}

void PointSW::curves(DLIList<Curve*> &curves)
{
	if (m_pSWVertex)
		SWQueryEngine::instance()->curves(m_pSWPart, m_pSWVertex, curves);
}

void PointSW::coedgesms(DLIList<CoEdgeSM*> &coedgesms)
{
	if (m_pSWVertex)
		SWQueryEngine::instance()->coedgesms(m_pSWPart, m_pSWVertex, coedgesms);
}

void PointSW::points(DLIList<Point*> &points)
{
	if (m_pSWVertex)
        SWQueryEngine::instance()->points(m_pSWPart, m_pSWVertex, points);
}

void PointSW::set_closed_edge_ptr(IEdge *edge)
{
    // the vertex should not reference a closed edge if it already reference a
    // SolidWorks vertex
    if (m_pSWVertex)
    {
        assert(false);
        return;
    }

	if (edge == m_pClosedEdge)
		return;

	if (m_pClosedEdge)
	{
		m_pSWPart->removePseudoVertex(m_pClosedEdge);
        m_pClosedEdge->Release();
	}

	m_pClosedEdge = edge;

	if (m_pClosedEdge)
	{
		m_pSWPart->addPseudoVertex(this, m_pClosedEdge);
        m_pClosedEdge->AddRef();
	}
}

void PointSW::get_parents_virt(DLIList<TopologyBridge*> &parents )
{
    IEdge         *pSWEdge = NULL;
    IEnumEdges    *pEdgeEnum = NULL;
    long nFetched = 0;

    if (m_pSWVertex)
    {
      HRESULT hr = m_pSWVertex->EnumEdges(&pEdgeEnum);
      if (SUCCEEDED(hr) && pEdgeEnum)
      {
          while (S_OK == pEdgeEnum->Next(1, &pSWEdge, &nFetched))
          {
            TopologyBridge *owner = m_pSWPart->cubit_owner(pSWEdge);
            if (owner)
              parents.append(owner);
          }
          pEdgeEnum->Release();
      }
    }
    else if (m_pClosedEdge)
    {
      TopologyBridge *owner = m_pSWPart->cubit_owner(m_pClosedEdge);
      if (owner)
        parents.append(owner);
    }
}

void PointSW::get_children_virt(DLIList<TopologyBridge*> &children )
{
  // no children to add
}

