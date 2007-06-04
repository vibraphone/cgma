#ifndef __SWIMPORT_HPP__
#define __SWIMPORT_HPP__

//#include "amapp.h"
#import "C:\Program Files\SolidWorks\sldworks.tlb" raw_interfaces_only, raw_native_types, no_namespace, named_guids     //Import the SolidWorks type library
#include "swconst.h"
#include "DLIList.hpp"

class SWPart;
class BodySM;
class Surface;
class Curve;
class Point;
class TopologyBridge;

class SWImport
{
public:
    SWImport();
    ~SWImport();

    CubitStatus import_solid_model(IModelDoc2 *pSWModelDoc,
                                   DLIList<TopologyBridge*>& imported_entities);

    IModelDoc2 *document();

private:
    void addPart(SWPart* &pSWPart);
    void addPart(IComponent2 *pSWComp, SWPart* &pSWPart);
    void updateComponentAndPartLists(IComponent2 *pSWComp, SWPart* &pNewSWPart);
    void findComponentPart(IComponent2 *pLeafComponent, DLIList<SWPart*> &oldPartList,
                           SWPart* &pSWPart);

    CubitStatus setDocument(IModelDoc2 *pSWDoc);

    BodySM *populate_topology_bridges(SWPart *pSWPart, IBody2 *pSWBody);
    Surface *populate_topology_bridges(SWPart *pSWPart, IFace2 *pSWFace);
    Curve *populate_topology_bridges(SWPart *pSWPart, IEdge *pSWEdge);
    Point *populate_topology_bridges(SWPart *pSWPart, IVertex *pSWVertex);

private:
    DLIList<SWPart*> m_SWPartList;
    DLIList<IComponent2 *> m_SWComponentList;

    enum swDocumentTypes_e m_SWDocType; // part/assembly

    IModelDoc2 *m_pSWModelDoc;
};
#endif // __SWIMPORT_HPP__