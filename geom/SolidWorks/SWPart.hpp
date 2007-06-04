#ifndef __SWPart__
#define __SWPart__

#pragma warning(disable:4786)
#include <map>
#include <vector>

//#include "amapp.h"
#import "C:\Program Files\SolidWorks\sldworks.tlb" raw_interfaces_only, raw_native_types, no_namespace, named_guids     //Import the SolidWorks type library
#include "DLIList.hpp"
#include "CubitTransformMatrix.hpp"
#include "CubitString.hpp"
#include "CubitSimpleAttrib.hpp"
#include "SWRefKey.h"

class SWImport;
class LumpSM;
class BodySM;
class ShellSM;
class TopologyBridge;
class CubitSimpleAttrib;

typedef std::map<IFace2*, TopologyBridge*> FACEMAP;
typedef std::map<ILoop2*, TopologyBridge*> LOOPMAP;
typedef std::map<ICoEdge*, TopologyBridge*> COEDGEMAP;
typedef std::map<IEdge*, TopologyBridge*> EDGEMAP;
typedef std::map<IVertex*, TopologyBridge*> VERTEXMAP;
typedef std::multimap<IUnknown *, CubitSimpleAttrib*> ATTRIBMAP;
typedef std::map<IUnknown *, SWRefKey> KEYMAP;

class SWPart
{
public:
    SWPart();
    ~SWPart();

//    void setBody(BodySM *pBody);
    BodySM *getBody();
//    void setLump(LumpSM *pLump);
//    LumpSM *getLump();
//    void setShell(ShellSM *pShell);
//    ShellSM* getShell();

    // functions for managing cubit owner maps
    TopologyBridge* cubit_owner(IFace2 *entity);
    TopologyBridge* cubit_owner(ILoop2 *entity);
    TopologyBridge* cubit_owner(ICoEdge *entity);
    TopologyBridge* cubit_owner(IEdge *entity);
    TopologyBridge* cubit_owner(IVertex *entity);
    TopologyBridge* cubit_owner(IUnknown *entity);
    void cubit_owner(DLIList<IUnknown *> &entity_list,
                     DLIList<TopologyBridge*> &tb_list);

    void set_cubit_owner(IUnknown *entity, TopologyBridge *cubit_entity);
    void remove_cubit_owner(IUnknown *entity);


    void lumps(DLIList<TopologyBridge*> &tb_list);
    void shells(DLIList<TopologyBridge*> &tb_list);

    // map for pseudo vertices - cubit vertex mapped to a closed SW edge
    void addPseudoVertex(TopologyBridge* pPseudoVertex, IEdge *pClosedSWEdge);
    TopologyBridge *findPseudoVertex(IEdge *pClosedSWEdge);
    void removePseudoVertex(IEdge *pClosedSWEdge);

    // manage and use part instance transforms
    HRESULT setTransforms(IComponent2 *pSWComponent);
    CubitTransformMatrix const & getPartToAssemblyTransform();
    CubitTransformMatrix const & getAssemblyToPartTransform();

    void transformPartToAssembly(CubitVector &partVect, CubitVector &assemblyVect);
    void transformAssemblyToPart(CubitVector &assemblyVect, CubitVector &partVect);

    // manage attributes associated with this part
    void append_simple_attribute_virt(IUnknown *pSWEntity, CubitSimpleAttrib* attrib_ptr);
    void remove_simple_attribute_virt(IUnknown *pSWEntity, CubitSimpleAttrib* attrib_ptr);
    void remove_all_simple_attribute_virt(IUnknown *pSWEntity);
    CubitStatus get_simple_attribute(IUnknown *pSWEntity,
                                     DLIList<CubitSimpleAttrib*>& cubit_simple_attrib_list);
    CubitStatus get_simple_attribute(IUnknown *pSWEntity, const CubitString& name,
                                     DLIList<CubitSimpleAttrib*>& cubit_simple_attrib_list);

    // special handling for body and lump attributes
    void append_simple_attribute_virt(TopologyBridge* pBridge, CubitSimpleAttrib* attrib_ptr);
    void remove_simple_attribute_virt(TopologyBridge* pBridge, CubitSimpleAttrib* attrib_ptr);
    void remove_all_simple_attribute_virt(TopologyBridge* pBridge);
    CubitStatus get_simple_attribute(TopologyBridge* pBridge, const CubitString& name,
                                     DLIList<CubitSimpleAttrib*>& cubit_simple_attrib_list);
    CubitStatus get_simple_attribute(TopologyBridge* pBridge,
                                     DLIList<CubitSimpleAttrib*>& cubit_simple_attrib_list);

    void updateAttribs();
    void setImport(SWImport *pImport);

    void setName(IComponent2 *pSWComp);
    CubitString getName();

private:

    SWImport *m_pImport;

	CubitString m_csPartName; // part instance name for assembly

// Ownership maps and data - cubit owner data and corresponding SW entity

    // special case for body, lump, and shell
    // SolidWorks has a body entity, but no lump.  Since there is only one body and
    // one lump per part, let both the body and the lump own the SolidWorks body.
    IBody2 *m_pSWBody;
    TopologyBridge* m_body;
    TopologyBridge* m_lump;

    // SolidWorks has no shell -- for now I am assuming a single shell per part, so
    // the shell also maps to the SolidWorks body
    TopologyBridge* m_shell;

    // owner map for faces
	FACEMAP m_FaceMap;
	FACEMAP::iterator m_FaceMapIter;

    // owner map for loops
	LOOPMAP m_LoopMap;
	LOOPMAP::iterator m_LoopMapIter;

    // owner map for coedges
	COEDGEMAP m_CoEdgeMap;
	COEDGEMAP::iterator m_CoEdgeMapIter;

    //owner map for edges
	EDGEMAP m_EdgeMap;
	EDGEMAP::iterator m_EdgeMapIter;

    // owner map for vertices
	VERTEXMAP m_VertexMap;
	VERTEXMAP::iterator m_VertexMapIter;

    // special case map for vertices on closed edges - SolidWorks does not have a
    // vertex on a closed edge, so I create a cubit vertex that points to a SolidWorks
    // edge
	EDGEMAP m_PseudoVertexMap;

// Associativity maps

    // special case lists for body and volume - since SolidWorks only has a body pointer,
    // store body and volume attributes in separate lists.
    std::vector<CubitSimpleAttrib*> m_BodyAttribs;
    std::vector<CubitSimpleAttrib*> m_LumpAttribs;

	// multimap for associativity - maps multiple attributes to an LPUNKNOWN (SW entity)
    // this map holds attributes for all faces, edges, and vertices
	ATTRIBMAP m_AttribMap;
	ATTRIBMAP::iterator m_AttribMapIter;

    // multimap for any dangling attributes - could not be mapped to an entity in the
    // current brep
    ATTRIBMAP m_DanglingAttribMap;


    KEYMAP m_KeyMap;
    KEYMAP::iterator m_KeyMapIter;

// transformations
    CubitTransformMatrix m_PartToAsm; // part to assembly transformation
    CubitTransformMatrix m_AsmToPart; // assembly to part transformation
};
#endif // __SWPart__
