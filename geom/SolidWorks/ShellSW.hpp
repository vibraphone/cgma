//-------------------------------------------------------------------------
// Filename      : ShellSW.hpp
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

#ifndef SHELL_SW_HPP
#define SHELL_SW_HPP


#include "CubitDefines.h"
#include "CubitEntity.hpp"
#include "ShellSM.hpp"
#include "SWPart.hpp"

template <class X> class DLIList;


class ShellSW : public ShellSM
{
public:
  
  ShellSW(SWPart *pPart);
  
  virtual ~ShellSW() ;
    //- Destructor.
  
  IBody2 *get_BODY_ptr() const;
    //R LPBODY
    //R- Pointer to an SW body
    // Returns a pointer to the SW body associated with this ShellSW
    // object.

  void set_BODY_ptr(IBody2 *shell);
    //I shellDisp - The SHELL (body) to be associated with this object.

  virtual GeometryQueryEngine* get_geometry_query_engine() const;
    //R GeometryQueryEngine*
    //R- A pointer to the geometric modeling engine associated with
    //R- the object.
    //- This function returns a pointer to the geometric modeling engine
    //- associated with the object.
  
  virtual void append_simple_attribute_virt(CubitSimpleAttrib*);
    //R void
    //I 
    //I- 
    //I- that is to be appended to this OSME object.
    //- The purpose of this function is to append a 
    //- attribute to the OSME. The  is attached to each of the 
    //- underlying solid model entities this one points to.
  
  virtual void remove_simple_attribute_virt(CubitSimpleAttrib*);
    //R void
    //I CubitSimpleAttrib*
    //I- A reference to a CubitSimpleAttrib object which is the object
    //I- that is to be removed to this OSME object.
    //- The purpose of this function is to remove a simple
    //- attribute from the OSME. The attribute is attached to each of the
    //- underlying solid model entities this one points to.
  
  virtual void remove_all_simple_attribute_virt();
    //R void
    //I-
    //- The purpose of this function is to remove all simple
    //- attributes from the OSME. 
  
  virtual CubitStatus get_simple_attribute(const CubitString& name,
                                           DLIList<CubitSimpleAttrib*>&);
  virtual CubitStatus get_simple_attribute(DLIList<CubitSimpleAttrib*>&);
    //R CubitSimpleAttrib*
    //R- the returned cubit simple attribute.
    //- The purpose of this function is to get the attributes
    //- of the geometry entity. The name is attached to the underlying solid
    //- model entity(ies) this one points to.
    //- MJP Note:
    //- This is the code that implements the requirement that names
    //- of VGI Entities propagate across solid model boolean
    //- operations.  The success of this relies, of course, on the underlying
    //- solid modeler being able to propagate attributes across
    //- such operations on its entities. If it cannot, then "names"
    //- of VGI entities will not propagate.
  
  void bodysms(DLIList<BodySM*> &bodies);
  void lumps(DLIList<Lump*> &lumps);
  void shellsms(DLIList<ShellSM*> &shellsms);
  void surfaces(DLIList<Surface*> &surfaces);
  void loopsms(DLIList<LoopSM*> &loopsms);
  void curves(DLIList<Curve*> &curves);
  void coedgesms(DLIList<CoEdgeSM*> &coedgesms);
  void points(DLIList<Point*> &points);
    //- topology traversal of TB's; need to implement at this level 'cuz
    //- don't know how many SW entities per TB entity


  virtual void get_parents_virt(DLIList<TopologyBridge*> &parents );
  virtual void get_children_virt(DLIList<TopologyBridge*> &children );


protected: 

private:
	IBody2 *m_pSWBody;
    SWPart *m_pSWPart;
};


#endif

