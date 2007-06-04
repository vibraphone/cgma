//-------------------------------------------------------------------------
// Filename      : LumpSW.hpp
//
// Purpose       : 
//
// Special Notes :
//
// Creator       : Joel Kopp, Xuechen Liu
//
// Creation Date : 8/23/00
//
// Owner         : Byron Hanks
//-------------------------------------------------------------------------

#ifndef LUMP_SW_HPP
#define LUMP_SW_HPP

#include "CubitDefines.h"
#include "CubitEntity.hpp"
#include "LumpSM.hpp"
#include "SWPart.hpp"

struct IBody2;
class TopologyEntity;

class LumpSW : public LumpSM
{
public:
  
  LumpSW(SWPart *pPart);
  
  virtual ~LumpSW();
    //- The destructor
  
  IBody2 *get_BODY_ptr() const;
  void set_BODY_ptr(IBody2 *lump);
    // Get/set the LUMP associated with this object.

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
  
  virtual CubitBox bounding_box() const ;
  
//  virtual GeometryQueryEngine* 
//  get_geometry_query_engine() const;
    //R GeometryQueryEngine*
    //R- A pointer to the geometric modeling engine associated with
    //R- the object.
    //- This function returns a pointer to the geometric modeling engine
    //- associated with the object.
  
  virtual CubitStatus merge( GeometryEntity* /*GEPtr*/)
    {
      PRINT_ERROR("BUG: In LumpSW::merge\n"
                  "     This function should not be called at all\n"
                  "  This is a Bug -- please report it!\n");
      assert(0);
      return CUBIT_FAILURE;
    }
    //R CubitStatus
    //R- CUBIT_SUCCESS/CUBIT_FAILURE
    //I GEPtr
    //I- A pointer to a GeometryEntity object which *should* be of type
    //I- LumpSW.
    //- This function merges "this" LumpSW object with the input GEPtr
    //- object.  Checks the type of the input object to make sure it is
    //- another LumpSW object.
  
  virtual TopologyEntity* unmerge(DLIList<RefVolume*>)
    {
      PRINT_ERROR( "BUG: In LumpSW::unmerge\n"
                   "     This function should not be called\n"
                   "  This is a Bug -- please report it!\n" );
      assert(0);
      return (TopologyEntity*)NULL;
    }
  
  virtual double measure();
    //R double
    //R- The numeric value of the measure (its units depend on the dimension
    //R- of the RefEntity being "measured")
    //- A generic geometric extent function.
    //- Returns volume for Lump, area for Surface, length for Curve and 
    //- 1.0 for Point
  
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

  virtual CubitStatus mass_properties( CubitVector &centroid, double &volume );

protected: 
  
private:
	IBody2 *m_pSWBody; // SolidWorks doesn't have a lump object, so just store the body pointer
    SWPart *m_pSWPart;
} ;

#endif

