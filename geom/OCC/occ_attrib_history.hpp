//-------------------------------------------------------------------------
// Filename      : occ_attrib_history.hpp
//
// Purpose       : Attributes needed to track subdivisions, 
//                 merges, copy, and geometry-changing events
//
//-------------------------------------------------------------------------

#ifndef OCC_ATTRIB_HISTORY_HPP
#define OCC_ATTRIB_HISTORY_HPP

// ********** BEGIN OCC INCLUDES             **********
// ********** END OCC INCLUDES               **********

// ********** BEGIN CUBIT INCLUDES            **********
#include <set>
// ********** END CUBIT INCLUDES              **********

// ********** BEGIN MACRO DEFINITIONS         **********
extern int OCC_ATTRIB_HISTORY_TYPE;
#define OCC_ATTRIB_HISTORY_LEVEL (ATTRIB_SNL_LEVEL + 1)
// ********** END MACRO DEFINITIONS           **********

// ********** BEGIN FORWARD DECLARATIONS      **********

// ********** END FORWARD DECLARATIONS        **********

class OCCHistory;
class TopoDS_Shape;

class OCC_ATTRIB_HISTORY 
{
public:
  
  OCC_ATTRIB_HISTORY( TopoDS_Shape* entity = NULL, 
                      OCCHistory* occ_history = NULL ); 

  static OCC_ATTRIB_HISTORY* get_history_attrib( TopoDS_Shape *occ_entity,
                                            bool create_if_necessary = false,
                                            OCCHistory *occ_history = NULL ); 
  
  //-------------------------
  void split_owner( TopoDS_Shape *entity);
  
  void merge_owner( TopoDS_Shape *entity, CubitBoolean delete_this);
  
  void trans_owner(  );

  void to_tolerant_owner( TopoDS_Shape *tol_ent );

  void copy_owner( TopoDS_Shape *copy_ent );

  void replace_owner( TopoDS_Shape *other_entity, 
                      CubitBoolean replace_owner );

  void lop_change_owner();

  void replace_owner_geometry( TopoDS_Shape *new_geom );

  void reverse_owner(); 

  //-------------------------

  std::set<int> get_tracking_ids();
  void add_tracking_id( int id ); 
  void add_occ_history( OCCHistory *occ_history );
  
  static void remove_all_attribs();

  //OCC_ATTRIB_FUNCTIONS(OCC_ATTRIB_HISTORY, NONE)
    
  private:
    //the history object that is being added to
    OCCHistory *occHistory; 

    //when an entity gets copied, the resultant copy entity needs to know 
    //the entity from whence it was copied.  This ATTRIB_HISTORY pointer 
    //facilitates mapping the original TopoDS_Shape to the resultant copied 
    //TopoDS_Shape 
    OCC_ATTRIB_HISTORY *fromHistoryAttrib;
    
    std::set<int> trackingIds;
    static std::set<OCC_ATTRIB_HISTORY*> allHistoryAttribs; 
    
    static int num_attribs;
    static bool addOrRemoveFromList; 
};

#endif


