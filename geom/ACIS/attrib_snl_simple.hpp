//- Class:       attrib_snl_simple.hpp
//- Owner:       Greg Nielson
//- Description: The ACIS attribute class used to store attribute information
//-              on ACIS (.sat) files.
//- Checked by:
//- Version: $Id:

#if !defined( ATTRIB_SNL_SIMPLE_CLASS )
#define ATTRIB_SNL_SIMPLE_CLASS

#include "attrib_snl.hpp"
#include <string.h>
#include "DLIList.hpp"
#include "CubitString.hpp"
#include "AcisTypes.h"

class CubitSimpleAttrib;

extern int ATTRIB_SNL_SIMPLE_TYPE;
#define ATTRIB_SNL_SIMPLE_LEVEL (ATTRIB_SNL_LEVEL + 1)

class ATTRIB_SNL_SIMPLE: public ATTRIB_SNL {
  friend class BodyACIS;
  friend class CoEdgeACIS;
  friend class LoopACIS;
  friend class ShellACIS;
  
public:
  
  class AttribData {
    public:
    
      AttribData(CubitSimpleAttrib*);

      ~AttribData();


      inline void inc_use_count() { useCount++; }

      void dec_use_count();


      CubitSimpleAttrib* make_CSA() const;

      CubitBoolean equivalent( CubitSimpleAttrib* ) const;


      inline CubitString name() const
        { return numStrings ? stringData[0] : CubitString(); }
        
      void save() const;
      
      static AttribData* read_old_attrib( const char* name );
      static AttribData* read_new_attrib();
   
    private: 

      AttribData(int num_strings, CubitString* strings,
                 int num_ints, int* ints,
                 int num_reals, double* reals );
    
      int useCount;
      int numStrings;
      int numInts;
      int numReals;

      CubitString* stringData;
      int* intData;
      double* realData;
  
  };
  
private:
   AttribData *data;
    
public:

  inline ATTRIB_SNL_SIMPLE( ENTITY* owner = NULL) 
    : ATTRIB_SNL(owner), data(0) {}
  
  inline ATTRIB_SNL_SIMPLE( ENTITY* owner, const ATTRIB_SNL_SIMPLE& from)
    : ATTRIB_SNL(owner) { data = from.data; data->inc_use_count(); }

  inline ATTRIB_SNL_SIMPLE( ENTITY* owner, CubitSimpleAttrib* csa )
    : ATTRIB_SNL(owner) { data = new AttribData(csa); }
  
  inline CubitString attribute_name() const 
    { return data->name(); }                    
                     
  inline CubitSimpleAttrib* get_CSA() const
    { return data->make_CSA(); }
  
// Functions called to aid attribute migration during modeling
// (boolean) operations.

  virtual void split_owner( ENTITY * );
//- The owner of this attribute is about to be split in two - the
//- argument is the new piece.  

  virtual void merge_owner( ENTITY *, // "other entity"
                             logical  // deleting_owner
                          );
//- The owner of this attribute is about to be merged with the 
//- given entity.  The logical argument is TRUE if the owner is to 
//- be deleted in the merge.

  virtual void trans_owner( SPAtransf const&);
  
  //-callback to notify this attribute that it's owner is being copied.
  virtual void copy_owner( ENTITY *copy_ent );

  inline CubitBoolean equivalent(CubitSimpleAttrib *csa_ptr) const
    { return data->equivalent(csa_ptr); }
    //- return true if csa_ptr stores the same data as this attrib

  ATTRIB_FUNCTIONS( ATTRIB_SNL_SIMPLE, NONE )
};

#endif
