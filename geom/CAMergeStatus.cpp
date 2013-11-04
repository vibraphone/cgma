//-------------------------------------------------------------------------
// Filename      : CAMergeStatus.cpp
//
// Purpose       : Save RefEntity::autoMergeStatus
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 06/06/02
//-------------------------------------------------------------------------

#include "RefEntity.hpp"
#include "CAMergeStatus.hpp"
#include "CADefines.hpp"

CubitAttrib* CAMergeStatus_creator(RefEntity* entity, const CubitSimpleAttrib &p_csa)
{
  return new CAMergeStatus(entity, p_csa);
}

CAMergeStatus::CAMergeStatus( RefEntity* owner, const CubitSimpleAttrib& csa )
  : CubitAttrib( owner )
{ 
  status = AUTO_MERGE_AUTO;
  if(!csa.isEmpty())
  {
    assert( csa.int_data_list().size() == 1 );
    int i = csa.int_data_list()[0];
    assert( i == 0 || i == 1 || i == 2 );
    status = (AutoMergeStatus)i;
  }
}

CAMergeStatus::~CAMergeStatus() 
{ }

CubitStatus CAMergeStatus::actuate() 
{
  if( hasActuated )
    return CUBIT_SUCCESS;
  
  if( !attribOwnerEntity )
    return CUBIT_FAILURE;
    
  DLIList<CubitAttrib*> att_list;
  attribOwnerEntity->find_cubit_attrib_type(CA_MERGE_PARTNER, att_list);
  if( att_list.size() )
    return CUBIT_FAILURE;
  
  assert( status == AUTO_MERGE_ON || status == AUTO_MERGE_OFF );  
  attribOwnerEntity->is_mergeable(status);
  
  deleteAttrib = CUBIT_FALSE;
  hasActuated = CUBIT_TRUE;

  return CUBIT_SUCCESS;
}

CubitStatus CAMergeStatus::update()
{
  if( hasUpdated )
    return CUBIT_SUCCESS;
    
  if( !attribOwnerEntity )
    return CUBIT_FAILURE;
  
  status = attribOwnerEntity->merge_status();
  if( status == AUTO_MERGE_AUTO )
    delete_attrib(CUBIT_TRUE);
  else
    hasUpdated = CUBIT_TRUE;
  
  return CUBIT_SUCCESS;
}

CubitStatus CAMergeStatus::reset()
{
  return CUBIT_SUCCESS;
}

CubitSimpleAttrib CAMergeStatus::cubit_simple_attrib()
{
  if( deleteAttrib )
    return CubitSimpleAttrib();
  
  assert( status != AUTO_MERGE_AUTO );
  CubitSimpleAttrib* result = 0;

  std::vector<CubitString> string_list;
  string_list.push_back(att_internal_name());
  
  int int_data = (int)status;
  std::vector<int> int_list;
  int_list.push_back( int_data );
  
  return CubitSimpleAttrib( &string_list, 0, &int_list );
}

int CAMergeStatus::int_attrib_type()
{
  return CA_MERGE_STATUS;
}

