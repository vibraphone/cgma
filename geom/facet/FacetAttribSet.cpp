//-------------------------------------------------------------------------
// Filename      : FacetAttribSet.cpp
//
// Purpose       : Common attribute functionality for MBG
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 03/01/03
//-------------------------------------------------------------------------

#include "FacetAttribSet.hpp"
#include "FacetAttrib.hpp"
#include "CubitSimpleAttrib.hpp"
#include "CubitFileIOWrapper.hpp"

void FacetAttribSet::append_attribute( const CubitSimpleAttrib& csa )
{
  FacetAttrib* new_attrib = new FacetAttrib(csa);
  new_attrib->listNext = listHead;
  listHead = new_attrib;
}

void FacetAttribSet::remove_attribute( const CubitSimpleAttrib& csa )
{
  if( !listHead )
    return;
    
  FacetAttrib* attrib = 0;
  if ( listHead->equals(csa) )
  {
    attrib = listHead;
    listHead = listHead->listNext;
    delete attrib;
    return;
  }
  
  for ( FacetAttrib* prev = listHead; prev->listNext; prev = prev->listNext )
  {
    if( prev->listNext->equals(csa) )
    {
      attrib = prev->listNext;
      prev->listNext = attrib->listNext;
      delete attrib;
      return;
    }
  }
}

void FacetAttribSet::remove_all_attributes()
{
  while( listHead )
  {
    FacetAttrib* dead = listHead;
    listHead = dead->listNext;
    delete dead;
  }
}

CubitStatus FacetAttribSet::get_attributes( DLIList<CubitSimpleAttrib>& list ) const
{
  for( FacetAttrib* attrib = listHead; attrib; attrib = attrib->listNext )
    list.append( attrib->get_CSA() );
  return CUBIT_SUCCESS;
}

CubitStatus FacetAttribSet::get_attributes( const CubitString& name,
                                    DLIList<CubitSimpleAttrib>& list ) const
{
  for( FacetAttrib* attrib = listHead; attrib; attrib = attrib->listNext )
    if( attrib->name() == name )
      list.append( attrib->get_CSA() );
  return CUBIT_SUCCESS;
}

CubitStatus FacetAttribSet::save_attributes( FILE* file_ptr ) const
{
  FacetAttrib *curr_attrib;
  CubitStatus status = CUBIT_SUCCESS;
  
  //save # attribs
  unsigned int size = attribute_count();
  NCubitFile::CIOWrapper wrapper( file_ptr );
  wrapper.Write( &size, 1 ); 

  //save each attrib
  for( curr_attrib = listHead; curr_attrib; curr_attrib = curr_attrib->listNext )
    if( !curr_attrib->save(file_ptr) )
      status = CUBIT_FAILURE;

  return status;
}
  
CubitStatus FacetAttribSet::restore_attributes( FILE* file_ptr, unsigned endian )
{
  FacetAttrib *curr_attrib;
  
  //Read # attribs
  unsigned int size;
  NCubitFile::CIOWrapper wrapper( endian, file_ptr );
  wrapper.Read( &size, 1 ); 

  for (unsigned i = 0; i < size; i++)
  {
    curr_attrib = FacetAttrib::restore( file_ptr, endian);  
    if (!curr_attrib)
    {
        // file corrupt?  don't try to read any more
      return CUBIT_FAILURE;
    }
    
    curr_attrib->listNext = listHead;
    listHead = curr_attrib;
  }

  return CUBIT_SUCCESS;
}


int FacetAttribSet::attribute_count() const
{
  int count = 0;
  for( FacetAttrib* attrib = listHead; attrib; attrib = attrib->listNext )
    count++;
  return count;
}

