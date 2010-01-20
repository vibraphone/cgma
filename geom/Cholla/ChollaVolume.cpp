//- Class:       ChollaVolume
//- Description: Temporary class for constructing the facet-based geometry
//-              
//- Owner:       Steven J. Owen
//- Checked by:
//- Version:
#include "ChollaVolume.hpp"
#include "ChollaSurface.hpp"


//===============================================================================
//Function:  ChollaVolume (PUBLIC) (constructor)
//===============================================================================
ChollaVolume::ChollaVolume(int block_id)
{
  static int count = 100;
  id = count++;
  myVolume = NULL;
  blockId = block_id;
}

//===============================================================================
//Function:  ~ChollaVolume (PUBLIC) (destructor)
//===============================================================================
ChollaVolume::~ChollaVolume()
{
}



//EOF

