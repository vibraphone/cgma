/*******************************************************************************
    COPYRIGHT 2002 CATERPILLAR INC.  ALL RIGHTS RESERVED

    This program is the property of Caterpillar Inc., includes Caterpillar's
    confidential and trade secret information, and is maintained as an
    unpublished copyrighted work.  It is not to be copied or used by others
    except under license from Caterpillar.  This program is also protected as an
    unpublished work in accordance with the copyright act of 1976.  In the event
    of either inadvertent or deliberate publication, Caterpillar Inc. intends to
    maintain copyright protection for this work under the relevant copyright
    laws pertaining to published works.  The inclusion of a copyright notice
    hereon is precautionary only, and does not imply publication or disclosure.

 
 Filename      : CubitFileSimModel.cpp

 Purpose       : Implements the reading and writing functionality for a Simulation
                 model section of a Cubit (*.cub) format file.
           
 Special Notes :

 Creator       : Andrew Rout

 Creation Date : 02/20/10

 Owner         : Andrew Rout

*******************************************************************************/

#include "CubitFileSimModel.hpp"
#include "CubitFileIOWrapper.hpp"
#include "CubitFile.hpp"

using namespace NCubitFile;

///////////////////////////////////////////////////////////////////////////////
// CSimModel
///////////////////////////////////////////////////////////////////////////////

// The number of 32 bit words contained int each of the stored structures.
const UnsignedInt32 CSimModel::mintSizeOfSimModelHeader =
    sizeof(SCubitFileSimModelHeader) / sizeof(UnsignedInt32);
const UnsignedInt32 CSimModel::mintSizeOfSimModelHeader2 =
    sizeof(SCubitFileSimModelHeader2) / sizeof(UnsignedInt32);
const UnsignedInt32 CSimModel::mintSizeOfBCEntry =
    sizeof(SCubitFileBCEntry) / sizeof(UnsignedInt32);
const UnsignedInt32 CSimModel::mintSizeOfICEntry =
    sizeof(SCubitFileICEntry) / sizeof(UnsignedInt32);
const UnsignedInt32 CSimModel::mintSizeOfBCSetEntry =
    sizeof(SCubitFileBCSetEntry) / sizeof(UnsignedInt32);
const UnsignedInt32 CSimModel::mintSizeOfMaterialEntry =
    sizeof(SCubitFileMaterialEntry) / sizeof(UnsignedInt32);
const UnsignedInt32 CSimModel::mintSizeOfAmplitudeEntry =
    sizeof(SCubitFileAmplitudeEntry) / sizeof(UnsignedInt32);
const UnsignedInt32 CSimModel::mintSizeOfConstraintEntry =
    sizeof(SCubitFileConstraintEntry) / sizeof(UnsignedInt32);

CSimModel::CSimModel()
{
    mpReadFile = mpWriteFile = NULL;
    mpaBCs = NULL;
    mpaICs = NULL;
    mpaBCSets = NULL;
    mpaMaterials = NULL;
    mpaAmplitudes = NULL;
    mpaConstraints = NULL;

    mintSimModelOffset = 0;

    memset(&mSimModel, 0, sizeof(SCubitFileSimModelHeader));
    mSimModel.mintSimModelEndian = CCubitFile::mintNativeEndian;

    memset(&mSimModel2, 0, sizeof(SCubitFileSimModelHeader2));

    memset(&mBCSetBuff, 0, sizeof(SBCSetReturnBuffer));
    memset(&mMaterialBuff, 0, sizeof(SMaterialReturnBuffer));
    memset(&mConstraintBuff, 0, sizeof(SConstraintReturnBuffer));

}

CSimModel::~CSimModel()
{
    if(mpaBCs)
        delete [] mpaBCs;
    if(mpaICs)
        delete [] mpaICs;
    if(mpaBCSets)
        delete [] mpaBCSets;
    if(mpaMaterials)
        delete [] mpaMaterials;
    if(mpaAmplitudes)
        delete [] mpaAmplitudes;
    if(mpaConstraints)
        delete [] mpaConstraints;

    if(mBCSetBuff.mpaBCSetRestraintData)
        delete [] mBCSetBuff.mpaBCSetRestraintData;
    if(mBCSetBuff.mpaintRestraintMemberIDs)
        delete [] mBCSetBuff.mpaintRestraintMemberIDs;

    if(mBCSetBuff.mpaBCSetLoadData)
        delete [] mBCSetBuff.mpaBCSetLoadData;
    if(mBCSetBuff.mpaintLoadMemberIDs)
        delete [] mBCSetBuff.mpaintLoadMemberIDs;

    if(mBCSetBuff.mpaBCSetContactPairData)
        delete [] mBCSetBuff.mpaBCSetContactPairData;
    if(mBCSetBuff.mpaintContactPairMemberIDs)
        delete [] mBCSetBuff.mpaintContactPairMemberIDs;

    if(mMaterialBuff.mpaMaterialData)
        delete [] mMaterialBuff.mpaMaterialData;
    if(mMaterialBuff.mpadblData)
        delete [] mMaterialBuff.mpadblData;
    
    if(mConstraintBuff.mpaintIndependentIDs)
        delete [] mConstraintBuff.mpaintIndependentIDs;
    if(mConstraintBuff.mpaintDependentIDs)
        delete [] mConstraintBuff.mpaintDependentIDs;
    if(mConstraintBuff.mpaIndependentData)
        delete [] mConstraintBuff.mpaIndependentData;
    if(mConstraintBuff.mpaDependentData)
        delete [] mConstraintBuff.mpaDependentData;

}


///////////////////////////////////////////////////////////////////////////////
// Write Functions
///////////////////////////////////////////////////////////////////////////////

UnsignedInt32 CSimModel::InitWrite(FILE* xpFile,
                                   UnsignedInt32 xintBCCount,
                                   UnsignedInt32 xintICCount,
                                   UnsignedInt32 xintBCSetCount,
                                   UnsignedInt32 xintMaterialCount,
                                   UnsignedInt32 xintAmplitudeCount,
                                   UnsignedInt32 xintConstraintCount)
{
    if(mpWriteFile)  throw CCubitFile::eFileWriteError;

    mpWriteFile = xpFile;
    CIOWrapper lIO(mpWriteFile);

    // Write out the Sim model header to reserve its position and size in the
    // file.
    mintSimModelOffset = lIO.BeginWriteBlock();
    lIO.Write((UnsignedInt32*)&mSimModel, mintSizeOfSimModelHeader);
    mSimModel.mintSimModelLength = lIO.EndWriteBlock();

    mSimModel.mintFutureTableOffset = lIO.BeginWriteBlock(mintSimModelOffset);
    lIO.Write((UnsignedInt32*)&mSimModel2, mintSizeOfSimModelHeader2);
    mSimModel.mintSimModelLength += lIO.EndWriteBlock();

    mSimModel.mintBCCount = xintBCCount;
    if(xintBCCount) {
        // Create a geometry entity array for storing ownership statistics and
        // initially blank it.
        mpaBCs = new SCubitFileBCEntry[xintBCCount];
        memset(mpaBCs, 0, sizeof(SCubitFileBCEntry) * xintBCCount);

        // Write the blank geometry table to the file to reserve its position
        // and size in the file.
        mSimModel.mintBCTableOffset = lIO.BeginWriteBlock(mintSimModelOffset);
        lIO.Write((UnsignedInt32*)mpaBCs,
            xintBCCount * mintSizeOfBCEntry);
        mSimModel.mintSimModelLength += lIO.EndWriteBlock();
    }

    mSimModel.mintICCount = xintICCount;
    if(xintICCount) {
        // Create a group array for storing ownership statistics and
        // initially blank it.
        mpaICs = new SCubitFileICEntry[xintICCount];
        memset(mpaICs, 0, sizeof(SCubitFileICEntry) * xintICCount);

        // Write the blank group table to the file to reserve its position
        // and size in the file.
        mSimModel.mintICTableOffset = lIO.BeginWriteBlock(mintSimModelOffset);
        lIO.Write((UnsignedInt32*)mpaICs,
            xintICCount * mintSizeOfICEntry);
        mSimModel.mintSimModelLength += lIO.EndWriteBlock();
    }

    mSimModel.mintBCSetCount = xintBCSetCount;
    if(xintBCSetCount) {
        // Create a block array for storing ownership statistics and
        // initially blank it.
        mpaBCSets = new SCubitFileBCSetEntry[xintBCSetCount];
        memset(mpaBCSets, 0, sizeof(SCubitFileBCSetEntry) * xintBCSetCount);

        // Write the blank block table to the file to reserve its position
        // and size in the file.
        mSimModel.mintBCSetTableOffset = lIO.BeginWriteBlock(mintSimModelOffset);
        lIO.Write((UnsignedInt32*)mpaBCSets,
            xintBCSetCount * mintSizeOfBCSetEntry);
        mSimModel.mintSimModelLength += lIO.EndWriteBlock();
    }

    mSimModel.mintMaterialCount = xintMaterialCount;
    if(xintMaterialCount) {
        // Create a node set array for storing ownership statistics and
        // initially blank it.
        mpaMaterials = new SCubitFileMaterialEntry[xintMaterialCount];
        memset(mpaMaterials, 0, sizeof(SCubitFileMaterialEntry) * xintMaterialCount);

        // Write the blank geometry table to the file to reserve its position
        // and size in the file.
        mSimModel.mintMaterialTableOffset = lIO.BeginWriteBlock(mintSimModelOffset);
        lIO.Write((UnsignedInt32*)mpaMaterials,
            xintMaterialCount * mintSizeOfMaterialEntry);
        mSimModel.mintSimModelLength += lIO.EndWriteBlock();
    }

    mSimModel.mintAmplitudeCount = xintAmplitudeCount;
    if(xintAmplitudeCount) {
        // Create a SideSet array for storing ownership statistics and
        // initially blank it.
        mpaAmplitudes = new SCubitFileAmplitudeEntry[xintAmplitudeCount];
        memset(mpaAmplitudes, 0, sizeof(SCubitFileAmplitudeEntry) * xintAmplitudeCount);

        // Write the blank geometry table to the file to reserve its position
        // and size in the file.
        mSimModel.mintAmplitudeTableOffset = lIO.BeginWriteBlock(mintSimModelOffset);
        lIO.Write((UnsignedInt32*)mpaAmplitudes,
            xintAmplitudeCount * mintSizeOfAmplitudeEntry);
        mSimModel.mintSimModelLength += lIO.EndWriteBlock();
    }

    mSimModel2.mintConstraintCount = xintConstraintCount;
    if(xintConstraintCount) {
        // Create a SideSet array for storing ownership statistics and
        // initially blank it.
        mpaConstraints = new SCubitFileConstraintEntry[xintConstraintCount];
        memset(mpaConstraints, 0, sizeof(SCubitFileConstraintEntry) * xintConstraintCount);

        // Write the blank geometry table to the file to reserve its position
        // and size in the file.
        mSimModel2.mintConstraintTableOffset = lIO.BeginWriteBlock(mintSimModelOffset);
        lIO.Write((UnsignedInt32*)mpaConstraints,
            xintConstraintCount * mintSizeOfConstraintEntry);
        mSimModel.mintSimModelLength += lIO.EndWriteBlock();
    }

    return mintSimModelOffset;
}

void CSimModel::WriteBCSet(UnsignedInt32 xintIndex,
                          UnsignedInt32 xintBCSetID,
                          UnsignedInt32 xintBCSetUniqueID,
                          UnsignedInt32 xintBCSetAnalysisType,
                          UnsignedInt32 xintRestraintTypesCount,
                          UnsignedInt32 xintLoadTypesCount,
                          UnsignedInt32 xintContactPairTypesCount,
                          SBCSetData* xpaBCSetRestraintData,
                          SBCSetData* xpaBCSetLoadData,
                          SBCSetData* xpaBCSetContactPairData)
{
    if(!mpWriteFile)
        throw CCubitFile::eFileWriteError;
    if(xintIndex >= mSimModel.mintBCSetCount)
        throw CCubitFile::eNotFound;
    if(!mpaBCSets)
        throw CCubitFile::eOrderError;
    if(mpaBCSets[xintIndex].mintBCSetLength) // something that won't be zero
        throw CCubitFile::eDuplicateWrite;

    if(xintRestraintTypesCount)
    {  
        // An empty bcset is valid, but an incompletely defined one is not.
        if(!xpaBCSetRestraintData)
            throw CCubitFile::ePassedNullPointer;
        for(UnsignedInt32 lintRestraints = 0; lintRestraints < xintRestraintTypesCount; lintRestraints++) {
            if(xpaBCSetRestraintData[lintRestraints].mintMemberCount &&
                !xpaBCSetRestraintData[lintRestraints].mpaMemberIDs)
                throw CCubitFile::ePassedNullPointer;
        }
    }
    if(xintLoadTypesCount)
    {  
        // An empty bcset is valid, but an incompletely defined one is not.
        if(!xpaBCSetRestraintData)
            throw CCubitFile::ePassedNullPointer;
        for(UnsignedInt32 lintLoads = 0; lintLoads < xintLoadTypesCount; lintLoads++) {
            if(xpaBCSetLoadData[lintLoads].mintMemberCount &&
                !xpaBCSetLoadData[lintLoads].mpaMemberIDs)
                throw CCubitFile::ePassedNullPointer;
        }
    }
    if(xintContactPairTypesCount)
    {  
        // An empty bcset is valid, but an incompletely defined one is not.
        if(!xpaBCSetContactPairData)
            throw CCubitFile::ePassedNullPointer;
        for(UnsignedInt32 lintContactPairs = 0; lintContactPairs < xintContactPairTypesCount; lintContactPairs++) {
            if(xpaBCSetContactPairData[lintContactPairs].mintMemberCount &&
                !xpaBCSetContactPairData[lintContactPairs].mpaMemberIDs)
                throw CCubitFile::ePassedNullPointer;
        }
    }

    // BCSet header written out elsewhere.  Here, just populate the header struct
    mpaBCSets[xintIndex].mintBCSetID = xintBCSetID;
    mpaBCSets[xintIndex].mintBCSetUniqueID = xintBCSetUniqueID;
    mpaBCSets[xintIndex].mintBCSetAnalysisType = xintBCSetAnalysisType;
    mpaBCSets[xintIndex].mintRestraintTypesCount = xintRestraintTypesCount;
    mpaBCSets[xintIndex].mintLoadTypesCount = xintLoadTypesCount;
    mpaBCSets[xintIndex].mintContactPairTypesCount = xintContactPairTypesCount;
    
    CIOWrapper* lpIO = new CIOWrapper(mpWriteFile);

    // write restraints
    if(xintRestraintTypesCount)
    {        
        mpaBCSets[xintIndex].mintRestraintsOffset =
            lpIO->BeginWriteBlock(mintSimModelOffset);

        for(UnsignedInt32 lintBCSet = 0; lintBCSet < xintRestraintTypesCount; lintBCSet++)
        {
            if(!xpaBCSetRestraintData[lintBCSet].mintMemberCount)
            {
                mpaBCSets[xintIndex].mintRestraintTypesCount--;
                continue;
            }
            mpaBCSets[xintIndex].mintRestraintMembersCount +=
                xpaBCSetRestraintData[lintBCSet].mintMemberCount;
            
            lpIO->Write(&xpaBCSetRestraintData[lintBCSet].mintMemberType, 1);
            lpIO->Write(&xpaBCSetRestraintData[lintBCSet].mintMemberCount, 1);
            lpIO->Write(xpaBCSetRestraintData[lintBCSet].mpaMemberIDs,
                xpaBCSetRestraintData[lintBCSet].mintMemberCount);
        }

        mpaBCSets[xintIndex].mintBCSetLength += lpIO->EndWriteBlock();
    }
    else
    {
        // BCSet has no restraints in the file.
        mpaBCSets[xintIndex].mintRestraintsOffset = 0;
    }

    //write loads
    if(xintLoadTypesCount)
    {        
        mpaBCSets[xintIndex].mintLoadsOffset =
            lpIO->BeginWriteBlock(mintSimModelOffset);

        for(UnsignedInt32 lintBCSet = 0; lintBCSet < xintLoadTypesCount; lintBCSet++)
        {
            if(!xpaBCSetLoadData[lintBCSet].mintMemberCount)
            {
                mpaBCSets[xintIndex].mintLoadTypesCount--;
                continue;
            }
            mpaBCSets[xintIndex].mintLoadMembersCount +=
                xpaBCSetLoadData[lintBCSet].mintMemberCount;
            
            lpIO->Write(&xpaBCSetLoadData[lintBCSet].mintMemberType, 1);
            lpIO->Write(&xpaBCSetLoadData[lintBCSet].mintMemberCount, 1);
            lpIO->Write(xpaBCSetLoadData[lintBCSet].mpaMemberIDs,
                xpaBCSetLoadData[lintBCSet].mintMemberCount);
        }

        mpaBCSets[xintIndex].mintBCSetLength += lpIO->EndWriteBlock();
    }
    else
    {
        // BCSet has no Loads in the file.
        mpaBCSets[xintIndex].mintLoadsOffset = 0;
    }

    // write contact pairs
    if(xintContactPairTypesCount)
    {        
        mpaBCSets[xintIndex].mintContactPairsOffset =
            lpIO->BeginWriteBlock(mintSimModelOffset);

        for(UnsignedInt32 lintBCSet = 0; lintBCSet < xintContactPairTypesCount; lintBCSet++)
        {
            if(!xpaBCSetContactPairData[lintBCSet].mintMemberCount)
            {
                mpaBCSets[xintIndex].mintContactPairTypesCount--;
                continue;
            }
            mpaBCSets[xintIndex].mintContactPairMembersCount +=
                xpaBCSetContactPairData[lintBCSet].mintMemberCount;
            
            lpIO->Write(&xpaBCSetContactPairData[lintBCSet].mintMemberType, 1);
            lpIO->Write(&xpaBCSetContactPairData[lintBCSet].mintMemberCount, 1);
            lpIO->Write(xpaBCSetContactPairData[lintBCSet].mpaMemberIDs,
                xpaBCSetContactPairData[lintBCSet].mintMemberCount);
        }

        mpaBCSets[xintIndex].mintBCSetLength += lpIO->EndWriteBlock();
        
    }
    else
    {
        // BCSet has no ContactPairs in the file.
        mpaBCSets[xintIndex].mintContactPairsOffset = 0;
    }
    
    mSimModel.mintSimModelLength += mpaBCSets[xintIndex].mintBCSetLength;
    delete lpIO;
}

void CSimModel::WriteMaterial(UnsignedInt32 xintIndex,
                              UnsignedInt32 xintMaterialID,
                              UnsignedInt32 xintMaterialUniqueID,
                              UnsignedInt32 xintPropertiesCount,
                              SMaterialData* xpaMaterialData)
{
    if(!mpWriteFile)
        throw CCubitFile::eFileWriteError;
    if(xintIndex >= mSimModel.mintMaterialCount)
        throw CCubitFile::eNotFound;
    if(!mpaBCSets)
        throw CCubitFile::eOrderError;
    if(mpaMaterials[xintIndex].mintMaterialLength) // something that won't be zero
        throw CCubitFile::eDuplicateWrite;

    if(xintPropertiesCount)
    {  
        // An empty material is valid, but an incompletely defined one is not.
        if(!xpaMaterialData)
            throw CCubitFile::ePassedNullPointer;
        for(UnsignedInt32 lintMaterialData = 0; lintMaterialData < xintPropertiesCount; lintMaterialData++) {
            if(xpaMaterialData[lintMaterialData].mintMemberRows &&
                !xpaMaterialData[lintMaterialData].mpadblMemberData)
                throw CCubitFile::ePassedNullPointer;
        }
    }

    // Material header written out elsewhere.  Here, just populate the header struct
    mpaMaterials[xintIndex].mintMaterialID = xintMaterialID;
    mpaMaterials[xintIndex].mintMaterialUniqueID = xintMaterialUniqueID;
    mpaMaterials[xintIndex].mintMaterialPropertiesCount = xintPropertiesCount;
    
    CIOWrapper* lpIO = new CIOWrapper(mpWriteFile);

    // write material data
    if(xintPropertiesCount)
    {        
        mpaMaterials[xintIndex].mintMaterialPropertiesOffset =
            lpIO->BeginWriteBlock(mintSimModelOffset);

        for(UnsignedInt32 lintMaterial = 0; lintMaterial < xintPropertiesCount; lintMaterial++)
        {
            if(!xpaMaterialData[lintMaterial].mintMemberRows)
            {
                mpaMaterials[xintIndex].mintMaterialPropertiesCount--;
                continue;
            }
             mpaMaterials[xintIndex].mintPropertyDataCount +=
                (xpaMaterialData[lintMaterial].mintMemberRows *
                 xpaMaterialData[lintMaterial].mintMemberColumns);
            
            lpIO->Write(&xpaMaterialData[lintMaterial].mintMemberType, 1);
            lpIO->Write(&xpaMaterialData[lintMaterial].mintMemberRows, 1);
            lpIO->Write(&xpaMaterialData[lintMaterial].mintMemberColumns, 1);
            lpIO->Write(xpaMaterialData[lintMaterial].mpadblMemberData,
                        xpaMaterialData[lintMaterial].mintMemberRows *
                            xpaMaterialData[lintMaterial].mintMemberColumns);
            
        }

        mpaMaterials[xintIndex].mintMaterialLength += lpIO->EndWriteBlock();
    }
    else
    {
        // Material has no data in the file.
        mpaMaterials[xintIndex].mintMaterialPropertiesOffset = 0;
    }
    
    mSimModel.mintSimModelLength += mpaMaterials[xintIndex].mintMaterialLength;
    delete lpIO;
}

void CSimModel::WriteConstraint(UnsignedInt32 xintIndex,
        UnsignedInt32 xintConstraintID, UnsignedInt32 xintConstraintUniqueID,
        UnsignedInt32 xintConstraintType,
        UnsignedInt32 xintIndependentTypeCount, SConstraintData* xpaIndependentData,
        UnsignedInt32 xintDependentTypeCount, SConstraintData* xpaDependentData)
{
    if(!mpWriteFile)
        throw CCubitFile::eFileWriteError;
    if(xintIndex >= mSimModel2.mintConstraintCount)
        throw CCubitFile::eNotFound;
    if(!mpaBCSets)
        throw CCubitFile::eOrderError;
    if(mpaConstraints[xintIndex].mintConstraintLength) // something that won't be zero
        throw CCubitFile::eDuplicateWrite;

    if(xintDependentTypeCount)
    {  
        // An empty Constraint is valid, but an incompletely defined one is not.
        if(!xpaDependentData)
            throw CCubitFile::ePassedNullPointer;
        for(int lintConstraintData = 0; lintConstraintData < (int)xintDependentTypeCount; lintConstraintData++) {
            if(xpaDependentData[lintConstraintData].mintMemberCount &&
                !xpaDependentData[lintConstraintData].mpaMemberIDs)
                throw CCubitFile::ePassedNullPointer;
        }
    }
    if(xintIndependentTypeCount)
    {  
        // An empty Constraint is valid, but an incompletely defined one is not.
        if(!xpaIndependentData)
            throw CCubitFile::ePassedNullPointer;
        for(int lintConstraintData = 0; lintConstraintData < (int)xintIndependentTypeCount; lintConstraintData++) {
            if(xpaIndependentData[lintConstraintData].mintMemberCount &&
                !xpaIndependentData[lintConstraintData].mpaMemberIDs)
                throw CCubitFile::ePassedNullPointer;
        }
    }

    // Constraint header written out elsewhere.  Here, just populate the header struct
    mpaConstraints[xintIndex].mintConstraintID = xintConstraintID;
    mpaConstraints[xintIndex].mintConstraintUniqueID = xintConstraintUniqueID;
    mpaConstraints[xintIndex].mintConstraintType = xintConstraintType;
    mpaConstraints[xintIndex].mintIndependentTypeCount = xintIndependentTypeCount;
    mpaConstraints[xintIndex].mintDependentTypeCount = xintDependentTypeCount;
    
    CIOWrapper* lpIO = new CIOWrapper(mpWriteFile);

    // write Independent Constraint data
    if(xintIndependentTypeCount)
    {        
        mpaConstraints[xintIndex].mintIndependentDataOffset =
            lpIO->BeginWriteBlock(mintSimModelOffset);

        for(int lintConstraint = 0; lintConstraint < (int)xintIndependentTypeCount; lintConstraint++)
        {
            if(!xpaIndependentData[lintConstraint].mintMemberCount)
            {
                mpaConstraints[xintIndex].mintIndependentTypeCount--;
                continue;
            }

            lpIO->Write(&xpaIndependentData[lintConstraint].mintMemberType, 1);
            lpIO->Write(&xpaIndependentData[lintConstraint].mintMemberCount, 1);
            lpIO->Write(xpaIndependentData[lintConstraint].mpaMemberIDs,
                        xpaIndependentData[lintConstraint].mintMemberCount);
        }

        mpaConstraints[xintIndex].mintConstraintLength += lpIO->EndWriteBlock();
    }
    else
    {
        // Constraint has no data in the file.
        mpaConstraints[xintIndex].mintIndependentDataOffset = 0;
    }
    
    // write Dependent Constraint data
    if(xintDependentTypeCount)
    {        
        mpaConstraints[xintIndex].mintDependentDataOffset =
            lpIO->BeginWriteBlock(mintSimModelOffset);

        for(int lintConstraint = 0; lintConstraint < (int)xintDependentTypeCount; lintConstraint++)
        {
            if(!xpaDependentData[lintConstraint].mintMemberCount)
            {
                mpaConstraints[xintIndex].mintIndependentTypeCount--;
                continue;
            }

            lpIO->Write(&xpaDependentData[lintConstraint].mintMemberType, 1);
            lpIO->Write(&xpaDependentData[lintConstraint].mintMemberCount, 1);
            lpIO->Write(xpaDependentData[lintConstraint].mpaMemberIDs,
                        xpaDependentData[lintConstraint].mintMemberCount);
        }

        mpaConstraints[xintIndex].mintConstraintLength += lpIO->EndWriteBlock();
    }
    else
    {
        // Constraint has no data in the file.
        mpaConstraints[xintIndex].mintDependentDataOffset = 0;
    }
    
    mSimModel.mintSimModelLength += mpaConstraints[xintIndex].mintConstraintLength;

    delete lpIO;
}

UnsignedInt32 CSimModel::EndWrite()
{
    CIOWrapper lIO(mpWriteFile);

    if(mSimModel.mintBCCount) {
        lIO.BeginRewriteBlock(mintSimModelOffset, mSimModel.mintBCTableOffset);
        lIO.Write((UnsignedInt32*)mpaBCs,
            mSimModel.mintBCCount * mintSizeOfBCEntry);
        lIO.EndWriteBlock();
    }

    if(mSimModel.mintICCount) {
        lIO.BeginRewriteBlock(mintSimModelOffset, mSimModel.mintICTableOffset);
        lIO.Write((UnsignedInt32*)mpaICs,
            mSimModel.mintICCount * mintSizeOfICEntry);
        lIO.EndWriteBlock();
    }

    if(mSimModel.mintBCSetCount) {
        lIO.BeginRewriteBlock(mintSimModelOffset, mSimModel.mintBCSetTableOffset);
        lIO.Write((UnsignedInt32*)mpaBCSets,
            mSimModel.mintBCSetCount * mintSizeOfBCSetEntry);
        lIO.EndWriteBlock();
    }

    if(mSimModel.mintMaterialCount) {
        lIO.BeginRewriteBlock(mintSimModelOffset, mSimModel.mintMaterialTableOffset);
        lIO.Write((UnsignedInt32*)mpaMaterials,
            mSimModel.mintMaterialCount * mintSizeOfMaterialEntry);
        lIO.EndWriteBlock();
    }

    if(mSimModel.mintAmplitudeCount) {
        lIO.BeginRewriteBlock(mintSimModelOffset, mSimModel.mintAmplitudeTableOffset);
        lIO.Write((UnsignedInt32*)mpaAmplitudes,
            mSimModel.mintAmplitudeCount * mintSizeOfAmplitudeEntry);
        lIO.EndWriteBlock();
    }

    if(mSimModel2.mintConstraintCount) {
        lIO.BeginRewriteBlock(mintSimModelOffset, mSimModel2.mintConstraintTableOffset);
        lIO.Write((UnsignedInt32*)mpaConstraints,
            mSimModel2.mintConstraintCount * mintSizeOfConstraintEntry);
        lIO.EndWriteBlock();
    }

    UnsignedInt32 lintMetaDataLength;
    mBCMetaData.WriteMetaData(mpWriteFile, mSimModel.mintBCMetaDataOffset,
        lintMetaDataLength, mintSimModelOffset);
    mSimModel.mintSimModelLength += lintMetaDataLength;
    mICMetaData.WriteMetaData(mpWriteFile, mSimModel.mintICMetaDataOffset,
        lintMetaDataLength, mintSimModelOffset);
    mSimModel.mintSimModelLength += lintMetaDataLength;
    mBCSetMetaData.WriteMetaData(mpWriteFile, mSimModel.mintBCSetMetaDataOffset,
        lintMetaDataLength, mintSimModelOffset);
    mSimModel.mintSimModelLength += lintMetaDataLength;
    mMaterialMetaData.WriteMetaData(mpWriteFile, mSimModel.mintMaterialMetaDataOffset,
        lintMetaDataLength, mintSimModelOffset);
    mSimModel.mintSimModelLength += lintMetaDataLength;
    mAmplitudeMetaData.WriteMetaData(mpWriteFile, mSimModel.mintAmplitudeMetaDataOffset,
        lintMetaDataLength, mintSimModelOffset);
    mSimModel.mintSimModelLength += lintMetaDataLength;
    mConstraintMetaData.WriteMetaData(mpWriteFile, mSimModel2.mintConstraintMetaDataOffset,
        lintMetaDataLength, mintSimModelOffset);
    mSimModel.mintSimModelLength += lintMetaDataLength;

    lIO.BeginRewriteBlock(mintSimModelOffset, 0);
    lIO.Write((UnsignedInt32*)&mSimModel, mintSizeOfSimModelHeader);
    lIO.EndWriteBlock();

    lIO.BeginRewriteBlock(mintSimModelOffset, mSimModel.mintFutureTableOffset);
    lIO.Write((UnsignedInt32*)&mSimModel2, mintSizeOfSimModelHeader2);
    lIO.EndWriteBlock();

    mpWriteFile = NULL;
    return mSimModel.mintSimModelLength;
}


///////////////////////////////////////////////////////////////////////////////
// Read Functions
///////////////////////////////////////////////////////////////////////////////

void CSimModel::InitRead(FILE* xpFile,
                         UnsignedInt32& xintAbsoluteOffset,
                         UnsignedInt32& xintBCCount,
                         UnsignedInt32& xintICCount,
                         UnsignedInt32& xintBCSetCount,
                         UnsignedInt32& xintMaterialCount,
                         UnsignedInt32& xintAmplitudeCount,
                         UnsignedInt32& xintConstraintCount)
{
    if(mpReadFile)
        throw CCubitFile::eOrderError;

    mpReadFile = xpFile;
    mintSimModelOffset = xintAbsoluteOffset;
    CIOWrapper lIO(mpReadFile, xintAbsoluteOffset, 0);

    lIO.BeginReadBlock(mintSimModelOffset, 0);
    lIO.Read((UnsignedInt32*)&mSimModel, mintSizeOfSimModelHeader);
    lIO.EndReadBlock();

    lIO.BeginReadBlock(mintSimModelOffset, mSimModel.mintFutureTableOffset);
    lIO.Read((UnsignedInt32*)&mSimModel2, mintSizeOfSimModelHeader2);
    lIO.EndReadBlock();

    xintBCCount = mSimModel.mintBCCount;
    xintICCount = mSimModel.mintICCount;
    xintBCSetCount = mSimModel.mintBCSetCount;
    xintMaterialCount = mSimModel.mintMaterialCount;
    xintAmplitudeCount = mSimModel.mintAmplitudeCount;
    xintConstraintCount = mSimModel2.mintConstraintCount;

    // Restore the Boundary Condition definition header table if there is one.
    if(mSimModel.mintBCCount) {
        mpaBCs = new SCubitFileBCEntry[mSimModel.mintBCCount];

        lIO.BeginReadBlock(mintSimModelOffset, mSimModel.mintBCTableOffset);
        lIO.Read((UnsignedInt32*)mpaBCs,
            mSimModel.mintBCCount * mintSizeOfBCEntry);
        lIO.EndReadBlock();
    }

    // Restore the Initial Condition definition header table if there is one.
    if(mSimModel.mintICCount) {
        mpaICs = new SCubitFileICEntry[mSimModel.mintICCount];

        lIO.BeginReadBlock(mintSimModelOffset, mSimModel.mintICTableOffset);
        lIO.Read((UnsignedInt32*)mpaICs,
            mSimModel.mintICCount * mintSizeOfICEntry);
        lIO.EndReadBlock();
    }

    // Restore the BCSet definition header table if there is one.
    if(mSimModel.mintBCSetCount) {
        mpaBCSets = new SCubitFileBCSetEntry[mSimModel.mintBCSetCount];

        lIO.BeginReadBlock(mintSimModelOffset, mSimModel.mintBCSetTableOffset);
        lIO.Read((UnsignedInt32*)mpaBCSets,
            mSimModel.mintBCSetCount * mintSizeOfBCSetEntry);
        lIO.EndReadBlock();
    }

    // Restore the Material definition header table if there is one.
    if(mSimModel.mintMaterialCount) {
        mpaMaterials = new SCubitFileMaterialEntry[mSimModel.mintMaterialCount];

        lIO.BeginReadBlock(mintSimModelOffset, mSimModel.mintMaterialTableOffset);
        lIO.Read((UnsignedInt32*)mpaMaterials,
            mSimModel.mintMaterialCount * mintSizeOfMaterialEntry);
        lIO.EndReadBlock();
    }

    // Restore the Amplitude definition header table if there is one.
    if(mSimModel.mintAmplitudeCount) {
        mpaAmplitudes = new SCubitFileAmplitudeEntry[mSimModel.mintAmplitudeCount];

        lIO.BeginReadBlock(mintSimModelOffset, mSimModel.mintAmplitudeTableOffset);
        lIO.Read((UnsignedInt32*)mpaAmplitudes,
            mSimModel.mintAmplitudeCount * mintSizeOfAmplitudeEntry);
        lIO.EndReadBlock();
    }

    // Restore the Amplitude definition header table if there is one.
    if(mSimModel2.mintConstraintCount) {
        mpaConstraints = new SCubitFileConstraintEntry[mSimModel2.mintConstraintCount];

        lIO.BeginReadBlock(mintSimModelOffset, mSimModel2.mintConstraintTableOffset);
        lIO.Read((UnsignedInt32*)mpaConstraints,
            mSimModel2.mintConstraintCount * mintSizeOfConstraintEntry);
        lIO.EndReadBlock();
    }

    mBCMetaData.ReadMetaData(mpReadFile, mintSimModelOffset,
        mSimModel.mintBCMetaDataOffset, mSimModel.mintSimModelEndian);
    mICMetaData.ReadMetaData(mpReadFile, mintSimModelOffset,
        mSimModel.mintICMetaDataOffset, mSimModel.mintSimModelEndian);
    mBCSetMetaData.ReadMetaData(mpReadFile, mintSimModelOffset,
        mSimModel.mintBCSetMetaDataOffset, mSimModel.mintSimModelEndian);
    mMaterialMetaData.ReadMetaData(mpReadFile, mintSimModelOffset,
        mSimModel.mintMaterialMetaDataOffset, mSimModel.mintSimModelEndian);
    mAmplitudeMetaData.ReadMetaData(mpReadFile, mintSimModelOffset,
        mSimModel.mintAmplitudeMetaDataOffset, mSimModel.mintSimModelEndian);
    mConstraintMetaData.ReadMetaData(mpReadFile, mintSimModelOffset,
        mSimModel2.mintConstraintMetaDataOffset, mSimModel.mintSimModelEndian);
}

void CSimModel::ReadBCSet(UnsignedInt32 xintIndex,
                          UnsignedInt32& xintBCSetID,
                          UnsignedInt32& xintBCSetUniqueID,
                          UnsignedInt32& xintBCSetAnalysisType,
                          UnsignedInt32& xintRestraintTypesCount,
                          UnsignedInt32& xintLoadTypesCount,
                          UnsignedInt32& xintContactPairTypesCount,
                          SBCSetData*& xpaBCSetRestraintData,
                          SBCSetData*& xpaBCSetLoadData,
                          SBCSetData*& xpaBCSetContactPairData)
{
    if(!mpReadFile)
        throw CCubitFile::eFileReadError;
    if(xintIndex >= mSimModel.mintBCSetCount)
        throw CCubitFile::eNotFound;
    if(!mpaBCSets)
        throw CCubitFile::eOrderError;

    xintBCSetID = mpaBCSets[xintIndex].mintBCSetID;
    xintBCSetUniqueID = mpaBCSets[xintIndex].mintBCSetUniqueID;
    xintBCSetAnalysisType = mpaBCSets[xintIndex].mintBCSetAnalysisType;
    xintRestraintTypesCount = mpaBCSets[xintIndex].mintRestraintTypesCount;
    xintLoadTypesCount = mpaBCSets[xintIndex].mintLoadTypesCount;
    xintContactPairTypesCount = mpaBCSets[xintIndex].mintContactPairTypesCount;

    // read Restraints
    if(xintRestraintTypesCount) 
    {
        // Resize the node set return buffer if necessary and then set the return
        // pointers to the buffer.
        xpaBCSetRestraintData = AdjustBuffer(xintRestraintTypesCount,
            mBCSetBuff.mintNumRestraintTypes, mBCSetBuff.mpaBCSetRestraintData);
        UnsignedInt32* lpIDs =
            AdjustBuffer(mpaBCSets[xintIndex].mintRestraintMembersCount,
            mBCSetBuff.mintNumRestraintMembers, mBCSetBuff.mpaintRestraintMemberIDs);

        // Read the node set from the file.
        UnsignedInt32 lintNumMembers, lintTotalMembers = 0;
        CIOWrapper* lpIO = new CIOWrapper(mpReadFile, mSimModel.mintSimModelEndian);
        lpIO->BeginReadBlock(mintSimModelOffset,
            mpaBCSets[xintIndex].mintRestraintsOffset);
        long start_location = lpIO->GetLocation();
        for(UnsignedInt32 lintType = 0; lintType < xintRestraintTypesCount; lintType++) {
            lpIO->Read(&xpaBCSetRestraintData[lintType].mintMemberType, 1);
            lpIO->Read(&xpaBCSetRestraintData[lintType].mintMemberCount, 1);

            xpaBCSetRestraintData[lintType].mpaMemberIDs = lpIDs;
            lintNumMembers = xpaBCSetRestraintData[lintType].mintMemberCount;
            // Make sure the total number of node set members does not exceed
            // what was specified in the node set table entry.
            lintTotalMembers += lintNumMembers;
            if(lintTotalMembers > mpaBCSets[xintIndex].mintRestraintMembersCount)
                throw CCubitFile::eFileReadError;

            lpIO->Read(xpaBCSetRestraintData[lintType].mpaMemberIDs, lintNumMembers);

            lpIDs = &lpIDs[lintNumMembers];
        }

        lpIO->EndReadBlock();
        delete lpIO;
    }
    else
        xpaBCSetRestraintData = NULL;

    // read Loads
    if(xintLoadTypesCount) 
    {
        // Resize the node set return buffer if necessary and then set the return
        // pointers to the buffer.
        xpaBCSetLoadData = AdjustBuffer(xintLoadTypesCount,
            mBCSetBuff.mintNumLoadTypes, mBCSetBuff.mpaBCSetLoadData);
        UnsignedInt32* lpIDs =
            AdjustBuffer(mpaBCSets[xintIndex].mintLoadMembersCount,
            mBCSetBuff.mintNumLoadMembers, mBCSetBuff.mpaintLoadMemberIDs);

        // Read the node set from the file.
        UnsignedInt32 lintNumMembers, lintTotalMembers = 0;
        CIOWrapper* lpIO = new CIOWrapper(mpReadFile, mSimModel.mintSimModelEndian);
        lpIO->BeginReadBlock(mintSimModelOffset,
            mpaBCSets[xintIndex].mintLoadsOffset);
        long start_location = lpIO->GetLocation();
        for(UnsignedInt32 lintType = 0; lintType < xintLoadTypesCount; lintType++) {
            lpIO->Read(&xpaBCSetLoadData[lintType].mintMemberType, 1);
            lpIO->Read(&xpaBCSetLoadData[lintType].mintMemberCount, 1);

            xpaBCSetLoadData[lintType].mpaMemberIDs = lpIDs;
            lintNumMembers = xpaBCSetLoadData[lintType].mintMemberCount;
            // Make sure the total number of node set members does not exceed
            // what was specified in the node set table entry.
            lintTotalMembers += lintNumMembers;
            if(lintTotalMembers > mpaBCSets[xintIndex].mintLoadMembersCount)
                throw CCubitFile::eFileReadError;

            lpIO->Read(xpaBCSetLoadData[lintType].mpaMemberIDs, lintNumMembers);

            lpIDs = &lpIDs[lintNumMembers];
        }

        lpIO->EndReadBlock();
        delete lpIO;
    }
    else
        xpaBCSetLoadData = NULL;
    
    // read ContactPairs
    if(xintContactPairTypesCount) 
    {
        // Resize the node set return buffer if necessary and then set the return
        // pointers to the buffer.
        xpaBCSetContactPairData = AdjustBuffer(xintContactPairTypesCount,
            mBCSetBuff.mintNumContactPairTypes, mBCSetBuff.mpaBCSetContactPairData);
        UnsignedInt32* lpIDs =
            AdjustBuffer(mpaBCSets[xintIndex].mintContactPairMembersCount,
            mBCSetBuff.mintNumContactPairMembers, mBCSetBuff.mpaintContactPairMemberIDs);

        // Read the node set from the file.
        UnsignedInt32 lintNumMembers, lintTotalMembers = 0;
        CIOWrapper* lpIO = new CIOWrapper(mpReadFile, mSimModel.mintSimModelEndian);
        lpIO->BeginReadBlock(mintSimModelOffset,
            mpaBCSets[xintIndex].mintContactPairsOffset);
        long start_location = lpIO->GetLocation();
        for(UnsignedInt32 lintType = 0; lintType < xintContactPairTypesCount; lintType++) {
            lpIO->Read(&xpaBCSetContactPairData[lintType].mintMemberType, 1);
            lpIO->Read(&xpaBCSetContactPairData[lintType].mintMemberCount, 1);

            xpaBCSetContactPairData[lintType].mpaMemberIDs = lpIDs;
            lintNumMembers = xpaBCSetContactPairData[lintType].mintMemberCount;
            // Make sure the total number of node set members does not exceed
            // what was specified in the node set table entry.
            lintTotalMembers += lintNumMembers;
            if(lintTotalMembers > mpaBCSets[xintIndex].mintContactPairMembersCount)
                throw CCubitFile::eFileReadError;

            lpIO->Read(xpaBCSetContactPairData[lintType].mpaMemberIDs, lintNumMembers);

            lpIDs = &lpIDs[lintNumMembers];
        }

        lpIO->EndReadBlock();
        delete lpIO;
    }
    else
        xpaBCSetContactPairData = NULL;

}

void CSimModel::ReadMaterial(UnsignedInt32 xintIndex,
                             UnsignedInt32& xintMaterialID,
                             UnsignedInt32& xintMaterialUniqueID,
                             UnsignedInt32& xintPropertiesCount,
                             SMaterialData*& xpaMaterialData)
{
    if(!mpReadFile)
        throw CCubitFile::eFileReadError;
    if(xintIndex >= mSimModel.mintMaterialCount)
        throw CCubitFile::eNotFound;
    if(!mpaMaterials)
        throw CCubitFile::eOrderError;

    xintMaterialID = mpaMaterials[xintIndex].mintMaterialID;
    xintMaterialUniqueID = mpaMaterials[xintIndex].mintMaterialUniqueID;
    xintPropertiesCount = mpaMaterials[xintIndex].mintMaterialPropertiesCount;

    // read material data
    if(xintPropertiesCount) 
    {
        // Resize the material return buffer if necessary and then set the return
        // pointers to the buffer.
        xpaMaterialData = AdjustBuffer(xintPropertiesCount,
            mMaterialBuff.mintNumDataTypes, mMaterialBuff.mpaMaterialData);
        
        double* lpData =
            AdjustBuffer(mpaMaterials[xintIndex].mintPropertyDataCount,
                mMaterialBuff.mintNumDataMembers, mMaterialBuff.mpadblData);


        // Read the material property data from the file.
        UnsignedInt32 lintNumMembers, lintTotalMembers = 0;
        CIOWrapper* lpIO = new CIOWrapper(mpReadFile, mSimModel.mintSimModelEndian);
        lpIO->BeginReadBlock(mintSimModelOffset,
            mpaMaterials[xintIndex].mintMaterialPropertiesOffset);

        long start_location = lpIO->GetLocation();
        for(UnsignedInt32 lintType = 0; lintType < xintPropertiesCount; lintType++) {
            lpIO->Read(&xpaMaterialData[lintType].mintMemberType, 1);
            lpIO->Read(&xpaMaterialData[lintType].mintMemberRows, 1);
            lpIO->Read(&xpaMaterialData[lintType].mintMemberColumns, 1);

            xpaMaterialData[lintType].mpadblMemberData = lpData;

            lintNumMembers = xpaMaterialData[lintType].mintMemberRows *
                xpaMaterialData[lintType].mintMemberColumns;

            // Make sure the total number of material members does not exceed
            // what was specified in the material table entry.
            lintTotalMembers += lintNumMembers;
            if(lintTotalMembers > mpaMaterials[xintIndex].mintPropertyDataCount)
                throw CCubitFile::eFileReadError;
            
            lpIO->Read(xpaMaterialData[lintType].mpadblMemberData,
                       xpaMaterialData[lintType].mintMemberRows *
                        xpaMaterialData[lintType].mintMemberColumns);

            lpData = &lpData[lintNumMembers];
        }

        lpIO->EndReadBlock();
        delete lpIO;
    }
    else
        xpaMaterialData = NULL;
}

void CSimModel::ReadConstraint(UnsignedInt32 xintIndex,
                               UnsignedInt32& xintConstraintID, UnsignedInt32& xintConstraintUniqueID,
                               UnsignedInt32& xintConstraintType,
                               UnsignedInt32& xintIndependentTypeCount, SConstraintData*& xpaIndependentData,
                               UnsignedInt32& xintDependentTypeCount, SConstraintData*& xpaDependentData)
{
    if(!mpReadFile)
        throw CCubitFile::eFileReadError;
    if(xintIndex >= mSimModel2.mintConstraintCount)
        throw CCubitFile::eNotFound;
    if(!mpaConstraints)
        throw CCubitFile::eOrderError;

    xintConstraintID = mpaConstraints[xintIndex].mintConstraintID;
    xintConstraintUniqueID = mpaConstraints[xintIndex].mintConstraintUniqueID;
    xintConstraintType = mpaConstraints[xintIndex].mintConstraintType;
    xintIndependentTypeCount = mpaConstraints[xintIndex].mintIndependentTypeCount;
    xintDependentTypeCount = mpaConstraints[xintIndex].mintDependentTypeCount;

    // read Dependent Constraint data
    if(xintDependentTypeCount) 
    {
        // Resize the Constraint return buffer if necessary and then set the return
        // pointers to the buffer.
        xpaDependentData = AdjustBuffer(xintDependentTypeCount,
            mConstraintBuff.mintNumDependentTypes, mConstraintBuff.mpaDependentData);
        UnsignedInt32* lpIDs =
            AdjustBuffer(mpaConstraints[xintIndex].mintDependentTypeCount,
            mConstraintBuff.mintNumDependentMembers, mConstraintBuff.mpaintDependentIDs);

        // Read the Constraint property data from the file.
        UnsignedInt32 lintNumMembers, lintTotalMembers = 0;
        CIOWrapper* lpIO = new CIOWrapper(mpReadFile, mSimModel.mintSimModelEndian);
        lpIO->BeginReadBlock(mintSimModelOffset,
            mpaConstraints[xintIndex].mintDependentDataOffset);

        long start_location = lpIO->GetLocation();
        for(UnsignedInt32 lintType = 0; lintType < xintDependentTypeCount; lintType++)
        {
            lpIO->Read(&xpaDependentData[lintType].mintMemberType, 1);
            lpIO->Read(&xpaDependentData[lintType].mintMemberCount, 1);

            xpaDependentData[lintType].mpaMemberIDs = lpIDs;
            lintNumMembers = xpaDependentData[lintType].mintMemberCount;
            // Make sure the total number of node set members does not exceed
            // what was specified in the node set table entry.
            lintTotalMembers += lintNumMembers;
            if(lintTotalMembers > mpaConstraints[xintIndex].mintDependentTypeCount)
                throw CCubitFile::eFileReadError;

            lpIO->Read(xpaDependentData[lintType].mpaMemberIDs, lintNumMembers);

            lpIDs = &lpIDs[lintNumMembers];
        }

        lpIO->EndReadBlock();
        delete lpIO;
    }
    else
        xpaDependentData = NULL;

    // read Independent Constraint data
    if(xintIndependentTypeCount) 
    {
        // Resize the Constraint return buffer if necessary and then set the return
        // pointers to the buffer.
        xpaIndependentData = AdjustBuffer(xintIndependentTypeCount,
            mConstraintBuff.mintNumIndependentTypes, mConstraintBuff.mpaIndependentData);
        UnsignedInt32* lpIDs =
            AdjustBuffer(mpaConstraints[xintIndex].mintIndependentTypeCount,
            mConstraintBuff.mintNumIndependentMembers, mConstraintBuff.mpaintIndependentIDs);

        // Read the Constraint property data from the file.
        UnsignedInt32 lintNumMembers, lintTotalMembers = 0;
        CIOWrapper* lpIO = new CIOWrapper(mpReadFile, mSimModel.mintSimModelEndian);
        lpIO->BeginReadBlock(mintSimModelOffset,
            mpaConstraints[xintIndex].mintIndependentDataOffset);

        long start_location = lpIO->GetLocation();
        for(UnsignedInt32 lintType = 0; lintType < xintIndependentTypeCount; lintType++)
        {
            lpIO->Read(&xpaIndependentData[lintType].mintMemberType, 1);
            lpIO->Read(&xpaIndependentData[lintType].mintMemberCount, 1);

            xpaIndependentData[lintType].mpaMemberIDs = lpIDs;
            lintNumMembers = xpaIndependentData[lintType].mintMemberCount;
            // Make sure the total number of node set members does not exceed
            // what was specified in the node set table entry.
            lintTotalMembers += lintNumMembers;
            if(lintTotalMembers > mpaConstraints[xintIndex].mintIndependentTypeCount)
                throw CCubitFile::eFileReadError;

            lpIO->Read(xpaIndependentData[lintType].mpaMemberIDs, lintNumMembers);

            lpIDs = &lpIDs[lintNumMembers];
        }

        lpIO->EndReadBlock();
        delete lpIO;
    }
    else
        xpaIndependentData = NULL;

}

void CSimModel::EndRead()
{
    mpReadFile = NULL;
}

///////////////////////////////////////////////////////////////////////////////
// Meta-Data Functions
///////////////////////////////////////////////////////////////////////////////

CMetaData& CSimModel::GetBCMetaData()
{
    return mBCMetaData;
}

CMetaData& CSimModel::GetICMetaData()
{
    return mICMetaData;
}

CMetaData& CSimModel::GetBCSetMetaData()
{
    return mBCSetMetaData;
}

CMetaData& CSimModel::GetMaterialMetaData()
{
    return mMaterialMetaData;
}

CMetaData& CSimModel::GetAmplitudeMetaData()
{
    return mAmplitudeMetaData;
}

CMetaData& CSimModel::GetConstraintMetaData()
{
    return mConstraintMetaData;
}

