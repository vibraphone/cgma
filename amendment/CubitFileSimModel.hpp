/*******************************************************************************
    COPYRIGHT 2010 CATERPILLAR INC.  ALL RIGHTS RESERVED

    This program is the property of Caterpillar Inc., includes Caterpillar's
    confidential and trade secret information, and is maintained as an
    unpublished copyrighted work.  It is not to be copied or used by others
    except under license from Caterpillar.  This program is also protected as an
    unpublished work in accordance with the copyright act of 1976.  In the event
    of either inadvertent or deliberate publication, Caterpillar Inc. intends to
    maintain copyright protection for this work under the relevant copyright
    laws pertaining to published works.  The inclusion of a copyright notice
    hereon is precautionary only, and does not imply publication or disclosure.

 
 Filename      : CubitFileSimModel.hpp

 Purpose       : Defines an interface for the reading and writing functionality
                 for a Simulation model section of a Cubit (*.cub) format file.
           
 Special Notes :

 Creator       : Andrew Rout

 Creation Date : 02/20/10

 Owner         : Andrew Rout

*******************************************************************************/

#ifndef CubitFileSimModel_HPP
#define CubitFileSimModel_HPP

#include "CubitFileMetaData.hpp"
#include "CubitUtilConfigure.h"
#include <vector>

namespace NCubitFile {

class CUBIT_UTIL_EXPORT CSimModel {
public:
    CSimModel();
    virtual ~CSimModel();
    
    UnsignedInt32 InitWrite(FILE* xpFile,
        UnsignedInt32 xintBCCount, UnsignedInt32 xintICCount,
        UnsignedInt32 xintBCSetCount, UnsignedInt32 xintMaterialCount,
        UnsignedInt32 xintAmplitudeCount, UnsignedInt32 xintConstraintCount);
    void WriteBCSet(UnsignedInt32 xintIndex,
        UnsignedInt32 xintBCSetID, UnsignedInt32 xintBCSetUniqueID,
        UnsignedInt32 xintBCSetAnalysisType, UnsignedInt32 xintRestraintTypesCount,
        UnsignedInt32 xintLoadTypesCount, UnsignedInt32 xintContactPairTypesCount,
        SBCSetData* xpaBCSetRestraintData, SBCSetData* xpaBCSetLoadData,
        SBCSetData* xpaBCSetContactPairData);
    void WriteMaterial(UnsignedInt32 xintIndex,
        UnsignedInt32 xintMaterialID, UnsignedInt32 xintMaterialUniqueID,
        UnsignedInt32 xintPropertiesCount, SMaterialData* xpaMaterialData);
    void WriteConstraint(UnsignedInt32 xintIndex,
        UnsignedInt32 xintConstraintID, UnsignedInt32 xintConstraintUniqueID,
        UnsignedInt32 xintConstraintType,
        UnsignedInt32 xintIndependentTypeCount, SConstraintData* xpaIndependentData,
        UnsignedInt32 xintDependentTypeCount, SConstraintData* xpaDependentData);
    UnsignedInt32 EndWrite();
    
    void InitRead(FILE* xpFile, UnsignedInt32& xintAbsoluteOffset,
        UnsignedInt32& xintBCCount, UnsignedInt32& xintICCount,
        UnsignedInt32& xintBCSetCount, UnsignedInt32& xintMaterialCount,
        UnsignedInt32& xintAmplitudeCount, UnsignedInt32& xintConstraintCount);
    void ReadBCSet(UnsignedInt32 xintIndex,
        UnsignedInt32& xintBCSetID, UnsignedInt32& xintBCSetUniqueID,
        UnsignedInt32& xintBCSetAnalysisType, UnsignedInt32& xintRestraintTypesCount,
        UnsignedInt32& xintLoadTypesCount, UnsignedInt32& xintContactPairTypesCount,
        SBCSetData*& xpaBCSetRestraintData, SBCSetData*& xpaBCSetLoadData,
        SBCSetData*& xpaBCSetContactPairData);
    void ReadMaterial(UnsignedInt32 xintIndex,
        UnsignedInt32& xintMaterialID, UnsignedInt32& xintMaterialUniqueID,
        UnsignedInt32& xintPropertiesCount, SMaterialData*& xpaMaterialData);
    void ReadConstraint(UnsignedInt32 xintIndex,
        UnsignedInt32& xintConstraintID, UnsignedInt32& xintConstraintUniqueID,
        UnsignedInt32& xintConstraintType,
        UnsignedInt32& xintIndependentTypeCount, SConstraintData*& xpaIndependentData,
        UnsignedInt32& xintDependentTypeCount, SConstraintData*& xpaDependentData);
    void EndRead();
    
    CMetaData& GetBCMetaData();
    CMetaData& GetICMetaData();
    CMetaData& GetBCSetMetaData();
    CMetaData& GetMaterialMetaData();
    CMetaData& GetAmplitudeMetaData();
    CMetaData& GetConstraintMetaData();
    
private:
    FILE* mpReadFile;
    FILE* mpWriteFile;
    UnsignedInt32 mintSimModelOffset;
    CMetaData mBCMetaData;
    CMetaData mICMetaData;
    CMetaData mBCSetMetaData;
    CMetaData mMaterialMetaData;
    CMetaData mAmplitudeMetaData;
    CMetaData mConstraintMetaData;

    // Data storage structures:
    // CAUTION: These structures must be 64 bit word aligned!!!
    struct SCubitFileSimModelHeader {
        UnsignedInt32 mintSimModelEndian;
        UnsignedInt32 mintSimModelSchema;
        UnsignedInt32 mintSimModelCompress;
        UnsignedInt32 mintSimModelLength;
        UnsignedInt32 mintBCCount;
        UnsignedInt32 mintBCTableOffset;
        UnsignedInt32 mintBCMetaDataOffset;
        UnsignedInt32 mintICCount;
        UnsignedInt32 mintICTableOffset;
        UnsignedInt32 mintICMetaDataOffset;
        UnsignedInt32 mintBCSetCount;
        UnsignedInt32 mintBCSetTableOffset;
        UnsignedInt32 mintBCSetMetaDataOffset;
        UnsignedInt32 mintMaterialCount;
        UnsignedInt32 mintMaterialTableOffset;
        UnsignedInt32 mintMaterialMetaDataOffset;
        UnsignedInt32 mintAmplitudeCount;
        UnsignedInt32 mintAmplitudeTableOffset;
        UnsignedInt32 mintAmplitudeMetaDataOffset;
        UnsignedInt32 mintFutureTableOffset; // for forwards compatibility
        //UnsignedInt32 mintPadFor64bitOS; // only if odd number of items in struct
    } mSimModel;
    static const UnsignedInt32 mintSizeOfSimModelHeader;
    //
    // Continuation from previous table
    // CAUTION: These structures must be 64 bit word aligned!!!
    struct SCubitFileSimModelHeader2 {
        UnsignedInt32 mintConstraintCount;
        UnsignedInt32 mintConstraintTableOffset;
        UnsignedInt32 mintConstraintMetaDataOffset;
        UnsignedInt32 mintFutureTableOffset; // for forwards compatibility
        //UnsignedInt32 mintPadFor64bitOS; // only if odd number of items in struct
    } mSimModel2;
    static const UnsignedInt32 mintSizeOfSimModelHeader2;

    struct SCubitFileBCSetEntry {
        UnsignedInt32 mintBCSetID;
        UnsignedInt32 mintBCSetUniqueID;
        UnsignedInt32 mintBCSetAnalysisType;
        UnsignedInt32 mintRestraintTypesCount;
        UnsignedInt32 mintRestraintMembersCount;
        UnsignedInt32 mintRestraintsOffset;
        UnsignedInt32 mintLoadTypesCount;
        UnsignedInt32 mintLoadMembersCount;
        UnsignedInt32 mintLoadsOffset;
        UnsignedInt32 mintContactPairTypesCount;
        UnsignedInt32 mintContactPairMembersCount;
        UnsignedInt32 mintContactPairsOffset;
        UnsignedInt32 mintBCSetLength;
        UnsignedInt32 mintPadFor64bitOS;
    } *mpaBCSets;
    static const UnsignedInt32 mintSizeOfBCSetEntry;

    struct SCubitFileMaterialEntry {
        UnsignedInt32 mintMaterialID;
        UnsignedInt32 mintMaterialUniqueID;
        UnsignedInt32 mintMaterialPropertiesOffset;
        UnsignedInt32 mintMaterialPropertiesCount;
        UnsignedInt32 mintPropertyDataCount; // total num of all property data
        UnsignedInt32 mintMaterialLength;
    } *mpaMaterials;
    static const UnsignedInt32 mintSizeOfMaterialEntry;

    struct SCubitFileConstraintEntry {
        UnsignedInt32 mintConstraintID;
        UnsignedInt32 mintConstraintUniqueID;
        UnsignedInt32 mintConstraintType;
        UnsignedInt32 mintIndependentTypeCount; // usually 1 independent node/vertex
        UnsignedInt32 mintIndependentDataOffset;
        UnsignedInt32 mintDependentTypeCount;
        UnsignedInt32 mintDependentDataOffset;
        UnsignedInt32 mintConstraintLength;
        //UnsignedInt32 mintPadFor64bitOS;
    } *mpaConstraints;
    static const UnsignedInt32 mintSizeOfConstraintEntry;

    // These structs aren't used yet (12.1), can modify them until they are first used.
    // Then, they must be locked in.
    struct SCubitFileBCEntry {
        UnsignedInt32 mintBCID;
        UnsignedInt32 mintBCUniqueID;
        UnsignedInt32 mintBCType;
        UnsignedInt32 mintDataLength;
        UnsignedInt32 mintDataOffset;
        UnsignedInt32 mintBCLength;
        //UnsignedInt32 mintPadFor64bitOS;
    } *mpaBCs;
    static const UnsignedInt32 mintSizeOfBCEntry;

    struct SCubitFileICEntry {
        UnsignedInt32 mintICID;
        UnsignedInt32 mintICUniqueID;
        UnsignedInt32 mintICType;
        UnsignedInt32 mintDataLength;
        UnsignedInt32 mintDataOffset;
        UnsignedInt32 mintICCLength;
        //UnsignedInt32 mintPadFor64bitOS;
    } *mpaICs;
    static const UnsignedInt32 mintSizeOfICEntry;

    struct SCubitFileAmplitudeEntry {
        UnsignedInt32 mintAmplitudeID;
        UnsignedInt32 mintAmplitudeUniqueID;
        UnsignedInt32 mintDataLength;
        UnsignedInt32 mintDataOffset;
        UnsignedInt32 mintAmplitudeLength;
        UnsignedInt32 mintPadFor64bitOS;
    } *mpaAmplitudes;
    static const UnsignedInt32 mintSizeOfAmplitudeEntry;
    
    // Buffers for read operations.
    struct SBCSetReturnBuffer {
        //restraints
        UnsignedInt32 mintNumRestraintTypes;
        UnsignedInt32 mintNumRestraintMembers;
        SBCSetData* mpaBCSetRestraintData;
        UnsignedInt32* mpaintRestraintMemberIDs;
        //loads
        UnsignedInt32 mintNumLoadTypes;
        UnsignedInt32 mintNumLoadMembers;
        SBCSetData* mpaBCSetLoadData;
        UnsignedInt32* mpaintLoadMemberIDs;
        //contact pairs
        UnsignedInt32 mintNumContactPairTypes;
        UnsignedInt32 mintNumContactPairMembers;
        SBCSetData* mpaBCSetContactPairData;
        UnsignedInt32* mpaintContactPairMemberIDs;
    } mBCSetBuff;

    struct SMaterialReturnBuffer {
        UnsignedInt32 mintNumDataTypes;
        UnsignedInt32 mintNumDataMembers; //i.e. current size of buffer
        SMaterialData* mpaMaterialData;
        double* mpadblData;
    } mMaterialBuff;

    struct SConstraintReturnBuffer {
        UnsignedInt32 mintNumIndependentTypes;
        UnsignedInt32 mintNumIndependentMembers;
        SConstraintData* mpaIndependentData;
        UnsignedInt32* mpaintIndependentIDs;
        //--
        UnsignedInt32 mintNumDependentTypes;
        UnsignedInt32 mintNumDependentMembers;
        SConstraintData* mpaDependentData;
        UnsignedInt32* mpaintDependentIDs;
    } mConstraintBuff;


    template <class TBuffer>
    TBuffer* AdjustBuffer(UnsignedInt32 xintRequiredSize,
        UnsignedInt32& xintActualSize, TBuffer*& xpaBuffer)
    {
        if(!xintRequiredSize)  return NULL;  // Nothing requested, return nothing.

        if(xintActualSize < xintRequiredSize) {
            if(xpaBuffer)
                delete [] xpaBuffer;
            xintActualSize = xintRequiredSize;
            xpaBuffer = new TBuffer[xintActualSize];
            if(!xpaBuffer) {
                xintActualSize = 0;
                throw CCubitFile::eMemoryError;
            }
        }
        return xpaBuffer;
    }
};

} // namespace NCubitFile

#endif

