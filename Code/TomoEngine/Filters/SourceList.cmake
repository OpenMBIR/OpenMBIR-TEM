#--////////////////////////////////////////////////////////////////////////////
#-- Copyright (c) 2009, Michael A. Jackson. BlueQuartz Software
#-- All rights reserved.
#-- BSD License: http://www.opensource.org/licenses/bsd-license.html
#-- This code was partly written under US Air Force Contract FA8650-07-D-5800
#--////////////////////////////////////////////////////////////////////////////
set (TomoEngine_Filters_SRCS
    ${TomoEngine_SOURCE_DIR}/Filters/CalculateAMatrixColumn.cpp
    ${TomoEngine_SOURCE_DIR}/Filters/CostData.cpp
    ${TomoEngine_SOURCE_DIR}/Filters/MRCSinogramInitializer.cpp
    ${TomoEngine_SOURCE_DIR}/Filters/RawSinogramInitializer.cpp
    ${TomoEngine_SOURCE_DIR}/Filters/GainsOffsetsReader.cpp
    ${TomoEngine_SOURCE_DIR}/Filters/InitialReconstructionInitializer.cpp
    ${TomoEngine_SOURCE_DIR}/Filters/InitialReconstructionBinReader.cpp
    ${TomoEngine_SOURCE_DIR}/Filters/TomoFilter.cpp
    ${TomoEngine_SOURCE_DIR}/Filters/ComputeInitialOffsets.cpp
    ${TomoEngine_SOURCE_DIR}/Filters/DetectorResponse.cpp
    ${TomoEngine_SOURCE_DIR}/Filters/DetectorResponseWriter.cpp
    ${TomoEngine_SOURCE_DIR}/Filters/TargetGainSigmaXEstimation.cpp
)

set (TomoEngine_Filters_HDRS
    ${TomoEngine_SOURCE_DIR}/Filters/CalculateAMatrixColumn.h
    ${TomoEngine_SOURCE_DIR}/Filters/CostData.h
    ${TomoEngine_SOURCE_DIR}/Filters/MRCSinogramInitializer.h
    ${TomoEngine_SOURCE_DIR}/Filters/RawSinogramInitializer.h
    ${TomoEngine_SOURCE_DIR}/Filters/GainsOffsetsReader.h
    ${TomoEngine_SOURCE_DIR}/Filters/InitialReconstructionInitializer.h
    ${TomoEngine_SOURCE_DIR}/Filters/InitialReconstructionBinReader.h
    ${TomoEngine_SOURCE_DIR}/Filters/TomoFilter.h
    ${TomoEngine_SOURCE_DIR}/Filters/ComputeInitialOffsets.h 
    ${TomoEngine_SOURCE_DIR}/Filters/DetectorResponse.h
    ${TomoEngine_SOURCE_DIR}/Filters/DetectorResponseWriter.h
    ${TomoEngine_SOURCE_DIR}/Filters/TargetGainSigmaXEstimation.h
)

cmp_IDE_SOURCE_PROPERTIES( "TomoEngine/Filters" "${TomoEngine_Filters_HDRS}" "${TomoEngine_Filters_SRCS}" "${CMP_INSTALL_FILES}")
