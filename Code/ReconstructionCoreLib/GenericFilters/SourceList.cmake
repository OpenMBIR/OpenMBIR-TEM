#--////////////////////////////////////////////////////////////////////////////
#-- Copyright (c) 2009, Michael A. Jackson. BlueQuartz Software
#-- All rights reserved.
#-- BSD License: http://www.opensource.org/licenses/bsd-license.html
#-- This code was partly written under US Air Force Contract FA8650-07-D-5800
#--////////////////////////////////////////////////////////////////////////////
set (ReconstructionCoreLib_GenericFilters_SRCS
    ${ReconstructionCoreLib_SOURCE_DIR}/GenericFilters/CalculateAMatrixColumn.cpp
    ${ReconstructionCoreLib_SOURCE_DIR}/GenericFilters/CostData.cpp
    ${ReconstructionCoreLib_SOURCE_DIR}/GenericFilters/MRCSinogramInitializer.cpp
    ${ReconstructionCoreLib_SOURCE_DIR}/GenericFilters/RawSinogramInitializer.cpp
    
    ${ReconstructionCoreLib_SOURCE_DIR}/GenericFilters/InitialReconstructionInitializer.cpp
    ${ReconstructionCoreLib_SOURCE_DIR}/GenericFilters/InitialReconstructionBinReader.cpp
    ${ReconstructionCoreLib_SOURCE_DIR}/GenericFilters/TomoFilter.cpp
    ${ReconstructionCoreLib_SOURCE_DIR}/GenericFilters/ComputeInitialOffsets.cpp
    ${ReconstructionCoreLib_SOURCE_DIR}/GenericFilters/DetectorResponse.cpp

    ${ReconstructionCoreLib_SOURCE_DIR}/GenericFilters/SigmaXEstimation.cpp
)

set (ReconstructionCoreLib_GenericFilters_HDRS
    ${ReconstructionCoreLib_SOURCE_DIR}/GenericFilters/CalculateAMatrixColumn.h
    ${ReconstructionCoreLib_SOURCE_DIR}/GenericFilters/CostData.h
    ${ReconstructionCoreLib_SOURCE_DIR}/GenericFilters/MRCSinogramInitializer.h
    ${ReconstructionCoreLib_SOURCE_DIR}/GenericFilters/RawSinogramInitializer.h
    ${ReconstructionCoreLib_SOURCE_DIR}/GenericFilters/InitialReconstructionInitializer.h
    ${ReconstructionCoreLib_SOURCE_DIR}/GenericFilters/InitialReconstructionBinReader.h
    ${ReconstructionCoreLib_SOURCE_DIR}/GenericFilters/TomoFilter.h
    ${ReconstructionCoreLib_SOURCE_DIR}/GenericFilters/ComputeInitialOffsets.h 
    ${ReconstructionCoreLib_SOURCE_DIR}/GenericFilters/DetectorResponse.h
    ${ReconstructionCoreLib_SOURCE_DIR}/GenericFilters/SigmaXEstimation.h
)

cmp_IDE_SOURCE_PROPERTIES( "ReconstructionCoreLib/GenericFilters" "${ReconstructionCoreLib_GenericFilters_HDRS}" "${ReconstructionCoreLib_GenericFilters_SRCS}" "${CMP_INSTALL_FILES}")
