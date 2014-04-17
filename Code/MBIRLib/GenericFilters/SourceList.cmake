#--////////////////////////////////////////////////////////////////////////////
#-- Copyright (c) 2009, Michael A. Jackson. BlueQuartz Software
#-- All rights reserved.
#-- BSD License: http://www.opensource.org/licenses/bsd-license.html
#-- This code was partly written under US Air Force Contract FA8650-07-D-5800
#--////////////////////////////////////////////////////////////////////////////
set (MBIRLib_GenericFilters_SRCS
    ${MBIRLib_SOURCE_DIR}/GenericFilters/CostData.cpp
    ${MBIRLib_SOURCE_DIR}/GenericFilters/MRCSinogramInitializer.cpp
    ${MBIRLib_SOURCE_DIR}/GenericFilters/RawSinogramInitializer.cpp
    ${MBIRLib_SOURCE_DIR}/GenericFilters/InitialReconstructionInitializer.cpp
    ${MBIRLib_SOURCE_DIR}/GenericFilters/InitialReconstructionBinReader.cpp
    ${MBIRLib_SOURCE_DIR}/GenericFilters/TomoFilter.cpp
    ${MBIRLib_SOURCE_DIR}/GenericFilters/DetectorResponse.cpp
    ${MBIRLib_SOURCE_DIR}/GenericFilters/SigmaXEstimation.cpp
    ${MBIRLib_SOURCE_DIR}/GenericFilters/BackgroundCalculation.cpp
    ${MBIRLib_SOURCE_DIR}/GenericFilters/DetectorParameters.cpp
)

set (MBIRLib_GenericFilters_HDRS
    ${MBIRLib_SOURCE_DIR}/GenericFilters/CostData.h
    ${MBIRLib_SOURCE_DIR}/GenericFilters/MRCSinogramInitializer.h
    ${MBIRLib_SOURCE_DIR}/GenericFilters/RawSinogramInitializer.h
    ${MBIRLib_SOURCE_DIR}/GenericFilters/InitialReconstructionInitializer.h
    ${MBIRLib_SOURCE_DIR}/GenericFilters/InitialReconstructionBinReader.h
    ${MBIRLib_SOURCE_DIR}/GenericFilters/TomoFilter.h
    ${MBIRLib_SOURCE_DIR}/GenericFilters/DetectorResponse.h
    ${MBIRLib_SOURCE_DIR}/GenericFilters/SigmaXEstimation.h
    ${MBIRLib_SOURCE_DIR}/GenericFilters/BackgroundCalculation.h
    ${MBIRLib_SOURCE_DIR}/GenericFilters/DetectorParameters.h
)

cmp_IDE_SOURCE_PROPERTIES( "MBIRLib/GenericFilters" "${MBIRLib_GenericFilters_HDRS}" "${MBIRLib_GenericFilters_SRCS}" "${CMP_INSTALL_FILES}")
