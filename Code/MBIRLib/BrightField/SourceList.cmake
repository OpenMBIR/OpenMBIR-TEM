#--////////////////////////////////////////////////////////////////////////////
#-- Copyright (c) 2009, Michael A. Jackson. BlueQuartz Software
#-- All rights reserved.
#-- BSD License: http://www.opensource.org/licenses/bsd-license.html
#-- This code was partly written under US Air Force Contract FA8650-07-D-5800
#--////////////////////////////////////////////////////////////////////////////

set (MBIRLib_BrightField_HDRS
    ${MBIRLib_SOURCE_DIR}/BrightField/BFForwardModel.h
    ${MBIRLib_SOURCE_DIR}/BrightField/BFAMatrixCol.h
    ${MBIRLib_SOURCE_DIR}/BrightField/BFForwardProject.h
    ${MBIRLib_SOURCE_DIR}/BrightField/BFQGGMRFPriorModel.h
    ${MBIRLib_SOURCE_DIR}/BrightField/BFUpdateYSlice.h
    ${MBIRLib_SOURCE_DIR}/BrightField/BFDetectorParameters.h
)

set (MBIRLib_BrightField_SRCS
    ${MBIRLib_SOURCE_DIR}/BrightField/BFForwardModel.cpp
    ${MBIRLib_SOURCE_DIR}/BrightField/BFAMatrixCol.cpp
    ${MBIRLib_SOURCE_DIR}/BrightField/BFForwardProject.cpp
    ${MBIRLib_SOURCE_DIR}/BrightField/BFQGGMRFPriorModel.cpp
    ${MBIRLib_SOURCE_DIR}/BrightField/BFUpdateYSlice.cpp
    ${MBIRLib_SOURCE_DIR}/BrightField/BFDetectorParameters.cpp
)

cmp_IDE_SOURCE_PROPERTIES( "MBIRLib/BrightField" "${MBIRLib_BrightField_HDRS}" "${MBIRLib_BrightField_SRCS}" "${CMP_INSTALL_FILES}")


set (MBIRLib_BF_FILTERS_HDRS
    ${MBIRLib_SOURCE_DIR}/BrightField/Filters/CalculateAMatrixColumn.h
    ${MBIRLib_SOURCE_DIR}/BrightField/Filters/ComputeInitialOffsets.h 
    ${MBIRLib_SOURCE_DIR}/BrightField/Filters/GainsOffsetsReader.h
    ${MBIRLib_SOURCE_DIR}/BrightField/Filters/NuisanceParamWriter.h
    ${MBIRLib_SOURCE_DIR}/BrightField/Filters/NuisanceParamReader.h
    ${MBIRLib_SOURCE_DIR}/BrightField/Filters/SinogramBinWriter.h
)

set (MBIRLib_BF_FILTERS_SRCS
    ${MBIRLib_SOURCE_DIR}/BrightField/Filters/CalculateAMatrixColumn.cpp
    ${MBIRLib_SOURCE_DIR}/BrightField/Filters/ComputeInitialOffsets.cpp
    ${MBIRLib_SOURCE_DIR}/BrightField/Filters/GainsOffsetsReader.cpp
    ${MBIRLib_SOURCE_DIR}/BrightField/Filters/NuisanceParamWriter.cpp
    ${MBIRLib_SOURCE_DIR}/BrightField/Filters/NuisanceParamReader.cpp
    ${MBIRLib_SOURCE_DIR}/BrightField/Filters/SinogramBinWriter.cpp
)

cmp_IDE_SOURCE_PROPERTIES( "MBIRLib/BrightField/Filters" "${MBIRLib_BF_FILTERS_HDRS}" "${MBIRLib_BF_FILTERS_SRCS}" "${CMP_INSTALL_FILES}")

set(MBIRLib_BrightField_HDRS
    ${MBIRLib_BrightField_HDRS}
    ${MBIRLib_BF_FILTERS_HDRS}
    )
    
set(MBIRLib_BrightField_SRCS
    ${MBIRLib_BrightField_SRCS}
    ${MBIRLib_BF_FILTERS_SRCS})
    