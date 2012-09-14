#--////////////////////////////////////////////////////////////////////////////
#-- Copyright (c) 2009, Michael A. Jackson. BlueQuartz Software
#-- All rights reserved.
#-- BSD License: http://www.opensource.org/licenses/bsd-license.html
#-- This code was partly written under US Air Force Contract FA8650-07-D-5800
#--////////////////////////////////////////////////////////////////////////////
set (MBIRLib_HAADF_SRCS

    ${MBIRLib_SOURCE_DIR}/HAADF/HAADFForwardModel.cpp
    ${MBIRLib_SOURCE_DIR}/HAADF/HAADFAMatrixCol.cpp
    ${MBIRLib_SOURCE_DIR}/HAADF/Filters/CalculateAMatrixColumn.cpp
    ${MBIRLib_SOURCE_DIR}/HAADF/Filters/ComputeInitialOffsets.cpp
    ${MBIRLib_SOURCE_DIR}/HAADF/Filters/GainsOffsetsReader.cpp
    ${MBIRLib_SOURCE_DIR}/HAADF/Filters/NuisanceParamWriter.cpp
    ${MBIRLib_SOURCE_DIR}/HAADF/Filters/NuisanceParamReader.cpp
    ${MBIRLib_SOURCE_DIR}/HAADF/Filters/SinogramBinWriter.cpp
    ${MBIRLib_SOURCE_DIR}/HAADF/ForwardProject.cpp
    ${MBIRLib_SOURCE_DIR}/HAADF/QGGMRFPriorModel.cpp
    ${MBIRLib_SOURCE_DIR}/HAADF/UpdateYSlice.cpp
    ${MBIRLib_SOURCE_DIR}/HAADF/HAADFDetectorParameters.cpp
)

set (MBIRLib_HAADF_HDRS
    ${MBIRLib_SOURCE_DIR}/HAADF/Filters/CalculateAMatrixColumn.h
    ${MBIRLib_SOURCE_DIR}/HAADF/Filters/ComputeInitialOffsets.h 
    ${MBIRLib_SOURCE_DIR}/HAADF/Filters/GainsOffsetsReader.h
    ${MBIRLib_SOURCE_DIR}/HAADF/Filters/NuisanceParamWriter.h
    ${MBIRLib_SOURCE_DIR}/HAADF/Filters/NuisanceParamReader.h
    ${MBIRLib_SOURCE_DIR}/HAADF/Filters/SinogramBinWriter.h
    ${MBIRLib_SOURCE_DIR}/HAADF/HAADFForwardModel.h
    ${MBIRLib_SOURCE_DIR}/HAADF/HAADFAMatrixCol.h
    ${MBIRLib_SOURCE_DIR}/HAADF/ForwardProject.h
    ${MBIRLib_SOURCE_DIR}/HAADF/QGGMRFPriorModel.h
    ${MBIRLib_SOURCE_DIR}/HAADF/UpdateYSlice.h
    ${MBIRLib_SOURCE_DIR}/HAADF/HAADFDetectorParameters.h
)

cmp_IDE_SOURCE_PROPERTIES( "MBIRLib/HAADF" "${MBIRLib_HAADF_HDRS}" "${MBIRLib_HAADF_SRCS}" "${CMP_INSTALL_FILES}")
