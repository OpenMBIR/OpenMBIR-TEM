#--////////////////////////////////////////////////////////////////////////////
#-- Copyright (c) 2009, Michael A. Jackson. BlueQuartz Software
#-- All rights reserved.
#-- BSD License: http://www.opensource.org/licenses/bsd-license.html
#-- This code was partly written under US Air Force Contract FA8650-07-D-5800
#--////////////////////////////////////////////////////////////////////////////

set (MBIRLib_BrightField_HDRS
    ${MBIRLib_SOURCE_DIR}/BrightField/BFForwardModel.h
    ${MBIRLib_SOURCE_DIR}/BrightField/BFForwardProject.h
    ${MBIRLib_SOURCE_DIR}/BrightField/BF_QGGMRFPriorModel.h
    ${MBIRLib_SOURCE_DIR}/BrightField/BFUpdateYSlice.h
    ${MBIRLib_SOURCE_DIR}/BrightField/BFReconstructionEngine.h
    ${MBIRLib_SOURCE_DIR}/BrightField/BFMultiResolutionReconstruction.h
    ${MBIRLib_SOURCE_DIR}/BrightField/BFConstants.h
)

set (MBIRLib_BrightField_SRCS
    ${MBIRLib_SOURCE_DIR}/BrightField/BFForwardModel.cpp
    ${MBIRLib_SOURCE_DIR}/BrightField/BFForwardProject.cpp
    ${MBIRLib_SOURCE_DIR}/BrightField/BF_QGGMRFPriorModel.cpp
    ${MBIRLib_SOURCE_DIR}/BrightField/BFUpdateYSlice.cpp
    ${MBIRLib_SOURCE_DIR}/BrightField/BFReconstructionEngine.cpp
    ${MBIRLib_SOURCE_DIR}/BrightField/BFReconstructionEngine_Extra.cpp
    ${MBIRLib_SOURCE_DIR}/BrightField/BFMultiResolutionReconstruction.cpp
)

set_source_files_properties( ${MBIRLib_SOURCE_DIR}/BrightField/BFReconstructionEngine_Extra.cpp
                             PROPERTIES HEADER_FILE_ONLY TRUE)

cmp_IDE_SOURCE_PROPERTIES( "MBIRLib/BrightField" "${MBIRLib_BrightField_HDRS}" "${MBIRLib_BrightField_SRCS}" "${CMP_INSTALL_FILES}")
