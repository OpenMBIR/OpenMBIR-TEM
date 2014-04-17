#--////////////////////////////////////////////////////////////////////////////
#--
#--  Copyright (c) 2011, Michael A. Jackson. BlueQuartz Software

#--  All rights reserved.
#--  BSD License: http://www.opensource.org/licenses/bsd-license.html
#--
#-- This code was partly written under US Air Force Contract FA8650-07-D-5800
#--
#--////////////////////////////////////////////////////////////////////////////


set (MBIRLib_HAADF_SRCS
  ${MBIRLib_SOURCE_DIR}/HAADF/HAADF_QGGMRFPriorModel.cpp
  ${MBIRLib_SOURCE_DIR}/HAADF/HAADF_MultiResolutionReconstruction.cpp
  ${MBIRLib_SOURCE_DIR}/HAADF/HAADF_ForwardProject.cpp
  ${MBIRLib_SOURCE_DIR}/HAADF/HAADF_ReconstructionEngine.cpp
  ${MBIRLib_SOURCE_DIR}/HAADF/HAADF_ReconstructionEngine_UpdateVoxels.cpp
  ${MBIRLib_SOURCE_DIR}/HAADF/HAADF_ReconstructionEngine_Extra.cpp
)

set (MBIRLib_HAADF_HDRS
    ${MBIRLib_SOURCE_DIR}/HAADF/HAADFConstants.h
    ${MBIRLib_SOURCE_DIR}/HAADF/HAADF_QGGMRFPriorModel.h
    ${MBIRLib_SOURCE_DIR}/HAADF/HAADF_MultiResolutionReconstruction.h
    ${MBIRLib_SOURCE_DIR}/HAADF/HAADF_ForwardProject.h
    ${MBIRLib_SOURCE_DIR}/HAADF/HAADF_ReconstructionEngine.h
)

set_source_files_properties( ${MBIRLib_SOURCE_DIR}/HAADF/HAADF_ReconstructionEngine_UpdateVoxels.cpp
                             ${MBIRLib_SOURCE_DIR}/HAADF/HAADF_ReconstructionEngine_Extra.cpp
                             PROPERTIES HEADER_FILE_ONLY TRUE)

cmp_IDE_SOURCE_PROPERTIES( "MBIRLib/HAADF" "${MBIRLib_HAADF_HDRS}" "${MBIRLib_HAADF_SRCS}" "${CMP_INSTALL_FILES}")
