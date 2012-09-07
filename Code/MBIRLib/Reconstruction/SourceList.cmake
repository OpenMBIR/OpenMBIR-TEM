#--////////////////////////////////////////////////////////////////////////////
#--
#--  Copyright (c) 2011, Michael A. Jackson. BlueQuartz Software

#--  All rights reserved.
#--  BSD License: http://www.opensource.org/licenses/bsd-license.html
#--
#-- This code was partly written under US Air Force Contract FA8650-07-D-5800
#--
#--////////////////////////////////////////////////////////////////////////////


set (MBIRLib_Reconstruction_SRCS
    ${MBIRLib_SOURCE_DIR}/Reconstruction/ReconstructionEngine.cpp
    ${MBIRLib_SOURCE_DIR}/Reconstruction/ReconstructionEngine_UpdateVoxels.cpp
    ${MBIRLib_SOURCE_DIR}/Reconstruction/ReconstructionEngine_Extra.cpp
    ${MBIRLib_SOURCE_DIR}/Reconstruction/ReconstructionInputs.cpp
    ${MBIRLib_SOURCE_DIR}/Reconstruction/MultiResolutionReconstruction.cpp
)

set (MBIRLib_Reconstruction_HDRS
    ${MBIRLib_SOURCE_DIR}/Reconstruction/ReconstructionConstants.h
    ${MBIRLib_SOURCE_DIR}/Reconstruction/ReconstructionEngine.h
    ${MBIRLib_SOURCE_DIR}/Reconstruction/ReconstructionInputs.h
    ${MBIRLib_SOURCE_DIR}/Reconstruction/ReconstructionStructures.h
    ${MBIRLib_SOURCE_DIR}/Reconstruction/MultiResolutionReconstruction.h
)
set_source_files_properties( ${MBIRLib_SOURCE_DIR}/Reconstruction/ReconstructionEngine_Extra.cpp
                             PROPERTIES HEADER_FILE_ONLY TRUE)
cmp_IDE_SOURCE_PROPERTIES( "MBIRLib/Reconstruction" "${MBIRLib_Reconstruction_HDRS}" "${MBIRLib_Reconstruction_SRCS}" "${CMP_INSTALL_FILES}")
