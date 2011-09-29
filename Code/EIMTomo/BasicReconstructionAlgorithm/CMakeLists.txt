#--////////////////////////////////////////////////////////////////////////////
#--
#--  Copyright (c) 2011, Michael A. Jackson. BlueQuartz Software

#--  All rights reserved.
#--  BSD License: http://www.opensource.org/licenses/bsd-license.html
#--
#-- This code was partly written under US Air Force Contract FA8650-07-D-5800
#--
#--////////////////////////////////////////////////////////////////////////////
project (BasicReconstructionAlgorithm)
cmake_minimum_required(VERSION 2.8.3)


set(BasicReconstructionAlgorithm_HDRS
  ${PROJECT_CODE_DIR}/EIMTomo/BasicReconstructionAlgorithm/ComputationEngineOptimized.h
  ${PROJECT_CODE_DIR}/EIMTomo/BasicReconstructionAlgorithm/ComputationInputs.h
)
set(BasicReconstructionAlgorithm_SRCS
  ${PROJECT_CODE_DIR}/EIMTomo/BasicReconstructionAlgorithm/ComputationEngineOptimized.c
  ${PROJECT_CODE_DIR}/EIMTomo/BasicReconstructionAlgorithm/ComputationInputs.c
  ${PROJECT_CODE_DIR}/EIMTomo/BasicReconstructionAlgorithm/Main.c
)

cmp_IDE_SOURCE_PROPERTIES( "EIMTomo/BasicReconstructionAlgorithm" "${BasicReconstructionAlgorithm_HDRS}" "${BasicReconstructionAlgorithm_SRCS}" "${CMP_INSTALL_FILES}")

# Include the source files for Common source files and the Mersenne Twister RNG
include (${PROJECT_CODE_DIR}/EIMTomo/common/SourceList.cmake)
include (${PROJECT_CODE_DIR}/EIMTomo/mt/SourceList.cmake)

set(GENERATED_HDRS
    ${CMP_HEADER_DIR}/${CMP_EXTRA_CONFIGURATION_FILE}
    ${CMP_HEADER_DIR}/${CMP_CONFIGURATION_FILE_NAME}
    ${CMP_HEADER_DIR}/${CMP_TYPES_FILE_NAME}
    ${CMP_HEADER_DIR}/${CMP_VERSION_HEADER_FILE_NAME}
)

cmp_IDE_GENERATED_PROPERTIES("${Project_Name}/Generated" "${GENERATED_HDRS}" "")

set(BasicReconstructionAlgorithm_PROJECT_FILES
  ${BasicReconstructionAlgorithm_HDRS}
  ${BasicReconstructionAlgorithm_SRCS}
  ${EIMTomo_Common_SRCS}
  ${EIMTomo_Common_HDRS}
  ${Mersenne_Twister_SRCS}
  ${Mersenne_Twister_HDRS}
  ${GENERATED_HDRS}
)

#add some additional include directories:
include_directories(${PROJECT_CODE_DIR})

add_executable(BasicReconstructionAlgorithm ${BasicReconstructionAlgorithm_PROJECT_FILES} )

