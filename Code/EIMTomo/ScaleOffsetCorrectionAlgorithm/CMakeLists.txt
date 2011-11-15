#--////////////////////////////////////////////////////////////////////////////
#--
#--  Copyright (c) 2011, Michael A. Jackson. BlueQuartz Software

#--  All rights reserved.
#--  BSD License: http://www.opensource.org/licenses/bsd-license.html
#--
#-- This code was partly written under US Air Force Contract FA8650-07-D-5800
#--
#--////////////////////////////////////////////////////////////////////////////
project (ScaleOffsetCorrectionAlgorithm)
cmake_minimum_required(VERSION 2.8.3)


set(ScaleOffsetCorrectionAlgorithm_HDRS
  ${PROJECT_CODE_DIR}/EIMTomo/ScaleOffsetCorrectionAlgorithm/ScaleOffsetCorrectionEngine.h
  ${PROJECT_CODE_DIR}/EIMTomo/ScaleOffsetCorrectionAlgorithm/ScaleOffsetCorrectionInputs.h
  ${PROJECT_CODE_DIR}/EIMTomo/ScaleOffsetCorrectionAlgorithm/ScaleOffsetStructures.h
)
set(ScaleOffsetCorrectionAlgorithm_SRCS
  ${PROJECT_CODE_DIR}/EIMTomo/ScaleOffsetCorrectionAlgorithm/ScaleOffsetCorrectionEngine.cpp
  ${PROJECT_CODE_DIR}/EIMTomo/ScaleOffsetCorrectionAlgorithm/ScaleOffsetCorrectionInputs.cpp
  ${PROJECT_CODE_DIR}/EIMTomo/ScaleOffsetCorrectionAlgorithm/ScaleOffsetCorrectionMain.cpp
)

cmp_IDE_SOURCE_PROPERTIES( "EIMTomo/ScaleOffsetCorrectionAlgorithm" "${ScaleOffsetCorrectionAlgorithm_HDRS}" "${ScaleOffsetCorrectionAlgorithm_SRCS}" "${CMP_INSTALL_FILES}")

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

set(ScaleOffsetCorrectionAlgorithm_PROJECT_FILES
  ${ScaleOffsetCorrectionAlgorithm_HDRS}
  ${ScaleOffsetCorrectionAlgorithm_SRCS}
  ${EIMTomo_Common_SRCS}
  ${EIMTomo_Common_HDRS}
  ${Mersenne_Twister_SRCS}
  ${Mersenne_Twister_HDRS}
  ${GENERATED_HDRS}
)

if ( ${EIMTomo_USE_MERSENNE_TWISTER} )
    set(ScaleOffsetCorrectionAlgorithm_PROJECT_FILES
      ${ScaleOffsetCorrectionAlgorithm_PROJECT_FILES}
      ${Mersenne_Twister_SRCS}
      ${Mersenne_Twister_HDRS}
    )
endif()

#add some additional include directories:
include_directories(${PROJECT_CODE_DIR})

add_executable(ScaleOffsetCorrectionAlgorithm ${ScaleOffsetCorrectionAlgorithm_PROJECT_FILES} )
target_link_libraries(ScaleOffsetCorrectionAlgorithm MXA)

configure_file(${EIMTomo_SOURCE_DIR}/Resources/InputFiles/RunningTheCode.sh.in
               ${EIMTomo_BINARY_DIR}/ScaleOffsetCorrectionAlgorithm.sh)
