#--////////////////////////////////////////////////////////////////////////////
#--
#--  Copyright (c) 2011, Michael A. Jackson. BlueQuartz Software

#--  All rights reserved.
#--  BSD License: http://www.opensource.org/licenses/bsd-license.html
#--
#-- This code was partly written under US Air Force Contract FA8650-07-D-5800
#--
#--////////////////////////////////////////////////////////////////////////////
project (ScaleOffsetMotionCorrectionAlgorithm)
cmake_minimum_required(VERSION 2.8.3)


set(ScaleOffsetMotionCorrectionAlgorithm_HDRS
  ${PROJECT_CODE_DIR}/EIMTomo/ScaleOffsetMotionCorrectionAlgorithm/ScaleOffsetMotionCorrectionEngine.h
  ${PROJECT_CODE_DIR}/EIMTomo/ScaleOffsetMotionCorrectionAlgorithm/ScaleOffsetMotionCorrectionParser.h
  ${PROJECT_CODE_DIR}/EIMTomo/ScaleOffsetMotionCorrectionAlgorithm/ScaleOffsetMotionStructures.h
  ${PROJECT_CODE_DIR}/EIMTomo/ScaleOffsetMotionCorrectionAlgorithm/ScaleOffsetMotionCorrectionConstants.h
)
set(ScaleOffsetMotionCorrectionAlgorithm_SRCS
  ${PROJECT_CODE_DIR}/EIMTomo/ScaleOffsetMotionCorrectionAlgorithm/ScaleOffsetMotionCorrectionEngine.cpp
  ${PROJECT_CODE_DIR}/EIMTomo/ScaleOffsetMotionCorrectionAlgorithm/ScaleOffsetMotionCorrectionParser.cpp
  ${PROJECT_CODE_DIR}/EIMTomo/ScaleOffsetMotionCorrectionAlgorithm/ScaleOffsetMotionCorrectionMain.cpp
)

cmp_IDE_SOURCE_PROPERTIES( "EIMTomo/ScaleOffsetMotionCorrectionAlgorithm" "${ScaleOffsetMotionCorrectionAlgorithm_HDRS}" "${ScaleOffsetMotionCorrectionAlgorithm_SRCS}" "${CMP_INSTALL_FILES}")

# Include the source files for Common source files and the Mersenne Twister RNG
#include (${PROJECT_CODE_DIR}/EIMTomo/common/SourceList.cmake)
#include (${PROJECT_CODE_DIR}/EIMTomo/mt/SourceList.cmake)
#include (${PROJECT_CODE_DIR}/EIMTomo/IO/SourceList.cmake)

set(GENERATED_HDRS
    ${CMP_HEADER_DIR}/${CMP_EXTRA_CONFIGURATION_FILE}
    ${CMP_HEADER_DIR}/${CMP_CONFIGURATION_FILE_NAME}
    ${CMP_HEADER_DIR}/${CMP_TYPES_FILE_NAME}
    ${CMP_HEADER_DIR}/${CMP_VERSION_HEADER_FILE_NAME}
)

cmp_IDE_GENERATED_PROPERTIES("${Project_Name}/Generated" "${GENERATED_HDRS}" "")

set(ScaleOffsetMotionCorrectionAlgorithm_PROJECT_FILES
  ${ScaleOffsetMotionCorrectionAlgorithm_HDRS}
  ${ScaleOffsetMotionCorrectionAlgorithm_SRCS}
  ${GENERATED_HDRS}
)

if ( ${EIMTomo_USE_MERSENNE_TWISTER} )
    set(ScaleOffsetMotionCorrectionAlgorithm_PROJECT_FILES
      ${ScaleOffsetMotionCorrectionAlgorithm_PROJECT_FILES}
    )
endif()

#add some additional include directories:
include_directories(${PROJECT_CODE_DIR})

add_executable(ScaleOffsetMotionCorrectionAlgorithm ${ScaleOffsetMotionCorrectionAlgorithm_PROJECT_FILES} )
target_link_libraries(ScaleOffsetMotionCorrectionAlgorithm MXA TomoEngine)

configure_file(${EIMTomo_SOURCE_DIR}/Resources/InputFiles/RunningTheCode.sh.in
               ${EIMTomo_BINARY_DIR}/ScaleOffsetMotionCorrectionAlgorithm.sh)
