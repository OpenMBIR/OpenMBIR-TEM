#--////////////////////////////////////////////////////////////////////////////
#--
#--  Copyright (c) 2011, Michael A. Jackson. BlueQuartz Software

#--  All rights reserved.
#--  BSD License: http://www.opensource.org/licenses/bsd-license.html
#--
#-- This code was partly written under US Air Force Contract FA8650-07-D-5800
#--
#--////////////////////////////////////////////////////////////////////////////
project (MotionCorrection)
cmake_minimum_required(VERSION 2.8.3)

#Define the headers for this program

set(MotionCorrection_HDRS
  ${PROJECT_CODE_DIR}/EIMTomo/MotionCorrectionAlgorithm/MotionCorrectionEngine.h
  ${PROJECT_CODE_DIR}/EIMTomo/MotionCorrectionAlgorithm/MotionCorrectionInputs.h
  ${PROJECT_CODE_DIR}/EIMTomo/MotionCorrectionAlgorithm/MotionCorrectionStructures.h
)

#Define the source files for this program
set(MotionCorrection_SRCS
  ${PROJECT_CODE_DIR}/EIMTomo/MotionCorrectionAlgorithm/MotionCorrectionEngine.cpp
  ${PROJECT_CODE_DIR}/EIMTomo/MotionCorrectionAlgorithm/MotionCorrectionInputs.cpp
  ${PROJECT_CODE_DIR}/EIMTomo/MotionCorrectionAlgorithm/MotionCorrectionMain.cpp
)

cmp_IDE_SOURCE_PROPERTIES( "EIMTomo/MotionCorrection" "${MotionCorrection_HDRS}" "${MotionCorrection_SRCS}" "${CMP_INSTALL_FILES}")


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

# Gather all the source files into a single CMake Variable
set(MotionCorrection_PROJECT_FILES
  ${MotionCorrection_HDRS}
  ${MotionCorrection_SRCS}
  ${EIMTomo_Common_SRCS}
  ${EIMTomo_Common_HDRS}
  ${Mersenne_Twister_SRCS}
  ${Mersenne_Twister_HDRS}
  ${GENERATED_HDRS}
)

#Add Some Definitions



#add some additional include directories:
include_directories(${PROJECT_CODE_DIR})

# Create the executable.
add_executable(MotionCorrection ${MotionCorrection_PROJECT_FILES} )

