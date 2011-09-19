#--////////////////////////////////////////////////////////////////////////////
#--
#--  Copyright (c) 2011, Michael A. Jackson. BlueQuartz Software

#--  All rights reserved.
#--  BSD License: http://www.opensource.org/licenses/bsd-license.html
#--
#-- This code was partly written under US Air Force Contract FA8650-07-D-5800
#--
#--////////////////////////////////////////////////////////////////////////////
project (JointEstimation)
cmake_minimum_required(VERSION 2.8.3)

#Define the headers for this program

set(JointEstimation_HDRS
  ${PROJECT_CODE_DIR}/EIMTomo/JointEstimationAlgorithm/ComputationEngineMotionCorrection.h
  ${PROJECT_CODE_DIR}/EIMTomo/JointEstimationAlgorithm/ComputationInputsMotionCorrection.h
)

#Define the source files for this program
set(JointEstimation_SRCS
  ${PROJECT_CODE_DIR}/EIMTomo/JointEstimationAlgorithm/ComputationEngineMotionCorrection.c
  ${PROJECT_CODE_DIR}/EIMTomo/JointEstimationAlgorithm/ComputationInputsMotionCorrection.c
  ${PROJECT_CODE_DIR}/EIMTomo/JointEstimationAlgorithm/JointEstimationMain.c
)

cmp_IDE_SOURCE_PROPERTIES( "EIMTomo/JointEstimation" "${JointEstimation_HDRS}" "${JointEstimation_SRCS}" "${CMP_INSTALL_FILES}")


# Include the source files for Common source files and the Mersenne Twister RNG
include (${PROJECT_CODE_DIR}/EIMTomo/Common/SourceList.cmake)
include (${PROJECT_CODE_DIR}/EIMTomo/mt/SourceList.cmake)


# Gather all the source files into a single CMake Variable
set(JointEstimation_PROJECT_FILES
  ${JointEstimation_HDRS}
  ${JointEstimation_SRCS}
  ${EIMTomo_Common_SRCS}
  ${EIMTomo_Common_HDRS}
  ${Mersenne_Twister_SRCS}
  ${Mersenne_Twister_HDRS}
)

#add some additional include directories:
include_directories(${PROJECT_CODE_DIR})

# Create the executable.
add_executable(JointEstimation ${JointEstimation_PROJECT_FILES} )

