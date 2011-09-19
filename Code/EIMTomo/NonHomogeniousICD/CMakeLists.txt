#--////////////////////////////////////////////////////////////////////////////
#--
#--  Copyright (c) 2011, Michael A. Jackson. BlueQuartz Software

#--  All rights reserved.
#--  BSD License: http://www.opensource.org/licenses/bsd-license.html
#--
#-- This code was partly written under US Air Force Contract FA8650-07-D-5800
#--
#--////////////////////////////////////////////////////////////////////////////
project (NonHomogeniousICD)
cmake_minimum_required(VERSION 2.8.3)


set(NonHomogeniousICD_HDRS
  ${PROJECT_CODE_DIR}/EIMTomo/NonHomogeniousICD/ComputationEngineOptimizedNHICD.h
  ${PROJECT_CODE_DIR}/EIMTomo/NonHomogeniousICD/ComputationInputsNHICD.h
)
set(NonHomogeniousICD_SRCS
  ${PROJECT_CODE_DIR}/EIMTomo/NonHomogeniousICD/ComputationEngineOptimizedNHICD.c
  ${PROJECT_CODE_DIR}/EIMTomo/NonHomogeniousICD/ComputationInputsNHICD.c
  ${PROJECT_CODE_DIR}/EIMTomo/NonHomogeniousICD/NHICDMain.c
)

cmp_IDE_SOURCE_PROPERTIES( "EIMTomo/NHICD" "${NonHomogeniousICD_HDRS}" "${NonHomogeniousICD_SRCS}" "${CMP_INSTALL_FILES}")

# Include the source files for Common source files and the Mersenne Twister RNG
include (${PROJECT_CODE_DIR}/EIMTomo/Common/SourceList.cmake)
include (${PROJECT_CODE_DIR}/EIMTomo/mt/SourceList.cmake)

set(NonHomogeniousICD_PROJECT_FILES
  ${NonHomogeniousICD_HDRS}
  ${NonHomogeniousICD_SRCS}
    ${EIMTomo_Common_SRCS}
  ${EIMTomo_Common_HDRS}
  ${Mersenne_Twister_SRCS}
  ${Mersenne_Twister_HDRS}
)

#add some additional include directories:
include_directories(${PROJECT_CODE_DIR})

add_executable(NHICD ${NonHomogeniousICD_PROJECT_FILES} )

