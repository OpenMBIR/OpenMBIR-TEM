#--////////////////////////////////////////////////////////////////////////////
#-- Copyright (c) 2011, Michael A. Jackson. BlueQuartz Software
#-- All rights reserved.
#-- BSD License: http://www.opensource.org/licenses/bsd-license.html
#-- This code was partly written under US Air Force Contract FA8650-07-D-5800
#--////////////////////////////////////////////////////////////////////////////
set (TomoEngine_SOC_SRCS
    ${TomoEngine_SOURCE_DIR}/SOC/SOCEngine.cpp
    ${TomoEngine_SOURCE_DIR}/SOC/SOCInputs.cpp
)

set (TomoEngine_SOC_HDRS
    ${TomoEngine_SOURCE_DIR}/SOC/SOCEngine.h
    ${TomoEngine_SOURCE_DIR}/SOC/SOCInputs.h
)
cmp_IDE_SOURCE_PROPERTIES( "SOC" "${TomoEngine_SOC_HDRS}" "${TomoEngine_SOC_SRCS}" "${CMP_INSTALL_FILES}")

