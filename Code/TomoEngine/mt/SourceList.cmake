#--////////////////////////////////////////////////////////////////////////////
#-- Copyright (c) 2009, Michael A. Jackson. BlueQuartz Software
#-- All rights reserved.
#-- BSD License: http://www.opensource.org/licenses/bsd-license.html
#-- This code was partly written under US Air Force Contract FA8650-07-D-5800
#--////////////////////////////////////////////////////////////////////////////
set (TomoEngine_mt_SRCS
    ${TomoEngine_SOURCE_DIR}/mt/mt19937ar.c
)

set (TomoEngine_mt_HDRS
    ${TomoEngine_SOURCE_DIR}/mt/mt19937ar.h
)
cmp_IDE_SOURCE_PROPERTIES( "mt" "${TomoEngine_mt_HDRS}" "${TomoEngine_mt_SRCS}" "${CMP_INSTALL_FILES}")
