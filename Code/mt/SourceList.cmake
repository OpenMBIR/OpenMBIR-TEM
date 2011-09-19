#--////////////////////////////////////////////////////////////////////////////
#-- Copyright (c) 2009, Michael A. Jackson. BlueQuartz Software
#-- All rights reserved.
#-- BSD License: http://www.opensource.org/licenses/bsd-license.html
#-- This code was partly written under US Air Force Contract FA8650-07-D-5800
#--////////////////////////////////////////////////////////////////////////////
set (Mersenne_Twister_SRCS
    ${PROJECT_CODE_DIR}/mt/mt19937ar.c
)

set (Mersenne_Twister_HDRS
    ${PROJECT_CODE_DIR}/mt/mt19937ar.h
)
#cmp_IDE_SOURCE_PROPERTIES( "mt" "${Mersenne_Twister_HDRS}" "${Mersenne_Twister_SRCS}" "${CMP_INSTALL_FILES}")
