#--////////////////////////////////////////////////////////////////////////////
#-- Copyright (c) 2009, Michael A. Jackson. BlueQuartz Software
#-- All rights reserved.
#-- BSD License: http://www.opensource.org/licenses/bsd-license.html
#-- This code was partly written under US Air Force Contract FA8650-07-D-5800
#--////////////////////////////////////////////////////////////////////////////
set (EIMTomo_Common_SRCS
    ${PROJECT_CODE_DIR}/EIMTomo/common/allocate.c
  #  ${PROJECT_CODE_DIR}/EIMTomo/common/randlib.c
)

set (EIMTomo_Common_HDRS
    ${PROJECT_CODE_DIR}/EIMTomo/common/allocate.h
  #  ${PROJECT_CODE_DIR}/EIMTomo/common/randlib.h
)
cmp_IDE_SOURCE_PROPERTIES( "EIMTomo/common" "${EIMTomo_Common_HDRS}" "${EIMTomo_Common_SRCS}" "${CMP_INSTALL_FILES}")
