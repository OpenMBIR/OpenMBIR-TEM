#--////////////////////////////////////////////////////////////////////////////
#-- Copyright (c) 2009, Michael A. Jackson. BlueQuartz Software
#-- All rights reserved.
#-- BSD License: http://www.opensource.org/licenses/bsd-license.html
#-- This code was partly written under US Air Force Contract FA8650-07-D-5800
#--////////////////////////////////////////////////////////////////////////////
set (TomoEngine_Common_SRCS
    ${TomoEngine_SOURCE_DIR}/Common/allocate.c
    ${TomoEngine_SOURCE_DIR}/Common/EIMTime.c
    ${TomoEngine_SOURCE_DIR}/Common/EIMImage.cpp
    ${TomoEngine_SOURCE_DIR}/Common/AbstractFilter.cpp
    ${TomoEngine_SOURCE_DIR}/Common/FilterPipeline.cpp
    ${TomoEngine_SOURCE_DIR}/Common/Observer.cpp
    ${TomoEngine_SOURCE_DIR}/Common/Observable.cpp
)

set (TomoEngine_Common_HDRS
    ${TomoEngine_SOURCE_DIR}/Common/allocate.h
    ${TomoEngine_SOURCE_DIR}/Common/TomoEngineDLLExport.h
    ${TomoEngine_SOURCE_DIR}/Common/MSVCDefines.h
    ${TomoEngine_SOURCE_DIR}/Common/EIMImage.h
    ${TomoEngine_SOURCE_DIR}/Common/EIMTime.h
    ${TomoEngine_SOURCE_DIR}/Common/EIMMath.h
    ${TomoEngine_SOURCE_DIR}/Common/AbstractFilter.h
    ${TomoEngine_SOURCE_DIR}/Common/FilterPipeline.h
    ${TomoEngine_SOURCE_DIR}/Common/Observer.h
    ${TomoEngine_SOURCE_DIR}/Common/Observable.h
    ${TomoEngine_SOURCE_DIR}/Common/CE_ConstraintEquation.hpp
    ${TomoEngine_SOURCE_DIR}/Common/DerivOfCostFunc.hpp
    ${TomoEngine_SOURCE_DIR}/Common/TomoArray.hpp
)

cmp_IDE_SOURCE_PROPERTIES( "TomoEngine/Common" "${TomoEngine_Common_HDRS}" "${TomoEngine_Common_SRCS}" "${CMP_INSTALL_FILES}")
