#--////////////////////////////////////////////////////////////////////////////
#-- Copyright (c) 2009, Michael A. Jackson. BlueQuartz Software
#-- All rights reserved.
#-- BSD License: http://www.opensource.org/licenses/bsd-license.html
#-- This code was partly written under US Air Force Contract FA8650-07-D-5800
#--////////////////////////////////////////////////////////////////////////////
set (MBIRLib_Common_SRCS
    ${MBIRLib_SOURCE_DIR}/Common/allocate.c
    ${MBIRLib_SOURCE_DIR}/Common/EIMTime.c
    ${MBIRLib_SOURCE_DIR}/Common/EIMImage.cpp
    ${MBIRLib_SOURCE_DIR}/Common/AbstractFilter.cpp
    ${MBIRLib_SOURCE_DIR}/Common/FilterPipeline.cpp
    ${MBIRLib_SOURCE_DIR}/Common/Observer.cpp
    ${MBIRLib_SOURCE_DIR}/Common/Observable.cpp
    ${MBIRLib_SOURCE_DIR}/Common/VoxelUpdateList.cpp
)

set (MBIRLib_Common_HDRS
    ${MBIRLib_SOURCE_DIR}/Common/allocate.h
    ${MBIRLib_SOURCE_DIR}/Common/MBIRLibDLLExport.h
    ${MBIRLib_SOURCE_DIR}/Common/MSVCDefines.h
    ${MBIRLib_SOURCE_DIR}/Common/EIMImage.h
    ${MBIRLib_SOURCE_DIR}/Common/EIMTime.h
    ${MBIRLib_SOURCE_DIR}/Common/EIMMath.h
    ${MBIRLib_SOURCE_DIR}/Common/AbstractFilter.h
    ${MBIRLib_SOURCE_DIR}/Common/FilterPipeline.h
    ${MBIRLib_SOURCE_DIR}/Common/Observer.h
    ${MBIRLib_SOURCE_DIR}/Common/Observable.h
    ${MBIRLib_SOURCE_DIR}/Common/CE_ConstraintEquation.hpp
    ${MBIRLib_SOURCE_DIR}/Common/DerivOfCostFunc.hpp
    ${MBIRLib_SOURCE_DIR}/Common/TomoArray.hpp
    ${MBIRLib_SOURCE_DIR}/Common/VoxelUpdateList.h
)

cmp_IDE_SOURCE_PROPERTIES( "MBIRLib/Common" "${MBIRLib_Common_HDRS}" "${MBIRLib_Common_SRCS}" "${CMP_INSTALL_FILES}")
