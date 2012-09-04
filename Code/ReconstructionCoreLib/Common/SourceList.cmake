#--////////////////////////////////////////////////////////////////////////////
#-- Copyright (c) 2009, Michael A. Jackson. BlueQuartz Software
#-- All rights reserved.
#-- BSD License: http://www.opensource.org/licenses/bsd-license.html
#-- This code was partly written under US Air Force Contract FA8650-07-D-5800
#--////////////////////////////////////////////////////////////////////////////
set (ReconstructionCoreLib_Common_SRCS
    ${ReconstructionCoreLib_SOURCE_DIR}/Common/allocate.c
    ${ReconstructionCoreLib_SOURCE_DIR}/Common/EIMTime.c
    ${ReconstructionCoreLib_SOURCE_DIR}/Common/EIMImage.cpp
    ${ReconstructionCoreLib_SOURCE_DIR}/Common/AbstractFilter.cpp
    ${ReconstructionCoreLib_SOURCE_DIR}/Common/FilterPipeline.cpp
    ${ReconstructionCoreLib_SOURCE_DIR}/Common/Observer.cpp
    ${ReconstructionCoreLib_SOURCE_DIR}/Common/Observable.cpp
)

set (ReconstructionCoreLib_Common_HDRS
    ${ReconstructionCoreLib_SOURCE_DIR}/Common/allocate.h
    ${ReconstructionCoreLib_SOURCE_DIR}/Common/ReconstructionCoreLibDLLExport.h
    ${ReconstructionCoreLib_SOURCE_DIR}/Common/MSVCDefines.h
    ${ReconstructionCoreLib_SOURCE_DIR}/Common/EIMImage.h
    ${ReconstructionCoreLib_SOURCE_DIR}/Common/EIMTime.h
    ${ReconstructionCoreLib_SOURCE_DIR}/Common/EIMMath.h
    ${ReconstructionCoreLib_SOURCE_DIR}/Common/AbstractFilter.h
    ${ReconstructionCoreLib_SOURCE_DIR}/Common/FilterPipeline.h
    ${ReconstructionCoreLib_SOURCE_DIR}/Common/Observer.h
    ${ReconstructionCoreLib_SOURCE_DIR}/Common/Observable.h
    ${ReconstructionCoreLib_SOURCE_DIR}/Common/CE_ConstraintEquation.hpp
    ${ReconstructionCoreLib_SOURCE_DIR}/Common/DerivOfCostFunc.hpp
    ${ReconstructionCoreLib_SOURCE_DIR}/Common/TomoArray.hpp
)

cmp_IDE_SOURCE_PROPERTIES( "ReconstructionCoreLib/Common" "${ReconstructionCoreLib_Common_HDRS}" "${ReconstructionCoreLib_Common_SRCS}" "${CMP_INSTALL_FILES}")
