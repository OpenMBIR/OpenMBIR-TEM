#--////////////////////////////////////////////////////////////////////////////
#-- Copyright (c) 2009, Michael A. Jackson. BlueQuartz Software
#-- All rights reserved.
#-- BSD License: http://www.opensource.org/licenses/bsd-license.html
#-- This code was partly written under US Air Force Contract FA8650-07-D-5800
#--////////////////////////////////////////////////////////////////////////////
set (TomoEngine_Common_SRCS
    ${TomoEngine_SOURCE_DIR}/Common/allocate.c
    ${TomoEngine_SOURCE_DIR}/Common/EIMTime.c
    ${TomoEngine_SOURCE_DIR}/Common/AbstractPipeline.cpp
    ${TomoEngine_SOURCE_DIR}/Common/AbstractFilter.cpp
    ${TomoEngine_SOURCE_DIR}/Common/Observer.cpp

    ${TomoEngine_SOURCE_DIR}/Common/Observable.cpp
    ${TomoEngine_SOURCE_DIR}/Common/CalculateAMatrixColumn.cpp   
    ${TomoEngine_SOURCE_DIR}/Common/MRCSinogramInitializer.cpp
    ${TomoEngine_SOURCE_DIR}/Common/RawSinogramInitializer.cpp
    ${TomoEngine_SOURCE_DIR}/Common/GainsOffsetsReader.cpp
    ${TomoEngine_SOURCE_DIR}/Common/InitialReconstructionInitializer.cpp
    ${TomoEngine_SOURCE_DIR}/Common/InitialReconstructionBinReader.cpp
    ${TomoEngine_SOURCE_DIR}/Common/TomoFilter.cpp
    ${TomoEngine_SOURCE_DIR}/Common/ComputeGainsOffsets.cpp
    ${TomoEngine_SOURCE_DIR}/Common/DetectorResponse.cpp
    ${TomoEngine_SOURCE_DIR}/Common/DetectorResponseWriter.cpp
)

set (TomoEngine_Common_HDRS
    ${TomoEngine_SOURCE_DIR}/Common/allocate.h
    ${TomoEngine_SOURCE_DIR}/Common/TomoEngineDLLExport.h
    ${TomoEngine_SOURCE_DIR}/Common/MSVCDefines.h
    ${TomoEngine_SOURCE_DIR}/Common/EIMTime.h
    ${TomoEngine_SOURCE_DIR}/Common/EIMMath.h
    ${TomoEngine_SOURCE_DIR}/Common/AbstractPipeline.h
    ${TomoEngine_SOURCE_DIR}/Common/AbstractFilter.h
    ${TomoEngine_SOURCE_DIR}/Common/Observer.h

    ${TomoEngine_SOURCE_DIR}/Common/Observable.h
    ${TomoEngine_SOURCE_DIR}/Common/CalculateAMatrixColumn.h
    ${TomoEngine_SOURCE_DIR}/Common/MRCSinogramInitializer.h
    ${TomoEngine_SOURCE_DIR}/Common/RawSinogramInitializer.h
    ${TomoEngine_SOURCE_DIR}/Common/GainsOffsetsReader.h
    ${TomoEngine_SOURCE_DIR}/Common/InitialReconstructionInitializer.h
    ${TomoEngine_SOURCE_DIR}/Common/InitialReconstructionBinReader.h
    ${TomoEngine_SOURCE_DIR}/Common/TomoFilter.h
    ${TomoEngine_SOURCE_DIR}/Common/ComputeGainsOffsets.h 
    ${TomoEngine_SOURCE_DIR}/Common/CE_ConstraintEquation.hpp
    ${TomoEngine_SOURCE_DIR}/Common/DerivOfCostFunc.hpp
    ${TomoEngine_SOURCE_DIR}/Common/TomoArray.hpp
    ${TomoEngine_SOURCE_DIR}/Common/DetectorResponse.h
    ${TomoEngine_SOURCE_DIR}/Common/DetectorResponseWriter.h
)

cmp_IDE_SOURCE_PROPERTIES( "Common" "${TomoEngine_Common_HDRS}" "${TomoEngine_Common_SRCS}" "${CMP_INSTALL_FILES}")
