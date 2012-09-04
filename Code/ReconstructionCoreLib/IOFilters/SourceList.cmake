#--////////////////////////////////////////////////////////////////////////////
#-- Copyright (c) 2011, Michael A. Jackson. BlueQuartz Software
#-- All rights reserved.
#-- BSD License: http://www.opensource.org/licenses/bsd-license.html
#-- This code was partly written under US Air Force Contract FA8650-07-D-5800
#--////////////////////////////////////////////////////////////////////////////
set (ReconstructionCoreLib_IOFilters_SRCS
    ${ReconstructionCoreLib_SOURCE_DIR}/IOFilters/AvizoUniformCoordinateWriter.cpp
    ${ReconstructionCoreLib_SOURCE_DIR}/IOFilters/DetectorResponseWriter.cpp
    ${ReconstructionCoreLib_SOURCE_DIR}/IOFilters/GainsOffsetsReader.cpp
    ${ReconstructionCoreLib_SOURCE_DIR}/IOFilters/RawGeometryWriter.cpp
    ${ReconstructionCoreLib_SOURCE_DIR}/IOFilters/MRCReader.cpp
    ${ReconstructionCoreLib_SOURCE_DIR}/IOFilters/MRCWriter.cpp
    ${ReconstructionCoreLib_SOURCE_DIR}/IOFilters/NuisanceParamWriter.cpp
    ${ReconstructionCoreLib_SOURCE_DIR}/IOFilters/NuisanceParamReader.cpp
    ${ReconstructionCoreLib_SOURCE_DIR}/IOFilters/SinogramBinWriter.cpp
    )

set (ReconstructionCoreLib_IOFilters_HDRS
    ${ReconstructionCoreLib_SOURCE_DIR}/IOFilters/AvizoUniformCoordinateWriter.h
    ${ReconstructionCoreLib_SOURCE_DIR}/IOFilters/DetectorResponseWriter.h
    ${ReconstructionCoreLib_SOURCE_DIR}/IOFilters/GainsOffsetsReader.h
    ${ReconstructionCoreLib_SOURCE_DIR}/IOFilters/VTKFileWriters.hpp
    ${ReconstructionCoreLib_SOURCE_DIR}/IOFilters/VTKWriterMacros.h
    ${ReconstructionCoreLib_SOURCE_DIR}/IOFilters/RawGeometryWriter.h
    ${ReconstructionCoreLib_SOURCE_DIR}/IOFilters/MRCReader.h
    ${ReconstructionCoreLib_SOURCE_DIR}/IOFilters/MRCWriter.h
    ${ReconstructionCoreLib_SOURCE_DIR}/IOFilters/MRCHeader.h
    ${ReconstructionCoreLib_SOURCE_DIR}/IOFilters/NuisanceParamWriter.h
    ${ReconstructionCoreLib_SOURCE_DIR}/IOFilters/NuisanceParamReader.h
    ${ReconstructionCoreLib_SOURCE_DIR}/IOFilters/SinogramBinWriter.h
)
cmp_IDE_SOURCE_PROPERTIES( "ReconstructionCoreLib/IOFilters" "${ReconstructionCoreLib_IOFilters_HDRS}" "${ReconstructionCoreLib_IOFilters_SRCS}" "${CMP_INSTALL_FILES}")

