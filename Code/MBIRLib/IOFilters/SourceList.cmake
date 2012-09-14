#--////////////////////////////////////////////////////////////////////////////
#-- Copyright (c) 2011, Michael A. Jackson. BlueQuartz Software
#-- All rights reserved.
#-- BSD License: http://www.opensource.org/licenses/bsd-license.html
#-- This code was partly written under US Air Force Contract FA8650-07-D-5800
#--////////////////////////////////////////////////////////////////////////////
set (MBIRLib_IOFilters_SRCS
    ${MBIRLib_SOURCE_DIR}/IOFilters/AvizoUniformCoordinateWriter.cpp
    ${MBIRLib_SOURCE_DIR}/IOFilters/DetectorResponseWriter.cpp
    ${MBIRLib_SOURCE_DIR}/IOFilters/RawGeometryWriter.cpp
    ${MBIRLib_SOURCE_DIR}/IOFilters/MRCReader.cpp
    ${MBIRLib_SOURCE_DIR}/IOFilters/MRCWriter.cpp
    )

set (MBIRLib_IOFilters_HDRS
    ${MBIRLib_SOURCE_DIR}/IOFilters/AvizoUniformCoordinateWriter.h
    ${MBIRLib_SOURCE_DIR}/IOFilters/DetectorResponseWriter.h
    ${MBIRLib_SOURCE_DIR}/IOFilters/VTKFileWriters.hpp
    ${MBIRLib_SOURCE_DIR}/IOFilters/VTKWriterMacros.h
    ${MBIRLib_SOURCE_DIR}/IOFilters/RawGeometryWriter.h
    ${MBIRLib_SOURCE_DIR}/IOFilters/MRCReader.h
    ${MBIRLib_SOURCE_DIR}/IOFilters/MRCWriter.h
    ${MBIRLib_SOURCE_DIR}/IOFilters/MRCHeader.h

)

cmp_IDE_SOURCE_PROPERTIES( "MBIRLib/IOFilters" "${MBIRLib_IOFilters_HDRS}" "${MBIRLib_IOFilters_SRCS}" "${CMP_INSTALL_FILES}")

