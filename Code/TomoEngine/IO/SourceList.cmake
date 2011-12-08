#--////////////////////////////////////////////////////////////////////////////
#-- Copyright (c) 2011, Michael A. Jackson. BlueQuartz Software
#-- All rights reserved.
#-- BSD License: http://www.opensource.org/licenses/bsd-license.html
#-- This code was partly written under US Air Force Contract FA8650-07-D-5800
#--////////////////////////////////////////////////////////////////////////////
set (TomoEngine_IO_SRCS
    ${TomoEngine_SOURCE_DIR}/IO/RawGeometryWriter.cpp
    ${TomoEngine_SOURCE_DIR}/IO/MRCReader.cpp
)

set (TomoEngine_IO_HDRS
    ${TomoEngine_SOURCE_DIR}/IO/VTKFileWriters.hpp
    ${TomoEngine_SOURCE_DIR}/IO/VTKWriterMacros.h
    ${TomoEngine_SOURCE_DIR}/IO/RawGeometryWriter.h
    ${TomoEngine_SOURCE_DIR}/IO/MRCReader.h
    ${TomoEngine_SOURCE_DIR}/IO/MRCHeader.h
)
cmp_IDE_SOURCE_PROPERTIES( "IO" "${TomoEngine_IO_HDRS}" "${TomoEngine_IO_SRCS}" "${CMP_INSTALL_FILES}")
