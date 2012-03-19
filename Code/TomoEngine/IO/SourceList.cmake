#--////////////////////////////////////////////////////////////////////////////
#-- Copyright (c) 2011, Michael A. Jackson. BlueQuartz Software
#-- All rights reserved.
#-- BSD License: http://www.opensource.org/licenses/bsd-license.html
#-- This code was partly written under US Air Force Contract FA8650-07-D-5800
#--////////////////////////////////////////////////////////////////////////////
set (TomoEngine_IO_SRCS
    ${TomoEngine_SOURCE_DIR}/IO/RawGeometryWriter.cpp
    ${TomoEngine_SOURCE_DIR}/IO/MRCReader.cpp
    ${TomoEngine_SOURCE_DIR}/IO/MRCWriter.cpp
    ${TomoEngine_SOURCE_DIR}/IO/NuisanceParamWriter.cpp
    ${TomoEngine_SOURCE_DIR}/IO/NuisanceParamReader.cpp
    ${TomoEngine_SOURCE_DIR}/IO/SinogramBinWriter.cpp
    )

set (TomoEngine_IO_HDRS
    ${TomoEngine_SOURCE_DIR}/IO/VTKFileWriters.hpp
    ${TomoEngine_SOURCE_DIR}/IO/VTKWriterMacros.h
    ${TomoEngine_SOURCE_DIR}/IO/RawGeometryWriter.h
    ${TomoEngine_SOURCE_DIR}/IO/MRCReader.h
    ${TomoEngine_SOURCE_DIR}/IO/MRCWriter.h
    ${TomoEngine_SOURCE_DIR}/IO/MRCHeader.h
    ${TomoEngine_SOURCE_DIR}/IO/NuisanceParamWriter.h
    ${TomoEngine_SOURCE_DIR}/IO/NuisanceParamReader.h
    ${TomoEngine_SOURCE_DIR}/IO/SinogramBinWriter.h
)
cmp_IDE_SOURCE_PROPERTIES( "IO" "${TomoEngine_IO_HDRS}" "${TomoEngine_IO_SRCS}" "${CMP_INSTALL_FILES}")

