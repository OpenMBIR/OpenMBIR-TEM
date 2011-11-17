#--////////////////////////////////////////////////////////////////////////////
#-- Copyright (c) 2011, Michael A. Jackson. BlueQuartz Software
#-- All rights reserved.
#-- BSD License: http://www.opensource.org/licenses/bsd-license.html
#-- This code was partly written under US Air Force Contract FA8650-07-D-5800
#--////////////////////////////////////////////////////////////////////////////
set (EIMTomo_IO_SRCS
    ${PROJECT_CODE_DIR}/EIMTomo/IO/RawGeometryWriter.cpp
)

set (EIMTomo_IO_HDRS
    ${PROJECT_CODE_DIR}/EIMTomo/IO/VTKFileWriters.hpp
    ${PROJECT_CODE_DIR}/EIMTomo/IO/VTKWriterMacros.h
    ${PROJECT_CODE_DIR}/EIMTomo/IO/RawGeometryWriter.h
)
cmp_IDE_SOURCE_PROPERTIES( "EIMTomo/IO" "${EIMTomo_IO_HDRS}" "${EIMTomo_IO_SRCS}" "${CMP_INSTALL_FILES}")
