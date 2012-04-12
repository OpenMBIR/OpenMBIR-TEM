#--////////////////////////////////////////////////////////////////////////////
#-- Copyright (c) 2009, Michael A. Jackson. BlueQuartz Software
#-- All rights reserved.
#-- BSD License: http://www.opensource.org/licenses/bsd-license.html
#-- This code was partly written under US Air Force Contract FA8650-07-D-5800
#--////////////////////////////////////////////////////////////////////////////
set (TomoEngine_Tiff_SRCS
    ${TomoEngine_SOURCE_DIR}/Tiff/TiffUtilities.cpp
)

set (TomoEngine_Tiff_HDRS
    ${TomoEngine_SOURCE_DIR}/Tiff/TiffUtilities.h
)
cmp_IDE_SOURCE_PROPERTIES( "TomoEngine/Tiff" "${TomoEngine_Tiff_HDRS}" "${TomoEngine_Tiff_SRCS}" "${PROJECT_INSTALL_HEADERS}")

add_executable(tiffTest16Bit ${TomoEngine_SOURCE_DIR}/Tiff/TiffTest16Bit.cpp ${TomoEngine_Tiff_SRCS})
target_link_libraries(tiffTest16Bit ${TIFF_LIBRARIES})




