#--////////////////////////////////////////////////////////////////////////////
#-- Copyright (c) 2011, Michael A. Jackson. BlueQuartz Software
#-- All rights reserved.
#-- BSD License: http://www.opensource.org/licenses/bsd-license.html
#-- This code was partly written under US Air Force Contract FA8650-07-D-5800
#--////////////////////////////////////////////////////////////////////////////

project(MRCSubset)
cmake_minimum_required(VERSION 2.8.5)


# --------------------------------------------------------------------
# Setup the install rules for the various platforms
set(install_dir "tools")
if (WIN32)
    set (install_dir ".")
endif()

set(MRCSubset_SRCS
    ${MRCSubset_SOURCE_DIR}/MRCSubset.cpp
    )

BuildToolBundle(TARGET MRCSubset
            SOURCES ${MRCSubset_SRCS}
            DEBUG_EXTENSION ${EXE_DEBUG_EXTENSION}
            VERSION_MAJOR
            VERSION_MAJOR ${OpenMBIR_VER_MAJOR}
            VERSION_MINOR ${OpenMBIR_VER_MINOR}
            VERSION_PATCH ${OpenMBIR_VER_PATCH}
            BINARY_DIR    ${PROJECT_BINARY_DIR}/MRCSubset
            COMPONENT     Applications
            INSTALL_DEST  ${install_dir}
            LINK_LIBRARIES MXA MBIRLib 
            LIB_SEARCH_DIRS ${CMAKE_LIBRARY_OUTPUT_DIRECTORY} ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}
)
