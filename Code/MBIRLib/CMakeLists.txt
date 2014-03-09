#--////////////////////////////////////////////////////////////////////////////
#--
#--  Copyright (c) 2011, Michael A. Jackson. BlueQuartz Software

#--  All rights reserved.
#--  BSD License: http://www.opensource.org/licenses/bsd-license.html
#--
#-- This code was partly written under US Air Force Contract FA8650-07-D-5800
#--
#--////////////////////////////////////////////////////////////////////////////

project(MBIRLib)
cmake_minimum_required(VERSION 2.8.6)
SET (MBIRLib_BUILT_AS_DYNAMIC_LIB)
IF (BUILD_SHARED_LIBS)
    SET (MBIRLib_BUILT_AS_DYNAMIC_LIB 1)
    IF (WIN32)
        ADD_DEFINITIONS("-DMBIRLib_BUILT_AS_DYNAMIC_LIB")
    ENDIF (WIN32)
ENDIF (BUILD_SHARED_LIBS)

SET (PROJECT_CODE_DIR ${PROJECT_SOURCE_DIR}/Source)
if (NOT PROJECT_RESOURCES_DIR)
SET (PROJECT_RESOURCES_DIR ${MBIRLib_SOURCE_DIR}/Resources)
endif()
if (NOT DEFINED CMP_SOURCE_DIR)
set(CMP_SOURCE_DIR ${MBIRLib_SOURCE_DIR}/Support/cmp)
endif()

SET (PROJECT_PREFIX "MBIRLib" CACHE STRING "The Prefix to be used for Preprocessor definitions")
set (MBIRLib_INSTALL_FILES 0)


#-------------------------------------------------------------------------------
# Include the cmp project to perform all the basic configuration tests for the 
# project.
#------------------------------------------------------------------------------- 
set (CMP_HEADER_DIR ${MBIRLib_BINARY_DIR})
set (CMP_CONFIGURATION_FILE_NAME "MBIRLibConfiguration.h")
set (CMP_TYPES_FILE_NAME "MBIRLibTypes.h")
set (CMP_VERSION_HEADER_FILE_NAME "MBIRLibVersion.h")
set (CMP_VERSION_SOURCE_FILE_NAME "MBIRLibVersion.cpp")
set (CMP_TOP_HEADER_FILE "MBIRLib.h")
set (CMP_PROJECT_NAMESPACE "MBIRLib")
set (CMP_PROJECT_NAME "MBIRLib")
set (CMP_ENABLE_PLUGINS "0")
set (CMP_LIB_SEARCH_DIRS "")
set (CMP_GENERATE_VERSION_STRING 1)
set (CMP_INSTALL_FILES 0)
INCLUDE (${CMP_SOURCE_DIR}/cmpProject.cmake)


# --------------------------------------------------------------------
# Generate a Header file with Compile Version variables
# --------------------------------------------------------------------
cmpVersionStringsFromGit(GENERATED_HEADER_FILE_PATH ${CMP_HEADER_DIR}/MBIRLibVersion.h
                         GENERATED_SOURCE_FILE_PATH ${CMP_HEADER_DIR}/MBIRLibVersion.cpp
                         NAMESPACE "${CMP_PROJECT_NAMESPACE}"
                         cmpProjectName "${CMP_PROJECT_NAME}" )

# --------------------------------------------------------------------
# Generate a Header file with Compile Version variables
# --------------------------------------------------------------------
configure_file(${MBIRLib_SOURCE_DIR}/MBIRLibConfiguration.h.in
               ${MBIRLib_BINARY_DIR}/${CMP_TOP_HEADER_FILE})   

cmp_IDE_GENERATED_PROPERTIES("Generated" 
                             "${CMP_HEADER_DIR}/${CMP_TOP_HEADER_FILE};${CMP_HEADER_DIR}/${CMP_TYPES_FILE_NAME};${CMP_HEADER_DIR}/${CMP_CONFIGURATION_FILE_NAME};${CMP_HEADER_DIR}/${CMP_VERSION_HEADER_FILE_NAME}"
                             "${CMP_HEADER_DIR}/${CMP_VERSION_SOURCE_FILE_NAME}")
set(Project_SRCS 
    ${CMP_HEADER_DIR}/${CMP_TOP_HEADER_FILE}
    ${CMP_HEADER_DIR}/${CMP_TYPES_FILE_NAME}
    ${CMP_HEADER_DIR}/${CMP_CONFIGURATION_FILE_NAME}
    ${CMP_HEADER_DIR}/${CMP_VERSION_HEADER_FILE_NAME}
    ${CMP_HEADER_DIR}/${CMP_VERSION_SOURCE_FILE_NAME}
)
                             
# --------------------------------------------------------------------
# Gather all the source files that will be put into the library

Macro(CMP_CompileDir DirName)
    include(${MBIRLib_SOURCE_DIR}/${DirName}/SourceList.cmake)
    set (Project_SRCS ${Project_SRCS} ${MBIRLib_${DirName}_SRCS} ${MBIRLib_${DirName}_HDRS})
endMacro()


CMP_CompileDir(Common)
CMP_CompileDir(GenericFilters)
CMP_CompileDir(IOFilters)
CMP_CompileDir(HAADF)
CMP_CompileDir(Reconstruction)

include_directories(${MBIRLib_SOURCE_DIR}/..)            
include_directories(${MBIRLib_SOURCE_DIR})
include_directories(${MBIRLib_BINARY_DIR})

# --------------------------------------------------------------------
# Add the MBIRLib Library
add_library(MBIRLib  ${LIB_TYPE} ${Project_SRCS})
if (OpenMBIR_USE_PARALLEL_ALGORITHMS)
    target_link_libraries(MBIRLib MXA ${TBB_LIBRARIES})
else()
    target_link_libraries(MBIRLib MXA)
endif()
LibraryProperties( MBIRLib  ${EXE_DEBUG_EXTENSION} )

set(install_dir "tools")
set(lib_install_dir "lib")
if (WIN32)
    set (install_dir ".")
    set(lib_install_dir "")
endif()


if (BUILD_SHARED_LIBS AND NOT APPLE)
    INSTALL(TARGETS ${PROJECT_NAME}  
        COMPONENT Applications
        RUNTIME DESTINATION ${install_dir}
        LIBRARY DESTINATION ${lib_install_dir}
        ARCHIVE DESTINATION ${lib_install_dir}       
        BUNDLE DESTINATION "."
    ) 
endif() 
