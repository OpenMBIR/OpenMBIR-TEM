#--////////////////////////////////////////////////////////////////////////////
#--
#--  Copyright (c) 2011, Michael A. Jackson. BlueQuartz Software

#--  All rights reserved.
#--  BSD License: http://www.opensource.org/licenses/bsd-license.html
#--
#-- This code was partly written under US Air Force Contract FA8650-07-D-5800
#--
#--////////////////////////////////////////////////////////////////////////////

project(HaadfMbirLib)
cmake_minimum_required(VERSION 2.8.6)
SET (HaadfMbirLib_BUILT_AS_DYNAMIC_LIB)
IF (BUILD_SHARED_LIBS)
    SET (HaadfMbirLib_BUILT_AS_DYNAMIC_LIB 1)
    IF (WIN32)
        ADD_DEFINITIONS("-DHaadfMbirLib_BUILT_AS_DYNAMIC_LIB")
    ENDIF (WIN32)
ENDIF (BUILD_SHARED_LIBS)

SET (PROJECT_CODE_DIR ${PROJECT_SOURCE_DIR}/Source)
if (NOT PROJECT_RESOURCES_DIR)
SET (PROJECT_RESOURCES_DIR ${HaadfMbirLib_SOURCE_DIR}/Resources)
endif()
if (NOT DEFINED CMP_SOURCE_DIR)
set(CMP_SOURCE_DIR ${HaadfMbirLib_SOURCE_DIR}/Support/cmp)
endif()

SET (PROJECT_PREFIX "HaadfMbirLib" CACHE STRING "The Prefix to be used for Preprocessor definitions")
set (HaadfMbirLib_INSTALL_FILES 0)


#-------------------------------------------------------------------------------
# Include the cmp project to perform all the basic configuration tests for the 
# project.
#------------------------------------------------------------------------------- 
set (CMP_HEADER_DIR ${HaadfMbirLib_BINARY_DIR})
set (CMP_CONFIGURATION_FILE_NAME "HaadfMbirLibConfiguration.h")
set (CMP_TYPES_FILE_NAME "HaadfMbirLibTypes.h")
set (CMP_VERSION_HEADER_FILE_NAME "HaadfMbirLibVersion.h")
set (CMP_TOP_HEADER_FILE "HaadfMbirLib.h")
set (CMP_PROJECT_NAMESPACE "HaadfMbirLib")
set (CMP_PROJECT_NAME "HaadfMbirLib")
set (CMP_ENABLE_PLUGINS "0")
set (CMP_LIB_SEARCH_DIRS "")
set (CMP_GENERATE_VERSION_STRING 1)
set (CMP_INSTALL_FILES 0)
INCLUDE (${CMP_SOURCE_DIR}/cmpProject.cmake)


# --------------------------------------------------------------------
# Generate a Header file with Compile Version variables
# --------------------------------------------------------------------
configure_file(${HaadfMbirLib_SOURCE_DIR}/HaadfMbirLibConfiguration.h.in
               ${HaadfMbirLib_BINARY_DIR}/${CMP_TOP_HEADER_FILE})   

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



set (HaadfMbirLib_SRCS
    ${HaadfMbirLib_SOURCE_DIR}/ForwardProject.cpp
    ${HaadfMbirLib_SOURCE_DIR}/HaadfMbirEngine.cpp
    ${HaadfMbirLib_SOURCE_DIR}/HaadfMbirEngine_UpdateVoxels.cpp
    ${HaadfMbirLib_SOURCE_DIR}/HaadfMbirEngine_Extra.cpp
    ${HaadfMbirLib_SOURCE_DIR}/HaadfMbirInputs.cpp
    ${HaadfMbirLib_SOURCE_DIR}/MultiResolutionHaadfMbir.cpp
    ${HaadfMbirLib_SOURCE_DIR}/QGGMRF_Functions.cpp
)

set (HaadfMbirLib_HDRS
    ${HaadfMbirLib_SOURCE_DIR}/ForwardProject.h
    ${HaadfMbirLib_SOURCE_DIR}/HaadfMbirConstants.h
    ${HaadfMbirLib_SOURCE_DIR}/HaadfMbirEngine.h
    ${HaadfMbirLib_SOURCE_DIR}/HaadfMbirInputs.h
    ${HaadfMbirLib_SOURCE_DIR}/HaadfMbirStructures.h
    ${HaadfMbirLib_SOURCE_DIR}/MultiResolutionHaadfMbir.h
    ${HaadfMbirLib_SOURCE_DIR}/QGGMRF_Functions.h
)
set_source_files_properties( ${HaadfMbirLib_SOURCE_DIR}/HaadfMbirEngine_UpdateVoxels.cpp
                             ${HaadfMbirLib_SOURCE_DIR}/HaadfMbirEngine_Extra.cpp
                             PROPERTIES HEADER_FILE_ONLY TRUE)
cmp_IDE_SOURCE_PROPERTIES( "HaadfMbirLib" "${HaadfMbirLib_HDRS}" "${HaadfMbirLib_SRCS}" "${CMP_INSTALL_FILES}")


# --------------------------------------------------------------------
# Add the ReconstructionCoreLib Library
add_library(HaadfMbirLib  ${LIB_TYPE} ${HaadfMbirLib_SRCS} ${HaadfMbirLib_HDRS} )
if (OpenMBIR_USE_MultiThreading)
    target_link_libraries(HaadfMbirLib MXA ${TBB_LIBRARIES})
else()
    target_link_libraries(HaadfMbirLib MXA)
endif()
LibraryProperties( HaadfMbirLib  ${EXE_DEBUG_EXTENSION} )

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

