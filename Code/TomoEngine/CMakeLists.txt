

project(TomoEngine)
cmake_minimum_required(VERSION 2.8)
SET (TomoEngine_BUILT_AS_DYNAMIC_LIB)
IF (BUILD_SHARED_LIBS)
    SET (TomoEngine_BUILT_AS_DYNAMIC_LIB 1)
    IF (WIN32)
        ADD_DEFINITIONS("-DTomoEngine_BUILT_AS_DYNAMIC_LIB")
    ENDIF (WIN32)
ENDIF (BUILD_SHARED_LIBS)

SET (PROJECT_CODE_DIR ${PROJECT_SOURCE_DIR}/Source)
if (NOT PROJECT_RESOURCES_DIR)
SET (PROJECT_RESOURCES_DIR ${TomoEngine_SOURCE_DIR}/Resources)
endif()
if (NOT DEFINED CMP_SOURCE_DIR)
set(CMP_SOURCE_DIR ${TomoEngine_SOURCE_DIR}/Support/cmp)
endif()

SET (PROJECT_PREFIX "TomoEngine" CACHE STRING "The Prefix to be used for Preprocessor definitions")
set (TomoEngine_INSTALL_FILES 0)
#-------------------------------------------------------------------------------
# Include the cmp project to perform all the basic configuration tests for the 
# project.
#------------------------------------------------------------------------------- 
set (CMP_HEADER_DIR ${TomoEngine_BINARY_DIR})
set (CMP_CONFIGURATION_FILE_NAME "TomoEngineConfiguration.h")
set (CMP_TYPES_FILE_NAME "TomoEngineTypes.h")
set (CMP_VERSION_HEADER_FILE_NAME "TomoEngineVersion.h")
set (CMP_TOP_HEADER_FILE "TomoEngine.h")
set (CMP_PROJECT_NAMESPACE "TomoEngine")
set (CMP_PROJECT_NAME "TomoEngine")
set (CMP_ENABLE_PLUGINS "0")
set (CMP_LIB_SEARCH_DIRS "")
INCLUDE (${CMP_SOURCE_DIR}/cmpProject.cmake)

# --------------------------------------------------------------------
# Generate a Header file with Compile Version variables
# --------------------------------------------------------------------
set (VERSION_GEN_NAME "TomoEngine")
set (VERSION_GEN_NAMESPACE "TomoEngine")
set (PROJECT_PREFIX "TomoEngine")
set (VERSION_GEN_COMPLETE  "1.0.0")
set (VERSION_GEN_VER_MAJOR "1")
set (VERSION_GEN_VER_MINOR "0")
set (VERSION_GEN_VER_PATCH "0")

configure_file(${CMP_CONFIGURED_FILES_SOURCE_DIR}/cmpVersion.h.in   
               ${TomoEngine_BINARY_DIR}/TomoEngineVersion.h  )

configure_file(${TomoEngine_SOURCE_DIR}/TomoEngineConfiguration.h.in
               ${TomoEngine_BINARY_DIR}/${CMP_TOP_HEADER_FILE})   

include_directories(${TomoEngine_SOURCE_DIR}/..)            
include_directories(${TomoEngine_SOURCE_DIR})
include_directories(${TomoEngine_BINARY_DIR})


# --------------------------------------------------------------------
# Gather all the source files that will be put into the library
set(Project_SRCS "")
Macro(D3DCompileDir DirName)
    include(${TomoEngine_SOURCE_DIR}/${DirName}/SourceList.cmake)
    set (Project_SRCS ${Project_SRCS} ${TomoEngine_${DirName}_SRCS} ${TomoEngine_${DirName}_HDRS})
endMacro()


D3DCompileDir(Common)
D3DCompileDir(IO)
D3DCompileDir(mt)
D3DCompileDir(SOC)



add_library(${PROJECT_NAME}  ${LIB_TYPE} ${Project_SRCS})
target_link_libraries(${PROJECT_NAME} MXA )


LibraryProperties( ${PROJECT_NAME}  ${EXE_DEBUG_EXTENSION} )

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
