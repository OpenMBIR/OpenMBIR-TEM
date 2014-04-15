#--////////////////////////////////////////////////////////////////////////////
#--
#--  Copyright (c) 2011, Michael A. Jackson. BlueQuartz Software

#--  All rights reserved.
#--  BSD License: http://www.opensource.org/licenses/bsd-license.html
#--
#-- This code was partly written under US Air Force Contract FA8650-07-D-5800
#--
#--////////////////////////////////////////////////////////////////////////////
project (OpenMBIR)
cmake_minimum_required(VERSION 2.8.9)

# We setup all the libraries and exectuables to all be placed into a single directory
# ---------- Setup Library output Directories -------------------------
SET (CMAKE_LIBRARY_OUTPUT_DIRECTORY
  ${PROJECT_BINARY_DIR}/Bin
  CACHE PATH
  "Single Directory for all Libraries"
  )

# --------- Setup the Executable output Directory -------------
SET (CMAKE_RUNTIME_OUTPUT_DIRECTORY
  ${PROJECT_BINARY_DIR}/Bin
  CACHE PATH
  "Single Directory for all Executables."
  )

# --------- Setup the Executable output Directory -------------
SET (CMAKE_ARCHIVE_OUTPUT_DIRECTORY
  ${PROJECT_BINARY_DIR}/Bin
  CACHE PATH
  "Single Directory for all static libraries."
  )


if (NOT DEFINED CMP_SOURCE_DIR)
set(CMP_SOURCE_DIR ${OpenMBIR_SOURCE_DIR}/Support/cmp)
endif()

if(UNIX AND NOT APPLE)
    set(OpenMBIR_LINUX_SYSTEM ${CMAKE_SYSTEM_NAME} CACHE STRING "The Name of the linux distrubution for packaging")
endif()


SET (PROJECT_CODE_DIR ${PROJECT_SOURCE_DIR}/Code)
SET (PROJECT_RESOURCES_DIR ${PROJECT_SOURCE_DIR}/Resources)
SET (PROJECT_PREFIX "OpenMBIR" CACHE STRING "The Prefix to be used for Preprocessor definitions")
SET (EXE_DEBUG_EXTENSION "_debug")
SET (PROJECT_INSTALL_HEADERS "0")
SET (PROJECT_INSTALL_EXECUTABLES "0")

#-------------------------------------------------------------------------------
# Include the cmp project to perform all the basic configuration tests for the 
# project.
#------------------------------------------------------------------------------- 
set (CMP_HEADER_DIR ${OpenMBIR_BINARY_DIR}/DONOTUSE)
set (CMP_INSTALL_FILES 1)
set (CMP_CONFIGURED_FILES_SOURCE_DIR ${CMP_SOURCE_DIR}/ConfiguredFiles)
set (CMP_MODULES_SOURCE_DIR ${CMP_SOURCE_DIR}/Modules CACHE INTERNAL "")
include( ${CMP_SOURCE_DIR}/cmpCMakeMacros.cmake)


if(UNIX AND NOT APPLE)
    set(OpenMBIR_LINUX_SYSTEM ${CMAKE_SYSTEM_NAME} CACHE STRING "The Name of the linux distrubution for packaging")
endif()


# --------------------------------------------------------------------
# In order to increment the version of DREAM3D simply adjust these variables.
# we are using a 3 part Version scheme commonly found in software development. The
# example would be Version 2.1.6543 which translates to Version [major].[minor].[patch].
# The Major version is set here in the CMake file. We will use Git
# to figure out the "patch" leve of the version by using "git describe" which will
# calculate how many commits since the last "tag" was applied to the repository.
# this basically gives us an always incrementing number for the patch level
# and we can set the Major and Minor versions as we see fit.
set(CMP_VERSION_MAJOR "v1")

# --------------------------------------------------------------------
# Over ride CMake's built in module directory by prepending cmp's module
# directory first
SET (CMAKE_MODULE_PATH ${CMP_MODULES_SOURCE_DIR} ${CMAKE_MODULE_PATH})


# We need to add the binary directory to our list of include directories so the source files
# can find the generated header files
include_directories(${OpenMBIR_BINARY_DIR})
include_directories(${PROJECT_CODE_DIR})

# --------------------------------------------------------------------
# Generate our ReadMe and License Files
configure_file(${PROJECT_RESOURCES_DIR}/OpenMBIR/License.txt.in
                ${PROJECT_BINARY_DIR}/License.txt )
configure_file(${PROJECT_RESOURCES_DIR}/OpenMBIR/ReadMe.txt.in
                ${PROJECT_BINARY_DIR}/ReadMe.txt )
                
# --------------------------------------------------------------------
# Build shared libraries
OPTION (BUILD_SHARED_LIBS "Build Shared Libraries" OFF)
set (OpenMBIR_BUILT_AS_DYNAMIC_LIB )
SET (LIB_TYPE STATIC)
IF (BUILD_SHARED_LIBS)
  SET (LIB_TYPE SHARED)
  SET (OpenMBIR_BUILT_AS_DYNAMIC_LIB 1)
    IF (WIN32)
    ADD_DEFINITIONS("-DOpenMBIR_BUILT_AS_DYNAMIC_LIB")
  ENDIF (WIN32)
ENDIF (BUILD_SHARED_LIBS)

set (OpenMBIR_HDF5_SUPPORT 0)

# ---------- Find Boost Headers/Libraries -----------------------
#SET (Boost_FIND_REQUIRED FALSE)
SET (Boost_DEBUG FALSE)
set (Boost_USE_MULTITHREADED TRUE)
set (Boost_USE_STATIC_LIBS TRUE)
SET (Boost_ADDITIONAL_VERSIONS "1.55.0" "1.47.0" "1.46.0" "1.44.0" "1.44" "1.41" "1.41.0" )
set (Boost_FIND_COMPONENTS "")
set (Boost_NO_BOOST_CMAKE 1)
FIND_PACKAGE(Boost)
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})
LINK_DIRECTORIES(${Boost_LIBRARY_DIRS})

# --------------------------------------------------------------------
# Should we use Intel Threading Building Blocks
# --------------------------------------------------------------------
set(OpenMBIR_USE_PARALLEL_ALGORITHMS 0)
option(OpenMBIR_USE_MultiThreading "Use MultiThreaded Algorithms where possible" ON)
if (OpenMBIR_USE_MultiThreading)
    find_package(TBB)
    if (TBB_FOUND)
        set(OpenMBIR_USE_PARALLEL_ALGORITHMS 1)
        include_directories(${TBB_INCLUDE_DIRS} )
        CMP_COPY_DEPENDENT_LIBRARIES("tbb;tbb_malloc")
        # Append the locations of the TBB libraries to our Library Search Paths
        list(APPEND CMP_LIB_SEARCH_DIRS ${TBB_BINARY_DIR} ${TBBG_LIBRARY_DIR} )
        # Create an install rule to copy the dependent shared library over to the installed location
        CMP_LIBRARIES_INSTALL_RULES("tbb;tbb_malloc" ".")
    else()
        message(FATAL_ERROR "The Intel Threading Building Blocks library is needed to enable the multithreaded algorithms. Please make sure it is installed properly")
    endif()
endif()

# --------------------------------------------------------------------
# Look for Lib Tiff
# --------------------------------------------------------------------
set (OpenMBIR_TIFF_SUPPORT "0")
OPTION (OpenMBIR_USE_TIFF "Add libTiff Support" OFF)
if (OpenMBIR_USE_TIFF)
    set(CMAKE_MODULE_PATH ${CMP_MODULES_SOURCE_DIR} ${CMAKE_MODULE_PATH} )
    FIND_PACKAGE(Tiff)
    IF (TIFF_FOUND)
        set (OpenMBIR_TIFF_SUPPORT "1")
        INCLUDE_DIRECTORIES(${TIFF_INCLUDE_DIRS})
        # Copy DLL libraries into the correct build location on MSVC platforms
        if (HAVE_TIFF_DLL)
            CMP_COPY_DEPENDENT_LIBRARIES ("tiff")
            CMP_LIBRARIES_INSTALL_RULES ("tiff" "." )
        endif()  
    else()
       message(STATUS "Tiff was asked to be used but was not located on your system. Tiff Support will be disabled")
       set(OpenMBIR_USE_TIFF "OFF")
       set (OpenMBIR_TIFF_SUPPORT "0")
    ENDIF ()
endif()

# --------------------------------------------------------------------
# Add in some compiler definitions
# --------------------------------------------------------------------
IF ( CMAKE_BUILD_TYPE MATCHES Debug )
  ADD_DEFINITIONS(-DDEBUG)
  IF (NOT WIN32)
     #Doing this with NMake files will overwhelm the output
    ADD_DEFINITIONS(-Wall)
  ENDIF()
ENDIF ( CMAKE_BUILD_TYPE MATCHES Debug )

# --------------------------------------------------------------------
# If we are using GCC, make the compiler messages on a single line
IF(CMAKE_COMPILER_IS_GNUCC)
  SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fmessage-length=0")
ENDIF(CMAKE_COMPILER_IS_GNUCC)
IF(CMAKE_COMPILER_IS_GNUCXX)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fmessage-length=0")
ENDIF(CMAKE_COMPILER_IS_GNUCXX)

if (MSVC)
  add_definitions(-D_CRT_SECURE_NO_WARNINGS)
endif()

# --------------------------------------------------------------------
# Generate the License File Headers
include(${PROJECT_SOURCE_DIR}/License/License.cmake)

# Add the MXA as a library project
add_subdirectory( ${PROJECT_CODE_DIR}/MXA ${PROJECT_BINARY_DIR}/MXA)


add_subdirectory( ${PROJECT_CODE_DIR}/MBIRLib ${PROJECT_BINARY_DIR}/MBIRLib)

add_subdirectory( ${PROJECT_CODE_DIR}/Applications/BFReconstruction ${PROJECT_BINARY_DIR}/Applications/BFReconstruction)

add_subdirectory( ${PROJECT_CODE_DIR}/Applications/MRCSubset ${PROJECT_BINARY_DIR}/Applications/MRCSubset)

# --------------------------------------------------------------------
# Macro that takes a subdirectory and creates a cmake option variable
# which allows developers to turn off those sub projects if needed.
# --------------------------------------------------------------------
macro(Add_Algorithm name)
    option( OpenMBIR_BUILD_${name} "" ON )
    if ( ${OpenMBIR_BUILD_${name}} )
        add_subdirectory( ${PROJECT_CODE_DIR}/OpenMBIR/${name})
    endif()
endmacro()


add_subdirectory(${PROJECT_CODE_DIR}/Test ${PROJECT_BINARY_DIR}/Test)

# --------------------------------------------------------------------
# Build the Qt Based GUI
option(OpenMBIR_BUILD_GUI "Build the Qt Based Gui Application" ON)
if (OpenMBIR_BUILD_GUI)
# ------------------------------------------------------------------------------
#  Qt 4 Section
#  if we are going to build Qt based applications then we need to find Qt 4 first
# ------------------------------------------------------------------------------
# by default only QtCore and QtGui modules are enabled
# other modules must be enabled like this:  
    IF (WIN32)
        SET (QT_USE_QTMAIN TRUE)
    ENDIF (WIN32)

   # this command finds Qt4 libraries and sets all required variables
   # note that it's Qt4, not QT4 or qt4
    FIND_PACKAGE( Qt4 REQUIRED )
  #  message(STATUS "QT_INCLUDE_DIR: ${QT_INCLUDE_DIR}")
  #  message(STATUS "QT_USE_FRAMEWORKS: ${QT_USE_FRAMEWORKS}")
    IF (QT4_FOUND)
        CMP_COPY_QT4_RUNTIME_LIBRARIES( "QtCore;QtGui")
        CMP_QT_LIBRARIES_INSTALL_RULES("QtCore;QtGui" ".")
    endif()
 #   message(STATUS "QT_INCLUDE_DIR: ${QT_INCLUDE_DIR}")
    # add some useful macros and variables
    # (QT_USE_FILE is a variable defined by FIND_PACKAGE( Qt4 ) that 
    # contains a path to CMake script)
    INCLUDE( ${QT_USE_FILE} )
 #   message(STATUS "QT_INCLUDE_DIR: ${QT_INCLUDE_DIR}")
 #   message(STATUS "QT_USE_FRAMEWORKS: ${QT_USE_FRAMEWORKS}")
    # Append the locations of the Qt libraries to our Library Search Paths
    list(APPEND CMP_LIB_SEARCH_DIRS ${QT_BINARY_DIR} ${QT_LIBRARY_DIR} )
   
    # Build the static QtSupport library - which MUST be a static library or the linking
    # will not work due to signals/slots mechanism.
    add_subdirectory(${OpenMBIR_SOURCE_DIR}/Code/QtSupport ${PROJECT_BINARY_DIR}/QtSupport)
   
    add_subdirectory(${OpenMBIR_SOURCE_DIR}/Code/Applications/BrightFieldGui ${PROJECT_BINARY_DIR}/Applications/BrightFieldGui)

endif()


# Should we build the documentation using DOxygen
include (${PROJECT_SOURCE_DIR}/Documentation/CMakeLists.txt)


# This should be the last line in this file:
include(${PROJECT_RESOURCES_DIR}/CPack/PackageProject.cmake)



