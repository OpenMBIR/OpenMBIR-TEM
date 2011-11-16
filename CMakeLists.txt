#--////////////////////////////////////////////////////////////////////////////
#--
#--  Copyright (c) 2011, Michael A. Jackson. BlueQuartz Software

#--  All rights reserved.
#--  BSD License: http://www.opensource.org/licenses/bsd-license.html
#--
#-- This code was partly written under US Air Force Contract FA8650-07-D-5800
#--
#--////////////////////////////////////////////////////////////////////////////
project (EIMTomo)
cmake_minimum_required(VERSION 2.8.5)

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
set(CMP_SOURCE_DIR ${EIMTomo_SOURCE_DIR}/Support/cmp)
endif()



SET (PROJECT_CODE_DIR ${PROJECT_SOURCE_DIR}/Code)
SET (PROJECT_RESOURCES_DIR ${PROJECT_SOURCE_DIR}/Resources)
SET (PROJECT_PREFIX "EIMTOMO" CACHE STRING "The 3 Letter Prefix to be used for Preprocessor definitions")
SET (EXE_DEBUG_EXTENSION "_debug")
SET (PROJECT_INSTALL_HEADERS "1")
SET (PROJECT_INSTALL_EXECUTABLES "1")

#-------------------------------------------------------------------------------
# Include the cmp project to perform all the basic configuration tests for the 
# project.
#------------------------------------------------------------------------------- 
set (CMP_HEADER_DIR ${EIMTomo_BINARY_DIR}/EIMTomo)
set (CMP_CONFIGURATION_FILE_NAME "EIMTomoConfiguration.h")
set (CMP_TYPES_FILE_NAME "EIMTomoTypes.h")
set (CMP_VERSION_HEADER_FILE_NAME "EIMTomoVersion.h")
set (CMP_TOP_HEADER_FILE "EIMTomo.h")
set (CMP_PROJECT_NAMESPACE "${PROJECT_NAME}")
set (CMP_PROJECT_NAME "${PROJECT_NAME}")
set (CMP_ENABLE_PLUGINS "0")
set (CMP_LIB_SEARCH_DIRS "")
set (CMP_GENERATE_VERSION_STRING 1)
set (CMP_INSTALL_FILES 1)
INCLUDE (${CMP_SOURCE_DIR}/cmpProject.cmake)

# We need to add the binary directory to our list of include directories so the source files
# can find the generated header files
include_directories(${EIMTomo_BINARY_DIR})

# --------------------------------------------------------------------
# Generate our ReadMe and License Files
configure_file(${PROJECT_RESOURCES_DIR}/EIMTomo/License.txt.in
                ${PROJECT_BINARY_DIR}/License.txt )
configure_file(${PROJECT_RESOURCES_DIR}/EIMTomo/ReadMe.txt.in
                ${PROJECT_BINARY_DIR}/ReadMe.txt )
                
# --------------------------------------------------------------------
# Build shared libraries
OPTION (BUILD_SHARED_LIBS "Build Shared Libraries" OFF)
set (EIMTomo_BUILT_AS_DYNAMIC_LIB )
SET (LIB_TYPE STATIC)
IF (BUILD_SHARED_LIBS)
  SET (LIB_TYPE SHARED)
  SET (EIMTomo_BUILT_AS_DYNAMIC_LIB 1)
    IF (WIN32)
    ADD_DEFINITIONS("-DEIMTomo_BUILT_AS_DYNAMIC_LIB")
  ENDIF (WIN32)
ENDIF (BUILD_SHARED_LIBS)

OPTION (EIMTomo_TIFF_SUPPORT "Add libTiff Support" ON)
OPTION (EIMTomo_HDF5_SUPPORT "Add HDF5 Support" OFF)
OPTION (EIMTomo_USE_MERSENNE_TWISTER "Use the Mersenne Twister RNG" ON)

# --------------------------------------------------------------------
# Should we use Intel Threading Building Blocks
# --------------------------------------------------------------------
option(EIMTomo_USE_PARALLEL_ALGORITHMS "Use MultiThreaded Algorithms where possible" OFF)
if (EIMTomo_USE_PARALLEL_ALGORITHMS)
    find_package(TBB)
    if (TBB_FOUND)
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
# Now configure this extra file with some additional definitions
configure_file(${EIMTomo_SOURCE_DIR}/Resources/EIMTomo/EIMTomoConfiguration.h.in
               ${EIMTomo_BINARY_DIR}/EIMTomo/${CMP_TOP_HEADER_FILE})
               
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

# Add the MXA as a library project
add_subdirectory( ${PROJECT_CODE_DIR}/MXA ${PROJECT_BINARY_DIR}/MXA)

# --------------------------------------------------------------------
# Macro that takes a subdirectory and creates a cmake option variable
# which allows developers to turn off those sub projects if needed.
# --------------------------------------------------------------------
macro(Add_Algorithm name)
    option( EIMTomo_BUILD_${name} "" ON )
    if ( ${EIMTomo_BUILD_${name}} )
        add_subdirectory( ${PROJECT_CODE_DIR}/EIMTomo/${name})
    endif()
endmacro()
# Add our executables

Add_Algorithm(JointEstimationAlgorithm)
Add_Algorithm(NonHomogeniousICD)
Add_Algorithm(BasicReconstructionAlgorithm)
Add_Algorithm(ScaleOffsetCorrectionAlgorithm)



