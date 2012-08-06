
option (EIMTomo_USE_TIFF "Use libTiff for some image output" ON)
if (EIMTomo_USE_TIFF)
 # Prepend our own "FindTiff" code in front of what CMake itself provides
    set(CMAKE_MODULE_PATH ${CMP_MODULES_SOURCE_DIR} ${CMAKE_MODULE_PATH} )
    FIND_PACKAGE(Tiff)
    IF (NOT TIFF_FOUND)
          MESSAGE (FATAL_ERROR " Tiff Library was NOT Found and is needed.")
    ENDIF (NOT TIFF_FOUND)
    INCLUDE_DIRECTORIES(${TIFF_INCLUDE_DIRS})
    
include_directories(${PROJECT_CODE_DIR})
add_executable(mrctest ${PROJECT_CODE_DIR}/Test/MRCReaderTest.cpp)
target_link_libraries(mrctest TomoEngine ${TIFF_LIBRARIES})


endif()


# --------------------------------------------------------------------
# 
# --------------------------------------------------------------------
add_executable(WrappedArrayTest WrappedArrayTest.cpp)
target_link_libraries(WrappedArrayTest MXA TomoEngine)


# --------------------------------------------------------------------
# 
# --------------------------------------------------------------------
if (1)
    add_executable(MRCSubset MRCSubset.cpp)
    target_link_libraries(MRCSubset MXA TomoEngine)
else()
    BuildToolBundle(TARGET MRCSubset
                SOURCES MRCSubset.cpp
                DEBUG_EXTENSION ${EXE_DEBUG_EXTENSION}
                VERSION_MAJOR
                VERSION_MAJOR "1"
                VERSION_MINOR "0"
                VERSION_PATCH "0"
                BINARY_DIR    ${PROJECT_BINARY_DIR}/Test
                COMPONENT     Applications
                INSTALL_DEST  "."
                LINK_LIBRARIES MXA TomoEngine
                LIB_SEARCH_DIRS ${CMAKE_LIBRARY_OUTPUT_DIRECTORY} ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}
)

endif()

