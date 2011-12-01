
option (EIMTomo_USE_TIFF "Use libTiff for some image output" OFF)
if (EIMTomo_USE_TIFF)
 # Prepend our own "FindTiff" code in front of what CMake itself provides
    set(CMAKE_MODULE_PATH ${CMP_MODULES_SOURCE_DIR} ${CMAKE_MODULE_PATH} )
    FIND_PACKAGE(Tiff)
    IF (NOT TIFF_FOUND)
          MESSAGE (FATAL_ERROR " Tiff Library was NOT Found and is needed.")
    ENDIF (NOT TIFF_FOUND)
    INCLUDE_DIRECTORIES(${TIFF_INCLUDE_DIRS})
    
include_directories(${PROJECT_CODE_DIR})
add_executable(mrctest ${PROJECT_CODE_DIR}/EIMTomo/IO/MRCReader.cpp 
                        ${PROJECT_CODE_DIR}/EIMTomo/Test/MRCReaderTest.cpp)
target_link_libraries(mrctest MXA ${TIFF_LIBRARIES})


endif()

