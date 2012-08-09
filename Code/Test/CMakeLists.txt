

if (EIMTomo_USE_TIFF)  
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

