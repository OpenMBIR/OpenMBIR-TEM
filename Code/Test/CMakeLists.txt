

if (OpenMBIR_USE_TIFF)  
    include_directories(${PROJECT_CODE_DIR})
    add_executable(mrctest ${PROJECT_CODE_DIR}/Test/MRCReaderTest.cpp)
    target_link_libraries(mrctest MXA ReconstructionCoreLib HaadfMbirLib ${TIFF_LIBRARIES})
endif()


# --------------------------------------------------------------------
# 
# --------------------------------------------------------------------
add_executable(WrappedArrayTest WrappedArrayTest.cpp)
target_link_libraries(WrappedArrayTest MXA ReconstructionCoreLib HaadfMbirLib)


