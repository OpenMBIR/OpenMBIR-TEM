
include_directories(${PROJECT_CODE_DIR})
add_executable(mrctest ${PROJECT_CODE_DIR}/EIMTomo/IO/MRCReader.cpp ${PROJECT_CODE_DIR}/EIMTomo/Test/MRCReaderTest.cpp)
target_link_libraries(mrctest MXA)


