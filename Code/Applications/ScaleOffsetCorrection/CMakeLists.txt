project(ScaleOffsetCorrection)
cmake_minimum_required(VERSION 2.8.5)


set(ScaleOffsetCorrection_PROJECT_FILES
    ${ScaleOffsetCorrection_SOURCE_DIR}/SOCMain.cpp
    ${ScaleOffsetCorrection_SOURCE_DIR}/SOCArgsParser.cpp
)


add_executable(soc ${ScaleOffsetCorrection_PROJECT_FILES} )
target_link_libraries(soc TomoEngine)