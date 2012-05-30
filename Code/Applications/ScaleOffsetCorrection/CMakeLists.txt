project(ScaleOffsetCorrection)
cmake_minimum_required(VERSION 2.8.5)
set (EXECUTABLE_NAME soc)

set(BiSphere_Input_File "" CACHE PATH "")

set(ScaleOffsetCorrection_PROJECT_FILES
    ${ScaleOffsetCorrection_SOURCE_DIR}/SOCMain.cpp
    ${ScaleOffsetCorrection_SOURCE_DIR}/SOCArgsParser.cpp
)


add_executable(${EXECUTABLE_NAME} ${ScaleOffsetCorrection_PROJECT_FILES} )
target_link_libraries(${EXECUTABLE_NAME} TomoEngine)


configure_file(${EIMTomo_SOURCE_DIR}/Resources/BiSphereTest/BiSphereTest.sh.in
               ${EIMTomo_BINARY_DIR}/BiSphereTest.sh
               )
        

