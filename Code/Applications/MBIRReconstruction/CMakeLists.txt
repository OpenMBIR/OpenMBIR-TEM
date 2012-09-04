project(MbirReconstruction)

cmake_minimum_required(VERSION 2.8.5)

set (EXECUTABLE_NAME MbirReconstruction)

set(BiSphere_Input_File "" CACHE PATH "")

set(MbirReconstruction_PROJECT_FILES
    ${MbirReconstruction_SOURCE_DIR}/main.cpp
    ${MbirReconstruction_SOURCE_DIR}/MBIRReconstructionArgsParser.cpp
    ${MbirReconstruction_SOURCE_DIR}/MBIRReconstructionArgsParser.h
)


add_executable(${EXECUTABLE_NAME} ${MbirReconstruction_PROJECT_FILES} )
target_link_libraries(${EXECUTABLE_NAME} TomoEngine)


configure_file(${EIMTomo_SOURCE_DIR}/Resources/BiSphereTest/BiSphereTest.sh.in
               ${EIMTomo_BINARY_DIR}/BiSphereTest.sh
               )
        

