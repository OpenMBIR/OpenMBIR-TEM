

project(Reconstruct)





set(Reconstruct_HDRS
  ${PROJECT_CODE_DIR}/allocate.h
  ${PROJECT_CODE_DIR}/ComputationEngineOptimized.h
  ${PROJECT_CODE_DIR}/ComputationInputs.h
)
set(Reconstruct_SRCS
  ${PROJECT_CODE_DIR}/allocate.c
  ${PROJECT_CODE_DIR}/ComputationEngineOptimized.c
  ${PROJECT_CODE_DIR}/ComputationInputs.c
  ${PROJECT_CODE_DIR}/Main.c
)


set(Reconstruct_PROJECT_FILES
  ${Reconstruct_HDRS}
  ${Reconstruct_SRCS}
)

add_executable(reconstruct ${Reconstruct_PROJECT_FILES} )


