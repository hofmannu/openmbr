cmake_minimum_required(VERSION 3.12)

# set project name
project(FluenceGUI 	VERSION 1.0
	DESCRIPTION "Towards a GUI for Monte Carlo simulations"
	LANGUAGES CUDA CXX)

# prepare for cuda compilation
enable_language(CUDA)
find_package(CUDA REQUIRED QUIET)
INCLUDE(FindCUDA)
INCLUDE_DIRECTORIES(/opt/cuda/include include)

set(CUDA_ATTACH_VS_BUILD_RULE_TO_CUDA_FILE OFF)
set(BUILD_SHARED_LIBS OFF)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON )
set(
    CUDA_NVCC_FLAGS
		${CUDA_NVCC_FLAGS};
		-x cu -rdc=true -g
)

add_subdirectory(src/)
add_subdirectory(lib/)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED true)
set(OpenGL_GL_PREFERENCE GLVND)
cmake_policy(SET CMP0079 NEW)
add_compile_options("-g3 -ggdb")

find_library(TERM_LIB curses)

add_executable(main src/main.cu)
set_property(TARGET main PROPERTY CUDA_SEPARABLE_COMPILATION ON)
set_property(TARGET main PROPERTY CUDA_ARCHITECTURES 35 50)
target_link_libraries(main PUBLIC
	MC
	Interface)

install(TARGETS main DESTINATION bin)
