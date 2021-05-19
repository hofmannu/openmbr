cmake_minimum_required(VERSION 3.12)

# set project name
project(MB-GUI 	VERSION 1.0
	DESCRIPTION "Towards a GUI for model based AR reconstructions"
	LANGUAGES CUDA CXX)

# prepare for cuda compilation
enable_language(CUDA)
find_package(CUDA REQUIRED QUIET)
INCLUDE(FindCUDA)
INCLUDE_DIRECTORIES(/opt/cuda/include include)

set(CUDA_ATTACH_VS_BUILD_RULE_TO_CUDA_FILE OFF)
set(BUILD_SHARED_LIBS OFF)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON )
# set(CMAKE_CUDA_FLAGS "${CUDA_CUDA_FLAGS} -rdc=true -ftz=false -prec-div=true -prec-sqrt=true")
set(CMAKE_CUDA_FLAGS "${CUDA_CUDA_FLAGS} -rdc=true -Xptxas -v -O2")

set(
    CUDA_NVCC_FLAGS
		${CUDA_NVCC_FLAGS};
		-x cu -rdc=true -ftz=false -prec-div=true -prec-sqrt=true -g
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
set_property(TARGET main PROPERTY CUDA_ARCHITECTURES 75)
target_link_libraries(main PUBLIC 
	MbRecon TransducerProperties "${TERM_LIB}")

add_executable(gui src/gui.cu)
set_property(TARGET gui PROPERTY CUDA_SEPARABLE_COMPILATION ON)
set_property(TARGET gui PROPERTY CUDA_ARCHITECTURES 75)
target_link_libraries(gui PUBLIC 
	Interface
	"${TERM_LIB}")

# add_executable(sim src/sim.cu)
# set_property(TARGET sim PROPERTY CUDA_SEPARABLE_COMPILATION ON)
# set_property(TARGET sim PROPERTY CUDA_ARCHITECTURES 75)
# target_link_libraries(sim PUBLIC 
# 	Simulation
# 	"${TERM_LIB}")

add_executable(test_kernelspeed unit_tests/test_kernelspeed.cu)
set_property(TARGET test_kernelspeed PROPERTY CUDA_SEPARABLE_COMPILATION ON)
set_property(TARGET test_kernelspeed PROPERTY CUDA_ARCHITECTURES 75)
target_link_libraries(test_kernelspeed PUBLIC 
	MbRecon
	TransducerProperties
	Volume
	ReconSettings
	"${TERM_LIB}")


target_include_directories(main PUBLIC
	"${PROJECT_BINARY_DIR}"
	"${PROJECT_SOURCE_DIR}"
	)

install(TARGETS main DESTINATION bin)
