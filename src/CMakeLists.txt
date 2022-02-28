set(OpenGL_GL_PREFERENCE GLVND)
cmake_policy(SET CMP0079 NEW)

# prepare for cuda compilation
enable_language(CUDA)
find_package(CUDA REQUIRED QUIET)
INCLUDE(FindCUDA)
INCLUDE_DIRECTORIES(/opt/cuda/include include)
set(CUDA_SEPARABLE_COMPILATION ON)
set(CUDA_ATTACH_VS_BUILD_RULE_TO_CUDA_FILE OFF)
set(BUILD_SHARED_LIBS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
# set(CMAKE_CUDA_FLAGS "${CUDA_CUDA_FLAGS} -rdc=true -ftz=false -prec-div=true -prec-sqrt=true")
set(CMAKE_CUDA_FLAGS "${CUDA_CUDA_FLAGS} -rdc=true --ptxas-options=-v -O2 -use_fast_math")


find_library(H5CPP_LIB hdf5_cpp)
find_library(H5_LIB hdf5)

add_library(Reconsett reconsett.cpp)

add_library(Model model.cpp)
target_link_libraries(Model PUBLIC
	"${H5CPP_LIB}" 
	"${H5_LIB}" 
	Volume)

add_library(Recon recon.cpp)
target_link_libraries(Recon PUBLIC
	Reconsett
	Model
	Volume
	MbrTranssym)
set_property(TARGET Recon PROPERTY CUDA_SEPARABLE_COMPILATION ON)

add_library(ColorMapper color_mapper.cpp)

add_library(Gui gui.cpp)
target_link_libraries(Gui PUBLIC
	ImGuiPlot 
	GLEW 
	ImGuiImplSdl 
	ImGuiImplOpenGL3
	ImGuiFileDialog
	Recon
	Volume
	ColorMapper
	)
set_property(TARGET Gui PROPERTY CUDA_SEPARABLE_COMPILATION ON)


add_library(ImGuiFileDialog ImGuiFileDialog.cpp)
target_link_libraries(ImGuiFileDialog PUBLIC ImGui)

add_library(ImGuiPlot imgui_plot.cpp)

add_library(ImGuiImplOpenGL3 imgui_impl_opengl3.cpp)
target_link_libraries(ImGuiImplOpenGL3 PUBLIC ImGui OpenGL)

find_package(SDL2 REQUIRED)
add_library(ImGuiImplSdl imgui_impl_sdl.cpp)
target_link_libraries(ImGuiImplSdl PUBLIC ImGui SDL2 OpenGL)

find_package(OpenGL REQUIRED)
include_directories(${OPENGL_INCLUDE_DIRS})

find_package(glfw3 REQUIRED)
include_directories(${GLFW3_INCLUDE_DIR})
link_libraries(${GLFW_LIBRARY_DIRS})

add_library(MbrTranssym mbr_transsym.cu)
set_property(TARGET MbrTranssym PROPERTY CUDA_SEPARABLE_COMPILATION ON)
set_property(TARGET MbrTranssym PROPERTY CUDA_ARCHITECTURES 75)