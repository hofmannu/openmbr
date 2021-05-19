
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
set(CUDA_NVCC_FLAGS "${CUDA_CUDA_FLAGS} -rdc=true -O3 -gencode=arch=compute_75,code=compute_75")


# set(CMAKE_CUDA_FLAGS "${CUDA_CUDA_FLAGS} -rdc=true
# 	-gencode=arch=compute_50,code=sm_50
#   -gencode=arch=compute_52,code=sm_52
#   -gencode=arch=compute_60,code=sm_60
#   -gencode=arch=compute_61,code=sm_61
#   -gencode=arch=compute_70,code=sm_70
#   -gencode=arch=compute_75,code=sm_75
#   -gencode=arch=compute_75,code=compute_75
#   -O2")

find_library(H5CPP_LIB hdf5_cpp)
find_library(H5_LIB hdf5)

add_library(FieldProperties fieldProperties.cpp)

add_library(BasicMathOp basicMathOp.cpp)

add_library(BaseClass baseClass.cpp)

# class representig our reconstruction settings
add_library(ReconSettings reconSettings.cpp)


add_library(Noise noise.cpp)

add_library(Volume volume.cpp)
target_link_libraries(Volume PUBLIC 
	BasicMathOp
	"${H5CPP_LIB}" "${H5_LIB}")

# class which performs modeling of our ultrasound transducer
add_library(TransModel transModel.cu)
set_property(TARGET TransModel PROPERTY CUDA_SEPARABLE_COMPILATION ON)
set_property(TARGET TransModel PROPERTY CUDA_ARCHITECTURES 75)
target_link_libraries(TransModel PUBLIC
	VtkWriter
	TransducerProperties
	FieldProperties
	ModelMatrix
	GaussianBeam)

add_library(TransducerProperties transducerProperties.cpp)
target_link_libraries(TransducerProperties PUBLIC 
	BaseClass
	"${H5CPP_LIB}" "${H5_LIB}")

add_library(ModelMatrix modelMatrix.cpp)

# Class for model based reconstruction
add_library(MbRecon mbrecon.cu)
set_property(TARGET MbRecon PROPERTY CUDA_SEPARABLE_COMPILATION ON)
set_property(TARGET MbRecon PROPERTY CUDA_ARCHITECTURES 75)

target_link_libraries(MbRecon PUBLIC
	CudaErrorCheck
	BaseClass
	ReconSettings
	BasicMathOp
	VtkWriter
	TransModel
	Volume
	ModelMatrix
	GaussianBeam)

# add_library(Simulation simulation.cu)
# set_property(TARGET Simulation PROPERTY CUDA_SEPARABLE_COMPILATION ON)
# set_property(TARGET Simulation PROPERTY CUDA_ARCHITECTURES 35 50 52 60 61 70 72 75)
# target_link_libraries(Simulation PUBLIC
# 	CudaErrorCheck
# 	BaseClass
# 	BasicMathOp
# 	VtkWriter
# 	TransModel
# 	Volume
# 	ModelMatrix
# 	ImGuiPlot 
# 	GLEW 
# 	ImGuiImplSdl 
# 	ImGuiImplOpenGL3
# 	ImGuiFileDialog
# 	MbRecon
# 	ColorMapper
# 	Noise
# 	MC)

add_library(CudaErrorCheck cudaErrorCheck.cu)
set_property(TARGET CudaErrorCheck PROPERTY CUDA_ARCHITECTURES 75)
set_property(TARGET CudaErrorCheck PROPERTY CUDA_SEPARABLE_COMPILATION ON)

add_library(ColorMapper color_mapper.cu)
set_property(TARGET ColorMapper PROPERTY CUDA_SEPARABLE_COMPILATION ON)
set_property(TARGET ColorMapper PROPERTY CUDA_ARCHITECTURES 75)

add_library(Interface interface.cu)
set_property(TARGET Interface PROPERTY CUDA_SEPARABLE_COMPILATION ON)
set_property(TARGET Interface PROPERTY CUDA_ARCHITECTURES 75)
target_link_libraries(Interface PUBLIC
		ReconSettings
		BaseClass
		MbRecon
		TransModel
		TransducerProperties
 		ImGuiPlot 
		GLEW 
		ImGuiImplSdl 
		ImGuiImplOpenGL3
		ImGuiFileDialog
		ColorMapper
		GaussianBeam
		)

# modified imgui related libraries
add_library(ImGuiFileDialog ImGuiFileDialog.cpp)
target_link_libraries(ImGuiFileDialog PUBLIC ImGui)
add_library(ImGuiPlot imgui_plot.cpp)

add_library(ImGuiImplOpenGL3 imgui_impl_opengl3.cpp)
target_link_libraries(ImGuiImplOpenGL3 PUBLIC ImGui OpenGL)

add_library(GaussianBeam gaussian_beam.cpp)

find_package(SDL2 REQUIRED)
add_library(ImGuiImplSdl imgui_impl_sdl.cpp)
target_link_libraries(ImGuiImplSdl PUBLIC ImGui SDL2 OpenGL)

find_package(OpenGL REQUIRED)
include_directories(${OPENGL_INCLUDE_DIRS})

find_package(glfw3 REQUIRED)
include_directories(${GLFW3_INCLUDE_DIR})
link_libraries(${GLFW_LIBRARY_DIRS})

