set(OpenGL_GL_PREFERENCE GLVND)
cmake_policy(SET CMP0079 NEW)

# clean imgui libraries
add_library(ImGui imgui/imgui.cpp)
target_link_libraries(ImGui PUBLIC ImGuiWidgets) 

add_library(ImGuiWidgets imgui/imgui_widgets.cpp)
target_link_libraries(ImGuiWidgets PUBLIC ImGuiDraw)

add_library(ImGuiDraw imgui/imgui_draw.cpp)

# add_library(FiberProperties arillumination/GUI/src/fiberProperties.cpp)
# set_property(TARGET FiberProperties PROPERTY CUDA_SEPARABLE_COMPILATION ON)

# add_library(McFieldProperties arillumination/GUI/src/mcFieldProperties.cpp)
# set_property(TARGET McFieldProperties PROPERTY CUDA_SEPARABLE_COMPILATION ON)

# add_library(OptProperties arillumination/GUI/src/optProperties.cpp)
# set_property(TARGET OptProperties PROPERTY CUDA_SEPARABLE_COMPILATION ON)

# add_library(SimProperties arillumination/GUI/src/simProperties.cpp)
# set_property(TARGET SimProperties PROPERTY CUDA_SEPARABLE_COMPILATION ON)

add_library(VtkWriter vtkwriter/vtkwriter.cpp)

# add_library(MC arillumination/GUI/src/mc.cu)
# set_property(TARGET MC PROPERTY CUDA_SEPARABLE_COMPILATION ON)
# set_property(TARGET MC PROPERTY CUDA_ARCHITECTURES 35 50)
# target_link_libraries(MC PUBLIC
# 	FiberProperties
# 	McFieldProperties
# 	OptProperties
# 	SimProperties
# 	VtkWriter)