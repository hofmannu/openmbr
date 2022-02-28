
add_executable(UtestKernelTrans utest_kernelTrans.cpp)
target_link_libraries(UtestKernelTrans PUBLIC MbrTranssym)
set_property(TARGET UtestKernelTrans PROPERTY CUDA_SEPARABLE_COMPILATION ON)

add_executable(UtestKernelFwd utest_kernelFwd.cpp)
target_link_libraries(UtestKernelFwd PUBLIC MbrTranssym)
set_property(TARGET UtestKernelFwd PROPERTY CUDA_SEPARABLE_COMPILATION ON)

add_executable(UtestSettings utest_reconsettings.cpp)
target_link_libraries(UtestSettings PUBLIC Reconsett)

add_executable(UtestModel utest_model.cpp)
target_link_libraries(UtestModel PUBLIC Model)

add_executable(UtestRecon utest_recon.cpp)
target_link_libraries(UtestRecon PUBLIC Recon Model Reconsett)
set_property(TARGET UtestRecon PROPERTY CUDA_SEPARABLE_COMPILATION ON)