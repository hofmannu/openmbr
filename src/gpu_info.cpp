#include "gpu_info.h"


gpu_info::gpu_info()
{

	cudaGetDeviceCount(&nGpus);
	for (int i = 0; i < nGpus; i++) {
    cudaDeviceProp currProp;
    cudaGetDeviceProperties(&currProp, i);
    props.push_back(currProp);
  }

	return;
}