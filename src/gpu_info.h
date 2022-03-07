#ifndef GPU_INFO_H
#define GPU_INFO_H

#include <stdio.h> 
#include <cuda_runtime.h>
#include <cuda.h>
#include <vector>

class gpu_info
{
private:
	int nGpus = 0; // numbe of GPUs detected
	std::vector<cudaDeviceProp> props;
	// name: device name
	// memoryclockrate: printf("Device Number: %d\n", i);
  // printf("  Device name: %s\n", prop.name);
  // printf("  Memory Clock Rate (KHz): %d\n",
  //        prop.memoryClockRate);
  // printf("  Memory Bus Width (bits): %d\n",
  //        prop.memoryBusWidth);
  // printf("  Peak Memory Bandwidth (GB/s): %f\n\n",
  //        2.0*prop.memoryClockRate*(prop.memoryBusWidth/8)/1.0e6);
public:
	gpu_info();
	const cudaDeviceProp* get_props(const int iGpu) const {return &props[iGpu];};
	int get_nGpus() const {return nGpus;};
};

#endif