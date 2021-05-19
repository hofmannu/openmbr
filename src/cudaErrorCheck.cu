#include "cudaErrorCheck.cuh"

void cudaErrorCheck::checkCudaReturn(const cudaError_t err, const std::string errMsg)
{
	if (err != cudaSuccess){
		printf("CUDA error string: ");
		printf(cudaGetErrorString(err));
		printf("\n");
		printf(errMsg.c_str());
		printf("\n");
		throw "cudaError";
	}
	return;
}	