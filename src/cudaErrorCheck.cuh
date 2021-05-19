/*
	error checking toolbox
	Author: Urs Hofmann
	Mail: hofmannu@ethz.ch
	Date: 20.05.2021
*/


#ifndef CUDAERRORCHECK_H
#define CUDAERRORCHECK_H

#include <cuda.h>
#include <cuda_runtime.h>
#include <iostream>

class cudaErrorCheck
{
public:
	void checkCudaReturn(const cudaError_t err, const std::string errMsg);
};

#endif