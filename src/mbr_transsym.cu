/*
	kernel for the backwards transformation from a given signal matrix into absorption (transpose)
	Author: Urs Hofmann
	Mail: hofmannu@ethz.ch
	Date: 25.11.2020

	Changelog:
		- moved to temporary variable in code to speed up process
		- flipped access direction for model matrix
		- switched to fmaf for increased accuracy
		- moved kernel definition to a separate file
		- moved forward kernel to shared memory
*/

#include "mbr_transsym.cuh"
#include "mbr_transsym_kernel.cu"

// class constructor and destructor
mbr_transsym::~mbr_transsym()
{

}

// starts the execution of the transpose kernel multiplication
void mbr_transsym::run_trans()
{
	// allocate memory on GPU
	cudaErr = cudaMalloc((void**) &modelMat_dev, get_nBytesModel());
	cudaErrCheck("Could not allocate memory for model matrix");

	cudaErr = cudaMalloc( (void**) &absMat_dev, get_nBytesAbs());
	cudaErrCheck("Could not allocate memory for absorber matrix");

	cudaErr = cudaMalloc( (void**) &sigMat_dev, get_nBytesSignal());
	cudaErrCheck("Could not allocate memory for signal matrix");

	cudaErr = cudaMemcpy(sigMat_dev, sigMat, get_nBytesSignal(), cudaMemcpyHostToDevice);
	cudaErrCheck("Could not copy memory for signal matrix");

	cudaErr = cudaMemcpy(modelMat_dev, modelMat, get_nBytesModel(), cudaMemcpyHostToDevice);
	cudaErrCheck("Could not copy memory for model to GPU");

	// define input arguments for kernel execution and copy to GPU
	tran_kernel_arg inArgs_h;
	inArgs_h.nxScan = nxScan;
	inArgs_h.nyScan = nyScan;
	inArgs_h.nz = nz;
	inArgs_h.nt = nt;
	inArgs_h.nxModel = nxModel;
	inArgs_h.nyModel = nyModel;
	inArgs_h.sigMat = sigMat_dev;
	inArgs_h.modelMat = modelMat_dev;

	cudaErr = cudaMemcpyToSymbol(inArgs_tran_d, &inArgs_h, sizeof(tran_kernel_arg));
	cudaErrCheck("Could not copy symbol to GPU");

	// define block and grid size for kernel
	dim3 blockSize(64, 1, 1);
	dim3 gridSize(
		(nxScan + blockSize.x - 1) / blockSize.x,
		(nyScan + blockSize.y - 1) / blockSize.y,
		(nz + blockSize.z - 1) / blockSize.z);

	// calculate amount of shared memory required for this operation
	const uint64_t nShared = (nxModel + nxModel + blockSize.x - 1) * sizeof(float);

	// run kernel
	const auto tStart = std::chrono::system_clock::now();
	gpuMBTranslationalTranKernel<<< gridSize , blockSize, nShared>>>(absMat_dev);

	cudaDeviceSynchronize(); // wait for kernel to finish
	const auto tStop = std::chrono::system_clock::now();
	const auto duration = std::chrono::duration_cast<std::chrono::milliseconds>
		(tStop - tStart);
	tExec = duration.count();

	cudaErr = cudaGetLastError();
	cudaErrCheck("Error during kernel execution");
	
	// copy matrix back from GPU
	cudaErr = cudaMemcpy(absMat, absMat_dev, get_nBytesAbs(), cudaMemcpyDeviceToHost);
	cudaErrCheck("Error while copying data back from device");
	
	// free device memory
	cudaFree(modelMat_dev);
	cudaFree(absMat_dev);
	cudaFree(sigMat_dev);
	return;
}

void mbr_transsym::run_fwd()
{

	// allocate memory on GPU
	cudaErr = cudaMalloc((void**) &modelMat_dev, get_nBytesModel());
	cudaErrCheck("Could not allocate memory for model matrix");

	cudaErr = cudaMalloc( (void**) &absMat_dev, get_nBytesAbs());
	cudaErrCheck("Could not allocate memory for absorber matrix");

	cudaErr = cudaMalloc( (void**) &sigMat_dev, get_nBytesSignal());
	cudaErrCheck("Could not allocate memory for signal matrix");

	cudaErr = cudaMemcpy(absMat_dev, absMat, get_nBytesAbs(), cudaMemcpyHostToDevice);
	cudaErrCheck("Could not copy memory for signal matrix");

	cudaErr = cudaMemcpy(modelMat_dev, modelMat, get_nBytesModel(), cudaMemcpyHostToDevice);
	cudaErrCheck("Could not copy memory for model to GPU");

	fwd_kernel_arg inArgs_h;
	inArgs_h.nxScan = nxScan;
	inArgs_h.nyScan = nyScan;
	inArgs_h.nz = nz;
	inArgs_h.nt = nt;
	inArgs_h.nxModel = nxModel;
	inArgs_h.nyModel = nyModel;
	inArgs_h.absMat = absMat_dev;
	inArgs_h.modelMat = modelMat_dev;

	cudaErr = cudaMemcpyToSymbol(inArgs_fwd_d, &inArgs_h, sizeof(fwd_kernel_arg));
	cudaErrCheck("Could not copy symbol to GPU");

	// define block and grid size for kernel
	dim3 blockSize(64, 1, 1);
	dim3 gridSize(
		(nxScan + blockSize.x - 1) / blockSize.x,
		(nyScan + blockSize.y - 1) / blockSize.y,
		(nt + blockSize.z - 1) / blockSize.z);

	// calculate amount of shared memory required for this operation
	const uint64_t nShared = (nxModel + nxModel + blockSize.x - 1) * sizeof(float);

	// run kernel
	const auto tStart = std::chrono::system_clock::now();
	gpuMBTranslationalFwdKernel<<< gridSize , blockSize, nShared>>>(sigMat_dev);

	cudaDeviceSynchronize(); // wait for kernel to finish
	const auto tStop = std::chrono::system_clock::now();
	const auto duration = std::chrono::duration_cast<std::chrono::milliseconds>
		(tStop - tStart);
	tExec = duration.count();

	cudaErr = cudaGetLastError();
	cudaErrCheck("Error during kernel execution");
	
	// copy matrix back from GPU
	cudaErr = cudaMemcpy(sigMat, sigMat_dev, get_nBytesSignal(), cudaMemcpyDeviceToHost);
	cudaErrCheck("Error while copying data back from device");
	
	// free device memory
	cudaFree(modelMat_dev);
	cudaFree(absMat_dev);
	cudaFree(sigMat_dev);
	return;
}

// define size of signal matrix
void mbr_transsym::set_sizeSigMat(
	const uint64_t _nxScan, const uint64_t _nyScan, const uint64_t _ntScan)
{
	nxScan = _nxScan;
	nyScan = _nyScan;
	if (nt == 0)
	{
		nt = _ntScan;
	}
	else
	{
		if (nt != _ntScan)
		{
			printf("Size of model and data needs to align along time.\n");
			throw "InvalidValue";
		}
	}
	return;
}

// define size of absorber matrix
void mbr_transsym::set_sizeAbsMat(
	const uint64_t _nxAbs, const uint64_t _nyAbs, const uint64_t _nzAbs)
{
	nxScan = _nxAbs;
	nyScan = _nxAbs;
	if (nz == 0)
	{
		nz = _nzAbs;
	}
	else
	{
		if (nz != _nzAbs)
		{
			printf("Size of model and data needs to align along z.\n");
			throw "InvalidValue";
		}
	}
	return;
}

// define size of model matrix
void mbr_transsym::set_sizeModel(
		const uint64_t _nxModel, const uint64_t _nyModel, 
		const uint64_t _ntModel, const uint64_t _nz)
{
	if ((_nxModel % 2) == 0)
	{
		printf("Model matrix must be uneven along x");
		throw "InvalidValue";
	}
	else
	{
		nxModel = _nxModel;
	}

	if ((_nyModel % 2) == 0)
	{
		printf("Model matrix must be uneven along y");
		throw "InvalidValue";
	}
	else
	{
		nyModel = _nyModel;
	}

	if (nt == 0)
	{
		nt = _ntModel;
	}
	else
	{
		if (nt != _ntModel)
		{
			printf("Size of model and data needs to align along time.\n");
			throw "InvalidValue";
		}
	}

	nz = _nz;

	return;
}

// helper function to check CUDA errors
void mbr_transsym::cudaErrCheck(const char* errMsg)
{
	if (cudaErr != cudaSuccess)
	{
		printf("Some error occured: %s, %s\n", errMsg, cudaGetErrorString(cudaErr));
		throw "CudaError";
	}
	return;
}

// prints some basic information about the kernel prior execution
void mbr_transsym::print_info()
{
	printf("Signal size [x, y, t]: %d, %d, %d\n",
		(int) nxScan, (int) nyScan, (int) nt);
	printf("Absorber size [x, y, z]: %d, %d, %d\n", 
		(int) nxScan, (int) nyScan, (int) nz);
	printf("Model size [x, y, t, z]: %d, %d, %d, %d\n", 
		(int) nxModel, (int) nyModel, (int) nt, (int) nz);
	return;
}