/*
	kernel for the backwards transformation from a given signal matrix into absorption (transpose)
	Author: Urs Hofmann
	Mail: hofmannu@ethz.ch
	Date: 25.11.2020

	Changelog:
		- moved to temporary variable in code to speed up process
		- flipped access direction for model matrix
		- switched to fmaf for increased accuracy
*/

#include "mbr_transsym.cuh"

__constant__ tran_kernel_arg inArgs_tran_d[1];
__constant__ fwd_kernel_arg inArgs_fwd_d[1];

__global__ void gpuMBTranslationalTranKernel(float* absMat)
{
	// get three dimension index of current output voxel
	const uint64_t xIm = blockIdx.x * blockDim.x + threadIdx.x;
	const uint64_t yIm = blockIdx.y * blockDim.y + threadIdx.y;
	const uint64_t zIm = blockIdx.z * blockDim.z + threadIdx.z;

	if (
		(zIm < inArgs_tran_d->nz) &&
		(xIm < inArgs_tran_d->nxScan) &&
		(yIm < inArgs_tran_d->nyScan)) // if we are in range
	{
		float temp = 0; // gonna abuse this boy for temp saving
		const uint64_t xRange = (inArgs_tran_d->nxModel - 1) / 2; // range in x direction (+/-)
		const uint64_t yRange = (inArgs_tran_d->nyModel - 1) / 2; // range in y direction (+/-)

		// calculate boundaries of xAbs and yAbs
		const uint64_t xMin = (xIm < xRange) ? 0 : (xIm - xRange);
		const uint64_t yMin = (yIm < yRange) ? 0 : (yIm - yRange);
		const uint64_t xMax = ((xIm + xRange) >= inArgs_tran_d->nxScan) ? 
			(inArgs_tran_d->nxScan - 1) : (xIm + xRange);
		const uint64_t yMax = ((yIm + yRange) >= inArgs_tran_d->nyScan) ? 
			(inArgs_tran_d->nyScan - 1) : (yIm + yRange);

		uint64_t sigIdx, modelIdx; // position of model along x and y
		uint64_t xMod, yMod;

		// run over neighbouring a scans
		#pragma unroll 
		for (uint64_t it = 0; it < inArgs_tran_d->nt; it++)
		{
			yMod = 0;
			for (uint64_t yAbs = yMin; yAbs <= yMax; yAbs++)
			{
				xMod = 0;
	 			#pragma unroll
				for (uint64_t xAbs = xMin; xAbs <= xMax; xAbs++)
				{
					// load next line to shared memory
					sigIdx = xAbs + inArgs_tran_d->nxScan * (yAbs + inArgs_tran_d->nyScan * it);
					modelIdx = xMod + inArgs_tran_d->nxModel * (yMod + inArgs_tran_d->nyModel * 
						(it + inArgs_tran_d->nt * zIm)); 
					temp = fmaf(inArgs_tran_d->sigMat[sigIdx], inArgs_tran_d->modelMat[modelIdx], temp);
					xMod++;
				}
				yMod++;
			}
		}

		// find initial absorber point and set it to temp [order: x, y, z]
		const uint64_t outputIdx = xIm + inArgs_tran_d->nxScan * (yIm + zIm * inArgs_tran_d->nyScan);
		absMat[outputIdx] = temp;
	}

	return;
}


// this kernel needs to have x and y flipped
__global__ void gpuMBTranslationalFwdKernel(float* sigMat)
{
	// get three dimension index of current output voxel
	const uint64_t xIm = blockIdx.x * blockDim.x + threadIdx.x;
	const uint64_t yIm = blockIdx.y * blockDim.y + threadIdx.y;
	const uint64_t tIm = blockIdx.z * blockDim.z + threadIdx.z;
	const uint64_t linThread = threadIdx.x + blockDim.x * (threadIdx.y + blockDim.y * threadIdx.z);
	const uint64_t nThreads = blockDim.x * blockDim.y * blockDim.z;
	
	float temp = 0; // gonna abuse this boy for temp saving
	const uint64_t xRange = (inArgs_fwd_d->nxModel - 1) / 2; // range in x direction (+/-)
	const uint64_t yRange = (inArgs_fwd_d->nyModel - 1) / 2; // range in y direction (+/-)

	// calculate boundaries of xAbs and yAbs
	const uint64_t yMin = (yIm < yRange) ? 0 : (yIm - yRange);
	const uint64_t yMax = ((yIm + yRange) >= inArgs_fwd_d->nyScan) ? 
		(inArgs_fwd_d->nyScan - 1) : (yIm + yRange);

	uint64_t absIdx, modelIdx; // position of model along x and y
	uint64_t xModLocal, yMod;

	const uint64_t nxSpan = inArgs_fwd_d->nxModel + blockDim.x - 1;

	// get pointers for shared memory
	extern __shared__ float sModVec[];
	float* absMatVec = &sModVec[inArgs_fwd_d->nxModel]; 

	// run over neighbouring a scans
	#pragma unroll
	for (uint64_t iz = 0; iz < inArgs_fwd_d->nz; iz++)
	{
		yMod = inArgs_fwd_d->nyModel - 1;
		for (uint64_t yAbs = yMin; yAbs <= yMax; yAbs++)
		{
			// load next vector along x from model file
			for (uint64_t xMod = 0; xMod < inArgs_fwd_d->nxModel; xMod = xMod + nThreads)
			{
				const uint64_t xModThread = xMod + linThread;
				if (xModThread < inArgs_fwd_d->nxModel)
				{
					sModVec[xModThread] = inArgs_fwd_d->modelMat[
						xModThread + inArgs_fwd_d->nxModel * (yMod + inArgs_fwd_d->nyModel * 
						(tIm + inArgs_fwd_d->nt * iz))]; 
				}
			}
			__syncthreads();

			// load next vector along y from model file
			for (uint64_t xSpanAbs = 0; xSpanAbs < nxSpan; xSpanAbs = xSpanAbs + nThreads)
			{
				const int64_t xAbsThread = xSpanAbs + linThread;
				if (xAbsThread < nxSpan)
				{
					// overall x index in absorber matrix
					const int64_t xIdxAbs = 
						(int64_t) blockIdx.x * blockDim.x - 
						(int64_t) xRange + xAbsThread;
					if ((xIdxAbs >= 0) && (xIdxAbs < inArgs_fwd_d->nxScan))
					{
						absIdx = xIdxAbs + inArgs_fwd_d->nxScan * (yAbs + inArgs_fwd_d->nyScan * iz);
						absMatVec[xAbsThread] = inArgs_fwd_d->absMat[absIdx]; 
					}
					else
					{
						absMatVec[xAbsThread] = 0.0f;
					}
				}
			}

			__syncthreads();

 			#pragma unroll
			for (int ix = 0; ix < inArgs_fwd_d->nxModel; ix++)
			{
				temp = fmaf(absMatVec[ix + threadIdx.x], sModVec[xModLocal - 1 - ix], temp);
			}
			yMod--;
		}
	}

	// push calculated value to output volume [ix, iy, it] if we are in range
	if ((tIm < inArgs_fwd_d->nt) &&	(xIm < inArgs_fwd_d->nxScan) &&	(yIm < inArgs_fwd_d->nyScan)) 
	{
		const int outputIdx = xIm + inArgs_fwd_d->nxScan * (yIm + tIm * inArgs_fwd_d->nyScan);
		sigMat[outputIdx] = temp;
	}

	return;
}

mbr_transsym::~mbr_transsym()
{

}

void mbr_transsym::run_trans()
{

	printf("Dimensions:\n");
	printf(" - nxScan, nyScan: %lu x %lu \n", nxScan, nyScan);
	printf(" - nxModel, nyModel: %lu x %lu \n", nxModel, nyModel);
	printf(" - nt: %lu \n", nt);
	printf(" - nz: %lu \n", nz);

	printf("Array sizes");

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
	dim3 blockSize(16, 16, 1);
	dim3 gridSize(
		(nxScan + blockSize.x - 1) / blockSize.x,
		(nyScan + blockSize.y - 1) / blockSize.y,
		(nz + blockSize.z - 1) / blockSize.z);

	// run kernel
	const auto tStart = std::chrono::system_clock::now();
	gpuMBTranslationalTranKernel<<< gridSize , blockSize>>>(absMat_dev);

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

	uint64_t nShared = (nxModel + nxModel + blockSize.x - 1) * sizeof(float);

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