
__constant__ tran_kernel_arg inArgs_tran_d[1];
__constant__ fwd_kernel_arg inArgs_fwd_d[1];

__global__ void gpuMBTranslationalTranKernel(float* absMat)
{
	// get three dimension index of current output voxel
	const uint64_t xIm = blockIdx.x * blockDim.x + threadIdx.x;
	const uint64_t yIm = blockIdx.y * blockDim.y + threadIdx.y;
	const uint64_t zIm = blockIdx.z * blockDim.z + threadIdx.z;
	const uint64_t linThread = 
		threadIdx.x + blockDim.x * (threadIdx.y + blockDim.y * threadIdx.z);
	const uint64_t nThreads = blockDim.x * blockDim.y * blockDim.z;


	float temp = 0; // gonna abuse this boy for temp saving
	const uint64_t xRange = (inArgs_tran_d->nxModel - 1) / 2; // range in x direction (+/-)
	const uint64_t yRange = (inArgs_tran_d->nyModel - 1) / 2; // range in y direction (+/-)

	// calculate boundaries of xAbs and yAbs
	const uint64_t yMin = (yIm < yRange) ? 0 : (yIm - yRange);
	const uint64_t yMax = ((yIm + yRange) >= inArgs_tran_d->nyScan) ? 
		(inArgs_tran_d->nyScan - 1) : (yIm + yRange);
	const uint64_t yStart = (yIm < yRange) ? (yRange - yIm) : 0;

	uint64_t sigIdx, modelIdx; // position of model along x and y
	uint64_t xModLocal, yMod;

	const uint64_t nxSpan = inArgs_tran_d->nxModel + blockDim.x - 1;

	// get pointers for shared memory
	extern __shared__ float sModVec[];
	float* sigMatVec = &sModVec[inArgs_tran_d->nxModel]; 

	#pragma unroll 
	for (uint64_t it = 0; it < inArgs_tran_d->nt; it++)
	{
		yMod = yStart;
		for (uint64_t yAbs = yMin; yAbs <= yMax; yAbs++)
		{
			// load next vector along x from model file
			for (uint64_t xMod = 0; xMod < inArgs_tran_d->nxModel; xMod = xMod + nThreads)
			{
				const uint64_t xModThread = xMod + linThread;
				if (xModThread < inArgs_tran_d->nxModel)
				{
					sModVec[xModThread] = inArgs_tran_d->modelMat[
						xModThread + inArgs_tran_d->nxModel * (yMod + inArgs_tran_d->nyModel * 
						(it + inArgs_tran_d->nt * zIm))]; 
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
					const int64_t xIdxSig =	(int64_t) blockIdx.x * blockDim.x - 
						(int64_t) xRange + xAbsThread;
					if ((xIdxSig >= 0) && (xIdxSig < inArgs_tran_d->nxScan))
					{
						sigIdx = xIdxSig + inArgs_tran_d->nxScan * (yAbs + inArgs_tran_d->nyScan * it);
						sigMatVec[xAbsThread] = inArgs_tran_d->sigMat[sigIdx]; 
					}
					else
					{
						sigMatVec[xAbsThread] = 0.0f;
					}
				}
			}

			__syncthreads();

 			#pragma unroll
			for (int ix = 0; ix < inArgs_tran_d->nxModel; ix++)
			{
				temp = fmaf(sigMatVec[ix + threadIdx.x], sModVec[ix], temp);
			}

			yMod++;
		}
	}



	// thread writes only to output volume if we are in range
	if ((zIm < inArgs_tran_d->nz) && (xIm < inArgs_tran_d->nxScan) && (yIm < inArgs_tran_d->nyScan)) // if we are in range
	{
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
	const uint64_t linThread = 
		threadIdx.x + blockDim.x * (threadIdx.y + blockDim.y * threadIdx.z);
	const uint64_t nThreads = blockDim.x * blockDim.y * blockDim.z;
	
	float temp = 0; // gonna abuse this boy for temp saving
	const uint64_t xRange = (inArgs_fwd_d->nxModel - 1) / 2; // range in x direction (+/-)
	const uint64_t yRange = (inArgs_fwd_d->nyModel - 1) / 2; // range in y direction (+/-)

	// calculate boundaries of xAbs and yAbs
	const uint64_t yMin = (yIm < yRange) ? 0 : (yIm - yRange);
	const uint64_t yMax = ((yIm + yRange) >= inArgs_fwd_d->nyScan) ? 
		(inArgs_fwd_d->nyScan - 1) : (yIm + yRange);
	const uint64_t yStart = (yIm < yRange) ? 
		(inArgs_fwd_d->nyModel - 1 - (yRange - yIm)) : 
		inArgs_fwd_d->nyModel - 1;

	uint64_t absIdx, modelIdx; // position of model along x and y
	uint64_t xModLocal, yMod;

	const uint64_t nxSpan = inArgs_fwd_d->nxModel + blockDim.x - 1;

	// get pointers for shared memory
	extern __shared__ float sModVec[];
	float* absMatVec = &sModVec[inArgs_fwd_d->nxModel]; 

	#pragma unroll
	for (uint64_t iz = 0; iz < inArgs_fwd_d->nz; iz++)
	{
		yMod = yStart;
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
