#include <mex.h>
#include <cuda.h>

struct constInArgs{
	float* surface;
	float* fluenceField;
	float dr; // resolution of fluence field in radial direction
	float dx;
	float dy;
	unsigned int nx; // number of elements in x direction of surface
	unsigned int ny; // number of elements in y direction of surface
	unsigned int nz; // number of elements in z direction of volume
	unsigned int nr; // number of elements in radial direction of fluence field
};

__global__ void Conv_Fluence(float* fluenceVolume, const constInArgs* inArgs){
	const unsigned int izOut = blockIdx.x * blockDim.x + threadIdx.x;
	const unsigned int ixOut = blockIdx.y * blockDim.y + threadIdx.y;
	const unsigned int iyOut = blockIdx.z * blockDim.z + threadIdx.z;
	
	if((izOut < inArgs->nz) && (ixOut < inArgs->nx) && (iyOut < inArgs->ny)){
		const unsigned int outputIdx = izOut + inArgs->nz * (ixOut + iyOut * inArgs->nx);
		fluenceVolume[outputIdx] = 0;
		
		const float xPos = (float) ixOut * inArgs->dx;
		const float yPos = (float) iyOut * inArgs->dy;
		const float rRange = (float) inArgs->nr * inArgs->dr;

		const float xMin = xPos - rRange;
		int xIdxMin = xMin / inArgs->dx + 0.5;
		if (xIdxMin < 0)
			xIdxMin = 0;
		const float xMax = xPos + rRange;
		int xIdxMax = xMax / inArgs->dx + 0.5;
		if (xIdxMax >= inArgs->nx)
			xIdxMax = inArgs->nx - 1;

		const float yMin = yPos - rRange;
		int yIdxMin = yMin / inArgs->dy + 0.5;
		if (yIdxMin < 0)
			yIdxMin = 0;
		const float yMax = yPos + rRange;
		int yIdxMax = yMax / inArgs->dy + 0.5;
		if (yIdxMax >= inArgs->ny)
			yIdxMax = inArgs->ny - 1;	

		float xDist, yDist, rDist;
		int rIdx, fluenceIdx, surfIdx;
		for (int ix = xIdxMin; ix <= xIdxMax; ix++){
			xDist = (float) ix * inArgs->dx - xPos;
			for (int iy = yIdxMin; iy <= yIdxMax; iy++){
				yDist = (float) iy * inArgs->dy - yPos;
				// get radial distance
				rDist = sqrt(xDist * xDist + yDist * yDist);
				rIdx = rDist / inArgs->dr + 0.5;
				if (rIdx < inArgs->nr){
					// get surface index
					surfIdx = inArgs->surface[ix + iy * inArgs->nx];	
					fluenceIdx = izOut + rIdx * inArgs->nz + surfIdx * inArgs->nz * inArgs->nr;
					fluenceVolume[outputIdx] += inArgs->fluenceField[fluenceIdx];
				}
			}
		}
	}

	return;
}

void mexFunction(
	int nlhs,
	mxArray *plhs[],
	int nrhs,
	const mxArray *prhs[]
		)
{

	constInArgs inputArguments;

	const float * surface = (float*) mxGetData(prhs[0]);
	inputArguments.nx = (unsigned int) mxGetDimensions(prhs[0])[0];
	inputArguments.ny = (unsigned int) mxGetDimensions(prhs[0])[1];
	
	const float * fluenceField = (float*) mxGetData(prhs[1]);	
	inputArguments.nz = (unsigned int) mxGetDimensions(prhs[1])[0];
	inputArguments.nr = (unsigned int) mxGetDimensions(prhs[1])[1];

	const float * res = (float*) mxGetData(prhs[2]);	
	inputArguments.dr = res[0];
	inputArguments.dx = res[1];
	inputArguments.dy = res[2];

	// prepare and allocate array for output volume in matlab
	const mwSize dims[3] = {inputArguments.nz, inputArguments.nx, inputArguments.ny};
	mxArray* volArray = mxCreateNumericArray(3, dims, mxSINGLE_CLASS, mxREAL);
	plhs[0] = volArray;
	float* outData = (float*) mxGetData(volArray);

	constInArgs* inputArguments_dev;
// allocate memory on gpu
	float* fluenceVolume_dev;
	bool m1 = (cudaSuccess != cudaMalloc( (void**) &inputArguments.surface,
		inputArguments.nx * inputArguments.ny * sizeof(float) )); 
	m1 |= (cudaSuccess != cudaMalloc( (void**) &fluenceVolume_dev, 
		inputArguments.nx * inputArguments.ny * inputArguments.nz * sizeof(float)));
	m1 |= (cudaSuccess != cudaMalloc( (void**) &inputArguments_dev, 
		sizeof(constInArgs) ));
	m1 |= (cudaSuccess != cudaMalloc( (void**) &inputArguments.fluenceField,
		inputArguments.nz * inputArguments.nz * inputArguments.nr * sizeof(float) ));

	if (m1){
		mexErrMsgTxt("Could not allocate memory on device");
	}else{
		// copy data over to GPU
		bool cpy1 = (cudaSuccess != cudaMemcpy(inputArguments_dev, &inputArguments,
			sizeof(constInArgs), cudaMemcpyHostToDevice ));
		cpy1 |= (cudaSuccess != cudaMemcpy(inputArguments.surface, surface,
			inputArguments.nx * inputArguments.ny * sizeof(float), cudaMemcpyHostToDevice));
		cpy1 |= (cudaSuccess != cudaMemcpy(inputArguments.fluenceField, fluenceField,
			inputArguments.nz * inputArguments.nz * inputArguments.nr * sizeof(float), 
			cudaMemcpyHostToDevice));
		if (cpy1){
			mexErrMsgTxt("Could not copy arrays to device");
		}else{
			// start kernel execution
			const dim3 blockSize(32, 4, 4); // z, x, y
			const dim3 gridSize(
				(inputArguments.nz + blockSize.x - 1) / blockSize.x,
				(inputArguments.nx + blockSize.y - 1) / blockSize.y,
				(inputArguments.ny + blockSize.z - 1) / blockSize.z);

			Conv_Fluence<<< gridSize, blockSize>>>(fluenceVolume_dev, inputArguments_dev);
			cudaDeviceSynchronize();
			
			// check if kernel execution was successfull 
			cudaError_t err = cudaGetLastError();
			if (err != cudaSuccess){
				mexErrMsgTxt(cudaGetErrorString(err));
			}else{
				// copy data back from GPU and push to MATLAB
				err = cudaMemcpy(outData, fluenceVolume_dev, 
					inputArguments.nz * inputArguments.nz * inputArguments.nr * sizeof(float),
					cudaMemcpyDeviceToHost);
				if (err != cudaSuccess){
					mexErrMsgTxt(cudaGetErrorString(err));
				}
			}
		}
		// free device memory
		cudaFree(inputArguments.surface);
		cudaFree(fluenceVolume_dev);
		cudaFree(inputArguments_dev);
	
		return;	
	}


	return;
};
