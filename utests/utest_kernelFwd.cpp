#include "../src/mbr_transsym.cuh"

int main()
{
	srand(1);

	const uint64_t nt = 301;
	const uint64_t nz = 201;
	const uint64_t nxModel = 11;
	const uint64_t nx = 201;
	const uint64_t ny = 101;

	mbr_transsym myKernel;

	// define the size of all matrixes
	myKernel.set_sizeAbsMat(nx, ny, nz);
	myKernel.set_sizeSigMat(nx, ny, nt);
	myKernel.set_sizeModel(nxModel, nxModel, nt, nz);
	
	// allocate memory for test datasets 
	float* sigMat = new float[myKernel.get_nElementsSig()];

	float* absMat = new float[myKernel.get_nElementsAbs()];
	for (uint64_t iElem = 0; iElem < myKernel.get_nElementsAbs(); iElem++)
	{
		absMat[iElem] = 1;
		absMat[iElem] = ((float) rand()) / ((float) RAND_MAX) - 0.5f;
	}

	float* modelMat = new float[myKernel.get_nElementsMod()];
	for (uint64_t iElem = 0; iElem < myKernel.get_nElementsMod(); iElem++)
		modelMat[iElem] = 1;
	
	// push pointers over to kernel
	myKernel.set_sigMat(sigMat);
	myKernel.set_modelMat(modelMat);
	myKernel.set_absMat(absMat);

	myKernel.print_info();
	myKernel.run_fwd();

	printf("Kernel execution took %f s\n", 
		((float) myKernel.get_tExec()) / 1000.0);

	uint64_t xMod, yMod;
	// calculate result for a single output position on the CPU and validate
	const uint64_t range = (nxModel - 1) / 2;
	const uint64_t idxOut[3] = {nx / 2, ny / 2, nt / 2}; // x, y, t
	float cpuVal = 0;
	for (uint64_t iz = 0; iz < nz; iz++)
	{
		yMod = nxModel - 1;
		for (uint64_t iy = 0; iy < nxModel; iy++)
		{
			const uint64_t yAbs = idxOut[1] - range + iy; 
			xMod = nxModel - 1;
			for (uint64_t ix = 0; ix < nxModel; ix++)
			{
				const uint64_t xAbs = idxOut[0] - range + ix;
				const float currAbs = absMat[xAbs + nx * (yAbs + ny * iz)]; 
				const float currModel = modelMat[xMod + nxModel * (yMod + nxModel * (idxOut[2] + nt * iz))];
				cpuVal = cpuVal + currAbs * currModel;
				xMod--;
			}
			yMod--;
		}
	}
	
	const uint64_t idxSigMat = idxOut[0] + nx * (idxOut[1] + ny * idxOut[2]); 
	printf("CPU calculation resulted in: %.4f\n", cpuVal);
	printf("GPU calculation resulted in: %.4f\n", sigMat[idxSigMat]);
	// printf("Value should be %d\n", (int) nxModel * nxModel * nz);
	const float error = abs(cpuVal - sigMat[idxSigMat]);
	if (error > 1e-6)
	{
		printf("There seems to be quite some difference between GPU and CPU\n");
		throw "InvalidResult";
	}

	// free host memory
	delete[] modelMat;
	delete[] sigMat;
	delete[] absMat;

	return 0;
}