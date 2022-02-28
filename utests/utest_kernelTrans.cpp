#include "../src/mbr_transsym.cuh"
#include <cmath>

int main()
{
	srand(1);

	const uint64_t nt = 111;
	const uint64_t nz = 101;
	const uint64_t nxModel = 31;
	const uint64_t nx = 501;
	const uint64_t ny = 201;

	mbr_transsym myKernel;

	myKernel.set_sizeSigMat(nx, ny, nt);
	myKernel.set_sizeModel(nxModel, nxModel, nt, nz);
	
	// allocate memory for test datasets 
	float* sigMat = new float[myKernel.get_nElementsSig()];
	for (uint64_t iElem = 0; iElem < myKernel.get_nElementsSig(); iElem++)
		sigMat[iElem] = ((float)rand()) / ((float) RAND_MAX) - 0.5f;

	float* absMat = new float[myKernel.get_nElementsAbs()];
	float* modelMat = new float[myKernel.get_nElementsMod()];
	for (uint64_t iElem = 0; iElem < myKernel.get_nElementsMod(); iElem++)
		modelMat[iElem] = 1;
	
	myKernel.set_sigMat(sigMat);
	myKernel.set_modelMat(modelMat);
	myKernel.set_absMat(absMat);

	myKernel.print_info();
	myKernel.run_trans();

	printf("Kernel execution took %f s\n", ((float) myKernel.get_tExec()) / 1000.0);


	// calculate result for a single output position on the CPU and validate
	const uint64_t range = (nxModel - 1) / 2;
	const uint64_t idxOut[3] = {nx / 2, ny / 2, nz / 2}; // x, y, z
	float cpuVal = 0;
	for (uint64_t it = 0; it < nt; it++)
	{
		for (uint64_t iy = 0; iy < nxModel; iy++)
		{
			const uint64_t yAbs = idxOut[1] - range + iy; 
			for (uint64_t ix = 0; ix < nxModel; ix++)
			{
				const uint64_t xAbs = idxOut[0] - range + ix;
				const float currSig = sigMat[xAbs + nx * (yAbs + ny * it)]; 
				const float currModel = modelMat[ix + nxModel * (iy + nxModel * (it + nt * idxOut[2]))];
				cpuVal = cpuVal + currSig * currModel;
			}
		}
	}
	
	const uint64_t idxAbsMat = idxOut[0] + nx * (idxOut[1] + ny * idxOut[2]); 
	printf("CPU calculation resulted in: %.4f\n", cpuVal);
	printf("GPU calculation resulted in: %.4f\n", absMat[idxAbsMat]);
	const float error = abs(cpuVal - absMat[idxAbsMat]);
	if (error > 1e-6)
	{
		printf("Some discrepancies are inacceptable\n");
		throw "InvalidResult";
	}
	
	// free host memory
	delete[] modelMat;
	delete[] sigMat;
	delete[] absMat;

	return 0;
}