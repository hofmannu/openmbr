#include <mex.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fstream>
#include <time.h>

#define THREADS_PER_BLOCK 1024

#include "calc_r.cpp"
#include "mc_toolbox.cpp"

using namespace std;

// main procedure
// nlhs: number of left hand side arguments (out)
// plhs: pointer to left hand side arguments
// same for right side (input)
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){

	// vector tissue conatins information like scattering and absorption coefficients
	float* tissue_cp = (float *) mxGetData(prhs[0]);
	float mu_a = tissue_cp[0]; // absorption coefficient in 1/m
	float mu_s = tissue_cp[1]; // scattering coefficient in 1/m
	float g = tissue_cp[2]; // anisotropy of tissue
	float nTissue = tissue_cp[3]; // optical index of tissue

	// vector fiber contains important fiber properties
	float* fiber_cp = (float *) mxGetData(prhs[1]);	
	float rCore = fiber_cp[0]; // core radius of multimode fiber [m]
	float na = fiber_cp[1]; // numerical aperture of fiber 

	float* field_cp = (float *) mxGetData(prhs[2]);
	float rMax = field_cp[0]; // maximum simulated radial distance from acoutic axis [m]
	float zEnd = field_cp[1]; // distance between lower field end and fiber output [m]
	float zSurface = field_cp[2]; // distance between surface and fiber output [m] 
	float drr = field_cp[3]; // spatial resolution of output field in radial direction [m] 
	float drz = field_cp[4]; // spaital resolution of output field in axial direction [m]

	// vector simulation determines important simulation properties
	float* sim_cp = (float *) mxGetData(prhs[3]);
	const long nPPerThread = sim_cp[0]; // number of photons per thread 
	const long nPhotons = sim_cp[1]; // overall number of simulated photons 
	
	// calculate number of blocks and true number of photons	
	const long nBlocks = nPhotons / (THREADS_PER_BLOCK * nPPerThread);
	printf("True number of simulated photons: %lu\n", 
		nBlocks * THREADS_PER_BLOCK * nPPerThread);	

	float nWater = 1.33;

	float theta = asin(na / nWater); // opening angle of beam (one sided)
	float rSpot = rCore + tan(theta) * zSurface; // spot diameter on surface

	float crit_angle = asin(nWater / nTissue); // total ref when water to tissue

	// calculate reflectance index matrix
	unsigned int nR = 1000; // number of substeps simulated for spectral reflectance
	float dalpha = crit_angle / ((float) nR); // stepsize as function of critical angle
	float* R = new float[nR]; // allocate memory for reflectance vector from tissue to water
	calc_r(R, nR, dalpha, nWater, nTissue); // call function

	crit_angle = cos(crit_angle); // crit_angle needs to be cosine for simulation

	// generate field and set everything to 0
	int nr = (rMax / drr) + 1; // number of bins in radial direction
	int nz = (zEnd - zSurface) / drz + 1; // number of bins in axial direction
	float * heat = new float[nr * nz]; // output matrix
	for (int index = 0; index < (nr * nz); index++) // set output matrix to 0
		heat[index] = 0;
	printf("Field size: nr = %d, nz = %d \n", nr, nz);

	float albedo = mu_s / (mu_s + mu_a);

	// allocate memory on device for absorption map
	float* heat_dev;
	bool m1 = (cudaSuccess != cudaMalloc( (void**)&heat_dev, nr * nz * sizeof(float) ));
	
	// Allocate memory on device for reflectance vector
	float* R_dev;
	m1 |= (cudaSuccess != cudaMalloc( (void**)&R_dev, nR * sizeof(float) )); 
	
	if (m1){
		printf("Could not allocate memory on card\n");
		delete[] heat;
		delete[] R;
		return;
	}	

	bool cpy1 = (cudaSuccess != cudaMemcpy(heat_dev, heat, nr * nz * sizeof(float), cudaMemcpyHostToDevice));
	cpy1 |= (cudaSuccess != cudaMemcpy(R_dev, R, nR * sizeof(float), cudaMemcpyHostToDevice));
	if (cpy1){
		printf("Could not copy memory to card\n");	
		return;
	}else{

		// generate a random number generator for each worker
		curandState *devStates;
		cudaMalloc( (void **)&devStates, THREADS_PER_BLOCK * nBlocks * sizeof(curandState) );

		// start stopwatch
		clock_t begin = clock(); 
		printf("Starting simulation...\n");
		
		// run on cuda!!!
		simPhoton<<<nBlocks, THREADS_PER_BLOCK >>>(
			heat_dev, // output matrix into which we write our absorption
			R_dev, // 1d vector containing reflectance at different angles 
			dalpha, // incremental steps for angular reflection calculation [rad]
			rSpot, // size of spot at surface [m]
			zSurface, // distance between surface and fiber [m]
			mu_a, // absorption coefficient
			mu_s, // scattering coefficient
			g, // anisotropy coefficient
			crit_angle, // critical angle for backreflection
			drr, // resolution in radial direction [m]
			drz, // resolution of output field in axial direction [m]
			nr, // number of field steps in radial direction
			nz, // number of field steps in axial direction
			albedo, // albedo operator
			nWater, // refractive index of water
			nTissue, // refractive index of tissue
			nPPerThread, // number of photons simulated in each thread
			devStates); // random number generators

		cudaDeviceSynchronize();	
		
		cudaError_t err = cudaGetLastError();
		if (err != cudaSuccess){
			printf(cudaGetErrorString(err));
		}else{	

			// copy data back from gpu
			err = cudaMemcpy(heat, heat_dev, 
				nr * nz * sizeof(float), cudaMemcpyDeviceToHost);
			if (err != cudaSuccess){
				printf("Could not copy memory back from GPU: ");
				printf(cudaGetErrorString(err));
			}else{

			printf("Compensating for dimension shift...\n");
			clock_t end = clock();
			// compensate for increase in intensity with increasing radius
			float areaOuter, areaInner, area, factorR;
			for (int ir = 0; ir < nr; ir++){
				areaOuter = ((float) ir + 1.0) * ((float) ir + 1.0);
				areaInner = ((float) ir) * ((float) ir);
				area = areaOuter - areaInner; 
				factorR = 1 / area;
				for (int iz = 0; iz < nz; iz++)
					heat[ir + iz * nr] *= factorR; 
			}

			double time_spent = (end - begin) / (double) CLOCKS_PER_SEC;
			printf("Time elapsed in [s]: %f \n", time_spent);
			printf("Time per photon package in [ns]: %f \n", time_spent / (double) nPhotons * 1e9);
			printf("Photon packages per ms: %f \n", (double) nPhotons / (time_spent * 1000));

			plhs[0] = mxCreateNumericMatrix(nr, nz, mxSINGLE_CLASS, mxREAL);
			memcpy(mxGetPr(plhs[0]), heat, nr * nz * sizeof(float));
			}	
		}
		cudaFree(devStates);
	}
	// free host and device memory
	cudaFree(heat_dev);
	delete[] heat;
	cudaFree(R_dev);
	delete[] R;	

	return;	

}

