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
int main(){

	// start stopwatch
	clock_t begin = clock(); 

	// overall number of photons
	const long nPhotons = 1e9; 
	
	// number of photons in each thread, if this number is too high, gpu cant finish
	const long nPPerThread = 5e3; 

	// calculate number of blocks and true number of photons	
	const long nBlocks = nPhotons / (THREADS_PER_BLOCK * nPPerThread);
	printf("True number of simulated photons: %i\n", nBlocks * THREADS_PER_BLOCK * nPPerThread);	

	// definition of scattering and absorption coefficient
	float mu_a = 200; // absorption coefficient in 1/m
	float mu_s = 1000;	// reduced scattering coefficient in 1/m

	// fiber and field properties
	float zSurface = 5e-3; // distance between surface and fiber output [m] 
	float zEnd = 11e-3; // distance between lower field end and fiber output [m]
	float rMax = 7e-3; // maximum simulated radial distance from acoutic axis [m]

	float rCore = 100e-6; // core radius of multimode fiber [m]
	float na = 0.2; // numerical aperture of fiber 
	float nWater = 1.33;
	float dr = 10e-6; // spatial resolution of output field in [m] 
	float nTissue = 1.41;
	float g = 0.7; // anisotropy coefficeint

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
	int nr = (rMax / dr) + 1; // number of bins in radial direction
	int nz = (zEnd - zSurface) / dr + 1; // number of bins in axial direction
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
		delete heat;
		delete R;
		return 0;
	}	

	bool cpy1 = (cudaSuccess != cudaMemcpy(heat_dev, heat, nr * nz * sizeof(float), cudaMemcpyHostToDevice));
	cpy1 |= (cudaSuccess != cudaMemcpy(R_dev, R, nR * sizeof(float), cudaMemcpyHostToDevice));
	if (cpy1){
		printf("Could not copy memory to card\n");	
		cudaFree(heat_dev);
		cudaFree(R_dev);
		delete heat;
		delete R;
		return 0;
	}

	// generate a random number generator for each worker
	curandState *devStates;
	cudaMalloc( (void **)&devStates, THREADS_PER_BLOCK * nBlocks * sizeof(curandState) );

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
		dr, // resolution in both axial and radial direction [m]
		nr, // number of field steps in radial direction
		nz, // number of field steps in axial direction
		albedo, // albedo operator
		nWater, // refractive index of water
		nTissue, // refractive index of tissue
		nPPerThread, // number of photons simulated in each thread
		devStates); // random number generators
	
	cudaError_t err = cudaGetLastError();
	if (err != cudaSuccess){
		printf(cudaGetErrorString(err));
		cudaFree(heat_dev);
		cudaFree(R_dev);
		delete heat;
		delete R;
		return 0;
	}	

	// copy data back from gpu
	err = cudaMemcpy(heat, heat_dev, nr * nz * sizeof(float), cudaMemcpyDeviceToHost);
	if (err != cudaSuccess){
		printf("Could not copy memory back from GPU: ");
		printf(cudaGetErrorString(err));
		cudaFree(heat_dev);
		cudaFree(R_dev);
		delete R;
		delete heat;
		return 0;
	}

	// compensate for increase in intensity with increasing radius
	for (int ir = 0; ir < nr; ir++){
		float factorR = 1 / ((float) ir + 0.5);
		for (int iz = 0; iz < nz; iz++){
			heat[ir + iz * nr] = heat[ir + iz * nr] * factorR; 
		}
	}

	printf("Writing to output file...\n");
	ofstream outfile("outputFile.bin", ofstream::out | ofstream::binary);
	long nByte = sizeof(float) * nr * nz;
	outfile.write(reinterpret_cast<char * > (heat), nByte);
	outfile.close();

	cudaFree(heat_dev);
	delete heat;
	cudaFree(R_dev);
	delete R;	

	clock_t end = clock();
	double time_spent = (end - begin) / (double) CLOCKS_PER_SEC;
	printf("Time elapsed in [s]: %f \n", time_spent);
	printf("Time per photon package in [ns]: %f \n", time_spent / (double) nPhotons * 1e9);

	return 0;
}

