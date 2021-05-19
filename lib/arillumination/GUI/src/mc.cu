#include "mc.cuh"

// scatters the photon into a new direction depending on currently employed 

// perfrom quaternion rotation
//   - r is normalized vector of rotation
//   - alpha is angle
//   - xyz is point to rotate, will be modified directly

// default structure of a quaternion
__device__ struct quaternion
{
	float x, y, z, w;
};

// get conjugate of it
__device__ static quaternion conjugate(quaternion quat)
{
  quat.x = -quat.x;
  quat.y = -quat.y;
  quat.z = -quat.z;
	// quat.w stays the same!
  return quat;
}

// multiply them
__device__ static quaternion mult(const quaternion& A, const quaternion& B)
{
  quaternion C;

  C.x = __fmul_rn(A.w, B.x) + __fmul_rn(A.x, B.w) + __fmul_rn(A.y, B.z) - __fmul_rn(A.z, B.y);
  C.y = __fmul_rn(A.w, B.y) - __fmul_rn(A.x, B.z) + __fmul_rn(A.y, B.w) + __fmul_rn(A.z, B.x);
  C.z = __fmul_rn(A.w, B.z) + __fmul_rn(A.x, B.y) - __fmul_rn(A.y, B.x) + __fmul_rn(A.z, B.w);
  C.w = __fmul_rn(A.w, B.w) - __fmul_rn(A.x, B.x) - __fmul_rn(A.y, B.y) - __fmul_rn(A.z, B.z);
  return C;
}

// perform a quaternion based rotation of our photon
__device__ static void quaternion_rot(
	const float& rx, 
	const float& ry, 
	const float& rz, 
	const float& alpha, 
	float& x, 
	float& y, 
	float& z){

	quaternion temp, quat_view, result;

	const float alpha2 = __fdividef(alpha, 2.0);
	const float sinalpha2 = __sinf(alpha2);

	// this vector is already by its definition 1 if rx, ry, rz has length 1 as well
	temp.x = __fmul_rn(rx, sinalpha2); 
	temp.y = __fmul_rn(ry, sinalpha2);
	temp.z = __fmul_rn(rz, sinalpha2);
	temp.w = __cosf(alpha2);

	// feed current position to first quaternion, already normalized as well
	quat_view.x = x;
	quat_view.y = y;
	quat_view.z = z;
	quat_view.w = 0.0;

	result = mult(mult(temp, quat_view), conjugate(temp));
	
	x = result.x; // update x output value
	y = result.y; // update y output value
	z = result.z; // update z output value
	
	return;
}

__device__ static float getsigma(
	const float& g, 
	const float& gx2, 
	const float& g2, 
	curandState& state){

	float sigma;

	// scatter photon
	if (g == 0){
		sigma = acosf(__fmul_rn(2.0, curand_uniform(&state)) - 1.0); 
		// cos(sigma) = 2 * RND - 1
	}else{
		// g - anisotroy coefficient, g2 = g^2, gx2 = 2 * g
		const float temp = __fdividef((1 - g2), (1 - g + __fmul_rn(gx2, curand_uniform(&state))));
		sigma = acosf(__fdividef((1 + g2 - __fmul_rn(temp, temp)), gx2));
	}

	return sigma;
}

__device__ static void scatter(
	float& u, // movement component in x
	float& v, // movement component in y
	float& w, // movement component in z
	curandState& state,
	const constArgsIn* inArgs){ // random number generator

	// generate rotation angles for scattering
	const float sigma = getsigma(inArgs->g, inArgs->gx2, inArgs->g2, state);
	const float phi = __fmul_rn(__fmul_rn(2.0, M_PI), curand_uniform(&state));

	float bx, by, bz;
	if ((abs(u) >= abs(v)) && (abs(u) >= abs(w))){
		// photon is mostly moving in x, so rectangular vector lies in y z plane
		by = abs(__fdividef(u, __fsqrt_rn(1 + __fmul_rn(u, u) + 
			__fmul_rn(2.0, __fmul_rn(v, w)) )));
		bz = by;
		bx = -__fmul_rn( __fdividef(by, u),  (v + w));
	}
	else if ((abs(v) >= abs(u)) && (abs(v) >= abs(w))){
		// photon is mostly moving in y, so rectangular vector lies in x y plane
		bx = abs(__fdividef(v, __fsqrt_rn(1 + __fmul_rn(v, v) + 
			__fmul_rn(2.0, __fmul_rn(u, w)) )));
		bz = bx;
		by = -__fmul_rn(__fdividef(bx, v), (u + w));
	}
	else if ((abs(w) >= abs(u)) && (abs(w) >= abs(v))){
		// photon is mostly moving in z, so rectangular vector lies in x y plane
		bx = abs(__fdividef(w, __fsqrt_rn(1 + __fmul_rn(w, w) + 
			__fmul_rn(2.0, __fmul_rn(u, v)) )));
		by = bx;
		bz = -__fmul_rn(__fdividef(bx, w), (u + v));
	}

	// rotate vector b around u to generate actual deflection vector
	quaternion_rot(u, v, w, phi, bx, by, bz);

	// rotate u aroung rotated b by sigma
	quaternion_rot(bx, by, bz, sigma, u, v, w); 
	return;
}

__device__ static void launch(
	float& weight, // weight of photon
	float& x, float& y, float& z, // position of photon in [x, y, z] in [m]
	float& u, float& v, float& w, // directory of movement in [x, y, z]
	curandState& state,
	const constArgsIn* inArgs){ // random number generator

	weight = 1.0; // reset weight of current photon to 1
	z = 0.0; // move z position to 0 [m]

	float r;

	// generate random x and y positions until in spot size
	do{
		x = (curand_uniform(&state) * 2.0 - 1) * inArgs->rSpot;
		y = (curand_uniform(&state) * 2.0 - 1) * inArgs->rSpot;
		r = __fsqrt_rn(__fmul_rn(x, x) + __fmul_rn(y, y));
	}while(r > inArgs->rSpot);

	// calculate direction of photon
	const float vLength = __fsqrt_rn(__fmul_rn(x, x) + __fmul_rn(y, y) + 
			__fmul_rn(inArgs->zSurface, inArgs->zSurface));
	u = __fdividef(x, vLength);
	v = __fdividef(y, vLength);
	w = __fdividef(inArgs->zSurface, vLength);
	return;
}

__device__ static void absorb(
	float* heat_dev, // matrix containing absorption
	const float& r, // radial position of photon
	const float& z, // axial position of photon
	float& weight,
	curandState& state,
	const constArgsIn* inArgs){

	// calculate index in radial direction
	uint64_t ir = __fdividef(r, inArgs->dr); // calculate index in radial direction
	if (ir >= inArgs->nr)
		ir = inArgs->nr - 1;
	uint64_t iz = __fdividef(z, inArgs->dz); // calculate index in axial direction
	if (iz >= inArgs->nz)
		iz = inArgs->nz - 1;

	// assign absorption event to heat map and decrease weight
	heat_dev[ir + inArgs->nr * iz] += __fmul_rn(inArgs->albedo1, weight);
	weight = __fmul_rn(weight, inArgs->albedo);

	// play roulette with photon
	if (weight < 1e-4){
		if (curand_uniform(&state) > 0.1)
			weight = 0;
		else
			weight = __fmul_rn(weight, 10);
	}
	return;
}

__device__ static void move(
	float& x, // most recent x position of photon
	float& y, // most recent y position of photon
	float& z, // most recent z position of photon
	float& r, // most recent radial position of photon
	const float& u, // movement component in x
	const float& v, // movement component in y
	float& w, // movement component in z
	float& weight, // weight of photon
	const float* R_dev, // reflectance matrix
	curandState& state,
	const constArgsIn* inArgs){ // random number generator

	// probability distribution for travel time
	const float t = -__fmul_rn(__logf(curand_uniform(&state)), inArgs->mu_as);
			
	x = __fadd_rn(x, __fmul_rn(t, u)); // move photon in x
	y = __fadd_rn(y, __fmul_rn(t, v)); // move photon in y
	z = __fadd_rn(z, __fmul_rn(t, w)); // move photon in z
	r = __fsqrt_rn(__fmul_rn(x, x) + __fmul_rn(y, y)); // recalculate radial pos

	// If we are leaving our tissue block
	if (z < 0){
		// move photon back to tissue surface
		// reduce photon intensity by reflectance index if there is no total internal reflection
		if (w < inArgs->critAngle){ // critical angle is already cosined
			const int idxR = __fdividef(acosf(w), inArgs->dAlpha);
			weight = __fmul_rn(weight, R_dev[idxR]);
		}
		z = -z; // invert position in z
		w = -w; // invert velocity in z
	}
	return;	
}

// cuda kernel definition
__global__ void simPhoton
(
		float * heat_dev,
		const float * R_dev,
		const constArgsIn* inArgs,
		curandState *states)
{
	float weight; // ranges from 0 to 1 representing the weight of ph

	// generate starting position
	float x, y, z; // position of photon [m]
	float r; // current radial position of photon and its inverse [m]
	float u, v, w; // directivity of photon propagation in tissue
	
	// temp variables for scattering
	const int tid = threadIdx.x + blockDim.x * blockIdx.x; // thread index
	curand_init(1234, tid, 0, &states[tid]); // Initialize CURAND

	#pragma unroll	
	for (uint64_t iPhoton = 0; iPhoton < inArgs->nPPerThread; iPhoton++){	
		// launch new photon
		launch(weight, x, y, z, u, v, w, states[tid], inArgs);
		
		// here starts the intense loop where calculation speed is critical	
		while(weight > 0)
		{
			// move photon to next position
			move(x, y, z, r, u, v, w, weight, R_dev, states[tid], inArgs);
			absorb(heat_dev, r, z, weight, states[tid], inArgs);	
			scatter(u, v, w, states[tid], inArgs);
		}
	}
	return;
}


// class destructor
mc::~mc()
{
	if (isHeatAlloc)
	{
		delete[] heat;
		delete[] heat_log;
	}

	if (isHeatDevAlloc)
		cudaFree(heat_dev);

	if (isRAlloc)
		delete[] R;

	if (isRDevAlloc)
		cudaFree(R_dev);
}

// class constructor
mc::mc()
{

}


void mc::init_vars()
{
	cudaError_t err;
	// alloc memory for heat map and log version of it on CPU
	if (isHeatAlloc)
	{
		delete[] heat;
		delete[] heat_log;
	}

	heat = new float [field.get_nElements()];
	heat_log = new float [field.get_nElements()];

	isHeatAlloc = 1;
	for (uint64_t idx = 0; idx < field.get_nElements(); idx++)
	{
		heat[idx] = 0;
		heat_log[idx] = 0;
	}

	// allocate matrix for heat map on GPU
	if (isHeatDevAlloc)
		cudaFree(heat_dev);

	// allocate memory on device for absorption map
	err = cudaMalloc( (void**)&heat_dev, field.get_nElements() * sizeof(float) );
	isHeatDevAlloc = 1;
	if (err != cudaSuccess)
	{
		printf("[mc] Could not allocate required memory on card\n");
		printf("[mc] Size of requested array: nr = %d, nz = %d, nElements = %d\n", 
			field.get_nr(), field.get_nz(), field.get_nElements());
		printf(cudaGetErrorString(err));
		printf("\n");
		throw "MemoryAllocError";
	}
	
	// Allocate memory on device for reflectance vector
	err = cudaMalloc( (void**)&R_dev, nR * sizeof(float) ); 
	isRDevAlloc = 1;

	if (err != cudaSuccess)
	{
		printf("[mc] Could not allocate required memory on card\n");
		printf("[mc] Size of requested array: nr = %d, nz = %d, nElements = %d\n", 
			field.get_nr(), field.get_nz(), field.get_nElements());
		printf(cudaGetErrorString(err));
		printf("\n");
		throw "MemoryAllocError";
	}

	bool cpy1 = (cudaSuccess != cudaMemcpy(heat_dev, heat, 
		field.get_nElements() * sizeof(float), cudaMemcpyHostToDevice));
	cpy1 |= (cudaSuccess != cudaMemcpy(R_dev, R, 
		nR * sizeof(float), cudaMemcpyHostToDevice));

	if (cpy1){
		printf("Could not copy memory to card\n");	
		throw "MemoryCopyError";
	}

	return;
}

void mc::calc_reflectance(const float nWater, const float nTissue)
{
	// calculate critical reflectance angle
	critAngle = asin(nWater / nTissue); // total ref when water to tissue
	dAlpha = critAngle / ((float) nR);

	if (isRAlloc)
		delete[] R;

	R = new float [nR];
	isRAlloc = 1;

	float alpha; // angle of incident beam
	float beta; // angle of outgoing beam
	float rs; // reflection coefficient for parallel pol
	float rp; // reflection coefficient for vertically pol
	float minR = 1000;
	float maxR = 0;
	
	for (uint64_t ir = 0; ir < nR; ir++){
	
		// calculate incoming angle for current step
		alpha = dAlpha * ir;
		
		// use snells law to calculate outgoing angle
		beta = asin(nTissue / nWater * sin(alpha));

		rp = tan(alpha - beta) / tan(alpha + beta);
		rs = sin(alpha - beta) / sin(alpha + beta); 	
		if (alpha == 0)
			R[ir] = pow(((nTissue - nWater) / (nTissue + nWater)), 2);
		else if (alpha > critAngle)
			R[ir] = 1; // total reflection
		else
			R[ir] = 0.5 * (pow(rp, 2) + pow(rs, 2));

		// printf("%2.2f %2.2f %2.2f\n", alpha, beta, *RVector_cp);
		// check for min and max value
		if (R[ir] > maxR)
			maxR = *R;

		if (R[ir] < minR)
			minR = R[ir];
	}

	return;
}

void mc::run_sim()
{
	cudaError_t err;

	// random number generator
	curandState *devStates;
	err = cudaMalloc( (void **)&devStates, 
		sim.get_threadsPerBlock() * sim.get_nBlocks() * sizeof(curandState) );
	if (err != cudaSuccess)
	{
		printf("Could not allocate memory for random number generators on GPU.\n");
		throw "cudaMallocErr";
	}

	// start stopwatch
	clock_t begin = clock(); 

	constArgsIn inArgs;
	inArgs.nr = field.get_nr(); // number of elements in radial direction
	inArgs.nz = field.get_nz(); // number of elements in axial direction 
	inArgs.dr = field.get_dr(); // resolution in radial direction [m]
	inArgs.dz = field.get_dz(); // resolution in axial direction [m]
	inArgs.zSurface = field.get_zBound(); // depth of boundary from fiber output on [m]
	inArgs.zSurface1 = 1 / inArgs.zSurface;
	inArgs.mu_a = tissue.get_mua(); // absorption coefficient [1/m]
	inArgs.mu_s = tissue.get_mus(); // scattering coefficient [1/m]
	inArgs.mu_as = 1 / (inArgs.mu_a + inArgs.mu_s);
	inArgs.albedo = tissue.get_albedo();
	inArgs.albedo1 = 1 - inArgs.albedo;
	inArgs.g = tissue.get_g();
	inArgs.g2 = inArgs.g * inArgs.g;
	inArgs.gx2 = 2 * inArgs.g;			

	inArgs.rSpot = fiber.get_rSpot(nWater, field.get_zBound()); // radius of light spot on surface

	// refractive index of both media
	inArgs.nWater = nWater; // refractive index of coupling medium
	inArgs.nTissue = tissue.get_n(); // refractive index of tissue
	
	// reflectance stuff
	inArgs.critAngle = critAngle;
	inArgs.nR = nR;
	inArgs.Rdev = R_dev;
	inArgs.dAlpha = dAlpha; 
	inArgs.nPPerThread = sim.get_nPPerThread(); // number of photons simulated per thread
 
 	// print out some important information about the simulation if debugging mode
	if (flagDebug)
	{
		printf("[debug] Field dimensions: nr = %d, nz = %d\n", inArgs.nr, inArgs.nz);
		printf("[debug] Number of blocks: %d, number of threads per block%d\n", 
			sim.get_nBlocks(), sim.get_threadsPerBlock());
		printf("[debug] Number of photons per thread: %d", inArgs.nPPerThread);
		printf("[debug] mua = %f, mus = %f, albedo = %f\n", 
			inArgs.mu_a, inArgs.mu_s, inArgs.albedo);
	}

	// allocate memory on GPU for constant arguments and copy them over
	constArgsIn* inArgs_dev;
	err = cudaMalloc( (void**)&inArgs_dev, sizeof(constArgsIn) );
	if (err != cudaSuccess)
	{
		printf("Could not allocate memory on card for settings struct.\n");
		throw "cudaMallocErr";
	}

	cudaMemcpy(inArgs_dev, &inArgs, sizeof(constArgsIn), cudaMemcpyHostToDevice);
	if (err != cudaSuccess)
	{
		printf("Could not copy settings struct over to GPU.\n");
		throw "cudaMemcpyErr";
	}

	// start actual simulation
	simPhoton<<<sim.get_nBlocks(), sim.get_threadsPerBlock() >>>(
		heat_dev, // output matrix into which we write our absorption
		R_dev, // 1d vector containing reflectance at different angles 
		inArgs_dev,
		devStates); // random number generators
	cudaDeviceSynchronize();	

	err = cudaGetLastError();
	if (err != cudaSuccess){
		printf(cudaGetErrorString(err));
	}else{	

		// copy data back from gpu
		err = cudaMemcpy(heat, heat_dev, 
			field.get_nElements() * sizeof(float), cudaMemcpyDeviceToHost);
		if (err != cudaSuccess){
			printf("Could not copy memory back from GPU\n");
			throw "MemoryCopyError";
		}
		else
		{
			printf("Compensating for dimension shift...\n");
			clock_t end = clock();
			// compensate for increase in intensity with increasing radius
			float areaOuter, areaInner, area, factorR;
			for (uint64_t ir = 0; ir < field.get_nr(); ir++){
				areaOuter = ((float) ir + 1.0) * ((float) ir + 1.0);
				areaInner = ((float) ir) * ((float) ir);
				area = areaOuter - areaInner; 
				factorR = 1 / area;
				for (uint64_t iz = 0; iz < field.get_nz(); iz++)
					heat[ir + iz * field.get_nr()] *= factorR; 
			}

			double time_spent = (end - begin) / (double) CLOCKS_PER_SEC;
			printf("Time elapsed in [s]: %f \n", time_spent);
			printf("Time per photon package in [ns]: %f \n", time_spent / (double) sim.get_nPhotonsTrue() * 1e9);
			printf("Photon packages per ms: %f \n", (double) sim.get_nPhotonsTrue() / (time_spent * 1000));
		}	
	}
	cudaFree(devStates);
	
	// free host and device memory
	cudaFree(heat_dev);
	isHeatDevAlloc = 0;

	cudaFree(R_dev);
	isRDevAlloc = 0;

	cudaFree(inArgs_dev);

	delete[] R;
	isRAlloc = 0;

	calcMinMax(); // calculate minimum and maximum value in heat map
	if (flagDebug)
	{
		printf("[debug] maximum value in field: %f, minimum value in fiedl: %f\n", maxVal, minVal);
	}

	calcLogHeat();

	return;
}

void mc::calcMinMax()
{
	maxVal = heat[0];
	minVal = heat[0];

	for (uint64_t iElement = 0; iElement < (field.get_nr() * field.get_nz()); iElement++)
	{

		if (heat[iElement] > maxVal)
		{
			maxVal = heat[iElement];
			
		}

		if (heat[iElement] < minVal)
		{
			minVal = heat[iElement];
			
		}
	}

	maxValLog = log10(maxVal);

	if (minVal > 0)
		minValLog = log10(minVal);
	else
		minValLog = -38;

	return;
}

void mc::calcLogHeat()
{
	for (uint64_t iElement = 0; iElement < field.get_nElements(); iElement++)
	{
		if (heat[iElement] > 0)
			heat_log[iElement] = log10(heat[iElement]);
		else
			heat_log[iElement] = -38; // lower end of float accuracy
	}
	return;
}

void mc::run()
{
	calc_reflectance(nWater, tissue.get_n());
	init_vars();
	run_sim();
	return;
}

// exports fluence field as a vtk file after rotating it around its own axis
void mc::exportVtk(const string filePath)
{
	vtkwriter myWriter;

	const string title ("fluence");
	myWriter.set_title(title);
	const string type ("STRUCTURED_POINTS");
	myWriter.set_type(type);
	myWriter.set_outputPath(filePath);

	uint64_t nRVol = ((float) field.get_nr()) / sqrt(2.0); // leave last container out due to trash bin
	uint64_t nRField = field.get_nr();
	uint64_t nZ = field.get_nz() - 1;
	uint64_t nXY = nRVol * 2 - 1;

	float * dataSet = new float [nXY * nXY * nZ];
	// dataset indexing: iz + ix * nz + iy * nz * nx

	// fill up 3d array
	uint64_t idxOutput;
	uint64_t idxR;
	float xPos;
	float yPos;
	float rPos;

	if (flagDebug)
		printf("[debug] Converting slice into 3D volume\n");
	
	#pragma unroll
	for (uint64_t ix = 0; ix < nXY; ix++)
	{
		xPos = (float) ix - (float) (nRVol - 1);
		#pragma unroll
		for (uint64_t iy = 0; iy < nXY; iy++)
		{
			yPos = (float) iy - (float) (nRVol - 1);
			rPos = sqrt(xPos * xPos + yPos * yPos);
			idxR = round(rPos);
			if (idxR < nRField)
			{
				#pragma unroll
				for (uint64_t iz = 0; iz < nZ; iz++)
				{
					idxOutput = iz + ix * nZ + iy * nZ * nXY; // order: [iz, ix, iy]
					dataSet[idxOutput] = heat_log[idxR + iz * nRField];
				}
			}
			else // set to minimum value if we are outside the range
			{
				for (uint64_t iz = 0; iz < nZ; iz++)
				{
					idxOutput = iz + ix * nZ + iy * nZ * nXY; // order: [iz, ix, iy]
					dataSet[idxOutput] = minValLog;
				}
			}
		}
	}

	// cut out a quater of the fluence field
	for (uint64_t ix = 0; ix < (nRVol - 1); ix++)
	{
		for (uint64_t iy = 0; iy < (nRVol - 1); iy++)
		{
			for (uint64_t iz = 0; iz < nZ; iz++)
			{
				idxOutput = iz + ix * nZ + iy * nZ * nXY;
				dataSet[idxOutput] = minValLog; // set to minimum value
			}
		}
	}

	griddedData myData;
	myData.data = dataSet;

	myData.origin[0] = field.get_zBound(); // origin in z
	myData.origin[1] = -((float) nRVol) * field.get_dr();// origin in x
	myData.origin[2] = -((float) nRVol) * field.get_dr();// origin in y
	
	myData.res[0] = field.get_dz(); // resolution in z direction
	myData.res[1] = field.get_dr(); // resolution in x direction
	myData.res[2] = field.get_dr(); // resolution in y direction

	myData.dim[0] = nZ;
	myData.dim[1] = nXY;
	myData.dim[2] = nXY;

	myWriter.set_structuredPoints(&myData);
	myWriter.set_binary();

	if (flagDebug)
		printf("[debug] Writing data to file\n");

	myWriter.write();

	delete[] dataSet;

	return;
}