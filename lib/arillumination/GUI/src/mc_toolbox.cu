#include <cuda.h>
#include <curand.h>
#include <curand_kernel.h>
#include <math.h>
#include "structArgsIn.h"


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
__global__ void simPhoton(
		float * heat_dev,
		const float * R_dev,
		const constArgsIn* inArgs,
		curandState *states){
	// precalculated to save loop time

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
		while(weight > 0){
			// move photon to next position
			move(x, y, z, r, u, v, w, weight, R_dev, states[tid], inArgs);
			absorb(heat_dev, r, z, weight, states[tid], inArgs);	
			scatter(u, v, w, states[tid], inArgs);
		}
	}
	return;
}
