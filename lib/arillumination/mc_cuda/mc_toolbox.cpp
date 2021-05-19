#include <cuda.h>
#include <curand.h>
#include <curand_kernel.h>
#include <math.h>

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
	const float& g, // anisotropy factor g
	const float& g2, // pow(g, 2)
	const float& gx2, // g * 2
	curandState& state){ // random number generator

	// generate rotation angles for scattering
	const float sigma = getsigma(g, gx2, g2, state);
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
	const float& zSurface, // surface distance
	const float& rSpot, // radius of spot
	curandState& state){ // random number generator

	weight = 1.0; // reset weight of current photon to 1
	z = 0.0; // move z position to 0 [m]

	float r;

	// generate random x and y positions until in spot size
	do{
		x = (curand_uniform(&state) * 2.0 - 1) * rSpot;
		y = (curand_uniform(&state) * 2.0 - 1) * rSpot;
		r = __fsqrt_rn(__fmul_rn(x, x) + __fmul_rn(y, y));
	}while(r > rSpot);

	// calculate direction of photon
	const float vLength = __fsqrt_rn(__fmul_rn(x, x) + __fmul_rn(y, y) + 
			__fmul_rn(zSurface, zSurface));
	u = __fdividef(x, vLength);
	v = __fdividef(y, vLength);
	w = __fdividef(zSurface, vLength);
	return;
}

__device__ static void absorb(
	float* heat_dev, // matrix containing absorption
	const float& r, // radial position of photon
	const float& z, // axial position of photon
	const float& drr, const float& drz, // resolution in radial and axial dir
	float& weight,
	const float& albedo,
	const float& albedo1,
	const int& nr, // number of elements in radial direction
	const int& nz, // number of elements in axial direction
	curandState& state){

	// calculate index in radial direction
	unsigned int ir = __fdividef(r, drr); // calculate index in radial direction
	if (ir >= nr)
		ir = nr - 1;
	unsigned int iz = __fdividef(z, drz); // calculate index in axial direction
	if (iz >= nz)
		iz = nz - 1;

	// assign absorption event to heat map and decrease weight
	heat_dev[ir + nr * iz] += __fmul_rn(albedo1, weight);
	weight = __fmul_rn(weight, albedo);

	// play roulette with photon
	if (weight < 1e-4){
		if (curand_uniform(&state) > 0.1)
			weight = 0;
		else
			weight = __fmul_rn(weight, 10);
	}
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
	const float& dalpha, // step size of angular reflectance matrix
	const float& mu_as, // combined absorption and scattering coefficient 
	const float& crit_angle, // total reflection angle cosined
	curandState& state){ // random number generator

	// probability distribution for travel time
	const float t = -__fmul_rn(__logf(curand_uniform(&state)), mu_as);
			
	x = __fadd_rn(x, __fmul_rn(t, u)); // move photon in x
	y = __fadd_rn(y, __fmul_rn(t, v)); // move photon in y
	z = __fadd_rn(z, __fmul_rn(t, w)); // move photon in z
	r = __fsqrt_rn(__fmul_rn(x, x) + __fmul_rn(y, y)); // recalculate radial pos

	// If we are leaving our tissue block
	if (z < 0){
		// move photon back to tissue surface
		// reduce photon intensity by reflectance index if there is no total internal reflection
		if (w < crit_angle){ // critical angle is already cosined
			const int idxR = __fdividef(acosf(w), dalpha);
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
		const float dalpha,
		const float rSpot, // radius of light spot at surface [m]
		const float zSurface, // distance btw fiber output and surface [m]
		const float mu_a, // absorption coefficient [1/m]
		const float mu_s, // scattering coefficient [1/m]
		const float g, // anisotropy coefficient
		const float crit_angle, // critical reflection angle cosined
		const float drr, // resolution in radial and axial direction [m]
		const float drz,
		const int nr, // number of elements in radial direction
		const int nz, // number of elements in axial direction
		const float albedo,
		const float nWater, // refractive index of water
		const float nTissue, // refractive index of tissue
		const int nPPerThread, // number of photons simulated in each thread
		curandState *states){

	const float albedo1 = 1 - albedo; // precalculated to save loop time
	const float mu_as = __fdividef(1, __fadd_rn(mu_a, mu_s)); // precalculated to save loop time

	float weight; // ranges from 0 to 1 representing the weight of ph

	// generate starting position
	float x, y, z; // position of photon [m]
	float r; // current radial position of photon and its inverse [m]
	float u, v, w; // directivity of photon propagation in tissue
	
	// temp variables for scattering
	const float g2 = __fmul_rn(g, g); // g^2
	const float gx2 = __fmul_rn(2.0, g);
	const float zSurface1 = __fdividef(1.0, zSurface); // inverse of surface distance
	const int tid = threadIdx.x + blockDim.x * blockIdx.x; // thread index
	curand_init(1234, tid, 0, &states[tid]); // Initialize CURAND

	#pragma unroll	
	for (int iPhoton = 0; iPhoton < nPPerThread; iPhoton++){	
		// launch new photon
		launch(weight, x, y, z, u, v, w, zSurface, rSpot, states[tid]);
		
		// here starts the intense loop where calculation speed is critical	
		while(weight > 0){
			// move photon to next position
			move(x, y, z, r, u, v, w, weight, R_dev, dalpha, mu_as, crit_angle, states[tid]);
			absorb(heat_dev, r, z, drr, drz, weight, albedo, albedo1, nr, nz, states[tid]);	
			scatter(u, v, w, g, g2, gx2, states[tid]);
		}
	}
	return;
}
