
#ifndef CONSTARGSIN_H
#define CONSTARGSIN_H

struct constArgsIn
{
	// properties of output field
	uint64_t nr; // number of elements in field in lateral / radial direction
	uint64_t nz; // number of elements in field in axial direction
	float dr; // resolution of field in radial direction
	float dz; // resolution of field in axial direction
	float zSurface;
	float zSurface1;
	
	// optical properties of tissue
	float mu_a; // absorption coefficient of tissue [1/m]
	float mu_s; // scattering coefficient in tissue [1/m]
	float mu_as;
	float g; // anisotropy coefficient
	float g2;
	float gx2;
	float albedo; // albedo operator / value
	float albedo1;

	float rSpot; // spot radius on tissue surface

	// angular reflection stuff
	float critAngle;
	uint64_t nR;
	float* Rdev;
	float dAlpha;

	float nWater;
	float nTissue;

	uint64_t nPPerThread;


};

#endif