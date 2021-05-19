/*
	small helper class to build a gaussian beam model for optics of OR
	Author: Urs Hofmann
	Mail: hofmannu@ethz.ch
	Date: 19.05.2021
*/

#ifndef GAUSSIAN_BEAM_H
#define GAUSSIAN_BEAM_H

#include <cstdint>
#include <cmath>
#include <iostream>

using namespace std;

class gaussian_beam
{
private:
	// optical properties
	float na = 0.02; // numerical aperture of focusing lens
	float wavelength = 532; // wavelength of used light [nm]
	float n = 1.33; // refractive index of medium
	float I0 = 1;

	// field properties
	uint64_t nr;
	uint64_t nz;
	float z0; // starting point of field along z [m]
	float dr; // resolution in radial direction [m]
	float dz; // resolution in axial direction [m]

	bool flagVerbose = 1;

public:
	// set functions
	void define_r(const float _dr, const uint64_t _nr);
	void define_z(const float _dz, const float _z0, const uint64_t nz);
	
	void set_wavelength(const float _wavelength);
	void set_n(const float _n);
	void set_na(const float _na);
	void set_I0(const float _I0);

	float get_wavelength() const {return wavelength;};
	float get_n() const {return n;};
	float get_na() const {return na;};
	float get_I0() const {return I0;};

	float* get_pwavelength() {return &wavelength;};
	float* get_pna() {return &na;};
	float* get_pn() {return &n;};
	float* get_pI0() {return &I0;};

	void convolve_model(float* modelMatrix, const uint64_t nt);

};

#endif