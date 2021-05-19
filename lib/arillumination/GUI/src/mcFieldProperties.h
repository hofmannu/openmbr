#ifndef MCFIELDPROPERTIES_H
#define MCFIELDPROPERTIES_H

#include <cstdio>
#include <cstdint>

class mcFieldProperties
{
private:
	float rMax = 3e-3; // maximum simulated radial distance
	float zBound = 3e-3; // minimum distance between fiber output and tissue
	float zMax = 7e-3; // maximum simulated distance between fiber output and tissue
	float dz = 5e-6; // resolution in z direction
	float dr = 5e-6; // resolution in radial direction

	// dependent variables
	uint64_t nr;
	uint64_t nz;

	// internal calculation functions
	void calc_nr();
	void calc_nz();
public:
	// class constructors and destructors
	mcFieldProperties();

	// constant get functions of class
	float get_rMax() const;
	float get_zBound() const;
	float get_zMax() const;
	float get_dz() const;
	float get_dr() const;
	float get_zExtend() const {return (zMax - zBound);};
	uint64_t get_nr() const;
	uint64_t get_nz() const;
	uint64_t get_nElements() const;

	void set_rMax(const float _rMax);
	void set_zBound(const float _zBound);
	void set_zMax(const float _zMax);
	void set_dz(const float _dz);
	void set_dr(const float _dr);

};

#endif