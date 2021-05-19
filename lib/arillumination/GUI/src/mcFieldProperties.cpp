#include "mcFieldProperties.h"

mcFieldProperties::mcFieldProperties()
{
	// initialize default values
	calc_nr();
	calc_nz();
}

float mcFieldProperties::get_rMax() const {return rMax;}
float mcFieldProperties::get_zBound() const {return zBound;}
float mcFieldProperties::get_zMax() const {return zMax;}
float mcFieldProperties::get_dz() const {return dz;}
float mcFieldProperties::get_dr() const {return dr;}
uint64_t mcFieldProperties::get_nr() const {return nr;}
uint64_t mcFieldProperties::get_nz() const {return nz;}
uint64_t mcFieldProperties::get_nElements() const 
{
	return get_nr() * get_nz();
}

void mcFieldProperties::calc_nr()
{
	nr = rMax / dr + 1;
	return;
}

void mcFieldProperties::calc_nz()
{
	nz = (zMax - zBound) / dz + 1;
	return;
}

void mcFieldProperties::set_rMax(const float _rMax)
{
	rMax = _rMax;
	calc_nr();
	return;
}

void mcFieldProperties::set_zBound(const float _zBound)
{
	zBound = _zBound;
	calc_nz();
	return;
}

void mcFieldProperties::set_zMax(const float _zMax)
{
	zMax = _zMax;
	calc_nz();
	return;
}

void mcFieldProperties::set_dz(const float _dz)
{
	dz = _dz;
	calc_nz();
	return;
}

void mcFieldProperties::set_dr(const float _dr)
{
	dr = _dr;
	calc_nr();
	return;
}