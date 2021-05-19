#include "fieldProperties.h"

// constant get functions
float fieldProperties::get_fs() const {return fs;}
float fieldProperties::get_dt() const {return dt;}
float fieldProperties::get_dr() const {return dr;}
float fieldProperties::get_SOS() const {return SOS;}
unsigned int fieldProperties::get_nr() const {return nr;}
unsigned int fieldProperties::get_nz() const {return nz;}
unsigned int fieldProperties::get_nElements() const {return nElements;}
float fieldProperties::get_dSir() const {return dSir;}
float fieldProperties::get_sir0() const {return sir0;}
float fieldProperties::get_z0() const {return z0;}
float fieldProperties::get_dz() const {return dz;}

void fieldProperties::set_fs(const float _fs)
{
	if (_fs <= 0)
		throw "Sampling frequency must be bigger then 0";
	else
	{
		fs = _fs;
		dt = 1 / fs;
	}
	return;
}

void fieldProperties::set_dt(const float _dt)
{
	if (dt <= 0)
		throw "time sampling must be bigger then 0";
	else
	{
		dt = _dt;
		set_fs(1 / dt);
	}
	return;
}

void fieldProperties::set_dr(const float _dr)
{
	if (_dr <= 0)
	{
		dr = 1e-9; 
	}
	else
	{
		dr = _dr;
	}
	return;
}

void fieldProperties::set_dz(const float _dz)
{
	if (_dz <= 0)
		dz = 1e-9;
	else
		dz = _dz;

	return;
}


void fieldProperties::set_SOS(const float _SOS)
{
	if (_SOS <= 0)
		SOS = 1e-9;
	else
	{
		SOS = _SOS;
	}
	return;
}

void fieldProperties::set_nr(const unsigned int _nr)
{
	if (_nr <= 0)
	{
		nr = 1;
	}	
	else
	{
		nr = _nr;
		nElements = nr * nz;
	}
	return;
}


void fieldProperties::set_nz(const unsigned int _nz)
{
	if (_nz <= 0)
	{
		nz = 1;
	}
	else
	{
		nz = _nz;
		nElements = nr * nz;
	}
	return;
}

void fieldProperties::set_z0(const float _z0)
{
	z0 = _z0;
	return;
}

void fieldProperties::set_sir0(const float _sir0)
{
	sir0 = _sir0;
	return;
}

void fieldProperties::set_dSir(const float _dSir)
{
	dSir = _dSir;
	return;
}

void fieldProperties::print_properties()
{
	printf("*** Field properties ***\n");
	printf(" - Sampling frequency: %.1f MHz\n", fs * 1e-6);
	printf(" - dt: %.1f ns\n", dt * 1e9);
	printf(" - dr: %.1f microm\n", dr * 1e6);
	printf(" - dz: %.2f microm\n", dz * 1e6);
	printf(" - SOS: %.2f m/s\n", SOS);
	printf(" - nr: %d\n", nr);
	printf(" - nz: %d\n", nz);
	printf(" - z0: %f mm\n", z0 * 1e3);
	printf(" - dSir: %f microm\n", dSir * 1e6);
	printf(" - sir0: %f mm\n", sir0 * 1e3);
	return;
}