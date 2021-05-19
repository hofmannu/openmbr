/*
	File: fielProperties.h
	Author: Urs Hofmann
	Mail: hofmannu@ethz.ch
	Date: 27.07.2020
*/

#ifndef FIELDPROPERTIES_H
#define FIELDPROPERTIES_H
// structure holding properties of our output grid in spatial dimensions

#include <iostream>
#include <cstdio>

class fieldProperties{
	private:
		float fs = 250e6; // sampling frequency [Hz]
		float dt = 1 / 250e6; // time step [s]
		float dr = 1e-6; // radial spacing of field [m]
		float dz = 10e-6; // axial spacing of field [m]
		float SOS = 1490; // assumed speed of sound

		unsigned int nr = 401; // number of field elements in radial direction
		unsigned int nz = 401; // number of field elements in axial direction
		unsigned int nElements = 401 * 401;
		float z0 = -2e-3; // z position relative to focus of element iz = 0
		float dSir = 1490 / 250e6 / 4; // set resolution of spatial impulse response to small value
		float sir0; // 0 element of spatial impulse response
		//float sir0;
	public:

		void set_fs(const float _fs);
		float get_fs() const;
		float* get_pfs() {return &fs;};
		
		void set_dt(const float _dt);
		float get_dt() const;

		void set_dr(const float _dr);
		float get_dr() const;
		float* get_pdr() {return &dr;};
		
		void set_dz(const float _dz);
		float get_dz() const;
		
		void set_SOS(const float SOS);
		float get_SOS() const;

		void set_nr(const unsigned int _nr);
		unsigned int get_nr() const;

		void set_nz(const unsigned int _nz);
		unsigned int get_nz() const;

		unsigned int get_nElements() const;

		void set_z0(const float _z0);
		float get_z0() const;
		
		void set_dSir(const float _dSir);
		float get_dSir() const;

		void set_sir0(const float _sir0);
		float get_sir0() const;

		void print_properties();

};
#endif