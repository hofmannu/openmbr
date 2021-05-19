/*
	c++ class handling the settings for our model based reconstruction algorithm
	Author: Urs Hofmann
	Mail: hofmannu@biomed.ee.ethz.ch
	Date: 22.08.2020
*/

#ifndef RECONSETTINGS_H
#define RECONSETTINGS_H

#include <iostream>
#include <cmath>

using namespace std;

class reconSettings
{

private:
	float SOS = 1490; // speed of sound in water [m/s]
	float threshold = 0.005; // threshold used to determine start and stop idx of mod matrix
	uint8_t nIter = 20; // number of LSQR rounds
	unsigned int downsampling[3] = {1, 1, 1};

	float boundRecon[6] = {6e-3, 7.5e-3, 1e-3, 18e-3, 1e-3, 18e-3};
	// z0, z1, x0, x1, y0, y1
	// zRecon, xRecon and yRecon defined as part of boundRecon
	
	float drRecon[3] = {20e-6, 25e-6, 25e-6}; // resolution of output volume [z, x, y]
	float rRes = 1e-6; // resolution in radial direction (model matrix) [m]

	// regularization parameters
	bool flagRegularization = 1; // turn regularization on or off
	float lambdaReg = 100; // regularization parameter
	string regMethod = "tk";

	bool flagH5Export = 1;
	bool flagVTKExport = 1;
	string polarityHandling = "pos";

public:
	string get_polarityHandling() const;
	string get_regMethod() const;
	bool isRegMethod(const string _regMethod) const;
	
	bool get_flagRegularization() const;
	bool* get_pflagRegularization() {return &flagRegularization;};

	bool get_flagH5Export() const;
	bool get_flagVTKExport() const;

	// Speed of sound
	void set_SOS(const float _SOS);
	float get_SOS() const;
	float* get_pSOS() {return &SOS;};

	// regularization parameter
	void set_lambdaReg(const float _lambdaReg);
	float get_lambdaReg() const;
	float* get_plambdaReg() {return &lambdaReg;};

	// threshold
	void set_threshold(const float _threshold);
	float get_threshold() const;
	float* get_pthreshold() {return &threshold;};

	void set_nIter(const uint8_t _nIter);
	uint8_t get_nIter() const;

	// rRes
	void set_rRes(const float _rRes);
	float get_rRes() const;
	float* get_prRes() {return &rRes;};

	float* get_dr();
	float get_dr(const unsigned int _dim) const;
	float get_dr0() const;
	float get_dr1() const;
	float get_dr2() const;
	void set_drRecon(const uint8_t _dim, const float dr);
	void set_drRecon(const float dr0, const float dr1, const float dr2);

	// reconstruction boundaries along z dimension
	float* get_zRecon();
	float get_zLower() const;
	float get_zUpper() const;
	void set_zUpper(const float _zUpper);
	void set_zLower(const float _zLower);
	void set_zLim(const float _zLower, const float _zUpper);

	// reconstruction boundaries along x dimension
	float* get_xRecon();
	float get_xLower() const;
	float get_xUpper() const;
	void set_xUpper(const float _xUpper);
	void set_xLower(const float _xLower);
	void set_xLim(const float _xLower, const float _xUpper);

	// reconstruction boundaries aling y dimension
	float* get_yRecon();
	float get_yLower() const;
	float get_yUpper() const;
	void set_yLower(const float _yLower);
	void set_yUpper(const float _yUpper);
	void set_yLim(const float _yLower, const float _yUpper);

	unsigned int* get_downsampling();
	unsigned int get_downsampling(const unsigned int _dim) const;

	uint64_t get_dim0() const;
	uint64_t get_dim1() const;
	uint64_t get_dim2() const;
	uint64_t get_dim(const uint8_t _dim) const;

	void print_properties();

};

#endif
