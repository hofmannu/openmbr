
#include "vector3.h"
#include <cstdint>

#ifndef RECONSETT_H
#define RECONSETT_H

class reconsett
{
private:
	// vector3<int16_t> medfiltKernel = {1, 1, 1};
	float regStrength = 0.1f; // regularization strength
	vector3<float> dr = {10e-6f, 20e-6f, 20e-6f}; // resolution of output
	int16_t nIter = 10; // number of iteration
	// regularizer(1, :) char = 'Tikhonov'; % Tikhonov, L1, TV
	// cropZFlag(1, 1) = 0;

	// cropping of dataset in x, y, z
	vector3<float> startPos = {0.0f, 0.0f, -1e-3f};
	vector3<float> endPos = {10e-3f, 5e-3f, 1e-3f};
	float tCrop[2] = {0.0f, 20e-6f};

	float sos = 1495.0f; // speed of sound [m/s]

	// zRange(1, 2) single = [-1.25e-3, 1.55e-3];
	// xRange(1, 2) single; % if set to the same value: full range
	// yRange(1, 2) single; % if set to the same value: full range
	// xIdx(1, :) = []; 
	// yIdx(1, :) = [];
	// rRange(1, 1) single = 0.5e-3; % how much should we consider lateral range for model
	// flagFreqFiltModel(1, 1) logical = 1; % freq filter model same way as data
	// flagMedFiltModel(1, 1) logical = 1; % median filter model same way as data
	// flagFoundScanSett(1, 1) logical = 0; % did we find scanSett while loading?
	// flagFoundPreprocSett(1, 1) logical = 0; % did we find preprocSett while loading?
	// stepSize = 1; % reguilarization parameter for L1
	// gammaScale(1, 1) = 1;

public:
	void sortCropping();

	int16_t get_nIter() const {return nIter;};
	void set_nIter(const int16_t _nIter);

	float get_regStrength() const {return regStrength;};
	void set_regStrength(const float _regStrength);

	vector3<float> get_dr() const {return dr;};
	void set_dr(const vector3<float> _dr);

	vector3<float> get_startPos() const {return startPos;};
	void set_startPos(const vector3<float> _startPos);

	vector3<float> get_endPos() const {return endPos;};
	void set_endPos(const vector3<float> _endPos);


	float get_tMin() const {return tCrop[0];};
	float get_tMax() const {return tCrop[1];};

	float get_xMin() const {return startPos.x;};
	float get_xMax() const {return endPos.x;};
	float get_yMin() const {return startPos.y;};
	float get_yMax() const {return endPos.y;};

	float get_sos() const {return sos;};

	void set_tMin(const float tMin);
	void set_tMax(const float tMax);

	void set_xMin(const float xMin);
	void set_xMax(const float xMax);
	void set_yMin(const float yMin);
	void set_yMax(const float yMax);

	void set_cropT(const float tMin, const float tMax);
	void set_cropX(const float xMin, const float xMax);
	void set_cropY(const float yMin, const float yMax);

};

#endif