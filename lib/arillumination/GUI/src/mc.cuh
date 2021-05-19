#ifndef MC_H
#define MC_H

#include "fiberProperties.h"
#include "mcFieldProperties.h"
#include "optProperties.h"
#include "simProperties.h"
#include <cuda.h>
#include <cuda_runtime.h>
#include <curand.h>
#include <curand_kernel.h>
#include <math.h>

#include "structArgsIn.h"
#include "../lib/vtkwriter/vtkwriter.h"

class mc
{
private:
	fiberProperties fiber;
	optProperties tissue;
	optProperties water;
	mcFieldProperties field;
	simProperties sim;
 
	float* heat; // normal heat map indexing: [ir + iz * nr]
	float* heat_log; // logarthmic scaled heat map [ir + iz * nr]
	bool isHeatAlloc = 0;
	float maxVal = 0; // maximum value in fluence map
	float minVal = 0; // minimum value in fluence map
	float maxValLog = 0;
	float minValLog = 0;

	float* heat_dev;
	bool isHeatDevAlloc = 0;

	// all the stuff we need for reflectance, should go in alonme class at some point
	uint64_t nR = 1000;
	float critAngle;
	float dAlpha;

	float* R;
	bool isRAlloc = 0;

	float* R_dev;
	bool isRDevAlloc = 0;

	void init_vars(); // initialize and allocate variables
	void calc_reflectance(const float nWater, const float nTissue); // calculates the tissues interal reflectance vector
	void run_sim(); // runs the actual simulation
	void calcMinMax();
	void calcLogHeat();

	float nWater = 1.33;


	bool flagDebug = 1; // give more debugging related output in terminal
public:

	// class constructor and class destructors
	mc();
	~mc();

	fiberProperties* get_pfiber() {return &fiber;};
	optProperties* get_ptissue() {return &tissue;};
	optProperties* get_pwater() {return &water;};
	mcFieldProperties* get_pfield() {return &field;};
	simProperties* get_psim() {return &sim;};	
	float* get_pheat() {return heat;};
	float* get_pheat_log() {return heat_log;};

	float get_maxVal() const {return maxVal;};
	float get_minVal() const {return minVal;};
	float get_maxValLog() const {return maxValLog;};
	float get_minValLog() const {return minValLog;};

	void run();
	void exportVtk(const string filePath);
};

#endif