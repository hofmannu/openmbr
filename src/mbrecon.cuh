#ifndef MBRECON_H
#define MBRECON_H

#include <cstring>
#include <cmath>
#include <fstream>
#include "reconSettings.h"
#include "volume.h"
#include "transModel.cuh"
#include "fieldProperties.h"
#include "baseClass.h"
#include "basicMathOp.h"
#include "cudaErrorCheck.cuh"
#include "griddedData.cuh"
#include "structs.h"
#include "transducerProperties.h"
#include "../lib/vtkwriter/vtkwriter.h"
#include "gaussian_beam.h"
#include <time.h>

class mbrecon : public baseClass, public basicMathOp, public cudaErrorCheck
{
private:

	reconSettings sett; // reconstruction settings
	volume preprocVol; // preprocessed input volume
	volume reconVol; // reconstructed volume
	transModel transMod; // transducer modelling class
	transducerProperties transProp; // class holding transducer properties
	fieldProperties field;
	string inputPath = "/home/hofmannu/inputData.h5";
	modelMatrix* model; // pointer to model matrix (will be mapped to trans->model)
	griddedData croppedVol; // signal matrix stored as [it, ix, iy]
	griddedData estimAbs; // absorber map stored as [iz, ix, iy]
	
	// float* model; // transducer model
	float* Av; // estimated signal matrix [it, ix, iy] 
	float* p; // same size as signal matrix [it, ix, iy]
	
	uint64_t nr; // number of radial elements in model matrix
	float rMax = 0; // maximum radial distance to take into account [m]

	bool isRawDataLoaded = 0; // is raw data loaded?
	bool isDataCropped = 0;
	bool isRecon = 0; // is reconstruction performed
	bool flagDebug = 1; // enables additional output
	bool flagFluenceComp = 1; // should we include fluence compensation or not?

	// private functions
	void calc_nr(); // returns the number of elements we need in radial direction
	void build_model(); // build transdcuer model
	void getCroppedVol(); // crops input volume to region of interest
	void prepareOutputVolume();
	void forward(float* signalMatrix, const float* absorberMatrix); // forward matrix multiplication
	void transpose(float* absorberMatrix, const float* signalMatrix); // transpose matrix multiplication

public:
	mbrecon(); // class constructor
	~mbrecon(); // class destructor
	void load_data(); // load datasets from file with default path
	void load_data(const string _inputPath); // load dataset from file with specified path
	void recon(); // run iterative lsqr procedure
	void recon(uint8_t* iterator);
	void simulate(); 
	void export_data(); // export datasets
	void export_vtk();
	void print_recon_info(); // print information about reconstruction
	void set_transProp(const transducerProperties _transProp);
	void set_sett(const reconSettings _sett);
	bool get_isRecon() const;
	volume* get_ppreprocVol();
	volume* get_preconVol();
	transducerProperties* get_ptransProp();
	transModel* get_ptransModel();
	reconSettings* get_psett();
};

// structs for arguments passed to transpose kernel
struct constArgsInTrans{
	uint64_t nr; // number of elements in radial direction for model matrix
	griddedData* outputGrid; // definition of output grid (absorber map)
 	griddedData* inputGrid; // definition of input grid (signal matrix)
	float* sigMat; // signal matrix [it, ix, iy]
	uint64_t* modelStIdx; // model matrix start and stop index
	float* modelMat; // model matrix
	float dr; // resolution of model matrix in radial direction [m]
	float* surface; // surface representation [ix, iy] in index
	float* fluenceModel; // fluence model [iz, ix, iy]
};

struct constArgKTrans{
	int ir; // radial index of currently reconstructed line
	int xoff; // x offset in signal matrix
	int yoff; // y offset in signal matrix
	unsigned int xIm; // position of output matrix index in x
	unsigned int yIm; // position of output matrix index in y
};

// structs for arguments passed to forward kernel
struct constArgsInFwd{
	uint64_t nr; // number of elements in radial direction for model matrix
	griddedData* outputGrid; // definition of output grid (absorber map)
 	griddedData* inputGrid; // definition of input grid (signal matrix)
	float* estimatedAbs; // signal matrix [it, ix, iy]
	uint64_t* modelStIdx; // model matrix start and stop index
	float* modelMat; // model matrix
	float dr; // resolution of model matrix in radial direction [m]
	float* surface; // surface representation [ix, iy] in index
	float* fluenceModel; // fluence model [iz, ix, iy]
};

struct constArgKFwd{
	int ir; // radial index of currently reconstructed line
	int xoff; // x offset in signal matrix
	int yoff; // y offset in signal matrix
	uint64_t voloff; // volume offset
	uint64_t xIm; // position of output matrix index in x
	uint64_t yIm; // position of output matrix index in y
	uint64_t xAScan;
	uint64_t yAScan;
};

#endif
