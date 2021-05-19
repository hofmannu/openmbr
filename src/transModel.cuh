/*
	c++ class respresenting an ultrasound transducer
	Author: Urs Hofmann
	Mail: hofmannu@biomed.ee.ethz.ch
	Date: 22.07.2020
*/

#ifndef TRANSMODEL_H
#define TRANSMODEL_H

#include <cmath>
#include <iostream>
#include <H5Cpp.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <time.h>
#include <cstdio>
#include <cstring>
#include <sstream>
#include <fstream> // used to export variables into files to lateron debug them
#include <thread>

// include our own libraries
#include "fieldProperties.h"
#include "structs.h"
#include "../lib/vtkwriter/vtkwriter.h"
#include "baseClass.h"
#include "modelMatrix.h"
#include "transducerProperties.h"
#include "gaussian_beam.h"

using namespace std;

class transModel: public baseClass
{

private:
	bool flagDebug = 1;

	// properties related to impulse response
	bool isImpRespDefined = 0; // flag if we loaded imp response
	float* impResp; // array holding the electrical impulse response
	float* tImp; // time vector of impulse response
	uint64_t nElem; // number of elements in time domain of impulse response
	float fSampling; // sampling frequency of impulse response [Hz]
	float dt; // incremental time step size (1 / fs) [s]
	float tOffset;

	// modelling parameters
	uint64_t nElements = 5e5; // number of discrete elements used for transducer modeling
	uint64_t nPoints = 0;

	fieldProperties fieldProp; 
	bool flagExportDiscret = 1; // wether transducer discretization should be exported
	string pathExportDiscret = "/home/hofmannu/transducerDiscretization.vtk";
	
	// transducer element description
	subElement* transducerElements; // transducer elements
	point_grid* pGrid;
	idx3* triangleIdR;
	bool isTriangleIdRAlloc = 0;
	bool isPGridAlloc = 0;
	bool isTransdElemAlloc = 0; // flag if allocation happened

	modelMatrix model; // model matrix (3d)
	float starttime; // starttime of first model matrix array

	void checkCudaReturn(const cudaError_t err, const string errMsg);
	void vprintf(const string txtMsg);

	bool flagVerbose = 1;
	bool isModelBuilt = 0;
	bool isRunning = 0;
	transducerProperties prop;
	gaussian_beam gauss;

	std::thread modelBuildingThread;

public:

	~transModel();

	// constant get functions
	float get_tOffset() const {return tOffset;};
	float get_fSampling() const {return fSampling;};
	uint64_t getNtModel() const;
	float getStarttime() const;
	float* get_pimpResp() {return impResp;};
	float* get_ptImp() {return tImp;};
	void resample(const float _fSampling);

	uint64_t getNElements() const {return nElements;};
	void setNElements(const uint64_t _nElements);
	
	uint64_t getNElem() const {return nElem;};
	uint64_t* get_pnElements() {return &nElem;};

	// set and get functions for field properties
	fieldProperties get_fieldProp() const {return fieldProp;};
	fieldProperties* get_pfieldProp() {return &fieldProp;};
	void setFieldProperties(const fieldProperties _fieldProp);

	// model matrix related function

	// impulse response related get functions
	float* getImpResp();

	void setTransProp(const transducerProperties _prop);

	void loadImpulseResponse();
	void loadImpulseResponse(const string filePath);

	// functions to build the transducer model on the GPU
	void buildGPUModel();
	float spherical_triangle_area(double *r1, double *r2, double *r3, const float R);
	void discretizeTransducer();
	modelMatrix* getModel();

	// opticla modelling parameters
	gaussian_beam* get_pgauss() {return &gauss;};

	// not implemented yet
	void exportElementsVtk(); // function to export discretized elements for paraview

};

#endif