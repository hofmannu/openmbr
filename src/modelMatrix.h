/*
	c++ class representing the model matrix
	Author: Urs Hofmann
	Mail: hofmannu@ethz.ch
	Date: 16.08.2020
*/

#ifndef MODELMATRIX_H
#define MODELMATRIX_H

#include "basicMathOp.h"
#include "baseClass.h"
#include <cstdio>
#include "../lib/vtkwriter/vtkwriter.h"
#include "griddedData.cuh"
#include <iostream>
#include <H5Cpp.h>

using namespace std;

class modelMatrix : public basicMathOp, public baseClass
{
private:
	uint64_t nElements;	// number of elements in model matrix
	uint64_t nt; // number of elements in time domain for model matrix
	uint64_t nz; // number of elements in z domain for model and sir
	uint64_t nr; // number of elements in radial domain for model and sir
	
	uint64_t nElementsSir; // number of elements in spatial impulse response
	uint64_t nzSir; // number of elements in SIR along 1st dim
	
	float maxVal; // maximum value in model matrix
	float maxValSir; // maximum value in spatial impulse response
	
	bool isModelMatrixAlloc = 0; // is model matrix matrix allocated
	bool isSirAlloc = 0; // is spatial impulse response matrix allocated
	bool isStartStopIdxAlloc = 0;
	bool isStartStopIdxPermutedAlloc = 0;
	float dr; // resolution of model matrix in radial direction [m]
	float dz; // resolution of model matrix in axial direction [m]
	float z0; // zero point of model matrix relative to focal plane [m]
	float tRes; // temporal resolution of model matrix
	float t0; // time offset of model matrix
	float sir0; // offset of spatial impulse response matrix along 1st dim
	float sirRes; // resolution of spatial impulse response along 1st dim

	float* sensField; // sensitivity field as max abs along t
	bool isSensFieldAlloc = 0;

	float* modSlice; // needs to have order ir, it
	bool isModSliceAlloc = 0;
	uint64_t lastSlice = 0;
	float maxValSlice = 1.0; // maximum value in slice
	float minValSlice = -1.0; // minimum value in slice

	uint64_t modVecR = 0; // last radial idx of model vector
	uint64_t modVecZ = 0; // last axial idx of model vector
	
	float modVecAbsMax = 1;
public:
	float* data; // idx order standard [it, ir, iz] transpose calc
	float* dataPermuted; // idx order permuted [iz, ir, it] forward calc
	uint64_t* startStopIdx; // indexing order 
	uint64_t* startStopIdxPermuted; // indexing order (forward calc)

	float* sir; // stored spatial impulse response [isir, ir, iz]

	float threshold = 0.01; // defines start and stop index as percentage of max

	// constructor / destructor
	modelMatrix();
	~modelMatrix();

	// constant get functions
	uint64_t getNElements() const {return nElements;};
	uint64_t getNElementsSir() const {return nElementsSir;};
	uint64_t get_nt() const {return nt;};
	uint64_t get_nz() const {return nz;};
	uint64_t get_nr() const {return nr;};
	uint64_t getNzSir() const {return nzSir;};
	float get_rDist() const {return (float) nr * dr;};
	float get_zDist() const {return (float) nz * dz;};
	float get_zMin() const {return z0;};
	float get_zMax() const {return z0 + (float) nz * dz;};
	float get_rMax() const {return (float) nr * dr;};
	float get_t0() const {return t0;};
	float get_maxValSlice() const {return maxValSlice;};
	float get_minValSlice() const {return minValSlice;};
	float* get_psensField() {return sensField;};
	float getMaxVal() const {return maxVal;};
	float get_dr() const {return dr;};
	float get_dz() const {return dz;}; 
	float get_z0() const {return z0;};
	float* get_data() {return data;};

	// set functions
	void setNtModel(uint64_t _ntModel);
	void setNzModel(uint64_t _nzModel);
	void setNrModel(uint64_t _nrModel);
	void setNzSir(uint64_t _nzSir);
	void setRRes(const float _rRes);
	void setZRes(const float _zRes);
	void setZ0(const float _z0);
	void setT0(const float _t0);
	void setTRes(const float _tRes);
	void setSir0(const float _sir0);
	void setSirRes(const float _sirRes);

	void allocate();
	void allocate_sir();
	void permute();

	// finding start and stop indices for model multiplication
	void getStartStopIdx(); // runs start stop idx for both
	void getStartStopIdxNonPermuted();
	void getStartStopIdxPermuted();

	// simple model matrix manipulation functions
	void normalize_matrix();
	void scale(const float _maxAbsVal);

	void exportSensFieldVtk() const;
	void exportSensFieldVtk(const string filePath) const; // exports the sensitivity field to a vtk file
	void saveToFile(const string filePath) const;

	// not implemented yet
	void readFromFile();

	void calcSensField();
	float* get_pmodSlice(const float zLevel); // get slice at z position
	float* get_pmodSlice(const uint64_t iz); // get slice at z index

	float* get_modelVec(const uint64_t ir, const uint64_t iz);
	float* get_modelVec(const float rPos, const float zPos);
	float get_modVecAbsMax() const {return modVecAbsMax;};
};

#endif