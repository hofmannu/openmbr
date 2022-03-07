/*
	File: model.h
	Author: Urs Hofmann
	Mail: mail@hofmannu.org
	Date: 28.02.2022
*/

#ifndef MODEL_H
#define MODEL_H

#include <cstdint>
#include <string>
#include <H5Cpp.h>
#include <cstdlib>
#include <cstring>
#include "../lib/CVolume/src/volume.h"

using namespace std;

class model
{
private:

	float* modelData; // should have order: x, y, t, z
	bool isModelDataAlloc = 0;

	volume sensField; // data order: [x, y, z]

	float res[4] = {4e-9, 20e-6, 20e-6, 20e-6}; // resolution of model [dt, dx, dy, dz]
	float origin[4] = {0.0, -1e-3, -1e-3, -1e-3}; // origin of model vectot [t0, x0, y0, z0]
	uint64_t dim[4] = {200, 101, 101, 101}; // size of model matrix: [nt, nx, ny, nz]
	bool isLoaded = 0; // is model loaded?
	float minVal = 0;
	float maxVal = 0;

	string filePath; // from where do we load the model
	string folderPath; // path where our file is located
	string fileName; // name of the file in the target location

	void alloc_model();
	void set_filePath(const string _filePath);

public:

	// constructor and destructor
	~model();

	// operator overloads
	model& operator = (model& modelIn);

	void crop(const uint64_t* startIdx, const uint64_t* stopIdx);

	void load_from_file(); // use existing filepath to load model
	void load_from_file(const string _filePath);

	void calc_sensField();
	void calc_minmax();

	float* get_pdata() {return modelData;};
	inline float get_val(
		const uint64_t idxX, const uint64_t idxY, 
		const uint64_t idxT, const uint64_t idxZ) const;

	// boring set dunctions
	void set_dim(const uint64_t nt, const uint64_t nx, const uint64_t ny, const uint64_t nz);
	void set_origin(const uint8_t iDim, const float originIn);
	void set_res(const uint8_t iDim, const float resIn);

	// boring get functions
	float get_res(const uint8_t iDim) const {return res[iDim];};
	float get_origin(const uint8_t iDim) const {return origin[iDim];};
	uint64_t get_dim(const uint8_t iDim) const {return dim[iDim];};

	// get functions for resolution of model
	inline float get_dt() const {return res[0];};
	inline float get_dx() const {return res[1];};
	inline float get_dy() const {return res[2];};
	inline float get_dz() const {return res[3];};

	// get funciton for origin of model
	inline float get_t0() const {return origin[0];};
	inline float get_x0() const {return origin[1];};
	inline float get_y0() const {return origin[2];};
	inline float get_z0() const {return origin[3];};

	uint64_t get_tIdx(const float tPos) const;

	inline float get_tEnd() const {return origin[0] + (((float) dim[0] - 1) * res[0]);};
	inline float get_xEnd() const {return origin[1] + (((float) dim[1] - 1) * res[1]);};
	inline float get_yEnd() const {return origin[2] + (((float) dim[2] - 1) * res[2]);};
	inline float get_zEnd() const {return origin[3] + (((float) dim[3] - 1) * res[3]);};

	inline uint64_t get_nt() const {return dim[0];};
	inline uint64_t get_nx() const {return dim[1];};
	inline uint64_t get_ny() const {return dim[2];};
	inline uint64_t get_nz() const {return dim[3];};

	inline uint64_t get_nElements() const {return dim[0] * dim[1] * dim[2] * dim[3];}; 

	inline float get_minVal() const {return minVal;};
	inline float get_maxVal() const {return maxVal;};

	bool get_isLoaded() const {return isLoaded;};
	
	const char* get_filePath() const {return filePath.c_str();};
	const char* get_folderPath() const {return folderPath.c_str();};
	const char* get_fileName() const {return fileName.c_str();};
	
	float get_memory() const; // return the memory of the matrix in bytes

	void print_information();

	volume* get_psensField() {return &sensField;};
};

#endif