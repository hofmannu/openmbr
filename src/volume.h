/*
	c++ based implementation of a volume handling class
	Author: Urs Hofmann
	Mail: hofmannU@ethz.ch
	Date: 05.08.2020

	Changelog:
		hofmannu - 14.08.2020 - added const to get functions
		hofmannu - 14.08.2020 - added #ifndef loop

	Todo:
		hofmannu - add paraview export functionality
*/

#ifndef VOLUME_H
#define VOLUME_H

#include <string>
#include "baseClass.h"
#include <H5Cpp.h>
#include "basicMathOp.h"
#include "../lib/vtkwriter/vtkwriter.h"
#include "griddedData.cuh"

using namespace std;

class volume : public baseClass, public basicMathOp
{
private:
	bool isMemAlloc = 0; // did we allocate memory
	uint64_t dim[3]; // dimensionailty of volume
	uint64_t nElements; // overall number of elements in 
	float origin[3]; // origin of volume
	float res[3]; // resolution of volume
	float* data; // matrix containing data
	float minVal;
	float maxVal;
	float maxAbsVal;

	// z slice array which will be only updated if new slice is requested
	float* sliceZ;
	uint64_t lastSliceZ = 0;
	float* sliceX;
	uint64_t lastSliceX = 0;
	float* sliceY;
	uint64_t lastSliceY = 0;

public:
	volume(); // class constructor
	~volume(); // class destructor

	// important properties
	float get_minVal() const {return minVal;};
	float get_maxVal() const {return maxVal;};
	float get_maxAbsVal() const {return maxAbsVal;};

	// functions to act on volume
	void multiply(const float factor); // multiplies all elements in volume with factor

	// dimensionality of the dataset	
	void setDim(const uint64_t dim0, const uint64_t dim1, const uint64_t dim2);
	void setDim(const uint64_t* _dim);
	void setDim(const uint8_t _dim, const uint64_t newDim);
	uint64_t getDim(const uint8_t _dim) const; // returns dimension along axis

	// define resolution of the dataset
	void setRes(const float* dx);
	void setRes(const float dx0, const float dx1, const float dx2);
	void setRes(const uint8_t _dim, const float _res);

	float getRes(const uint8_t _dim) const;
	float getRes0() const;
	float getRes1() const;
	float getRes2() const;

	// define origin of the dataset
	void setOrigin(const float* _origin);
	void setOrigin(const uint8_t _dim, const float _origin);
	void setOrigin(const float origin0, const float origin1, const float origin2);
	float getOrigin(const uint8_t _dim);

	// different ways to define matrix elements
	void setValue(const float value); // set whole array to certain value
	void setValue(const unsigned int x0, const unsigned int x1, const unsigned int x2, const float value);
	void setValue(const uint64_t* pos, const float value);


	float getValue(const unsigned int x0, const unsigned int x1, const unsigned int x2) const;
	float getValue(const unsigned int iElem) const;
	float getValue(const unsigned int* pos) const;

	// get position along axis
	float getPos0(const uint64_t idx0) const;
	float getPos1(const uint64_t idx1) const;
	float getPos2(const uint64_t idx2) const;
	float getPos(const uint64_t idx, const uint8_t iDim) const;

	float getCenterPos(const uint8_t _dim); // returns center position along dimension

	// get index of a certain position
	uint64_t getIdx0(const float pos0) const;
	uint64_t getIdx1(const float pos1) const;
	uint64_t getIdx2(const float pos2) const;
	uint64_t getIdx(const float pos, const uint8_t iDim) const;

	float getLength(const uint8_t _dim); 
	// returns length of dataset along a certain dimension

	uint64_t get_nElements() const;
	void allocMemory();

	void readFromFile(const string filePath); // read from h5 file
	void saveToFile(const string filePath) const;

	void printInformation() const;

	float getMinPos(const unsigned int _dim) const;
	float getMaxPos(const unsigned int _dim) const;
	float getRangeLimitedPos(const float pos, const unsigned int _dim) const;

	void getCroppedVolume(
		float* vol, // array containing cropped volume
		const uint64_t start0, const uint64_t stop0,
		const uint64_t start1, const uint64_t stop1,
		const uint64_t start2, const uint64_t stop2) const;

	void getCroppedVolume(
		float* vol, // array containing cropped volume 
		const uint64_t* startIdx, 
		const uint64_t* stopIdx) const;

	void exportVtk(const string filePath);

	void calcMinMax();
	float getMinVal() const {return minVal;};
	float getMaxVal() const {return maxVal;};

	float* get_pdata() {return data;};
	void set_pdata(float* _data);

	// get slices of volume
	float* get_psliceZ(const uint64_t zLevel);
	float* get_psliceX(const uint64_t xLevel);
	float* get_psliceY(const uint64_t yLevel);
	float* get_psliceZ(const float zPos);
	float* get_psliceX(const float xPos);
	float* get_psliceY(const float yPos);
};

#endif
