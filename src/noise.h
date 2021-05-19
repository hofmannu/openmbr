/*
	short object oriented wrapper to add noise to datasets
	Author: Urs Hofmann
	Mail: hofmannu@biomed.ee.ethz.ch
	Date: 09.10.2020
*/

#ifndef NOISE_H
#define NOISE_H

#include <cstring>
#include <iostream>
#include <cstdio>
#include <random>

class noise
{

public:
	// cosntructor / destructor
	noise();
	~noise();

	double* get_pnoiseLevel() {return &noiseLevel;};
	float get_noiseLevel() const {return noiseLevel;};

	void addNoise(float* array, const uint64_t nElements);
	void addGaussianNoise(float* array, const uint64_t nElements);

	void set_noiseLevel(const float _noiseLevel);

private:
	
	uint8_t noiseType = 0;
	// 0: gaussian noise

	double mean = 0;
	double noiseLevel = 1;

};

#endif