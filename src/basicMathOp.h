/*
	File: basicMathOp.h
	Author: Urs Hofmann
	Mail: hofmannu@ethz.ch
	Date: 05.08.2020

	Description: Implementation of basic mathematical functions applied over arrays
*/

#ifndef BASICMATHOP_H
#define BASICMATHOP_H

#include <cmath>
#include <iostream>
#include <cstdlib>

using namespace std;

class basicMathOp
{
	public:
		void handlePolarity(float* _array, const unsigned int _nElements, const string _method) const; 

		float getNorm(const float* _array, const unsigned int _nElements) const;
		
		// MULTIPLICATION
		// _array = _array * factor
		void multiply(float* _array, const float _factor, const unsigned int _nElements) const;
		// _arrayA = _arrayA * _arrayB
		void multiply(float* _arrayA, const float*  _arrayB, const unsigned int _nElements) const;
		// _arrayA = _arrayB * arrayC
		void multiply(float* _arrayA, const float* _arrayB, const float* _arrayC, const unsigned int _nElements) const;
		// _arrayA = _arrayB * _factor
		void multiply(float* _arrayA, const float* _arrayB, const float factor, const unsigned int _nElements) const;

		// DIVISIONS
		// _array = _array / _factor
		void divide(float* _array, const float _factor, const unsigned int _nElements) const;
		// _arrayA = _arrayA / _arrayB
		void divide(float* _arrayA, const float* _arrayB, const unsigned int _nElements) const;
		// _arrayA = _arrayB / _arrayC
		void divide(float* _arrayA, const float* _arrayB, const float* _arrayC, const unsigned int _nElements) const;
		// _arrayA = _arrayB / _factor
		void divide(float* _arrayA, const float* _arrayB, const float factor, const unsigned int _nElements) const;

		void normalize(float* _array, const uint64_t _nElements) const;
		
		// SUBSTRACTIONS
		// arrayA = arrayA - arrayB
		void substract(float* _arrayA, const float* _arrayB, const unsigned int _nElements) const; 
		// arrayA = arrayA - value 
		void substract(float* _arrayA, const float _value, const unsigned int _nElements) const;
		// arrayA = arrayB - arrayC
		void substract(float* _arrayA, const float* _arrayB, const float* _arrayC, const unsigned int _nElements) const;

		// ADDITIONS
		// arrayA = arrayA + arrayB
		void add(float* _arrayA, const float* _arrayB, const unsigned int _nElements) const;
		// arrayA = arrayA + value
		void add(float* _arrayA, const float _value, const unsigned int _nElements) const;
		// arrayA = arrayB + arrayC
		void add(float* _arrayA, const float* _arrayB, const float* _arrayC, const unsigned int _nElements) const;

		// set an array to equal elements for its full length _arrayOut = _arrayIn
		void assign(float* _arrayOut, const float* _arrayIn, const unsigned int _nElements) const;

		// set all elements in array to random values
		void assignRand(float* _array, const uint64_t _nElements) const;

		float getMaxAbs(const float* _array, const unsigned int _nElements) const;
		// returns maximum absolute value in array
		float getMin(const float* _array, const uint64_t _nElements);
		float getMax(const float* _array, const uint64_t _nElements);


};

#endif