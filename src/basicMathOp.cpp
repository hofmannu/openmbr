#include "basicMathOp.h"

void basicMathOp::handlePolarity(float *_array, const unsigned int _nElements, const string _method) const
{
	if (_method.compare("abs"))
	{
		#pragma unroll
		for (unsigned int iElement = 0; iElement < _nElements; iElement++)
			_array[iElement] = fabs(_array[iElement]);			
	}
	else if (_method.compare("pos"))
	{
		#pragma unroll
		for (unsigned int iElement = 0; iElement < _nElements; iElement++)
		{
			if (_array[iElement] < 0)
				_array[iElement] = 0;
		}
	}
	else if (_method.compare("neg"))
	{
		#pragma unroll
		for (unsigned int iElement = 0; iElement < _nElements; iElement++)
		{
			_array[iElement] = (_array[iElement] < 0) ? -_array[iElement] : 0;
		}	
	}
	else if (_method.compare("none"))
	{	
		// nothing to do in this case
	}
	else
	{
		throw "Invalid argument passed to function";
	}

	return;
}

float basicMathOp::getNorm(const float* _array, const unsigned int _nElements) const
{
	float norm = 0;
	#pragma unroll
	for (unsigned int iElement = 0; iElement < _nElements; iElement++)
		norm += (_array[iElement] * _array[iElement]);
	
	norm = pow(norm, 0.5);
	return norm;
}

void basicMathOp::multiply(float *_array, const float _factor, const unsigned int _nElements) const
{
	#pragma unroll
	for (unsigned int iElement = 0; iElement < _nElements; iElement++)
		_array[iElement] *= _factor;

	return;
}

void basicMathOp::multiply(float* _arrayA, const float* _arrayB, const unsigned int _nElements) const
{
	#pragma unroll
	for (unsigned int iElement = 0; iElement < _nElements; iElement++)
		_arrayA[iElement] *= _arrayB[iElement];

	return;
}

void basicMathOp::multiply(float* _arrayA, const float* _arrayB, const float* _arrayC, const unsigned int _nElements) const
{
	#pragma unroll
	for (unsigned int iElement = 0; iElement < _nElements; iElement++)
		_arrayA[iElement] = _arrayB[iElement] + _arrayC[iElement];

	return;
}

void basicMathOp::multiply(float* _arrayA, const float* _arrayB, const float _factor, const unsigned int _nElements) const
{
	#pragma unroll
	for (unsigned int iElement = 0; iElement < _nElements; iElement++)
		_arrayA[iElement] = _arrayB[iElement] * _factor;

	return;
}

void basicMathOp::divide(float* _array, const float _factor, const unsigned int _nElements) const
{
	float _invFactor = 1 / _factor;
	multiply(_array, _invFactor, _nElements);
	return;
}

void basicMathOp::divide(float* _arrayA, const float* _arrayB, const unsigned int _nElements) const
{
	#pragma unroll
	for (unsigned int iElement = 0; iElement < _nElements; iElement++)
		_arrayA[iElement] /= _arrayB[iElement];
	
	return;
}

void basicMathOp::divide(float *_arrayA, const float *_arrayB, const float *_arrayC, const unsigned int _nElements) const
{
	#pragma unroll
	for (unsigned int iElement = 0; iElement < _nElements; iElement++)
		_arrayA[iElement] = _arrayB[iElement] / _arrayC[iElement];

	return;
}

void basicMathOp::divide(float* _arrayA, const float* _arrayB, const float _factor, const unsigned int _nElements) const
{
	#pragma unroll
	for (unsigned int iElement = 0; iElement < _nElements; iElement++)
		_arrayA[iElement] = _arrayB[iElement] / _factor;

	return;
}

void basicMathOp::normalize(float *_array, const uint64_t _nElements) const
{
	const float _norm = getNorm(_array, _nElements);
	divide(_array, _norm, _nElements);
	return;
}

// _arrayA = _arrayA - arrayB
void basicMathOp::substract(float* _arrayA, const float* _arrayB, const unsigned int _nElements) const
{
	#pragma unroll
	for (unsigned int iElement = 0; iElement < _nElements; iElement++)
		_arrayA[iElement] -= _arrayB[iElement];

	return;
}

// _arrayA = _arrayA - value 
void basicMathOp::substract(float* _arrayA, const float _value, const unsigned int _nElements) const
{
	#pragma unroll
	for (unsigned int iElement = 0; iElement < _nElements; iElement++)
		_arrayA[iElement] -= _value;

	return;
}

void basicMathOp::substract(float* _arrayA, const float* _arrayB, const float* _arrayC, const unsigned int _nElements) const
{
	#pragma unroll
	for (unsigned int iElement = 0; iElement < _nElements; iElement++)
		_arrayA[iElement] = _arrayB[iElement] - _arrayC[iElement];
	return;
}


// arrayA = arrayA + arrayB
void basicMathOp::add(float* _arrayA, const float* _arrayB, const unsigned int _nElements) const
{
	#pragma unroll
	for (unsigned int iElement = 0; iElement < _nElements; iElement++)
		_arrayA[iElement] += _arrayB[iElement];

	return;
}

// arrayA = arrayA + value
void basicMathOp::add(float* _arrayA, const float _value, const unsigned int _nElements) const
{
	#pragma unroll
	for (unsigned int iElement = 0; iElement < _nElements; iElement++)
		_arrayA[iElement] += _value;

	return;
}

// arrayA = arrayB + arrayC
void basicMathOp::add(float *_arrayA, const float *_arrayB, const float *_arrayC, const unsigned int _nElements) const
{
	#pragma unroll
	for (unsigned int iElement = 0; iElement < _nElements; iElement++)
		_arrayA[iElement] = _arrayB[iElement] + _arrayC[iElement];

	return;
}

void basicMathOp::assign(float* _arrayOut, const float* _arrayIn, const unsigned int _nElements) const
{
	#pragma unroll
	for (unsigned int iElement = 0; iElement < _nElements; iElement++)
		_arrayOut[iElement] = _arrayIn[iElement];

	return;

}

void basicMathOp::assignRand(float* _array, const uint64_t _nElements) const
{
	#pragma unroll
	for (uint64_t idx = 0; idx < _nElements; idx++)
	{
		_array[idx] = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
	}
	return;
}

float basicMathOp::getMaxAbs(const float *_array, const unsigned int _nElements) const
{
	float maxAbsVal = 0;
	#pragma unroll
	for (unsigned int iElement = 0; iElement < _nElements; iElement++){
		if (fabs(_array[iElement]) > maxAbsVal){
			maxAbsVal = fabs(_array[iElement]);
		}
	}
	return maxAbsVal;	
}

float basicMathOp::getMin(const float* _array, const uint64_t _nElements)
{
	float minVal = _array[0];
	for (uint64_t iElement = 0; iElement < _nElements; iElement++)
	{
		if (_array[iElement] < minVal)
			minVal = _array[iElement]; 
	}
	return minVal;
}

float basicMathOp::getMax(const float* _array, const uint64_t _nElements)
{
	float maxVal = _array[0];
	for (uint64_t iElement = 0; iElement < _nElements; iElement++)
	{
		if (_array[iElement] > maxVal)
			maxVal = _array[iElement];
	}
	return maxVal;
}