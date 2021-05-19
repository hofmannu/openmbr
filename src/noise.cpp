#include "noise.h"

// class construtor and destructor
noise::noise()
{

}

noise::~noise()
{
	
}

// define the amplitude of the noise
void noise::set_noiseLevel(const float _noiseLevel)
{
	noiseLevel = _noiseLevel;
	return;
}

// add noise on top of array
void noise::addNoise(float* array, const uint64_t nElements)
{

	// switch case depending on noise type
	if (noiseType == 0)
	{
		addGaussianNoise(array, nElements);
	}
	else
	{
		printf("Invalid noise tye defined\n");
		throw "invalidValue";
	}


	return;
}

// add gaussian noise on top of array
void noise::addGaussianNoise(float* array, const uint64_t nElements)
{

  std::default_random_engine generator;
  std::normal_distribution<double> dist(mean, noiseLevel);
	// case gaussian noise
	for (uint64_t iElement = 0; iElement < nElements; iElement++)
	{
		array[iElement] = array[iElement] + (float) dist(generator);
	}

	return;
}