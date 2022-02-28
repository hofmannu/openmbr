#include "../src/reconsett.h"
#include <catch2/catch.hpp>

int main()
{

	reconsett sett;

	// check if regularization strength works fine
	sett.set_regStrength(0.1f);
	if (sett.get_regStrength() != 0.1f)
	{
		printf("Could not return the correct reg strength\n");
		throw "InvalidValue";
	}

	// check if number of iterations works fine
	sett.set_nIter(100);
	if (sett.get_nIter() != 100)
	{
		printf("Could not retrieve the correct number of iterations\n");
		throw "Invalid value";
	}

	// check all cropping
	sett.set_tMin(0.123f);
	if (sett.get_tMin() != 0.123f)
	{
		printf("Could not define minimum time point\n");
		throw "Invalid value";
	}

	// check all cropping
	sett.set_tMax(0.323f);
	if (sett.get_tMax() != 0.323f)
	{
		printf("Could not define minimum time point\n");
		throw "Invalid value";
	}

	return 0;
}