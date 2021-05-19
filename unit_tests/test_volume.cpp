#include <stdio.h>
#include "../src/volume.h"

int main()
{
	printf("Running unit test for testVolume...\n");

	volume testVolume;
	testVolume.setRes(0.1, 0.1, 0.1); // define resolution of testVolume
	testVolume.setOrigin(0.5, 0, 0); // define origin of testVolume
	testVolume.setDim(1024, 1024, 1024); // define dimension of test volume

	printf("Allocating memory...\n");
	testVolume.allocMemory();

	printf("Assigning constant value to memory...\n");
	testVolume.setValue(1);
	testVolume.setValue(0, 0, 0, 2);

	printf("Return value of pos (0, 0, 0): %1.0f \n", testVolume.getValue(0, 0, 0));

	printf("Multiplying volume by 2...\n");
	testVolume.multiply(2);

	printf("Return value of pos (0, 0, 0): %1.0f \n", testVolume.getValue(0, 0, 0));

	printf("Saving dataset to file...\n");
	testVolume.saveToFile("/home/hofmannu/test.h5");

	printf("Loading dataset from file...\n");
	volume testVolume2;
	testVolume2.readFromFile("/home/hofmannu/test.h5");

	printf("Return value of pos (0, 0, 0): %1.0f \n", testVolume.getValue(0, 0, 0));
	return 0;
}