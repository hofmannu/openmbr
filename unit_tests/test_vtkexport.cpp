/*
	text_vtkexport.cpp
	unit test file to check the functionality of the vtkExport class

	Author: Urs Hofmann
	Mail: hofmannu@ethz.ch
	Date: 20.08.2020
*/

#include <cmath>
#include "../src/vtkWriter.h"
#include "../src/griddedData.cuh"
#include <iostream>

using namespace std;

int main()
{

	vtkWriter myWriter;
	griddedData testData;
	
	const unsigned int nX = 400;
	const unsigned int nY = 400;
	const unsigned int nZ = 200;

	testData.dim[0] = nZ; // z direction for our images
	testData.dim[1] = nX; // x direction
	testData.dim[2] = nY; // y direction

			
	testData.nElements = 1;
	for (unsigned char iDim = 0; iDim < 3; iDim++)
	{
		testData.res[iDim] = 1;
		testData.origin[iDim] = 0;
		testData.nElements *= testData.dim[iDim];
	}

	testData.data = new float [testData.nElements];

	// create something like a cone
	float dx, dy, dr;
	unsigned int idxVol;
	for (unsigned int ix = 0; ix < nX; ix++)
	{
		dx = ((float) ix) - ((float) nX) / 2.0;
		for (unsigned int iy = 0; iy < nY; iy++)
		{
			dy = ((float) iy) - ((float) nY) / 2.0;
		 	dr = sqrt(dx * dx + dy * dy);
		 	for (unsigned int iz = 0; iz < nZ; iz++)
		 	{
		 		idxVol = iz + ix * nZ + iy * nZ * nX;
		 		testData.data[idxVol] = dr;
		 	}
		}
	}

	string title ("reconVol");
	myWriter.set_title(title);

	string type ("STRUCTURED_POINTS");
	myWriter.set_type(type);

	string outputPath ("/home/hofmannu/unitTestVtkExport.vtk");
	myWriter.set_outputPath(outputPath);

	myWriter.set_structuredPoints(&testData);

	myWriter.set_binary();

	myWriter.write();

	delete[] testData.data;

	return 0;
}