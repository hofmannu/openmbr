#include <string>
#include <cstdio>
#include <iostream>
#include <fstream>
#include "griddedData.h"

using namespace std;

#ifndef VTKWRITER_H
#define VTKWRITER_H

struct unstructuredGrid{
	unsigned int nPoints;
	float* x;
	float* y;
	float* z;
};

struct polydata{
	unsigned int nPoints;
	float* xPoints;
	float* yPoints;
	float* zPoints;
	unsigned int nPolygons;
	unsigned int dimPolygon;
	unsigned int* idPolygons; // holding ids of points 
	// [iCorner + dimPolygon * iPolygon]
};


class vtkwriter{
	public:
		vtkwriter();

		void write();
	
		void set_unstructuredGrid(unstructuredGrid* _unstGrid);
		void set_polydata(polydata* _polDat);
		void set_structuredPoints(griddedData* _structPoints);

		void set_title(const string _title);
		void set_type(const string _type);
		void set_dataType(const string _dataType);
		void set_outputPath(const string _outputPath);
		void set_title(const char* _title, unsigned long length);
		void set_type(const char* _type, unsigned long length);
		void set_dataType(const char* _dataType, unsigned long length);
		void set_outputPath(const char* _outputPath, unsigned long length);
		void set_ascii();
		void set_binary();
	
	private:
		unstructuredGrid* unstrGrid;
		polydata* polData;
		griddedData* structPoints;
		string title = "Unnamed dataset";
		bool flagBinary = 0;
		string dataType = "ASCII";
		string type = "STRUCTURED_POINTS";
		string outputPath = "/home/hofmannu/test.vtk";
		unsigned char idType = 0;
		bool flagVerbose = 1;
		bool flagLittleEndian = 0; // if system is big or little endian
		// STRUCTURED_POINTS --> idType = 0
		// STRUCTURED_GRID --> idType = 1
		// UNSTRUCTURED_GRID --> idType = 2
		// POLYDATA --> idType = 3
		// RECTILINEAR_GRID --> idType = 4
		// FIELD --> idType = 5
};

#endif
