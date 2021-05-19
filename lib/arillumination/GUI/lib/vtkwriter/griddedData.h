/*
	c++ class representing a simple gridded dataset
	Author: Urs Hofmann
	Mail: hofmannu@ethz.ch
	Date: 14.08.2020
*/

#ifndef GRIDDEDDATA_H
#define GRIDDEDDATA_H

#include <cstdio>
#include <cstring>
#include <fstream>

using namespace std;

class griddedData{
	private:
	public:
		float ires[3]; // inverted resolution
		float * data; // data pointer
		unsigned int dim[3]; // dimensions of volume [z, x, y]
		unsigned int nElements;
		float res[3]; // resolution of volume
		float origin[3]; // position of voxel [0,0,0]
		void set_origin(const float* _origin);
		void set_res(const float* _res);
		void set_dim(const unsigned int* _dim);
		// unsigned int get_nElements() {return dim[0] * dim[1] * dim[2];};
		void print_info();
};

#endif