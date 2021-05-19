#include "griddedData.h"

// define resolution of dataset
void griddedData::set_res(const float* _res){
	for (unsigned char iDim = 0; iDim < 3; iDim++){
		if (_res[iDim] <= 0)
			throw "Invalid resolution for grid";
		res[iDim] = _res[iDim];
		ires[iDim] = 1 / _res[iDim];
	}
	return;
}

// set dimension of dataset
void griddedData::set_dim(const unsigned int* _dim){
	for (unsigned char iDim = 0; iDim < 3; iDim++)
		dim[iDim] = _dim[iDim];
	return;
}

// define position of element [0, 0, 0] in space
void griddedData::set_origin(const float* _origin){
	for (unsigned char iDim = 0; iDim < 3; iDim++)
		origin[iDim] = _origin[iDim];
	return;
}

// print information about dataset
void griddedData::print_info(){
	printf(" - Size:       %i, %i, %i\n", dim[0], dim[1], dim[2]);
	printf(" - Resolution: %f, %f, %f\n", res[0], res[1], res[2]);
	printf(" - Offset:     %f, %f, %f\n", origin[0], origin[1], origin[2]);
	return;
}