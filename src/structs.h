#include <cmath>

// #ifndef TRANSDUCERPROPERTIES_H
// #define TRANSDUCERPROPERTIES_H
// 	// gemoetrical properties of our transducer
// 	struct transducerProperties{
// 		float focalDistance;
// 		float rAperture;
// 		float rHole;
// 		unsigned int nSigma; // number of elements for spherical cosy
// 		unsigned int nPhi; // number of elements for spherical cosy
// 		unsigned int nElements; // number of elements for Kurihara approach
// 	};
// #endif

#ifndef FIELDSTRUCT_H
#define FIELDSTRUCT_H
struct fieldStruct
{
	float fs;
	float dt; // temporal resolution of model matrix
	float dr; // resolution of field in radial direction
	float dz; // resolution of field in axial direction
	float SOS; // assumed speed of soud
	unsigned int nr; // number of elements in radial direction
	unsigned int nz; // number of elements in axial direction
	unsigned int nElements;
	float z0; // starting position of field in z
	float dSir; // resolution of field in spatial impulse response along z 
	float sir0; // starting indeo of spatial impulse response in z
	float t0; // start time of model matrix
};
#endif

#ifndef POS_H
#define POS_H
struct pos{
	float x;
	float y;
	float z;
};
#endif

#ifndef IDX3
#define IDX3
class idx3{
	public:
		unsigned int id0; 
		unsigned int id1; 
		unsigned int id2;
		
};
#endif

#ifndef POINT_GRID
#define POINT_GRID
class point_grid{
	public:
	float phi;
	float theta2;
	float x;
	float y;
	float z;
	
	// calculates the element position based on angles and radius
	void calc_pos(const float radius){
		x = radius * cos(phi) * sin(theta2);
		y = radius * sin(phi) * sin(theta2);
		z = -radius * cos(theta2);
		return;
	};
};
#endif

#ifndef SUBELEMENT_H
#define SUBELEMENT_H
// transducer subelement
struct subElement{
	pos centerPos; // 3d position of transducer element relative to focus
	float area; // physical area of transducer element
	float rDist; // radial distance from acoustic axis
};
#endif
