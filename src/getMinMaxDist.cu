
#include <cuda.h>
#include "structs.h"


__global__ void getMinMaxDist(
		float* minDist,
		float* maxDist,
		const subElement* transducerElements, 
		const fieldProperties* field,
		unsigned int nElem)
{

	const unsigned int iR = blockIdx.x * blockDim.x + threadIdx.x;
	const unsigned int iZ = blockIdx.y * blockDim.y + threadIdx.y;

	if ((iR < field->nR) && (iZ < field->nZ)){
		const unsigned int outputIdx = iR + iZ * field->nZ;
	 	minDist[outputIdx] = 1e5;
		maxDist[outputIdx] = 0;

		const float rPos = ((float) iR) * field->dr;
		const float zPos = ((float) iZ) * field->dz + field->z0;

		float deltaX, deltaY, deltaZ, dist;

		for (unsigned int iElem = 0; iElem < nElem; iElem++){
			deltaX = transducerElements[iElem].centerPos.x - rPos;
			deltaY = transducerElements[iElem].centerPos.y;
			deltaZ = transducerElements[iElem].centerPos.z - zPos;
			dist = sqrt(deltaX * deltaX + deltaY * deltaY + deltaZ * deltaZ);
			
			if (dist > maxDist[outputIdx])
				maxDist[outputIdx] = dist;

			if (dist < minDist[outputIdx])
				minDist[outputIdx] = dist;
		}
	
	}

	return;
}
