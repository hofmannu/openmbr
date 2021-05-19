#include "transModel.cuh"

// destructor
transModel::~transModel()
{
	// free all memory which was allocated
	if (isImpRespDefined)
	{
		delete[] impResp;
		delete[] tImp;
	}

	if (isTransdElemAlloc)
		delete[] transducerElements;

	if (isPGridAlloc)
		delete[] pGrid;

	if (isTriangleIdRAlloc)
		delete[] triangleIdR;

}

// constant get functions
uint64_t transModel::getNtModel() const {return model.get_nt();}
float transModel::getStarttime() const {return starttime;}

// non constant get functions
float* transModel::getImpResp() {return impResp;}

// set functions
void transModel::setNElements(const uint64_t _nElements)
{
	if (_nElements <= 0)
		throw "nElements must be bigger then 0";
	else
		nElements = _nElements;
	return;
}

void transModel::setFieldProperties(const fieldProperties _fieldProp)
{
	fieldProp = _fieldProp;
	return;
}

// set properties of transducer like focal distance etc.
void transModel::setTransProp(const transducerProperties _prop)
{
	prop = _prop;
	return;
}

void transModel::checkCudaReturn(const cudaError_t err, const string errMsg)
{
	if (err != cudaSuccess){
		printf("CUDA error string: ");
		printf(cudaGetErrorString(err));
		printf("\n");
		printf(errMsg.c_str());
		printf("\n");
		throw "cudaError";
	}
	return;
}

void transModel::vprintf(const string txtMsg)
{
	if (flagVerbose){
		printf(txtMsg.c_str());
		fflush(stdout);
	}
	return;
}

// if nothing is passed, generate file name automatically
void transModel::loadImpulseResponse()
{
	const string filePath = "transFiles/" + prop.getName() + "_imp.h5";
	loadImpulseResponse(filePath);
	return;
}

void transModel::loadImpulseResponse(const string filePath)
{
	// structure in file
	// fSampling 		single 		sampling frequency of signal [Hz]
	// nSgnl									number of elements in signal
	// sgnl 				single		array containing the actual signal

	H5::H5File file(filePath, H5F_ACC_RDONLY);

	// load sampling frequency from file
	H5::DataSet fsDataset = file.openDataSet("fSampling");
	const hsize_t col_dims = 1;
	H5::DataSpace mspaceFs (1, &col_dims); 
	H5::DataSpace filespace = fsDataset.getSpace();
	fsDataset.read(&fSampling, H5::PredType::NATIVE_FLOAT, mspaceFs, filespace);
	dt = 1 / fSampling;	
	// printf("Sampling frequency read: %f\n", fSampling);

	// load length of impulse signal from file
	H5::DataSet nsDataset = file.openDataSet("nSgnl");
	H5::DataSpace mspaceNs (1, &col_dims);
	filespace = nsDataset.getSpace();
	nsDataset.read(&nElem, H5::PredType::NATIVE_UINT64, mspaceNs, filespace);
	printf("Number of samples in signal: %d\n", nElem);

	// load time signal from file
	impResp = new float[nElem];
	tImp = new float[nElem];
	H5::DataSet sgnlDataset = file.openDataSet("sgnl");
	const hsize_t col_sgnl = nElem;
	H5::DataSpace mspaceSgnl (1, &col_sgnl);
	filespace = sgnlDataset.getSpace();
	sgnlDataset.read(impResp, H5::PredType::NATIVE_FLOAT, mspaceSgnl, filespace);
	// printf("Number of samples in signal: %f\n", impResp[0]);

	// extract max and min value and their position
	float maxVal = 0; 
	float minVal = 0; 
	unsigned int posMin, posMax;
	for (unsigned int iElem = 0; iElem < nElem; iElem++){
		if (impResp[iElem] > maxVal){
			maxVal = impResp[iElem];
			posMax = iElem;
		}
		if (impResp[iElem] < minVal){
			minVal = impResp[iElem];
			posMin = iElem;
		}
	}

	tOffset = ((float) (posMax + posMin)) * 0.5 * dt;
	// printf("offset in time vector: %.0f ns\n", tOffset * 1e9);

	// normalize vector to have a maximum amplitude of 1
	if (abs(minVal) > abs(maxVal))
		maxVal = abs(minVal);

	for (unsigned int iElem = 0; iElem < nElem; iElem++)
		impResp[iElem] /= maxVal;

	// create corresponding time vector
	for (uint32_t iElem = 0; iElem < nElem; iElem++)
		tImp[iElem] = tOffset + (float) iElem / fSampling;
	// double check what max abs of vector is
	// float maxAbsVal = 0;
	// for (unsigned int iElem = 0; iElem < nElem; iElem++)
	// {
	// 	if (abs(impResp[iElem]) > maxAbsVal)
	// 		maxAbsVal = impResp[iElem];
	// }
	// printf("Maximum detected value after normalization: %f\n", maxAbsVal);

	file.close();
}

inline float dot(float* vec1, float* vec2, unsigned int nElem)
{
	float dotProduct = 0;
	for (unsigned int iElem = 0; iElem < nElem; iElem++)
		dotProduct += vec1[iElem] * vec2[iElem];

	return dotProduct;
}

inline double dot(double* vec1, double* vec2, unsigned int nElem)
{
	double dotProduct = 0;
	for (unsigned int iElem = 0; iElem < nElem; iElem++)
		dotProduct += vec1[iElem] * vec2[iElem];

	return dotProduct;
}

float transModel::spherical_triangle_area(
	double *r1, double *r2, double *r3, const float R)
{
	// calculates the area of the spherical triangle given
	// by the cartesian vectors r1, r2 and r3 of shape 3x0
	// ad radius R. */
								  
	const double aux1 = 1 / ((double) R);
	
	// Normalising vectors
	#pragma unroll
	for (uint8_t n = 0; n < 3; n++){
		r1[n] = r1[n] * aux1;
		r2[n] = r2[n] * aux1;
		r3[n] = r3[n] * aux1;
	}

	// Calculating the cosines and sines
	const double cosa = dot(r2, r3, 3); 
	const double sina = sqrt(1 - cosa * cosa);
	const double cosb = dot(r1, r3, 3); 
	const double sinb = sqrt(1 - cosb * cosb);
	const double cosc = dot(r1, r2, 3); 
	const double sinc = sqrt(1 - cosc * cosc);
	                                  
	// Calculating angles
	const double A = acos((cosa - cosb * cosc) / (sinb * sinc));
	const double B = acos((cosb - cosa * cosc) / (sina * sinc));
	const double C = acos((cosc - cosa * cosb) / (sina * sinb));
	                                          

	const double dA = R * R * (A + B + C - M_PI);
	const float dAfloat = (float) dA;            

	return dAfloat;
}

void transModel::discretizeTransducer()
{
	// calculate number of points required for grid
	unsigned int nAzi = round(sqrt((float) nElements / 4)); // number of discret steps in azi dir
	nElements = nAzi * nAzi * 4; // update actual number of elements
	nPoints = 2 * nAzi * nAzi + 2 * nAzi + 1; // number of points in grid

	// define arrays needed for point calculations
	if (isPGridAlloc)
		delete[] pGrid;

	pGrid = new point_grid [nPoints];
	isPGridAlloc = 1;

	idx3* triangleId = new idx3 [nElements * 3];

	// if already allocated, first free existing memory
	subElement* subElems = new subElement[nElements];

	// get thetaMin and thetaMax based on transducer properties
	const float thetaMin = prop.getThetaHole();
	const float thetaMax = prop.getTheta();;

	unsigned int count = 1;
	pGrid[0].theta2 = thetaMin; 
	pGrid[0].phi = 0;
	#pragma unroll
	for (unsigned int iAzi = 0; iAzi <= nAzi; iAzi++){
		#pragma unroll
		for (unsigned int iPhi = 0; iPhi < (4 * iAzi); iPhi++){
			pGrid[count].theta2 = thetaMax * ((float) iAzi) / ((float) nAzi);
			pGrid[count].phi = M_PI * 0.5 * ((float) iPhi) / ((float) iAzi);
			count++;
		}
	}

	#pragma unroll
	for (unsigned int iPoint = 0; iPoint < nPoints; iPoint++)
		pGrid[iPoint].calc_pos(prop.getFocalDistance());

	int upDown, du, dl; // triangle direction and increment flag
	
	// zero level manually
	triangleId[0].id0 = 0; triangleId[0].id1 = 1; triangleId[0].id2 = 2;
	triangleId[1].id0 = 0; triangleId[1].id1 = 2; triangleId[1].id2 = 3;
	triangleId[2].id0 = 0; triangleId[2].id1 = 3; triangleId[2].id2 = 4;
	triangleId[3].id0 = 0; triangleId[3].id1 = 4; triangleId[3].id2 = 1;
	count=4;
	for (int iAzi = 1; iAzi < nAzi; iAzi++) { // number of levels                 
		upDown = -1;                                                                  
		du = 1; dl = 1;                                                                 
		for (int m = 1; m < 4 * (2 * (iAzi + 1) - 1) + 1; m++) { // triangles in one level  
			if (m == 4 * (2 * (iAzi + 1) - 1)){             
				//print 'hey - '+ repr(count)                                           
				triangleId[count].id0 = 2 * iAzi * iAzi - 2 * iAzi + 1;   
				triangleId[count].id1 = 2 * (iAzi + 1) * (iAzi + 1) - 2 * (iAzi + 1) + dl;
				triangleId[count].id2 = 2 * (iAzi + 1) * (iAzi + 1) - 2 * (iAzi + 1) + 1; 
				dl++;                                                                  
			}else{                                                                    
				if (m==4*(2*(iAzi+1)-1)-1){                                                
					//print 'hey + '+ repr(count)                                           
					triangleId[count].id0 = 2 * iAzi * iAzi - 2 * iAzi + du;  
					triangleId[count].id1 = 2 * (iAzi + 1) * (iAzi + 1) - 2 * (iAzi + 1) + dl;
					triangleId[count].id2 = 2 * iAzi * iAzi - 2 * iAzi + 1; 
					du++;                                                                  
				}else{                                                                  
					if (upDown<0) {                                                       
					//print '- '+ repr(count)                                           
					triangleId[count].id0 = 2 * iAzi * iAzi - 2 * iAzi + du;
					triangleId[count].id1 = 2 * (iAzi + 1) * (iAzi + 1) - 2 * (iAzi + 1) + dl;
					triangleId[count].id2 = 2 * (iAzi + 1) * (iAzi + 1) - 2 * (iAzi + 1) + dl + 1;
					dl++;                                                              
					}else{                                                                
						//print '+ '+ repr(count)                                           
						triangleId[count].id0 = 2 * iAzi * iAzi - 2 * iAzi + du;
						triangleId[count].id1 = 2 * (iAzi + 1) * (iAzi + 1) - 2 * (iAzi + 1) + dl;
						triangleId[count].id2 = 2 * iAzi * iAzi - 2 * iAzi + du + 1; 
						du++;                                                              
					}                                                                     
				}                                                                       
			}                                                                         
			if ((m == 2 * (iAzi + 1) - 1) || (m == 2 * (2 * (iAzi + 1) - 1)) || (m == 3 * (2 * (iAzi + 1) - 1)) ) { // corners
				upDown = -1;
			}                                                                         
			else {                                                                    
				upDown = -upDown;
			}                                                                         
			count++;                                                                 
		}                                                                           
	}  
	
	double *aux0 = new double [3];                                                  
	double *aux1 = new double [3];                                                  
	double *aux2 = new double [3];                                                  
	
	// float phiAv, thetaAv;

	float xTemp, yTemp, zTemp, norm;

	#pragma unroll	
	for (int iElem = 0; iElem < nElements; iElem++){
		
		// calculate element area
		// convert all distances from m to microm
		aux0[0] = (double) pGrid[triangleId[iElem].id0].x * 1e6;
		aux0[1] = (double) pGrid[triangleId[iElem].id0].y * 1e6; 
		aux0[2] = (double) pGrid[triangleId[iElem].id0].z * 1e6;
		
		aux1[0] = (double) pGrid[triangleId[iElem].id1].x * 1e6; 
		aux1[1] = (double) pGrid[triangleId[iElem].id1].y * 1e6; 
		aux1[2] = (double) pGrid[triangleId[iElem].id1].z * 1e6;
		
		aux2[0] = (double) pGrid[triangleId[iElem].id2].x * 1e6; 
		aux2[1] = (double) pGrid[triangleId[iElem].id2].y * 1e6; 
		aux2[2] = (double) pGrid[triangleId[iElem].id2].z * 1e6;
		
		subElems[iElem].area = spherical_triangle_area(
			aux0, aux1, aux2, prop.getFocalDistance() * 1e6);// in (µm^2)
		// printf("Area of element %d is %f\n", iElem, subElems[iElem].area);

		// calculate x, y, and z position based on averaging x, y, and z
		xTemp = (pGrid[triangleId[iElem].id0].x +
				pGrid[triangleId[iElem].id1].x + 
				pGrid[triangleId[iElem].id2].x) / 3;         
		yTemp = (pGrid[triangleId[iElem].id0].y +
				pGrid[triangleId[iElem].id1].y +
				pGrid[triangleId[iElem].id2].y) / 3; 
		zTemp = (pGrid[triangleId[iElem].id0].z +
				pGrid[triangleId[iElem].id1].z + 
				pGrid[triangleId[iElem].id2].z) / 3;         
	
		// calculate direction vector norm
		norm = sqrt(xTemp * xTemp + yTemp * yTemp + zTemp * zTemp);

		// get center position back 	from spherical coordinate system
		subElems[iElem].centerPos.x = prop.getFocalDistance() * xTemp / norm;
		subElems[iElem].centerPos.y = prop.getFocalDistance() * yTemp / norm;
		subElems[iElem].centerPos.z = prop.getFocalDistance() * zTemp / norm;
		
		// calculate radial distance from acoustic axis for later sorting
		subElems[iElem].rDist = sqrt(
				subElems[iElem].centerPos.x * subElems[iElem].centerPos.x + 
				subElems[iElem].centerPos.y * subElems[iElem].centerPos.y);
	}                                                                             

	// sort out all transducer elements which don't have correct radius
	unsigned int nElementsR = 0;
	#pragma unroll
	for(unsigned int iElem = 0; iElem < nElements; iElem++){
		if (subElems[iElem].rDist >= prop.getRHole())
			nElementsR++;
	}

	// allocate memory for transducer elements and set flag to 1
	if (isTransdElemAlloc)
		delete[] transducerElements;
	
	transducerElements = new subElement [nElementsR];
	isTransdElemAlloc = 1;

	if (isTriangleIdRAlloc)
		delete[] triangleIdR;

	triangleIdR = new idx3 [nElementsR];
	isTriangleIdRAlloc = 1;

	unsigned int counterElem = 0;
	#pragma unroll
	for (unsigned int iElem = 0; iElem < nElements; iElem++){
		if (subElems[iElem].rDist >= prop.getRHole()){
			transducerElements[counterElem] = subElems[iElem];
			triangleIdR[counterElem] = triangleId[iElem];
			counterElem++;
		}
	}
	
	nElements = nElementsR; // update true number of elements in class
	
	// delete all dynamically allocated arrays
	delete[] triangleId;
	delete[] aux0;
	delete[] aux1;
	delete[] aux2;
	delete[] subElems;
	
	// return pointer to subelements
	return;
}

void transModel::exportElementsVtk()
{
	polydata polData; // only used for export?

	// define points of our grid
	polData.nPoints = nPoints;
	polData.xPoints = new float [nPoints];
	polData.yPoints = new float [nPoints];
	polData.zPoints = new float [nPoints];
	#pragma unroll 
	for (unsigned int iPoint = 0; iPoint < nPoints; iPoint++){
		polData.xPoints[iPoint] = pGrid[iPoint].z;
		polData.yPoints[iPoint] = pGrid[iPoint].x;
		polData.zPoints[iPoint] = pGrid[iPoint].y;
	}

	// define elements
	polData.dimPolygon = 3;
	polData.nPolygons = nElements;
	polData.idPolygons = new unsigned int [nElements * 3];
	#pragma unroll
	for (unsigned int iElem = 0; iElem < nElements; iElem++){
		polData.idPolygons[iElem * 3] = triangleIdR[iElem].id2;
		polData.idPolygons[iElem * 3 + 1] = triangleIdR[iElem].id0;
		polData.idPolygons[iElem * 3 + 2] = triangleIdR[iElem].id1;
	}

	vtkwriter outputter;
		
	outputter.set_polydata(&polData);

	string dataType ("ASCII");
	outputter.set_dataType(dataType); //, sizeof(dataType));
	
	string type ("POLYDATA");
	outputter.set_type(type);
	
	string title ("Transducer Discretization Grid");
	outputter.set_title(title);
	
	string outputPath ("/home/hofmannu/transducerDiscretization.vtk");
	outputter.set_outputPath(outputPath);

	outputter.write();

	delete[] polData.xPoints;
	delete[] polData.yPoints;
	delete[] polData.zPoints;
	delete[] polData.idPolygons;
	return;
}

#ifndef CONSTGETSIRARGS_H
#define CONSTGETSIRARGS_H

struct constGetSirArgs{
	subElement* trans;
	fieldStruct* field;
	uint64_t nElem;
	uint64_t nzSir;
};

#endif

__global__ void getSIR(
	float* sir_dev, // izsir, irfield, izfield
	const constGetSirArgs* inArgsConst
		)
{

	const unsigned int iR = blockIdx.x * blockDim.x + threadIdx.x;
	const unsigned int iZ = blockIdx.y * blockDim.y + threadIdx.y;

	if ((iR < inArgsConst->field->nr) && (iZ < inArgsConst->field->nz)){
		// reset everything in here to zeros before we sum up
		const unsigned int sirOffset = inArgsConst->nzSir * (iR + iZ * inArgsConst->field->nr);
		#pragma unroll
		for (unsigned int iSir = 0; iSir < inArgsConst->nzSir; iSir++)
			sir_dev[iSir + sirOffset] = 0;
		
		const float rField = __fmul_rn((float) iR, inArgsConst->field->dr);
		const float zField = __fadd_rn(__fmul_rn((float) iZ, inArgsConst->field->dz), 
			inArgsConst->field->z0);
		// printf("Field coordinates: %f, %f\n", rField, zField);
		float deltaX, deltaZ, ratio;
		unsigned int distI, index;

		# pragma unroll
		for (unsigned int iElem = 0; iElem < inArgsConst->nElem; iElem++){
			// calculate distance between field point and subelement
			// deltaX = __fsub_rn(rField, trans[iElem].centerPos.x);
			// deltaY = -trans[iElem].centerPos.y;
			// deltaZ = __fsub_rn(zField, trans[iElem].centerPos.z);
			// dist = __fsqrt_rn(
			// 		__fmul_rn(deltaX, deltaX) + 
			// 		__fmul_rn(deltaY, deltaY) + 
			// 		__fmul_rn(deltaZ, deltaZ));

			deltaX = __fsub_rn(rField, inArgsConst->trans[iElem].centerPos.x);
			// deltaY = -trans[iElem].centerPos.y;
			deltaZ = __fsub_rn(zField, inArgsConst->trans[iElem].centerPos.z);
			
			// we reuse deltaX here as dist to save memory
			deltaX = __fsqrt_rn(
					__fmul_rn(deltaX, deltaX) + 
					__fmul_rn(-inArgsConst->trans[iElem].centerPos.y, -inArgsConst->trans[iElem].centerPos.y) + 
					__fmul_rn(deltaZ, deltaZ));
			
			// convert distance in index
			deltaX = __fdividef(__fsub_rn(deltaX, inArgsConst->field->sir0), inArgsConst->field->dSir);
			
			// convert distance into index and interpolate
			distI = (unsigned int) deltaX;
			index = distI + sirOffset;
			ratio = __fsub_rn(deltaX, (float) distI);
			
			sir_dev[index] = __fadd_rn(sir_dev[index], 
					__fmul_rn(inArgsConst->trans[iElem].area, __fsub_rn(1, ratio))); 
			sir_dev[index + 1] = __fadd_rn(sir_dev[index + 1], 
					__fmul_rn(inArgsConst->trans[iElem].area, ratio)); 
		}
	}

	return;
};

// kernel is parallelized over full model matrix (3D)
__global__ void convolveSIR(
		float* mMatrix, // model matrix [it, ir, iz]
		const float* sir, // spatial impulse response [it, ir, iz]
		const float* exSig, // signal vector of length nT
		fieldStruct* field,
		const uint64_t nzSir, // number of elements along z in spatial impulse response
		const uint64_t nTEx, // number of elements in signal vector
		const uint64_t ntModel, // size of model matrix 1st dimension
		const float t0Ex // origin of excitation signal in s
		)
{

	// get current position in output volume
	const unsigned int iT = blockIdx.x * blockDim.x + threadIdx.x;
	const unsigned int iR = blockIdx.y * blockDim.y + threadIdx.y;
	const unsigned int iZ = blockIdx.z * blockDim.z + threadIdx.z;

	if ((iR < field->nr) && (iZ < field->nz) && (iT < ntModel)){
		// output index in model matrix
		const uint64_t iModelMatrix = iT + ntModel * (iR + iZ * field->nr);

		mMatrix[iModelMatrix] = 0; // reset output value to 0

		const float tVoxel = (float) iT * field->dt + field->t0; // time point of current matrix element
		const float tGesEx = (float) nTEx * field->dt; // overall length of excitation signal
		const float tMin = tVoxel - t0Ex; // lower end of convolution time frame 
		const float tMax = tMin + tGesEx; // upper end of convolution time frame

		const float sirMax = field->sir0 + ((float) nzSir) * field->dSir;
		const float zMin = ((tMin * field->SOS) < field->sir0) ? field->sir0 : (tMin * field->SOS);
		const float zMax = ((tMax * field->SOS) > sirMax) ? sirMax : (tMax * field->SOS);

		// convert from distance into start and stop index in z
		const uint64_t startIdx = (zMin - field->sir0) / field->dSir;
		const uint64_t stopIdx = (zMax - field->sir0) / field->dSir;

		const uint64_t sirOffset = nzSir * (iR + field->nr * iZ);

		#pragma unroll
		for (uint64_t idxSir = startIdx; idxSir <= stopIdx; idxSir++){
			float tSir = field->sir0 + (float) idxSir * field->dSir; // until now position
			tSir = tSir / field->SOS; // convert from position into time
			float deltaT = tVoxel - tSir;	// find position in excitation vector

			float idxEx = (t0Ex + deltaT) / field->dt; // accurate index in excitation signal
			if ((idxEx >= 0) && (idxEx <= (nTEx - 1)))
			{
				uint64_t idxExInt = idxEx; // cut down to integer value
				float ratio = idxEx - idxExInt; // ratio between upper and lower part [0...1]
				// do linear interpolation of excitation signal
				float exSigInterp = (1 - ratio) * exSig[idxExInt] + ratio * exSig[idxExInt + 1];

				mMatrix[iModelMatrix] = __fadd_rn(mMatrix[iModelMatrix], 
				 	__fmul_rn(sir[sirOffset + idxSir], exSigInterp));
			}
		}
	}

	return;
}

// main function used to build our transducer model
void transModel::buildGPUModel()
{
	cudaError_t err; // error of cuda execution
	time_t timer; // timer for execution timing
	
	// copy all properties over into struct
	fieldStruct field;
	field.fs = fieldProp.get_fs();
	field.dt = fieldProp.get_dt();
	field.dr = fieldProp.get_dr();
	field.dz = fieldProp.get_dz();
	field.SOS = fieldProp.get_SOS();
	field.nr = fieldProp.get_nr();
	field.nz = fieldProp.get_nz();
	field.nElements = fieldProp.get_nElements();
	field.z0 = fieldProp.get_z0();
	field.dSir = fieldProp.get_dSir();
	// sir0 and t0 are calculated later

	model.setNrModel(field.nr);
	model.setNzModel(field.nz);
	model.setTRes(field.dt);

	// get excitation signal for convolution from matlab
	const float* exSig = getImpResp();
	const unsigned int nExSig = getNElem(); 
	if (nExSig == 1)
		printf("Warning: length of excitation signal is 1");

	subElement* transducerElements_dev;
	transducerProperties* trans_dev;

	discretizeTransducer(); // performing transducer discretization

	if (flagVerbose)
		printf("Number of transducer elements: %d \n", nElements);

	// get transducer properties from matlab
	transducerProperties trans = prop;

	// allocate memory for transducer properties
	err = cudaMalloc( (void**) &trans_dev, sizeof(transducerProperties));
	checkCudaReturn(err, "Could not allocate memory for transducer properties on device");
	
	// allocate memory for transducer elements
	err =  cudaMalloc( (void**) &transducerElements_dev, nElements * sizeof(subElement));
	checkCudaReturn(err, "Could not allocate memory for transducer subelements on device");
 	
 	// copy transducer properties over to device
	err = cudaMemcpy(trans_dev, &trans,	sizeof(transducerProperties), 
		cudaMemcpyHostToDevice);
	checkCudaReturn(err, "Could not copy transducer properties to device");

	// copy subelements to device
	err = cudaMemcpy(transducerElements_dev, transducerElements,
		nElements * sizeof(subElement), cudaMemcpyHostToDevice);
	checkCudaReturn(err, "Could not copy subelements to device");

	// get maximum and minimum distance of elements to create SIR matrix
	// vprintf("Finding maximum and minimum in SIR matrix... \n");
	
	field.sir0 = 1e5;
	float maxDist = 0;
	float rDist, deltaX, deltaY, deltaZ;
	// note: we only need to check the boundaries of the field and not the center
	#pragma unroll
	for (unsigned int iz = 0; iz < field.nz; iz = iz + (field.nz - 1)){
		#pragma unroll
		for (unsigned int ir = 0; ir < field.nr; ir++){
			#pragma unroll
			for (unsigned int iElem = 0; iElem < nElements; iElem++){
				deltaY = transducerElements[iElem].centerPos.y;
				deltaX = transducerElements[iElem].centerPos.x - ((float) ir * field.dr);
				deltaZ = transducerElements[iElem].centerPos.z - (field.z0 + field.dz * 
					((float) iz));
				
				// calculate absolute distance between field point and element
				rDist = sqrt(deltaX * deltaX + deltaY * deltaY + deltaZ * deltaZ);

				if (rDist > maxDist)
					maxDist = rDist;

				if (rDist < field.sir0)
					field.sir0 = rDist;
			}
		}
	}
	field.t0 = field.sir0 / field.SOS;
	if (flagVerbose) 
		printf("Min/max distance in field %f / %f \n", field.sir0, maxDist);

	model.setNzSir((maxDist - field.sir0) / field.dSir + 1);
	
	float* sir_dev; // poitner to memory of sir on device
	fieldStruct* field_dev;
	bool m2 = (cudaSuccess != cudaMalloc( (void**) &sir_dev, 
		model.getNElementsSir() * sizeof(float)));

	// generate field properties on GPU and copy over
	m2 |= (cudaSuccess != cudaMalloc( (void**) &field_dev, sizeof(fieldStruct)));
	m2 |= (cudaSuccess != cudaMemcpy(field_dev, &field, sizeof(fieldStruct), 
		cudaMemcpyHostToDevice)); 	
	

	if (m2){
		printf("Could not allocate memory for field or spatial impulse resposne\n");
		printf("Requested dimensions of SIR: nSir = %d, nr = %d, nz = %d\n", 
			model.getNzSir(), field.nr, field.nz);
		printf("Required memory for this: %f Gb\n", (float) model.getNElementsSir() * 4 / 1024 / 1024 / 1024);
		throw "CudaMemAlloc";
	}else{
		// define thread and grid size for 2D spatial impulse response calculation
		dim3 blockSizeSir(4, 4);
		dim3 gridSizeSir((field.nr + blockSizeSir.x - 1) / blockSizeSir.x,
				(field.nz + blockSizeSir.y - 1) / blockSizeSir.y);

		// calculate spatial impulse response
		vprintf("Calculating spatial impulse response... ");

		constGetSirArgs constArgs;
		constArgs.trans = transducerElements_dev;
		constArgs.field = field_dev;
		constArgs.nElem = nElements;
		constArgs.nzSir = model.getNzSir();

		constGetSirArgs* constArgs_dev;
		err = cudaMalloc( (void**) &constArgs_dev, sizeof(constGetSirArgs));
		checkCudaReturn(err, "Could not allocate memory for constant arguments");
		
		err = cudaMemcpy(constArgs_dev, &constArgs, sizeof(constGetSirArgs), 
			cudaMemcpyHostToDevice);
		checkCudaReturn(err, "Could copy constant arguments");

		timer = time(NULL); // start stopwatch here

		getSIR<<< gridSizeSir, blockSizeSir >>>(
				sir_dev, // spatial impulse response [iz/it, irfield, izfield] 
				constArgs_dev); // number of elements along z for sir

		cudaDeviceSynchronize();

		err = cudaGetLastError();
		if (err != cudaSuccess){
			string errTxt = cudaGetErrorString(err);
			printf(errTxt.c_str());
			printf("Failed to get spatial impulse response\n");
			throw "cudaError";
		}else{
			double seconds = difftime(time(NULL), timer);
			if (flagVerbose)
				printf("done after %.0f s!\n", seconds);

			// allocate memory for output matrix
			model.allocate_sir();

			// copy matrix back from GPU into MATLAB workspace
			err = cudaMemcpy(model.sir, sir_dev, 
				model.getNElementsSir() * sizeof(float), cudaMemcpyDeviceToHost);
			if (err != cudaSuccess){
				string errTxt = cudaGetErrorString(err); 
				printf(errTxt.c_str());
				printf("Failed to copy spatial impulse response back\n");
				throw "cudaError";
			}else{
				// allocate memory for model matrix
				float* model_dev;	// structure of model matrix [it, ir, iz]
				// printf("range sir = %f\n ", maxDist - field.sir0);
				// printf("t range sir = %f\n", (maxDist - field.sir0) / fieldProp.get_SOS() * fieldProp.get_fs());
				uint64_t nTModel = (maxDist - field.sir0) / fieldProp.get_SOS() * fieldProp.get_fs() 
				 	+ nExSig - 1;
				
				if (flagDebug)
				{
					printf("Number of time elements in model: %d\n", nTModel);
					printf("maxDist: %f\n", maxDist);
					printf("sir0: %f\n", field.sir0);
					printf("SOS: %f\n", fieldProp.get_SOS());
					printf("Sampling frequency: %f\n", fieldProp.get_fs());
				}

				model.setNtModel(nTModel);

				// allocate memory for model matrix of GPU
				err = cudaMalloc( (void**) &model_dev, model.get_nt() * field.nr * field.nz *
					sizeof(float) );
				checkCudaReturn(err, "Could not allocate mdoel matrix");
				
				float* exSig_dev;
				err = cudaMalloc( (void**) &exSig_dev, nExSig * sizeof(float) );
				checkCudaReturn(err, "Could not allocate excitation signal");

				err = cudaMemcpy(exSig_dev, exSig, nExSig * sizeof(float), cudaMemcpyHostToDevice);
				checkCudaReturn(err, "Could not copy excitation signal");

				// define 3D block and grid size for time domain convolution
				dim3 blockSizeConv(2, 2, 4);
				dim3 gridSizeConv(
						(model.get_nt() + blockSizeConv.x - 1) / blockSizeConv.x, 
						(field.nr + blockSizeConv.y - 1) / blockSizeConv.y,
						(field.nz + blockSizeConv.z - 1) / blockSizeConv.z);

				// perform convolution
				vprintf("Convolving SIR with impulse response... ");
				
				timer = time(NULL); // start stopwatch here

				convolveSIR<<< gridSizeConv, blockSizeConv>>>(
						model_dev, // model matrix [it, ir, iz]
						sir_dev, // spatial impulse response [izSir, ir, iz]
						exSig_dev, // excitation signal vector
						field_dev, // properties of output field
						model.getNzSir(), // number of z elements in sir 
						nExSig, // length of exSig_dev
						model.get_nt(), // number of timepoints in model matrix
						tOffset); // temporal offset in excitation signal [s]	

				cudaDeviceSynchronize();
				err = cudaGetLastError();
				if (err != cudaSuccess){
					throw cudaGetErrorString(err);	
				}else{

					seconds = difftime(time(NULL), timer);
					if (flagVerbose)
						printf("done after %.0f s!\n", seconds);
				
					model.allocate(); // allocate memory for model matrix

					// save start time
					starttime = field.sir0 / field.SOS;
					model.setT0(starttime);

					// copy model matrix back from GPU
					err = cudaMemcpy(model.data, model_dev, 
						model.getNElements() * sizeof(float), cudaMemcpyDeviceToHost);
					checkCudaReturn(err, "Could not copy model back from device");					
				}
				cudaFree(model_dev); // free memory for model matrix on device
				cudaFree(exSig_dev); // free memory for excitation signal on device
			}
		}
	}

	model.setRRes(fieldProp.get_dr());
	model.setZRes(fieldProp.get_dz());
	model.setZ0(fieldProp.get_z0());

	cudaFree(sir_dev);
	cudaFree(field_dev);
	cudaFree(transducerElements_dev);
	cudaFree(trans_dev);
	
	return;
}

modelMatrix* transModel::getModel()
{
	return &model;
}
