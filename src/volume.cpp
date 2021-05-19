#include "volume.h"

volume::volume()
{
	className = "volume";
}

volume::~volume()
{

	// if memory was allocated in data, release now
	if (isMemAlloc)
	{
		delete[] data;
		delete[] sliceZ;
		delete[] sliceX;
		delete[] sliceY;
	}
}

float* volume::get_psliceZ(const uint64_t zLevel)
{
	if (zLevel != lastSliceZ)
	{
		lastSliceZ = zLevel;
		// update slice
		for (uint64_t ix = 0; ix < dim[1]; ix++)
		{
			for (uint64_t iy = 0; iy < dim[2]; iy++)
			{
				const uint64_t sliceIdx = ix + iy * dim[1];
				const uint64_t volIdx = zLevel + ix * dim[0] + iy
					* dim[0] * dim[1];
				sliceZ[sliceIdx] = data[volIdx];
			}
		}
	}

	return sliceZ;
}

float* volume::get_psliceX(const uint64_t xLevel)
{
	if (xLevel != lastSliceX)
	{
		lastSliceX = xLevel;
		// update slice
		for (uint64_t iz = 0; iz < dim[0]; iz++)
		{
			for (uint64_t iy = 0; iy < dim[2]; iy++)
			{
				const uint64_t sliceIdx = iz + iy * dim[0];
				const uint64_t volIdx = iz + xLevel * dim[0] + iy
					* dim[0] * dim[1];
				sliceX[sliceIdx] = data[volIdx];
			}
		}
	}

	return sliceX;
}

float* volume::get_psliceY(const uint64_t yLevel)
{
	if (yLevel != lastSliceY)
	{
		lastSliceY = yLevel;
		// update slice
		for (uint64_t iz = 0; iz < dim[0]; iz++)
		{
			for (uint64_t ix = 0; ix < dim[1]; ix++)
			{
				const uint64_t sliceIdx = ix + iz * dim[1];
				const uint64_t volIdx = iz + ix * dim[0] + yLevel
					* dim[0] * dim[1];
				sliceY[sliceIdx] = data[volIdx];
			}
		}
	}

	return sliceY;
}

// get slice of volume at position
float* volume::get_psliceZ(const float zPos)
{
	uint64_t zIdx = getIdx(zPos, 0);
	return get_psliceZ(zIdx);
}

float* volume::get_psliceX(const float xPos)
{
	uint64_t xIdx = getIdx(xPos, 1);
	return get_psliceX(xIdx);
}

float* volume::get_psliceY(const float yPos)
{
	uint64_t yIdx = getIdx(yPos, 2);
	return get_psliceY(yIdx);
}

float volume::getLength(const uint8_t _dim)
{
	return (float) dim[_dim] * res[_dim];
}

void volume::allocMemory()
{
	if (isMemAlloc)
	{
		delete[] data;
		delete[] sliceX;
		delete[] sliceY;
		delete[] sliceZ;
	}

	isMemAlloc = 1;
	unsigned int nElements = dim[0] * dim[1] * dim[2];
	data = new float[nElements];
	sliceZ = new float [dim[1] * dim[2]];
	sliceX = new float [dim[0] * dim[2]];
	sliceY = new float [dim[0] * dim[1]];
	return;
}

// scales whole data array by a number
void volume::multiply(const float factor)
{
	unsigned int nElements = dim[0] * dim[1] * dim[2];
	for (unsigned int iElement = 0; iElement < nElements; iElement++)
		data[iElement] *= factor;

	return;
}

// define dimensions of dataset
void volume::setDim(const uint64_t dim0, const uint64_t dim1, const uint64_t dim2)
{
	dim[0] = dim0;
	dim[1] = dim1;
	dim[2] = dim2;
	nElements = dim[0] * dim[1] * dim[2];
	return;
}

void volume::setDim(const uint64_t* _dim)
{
	dim[0] = _dim[0];
	dim[1] = _dim[1];
	dim[2] = _dim[2];
	nElements = dim[0] * dim[1] * dim[2];
	return;
}

void volume::setDim(const uint8_t _dim, const uint64_t newDim)
{
	dim[_dim] = newDim;
	nElements = dim[0] * dim[1] * dim[2];
	return;
}

uint64_t volume::getDim(const uint8_t _dim) const
{
	return dim[_dim];
}

// set origin of volumetrix dataset
void volume::setOrigin(const float* _origin)
{
	origin[0] = _origin[0];
	origin[1] = _origin[1];
	origin[2] = _origin[2];
	return;
}

void volume::setOrigin(const float origin0, const float origin1, const float origin2)
{
	origin[0] = origin0;
	origin[1] = origin1;
	origin[2] = origin2;
	return;
}

void volume::setOrigin(const uint8_t _dim, const float _origin)
{
	origin[_dim] = _origin;
	return;
}

float volume::getOrigin(const uint8_t _dim)
{
	return origin[_dim];
}

// set resolution of volumetric dataset
void volume::setRes(const float* dx)
{
	// res[0] = dx[0];
	// res[1] = dx[1];
	// res[2] = dx[2];
	for (uint8_t iDim = 0; iDim < 3; iDim++)
		setRes(iDim, dx[iDim]);
	return; 
}

void volume::setRes(const float dx0, const float dx1, const float dx2)
{
	// res[0] = dx0;
	// res[1] = dx1;
	// res[2] = dx2;
	setRes(0, dx0);
	setRes(1, dx1);
	setRes(2, dx2);
	return;
}

void volume::setRes(const uint8_t _dim, const float _res)
{
	if (_res <= 0)
	{
		printf("Resolution along %d axis needs to be bigger then 0\n");
		throw "invalidValue";
	}
	res[_dim] = _res;
	return;
}

float volume::getRes(const uint8_t _dim) const {return res[_dim];}
float volume::getRes0() const {return res[0];}
float volume::getRes1() const {return res[1];}
float volume::getRes2() const {return res[2];}

// sets whole array to a certain value
void volume::setValue(const float value)
{
	unsigned int nElements = dim[0] * dim[1] * dim[2];
	for (unsigned int iElement = 0; iElement < nElements; iElement++)
		data[iElement] = value;

	return;
}

// set only one specific value in volume defined by index
void volume::setValue(
	const unsigned int x0, const unsigned int x1, const unsigned int x2, const float value)
{
	unsigned int index = x0 + dim[0] * (x1 + x2 * dim[1]);
	data[index] = value;
	return;
}

// set only one specific value in volume
void volume::setValue(const uint64_t* pos, const float value)
{
	uint64_t index = pos[0] + dim[0] * (pos[1] + pos[2] * dim[1]);
	data[index] = value;
	return;
}

float volume::getValue(const unsigned int x0, const unsigned int x1, const unsigned int x2) const
{
	const unsigned int idx = x0 + dim[0] * (x1 + x2 * dim[2]);
	return data[idx];
}

float volume::getValue(const unsigned int iElem) const {return data[iElem];}

float volume::getValue(const unsigned int* pos) const
{
	const unsigned int idx = pos[0] + dim[0] * (pos[1] + pos[2] * dim[2]);
	return data[idx];
}

// get position along axis
float volume::getPos0(const uint64_t idx0) const
{
	return origin[0] + (float) idx0 * res[0];
}

float volume::getPos1(const uint64_t idx1) const
{
	return origin[1] + (float) idx1 * res[1];
}

float volume::getPos2(const uint64_t idx2) const
{
	return origin[2] + (float) idx2 * res[2];
}

float volume::getPos(const uint64_t idx, const uint8_t iDim) const
{
	return origin[iDim] + (float) idx * res[iDim];
}

float volume::getCenterPos(const uint8_t _dim)
{
	const float centerPos = origin[_dim] + ((float) dim[_dim] * res[_dim]) / 2;
	return centerPos;
}

// get index along a certain dimension
uint64_t volume::getIdx0(const float pos0) const {return getIdx(pos0, 0);}
uint64_t volume::getIdx1(const float pos1) const {return getIdx(pos1, 1);}
uint64_t volume::getIdx2(const float pos2) const {return getIdx(pos2, 2);}

uint64_t volume::getIdx(const float pos, const uint8_t iDim) const
{
	if (pos < origin[iDim]){
		printf("[volume] hitting lower boundary of volume along dim %d\n", iDim);
		return 0;
	}
	else{
		const unsigned int outputIdx = (pos - origin[iDim]) / res[iDim] + 0.5;
		if (outputIdx > (dim[iDim] - 1)){
			printf("[volume] hitting upper boundary of volume along dim %d\n", iDim);
			return dim[iDim] - 1;
		}
		else
			return outputIdx;
	}
}

uint64_t volume::get_nElements() const
{
	uint64_t nElements = dim[0] * dim[1] * dim[2];
	return nElements;
}

// not implemented yet
void volume::saveToFile(const string filePath) const
{
	vprintf("Saving data to file", 1);
	H5::H5File file(filePath, H5F_ACC_TRUNC);

	// write resolutiion to file
	//printf("Write resolution to file...\n");
	const hsize_t col_dims = 3;
	H5::DataSpace mspaceRes(1, &col_dims);
	H5::DataSet resDataset = file.createDataSet(
		"dr", H5::PredType::NATIVE_FLOAT, mspaceRes);
	resDataset.write(res, H5::PredType::NATIVE_FLOAT);
	resDataset.close();

	// write origin to file
	//printf("Write origin to file...\n");
	H5::DataSpace mspaceOrigin(1, &col_dims);
	H5::DataSet originDataset = file.createDataSet(
		"origin", H5::PredType::NATIVE_FLOAT, mspaceOrigin);
	originDataset.write(origin, H5::PredType::NATIVE_FLOAT);
	originDataset.close();

	// write dimension to file
	//printf("Write dimensions to file...\n");
	H5::DataSpace mspaceDim(1, &col_dims);
	H5::DataSet dimDataset = file.createDataSet(
		"dim", H5::PredType::NATIVE_UINT, mspaceDim);
	dimDataset.write(dim, H5::PredType::NATIVE_UINT64);
	dimDataset.close();

	// w4rite actual datamatrix to file
	const hsize_t col_data = nElements;
	H5::DataSpace mspaceData(1, &col_data);
	H5::DataSet dataDataset = file.createDataSet(
		"vol", H5::PredType::NATIVE_FLOAT, mspaceData);
	dataDataset.write(data, H5::PredType::NATIVE_FLOAT);
	dataDataset.close();

	file.close();
	vprintf("done!\n", 0);
	return;
}

void volume::readFromFile(const string filePath)
{
	vprintf("Reading from file... ", 1);
	H5::H5File file(filePath, H5F_ACC_RDONLY); // open dataset as read only

	// load resolution from file
	// printf("Loading resolution from file...\n");
	H5::DataSet resDataset = file.openDataSet("dr"); // init dataset for res 
	const hsize_t col_dims = 3;
	H5::DataSpace mspaceRes (1, &col_dims); 
	H5::DataSpace filespace = resDataset.getSpace();
	resDataset.read(res, H5::PredType::NATIVE_FLOAT, mspaceRes, filespace); 	

	// load origin from file
	// printf("Loading origin from file...\n");
	H5::DataSet originDataset = file.openDataSet("origin"); // dataset for origin
	H5::DataSpace mspaceOrigin (1, &col_dims); 
	filespace = originDataset.getSpace();
	originDataset.read(origin, H5::PredType::NATIVE_FLOAT, mspaceOrigin, filespace);

	// load number of dimensions from file
	// printf("Loading dimnesions from file...\n");
	H5::DataSet dimDataset = file.openDataSet("dim"); // init dimension reader
	H5::DataSpace mspaceDim (1, &col_dims);
	filespace = dimDataset.getSpace();
	dimDataset.read(dim, H5::PredType::NATIVE_UINT64, mspaceDim, filespace);
	nElements = dim[0] * dim[1] * dim[2];

	// read actual datamatrix
	// printf("Loading data from file...\n");
	H5::DataSet dataDataset = file.openDataSet("vol"); // read actual dataset
	const hsize_t col_dims_data = nElements; 
	H5::DataSpace mspaceData (1, &col_dims_data);	
	filespace = dataDataset.getSpace();
	allocMemory();
	dataDataset.read(data, H5::PredType::NATIVE_FLOAT, mspaceData, filespace);

	isMemAlloc = 1;

	file.close();
	vprintf("done!", 0);
	return;
}

// prints all the required information about the loaded volume
void volume::printInformation() const
{
	printf("Volumetric dataset properties: \n");
	printf(" - origin:        %.2e, %.2e, %.2e \n", origin[0], origin[1], origin[2]);
	printf(" - resolution:    %.2e, %.2e, %.2e \n", res[0], res[1], res[2]);
	printf(" - dimensions:    %d, %d, %d \n", dim[0], dim[1], dim[2]);
	printf(" - first element: %f \n - last element:  %f\n", data[0], data[nElements - 1]);
	return;
}

float volume::getMinPos(const unsigned int _dim) const {return origin[_dim];}
float volume::getMaxPos(const unsigned int _dim) const
{
	return origin[_dim] + res[_dim] * ((float) dim[_dim] - 1.0);
}

float volume::getRangeLimitedPos(const float pos, const unsigned int _dim) const
{
	float posOut = pos;
	if (posOut < getMinPos(_dim))
		posOut = getMinPos(_dim);

	if (posOut > getMaxPos(_dim))
		posOut = getMaxPos(_dim);

	return posOut;
}

// get cropped volume, start and stopped passed as array
void volume::getCroppedVolume(
	float* vol, const uint64_t *startIdx, const uint64_t *stopIdx) const
{
 
	const uint64_t start0 = startIdx[0];
	const uint64_t stop0 = stopIdx[0];
	const uint64_t start1 = startIdx[1];
	const uint64_t stop1 = stopIdx[1];
	const uint64_t start2 = startIdx[2];
	const uint64_t stop2 = stopIdx[2];

 getCroppedVolume(vol, start0, stop0, start1, stop1, start2, stop2);
 return;
}

void volume::getCroppedVolume(
	float* vol, // array pointing to volume
	const uint64_t start0, const uint64_t stop0,
	const uint64_t start1, const uint64_t stop1,
	const uint64_t start2, const uint64_t stop2) const
{

	// check that stop index is always bigger then start index
	if (start0 >= stop0)
	{
		printf("Stop index must be bigger then start index along dim0\n");
		throw "invalidVal";
	}

	if (start1 >= stop1)
	{
		printf("Stop index must be bigger then start index along dim1\n");
		throw "invalidVal";
	}

	if (start2 >= stop2)
	{
		printf("Stop index must be bigger then start index along dim2\n");
		throw "invalidVal";
	}

	// calc number of elements in each dimension
	const uint64_t n0 = stop0 - start0 + 1; // n along t/z
	const uint64_t n1 = stop1 - start1 + 1; // n alonmg x
	const uint64_t n2 = stop2 - start1 + 1; // n along y
	const uint64_t nElements = n0 * n1 * n2;

	// set all values in volume for now to 0
	for (uint64_t idx = 0; idx < nElements; idx++)
		vol[idx] = 0;

	// check if all stop indices are within range
	const uint64_t stopTrue0 = (stop0 >= dim[0]) ? (dim[0] - 1) : stop0;
	if (stopTrue0 != stop0)
		printf("Cropping not doable in dim0\n");

	const uint64_t stopTrue1 = (stop1 >= dim[1]) ? (dim[1] - 1) : stop1;
	if (stopTrue1 != stop1)
		printf("Cropping not doable in dim1\n");
	
	const uint64_t stopTrue2 = (stop2 >= dim[2]) ? (dim[2] - 1) : stop2;
	if (stopTrue2 != stop2)
		printf("Cropping not doable in dim2\n");

	// index in original data matrix
	// idx = i0 + n0 * i1 + n0 * n1 * i2
	uint64_t idx; // index in original data matrix
	uint64_t offset2, offset1; // offsets in original data matrix
	
	// index in nnewly created data matrix
	// idxOut = j0 + m0 * j1 + m0 * m1 * j2

	// original --> new
	// i0 --> i0 - start0
	// i1 --> i1 - start1
	// i2 --> i2 - start2
	
	uint64_t idxOut;
	uint64_t offset2Out, offset1Out; // offset in output data matrix
	#pragma unroll
	for (uint64_t i2 = start2; i2 <= stopTrue2; i2++){
		offset2 = i2 * dim[0] * dim[1];
		offset2Out = (i2 - start2) * n0 * n1;
		#pragma unroll
		for (uint64_t i1 = start1; i1 <= stopTrue1; i1++){
			offset1 = i1 * dim[0];
			offset1Out = (i1 - start1) * n0;
			#pragma unroll
			for (unsigned int i0 = start0; i0 < stopTrue0; i0++){
				idx = offset2 + offset1 + i0;
				idxOut = offset2Out + offset1Out + (i0 - start0);
				vol[idxOut] = data[idx];
			}
		}
	}

	return;
}

// calculates maximum and minimum value in matrix
void volume::calcMinMax()
{
	minVal = getMin(data, nElements);
	maxVal = getMax(data, nElements);

	if (abs(minVal) > abs(maxVal))
	{
		maxAbsVal = abs(minVal);
	}
	else
	{
		maxAbsVal = abs(maxVal);
	}
	return;
}

void volume::exportVtk(const string filePath)
{
	// create copy of datset before messing with polarity
	float* dataCopy = new float [nElements];
	assign(dataCopy, data, nElements);
	
	// const string polarityHandling (sett.get_polarityHandling());
	
	// handlePolarity(dataCopy, estimAbs.nElements, polarityHandling);

	vtkwriter outputter; // prepare output pipeline
	const string title ("reconVol"); // generate title
	outputter.set_title(title); // define title of outut volume
	const string type ("STRUCTURED_POINTS"); 
	outputter.set_type(type); // define output type
	const string outputPath (filePath);
	outputter.set_outputPath(outputPath);

	float * dataOld = data;
	data = dataCopy;

	griddedData myData;
	myData.data = data;
	for (uint8_t iDim = 0; iDim < 3; iDim++)
	{
		myData.origin[iDim] = origin[iDim];
		myData.res[iDim] = res[iDim];
		myData.dim[iDim] = dim[iDim];
	}

	outputter.set_structuredPoints(&myData);
	outputter.set_binary();
	outputter.write();
	delete[] dataCopy;
	data = dataOld;
	return;
}

void volume::set_pdata(float* _data)
{
	data = _data;
	return;
}