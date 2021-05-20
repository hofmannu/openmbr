#include "modelMatrix.h"

modelMatrix::modelMatrix()
{
	className = "modelMatrix";
}

modelMatrix::~modelMatrix()
{
	// free memory if required
	if (isModelMatrixAlloc){
		delete[] data;
		delete[] dataPermuted;
	}

	if (isStartStopIdxAlloc)
		delete[] startStopIdx;

	if (isStartStopIdxPermutedAlloc)
		delete[] startStopIdxPermuted;

	if (isSensFieldAlloc)
		delete[] sensField;

	if (isModSliceAlloc)
		delete[] modSlice;
}

// set functions
void modelMatrix::setNtModel(uint64_t _ntModel)
{
	if (_ntModel == 0){
		printf("Model matrix must have size bigger then 0 along t\n");
		throw "sizeErr";
	}else{
		nt = _ntModel;
		nElements = nt * nr * nz;
	}
	return;
}

void modelMatrix::setNzModel(uint64_t _nzModel)
{
	if (_nzModel == 0){
		printf("Model matrix must have size bigger then 0 along z\n");
		throw "sizeErr";
	}else{
		nz = _nzModel;
		nElements = nt * nr * nz;
		nElementsSir = nzSir * nr * nz;
	}
	return;
}

void modelMatrix::setNrModel(uint64_t _nrModel)
{
	if (_nrModel == 0){
		printf("Model matrix must have size bigger then 0 along t\n");
		throw "sizeErr";
	}else{
		nr = _nrModel;
		nElements = nt * nr * nz;
		nElementsSir = nzSir * nr * nz;
	}
	return;
}

void modelMatrix::setNzSir(uint64_t _nzSir)
{
	if (_nzSir == 0){
		printf("Sir matrix must have size bigger then 0 along 1st dim\n");
		throw "sizeErr";
	}else{
		nzSir = _nzSir;
		nElementsSir = nzSir * nr * nz;
	}
	return;
}

void modelMatrix::setRRes(const float _rRes)
{
	dr = _rRes;
	return;
}

void modelMatrix::setZRes(const float _zRes)
{
	dz = _zRes;
	return;
}

void modelMatrix::setZ0(const float _z0)
{
	z0 = _z0;
	return;
}

void modelMatrix::setT0(const float _t0)
{
	t0 = _t0;
	return;
}

void modelMatrix::setTRes(const float _tRes)
{
	tRes = _tRes;
	return;
}

void modelMatrix::setSir0(const float _sir0)
{
	sir0 = _sir0;
	return;
}

void modelMatrix::setSirRes(const float _sirRes)
{
	sirRes = _sirRes;
	return;
}

void modelMatrix::allocate()
{
	if (isModelMatrixAlloc)
	{
		delete[] data;
		delete[] dataPermuted;
	}

	data = new float [nElements];
	dataPermuted = new float [nElements];

	if (isModSliceAlloc)
		delete[] modSlice;

	modSlice = new float [nr * nt];
	isModSliceAlloc = 1;

	isModelMatrixAlloc = 1;
	return;
}

void modelMatrix::allocate_sir()
{
	if (isSirAlloc)
	{
		delete[] sir;
	}

	sir = new float [nElementsSir];
	isSirAlloc = 1;

	return;
}

// generates the dataPermuted [iz, ir, it] out of data [it, ir, iz]
void modelMatrix::permute()
{
	uint64_t idx, idxPermuted;
	#pragma unroll
	for (uint64_t iz = 0; iz < nz; iz++)
	{
		#pragma unroll
		for (uint64_t  ir = 0; ir < nr; ir++)
		{
			#pragma unroll
			for (uint64_t it = 0; it < nt; it++)
			{
				idxPermuted = iz + nz * (ir + nr * it);
				idx = it + nt * (ir + nr * iz);
				dataPermuted[idxPermuted] = data[idx];
			}
		}
	}
	return;
}

void modelMatrix::getStartStopIdx()
{
	getStartStopIdxPermuted(); // for forward kernel (absorber to signal)
	getStartStopIdxNonPermuted(); // for transpose kernel (signal to absorber)
	return;
}


void modelMatrix::getStartStopIdxNonPermuted()
{
	// model matrix indexing: iz, ir, it
	// iz + nzModel * (ir + it * nrModel)

	// startStopIndex mapping
	// ir + nr * it + (1/0)
	if (isStartStopIdxAlloc)
		delete[] startStopIdx;

	startStopIdx = new uint64_t [nr * nz * 2];
	isStartStopIdxAlloc = 1;

	bool foundStart; // bool representing if first value above threshold found
		bool foundGlobalPeak;
		int mmIdx; // current index in model matrix
		int mmStIdx; // current index in model start/stop index matrix
		float maxAbsVal;
		float globalPeak;
		int counterEmptyModel = 0;

		#pragma unroll	
		for (uint64_t iz = 0; iz < nz; iz++){
			// determine absolute maximum value in model matrix
			globalPeak = 0;
			#pragma unroll
			for (uint64_t ir = 0; ir < nr; ir++){
				#pragma unroll
				for (uint64_t it = 0; it < nt; it++){
					mmIdx = it + nt * (ir + iz * nr);	
					if (fabs(data[mmIdx]) > globalPeak)
						globalPeak = fabs(data[mmIdx]);
				}
			}	
			globalPeak *= threshold;

			#pragma unroll
			for (uint64_t ir = 0; ir < nr; ir++){
				// find max value of currently watched model vector
				maxAbsVal = 0;
				#pragma unroll
				for (uint64_t it = 0; it < nt; it++){
					mmIdx = it + nt * (ir + iz * nr);
					if (fabs(data[mmIdx]) > maxAbsVal)
						maxAbsVal = fabs(data[mmIdx]);
				}
				
				// determine threshold based on maximum found value
				maxAbsVal *= threshold;
				foundStart = 0; // reset found start flag
				foundGlobalPeak = 0;
				mmStIdx = 2 * (ir + iz * nr);
				#pragma unroll
				for (uint64_t it = 0; it < nt; it++){
					mmIdx = it + nt * (ir + iz * nr);	
					if ((fabs(data[mmIdx]) > maxAbsVal) && (!foundStart)){
						foundStart = 1;
						startStopIdx[mmStIdx] = it;
					}
					if (fabs(data[mmIdx]) > maxAbsVal)
						startStopIdx[mmStIdx + 1] = it;

					if (fabs(data[mmIdx]) > globalPeak)
						foundGlobalPeak = 1;
				}

				// check if there was no starting point at all, that means no need to
				// reconstruct anything --> put start and stop to 0
				if (!foundStart){
					startStopIdx[mmStIdx] = 0;
					startStopIdx[mmStIdx + 1] = 0;
					counterEmptyModel++;
				}else if(!foundGlobalPeak){ // if whole model vector is below plane thres
					startStopIdx[mmStIdx] = 0;
					startStopIdx[mmStIdx + 1] = 0;
					counterEmptyModel++;
				}
			}
		}
	float precEmpty = (float) counterEmptyModel / ((float) nr * nz);
	printf("Percent empty vectors: %.1f\n", precEmpty * 100);
	return;
}

void modelMatrix::getStartStopIdxPermuted()
{
	// data order [it, ir, iz]
	// find for forward calculation
	if (isStartStopIdxPermutedAlloc)
		delete[] startStopIdxPermuted;

	startStopIdxPermuted = new uint64_t [nt * nr * 2];
	isStartStopIdxPermutedAlloc = 1;

	bool foundStart; // bool representing if first value above threshold found
		bool foundGlobalPeak;
		int mmIdx; // current index in model matrix
		int mmStIdx; // current index in model start/stop index matrix
		float maxAbsVal; // maximum absolute value of current model vector
		float globalPeak; // maximum absolute value of current time plane
		int counterEmptyModel = 0;	
		#pragma unroll
		for (uint64_t it = 0; it < nt; it++){
			// determine absolute maximum value in model matrix plane
			globalPeak = 0;
			#pragma unroll
			for (uint64_t ir = 0; ir < nr; ir++){
				#pragma unroll
				for (uint64_t iz = 0; iz < nz; iz++){
					mmIdx = iz + nz * (ir + it * nr);
					if(fabs(dataPermuted[mmIdx]) > globalPeak)
						globalPeak = fabs(dataPermuted[mmIdx]);
				}
			}
			globalPeak *= threshold;	

			#pragma unroll
			for (uint64_t ir = 0; ir < nr; ir++){
				// find max value of currently watched model vector
				maxAbsVal = 0;
				#pragma unroll
				for (uint64_t iz = 0; iz < nz; iz++){
					mmIdx = iz + nz * (ir + it * nr);
					if (fabs(dataPermuted[mmIdx]) > maxAbsVal)
						maxAbsVal = fabs(dataPermuted[mmIdx]);
				}
				
				// determine threshold based on maximum found value
				maxAbsVal *= threshold;
				foundStart = 0; // reset found start flag
				foundGlobalPeak = 0;
				mmStIdx = 2 * (ir + it * nr);
				#pragma unroll
				for (uint64_t iz = 0; iz < nz; iz++){
					mmIdx = iz + nz * (ir + it * nr);	
					if ((fabs(dataPermuted[mmIdx]) > maxAbsVal) && (!foundStart)){
						foundStart = 1;
						startStopIdxPermuted[mmStIdx] = iz;
					}
					if (fabs(dataPermuted[mmIdx]) > maxAbsVal)
						startStopIdxPermuted[mmStIdx + 1] = iz;

					// check if element is above global peak * threshold
					if (fabs(dataPermuted[mmIdx]) > globalPeak)
						foundGlobalPeak = 1;
				}

				// if we did not fine a start at all
				if (!foundStart){
					startStopIdxPermuted[mmStIdx] = 0;
					startStopIdxPermuted[mmStIdx + 1] = 0;
					counterEmptyModel++;
				}else if (!foundGlobalPeak){ // if whole vector is below global threshold
					startStopIdxPermuted[mmStIdx] = 0;
					startStopIdxPermuted[mmStIdx + 1] = 0;
					counterEmptyModel++;
				}
			}
		}
	float precEmpty = (float) counterEmptyModel / ((float) nr * nt);
	printf("Percent empty vectors: %.1f\n", precEmpty * 100);
	return;
}

void modelMatrix::exportSensFieldVtk() const
{
	const string filePath = "~/sensField.vtk";
	exportSensFieldVtk(filePath);
	return;
}

// exports the sensitivity field as a 3D matrix for paraview
void modelMatrix::exportSensFieldVtk(const string filePath) const
{

	// calculate maximum along time dimension of each signal and store in mip
	float* mip = new float [nz * nr];
	for (unsigned int iz = 0; iz < nz; iz ++)
	{
		for (unsigned int ir = 0; ir < nr; ir++)
		{
			float maxVal = 0;
			for (unsigned int it = 0; it < nt; it++)
			{
				float currVal = data[it + nt * (ir + nr * iz)];
				if (maxVal < fabs(currVal))
					maxVal = fabs(currVal);
			}
			mip[iz + ir * nz] = maxVal;
		}
	}

	// generate volumetric representation of the field and assign
	griddedData exportVol;
	unsigned int nXY = 2 * nr - 1;
	exportVol.data = new float [nXY * nXY * nz];
	unsigned int ir;
	int xPos, yPos;
	#pragma unroll
	for (unsigned int ix = 0; ix < nXY; ix++)
	{
		xPos = (int) ix - nr;
		#pragma unroll
		for (unsigned int iy = 0; iy < nXY; iy++)
		{
			yPos = (int) iy - nr;
			ir = sqrt((float) (xPos * xPos + yPos * yPos)) + 0.5;
			#pragma unroll
			for (unsigned int iz = 0; iz < nz; iz++)
			{
				if (ir < nr){
					exportVol.data[iz + nz * (ix + iy * nXY)] = 
						mip[iz + ir * nz];
				}
				else{
					exportVol.data[iz + nz * (ix + iy * nXY)] = 0;
				}
			}
		}
	}

	// set a quater to 0 again to look inside
	#pragma unroll
	for (unsigned int ix = 0; ix < (nr - 1); ix++)
	{
		#pragma unroll
		for (unsigned int iy = 0; iy < (nr - 1); iy++)
		{
			#pragma unroll
			for (unsigned int iz = 0; iz < nz; iz++)
			{
					exportVol.data[iz + nz * (ix + iy * nXY)] = 0;
			}
		}
	}

	// fill other information into gridded data
	exportVol.dim[0] = nz;
	exportVol.dim[1] = nXY;
	exportVol.dim[2] = nXY;

	exportVol.res[0] = dz;
	exportVol.res[1] = dr;
	exportVol.res[2] = dr;

	exportVol.origin[0] = z0;
	exportVol.origin[1] = -dr * ((float) nr - 1);
	exportVol.origin[2] = -dr * ((float) nr - 1);

	exportVol.nElements = nz * nXY * nXY;

	vtkwriter myWriter;
	string title ("modelMatrix");
	myWriter.set_title(title);

	string type ("STRUCTURED_POINTS");
	myWriter.set_type(type);

	string outputPath(filePath);
	myWriter.set_outputPath(outputPath);

	myWriter.set_structuredPoints(&exportVol);

	myWriter.set_binary();

	myWriter.write();

	delete[] mip;
	delete[] exportVol.data;
	return;
}

// saves the model matrix we just built including all important properties to file
void modelMatrix::saveToFile(const string filePath) const
{
	// properties we want to save
	// 	ntModel, nrModel, nzModel
	//  rRes
	//	zRes
	//	z0
	//	data

	vprintf("Saving data to file... ", 1);
	H5::H5File file(filePath, H5F_ACC_TRUNC);

	const hsize_t col_dims = 1;
	// ntModel
	H5::DataSpace mspaceNTModel(1, &col_dims);
	H5::DataSet ntModelDataset = file.createDataSet(
		"ntModel", H5::PredType::NATIVE_UINT, mspaceNTModel);
	ntModelDataset.write(&nt, H5::PredType::NATIVE_UINT);
	ntModelDataset.close();

	// nrModel
	H5::DataSpace mspaceNRModel(1, &col_dims);
	H5::DataSet nrModelDataset = file.createDataSet(
		"nrModel", H5::PredType::NATIVE_UINT, mspaceNRModel);
	nrModelDataset.write(&nr, H5::PredType::NATIVE_UINT);
	nrModelDataset.close();

	// nzModel
	H5::DataSpace mspaceNZModel(1, &col_dims);
	H5::DataSet nzModelDataset = file.createDataSet(
		"nzModel", H5::PredType::NATIVE_UINT, mspaceNZModel);
	nzModelDataset.write(&nz, H5::PredType::NATIVE_UINT);
	nzModelDataset.close();

	// nZSir
	H5::DataSpace mspaceNZsir(1, &col_dims);
	H5::DataSet nzSirDataset = file.createDataSet(
		"nzSir", H5::PredType::NATIVE_UINT, mspaceNZsir);
	nzSirDataset.write(&nzSir, H5::PredType::NATIVE_UINT);
	nzSirDataset.close();

	// rRes
	H5::DataSpace mspaceRRes(1, &col_dims);
	H5::DataSet rResDataset = file.createDataSet(
		"rRes", H5::PredType::NATIVE_FLOAT, mspaceRRes);
	rResDataset.write(&dr, H5::PredType::NATIVE_FLOAT);
	rResDataset.close();

	// zRes
	H5::DataSpace mspaceZRes(1, &col_dims);
	H5::DataSet zResDataset = file.createDataSet(
		"zRes", H5::PredType::NATIVE_FLOAT, mspaceZRes);
	zResDataset.write(&dz, H5::PredType::NATIVE_FLOAT);
	zResDataset.close();

	// tRes
	H5::DataSpace mspaceTRes(1, &col_dims);
	H5::DataSet tResDataset = file.createDataSet(
		"tRes", H5::PredType::NATIVE_FLOAT, mspaceTRes);
	tResDataset.write(&tRes, H5::PredType::NATIVE_FLOAT);
	tResDataset.close();

	// t0
	H5::DataSpace mspacet0(1, &col_dims);
	H5::DataSet t0Dataset = file.createDataSet(
		"t0", H5::PredType::NATIVE_FLOAT, mspacet0);
	t0Dataset.write(&t0, H5::PredType::NATIVE_FLOAT);
	t0Dataset.close();

	// z0
	H5::DataSpace mspacez0(1, &col_dims);
	H5::DataSet z0Dataset = file.createDataSet(
		"z0", H5::PredType::NATIVE_FLOAT, mspacez0);
	z0Dataset.write(&z0, H5::PredType::NATIVE_FLOAT);
	z0Dataset.close();

	// data (actual model matrix)
	const hsize_t col_data = nElements;
	H5::DataSpace mspaceData(1, &col_data);
	H5::DataSet dataDataset = file.createDataSet(
		"data", H5::PredType::NATIVE_FLOAT, mspaceData);
	dataDataset.write(data, H5::PredType::NATIVE_FLOAT);
	dataDataset.close();

	// also export spatial impulse response [isir, ir, iz]
	const hsize_t col_datasir = nr * nz * nzSir;
	H5::DataSpace mspaceSir(1, &col_datasir);
	H5::DataSet sirDataset = file.createDataSet(
		"sir", H5::PredType::NATIVE_FLOAT, mspaceSir);
	sirDataset.write(sir, H5::PredType::NATIVE_FLOAT);
	sirDataset.close();

	file.close();
	vprintf("done!\n", 0);
	return;
}

// divides each element of the matrix by the norm
void modelMatrix::normalize_matrix()
{
	normalize(data, nElements);
	return;
}

// scales the whole model matrix to a new _maxAbsVal
void modelMatrix::scale(const float _maxAbsVal)
{
	float maxAbsVal = 0;
	for (uint64_t iElement = 0; iElement < nElements; iElement++)
	{
		if (abs(data[iElement]) > maxAbsVal)
			maxAbsVal = abs(data[iElement]);
	}

	for (uint64_t iElement = 0; iElement < nElements; iElement++)
		data[iElement] = data[iElement] / maxAbsVal * _maxAbsVal;
	
	return;
}

void modelMatrix::calcSensField()
{
	if (isSensFieldAlloc)
		delete[] sensField;
	else	
		isSensFieldAlloc = 1;
	
	sensField = new float [nz * nr];

	maxVal = 0;
	uint64_t idx;
	for (uint64_t ir = 0; ir < nr; ir++)
	{
		for (uint64_t iz = 0; iz < nz; iz++)
		{
			sensField[ir + iz * nr] = 0;
			for (uint64_t it = 0; it < nt; it++)
			{
				idx = it + ir * nt + iz * nr * nt;
				if (abs(data[idx]) > sensField[ir + iz * nr])
				{
					sensField[ir + iz * nr] = abs(data[idx]);
				}

			}
			if (sensField[ir + iz * nr] > maxVal)
				maxVal = sensField[ir + iz * nr];
		}
	}

	return;
}

float* modelMatrix::get_pmodSlice(const uint64_t iz)
{
	uint64_t izCorr = iz;
	if (iz >= nz)
		izCorr = iz - 1;

	if (lastSlice != izCorr)
	{
		cout << "generating new model slive "<< endl;
		// erset maxValSlice and minValSlice to 0
		maxValSlice = 0.0;
		minValSlice = 0.0;

		for (uint64_t ir = 0; ir < nr; ir++)
		{
			for (uint64_t it = 0; it < nt; it++)
			{
				modSlice[ir + it * nr] = data[it + nt * (ir + izCorr * nr)];
				if (modSlice[ir + it * nr] < minValSlice)
					minValSlice = modSlice[ir + it * nr];

				if (modSlice[ir + it * nr] > maxValSlice)
					maxValSlice = modSlice[ir + it * nr];

			}
		}
		lastSlice = izCorr;
	}
	return modSlice;
}

float* modelMatrix::get_pmodSlice(const float zLevel)
{
	float zCorr = zLevel;
	if (zCorr < get_zMin())
		zCorr = get_zMin();

	if (zCorr > get_zMax())
		zCorr = get_zMax();

	uint64_t zIdx = (zCorr - z0) / dz + 0.5;
	return get_pmodSlice(zIdx);
}

float *modelMatrix::get_modelVec(const uint64_t ir, const uint64_t iz)
{
	
	uint64_t offsetData = nt * (ir + iz * nr);
	if ((modVecR != ir) || (modVecZ != iz))
	{
		cout << "generating new model vector plot" << endl;
		modVecR = ir; // update model vector position to new value
		modVecZ = iz; // update model vector position to new value
		// update max and min value
		modVecAbsMax = 0;
		for (uint64_t it = 0; it < nt; it++)
		{
			if (abs(data[offsetData + it]) > modVecAbsMax)
				modVecAbsMax = abs(data[offsetData + it]);
		}
	}

	return (data + offsetData);
}

float *modelMatrix::get_modelVec(const float rPos, const float zPos)
{
	// convert radial position into index
	uint64_t rIdx = (rPos / dr) + 0.5;
	if (rIdx >= nr)
		rIdx = nr - 1;
			
	// covnert axial position into index
	uint64_t zIdx = ((zPos - z0) / dz) + 0.5;
	if (zIdx >= nz)
		zIdx = nz - 1;

	return get_modelVec(rIdx, zIdx);
}