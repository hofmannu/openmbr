#include "model.h"

model::~model()
{
	if (isModelDataAlloc)
		delete[] modelData;


}

// here go all out operator overloaders
model& model::operator = (model& modelIn)
{
	this->set_dim(modelIn.get_nt(), modelIn.get_nx(), modelIn.get_ny(), modelIn.get_nz());
	this->alloc_model();
	memcpy(this->modelData, modelIn.get_pdata(), get_nElements() * sizeof(float));

	// now we push origin and resolution over
	for (uint8_t iDim = 0; iDim < 4; iDim++)
	{
		this->set_res(iDim, modelIn.get_res(iDim));
		this->set_origin(iDim, modelIn.get_origin(iDim));
	}
	return *this;
}

// indexing of arrays: [t, x, y, z]
void model::crop(const uint64_t* startIdx, const uint64_t* stopIdx)
{
	// first check if everything we have looks valid in all 4 dimensions
	uint64_t newDim[4];
	#pragma unroll
	for (uint8_t iDim = 0; iDim < 4; iDim++)
	{
		if (startIdx[iDim] > stopIdx[iDim])
		{
			printf("Stop index needs to be larger than start index\n");
			throw "InvalidValue";
		}

		if (stopIdx[iDim] >= dim[iDim])
		{
			printf("Hitting upper boundary of model matrix along dimension %d (%lu of %lu)\n",
				iDim, stopIdx[iDim], dim[iDim] - 1);
			throw "InvalidValue";
		}
		newDim[iDim] = stopIdx[iDim] - startIdx[iDim] + 1;
	}

	const uint64_t newElements = newDim[0] * newDim[1] * newDim[2] * newDim[3];
	float *newData = new float[newElements];

	for (uint64_t iz = 0; iz < newDim[3]; iz++) 
	{
		const uint64_t zOld = iz + startIdx[3];
		for (uint64_t it = 0; it < newDim[0]; it++) 
		{
			const uint64_t tOld = it + startIdx[0];
			for (uint64_t iy = 0; iy < newDim[2]; iy++)
			{
				const uint64_t yOld = iy + startIdx[2];
				const uint64_t newIdx = newDim[1] * (iy + newDim[2] * (it + newDim[0] * iz));
				const uint64_t oldIdx = dim[1] * (yOld + dim[2] * (tOld + dim[0] * zOld));
				memcpy(&newData[newIdx], &modelData[oldIdx], newDim[1] * sizeof(float));
			}	
		}
	}
	
	#pragma unroll
	for (uint8_t iDim = 0; iDim < 4; iDim++)
	{
		dim[iDim] = stopIdx[iDim] - startIdx[iDim] + 1;
		origin[iDim] = origin[iDim] + res[iDim] * ((float) startIdx[iDim]);
	}

	delete[] modelData;
	modelData = newData;

	return;
}


// here go all set and get functions
void model::set_dim(const uint64_t nt, const uint64_t nx, const uint64_t ny, const uint64_t nz)
{
	dim[0] = nt;
	dim[1] = nx;
	dim[2] = ny;
	dim[3] = nz;
	alloc_model();
	return;
}

void model::set_origin(const uint8_t iDim, const float originIn)
{
	origin[iDim] = originIn;
	return;
}

void model::set_res(const uint8_t iDim, const float resIn)
{
	res[iDim] = resIn;
	return;
}

uint64_t model::get_tIdx(const float tPos) const
{	
	if (tPos < get_t0())
	{
		printf("The position you try to crop is outside of the model range");
		throw "InvalidValue";
	}

	const float deltaT = tPos - get_t0();
	const float tIdx = deltaT / get_dt();
	return ((uint64_t) (tIdx + 0.5f));
}


// here come all other fancy functions

void model::calc_sensField()
{
	sensField.set_dim(get_nx(), get_ny(), get_nz());
	sensField.alloc_memory();
	for (uint64_t iz = 0; iz < get_nz(); iz++)
	{
		for (uint64_t iy = 0; iy < get_ny(); iy++)
		{
			for (uint64_t ix = 0; ix < get_nx(); ix++)
			{
				float tempAbsMax = 0;
				for (uint64_t it = 0; it < get_nt(); it++)
				{
					const float currVal = fabsf(get_val(ix, iy, it, iz));
					if (currVal > tempAbsMax)
					{
						tempAbsMax = currVal;
					}
				}
				sensField.set_value(ix, iy, iz, tempAbsMax);
			}
		}
	}
	sensField.calcMinMax();
	return;
}



void model::alloc_model()
{
	if (isModelDataAlloc)
	{
		delete[] modelData;
	}	
	modelData = new float[get_nElements()];
	isModelDataAlloc = 1;
	return;
}

// load model matrix from an h5 file, careful with data order here!!!
void model::load_from_file(const string filePath)
{
	H5::H5File file(filePath, H5F_ACC_RDONLY);

	H5::DataSet resDataset = file.openDataSet("dr"); // init dataset for res 
	const hsize_t col_dims = 4;
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
	
	alloc_model();

	H5::DataSet modelDataset = file.openDataSet("modelMatrix");
	const hsize_t col_model = get_nElements(); 
	H5::DataSpace mspaceModel (1, &col_model);
	filespace = modelDataset.getSpace();
	modelDataset.read(modelData, H5::PredType::NATIVE_FLOAT, mspaceModel, filespace); 	
	
	isLoaded = 1;

	return;
}

inline float model::get_val(
	const uint64_t idxX, const uint64_t idxY, 
	const uint64_t idxT, const uint64_t idxZ) const
{
	return modelData[idxX + get_nx() * (idxY + get_ny() * (idxT + get_nt() * idxZ))];
}

// a small helper function to print some basic information about the model
void model::print_information()
{
	printf("Model information\n");
	// printf("File path: %s\n", );
	printf("Resolution [t, z, y, x]: %.1f ns, %.1f microm, %.1f microm, %.1f microm\n",
		res[0]  * 1e9, res[1] * 1e6, res[2] * 1e6, res[3] * 1e6);
	printf("Origin [t, z, y, x]: %.1f ns, %.1f mm, %.1f mm, %.1f mm\n",
		origin[0]  * 1e9, origin[1] * 1e3, origin[2] * 1e3, origin[3] * 1e3);
	printf("Dimensions [t, z, y, x]: %lu, %lu, %lu, %lu\n",
		dim[0], dim[1], dim[2], dim[3]);

	return;
}