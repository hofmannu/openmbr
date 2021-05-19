#include "mbrecon.cuh"

// class constructor
mbrecon::mbrecon()
{
	className = "mbrecon";
}

// class destructor
mbrecon::~mbrecon()
{

}

bool mbrecon::get_isRecon() const {return isRecon;}
volume* mbrecon::get_ppreprocVol() {return &preprocVol;}
volume* mbrecon::get_preconVol() {return &reconVol;}
transducerProperties* mbrecon::get_ptransProp() {return &transProp;};
transModel* mbrecon::get_ptransModel() {return &transMod;};
reconSettings* mbrecon::get_psett() {return &sett;};

// let transducer properties be overwritten from outside
void mbrecon::set_transProp(const transducerProperties _transProp)
{
	transProp = _transProp;
	return;
}

// define reconstruction settings
void mbrecon::set_sett(const reconSettings _sett)
{
	sett = _sett;
	return;
}

// load input data from hard drive
void mbrecon::load_data()
{
	load_data(inputPath);
	return;
}

// loads data from a file defined by _inputPath
void mbrecon::load_data(const string _inputPath)
{
	// check if we have a valid file
	ifstream f(_inputPath.c_str());
	if (f.good())
		inputPath = _inputPath;
	else{
		vprintf("Path is not pointing to a file.\n", 1);
		throw "Path is pointing to an invalid file";
	}

	// load volumetric dataPermutedaset from file
	preprocVol.readFromFile(_inputPath);
	preprocVol.printInformation();
	isRawDataLoaded = 1;

	// load transducer name from file
	// TODO needs to be adapted to dynamic loading
	// transProp.setName("russian_01");
	// transProp.readFromFile();
	// trans.printProperties();
	
}

// calculate number of elements required in radial direction
void mbrecon::calc_nr()
{
	float* zRecon = sett.get_zRecon();
	float upperDist = abs(zRecon[0] - transProp.getFocalDistance());
	float lowerDist = abs(zRecon[1] - transProp.getFocalDistance());
	float maxDist = upperDist > lowerDist ? upperDist : lowerDist;

	rMax = maxDist * tan(transProp.getTheta()); // maximum radial distance
	nr = rMax / sett.get_rRes() + 1;

	if (flagDebug)
	{
		printf("[debug] nr = %d\n", nr);
		printf("[debug] rMax [mm] = %.3f\n", rMax * 1e3);
	}

	return;
}

// build GPU model 
void mbrecon::build_model()
{
	bool flagDebug = 1;

	// define properties of simulation field
	float dt = preprocVol.getRes0() * ((float) sett.get_downsampling(0));
	field.set_dt(dt); // define temporal spacing and sampling freq
	field.set_dr(sett.get_rRes()); // get radial spacing of model
	field.set_dz(sett.get_dr0()); // get axial spacing of model
	field.set_SOS(sett.get_SOS()); // get speed of sound
	field.set_nr(nr); // get number of elements in rad direction
	field.set_nz(sett.get_dim0()); // get number of elemets in axial direction
	field.set_z0(sett.get_zLower() - transProp.getFocalDistance()); // get starting point in axial direction
	
	// print additional information if we are currently in debugging mode
	if (flagDebug)
	{
		field.print_properties();
		sett.print_properties();
		// reconVol.printInformation();
	}
	
	transMod.setFieldProperties(field); // push field properties over to class
	transMod.setTransProp(transProp); // push transducer properties over to class
	//transMod.loadImpulseResponse(); // load impulse response from file
	transMod.buildGPUModel(); // start actual model building process
	model = transMod.getModel(); // return model from model builder
	
	// convolve model with optical model
	gaussian_beam* gauss = transMod.get_pgauss();
	gauss->define_r(sett.get_rRes(), (uint64_t) nr);
	gauss->define_z(model->get_dz(), model->get_z0(), model->get_nz());
	gauss->convolve_model(model->get_data(), model->get_nt());
	
	model->scale(1.0); // scale whole model to range [-1 ... 1]
	model->permute(); // generate permuted version of model matrix
	model->getStartStopIdx(); // generate start and stop indices
	
	model->exportSensFieldVtk(); // save sensitivity field as vtk to file
	return;
}

// extract cropped volume from signal matrix
void mbrecon::getCroppedVol()
{
	// check if rMax was already defined
	if (rMax == 0)
	{
		printf("Invalid radial range defined for reconstruction, run calc_nr() first\n");
		throw "invalidVal";
	}

	uint64_t startIdx[3], stopIdx[3];
	
	// determine temporal start and stop index in signal matrix
	startIdx[0] = preprocVol.getIdx0(transMod.getStarttime());
	stopIdx[0] = startIdx[0] + transMod.getNtModel() - 1;

	// define spatial start and stop index in signal matrix
	startIdx[1] = preprocVol.getIdx1(sett.get_xLower() - rMax); 
	stopIdx[1] = preprocVol.getIdx1(sett.get_xUpper() + rMax);
	
	// define spatial start and stop index in signal matrix
	startIdx[2] = preprocVol.getIdx2(sett.get_yLower() - rMax);
	stopIdx[2] = preprocVol.getIdx2(sett.get_yUpper() + rMax);
	
	croppedVol.nElements = 1;
	for (uint8_t iDim = 0; iDim < 3; iDim++)
	{
		croppedVol.dim[iDim] = stopIdx[iDim] - startIdx[iDim] + 1;
		croppedVol.nElements = croppedVol.nElements * croppedVol.dim[iDim];
	}
	// allocate memory for new array
	croppedVol.data = new float [croppedVol.nElements];
	preprocVol.getCroppedVolume(croppedVol.data, &startIdx[0], &stopIdx[0]);
	
	// push dimensionality over to gridded data
	for (unsigned char iDim = 0; iDim < 3; iDim++)
	{
		croppedVol.origin[iDim] = preprocVol.getPos(startIdx[iDim], iDim);
		croppedVol.res[iDim] = preprocVol.getRes(iDim) * ((float) sett.get_downsampling(iDim));
		croppedVol.ires[iDim] = 1 / croppedVol.res[iDim];
	}

	if (flagDebug)
	{
		printf("[debug] *** Cropped volume information ***\n");
		printf("[debug] - startIdx: %d, %d, %d\n", startIdx[0], startIdx[1], startIdx[2]);
		printf("[debug] - stopIdx: %d, %d, %d\n", stopIdx[0], stopIdx[1], stopIdx[2]);
		printf("[debug] - resolution [ns/microm]: %f, %f, %f\n", 
			croppedVol.res[0] * 1e9, croppedVol.res[1] * 1e6, croppedVol.res[2] * 1e6);
		printf("[debug] - dimensions: %d, %d, %d\n", 
			croppedVol.dim[0], croppedVol.dim[1], croppedVol.dim[2]);
		printf("[debug] - origin [micros, millim]: %f, %f, %f\n\n",
			croppedVol.origin[0] * 1e6, croppedVol.origin[1] * 1e3, croppedVol.origin[2] * 1e3);
		
		// save croppedVol for paraview visualization
		// keep in mind that axis will be flipped:
		vtkwriter outputter;
		const string title ("croppedVol"); // generate title
		outputter.set_title(title); // define title of outut volume
		const string type ("STRUCTURED_POINTS"); 
		outputter.set_type(type); // define output type
		const string outputPath ("/home/hofmannu/croppedVol.vtk");
		outputter.set_outputPath(outputPath);

		// convert time domain into z position
		croppedVol.res[0] = croppedVol.res[0] * 1490.0; 
		croppedVol.origin[0] = croppedVol.origin[0] * 1490.0;

		outputter.set_structuredPoints(&croppedVol);
		outputter.set_binary();
		outputter.write();
		
		croppedVol.res[0] = croppedVol.res[0] / 1490.0;
		croppedVol.origin[0] = croppedVol.origin[0] / 1490.0;
	}
	
	isDataCropped = 1;
	return;
}

// prepares the estimated absorber matrix
void mbrecon::prepareOutputVolume()
{
	// define size of output volume
	estimAbs.nElements = 1;
	for (unsigned char iDim = 0; iDim < 3; iDim++){
		estimAbs.res[iDim] = sett.get_dr(iDim);
		estimAbs.ires[iDim] = 1 / sett.get_dr(iDim);
		estimAbs.dim[iDim] = sett.get_dim(iDim);
		estimAbs.nElements *= estimAbs.dim[iDim];
	}

	estimAbs.origin[0] = sett.get_zLower();
	estimAbs.origin[1] = sett.get_xLower();
	estimAbs.origin[2] = sett.get_yLower();

	// push information over to volumetric dataset
	for (uint8_t iDim = 0; iDim < 3; iDim++)
	{
		reconVol.setRes(iDim, estimAbs.res[iDim]);
		reconVol.setOrigin(iDim, estimAbs.origin[iDim]);
		reconVol.setDim(iDim, estimAbs.dim[iDim]);
	}
	reconVol.allocMemory();

	// allocate memory and fill with initial 0s
	estimAbs.data = reconVol.get_pdata();
	for (unsigned int iElement = 0; iElement < estimAbs.nElements; iElement++)
		estimAbs.data[iElement] = 0;

	if (flagDebug)
	{
		reconVol.printInformation();
	}
	return;	
}

// export volume for paraview
void mbrecon::export_vtk()
{
	// create copy of datset before messing with polarity
	float* dataCopy = new float [estimAbs.nElements];
	assign(dataCopy, estimAbs.data, estimAbs.nElements);
	const string polarityHandling (sett.get_polarityHandling());
	handlePolarity(estimAbs.data, estimAbs.nElements, polarityHandling);

	vtkwriter outputter; // prepare output pipeline
	const string title ("reconVol"); // generate title
	outputter.set_title(title); // define title of outut volume
	const string type ("STRUCTURED_POINTS"); 
	outputter.set_type(type); // define output type
	const string outputPath ("/home/hofmannu/testMbrecon.vtk");
	outputter.set_outputPath(outputPath);

	// float * dataOld = estimAbs.data;
	outputter.set_structuredPoints(&estimAbs);
	outputter.set_binary();
	outputter.write();
	
	// get back untreated bipolar dataset 
	assign(estimAbs.data, dataCopy, estimAbs.nElements);
	delete[] dataCopy;
	return;
}

void mbrecon::recon()
{
	uint8_t emptyIterator = 0;
	recon(&emptyIterator);
	return;
}

void mbrecon::simulate()
{

	// prepare estimAbs from reconVol
	estimAbs.nElements = 1;
	for (uint8_t iDim = 0; iDim < 3; iDim++)
	{
		estimAbs.res[iDim] = reconVol.getRes(iDim);
		estimAbs.ires[iDim] = 1 / reconVol.getRes(iDim);
		estimAbs.dim[iDim] = reconVol.getDim(iDim);
		estimAbs.origin[iDim] = reconVol.getOrigin(iDim); 
		estimAbs.nElements *= estimAbs.dim[iDim];
	}

	// prepare signal matrix depending on reconVol
	croppedVol.nElements = 1;
	for (uint8_t iDim = 1; iDim < 3; iDim++)
	{
		croppedVol.res[iDim] = reconVol.getRes(iDim);
		croppedVol.ires[iDim] = 1 / reconVol.getRes(iDim);
		croppedVol.dim[iDim] = reconVol.getDim(iDim);
		croppedVol.origin[iDim] = reconVol.getOrigin(iDim);
		croppedVol.nElements *= croppedVol.dim[iDim];
		preprocVol.setDim(iDim, croppedVol.dim[iDim]);
		preprocVol.setRes(iDim, croppedVol.res[iDim]);
		preprocVol.setOrigin(iDim, croppedVol.origin[iDim]);
	}

	croppedVol.res[0] = preprocVol.getRes(1);
	croppedVol.ires[0] = 1 / preprocVol.getRes(1);
	preprocVol.setRes(0, 1.0 / 250e6);
	
	calc_nr();
	build_model(); // build transducer model
	
	croppedVol.dim[0] = model->get_nt(); // push over number of time points
	croppedVol.nElements *= croppedVol.dim[0];
	
	preprocVol.setOrigin(0, model->get_t0());
	preprocVol.setDim(0, croppedVol.dim[0]);

	preprocVol.allocMemory();

	// simulate forward model
	forward(preprocVol.get_pdata(), reconVol.get_pdata());

	preprocVol.calcMinMax();
	// add noise to array if wanted
	return;
}

// actual iterative reconstruction procedure
void mbrecon::recon(uint8_t* iterator)
{
	unsigned int* ds = sett.get_downsampling(); // downsampling vector for readability
	calc_nr(); // calculate max radial distance to take into account
	build_model(); // build transducer model
	getCroppedVol(); // crop input volume to requried range
	prepareOutputVolume(); // prepares the output volume

	float* r = new float [estimAbs.nElements]; // r has dimension [nz, nx, ny]
	float* v = new float [estimAbs.nElements]; // r has dimension [nz, nx, ny]
	float* p = new float [croppedVol.nElements]; // p has dimension [nt, nx, ny]
	float* w = new float [estimAbs.nElements];
	float* wWeighted = new float [estimAbs.nElements];

	// normalize signal matrix
	const float beta = getNorm(croppedVol.data, croppedVol.nElements);
	
	if (flagDebug)
		printf("beta = %f\n", beta);
	
	divide(croppedVol.data, beta, croppedVol.nElements);

	unsigned int rIdx = rand() % croppedVol.nElements;
	
	if (flagDebug)
		printf("Random element in matrix v / w: %f\n", croppedVol.data[rIdx]);
	
	// r = A' * croppedVol.data
	transpose(r, croppedVol.data);
	// inserted for debugging
	// assign(estimAbs.data, r, estimAbs.nElements);

	rIdx = rand() % estimAbs.nElements;
	if (flagDebug)
		printf("Random element in matrix r: %f\n", r[rIdx]);
	
	float alpha = getNorm(r, estimAbs.nElements);
	if (flagDebug)
		printf("alpha = %f\n", alpha);
	
	divide(v, r, alpha, estimAbs.nElements); // v = r / alpa

	float phi_bar = beta;
	float rho_bar = alpha;

	assign(w, v, estimAbs.nElements); // w = v

	rIdx = rand() % estimAbs.nElements; // get random element of estimAbs
	if (flagDebug)
		printf("Random element in matrix v / w: %f\n", v[rIdx]);

	float norm_p;

	// sett.get_nIter()
	for (uint8_t iIter = 0; iIter < sett.get_nIter(); iIter++){
		
		*iterator = iIter;
		printf("Running iteration %d of %d...\n", iIter + 1, sett.get_nIter());
		
		// bidiagonalization 
		// p = A * v - alpha * croppedVol.data
		forward(p, v); // A * v
		rIdx = rand() % croppedVol.nElements; 
		printf("Random element in matrix p: %f\n", p[rIdx]);
		multiply(croppedVol.data, alpha, croppedVol.nElements); // alpha * croppedVol.data
		substract(p, croppedVol.data, croppedVol.nElements); // p - alpha * croppedVol.data

		// determine norm_p depending on reg method
		if (sett.get_flagRegularization())
		{
			float norm_p_ori = getNorm(p, croppedVol.nElements);
			norm_p_ori = pow(norm_p_ori, 2);
			multiply(wWeighted, v, sett.get_lambdaReg(), estimAbs.nElements); // wWeighted = v * lambda
			float norm_p_reg = getNorm(wWeighted, estimAbs.nElements); // norm(wWeighted)
			norm_p_reg = pow(norm_p_reg, 2); 
			norm_p = pow(norm_p_ori + norm_p_reg, 0.5);
		}
		else
		{
			norm_p = getNorm(p, croppedVol.nElements);
		}

		// croppedVol.data = p / beta
		divide(croppedVol.data, p, norm_p, croppedVol.nElements);

		// r = A' * croppedVol.data - beta * v;
		transpose(r, croppedVol.data); // r = A' * croppedVol.data

		rIdx = rand() % estimAbs.nElements;
		printf("Random element in matrix r: %f\n", r[rIdx]);

		if (sett.get_flagRegularization())
		{
			multiply(wWeighted, v, pow(sett.get_lambdaReg(), 2) / norm_p, estimAbs.nElements);
			// wWeighted = v * lambda^2 / norm_p
			add(r, wWeighted, estimAbs.nElements); // r = r + wWeighted 
		}

		// r = ATu - norm_p .* v
		multiply(v, beta, estimAbs.nElements); // v = v * beta
		substract(r, v, estimAbs.nElements); // r = r - v

		alpha = getNorm(r, estimAbs.nElements); // alpha = norm(r);
		// v = r / alpha;
		divide(v, r, alpha, estimAbs.nElements);
		
		// % ------------- orthogonal transformation ---------------
		// rrho = norm([rho_bar, beta]);
		float rrho = rho_bar * rho_bar + beta * beta;
		rrho = pow(rrho, 0.5);
		float c1 = rho_bar / rrho;
		float s1 = norm_p / rrho;
		float theta = s1 * alpha;
		rho_bar = -c1 * alpha;
		float phi = c1 * phi_bar;
		phi_bar = s1 * phi_bar;

		// update solution and search direction
		// x = x + (phi / rrho) * w
		multiply(wWeighted, w, phi / rrho, estimAbs.nElements); // wWeighted = w * phi / rrho
		add(estimAbs.data, wWeighted, estimAbs.nElements); // estimAbs.data = estim.data + wWeighted

		// w = v - (theta / rrho) * w;
		multiply(wWeighted, w, theta / rrho, estimAbs.nElements);
		substract(w, v, wWeighted, estimAbs.nElements);
		
		// export vtk after each iteration for preview purpose
		export_vtk();
	}
  isRecon = 1; // indicate that reconstruction is done

  	// free all temorarily used variables
	delete[] croppedVol.data;
	delete[] wWeighted;
	delete[] w;
	delete[] r;
	delete[] p;
	delete[] v;

	return;
}


// subkernel running for one reconstruction column along z
__global__ void kernelTVec(
		float* calcSig, // caluclated signal matrix (output volume) it, ix, iy
		const constArgsInFwd* inArgs, // constant global input arguments
		const constArgKFwd inArgsT // subkernel specific input arguments
	)
{
	const uint64_t tIm = threadIdx.x + blockIdx.x * blockDim.x; // it output matrix
	if (tIm < inArgs->outputGrid->dim[0]){
		// index of currently time signal in calcSig
		const uint64_t outputIdx = tIm + inArgs->outputGrid->dim[0] * 
			(inArgsT.xIm + inArgsT.yIm * inArgs->outputGrid->dim[1]);		
		const uint64_t mmIdx = inArgs->inputGrid->dim[0] * (inArgsT.ir + tIm * inArgs->nr);
		const uint64_t idxMStd = 2 * (inArgsT.ir + tIm * inArgs->nr);
		const uint64_t startIdx =  inArgs->modelStIdx[idxMStd];
		const uint64_t stopIdx = inArgs->modelStIdx[idxMStd + 1];
		
		// use temporary variable here to increase speed and decrease memory access
		float temp = 0;

		#pragma unroll
		for (uint64_t iz = startIdx; iz <= stopIdx; iz++){
		//for (unsigned int iz = 0; iz < inArgs->inputGrid->dim[0]; iz++){
			// const unsigned int iFluenceIdx = iSurface + inArgs->inputGrid->dim[0] * (
			// 	inArgsT.ir + iz * inArgs->nr);
			// 
			// calcSig = calcSig + abs * model * fluence
			// calcSig[outputIdx] = __fadd_rn(calcSig[outputIdx], 
			// 	__fmul_rn(inArgs->estimatedAbs[iVolIdx + iz], 
			// 	__fmul_rn(inArgs->modelMat[mmIdx + iz], 
			// 		inArgs->fluenceModel[iFluenceIdx])));

			// calcSig = calcSig + abs * model
			temp = __fadd_rn(temp, 
				__fmul_rn(inArgs->estimatedAbs[inArgsT.voloff + iz],	inArgs->modelMat[mmIdx + iz]));
		}
		
		// push temporary variable to output volume
		calcSig[outputIdx] = temp;
	}
	return;
}


__global__ void modelBasedFwd(
	float* calcSig, // [itOut + ixOut * ntOut + iyOut * ntOut * nxOut]
	const constArgsInFwd* inArgs
)
{
	constArgKFwd inArgsT; // struct to hold subthread information
	inArgsT.xIm = blockIdx.x * blockDim.x + threadIdx.x; // x index of recon voxel
	inArgsT.yIm = blockIdx.y * blockDim.y + threadIdx.y; // y index of recon voxel
	if ( // check if image voxel actually exists in volume (kernel overlap)
		(inArgsT.xIm < inArgs->outputGrid->dim[1]) && (inArgsT.yIm < inArgs->outputGrid->dim[2])){

		int outputIdx = inArgs->outputGrid->dim[0] * (inArgsT.xIm + inArgsT.yIm * inArgs->outputGrid->dim[1]);
		for(unsigned int tIm = 0; tIm < inArgs->outputGrid->dim[0]; tIm++) 
			calcSig[outputIdx + tIm] = 0; // set value to 0 before summing up
	
		// get reconstructed point x and y position	
		const float xPos = inArgs->outputGrid->origin[1] + 
			__fmul_rn(((float) inArgsT.xIm), inArgs->outputGrid->res[1]);
		const float yPos = inArgs->outputGrid->origin[2] + 
			__fmul_rn(((float) inArgsT.yIm), inArgs->outputGrid->res[2]);

		// calc recon range in true sized input matrix & convert to index range
		const float rRange = __fmul_rn(inArgs->dr, ((float) inArgs->nr));
		const float xmin = xPos - rRange;
		int xminI = __float2int_ru((xmin - inArgs->inputGrid->origin[1]) * inArgs->inputGrid->ires[1]); // round up
		if (xminI < 0)
			xminI = 0;
		
		const float xmax = __fadd_rn(xPos, rRange);
		int xmaxI = __float2int_rd((xmax - inArgs->inputGrid->origin[1]) * inArgs->inputGrid->ires[1]); // round down
		if (xmaxI >= inArgs->inputGrid->dim[1])
			xmaxI = inArgs->inputGrid->dim[1] - 1;

		const float ymin = yPos - rRange;
		int yminI = __float2int_ru((ymin - inArgs->inputGrid->origin[2]) * inArgs->inputGrid->ires[2]); // round up
		if (yminI < 0)
			yminI = 0;

		const float ymax = __fadd_rn(yPos, rRange);
		int ymaxI = __float2int_rd((ymax - inArgs->inputGrid->origin[2]) * inArgs->inputGrid->ires[2]); // round down
		if (ymaxI >= inArgs->inputGrid->dim[2])
			ymaxI = inArgs->inputGrid->dim[2] - 1;
		
		float rDist, deltaX, deltaY;
		
		// starting y index offset of estimatedAbs at beginning of for loop
		const dim3 blockSize(32);
		const dim3 gridSize(
			(inArgs->outputGrid->dim[0] + blockSize.x - 1) / blockSize.x);

		inArgsT.yoff = yminI * inArgs->inputGrid->dim[0] * inArgs->inputGrid->dim[1]; 
		#pragma unroll 4
		for (int yidx = yminI; yidx < ymaxI; yidx++){
			deltaY = yPos - __fadd_rn(inArgs->inputGrid->origin[2], 
					__fmul_rn((float) yidx, inArgs->inputGrid->res[2]));
			inArgsT.yAScan = yidx;
			inArgsT.xoff = xminI * inArgs->inputGrid->dim[0];
			#pragma unroll 4
			for (int xidx = xminI; xidx < xmaxI; xidx++){
				deltaX = xPos - __fadd_rn(inArgs->inputGrid->origin[1], 
						__fmul_rn((float) xidx, inArgs->inputGrid->res[1]));
				rDist = __fsqrt_rn(__fadd_rn(
							__fmul_rn(deltaX, deltaX),
						 	__fmul_rn(deltaY, deltaY)));
				inArgsT.ir = __float2int_rn(__fdividef(rDist, inArgs->dr));
				inArgsT.xAScan = xidx;	
				inArgsT.voloff = inArgsT.xoff + inArgsT.yoff;
				// crop rectangularly defined area to circle
				if (inArgsT.ir < inArgs->nr){
					kernelTVec<<<gridSize, blockSize>>>(calcSig, inArgs, inArgsT);
					cudaDeviceSynchronize();	
				}
				inArgsT.xoff += inArgs->inputGrid->dim[0]; 
				// move x off to next value
			}
			inArgsT.yoff += inArgs->inputGrid->dim[0] * inArgs->inputGrid->dim[1]; 
			// move y off to next value
		}
	}
	return;
}

// forward model calculation converting an estimated absorber map into a signal
// input: 
//		- signalMatrix (pointer to output data)
//		- absorberMatrix (pointer to estimated absorber map)

void mbrecon::forward(float* signalMatrix, const float* absorberMatrix)
{

	// this motherfucker requires implementation!!!
	cudaError_t err;

	// allocate memory on GPU and fill struct with pointers
	constArgsInFwd* inArg_dev;
	constArgsInFwd inArg;
	float* calcSig_dev;

	err = cudaMalloc( (void**) &inArg.estimatedAbs, estimAbs.nElements * 
		sizeof(float));
	checkCudaReturn(err, "Could not allocate absorber map on GPU");
	err = cudaMalloc( (void**) &calcSig_dev, croppedVol.nElements * sizeof(float));
	checkCudaReturn(err, "Could not allocate memory for signal matrix on GPU");
	err = cudaMalloc( (void**) &inArg.modelMat,	model->getNElements() * sizeof(float));
	checkCudaReturn(err, "Could not allocate memory for model matrix on GPU");
	err = cudaMalloc( (void**) &inArg.inputGrid, sizeof(griddedData));
	checkCudaReturn(err, "Could not allocate memory for input grid on GPU");
	err = cudaMalloc( (void**) &inArg.outputGrid, sizeof(griddedData));
	checkCudaReturn(err, "Could not allocate memory for output grid on GPU");
	err = cudaMalloc( (void**) &inArg.modelStIdx,	2 * nr * croppedVol.dim[0] * 
		sizeof(uint64_t));
	checkCudaReturn(err, "Could not allocate memory for start and stop index on card");
	inArg.dr = sett.get_rRes();
	inArg.nr = nr;

	err = cudaMalloc( (void**) &inArg_dev, sizeof(constArgsInFwd));
	checkCudaReturn(err, "Could not allocate memory for input struct on card");

	err = cudaMemcpy(inArg.estimatedAbs, absorberMatrix, estimAbs.nElements *
		sizeof(float), cudaMemcpyHostToDevice);
	checkCudaReturn(err, "Could not copy absorber map to device");
	err = cudaMemcpy(inArg.modelMat, model->dataPermuted, model->getNElements() * sizeof(float), 
		cudaMemcpyHostToDevice);
	checkCudaReturn(err, "Could not copy model matrix to device");
	err = cudaMemcpy(inArg.inputGrid, &estimAbs, sizeof(griddedData), 
		cudaMemcpyHostToDevice);
	checkCudaReturn(err, "Could not copy input grid to device");
	err = cudaMemcpy(inArg.modelStIdx, model->startStopIdxPermuted, 
		2 * nr * croppedVol.dim[0] * sizeof(uint64_t), 
		cudaMemcpyHostToDevice);
	checkCudaReturn(err, "Could not copy start and stop indices to device");
	err =  cudaMemcpy(inArg.outputGrid, &croppedVol,	sizeof(griddedData), 
		cudaMemcpyHostToDevice);
	checkCudaReturn(err, "Could not copy output grid to device");
	err = cudaMemcpy(inArg_dev, &inArg, sizeof(constArgsInFwd), cudaMemcpyHostToDevice);
	checkCudaReturn(err, "Could not copy input struct to device");

	if (flagDebug)
	{
		printf("[debug] Dimensions of absorber volume (nz x nx x ny): %d x %d x %d\n", 
			estimAbs.dim[0], estimAbs.dim[1], estimAbs.dim[2]);
		printf("[debug] Dimensions of signal matrix (nt x nx x ny): %d x %d x %d\n", 
			croppedVol.dim[0], croppedVol.dim[1], croppedVol.dim[2]);
		printf("[debug] Dimensions of model matrix (nz x nr x nt): %d x %d x %d\n", 
			model->get_nz(), model->get_nr(), model->get_nt());
	}
	// defining thread and block size as number of a scan positions
	// this block size is highly critical!!!
	dim3 blockSize(8, 1);
	// 1, 1 --> 
	// 2, 2 --> 1.39 s
	// 3, 3 --> 
	// 4, 1 --> 1.37 s
	// 4, 4 --> 3.57 s
	// 8, 1 --> 1.07 s
	// 8, 2 --> 3.72 s
	// 8, 8 -->
	// 16, 1 --> 3.56 s
	// 32, 1 --> 6.50 s
	// 32, 4 --> illegal memory access

	dim3 gridSize(
		(croppedVol.dim[1] + blockSize.x - 1) / blockSize.x,
		(croppedVol.dim[2] + blockSize.y - 1) / blockSize.y);
	// printf("%d, %d\n", gridSize.x, gridSize.y);
	
	vprintf("Starting forward kernel...\n", 1); fflush(stdout);
	clock_t tStart = clock();
	modelBasedFwd<<< gridSize , blockSize >>>(
		calcSig_dev, // calculated signal (output)
		inArg_dev // constant input arguments
	);
	cudaDeviceSynchronize();
	vprintf("Finished kernel execution...\n", 1);
	printf("Time taken: %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);
	err = cudaGetLastError();
	checkCudaReturn(err, "Error during kernel execution");
		
	vprintf("Copy data back from GPU...\n"); // copy matrix back from GPU
	err = cudaMemcpy(signalMatrix, calcSig_dev, croppedVol.nElements * sizeof(float), 
		cudaMemcpyDeviceToHost);
	checkCudaReturn(err, "Could not copy data back from GPU");
	
	// free device and host memory
	cudaFree(inArg.modelStIdx);
	cudaFree(calcSig_dev);
	cudaFree(inArg.inputGrid);
	cudaFree(inArg.modelMat);
	cudaFree(inArg.outputGrid);
	cudaFree(inArg.estimatedAbs);
	cudaFree(inArg_dev);
	
	return;
}

__global__ void kernelZLoopTrans(
		float* outputVol, // [izOut + nzOut * (ixOut + iyOut * nxOut)]
		const constArgsInTrans* inArgs,	
		const constArgKTrans inArgsK
	)
{

	// nt 		--> inArgs->inputGrid->dim[0]
	// nxOut	--> inArgs->outputGrid->dim[1]
	// nz 		--> inArgs->outputGrid->dim[0]
	// zIm 		--> index in output matrix in z direction

	const uint64_t zIm = threadIdx.x + blockIdx.x * blockDim.x;
	if (zIm < inArgs->outputGrid->dim[0])
	{
		// calculate index in output volume
		// iz + ixOut * nz + iyOut * nz * nxOut
		const uint64_t outputIdx = zIm + inArgs->outputGrid->dim[0] * 
			(inArgsK.xIm + inArgsK.yIm * inArgs->outputGrid->dim[1]);  	
		
		// constant offset in model matrix 
		// modMatOff = nt * ir + nnt * nr * iz
		const uint64_t modMatOff = inArgs->inputGrid->dim[0] * (inArgsK.ir + zIm * inArgs->nr);
		
		// constant offset in signal matrix
		// sigMatOff = nt * ix + nt * nx * iy
		const unsigned int sigMatOff = inArgsK.yoff + inArgsK.xoff; // const offset sigMat
		const uint64_t startIdx = inArgs->modelStIdx[2 * (inArgsK.ir + zIm * inArgs->nr)];
		const uint64_t stopIdx = inArgs->modelStIdx[2 * (inArgsK.ir + zIm * inArgs->nr) + 1];

		float temp = 0;

		#pragma unroll 4
		// for (uint64_t it = 0; it < inArgs->inputGrid->dim[0]; it++){
		for (uint64_t it = startIdx; it <= stopIdx; it++){
			temp = __fadd_rn(temp, 
			__fmul_rn(inArgs->sigMat[sigMatOff + it], inArgs->modelMat[modMatOff + it]));
		}

		outputVol[outputIdx] = temp;
		
		// perform fluence correction
		// const unsigned int iSurface = inArgs->surface[inArgsK.xIm + 
		// 	inArgsK.yIm * inArgs->outputGrid->dim[1]];
		// const unsigned int iFluenceIdx = iSurface + 
		// 	inArgs->outputGrid->dim[0] * (inArgsK.ir + zIm * inArgs->nr);
	
		// outputVol[outputIdx] = __fadd_rn(outputVol[outputIdx], __fmul_rn(
		// 		temp, inArgs->fluenceModel[iFluenceIdx]));
	}
	return;
};

__global__ void modelBasedTrans(
	float* outputVol, // [izOut + nzOut * (ixOut + iyOut * nxOut)]
	const constArgsInTrans* inArgs
)
{

	constArgKTrans inArgsK; // struct containing constant kernel arguments
	inArgsK.xIm = blockIdx.x * blockDim.x + threadIdx.x; // x index in output vol
	inArgsK.yIm = blockIdx.y * blockDim.y + threadIdx.y; // y index in output vol

	if ( // check if image voxel actually exists in volume (kernel overlap)
		(inArgsK.xIm < inArgs->outputGrid->dim[1]) && 
		(inArgsK.yIm < inArgs->outputGrid->dim[2])){


		// set current voxel column in output matrix to 0
		const int outputIdx = inArgs->outputGrid->dim[0] * (inArgsK.xIm +	inArgsK.yIm * inArgs->outputGrid->dim[1]);
		for (unsigned int zIm = 0; zIm < inArgs->outputGrid->dim[0]; zIm++)
			outputVol[outputIdx + zIm] = 0;

		// calculate position of current voxel column in COSY of output grid
		const float xPosVox = inArgs->outputGrid->origin[1] + 
			__fmul_rn((float) inArgsK.xIm, inArgs->outputGrid->res[1]);
		const float yPosVox = inArgs->outputGrid->origin[2] + 
			__fmul_rn((float) inArgsK.yIm, inArgs->outputGrid->res[2]);

		// calc reconstruction range in true sized input matrix & convert
		// to index range of inputGrid
		const float rRange = __fmul_rn(inArgs->dr, (float) inArgs->nr);
		
		// calculate minimum index in x position of input grid
		const float xmin = xPosVox - rRange;
		int xminI = __float2int_ru((xmin - inArgs->inputGrid->origin[1]) * 
				inArgs->inputGrid->ires[1]); // round up
		if (xminI < 0)
			xminI = 0;
		
		// calculate maximum index in x position of input grid
		const float xmax = xPosVox + rRange;
		int xmaxI = __float2int_rd((xmax - inArgs->inputGrid->origin[1]) * 
				inArgs->inputGrid->ires[1]); // round down
		if (xmaxI >= inArgs->inputGrid->dim[1])
			xmaxI = inArgs->inputGrid->dim[1] - 1;

		// calculate minimum index in y position of input grid
		const float ymin = yPosVox - rRange;
		int yminI = __float2int_ru((ymin - inArgs->inputGrid->origin[2]) * 
				inArgs->inputGrid->ires[2]); // round up
		if (yminI < 0)
			yminI = 0;

		// calculate maximum index in y position of input grid
		const float ymax = yPosVox + rRange;
		int ymaxI = __float2int_rd((ymax - inArgs->inputGrid->origin[2]) * 
			inArgs->inputGrid->ires[2]); // round down
		if (ymaxI >= inArgs->inputGrid->dim[2])
			ymaxI = inArgs->inputGrid->dim[2] - 1;

		
		// check for range validity
		if (xminI > xmaxI)
		{
			printf("Corrupt x range: %d ... %d\n", xminI, xmaxI);
		}		

		if (yminI > ymaxI)
		{
			printf("Corrupt y range: %d ... %d\n", yminI, ymaxI);
		}


		// run over each a scan
		float rDist; // radial distance between A scan and voxel
		float deltaX; // x distance between A scan and voxel
		float deltaY; // y distance between A scan and voxel
		float xPosAScan;
		float yPosAScan;

		const dim3 blockSize(32);
		const dim3 gridSize((inArgs->outputGrid->dim[0] + blockSize.x - 1) / blockSize.x);

		// run along neighbouring elements in y direction (inputGrid)
		#pragma unroll
		for (int yidx = yminI; yidx < ymaxI; yidx++){
			// index offset for input volume along y
			// yOffset = yAScanIdx * nt * nx
			inArgsK.yoff = yidx * inArgs->inputGrid->dim[0] * inArgs->inputGrid->dim[1]; 
			
			yPosAScan = inArgs->inputGrid->origin[2] + ((float) yidx) * inArgs->inputGrid->res[2];
			deltaY = yPosVox - yPosAScan;
			
			// run along neighbouring elements in x direction (inputGrid)
			#pragma unroll
			for (int xidx = xminI; xidx < xmaxI; xidx++){
				// index offset for input volume along x
				// xOffset = xAScanIdx * nt
				inArgsK.xoff = xidx * inArgs->inputGrid->dim[0];
				
				xPosAScan = inArgs->inputGrid->origin[1] + ((float) xidx) * inArgs->inputGrid->res[1];
				deltaX = xPosVox - xPosAScan;

				// calculate radial distance
				rDist = __fsqrt_rn(__fadd_rn(
						__fmul_rn(deltaX, deltaX), __fmul_rn(deltaY, deltaY)));
				// printf("%f\n", rDist);

				// convert radial distance into index
				inArgsK.ir = rDist / inArgs->dr + 0.5;
				
				// multiply over full z vector
				if (inArgsK.ir < inArgs->nr)
				{
					// inArgsK.ir = inArgs->nr - inArgsK.ir - 1;
					// here we could implement dynamic parallelism over z positions
					kernelZLoopTrans<<<gridSize, blockSize>>>(
						outputVol, // matrix containing reconstructed absorber map
						inArgs, // global constant arguments
						inArgsK // local constant arguments
					);
					cudaDeviceSynchronize();

				}
			}
		}
	}
	return;
}

// transpose operation takes an input matrix of size dimCropped and remodels it into the
// shape of our absorber matrix [iz, ix,  iy]
void mbrecon::transpose(float* absorberMatrix, const float* sigMat)
{
	cudaError_t err; // error container for cuda
	
	float* outputVol_dev; // array containing output volume on GPU [iz, ix, iy]
	constArgsInTrans* inArgs_dev; // constant arguments struct for device
	constArgsInTrans inArgs; // constant argument struct on host

	vprintf("Allocating memory on GPU and copy things over... ", 1);
	
	// signal matrix [it, ix, iy]
	err =  cudaMalloc( (void**) &inArgs.sigMat,	croppedVol.nElements * sizeof(float));
	checkCudaReturn(err, "Could not allocate signal matrix");
	
	err = cudaMemcpy(inArgs.sigMat, sigMat, croppedVol.nElements * sizeof(float), 
		cudaMemcpyHostToDevice);
	checkCudaReturn(err, "Could not copy signal matrix to device");
	
	// output matrix [iz, ix, iy]
	err = cudaMalloc( (void**) &outputVol_dev, estimAbs.nElements * sizeof(float));
	checkCudaReturn(err, "Could not allocate absorber matrix");
	
	// memory for model matrix [it, ir, iz]
	err = cudaMalloc( (void**) &inArgs.modelMat, model->getNElements() * sizeof(float));
	checkCudaReturn(err, "Could not allocate model matrix");
	
	err = cudaMemcpy(inArgs.modelMat, model->data, model->getNElements() * sizeof(float), 
		cudaMemcpyHostToDevice);
	checkCudaReturn(err, "Could not copy model matrix to device");
	
	// input grid definition based of croppedVol
	err = cudaMalloc( (void**) &inArgs.inputGrid, sizeof(griddedData));
	checkCudaReturn(err, "Could not allocate input grid");
	
	err = cudaMemcpy(inArgs.inputGrid, &croppedVol, sizeof(griddedData), 
		cudaMemcpyHostToDevice);
	checkCudaReturn(err, "Could not copy input grid to device");
	
	// output grid definition based on estimated absorbers
	err = cudaMalloc( (void**) &inArgs.outputGrid, sizeof(griddedData));
	checkCudaReturn(err, "Could not allocate output grid");
	
	err = cudaMemcpy(inArgs.outputGrid, &estimAbs, sizeof(griddedData), 
		cudaMemcpyHostToDevice);
	checkCudaReturn(err, "Could not copy output grid to device");
	
	// start and stop index in model matrix
	err = cudaMalloc( (void**) &inArgs.modelStIdx, 2 * estimAbs.dim[0] * nr * 
		sizeof(uint64_t) );
	checkCudaReturn(err, "Could not allocate model start index");

	err = cudaMemcpy(inArgs.modelStIdx, model->startStopIdx, 2 * estimAbs.dim[0] * nr * 
		sizeof(uint64_t), cudaMemcpyHostToDevice);
	checkCudaReturn(err, "Could not copy model start index to device");
	
	inArgs.nr = nr; // push number of elements in radial direction over
	inArgs.dr = sett.get_rRes(); // push resolution in radial direction (both fluence and model)

	// allocate memory for struct containing everything
	err = cudaMalloc( (void**) &inArgs_dev, sizeof(constArgsInTrans));
	checkCudaReturn(err, "Could not allocate const arguments");
	err = cudaMemcpy(inArgs_dev, &inArgs, sizeof(constArgsInTrans), cudaMemcpyHostToDevice);
	checkCudaReturn(err, "Could not copy const arguments to device");

	vprintf("done!\n");

	// generate processing grid
	dim3 blockSize(8, 1); // prev 4, 4, 8
	dim3 gridSize(
		(estimAbs.dim[1] + blockSize.x - 1) / blockSize.x,
		(estimAbs.dim[2] + blockSize.y - 1) / blockSize.y);
	
	if (flagDebug)
	{
		printf("[debug] Signal matrix size (input): nt = %d, nx = %d, ny = %d\n", 
			croppedVol.dim[0], croppedVol.dim[1], croppedVol.dim[2]);
		printf("[debug] Model matrix size: nt = %d, nr = %d, nz = %d\n", 
			model->get_nt(), model->get_nr(), model->get_nz());
		printf("[debug] Estimated absorber size (output): nz = %d, nx = %d, ny = %d\n", 
			estimAbs.dim[0], estimAbs.dim[1], estimAbs.dim[2]);
		printf("[debug] Kernel grid size: %d x %d\n", gridSize.x, gridSize.y);
		printf("[debug] nr = %d, dr = %f\n", inArgs.nr, inArgs.dr);
	}

	// check for potential dimension mismatchs
	if (croppedVol.dim[0] != model->get_nt())
	{
		printf("Dimension mismatch between model matrix and signal matrix\n");
		throw "dimError";
	}

	if (estimAbs.dim[0] != model->get_nz())
	{
		printf("Dimension mismatch between model matrix and absorber matrix\n");
		throw "dimError";
	}		

	vprintf("Starting transpose kernel... ", 1); fflush(stdout);
	clock_t tStart = clock();
	modelBasedTrans<<< gridSize, blockSize >>>(outputVol_dev, inArgs_dev);
	
	cudaDeviceSynchronize(); // wait for kernel to finish
	vprintf("done!\n", 0);
	printf("Time taken: %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);
	
	err = cudaGetLastError();
	checkCudaReturn(err, "Error during kernel execution");
	
	err = cudaMemcpy(absorberMatrix, outputVol_dev, estimAbs.nElements * sizeof(float), 
		cudaMemcpyDeviceToHost);
	checkCudaReturn(err, "Could not copy back memory from device");

	// free device and host memory
	cudaFree(inArgs.modelStIdx);
	cudaFree(inArgs.inputGrid);
	cudaFree(inArgs.outputGrid);
	cudaFree(inArgs.sigMat);
	cudaFree(inArgs.modelMat);
	cudaFree(outputVol_dev);
	cudaFree(inArgs_dev);
	return;
}

void mbrecon::export_data()
{
	// if wanted export as *.vtk
	if (sett.get_flagVTKExport())
	{

	}

	// if wanted export as *.h5
	if (sett.get_flagH5Export())
	{

	}

	return;
}
