#include "recon.h"


void recon::crop()
{
	statusVerbal = "Cropping signal matrix...";
	// for now we only "crop" along time direction
	const uint64_t tStartVol = sigMat.get_idx(sett.get_tMin(), 2);
	const uint64_t tStopVol = sigMat.get_idx(sett.get_tMax(), 2);
	const uint64_t nt = tStopVol - tStartVol + 1;
	
	const uint64_t startIdxSig[3] = {0, 0, tStartVol};
	const uint64_t stopIdxSig[3] = {sigMat.get_dim(0) - 1, sigMat.get_dim(1) - 1, tStopVol};

	// crop signal matrix to correct t range
	croppedSigMat = sigMat;
	croppedSigMat.crop(startIdxSig, stopIdxSig);

	statusVerbal = "Cropping model matrix...";
	// order t, x, y, z
	const uint64_t tStartMod = mod.get_tIdx(sett.get_tMin());
	const uint64_t tStopMod = tStartMod + nt - 1;
	const uint64_t startIdxMod[4] = {tStartMod, 0, 0, 0};
	const uint64_t stopIdxMod[4] = {tStopMod, mod.get_nx() - 1, mod.get_ny() - 1, mod.get_nz() - 1};

	// crop model matrix to correct range
	croppedMod = mod;
	croppedMod.crop(startIdxMod, stopIdxMod);

	return;
}

void recon::lsqr()
{
	statusVerbal = "Preparing variables for LSQR...";
	// here we load our mbr_transsym with all the important details
	absMat.set_dim(croppedSigMat.get_dim(0), croppedSigMat.get_dim(1), croppedMod.get_nz());
	absMat.alloc_memory();

	kernel.set_sizeAbsMat(absMat.get_dim(0), absMat.get_dim(1), absMat.get_dim(2)); // nx, ny, nz
	kernel.set_sizeSigMat(
		croppedSigMat.get_dim(0), croppedSigMat.get_dim(1), croppedSigMat.get_dim(2)); // nx, ny, nt
	kernel.set_sizeModel(croppedMod.get_nx(), croppedMod.get_ny(), croppedMod.get_nt(), croppedMod.get_nz());
	kernel.set_modelMat(croppedMod.get_pdata());

	volume r(absMat.get_dim(0), absMat.get_dim(1), absMat.get_dim(2));
	
	// prepare reconstruction matrix including origin and resolution
	volume recon(absMat.get_dim(0), absMat.get_dim(1), absMat.get_dim(2));
	recon.set_res(sigMat.get_res(0), sigMat.get_res(1), croppedMod.get_dz());
	recon.set_origin(0.0f, 0.0f, croppedMod.get_z0());

	volume v(absMat.get_dim(0), absMat.get_dim(1), absMat.get_dim(2));
	volume w(absMat.get_dim(0), absMat.get_dim(1), absMat.get_dim(2));
	volume wWeighted(absMat.get_dim(0), absMat.get_dim(1), absMat.get_dim(2));
	// float* r = new float [estimAbs.nElements]; // r has dimension [nx, ny, nz]
	volume p(croppedSigMat.get_dim(0), croppedSigMat.get_dim(1), croppedSigMat.get_dim(2));
	volume u(croppedSigMat.get_dim(0), croppedSigMat.get_dim(1), croppedSigMat.get_dim(2));
	
	// normalize signal matrix
	u = croppedSigMat;
	float beta = u.get_norm();
	u /= beta;
	
	// r = A' * croppedVol.data
	statusVerbal = "Running first transpose multiplication...";
	kernel.set_sigMat(u.get_pdata()); // set input pointer
	kernel.set_absMat(r.get_pdata()); // output pointer
	kernel.run_trans(); // run transpose multiplication

	recon = 0.0f; // set absorbance matrix initially to all 0
	float alpha = r.get_norm();
	
	v = r; v /= alpha; // v = r / alpa

	float phi_bar = beta;
	float rho_bar = alpha;

	w = v;

	for (uint8_t iIter = 0; iIter < sett.get_nIter(); iIter++)
	{
		statusVerbal = "Running iteration " + std::to_string(iIter + 1) + " of " + 
			std::to_string(sett.get_nIter());
		
		// bidiagonalization 
		// p = A * v - alpha * u
		// p = A * v;
		kernel.set_absMat(v.get_pdata()); // set input pointer
		kernel.set_sigMat(p.get_pdata()); // set output pointer
		kernel.run_fwd(); // afterwards: p = A * v

		u *= alpha; // u = u * alpha
		p -= u; // p = p - u (--> p = A * v - u * alpha)

		// determine beta
		// beta = sqrt(beta^2 + (norm(lambda .* v(:)))^2);
		// note: norm(lambda .* v(:)) is the same as lambda * norm(v(:))
		beta = p.get_norm();
		// const float vNorm = v.get_norm();
		// const float vNormWeighted = vNorm * sett.get_regStrength();
		// beta = powf(beta * beta + vNormWeighted * vNormWeighted, 0.5f);
		
		u = p; 	u /= beta; // u = p / beta

		// r = A' * u - beta * v;
		kernel.set_sigMat(u.get_pdata()); // set input pointer
		kernel.set_absMat(r.get_pdata()); // set output pointer
		kernel.run_trans(); // r = A' * u

		wWeighted = v;
		wWeighted *= (-1.0f * beta);
		// wWeighted = v * lambda^2 / norm_p

		r += wWeighted; // r = r + wWeighted
		
		// r = ATu - norm_p .* v
		// v *= beta; // v = v * norm_p
		// r -= v; // r = r - v

		alpha = r.get_norm(); // alpha = norm(r);
		v = r; v /= alpha; // v = r / alpha
		
		// % ------------- orthogonal transformation ---------------
		// rrho = norm([rho_bar, beta]);
		const float rrho = powf(rho_bar * rho_bar + beta * beta, 0.5f);
		const float c1 = rho_bar / rrho;
		const float s1 = beta / rrho;
		const float theta = s1 * alpha;
		rho_bar = -c1 * alpha;
		const float phi = c1 * phi_bar;
		phi_bar = s1 * phi_bar;

		printf("Phi bar after iteration %d: %f\n", iIter + 1, phi_bar);
		phi_bar_vec.push_back(phi_bar);

		// update solution and search direction
		// x = x + (phi / rrho) * w
		wWeighted = w;
		wWeighted *= (phi / rrho); // wWeighted = w * phi / rrho
		recon += wWeighted;

		// w = v - (theta / rrho) * w;
		wWeighted = w;
		wWeighted *= (theta / rrho);
		w = v; w -= wWeighted;
		
		if (flagAbort)
			break;

		// update absMat after each iteration to preview stuff in GUI
  	absMat = recon; 
  	absMat.calcMinMax();

	}

	return;
}

// main reconstruction loop. this will later allow distinguishing between different
// inversion schemes
void recon::reconstruct()
{
	flagAbort = 0;
	isRunning = 1;
	phi_bar_vec.clear();

	tStart = std::chrono::system_clock::now();
	crop();
	lsqr();
	absMat.calcMinMax();
	tEnd = std::chrono::system_clock::now();
	reconTime = duration_cast<seconds>(tEnd - tStart).count();
	
	statusVerbal = "reconstruction done";
	isRunning = 0;
	isRecon = 1;

	return;
}

// this should run our reconstruction in a detached thread
std::thread recon::reconstruct2thread()
{
	isRunning = 1;
	return std::thread([=] { reconstruct();} );
}

void recon::set_flagAbort(const bool _flagAbort)
{
	flagAbort = _flagAbort;
	return;
}