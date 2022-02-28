#include <cuda.h>
#include <cuda_runtime.h>
#include <cstdint>
#include <iostream>
#include <chrono>

using namespace std;

#ifndef TRAN_KERNEL_AG_H
#define TRAN_KERNEL_AG_H

struct tran_kernel_arg
{
	uint64_t nxScan;
	uint64_t nyScan;
	uint64_t nz;
	uint64_t nt;
	uint64_t nxModel;
	uint64_t nyModel;
	float* sigMat;
	float* modelMat; 
};

#endif

#ifndef FWD_KERNEL_AG_H
#define FWD_KERNEL_AG_H

struct fwd_kernel_arg
{
	uint64_t nxScan;
	uint64_t nyScan;
	uint64_t nz;
	uint64_t nt;
	uint64_t nxModel;
	uint64_t nyModel;
	float* absMat;
	float* modelMat; 
};

#endif

#ifndef MBR_TRANSSYM_H
#define MBR_TRANSSYM_H

class mbr_transsym
{
private:
	uint64_t nxScan = 0; // 0 means not set
	uint64_t nyScan = 0; // 0 means not set
	uint64_t nz = 0; // 0 means not set
	uint64_t nt = 0; // 0 means not set
	uint64_t nxModel = 0; // 0 means not set
	uint64_t nyModel = 0; // 0 means not set

	uint64_t tExec = 0; // execution time of last kernel run

	// pointers to host arrays
	float* modelMat; // not owned, order: x, y, t, z
	float* sigMat; // not owned
	float* absMat; // not owned

	// pointer to device arrays
	float* modelMat_dev;
	float* absMat_dev;
	float* sigMat_dev;

	cudaError_t cudaErr;
	void cudaErrCheck(const char* errMsg);
	uint64_t get_nBytesModel() const {return nxModel * nyModel * nz * nt * sizeof(float);};
	uint64_t get_nBytesSignal() const {return nt * nxScan * nyScan * sizeof(float);};
	uint64_t get_nBytesAbs() const {return  nz * nxScan * nyScan * sizeof(float);};

public:
	~mbr_transsym();

	void set_sizeSigMat(const uint64_t _nxScan, const uint64_t _nyScan, const uint64_t _ntScan);
	void set_sizeAbsMat(const uint64_t _nxAbs, const uint64_t _nyAbs, const uint64_t _nzAbs);
	void set_sizeModel(
		const uint64_t _nxModel, const uint64_t _nyModel, 
		const uint64_t _ntModel, const uint64_t _nz);
	void set_modelMat(float* _modelMat) {modelMat = _modelMat; return;};
	void set_sigMat(float* _sigMat) {sigMat = _sigMat; return;};
	void set_absMat(float* _absMat) {absMat = _absMat; return;};

	void print_info();

	// return number of elements in matrixes
	uint64_t get_nElementsAbs() const {return nxScan * nyScan * nz;};
	uint64_t get_nElementsSig() const {return nxScan * nyScan * nt;};
	uint64_t get_nElementsMod() const {return nxModel * nyModel * nt * nz;};
	
	// return size of stuff
	uint64_t get_nz() const {return nz;};
	uint64_t get_nt() const {return nt;};
	uint64_t get_nxScan() const {return nxScan;};
	uint64_t get_nyScan() const {return nyScan;};
	uint64_t get_nxModel() const {return nxModel;};
	uint64_t get_nyModel() const {return nyModel;};

	uint64_t get_tExec() const {return tExec;};

	void run_trans();
	void run_fwd();
};

#endif