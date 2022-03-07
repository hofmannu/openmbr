/* 
	File: recon.h
	Author: urs Hofmann
	Mail: mail@hofmannu.org
	Date: 19.02.2022
*/

#include "reconsett.h"
#include "model.h"
#include "../lib/CVolume/src/volume.h"
#include "mbr_transsym.cuh"
#include <string>
#include <stdio.h>
#include <math.h>
#include <vector>

#ifndef RECON_H
#define RECON_H

using namespace std::chrono;

class recon
{
private:

	volume sigMat;
	volume croppedSigMat;
	volume absMat;
	reconsett sett;
	model mod;
	model croppedMod;
	mbr_transsym kernel;

	bool isRecon = 0; // incicator if we are done reconstructing
	bool isRunning = 0; // indicator if a reconstruction is currently runningS
	bool flagAbort = 0; // did we receive an abort request

	uint8_t* iterator; // which iteration is the reconstruction running at?

	// reconstruction time
	time_point<system_clock> tStart; // start time of reconstruction
	time_point<system_clock> tEnd; // end time of reconstruction
	double tRemain = 0; // remaining time for execution
	float reconTime = 0; // time required for last reconstruction

	string statusVerbal = "not started";
	vector<float> phi_bar_vec;

public:
	reconsett* get_psett() {return &sett;}; // return pointer to settings
	model* get_pmodel() {return &mod;}; // return pointer to model matrix
	volume* get_psigMat() {return &sigMat;}; // return pointer to volume containing signals
	volume* get_pabsMat() {return &absMat;}; // return pointer to volume containing absorption

	std::thread reconstruct2thread(); // return reconstruction as a detatched thread
	void reconstruct();
	void lsqr();
	void crop(); // crop both signal and model along t

	// some status indicators
	bool get_isRecon() const {return isRecon;};
	bool get_isRunning() const {return isRunning;};
	void set_flagAbort(const bool _flagAbort); 
	bool get_flagAbort() const {return flagAbort;}; 
	int get_iIter() const {return phi_bar_vec.size() + 1;};
	float get_relConvergence() const {return (phi_bar_vec.back() / phi_bar_vec.front() * 100.0f);}

	time_point<system_clock> get_tStart() const {return tStart;};

	float get_reconTime() const {return reconTime;};
	const char* get_statusVerbal() const {return statusVerbal.c_str();};

};

#endif