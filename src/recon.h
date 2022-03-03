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
	uint8_t* iterator; // which iteration is the reconstruction running at?

	time_point<system_clock> tStart; // start time of reconstruction
	time_point<system_clock> tEnd; // end time of reconstruction
	double tRemain = 0; // remaining time for execution
	float reconTime = 0; // time required for last reconstruction
	string statusVerbal = "not started";

public:
	reconsett* get_psett() {return &sett;};
	model* get_pmodel() {return &mod;};
	volume* get_psigMat() {return &sigMat;};
	volume* get_pabsMat() {return &absMat;};

	std::thread reconstruct2thread();
	void reconstruct();
	void lsqr();
	void crop();

	bool get_isRecon() const {return isRecon;};
	bool get_isRunning() const {return isRunning;};

	time_point<system_clock> get_tStart() const {return tStart;};

	float get_reconTime() const {return reconTime;};
	const char* get_statusVerbal() const {return statusVerbal.c_str();};
};

#endif