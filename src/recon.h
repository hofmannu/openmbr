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

#ifndef RECON_H
#define RECON_H

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
	uint8_t* iterator; // which iteration is the reconstruction running at?

public:
	reconsett* get_psett() {return &sett;};
	model* get_pmodel() {return &mod;};
	volume* get_psigMat() {return &sigMat;};
	volume* get_pabsMat() {return &absMat;};

	void lsqr();
	void crop();

};

#endif