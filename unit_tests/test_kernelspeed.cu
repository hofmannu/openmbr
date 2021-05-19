/* 
	small standalone script to test the speed of the forward and transpose kernel
	Author: Urs Hofmann
	Mail: hofmannu@ethz.ch
	Date: 16.10.2020

	Starting point for optimization of forward kernel: 9.06 s
	Using temporary variable for calcSig: 5.59 s
	Serializing for loops 
*/

#include "../src/reconSettings.h"
#include "../src/mbrecon.cuh"
#include "../src/transducerProperties.h"
#include "../src/volume.h"

using namespace std;

int main()
{

	mbrecon recon;

	// reconstruction settings
	reconSettings sett;
	sett.set_SOS(1490.0);
	sett.set_drRecon(25e-6, 25e-6, 25e-6);
	sett.set_xLim(-2e-3, 2e-3);
	sett.set_yLim(-2e-3, 2e-3);
	sett.set_zLim(6.5e-3, 7.5e-3);
	sett.set_rRes(1e-6);
	sett.set_threshold(0.0);
	recon.set_sett(sett);

	// transducer properties
	transducerProperties* trans;
	trans = recon.get_ptransProp();
	trans->setFocalDistance(7e-3);
	trans->setRAperture(3e-3);
	trans->setRHole(0.0);
	trans->setCentralFrequency(60e6);
	trans->setName("russian_johanna");

	// absorber matrix
	volume* absorberMatrix;
	absorberMatrix = recon.get_preconVol();

	// set default values for absorber volume
	for (uint8_t iDim = 0; iDim < 3; iDim++)
	{
		absorberMatrix->setRes(iDim, sett.get_dr(iDim));
		absorberMatrix->setDim(iDim, sett.get_dim(iDim));
	}
	
	absorberMatrix->setOrigin(0, sett.get_zLower());
	absorberMatrix->setOrigin(1, sett.get_xLower());
	absorberMatrix->setOrigin(2, sett.get_yLower());

	absorberMatrix->allocMemory();
	absorberMatrix->assignRand(absorberMatrix->get_pdata(), 
		absorberMatrix->get_nElements());

	// forward kernel is tested with simulate
	recon.simulate();

	return 0;
}