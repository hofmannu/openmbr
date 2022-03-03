#include "../src/reconsett.h"
#include "../src/model.h"
#include "../src/recon.h"


int main()
{
	// load model matrix from file
	recon rec;

	reconsett* sett = rec.get_psett();
	sett->set_regStrength(5.0f);
	sett->set_nIter(1);
	sett->set_tMin(4e-6f);
	sett->set_tMax(7.34e-6f);

	model* mod = rec.get_pmodel();
	mod->load_from_file("/mnt/hofmannu/aroam/99_CbrPaper/models/exp_model.h5");

	volume* sigMat = rec.get_psigMat();
	sigMat->readFromFile("/mnt/hofmannu/aroam/06_HumanSkin/2021-10-20_willWrist/011_will_wrist_532_preproc.h5");

	rec.reconstruct();

	volume* absMat = rec.get_pabsMat();
	absMat->saveToFile("/mnt/hofmannu/aroam/06_HumanSkin/2021-10-20_willWrist/011_will_wrist_532_recon.h5");

	return 0;
}