#include "../src/reconsett.h"
#include "../src/model.h"
#include "../src/recon.h"


int main()
{

	// load model matrix from file
	recon rec;

	reconsett* sett = rec.get_psett();
	sett->set_regStrength(5.0f);
	sett->set_nIter(10);
	sett->set_tMin(4e-6f);
	sett->set_tMax(7.34e-6f);

	model* mod = rec.get_pmodel();
	mod->load_from_file("/mnt/hofmannu/aroam/99_CbrPaper/models/exp_model.h5");
	// mod->print_information();

	volume* sigMat = rec.get_psigMat();
	sigMat->readFromFile("/mnt/hofmannu/aroam/06_HumanSkin/2021-10-20_willWrist/011_will_wrist_532_preproc.h5");

	rec.crop();
	rec.lsqr();

	return 0;
}