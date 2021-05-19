/*
	small test script which runs execution of model based reconstruction on a test dataset
	Author: Urs Hofmann
	Mail: hofmannu@biomed.ee.ethz.ch
	Date: 28.08.2020
*/


#include "../src/mbrecon.cuh"

using namespace std;

int main()
{
	
	mbrecon myRecon;
	myRecon.load_data("/home/hofmannu/pavelsRabbit.h5");
	myRecon.recon();
	myRecon.export_vtk();

	return 0;
}
