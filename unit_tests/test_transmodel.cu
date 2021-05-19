#include <stdio.h>
#include "../src/transModel.cuh" 
#include "../src/transducerProperties.h"
#include "../src/fieldProperties.h"

using namespace std;

int main()
{

	printf("Running unit test for transducer model class...\n");

	// transModel testTrans;
	// testTrans.setName("testTransducer");
	// testTrans.setFocalDistance(7e-3);
	// testTrans.setRAperture(3.2e-3);
	// testTrans.setRHole(0.45e-3);
	// testTrans.setCentralFrequency(60e6);
	// testTrans.saveToFile();

	// transModel testTrans2;
	// testTrans2.setName("testTransducer");

	transducerProperties transProp;
	transProp.setName("russian_01");
	transProp.readFromFile();
	transProp.printProperties();

	fieldProperties fieldProp;
	fieldProp.set_dr(10e-6);
	fieldProp.set_dz(10e-6);
	fieldProp.set_nr(3e-3 / 10e-6);
	fieldProp.set_z0(-3e-3);
	fieldProp.set_SOS(1495);
	fieldProp.set_dSir(1495 / 250e6 / 4);
	unsigned int nz = (3e-3 / 10e-6) * 2.0 + 0.5;
	fieldProp.set_nz(nz);

	transModel testTran3;
	testTran3.setFieldProperties(fieldProp);
	testTran3.setTransProp(transProp);
	testTran3.setNElements(1e5);

	testTran3.loadImpulseResponse();
	testTran3.buildGPUModel();

	modelMatrix* myModel = testTran3.getModel();
	myModel->exportSensFieldVtk();
	myModel->saveToFile("/home/hofmannu/modelExample.h5");

	return 0;
}