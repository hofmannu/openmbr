/* 
	builds a transducer model and saves it to the disk as an h5 file
	Author: Urs Hofmann
	Mail: hofmannu@biomed.ee.ethz.ch
	Date: 27.08.2020
*/

#include "transModel.cuh"
#include "modelMatrix.h"
#include "transducerProperties.h"
#include "fieldProperties.h"
#include <iostream>

using namespace std;

int main()
{

	int selection = 1;
	transModel model;
	transducerProperties prop;
	prop.setName("russian_01");
	prop.readFromFile();
	string transducerName;
	fieldProperties field;

	do{
		cout << ">>> Welcome to model building helper script <<<" << endl;
		cout << " [1] Select transducer" << endl;
		cout << " [2] Build model" << endl;
		cout << " [3] Save model to file" << endl;
		cout << " [0] Exit " << endl;
		cin >> selection;

		if (selection == 0)
		{

		}
		else if (selection == 1)
		{
			cout << "Enter transducer name: ";
			cin >> transducerName;
			prop.setName(transducerName);
			prop.readFromFile();
		}
		else if (selection == 2)
		{
			model.setTransProp(prop);
			model.loadImpulseResponse();
			model.buildGPUModel();
		}
		else if (selection == 3)
		{
			// save model matrix to file
			modelMatrix* modmatrix;
			modmatrix = model.getModel();
			modmatrix->saveToFile("/home/hofmannu/modMatrix2.h5"); 
		}

	} while(selection != 0);


	return 0;
}