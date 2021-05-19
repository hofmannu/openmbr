/*
	allows definition and saving of a new transducer to a file
	Author:  Urs Hofmann
	Mail: hofmannu@biomed.ee.ethz.ch
	Date: 28.08.2020
*/

#include "transducerProperties.h"
#include <iostream>

using namespace std;

int main()
{

	cout << ">>> Transducer definition script <<<" << endl;
	transducerProperties myTran;

	int happy = 0;
	do{
		myTran.defineNew();
		myTran.printProperties();
		cout << "Correct result? [0 / 1]" << endl;
		cin >> happy;
	}while(happy == 0);
	
	myTran.saveToFile();
	return 0;
}
