/*
	c++ class representing properties of a spherically focused transducer
	Author: Urs Hofmann
	Mail: hofmannu@ethz.ch
	Date: 14.08.2020
*/

#ifndef TRANSDUCERPROPERTIES_H
#define TRANSDUCERPROPERTIES_H

#include <cmath>
#include <iostream>
#include <H5Cpp.h>
#include <fstream>
#include "baseClass.h"

using namespace std;

class transducerProperties : public baseClass
{
private:
	string name = "russian_pavel"; // unique transducer identifier / name
	
	// path to impulse response file
	string pathImpResp;
	bool isPathImpRespDefined = 0;
	
	float focalDistance = 7e-3; // focal distance of transducer [m]
	float rAperture = 3.2e-3; // radius of aperture of transducer [m]
	float centralFrequency = 60e6; // central frequency of transducer [Hz]
	float rHole = 0.5e-3; // hole radius of transducer [m]
	// constant get functions
public:
	// static get functions
	float getFocalDistance() const;
	float* get_pfocalDistance() {return &focalDistance;};

	float getRAperture() const;
	float* get_prAperture() {return &rAperture;};

	float getCentralFrequency() const;
	float* get_pcentralFrequency() {return &centralFrequency;};

	float getRHole() const;
	float* get_prHole() {return &rHole;};
	string getName() const;

	// dependent get functions
	float getTheta() const;
	float getThetaHole() const;

	// set functions
	void setFocalDistance(const float _focalDistance);
	void setRAperture(const float _rAperture);
	void setCentralFrequency(const float _centralFrequency);
	void setRHole(const float _rHole);
	void setName(const string _name);

	// save properties / load properties to / from file
	void readFromFile();
	void readFromFile(const string filePath);
	void saveToFile();
	void saveToFile(const string filePath);

	void printProperties();
	void defineNew();
};

#endif
