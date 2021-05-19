#ifndef FIBERPROPERTIES_H
#define FIBERPROPERTIES_H

#include <cstring>
#include <iostream>
#include <cmath>

using namespace std;

class fiberProperties
{
private:
	float numAp = 0.1; // numerical aperture of fiber
	float dCore = 0.2e-3; // core diameter of fiber
	float rCore = 0.1e-3;
	string name;
public:
	
	fiberProperties();

	void set_numAp(const float _numAp);
	void set_dCore(const float _dCore);
	void set_rCore(const float _rCore);
	void set_name(const string _name);

	float get_numAp() const;
	float get_dCore() const;
	float get_rCore() const;
	float get_theta(const float nMedium) const; // one sided opening angle
	float get_rSpot(const float nMedium, const float dist) const;
	string get_name() const;

	// not implemented yet
	void saveToFile();
	void saveToFile(const string fileName);
	void readFromFile();
	void readFromFile(const string fileName);
};

#endif