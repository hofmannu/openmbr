#include "fiberProperties.h"

fiberProperties::fiberProperties()
{
		name = "default";

}

float fiberProperties::get_numAp() const {return numAp;}
float fiberProperties::get_dCore() const {return dCore;}
float fiberProperties::get_rCore() const {return rCore;}
string fiberProperties::get_name() const {return name;}

void fiberProperties::set_numAp(const float _numAp)
{
	numAp = _numAp;
	return;
}

void fiberProperties::set_dCore(const float _dCore)
{
	dCore = _dCore;
	rCore = _dCore / 2;
	return;
}

void fiberProperties::set_rCore(const float _rCore)
{
	dCore = 2 * _rCore;
	rCore = _rCore;
	return;
}

void fiberProperties::set_name(const string _name)
{
	name = _name;
	return;
}

float fiberProperties::get_theta(const float nMedium) const
{
	float theta = asin(numAp / nMedium);
	return theta;
}

float fiberProperties::get_rSpot(const float nMedium, const float dist) const
{
	float rSpot = rCore + tan(get_theta(nMedium)) * dist;
	return rSpot;
}