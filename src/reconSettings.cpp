#include "reconSettings.h"

string reconSettings::get_polarityHandling() const {return polarityHandling;}
string reconSettings::get_regMethod() const {return regMethod;}
bool reconSettings::get_flagRegularization() const {return flagRegularization;}
float reconSettings::get_SOS() const {return SOS;}
float reconSettings::get_lambdaReg() const {return lambdaReg;}

bool reconSettings::isRegMethod(const string _regMethod) const
{
	bool result = _regMethod.compare(regMethod) ? 1 : 0;
	return result;
}


void reconSettings::set_SOS(const float _SOS)
{
	if (_SOS <= 0)
		throw "Speed of sound must be bigger then 0";
	else
		SOS = _SOS;

	return;
}


void reconSettings::set_lambdaReg(const float _lambdaReg)
{
	lambdaReg = _lambdaReg;
	return;
}


void reconSettings::set_threshold(const float _threshold)
{
	if (_threshold < 0)
	{
		printf("Threshold must be bigger or equal 0");
		throw "invalidValue";
	}
	else
	{
		threshold = _threshold;
	}

	return;
}

float reconSettings::get_threshold() const {return threshold;}

void reconSettings::set_nIter(const uint8_t _nIter)
{
	if (_nIter <= 0)
		throw "Number of iterations must be bigger then 0";
	else
		nIter = _nIter;
}

// define boundaries along z dimension
void reconSettings::set_zLower(const float _zLower)
{
	boundRecon[0] = _zLower;
	return;
}

void reconSettings::set_zUpper(const float _zUpper)
{
	boundRecon[1] = _zUpper;
	return;
}

void reconSettings::set_zLim(const float _zLower, const float _zUpper)
{
	set_zLower(_zLower);
	set_zUpper(_zUpper);
	return;
}

// define boudaries aling x dimension
void reconSettings::set_xLower(const float _xLower)
{
	boundRecon[2] = _xLower;
	return;
}

void reconSettings::set_xUpper(const float _xUpper)
{
	boundRecon[3] = _xUpper;
	return;
}

void reconSettings::set_xLim(const float _xLower, const float _xUpper)
{
	set_xLower(_xLower);
	set_xUpper(_xUpper);
	return;
}

// define reconstruction boundaries along y dimension
void reconSettings::set_yLower(const float _yLower)
{
	boundRecon[4] = _yLower;
	return;
}

void reconSettings::set_yUpper(const float _yUpper)
{
	boundRecon[5] = _yUpper;
	return;
}

void reconSettings::set_yLim(const float _yLower, const float _yUpper)
{
	set_yLower(_yLower);
	set_yUpper(_yUpper);
	return;
}

uint8_t reconSettings::get_nIter() const {return nIter;}
float* reconSettings::get_zRecon() {return &boundRecon[0];}
float reconSettings::get_zLower() const {return boundRecon[0];}
float reconSettings::get_zUpper() const {return boundRecon[1];}
float* reconSettings::get_xRecon() {return &boundRecon[2];}
float reconSettings::get_xLower() const {return boundRecon[2];}
float reconSettings::get_xUpper() const {return boundRecon[3];}
float* reconSettings::get_yRecon() {return &boundRecon[4];}
float reconSettings::get_yLower() const {return boundRecon[4];}
float reconSettings::get_yUpper() const {return boundRecon[5];}
bool reconSettings::get_flagVTKExport() const {return flagVTKExport;}
bool reconSettings::get_flagH5Export() const {return flagH5Export;}

unsigned int* reconSettings::get_downsampling() {return &downsampling[0];}
unsigned int reconSettings::get_downsampling(const unsigned int _dim) const 
{
	if ((_dim < 0) || (_dim > 2))
		throw "We only offer three dimensions until now";
	
	return downsampling[_dim];
}

// define resolution of model matrix along r
void reconSettings::set_rRes(const float _rRes)
{
	if (_rRes <= 0)
		throw "radial resolution must be bigger then 0";
	else
		rRes = _rRes;

	return;
}

// set resolution of output volume along certain dimension
void reconSettings::set_drRecon(const uint8_t _dim, const float dr)
{
	if (dr > 0)
	{
		drRecon[_dim] = dr;
	}
	else
	{
		printf("Resolution of reconstructed volume must be bigger 0\n");
		throw "invalidValue";
	}
	return;
}

// define resolution of reconstruction volume in a single go
void reconSettings::set_drRecon(const float dr0, const float dr1, const float dr2)
{
	set_drRecon(0, dr0);
	set_drRecon(1, dr1);
	set_drRecon(2, dr2);
	return;
}

float reconSettings::get_rRes() const {return rRes;}

float reconSettings::get_dr(const unsigned int _dim) const {return drRecon[_dim];}
float reconSettings::get_dr0() const {return drRecon[0];}
float reconSettings::get_dr1() const {return drRecon[1];}
float reconSettings::get_dr2() const {return drRecon[2];}
float* reconSettings::get_dr() {return &drRecon[0];}

uint64_t reconSettings::get_dim0() const {return get_dim(0);}
uint64_t reconSettings::get_dim1() const {return get_dim(1);}
uint64_t reconSettings::get_dim2() const {return get_dim(2);}

uint64_t reconSettings::get_dim(const uint8_t _dim) const
{
	const float delta = boundRecon[_dim * 2 + 1] - boundRecon[_dim * 2];
	// check if boundaries where defined correclty
	if (delta < 0)
	{
		printf("The upper boundary must be bigger then the lower along dfim %d\n", _dim);
		throw "invalidValue";
	}
	const unsigned int dim = round(delta / drRecon[_dim]) + 1;
	return dim;
}

void reconSettings::print_properties()
{
	printf("*** Recon settings ***\n");
	printf(" - SOS: %f m/s\n", SOS);
	printf(" - Threshold: %f\n", threshold);
	printf(" - nIter: %d\n", nIter);
	printf(" - zRange: %f ... %f\n", boundRecon[0], boundRecon[1]);
	printf(" - xRange: %f ... %f\n", boundRecon[2], boundRecon[3]);
	printf(" - yRange: %f ... %f\n", boundRecon[4], boundRecon[5]);
	printf(" - Resolution: %f x %f x %f\n", drRecon[0], drRecon[1], drRecon[2]);
	printf(" - Dim: %d x %d x %d\n", get_dim(0), get_dim(1), get_dim(2));
	printf(" - rRes: %f\n", rRes);
	return;
}