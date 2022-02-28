#include "reconsett.h"

void reconsett::set_nIter(const int16_t _nIter)
{
	if (_nIter < 1)
	{
		printf("Number of iterations needs to be at least 1");
		throw "InvalidValue";
	}
	nIter = _nIter;
	return;
}

void reconsett::set_regStrength(const float _regStrength)
{
	if (_regStrength < 0)
	{
		printf("Regularization strength needs to be at least 0");
		throw "InvalidValue";
	}
	return;
}

void reconsett::set_dr(const vector3<float> _dr)
{
	if ((dr.x < 0.0) || (dr.y < 0.0) || (dr.z < 0.0))
	{
		printf("Resolution must be bigger then one along all dimensions");
		throw "InvalidValue";
	}
	dr = _dr;
	return;
}

void reconsett::set_startPos(const vector3<float> _startPos)
{
	startPos = _startPos;
	return;
}

void reconsett::set_endPos(const vector3<float> _endPos)
{
	endPos = _endPos;
	return;
}

void reconsett::set_tMin(const float tMin)
{
	tCrop[0] = tMin;
	return;
}
void reconsett::set_tMax(const float tMax)
{
	tCrop[1] = tMax;
	return;
}
void reconsett::set_xMin(const float xMin)
{
	startPos.x = xMin;
	return;
}
void reconsett::set_xMax(const float xMax)
{
	endPos.x = xMax;
	return;
}

void reconsett::set_yMin(const float yMin)
{
	startPos.y = yMin;
	return;
}

void reconsett::set_yMax(const float yMax)
{
	endPos.y = yMax;
	return;
}

void reconsett::set_cropT(const float tMin, const float tMax)
{
	tCrop[0] = tMin;
	tCrop[1] = tMax;
	return;
}

void reconsett::set_cropX(const float xMin, const float xMax)
{
	startPos.x = xMin;
	endPos.x = xMax;
	return;
}

void reconsett::set_cropY(const float yMin, const float yMax)
{
	startPos.y = yMin;
	endPos.y = yMax;
	return;
}

