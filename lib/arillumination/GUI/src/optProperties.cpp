#include "optProperties.h"

optProperties::optProperties()
{
	calc_albedo();
}

float optProperties::get_mua() const {return mua;}
float optProperties::get_mus() const {return mus;}
float optProperties::get_g() const {return g;}
float optProperties::get_n() const {return n;}
float optProperties::get_albedo() const {return albedo;}


void optProperties::set_mua(const float _mua)
{
	mua = _mua;
	calc_albedo();
	return;
}

void optProperties::set_mus(const float _mus)
{
	mus = _mus;
	calc_albedo();
	return;
}

void optProperties::set_g(const float _g)
{
	g = _g;
	return;
}

void optProperties::set_n(const float _n)
{
	n = _n;
	return;
}

void optProperties::calc_albedo()
{
	albedo = mus / (mus + mua);
	return;
}