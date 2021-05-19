#ifndef OPTPROPERTIES_H
#define OPTPROPERTIES_H

class optProperties
{
private:
	float mua = 0.3 * 1000; // absorption coefficient [1/m]
	float mus  = 4.2 * 1000; // scattering coefficient [1/m]
	float g = 0.8; // anisotropy of medium
	float n = 1.41; // optical index

	// dependent properties (no set function)
	float albedo; 

	void calc_albedo();
public:
	optProperties();

	float get_mua() const;
	float get_mus() const;
	float get_g() const;
	float get_n() const;
	float get_albedo() const;

	float* get_pg() {return &g;};
	float* get_pn() {return &n;};

	void set_mua(const float _mua);
	void set_mus(const float _mus);
	void set_g(const float _g);
	void set_n(const float _n);
};

#endif