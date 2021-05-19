	// % build threedimensional gaussian beam
	// I = zeros(om.nX, om.nY, om.nZ, 'single');
	
	// % create radial distance field
	// rDist = sqrt(om.xVec.^2 + om.yVec'.^2);
	// rDist = repmat(rDist, [1, 1, om.nZ]);

	// % create axial distance field
	// zDist = om.zVec;
	// zDist = reshape(zDist, [1, 1, om.nZ]);
	// zDist = repmat(zDist, [om.nX, om.nY, 1]);

	// I0 = 1;

	// angSigma = asin(om.na / om.n);

	// w0 = om.wavelength / (pi * om.n * angSigma);

	// zr = pi * w0^2 * om.n / om.wavelength;

	// w = w0 * sqrt(1 + (zDist / zr).^2);
	// om.optModel = I0 * (w0 ./ w).^2 .* exp((-2 * rDist.^2) ./ w.^2);

	// % elongate optical model along time axis
	// optModelReshaped = reshape(om.optModel, [1, om.nX, om.nY, om.nZ]);
	// optModelReshaped = repmat(optModelReshaped, [om.nT, 1, 1, 1]);

	// % overlay acoustic and optical model
	// om.modelData = om.modelData .* optModelReshaped;

#include "gaussian_beam.h"

void gaussian_beam::convolve_model(float* modelMatrix, const uint64_t nt)
{
	// model matrix has order it, ir, iz
	if (flagVerbose)
		cout << "Convolving model with gaussian beam" << endl;

	float rPos, zPos;
	float fluence, w;
	uint64_t modelOffset;

	const float angSigma = asin(na / n);
	const float w0 = wavelength * 1e-9 / (M_PI * n * angSigma);
	const float zr = M_PI * pow(w0, 2) * n / (wavelength * 1e-9);

	for (uint64_t ir = 0; ir < nr; ir++)
	{
		rPos = ((float) ir) * dr;
		for (uint64_t iz = 0; iz < nz; iz++)
		{
			zPos = z0 + ((float) iz) * dz;

			w = w0 * sqrt(1 + pow(zPos / zr, 2));
			fluence = I0 * pow(w0 / w, 2) * exp((-2 * rPos * rPos) / pow(w, 2));
			
			// scale whole model with value
			modelOffset = nt * (ir + iz * nr);
			for (uint64_t it = 0; it < nt; it++)
			{
				modelMatrix[it + modelOffset] *= fluence; 
			}
		}
	}

	return;
}

void gaussian_beam::define_r(const float _dr, const uint64_t _nr)
{
	if (flagVerbose)
		cout << "Defining radial vector: " << _dr << ", " << _nr << endl;

	dr = _dr;
	nr = _nr;
	return;
}

void gaussian_beam::define_z(const float _dz, const float _z0, const uint64_t _nz)
{
	if (flagVerbose)
		cout << "Defining axial vector: " << _z0 << ", " << _dz << ", " << _nz << endl;

	nz = _nz;
	dz = _dz;
	z0 = _z0;
	return;
}

void gaussian_beam::set_wavelength(const float _wavelength)
{
	wavelength = _wavelength;
	return;
}

void gaussian_beam::set_n(const float _n)
{
	n = _n;
	return;
}

void gaussian_beam::set_na(const float _na)
{
	na = _na;
	return;
}

void gaussian_beam::set_I0(const float _I0)
{
	I0 = _I0;
	return;
}