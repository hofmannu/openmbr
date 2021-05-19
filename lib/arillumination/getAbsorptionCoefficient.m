% File: getAbsorptionCoefficient.m
% Author: Urs Hofmann
% Mail: hofmannu@biomed.ee.ethz.ch
% Date: 25.06.2019

% Description: Returns both absorption and scattering coefficient in skin.

function [mua, mus] = getAbsorptionCoefficient(lambda)

	if ((lambda >= 450) && (lambda <= 600))
		% data taken from Bosschaart2011
		lambdaV = 450:10:600;
		muaV = [0.72, 0.66, 0.58, 0.54, 0.43, 0.34, 0.29, 0.27, 0.33, 0.37, 0.34, 0.31, 0.33, 0.29, 0.18, 0.08];
		mua = interp1(lambdaV, muaV, double(lambda), 'linear') * 1e3;
		musV = [2.08, 2.03, 2.00, 1.96, 1.93, 1.90, 1.87, 1.83, 1.80, 1.77, 1.74, 1.71, 1.68, 1.65, 1.62, 1.60];
		mus = interp1(lambdaV, musV, double(lambda), 'linear') * 2.5e3;
	elseif (lambda == 1064)
		mua = 0.4e2;	
		% Bashkatov 2005
		mus = 1.1e12 * lambda^(-4) + 73.7 * lambda^(-0.22); % 1 / cm
		mus = mus * 100; % now in 1 / m
	else
		error('Did not specifiy wavelength yet');
	end



end
