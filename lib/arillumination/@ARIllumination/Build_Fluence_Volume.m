% File: Build_FLuence_VOlume.m @ ARIllumination
% Author: Urs Hofmann
% Mail: hofmannu@biomed.ee.ethz.ch
% Date: 07-Apr-2020

% Description: Uses the 

function Build_FLuence_Volume(ar, varargin)

	dx = 20e-6; % resolution of grid in x direction
	dy = 20e-6; % resolution of grid in y direction
	dz = 10e-6; % resolution of grid in z direction
	surface = zeros(500, 500, 'single');
	nz = 500;
	z0 = 3e-3;

	for iargin=1:2:(nargin - 1)
		switch varargin{iargin}
			
			otherwise
				error('Invalid input argument');
		end
	end

	nx = size(surface, 1);
	ny = size(surface, 2);

	fluenceVol = zeros(nz, nx, ny, 'single');

end
