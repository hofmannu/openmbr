% File: Build_Fluence_Field.m @ ARIllumination
% Author: Urs Hofmann
% Mail: hofmannu@biomed.ee.ethz.ch
% Date: 04.07.2019

% Input:
%  - rMax: maximum radius calculated [m]
%  - zSurface: surface layer steps calculated [m]

function Build_Fluence_Field(ai, varargin)

	% default settings
	rMax = 10e-3;
	zSurface = 4e-3:0.1e-3:15e-3;
	nPhotons = 5e8;
	flagDisplay = 1;

	for (iargin = 1:2:(nargin-1))
		switch varargin{iargin}
			case 'rMax'
				rMax = varargin{iargin + 1};
			case 'zSurface'
				zSurface = varargin{iargin + 1};
			case 'wavelength'
				ai.wavelength = varargin{iargin + 1};
			case 'nPhotons'
				nPhotons = varargin{iargin + 1};
			case 'flagDisplay'
				flagDisplay = varargin{iargin + 1};
			otherwise
				error('Invalid argument passed to function');
		end
	end

	% zSurface: vector containing surface steps to simulate
	nLayers = length(zSurface);
	nZ = int16((zSurface(end) - zSurface(1)) / ai.drSim) + 1;
	nR = floor(rMax / ai.drSim) + 1;
	rVec = (0:(nR - 1)) * ai.drSim;
 	zVec = zSurface(1) + single(0:(nZ - 1)) * ai.drSim;	

	% get scattering and absorption coefficient from lookup table
	[mua, mus] = getAbsorptionCoefficient(ai.wavelength);

	% tissue
	%   absorption coefficient in 1/m
	%   scattering coefficient in 1/m
	%   anisotropy factor
	% 	optical refractive index
	tissue = single([mua, mus, 0.8, 1.41]);

	% fiber
	% 	core radius [m]
	% 	numerical aperture of fiber
	fiber = single([ai.rCore, ai.NA]);

	% field
	% 	maximum simulated radial distance [m]
	%   maximum simulated axial distance [m]
	% 	minimum simulated axial distance (equal to surface layer) [m]
	%   radial resolution [m]
	%   axial resolution [m]
	field = single([rMax, zSurface(end), zSurface(1), ai.drSim, ai.drSim]);

	% sim
	% 	number of photons in each thread
	%   overall number of simulated photons
	sim = single([100, nPhotons]);

	ai.fluenceField = zeros(nLayers, nR, nZ, 'single');
	
	if (flagDisplay)
		figure();
		% plot fluence field
		ax1 = subplot(2, 2, [1, 3]);
		fieldViz = imagesc(rVec, zVec, squeeze(ai.fluenceField(1, :, :))');
		xlabel('Radial distance [m]');
		ylabel('Axial distance [m]');
		colormap(ax1, bone(1024));
		colorbar;
		title('Fluence field');

		% plot cut through central line
		ax2 = subplot(2, 2, 2);
		hold on;
		for i = 1:50:size(ai.fluenceField, 2)
			axialFluence{i} = plot(zVec, squeeze(ai.fluenceField(1, i, :)));
		end
		hold off;
		xlabel('Axial distance [m]');
		ylabel('Light intensity');
		axis tight;
		grid on;

		% plot cuts in radial direction at different depths
		ax4 = subplot(2, 2, 4);
		hold on;
		for i = 1:50:size(ai.fluenceField, 3)
			radialFluence{i} = plot(rVec, squeeze(ai.fluenceField(1, :, i)));
		end
		hold off;
		xlabel('Radial distance [m]');
		ylabel('Light intesnity');
		axis tight;
		grid on;
	end

	for iLayer = 1:nLayers
		fprintf('--> Simulating layer %d of %d <--\n', iLayer, nLayers);
		field(3) = single(zSurface(iLayer)); % redefine surface starting layer
		tempField = m_mc(tissue, fiber, field, sim);
		tempField(1, 1) = tempField(2, 1);
		idx = size(tempField, 2);
		idxStart = nZ - idx + 1;
		ai.fluenceField(iLayer, :, idxStart:end) = tempField;
	
		% fill empty void above surface
		rSurface = ai.rCore + tan(ai.theta) * zVec(idxStart);
		[~, rMaxIdx] = min(abs(rVec - rSurface));
		fluence = sum(ai.fluenceField(iLayer, 1:rMaxIdx, idxStart)) / rMaxIdx;
		for (i = 1:(idxStart-1))
			rMax = ai.rCore + tan(ai.theta) * zVec(i);
			[~, rMaxIdx] = min(abs(rVec - rMax));
			ratio = rSurface^2 / rMax^2;
			ai.fluenceField(iLayer, 1:rMaxIdx, i) = ratio * fluence;
		end

		% update plot to latest simulation result
		if (flagDisplay)
			set(fieldViz, 'cdata', log(squeeze(ai.fluenceField(iLayer, :, :))'));
			for i=1:50:size(ai.fluenceField, 2)
				set(axialFluence{i}, 'ydata', squeeze(ai.fluenceField(iLayer, i, :)));
			end
			for i=1:50:size(ai.fluenceField, 3)
				set(radialFluence{i}, 'ydata', squeeze(ai.fluenceField(iLayer, :, i)));
			end
			drawnow();
		end
	end

	% save calculated fluence field and all important properties to file
	fluenceField = ai.fluenceField;
	save(ai.fieldPath, 'fluenceField', '-nocompression', '-v7.3');
	save(ai.fieldPath, 'rMax', '-append');
	save(ai.fieldPath, 'zSurface', '-append');

	dr = ai.drSim;
	save(ai.fieldPath, 'dr', '-append');

	NA = ai.NA;
	save(ai.fieldPath, 'NA', '-append');

	dCore = ai.dCore;
	save(ai.fieldPath, 'dCore', '-append');

	wavelength = ai.wavelength;
	save(ai.fieldPath, 'wavelength', '-append');

end
