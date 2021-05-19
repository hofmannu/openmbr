% File: Compensate_1D.m @ ARIllumination
% Author: Urs Hofmann
% Mail: hofmannu@biomed.ee.ethz.ch
% Date: 25.06.2019

% Description: Function compensates fluence in 1D for each A scan which is considered to be a along the first dimension of inputVol

% Inputs:
%   - inputVol: 3D volume of measured data
%   - surface: 2D matrix describing the skin surface
%   - z: vector correlating z distance to index
% 	- lambda: wavelength

function inputVol = Compensate_1D(ari, inputVol, surface, z, lambda)

	% check first if dimensions are correct
	if (size(inputVol, 2) ~= size(surface, 1))
		error('2nd dim of inputVol and 1st of surface must be the same');
	end
	if (size(inputVol, 3) ~= size(surface, 2))
		error('3rd dim of inputVol and 2nd of surface must be the same');
	end
	if (length(z) ~= size(inputVol, 1))
		error('1st dim of inputVol and length of z must be the same');
	end

	[nZ, nX, nY] = size(inputVol); % get size of input vol
	nAScan = nX * nY; % number of a scans
	idxPlot = randi(nAScan); % which index should we plot
	inputVol = reshape(inputVol, [nZ, nAScan]); % reshape vol in a scan format
	surface = reshape(surface, [1, nAScan]); % reshape surface in a scan format
	ari.Load_Fluence_Field(); % load fluence map for fiber
	minSurface = min(surface); % min, max surface distance
	maxSurface = max(surface);
	dz = z(2) - z(1);
	zSurface = single(minSurface:dz:maxSurface);
	rVector = single(0); % we are only interested in central axis

	% get fluence field (returned as [surfaceLayer, r, z])
	field = ari.Get_Fluence_Field(rVector, z, zSurface);

	fprintf('[ARIllumination] Compensating fluence... ');
	for iAScan = 1:nAScan
		[~, minIdx] = min(abs(zSurface - surface(iAScan))); % find idx of surf
		% backup our old vector
		if (iAScan == idxPlot)
			backup = squeeze(inputVol(:, iAScan)); 
			backupFluence = squeeze(field(minIdx, 1, :));
		end
		inputVol(:, iAScan) = inputVol(:, iAScan) ./ squeeze(field(minIdx, 1, :));
			% perform actual compensation
		
		idxNoLight = (field(minIdx, 1, :) <= 0.005); % find no light areas
		inputVol(idxNoLight, iAScan) = 0; % set to 0
	end

	if ari.flagDisplay
		figure();
		hold on
		yyaxis right
		plot(z, backupFluence / max(abs(backupFluence)), '-');
		ylabel('Relative fluence');
		yyaxis left
		plot(z, backup, '--');
		plot(z, squeeze(inputVol(:, idxPlot)) / max(abs(inputVol(:, idxPlot))), '-');
		ylabel('Signal intensity');
		axis tight;
	end

	inputVol = reshape(inputVol, [nZ, nX, nY]); % reshape to original format
	fprintf('done!\n');

end
