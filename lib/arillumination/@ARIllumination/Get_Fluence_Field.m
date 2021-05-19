% File: Get_Fluence_Field.m @ ARIllumination
% Author: Urs Hofmann
% Mail: hofmannu@biomed.ee.ethz.ch
% Date: 04.07.2019

% Output structure: isurface, ir, iz

function field = Get_Fluence_Field(ai, zSurfacesI, rVectorI, zVectorI)

	% function returns the field interpolated to the requested values
	% after interpolation from the monte carloed stuff
	ai.VPrintf('Interpolating fluence field... ', 1);
	tic 
	
	% Reduce z distance by fiber protusion
	% zVectorI = zVectorI - ai.protusion;

	% load fluence field if not already loaded
	ai.Load_Fluence_Field();

	% check if we are out of borders and if so, throw warning
	if (min(zSurfacesI) < min(ai.zSurface))
		txtMsg = ['Requested zSurface interp not in simulated range: ', ...
			num2str(min(ai.zSurface) * 1000, 3), ' mm min'];
		error(txtMsg);
	end

	if (max(zSurfacesI) > max(ai.zSurface))
		txtMsg = ['Requested zSurface interp not in simulated range: ', ...
			num2str(max(ai.zSurface) * 1000, 3), ' mm max'];
		error(txtMsg);
	end

	if (max(rVectorI) > max(ai.rSim))
		warning('Simulated rVector does not cover requested interp');
	end

	if (min(zVectorI) < min(ai.zSim))
		txtMsg = ['Requested zVector not in simulated range: ', ...
			num2str(min(ai.zSim) * 1000, 3), ' mm min'];
		error(txtMsg);
	end

	if (max(zVectorI) > max(ai.zSim))
		txtMsg = ['Requested zVector not in simulated range: ', ...
			num2str(max(ai.zSim) * 1000, 3), 'mm max'];
		error(txtMsg);
	end

	oldGrid{1} = ai.zSurface;
	oldGrid{2} = ai.rSim;
	oldGrid{3} = ai.zSim;

	newGrid{1} = zSurfacesI;
	newGrid{2} = rVectorI;
	newGrid{3} = zVectorI;
	
	F = griddedInterpolant(oldGrid, ai.fluenceField);

	% Performing interpolation
	field = F(newGrid); % structure: [iSurface, iR, iz] 

	tElapsed = toc;
	txtMsg = ['done after ', num2str(tElapsed, 2), ' sec!\n'];
	ai.VPrintf(txtMsg, 0);

	if any(isnan(field(:)))
		error('Fluence field contains NaNs');
	end

	if any(~isfinite(field(:)))
		error('Fluence field contains infs');
	end

	% TODO normalization should be replaced witch actual fluence calculated based on PD
	field = field / max(abs(field(:)));

	% plot interpolated field to double check if error occured
	if (ai.flagDisplay)
		figure('Name', 'Interpolated fluence field');
		
		plotLayer = round(size(field, 1) / 2);
		fieldLayer = squeeze(field(plotLayer, :, :));
		nCuts = 10;
		nR = size(field, 2);
		nZ = size(field, 3);
		dR = floor(nR / nCuts);
		dZ = floor(nZ / nCuts);

		ax1 = subplot(2, 2, [1, 3]);
		imagesc(rVectorI, zVectorI, fieldLayer');
		colormap(ax1, bone(1024));
		colorbar;
		xlabel('Radial distance [m]');
		ylabel('Axial distance');
		title('Cross section through fluence field');

		ax2 = subplot(2, 2, 2);

		plot(zVectorI, fieldLayer(1:dR:end, :));
		title('Axial cuts');
		xlabel('Axial distance [m]');
		ylabel('Light intensity');
		axis tight;
		grid on;

		ax3 = subplot(2, 2, 4);

		plot(rVectorI, fieldLayer(:, 1:dZ:end));
		title('Cross sections');
		xlabel('Radial distance [m]');
		ylabel('Light intensity');
		axis tight;
		grid on;

	end
end
