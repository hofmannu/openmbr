% File: test_cuda.m
% Author: Urs Hofmann
% Mail: hofmannu@biomed.ee.ethz.ch
% Date: 24-02-2020

% Description: 
% Performs the unit test for our monte carlo based fluence field estiamtion

% mexcuda m_mc.cu % compile mc_cuda

tissue = single([0.3*1000, 4.2*1000, 0.8, 1.41]);
% [0] absorption coefficient [1/m]
% [1] scattering coefficient [1/m]
% [2] anisotroy factor
% [3] optical index of tissue

fiber = single([550e-6 / 2, 0.22]);
% [0] core radius of multimode fiber
% [1] numerical aperture

field = single([8e-3, 15e-3, 5e-3, 10e-6, 10e-6]);
% [0] maximum simulated radial distance [m]
% [1] distance between lower field end and fiber output [m]
% [2] distance between fiber output and skin surface [m]
% [3] spatial resolution of output field in radial direction [m]
% [4] spatial resolution of output field in axial direction [m]

sim = single([500, 1.25e9]);
% [0] number of photons per thread
% [1] overall number of photons to simulate

dSurf = 4e-3:0.1e-3:8e-3;
dSurf = single(dSurf);

for iSurf = 1:length(dSurf)
	field(3) = dSurf(iSurf);
	% run actual cuda code (there must be a file called m_mc.mexa64)
	tic
	output = m_mc(tissue, fiber, field, sim);
	toc
	% FIXME there is still a bug with the uppermost voxel which seems to disappear
	% in the moment we change the binning approach
	output(1, 1) = output(2, 1);
	
	% display result
	flagDisplay = 1;
	if flagDisplay
		figure();

		nR = size(output, 1); % number of elements in radial direction
		nZ = size(output, 2); % number of elements in axial direction

		r{iSurf} = ((0:(nR - 1)) + 0.5)* field(4); % radial vector
		z{iSurf} = ((0:(nZ - 1)) + 0.5)* field(5) + dSurf(iSurf); % axial vector
		% NOTE last element of vector in both radial and axial field contains escaped
		% photons

		% plot fluence field
		ax1 = subplot(2, 2, [1, 3]);
		imagesc(r{iSurf}, z{iSurf}, log(output'));
		colormap(bone(1024));
		c = colorbar;
		c.Label.String = 'Absorption';
		axis image
		xlabel('Radial distance of bins [m]');
		ylabel('Axial distance of bins [m]');

		ax2 = subplot(2, 2, 2);
		plot(z{iSurf}(1:end-1), output(1, (1:end-1)));
		xlabel('zDepth')
		ylabel('Photon count')
		title('Central depth fluence');
		axis tight
		grid on

		ax3 = subplot(2, 2, 4);
		plot(r{iSurf}(1:end-1), output(1:end-1, 1:50:end-100))
		xlabel('Lateral distance')
		ylabel('Photon count')
		title('Cross sections');
		axis tight
		grid on
		drawnow();

	end
	
	fluenceField{iSurf} = output;
end


flagSave = 1;
if flagSave
	save('/home/hofmannu/fluenceForPavel.mat', ...
		'tissue', 'fiber', 'field', 'sim', 'fluenceField', 'dSurf', 'r', 'z');
end
