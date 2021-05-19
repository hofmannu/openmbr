% File: ARIllumination.m @ ARIllumination
% Author: Urs Hofmann
% Mail: hofmannu@biomed.ee.ethz.ch
% Date: 27.04.2019

% Description: Calculates illumination pattern for acoustic resolution opto-
% acoustic microscopy.

% all dimensions in [m], [rad]

classdef ARIllumination < handle

  properties
 		name = 'FT200UMT';
		NA(1, 1) double {mustBeNumeric, mustBePositive} = 0.39; % numerical aperture of fiber
    dCore(1, 1) double {mustBeNumeric, mustBeNonnegative} = 200e-6; % core diameter of fiber
    wavelength(1, 1) double {mustBeNumeric, mustBeNonnegative} = 532; % wavelength of the used laser
    flagDisplay(1, 1) logical = 1; % determines if we should show any figures
		protusion(1, 1) double {mustBeNumeric} = 0e-3; % protusion of fiber through hole [m]
		flagVerbose(1, 1) logical = 1; % enables or suppresses verbose output
	end

  properties (Dependent = true)
    filePath(1, :) char; % path where fiber properties are stored e.g. NA, dCore ...
		fieldPath(1, :) char; % path where 3D calculated field (MC result) is stored
		theta(1, 1) double {mustBePositive}; % opening angle at output
    rCore(1, 1) double {mustBeNumeric, mustBeNonnegative}; % core radius of MM fiber
    dAxial(1, 1) double {mustBeNumeric, mustBeNonnegative}; % resolution in axial direction
    rSurface(1, 1) double {mustBeNonnegative}; % spot radius at skin surface
		rSim(1, :) double; % r vector of simulation
		zSim(1, :) double; % z vector of simulation
	end

	properties (SetAccess = private)
		% fluence field calculated using MC simulations
		fluenceField(:, :, :) double = [];
		% surface distance vector for individual fluence field steps
		zSurface(1, :) double;
	end

  properties (Hidden, SetAccess = private)
    nMedium = 1.33; % refractive index in coupling medium
  	drSim = 15e-6; % spacing used for simulations
	end

  methods

		function set.wavelength(ai, wavelength)
			if(ai.wavelength ~= wavelength)
				ai.fluenceField = []; 
				% if the new wavelength is not the same as the previous remove field
				ai.zSurface = [];
			end
			% Write wavelength to member variable
			ai.wavelength = wavelength;
		end

		% get radial vector of simulated field
		function rSim = get.rSim(ai)
			nR = size(ai.fluenceField, 2);
			rSim = 0 : ai.drSim: (nR-1) * ai.drSim;
		end

		% get axial vector of simulated field
		function zSim = get.zSim(ai)
			nZ = size(ai.fluenceField, 3);
			zSim = 0:ai.drSim:(nZ - 1) * ai.drSim;
			zSim = zSim + ai.zSurface(1);
		end

    % get angle at output
    function theta = get.theta(ai)
      theta = asin(ai.NA / ai.nMedium);
    end

		% get radius of core which depends on core diameter
    function rCore = get.rCore(ai)
      rCore = ai.dCore / 2;
    end

    % beam radius at skin surface
    function rSurface = get.rSurface(ai)
      rSurface = ai.rCore + tan(ai.theta) * (ai.zSurface - ai.zFiber);
    end

		% returns path where fiber properties are stored
		function filePath = get.filePath(ai)
			filePath = [get_path('arfiber'), ai.name, '.mat'];
		end	

		% returns path where field is stored
		function fieldPath = get.fieldPath(ai)
			fieldPath = [get_path('arfiber'), ...
				ai.name, '_', num2str(ai.wavelength), '_field.mat'];
		end
	
		function set.name(ai, name)
			% if name is different from currently implemented fiber
			% , try to read it from file
			if ~strcmp(ai.name, name)
				ai.name = name;
				ai.Read_From_File();
			end
		end

		% reads and writes fiber properties from or to file based on name
		Save_To_File(ai);
		Read_From_File(ai);
		Build_Fluence_Field(ai, varargin); % build and save fluence field
		Load_Fluence_Field(ai); % load fluence field from file
		field = Get_Fluence_Field(ai, rVectorI, zVectorI, zSurfacesI); % interpolate fluence field
		inputVol = Compensate_1D(ari, inputVol, surface, z, lambda);
		VPrintf(ai, txtMsg, flagName);
	end


end
