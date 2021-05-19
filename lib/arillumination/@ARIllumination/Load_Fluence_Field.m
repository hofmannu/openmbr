% File: Load_Fluence_Field.m @ ARIllumination
% Author: Urs Hofmann
% Mail: hofmannu@biomed.ee.ethz.ch
% Date: 04.07.2019

% Description: Loads the monte carlo based fluence field into class

function Load_Fluence_Field(ai)

	if isfile(ai.fieldPath)
		load(ai.fieldPath); % load fluence field file
		ai.fluenceField = fluenceField; % assign field matrix
		ai.drSim = dr; % assign resolution of simulation
		ai.zSurface = zSurface; % assing surface vector

		% check data integrity
		if (dCore ~= ai.dCore)
			error('Core diameters do not match.');
		end

		if (wavelength ~= ai.wavelength)
			error('Wavelengths do not match');
		end

		if (length(zSurface) ~= size(fluenceField, 1))
			error('Length of surface vector and fluence field have mismatch');
		end
	else
		errMsg = [ai.fieldPath, ' is not pointing to a valid file'];
		error(errMsg);
	end

end
