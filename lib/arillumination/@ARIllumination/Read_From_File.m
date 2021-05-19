% File: Read_From_File.m @ ARIllumination
% Author: Urs Hofmann
%	Mail: hofmannu@biomed.ee.ethz.ch
% Date: 25.06.2019

% Description: Reads fiber properties from file.

function Read_From_File(ai)
	
	if isfile(ai.filePath)			
		load(ai.filePath);

		ai.NA = NA; % numerical aperture of fiber
		ai.dCore = dCore; % core diameter of fiber
	else
		warning('Could not find file!');
	end
end
