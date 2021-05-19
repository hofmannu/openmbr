% File: Save_To_File.m @ ARIllumination
% Author: Urs Hofmann
% Mail: hofmannu@biomed.ee.ethz.ch
% Date: 25.06.2019

% Description: Save fiber properties to file

function Save_To_File(ai)
	NA = ai.NA;
	save(ai.filePath, 'NA');

	dCore = ai.dCore;
	save(ai.filePath, 'dCore', '-append');
end
