% File: VPrintf.m @ ARIllumination
% Author: Urs Hofmann
% Mail: hofmannu@biomed.ee.ethz.ch
% Dtae: 15.12.2019

% Description: Selective class output

function VPrintf(ai, txtMsg, flagName)

	if ai.flagVerbose
		if flagName
			txtMsg = ['[ARIlluminationm] ', txtMsg];
		end
		fprintf(txtMsg);
	end

end
