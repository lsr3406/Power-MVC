%% getMpcSteady: TODO
function [mpc] = getMpcSteady(str)

	mpc = [];
	if nargin == 1
		eval(['mpc = ', str, '_ss();']);
		return ;
	end

end