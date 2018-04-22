%% getMpcFault: TODO
function [mpc] = getMpcFault(str)

	mpc = [];
	if nargin == 1
		eval(['mpc = ', str, '_ft();']);
		return ;
	end

end