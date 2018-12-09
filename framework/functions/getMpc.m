%% getMpc: 
function [mpc] = getMpc(str)

	mpc = [];
	if nargin == 1
		eval(['mpc = ', str, '();']);
		return ;
	end

end