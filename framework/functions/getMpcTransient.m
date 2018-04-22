%% getMpcTransient: TODO getMpcTransient description
function [mpc] = getMpcTransient(str)

	mpc = [];
	if nargin == 1
		eval(['mpc = ', str, '_ts();']);
		return ;
	end	

end