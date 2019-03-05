%% getMpcSteady: TODO
function [mpc] = getMpcSteady(str)

	mpc = [];
	if nargin == 1

        if exist([str, '_ss'])
    		eval(['mpc = ', str, '_ss();']);
        elseif exist(str)
            eval(['mpc = ', str, '();']);
        end
	end

end