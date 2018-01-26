%% getPowerLoss: 计算线路有功和无功损耗
function [dS] = getPowerLoss(S, U, z)
	dS = S.*conj(S)./U.^2.*z;
end
