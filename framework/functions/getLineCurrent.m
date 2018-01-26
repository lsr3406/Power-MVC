%% getLineCurrent: 计算线路电流
function [I] = getLineCurrent(S, U)
	I = conj(S./sqrt(3)./U);
end
