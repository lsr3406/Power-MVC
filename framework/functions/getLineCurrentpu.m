%% getLineCurrentpu: 计算线路电流
function [I] = getLineCurrentpu(S, U)
	I = conj(S./U);
end
