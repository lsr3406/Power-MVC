%% getParallel: 返回阻抗的并联值
function [z] = getParallel(arr)
	z = 1./sum(1./arr);
end