%% conf: 基本配置信息
function [config] = conf()
	config.controller = 'SteadyState';
	config.method = 'test';
	% config.solver.method = 'NR';
	% config.solver.maxIteration = 50;
	% config.solver.epsilon = 1e-3;
end