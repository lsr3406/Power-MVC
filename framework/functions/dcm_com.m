% encoding: utf-8

%% dcm_com: 直流输电原始数据, 替代 14 节点系统 5-4 交流线路
function [mpc] = dcm_com()
	mpc.data = [
		% source sink Tr Ti Rcr Rci Rl mode sourceParam sinkParam
		5	4	0.96	0.94	0.028	0.069	0.334e-2	1	58.5	10/180*pi;
	];
end
