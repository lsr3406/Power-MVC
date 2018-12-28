% encoding: utf-8
% @author: siru
% @create on: 2018-09-29 14:43:39
% @update on: 2018-12-10 11:02:37
%% dcm_test: 直流输电原始数据, 测试版
function [mpc] = dcm_test()
	mpc.data = [
		% source sink Tr Ti Rcr Rci Rl mode sourceParam sinkParam
		2	3	1.00	1.00	0.2	0.2	0.5	1	30	15/180*pi;
	];
end
