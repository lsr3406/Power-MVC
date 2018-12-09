%% case9_ss: IEEE-9 暂态
function [mpc] = case9_ts()

	% 所有负荷目前按照纯阻抗处理, 不考虑有发电机, 电力电子设备等
	% 未直接考虑各阻尼绕组的影响
	% nodeId, TJ, Ra, Xd, X'd, Xq, X'q, T'd0, T'q0, D
	mpc.gen = [
		1	47.28	0	0.146	0.0608	0.0608	0.0608	8.96	0.31	0;
		2	12.80	0	0.8958	0.1198	0.1198	0.1198	6.00	0.535	0;
		3	6.02	0	1.3125	0.1813	0.1813	0.1813	8.59	0.600	0;
		% 1	47.28	0	0.146	0.0608	0.0969	0.0969	8.96	0.31	0;
		% 2	12.80	0	0.8958	0.1198	0.8645	0.1969	6.00	0.535	0;
		% 3	6.02	0	1.3125	0.1813	1.2578	0.2500	8.59	0.600	0;
	];

	% 电网操作
	mpc.operating = [];

	% 节点 7 发生三相短路 (线路 5-7 末端两相短路接地)
	mpc.operating(1).time = 0;
	mpc.operating(1).nid = 7;
	mpc.operating(1).ntype = 'f3';
	mpc.operating(1).zf = 1e-8;
	mpc.operating(1).zg = 1e-8;

	% 线路 5-7 退出运行
	mpc.operating(2).time = 0.0833;
	mpc.operating(2).fid = 5;
	mpc.operating(2).tid = 7;
	mpc.operating(2).btype = 'break';
	mpc.operating(2).nid = 7;
	mpc.operating(2).ntype = 'restore';
end