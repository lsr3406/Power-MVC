%% lgq_eg8_1_ts: lgq_eg8_1 暂态稳定
function [mpc] = lgq_eg8_1_ts()

	mpc.gen = [
		1	8.47	0	2.241	0.304	2.241	2.241	0	0	0;
		3	1e7		0	0.0001	0.0001	0.0001	0.0001	0	0	0;
	];

	% 电网操作
	mpc.operating = [];

	mpc.operating(1).time = 0;
	mpc.operating(1).nid = 2;
	mpc.operating(1).ntype = 'f3';
	mpc.operating(1).zf = 1e-7;
	mpc.operating(1).zg = 1e-7;

	mpc.operating(2).time = 0.056;
	mpc.operating(2).fid = 2;
	mpc.operating(2).tid = 3;
	mpc.operating(2).btype = 'break';
	mpc.operating(2).nid = 2;
	mpc.operating(2).ntype = 'restore';
end