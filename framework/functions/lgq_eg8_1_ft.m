%% lgq_eg8_1_ft: lgq_eg8_1 故障
function [mpc] = lgq_eg8_1_ft()

	% nid, Ra, Xd'1, Xd'2, Xd'0, Xn, Connection
	mpc.gen = {
		1	0	0.304	0.043	0.304	0	'Y';
	};

	% fid, tid, r, x1, x0, 
	mpc.line = [
		2	3	0	0.47	1.88;
		2	3	0	0.47	1.88;
	];

	% fid, tid, rt, Xt1, Xt0, Xn, ratio, angle, Connection
	mpc.transformer = {
		1	2	0	0.13	0.13	0	1	0	'Dyn11';
		4	3	0	0.108	0.108	0	1	0	'Dyn11';
	};

	mpc.fault = [];
	mpc.fault.nid = 2;
	mpc.fault.zf = 0;
	mpc.fault.zg = 0;
	mpc.fault.type = 'f2';
	mpc.fault.phase = 'ab';
end
