%% case9_ft: IEEE-9 故障
function [mpc] = case9_ft()
	
	% nid, Ra, Xd'1, Xd'2, Xd'0, Xn, Connection
	mpc.gen = {
		1	0	0.0608	0.0608	0.18	0	'Yn';
		2	0	0.1198	0.1198	0.36	0	'Yn';
		3	0	0.1813	0.1813	0.54	0	'Yn';
	};

	% fid, tid, r, x1, x0, 
	mpc.line = [
		4	5	0.010	0.085	0.09;
		4	6	0.017	0.092	0.09;
		5	7	0.032	0.161	0.16;
		6	9	0.039	0.170	0.17;
		7	8	0.0085	0.072	0.07;
		8	9	0.0119	0.1008	0.10;
	];

	% fid, tid, rt, Xt1, Xt0, Xn, ratio, angle, Connection
	mpc.transformer = {
		1	4	0	0.0576	0.14	0.085*1.5 - 0.1205	1	0	'Dyn11';
		2	7	0	0.0625	0.16	0.085*1.5 - 0.1205	1	0	'Dyn11';
		3	9	0	0.0586	0.15	0.085*1.5 - 0.1205	1	0	'Dyn11';
	};

	mpc.fault = [];
	mpc.fault.nid = 5;
	mpc.fault.zf = 0;
	mpc.fault.zg = 0;
	mpc.fault.type = 'f1';
	mpc.fault.phase = 'b';

end