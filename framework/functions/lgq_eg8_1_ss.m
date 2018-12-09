%% lgq_eg8_1_ss: lgq_eg8_1 稳态
function [mpc] = lgq_eg8_1_ss()
	
	mpc.baseMVA = 100;
	% Nodes(id,type,Pd,Qd,conductance,susceptance,area,mag0,ang0,baseKV,zone,Vmax,Vmin)
	mpc.bus = [

		1	3	0	0	0	0	1	1.00	0	16.5	1	1.1	0.9;
		2	1	0	0	0	0	1	1.00	0	18.0	1	1.1	0.9;
		3	1	0	0	0	0	1	1.00	0	13.8	1	1.1	0.9;
		4	1	110	22	0	0	1	1.00	0	230.0	1	1.1	0.9;
	];

	% Generator(nodeId,Pg,Qg,Qmax,Qmin,Vg,mBase,status,Pmax,Pmin)
	mpc.gen = [

		3	0	0	300	-300	1.00	100	1	250	-10;
	];

	% Branches(from, to, r, x, b, rateA, rateB, rateC, ratio, angle, status, angmin, angmax)
	mpc.branch = [

		1	2	0	0.05	0	0	0	0	1	0	1	-180	180;
		2	3	0	0.12	0	0	0	0	0	0	1	-180	180;
		% 2	3	0	0.12	0	0	0	0	0	0	1	-180	180;
		3	4	0.005	0.05	0	0	0	0	1	0	1	-180	180;
	];
end
