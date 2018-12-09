function mpc = case5_test
%CASE5  Power flow data for modified 5 bus, 5 gen case based on PJM 5-bus system
%   Please see CASEFORMAT for details on the case file format.
%
%   Based on data from ...
%     F.Li and R.Bo, "Small Test Systems for Power System Economic Studies",
%     Proceedings of the 2010 IEEE Power & Energy Society General Meeting

%   Created by Rui Bo in 2006, modified in 2010, 2014.
%   Distributed with permission.

%   MATPOWER

%% MATPOWER Case Format : Version 2
mpc.version = '2';

%%-----  Power Flow Data  -----%%
%% system MVA base
mpc.baseMVA = 100;

%% bus data
%	bus_i	type	Pd	Qd	Gs	Bs	area	Vm	Va	baseKV	zone	Vmax	Vmin
mpc.bus = [
	1	1	160	80	0	0	1	1.00	0	230	1	1.1	0.9;
	2	1	200	100	0	0	1	1.00	0	230	1	1.1	0.9;
	3	1	370	130	0	0	1	1.00	0	230	1	1.1	0.9;
	4	2	0	0	0	0	1	1.05	0	230	1	1.1	0.9;
	5	3	0	0	0	0	1	1.05	0	230	1	1.1	0.9;
];

%% generator data
%	bus	Pg	Qg	Qmax	Qmin	Vg	mBase	status	Pmax	Pmin	Pc1	Pc2	Qc1min	Qc1max	Qc2min	Qc2max	ramp_agc	ramp_10	ramp_30	ramp_q	apf
mpc.gen = [
	4	500	0	300	-300	1.05	100	1	800	100	0	0	0	0	0	0	0	0	0	0	0;
	5	0	0	500	-210	1.05	100	1	800	100	0	0	0	0	0	0	0	0	0	0	0;
];

%% branch data
%	fbus	tbus	r	x	b	rateA	rateB	rateC	ratio	angle	status	angmin	angmax
mpc.branch = [
	1	2	0.04	0.25	0.50	200	0	0	0.00	0	1	-360	360;
	1	3	0.10	0.35	0.00	65	0	0	0.00	0	1	-360	360;
	2	3	0.08	0.30	0.50	200	0	0	0.00	0	1	-360	360;
	2	4	0.00	0.015	0.00	600	0	0	1.05	0	1	-360	360;
	3	5	0.00	0.03	0.00	500	0	0	1.05	0	1	-360	360;
];

%%-----  OPF Data  -----%%
%% generator cost data
%	1	startup	shutdown	n	x1	y1	...	xn	yn
%	2	startup	shutdown	n	c(n-1)	...	c0
mpc.gencost = [
	2	0	0	3	50.4395.*1e-4	200.4335.*1e-2	1200.5485;
	2	0	0	3	200.55.*1e-4	500.7460.*1e-2	1857.2010;
];

