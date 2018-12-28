%% SteadyState
classdef SteadyState < handle
	properties
		mBase;

		%% 下面这 4 个表示电网中一些元素的状态, 它们可能经常改变.
		bus;		% 节点状态, [id, type, mag, ang, Pg, Qg, Pc, Qc, Pd, Qd, Pis, Qis, Pout, Qout, dP, dQ, Pmax, Qmax, Pmin, Qmin, Vmax, Vmin]
		gen;	% 发电机状态, [id, nid, Pg, Qg, status]
		branch;	% 线路状态,	[id, fid, tid, ratio, Pij, Qij, Pji, Qji, dP, dQ, status]
		gencost;

		NAM;		% 节点导纳矩阵

		FDM1;		% 快速分解法修正方程第一系数矩阵
		FDM2;		% 快速分解法修正方程第二系数矩阵
		%% 雅可比矩阵在求解过程中一直改变, 不直接当做属性存储

		itlog;		% 和电力系统潮流求解相关的记录

		% 普通潮流业务
		pf;
		he;

		shuntFlag;

		% 最优潮流业务
		opf;
		% opf.n_eq: 等式约束的数量 (节点数 * 2)
		% opf.n_state: 自变量的数量 (节点数 + 2 * 发电机数)
		% opf.n_ieq: 不等式约束的数量 (节点数 + 2 * 发电机数 + 2 * 线路数)
		% opf.x: 自变量对应的列向量 (n_state * 1)
		% 	[发电机有功, 发电机无功, 节点相角, 节点电压]
		% opf.y: 等式约束对应的拉格朗日乘子列向量 (n_eq * 1)
		% 	[潮流方程有功乘子, 潮流方程无功乘子]
		% opf.z: 不等式(下限)约束对应的拉格朗日乘子列向量 (n_ieq * 1)
		% 	[发电机有功约束乘子, 发电机无功约束乘子, 节点电压约束乘子, 线路始端有功约束乘子, 线路末端有功约束乘子]
		% opf.w: 不等式(上限)约束对应的拉格朗日乘子列向量 (n_ieq * 1)
		% 	[发电机有功约束乘子, 发电机无功约束乘子, 节点电压约束乘子, 线路始端有功约束乘子, 线路末端有功约束乘子]
		% opf.l: 不等式(下限)约束对应的松弛变量列向量 (n_ieq * 1)
		% 	[发电机有功约束松弛变量, 发电机无功约束松弛变量, 节点电压约束松弛变量, 线路有功约束松弛变量, 线路有功约束松弛变量]
		% opf.u: 不等式(下限)约束对应的松弛变量列向量 (n_ieq * 1)
		% 	[发电机有功约束松弛变量, 发电机无功约束松弛变量, 节点电压约束松弛变量, 线路有功约束松弛变量, 线路有功约束松弛变量]
		% opf.g_min: 不等式约束对应的下限 (n_ieq * 1)
		% 	[发电机有功下限, 发电机无功下限, 节点电压下限, 线路有功上限相反数, 线路有功上限相反数]
		% opf.g_max: 不等式约束对应的上限 (n_ieq * 1)
		% 	[发电机有功上限, 发电机无功上限, 节点电压上限, 线路有功上限, 线路有功上限]
		% opf.epsilon: 内点法求解精度, 默认 1e-6
		% opf.n_iters_max: 最大迭代次数, 默认 50
		% opf.sigma: 向心参数, 默认 0.01
		% opf.n_iters: 迭代次数
		% opf.gap: 对偶间隙
		% opf.mu: 障碍参数
		% opf.L_inv_Z: 临时变量, z/l 的元素组成的对角阵 (n_ieq * n_ieq)
		% opf.U_inv_W: 临时变量, w/u 的元素组成的对角阵 (n_ieq * n_ieq)
		% opf.Yed: 临时变量, 论文中有详细说明 (节点数 * 节点数)
		% opf.Yed_Re: opf.Yed 的实部 (节点数 * 节点数)
		% opf.Yed_Im: opf.Yed 的虚部 (节点数 * 节点数)
		% opf.f: 总费用
		% opf.nable_F: 目标函数的梯度 (n_state * 1)
		% opf.nable2_F: 目标函数的海森矩阵 (n_state * n_state)
		% opf.nable_G: 不等式约束方程组的雅可比矩阵 (n_ieq * n_state)
		% opf.nable2_G_c: 不等式约束方程组的海森矩阵与相应乘子相乘的结果 (n_state * n_state)
		% opf.nable_H: 等式约束对应的雅可比矩阵 (n_eq * n_state)
		% opf.nable2_H_y: 等式约束对应的海森矩阵与相应乘子相乘的结果 (n_state * n_state)
		% opf.valueMatrix: 修正方程等号右侧的列向量 ( (n_state + n_eq + 4 * n_ieq) * 1 )
		% opf.h: 等式约束的函数值 (n_eq * 1)
		% 	[节点有功功率误差, 节点无功功率误差]
		% opf.Ly: 同上, 等式约束的函数值 (n_eq * 1)
		% opf.g: 不等式约束对应的实际变量 (n_ieq * 1)
		% 	[发电机有功, 发电机无功, 节点电压, 线路有功, 线路有功]
		% opf.Lz: 不等式约束对应的函数值 (实际值 - 松弛变量 - 下限) (n_ieq * 1)
		% 	[发电机有功函数值, 发电机无功函数值, 节点电压函数值, 线路有功函数值, 线路有功函数值]
		% opf.Lw: 不等式约束对应的函数值 (实际值 + 松弛变量 - 上限) (n_ieq * 1)
		% 	[发电机有功函数值, 发电机无功函数值, 节点电压函数值, 线路有功函数值, 线路有功函数值]
		% opf.Ll_mu: 衡量对偶间隙的一个变量 (n_ieq * 1)
		% opf.Lu_mu: 衡量对偶间隙的一个变量 (n_ieq * 1)
		% opf.Lx: 拉格朗日函数对自变量的梯度列向量 (n_state * 1)
		% opf.Lx_d: 对 opf.Lx 作了调整后的列向量
		% opf.parameterMatrix: 修正方程系数矩阵 ( (n_state + n_eq + 4 * n_ieq) * (n_state + n_eq + 4 * n_ieq) )
		% opf.H: 拉格朗日函数的海森矩阵 (n_state * n_state)
		% opf.H_d: opf.H 作了调整后的矩阵 (n_state * n_state)
		% opf.correction: 修正量 ( (n_state + n_eq + 4 * n_ieq) * 1 )
		% opf.delta_l: 不等式(下限)约束松弛变量修正量 (n_ieq * 1)
		% opf.delta_u: 不等式(上限)约束松弛变量修正量 (n_ieq * 1)
		% opf.delta_z: 不等式(下限)约束拉格朗日乘子修正量 (n_ieq * 1)
		% opf.delta_w: 不等式(下限)约束拉格朗日乘子修正量 (n_ieq * 1)
		% opf.alpha_p: 自变量和松弛变量步长, 最大 0.9995
		% opf.alpha_d: 拉格朗日乘子步长, 最大 0.9995
		% opf.delta_x: 自变量修正量 (n_state * 1)
		% opf.delta_y: 等式约束修正量 (n_eq * 1)

		% hvdc
		hvdc;

	end
	properties (Dependent)
		n_bus;
		n_branch;
		n_gen;

		pq;		% pq 节点索引
		n_pq;		% pq 节点数量
		pv;		% pv 节点索引
		n_pv;		% pv 节点数量
		slack;		% 平衡节点索引
		n_slack;		% 平衡节点索引
		
		load;		% 负荷节点索引
		n_load;		% 负荷节点数量

		pqGens;			% 处于 pq 节点的发电机索引
		pvGens;			% 处于 pv 节点的发电机索引
		slackGens;		% 处于平衡节点的发电机索引

		line;			% 普通线路索引(id)
		n_line;			% 普通线路数量
		trans;	% 变压器索引(id)
		n_trans;	% 变压器数量
	end
	methods

		%% 下面这些 getter 没有缓存, 之后将进行改进

		%% get.pq: 获取 pq 节点的索引
		function [pq] = get.pq(self)
			pq = find(self.bus.type == 1);
		end
		%% get.n_pq: 获取 pq 节点的数量
		function [n_pq] = get.n_pq(self)
			n_pq = sum(self.bus.type == 1);
		end
		%% get.pv: 获取 pv 节点的索引
		function [pv] = get.pv(self)
			pv = find(self.bus.type == 2);
		end
		%% get.n_pv: 获取 pv 节点的数量
		function [n_pv] = get.n_pv(self)
			n_pv = sum(self.bus.type == 2);
		end
		%% get.slack: 获取 ref 节点的索引
		function [slack] = get.slack(self)
			slack = find(self.bus.type == 3);
		end
		%% get.n_slack: 获取 ref 节点的数量
		function [n_slack] = get.n_slack(self)
			n_slack = sum(self.bus.type == 3);
		end
		%% get.load: 获取负荷节点的索引
		function [load] = get.load(self)
			load = find(self.bus.Pd > 0);
		end
		%% get.n_load: 获取负荷节点的数量
		function [n_load] = get.n_load(self)
			n_load = sum(self.bus.Pd > 0);
		end
		%% get.n_bus: 获取节点的数量
		function [n_bus] = get.n_bus(self)
			n_bus = length(self.bus.id);
		end
		
		%% get.n_gen: 获取节点的数量
		function [n_gen] = get.n_gen(self)
			n_gen = length(self.gen.id);
		end

		%% get.pqGens: 小型发电机数量
		function [pqGens] = get.pqGens(self)
			pqGens = sum(self.gen.type == 1);
		end
		%% get.pvGens: 中型发电机数量
		function [pvGens] = get.pvGens(self)
			pvGens = sum(self.gen.type == 2);
		end
		%% get.slackGens: 大型发电机数量
		function [slackGens] = get.slackGens(self)
			slackGens = sum(self.gen.type == 3);
		end

		%% get.line: 线路索引
		function [line] = get.line(self)
			line = find(self.branch.ratio==0);
		end
		%% get.n_line: 线路数量
		function [n_line] = get.n_line(self)
			n_line = sum(self.branch.ratio==0);
		end
		%% get.trans: 变压器索引
		function [trans] = get.trans(self)
			trans = find(self.branch.ratio~=0);
		end
		%% get.n_trans: 变压器数量
		function [n_trans] = get.n_trans(self)
			n_trans = sum(self.branch.ratio~=0);
		end
		%% get.n_branch: 变压器数量
		function [n_branch] = get.n_branch(self)
			n_branch = length(self.branch.id);
		end

		%% get.isLine: 返回该支路的类型, 是否为普通线路
		function [isLine] = isLine(self, id)
			% 支路的 id 也是 index
			isLine = self.branch.ratio(id) == 0;
		end

		%% get.Transformer: 返回该支路的类型, 是否为普通线路
		function [isTransformer] = isTransformer(self, id)
			% 支路的 id 也是 index
			isTransformer = self.branch.ratio(id) ~= 0;
		end

		%% init: 电力系统稳态模型初始化
		function init(self, mpc)

			if ~isstruct(mpc) && regexp(mpc, '^\w+$')
				eval(['mpc = ', mpc, '();']);
			end

			mpc.bus = sortrows(mpc.bus);
			self.mBase = mpc.baseMVA;
			self.bus.id = mpc.bus(:, 1);
			busLength = length(self.bus.id);
			self.bus.type = mpc.bus(:, 2);
			self.bus.g = mpc.bus(:,5)./self.mBase;
			self.bus.b = mpc.bus(:,6)./self.mBase;
			self.bus.mag = mpc.bus(:,8);
			self.bus.ang = mpc.bus(:,9).*pi./180;
			self.bus.Pg = zeros(busLength,1);
			self.bus.Qg = zeros(busLength,1);
			self.bus.Pg0 = zeros(busLength,1);
			self.bus.Qg0 = zeros(busLength,1);
			self.bus.Pd = mpc.bus(:, 3)./mpc.baseMVA;
			self.bus.Qd = mpc.bus(:, 4)./mpc.baseMVA;
			self.bus.Pdc = zeros(busLength,1);
			self.bus.Qdc = zeros(busLength,1);
			self.bus.Pc = zeros(busLength,1);
			self.bus.Qc = zeros(busLength,1);
			self.bus.Pis = zeros(busLength,1);
			self.bus.Qis = zeros(busLength,1);
			self.bus.Pout = zeros(busLength,1);
			self.bus.Qout = zeros(busLength,1);
			self.bus.dP = zeros(busLength,1);
			self.bus.dQ = zeros(busLength,1);
			self.bus.Qmax = zeros(busLength,1);
			self.bus.Qmin = zeros(busLength,1);
			self.bus.Vmax = mpc.bus(:, 12);
			self.bus.Vmin = mpc.bus(:, 13);

			genSize = size(mpc.gen);
			genLength = genSize(1);
			self.gen.nid = mpc.gen(:, 1);
			self.gen.id = (1:genLength)';
			self.gen.Pg = mpc.gen(:, 2)./mpc.gen(:, 7);	% Pg./mBase
			self.gen.Qg = mpc.gen(:, 3)./mpc.gen(:, 7);	% Qg./mBase
			self.gen.mBase = mpc.gen(:,7);
			self.gen.Qmax = mpc.gen(:,4)./self.mBase;
			self.gen.Qmin = mpc.gen(:,5)./self.mBase;
			self.gen.votage = mpc.gen(:,6);
			self.gen.status = mpc.gen(:,8);
			self.gen.Pmax = mpc.gen(:,9)./self.mBase;
			self.gen.Pmin = mpc.gen(:,10)./self.mBase;
			self.gen.type = zeros(genLength, 1);
			for k = 1:genLength
				self.gen.type(k) = self.bus.type(self.gen.nid(k) == self.bus.id);
			end

			branchSize = size(mpc.branch);
			branchLength = branchSize(1);
			self.branch.fid = mpc.branch(:, 1);
			self.branch.tid = mpc.branch(:, 2);
			self.branch.id = (1:branchLength)';
			self.branch.r = mpc.branch(:, 3);
			self.branch.x = mpc.branch(:, 4);
			self.branch.g = zeros(branchLength, 1);	% 目前我们不支持添加线路对地导纳
			self.branch.b = mpc.branch(:, 5);
			self.branch.ratio = mpc.branch(:, 9);
			self.branch.angle = mpc.branch(:, 10).*pi./360;
			self.branch.Pij = zeros(branchLength, 1);
			self.branch.Qij = zeros(branchLength, 1);
			self.branch.Pji = zeros(branchLength, 1);
			self.branch.Qji = zeros(branchLength, 1);
			self.branch.dP = zeros(branchLength, 1);
			self.branch.dQ = zeros(branchLength, 1);
			self.branch.status = mpc.branch(:, 11);
			self.branch.angmin = mpc.branch(:, 12).*pi./360;
			self.branch.angmax = mpc.branch(:, 13).*pi./360;

			self.branch.Pmax = mpc.branch(:, 6) ./ self.mBase;

			% 最优潮流
			if isfield(mpc, 'gencost')

				assert(size(mpc.gencost, 1) == genLength);

				assert(all(mpc.gencost(:, 1) == 2));
				% assert(all(mpc.gencost(:, 2) == 0));
				% assert(all(mpc.gencost(:, 3) == 0));
				assert(all(mpc.gencost(:, 4) <= 3));
				assert(all(mpc.gencost(:, 4) == mpc.gencost(1, 4)));

				self.gencost.model = mpc.gencost(:, 1);
				self.gencost.startup = mpc.gencost(:, 2);
				self.gencost.shutdown = mpc.gencost(:, 3);
				self.gencost.ncost = mpc.gencost(:, 4);

				self.gencost.costParams = mpc.gencost(:, 5:end) .* logspace((self.gencost.ncost(1)-1)*2, 0, self.gencost.ncost(1));

			end

			self.itlog.mag = [];
			self.itlog.ang = [];
			self.itlog.dP = [];
			self.itlog.dQ = [];

			self.NAM = [];
			self.FDM1 = [];
			self.FDM2 = [];

		end

		%% hvdcInit: hvdc 稳态分析模型初始化
		function hvdcInit(self, hvdc)

			assert(nargin == 2);
			self.hvdc = hvdc;

			self.hvdc.setBusIndex(self.bus);
		end

		%% getBusData: 生成节点参数, 主要是带有独立导纳设备的节点参数
		% 这里将节点上的并联电容视为恒阻抗模型,并将其归算至节点导纳矩阵
		function [busData] = getBusData(self)
			busData = [self.bus.id, self.bus.g, self.bus.b];
		end
		%% getBranchData: 生成线路参数, 包括线路和变压器
		function [branchData] = getBranchData(self)
			% 起始节点	终止节点	线路电阻	线路电抗	线路对地电导	线路对地电纳	变比
			branchData = [self.branch.fid, self.branch.tid, self.branch.r, self.branch.x, self.branch.g, self.branch.b, self.branch.ratio, self.branch.angle];
		end

		%% NAMInit: 计算节点导纳矩阵, 结果存放于电力网相应属性中
		% @return void
		function NAMInit(self)
			self.NAM = Model.NAM();
			self.NAM.generate(self.getBusData(), self.getBranchData());
		end
		
		%% votageInit: 获取迭代初始值(考虑对节点的设置及发电机的设置), 结果存放于 bus 字段相应属性中
		% @return void
		function votageInit(self, solver)

			startMode = '';

			if isstruct(solver) && isfield(solver, 'start') && regexp(solver.start, '^\w+$')
				startMode = solver.start;
			end

			switch startMode
				% 0-1 启动
				case {'flat', 'Flat', 'FLAT', 'worm', 'Worm', 'WORM'}
					self.bus.mag(self.pq) = 1;
					self.bus.ang([self.pq; self.pv]) = 0;
				% 直流潮流启动
				case {'dc', 'DC'}
					self.bus.mag(self.pq) = 1;
					self.setFDM('FD');
					self.planUpdate();
					self.bus.ang(self.bus.type ~= 3) = self.FDM1 \ self.bus.Pis(self.bus.type ~= 3);
				otherwise
					self.bus.mag(self.pq) = 1;
					self.bus.ang([self.pq; self.pv]) = 0;
			end

			self.bus.mag0 = self.bus.mag;
			self.bus.ang0 = self.bus.ang;
		end

		%% planInit: 初始化,用于计算各节点功率计划值,并根据发电机的情况计算出该节点可发出的的最大无功功率, 结果存放于 bus 字段相应属性中
		% @return void
		function planInit(self)

			self.bus.Pg = zeros(size(self.bus.Pg));
			self.bus.Qg = zeros(size(self.bus.Qg));

			for k = 1:length(self.gen.id)
				index = find(self.bus.id == self.gen.nid(k));
				% 这里对各节点的功率计划的赋值有一大部分(PQ节点)是没有意义的
				self.bus.Pg(index) = self.bus.Pg(index) + self.gen.Pg(k);
				self.bus.Qg(index) = self.bus.Qg(index) + self.gen.Qg(k);
				% self.bus.Pmin(index) = self.bus.Pmin(index) + self.gen.Pmin(k);
				% self.bus.Pmax(index) = self.bus.Pmax(index) + self.gen.Pmax(k);
				% self.bus.Qmin(index) = self.bus.Qmin(index) + self.gen.Qmin(k);
				% self.bus.Qmax(index) = self.bus.Qmax(index) + self.gen.Qmax(k);
			end

			self.bus.Pg0 = self.bus.Pg;
			self.bus.Qg0 = self.bus.Qg;
		end

		%% planUpdate: 更新各节点功率计划值, 结果存放于 bus 字段相应属性中
		% @return void
		function planUpdate(self)
			self.bus.Pis = self.bus.Pg + self.bus.Pdc - self.bus.Pd;
			self.bus.Qis = self.bus.Qg + self.bus.Qdc - self.bus.Qd;
		end

		%% setPowerOutflow: 计算从 nid 节点注入电网的潮流, 结果存放于 bus 字段相应属性中
		% @param  int32   nid   节点 id
		% @return void
		function setPowerOutflow(self, nid)

			if nargin == 1
				votage = self.bus.mag .* exp(1i .* self.bus.ang);
				Sf = votage .* (conj(self.NAM.value * votage));
				self.bus.Pout = real(Sf);	
				self.bus.Qout = imag(Sf);
				return;
			end

			for k = 1:length(nid)
				index = find(self.bus.id == nid(k));
				Sf = conj(self.NAM.value(index, :))*(self.bus.mag.*exp(1i.*(self.bus.ang(index)-self.bus.ang))).*self.bus.mag(index);
				self.bus.Pout(k) = real(Sf);
				self.bus.Qout(k) = imag(Sf);
			end
		end

		%% setPowerUnbalance: 计算功率不平衡量, 仅用于牛拉法和 PQ 分解法计算, 返回变量均涵盖所有节点并已分类, 结果存放于 bus 字段相应属性中
		% @return void
		function setPowerUnbalance(self)
			self.bus.dP = self.bus.Pis - self.bus.Pout;
			self.bus.dQ = self.bus.Qis - self.bus.Qout;
		end

		%% setIterationData: 每次使用潮流方程之后调用此方法获取迭代信息, 结果存放于 bus 字段相应属性中
		% @return void
		function setIterationData(self)
			% TODO 优化内存使用方式
			self.itlog.mag = [self.itlog.mag, self.bus.mag];
			self.itlog.ang = [self.itlog.ang, self.bus.ang];
			self.itlog.dP = [self.itlog.dP, self.bus.dP];
			self.itlog.dQ = [self.itlog.dQ, self.bus.dQ];
		end

		%% checkConverged: 判断是否收敛, 同时检查是否结束循环
		% @param  struct  solver         与判断收敛相关的设置
		% --attr  double  epsilon        允许的功率误差
		% --attr  int32   n_iters_max   最大迭代次数
		% @return int32   status       判断结果, 0 表示未收敛, 1 表示已收敛
		function [status] = checkConverged(self, solver, it)
			% 收敛判据
			if solver.epsilon >= max(abs([self.bus.dQ; self.bus.dP]))
				status = 1;	% 收敛
			elseif it >= solver.n_iters_max
				warning('Exception 101: Number of iterations exceeds the limit');
				status = 101;	% 迭代次数超出上限
			% elseif any(self.bus.mag > 2.5) || any(self.bus.mag < 0.4)
			% 	warning('Exception 102: Some node have abnormal voltages and iteration is terminated');
			% 	status = 102;	% 部分节点电压不正常
			else
				status = 0;	% 未收敛, 继续计算
			end
		end

		%% getJacobian: 在当前网络状态下获取雅可比矩阵,用于牛顿法的计算
		% @return n*n-double jacobianMatrix   雅克比矩阵
		function [jacobianMatrix] = getJacobian(self)

			% 为方便构造雅可比矩阵,现分别对 PQ 节点及 PV 节点提取出节点电压及节点导纳, 并计算节点之间的相角矩阵
			% 本方法中, 变量名以 1 结尾的均与雅可比矩阵中的 H N M 有关, 变量名以 2 结尾的均与雅可比矩阵中的 N M L 有关
			NAM1 = self.NAM.value([self.pq; self.pv], [self.pq; self.pv]);
			NAM2 = self.NAM.value([self.pq], [self.pq]);

			mag1 = self.bus.mag([self.pq; self.pv]);
			ang1 = self.bus.ang([self.pq; self.pv]);
			mag2 = self.bus.mag([self.pq]);
			ang2 = self.bus.ang([self.pq]);

			nodeAngle = ang1 * ones(1, length(ang1));
			nodeAngle = nodeAngle - nodeAngle';

			% 开始计算雅可比矩阵,创建初始状态. jacPrimary 是一个复数矩阵, 雅克比矩阵的元素主要来自这个矩阵
			jacPrimary = mag1 * (mag1') .* conj(NAM1) .* exp(1i .* nodeAngle);
			% 这个结果的对角元无实际意义, 去掉对角元
			jacPrimary = jacPrimary - diag(diag(jacPrimary));

			% 先计算各个分块矩阵的非对角元
			% 需要注意的是, 下面的两个索引表示在雅可比矩阵的初始状态中 PQPV 节点及 PQ 节点所在的位置, 而 self.pq 表示在节点对象中 PQ 节点所在的位置
			pqpv = 1:(self.n_pq + self.n_pv);
			pq = 1:self.n_pq;

			% 雅可比矩阵的非对角元均可以在初始状态矩阵中提取
			Hij = imag(jacPrimary([pqpv], [pqpv]));	% Hij 表示系统中有功潮流差与相角差的数量关系
			Nij = real(jacPrimary([pqpv], [pq]));	% Nij 表示系统中有功潮流差与电压差的相对值的数量关系
			Mij = -real(jacPrimary([pq], [pqpv]));	% Mij 表示系统中无功潮流差与相角差的数量关系
			Lij = imag(jacPrimary([pq], [pq]));	% Lij 表示系统中无功潮流差与电压差的相对值的数量关系

			% 再计算各个分块矩阵的对角元. 求导的原因, 使得对角元需要单独计算
			Hii = -self.bus.Qout([self.pq; self.pv]) - mag1.^2.*imag(diag(NAM1));
			Nii = self.bus.Pout(self.pq) + mag2.^2.*real(diag(NAM2));
			Mii = self.bus.Pout(self.pq) - mag2.^2.*real(diag(NAM2));
			Lii = self.bus.Qout(self.pq) - mag2.^2.*imag(diag(NAM2));

			Hij = Hij + diag(Hii);
			Nij(pq, pq) = Nij(pq, pq) + diag(Nii);
			Mij(pq, pq) = Mij(pq, pq) + diag(Mii);
			Lij = Lij + diag(Lii);

			% 求解普通潮流用的雅可比矩阵
			jacobianMatrix = [Hij, Nij; Mij, Lij;];
		end

		%% getOpfJacobian: 获取真正意义的潮流方程雅可比矩阵
		% @return n*n-double jacobianMatrix   雅克比矩阵
		function [jacobianMatrix] = getOpfJacobian(self)

			n_bus = self.n_bus;

			% 开始计算雅可比矩阵,创建初始状态. jacPrimary 是一个复数矩阵, 雅克比矩阵的元素主要来自这个矩阵
			jacPrimary = self.bus.mag .* self.opf.Yed .* self.bus.mag';
			% 这个结果的对角元无实际意义, 去掉对角元
			diag_index = sub2ind(size(jacPrimary), 1:n_bus, 1:n_bus);
			% jacPrimary(diag_index) = 0;

			% 雅可比矩阵的非对角元均可以在初始状态矩阵中提取
			Hij = imag(jacPrimary);	% Hij 表示系统中有功潮流差与相角差的数量关系
			Nij = real(jacPrimary);	% Nij 表示系统中有功潮流差与电压差的相对值的数量关系
			Mij = -Nij;	% Mij 表示系统中无功潮流差与相角差的数量关系
			Lij = Hij;	% Lij 表示系统中无功潮流差与电压差的相对值的数量关系

			% 再计算各个分块矩阵的对角元. 求导的原因, 使得对角元需要单独计算
			Hij(diag_index) = -self.bus.Qout - self.bus.mag .^ 2 .* imag(diag(self.NAM.value));
			Nij(diag_index) = self.bus.Pout + self.bus.mag .^ 2 .* real(diag(self.NAM.value));
			Mij(diag_index) = self.bus.Pout - self.bus.mag .^ 2 .* real(diag(self.NAM.value));
			Lij(diag_index) = self.bus.Qout - self.bus.mag .^ 2 .* imag(diag(self.NAM.value));

			% 求解普通潮流用的雅可比矩阵
			jacobianMatrix = [Hij, Nij ./ self.bus.mag'; Mij, Lij ./ self.bus.mag';]';
		end

		%% calFDM: 返回 PQ 分解法的 B 矩阵
		% @param  string     fdType    PQ 分解法的系数矩阵类型(FD, FDBX, FDXB)
		% @param  array      busType 矩阵中包含的节点类型, 要得到 B' 时传 [1, 2], 要得到 B'' 时传 [1]
		% @return n*n double 系数矩阵
		function [matrix] = calFDM(self, fdType, busType)

			assert(ismember(1, busType) && ~ismember(3, busType));

			matrix = sparse(zeros(self.n_pq + self.n_pv + self.n_slack));
			switch fdType
				case 'r'	% 考虑线路电阻
					for k = 1:length(self.branch.id)
						fi = find(self.bus.id == self.branch.fid(k));
						ti = find(self.bus.id == self.branch.tid(k));
						b = self.branch.x(k)./(self.branch.r(k).^2 + self.branch.x(k).^2);
						matrix([fi, ti], [fi, ti]) = matrix([fi, ti], [fi, ti]) + [b, -b; -b, b];
					end
				case 'b'	% 考虑节点电纳
					for k = 1:length(self.branch.id)
						fi = find(self.bus.id == self.branch.fid(k));
						ti = find(self.bus.id == self.branch.tid(k));
						b = 1./self.branch.x(k);
						matrix([fi, ti], [fi, ti]) = matrix([fi, ti], [fi, ti]) + [b, -b; -b, b];
					end
					matrix = matrix - diag(self.bus.b);
				otherwise
					error('illegal fast duplicate type in calFDM');
			end

			index = find(self.bus.type == 1);
			if ismember(2, busType)
				index = [index; find(self.bus.type == 2)];
			end

			matrix = -matrix(index, index);
		end

		%% setFDM: 计算 PQ 分解法的系数矩阵
		% @param  string     fdType    PQ 分解法的系数矩阵类型(FD, FDBX, FDXB)
		% @return  void
		function setFDM(self, fdType)
			switch fdType
				case 'FD'
					assert(length(self.NAM.value) ~= 0);
					self.FDM1 = imag(self.NAM.value([self.pq; self.pv],[self.pq; self.pv]));
					self.FDM2 = imag(self.NAM.value([self.pq], [self.pq]));
				case 'FDBX'
					self.FDM1 = self.calFDM('r', [1, 2]);	% 考虑线路电阻
					self.FDM2 = self.calFDM('b', [1]);		% 考虑节点电纳
				case 'FDXB'
					self.FDM1 = self.calFDM('b', [1, 2]);
					self.FDM2 = self.calFDM('r', [1]);
				otherwise
					error('illegal fast duplicate type in setFDM');
			end
		end

		%% getVotageAndAngleCorrection: 修正方程, 用于牛顿法
		%　@param  n*n-double  jacobianMatrix   雅克比矩阵
		% @return double      dTheta           相角修正量, 用于迭代方程
		% @return double      dUpn             电压修正量标幺值, 用于迭代方程
		function [dTheta,dUpn] = getVotageAndAngleCorrection(self, jacobianMatrix, standardFlag)

			dP = self.bus.dP([self.pq; self.pv]);
			dQ = self.bus.dQ([self.pq]);

			dPara = jacobianMatrix \ [dP; dQ];

			dTheta = dPara(1:(self.n_pq + self.n_pv));	% PQ 及 PV 节点相角的修正量
			dUpn = dPara((self.n_pq + self.n_pv + 1):end);	% PQ 节点电压相对值的修正量

			if nargin == 3 && standardFlag == true
				% 将电压和相角误差规范化, 即填充至大向量中, 方便迭代方程直接使用
				dUpn_t = zeros(self.n_pq + self.n_pv + self.n_slack, 1);	% 初始化相角修正量
				dUpn_t([self.pq]) = dUpn;	% 代入 PQ 及 PV 节点相角的修正量
				dUpn = dUpn_t;
				dTheta_t = zeros(self.n_pq + self.n_pv + self.n_slack, 1);
				dTheta_t([self.pq; self.pv]) = dTheta;
				dTheta = dTheta_t;
			end
		end

		%% getMagCorrection: 计算电压修正量, 用于 PQ 分解法
		function [dUpn] = getMagCorrection(self, standardFlag)
			if nargin == 2 && standardFlag == true
				dUpn = zeros(self.n_pq + self.n_pv + self.n_slack, 1);
				dUpn([self.pq]) = self.FDM2 \ (self.bus.dQ([self.pq])./self.bus.mag([self.pq]));
			else
				dUpn = self.FDM2 \ (self.bus.dQ([self.pq])./self.bus.mag([self.pq]));
			end
		end

		%% getAngCorrection: 计算相角修正量, 用于 PQ 分解法
		function [dTheta] = getAngCorrection(self, standardFlag)
			if nargin == 2 && standardFlag == true
				dTheta = zeros(self.n_pq + self.n_pv + self.n_slack, 1);
				dTheta([self.pq; self.pv]) = (self.FDM1 \ (self.bus.dP([self.pq; self.pv])./self.bus.mag([self.pq; self.pv])))./self.bus.mag([self.pq; self.pv]);
			else
				dTheta = (self.FDM1 \ (self.bus.dP([self.pq; self.pv])./self.bus.mag([self.pq; self.pv])))./self.bus.mag([self.pq; self.pv]);
			end
		end

		%% setNewVotage: 获取新的节点电压及相角, 用于牛顿法. 更新电力网 bus 字段的属性 mag 和 ang
		% @param  double dTheta    相角修正量, 来自修正方程
		% @param  double dUpn      电压修正量标幺值, 来自修正方程
		% @return void
		function setNewVotage(self, dTheta, dUpn, standardFlag)
			if nargin == 4 && standardFlag == true
				self.bus.ang = self.bus.ang + (dTheta);
				self.bus.mag = self.bus.mag.*(1+dUpn);
			else
				self.bus.ang([self.pq; self.pv]) = self.bus.ang([self.pq; self.pv]) + (dTheta);
				self.bus.mag([self.pq]) = self.bus.mag([self.pq]).*(1+dUpn);
			end
		end

		%% setGeneratorPower: 计算各节点的实际发电量 (仅计算了 pv 节点的无功和平衡节点的有功无功)
		%  结果存放于电力网 node 字段的属性 Pg 和 Qg
		% @return void
		function setGeneratorPower(self)
			% 每个节点发电机发出的功率都等于当地需求和收敛时潮流方程计算的注入电网的功率的和
			self.bus.Pg(self.slack) = self.bus.Pout(self.slack) - self.bus.Pdc(self.slack) + self.bus.Pd(self.slack);
			self.bus.Qg([self.pv; self.slack]) = self.bus.Qout([self.pv; self.slack]) - self.bus.Qdc([self.pv; self.slack]) + self.bus.Qd([self.pv; self.slack]);
		end

		%% getLinePower: 计算普通线路的潮流
		% @param  int32   index    线路 id(index)
		% @return complex Sij      从起始节点到终止节点的复功率, 标幺值
		% @return complex Sji      从终止节点到起始节点的复功率, 标幺值
		% @return complex dS       线路中的损耗, 但线路充电电容的影响也考虑了进去. 标幺值
		function [Sij,Sji,dS] = getLinePower(self, index, shuntFlag)
			%　计算三个导纳值，计算功率时会用到
			if shuntFlag
				conjYi0 = conj((self.branch.g(index) + self.branch.b(index) * 1i ) ./ 2 );
				conjYj0 = conjYi0;
			else
				conjYi0 = 0;
				conjYj0 = 0;
			end
			conjYij = 1 ./ conj((self.branch.r(index) + self.branch.x(index) * 1i));

			fid = find(self.bus.id == self.branch.fid(index));	% 得到线路始末端节点的id
			tid = find(self.bus.id == self.branch.tid(index));	% 得到线路始末端节点的id

			% 使用公式计算
			Sij = (self.bus.mag(fid)).^2.*(conjYi0 + conjYij) - self.bus.mag(fid).*self.bus.mag(tid).*conjYij.*exp((self.bus.ang(fid)-self.bus.ang(tid)).*i);
			Sji = (self.bus.mag(tid)).^2.*(conjYj0 + conjYij) - self.bus.mag(fid).*self.bus.mag(tid).*conjYij.*exp((self.bus.ang(tid)-self.bus.ang(fid)).*i);
			dS = Sij + Sji;
		end

		%% getTransformerPower: 计算变压器的潮流
		% @param  int32   index    线路 id(index)
		% @return complex Sij      从起始节点到终止节点的复功率, 标幺值
		% @return complex Sji      从终止节点到起始节点的复功率, 标幺值
		% @return complex dS       线路中的损耗, 但线路充电电容的影响也考虑了进去. 标幺值
		function [Sij,Sji,dS] = getTransformerPower(self, index, shuntFlag)
			%　首先确定变压器　pi 型等效电路的三个参数
			[trZ, trY1, trY2] = gamma2pi((self.branch.r(index)+self.branch.x(index)*1i), (self.branch.g(index)+self.branch.b(index)*1i), self.branch.ratio(index));
			%　计算三个导纳值，计算功率时会用到
			if shuntFlag
				conjYi0 = conj(trY1) - (self.branch.g(index)+self.branch.b(index)*1i)./2;
				conjYj0 = conj(trY2) - (self.branch.g(index)+self.branch.b(index)*1i)./2;
			else
				conjYi0 = 0;
				conjYj0 = 0;
			end
			conjYij = conj(1./trZ);

			fid = find(self.bus.id == self.branch.fid(index));	% 得到线路始始端节点的id
			tid = find(self.bus.id == self.branch.tid(index));	% 得到线路始末端节点的id
			
			% 使用公式计算
			Sij = (self.bus.mag(fid)).^2.*(conjYi0 + conjYij) - self.bus.mag(fid).*self.bus.mag(tid).*conjYij.*exp((self.bus.ang(fid)-self.bus.ang(tid)).*1i);
			Sji = (self.bus.mag(tid)).^2.*(conjYj0 + conjYij) - self.bus.mag(fid).*self.bus.mag(tid).*conjYij.*exp((self.bus.ang(tid)-self.bus.ang(fid)).*1i);
			dS = Sij + Sji;
		end
		
		%% setBranchPower: 计算线路功率及功率损耗
		%  该方法在潮流迭代结束之后再使用, 用于求解收敛时电力网中各线路传输的功率, 求解完成后更新电力网 branch 字段的相应属性
		% @return void
		function setBranchPower(self, shuntFlag)

			if nargin == 1
				shuntFlag = true;
			end

			Sij = zeros(length(self.branch.id), 1);
			Sji = zeros(length(self.branch.id), 1);
			dS = zeros(length(self.branch.id), 1);
			for k = 1:length(self.branch.id)
				if self.isLine(k)
					[Sij(k), Sji(k), dS(k)] = self.getLinePower(k, shuntFlag);
				else
					[Sij(k), Sji(k), dS(k)] = self.getTransformerPower(k, shuntFlag);
				end	
			end
			self.branch.Pij = real(Sij);
			self.branch.Qij = imag(Sij);
			self.branch.Pji = real(Sji);
			self.branch.Qji = imag(Sji);
			self.branch.dP = real(dS);
			self.branch.dQ = imag(dS);
		end

		%% setConpensatorPower: 计算无功补偿器电量
		% @return void
		function setConpensatorPower(self)
			self.bus.Pc = self.bus.mag .^ 2 .* self.bus.g;
			self.bus.Qc = self.bus.mag .^ 2 .* self.bus.b;
		end

		%% setVotageRealAndImag: 将节点电压转化为直角坐标形式, 存储在 bus.Ux, bus.Uy 中
		% @return void
		function setVotageRealAndImag(self)
			self.bus.Ux = self.bus.mag .* cos(self.bus.ang);
			self.bus.Uy = self.bus.mag .* sin(self.bus.ang);
		end

		%% NRIteration: 牛顿法迭代程序核心部分
		% @param    struct  solver	        对迭代程序的基本设置
		% --attr    string  method            求解方法, NR FDBX FDXB
		% --attr    int32   n_iters_max      最大迭代次数
		% --attr    double  epsilon           收敛判据, 功率不平衡量标幺
		% --attr    string  start            
		% --attr    string  documentName      如果是文本输出, 需设置文本计算报告 src
		% @return   struct  result          迭代结果
		% --attr    int32   status            错误码, 为 0 表示正常
		% --attr    int32   it                迭代次数
		function [result] = NRIteration(self, solver)

			% 从这里开始进入循环,跳出循环的条件为迭代次数超过最大次数或潮流不平衡量的无穷范数小于给定值
			result.it = 0;
			while 1
				self.planUpdate();	% 确定各节点功率需求
				self.setPowerOutflow();	% 潮流方程
				self.setPowerUnbalance();	% 误差方程

				self.setGeneratorPower();
				self.setIterationData();	% 存储迭代信息,包括电压,相角,有功不平衡量,无功不平衡量
				
				% 判收敛以及是否迭代次数上限, 返回一个状态码, 得到非零值直接返回
				result.status = self.checkConverged(solver, result.it);
				if(result.status ~= 0)	% result.status 非 0 表示收敛或错误
					return;
				end

				J = self.getJacobian();	% ❀ 将各节点电压的初值代入修正方程,求解系数矩阵.
				[dTheta, dUpn] = self.getVotageAndAngleCorrection(J);	% 修正方程
				self.setNewVotage(dTheta, dUpn);	% 迭代方程

				result.it = result.it + 1;
				% ❀ 运用各节点电压的新值进行下一步的迭代.
			end
		end

		%% PQIteration: PQ 分解法迭代程序核心部分
		% @param    struct  solver	        对迭代程序的基本设置
		% --attr    string  method            求解方法, FD FDBX FDXB
		% --attr    int32   n_iters_max      最大迭代次数
		% --attr    double  epsilon           收敛判据, 功率不平衡量标幺
		% --attr    string  start             启动方式, default 为按发电机端电压起动
		% --attr    string  documentName      如果是文本输出, 需设置文本计算报告 src
		% @return   struct  result          迭代结果
		% --attr    int32   status            错误码, 为 0 表示正常
		% --attr    int32   it                迭代次数
		function [result] = PQIteration(self, solver)
			% 计算 PQ 分解法的两个系数矩阵
			self.setFDM(solver.method);
			% 从这里开始进入循环, 跳出循环的条件为迭代次数超过最大次数或潮流不平衡量的一范数小于给定值
			result.it = 0;

			while 1
				self.planUpdate();	% 确定各节点功率需求
				self.setPowerOutflow();	% 潮流方程
				self.setPowerUnbalance();	% 误差方程

				self.setGeneratorPower();
				self.setIterationData();	% 存储迭代信息,包括电压,相角,有功不平衡量,无功不平衡量

				% 判收敛以及是否迭代次数上限, 返回一个状态码, 得到非零值直接返回
				result.status = self.checkConverged(solver, result.it);
				if(result.status ~= 0)	% result.status 非 0 表示收敛或错误
					return;
				end

				dUpn = self.getMagCorrection();
				self.setNewVotage(0, -dUpn);	% 迭代方程

				dTheta = self.getAngCorrection();
				self.setNewVotage(-dTheta, 0);	% 迭代方程
				
				result.it = result.it + 1;
				% ❀ 运用各节点电压的新值进行下一步的迭代.
			end
		end
		
		%% checkReactivePower: 考察各节点的发电机无功功率是否越限, 并采取相应措施(节点类型转化)
		function [output] = checkReactivePower(self)

			% 下面两句 返回的 busIndex 包括所有类型的节点, 但无需剔除不必要的节点
			lack = find(self.bus.Qmax < self.bus.Qg);
			excess = find(self.bus.Qmin > self.bus.Qg);

			% 只保留 PV 节点
			lack(find(self.bus.type(lack) ~=2 )) = [];
			excess(find(self.bus.type(excess) ~=2 )) = [];

			if isempty([lack; excess])
				output = false;
				return;
			end
			output = true;

			self.convertToPQNode(lack, excess);

			for k = 1:length(lack)
				fprintf('%s %d %s\n', '	  Node ', self.bus.id(lack(k)),' has insufficient reactive power and has been converted to a PQ node');
			end
			for k = 1:length(excess)
				fprintf('%s %d %s\n', '	  Node ', self.bus.id(excess(k)),' has excess reactive power and has been converted to a PQ node');
			end
		end

		%% convertToPQNode: 解决某些节点无功不足或过剩的情况
		function convertToPQNode(self, lack, excess)

			self.bus.Qg(lack) = self.bus.Qmax(lack);
			self.bus.Qg(excess) = self.bus.Qmin(excess);

			% 修改节点类型,将其转化为 PQ 节点,并更新节点对象相应的的 PQ PV 节点的索引及数量属性
			self.bus.type([lack; excess]) = 1;

		end

		%% setSourceFromHVDC: 根据直流输电系统的功率设置交流系统节点注入功率
		function setSourceFromHVDC(self)

			n_hvdc = length(self.hvdc.f);

			self.bus.Pg = self.bus.Pg0;
			self.bus.Pg = self.bus.Pg0;

			self.bus.Pdc = zeros(self.n_bus, 1);
			self.bus.Qdc = zeros(self.n_bus, 1);
			for k = 1:n_hvdc
				self.bus.Pdc(self.hvdc.f(k)) = self.bus.Pdc(self.hvdc.f(k)) - self.hvdc.Pr(k);
				self.bus.Pdc(self.hvdc.t(k)) = self.bus.Pdc(self.hvdc.t(k)) + self.hvdc.Pi(k);
				self.bus.Qdc(self.hvdc.f(k)) = self.bus.Qdc(self.hvdc.f(k)) - self.hvdc.Qr(k);
				self.bus.Qdc(self.hvdc.t(k)) = self.bus.Qdc(self.hvdc.t(k)) - self.hvdc.Qi(k);
			end
		end

		%% solvePowerFlow: 电力系统潮流计算
		function [result] = solvePowerFlow(self, solver)
		    assert(nargin == 2);

		    % ❀ 生成节点导纳矩阵.
		    self.NAMInit();
		    % ❀ 设置节点电压的初值(电压,相角).
		    self.votageInit(solver);
		    % ❀ 将各节点电压的初值代入潮流方程,并求解有功无功不平衡量
		    self.planInit();    % 计算各节点功率的计划值
		    
	        if isobject(self.hvdc)
	            self.hvdc.render(self.bus.mag(self.hvdc.f), self.bus.mag(self.hvdc.t));
	            self.setSourceFromHVDC();
	        end
		    % 这个循环是与检测无功功率或与直流输电有关的, 与具体的求解无关
		    while 1


		        % ❀ 进入迭代. 得到的 result 包含是否收敛, 迭代次数等信息
		        switch solver.method
		            case 'NR'   % 牛顿法
		                result = self.NRIteration(solver);
		            case {'FD', 'FDBX', 'FDXB'} % PQ 分解法
		                result = self.PQIteration(solver);
		            otherwise
		                error('illegal solver');
		        end

		        % 未收敛时直接返回
		        if(result.status ~= 1)
		            return;
		        end

		        self.setGeneratorPower();   % 根据迭代结果计算发电机功率, 结果记录到 bus 属性中
		        
		        % if isfield(solver, 'checkReactivePower') && solver.checkReactivePower == true && self.checkReactivePower()
		        %     continue;
		        % end

		        if isobject(self.hvdc) && ~self.hvdc.render(self.bus.mag(self.hvdc.f), self.bus.mag(self.hvdc.t))
		        	self.setSourceFromHVDC();
		        	disp([self.bus.mag(4:5)', self.hvdc.Pr, self.hvdc.Pi, self.hvdc.Qr, self.hvdc.Qi]);
		            continue;
		        end
		        break;
		    end

		    self.setBranchPower(true);  % 根据迭代结果计算线路功率及损耗, 结果记录到 branch 属性中
		    self.setConpensatorPower(); % 根据迭代结果计算各节点对地电容功率, 结果记录到 bus 属性中
		end

		% addGenImpedanceToNodes: 将发电机的暂态电抗加至各节点, 为计算暂态参数做准备
		% 不是给潮流分析留的接口
		function addGenImpedanceToNodes(self, gen)
			if nargin == 1
				return;
			end

			index = getIndex(self.bus.id, gen.nid);
			yg = 1./(gen.Ra + gen.Xd1.*1i);
			
			self.bus.g(index) = self.bus.g(index) + real(yg);
			self.bus.b(index) = self.bus.b(index) + imag(yg);
			self.NAM.addAdmittance(index, yg);
		end

		% addMoterImpedanceToBus: 将异步电动机的暂态电抗加至各节点, 为计算暂态参数做准备
		% 不是给潮流分析留的接口
		function addMoterImpedanceToBus(self, moter)
			% TODO addMoterImpedanceToBus
			if nargin == 1
				return;
			end
		end

		% addLoadImpedanceToBus: 将负荷的电抗加至节点, 为计算暂态参数做准备
		% 不是给潮流分析留的接口
		function addLoadImpedanceToBus(self)
			self.bus.g = self.bus.g + self.bus.Pd./self.bus.mag.^2;
			self.bus.b = self.bus.b - self.bus.Qd./self.bus.mag.^2;
			self.NAM.addAdmittance(1:self.n_bus, (self.bus.Pd - self.bus.Qd.*1i)./self.bus.mag.^2);
		end

		%% changeLoadCapacity: 改变所有负荷的容量
		function changeLoadCapacity(self, rate)
			self.bus.Pd = self.bus.Pd.*rate;
			self.bus.Qd = self.bus.Qd.*rate;
		end

		%% changeLoadCapacity: 改变所有负荷有功
		function changeLoadActiveCapacity(self, rate)
			self.bus.Pd = self.bus.Pd.*rate;
		end

		%% changeLoadFactorPerDeg: 改变所有负荷的功率因数
		function changeLoadFactorPerDeg(self, hysteresisDeg)
			PM = [
				cos(hysteresisDeg./180.*pi), -sin(hysteresisDeg./180.*pi);
				sin(hysteresisDeg./180.*pi), cos(hysteresisDeg./180.*pi);
			];
			[self.bus.Pd, self.bus.Qd] = PM*[self.bus.Pd, self.bus.Qd];
		end

		%% getLoadFactor: 计算负荷的功率因数, 并返回负荷的性质(感性0, 容性1)
		function [loadFactor, nature] = getLoadFactor(self)
			loadFactor = ones(self.n_bus, 1);
			nature = zeros(self.n_bus, 1);
			loadFactor(self.load) = self.bus.Pd(self.load)./sqrt((self.bus.Pd(self.load) + self.bus.g(self.load)).^2+(self.bus.Qd(self.load) - self.bus.b(self.load)).^2);
			nature(self.load) = (self.bus.Qd(self.load) - self.bus.b(self.load) < 0);
		end

		%% compensateReactivePower: 对所有负荷使用调相机进行无功补偿,用于测试系统对负荷功率因数的反应.
		% TODO add bus.b
		function compensateReactivePowerByCompensator(self, LF_obj)
			tanphi_t = tan(acos(LF_obj));
			[LF_org, LF_nat] = self.getLoadFactor();
			% 不能使用 for k = self.load
			for k = self.load'
				if( (LF_org(k) < LF_obj) & (LF_nat(k) == 0) )
					self.bus.Qd(k) = (self.bus.Pd(k) + self.bus.g(k)).*tanphi_t + self.bus.b(k);
				end
			end
		end

		%% compensateReactivePowerByCapacitance: 对所有负荷使用电容进行无功补偿,用于测试系统对负荷功率因数的反应.
		% TODO add bus.b
		function compensateReactivePowerByCapacitance(self, LF_obj)
			tanphi_t = tan(acos(LF_obj));
			[LF_org,LF_nat] = self.getLoadFactor();
			% 不能使用 for k = self.load
			for k = self.load'
				if( (LF_org(k) < LF_obj) & (LF_nat(k) == 0) )
					comp_t = self.bus.Qd(k) - self.bus.b(k) - (self.bus.Pd(k) + self.bus.g(k)).*tanphi_t;	% 待补偿量
					if 0 < comp_t
						self.bus.b(k) = self.bus.b(k) + comp_t;
					end
				end
			end
		end

		%% solveOptimalPowerFlow: 最优潮流
		function [result] = solveOptimalPowerFlow(self, solver)

			assert(isstruct(solver));

			%% 初始化
			self.opfInit(solver);

			self.opf.n_iters = 0;
			while true

				%% 计算互补间隙
				self.calcGap();

				result.status = self.opfConverged();

				%% 收敛判据(互补间隙满足要求)
				if result.status ~= 0
					break;
				end

				%% 计算扰动因子 μ
				self.opf.mu = self.opf.sigma .* self.opf.gap ./ (2 .* self.opf.n_ieq);
				assert(self.opf.mu > 0);

				%% 求解修正方程
				self.calcOpfCorrection();

				%% 计算步长
				self.calcOpfStep();

				%% 修正
				self.setNewOpfParameters();
				self.opf.n_iters = self.opf.n_iters + 1;
			end

			self.planUpdate();
			self.setPowerOutflow();
			self.setPowerUnbalance();
			self.setBranchPower(true);

		end

		%% opfInit: 最优潮流计算初始化
		function opfInit(self, solver)
			% 各种约束和变量的个数
			self.opf.n_eq = 2 .* self.n_bus;
			self.opf.n_state = 2 .* self.n_gen + self.opf.n_eq;
			% self.opf.n_ieq = 2 .* self.n_gen + self.n_bus + self.n_branch;
			self.opf.n_ieq = 2 .* self.n_gen + self.n_bus + 2 .* self.n_branch;

			% ❀ 生成节点导纳矩阵.
			self.NAMInit();
			% ❀ 设置节点电压的初值(电压,相角).
			self.votageInit('');

			%% 变量的初始化方案参考教材
			self.gen.Pg = (self.gen.Pmax + self.gen.Pmin) ./ 2;
			self.gen.Qg = (self.gen.Qmax + self.gen.Qmin) ./ 2;
			self.opf.x = [self.gen.Pg; self.gen.Qg; self.bus.ang; self.bus.mag];

			% ❀ 将各节点电压的初值代入潮流方程,并求解有功无功不平衡量
			self.planInit();	% 计算各节点功率的计划值

			self.opf.y = [1e-10 .* ones(self.n_bus, 1); -1e-10 .* ones(self.n_bus, 1)];

			self.opf.z = ones(self.opf.n_ieq, 1);
			self.opf.w = -0.5 .* self.opf.z;
			self.opf.l = self.opf.z;
			self.opf.u = self.opf.z;

			self.opf.g_min = [self.gen.Pmin; self.gen.Qmin; self.bus.Vmin; -self.branch.Pmax; -self.branch.Pmax];
			self.opf.g_max = [self.gen.Pmax; self.gen.Qmax; self.bus.Vmax; self.branch.Pmax; self.branch.Pmax];

			%%
			self.opf.epsilon = 1e-6;
			self.opf.n_iters_max = 50;
			self.opf.sigma = 0.01;	% 书本给出的经验值

			%% 
			if nargin ~= 2 || ~isstruct(solver)
				return;
			end
			if isfield(solver, 'epsilon')
				self.opf.epsilon = solver.epsilon;
			end
			if isfield(solver, 'n_iters_max')
				self.opf.n_iters_max = solver.n_iters_max;
			end
			if isfield(solver, 'sigma')
				self.opf.sigma = solver.sigma;
			end

			self.itlog.gap = [];
			self.itlog.ang = [self.bus.ang];
			self.itlog.mag = [self.bus.mag];
			self.itlog.Pg = [self.gen.Pg];
			self.itlog.Qg = [self.gen.Qg];
			self.itlog.Pij = [self.branch.Pij];
			self.itlog.Pji = [self.branch.Pji];

			self.itlog.Lx = [];
			self.itlog.Ly = [];
			self.itlog.Lz = [];
			self.itlog.Lw = [];
		end

		%% calcGap: 计算互补间隙
		function calcGap(self)
			self.opf.gap = self.opf.l' * self.opf.z - self.opf.u' * self.opf.w;
			self.itlog.gap = [self.itlog.gap, self.opf.gap];
		end

		%% opfConverged: 判断最优潮流求解过程是否收敛
		function [status] = opfConverged(self)

			status = 0;
			if abs(self.opf.gap) < self.opf.epsilon
				status = 1;
				return;
			end
			if self.opf.n_iters >= self.opf.n_iters_max
				status = 101;
				return;
			end
		end

		%% calcOpfCorrection: 计算修正方程
		function calcOpfCorrection(self)

			% 更新功率
			self.planUpdate();
			self.setPowerOutflow();
			self.setPowerUnbalance();
			self.setBranchPower(true);

			self.opf.L_inv_Z = diag(self.opf.z ./ self.opf.l);
			self.opf.U_inv_W = diag(self.opf.w ./ self.opf.u);
			self.calcYed();

			self.calcNableF();
			self.calcNableG();
			self.calcNableH();

			% 计算右侧结果向量
			self.opf.valueMatrix = [];

			% Ly
			self.opf.h = [self.bus.dP; self.bus.dQ];
			self.opf.Ly = self.opf.h;

			% Lz, Lw
			% self.opf.g = [self.gen.Pg; self.gen.Qg; self.bus.mag; self.branch.Pij];
			self.opf.g = [self.gen.Pg; self.gen.Qg; self.bus.mag; self.branch.Pij; self.branch.Pji];
			self.opf.Lz = self.opf.g - self.opf.l - self.opf.g_min;
			self.opf.Lw = self.opf.g + self.opf.u - self.opf.g_max;

			% Ll_mu, Lu_mu
			self.opf.Ll_mu = self.opf.l .* self.opf.z - self.opf.mu;
			self.opf.Lu_mu = self.opf.u .* self.opf.w + self.opf.mu;

			% Lx, Lx_d
			self.opf.Lx = self.opf.nable_F - self.opf.nable_H * self.opf.y - self.opf.nable_G * (self.opf.z + self.opf.w);
			self.opf.Lx_d = self.opf.Lx + self.opf.nable_G * ((self.opf.Ll_mu + self.opf.z .* self.opf.Lz) ./ self.opf.l + (self.opf.Lu_mu - self.opf.w .* self.opf.Lw) ./ self.opf.u);

			self.opf.valueMatrix = [
				-self.opf.Ll_mu ./ self.opf.l;
				self.opf.Lz;
				-self.opf.Lu_mu ./ self.opf.u;
				-self.opf.Lw;
				self.opf.Lx_d;
				-self.opf.Ly;
			];

			% 计算系数矩阵
			self.opf.parameterMatrix = [];

			self.opf.H = (-self.opf.nable2_F + self.opf.nable2_H_y + self.opf.nable2_G_c);
			self.opf.H_d = self.opf.H - self.opf.nable_G * diag(self.opf.z ./ self.opf.l - self.opf.w ./ self.opf.u) * self.opf.nable_G';

			self.opf.parameterMatrix = [
				eye(self.opf.n_ieq), self.opf.L_inv_Z, zeros(self.opf.n_ieq), zeros(self.opf.n_ieq), zeros(self.opf.n_ieq, self.opf.n_state), zeros(self.opf.n_ieq, self.opf.n_eq);
				zeros(self.opf.n_ieq), eye(self.opf.n_ieq), zeros(self.opf.n_ieq), zeros(self.opf.n_ieq), -self.opf.nable_G', zeros(self.opf.n_ieq, self.opf.n_eq);
				zeros(self.opf.n_ieq), zeros(self.opf.n_ieq), eye(self.opf.n_ieq), self.opf.U_inv_W, zeros(self.opf.n_ieq, self.opf.n_state), zeros(self.opf.n_ieq, self.opf.n_eq);
				zeros(self.opf.n_ieq), zeros(self.opf.n_ieq), zeros(self.opf.n_ieq), eye(self.opf.n_ieq), self.opf.nable_G', zeros(self.opf.n_ieq, self.opf.n_eq);
				zeros(self.opf.n_state, self.opf.n_ieq), zeros(self.opf.n_state, self.opf.n_ieq), zeros(self.opf.n_state, self.opf.n_ieq), zeros(self.opf.n_state, self.opf.n_ieq), self.opf.H_d, self.opf.nable_H;
				zeros(self.opf.n_eq, self.opf.n_ieq), zeros(self.opf.n_eq, self.opf.n_ieq), zeros(self.opf.n_eq, self.opf.n_ieq), zeros(self.opf.n_eq, self.opf.n_ieq), self.opf.nable_H', zeros(self.opf.n_eq);
			];

			self.opf.correction = self.opf.parameterMatrix \ self.opf.valueMatrix;

		end

		%% calcYed: 计算一个矩阵, Y^exp(j*theta_ij)
		function calcYed(self)
			% [theta_i, theta_j] = meshgrid(self.bus.ang);
			self.opf.Yed = conj(self.NAM.get()) .* exp(1i .* (self.bus.ang - self.bus.ang'));
			self.opf.Yed_Re = real(self.opf.Yed);
			self.opf.Yed_Im = imag(self.opf.Yed);
		end

		%% calcNableF: 计算目标函数的梯度
		function calcNableF(self)

			self.opf.f = 0;
			self.opf.nable_F = zeros(self.opf.n_state, 1);
			self.opf.nable2_F = zeros(self.opf.n_state);
			index = sub2ind([self.opf.n_state, self.opf.n_state], 1:self.n_gen, 1:self.n_gen);
			
			for k = 1:self.n_gen

				param_d0 = self.gencost.costParams(k, :);
				param_d1 = polyder(param_d0);
				param_d2 = polyder(param_d1);
				
				self.opf.f = self.opf.f + polyval(param_d0, self.gen.Pg(k));
				self.opf.nable_F(k) = polyval(param_d1, self.gen.Pg(k));
				self.opf.nable2_F(index(k)) = polyval(param_d2, self.gen.Pg(k));
			end
		end

		%% calcNableG: 计算不等式约束对应的方程组左侧的梯度
		function calcNableG(self)
			
			n_bus = self.n_bus;
			n_gen = self.n_gen;
			n_branch = self.n_branch;
			%% 计算一阶梯度
			% partial1 = zeros(n_bus, n_branch);
			% partial2 = zeros(n_bus, n_branch);
			partial1 = zeros(n_bus, 2.*n_branch);
			partial2 = zeros(n_bus, 2.*n_branch);

			index_i = getIndex(self.bus.id, self.branch.fid);
			index_j = getIndex(self.bus.id, self.branch.tid);

			Sij_dummy = self.bus.mag .* self.opf.Yed .* self.bus.mag';
			Pij_dummy = real(Sij_dummy);
			Qij_dummy = imag(Sij_dummy);

			% 遍历线路, 后期有时间再向量化处理, 这里暂时先用 for 循环
			for k = 1:n_branch

				f = index_i(k);
				t = index_j(k);

				% from
				% partial1(index_i(k), k) = partial1(index_i(k), k) - Qij_dummy(index_i(k), index_j(k));
				% partial1(index_j(k), k) = partial1(index_j(k), k) + Qij_dummy(index_i(k), index_j(k));
				% partial2(index_i(k), k) = partial2(index_i(k), k) + Pij_dummy(index_i(k), index_j(k)) ./ self.bus.mag(index_i(k)) - 2 .* self.bus.mag(index_i(k)) .* real(self.NAM.value(index_i(k), index_j(k)));
				% partial2(index_j(k), k) = partial2(index_j(k), k) + Pij_dummy(index_i(k), index_j(k)) ./ self.bus.mag(index_j(k));
				% from
				partial1(f, k) = partial1(f, k) - Qij_dummy(f, t);
				partial1(t, k) = partial1(t, k) + Qij_dummy(f, t);
				partial2(f, k) = partial2(f, k) + Pij_dummy(f, t) ./ self.bus.mag(f) - 2 .* self.bus.mag(f) .* real(self.NAM.value(f, t));
				partial2(t, k) = partial2(t, k) + Pij_dummy(f, t) ./ self.bus.mag(t);
				% to
				partial1(t, k+n_branch) = partial1(t, k+n_branch) - Qij_dummy(t, f);
				partial1(f, k+n_branch) = partial1(f, k+n_branch) + Qij_dummy(t, f);
				partial2(t, k+n_branch) = partial2(t, k+n_branch) + Pij_dummy(t, f) ./ self.bus.mag(t) - 2 .* self.bus.mag(t) .* real(self.NAM.value(t, f));
				partial2(f, k+n_branch) = partial2(f, k+n_branch) + Pij_dummy(t, f) ./ self.bus.mag(f);
			end

			% 拼接
			self.opf.nable_G = [
				% eye(n_gen), zeros(n_gen), zeros(n_gen, n_bus), zeros(n_gen, n_branch);
				eye(n_gen), zeros(n_gen), zeros(n_gen, n_bus), zeros(n_gen, 2.*n_branch);
				% zeros(n_gen), eye(n_gen), zeros(n_gen, n_bus), zeros(n_gen, n_branch);
				zeros(n_gen), eye(n_gen), zeros(n_gen, n_bus), zeros(n_gen, 2.*n_branch);
				zeros(n_bus, n_gen), zeros(n_bus, n_gen), zeros(n_bus), partial1;
				zeros(n_bus, n_gen), zeros(n_bus, n_gen), eye(n_bus), partial2;
			];

			%% 计算二阶梯度
			part1 = zeros(n_bus);
			part2 = zeros(n_bus);
			part4 = zeros(n_bus);

			% c_branch = self.opf.z(end-n_branch+1:end) + self.opf.w(end-n_branch+1:end);
			c1_branch = self.opf.z(end-2.*n_branch+1:end-n_branch) + self.opf.w(end-2.*n_branch+1:end-n_branch);
			c2_branch = self.opf.z(end-n_branch+1:end) + self.opf.w(end-n_branch+1:end);

			% from = getIndex(self.bus.id, self.branch.fid);
			% to = getIndex(self.bus.id, self.branch.tid);

			% 这里写法非常不优雅, 时间紧迫, 以后有时间再改
			for k = 1:n_branch

				f = index_i(k);
				t = index_j(k);
				% c = c_branch(k);
				c1 = c1_branch(k);
				c2 = c2_branch(k);

				% part1
				% part1(f, f) = part1(f, f) - Pij_dummy(f, t) .* c;
				% part1(f, t) = part1(f, t) + Pij_dummy(f, t) .* c;
				% part1(t, f) = part1(t, f) + Pij_dummy(f, t) .* c;
				% part1(t, t) = part1(t, t) - Pij_dummy(f, t) .* c;	% debug
				part1(f, f) = part1(f, f) - Pij_dummy(f, t) .* c1 - Pij_dummy(t, f) .* c2;
				part1(f, t) = part1(f, t) + Pij_dummy(f, t) .* c1 + Pij_dummy(t, f) .* c2;
				part1(t, f) = part1(t, f) + Pij_dummy(f, t) .* c1 + Pij_dummy(t, f) .* c2;
				part1(t, t) = part1(t, t) - Pij_dummy(f, t) .* c1 - Pij_dummy(t, f) .* c2;	% debug

				% part2
				% part2(f, f) = part2(f, f) - Qij_dummy(f, t) ./ self.bus.mag(f) .* c;
				% part2(f, t) = part2(f, t) - Qij_dummy(f, t) ./ self.bus.mag(t) .* c;
				% part2(t, f) = part2(t, f) + Qij_dummy(f, t) ./ self.bus.mag(f) .* c;
				% part2(t, t) = part2(t, t) + Qij_dummy(f, t) ./ self.bus.mag(t) .* c;
				part2(f, f) = part2(f, f) - Qij_dummy(f, t) ./ self.bus.mag(f) .* c1 + Qij_dummy(t, f) ./ self.bus.mag(f) .* c2;
				part2(f, t) = part2(f, t) - Qij_dummy(f, t) ./ self.bus.mag(t) .* c1 + Qij_dummy(t, f) ./ self.bus.mag(t) .* c2;
				part2(t, f) = part2(t, f) + Qij_dummy(f, t) ./ self.bus.mag(f) .* c1 - Qij_dummy(t, f) ./ self.bus.mag(f) .* c2;
				part2(t, t) = part2(t, t) + Qij_dummy(f, t) ./ self.bus.mag(t) .* c1 - Qij_dummy(t, f) ./ self.bus.mag(t) .* c2;

				% part4
				% part4(f, f) = part4(f, f) - 2.* real(self.NAM.value(f, t)) .* c;
				% part4(f, t) = part4(f, t) + Pij_dummy(f, t) ./ (self.bus.mag(f) .* self.bus.mag(t)) .* c;
				% part4(t, f) = part4(t, f) + Pij_dummy(f, t) ./ (self.bus.mag(f) .* self.bus.mag(t)) .* c;
				% part4(t, t) = part4(t, t) + 0 .* c;
				part4(f, f) = part4(f, f) - 2.* real(self.NAM.value(f, t)) .* c1 + 0 .* c2;
				part4(f, t) = part4(f, t) + (Pij_dummy(f, t) .* c1 + Pij_dummy(t, f) .* c2) ./ (self.bus.mag(f) .* self.bus.mag(t));
				part4(t, f) = part4(t, f) + (Pij_dummy(f, t) .* c1 + Pij_dummy(t, f) .* c2) ./ (self.bus.mag(f) .* self.bus.mag(t));
				part4(t, t) = part4(t, t) + 0 .* c1 - 2.* real(self.NAM.value(t, f)) .* c2;
			end

			self.opf.nable2_G_c = [
				zeros(n_gen), zeros(n_gen), zeros(n_gen, n_bus), zeros(n_gen, n_bus);
				zeros(n_gen), zeros(n_gen), zeros(n_gen, n_bus), zeros(n_gen, n_bus);
				zeros(n_bus, n_gen), zeros(n_bus, n_gen), part1, part2;
				zeros(n_bus, n_gen), zeros(n_bus, n_gen), part2', part4;
			];
		end

		%% calcNableH: 计算等式约束对应的方程组左侧的梯度
		function calcNableH(self)

			n_bus = self.n_bus;
			n_gen = self.n_gen;

			partial1 = zeros(n_gen, n_bus);

			index = sub2ind([n_gen, n_bus], 1:n_gen, getIndex(self.bus.id, self.gen.nid));
			partial1(index) = 1;

			self.opf.nable_H = [
				partial1, zeros(n_gen, n_bus);
				zeros(n_gen, n_bus), partial1;
				-self.getOpfJacobian();
			];

			%% 计算二阶梯度
			yp = self.opf.y(1:end/2);
			yq = self.opf.y(end/2+1:end);
			Y_diag_re = real(diag(self.NAM.value));
			Y_diag_im = imag(diag(self.NAM.value));
			U = self.bus.mag;

			Uyp = U .* yp;
			Uyq = U .* yq;
			Uyp_diag = diag(Uyp);
			Uyq_diag = diag(Uyq);

			% 左上方矩阵
			part1 = -Uyp .* (self.opf.Yed_Re * U) - (Uyp' * self.opf.Yed_Re)' .* U - Uyq .* (self.opf.Yed_Im * U) - (Uyq' * self.opf.Yed_Im)' .* U + 2 .* U .^ 2 .* (yp .* Y_diag_re + yq .* Y_diag_im);
			% part1 = -Uyp .* (self.opf.Yed_Re * U) - (Uyp' * self.opf.Yed_Re)' .* U - Uyq .* (self.opf.Yed_Im * U) - (Uyq' * self.opf.Yed_Im)' .* U;
			part3 = (Uyp .* self.opf.Yed_Re + Uyq .* self.opf.Yed_Im) .* U';
			part3 = part3 + part3';

			% 右上方, 左下方
			part2 = -yp .* self.opf.Yed_Im * U + self.opf.Yed_Im' * Uyp + yq .* self.opf.Yed_Re * U - self.opf.Yed_Re' * Uyq;
			part4 = -Uyp .* self.opf.Yed_Im + Uyq .* self.opf.Yed_Re;
			part4 = part4 + part4';

			% 右下方
			part5 = 2 .* (yp .* diag(self.opf.Yed_Re) + yq .* diag(self.opf.Yed_Im));
			part7 = yp .* self.opf.Yed_Re + yq .* self.opf.Yed_Im;
			part7 = part7 + part7';

			% 拼接
			diag_index = sub2ind(size(part3), (1:n_bus)', (1:n_bus)');
			part3(diag_index) = part1;
			part4(diag_index) = part2;
			part7(diag_index) = part5;

			self.opf.nable2_H_y = [
				zeros(n_gen), zeros(n_gen), zeros(n_gen, n_bus), zeros(n_gen, n_bus);
				zeros(n_gen), zeros(n_gen), zeros(n_gen, n_bus), zeros(n_gen, n_bus);
				zeros(n_bus, n_gen), zeros(n_bus, n_gen), -part3, -part4;
				zeros(n_bus, n_gen), zeros(n_bus, n_gen), -part4', -part7;
			];
		end

		%% calcOpfStep: 计算步长
		function calcOpfStep(self)

			self.opf.delta_l = self.opf.correction(self.opf.n_ieq+1:2.*self.opf.n_ieq);
			self.opf.delta_u = self.opf.correction(3.*self.opf.n_ieq+1:4.*self.opf.n_ieq);
			self.opf.delta_z = self.opf.correction(1:self.opf.n_ieq);
			self.opf.delta_w = self.opf.correction(2.*self.opf.n_ieq+1:3.*self.opf.n_ieq);

			index_l = find(self.opf.delta_l < 0);
			index_u = find(self.opf.delta_u < 0);
			index_z = find(self.opf.delta_z < 0);
			index_w = find(self.opf.delta_w > 0);

			self.opf.alpha_p = 0.9995 .* min([-self.opf.l(index_l) ./ self.opf.delta_l(index_l); -self.opf.u(index_u) ./ self.opf.delta_u(index_u); 1]);
			self.opf.alpha_d = 0.9995 .* min([-self.opf.z(index_z) ./ self.opf.delta_z(index_z); -self.opf.w(index_w) ./ self.opf.delta_w(index_w); 1]);
		end

		%% setNewOpfParameters: 修正
		function setNewOpfParameters(self)

			self.opf.delta_x = self.opf.correction(end-self.opf.n_eq-self.opf.n_state+1:end-self.opf.n_eq);
			self.opf.delta_y = self.opf.correction(4.*self.opf.n_ieq+1:4.*self.opf.n_ieq+self.opf.n_eq);

			self.opf.x = self.opf.x + self.opf.alpha_p .* self.opf.delta_x;
			self.opf.y = self.opf.y + self.opf.alpha_d .* self.opf.delta_y;
			self.opf.z = self.opf.z + self.opf.alpha_d .* self.opf.delta_z;
			self.opf.w = self.opf.w + self.opf.alpha_d .* self.opf.delta_w;
			self.opf.l = self.opf.l + self.opf.alpha_p .* self.opf.delta_l;
			self.opf.u = self.opf.u + self.opf.alpha_p .* self.opf.delta_u;

			assert(all(self.opf.l > 0));
			assert(all(self.opf.u > 0));

			self.gen.Pg = self.opf.x(1:self.n_gen);
			self.gen.Qg = self.opf.x(self.n_gen+1:2*self.n_gen);
			self.bus.ang = self.opf.x(end-2*self.n_bus+1:end-self.n_bus);
			self.bus.mag = self.opf.x(end-self.n_bus+1:end);

			self.bus.ang = self.bus.ang - self.bus.ang(self.slack);

			self.planInit();

			self.itlog.ang = [self.itlog.ang, self.bus.ang];
			self.itlog.mag = [self.itlog.mag, self.bus.mag];
			self.itlog.Pg = [self.itlog.Pg, self.gen.Pg];
			self.itlog.Qg = [self.itlog.Qg, self.gen.Qg];
			self.itlog.Pij = [self.itlog.Pij, self.branch.Pij];
			self.itlog.Pji = [self.itlog.Pji, self.branch.Pji];

			self.itlog.Lx = [self.itlog.Lx, self.opf.Lx];
			self.itlog.Ly = [self.itlog.Ly, self.opf.Ly];
			self.itlog.Lz = [self.itlog.Lz, self.opf.Lz];
			self.itlog.Lw = [self.itlog.Lw, self.opf.Lw];

		end



	end
end

