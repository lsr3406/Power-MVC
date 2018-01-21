%% SteadyState
classdef SteadyState < handle
	properties
		% 下面这三个表示电网中一些元素的状态, 它们可能经常改变.
		nodes;		% 节点状态, [id, type, mag, ang, Pg, Qg, Pc, Qc, Pd, Qd, Pis, Qis, Pout, Qout, dP, dQ, Pmax, Qmax, Pmin, Qmin, Vmax, Vmin]
		generator;	% 发电机状态, [id, nid, Pg, Qg, status]
		branches;	% 线路状态,	[id, fid, tid, ratio, Pij, Qij, Pji, Qji, dP, dQ, status]

		% 下面存放的是一些中间计算结果, 由于这些矩阵较大, 调用较多, 但不经常改变, 因此当做属性存储
		NAM;		% 节点导纳矩阵
		FDM1;		% 快速分解法修正方程第一系数矩阵
		FDM2;		% 快速分解法修正方程第二系数矩阵
		% 雅可比矩阵在求解过程中一直改变, 不直接当做属性存储

		itlog;
	end
	properties (Dependent)
		pqIndex;
		pqCount;
		pvIndex;
		pvCount;
		refIndex;
		refCount;
		loadIndex;
		loadCount;
	end
	methods

		%% get.pqIndex: 获取 pq 节点的索引
		function [pqIndex] = get.pqIndex(self)
			pqIndex = find(self.nodes.type == 1);
		end
		%% get.pqCount: 获取 pq 节点的数量
		function [pqCount] = get.pqCount(self)
			pqCount = sum(self.nodes.type == 1);
		end
		%% get.pvIndex: 获取 pv 节点的索引
		function [pvIndex] = get.pvIndex(self)
			pvIndex = find(self.nodes.type == 2);
		end
		%% get.pvCount: 获取 pv 节点的数量
		function [pvCount] = get.pvCount(self)
			pvCount = sum(self.nodes.type == 2);
		end
		%% get.refIndex: 获取 ref 节点的索引
		function [refIndex] = get.refIndex(self)
			refIndex = find(self.nodes.type == 3);
		end
		%% get.refCount: 获取 ref 节点的数量
		function [refCount] = get.refCount(self)
			refCount = sum(self.nodes.type == 3);
		end
		%% get.loadIndex: 获取负荷节点的索引
		function [loadIndex] = get.loadIndex(self)
			loadIndex = find(self.nodes.Pd > 0);
		end
		%% get.loadCount: 获取负荷节点的数量
		function [loadCount] = get.loadCount(self)
			loadCount = sum(self.nodes.Pd > 0);
		end

		%% SteadyState: 电力系统稳态模型构造器
		% function [self] = SteadyState(powerSystemConfig)
		% end

		%% init: 电力系统稳态模型初始化
		function init(self, mpc, nodes, gens, branches)
			nodesLength = length(nodes.id);
			gensLength = length(gens.id);
			branchesLength = length(branches.id);

			self.nodes.id = mpc.bus(:, 1);
			self.nodes.type = mpc.bus(:, 2);
			self.nodes.mag = ones(nodesLength, 1);
			self.nodes.ang = zeros(nodesLength, 1);
			self.nodes.Pg = zeros(nodesLength,1);
			self.nodes.Qg = zeros(nodesLength,1);
			self.nodes.Pd = mpc.bus(:, 3)./mpc.baseMVA;
			self.nodes.Qd = mpc.bus(:, 4)./mpc.baseMVA;
			self.nodes.Pc = zeros(nodesLength,1);
			self.nodes.Qc = zeros(nodesLength,1);
			self.nodes.Pis = zeros(nodesLength,1);
			self.nodes.Qis = zeros(nodesLength,1);
			self.nodes.Pout = zeros(nodesLength,1);
			self.nodes.Qout = zeros(nodesLength,1);
			self.nodes.dP = zeros(nodesLength,1);
			self.nodes.dQ = zeros(nodesLength,1);
			self.nodes.Qmax = zeros(nodesLength,1);
			self.nodes.Qmin = zeros(nodesLength,1);
			self.nodes.Vmax = mpc.bus(:, 12);
			self.nodes.Vmin = mpc.bus(:, 13);

			self.generator.nId = mpc.gen(:, 1);
			self.generator.id = (1:length(self.generator.nId))';
			self.generator.Pg = mpc.gen(:, 2)./mpc.gen(:, 7);	% Pg./mBase
			self.generator.Qg = mpc.gen(:, 3)./mpc.gen(:, 7);	% Qg./mBase
			self.generator.status = zeros(gensLength, 1);

			self.branches.fid = mpc.branch(:, 1);
			self.branches.tid = mpc.branch(:, 2);
			self.branches.id = (1:length(self.branches.fid))';
			self.branches.ratio = mpc.branch(:, 9);
			self.branches.Pij = zeros(branchesLength, 1);
			self.branches.Qij = zeros(branchesLength, 1);
			self.branches.Pji = zeros(branchesLength, 1);
			self.branches.Qji = zeros(branchesLength, 1);
			self.branches.dP = zeros(branchesLength, 1);
			self.branches.dQ = zeros(branchesLength, 1);
			self.branches.status = ones(branchesLength, 1);

			self.itlog.mag = [];
			self.itlog.ang = [];
			self.itlog.dP = [];
			self.itlog.dQ = [];
		end

		%% NAMInit: 计算节点导纳矩阵
		function NAMInit(self, nodes, branches)
			lineData = branches.getLineData();
			nodeData = nodes.getNodeData();
			getAdmittanceMatrix(nodeData, lineData, self);
		end
		
		%% votageInit: 获取迭代初始值(考虑对节点的设置及发电机的设置)
		function votageInit(self,nodes,config)
			% if 0
			% 	for k = generator.id
			% 		index = find(self.id == generator.nid(k));
			% 		self.mag0(index) = generator.votage(k);
			% 	end
			% end
			if isstruct(config) && isfield(config, 'start') && strcmp(config.start, 'default')
				self.nodes.mag = nodes.mag0;
				self.nodes.ang = nodes.ang0;			
			else	% 0-1 启动
				self.nodes.mag = ones(length(self.nodes.id),1);
				self.nodes.ang = zeros(length(self.nodes.id),1);
			end
			% self.nodes.iterationData_mag = [];
			% self.nodes.iterationData_ang = [];
			% self.nodes.iterationData_dq = [];
			% self.nodes.iterationData_dp = [];
		end

		%% planInit: 初始化,用于计算各节点功率计划值,并根据发电机的情况计算出该节点可发出的的最大无功功率
		function planInit(self,nodes,generator)
			for k = generator.id
				index = find(nodes.id == generator.nid(k));
				% 这里对各节点的功率计划的赋值有一大部分(PQ节点)是没有意义的
				self.nodes.Pg(index) = self.nodes.Pg(index) + generator.Pg(k);
				self.nodes.Qg(index) = self.nodes.Qg(index) + generator.Qg(k);
				self.nodes.Qmin(index) = self.nodes.Qmin(index) + generator.Qmin(k);
				self.nodes.Qmax(index) = self.nodes.Qmax(index) + generator.Qmax(k);
			end
			% self.nodes.Pd = nodes.Pd;
			% self.nodes.Qd = nodes.Qd;
			% self.nodes.Pis = self.nodes.Pg - self.nodes.Pd;
			% self.nodes.Qis = self.nodes.Qg - self.nodes.Qd;
		end

		%% planUpdate: 更新各节点功率计划值
		function planUpdate(self)
			self.nodes.Pis = self.nodes.Pg - self.nodes.Pd;
			self.nodes.Qis = self.nodes.Qg - self.nodes.Qd;
		end

		%% getPowerOutflow: 计算从 nid 节点注入电网的潮流
		function getPowerOutflow(self,nid)
			for k = 1:length(nid)
				index = find(self.nodes.id == nid(k));
				% index
				% self
				Sf = conj(self.NAM(index,:))*(self.nodes.mag.*exp(i.*(self.nodes.ang(index)-self.nodes.ang))).*self.nodes.mag(index);
				self.nodes.Pout(k) = real(Sf);
				self.nodes.Qout(k) = imag(Sf);
			end
		end

		%% getPowerUnbalance: 计算功率不平衡量,仅用于牛拉法和 PQ 分解法计算,返回变量均涵盖所有节点并已分类
		function getPowerUnbalance(self)
			self.nodes.dP = self.nodes.Pis - self.nodes.Pout;
			self.nodes.dQ = self.nodes.Qis - self.nodes.Qout;
		end

		%% getIterationData: 每次使用潮流方程之后调用此方法获取迭代信息
		% TODO 优化内存使用方式
		function getIterationData(self)
			self.itlog.mag = [self.itlog.mag, self.nodes.mag];
			self.itlog.ang = [self.itlog.ang, self.nodes.ang];
			self.itlog.dP = [self.itlog.dP, self.nodes.dP];
			self.itlog.dQ = [self.itlog.dQ, self.nodes.dQ];
		end

		%% checkConverged: 判断是否收敛, 同时检查是否结束循环
		function [status] = checkConverged(self, config, it)
			% 收敛判据
			if config.epsilon >= max(abs([self.nodes.dP([self.pqIndex;self.pvIndex]);self.nodes.dQ(self.pqIndex)]))
				status = 1;	% 收敛
			elseif it >= config.maxIteration
				status = 101;	% 迭代次数超出上限
			else
				status = 0;	% 未收敛, 继续计算
			end
		end

		%% getJacobian: 在当前网络状态下获取雅可比矩阵,用于牛顿法的计算
		function [jacobianMatrix] = getJacobian(self)

			% 为方便构造雅可比矩阵,现分别对 PQ 节点及 PV 节点提取出节点电压及节点导纳, 并计算节点之间的相角矩阵
			% 本方法中, 变量名以 1 结尾的均与雅可比矩阵中的 H N M 有关, 变量名以 2 结尾的均与雅可比矩阵中的 N M L 有关
			NAM1 = self.NAM([self.pqIndex; self.pvIndex], [self.pqIndex; self.pvIndex]);
			NAM2 = self.NAM([self.pqIndex], [self.pqIndex]);

			mag1 = self.nodes.mag([self.pqIndex; self.pvIndex]);
			ang1 = self.nodes.ang([self.pqIndex; self.pvIndex]);
			mag2 = self.nodes.mag([self.pqIndex]);
			ang2 = self.nodes.ang([self.pqIndex]);

			nodeAngle = ang1*ones(1,length(ang1));
			nodeAngle = nodeAngle - nodeAngle';

			% 开始计算雅可比矩阵,创建初始状态
			jacPrimary = mag1*(mag1').*conj(NAM1).*exp(i.*nodeAngle);
			% 这个结果的对角元无实际意义,下面去掉对角元
			jacPrimary = jacPrimary - diag(diag(jacPrimary));

			% 先计算各个分块矩阵的非对角元
			% 需要注意的是, 下面的两个索引表示在雅可比矩阵的初始状态中 PQPV 节点及 PQ 节点所在的位置, 而 self.pqIndex 表示在节点对象中 PQ 节点所在的位置
			pqpvIndex = 1:(self.pqCount + self.pvCount);
			pqIndex = 1:self.pqCount;

			% 雅可比矩阵的非对角元均可以在初始状态矩阵中提取

			Hij = imag(jacPrimary([pqpvIndex],[pqpvIndex]));	% Hij 表示系统中有功潮流差与相角差的数量关系
			Nij = real(jacPrimary([pqpvIndex],[pqIndex]));	% Nij 表示系统中有功潮流差与电压差的相对值的数量关系
			Mij = -real(jacPrimary([pqIndex],[pqpvIndex]));	% Mij 表示系统中无功潮流差与相角差的数量关系
			Lij = imag(jacPrimary([pqIndex],[pqIndex]));	% Lij 表示系统中无功潮流差与电压差的相对值的数量关系

			% 再计算各个分块矩阵的对角元. 求导的原因, 使得对角元需要单独计算
			Hii = -self.nodes.Qout([self.pqIndex; self.pvIndex]) - mag1.^2.*imag(diag(NAM1));
			Nii = self.nodes.Pout(self.pqIndex) + mag2.^2.*real(diag(NAM2));
			Jii = self.nodes.Pout(self.pqIndex) - mag2.^2.*real(diag(NAM2));
			Lii = self.nodes.Qout(self.pqIndex) - mag2.^2.*imag(diag(NAM2));

			Hij = Hij + diag(Hii);
			Nij(pqIndex,pqIndex) = Nij(pqIndex,pqIndex) + diag(Nii);
			Mij(pqIndex,pqIndex) = Mij(pqIndex,pqIndex) + diag(Jii);
			Lij = Lij + diag(Lii);

			% 在这里拼成雅可比矩阵
			jacobianMatrix = [Hij, Nij; Mij, Lij;];
		end

		%% getNewNodeParameter: 修正方程, 用于牛顿法
		function [dTheta,dUpn] = getVotageCorrection(self,jacobianMatrix)

			dP = self.nodes.dP([self.pqIndex; self.pvIndex]);
			dQ = self.nodes.dQ([self.pqIndex]);

			dPara = jacobianMatrix \ [dP; dQ];

			dTheta = dPara(1:(self.pqCount + self.pvCount));	% PQ 及 PV 节点相角的修正量
			dUpn = dPara((self.pqCount + self.pvCount + 1):end);	% PQ 节点电压相对值的修正量
			
			% 将电压和相角误差规范化, 即填充至大向量中, 方便迭代方程直接使用
			dUpn_t = zeros(self.pqCount + self.pvCount + self.refCount, 1);	% 初始化相角修正量
			dUpn_t([self.pqIndex]) = dUpn;	% 代入 PQ 及 PV 节点相角的修正量
			dUpn = dUpn_t;
			dTheta_t = zeros(self.pqCount + self.pvCount + self.refCount, 1);
			dTheta_t([self.pqIndex; self.pvIndex]) = dTheta;
			dTheta = dTheta_t;
		end

		%% getNewVotage: 获取新的节点电压及相角, 用于牛顿法
		function getNewVotage(self, dTheta, dUpn)
			self.nodes.ang = self.nodes.ang + (dTheta);
			self.nodes.mag = self.nodes.mag.*(1+dUpn);
		end

		%% getNodePowerGeneration: 计算各节点的实际发电量 (仅计算了 pv 节点的无功和平衡节点的有功无功)
		function getNodePowerGeneration(self)
			self.nodes.Pg(self.refIndex) = self.nodes.Pout(self.refIndex) + self.nodes.Pd(self.refIndex);
			self.nodes.Qg([self.pvIndex; self.refIndex]) = self.nodes.Pout([self.pvIndex; self.refIndex]) + self.nodes.Pd([self.pvIndex; self.refIndex]);
		end

		%% getLinePower: 计算普通线路的潮流
		function [Sij,Sji,dS] = getLinePower(self,branches,index)
			% matpower 默认在计算线路损耗时, 并没有考虑线路的充电电容
			conjYi0 = conj(branches.y(index)./2);
			% conjYi0 = 0;
			conjYj0 = conjYi0;
			conjYij = 1./conj(branches.z(index));
			fid = find(self.nodes.id == self.branches.fid(index));	% 得到线路始末端节点的id
			tid = find(self.nodes.id == self.branches.tid(index));	% 得到线路始末端节点的id
			Sij = (self.nodes.mag(fid)).^2.*(conjYi0 + conjYij) - self.nodes.mag(fid).*self.nodes.mag(tid).*conjYij.*exp((self.nodes.ang(fid)-self.nodes.ang(tid)).*i);
			Sji = (self.nodes.mag(tid)).^2.*(conjYj0 + conjYij) - self.nodes.mag(fid).*self.nodes.mag(tid).*conjYij.*exp((self.nodes.ang(tid)-self.nodes.ang(fid)).*i);
			dS = Sij + Sji;
		end

		%% getTransformerPower: 计算变压器的潮流
		function [Sij,Sji,dS] = getTransformerPower(self,branches,index)
			[trZ,trY1,trY2] = gamma2pi(branches.z(index), branches.y(index), branches.ratio(index));
			conjYi0 = conj(trY1) - branches.y(index)./2;
			conjYj0 = conj(trY2) - branches.y(index)./2;
			conjYij = conj(1./trZ);
			fid = find(self.nodes.id == self.branches.fid(index));	% 得到线路始末端节点的id
			tid = find(self.nodes.id == self.branches.tid(index));	% 得到线路始末端节点的id
			Sij = (self.nodes.mag(fid)).^2.*(conjYi0 + conjYij) - self.nodes.mag(fid).*self.nodes.mag(tid).*conjYij.*exp((self.nodes.ang(fid)-self.nodes.ang(tid)).*i);
			Sji = (self.nodes.mag(tid)).^2.*(conjYj0 + conjYij) - self.nodes.mag(fid).*self.nodes.mag(tid).*conjYij.*exp((self.nodes.ang(tid)-self.nodes.ang(fid)).*i);
			dS = Sij + Sji;
		end

		%% getBranchPower: 计算线路功率及功率损耗
		function getBranchPower(self, branches)
			Sij = zeros(length(self.branches.id), 1);
			Sji = zeros(length(self.branches.id), 1);
			dS = zeros(length(self.branches.id), 1);
			for k = self.branches.id'
				if (self.branches.ratio(k)==1) || (self.branches.ratio(k)==0)
					[Sij(k), Sji(k), dS(k)] = self.getLinePower(branches, k);
				else
					[Sij(k), Sji(k), dS(k)] = self.getTransformerPower(branches, k);
				end	
			end
			self.branches.Pij = real(Sij);
			self.branches.Qij = imag(Sij);
			self.branches.Pji = real(Sji);
			self.branches.Qji = imag(Sji);
			self.branches.dP = real(dS);
			self.branches.dQ = imag(dS);
		end

		%% getConpensatorPower: 计算无功补偿器电量
		function getConpensatorPower(self, nodes)
			self.nodes.Pc = self.nodes.mag.^2.*nodes.g;
			self.nodes.Qc = self.nodes.mag.^2.*nodes.b;
		end

		%% NRIteration: 牛顿法迭代程序核心部分
		function [result] = NRIteration(self, config)
			% 从这里开始进入循环,跳出循环的条件为迭代次数超过最大次数或潮流不平衡量的一范数小于给定值
			result.it = 0;
			while 1
				self.planUpdate();	% 确定各节点功率需求
				self.getPowerOutflow(self.nodes.id);	% 潮流方程
				self.getPowerUnbalance();	% 误差方程
				self.getIterationData();	% 存储迭代信息,包括电压,相角,有功不平衡量,无功不平衡量
				
				% 判收敛
				result.status = self.checkConverged(config, result.it);
				if(result.status)	% result.status 非 0 表示收敛或错误
					return;
				end

				J = self.getJacobian();	% ❀ 将各节点电压的初值代入修正方程,求解系数矩阵.
				[dTheta, dUpn] = self.getVotageCorrection(J);	% 修正方程
				self.getNewVotage(dTheta, dUpn);	% 迭代方程
				result.it = result.it + 1;
				% ❀ 运用各节点电压的新值进行下一步的迭代.
			end
		end

		%% solvePowerFlow: 电力系统潮流计算 (写好后将其放至 Model 层)
		function [result] = solvePowerFlow(self, nodes, generator, branches, config)

			% ❀ 生成节点导纳矩阵.
			self.NAMInit(nodes, branches);

			% ❀ 设置节点电压的初值(电压,相角).
			self.votageInit(nodes, config);

			% ❀ 将各节点电压的初值代入潮流方程,并求解有功无功不平衡量
			self.planInit(nodes, generator);	% 计算各节点功率的计划值

			switch config.method
				case 'NR'	% 牛顿法
					result = self.NRIteration(config);
				otherwise
					% TODO otherwise
			end
			if(result.status ~= 1)
				return;
			end
			self.getBranchPower(branches);
			self.getConpensatorPower(nodes);
		end
		
	end
end