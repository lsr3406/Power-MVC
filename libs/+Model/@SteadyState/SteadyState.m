%% SteadyState
classdef SteadyState < handle
	properties
		mBase;

		%% 下面这三个表示电网中一些元素的状态, 它们可能经常改变.
		nodes;		% 节点状态, [id, type, mag, ang, Pg, Qg, Pc, Qc, Pd, Qd, Pis, Qis, Pout, Qout, dP, dQ, Pmax, Qmax, Pmin, Qmin, Vmax, Vmin]
		generator;	% 发电机状态, [id, nid, Pg, Qg, status]
		branches;	% 线路状态,	[id, fid, tid, ratio, Pij, Qij, Pji, Qji, dP, dQ, status]

		%% 下面存放的是一些中间计算结果, 由于这些矩阵较大, 调用较多, 但不经常改变, 因此当做属性存储
		NAM;		% 节点导纳矩阵
		FDM1;		% 快速分解法修正方程第一系数矩阵
		FDM2;		% 快速分解法修正方程第二系数矩阵
		%% 雅可比矩阵在求解过程中一直改变, 不直接当做属性存储

		itlog;		% 和电力系统潮流求解相关的记录

		magStep;
		angStep;
	end
	properties (Dependent)
		pqIndex;		% pq 节点索引
		pqCount;		% pq 节点数量
		pvIndex;		% pv 节点索引
		pvCount;		% pv 节点数量
		refIndex;		% 平衡节点索引
		refCount;		% 平衡节点索引
		loadIndex;		% 负荷节点索引
		loadCount;		% 负荷节点数量

		pqGens;			% 处于 pq 节点的发电机索引
		pvGens;			% 处于 pv 节点的发电机索引
		refGens;		% 处于平衡节点的发电机索引

		lineIndex;			% 普通线路索引(id)
		lineCount;			% 普通线路数量
		transformerIndex;	% 变压器索引(id)
		transformerCount;	% 变压器数量
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

		%% get.pqGens: 小型发电机数量
		function [pqGens] = get.pqGens(self)
			pqGens = sum(self.generator.type == 1);
		end
		%% get.pvGens: 中型发电机数量
		function [pvGens] = get.pvGens(self)
			pvGens = sum(self.generator.type == 2);
		end
		%% get.refGens: 大型发电机数量
		function [refGens] = get.refGens(self)
			refGens = sum(self.generator.type == 3);
		end

		%% get.lineIndex: 线路索引
		function [lineIndex] = get.lineIndex(self)
			lineIndex = find(self.branches.ratio==0);
		end
		%% get.lineCount: 线路数量
		function [lineCount] = get.lineCount(self)
			lineCount = sum(self.branches.ratio==0);
		end
		%% get.transformerIndex: 变压器索引
		function [transformerIndex] = get.transformerIndex(self)
			transformerIndex = find(self.branches.ratio~=0);
		end
		%% get.transformerCount: 变压器数量
		function [transformerCount] = get.transformerCount(self)
			transformerCount = sum(self.branches.ratio~=0);
		end

		%% get.isLine: 返回该支路的类型, 是否为普通线路
		function [isLine] = isLine(self, id)
			% 支路的 id 也是 index
			isLine = self.branches.ratio(id) == 0;
		end

		%% get.isTransformer: 返回该支路的类型, 是否为普通线路
		function [isTransformer] = isTransformer(self, id)
			% 支路的 id 也是 index
			isTransformer = self.branches.ratio(id) ~= 0;
		end

		%% init: 电力系统稳态模型初始化
		function init(self, mpc)

			mpc.bus = sortrows(mpc.bus);
			self.mBase = mpc.baseMVA;
			self.nodes.id = mpc.bus(:, 1);
			nodesLength = length(self.nodes.id);
			self.nodes.type = mpc.bus(:, 2);
			self.nodes.g = mpc.bus(:,5)./self.mBase;
			self.nodes.b = mpc.bus(:,6)./self.mBase;
			self.nodes.mag0 = mpc.bus(:,8);
			self.nodes.ang0 = mpc.bus(:,9).*pi./180;
			self.nodes.mag = self.nodes.mag0;
			self.nodes.ang = self.nodes.ang0;
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

			self.generator.nid = mpc.gen(:, 1);
			gensLength = length(self.generator.nid);
			self.generator.id = (1:gensLength)';
			self.generator.Pg = mpc.gen(:, 2)./mpc.gen(:, 7);	% Pg./mBase
			self.generator.Qg = mpc.gen(:, 3)./mpc.gen(:, 7);	% Qg./mBase
			self.generator.mBase = mpc.gen(:,7);
			self.generator.Qmax = mpc.gen(:,4)./self.mBase;
			self.generator.Qmin = mpc.gen(:,5)./self.mBase;
			self.generator.votage = mpc.gen(:,6);
			self.generator.status = mpc.gen(:,8);
			self.generator.Pmax = mpc.gen(:,9)./self.mBase;
			self.generator.Pmin = mpc.gen(:,10)./self.mBase;
			self.generator.type = zeros(gensLength, 1);
			for k = 1:gensLength
				self.generator.type(k) = self.nodes.type(self.generator.nid(k) == self.nodes.id);
			end

			self.branches.fid = mpc.branch(:, 1);
			self.branches.tid = mpc.branch(:, 2);
			branchesLength = length(self.branches.fid);
			self.branches.id = (1:branchesLength)';
			self.branches.r = mpc.branch(:, 3);
			self.branches.x = mpc.branch(:, 4);
			self.branches.g = zeros(branchesLength, 1);	% 目前我们不支持添加线路对地导纳
			self.branches.b = mpc.branch(:, 5);
			self.branches.ratio = mpc.branch(:, 9);
			self.branches.angle = mpc.branch(:, 10).*pi./360;
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

			% Mar 21 2018	设立节点电压 x, y分量, 用于机网转换
			self.nodes.Ux = zeros(nodesLength, 1);
			self.nodes.Uy = zeros(nodesLength, 1);
		end

		%% getNodeData: 生成节点参数,主要是带有独立导纳设备的节点参数
		% 这里将节点上的并联电容视为恒阻抗模型,并将其归算至节点导纳矩阵
		function [nodeData] = getNodeData(self)
			nodeData = [self.nodes.id, self.nodes.g, self.nodes.b];
		end
		%% getLineData: 生成线路参数,包括线路和变压器
		function [lineData] = getLineData(self)
			% 起始节点	终止节点	线路电阻	线路电抗	线路对地电导	线路对地电纳	变比
			lineData = [self.branches.fid, self.branches.tid, self.branches.r, self.branches.x, self.branches.g, self.branches.b, self.branches.ratio, self.branches.angle];
		end

		%% NAMInit: 计算节点导纳矩阵, 结果存放于电力网相应属性中
		% @return void
		function NAMInit(self)
			self.NAM = Model.NAM();
			self.NAM.generate(self.getNodeData(), self.getLineData());
		end
		
		%% votageInit: 获取迭代初始值(考虑对节点的设置及发电机的设置), 结果存放于 nodes 字段相应属性中
		% @return void
		function votageInit(self, solver)

			if isstruct(solver) && isfield(solver, 'start') && strcmp(solver.start, 'default')
				% do nothing			
			else	% 0-1 启动
				self.nodes.mag(self.pqIndex) = 1.05;
				self.nodes.ang([self.pqIndex; self.pvIndex]) = 0;
			end
		end

		%% planInit: 初始化,用于计算各节点功率计划值,并根据发电机的情况计算出该节点可发出的的最大无功功率, 结果存放于 nodes 字段相应属性中
		% @return void
		function planInit(self)
			for k = 1:length(self.generator.id)
				index = find(self.nodes.id == self.generator.nid(k));
				% 这里对各节点的功率计划的赋值有一大部分(PQ节点)是没有意义的
				self.nodes.Pg(index) = self.nodes.Pg(index) + self.generator.Pg(k);
				self.nodes.Qg(index) = self.nodes.Qg(index) + self.generator.Qg(k);
				self.nodes.Qmin(index) = self.nodes.Qmin(index) + self.generator.Qmin(k);
				self.nodes.Qmax(index) = self.nodes.Qmax(index) + self.generator.Qmax(k);
			end
		end

		%% planUpdate: 更新各节点功率计划值, 结果存放于 nodes 字段相应属性中
		% @return void
		function planUpdate(self)
			self.nodes.Pis = self.nodes.Pg - self.nodes.Pd;
			self.nodes.Qis = self.nodes.Qg - self.nodes.Qd;
		end

		%% setPowerOutflow: 计算从 nid 节点注入电网的潮流, 结果存放于 nodes 字段相应属性中
		% @param  int32   nid   节点 id
		% @return void
		function setPowerOutflow(self, nid)
			for k = 1:length(nid)
				index = find(self.nodes.id == nid(k));
				% disp(index)
				Sf = conj(self.NAM.value(index, :))*(self.nodes.mag.*exp(i.*(self.nodes.ang(index)-self.nodes.ang))).*self.nodes.mag(index);
				self.nodes.Pout(k) = real(Sf);
				self.nodes.Qout(k) = imag(Sf);
			end
		end

		%% setPowerUnbalance: 计算功率不平衡量, 仅用于牛拉法和 PQ 分解法计算, 返回变量均涵盖所有节点并已分类, 结果存放于 nodes 字段相应属性中
		% @return void
		function setPowerUnbalance(self)
			self.nodes.dP = self.nodes.Pis - self.nodes.Pout;
			self.nodes.dQ = self.nodes.Qis - self.nodes.Qout;
		end

		%% setIterationData: 每次使用潮流方程之后调用此方法获取迭代信息, 结果存放于 nodes 字段相应属性中
		% @return void
		function setIterationData(self)
			% TODO 优化内存使用方式
			self.itlog.mag = [self.itlog.mag, self.nodes.mag];
			self.itlog.ang = [self.itlog.ang, self.nodes.ang];
			self.itlog.dP = [self.itlog.dP, self.nodes.dP];
			self.itlog.dQ = [self.itlog.dQ, self.nodes.dQ];
		end

		%% checkConverged: 判断是否收敛, 同时检查是否结束循环
		% @param  struct  solver         与判断收敛相关的设置
		% --attr  double  epsilon        允许的功率误差
		% --attr  int32   maxIteration   最大迭代次数
		% @return int32   status       判断结果, 0 表示未收敛, 1 表示已收敛
		function [status] = checkConverged(self, solver, it)
			% 收敛判据
			if solver.epsilon >= max(abs([self.nodes.dQ; self.nodes.dP]))
				status = 1;	% 收敛
			elseif it >= solver.maxIteration
				warning('Exception 101: Number of iterations exceeds the limit\n');
				status = 101;	% 迭代次数超出上限
			elseif any(self.nodes.mag > 1.5) || any(self.nodes.mag < 0.5)
				warning('Exception 102: Some node have abnormal voltages and iteration is terminated\n');
				status = 102;	% 部分节点电压不正常
			else
				self.magStep = min(max(1.0 + 0.1.*log10(abs(self.nodes.dQ([self.pqIndex]))), 0.05), 1.05);
				self.angStep = min(max(1.0 + 0.1.*log10(abs(self.nodes.dP([self.pqIndex; self.pvIndex]))), 0.05), 1.05);
				status = 0;	% 未收敛, 继续计算
			end
		end

		%% getJacobian: 在当前网络状态下获取雅可比矩阵,用于牛顿法的计算
		% @return n*n-double jacobianMatrix   雅克比矩阵
		function [jacobianMatrix] = getJacobian(self)

			% 为方便构造雅可比矩阵,现分别对 PQ 节点及 PV 节点提取出节点电压及节点导纳, 并计算节点之间的相角矩阵
			% 本方法中, 变量名以 1 结尾的均与雅可比矩阵中的 H N M 有关, 变量名以 2 结尾的均与雅可比矩阵中的 N M L 有关
			NAM1 = self.NAM.value([self.pqIndex; self.pvIndex], [self.pqIndex; self.pvIndex]);
			NAM2 = self.NAM.value([self.pqIndex], [self.pqIndex]);

			mag1 = self.nodes.mag([self.pqIndex; self.pvIndex]);
			ang1 = self.nodes.ang([self.pqIndex; self.pvIndex]);
			mag2 = self.nodes.mag([self.pqIndex]);
			ang2 = self.nodes.ang([self.pqIndex]);

			nodeAngle = ang1*ones(1,length(ang1));
			nodeAngle = nodeAngle - nodeAngle';

			% 开始计算雅可比矩阵,创建初始状态. jacPrimary 是一个复数矩阵, 雅克比矩阵的元素主要来自这个矩阵
			jacPrimary = mag1*(mag1').*conj(NAM1).*exp(i.*nodeAngle);
			% 这个结果的对角元无实际意义, 去掉对角元
			jacPrimary = jacPrimary - diag(diag(jacPrimary));

			% 先计算各个分块矩阵的非对角元
			% 需要注意的是, 下面的两个索引表示在雅可比矩阵的初始状态中 PQPV 节点及 PQ 节点所在的位置, 而 self.pqIndex 表示在节点对象中 PQ 节点所在的位置
			pqpv = 1:(self.pqCount + self.pvCount);
			pq = 1:self.pqCount;

			% 雅可比矩阵的非对角元均可以在初始状态矩阵中提取
			Hij = imag(jacPrimary([pqpv], [pqpv]));	% Hij 表示系统中有功潮流差与相角差的数量关系
			Nij = real(jacPrimary([pqpv], [pq]));	% Nij 表示系统中有功潮流差与电压差的相对值的数量关系
			Mij = -real(jacPrimary([pq], [pqpv]));	% Mij 表示系统中无功潮流差与相角差的数量关系
			Lij = imag(jacPrimary([pq], [pq]));	% Lij 表示系统中无功潮流差与电压差的相对值的数量关系

			% 再计算各个分块矩阵的对角元. 求导的原因, 使得对角元需要单独计算
			Hii = -self.nodes.Qout([self.pqIndex; self.pvIndex]) - mag1.^2.*imag(diag(NAM1));
			Nii = self.nodes.Pout(self.pqIndex) + mag2.^2.*real(diag(NAM2));
			Mii = self.nodes.Pout(self.pqIndex) - mag2.^2.*real(diag(NAM2));
			Lii = self.nodes.Qout(self.pqIndex) - mag2.^2.*imag(diag(NAM2));

			Hij = Hij + diag(Hii);
			Nij(pq, pq) = Nij(pq, pq) + diag(Nii);
			Mij(pq, pq) = Mij(pq, pq) + diag(Mii);
			Lij = Lij + diag(Lii);

			% 在这里拼成雅可比矩阵
			jacobianMatrix = [Hij, Nij; Mij, Lij;];
		end

		%% getFDM: 返回 PQ 分解法的 B 矩阵
		% @param  string     fdType    PQ 分解法的系数矩阵类型(FD, FDBX, FDXB)
		% @param  array      nodesType 矩阵中包含的节点类型, 要得到 B' 时传 [1, 2], 要得到 B'' 时传 [1]
		% @return n*n double 系数矩阵
		function [matrix] = getFDM(self, fdType, nodesType)

			assert(ismember(1, nodesType) && ~ismember(3, nodesType));

			matrix = sparse(zeros(self.pqCount + self.pvCount + self.refCount));
			switch fdType
				case 'r'	% 考虑线路电阻
					for k = 1:length(self.branches.id)
						fi = find(self.nodes.id == self.branches.fid(k));
						ti = find(self.nodes.id == self.branches.tid(k));
						b = self.branches.x(k)./(self.branches.r(k).^2 + self.branches.x(k).^2);
						matrix([fi, ti], [fi, ti]) = matrix([fi, ti], [fi, ti]) + [b, -b; -b, b];
					end
				case 'b'	% 考虑节点电纳
					for k = 1:length(self.branches.id)
						fi = find(self.nodes.id == self.branches.fid(k));
						ti = find(self.nodes.id == self.branches.tid(k));
						b = 1./self.branches.x(k);
						matrix([fi, ti], [fi, ti]) = matrix([fi, ti], [fi, ti]) + [b, -b; -b, b];
					end
					matrix = matrix - diag(self.nodes.b);
				otherwise
					error('illegal fast duplicate type in getFDM');
			end

			index = find(self.nodes.type == 1);
			if ismember(2, nodesType)
				index = [index; find(self.nodes.type == 2)];
			end

			matrix = -matrix(index, index);
		end

		%% setFDM: 计算 PQ 分解法的系数矩阵
		% @param  string     fdType    PQ 分解法的系数矩阵类型(FD, FDBX, FDXB)
		% @return  void
		function setFDM(self, fdType)
			switch fdType
				case 'FD'
					self.FDM1 = imag(self.NAM.value([self.pqIndex; self.pvIndex],[self.pqIndex; self.pvIndex]));
					self.FDM2 = imag(self.NAM.value([self.pqIndex], [self.pqIndex]));
				case 'FDBX'
					self.FDM1 = self.getFDM('r', [1, 2]);	% 考虑线路电阻
					self.FDM2 = self.getFDM('b', [1]);		% 考虑节点电纳
				case 'FDXB'
					self.FDM1 = self.getFDM('b', [1, 2]);
					self.FDM2 = self.getFDM('r', [1]);
				otherwise
					error('illegal fast duplicate type in setFDM');
			end
		end

		%% getVotageAndAngleCorrection: 修正方程, 用于牛顿法
		%　@param  n*n-double  jacobianMatrix   雅克比矩阵
		% @return double      dTheta           相角修正量, 用于迭代方程
		% @return double      dUpn             电压修正量标幺值, 用于迭代方程
		function [dTheta,dUpn] = getVotageAndAngleCorrection(self, jacobianMatrix, standardFlag)

			dP = self.nodes.dP([self.pqIndex; self.pvIndex]);
			dQ = self.nodes.dQ([self.pqIndex]);

			dPara = jacobianMatrix \ [dP; dQ];

			dTheta = dPara(1:(self.pqCount + self.pvCount));	% PQ 及 PV 节点相角的修正量
			dUpn = dPara((self.pqCount + self.pvCount + 1):end);	% PQ 节点电压相对值的修正量

			if nargin == 3 && standardFlag == true
				% 将电压和相角误差规范化, 即填充至大向量中, 方便迭代方程直接使用
				dUpn_t = zeros(self.pqCount + self.pvCount + self.refCount, 1);	% 初始化相角修正量
				dUpn_t([self.pqIndex]) = dUpn;	% 代入 PQ 及 PV 节点相角的修正量
				dUpn = dUpn_t;
				dTheta_t = zeros(self.pqCount + self.pvCount + self.refCount, 1);
				dTheta_t([self.pqIndex; self.pvIndex]) = dTheta;
				dTheta = dTheta_t;
			end
		end

		%% getMagCorrection: 计算电压修正量, 用于 PQ 分解法
		function [dUpn] = getMagCorrection(self, standardFlag)
			if nargin == 2 && standardFlag == true
				dUpn = zeros(self.pqCount + self.pvCount + self.refCount, 1);
				dUpn([self.pqIndex]) = self.FDM2 \ (self.nodes.dQ([self.pqIndex])./self.nodes.mag([self.pqIndex]));
			else
				dUpn = self.FDM2 \ (self.nodes.dQ([self.pqIndex])./self.nodes.mag([self.pqIndex]));
			end
		end

		%% getAngCorrection: 计算相角修正量, 用于 PQ 分解法
		function [dTheta] = getAngCorrection(self, standardFlag)
			if nargin == 2 && standardFlag == true
				dTheta = zeros(self.pqCount + self.pvCount + self.refCount, 1);
				dTheta([self.pqIndex; self.pvIndex]) = (self.FDM1 \ (self.nodes.dP([self.pqIndex; self.pvIndex])./self.nodes.mag([self.pqIndex; self.pvIndex])))./self.nodes.mag([self.pqIndex; self.pvIndex]);
			else
				dTheta = (self.FDM1 \ (self.nodes.dP([self.pqIndex; self.pvIndex])./self.nodes.mag([self.pqIndex; self.pvIndex])))./self.nodes.mag([self.pqIndex; self.pvIndex]);
			end
		end

		%% setNewVotage: 获取新的节点电压及相角, 用于牛顿法. 更新电力网 nodes 字段的属性 mag 和 ang
		% @param  double dTheta    相角修正量, 来自修正方程
		% @param  double dUpn      电压修正量标幺值, 来自修正方程
		% @return void
		function setNewVotage(self, dTheta, dUpn, standardFlag)
			if nargin == 4 && standardFlag == true
				self.nodes.ang = self.nodes.ang + (dTheta);
				self.nodes.mag = self.nodes.mag.*(1+dUpn);
			else
				self.nodes.ang([self.pqIndex; self.pvIndex]) = self.nodes.ang([self.pqIndex; self.pvIndex]) + (dTheta);
				self.nodes.mag([self.pqIndex]) = self.nodes.mag([self.pqIndex]).*(1+dUpn);
			end
		end

		%% setGeneratorPower: 计算各节点的实际发电量 (仅计算了 pv 节点的无功和平衡节点的有功无功)
		%  结果存放于电力网 node 字段的属性 Pg 和 Qg
		% @return void
		function setGeneratorPower(self)
			% 每个节点发电机发出的功率都等于当地需求和收敛时潮流方程计算的注入电网的功率的和
			self.nodes.Pg(self.refIndex) = self.nodes.Pout(self.refIndex) + self.nodes.Pd(self.refIndex);
			self.nodes.Qg([self.pvIndex; self.refIndex]) = self.nodes.Qout([self.pvIndex;  self.refIndex]) + self.nodes.Qd([self.pvIndex; self.refIndex]);
		end

		%% getLinePower: 计算普通线路的潮流
		% @param  int32   index    线路 id(index)
		% @return complex Sij      从起始节点到终止节点的复功率, 标幺值
		% @return complex Sji      从终止节点到起始节点的复功率, 标幺值
		% @return complex dS       线路中的损耗, 但线路充电电容的影响也考虑了进去. 标幺值
		function [Sij,Sji,dS] = getLinePower(self, index)
			%　计算三个导纳值，计算功率时会用到
			conjYi0 = conj((self.branches.g(index)+self.branches.b(index)*1i)./2);
			conjYj0 = conjYi0;
			conjYij = 1./conj((self.branches.r(index)+self.branches.x(index)*1i));

			fid = find(self.nodes.id == self.branches.fid(index));	% 得到线路始末端节点的id
			tid = find(self.nodes.id == self.branches.tid(index));	% 得到线路始末端节点的id

			% 使用公式计算
			Sij = (self.nodes.mag(fid)).^2.*(conjYi0 + conjYij) - self.nodes.mag(fid).*self.nodes.mag(tid).*conjYij.*exp((self.nodes.ang(fid)-self.nodes.ang(tid)).*i);
			Sji = (self.nodes.mag(tid)).^2.*(conjYj0 + conjYij) - self.nodes.mag(fid).*self.nodes.mag(tid).*conjYij.*exp((self.nodes.ang(tid)-self.nodes.ang(fid)).*i);
			dS = Sij + Sji;
		end

		%% getTransformerPower: 计算变压器的潮流
		% @param  int32   index    线路 id(index)
		% @return complex Sij      从起始节点到终止节点的复功率, 标幺值
		% @return complex Sji      从终止节点到起始节点的复功率, 标幺值
		% @return complex dS       线路中的损耗, 但线路充电电容的影响也考虑了进去. 标幺值
		function [Sij,Sji,dS] = getTransformerPower(self, index)
			%　首先确定变压器　pi 型等效电路的三个参数
			[trZ, trY1, trY2] = gamma2pi((self.branches.r(index)+self.branches.x(index)*1i), (self.branches.g(index)+self.branches.b(index)*1i), self.branches.ratio(index));
			%　计算三个导纳值，计算功率时会用到
			conjYi0 = conj(trY1) - (self.branches.g(index)+self.branches.b(index)*1i)./2;
			conjYj0 = conj(trY2) - (self.branches.g(index)+self.branches.b(index)*1i)./2;
			conjYij = conj(1./trZ);

			fid = find(self.nodes.id == self.branches.fid(index));	% 得到线路始始端节点的id
			tid = find(self.nodes.id == self.branches.tid(index));	% 得到线路始末端节点的id
			
			% 使用公式计算
			Sij = (self.nodes.mag(fid)).^2.*(conjYi0 + conjYij) - self.nodes.mag(fid).*self.nodes.mag(tid).*conjYij.*exp((self.nodes.ang(fid)-self.nodes.ang(tid)).*i);
			Sji = (self.nodes.mag(tid)).^2.*(conjYj0 + conjYij) - self.nodes.mag(fid).*self.nodes.mag(tid).*conjYij.*exp((self.nodes.ang(tid)-self.nodes.ang(fid)).*i);
			dS = Sij + Sji;
		end
		
		%% setBranchPower: 计算线路功率及功率损耗
		%  该方法在潮流迭代结束之后再使用, 用于求解收敛时电力网中各线路传输的功率, 求解完成后更新电力网 branches 字段的相应属性
		% @return void
		function setBranchPower(self)
			Sij = zeros(length(self.branches.id), 1);
			Sji = zeros(length(self.branches.id), 1);
			dS = zeros(length(self.branches.id), 1);
			for k = 1:length(self.branches.id)
				if self.isLine(k)
					[Sij(k), Sji(k), dS(k)] = self.getLinePower(k);
				else
					[Sij(k), Sji(k), dS(k)] = self.getTransformerPower(k);
				end	
			end
			self.branches.Pij = real(Sij);
			self.branches.Qij = imag(Sij);
			self.branches.Pji = real(Sji);
			self.branches.Qji = imag(Sji);
			self.branches.dP = real(dS);
			self.branches.dQ = imag(dS);
		end

		%% setConpensatorPower: 计算无功补偿器电量
		% @return void
		function setConpensatorPower(self)
			self.nodes.Pc = self.nodes.mag.^2.*self.nodes.g;
			self.nodes.Qc = self.nodes.mag.^2.*self.nodes.b;
		end

		%% setVotageRealAndImag: 将节点电压转化为直角坐标形式, 存储在 nodes.Ux, nodes.Uy 中
		% @return void
		function setVotageRealAndImag(self)
			self.nodes.Ux = self.nodes.mag.*cos(self.nodes.ang);
			self.nodes.Uy = self.nodes.mag.*sin(self.nodes.ang);
		end

		%% NRIteration: 牛顿法迭代程序核心部分
		% @param    struct  solver	        对迭代程序的基本设置
		% --attr    string  method            求解方法, NR FDBX FDXB
		% --attr    int32   maxIteration      最大迭代次数
		% --attr    double  epsilon           收敛判据, 功率不平衡量标幺
		% --attr    string  start             启动方式, default 为按发电机端电压起动
		% --attr    string  documentName      如果是文本输出, 需设置文本计算报告 src
		% @return   struct  result          迭代结果
		% --attr    int32   status            错误码, 为 0 表示正常
		% --attr    int32   it                迭代次数
		function [result] = NRIteration(self, solver)
			% 从这里开始进入循环,跳出循环的条件为迭代次数超过最大次数或潮流不平衡量的一范数小于给定值
			result.it = 0;

			while 1
				self.planUpdate();	% 确定各节点功率需求
				self.setPowerOutflow(self.nodes.id);	% 潮流方程
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
		% --attr    string  method            求解方法, NR FD FDBX FDXB
		% --attr    int32   maxIteration      最大迭代次数
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
				self.setPowerOutflow(self.nodes.id);	% 潮流方程
				self.setPowerUnbalance();	% 误差方程

				self.setGeneratorPower();
				self.setIterationData();	% 存储迭代信息,包括电压,相角,有功不平衡量,无功不平衡量

				% 判收敛以及是否迭代次数上限, 返回一个状态码, 得到非零值直接返回
				result.status = self.checkConverged(solver, result.it);
				if(result.status ~= 0)	% result.status 非 0 表示收敛或错误
					return;
				end

				dUpn = self.getMagCorrection();
				self.setNewVotage(0, -dUpn.*self.magStep);	% 迭代方程

				dTheta = self.getAngCorrection();
				self.setNewVotage(-dTheta.*self.angStep, 0);	% 迭代方程

				result.it = result.it + 1;
				% ❀ 运用各节点电压的新值进行下一步的迭代.
			end
		end
		
		function [output] = checkReactivePower(self)

			% 下面两句 返回的 nodesIndex 包括所有类型的节点, 但无需剔除不必要的节点
			lack = find(self.nodes.Qmax < self.nodes.Qg);
			excess = find(self.nodes.Qmin > self.nodes.Qg);

			% 只保留 PV 节点
			lack(find(self.nodes.type(lack) ~=2 )) = [];
			excess(find(self.nodes.type(excess) ~=2 )) = [];

			if isempty([lack; excess])
				output = false;
				return;
			end
			output = true;

			self.convertToPQNode(lack, excess);

			for k = 1:length(lack)
				fprintf('%s %d %s\n', '	  Node ', self.nodes.id(lack(k)),' has insufficient reactive power and has been converted to a PQ node');
			end
			for k = 1:length(excess)
				fprintf('%s %d %s\n', '	  Node ', self.nodes.id(excess(k)),' has excess reactive power and has been converted to a PQ node');
			end
		end

		%% convertToPQNode: 解决某些节点无功不足或过剩的情况
		function convertToPQNode(self, lack, excess)

			self.nodes.Qg(lack) = self.nodes.Qmax(lack);
			self.nodes.Qg(excess) = self.nodes.Qmin(excess);

			% 修改节点类型,将其转化为 PQ 节点,并更新节点对象相应的的 PQ PV 节点的索引及数量属性
			self.nodes.type([lack; excess]) = 1;

		end

		%% solvePowerFlow: 电力系统潮流计算
		function [result] = solvePowerFlow(self, solver)
			if nargin ~= 2
				error('illegal power flow solver');
			end

			% ❀ 生成节点导纳矩阵.
			self.NAMInit();
			% ❀ 设置节点电压的初值(电压,相角).
			self.votageInit(solver);
			% ❀ 将各节点电压的初值代入潮流方程,并求解有功无功不平衡量
			self.planInit();	% 计算各节点功率的计划值

			% 这个循环是与检测无功功率有关的, 与具体的求解无关
			while 1

				% ❀ 进入迭代. 得到的 result 包含是否收敛, 迭代次数等信息
				switch solver.method
					case 'NR'	% 牛顿法
						result = self.NRIteration(solver);
					case {'FD', 'FDBX', 'FDXB'}	% PQ 分解法
						result = self.PQIteration(solver);
					% case 'FDBX'	% PQ 分解法
					% 	result = self.PQIteration(solver);
					% case 'FDXB'	% PQ 分解法
					% 	result = self.PQIteration(solver);
					otherwise
						% TODO otherwise
						error('illegal solver');
				end

				% 未收敛时直接返回
				if(result.status ~= 1)
					return;
				end
				
				self.setGeneratorPower();	% 根据迭代结果计算发电机功率, 结果记录到 nodes 属性中
				if ~isfield(solver, 'checkReactivePower') || solver.checkReactivePower ~= true || ~self.checkReactivePower()
					break;
				end
				fprintf('continue\n');
			end

			self.setBranchPower();	% 根据迭代结果计算线路功率及损耗, 结果记录到 branches 属性中
			self.setConpensatorPower();	% 根据迭代结果计算各节点对地电容功率, 结果记录到 nodes 属性中
			% self.setVotageRealAndImag();	% 若需要作后续暂态分析, 可以做这一步的电压转化
		end

		% addGenImpedanceToNodes: 将发电机的暂态电抗加至各节点, 为计算暂态参数做准备
		% 不是给潮流分析留的接口
		function addGenImpedanceToNodes(self, generator)
			if nargin == 1
				return;
			end

			index = getIndex(self.nodes.id, generator.nid);
			yg = 1./(generator.Ra + generator.Xd1.*1i);
			
			self.nodes.g(index) = self.nodes.g(index) + real(yg);
			self.nodes.b(index) = self.nodes.b(index) + imag(yg);
			self.NAM.addAdmittance(index, yg);
		end

		% addMoterImpedanceToNodes: 将异步电动机的暂态电抗加至各节点, 为计算暂态参数做准备
		% 不是给潮流分析留的接口
		function addMoterImpedanceToNodes(self, moter)
			% TODO addMoterImpedanceToNodes
			if nargin == 1
				return;
			end
		end

		% addLoadImpedanceToNodes: 将负荷的电抗加至节点, 为计算暂态参数做准备
		% 不是给潮流分析留的接口
		function addLoadImpedanceToNodes(self)
			self.nodes.g = self.nodes.g + self.nodes.Pd./self.nodes.mag.^2;
			self.nodes.b = self.nodes.b - self.nodes.Qd./self.nodes.mag.^2;
			self.NAM.addAdmittance(1:length(self.nodes.id), (self.nodes.Pd - self.nodes.Qd.*1i)./self.nodes.mag.^2);
		end

		% getShortCircultCapacity: 计算单个节点的短路容量
		function [result] = getShortCircultCapacity(self, id)
			
			% 1. 根据发电机信息, 修改所有节点的接地参数;
			self.addLoadImpedanceToNodes();
			self.addGenImpedanceToNodes();
			% 2. 生成节点导纳矩阵;
			self.NAMInit();
			% 3. 解方程 Yz = ek;
			k = find(self.nodes.id == id);
			ek = zeros(length(self.nodes.id), 1);
			ek(k) = 1;
			z = self.NAM \ ek;
			% 4. 计算短路容量, S = 1/zk, S = 1.732*Un/zk
			result.scc = 1./z(k);
			result.z = z;
		end

		% getShortCircultCapacity: 计算所有节点的短路容量
		function [result] = getAllShortCircultCapacity(self)
			
			% 1. 根据发电机信息, 修改所有节点的接地参数;
			% self.addLoadImpedanceToNodes(nodes);
			self.addGenImpedanceToNodes();
			% 2. 生成节点导纳矩阵;
			self.NAMInit();
			% 3. 求逆;
			Z = inv(self.NAM);
			% 4. 计算短路容量, S = 1/zk, S = 1.732*Un/zk
			result.scc = 1./diag(Z);
			result.Z = Z;
		end

	end
end