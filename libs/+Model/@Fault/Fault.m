%% Fault
classdef Fault < handle

	properties

		ss;			% 电力系统稳态模型

		nodes;		% 节点故障参数
		generator;	% 发电机故障参数
		branches;	% 线路故障参数

		fault;		% 故障原始信息

		NAM1;	% 正序网络节点导纳矩阵
		NAM2;	% 负序网络节点导纳矩阵
		NAM0;	% 零序网络节点导纳矩阵

		basePhase;
		T;	% 120 变换矩阵
		T_inv;	% 120 反变换矩阵

		itlog;

	end

	properties (Dependent)
		lineIndex;
		lineCount;
		transformerIndex;
		transformerCount;
	end

	methods

		%% get.transformerIndex: 返回变压器的索引
		function [transformerIndex] = get.transformerIndex(self)
			transformerIndex = 1:length(self.branches.xn);
		end

		%% get.transformerCount: 返回变压器的数量
		function [transformerCount] = get.transformerCount(self)
			transformerCount = length(self.branches.xn);
		end

		%% get.lineIndex: 返回变压器的索引
		function [lineIndex] = get.lineIndex(self)
			lineIndex = (self.transformerCount + 1):length(self.branches.x1);
		end

		%% get.lineCount: 返回变压器的数量
		function [lineCount] = get.lineCount(self)
			lineCount = length(self.branches.x1) - length(self.branches.xn);
		end		

		%% init: 电力系统故障分析模型初始化
		function init(self, mpc)
			% 故障信息, 节点导纳矩阵初始化
			self.fault = mpc.fault;
			self.NAM1 = Model.NAM();
			self.NAM2 = Model.NAM();
			self.NAM0 = Model.NAM();

			% 分别初始化线路和变压器故障参数
			% lineCount = size(mpc.line)(1);
			transformerCount = size(mpc.transformer);
			transformerCount = transformerCount(1);

			% 下面这些字段都最终都是列向量或列细胞数组
			self.branches.fid = [[mpc.transformer{:, 1}]'; mpc.line(:, 1)];
			self.branches.tid = [[mpc.transformer{:, 2}]'; mpc.line(:, 2)];
			self.branches.r = [[mpc.transformer{:, 3}]'; mpc.line(:, 3)];
			self.branches.x1 = [[mpc.transformer{:, 4}]'; mpc.line(:, 4)];
			self.branches.x0 = [[mpc.transformer{:, 5}]'; mpc.line(:, 5)];
			self.branches.xn = [mpc.transformer{:, 6}]';		% 没有线路的
			self.branches.ratio = [mpc.transformer{:, 7}]';		% 没有线路的
			self.branches.angle = [mpc.transformer{:, 8}]';		% 没有线路的
			self.branches.group = {mpc.transformer{:, 9}}';	% cell array

			self.nodes.Ua = [];
			self.nodes.Ub = [];
			self.nodes.Uc = [];
			self.nodes.U1 = [];
			self.nodes.U2 = [];
			self.nodes.U0 = [];
			self.branches.Ia = [];
			self.branches.Ib = [];
			self.branches.Ic = [];
			self.branches.I1 = [];
			self.branches.I2 = [];
			self.branches.I0 = [];

			self.generator.nid = [mpc.gen{:, 1}]';
			self.generator.r = [mpc.gen{:, 2}]';
			self.generator.x1 = [mpc.gen{:, 3}]';
			self.generator.x2 = [mpc.gen{:, 4}]';
			self.generator.x0 = [mpc.gen{:, 5}]';
			self.generator.xn = [mpc.gen{:, 6}]';
			self.generator.group = {mpc.gen{:, 7}}';	% cell array
		end

		%% set120MatrixByBasePhase: 根据基本相返回 120 变换矩阵
		function [matrix] = set120MatrixByBasePhase(self, basePhase)
			a = -0.5 + sqrt(3)./2.*1i;
			c = -0.5 - sqrt(3)./2.*1i;
			switch basePhase
				case 'a'
					self.T = [
						1	1	1;
						c	a	1;
						a	c	1;
					];
				case 'b'
					self.T = [
						a	c	1;
						1	1	1;
						c	a	1;
					];
				case 'c'
					self.T = [
						c	a	1;
						a	c	1;
						1	1	1;
					];
				case 'none'
					self.T = [];
				otherwise
					error('illegal base phase in set120MatrixByBasePhase');
			end
		end

		%% set120InvMatrixByBasePhase: 根据基本相返回 120 反变换矩阵
		function set120InvMatrixByBasePhase(self, basePhase)
			e = 1./3;
			a = -1./6 + sqrt(3)./6.*1i;
			c = -1./6 - sqrt(3)./6.*1i;
			switch basePhase
				case 'a'
					self.T_inv = [
						e	a	c;
						e	c	a;
						e	e	e;
					];
				case 'b'
					self.T_inv = [
						c	e	a;
						a	e	c;
						e	e	e;
					];
				case 'c'
					self.T_inv = [
						a	c	e;
						c	a	e;
						e	e	e;
					];
				case 'none'
					self.T_inv = [];
				otherwise
					error('illegal base phase in set120InvMatrixByBasePhase');
			end
		end

		%% getBasePhase: 返回基本相
		function [basePhase] = getBasePhase(self)
			switch self.fault.type
				case 'f3'
					basePhase = 'none';
				case {'f2', 'f11'}
					if ~isempty(regexpi(self.fault.phase, '(ab|ba)', 'match'))
						basePhase = 'c';
					elseif ~isempty(regexpi(self.fault.phase, '(ac|ca)', 'match'))
						basePhase = 'b';
					elseif ~isempty(regexpi(self.fault.phase, '(bc|cb)', 'match'))
						basePhase = 'a';
					else
						error('illegal fault phase');
					end
				case 'f1'
					if ~isempty(regexpi(self.fault.phase, 'a', 'match'))
						basePhase = 'a';
					elseif ~isempty(regexpi(self.fault.phase, 'b', 'match'))
						basePhase = 'b';
					elseif ~isempty(regexpi(self.fault.phase, 'c', 'match'))
						basePhase = 'c';
					else
						error('illegal fault phase');
					end
				otherwise
					error('illegal fault type');
			end
		end

		%% checkTransformerGroup: 检查联结组的合法性, 并返回联结组信息
		function [group] = checkTransformerGroup(self, string)

			group = regexpi(string, '(?<from>Y|YN|Y0|D)(?<to>y|yn|y0|d)(?<number>\d+)', 'names');
			group.number = str2num(group.number);

			matchSize = size(group);
			if matchSize(2) ~= 1	% 匹配不到或匹配多个
				error(['illegal transformer connection group: ', string]);
			end

			if strcmp(group.from, 'D') && strcmp(group.to, 'd')	% 不能为 Dd*
				error('connection group cannot be Dd*');
			end
			if group.number < 0 || group.number > 12
				error(['connection group number cannot be ', group.number]);
			end
			if (strcmp(group.from, 'D') || strcmp(group.to, 'd')) && mod(group.number, 2) == 0
				error(['connection number should be odd in ', group.from, group.to, ' connection']);
			end
			% 至此联结组标号合法性检查完成
		end

		%% getZeroSequenceInfo(self, id, group): 返回变压器零序参数向节点导纳矩阵中添加的信息
		function [nodes, item] = getZeroSequenceInfo(self, id, group)
			% 不计零序电阻
			if strcmp(group.from, 'D') && ~strcmp(group.to, 'y')
				nodes = 'to';
				item = 1./((self.branches.x0(id) + 3.*self.branches.xn(id)).*1i);	% 向节点上添加的是导纳
				return ;
			end
			if strcmp(group.to, 'd') && ~strcmp(group.from, 'Y')
				nodes = 'from';
				item = 1./((self.branches.x0(id) + 3.*self.branches.xn(id)).*1i);	% 向节点上添加的是导纳
				return ;
			end
			if ~strcmp(group.from, 'Y') && ~strcmp(group.to, 'y')
				nodes = 'both';
				item = (self.branches.x0(id) + 3.*self.branches.xn(id)).*1i;	% 向线路上添加的是阻抗
				return ;
			end
			nodes = 'none';
			item = [];
		end

		%% addLineToNAM: 将单条线路添加到正序负序零序节点导纳矩阵
		function addLineToNAM(self, id, seq)

			fi = find(self.ss.nodes.id == self.branches.fid(id));
			ti = find(self.ss.nodes.id == self.branches.tid(id));
			% 只考虑线路的阻抗, 不考虑充电电容. 由于线路的零序电抗已经考虑了互感, 故直接添加
			if findstr(seq, '1')
				self.NAM1.addImpedance(fi, ti, self.branches.r(id) + self.branches.x1(id).*1i);
			end
			if findstr(seq, '2')
				self.NAM2.addImpedance(fi, ti, self.branches.r(id) + self.branches.x1(id).*1i);
			end
			if findstr(seq, '0')
				self.NAM0.addImpedance(fi, ti, self.branches.x0(id).*1i);	% 不计零序电阻
			end
		end
		
		%% addTransformerToNAM: 将单台变压器添加到正序负序零序节点导纳矩阵
		function addTransformerToNAM(self, id, seq)

			fi = find(self.ss.nodes.id == self.branches.fid(id));
			ti = find(self.ss.nodes.id == self.branches.tid(id));

			% 先取得联结组信息
			group = self.checkTransformerGroup(self.branches.group{id});
			dAngle = (group.number./6.*pi);	% 由于联结组带来的移相效果

			% 正序负序比较简单, 可以忽略励磁
			% 正序需要按照 angle 和联结组号进行移相, 从一次侧到二次侧, 每点滞后 30 度; 负序需要按照 angle 和联结组号进行移相, 从一次侧到二次侧, 每点超前 30 度
			if findstr(seq, '1')
				self.NAM1.addTransformer(fi, ti, self.branches.r(id) + self.branches.x1(id).*1i, 0, self.branches.ratio(id).*exp(1i.*(self.branches.angle(id))));
			end
			if findstr(seq, '2')
				self.NAM2.addTransformer(fi, ti, self.branches.r(id) + self.branches.x1(id).*1i, 0, self.branches.ratio(id).*exp(1i.*(self.branches.angle(id) + 2.*dAngle)));
			end

			% 变压器的零序等效电路参数与变压器铁芯结构有关, 我们忽略了三柱变压器对零序励磁电抗带来的影响, 认为所有变压器的励磁电抗都是无穷大
			% 变压器的零序等效电路结构与变压器的联结组及标号有关. 下面分别讨论
			% 零序不需要根据点数移相, 也不用考虑变压器本身的移相功能
			if findstr(seq, '0')
				[zsNodes, zsItem] = self.getZeroSequenceInfo(id, group);	% TODO 阻抗应分一半移相
				switch zsNodes
					case 'from'
						self.NAM0.addAdmittance(self.branches.fid(id), zsItem);	% 环流 + 开路
					case 'to'
						self.NAM0.addAdmittance(self.branches.tid(id), zsItem);	% 环流 + 开路
					case 'both'
						% 零序能通过两侧都接地的变压器, 在导纳矩阵中应放通
						self.NAM0.addBranch(fi, ti, zsItem, 0, self.branches.ratio(id).*exp(1i.*(self.branches.angle(id) + dAngle)));
					otherwise	% 'none'
						% do nothing
				end
			end

		end

		%% addGeneratorToNAM: 向节点导纳矩阵添加发电机
		function addGeneratorToNAM(self, seq)
			index = getIndex(self.ss.nodes.id, self.generator.nid);
			if findstr(seq, '1')
				self.NAM1.addAdmittance(index, 1./(self.generator.r + self.generator.x1.*i));
			end
			if findstr(seq, '2')
				self.NAM2.addAdmittance(index, 1./(self.generator.r + self.generator.x2.*i));
			end

			if findstr(seq, '0')
				zeroPass = regexp(self.generator.group, '^[Yy][Nn0]$');
				for k = 1:length(self.generator.nid)
					if zeroPass{k}
						self.NAM0.addAdmittance(index, 1./((self.generator.x0 + 3.*self.generator.xn).*1i));
					end
				end
			end
		end

		%% addLoadsToNAM: 向节点导纳矩阵添加负荷
		function addLoadsToNAM(self, seq)
			if findstr(seq, '1')
				index1 = 1:length(self.ss.nodes.id);	% 正序
				self.NAM1.addAdmittance(index1, self.ss.nodes.Pd - self.ss.nodes.Qd.*1i);
			end
			if findstr(seq, '2')
				index2 = find(self.ss.nodes.Pd ~= 0);	% 负序, 有负荷的地方
				self.NAM2.addAdmittance(index2, ones(length(index2), 1)./(0.19 + 0.36.*1i));	% 普通负荷负序阻抗典型值
			end
			
			% 普通负荷不考虑零序
		end

		%% NAMInit: 节点导纳矩阵初始化, 填入故障前参数
		function NAMInit(self, seq)

			% 线路与变压器
			for k = self.transformerIndex
				self.addTransformerToNAM(k, seq);
			end
			for k = self.lineIndex
				self.addLineToNAM(k, seq);
			end

			% 发电机
			self.addGeneratorToNAM(seq);
			% 负荷
			self.addLoadsToNAM(seq);

		end
		

		%% getPowerFlowResult: 计算故障前电网的潮流
		function getPowerFlowResult(self, mpcSteady)

			% 建立电力网稳态模型并初始化
			self.ss = Model.SteadyState();
			self.ss.init(mpcSteady);

			%% 设置求解器的基本信息
			solver.method = 'FDBX';	% 求解方法
			solver.maxIteration = 50;	% 最大迭代
			solver.epsilon = 1e-5;	% 收敛判据, 功率不平衡量标幺
			solver.start = 'default';	% 启动方式, default 为按发电机端电压起动

			%% 求解
			self.ss.solvePowerFlow(solver);
		end

		%% calcSumInfo: 
		function [za, k2, k0] = calcSumInfo(self, faultType, z2f, z0f)
			switch faultType
				case 'f3'
					za = self.fault.zf;
					k2 = 0;
					k0 = 0;
				case 'f2'
					if nargin < 3 || ~isnumeric(z2f)
						error('illegal negative sequence impedace');	
					end
					za = self.fault.zf + z2f;
					k2 = -1;
					k0 = 0;
				case 'f1'
					if nargin < 4 || ~isnumeric(z2f) || ~isnumeric(z0f)
						error('illegal negative or zero sequence impedace');	
					end
					za = self.fault.zf + z2f + z0f;
					k2 = 1;
					k0 = 1;
				case 'f11'
					if nargin < 4 || ~isnumeric(z2f) || ~isnumeric(z0f)
						error('illegal negative or zero sequence impedace');	
					end
					za = self.fault.zf + getParallel([z2f + self.fault.zf, z0f + self.fault.zf + 3.*self.fault.zg]);
					k2 = -(z0f + self.fault.zf + 3.*self.fault.zg)./(z2f + z0f + 2.*self.fault.zf + 3.*self.fault.zg);
					k0 = -(z2f + self.fault.zf)./(z2f + z0f + 2.*self.fault.zf + 3.*self.fault.zg);
				otherwise
					error('illegal fault type in calcSumInfo');
			end
		end

		%% calcNetCurrent: 计算全网电流 (正序)
		function calcNetCurrentPS(self, U1)

			self.branches.I1 = zeros(length(self.branches.fid), 1);
			% TODO remove for loop
			for k = self.transformerIndex
				fi = find(self.ss.nodes.id == self.branches.fid(k));
				ti = find(self.ss.nodes.id == self.branches.tid(k));
				self.branches.I1(k) = (self.nodes.U1(fi).*self.branches.ratio(k).*exp(self.branches.angle(k).*1i) - self.nodes.U1(ti)) ./ (self.branches.r(k) + self.branches.x1(k).*1i);
			end
			for k = self.lineIndex
				fi = find(self.ss.nodes.id == self.branches.fid(k));
				ti = find(self.ss.nodes.id == self.branches.tid(k));
				self.branches.I1(k) = (self.nodes.U1(fi) - self.nodes.U1(ti)) ./ (self.branches.r(k) + self.branches.x1(k).*1i);
			end
		end

		%% calcNetCurrent: 计算全网电流 (负序)
		function calcNetCurrentNS(self, U2)

			self.branches.I2 = zeros(length(self.branches.fid), 1);
			for k = self.transformerIndex
				fi = find(self.ss.nodes.id == self.branches.fid(k));
				ti = find(self.ss.nodes.id == self.branches.tid(k));
				self.branches.I2(k) = (self.nodes.U2(fi).*self.branches.ratio(k).*exp(self.branches.angle(k).*1i) - self.nodes.U2(ti)) ./ (self.branches.r(k) + self.branches.x1(k).*1i);
			end
			for k = self.lineIndex
				fi = find(self.ss.nodes.id == self.branches.fid(k));
				ti = find(self.ss.nodes.id == self.branches.tid(k));
				self.branches.I2(k) = (self.nodes.U2(fi) - self.nodes.U2(ti)) ./ (self.branches.r(k) + self.branches.x1(k).*1i);
			end
		end

		%% calcNetCurrent: 计算全网电流 (零序)
		function calcNetCurrentZS(self, U0)

			self.branches.I0 = zeros(length(self.branches.fid), 1);
			for k = self.transformerIndex
				fi = find(self.ss.nodes.id == self.branches.fid(k));
				ti = find(self.ss.nodes.id == self.branches.tid(k));
				self.branches.I0(k) = (self.nodes.U0(fi).*self.branches.ratio(k).*exp(self.branches.angle(k).*1i) - self.nodes.U0(ti)) ./ (self.branches.r(k) + self.branches.x0(k).*1i);
			end
			for k = self.lineIndex
				fi = find(self.ss.nodes.id == self.branches.fid(k));
				ti = find(self.ss.nodes.id == self.branches.tid(k));
				self.branches.I0(k) = (self.nodes.U0(fi) - self.nodes.U0(ti)) ./ (self.branches.r(k) + self.branches.x0(k).*1i);
			end
		end

		%% solveThreePhaseShortCircult: 三相短路
		function solveThreePhaseShortCircult(self)
			% 三相短路, 无需考虑负序与零序
			self.NAM1.init(length(self.ss.nodes.id));
			self.NAMInit('1');

			%　求解节点阻抗矩阵的第　ｋ　列 (目前仅实现节点的短路故障)
			if isfield(self.fault, 'nid')
				k = find(self.ss.nodes.id == self.fault.nid);
				ek = zeros(length(self.ss.nodes.id), 1);
				ek(k) = 1;
				z1 = self.NAM1.get() \ ek;
			end

			% 计算附加阻抗与电流系数
			[za, k2, k0] = self.calcSumInfo('f3');

			% 计算短路点正序电压与全网正序电压
			Uss = self.ss.nodes.mag.*exp(self.ss.nodes.ang.*1i);	% 稳态
			If1 = Uss(k) ./ (z1(k) + za);	% z1(k) 是自阻抗, 在这里表示正序阻抗

			self.nodes.U1 = Uss - If1 .* z1;

			% 计算全网正序电流
			self.calcNetCurrentPS(self.nodes.U1);

			% 计算全网电压与电流
			self.branches.Ia = self.branches.I1;
			self.branches.Ib = self.branches.I1.*exp(-2./3.*pi.*1i);
			self.branches.Ic = self.branches.I1.*exp(2./3.*pi.*1i);
			self.nodes.Ua = self.nodes.U1;
			self.nodes.Ub = self.nodes.U1.*exp(-2./3.*pi.*1i);
			self.nodes.Uc = self.nodes.U1.*exp(2./3.*pi.*1i);
		end

		%% solveTwoPhaseShortCircult: 两相短路
		function solveTwoPhaseShortCircult(self)
			% 两相短路, 无需考虑零序
			self.NAM1.init(length(self.ss.nodes.id));
			self.NAM2.init(length(self.ss.nodes.id));
			self.NAMInit('12');

			%　求解节点阻抗矩阵的第　ｋ　列 (目前仅实现节点的短路故障)
			if isfield(self.fault, 'nid')
				k = find(self.ss.nodes.id == self.fault.nid);
				ek = zeros(length(self.ss.nodes.id), 1);
				ek(k) = 1;
				z1 = self.NAM1.get() \ ek;
				z2 = self.NAM2.get() \ ek;
			end

			% 计算附加阻抗与电流系数
			[za, k2, k0] = self.calcSumInfo('f2', z2(k));

			% 计算短路点正序电压与全网正序电压
			Uss = self.ss.nodes.mag.*exp(self.ss.nodes.ang.*1i);	% 稳态

			If1 = Uss(k) ./ (z1(k) + za);	% z1(k) 是自阻抗, 在这里表示正序阻抗
			If2 = k2.*If1;

			Uf1 = If1 .* za;
			Uf2 = -If2 .* z2;

			self.nodes.U1 = Uss - z1 .* If1;
			self.nodes.U2 = - z2 .* If2;

			% 计算全网正序电流
			self.calcNetCurrentPS(self.nodes.U1);
			self.calcNetCurrentNS(self.nodes.U2);

			% 计算全网电压与电流
			I = self.T * [self.branches.I1, self.branches.I2, zeros(length(self.branches.I1), 1)]';	% 3 * n
			self.branches.Ia = I(1, :)';
			self.branches.Ib = I(2, :)';
			self.branches.Ic = I(3, :)';
			U = self.T * [self.nodes.U1, self.nodes.U2, zeros(length(self.nodes.U1), 1)]';	% 3 * n
			self.nodes.Ua = U(1, :)';
			self.nodes.Ub = U(2, :)';
			self.nodes.Uc = U(3, :)';
		end

		%% solveSinglePhaseShortCircult: 单相短路
		function solveSinglePhaseShortCircult(self)
			% 单相短路
			self.NAM1.init(length(self.ss.nodes.id));
			self.NAM2.init(length(self.ss.nodes.id));
			self.NAM0.init(length(self.ss.nodes.id));
			self.NAMInit('120');

			%　求解节点阻抗矩阵的第　ｋ　列 (目前仅实现节点的短路故障)
			if isfield(self.fault, 'nid')
				k = find(self.ss.nodes.id == self.fault.nid);
				ek = zeros(length(self.ss.nodes.id), 1);
				ek(k) = 1;
				z1 = self.NAM1.get() \ ek;
				z2 = self.NAM2.get() \ ek;
				z0 = self.NAM0.get() \ ek;
				% z0 = pinv(full(self.NAM0.get())) * ek;	% 零序网络有可能不能连通, 这里使用广义逆
			end

			% debug
			% disp(z1);
			% disp(z2);
			% disp(z0);

			% 计算附加阻抗与电流系数
			[za, k2, k0] = self.calcSumInfo('f1', z2(k), z0(k));

			% 计算短路点正序电压与全网正序电压
			Uss = self.ss.nodes.mag.*exp(self.ss.nodes.ang.*1i);	% 稳态

			If1 = Uss(k) ./ (z1(k) + za);	% z1(k) 是自阻抗, 在这里表示正序阻抗
			If2 = k2.*If1;
			If0 = k0.*If1;

			% Uf1 = If1 .* za;
			% Uf2 = -If2 .* z2;
			% Uf0 = -If0 .* z0;

			self.nodes.U1 = Uss - z1 .* If1;
			self.nodes.U2 = - z2 .* If2;
			self.nodes.U0 = - z0 .* If0;

			% 计算全网正序电流
			self.calcNetCurrentPS(self.nodes.U1);
			self.calcNetCurrentNS(self.nodes.U2);
			self.calcNetCurrentZS(self.nodes.U0);

			% 计算全网电压与电流
			I = self.T * [self.branches.I1, self.branches.I2, self.branches.I0]';	% 3 * n
			self.branches.Ia = I(1, :)';
			self.branches.Ib = I(2, :)';
			self.branches.Ic = I(3, :)';
			U = self.T * [self.nodes.U1, self.nodes.U2, self.nodes.U0]';	% 3 * n
			self.nodes.Ua = U(1, :)';
			self.nodes.Ub = U(2, :)';
			self.nodes.Uc = U(3, :)';
		end

		%% solveTwoPhaseShortCircultToGround: 两相短路接地
		function solveTwoPhaseShortCircultToGround(self)
			% 两相短路接地
			self.NAM1.init(length(self.ss.nodes.id));
			self.NAM2.init(length(self.ss.nodes.id));
			self.NAM0.init(length(self.ss.nodes.id));
			self.NAMInit('120');

			%　求解节点阻抗矩阵的第　ｋ　列 (目前仅实现节点的短路故障)
			if isfield(self.fault, 'nid')
				k = find(self.ss.nodes.id == self.fault.nid);
				ek = zeros(length(self.ss.nodes.id), 1);
				ek(k) = 1;
				z1 = self.NAM1.get() \ ek;
				z2 = self.NAM2.get() \ ek;
				z0 = self.NAM0.get() \ ek;
			end

			% 计算附加阻抗与电流系数
			[za, k2, k0] = self.calcSumInfo('f11', z2(k), z0(k));

			% 计算短路点正序电压与全网正序电压
			Uss = self.ss.nodes.mag.*exp(self.ss.nodes.ang.*1i);	% 稳态

			If1 = Uss(k) ./ (z1(k) + za);	% z1(k) 是自阻抗, 在这里表示正序阻抗
			If2 = k2.*If1;
			If0 = k0.*If1;

			Uf1 = If1 .* za;
			% Uf2 = -If2 .* z2;
			% Uf0 = -If0 .* z0;

			self.nodes.U1 = Uss - z1 .* If1;
			self.nodes.U2 = - z2 .* If2;
			self.nodes.U0 = - z0 .* If0;

			% 计算全网正序电流
			self.calcNetCurrentPS(self.nodes.U1);
			self.calcNetCurrentNS(self.nodes.U2);
			self.calcNetCurrentZS(self.nodes.U0);

			% 计算全网电压与电流
			I = self.T * [self.branches.I1, self.branches.I2, self.branches.I0]';	% 3 * n
			self.branches.Ia = I(1, :)';
			self.branches.Ib = I(2, :)';
			self.branches.Ic = I(3, :)';
			U = self.T * [self.nodes.U1, self.nodes.U2, self.nodes.U0]';	% 3 * n
			self.nodes.Ua = U(1, :)';
			self.nodes.Ub = U(2, :)';
			self.nodes.Uc = U(3, :)';
		end

		%% solveFault: 电力系统故障分析
		function solveFault(self, solver, mpcSteady)
			
			% 潮流计算
			self.getPowerFlowResult(mpcSteady);

			% 设置基本相
			self.basePhase = self.getBasePhase();
			self.set120MatrixByBasePhase(self.basePhase);

			% 确定故障信息并求解
			switch self.fault.type
				case 'f3'
					self.solveThreePhaseShortCircult();
				case 'f2'
					self.solveTwoPhaseShortCircult();
				case 'f1'
					self.solveSinglePhaseShortCircult();
				case 'f11'
					self.solveTwoPhaseShortCircultToGround();
				otherwise
					error('illegal fault type');
			end
		end

	end

end