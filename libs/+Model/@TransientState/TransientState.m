%% TransientState
classdef TransientState < handle
	properties

		ss;		% 电力系统稳态模型
		ft;		% 电力系统故障模型
		
		bus;		% 
		gen;	% 
		branch;	% 

		dt;			% 大干扰稳定分析时用到的步长, 0.001s 左右
		ct;			% 大干扰稳定分析时用到的离散时间, 以向量形式存储
		time;		% 
		operating;	% 电力系统发生状态变化的信息
		NAM;

		itlog;
		vot;	% 存储暂态过程中各节点的电压
		cur;	% 存储暂态过程中各节点的电流
		delta;	% 存储暂态过程中各发电机的功角
		omega;	% 存储暂态过程中各发电机的转速

	end

	methods

		%% init: 电力系统暂态分析模型属性初始化
		function init(self, mpc)

			self.operating = mpc.operating;
			self.NAM = Model.NAM();

			% 目前我们只用发电机暂态参数
			self.gen.nid = mpc.gen(:, 1);
			gensLength = length(self.gen.nid);
			self.gen.id = (1:gensLength)';

			self.gen.TJ = mpc.gen(:, 2);			% 转子时间常数
			self.gen.D = mpc.gen(:, 10);	% 阻尼

			% self.gen.order = 2;
			self.gen.Pm = zeros(gensLength, 1);		% 机械功率
			self.gen.Pe = zeros(gensLength, 1);		% 电磁功率
			self.gen.Udq = zeros(gensLength, 1);		% 定子电压(电机坐标系下)
			self.gen.Uxy = zeros(gensLength, 1);		% 定子电压(电网坐标系下)
			self.gen.Idq = zeros(gensLength, 1);		% 定子电流(电机坐标系下)
			self.gen.Ixy = zeros(gensLength, 1);		% 定子电流(电网坐标系下)
			% self.gen.Ux = zeros(gensLength, 1);		% 定子电压(电网坐标系下)
			% self.gen.Ix = zeros(gensLength, 1);		% 定子电流(电网坐标系下)
			% self.gen.Uy = zeros(gensLength, 1);		% 定子电压(电网坐标系下)
			% self.gen.Iy = zeros(gensLength, 1);		% 定子电流(电网坐标系下)
			
			self.gen.Ra = mpc.gen(:, 3);			% 电机参数, 电枢电阻
			self.gen.Xd = mpc.gen(:, 4);			% 电机参数, 直轴同步电抗
			self.gen.Xd1 = mpc.gen(:, 5);			% 电机参数, 直轴暂态电抗
			self.gen.Xd2 = zeros(gensLength, 1);	% TODO 电机参数, 直轴次暂态电抗
			self.gen.Xq = mpc.gen(:, 6);			% 电机参数, 交轴同步电抗
			self.gen.Xq1 = mpc.gen(:, 7);			% 电机参数, 交轴暂态电抗
			self.gen.Xq2 = zeros(gensLength, 1);	% TODO 电机参数, 交轴次暂态电抗

			self.gen.EQ = zeros(gensLength, 1);	% 虚构电势, 确定功角用
			self.gen.Eq = zeros(gensLength, 1);		% 交轴励磁电势
			self.gen.Eq1 = zeros(gensLength, 1);	% 交轴暂态电势
			self.gen.Eq2 = zeros(gensLength, 1);	% TODO 交轴次暂态电势
			self.gen.Ed = zeros(gensLength, 1);		% 直轴励磁电势
			self.gen.Ed1 = zeros(gensLength, 1);	% 直轴暂态电势
			self.gen.Ed2 = zeros(gensLength, 1);	% TODO 直轴次暂态电势

			self.gen.omegaBase = 376.8;
			self.gen.omega = ones(gensLength, 1);	% 转速标幺值
			self.gen.delta = zeros(gensLength, 1);	% 这个功角是发电机相对出口的功角, 不是相对平衡节点的

			self.gen.Td01 = mpc.gen(:, 8);			% 电机参数, 直轴开路暂态时间常数
			self.gen.Td02 = zeros(gensLength, 1);	% TODO 电机参数, 直轴开路次暂态时间常数
			self.gen.Tq01 = mpc.gen(:, 9);			% 电机参数, 交轴开路暂态时间常数
			self.gen.Tq02 = zeros(gensLength, 1);	% TODO 电机参数, 交轴开路次暂态时间常数

			self.delta = [];
			self.omega = [];
			self.vot = [];
			self.cur = [];

			self.itlog.Pe = [];
			self.itlog.Pm = [];

			self.itlog.bxgy = [];

			% 线路和节点的暂态参数列表
			self.bus.id = [];
			self.bus.trY = [];
			self.branch.fid = [];
			self.branch.tid = [];
			self.branch.status = [];	% 0 断开, 1 闭合

		end

		%% setTimeAndStep: 设置时间与步长(时间)
		function setTimeAndStep(self, time, dt)
			self.dt = dt;
			self.ct = 0:dt:time;
		end

		%% getPowerFlowResult: 计算大干扰前电网的潮流
		function getPowerFlowResult(self, mpcSteady)

			% 建立电力网稳态模型并初始化
			self.ss = Model.SteadyState();
			self.ss.init(mpcSteady);

			%% 设置求解器的基本信息
			solver.method = 'NR';	% 求解方法
			solver.n_iters_max = 6;	% 最大迭代
			solver.epsilon = 1e-5;	% 收敛判据, 功率不平衡量标幺
			solver.start = 'default';	% 启动方式, default 为按发电机端电压起动

			%% 求解
			self.ss.solvePowerFlow(solver);
		end

		%% setGeneratorCurrent: 计算发电机注入网络的电流
		function setGeneratorCurrent(self)
			index = getIndex(self.ss.bus.id, self.ss.gen.nid);
			Sxy_conj = self.ss.bus.Pg(index) - self.ss.bus.Qg(index) * 1i;
			self.gen.Uxy = self.ss.bus.mag(index).*exp(self.ss.bus.ang(index).*1i);
			self.gen.Ixy = Sxy_conj./conj(self.gen.Uxy);
		end

		%% setDummyPotential: 计算发电机虚构电势
		function setDummyPotential(self)
			self.gen.EQ = self.gen.Uxy + (self.gen.Ra + self.gen.Xq*1i).*self.gen.Ixy;
		end

		%% powerAngleInit: 计算发电机转子角度初始值
		function powerAngleInit(self)
			self.gen.delta = angle(self.gen.EQ);
		end

		%% setStatorStatus: 电压电流的机网转换
		function [Udq, Idq] = setStatorStatus(self, Uxy, Ixy, delta)
			% 下面的 self.gen.delta 是电机当前的功角, 注意这个接口只用于初值计算
			Udq = Uxy.*exp((pi/2 - delta).*1i);
			Idq = Ixy.*exp((pi/2 - delta).*1i);
		end

		%% setTransientPotential: 计算暂态电势
		function setTransientPotential(self, Udq, Idq)
			% 注意这个接口只用于初值计算
			self.gen.Eq = imag(Udq) + self.gen.Ra .* imag(Idq) + self.gen.Xd .* real(Idq);
			self.gen.Ed = real(Udq) + self.gen.Ra .* real(Idq) - self.gen.Xq .* imag(Idq);
			self.gen.Eq1 = imag(Udq) + self.gen.Ra .* imag(Idq) + self.gen.Xd1 .* real(Idq);
			self.gen.Ed1 = real(Udq) + self.gen.Ra .* real(Idq) - self.gen.Xq1 .* imag(Idq);
		end

		%% updateMagForce: Eq1, Ed1
		% function updateMagForce(self)
		% 	self.gen.Eq1 = imag(self.gen.Udq) + self.gen.Ra .* imag(self.gen.Idq) + self.gen.Xd1 .* real(self.gen.Idq);
		% 	self.gen.Ed1 = real(self.gen.Udq) + self.gen.Ra .* real(self.gen.Idq) - self.gen.Xq1 .* imag(self.gen.Idq);
		% end

		%% setMechancialPower: 计算初始电磁功率
		function setMechancialPower(self, Udq, Idq)
			nid = self.ss.gen.nid;

			self.gen.Pe = self.ss.bus.Pg(getIndex(self.ss.bus.id, nid)) + abs(self.gen.Ixy).^2 .* self.gen.Ra;
			self.gen.Pm = self.gen.Pe;
		end

		%% net2Machine: 机网转换
		function net2Machine(self)
			self.gen.Udq = self.gen.Uxy.*exp((pi/2 - self.gen.delta).*1i);
			self.gen.Idq = self.gen.Ixy.*exp((pi/2 - self.gen.delta).*1i);
		end

		%% getGeneratorAdmittance: 返回发电机在暂态过程中等效的导纳
		function [YG1] = getGeneratorAdmittance(self)
			den = self.gen.Ra.^2 + self.gen.Xd1.*self.gen.Xq1;
			G = self.gen.Ra ./ den;
			B = - 0.5.*(self.gen.Xd1 + self.gen.Xq1) ./ den;
			YG1 = G + B.*1i;
		end

		%% addGenImpedanceToBus: 迭代法计算网络方程，发电机入阵
		function addGenImpedanceToBus(self)
			YG1 = self.getGeneratorAdmittance();
			self.ss.NAM.addAdmittance(getIndex(self.ss.bus.id, self.gen.nid), YG1);
		end

		%% setTransientInitialValue: 大干扰稳定分析初值计算
		function setTransientInitialValue(self, mpcSteady)
			
			% 1. 电力系统潮流计算
			self.getPowerFlowResult(mpcSteady);

			% 2. 计算发电机注入网络的电流
			% 3. 计算发电机虚构电势
			% 4. 计算发电机转子角度初始值
			% 5. 电压电流的机网转换
			self.setGeneratorCurrent();
			self.setDummyPotential();
			self.powerAngleInit();
			[Udq, Idq] = self.setStatorStatus(self.gen.Uxy, self.gen.Ixy, self.gen.delta);	% 我们默认稳态运行时, 电角速度标幺值为 1

			% 6. 计算暂态电势
			% 7. 计算初始电磁功率
			self.setTransientPotential(Udq, Idq);
			self.setMechancialPower(Udq, Idq);

			% 8. 电网稳态模型中负荷按恒阻抗归入节点导纳矩阵

			self.ss.addLoadImpedanceToBus();
			% self.ss.addGenImpedanceToBus(self.gen);	% 发电机按经典模型计算时可以在这里归入节点导纳矩阵
			% self.addGenImpedanceToBus();
			% self.ss.addMoterImpedanceToBus();

			% debug 通过
			% disp('发电机内电势, 功角: ');
			% disp([self.gen.Eq, self.gen.Ed, self.gen.Eq1, self.gen.Ed1, self.gen.delta*180/pi]);
			% disp('发电机电压, 电流(网)')
			% disp([abs(self.gen.Uxy), angle(self.gen.Uxy)*180/pi]);
			% disp([abs(self.gen.Ixy), angle(self.gen.Ixy)*180/pi]);
			% disp('发电机电压, 电流(机)')
			% disp([abs(Udq), angle(Udq)*180/pi]);
			% disp([abs(Idq), angle(Idq)*180/pi]);
			% disp('机械功率: ')
			% disp([self.gen.Pm]);

		end
		
		%% setData: 记录计算过程中得到的节点电压, 线路电流, 发电机功角
		function setData(self, Uxy, Ixy)
			self.delta = [self.delta, self.gen.delta];
			self.omega = [self.omega, self.gen.omega];

			% debug
			self.itlog.Pe = [self.itlog.Pe, self.gen.Pe];
			self.itlog.Pm = [self.itlog.Pm, self.gen.Pm];

			self.vot = [self.vot, Uxy];
			self.cur = [self.cur, Ixy];
		end

		%% getStage: 返回当前时刻电网所处的阶段, 对应各个节点导纳矩阵
		function [stage] = getStage(self, time)
			if isempty(self.operating)
				stage = 0;
			else
				stage = [self.operating.time];
				stage = length(find(time >= stage));
			end
		end

		%% addTransientImpedanceToNodes: 向节点添加过渡电阻
		function addTransientImpedanceToNodes(self, id, trZ)
			index = find(self.bus.id == id);
			if isempty(index)
				index = length(self.bus.id) + 1;
				self.bus.id(index) = id;
				self.bus.trY(index) = 1./trZ;
			else
				self.bus.trY(index) = self.bus.trY(index) + 1./trZ;
			end
			disp('Transient impedance added');
		end

		%% removeTransientImpedanceFromNodes: 移除某个节点的过渡电阻
		function removeTransientImpedanceFromNodes(self, id)
			index = find(self.bus.id == id);
			assert(~isempty(index));

			self.bus.trY(index) = 0;	% 移除
			disp('Transient impedance removed');
		end

		%% breakBranch: 断开线路
		function breakBranch(self, fid, tid)

			index = intersect(find(self.branch.fid == fid), find(self.branch.tid == tid)) ;
			if(isempty(index))
				index = length(self.branch.fid) + 1;
				self.branch.fid(index) = fid;
				self.branch.tid(index) = tid;
			end
			% 直接断开
			self.branch.status(index) = 0;
			disp('Branch breaked');
		end

		%% restoreBranch: 重合闸线路
		function restoreBranch(self, fid, tid)
			index = find(fid == self.branch.fid && tid == self.branch.tid);
			assert(~isempty(index));
			% 直接连接
			self.branch.status(index) = 1;
			disp('Branch restored');
		end

		%% setTransientData: 设置该阶段的暂态信息
		function setTransientData(self, stage, solver)

			% 检查是否有与节点有关的操作
			if isfield(self.operating(stage), 'nid') && ~isempty(self.operating(stage).nid)
				switch self.operating(stage).ntype
					case 'f3'	% 三相短路
						% self.addTransientImpedanceToNodes(self.operating(stage).nid, self.operating(stage).zf);
						self.ft.fault.zf = self.operating(stage).zf;
						[za, k2, k0] = self.ft.calcSumInfo(self.operating(stage).ntype);
						self.addTransientImpedanceToNodes(self.operating(stage).nid, za);
					case 'f2'	% 两相短路
						self.ft.fault.zf = self.operating(stage).zf;
						self.ft.NAM2.init(length(self.ss.bus.id));
						self.ft.NAMInit('2');
						%　(目前仅实现节点的短路故障)
						k = find(self.ss.bus.id == self.operating(stage).nid);
						ek = zeros(length(self.ss.bus.id), 1);
						ek(k) = 1;
						z2 = self.ft.NAM2.get() \ ek;
						[za, k2, k0] = self.ft.calcSumInfo(self.operating(stage).ntype, z2(k));
						self.addTransientImpedanceToNodes(self.operating(stage).nid, za);
					case 'f1'	% 单相短路
						self.ft.fault.zf = self.operating(stage).zf;
						self.ft.NAM2.init(length(self.ss.bus.id));
						self.ft.NAM0.init(length(self.ss.bus.id));
						self.ft.NAMInit('20');
						%　(目前仅实现节点的短路故障)
						k = find(self.ss.bus.id == self.operating(stage).nid);
						ek = zeros(length(self.ss.bus.id), 1);
						ek(k) = 1;
						z2 = self.ft.NAM2.get() \ ek;
						z0 = self.ft.NAM0.get() \ ek;
						[za, k2, k0] = self.ft.calcSumInfo(self.operating(stage).ntype, z2(k), z0(k));
						self.addTransientImpedanceToNodes(self.operating(stage).nid, za);
					case 'f11'	% 两相短路接地
						self.ft.fault.zf = self.operating(stage).zf;
						self.ft.fault.zg = self.operating(stage).zg;
						self.ft.NAM2.init(length(self.ss.bus.id));
						self.ft.NAM0.init(length(self.ss.bus.id));
						self.ft.NAMInit('20');
						%　(目前仅实现节点的短路故障)
						k = find(self.ss.bus.id == self.operating(stage).nid);
						ek = zeros(length(self.ss.bus.id), 1);
						ek(k) = 1;
						z2 = self.ft.NAM2.get() \ ek;
						z0 = self.ft.NAM0.get() \ ek;
						[za, k2, k0] = self.ft.calcSumInfo(self.operating(stage).ntype, z2(k), z0(k));
						self.addTransientImpedanceToNodes(self.operating(stage).nid, za);
					case 'restore'	% 解除短路状态
						self.removeTransientImpedanceFromNodes(self.operating(stage).nid);
					otherwise
						error('illegal node operation type');
				end
			end
			% 检查是否有与线路有关的操作
			if isfield(self.operating(stage), 'fid') && ~isempty(self.operating(stage).fid)
				switch self.operating(stage).btype
					case 'break'
						self.breakBranch(self.operating(stage).fid, self.operating(stage).tid);
						disp((self.delta(1, end) - self.delta(2, end))*180/pi);
					case 'restore'
						self.restoreBranch(self.operating(stage).fid, self.operating(stage).tid);
					otherwise
						error('illegal branch operation type');
				end
			end
		end

		%% setTransientNAM: 设置该阶段的节点导纳矩阵
		function setTransientNAM(self)
			% 这个算法的时间复杂度不乐观
			self.NAM.set(self.ss.NAM);

			% 短路节点的处理
			index = getIndex(self.ss.bus.id, self.bus.id);
			self.NAM.addAdmittance(index, self.bus.trY);
			
			% 线路的处理
			for k = 1:length(self.branch.fid)
				index = intersect(find(self.ss.branch.fid == self.branch.fid(k)), find(self.ss.branch.tid == self.branch.tid(k)));
				if self.branch.status == 0
					self.NAM.removeLine(self.branch.fid(k), self.branch.tid(k), self.ss.branch.r(index(1)) + self.ss.branch.x(index(1))*1i, self.ss.branch.g(index(1)) + self.ss.branch.b(index(1))*1i);
				end
			end

		end

		%% param_gb: 网络方程求解需要用到的参数
		function [array] = param_gb(self, delta)
			den = self.gen.Ra.^2 + self.gen.Xd1.*self.gen.Xq1;
			gx = (self.gen.Ra.*sin(delta) - self.gen.Xd1.*cos(delta))./den;
			gy = (self.gen.Ra.*sin(delta) - self.gen.Xq1.*cos(delta))./den;
			bx = (self.gen.Ra.*cos(delta) + self.gen.Xq1.*sin(delta))./den;
			by = (-self.gen.Ra.*cos(delta) - self.gen.Xd1.*sin(delta))./den;
			array = [gx, bx, by, gy];
		end

		% param_GB: 网络方程求解需要用到的参数
		function [array] = param_GB(self, delta)
			den = self.gen.Ra.^2 + self.gen.Xd1.*self.gen.Xq1;
			Gx = (self.gen.Ra - 0.5.*(self.gen.Xd1-self.gen.Xq1).*sin(2.*delta))./den;
			Gy = (self.gen.Ra + 0.5.*(self.gen.Xd1-self.gen.Xq1).*sin(2.*delta))./den;
			Bx = (0.5.*(self.gen.Xd1+self.gen.Xq1) + 0.5.*(self.gen.Xd1-self.gen.Xq1).*cos(2.*delta))./den;
			By = (-0.5.*(self.gen.Xd1+self.gen.Xq1) + 0.5.*(self.gen.Xd1-self.gen.Xq1).*cos(2.*delta))./den;
			array = [Gx, Bx, By, Gy];
		end

		%% calcElectricMagneticPower: 计算电磁功率
		function [Pe] = calcElectricMagneticPower(self, Uxy, Ixy)
			% index = getIndex(self.ss.bus.id, self.gen.nid);
			Pe = real(Uxy .* conj(Ixy)) + abs(Ixy).^2 .* self.gen.Ra;
		end

		%% stepForward: 求解微分代数方程
		function [omega1, delta1] = stepForward(self, Uxy, Ixy)
			
			self.gen.Pe = self.calcElectricMagneticPower(Uxy, Ixy);
			omega1 = (self.gen.Pm - self.gen.Pe - (self.gen.omega-1).*self.gen.D) ./ self.gen.TJ;
			% self.gen.omega = self.gen.omega + omega1 .* self.dt;
			delta1 = self.gen.omegaBase .* (self.gen.omega + omega1 .* self.dt./2 - 1);
			% self.gen.delta = self.gen.delta + delta1 .* self.dt;
		end

		%% getBlockForNAM: 向节点导纳矩阵实矩阵中叠加一个分块
		function [blk] =  getBlockForNAM(self, pos, param_GB)
			odd = 2.*pos - 1;
			even = 2.*pos;
			blk = sparse(2.*length(self.ss.bus.id), 2.*length(self.ss.bus.id));
			blk(odd, odd) = diag(param_GB(:, 1));
			blk(odd, even) = diag(param_GB(:, 2));
			blk(even, odd) = diag(param_GB(:, 3));
			blk(even, even) = diag(param_GB(:, 4));
		end

		%% solveNetEquation: 求解网络方程
		function [Uxy, Ixy] = solveNetEquation(self, net, delta)
			index = getIndex(self.ss.bus.id, self.gen.nid);
			n = self.NAM.size();
			odd = (2.*(1:n) - 1)';
			even = (2.*(1:n))';
			param_GB = self.param_GB(delta);

			% 计算暂态电动势与注入电流
			I = zeros(2.*n);
			
			% 二阶经典
			% Exy = (self.gen.Ed1 + self.gen.Eq1.*1i).*exp(-(pi./2 - delta).*1i);
			% Eq1 恒定
			Exy = (self.gen.Eq1.*1i).*exp(-(pi./2 - delta).*1i);

			I(odd(index)) = param_GB(:, 1).*real(Exy) + param_GB(:, 2).*imag(Exy);
			I(even(index)) = param_GB(:, 3).*real(Exy) + param_GB(:, 4).*imag(Exy);

			% 根据发电机的运行参数，修改导纳矩阵
			blk = self.getBlockForNAM(index, param_GB);

			% 求解节点电压方程
			U = (self.NAM.real + blk) \ I;
			Uxy = U(odd) + U(even).*1i;

			Ut = Exy - Uxy(index);
			I(odd(index)) = param_GB(:, 1).*real(Ut) + param_GB(:, 2).*imag(Ut);
			I(even(index)) = param_GB(:, 3).*real(Ut) + param_GB(:, 4).*imag(Ut);
			Ixy = I(odd) + I(even).*1i;

		end

		%% solveStep: 求解一步
		function [stable] = solveStep(self, dae, net, t)
			
			index = getIndex(self.ss.bus.id, self.gen.nid);
			switch dae
				case {'euler', 'Euler'}
					% [Uxy, Ixy] = self.solveNetEquation(net, self.gen.delta);
					[omega1, delta1] = self.stepForward(self.gen.Uxy, self.gen.Ixy);
					[Uxy, Ixy] = self.solveNetEquation(net, self.gen.delta + delta1.*self.dt);
					[omega2, delta2] = self.stepForward(Uxy(index), Ixy(index));

					self.gen.omega = self.gen.omega + (omega1 + omega2)./2.*self.dt;
					self.gen.delta = self.gen.delta + (delta1 + delta2)./2.*self.dt;
					[Uxy, Ixy] = self.solveNetEquation(net, self.gen.delta);
					self.gen.Uxy = Uxy(index);
					self.gen.Ixy = Ixy(index);
				otherwise
					error('illegal dae method');
			end

			self.setData(Uxy, Ixy);

			if ~self.isStable()
				warning(['Transient instability occurs at ', num2str(t), 's.']);
				stable = false;
				return;
			end
			stable = true;
		end

		%% checkStage: 根据当前时间选择性更新网络状态
		function [stage] = checkStage(self, t, stagePrev, solver)
			self.time = t;
			stage = self.getStage(t);
			% 系统状态发生了变化
			if stage ~= stagePrev
				self.setTransientData(stage, solver);
				self.setTransientNAM();
				self.NAM.setReal();	% 不用直接法求网络方程可以去掉
				% 求解网络方程
				[Uxy, Ixy] = self.solveNetEquation(solver.net, self.gen.delta);
				self.gen.Uxy = Uxy(getIndex(self.ss.bus.id, self.gen.nid));
				self.gen.Ixy = Ixy(getIndex(self.ss.bus.id, self.gen.nid));
			end
		end

		%% isStable: 判断电力系统当前是否稳定
		function [stable] = isStable(self)
			stable = range(self.gen.delta) < 1.6.*pi;
		end

		%% solveLargeDisturbance: 电力系统大干扰分析
		function solveLargeDisturbance(self, solver, mpcSteady, mpcFault)

			% 1. 初值计算
			self.setTransientInitialValue(mpcSteady);
			% 判断是否有不对称故障，选择性初始化故障模型
			% if any(strcmp(self.operating(1).ntype, {'f2', 'f1', 'f11'}))
				self.ft = Model.Fault();
				self.ft.init(mpcFault);
				self.ft.ss = self.ss;
			% end

			% 2. 
			self.setTimeAndStep(solver.time, solver.dt);

			% 3. 网络方程的求解, 微分代数方程的数值积分, 系统状态改变时更新相应的节点导纳矩阵
			stagePrev = 0;
			for t = self.ct
				stagePrev = self.checkStage(t, stagePrev, solver);
				stable = self.solveStep(solver.dae, solver.net, t);
				
				if ~stable
					return;
				end
			end

		end

	end
end
