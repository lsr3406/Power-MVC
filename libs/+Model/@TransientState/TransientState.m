%% TransientState
classdef TransientState < handle
	properties

		ss;		% 电力系统稳态模型
		
		nodes;		% 
		generator;	% 
		branches;	% 

		dt;			% 大干扰稳定分析时用到的步长, 0.001s 左右
		ct;			% 大干扰稳定分析时用到的离散时间, 以向量形式存储
		time;		% 
		operating;	% 电力系统发生状态变化的信息
		NAM;		% 

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
			self.generator.nid = mpc.gen(:, 1);
			gensLength = length(self.generator.nid);
			self.generator.id = (1:gensLength)';

			self.generator.TJ = mpc.gen(:, 2);			% 转子时间常数
			self.generator.D = zeros(gensLength, 1);	% 阻尼

			self.generator.Pm = zeros(gensLength, 1);		% 机械功率
			self.generator.Pe = zeros(gensLength, 1);		% 电磁功率
			self.generator.Udq = zeros(gensLength, 1);		% 定子电压(电机坐标系下)
			self.generator.Uxy = zeros(gensLength, 1);		% 定子电压(电网坐标系下)
			self.generator.Idq = zeros(gensLength, 1);		% 定子电流(电机坐标系下)
			self.generator.Ixy = zeros(gensLength, 1);		% 定子电流(电网坐标系下)
			% self.generator.Ux = zeros(gensLength, 1);		% 定子电压(电网坐标系下)
			% self.generator.Ix = zeros(gensLength, 1);		% 定子电流(电网坐标系下)
			% self.generator.Uy = zeros(gensLength, 1);		% 定子电压(电网坐标系下)
			% self.generator.Iy = zeros(gensLength, 1);		% 定子电流(电网坐标系下)
			
			self.generator.Ra = mpc.gen(:, 3);			% 电机参数, 电枢电阻
			self.generator.Xd = mpc.gen(:, 4);			% 电机参数, 直轴同步电抗
			self.generator.Xd1 = mpc.gen(:, 5);			% 电机参数, 直轴暂态电抗
			self.generator.Xd2 = zeros(gensLength, 1);	% TODO 电机参数, 直轴次暂态电抗
			self.generator.Xq = mpc.gen(:, 6);			% 电机参数, 交轴同步电抗
			self.generator.Xq1 = mpc.gen(:, 7);			% 电机参数, 交轴暂态电抗
			self.generator.Xq2 = zeros(gensLength, 1);	% TODO 电机参数, 交轴次暂态电抗

			self.generator.EQ = zeros(gensLength, 1);	% 虚构电势, 确定功角用
			self.generator.Eq = zeros(gensLength, 1);		% 交轴励磁电势
			self.generator.Eq1 = zeros(gensLength, 1);	% 交轴暂态电势
			self.generator.Eq2 = zeros(gensLength, 1);	% TODO 交轴次暂态电势
			self.generator.Ed = zeros(gensLength, 1);		% 直轴励磁电势
			self.generator.Ed1 = zeros(gensLength, 1);	% 直轴暂态电势
			self.generator.Ed2 = zeros(gensLength, 1);	% TODO 直轴次暂态电势

			self.generator.omegaBase = 376.8;
			self.generator.omega = ones(gensLength, 1);	% 转速标幺值
			self.generator.delta = zeros(gensLength, 1);	% 这个功角是发电机相对出口的功角, 不是相对平衡节点的

			self.generator.Td01 = mpc.gen(:, 8);			% 电机参数, 直轴开路暂态时间常数
			self.generator.Td02 = zeros(gensLength, 1);	% TODO 电机参数, 直轴开路次暂态时间常数
			self.generator.Tq01 = mpc.gen(:, 9);			% 电机参数, 交轴开路暂态时间常数
			self.generator.Tq02 = zeros(gensLength, 1);	% TODO 电机参数, 交轴开路次暂态时间常数


			self.delta = [];
			self.omega = [];
			self.vot = [];
			self.cur = [];

			self.itlog.Pe = [];
			self.itlog.Pm = [];

			self.itlog.bxgy = [];

			% 线路和节点的暂态参数列表
			self.nodes.id = [];
			self.nodes.trY = [];
			self.branches.fid = [];
			self.branches.tid = [];
			self.branches.status = [];	% 0 断开, 1 闭合

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
			solver.method = 'FDBX';	% 求解方法
			solver.maxIteration = 30;	% 最大迭代
			solver.epsilon = 1e-5;	% 收敛判据, 功率不平衡量标幺
			solver.start = 'default';	% 启动方式, default 为按发电机端电压起动

			%% 求解
			self.ss.solvePowerFlow(solver);
		end

		%% setGeneratorCurrent: 计算发电机注入网络的电流
		function setGeneratorCurrent(self)
			index = getIndex(self.ss.nodes.id, self.ss.generator.nid);
			Sxy_conj = self.ss.nodes.Pg(index) - self.ss.nodes.Qg(index) * 1i;
			self.generator.Uxy = self.ss.nodes.mag(index).*exp(self.ss.nodes.ang(index).*1i);
			self.generator.Ixy = Sxy_conj./conj(self.generator.Uxy);
		end

		%% setImaginaryPotential: 计算发电机虚构电势
		function setImaginaryPotential(self)
			self.generator.EQ = self.generator.Uxy + (self.generator.Ra + self.generator.Xq*1i).*self.generator.Ixy;
		end

		%% powerAngleInit: 计算发电机转子角度初始值
		function powerAngleInit(self)
			self.generator.delta = angle(self.generator.EQ);
		end

		%% setStatorStatus: 电压电流的机网转换
		function [Udq, Idq] = setStatorStatus(self, Uxy, Ixy)
			% 下面的 self.generator.delta 是电机当前的功角, 注意这个接口只用于初值计算
			Udq = Uxy.*exp((pi/2 - self.generator.delta).*1i);
			Idq = Ixy.*exp((pi/2 - self.generator.delta).*1i);
		end

		%% setTransientPotential: 计算暂态电势
		function setTransientPotential(self, Udq, Idq)
			% 注意这个接口只用于初值计算
			self.generator.Eq = imag(Udq) + self.generator.Ra .* imag(Idq) + self.generator.Xd .* real(Idq);
			self.generator.Ed = real(Udq) + self.generator.Ra .* real(Idq) - self.generator.Xq .* imag(Idq);
			self.generator.Eq1 = imag(Udq) + self.generator.Ra .* imag(Idq) + self.generator.Xd1 .* real(Idq);
			self.generator.Ed1 = real(Udq) + self.generator.Ra .* real(Idq) - self.generator.Xq1 .* imag(Idq);
		end

		%% setMechancialPower: 计算初始电磁功率
		function setMechancialPower(self)
			nid = self.ss.generator.nid;

			self.generator.Pe = self.ss.nodes.Pg(getIndex(self.ss.nodes.id, nid)) + abs(self.generator.Ixy).^2 .* self.generator.Ra;
			self.generator.Pm = self.generator.Pe;
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
			self.setImaginaryPotential();
			self.powerAngleInit();
			[Udq, Idq] = self.setStatorStatus(self.generator.Uxy, self.generator.Ixy);	% 我们默认稳态运行时, 电角速度标幺值为 1

			% 6. 计算暂态电势
			% 7. 计算初始电磁功率
			self.setTransientPotential(Udq, Idq);
			self.setMechancialPower();

			% 8. 电网稳态模型中负荷按恒阻抗归入节点导纳矩阵

			self.ss.addLoadImpedanceToNodes();
			self.ss.addGenImpedanceToNodes(self.generator);	% 发电机按经典模型计算时可以在这里归入节点导纳矩阵
			% self.ss.addMoterImpedanceToNodes();

			% debug 通过
			% disp('发电机内电势, 功角: ');
			% disp([self.generator.Eq, self.generator.Ed, self.generator.Eq1, self.generator.Ed1, self.generator.delta*180/pi]);
			% disp('发电机电压, 电流(网)')
			% disp([abs(self.generator.Uxy), angle(self.generator.Uxy)*180/pi]);
			% disp([abs(self.generator.Ixy), angle(self.generator.Ixy)*180/pi]);
			% disp('发电机电压, 电流(机)')
			% disp([abs(Udq), angle(Udq)*180/pi]);
			% disp([abs(Idq), angle(Idq)*180/pi]);
			% disp('机械功率: ')
			% disp([self.generator.Pm]);

		end
		
		%% setData: 记录计算过程中得到的节点电压, 线路电流, 发电机功角
		function setData(self)
			self.delta = [self.delta, self.generator.delta];
			self.omega = [self.omega, self.generator.omega];

			% debug
			self.itlog.Pe = [self.itlog.Pe, self.generator.Pe];
			self.itlog.Pm = [self.itlog.Pm, self.generator.Pm];
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
			index = find(self.nodes.id == id);
			if isempty(index)
				index = length(self.nodes.id) + 1;
				self.nodes.id(index) = id;
				self.nodes.trY(index) = 1./trZ;
			else
				self.nodes.trY(index) = self.nodes.trY(index) + 1./trZ;
			end
			disp('Transient impedance added');
		end

		%% removeTransientImpedanceFromNodes: 移除某个节点的过渡电阻
		function removeTransientImpedanceFromNodes(self, id)
			index = find(self.nodes.id == id);
			assert(~isempty(index));

			self.nodes.trY(index) = 0;	% 移除
			disp('Transient impedance removed');
		end

		%% breakBranch: 断开线路
		function breakBranch(self, fid, tid)

			index = intersect(find(self.branches.fid == fid), find(self.branches.tid == tid)) ;
			if(isempty(index))
				index = length(self.branches.fid) + 1;
				self.branches.fid(index) = fid;
				self.branches.tid(index) = tid;
			end
			% 直接断开
			self.branches.status(index) = 0;
			disp('Branch breaked');
		end

		%% restoreBranch: 重合闸线路
		function restoreBranch(self, fid, tid)
			index = find(fid == self.branches.fid && tid == self.branches.tid);
			assert(~isempty(index));
			% 直接连接
			self.branches.status(index) = 1;
			disp('Branch restored');
		end

		%% setTransientData: 设置该阶段的暂态信息
		function setTransientData(self, stage, solver)

			% 检查是否有与节点有关的操作
			if isfield(self.operating(stage), 'nid') && ~isempty(self.operating(stage).nid)
				switch self.operating(stage).ntype
					case 'f3'	% 三相短路
						self.addTransientImpedanceToNodes(self.operating(stage).nid, self.operating(stage).trZ);
					case 'f2'	% 两相短路
						% TODO 两相短路
					case 'f1'	% 单相短路
						% TODO 单相短路
					case 'f11'	% 两相短路接地
						% TODO 两相短路接地
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
					case 'restore'
						self.restoreBranch(self.operating(stage).fid, self.operating(stage).tid);
					otherwise
						error('illegal branches operation type');
				end
			end
		end

		%% setTransientNAM: 设置该阶段的节点导纳矩阵
		function setTransientNAM(self)
			% 这个算法的时间复杂度不乐观
			self.NAM.set(self.ss.NAM);

			% 短路节点的处理
			index = getIndex(self.ss.nodes.id, self.nodes.id);
			self.NAM.addAdmittance(index, self.nodes.trY);
			
			% 线路的处理
			for k = 1:length(self.branches.fid)
				index = intersect(find(self.ss.branches.fid == self.branches.fid(k)), find(self.ss.branches.tid == self.branches.tid(k)));
				if self.branches.status == 0
					self.NAM.removeLine(self.branches.fid(k), self.branches.tid(k), self.ss.branches.r(index) + self.ss.branches.x(index)*1i, self.ss.branches.g(index) + self.ss.branches.b(index)*1i);
				end
			end
			
		end

		%% param_gb: 网络方程求解需要用到的参数
		function [gxby, bxgy] = param_gb(self)
			den = self.generator.Ra.^2 + self.generator.Xd1.*self.generator.Xq;
			gx = (self.generator.Ra.*sin(self.generator.delta) - self.generator.Xd1.*cos(self.generator.delta))./den;
			gy = (self.generator.Ra.*sin(self.generator.delta) - self.generator.Xq.*cos(self.generator.delta))./den;
			bx = (self.generator.Ra.*cos(self.generator.delta) + self.generator.Xq.*sin(self.generator.delta))./den;
			by = (-self.generator.Ra.*cos(self.generator.delta) - self.generator.Xd1.*sin(self.generator.delta))./den;
			gxby = gx + by .* 1i;
			bxgy = bx + gy .* 1i;
		end

		%% param_GB: 网络方程求解需要用到的参数
		function [GxBy, BxGy] = param_GB(self)
			den = self.generator.Ra.^2 + self.generator.Xd1.*self.generator.Xq;
			Gx = (self.generator.Ra - 0.5.*(self.generator.Xd1-self.generator.Xq).*sin(2*self.generator.delta))./den;
			Gy = (self.generator.Ra + 0.5.*(self.generator.Xd1-self.generator.Xq).*sin(2*self.generator.delta))./den;
			Bx = (self.generator.Xd1.*cos(self.generator.delta).^2 + self.generator.Xq.*sin(self.generator.delta).^2)./den;
			By = (-self.generator.Xd1.*sin(self.generator.delta).^2 - self.generator.Xq.*cos(self.generator.delta).^2)./den;
			GxBy = Gx + By .* 1i;
			BxGy = Bx + Gy .* 1i;
		end

		%% solveDiffEquation: 求解微分代数方程
		function solveDiffEquation(self)
			omega1 = (self.generator.Pm - self.generator.Pe) ./ self.generator.TJ;
			self.generator.omega = self.generator.omega + omega1 .* self.dt;
			delta1 = self.generator.omegaBase .* (self.generator.omega - 1);
			self.generator.delta = self.generator.delta + delta1 .* self.dt;
		end

		%% solveNetEquation: 求解网络方程
		function solveNetEquation(self, stage)
			
			% 1. 忽略阻尼绕组, 求解发电机的虚拟注入电流
			[gxby, bxgy] = self.param_gb();
			[GxBy, BxGy] = self.param_GB();

			Ixy1 = bxgy .* self.generator.Eq1;
			% Ixy1 = (gxby .* self.generator.Ed1 + bxgy .* self.generator.Eq1);

			% 2. 求解节点电压方程
			I = zeros(length(self.ss.nodes.id), 1);
			index = getIndex(self.ss.nodes.id, self.generator.nid);	% 发电机节点索引
			I(index) = Ixy1;

			U = self.NAM.get() \ (I);

			% 3. 求解发电机的注入电流
			self.generator.Uxy = U(index);
			self.generator.Ixy = Ixy1 - GxBy .* real(self.generator.Uxy) - BxGy .* imag(self.generator.Uxy);

			% 4. 计算电磁功率
			% self.generator.Pe = real(self.generator.Uxy .* conj(self.generator.Ixy)) + abs(self.generator.Ixy).^2 .* self.generator.Ra;
			self.generator.Pe = (abs(self.generator.Ed1 + 1i.*self.generator.Eq1) .* abs(self.generator.Uxy))./self.generator.Xd1.*sin(self.generator.delta - angle(self.generator.Uxy));

			% debug
			self.vot = [self.vot, U];
			self.cur = [self.cur, I];
		end
		
		%% solveLargeDisturbance: 电力系统大干扰分析
		function solveLargeDisturbance(self, solver, mpcSteady)

			% 1. 初值计算
			self.setTransientInitialValue(mpcSteady);

			% 2. 网络方程的求解, 微分代数方程的数值积分, 系统状态改变时更新相应的节点导纳矩阵
			self.setTimeAndStep(solver.time, solver.dt);
			stagePrev = 0;
			for t = self.ct

				self.time = t;
				stage = self.getStage(t);
				% 系统状态发生了变化
				if stage ~= stagePrev
					self.setTransientData(stage, solver);
					self.setTransientNAM();
					stagePrev = stage;
				end

				self.solveNetEquation(stage);
				self.solveDiffEquation();
				self.setData();

			end
		end

	end

end
