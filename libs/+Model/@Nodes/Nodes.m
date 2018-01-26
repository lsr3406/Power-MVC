%% 定义电网中节点的类
classdef Nodes < handle

	properties

		id;				%　节点编号
		type;			% 节点类型	1PQ	2PV	3平衡
		g;	% 分流电导(MW)(或导纳标么值)
		b;	% 分流电纳(MVar)(或导纳标么值)
		mag0;			% 起始电压
		ang0;			% 起始相角
		mag;			% 电压
		ang;			% 相角
		Pis;			% 有功功率的计划值
		Qis;			% 无功功率的计划值
		Pd;				% 有功负荷
		Qd;				% 无功负荷
		Pc;				% 有功补偿量
		Qc;				% 无功补偿量

		Pout;			% 流出节点的有功功率
		Qout;			% 流出节点的无功功率
		Pg;				% 收敛时节点发电机发出的有功
		Qg;				% 收敛时节点发电机发出的无功

		Qmin;
		Qmax;

		area;			%　
		baseKV;		% 基准电压(kV)
		zone;			%
		Vmax;			% 最大电压(标幺值)
		Vmin;			% 最小电压(标幺值)
		AdmittanceMatrix;		% 节点导纳
		
		% PQM_B1_inv;
		% PQM_B2_inv;
		PQM_B1;
		PQM_B2;

		mBase;

		iterationData_mag;	% 存放每次迭代产生的结果
		iterationData_ang;	% 存放每次迭代产生的结果
		iterationData_dp;
		iterationData_dq;

	end

	properties (Dependent)
		pqIndex;		% 各类节点的索引
		pqCount;		% 各类节点的数量
		pvIndex;		% 各类节点的索引
		pvCount;		% 各类节点的数量
		refIndex;		% 各类节点的索引
		refCount;		% 各类节点的数量

		loadIndex;
		loadCount;
    end

	methods

		%% get.pqIndex: 获取 pq 节点的索引
		function [pqIndex] = get.pqIndex(self)
			pqIndex = find(self.type == 1);
		end
		%% get.pqCount: 获取 pq 节点的数量
		function [pqCount] = get.pqCount(self)
			pqCount = sum(self.type == 1);
		end
		%% get.pvIndex: 获取 pv 节点的索引
		function [pvIndex] = get.pvIndex(self)
			pvIndex = find(self.type == 2);
		end
		%% get.pvCount: 获取 pv 节点的数量
		function [pvCount] = get.pvCount(self)
			pvCount = sum(self.type == 2);
		end
		%% get.refIndex: 获取 ref 节点的索引
		function [refIndex] = get.refIndex(self)
			refIndex = find(self.type == 3);
		end
		%% get.refCount: 获取 ref 节点的数量
		function [refCount] = get.refCount(self)
			refCount = sum(self.type == 3);
		end
		%% get.loadIndex: 获取负荷节点的索引
		function [loadIndex] = get.loadIndex(self)
			loadIndex = find(self.Pd > 0);
		end
		%% get.loadCount: 获取负荷节点的数量
		function [loadCount] = get.loadCount(self)
			loadCount = sum(self.Pd > 0);
		end

		%% Nodes: 节点构造方法
		function [self] = Nodes(bus,baseMVA)
			self.mBase = baseMVA;
			% check();
			self.id = bus(:,1);
			self.type = bus(:,2);
			self.Pd = bus(:,3)./self.mBase;
			self.Qd = bus(:,4)./self.mBase;
			self.g = bus(:,5)./self.mBase;
			self.b = bus(:,6)./self.mBase;
			self.area = bus(:,7);
			self.baseKV = bus(:,10);
			self.mag0 = bus(:,8);
			self.ang0 = bus(:,9).*pi./180;
			self.zone = bus(:,11);
			% self.Vmax = bus(:,12);
			% self.Vmin = bus(:,13);

			% self.Qmin = zeros(length(self.id),1);
			% self.Qmax = zeros(length(self.id),1);
			if(self.refCount ~= 1)
				die('平衡节点数量不为 1');
			end
		end

		% 在方法的定义中需要将对象本身当做第一个参数传入,在调用时不需要传入
		%% getNodeData: 生成节点参数,主要是带有独立导纳设备的节点参数
		% 这里将节点上的并联电容视为恒阻抗模型,并将其归算至节点导纳矩阵
		function [nodeData] = getNodeData(self)

			nodeData = [self.id,self.g,self.b];
		end

		%% generateAdmittanceMatrix: 生成节点导纳矩阵
		% 调用时参数只传支路对象
		function [AM] = generateAdmittanceMatrix(self,branches)
			lineData = branches.getLineData();
			nodeData = self.getNodeData();
			AM = self.getMatrix(nodeData,lineData);
		end

		%% getIterationInitialValue: 获取迭代初始值(考虑对节点的设置及发电机的设置)
		function [self] = getIterationInitialValue(self,generator)

			if 0
				for k = generator.id
					index = find(self.id == generator.nid(k));
					self.mag0(index) = generator.votage(k);
				end
			end

			self.mag = self.mag0;
			self.ang = self.ang0;

			% 0-1 启动
			if 0
				self.mag = ones(length(self.id),1);
				self.ang = zeros(length(self.id),1);
			end
			self.iterationData_mag = [];
			self.iterationData_ang = [];
			self.iterationData_dq = [];
			self.iterationData_dp = [];
		end

		%% getPlan: 初始化,用于计算各节点功率计划值,并根据发电机的情况计算出该节点可发出的的最大无功功率
		function [self] = getPlan(self,generator)

			self.Pg = zeros(length(self.id),1);
			self.Qg = zeros(length(self.id),1);

			for k = generator.id
				index = find(self.id == generator.nid(k));

				% 这里对各节点的功率计划的赋值有一大部分(PQ节点)是没有意义的
				self.Pg(index) = self.Pg(index) + generator.Pg(k);
				self.Qg(index) = self.Qg(index) + generator.Qg(k);

				self.Qmin(index) = self.Qmin(index) + generator.Qmin(k);
				self.Qmax(index) = self.Qmax(index) + generator.Qmax(k);
			end
			% self.Pis = self.Pg - self.Pd;
			% self.Qis = self.Qg - self.Qd;
		end

		%% updatePlan: 更新各节点功率计划值. 4 月 7 日修复 bug
		function [self] = updatePlan(self)
			self.Pis = self.Pg - self.Pd;
			self.Qis = self.Qg - self.Qd;
		end

		%% getPowerOutflow: 计算从 id 节点注入电网的潮流
		function [Pf,Qf] = getPowerOutflow(self,id)
			Pf = zeros(length(id),1);
			Qf = zeros(length(id),1);
			for k = 1:length(id)
				index_t = find(self.id == id(k));
				Sf = conj(self.AdmittanceMatrix(index_t,:))*(self.mag.*exp(i.*(self.ang(index_t)-self.ang))).*self.mag(index_t);
				Pf(k) = real(Sf);
				Qf(k) = imag(Sf);
			end
		end

		%% getPowerUnbalance: 计算功率不平衡量,仅用于牛拉法和PQ分解法计算,返回变量均涵盖所有节点并已分类
		function [outflowP,outflowQ,deltaP,deltaQ] = getPowerUnbalance(self)

			[outflowP,outflowQ] = self.getPowerOutflow(self.id);
			deltaP = self.Pis - outflowP;
			deltaQ = self.Qis - outflowQ;
		end

		%% getConpensatorPower: 计算无功补偿器电量
		function [Pc,Qc] = getConpensatorPower(self)
			Pc = self.mag.^2.*self.g;
			Qc = self.mag.^2.*self.b;
		end

		%% getIterationData: 每次使用潮流方程之后调用此方法获取迭代信息
		function [self] = getIterationData(self,dp,dq)

			self.iterationData_mag = [self.iterationData_mag,self.mag];
			self.iterationData_ang = [self.iterationData_ang,self.ang];
			self.iterationData_dp = [self.iterationData_dp,dp];
			self.iterationData_dq = [self.iterationData_dq,dq];
		end

		%% changeLoadCapacity: 改变所有负荷的容量
		function [self] = changeLoadCapacity(self,rate)
			self.Pd = self.Pd.*rate;
			self.Qd = self.Qd.*rate;
		end

		%% changeLoadCapacity: 改变所有负荷有功
		function [self] = changeLoadActiveCapacity(self,rate)
			self.Pd = self.Pd.*rate;
			% self.Qd = self.Qd.*rate;
		end

		%% changeLoadFactorPerDeg: 改变所有负荷的功率因数
		function [self] = changeLoadFactorPerDeg(self,hysteresis_deg)
			PM = [
				cos(hysteresis_deg./180.*pi),-sin(hysteresis_deg./180.*pi);
				sin(hysteresis_deg./180.*pi),cos(hysteresis_deg./180.*pi);
			];
			[self.Pd,self.Qd] = PM*[self.Pd,self.Qd];
		end

		%% getLoadFactor: 计算负荷的功率因数,并返回负荷的性质(感性0,容性1)
		function [loadFactor,nature] = getLoadFactor(self)
			loadFactor = ones(length(self.id),1);
			nature = zeros(length(self.id),1);
			index = self.loadIndex;
			loadFactor(index) = self.Pd(index)./sqrt(self.Pd(index).^2+self.Qd(index).^2);
			nature(index) = (self.Qd(index) < 0);
		end

		%% compensateReactivePower: 对所有负荷使用调相机进行无功补偿,用于测试系统对负荷功率因数的反应
		function [self] = compensateReactivePowerByCompensator(self,LF_obj)
			tenfai = tan(acos(LF_obj));
			[LF_org,LF_nat] = self.getLoadFactor();
			% 不能使用 for k = self.loadIndex
			for k = self.loadIndex'
				if( (LF_org(k) < LF_obj) & (LF_nat(k) == 0) )
					self.Qd(k) = self.Pd(k).*tenfai;
				end
			end
		end

		%% compensateReactivePowerByCapacitance: 对所有负荷使用电容进行无功补偿,用于测试系统对负荷功率因数的反应
		function [self] = compensateReactivePowerByCapacitance(self,LF_obj)
			tenfai = tan(acos(LF_obj));
			[LF_org,LF_nat] = self.getLoadFactor();
			% 不能使用 for k = self.loadIndex
			for k = self.loadIndex'
				if( (LF_org(k) < LF_obj) & (LF_nat(k) == 0) )
					comp = self.Qd(k) - self.Pd(k).*tenfai;	% 待补偿量
					if self.b(k) < comp
						self.b(k) = comp;
					end
				end
			end
		end

	end
end

