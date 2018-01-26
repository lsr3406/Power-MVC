%% 定义电网中发电机的类
classdef Generator < handle

	properties

		id;				%　编号
		type;
		nid;			% 发电机所连接的节点编号
		Pg;				% 有功发电计划值(仅 pq pv 节点有效)
		Qg;				% 无功发电计划值(仅 pq 节点有效, 若 pv 节点转化为 pq 节点, 这个字段将被填充为 Qmax 的值)
		Pmax;			% 有功发电最大值
		Pmin;			% 无功发电最大值
		Qmax;			% 有功发电最小值
		Qmin;			% 无功发电最小值
		mBase;			% 基准容量
		status;			% 
		votage;			% 电压设定值(若此节点转为PQ节点则该属性会改变)

		% 暂态分析用
		xd;
		xd1;
		xq;
	end
	properties (Dependent)
		pqCount;
		pvCount;		% 调压发电机数量
		refCount;		% 调频发电机数量
	end

	properties (SetAccess = private)
        % nothing
    end

	methods

		%% Generator: 发电机构造方法
		function [self] = Generator(gen,nodes)
			self.nid = gen(:,1);
			self.id = (1:length(self.nid));
			self.mBase = gen(:,7);
			self.Pg = gen(:,2)./self.mBase;
			self.Qg = gen(:,3)./self.mBase;
			self.Qmax = gen(:,4)./self.mBase;
			self.Qmin = gen(:,5)./self.mBase;
			self.votage = gen(:,6);
			self.status = gen(:,8);
			self.Pmax = gen(:,9)./self.mBase;
			self.Pmin = gen(:,10)./self.mBase;
			self.status = ones(length(self.id), 1);

			self.type = zeros(length(self.id), 1);
			for k = self.id
				self.type(k) = nodes.type(self.nid(k) == nodes.id);
			end

			% self.xd = zeros(length(self.id), 1);
			self.xd1 = gen(:,end);
			% self.xq = zeros(length(self.id), 1);
		end	

		%% get.pqCount: 小型发电机数量
		function [pqCount] = get.pqCount(self)
			pqCount = sum(self.type == 1);
		end
		%% get.pvCount: 中型发电机数量
		function [pvCount] = get.pvCount(self)
			pvCount = sum(self.type == 2);
		end
		%% get.refCount: 大型发电机数量
		function [refCount] = get.refCount(self)
			refCount = sum(self.type == 3);
		end

	end
end