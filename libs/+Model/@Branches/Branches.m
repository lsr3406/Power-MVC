%% 定义电网中线路及变压器的类
classdef Branches

	properties

		id;				%　节点编号
		fid;			% 起始节点
		tid;				% 终止节点
		r;		%　线路电阻
		x;		%　线路电抗
		g;	%　线路电导
		b;	%　线路电纳
		rateA			% 
		rateB			% 
		rateC			% 
		ratio;			% 变比(0 为线路)
		angle;			% 滞后
		status;			% 
		angmin;			% 
		angmax;			% 

		Pij;			% 线路有功潮流
		Pji;			% 线路有功潮流
		Qij;			% 线路无功潮流
		Qji;			% 线路无功潮流
		dP;				% 线损(有功)
		dQ;				% 线损(无功)
	end
	properties (Dependent)
		z;
		y;

		lineIndex;
		lineCount;
		transformerIndex;
		transformerCount;
	end

	methods
		%% get.z: 线路阻抗计算
		function [z] = get.z(self)
			z = self.r + self.x.*1i;
		end
		%% get.y: 线路导纳计算
		function [y] = get.y(self)
			y = self.g + self.b.*1i;
		end

		%% get.lineIndex: 线路索引
		function [lineIndex] = get.lineIndex(self)
			lineIndex = find(self.ratio==0);
		end
		%% get.lineCount: 线路数量
		function [lineCount] = get.lineCount(self)
			lineCount = sum(self.ratio==0);
		end
		%% get.transformerIndex: 变压器索引
		function [transformerIndex] = get.transformerIndex(self)
			transformerIndex = find(self.ratio~=0);
		end
		%% get.transformerCount: 变压器数量
		function [transformerCount] = get.transformerCount(self)
			transformerCount = sum(self.ratio~=0);
		end

		%% Branches: 支路构造方法
		function [self] = Branches(branch)
			self.fid = branch(:,1);
			self.tid = branch(:,2);
			self.id = (1:length(self.fid));
			self.r = branch(:,3);
			self.x = branch(:,4);
			self.g = zeros(length(self.id),1);	% 目前我们不支持添加线路对地导纳
			self.b = branch(:,5);
			self.rateA = branch(:,6);
			self.rateB = branch(:,7);
			self.rateC = branch(:,8);
			self.ratio = branch(:,9);
			self.angle = branch(:,10).*pi./360;
			self.status = branch(:,11);
			self.angmin = branch(:,12);
			self.angmax = branch(:,13);

			% self.Pij = zeros(1,length(self.id));
			% self.Pji = zeros(1,length(self.id));
			% self.Qij = zeros(1,length(self.id));
			% self.Qji = zeros(1,length(self.id));
			% self.dP = zeros(1,length(self.id));
			% self.dQ = zeros(1,length(self.id));
		end

		%% getLineData: 生成线路参数,包括线路和变压器
		function [lineData] = getLineData(self)
			% 起始节点	终止节点	线路电阻	线路电抗	线路对地电导	线路对地电纳	变比
			lineData = [self.fid,self.tid,self.r,self.x,self.g,self.b,self.ratio,self.angle];
		end

		%% changeLineParameter: 按比例改变线路参数
		function [self] = changeLineParameter(self,rr,rx,rb)
			if rr ~= 1 
				self.r = self.r.*rr;
			end
			if rx ~= 1 
				self.x = self.x.*rx;
			end
			if rb ~= 1 
				self.b = self.b.*rb;
			end
		end

		%% getRelatedNodesByBranche: 获取与支路相连的节点编号,id传列向量
		function [fid,tid] = getRelatedNodesByBranche(self,bid)
			index = find(self.id == bid);
			fid = self.from(index);
			tid = self.to(index);
		end

		%% getRelatedNodesByBranches: 获取与支路相连的节点编号,id传列向量
		% TODO 将代码缩减为 2 行
		function [fid,tid] = getRelatedNodesByBranches(self,bid)
			fid = zeros(length(bid),1);
			tid = zeros(length(bid),1);
			for k = 1:length(bid)
				[fid(k,1),tid(k,1)] = self.getRelatedNodesByBranche(bid(k));
			end
		end

		%% getRelatedBranchesByNode: 获取与节点相连的支路编号.branches_major是节点的下游线路
		% TODO 缩减代码为 2 行
		function [major,minor] = getRelatedBranchesByNode(self,nid)
			tfIndex = find(self.fid == nid);
			ttIndex = find(self.tid == nid);
			major = self.id(tfIndex);
			minor = self.id(ttIndex);
		end

		%% getRelatedBranchesByNodes: 获取与节点相连的支路编号.branches_major是节点的下游线路
		% TODO 缩减代码为 2 行
		function [major,minor] = getRelatedBranchesByNodes(self,nid)
			for k = 1:length(nid)
				[major{k,1},minor{k,1}] = self.getRelatedBranchesByNode(nid(k));
			end
		end

		%% getRelatedNodesByNode: 获取与节点直接相连的节点编号(列向量), id 传单个节点的 id, nodes_major是该节点的下游节点
		function [major,minor] = getRelatedNodesByNode(self,nid)
			[major,minor] = self.getRelatedBranchesByNode(nid);
			[mj1,major] = self.getRelatedNodesByBranches(major);
			[minor,mn2] = self.getRelatedNodesByBranches(minor);
			% minor = mj2;
			% major = mn1;
		end

		%% getRelatedNodesByNodes: 返回与节点直接相连的所有节点,major是该节点的下游节点
		% TODO 缩减代码为 2 行
		function [major,minor] = getRelatedNodesByNodes(self,nid)
			for k = 1:length(nid)
				[major{k,1},minor{k,1}] = self.getRelatedNodesByNode(nid(k));
			end
		end

		%% setRatio: 改变某些变压器变比
		function [self] = setRatio(self,bid,tap,to_flag)
			if nargin == 3
				to_flag = 0;
			end
			index = find(self.id == bid);
			ratio = self.ratio(index);
			if ratio == 0 || isempty(ratio)	% 不是变压器
				return;
			end
			if ~to_flag	% 改变一次侧变比
				self.ratio(index) = self.ratio(index)*(1 + tap);
			else	% 改变二次侧变比
				self.ratio(index) = self.ratio(index)/(1 + tap);
			end
		end
	end
end 

