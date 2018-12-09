% 节点导纳矩阵
% TODO 三角分解的设置
classdef NAM < handle
	
	properties
		value;
		real;
	end
	properties (Dependent)
		size;
	end

	methods

		%% get.size: 大小
		function [res] = get.size(self)
			res = size(self.value);
		end

		%% get: 返回节点导纳矩阵
		function [value] = get(self)
			value = self.value;
		end

		%% getReal: 返回节点导纳矩阵
		function [value] = getReal(self)
			value = self.real;
		end

		%% get: 设置节点导纳矩阵
		function set(self, nam)
			assert(nargin == 2 && ~isempty(nam))
			if isa(nam, 'Model.NAM')
				self.value = sparse(nam.value);
			else
				s = size(nam);
				assert(s(1) == s(2))
				self.value = nam;
			end
		end

		%% setReal: 设置节点导纳矩阵 2n 维实矩阵（内部转化）
		function setReal(self)
			n = self.size;
			even = 2.*(1:n);
			odd = even - 1;
			self.real = sparse(2*n);
			self.real(odd, odd) = real(self.value);
			self.real(odd, even) = -imag(self.value);
			self.real(even, odd) = imag(self.value);
			self.real(even, even) = real(self.value);
		end

		%% addAdmittance: 向系统中某节点添加接地导纳, 返回新的节点导纳矩阵
		% 允许 index 和 y 传向量
		function addAdmittance(self, index, y)
			self.value(index, index) = self.value(index, index) + diag(y);
		end

		%% addImpedance: 向系统中某两个节点之间添加一条支路, 并返回更新后的节点导纳矩阵.
		function addImpedance(self, from, to, z)
			self.value([from,to],[from,to]) = [
				1./z	-1./z;
				-1./z	1./z;
			] + self.value([from,to],[from,to]);
		end

		%% addLine: 向节点导纳矩阵中添加一个线路
		function addLine(self, from, to, z, y)
			self.value([from,to],[from,to]) = [
				1./z + y./2,	-1./z;
				-1./z,	1./z + y./2;
			] + self.value([from,to],[from,to]);
		end

		%% removeLine: 向节点导纳矩阵中移除一个线路
		function removeLine(self, from, to, z, y)
			self.value([from,to], [from,to]) = [
				-1./z - y./2,	1./z;
				1./z,	-1./z - y./2;
			] + self.value([from,to], [from,to]);
		end

		%% addTransformer: 向节点导纳矩阵中添加变压器
		function addTransformer(self, from, to, z, y, k)
			% k 可以是复数
			self.value([from, to], [from, to]) = [
				1./(abs(k).^2.*z) + y./2,	-1./(conj(k)*z);
				-1./(k*z),	1./z + y./2;
			] + self.value([from,to],[from,to]);
		end

		%% addBranch: 向节点导纳矩阵中添加线路或变压器
		function addBranch(self, from, to, z, y, k)
			% 这个方法将全部按照 matpower 对变压器的定义作处理
			if nargin < 5
				error('illegal branch parameter');
			end
			if nargin == 6 && k ~= 0 && k ~= 1
				self.value([from, to], [from, to]) = [
					1./(abs(k).^2.*z) + y./2,	-1./(conj(k)*z);
					-1./(k*z),	1./z + y./2;
				] + self.value([from,to],[from,to]);
			else
				self.value([from,to],[from,to]) = [
					1./z + y./2,	-1./z;
					-1./z,	1./z + y./2;
				] + self.value([from,to],[from,to]);
			end
		end

		%% init: 零矩阵
		function init(self, m, n)
			if nargin == 2
				self.value = sparse(m, m);
			elseif nargin == 3
				self.value = sparse(m, n);
			else
				error('illegal NAM size in NAM init');
			end
		end

		%% generate: 根据已知线路参数及节点参数计算节点导纳矩阵
		function generate(self, nodeData, lineData)

			% 任务: 根据已知信息求解得到系统的节点导纳矩阵.
			% 注: 本程序未考虑三绕组变压器

			% 线路参数, 包括线路和变压器, 直接来自稳态模型的属性
			lineDataSize = size(lineData);
			lineDataSize = lineDataSize(1);
			% lineData = reshape(lineData,7,lineDataSize)';

			lineFrom = lineData(:, 1);
			lineTo = lineData(:, 2);
			lineImpedance = lineData(:, 3) + lineData(:, 4).*1i;
			lineAdmittance = lineData(:, 5) + lineData(:, 6).*1i;
			ratio = lineData(:, 7).*exp(1i.*lineData(:, 8));

			% 节点参数, 主要是带有独立导纳设备的节点参数, 直接来自稳态模型的属性
			nodeData_size = size(nodeData);
			nodeData_size = nodeData_size(1);
			% nodeData = reshape(nodeData,3,nodeData_size)';

			nid = nodeData(:, 1);
			nodeAdmittance_dep = nodeData(:, 2) + nodeData(:, 3).*1i;

			% 假设系统中不存在变压器, 将线路阻抗和导纳填至节点导纳阵. 当系统中有变压器时, 将其转换成较精确的π形等效电路, 再求解节点导纳阵.
			% 节点导纳矩阵存储时, 是按照节点在原始数据的排列顺序存的

			nodeIndex = unique([lineTo, lineFrom]);
			self.value = sparse(length(nodeIndex),  length(nodeIndex));

			% 尽管线路在电压等级不太高时没有电晕损耗, 但考虑到变压器存在一定的有功和励磁损耗, 以及系统中某些节点会经阻抗接地, 在这里并不忽略节点的对地电导.
			% 根据线路信息解出系统所有元件的等效电路

			for k = 1:lineDataSize	% 遍历所有的线路
				fi = find(nid == lineFrom(k));
				ti = find(nid == lineTo(k));
				self.addBranch(fi, ti, lineImpedance(k), lineAdmittance(k), ratio(k));
			end

			% 遍历所有的节点
			self.addAdmittance(1:nodeData_size, nodeAdmittance_dep);
		end

		%% generateTrans: 生成一个秩为 n-1 的节点导纳矩阵, 去掉了所有的对地支路的影响
		function generateTrans(self, nodeData, lineData)

			lineDataSize = size(lineData);
			lineDataSize = lineDataSize(1);

			lineFrom = lineData(:, 1);
			lineTo = lineData(:, 2);
			lineImpedance = lineData(:, 3) + lineData(:, 4).*1i;
			lineAdmittance = lineData(:, 5) + lineData(:, 6).*1i;
			ratio = lineData(:, 7).*exp(i.*lineData(:, 8));

			nodeData_size = size(nodeData);
			nodeData_size = nodeData_size(1);

			nid = nodeData(:, 1);
			nodeAdmittance_dep = nodeData(:, 2) + nodeData(:, 3).*1i;

			nodeIndex = unique([lineTo, lineFrom]);
			self.value = sparse(length(nodeIndex),  length(nodeIndex));

			for k = 1:lineDataSize	% 遍历所有的线路
				fi = find(nid == lineFrom(k));
				ti = find(nid == lineTo(k));
				self.addBranch(fi, ti, lineImpedance(k), 0, ratio(k));
			end
		end

		%% generateShunt: 根据已知节点参数计算节点导纳矩阵
		function generateShunt(self, nodeData, lineData)

			lineDataSize = size(lineData);
			lineDataSize = lineDataSize(1);

			lineFrom = lineData(:, 1);
			lineTo = lineData(:, 2);
			lineImpedance = lineData(:, 3) + lineData(:, 4).*1i;
			lineAdmittance = lineData(:, 5) + lineData(:, 6).*1i;
			ratio = lineData(:, 7).*exp(i.*lineData(:, 8));

			nodeData_size = size(nodeData);
			nodeData_size = nodeData_size(1);

			nid = nodeData(:, 1);
			nodeAdmittance_dep = nodeData(:, 2) + nodeData(:, 3).*1i;

			nodeIndex = unique([lineTo, lineFrom]);
			self.value = sparse(length(nodeIndex),  length(nodeIndex));

			for k = 1:lineDataSize	% 遍历所有的线路
				fi = find(nid == lineFrom(k));
				ti = find(nid == lineTo(k));
				self.addBranch(fi, ti, inf, lineAdmittance(k), ratio(k));
			end

			self.addAdmittance(1:nodeData_size, nodeAdmittance_dep);
		end

	end

end