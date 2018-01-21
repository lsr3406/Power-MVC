%% getAdmittanceMatrix: 根据已知线路参数及节点参数计算节点导纳矩阵
function [AdmittanceMatrix] = getAdmittanceMatrix(nodes,steadyState,nodeData,lineData)

	% 任务: 根据已知信息求解得到系统的节点导纳矩阵.
	% 注: 本程序未考虑三绕组变压器

	% 线路参数, 包括线路和变压器
	lineData_size = size(lineData);
	lineData_size = lineData_size(1);
	% lineData = reshape(lineData,7,lineData_size)';

	lineIndex_start = lineData(:,1);
	lineIndex_end = lineData(:,2);
	lineImpedance = lineData(:,3) + lineData(:,4).*i;
	lineAdmittance = lineData(:,5) + lineData(:,6).*i;
	ratio = lineData(:,7).*exp(i.*lineData(:,8));

	% 节点参数, 主要是带有独立导纳设备的节点参数
	nodeData_size = size(nodeData);
	nodeData_size = nodeData_size(1);
	% nodeData = reshape(nodeData,3,nodeData_size)';

	nodeIndex_dep = nodeData(:,1);
	nodeAdmittance_dep = nodeData(:,2) + nodeData(:,3).*i;

% 假设系统中不存在变压器, 将线路阻抗和导纳填至节点导纳阵. 当系统中有变压器时, 将其转换成较精确的π形等效电路, 再求解节点导纳阵.

	nodeIndex = unique([lineIndex_end,lineIndex_start]);
	self.NAM = sparse(length(nodeIndex),length(nodeIndex));

	% 尽管线路在电压等级不太高时没有电晕损耗, 但考虑到变压器存在一定的有功和励磁损耗, 以及系统中某些节点会经阻抗接地, 在这里并不忽略节点的对地电导.
	% 根据线路信息解出系统所有元件的等效电路

	for k = 1:lineData_size	% 遍历所有的线路

		index_start_t = find(nodes.id == lineIndex_start(k));
		index_end_t = find(nodes.id == lineIndex_end(k));

		if (ratio(k) == 1 || ratio(k) == 0)	% 该线路是普通线路
			self.NAM = nodes.addImpedanceToNodes(self.NAM,index_start_t,index_end_t,lineImpedance(k));
			% 这里认为传入的线路导纳是整条线路上的导纳
			self.NAM = nodes.addAdmittanceToNode(self.NAM,index_start_t,lineAdmittance(k)./2);
			self.NAM = nodes.addAdmittanceToNode(self.NAM,index_end_t,lineAdmittance(k)./2);
		else	% 变压器
			% 首先将变压器转化成π形等效电路
			[trIndex_start, trIndex_end, trImpedance, trAdmittance1, trAdmittance2] = nodes.gamma2pi(index_start_t,index_end_t,lineImpedance(k),lineAdmittance(k),ratio(k));
			% 按照常规办法计算
			self.NAM = nodes.addImpedanceToNodes(self.NAM,trIndex_start,trIndex_end,trImpedance);
			% 下面这两行按照教材例题修改
			self.NAM = nodes.addAdmittanceToNode(self.NAM,trIndex_start,trAdmittance1);
			self.NAM = nodes.addAdmittanceToNode(self.NAM,trIndex_end,trAdmittance2);
		end
	end

	for k = 1:nodeData_size	% 遍历所有节点,计算导纳
		node_index_t = find(nodes.id == nodeIndex_dep(k));
		self.NAM = nodes.addAdmittanceToNode(self.NAM,node_index_t,nodeAdmittance_dep(k));
	end

end