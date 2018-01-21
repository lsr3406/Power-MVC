%% checkReactivePower: 每迭代一次都检查一下无功功率是否充足
function [nodes_o] = checkReactivePower(nodes,generator,branches)

	nodes_o = [];
	% flag = 0;

	PV_Qmax = nodes.Qmax;
	PV_Qmin = nodes.Qmin;
	PV_Qrel = nodes.Qout+nodes.Qd;

	% 下面两句 返回的 nodesIndex 包括所有类型的节点
	nodesIndex_lack = find(PV_Qmax < PV_Qrel);
	nodesIndex_excess = find(PV_Qmin > PV_Qrel);

	% 只保留 PV 节点
	nodesIndex_lack(find(nodes.type(nodesIndex_lack)~=2)) = [];
	nodesIndex_excess(find(nodes.type(nodesIndex_excess)~=2)) = [];

	if isempty([nodesIndex_lack;nodesIndex_excess])
		return;
	end
	% flag = 1;

	nodes_o = nodes.convertToPQNode(generator,branches,nodesIndex_lack,nodesIndex_excess);

	% fprintf('%s', 'lack: ');
	% fprintf(' %d ',nodesIndex_lack);
	% fprintf('\n');
	% fprintf('%s', 'excess: ');
	% fprintf(' %d ',nodesIndex_excess);
	% fprintf('\n');

	for k = 1:length(nodesIndex_lack)
		fprintf('%s %d %s\n', '	  节点',nodes.id(nodesIndex_lack(k)),'无功功率不足,已转化为 PQ 节点');
	end
	for k = 1:length(nodesIndex_excess)
		fprintf('%s %d %s\n', '	  节点',nodes.id(nodesIndex_excess(k)),'无功功率过剩,已转化为 PQ 节点');
	end

end