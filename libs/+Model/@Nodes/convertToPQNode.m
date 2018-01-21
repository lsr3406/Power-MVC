%% convertToPQNode: 解决某些节点无功不足或过剩的情况
function [nodes] = convertToPQNode(nodes,generator,branches,nodesIndex_lack,nodesIndex_excess)

	if ~isempty(nodesIndex_lack)

		% 处理无功不足的节点,首先将功率计划值设置为该节点能发出的最大无功功率
		% nodes.Qis(nodesIndex_lack') = nodes.Qis(nodesIndex_lack') + nodes.Qmax(nodesIndex_lack') - nodes.Qg(nodesIndex_lack');	% bug
		nodes.Qg(nodesIndex_lack') = nodes.Qmax(nodesIndex_lack');
		% generator.Qg(nodesIndex_lack') = generator.Qmax(nodesIndex_lack');
		nodes.Qis(nodesIndex_lack') = nodes.Qg(nodesIndex_lack') - nodes.Qd(nodesIndex_lack');

	end

	if ~isempty(nodesIndex_excess)

		% 处理无功过剩的节点,首先将功率计划值设置为该节点能发出的最小无功功率
		% nodes.Qis(nodesIndex_excess') = nodes.Qis(nodesIndex_excess') + nodes.Qmin(nodesIndex_excess') - nodes.Qg(nodesIndex_excess');	% bug
		nodes.Qg(nodesIndex_excess') = nodes.Qmin(nodesIndex_excess');
		% generator.Qg(nodesIndex_excess') = generator.Qmin(nodesIndex_excess');
		nodes.Qis(nodesIndex_excess') = nodes.Qg(nodesIndex_excess') - nodes.Qd(nodesIndex_excess');

	end

	% 修改节点类型,将其转化为 PQ 节点,并更新节点对象相应的的 PQ PV 节点的索引及数量属性

	nodes.type([nodesIndex_lack',nodesIndex_excess']) = 1;

	nodes.PQ_index = find(nodes.type == 1);
	nodes.PV_index = find(nodes.type == 2);
	nodes.PQ_count = length(nodes.PQ_index);
	nodes.PV_count = length(nodes.PV_index);

	% 最后修改PQ分解法中的系数矩阵,直接调用系数矩阵矩阵构造方法
	if ~isempty(nodes.PQM_B1)
		[nodes.PQM_B1,nodes.PQM_B2] = nodes.getPQM_matrix();
	end

end