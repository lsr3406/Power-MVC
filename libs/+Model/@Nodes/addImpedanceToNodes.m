%% addImpedanceToNodes: 向系统中某两个节点之间添加一条支路, 并返回更新后的节点导纳矩阵.
function [output_nodeMatrix] = addImpedanceToNodes(nodes,nodeMatrix, from, to, impedance)
	
	nodeMatrix([from,to],[from,to]) = [
		1./impedance	-1./impedance;
		-1./impedance	1./impedance;
	] + nodeMatrix([from,to],[from,to]);
	output_nodeMatrix = nodeMatrix;

end
