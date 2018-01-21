%% addAdmittanceToNode: 向系统中某节点添加接地导纳, 返回新的节点导纳矩阵
function [outputs_nodeMatrix] = addAdmittanceToNode(nodeMatrix, nodeIndex, admittance)
	nodeMatrix(nodeIndex,nodeIndex) = nodeMatrix(nodeIndex,nodeIndex) + admittance;
	outputs_nodeMatrix = nodeMatrix;
end