%% getPQM_matrix: 获取 PQ 分解法中要用的第一及第二系数矩阵
function [PQM_B1,PQM_B2] = getPQM_matrix(nodes)

	PQM_B1 = imag(nodes.AdmittanceMatrix([nodes.PQ_index',nodes.PV_index'],[nodes.PQ_index',nodes.PV_index']));
	PQM_B2 = imag(nodes.AdmittanceMatrix([nodes.PQ_index'],[nodes.PQ_index']));

	% PQM_B1_inv = inv(PQM_B1);
	
	% PQM_B2_inv = inv(PQM_B2);

end
