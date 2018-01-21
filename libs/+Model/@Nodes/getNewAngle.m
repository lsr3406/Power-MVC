%% getNewAngle: 在 PQ 分解法中用来迭代一次得到一个新的相角
function [PQM_t1,newAngle] = getNewAngle(nodes,deltaP,PQM_B1)

	newAngle = nodes.ang;

	PQM_x1 = deltaP([nodes.PQ_index',nodes.PV_index'])./nodes.mag([nodes.PQ_index',nodes.PV_index']);
	
	PQM_t1 = PQM_B1\PQM_x1;
	
	newAngle([nodes.PQ_index',nodes.PV_index']) = newAngle([nodes.PQ_index',nodes.PV_index']) - PQM_t1./nodes.mag([nodes.PQ_index',nodes.PV_index']).*0.98;

end