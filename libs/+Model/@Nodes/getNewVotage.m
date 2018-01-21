%% getNewVotage: 在 PQ 分解法中用来迭代一次得到一个新的相角
function [PQM_t2,newVotage] = getNewVotage(nodes,deltaQ,PQM_B2)
	
	newVotage = nodes.mag;

	PQM_x2 = deltaQ(nodes.PQ_index')./nodes.mag(nodes.PQ_index');

	PQM_t2 = PQM_B2\PQM_x2;

	newVotage([nodes.PQ_index']) = newVotage([nodes.PQ_index']) - 0.98.*PQM_t2;

end