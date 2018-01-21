%% getNewNodeParameter: 获取新的节点电压及相角,用于牛顿法
function [deltaTheta,deltaUpn,votage,angle] = getNewNodeParameter(nodes,deltaP,deltaQ,jacobianMatrix)

	% disp(size(jacobianMatrix))

	deltaP_t = deltaP([nodes.PQ_index;nodes.PV_index]);
	deltaQ_t = deltaQ([nodes.PQ_index]);

	deltaPara = jacobianMatrix\[deltaP_t;deltaQ_t];

	deltaTheta = deltaPara(1:(nodes.PQ_count+nodes.PV_count));	% PQ 及 PV 节点相角的修正量
	deltaUpn = deltaPara((nodes.PQ_count+nodes.PV_count+1):end);	% PQ 节点电压相对值的修正量
	
	% 下面这句用于将电压差值扩展至PV节点

	deltaUpn_t = zeros(nodes.PQ_count+nodes.PV_count+nodes.Balance_count,1);	% 初始化相角修正量
	deltaUpn_t([nodes.PQ_index']) = deltaUpn;	% 代入 PQ 及 PV 节点相角的修正量
	deltaUpn = deltaUpn_t;

	deltaTheta_t = zeros(nodes.PQ_count+nodes.PV_count+nodes.Balance_count,1);
	deltaTheta_t([nodes.PQ_index',nodes.PV_index']) = deltaTheta;
	deltaTheta = deltaTheta_t;


	votage = nodes.mag.*(1+deltaUpn);
	angle = nodes.ang + (deltaTheta);

end