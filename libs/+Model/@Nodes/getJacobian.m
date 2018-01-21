%% getJacobian: 在当前网络状态下获取雅可比矩阵,用于牛顿法的计算
function [jacobianMatrix] = getJacobian(nodes)

	% 为方便构造雅可比矩阵,现分别对PQ节点及PV节点提取出节点电压及节点导纳,并计算节点之间的相角矩阵

	AdmittanceMatrix_PQPV = nodes.AdmittanceMatrix([nodes.PQ_index',nodes.PV_index'],[nodes.PQ_index',nodes.PV_index']);
	AdmittanceMatrix_PQ = nodes.AdmittanceMatrix([nodes.PQ_index'],[nodes.PQ_index']);

	mag_PQPV = nodes.mag([nodes.PQ_index',nodes.PV_index']);
	ang_PQPV = nodes.ang([nodes.PQ_index',nodes.PV_index']);
	mag_PQ = nodes.mag([nodes.PQ_index']);
	ang_PQ = nodes.ang([nodes.PQ_index']);

	nodeAngle = ang_PQPV*ones(1,length(ang_PQPV));
	nodeAngle = nodeAngle - nodeAngle';

	% 开始计算雅可比矩阵,创建初始状态
	jac_primary = mag_PQPV*(mag_PQPV').*conj(AdmittanceMatrix_PQPV).*exp(i.*nodeAngle);
	% 这个结果的对角元无实际意义,下面去掉对角元
	jac_primary = jac_primary - diag(diag(jac_primary));

	% 先计算各个分块矩阵的非对角元
	% 需要注意的是,下面的两个索引表示在雅可比矩阵的初始状态中PQPV节点及PQ节点所在的位置,而 nodes.PQ_index 表示在节点对象中 PQ节点所在的位置
	PQPV_index_t = 1:(nodes.PQ_count+nodes.PV_count);
	PQ_index_t = 1:nodes.PQ_count;

	% 雅可比矩阵的非对角元均可以在初始状态矩阵中提取

	Hij = imag(jac_primary([PQPV_index_t],[PQPV_index_t]));	% Hij 表示系统中有功潮流差与相角差的数量关系
	Nij = real(jac_primary([PQPV_index_t],[PQ_index_t]));	% Nij 表示系统中有功潮流差与电压差的相对值的数量关系
	Jij = -real(jac_primary([PQ_index_t],[PQPV_index_t]));	% Jij 表示系统中无功潮流差与相角差的数量关系
	Lij = imag(jac_primary([PQ_index_t],[PQ_index_t]));	% Lij 表示系统中无功潮流差与电压差的相对值的数量关系

	% 再计算各个分块矩阵的对角元.求导的原因,使得对角元需要单独计算

	Hii_d = -nodes.Qout([nodes.PQ_index',nodes.PV_index']) - mag_PQPV.^2.*imag(diag(AdmittanceMatrix_PQPV));
	Nii_d = nodes.Pout(nodes.PQ_index') + mag_PQ.^2.*real(diag(AdmittanceMatrix_PQ));
	Jii_d = nodes.Pout(nodes.PQ_index') - mag_PQ.^2.*real(diag(AdmittanceMatrix_PQ));
	Lii_d = nodes.Qout(nodes.PQ_index') - mag_PQ.^2.*imag(diag(AdmittanceMatrix_PQ));

	Hij = Hij + diag(Hii_d);
	Nij(PQ_index_t,PQ_index_t) = Nij(PQ_index_t,PQ_index_t) + diag(Nii_d);
	Jij(PQ_index_t,PQ_index_t) = Jij(PQ_index_t,PQ_index_t) + diag(Jii_d);
	Lij = Lij + diag(Lii_d);

	% 在这里拼成雅可比矩阵
	jacobianMatrix = [
		Hij	Nij;
		Jij	Lij;
	];
end