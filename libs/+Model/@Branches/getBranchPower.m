%% getBranchPower: 计算线路功率及功率损耗
% TODO 修改写法
function [self] = getBranchPower(self,nodes)
	% BranchPower数据存储格式: P_f2t Q_f2t P_t2f Q_t2f dP dQ
	branchPower = zeros(length(self.id),6);
	for k = self.id
		if (self.ratio(k)==1)|(self.ratio(k)==0)
			[branchPower(k,1),branchPower(k,2),branchPower(k,3),branchPower(k,4),branchPower(k,5),branchPower(k,6)] = self.getLinePower(nodes,k);
		else
			[branchPower(k,1),branchPower(k,2),branchPower(k,3),branchPower(k,4),branchPower(k,5),branchPower(k,6)] = self.getTransformerPower(nodes,k);
		end	
	end
	self.Pij = real(branchPower(:,1));
	self.Qij = imag(branchPower(:,2));
	self.Pji = real(branchPower(:,3));
	self.Qji = imag(branchPower(:,4));
	self.dP = real(branchPower(:,5));
	self.dQ = imag(branchPower(:,6));
end