%% getTransformerPower: 计算变压器线路的损耗
function [Sij,Sji,dS] = getTransformerPower(self,nodes,index)
	[trFrom,trTo,trZ,trY1,trY2] = nodes.gamma2pi(self.fid(index),self.tid(index),(self.resistance(index)+i*self.reactance(index)),(self.conductance(index)+i.*self.susceptance(index)),self.ratio(index));
	% 下面这句应该去掉对励磁支路导纳的折半
	conjYi0 = (self.conductance(index) - self.susceptance(index).*i)./2 + conj(trY1);
	conjYj0 = conj(trY2);
	conjYij = conj(1./trZ);
	fid = find(nodes.id == self.fid(index));	% 得到线路始末端节点的id
	tid = find(nodes.id == self.tid(index));	% 得到线路始末端节点的id
	Sij = (nodes.mag(fid)).^2.*(conjYi0 + conjYij) - nodes.mag(fid).*nodes.mag(tid).*conjYij.*exp((nodes.ang(fid)-nodes.ang(tid)).*i);
	Sji = (nodes.mag(tid)).^2.*(conjYj0 + conjYij) - nodes.mag(fid).*nodes.mag(tid).*conjYij.*exp((nodes.ang(tid)-nodes.ang(fid)).*i);
	dS = Sij + Sji;
	% Pij = real(Sij);
	% Qij = imag(Sij);
	% Pji = real(Sji);
	% Qji = imag(Sji);
	% dP = real(dS);
	% dQ = imag(dS);
end