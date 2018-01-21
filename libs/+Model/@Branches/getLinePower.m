%% getLinePower: 计算普通线路的潮流
function [Sij,Sji,dS] = getLinePower(self,nodes,index)
	conjYi0 = (self.conductance(index) - self.susceptance(index).*i)./2;
	conjYj0 = (self.conductance(index) - self.susceptance(index).*i)./2;
	conjYij = conj(1./(self.resistance(index) + self.reactance(index).*i));
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