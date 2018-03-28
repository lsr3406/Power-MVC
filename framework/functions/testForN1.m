%% TODO delete
%% testForN1: 
function [result] = testForN1(ss)
	Sij = (ss.branches.Pij + ss.branches.Qij.*1i).*100.*1000;
	Sji = (ss.branches.Pji + ss.branches.Qji.*1i).*100.*1000;

	nfindex = zeros(length(ss.branches.id), 1);
	ntindex = zeros(length(ss.branches.id), 1);
	for k = 1:length(ss.branches.id)
		nfindex = find(ss.nodes.id == ss.branches.fid(k));
		ntindex = find(ss.nodes.id == ss.branches.tid(k));
	end
	
	Uf = ss.nodes.mag(nfindex).*220;
	Ut = ss.nodes.mag(ntindex).*220;
	disp([ss.nodes.id, ss.nodes.mag]);
	disp([ss.branches.id,ss.branches.fid,ss.branches.tid, abs(getLineCurrent(Sij, Uf)), abs(getLineCurrent(Sij, Ut))]);
end