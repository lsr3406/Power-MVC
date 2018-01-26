%% testForN1: 
function [result] = testForN1(steadyState)
	Sij = (steadyState.branches.Pij + steadyState.branches.Qij.*1i).*100.*1000;
	Sji = (steadyState.branches.Pji + steadyState.branches.Qji.*1i).*100.*1000;

	nfindex = zeros(length(steadyState.branches.id), 1);
	ntindex = zeros(length(steadyState.branches.id), 1);
	for k = 1:length(steadyState.branches.id)
		nfindex = find(steadyState.nodes.id == steadyState.branches.fid(k));
		ntindex = find(steadyState.nodes.id == steadyState.branches.tid(k));
	end
	
	Uf = steadyState.nodes.mag(nfindex).*220;
	Ut = steadyState.nodes.mag(ntindex).*220;
	disp([steadyState.nodes.id, steadyState.nodes.mag]);
	disp([steadyState.branches.id,steadyState.branches.fid,steadyState.branches.tid, abs(getLineCurrent(Sij, Uf)), abs(getLineCurrent(Sij, Ut))]);
end