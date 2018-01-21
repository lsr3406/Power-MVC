%% SteadyState
classdef SteadyState < handle
	properties

	end
	methods

		%% getMpc: 获取电力系统数据
		function [mpc] = getMpc(name)
			% TODO getMpc()
			mpc = [];
		end


		

		%% test: TODO test description
		function test(self, config)

			mpc = case14();

			steadyState = Model.SteadyState();
			nodes = Model.Nodes(sortrows(mpc.bus),mpc.baseMVA);
			generator = Model.Generator(mpc.gen,nodes);
			branches = Model.Branches(mpc.branch);
			steadyState.init(mpc, nodes, generator, branches);

			solverConfig.method = 'NR';
			solverConfig.maxIteration = 50;
			solverConfig.epsilon = 1e-5;
			solverConfig.start = 'default';
			solverConfig.documentName = 'reportNR.txt';

			result = steadyState.solvePowerFlow(nodes, generator, branches, solverConfig);

			viewModel = View.Plain();
			viewModel.getPowerFlowReport(steadyState, nodes, generator, branches, solverConfig, result);

			save('test.mat');
			
		end

	end
end