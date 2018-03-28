%% SteadyState
classdef SteadyState < handle
	properties

	end
	methods

		%% testPF: 潮流程序测试
		function testPF(self)

			% 建立电力网稳态模型并初始化
			ss = Model.SteadyState();
			ss.init(case_ieee30());

			%% 设置求解器的基本信息
			solver.method = 'FDBX';	% 求解方法
			solver.maxIteration = 50;	% 最大迭代
			solver.epsilon = 1e-6;	% 收敛判据, 功率不平衡量标幺
			solver.start = 'defaul';	% 启动方式, default 为按发电机端电压起动
			% solver.checkReactivePower = true;	% 

			%% 求解
			result = ss.solvePowerFlow(solver);

			%% 建立视图对象并生成计算报告, 这里以文本文件作为输出结果
			viewModel = View.Plain();
			solver.documentName = 'report.txt';	% 文本计算报告
			viewModel.getPowerFlowReport(ss, solver, result);

			save('test.mat');	% 留作测试
		end

		%% testSC: 短路容量计算测试
		function testSC(self, config)

			mpc = getMpcSteady();

			steadyState = Model.SteadyState();
			steadyState.init(mpc);

			% result = steadyState.getShortCircultCapacity(nodes, generator, branches, 1);
			result = steadyState.getAllShortCircultCapacity();

			save('test.mat');
		end

	end
end