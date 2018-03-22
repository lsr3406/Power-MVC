%% SteadyState
classdef SteadyState < handle
	properties

	end
	methods		

		%% testPF: 潮流程序测试
		function testPF(self)

			% 获取原始数据
			% mpc = getmpc();
			mpc = case89pegase();

			% 建立电力网稳态模型并初始化
			steadyState = Model.SteadyState();
			steadyState.init(mpc);

			%% 设置求解器的基本信息
			solverConfig.method = 'FD';	% 求解方法
			solverConfig.maxIteration = 20;	% 最大迭代
			solverConfig.epsilon = 1e-5;	% 收敛判据, 功率不平衡量标幺
			solverConfig.start = 'default';	% 启动方式, default 为按发电机端电压起动
			solverConfig.documentName = 'reportNR.txt';	% 文本计算报告

			%% 求解
			result = steadyState.solvePowerFlow(solverConfig);

			%% 建立视图对象并生成计算报告, 这里以文本文件作为输出结果
			viewModel = View.Plain();
			viewModel.getPowerFlowReport(steadyState, solverConfig, result);

			save('test.mat');	% 留作测试
		end

		%% testSC: 短路容量计算测试
		function testSC(self, config)

			mpc = getmpc();

			steadyState = Model.SteadyState();
			steadyState.init(mpc);

			% result = steadyState.getShortCircultCapacity(nodes, generator, branches, 1);
			result = steadyState.getAllShortCircultCapacity();

			save('test.mat');
		end

	end
end