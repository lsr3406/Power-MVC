%% 电力系统暂态分析控制器
classdef Fault < handle
	properties
		% nothing
	end

	methods
		%% testF: 大干扰稳定测试
		function testF(self)
		
			ft = Model.Fault();
			ft.init(getMpcFault());

			solver = [];

			ft.solveFault(solver, getMpcSteady());

			save('test');

		end
	end
end
