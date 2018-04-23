%% 电力系统暂态分析控制器
classdef TransientState < handle
	properties
		% nothing
	end

	methods
		%% testLD: 大干扰稳定测试
		function testLD(self)
		
			obj = 'case9';

			ts = Model.TransientState();
			ts.init(getMpcTransient(obj));

			solver.diffSolver = 'default';	% 改进的欧拉法
			solver.netSolver = 'default';		% 直接法
			solver.dt = 0.001;	% 步长
			solver.time = 2;

			ts.solveLargeDisturbance(solver, getMpcSteady(obj));
			self.drawWave(ts);

			save('test.mat');

		end

		%% drawWave: 画功角曲线
		function drawWave(self, ts)
			figure()
			subplot(231)
			hold on
			plot(ts.ct, (ts.delta(1, :)).*180./pi, 'lineWidth', 1.5)
			plot(ts.ct, (ts.delta(2, :)).*180./pi, 'lineWidth', 1.5)
			plot(ts.ct, (ts.delta(3, :)).*180./pi, 'lineWidth', 1.5)
			title('发电机功角')

			subplot(232)                 
			hold on
			plot(ts.ct, ts.itlog.Pe(1, :), 'b')
			plot(ts.ct, ts.itlog.Pe(2, :), 'r')
			plot(ts.ct, ts.itlog.Pe(3, :), 'y')
			plot(ts.ct, ts.itlog.Pm(1, :), 'b', 'lineWidth', 1.2)
			plot(ts.ct, ts.itlog.Pm(2, :), 'r', 'lineWidth', 1.2)
			plot(ts.ct, ts.itlog.Pm(3, :), 'y', 'lineWidth', 1.2)
			title('发电机有功')

			subplot(233)
			hold on
			plot(ts.ct, imag(ts.vot(1, :).*conj(ts.cur(1, :))))
			plot(ts.ct, imag(ts.vot(2, :).*conj(ts.cur(2, :))))
			plot(ts.ct, imag(ts.vot(3, :).*conj(ts.cur(3, :))))
			title('发电机无功')

			% plot(ts.ct, ts.omega, 'lineWidth', 1.5);

			subplot(234)                 
			hold on
			% plot(ts.ct, (ts.delta(1, :)).*180./pi, 'lineWidth', 1.5)
			plot(ts.ct, (ts.delta(2, :) - ts.delta(1, :)).*180./pi, 'lineWidth', 1.5)
			plot(ts.ct, (ts.delta(3, :) - ts.delta(1, :)).*180./pi, 'lineWidth', 1.5)
			title('发电机摇摆角')

			subplot(235)
			hold on
			for k = 1:9
				% plot(ts.ct, abs(ts.vot(k, :)))
				plot(real(ts.vot(k, :)), imag(ts.vot(k, :)), '-');
			end
			legend('node 1', 'node 2', 'node 3', 'node 4', 'node 5', 'node 6', 'node 7', 'node 8', 'node 9');
			title('电网电压')

			subplot(236)
			hold on
			for k = 1:9
				% plot(ts.ct, abs(ts.cur(k, :)))
				plot(real(ts.cur(k, :)), imag(ts.cur(k, :)), '-');
			end
			legend('node 1', 'node 2', 'node 3', 'node 4', 'node 5', 'node 6', 'node 7', 'node 8', 'node 9');
			title('电网电流')
		end
	end
end
