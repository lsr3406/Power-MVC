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

			solver.dae = 'euler';	% 改进的欧拉法
			solver.net = 'default';		% 迭代法
			solver.dt = 0.001;	% 步长
			solver.time = 2;

			ts.solveLargeDisturbance(solver, getMpcSteady(obj));
			save('test.mat');
			self.logReport(ts);
			self.drawWave(ts);
		end

		%% testBasic: 39 节点系统基本测试
		function testBasic(self)

			obj = 'case39';

			ts = Model.TransientState();
			ts.init(getMpcTransient(obj));

			solver.dae = 'euler';	% 改进的欧拉法
			solver.net = 'default';		% 迭代法
			solver.dt = 0.001;	% 步长
			solver.time = 2;

			ts.solveLargeDisturbance(solver, getMpcSteady(obj));

			self.logReport(ts);
			self.drawWaveForPaper(ts);
			save('test.mat');
		end

		%% testAdvanced: 9 节点系统
		function testAdvanced(self)
			 
			obj = 'case9';

			ts = Model.TransientState();
			ts.init(getMpcTransient(obj));

			solver.dae = 'euler';	% 改进的欧拉法
			solver.net = 'default';		% 迭代法
			solver.dt = 0.001;	% 步长
			solver.time = 2;

			ts.solveLargeDisturbance(solver, getMpcSteady(obj));

			% self.logReport(ts);
			self.drawWaveForPaper(ts);
			
			index = [1:20:81, 101:50:2001];
			disp([ts.ct(index)', ts.delta(:, index)'*180/pi, (ts.delta(2, index)-ts.delta(1, index))'*180/pi]);
			save('test.mat');
		end

		%% drawWaveForPaper: 画功角曲线
		function drawWaveForPaper(self, ts)
			ts.ct = ts.ct(1:length(ts.delta(1, :)));
			genIndex = getIndex(ts.generator.nid, ts.ss.nodes.id);

			figure();
			subplot(1, 2, 1);
				hold on;
				for k = 1:length(ts.generator.nid)
					plot(ts.ct, (ts.delta(k, :)).*180./pi, 'lineWidth', 1.5);
				end
				set(get(gca,'XLabel'),'String','时间(s)');
				set(get(gca,'YLabel'),'String','δ(deg)');
				title('发电机功角');
			% subplot(1, 2, 2);
			% 	hold on;
			% 	for k = 1:length(ts.generator.nid)
			% 		plot(ts.ct, (ts.omega(k, :)), 'lineWidth', 1.5);
			% 	end
			% 	set(get(gca,'XLabel'),'String','时间(s)');
			% 	set(get(gca,'YLabel'),'String','ω');
			% 	title('发电机转速');
			subplot(1, 2, 2);
				hold on;
				res = nchoosek(genIndex, 2);
				res = res(1:ceil(length(res)./2), :);
				for k = 1:length(res(:, 1))
					plot(ts.ct, (ts.delta(res(k, 2), :) - ts.delta(res(k, 1), :)).*180./pi, 'lineWidth', 1.5);
				end
				set(get(gca,'XLabel'),'String','时间(s)');
				set(get(gca,'YLabel'),'String','Δδ(deg)');
				axis([0, 2, 0, 200])	% 坐标范围	x:[0, 2]	y:[ 0, 200]
				title('发电机摇摆角');
				legend('gen 2-1', 'gen 3-1');
		end

		%% logReport: 打印主要的信息
		function logReport(self, ts)
			maxSwingAngle = max(max(ts.delta) - min(ts.delta));
			maxSwingIndex = find(max(ts.delta) - min(ts.delta) == maxSwingAngle);
			generator1 = find(ts.delta(:, maxSwingIndex) == max(ts.delta(:, maxSwingIndex)));
			generator2 = find(ts.delta(:, maxSwingIndex) == min(ts.delta(:, maxSwingIndex)));
			fprintf('%s%6.2f%s\n', '    Solved time: ', ts.ct(end) - ts.ct(1), 's');
			fprintf('%s%6.4f%s\n', '           Step: ', ts.dt, 's');
			fprintf('%s%6.2f%s\n', 'Max Swing Angle: ', maxSwingAngle.*180./pi, '(deg)');
			fprintf('%s%6.4f%s\n', '           Time: ', ts.ct(maxSwingIndex), 's');
			fprintf('%s%d%s%d%s\n', '      Generator:  ', generator1, ' - ', generator2);
			fprintf('\n');
			self.logSwingMesBetween(ts, 0.0, 0.2);
			self.logSwingMesBetween(ts, 0.2, 0.4);
			self.logSwingMesBetween(ts, 0.4, 0.6);
			self.logSwingMesBetween(ts, 0.6, 0.8);
			self.logSwingMesBetween(ts, 0.8, 1.0);
			self.logSwingMesBetween(ts, 1.0, 1.2);
			self.logSwingMesBetween(ts, 1.2, 1.4);
			self.logSwingMesBetween(ts, 1.4, 1.6);
			self.logSwingMesBetween(ts, 1.6, 1.8);
			self.logSwingMesBetween(ts, 1.8, 2.0);
		end

		%% logSwingMesBetween: 获取某个区间内的摇摆信息
		function logSwingMesBetween(self, ts, l, r)
			l = (find(ts.ct >= l, 1));
			r = (find(ts.ct >= r, 1));
			delta = ts.delta(:, l:r);
			maxSwingAngle = max(max(delta) - min(delta));
			maxSwingIndex = find(max(delta) - min(delta) == maxSwingAngle);
			generator1 = find(delta(:, maxSwingIndex) == max(delta(:, maxSwingIndex)));
			generator2 = find(delta(:, maxSwingIndex) == min(delta(:, maxSwingIndex)));

			fprintf('%s%6.2f%s%6.4f%s%6.4f%s\n', 'Max Swing Angle: ', maxSwingAngle.*180./pi, '(deg) (in ', ts.ct(l), 's and ', ts.ct(r), 's)');
			fprintf('%s%6.4f%s%6.4f%s%6.4f%s\n', '           Time: ', ts.ct(l + maxSwingIndex - 1), 's (in ', ts.ct(l), 's and ', ts.ct(r), 's)');
			fprintf('%s%d%s%d%s\n', '      Generator:  ', generator1, ' - ', generator2);
			fprintf('\n');
		end

		%% drawWave: 画功角曲线
		function drawWave(self, ts)

			ts.ct = ts.ct(1:length(ts.delta(1, :)));
			genIndex = getIndex(ts.generator.nid, ts.ss.nodes.id);

			figure();
			subplot(231);
			hold on;
			for k = 1:length(ts.generator.nid)
				plot(ts.ct, (ts.delta(k, :)).*180./pi, 'lineWidth', 1.5);
			end
			title('发电机功角');

			subplot(232) ;
			hold on;
			for k = 1:length(ts.generator.nid)
				plot(ts.ct, ts.itlog.Pe(k, :), 'lineWidth', 1.5);
			end
			title('发电机有功');

			subplot(233);
			hold on;
			for k = genIndex
				plot(ts.ct, imag(ts.vot(k, :).*conj(ts.cur(k, :))));
			end
			title('发电机无功');

			subplot(234);
			hold on;
			res = nchoosek(genIndex, 2);
			for k = 1:length(res(:, 1))
				plot(ts.ct, (ts.delta(res(k, 2), :) - ts.delta(res(k, 1), :)).*180./pi, 'lineWidth', 1.5);
			end
			title('发电机摇摆角');

			subplot(235);
			hold on;
			for k = 1:length(ts.vot(:, 1))
				plot(real(ts.vot(k, :)), imag(ts.vot(k, :)), '-');
			end
			title('电网电压');

			subplot(236);
			hold on;
			for k = 1:length(ts.cur(:, 1))
				plot(real(ts.cur(k, :)), imag(ts.cur(k, :)), '-');
			end
			title('电网电流');
		end
	end
end
