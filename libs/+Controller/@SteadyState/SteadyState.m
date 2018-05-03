%% SteadyState
classdef SteadyState < handle
	properties

	end
	methods

		%% testPF: 潮流程序测试
		function testPF(self)

			% 建立电力网稳态模型并初始化
			obj = 'case9';
			ss = Model.SteadyState();
			ss.init(getMpcSteady(obj));

			%% 设置求解器的基本信息
			solver.method = 'FDBX';	% 求解方法
			solver.maxIteration = 50;	% 最大迭代
			solver.epsilon = 1e-6;	% 收敛判据, 功率不平衡量标幺
			solver.start = '';	% 启动方式, default 为按发电机端电压起动
			% solver.step = true;
			% solver.checkReactivePower = true;	% 

			%% 求解
			result = ss.solvePowerFlow(solver);
			%% 画图
			self.drawItlog(ss.itlog);

			%% 建立视图对象并生成计算报告, 这里以文本文件作为输出结果
			viewModel = View.Plain();
			solver.documentName = ['report_powerflow_', obj, '.txt'];	% 文本计算报告
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

		%% testIEEE: ieee 算例总结
		function testIEEE(self)

			mpcList = {'case5', 'case9', 'case14', 'case30', 'case39', 'case57', 'case118', 'case145', 'case300'};
			solverMethodList = {'NR', 'FD', 'FDBX', 'FDXB'};

			%% 设置求解器的基本信息
			solver.maxIteration = 100;	% 最大迭代
			solver.epsilon = 1e-5;	% 收敛判据, 功率不平衡量标幺
			solver.start = '';	% 启动方式, default 为按发电机端电压起动
			ss = Model.SteadyState();

			res = {};	% 存放测试结果
			for k1 = 1:length(mpcList)
				fprintf('%s', ['test for ', mpcList{k1}, '... ']);
				res{k1, 1} = mpcList{k1};

				%% 求解
				for k2 = 1:length(solverMethodList)
					fprintf('%s', [solverMethodList{k2}, '.. ']);

					% 建立电力网稳态模型并初始化
					eval(['ss.init(', mpcList{k1}, '());']);
					solver.method = solverMethodList{k2};
					result = ss.solvePowerFlow(solver);
					if(result.status == 1)	% 收敛
						res{k1, k2 + 1} = result.it;
					else
						res{k1, k2 + 1} = result.status;
					end
				end

				fprintf('%s\n', 'ok!');
			end
			disp(res);
		end

		%% testIT_EPS: 画出迭代次数与最大误差的曲线
		function testIT_EPS(self)

			mpcList = {'case9', 'case_ieee30'};
			titleList = {'IEEE-9', 'IEEE-30'};
			solverMethodList = {'NR', 'FD', 'FDBX', 'FDXB'};
			colorList = {[0.95, 0.26, 0.21], [0.3, 0.68, 0.31], [0.13, 0.59, 0.95], [1, 0.75, 0.03]};

			%% 设置求解器的基本信息
			solver.maxIteration = 7;	% 最大迭代
			solver.epsilon = 1e-20;	% 收敛判据, 功率不平衡量标幺
			solver.start = '';	% 启动方式, default 为按发电机端电压起动
			ss = Model.SteadyState();

			epsilon = {};	% 存放测试结果
			for k1 = 1:length(mpcList)
				fprintf('%s', ['test for ', mpcList{k1}, '... ']);
				epsilon{k1} = [];

				%% 求解
				for k2 = 1:length(solverMethodList)
					fprintf('%s', [solverMethodList{k2}, '.. ']);

					% 建立电力网稳态模型并初始化
					eval(['ss.init(', mpcList{k1}, '());']);
					solver.method = solverMethodList{k2};
					result = ss.solvePowerFlow(solver);
					epsilon{k1} = [epsilon{k1}; max(max(ss.itlog.dP), max(ss.itlog.dQ))];
				end
			end

			figure(1);
			for k1 = 1:length(epsilon)
				subplot(1, 2, k1);
				hold on;
				grid on;
				
				for k2 = 1:4
					cx = 0:0.1:solver.maxIteration;
					cy = interp1(0:solver.maxIteration, log10(epsilon{k1}(k2, :)), cx,'cubic');	% 计算插值函数在 0:0.1:solver.maxIteration 处的值，0:solver.maxIteration, log10(epsilon{k1}) 是观测值
					plot(cx, cy, 'lineWidth', 1.5, 'Color', colorList{k2});
				end
				legend('牛顿法','原生 PQ 分解法','BX 法','XB 法');
				for k2 = 1:4
					plot(0:solver.maxIteration, log10(epsilon{k1}(k2, :)), 'o', 'Color', colorList{k2});
				end

				set(get(gca,'XLabel'),'String','迭代次数');
				set(get(gca,'YLabel'),'String','最大功率误差');
				title(titleList{k1});
			end
		end

		%% testConvert: PV 节点转化测试
		function testConvert(self)

			% 建立电力网稳态模型
			ss = Model.SteadyState();

			%% 设置求解器的基本信息
			solver.method = 'NR';	% 求解方法
			solver.maxIteration = 100;	% 最大迭代
			solver.epsilon = 1e-6;	% 收敛判据, 功率不平衡量标幺
			solver.start = '';	% 启动方式, default 为按发电机端电压起动

			ss.init(case39());
			ss.solvePowerFlow(solver);
			res1 = ss.itlog;
			
			solver.checkReactivePower = true;	% 
			ss.init(case39());
			ss.solvePowerFlow(solver);
			res2 = ss.itlog;

			save('test.mat');	% 留作测试
		end

		%% testComp: 无功补偿测试
		function testComp(self)

			pfList = [0.8, 0.85, 0.9, 0.95];

			%% 设置求解器的基本信息
			solver.method = 'NR';	% 求解方法
			solver.maxIteration = 50;	% 最大迭代
			solver.epsilon = 1e-6;	% 收敛判据, 功率不平衡量标幺
			solver.start = '';	% 启动方式, default 为按发电机端电压起动

			% 建立电力网稳态模型
			ss = Model.SteadyState();
			mpc = case14();

			res1_mag = [];
			res1_ang = [];
			res2_mag = [];
			res2_ang = [];
			for k = 1:length(pfList)
				ss.init(mpc);
				ss.compensateReactivePowerByCapacitance(pfList(k));
				ss.solvePowerFlow(solver);
				res1_mag = [res1_mag, ss.nodes.mag];
				res1_ang = [res1_ang, ss.nodes.ang];

				ss.init(mpc);
				ss.compensateReactivePowerByCompensator(pfList(k));
				ss.solvePowerFlow(solver);
				res2_mag = [res2_mag, ss.nodes.mag];
				res2_ang = [res2_ang, ss.nodes.ang];
			end

			figure(1);
			subplot(1, 2, 1);
			hold on;
			grid on;
			plot(pfList, res1_mag', 'lineWidth', 1.5);
			set(get(gca,'XLabel'),'String','负荷功率因数');
			set(get(gca,'YLabel'),'String','节点电压');
			title('节点电压的变化（电容补偿）');
			legend('节点 1', '节点 2', '节点 3', '节点 4', '节点 5', '节点 6', '节点 7', '节点 8', '节点 9', '节点 10', '节点 11', '节点 12', '节点 13', '节点 14');

			subplot(1, 2, 2);
			hold on;
			grid on;
			plot(pfList, res1_ang'.*180./pi, 'lineWidth', 1.5);
			set(get(gca,'XLabel'),'String','负荷功率因数');
			set(get(gca,'YLabel'),'String','节点相角');
			title('节点相角的变化（电容补偿）');
			legend('节点 1', '节点 2', '节点 3', '节点 4', '节点 5', '节点 6', '节点 7', '节点 8', '节点 9', '节点 10', '节点 11', '节点 12', '节点 13', '节点 14');

			figure(2);
			subplot(1, 2, 1);
			hold on;
			grid on;
			plot(pfList, res2_mag', 'lineWidth', 1.5);
			set(get(gca,'XLabel'),'String','负荷功率因数');
			set(get(gca,'YLabel'),'String','节点电压');
			title('节点电压的变化（调相机补偿）');
			legend('节点 1', '节点 2', '节点 3', '节点 4', '节点 5', '节点 6', '节点 7', '节点 8', '节点 9', '节点 10', '节点 11', '节点 12', '节点 13', '节点 14');

			subplot(1, 2, 2);
			hold on;
			grid on;
			plot(pfList, res2_ang'.*180./pi, 'lineWidth', 1.5);
			set(get(gca,'XLabel'),'String','负荷功率因数');
			set(get(gca,'YLabel'),'String','节点相角');
			title('节点相角的变化（调相机补偿）');
			legend('节点 1', '节点 2', '节点 3', '节点 4', '节点 5', '节点 6', '节点 7', '节点 8', '节点 9', '节点 10', '节点 11', '节点 12', '节点 13', '节点 14');

			save('test.mat');
		end
		
		%% drawItlog: 画出迭代过程中的参数变化
		function drawItlog(self, itlog)
			s = size(itlog.mag);
			it = 0:(s(2) - 1);

			cit = 0:0.1:(s(2) - 1);
			mag =  interp1(it, itlog.mag', cit, 'pchip');	% 计算插值函数在 cit 处的值，it, itlog.mag 是观测值
			ang =  interp1(it, itlog.ang', cit, 'pchip');	% 计算插值函数在 cit 处的值，it, itlog.mag 是观测值
			dP =  interp1(it, abs(itlog.dP)', cit, 'pchip');	% 计算插值函数在 cit 处的值，it, itlog.mag 是观测值
			dQ =  interp1(it, abs(itlog.dQ)', cit, 'pchip');	% 计算插值函数在 cit 处的值，it, itlog.mag 是观测值

			figure();
			subplot(1, 2, 1);
			hold on;
			grid on;
			plot(cit, mag', 'LineWidth', 1.5);
			set(get(gca,'XLabel'),'String','迭代次数');
			set(get(gca,'YLabel'),'String','电压（pu）');
			title('节点电压的变化');
			axis([0, (s(2) - 1), 0.95, 1.05])	% 坐标范围	x:[0, (s(2) - 1)]
			legend('节点 1', '节点 2', '节点 3', '节点 4', '节点 5', '节点 6', '节点 7', '节点 8', '节点 9');
			subplot(1, 2, 2);
			hold on;
			grid on;
			plot(cit, ang'.*180./pi, 'LineWidth', 1.5);
			set(get(gca,'XLabel'),'String','迭代次数');
			set(get(gca,'YLabel'),'String','相角（deg）');
			title('节点相角的变化');
			legend('节点 1', '节点 2', '节点 3', '节点 4', '节点 5', '节点 6', '节点 7', '节点 8', '节点 9');
			axis([0, (s(2) - 1), -20, 20])	% 坐标范围	x:[0, (s(2) - 1)]

			figure();
			subplot(1, 2, 1);
			hold on;
			grid on;
			plot(cit, log10(dP'), 'LineWidth', 1.5);
			set(get(gca,'XLabel'),'String','迭代次数');
			set(get(gca,'YLabel'),'String','有功功率不平衡量（10^x pu）');
			title('有功功率不平衡量的变化');
			legend('节点 1', '节点 2', '节点 3', '节点 4', '节点 5', '节点 6', '节点 7', '节点 8', '节点 9');
			axis([0, (s(2) - 1), -9, 1])	% 坐标范围	x:[0, (s(2) - 1)] 	y:[-9, 1]
			subplot(1, 2, 2);
			hold on;
			grid on;
			plot(cit, log10(dQ'), 'LineWidth', 1.5);
			set(get(gca,'XLabel'),'String','迭代次数');
			set(get(gca,'YLabel'),'String','无功功率不平衡量（10^x）');
			title('无功功率不平衡量的变化');
			legend('节点 1', '节点 2', '节点 3', '节点 4', '节点 5', '节点 6', '节点 7', '节点 8', '节点 9');
			axis([0, (s(2) - 1), -9, 1])	% 坐标范围	x:[0, (s(2) - 1)] 	y:[-9, 1]
		end

		%% test_ly: 
		function test_ly(self)

			%% 设置求解器的基本信息
			solver.method = 'NR';	% 求解方法
			solver.maxIteration = 10;	% 最大迭代
			solver.epsilon = 1e-5;	% 收敛判据, 功率不平衡量标幺

			% 建立电力网稳态模型并初始化
			obj = 'case30';
			ss = Model.SteadyState();

			%% 设置参数（有功，电阻）
			% p 是与负荷的有功功率相关的系数
			param_p = 0.8:0.1:1.2;
			param_q = 0.8:0.1:1.2;

			% node id = 30
			U = zeros(length(param_p), length(param_q));
			theta = zeros(length(param_p), length(param_q));
			for k1 = 1:length(param_p)
				% disp(['param_p = ', num2str(param_p(k1))]);
				for k2 = 1:length(param_q)
					% 每次测试前先初始化
					ss.init(getMpcSteady(obj));
					% 调参数
					ss.nodes.Pd = ss.nodes.Pd .* param_p(k1);
					ss.generator.Pg = ss.generator.Pg .* param_p(k1);
					ss.nodes.Qd = ss.nodes.Qd .* param_q(k2);
					ss.generator.Qg = ss.generator.Qg .* param_q(k2);

					%% 求解
					result = ss.solvePowerFlow(solver);
					U(k1, k2) = ss.nodes.mag(30);
					theta(k1, k2) = ss.nodes.ang(30).*180./pi;
				end

			end
			% save('ly.mat');	% 留作测试
			figure()
			subplot(1, 2, 1);
				hold on;
				grid on;
				contourf(param_p, param_q, U);
				set(get(gca,'XLabel'),'String','p');
				set(get(gca,'YLabel'),'String','q');
				title('电压');
				% set(get(gca,'XTickLabel'),{'xTick1','xTick2','xTick3'});
				% legend('leg1','leg2','leg3');
				% gtext('gt1','gt2','gt3');
			subplot(1, 2, 2);
				hold on;
				grid on;
				contourf(param_p, param_q, theta);
				set(get(gca,'XLabel'),'String','p');
				set(get(gca,'YLabel'),'String','q');
				title('相角');
				% set(get(gca,'XTickLabel'),{'xTick1','xTick2','xTick3'});
				% legend('leg1','leg2','leg3');
				% gtext('gt1','gt2','gt3');

		end
		
		%% test_yq: 
		function test_yq(self)

			step = 0.1;
			x = 1:-step:-1;
			p = 64.^x;
			q = p;
			p = 2;
			q = 2;
			ploss = zeros(length(p));	% n*n matrix
			qloss = zeros(length(p));	% n*n matrix
			votage1 = zeros(length(p));
			votage2 = zeros(length(p));
			votage3 = zeros(length(p));
			votage4 = zeros(length(p));

			% 建立电力网稳态模型并初始化
			objMain = 'case39';
			objApp = 'case34';

			ssMain = Model.SteadyState();
			ssApp1 = Model.SteadyState();
			ssApp2 = Model.SteadyState();
			ssApp3 = Model.SteadyState();
			ssApp4 = Model.SteadyState();

			%% 设置求解器的基本信息
			solver.method = 'NR';	% 求解方法
			solver.maxIteration = 10;	% 最大迭代
			solver.epsilon = 1e-4;	% 收敛判据, 功率不平衡量标幺
			solver.start = '';	% 启动方式, default 为按发电机端电压起动

			epsilon = 1e-5;

			for k1 = 1:length(p)
				for k2 = 1:length(q)
					% fprintf('%s', ['p = ', num2str(p(k1)), '  q = ', num2str(q(k2)), '  ']);

					ssMain.init(getMpcSteady(objMain));
					ssApp1.init(getMpcSteady(objApp));
					ssApp2.init(getMpcSteady(objApp));
					ssApp3.init(getMpcSteady(objApp));
					ssApp4.init(getMpcSteady(objApp));
					% 设置配电网功率
					ssApp1.nodes.Pd = ssApp1.nodes.Pd .* p(k1);
					ssApp1.nodes.Qd = ssApp1.nodes.Qd .* q(k2);
					ssApp2.nodes.Pd = ssApp2.nodes.Pd .* p(k1);
					ssApp2.nodes.Qd = ssApp2.nodes.Qd .* q(k2);
					ssApp3.nodes.Pd = ssApp3.nodes.Pd .* p(k1);
					ssApp3.nodes.Qd = ssApp3.nodes.Qd .* q(k2);
					ssApp4.nodes.Pd = ssApp4.nodes.Pd .* p(k1);
					ssApp4.nodes.Qd = ssApp4.nodes.Qd .* q(k2);

					it = 0;
					dp = 0;
					dq = 0;
					dp_prev = inf;
					dq_prev = inf;
					while it < 20 && max(abs(dp - dp_prev), abs(dq - dq_prev)) > epsilon
						dp_prev = dp;
						dq_prev = dq;

						% 添加网损
						ssMain.nodes.Pd(40:43) = ssMain.nodes.Pd(40:43) + dp - dp_prev;
						ssMain.nodes.Qd(40:43) = ssMain.nodes.Qd(40:43) + dq - dq_prev;
						% 先计算主网，获取配网平衡节点电压
						resMain = ssMain.solvePowerFlow(solver);
						if resMain.status ~= 1
							dp = nan;dq = nan;break;
						end

						swingVotage = ssMain.nodes.mag(40:43)

						% 设置配网平衡节点电压
						ssApp1.nodes.mag(28) = swingVotage(1);
						ssApp2.nodes.mag(28) = swingVotage(2);
						ssApp3.nodes.mag(28) = swingVotage(3);
						ssApp4.nodes.mag(28) = swingVotage(4);

						% 计算配网潮流
						resApp1 = ssApp1.solvePowerFlow(solver);
						if resApp1.status ~= 1
							dp = nan;dq = nan;break;
						end
						resApp2 = ssApp2.solvePowerFlow(solver);
						if resApp2.status ~= 1
							dp = nan;dq = nan;break;
						end
						resApp3 = ssApp3.solvePowerFlow(solver);
						if resApp3.status ~= 1
							dp = nan;dq = nan;break;
						end
						resApp4 = ssApp4.solvePowerFlow(solver);
						if resApp4.status ~= 1
							dp = nan;dq = nan;break;
						end

						% 获取配网损耗
						dp1 = sum(ssApp1.branches.dP);
						dq1 = sum(ssApp1.branches.dQ);
						dp2 = sum(ssApp2.branches.dP);
						dq2 = sum(ssApp2.branches.dQ);
						dp3 = sum(ssApp3.branches.dP);
						dq3 = sum(ssApp3.branches.dQ);
						dp4 = sum(ssApp4.branches.dP);
						dq4 = sum(ssApp4.branches.dQ);

						dp = sum([dp1, dp2, dp3, dp4]);
						dq = sum([dq1, dq2, dq3, dq4]);
						it = it + 1;
					end
					% fprintf('%s', ['k = ', num2str(it), '  ']);
					% fprintf('\n');
					% 记录 dp, dq 
					ploss(k1, k2) = dp.*100;
					qloss(k1, k2) = dq.*100;
					votage1(k1, k2) = ssApp1.nodes.mag(28);
					votage2(k1, k2) = ssApp2.nodes.mag(28);
					votage3(k1, k2) = ssApp3.nodes.mag(28);
					votage4(k1, k2) = ssApp4.nodes.mag(28);

					% 找边界
					% if(~isnan(dp))
					% 	break;
					% end
				end
				save('yq.mat')
			end



			%% 求解
			
		end
	end
end






