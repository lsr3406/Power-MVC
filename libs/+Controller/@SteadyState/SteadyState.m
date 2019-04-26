%% SteadyState
classdef SteadyState < handle
    properties

    end
    methods

        %% testOPFBasic: OPF
        function testOPFBasic(self)

            ss = Model.SteadyState();
            ss.init('case5_test');

            solver.n_iters_max = 50;  % 最大迭代次数
            solver.epsilon = 1e-6;     % 对偶间隙需满足的精度
            solver.sigma = 0.1;        % 中心参数

            ss.solveOptimalPowerFlow(solver);
            disp(ss.opf.f);

            self.drawOPFLog(ss);

        end

        %% testPF: 普通潮流测试
        function testPF(self)

            ss = Model.SteadyState();
            ss.init(getMpcSteady('case118'));

            solver.method = 'NR';
            solver.n_iters_max = 50;
            solver.epsilon = 1e-6;
            solver.start = 'flat';

            res = ss.solvePowerFlow(solver);

            % disp([real(ss.bus.mag.*exp(ss.bus.ang.*1i)), imag(ss.bus.mag.*exp(ss.bus.ang.*1i)), ss.bus.Pg, ss.bus.Qg, ss.bus.Pd, ss.bus.Qd]);
            % disp([ss.branch.fid, ss.branch.tid, ss.branch.Pij*100, ss.branch.Qij*100, ss.branch.Pji*100, ss.branch.Qji*100, ss.branch.dP*100, ss.branch.dQ*100]);
            save('testPF.mat')
        end

        %% testIEEE: ieee 算例总结
        function testIEEE(self)

            mpcList = {'case5', 'case9', 'case14', 'case30', 'case39', 'case57', 'case118', 'case145', 'case300'};
            solverMethodList = {'NR', 'FD', 'FDBX', 'FDXB'};

            %% 设置求解器的基本信息
            solver.n_iters_max = 100;  % 最大迭代
            solver.epsilon = 1e-5;  % 收敛判据, 功率不平衡量标幺
            solver.start = '';  % 启动方式, default 为按发电机端电压起动
            ss = Model.SteadyState();

            res = {};   % 存放测试结果
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
                    if(result.status == 1)  % 收敛
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
            solver.n_iters_max = 7;    % 最大迭代
            solver.epsilon = 1e-20; % 收敛判据, 功率不平衡量标幺
            solver.start = '';  % 启动方式, default 为按发电机端电压起动
            ss = Model.SteadyState();

            epsilon = {};   % 存放测试结果
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
                    cx = 0:0.1:solver.n_iters_max;
                    cy = interp1(0:solver.n_iters_max, log10(epsilon{k1}(k2, :)), cx,'cubic'); % 计算插值函数在 0:0.1:solver.n_iters_max 处的值，0:solver.n_iters_max, log10(epsilon{k1}) 是观测值
                    plot(cx, cy, 'lineWidth', 1.5, 'Color', colorList{k2});
                end
                legend('牛顿法','原生 PQ 分解法','BX 法','XB 法');
                for k2 = 1:4
                    plot(0:solver.n_iters_max, log10(epsilon{k1}(k2, :)), 'o', 'Color', colorList{k2});
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
            solver.method = 'NR';   % 求解方法
            solver.n_iters_max = 100;  % 最大迭代
            solver.epsilon = 1e-6;  % 收敛判据, 功率不平衡量标幺
            solver.start = '';  % 启动方式, default 为按发电机端电压起动

            ss.init(case39());
            ss.solvePowerFlow(solver);
            res1 = ss.itlog;
            
            solver.checkReactivePower = true;   % 
            ss.init(case39());
            ss.solvePowerFlow(solver);
            res2 = ss.itlog;

            save('test.mat');   % 留作测试
        end

        %% testComp: 无功补偿测试
        function testComp(self)

            pfList = [0.8, 0.85, 0.9, 0.95];

            %% 设置求解器的基本信息
            solver.method = 'NR';   % 求解方法
            solver.n_iters_max = 50;   % 最大迭代
            solver.epsilon = 1e-6;  % 收敛判据, 功率不平衡量标幺
            solver.start = '';  % 启动方式, default 为按发电机端电压起动

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
                res1_mag = [res1_mag, ss.bus.mag];
                res1_ang = [res1_ang, ss.bus.ang];

                ss.init(mpc);
                ss.compensateReactivePowerByCompensator(pfList(k));
                ss.solvePowerFlow(solver);
                res2_mag = [res2_mag, ss.bus.mag];
                res2_ang = [res2_ang, ss.bus.ang];
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
            mag =  interp1(it, itlog.mag', cit, 'pchip');   % 计算插值函数在 cit 处的值，it, itlog.mag 是观测值
            ang =  interp1(it, itlog.ang', cit, 'pchip');   % 计算插值函数在 cit 处的值，it, itlog.mag 是观测值
            dP =  interp1(it, abs(itlog.dP)', cit, 'pchip');    % 计算插值函数在 cit 处的值，it, itlog.mag 是观测值
            dQ =  interp1(it, abs(itlog.dQ)', cit, 'pchip');    % 计算插值函数在 cit 处的值，it, itlog.mag 是观测值

            figure();
            subplot(1, 2, 1);
            hold on;
            grid on;
            plot(cit, mag', 'LineWidth', 1.5);
            set(get(gca,'XLabel'),'String','迭代次数');
            set(get(gca,'YLabel'),'String','电压（pu）');
            title('节点电压的变化');
            axis([0, (s(2) - 1), 0.95, 1.05])   % 坐标范围  x:[0, (s(2) - 1)]
            legend('节点 1', '节点 2', '节点 3', '节点 4', '节点 5', '节点 6', '节点 7', '节点 8', '节点 9');
            subplot(1, 2, 2);
            hold on;
            grid on;
            plot(cit, ang'.*180./pi, 'LineWidth', 1.5);
            set(get(gca,'XLabel'),'String','迭代次数');
            set(get(gca,'YLabel'),'String','相角（deg）');
            title('节点相角的变化');
            legend('节点 1', '节点 2', '节点 3', '节点 4', '节点 5', '节点 6', '节点 7', '节点 8', '节点 9');
            axis([0, (s(2) - 1), -20, 20])  % 坐标范围  x:[0, (s(2) - 1)]

            figure();
            subplot(1, 2, 1);
            hold on;
            grid on;
            plot(cit, log10(dP'), 'LineWidth', 1.5);
            set(get(gca,'XLabel'),'String','迭代次数');
            set(get(gca,'YLabel'),'String','有功功率不平衡量（10^x pu）');
            title('有功功率不平衡量的变化');
            legend('节点 1', '节点 2', '节点 3', '节点 4', '节点 5', '节点 6', '节点 7', '节点 8', '节点 9');
            axis([0, (s(2) - 1), -9, 1])    % 坐标范围  x:[0, (s(2) - 1)]   y:[-9, 1]
            subplot(1, 2, 2);
            hold on;
            grid on;
            plot(cit, log10(dQ'), 'LineWidth', 1.5);
            set(get(gca,'XLabel'),'String','迭代次数');
            set(get(gca,'YLabel'),'String','无功功率不平衡量（10^x）');
            title('无功功率不平衡量的变化');
            legend('节点 1', '节点 2', '节点 3', '节点 4', '节点 5', '节点 6', '节点 7', '节点 8', '节点 9');
            axis([0, (s(2) - 1), -9, 1])    % 坐标范围  x:[0, (s(2) - 1)]   y:[-9, 1]
        end

        %% drawOPFLog: 
        function drawOPFLog(self, ss)
            %% 线条颜色
            colorMap = {[0.96, 0.26, 0.21], [0.25, 0.32, 0.71], [0, 0.59, 0.53], [1, 0.92, 0.23], [0.47, 0.33, 0.28], [0.91, 0.12, 0.39], [0.13, 0.59, 0.95], [0.3, 0.69, 0.31], [1, 0.76, 0.03], [0.62, 0.62, 0.62], [0.61, 0.15, 0.69], [0.01, 0.66, 0.96], [0.55, 0.76, 0.29], [1, 0.6, 0], [0.38, 0.49, 0.55], [0.4, 0.23, 0.72], [0, 0.74, 0.83], [0.8, 0.86, 0.22], [1, 0.34, 0.13]};

            %% 插值初始化
            sz = size(ss.itlog.mag);
            iters = 0:(sz(2)-1);
            cx = 0:0.1:(sz(2)-1);

            figure();
            subplot(121);
            hold on;
            grid on;
            for k = 1:ss.n_bus
                plot(iters, ss.itlog.ang(k, :) .* 180 ./ pi, '.' , 'MarkerSize', 10, 'Color', colorMap{k});
                cy = interp1(iters, ss.itlog.ang(k, :)' .* 180 ./ pi, cx, 'pchip');
                plot(cx, cy, 'Color', colorMap{k});
            end
            title('节点相角在迭代过程中的变化');
            set(get(gca,'XLabel'),'String','iters');
            set(get(gca,'YLabel'),'String','\theta(deg)');

            % figure();
            subplot(122);
            hold on;
            grid on;
            for k = 1:ss.n_bus
                plot(iters, ss.itlog.mag(k, :), '.' , 'MarkerSize', 10, 'Color', colorMap{k});
                cy = interp1(iters, ss.itlog.mag(k, :)', cx, 'pchip');
                plot(cx, cy, 'Color', colorMap{k});
            end
            title('节点电压在迭代过程中的变化');
            set(get(gca,'XLabel'),'String','iters');
            set(get(gca,'YLabel'),'String','V(pu)');
            
            figure();
            subplot(121);
            hold on;
            grid on;
            cy = interp1(iters, ss.itlog.gap, cx, 'pchip');
            semilogy(iters, log10(ss.itlog.gap), '*', 'Color', colorMap{1});
            semilogy(cx, log10(cy), 'Color', colorMap{1});
            title('对偶间隙在迭代过程中的变化');
            set(get(gca,'XLabel'),'String','iters');
            set(get(gca,'YLabel'),'String','Gap');

            % figure();
            subplot(122);
            hold on;
            grid on;
            for k = 1:ss.n_branch
                plot(iters, ss.itlog.Pij(k, :) .* 100, '.' , 'MarkerSize', 10, 'Color', colorMap{k});
                cy = interp1(iters, ss.itlog.Pij(k, :)' .* 100, cx, 'pchip');
                plot(cx, cy, 'Color', colorMap{k});
            end
            title('线路有功在迭代过程中的变化(仅始端)');
            set(get(gca,'XLabel'),'String','iters');
            set(get(gca,'YLabel'),'String','P_{ij}(MW)');

            figure();
            subplot(121);
            hold on;
            grid on;
            for k = 1:ss.n_gen
                plot(iters, ss.itlog.Pg(k, :) .* 100, '.' , 'MarkerSize', 10, 'Color', colorMap{k});
                cy = interp1(iters, ss.itlog.Pg(k, :)' .* 100, cx, 'pchip');
                plot(cx, cy, 'Color', colorMap{k});
            end
            title('发电机有功在迭代过程中的变化');
            set(get(gca,'XLabel'),'String','iters');
            set(get(gca,'YLabel'),'String','P_G(MW)');
            
            % figure();
            subplot(122);
            hold on;
            grid on;
            for k = 1:ss.n_gen
                plot(iters, ss.itlog.Qg(k, :) .* 100, '.' , 'MarkerSize', 10, 'Color', colorMap{k});
                cy = interp1(iters, ss.itlog.Qg(k, :)' .* 100, cx, 'pchip');
                plot(cx, cy, 'Color', colorMap{k});
            end
            title('发电机无功在迭代过程中的变化');
            set(get(gca,'XLabel'),'String','iters');
            set(get(gca,'YLabel'),'String','Q_G(MVar)');
        end

    end
end






